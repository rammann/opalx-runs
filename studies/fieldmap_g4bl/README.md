# G4beamline field maps: the 2D cylinder and 3D grid formats

Does OPALX track correctly through a field map written in either of the two
G4beamline formats it can read? The maps here are written from fields given in
closed form, so the right answer is known before anything is tracked.

## Short answer

Yes, for the magnetic field, the electric field and both at once. **68 of 68 checks
pass.** The two formats hold the same field and track it identically, to the last digit. Tracking follows the closed form to 5e-6
relative for the quadrupole and 2.5e-5 for the uniform solenoid field, and what
is left is the time step resolving the map's hard ends, not the reading.

The sharp results:

| what | result |
|------|--------|
| the same field in both formats | every element of the 4x4 identical, and the final phase space identical |
| a posed `FIELDMAP` against a `SOLENOID` placed by `ELEMEDGE` | agree to 9e-14 |
| `SCALE = 2`, `normB/current = 2`, and doubling every value in the file | agree to 7e-11 |
| a comment block, shuffled rows and nine columns against the tidy file | identical |
| `ZREVERSE = TRUE` against a map written out already reversed | identical |
| uniform dipole against the arc | bend angle to 1.1e-6 rad, offset to 2.2 um |
| energy gained crossing a uniform Ez, against q E dz | 4.7e-9 relative |
| the same map with its electric half switched off | reproduces the magnet-only map exactly |
| quadrupole matrix against the closed form | 5.2e-6 relative, symplectic to 1.6e-9 |

## Run it

```bash
./run_all.sh              # write the maps and inputs, track all 14, test
./run_all.sh --test-only  # just re-run the analysis
```

Takes about five minutes: 18 runs at roughly 17 s each, one rank, no space
charge. `run_all.sh` finds the binary itself and refuses to run if it has no
`FIELDMAP` element, because a stale one gives 14 identical parse errors and no
hint why.

**The executable target is `opalx_exe`, not `opalx`.** `make opalx` builds only
the static library and leaves whatever executable was there before, which is
how this study first ran against a binary three days older than its sources.

## Files

| file | what it does |
|------|--------------|
| `make_inputs.py` | writes every map, `<case>/<case>.in`, `parts.txt`, the per-case README, and `cases.json`. The only place the geometry and field strengths are chosen |
| `fmlib.py` | writes the two map formats, holds the closed-form matrices, and reads the tracking output back |
| `run_tests.py` | the 14 tests, the tolerance table, `results.txt` |
| `plot_tests.py` | the six figures in `plots/` (not tracked; the repo ignores `plots/`) |
| `run_all.sh` | one command for the lot |
| `maps/` | the generated maps (not tracked; `run_all.sh` rewrites them) |
| `<case>/` | input, particles, per-case README, and the tracking output (not tracked) |

## The two formats

Both are read from a plain text file named by `FMAPFN`, in millimetres and
Tesla, and neither is normalised: the values are absolute, so the multiplier on
the element (`SCALE` on `FIELDMAP`, `KS` on `SOLENOID`) is a plain factor and 1
reproduces the map as written.

**3D cartesian, `grid`** — read by `src/Fields/G4BL3DGrid.cpp`. A
`grid` line giving the corner, the spacing and the point counts, then a `data`
section with one row per point carrying its own coordinates and field. Rows are
placed by the coordinates they carry, so their order in the file does not
matter. Six columns, or nine with `Ex Ey Ez` as well.

The electric columns are read: they are MV/m in the file and V/m once loaded. The two
fields are scaled by separate pairs of header keys, exactly as G4beamline does it --
`normB` and `current` for the magnetic field, `normE` and `gradient` for the electric one
-- and by separate attributes on the element, `SCALE` and `ESCALE`. Storage for the
electric field is only allocated once a non-zero value turns up, so the many nine-column
maps whose `Ex Ey Ez` are all zero cost nothing.

**2D axisymmetric, `cylinder`** — read by
`src/Fields/G4BL2DMagnetoStatic.cpp`. A `cylinder` line, then a `Bz` block and
a `Br` block, each one row per z holding all the r values, with column 0 on the
axis. A missing `Br` block leaves Br at zero, as in G4beamline.

Either may start with a `param` line. Only `normB` and `current` are used, and
the file's own pair is folded in on load as `normB/current`, so the element's
multiplier is exactly the `current=` written on the `fieldmap` line that places
the map in the G4beamline input.

`ZREVERSE` turns a map round -- mirror in z, negate the longitudinal component,
leave the radial one -- and works only on the cylinder format.

## Which fields, and why those

A constant field, and a field linear in each coordinate, are reproduced by
trilinear interpolation with **no discretisation error at all**, and the same
holds for bilinear interpolation of a constant field in the cylinder format. So
the map carries the analytic field exactly, and any disagreement is a reader,
placement or integration problem rather than a sampling one. (This is the
argument `unit_tests/AbsBeamline/TestFieldmapVsAnalytic.cpp` makes.)

- **Uniform By**, 0.058218 T over 1 m: a 0.1 GeV electron follows a circular arc
  of radius 5.7588 m and leaves at 10 degrees.
- **Quadrupole**, By = g x and Bx = g y with g = 0.083816 T/m, so
  k1 = -0.25 m^-2. This field is divergence free and curl free, so its transfer
  matrix is symplectic and that is checked.
- **Uniform Bz**, 0.1 T over 1 m, which turns the transverse momentum through
  -0.29827 rad. Note this is **not** a focusing solenoid: a box of uniform Bz
  has no radial field at its ends, so a particle entering parallel to the axis
  carries straight on. It is not a divergence-free field either, which is why
  symplecticity is not checked on it. It is used because both formats carry it
  exactly and its motion has a closed form:

      x  = x0 + S x'0 + D y'0          C = cos(kL),  S = sin(kL)/k
      x' =        C x'0 + s y'0        s = sin(kL),  D = (1 - cos(kL))/k
      y  = y0 - D x'0 + S y'0
      y' =      - s x'0 + C y'0

- **Uniform Ez**, 1 MV/m over 1 m. A particle crossing it gains q E dz of energy, and
  because the field has no z dependence that only depends on the z it crossed -- so the
  same number holds when a magnetic field is bending the orbit at the same time. For this
  electron it is -1.0000 MeV.
- **An asymmetric solenoid profile** for the `ZREVERSE` cases, with the radial
  field from the usual first-order expansion. Turning it round gives a
  different field, which is what makes the comparison worth anything; test 13
  checks that it does.

## Cases

Every case puts the map's own z = 0 at lab z = 1 m, so the field runs from
0.5 m to 1.5 m, and tracks 13 particles: a reference plus a +/- step in each
coordinate, which gives the transfer matrix by centred differences.

| case | format | element | what it is for |
|------|--------|---------|----------------|
| `grid_dipole` | grid | FIELDMAP | uniform By against the arc |
| `grid_dipole_shift` | grid | FIELDMAP | the same map posed 0.25 m downstream |
| `grid_quad` | grid | FIELDMAP | quadrupole against the closed form |
| `grid_quad_2x` | grid | FIELDMAP | every value in the file doubled; the scaling reference |
| `grid_quad_scale2` | grid | FIELDMAP | `SCALE = 2` on the element |
| `grid_quad_param` | grid | FIELDMAP | `normB = 4`, `current = 2` in the file |
| `grid_quad_messy` | grid | FIELDMAP | comments, shuffled rows, nine columns |
| `grid_sol` | grid | FIELDMAP | uniform Bz against the closed form |
| `grid_sol_dt2` | grid | FIELDMAP | the same with twice the time step |
| `cyl_sol` | cylinder | FIELDMAP | the same uniform Bz, the other format |
| `cyl_sol_solenoid` | cylinder | SOLENOID | the same map, placed by `ELEMEDGE` |
| `grid_efield` | grid | FIELDMAP | a purely electric map, against q E dz |
| `grid_efield_escale2` | grid | FIELDMAP | the same with `ESCALE = 2` |
| `grid_eb` | grid | FIELDMAP | electric and magnetic in one map, both on |
| `grid_eb_bonly` | grid | FIELDMAP | the same map with `ESCALE = 0` |
| `cyl_ramp` | cylinder | FIELDMAP | the asymmetric profile as written |
| `cyl_ramp_zrev` | cylinder | FIELDMAP | the same map with `ZREVERSE = TRUE` |
| `cyl_ramp_mirror` | cylinder | FIELDMAP | that field written out already reversed |

## The 14 tests

The tests are of two kinds, and it is worth telling them apart.

Checks **against a closed form** ask whether OPALX moves a particle the way the
field in the map says it should. They cannot be exact, because the map has hard
ends and the time step only resolves where those ends are to within one step.

Checks **between two runs** ask whether two ways of saying the same thing agree.
Both sides have the same hard ends and the same step, so those cancel and the
agreement is at round-off. These are the sharp ones.

| # | kind | what it checks |
|---|------|----------------|
| 1 | closed form | \|p\| does not move: a magnetic field does no work |
| 2 | closed form | the field starts and stops where the placement says |
| 3 | closed form | uniform dipole: bend angle and offset against the arc |
| 4 | two runs | posing the map 0.25 m downstream moves the field and nothing else |
| 5 | closed form | quadrupole 4x4 against the thick-quadrupole matrix |
| 6 | closed form | that matrix is symplectic |
| 7 | closed form | uniform Bz 4x4 against the closed form, in both formats |
| 8 | two runs | **the two formats holding the same field track the same** |
| 9 | two runs | posed `FIELDMAP` against `SOLENOID` placed by `ELEMEDGE` |
| 10 | two runs | the three ways of asking for twice the field |
| 11 | two runs | the awkward file against the tidy one |
| 12 | closed form | where the residual comes from: run it again with twice the step |
| 13 | two runs | `ZREVERSE` against a map written out already reversed |
| 14 | both | the electric field: energy against q E dz, `ESCALE`, and that `ESCALE = 0` leaves the magnet-only map behind |

## Results (`results.txt`)

68 of 68. The numbers worth quoting:

| test | quantity | measured | tolerance |
|------|----------|----------|-----------|
| 1 | \|p\| conserved, worst case | 2.2e-13 | 1e-12 |
| 3 | bend angle | 1.1e-6 rad | 1e-4 |
| 3 | offset at the exit face | 2.2e-6 m | 2e-4 |
| 5 | quadrupole 4x4, relative | 5.2e-6 | 5e-5 |
| 6 | symplecticity | 1.6e-9 | 1e-7 |
| 7 | uniform Bz 4x4, relative | 2.5e-5 | 5e-5 |
| 8 | grid against cylinder | 0 (exactly) | 1e-8 |
| 9 | `FIELDMAP` against `SOLENOID` | 9.3e-14 | 1e-11 |
| 10 | the three scalings | 7.0e-11 | 1e-8 |
| 11 | awkward file against tidy | 0 (exactly) | 1e-11 |
| 13 | `ZREVERSE` against pre-reversed | 0 (exactly) | 1e-8 |
| 14 | energy gained, against q E dz | 4.7e-9 relative | 2e-4 |
| 14 | the same with the orbit also bent | 1.7e-8 relative | 2e-4 |
| 14 | `ESCALE = 0` against the magnet-only map | 0 (exactly) | 1e-8 |

### About the residual against the closed form

Test 12 runs `grid_sol` again with twice the time step. The residual goes from
2.5e-5 to 3.4e-5, so it does follow the step -- but not as a clean factor of
two, because part of it is *where* the two hard ends fall inside a step rather
than how big the step is. That part also shifts the drifts stripped off when
the matrix is built. Allowing the field to be taken as running less than one
step beyond its nominal ends (-18 um and +14 um, against a 30 um step) takes
`grid_sol` from 2.5e-5 down to 4.6e-6. So the residual is the hard ends, not
the reading -- which is what test 8 says directly, by getting exactly zero
between the two formats.

## Tolerances

In the `TOL` table at the top of `run_tests.py`, with the reason for each
written next to it. The two that matter: `matrix` at 5e-5 relative is about ten
times the measured residual against the closed form, and `identical` /
`equivalent` at 1e-11 / 1e-8 are for runs that should be doing the same
arithmetic and in most cases do it bit for bit.

## Notes for whoever hits these next

- **`make opalx` does not build the executable.** The target is `opalx_exe`.
  The library target relinks happily and leaves a stale binary in `src/opalx`.
- **A `FIELDMAP` element takes no `L` and no `ELEMEDGE`.** It is placed by an
  absolute lab pose, its length comes from the map header, and OPALX allows
  only one placement convention per beamline -- so a line containing one must
  pose every element in it.
- **`ELEMEDGE` and the pose mean the same point**: the lab position of the
  map's own z = 0, which is what G4beamline's `place z=` refers to. OPALX adds
  the map's first z to `ELEMEDGE`. Test 9 is what pins this down.
- **A bent orbit is longer than its extent along z.** The 10 degree arc here is
  1.0051 m long over a 1 m field, so reading the orbit at path length 1.5 m
  gives an answer 0.9 mm out. Read it against `ref_z`. Along path length the
  straight section after a bend rises as sin(theta), not tan(theta).
- **`FMAPSCALE` does not exist.** The attribute is `SCALE` on `FIELDMAP` and
  `KS` on `SOLENOID`. `studies/g4bl_compare/asr61_dipole/asr61_dipole.in` still
  uses `FMAPSCALE` on an `SBEND`, and `SBEND` never reads `FMAPFN` on this
  branch either, so that input no longer works. Not touched here.
- **A quadrupole's field on the axis is zero**, so the reference particle sees
  nothing and the field window cannot be found from the `.stat` file. The
  window comes from `cases.json` instead. Same for the asymmetric solenoid
  profile, which has tapered to 2e-4 of its peak by the map ends.
- **The upper faces are excluded** by both readers (`r < end`, never `<=`),
  because interpolation needs a whole cell.
- **A map can carry an electric field, and it is scaled separately.** `SCALE` is the
  magnetic multiplier (G4beamline's `current=`) and `ESCALE` the electric one (its
  `gradient=`). One attribute will not do for both: G4beamline scales them
  independently, and the MUH2 separator is a map that needs that.
- **The load message says "3D static"**, not which fields are present. `getInfo()` runs
  from the element's `initialise()`, before the data is read, and nothing in the header
  says whether the electric columns are there.
- **What these maps do not cover**: `extendX/Y/Z` symmetry declarations, the `points`
  section and a cylinder `data` block are all rejected by the readers by design, with a
  named reason, and so are the 2D cylinder format's `Er` and `Ez` blocks -- so an electric
  field can only come in through the 3D grid format. Those rejections are covered by
  `unit_tests/Fields/TestG4BL2DMagnetoStatic.cpp` and `TestG4BL3DGrid.cpp`, not here. Note
  that `printfield type=cylinder` in G4beamline writes a `data` block, so its output
  cannot be fed to OPALX as it stands; `printfield type=grid` can be, as long as no
  `extend` line is emitted.
- **Checked against a real map too**, outside this study:
  `G4BL/runs/MuH2/separator_MUH2_400mm_320kV_MagneticFieldmapV18.3.BLFieldMap`, 991,440
  points of combined electric and magnetic field. Its on-axis `Ex` of 1 MV/m comes back
  through OPALX to 1.5e-6 MV/m.
