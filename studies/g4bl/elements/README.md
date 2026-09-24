# g4bl_compare — checking OPALX field maps against G4beamline

One folder per element. Each holds a G4beamline input, an OPALX input describing the
same element, and the same particles going through both. The two codes read the **same
`.g4blmap` file** — no conversion, no derived map.

The point is attribution. Comparing whole beamlines only says "roughly agrees";
comparing one map in isolation says which map, and where in it, is wrong.

**Scope.** These cases cover the five `grid` maps muE4 uses. They do **not** cover
`wsx_total.g4blmap`, the sixth map and the only `cylinder` one — the only file read by
`G4BL2DMagnetoStatic` rather than `G4BL3DGrid`. `make_cases.py` excludes it by design and
its only comparison is the older `wsx_solenoid/` notebooks, which run at the step sizes
this README shows are not converged and compare the field on the axis only, so that map's
off-axis `Br` has never been compared with anything. Read every statement below as being
about the five grid maps.

```
g4bl_compare/
  cmplib.py        readers for both codes, the six coordinates, tolerances, the result table
  make_cases.py    writes every case's decks and particle files from one description
  cases.json       generated; one record per case
  run_all.sh       generate, run both codes, test, plot
  run_tests.py     every check -> results.txt
  plot_tests.py    -> plots/
  results.txt      generated; ends with an N/N count
  plots/           generated; nine figures

  asr61_dipole/    ASR61_300d, the first bend's main map
  asr61_300sm/     ASR61_300sm, moved onto the axis so a beam can reach it
  asr62_dipole/    ASR62, the second and third bends' main map
  qsm600_quad/     QSM600, the only nine-column file
  asr61_group/     ASR61_300d + both ASR61_300sm placements, as muE4 places them
  asr62_d2_group/  ASR62 + both ASR62_sm placements of the bend at 7504 mm
  asr62_d3_group/  the same for the bend at 12032 mm

  wsx_solenoid/          WSX, 13 particles, notebook      (older, still valid)
  wsx_solenoid_gauss/    WSX, 15000 particles, notebook   (older, still valid)
```

## Running it

```bash
# The CMake target is opalx_exe. `make opalx` builds only the static library and
# leaves whatever executable was there before.
cd /Users/rammann/Code/OPALX/build && make -j8 opalx_exe

cd /Users/rammann/Code/OPALX/opalx-runs/studies/g4bl_compare
./run_all.sh                  # everything
./run_all.sh --pair-only      # skip the 20000-particle stage
./run_all.sh --test-only      # re-run the analysis without re-tracking
./run_all.sh asr61_dipole     # one case
cat results.txt
```

Python is the conda base environment, `/opt/homebrew/Caskroom/miniconda/base/bin/python`;
the system `python3` has no h5py or numpy. One MPI rank only — which row a particle lands
on across ranks is not guaranteed. Both codes resolve relative paths from the working
directory, so each case runs from inside its own folder.

## What the study measures

Three stages per case, in the order in which a failure makes the next one meaningless.

**The field, with no tracking.** Both codes are asked for the field at the same points,
OPALX through `DUMPEMFIELDS` and G4beamline through `fieldntuple`. Sampled at the map's
own grid points, which reads the table back, and halfway between them, which makes both
codes interpolate. If the interpolation schemes differed, the second would fail while the
first passed.

**Nineteen particles.** The reference, then a small ± step in each of the six coordinates,
then a large ± step in x, y and dp/p. The small steps give a transfer matrix from each code
by centred differences; the large ones probe where the map stops being linear, which the
small steps cannot see. The set is exactly symmetric in ±, so its mean is exactly zero:
OPALX under `FROMFILE` takes its reference orbit from the bunch mean, and that is what
makes the mean coincide with G4beamline's declared reference particle.

**Twenty thousand Gaussian muons.** Matched by id, so the comparison is still per particle
and not only a comparison of distribution widths. A map that is read correctly near the
axis and wrongly further out shows up here as a difference that grows with amplitude, and
nowhere else.

## What it found

**The two codes agree on the field to the last digit G4beamline prints.** In all seven
cases the worst difference is exactly one printing step — 1e-7 T on the 0.086 T ASR61 map,
1e-6 T on the 0.215 T quadrupole — at the map's own grid points and halfway between them
alike. So both the table and the trilinear interpolation match.

**Tracking agrees to between 3.9e-8 and 9.4e-6 m** at the exit plane, over orbits that bend
through 41 degrees and reach 1.27 m off axis, and to between 1.2e-7 and 8.2e-6 rad in angle.

**OPALX conserves |p| far better than G4beamline can report it.** A static magnetic field
does no work, and OPALX holds |p| between 1.7e-14 and 6.2e-11 across the seven cases.
G4beamline's own figure is 1.5e-6 to 2.2e-6, which is its six-figure ASCII output rather
than its tracking.

**The transfer matrices agree to within the precision the comparison has.** This is now a
real check, and it was not at first. `row_floors()` originally took each row's floor from
the worst measured OPALX-minus-G4beamline disagreement in that coordinate; since the matrix
difference is built from those same differences, the ratio was bounded by about one by
construction and, against a tolerance of 3, **could not fail**. Replacing the OPALX exit
state with zeros, with noise, or with x doubled all still passed, so the "0.49 to 0.89 of
the floor" it printed carried no information. The floor is now derived from G4beamline's
printed precision alone, per coordinate, with no reference to the OPALX data: the real data
passes at 0.49 and those same corruptions fail by two to three orders of magnitude.

**With 20000 particles, every one is recorded at every plane in both codes**, and the rms
of the bunch agrees to the digit G4beamline prints. Nothing is lost, so no number here is
biased by a missing tail.

**The difference does not grow with amplitude.** This is the result the large bunch exists
to produce. For the ASR61 dipole the median per-particle difference at the exit is flat at
about 8 µm from 0 to 40 mm of starting radius — the map is read as well 40 mm off axis as
on it. The one case where it does grow is `qsm600_quad`, from 3e-9 m on axis to 3e-7 m at
30 mm, which is what a quadrupole must do: its field is proportional to displacement, so a
particle further out sees more field and accumulates more of everything, including the
difference. Figure `07_gauss_per_particle.png`.

Two things the 20000-particle stage measures that 19 particles cannot, and both are about
where a monitor samples rather than about the maps:

- **At a plane inside the field the angle disagrees by up to 2.4e-4 rad** while the
  position still agrees to 2e-5 m. Both codes interpolate the position onto the plane, but
  the momentum is reported from whichever tracking step the particle was on, and inside a
  field the momentum is changing with z. At 0.086 T and Bρ = 0.0934 T·m that angle is
  **0.26 mm of longitudinal sampling offset** — one to two OPALX steps. Outside the field
  there is nothing to change and the same planes agree to 2e-8 rad. `run_tests.py` reports
  the in-field planes separately and converts the angle into the length that explains it;
  for the quadrupole it cannot, because a quadrupole is zero on its own axis and the
  `.stat` file only carries the reference orbit.
- **The largest of 20000 draws is not comparable to the largest of 19.** The per-particle
  difference has a distribution: median 8.0e-6 rad, worst 3.7e-5. The 19-particle stage's
  worst is 8.2e-6, i.e. the median of the same distribution. Halving the time step moves
  the worst only from 3.7e-5 to 3.1e-5, so it is the tail and not the step. The
  20000-particle stage is therefore tested on its median and 99th percentile, which do not
  depend on how many particles were drawn, and the maximum is reported rather than tested.

**Adding a map of zeros changes nothing, measured rather than assumed.** `asr62_d3_group`
is `asr62_dipole` plus two placements of `asr62shim_280_sm_track.g4blmap`, and its exit
state comes out identical — exactly in G4beamline, to 1.5e-11 m in OPALX.

### The two cases at the real muE4 angle

`asr61_bisector` and `asr61_group_bisector` place the ASR61 magnet the way muE4 does:
between two `cornerarc` commands that turn the centreline 20 degrees before it and 20
degrees after, so the beam enters the map at 20 degrees to the map's own axis instead of
straight down it. The positions and rotations come from `mue4lib.walk()`, the same code
the whole-line study uses, and they reproduce that study's numbers exactly --
`X = 0.069900004`, `Z = 2.974500020`, `THETA = 0.349065850`.

Two things these cases do that the straight ones cannot.

**They check that OPALX puts the magnet where G4beamline puts it.** G4beamline says
`cornerarc`; OPALX needs an absolute position and rotation. Every straight case uses a
rotation of 0 or 180 degrees, so the conversion between the two was never tested. Asking
both codes for the field at the same 10074 points in the lab gives the same answer to
1e-7 T, with no tracking involved. If the conversion were wrong this would fail at once.

**They see about a hundred times more finely.** The beam stays within 35 mm of the
centreline instead of swinging out to 1.27 m, and G4beamline prints six digits, so the
smallest visible difference drops from 1e-5 m to 1e-7 m.

That extra resolution immediately showed that `DT = 1e-12` s and `maxStep = 0.1` mm are
**not** converged after all. They only looked converged because the straight cases could
not see below 1e-5 m. So each of these cases is run twice, at one step and at half it, and
`run_tests.py` reports how much the answer moved. Where the error falls in proportion to
the step -- halving the step halves it -- the step's contribution is removed and what is
left is tested. For `asr61_bisector` that is **0.028 µm** in position and **0.069 µrad** in
angle, against a 0.1 µm printing floor: the two codes agree as closely as this comparison
can measure.

**`asr61_group_bisector` has not converged, and this is not hidden.** It adds the two
`ASR61_300sm` placements, and the beam does reach them -- it gets to x = 433 mm in the
map's own frame, past the main map's 390 mm edge, and those files move the exit by 7.7 mm
in both codes. Crossing from one field file into another is where an integrator converges
slowest, and it shows:

| step | worst position | worst angle | matrix, in printing floors |
|---|---|---|---|
| 1e-12 s / 0.1 mm | 2.28 µm | 0.97 µrad | 26.2 |
| 5e-13 s / 0.05 mm | 1.20 µm | 0.58 µrad | 7.4 |
| 2.5e-13 s / 0.025 mm | 0.75 µm | 0.40 µrad | 5.4 |

Each halving still reduces the difference, so most of what is left is step size, but the
reduction is slowing and the curve has not flattened. **Whether any of it is a real
difference between the codes is not yet answered.** Four checks on this case fail as a
result, and they are left failing rather than given a looser limit: the honest statement is
that this geometry needs a smaller step than is affordable here, not that the codes agree.

The same applies to its 20000-particle stage, which holds the largest disagreement anywhere
in the study -- the mean x at the last plane differs by 36 µm on a 23 mm wide beam, and the
99th percentile of the per-particle angle difference is 91 µrad. Those are at the finest
step the stage runs at, and they have not been shown to be converged either.

### Both codes had to be converged first, and neither was

This is the part worth knowing before trusting any number above. At G4beamline's
`maxStep = 1` mm and OPALX's `DT = 1e-11` s, the settings this study's ancestors used, the
ASR61 dipole comes out **49 µm apart**. Almost none of that is the map:

| OPALX `DT` [s] | worst exit \|dx\| against G4beamline, by `maxStep` [mm] |
|---|---|
| | **1.0** — **0.5** — **0.25** — **0.1** |
| 1e-11 | 48.8 — 72.0 — 62.0 — 62.0 µm |
| 5e-12 | 46.3 — 38.1 — 38.1 — 38.1 µm |
| 2e-12 | 24.4 — 12.1 — 12.0 — 12.0 µm |
| 1e-12 | 31.7 — 12.1 — 11.4 — **7.5** µm |

OPALX's own answer moves 84, then 34, then 14 µm as the step halves, so at 1e-11 s it is
still some 30 µm from its own limit. G4beamline at 1 mm is 30 µm from its limit too. Once
both are converged the difference falls to 7.5 µm, which is G4beamline's ASCII resolution
on a 1.27 m coordinate. The decks therefore use `DT = 1e-12` s and `maxStep = 0.1` mm for
the 19-particle stage, and one step coarser for the 20000-particle stage, where what is
left is 12 µm against a bunch 10 mm across.

**A comparison can sit inside its tolerance and measure nothing.** Every tolerance in
`cmplib.TOL` is written against its measurement floor, and two of them cannot be constants:

- The **field** tolerance is three printing steps *of that case's own peak field*. A fixed
  2e-7 T passes the dipoles and fails the quadrupole while both agree perfectly, because
  six significant figures resolve 0.086 T to 1e-7 T but 0.215 T only to 1e-6 T.
- The **transfer matrix** tolerance is three times *that row's own* floor, measured as the
  worst disagreement in that coordinate divided by the difference step. The six rows do not
  share a floor: after a 41 degree bend the `delta` row's is twenty times the `x` row's.
  A single number for the whole matrix reports the `x` row's floor and then fails the
  `delta` row with it.

A relative test needs the same care. M[ζ ← x] of the ASR61 group is 3.6% of the largest
element in the matrix, so dividing a difference that *is* the floor by it reports 3.4% and
looks like a real disagreement. The relative check is therefore secondary, and only runs on
elements at least a hundred times the floor.

### Things about the maps worth recording

- **`asr62shim_280_sm_track.g4blmap` is identically zero** — all 52726 rows, confirmed by
  reading the file. It is placed four times in muE4 and contributes nothing, so
  `asr62_d2_group` and `asr62_d3_group` add no field over `asr62_dipole`. What they do test
  is that OPALX handles `rotation=Y180` and `rotation=Y180,Z180` as a rotation of the
  element — grid maps do not support `ZREVERSE`, so it has to be expressed that way — and
  `asr62_d3_group`'s exit state does come out equal to `asr62_dipole`'s.
- **`ASR61_300sm` never sees a beam on the axis.** Its box is x = 380…630 mm; in muE4 the
  orbit reaches it only after the main map has bent the beam outward. Standalone its
  position is shifted so the box straddles the axis. Both codes get the identical shift, so
  the comparison stays one to one, and its behaviour at the real muE4 position is covered
  by `asr61_group`.
- **Both codes drop the top face of a map box** when the map is placed without rotation,
  although the file holds values there — the ASR61 map has 0.033 T at x = 390 mm. The
  OPALX readers exclude it by `r < end`, never `<=`, because interpolating needs a whole
  cell above the point, and G4beamline does the same. Where a map is placed with a 180
  degree rotation they do **not** agree on the face: `sin(pi)` is 1.22e-16 rather than 0,
  so the sample lands a few 1e-16 m either side of it and each code falls on a different
  side. `asr61_group` measures 0.0026 T there and `asr62_d2_group` 0.105 T. The `edge`
  sample set reports what each code returned; it is never tested, because a face has no
  thickness and no particle spends path length on it.
- **A sample point a few 1e-14 m inside a face is resolved differently by the two codes**:
  OPALX interpolates and G4beamline returns zero. That is a knife edge 0.04 picometres
  wide, so no particle can land on it, but a sample grid can — a scan written as
  `-190 + 2*290` lands on 389.99999999999996 and reads as a 0.05 T disagreement that is
  entirely the last bit of a float. `run_tests.py` reports points on a box face separately
  instead of testing them, and `three_inside()` keeps the scans a whole step clear of both
  ends. A face is found by checking the point is on the boundary **and** inside the box in
  the other two axes; testing the face plane alone flags every sample in the ASR62 groups,
  because `ASR62_sm` has a boundary at y = 0 and the sample plane is y = 0.

### Things about the two codes worth recording

- **OPALX writes `DUMPEMFIELDS` output into `data/`**, not the working directory.
- `DUMPEMFIELDS` does see a `FIELDMAP` element. It runs from `ParallelTracker.cpp` just
  before the tracking loop, after the maps are read, and `studies/mue4_g4bl/README.md`
  listed this as untried.
- **A `FIELDMAP` element takes neither `ELEMEDGE` nor `L`**, and OPALX allows one placement
  convention per beamline, so every element in these decks is placed by absolute position
  and rotation. OPALX then sorts the line alphabetically by element name rather than by
  position, which is why the names carry a number in front.
- **Monitors inside the field do record.** Four of the six planes in the bend cases sit
  inside the map and all of them fill.
- **`SCALE` on `FIELDMAP` is G4beamline's deck-side `current=`, verbatim.** The G4beamline
  readers store absolute Tesla and do not normalise.
- **`CHARGE = 1` on the OPALX `BEAM` is mandatory.** `ParticleProperties` maps `MUON` to
  −1 and `Beam::execute` only consults that table when the attribute is absent. A charge
  sign error looks exactly like a map turned the wrong way round.
- **G4beamline's parameter is `maxStep`, capital S.** `maxstep` is a different, unused one.
- **`fieldmap` places no physical volume**, so without the vacuum `WORLDBOX` every track
  leaves the world immediately and the ntuples come out empty without complaint.
- **`start` must come before every `place`**, or g4bl aborts with `Invalid start`.
- **In `run_all.sh`, `grep -c` not `grep -q`.** Under `set -o pipefail`, `grep -q` exits on
  the first match, `strings` then dies of SIGPIPE, and the pipeline reports failure on a
  binary that is perfectly fine.
- macOS ships **bash 3.2**, which has no `mapfile` and errors on expanding an empty array
  under `set -u`.

## The older solenoid cases

`wsx_solenoid/` and `wsx_solenoid_gauss/` predate this harness and still work. They cover
the one `cylinder` map in muE4 and are driven by notebooks rather than `run_all.sh`; their
`ELEMEDGE` placement is still valid because a `SOLENOID` — unlike a `FIELDMAP` — accepts
it. `asr61_dipole/` used to be a third case of that kind, built on `SBEND` with `FMAPSCALE`;
neither that attribute nor `SBEND` + `FMAPFN` exists on this branch, so it has been rebuilt
as a `FIELDMAP` placed by position. Its recorded result, exit |dx| ≤ 50 µm, is reproduced
exactly, and it now records all 19 particles where the old case dropped two.
