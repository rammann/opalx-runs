# g4bl_compare — checking OPALX elements against G4beamline

One folder per element. Each holds a G4BL input, an OPALX input describing the same
element, and a notebook that pushes the **same particles** through both and compares
the outputs particle by particle.

The point is attribution. Comparing whole beamlines only says "roughly agrees";
comparing one element in isolation says which element, and which term of its map, is
wrong.

```
g4bl_compare/
  wsx_solenoid/              13 test particles -> per-particle diff + transfer matrix
    wsx_solenoid.g4bl
    wsx_solenoid.in
    compare.ipynb
  wsx_solenoid_gauss/        15000 Gaussian particles -> monitor plots, 8 planes
    wsx_solenoid_gauss.g4bl
    wsx_solenoid_gauss.in
    compare_gauss.ipynb
  toy_transfer_matrix.py     how the finite-difference transfer matrix works
```

Two ways of looking at the same solenoid. `wsx_solenoid/` is the strict one: few
particles, compared one by one, plus a 6x6 transfer matrix. `wsx_solenoid_gauss/` is
the visual one: a real bunch, phase-space plots at 8 planes from both codes side by
side, with sigma_x != sigma_y so the Larmor rotation is directly visible as the x-y
ellipse turning.

`parts.txt` (OPALX) and `beam.txt` (G4BL) are written by the first cell of the
notebook from one array, so the two codes cannot be fed different particles. They are
gitignored for the same reason run output is.

## Running a case

```bash
cd wsx_solenoid

# 1. first two cells of compare.ipynb -> parts.txt, beam.txt

# 2. G4BL
source /Users/rammann/Code/G4BL/G4beamline-3.08.app/Contents/root/bin/thisroot.sh
g4bl wsx_solenoid.g4bl > g4bl.log 2>&1

# 3. OPALX  (one rank only -- particle-to-row identity across ranks is not guaranteed)
mpirun -n 1 /Users/rammann/Code/OPALX/build/src/opalx wsx_solenoid.in --info 1 > run.log 2>&1

# 4. rest of compare.ipynb
```

The notebook has a cell that runs both for you if the binaries are in the usual places
(`build/src/opalx`, `opalx/build_serial/src/opalx`, `opalx/build/src/opalx`).

Python is the conda base env — `/opt/homebrew/Caskroom/miniconda/base/bin/python`. The
system `python3` has no h5py or numpy.

Both codes must run **from inside the case folder**: field map and particle file paths
are relative to the working directory.

## How a case is built

**Same particles.** The reference plus a ± step in each of the 6 phase-space
coordinates — 13 particles. Small enough to compare one by one instead of
statistically, and the ± pairs give a 6×6 transfer matrix from each code by centered
finite differences.

The set must stay **exactly symmetric in ±**. OPALX under `FROMFILE` takes its
reference orbit from the bunch *mean*, and that symmetry is what makes the mean
coincide with G4BL's `reference particle`. One asymmetric particle puts a silent
offset into every transverse comparison.

**Same frame.** The element is tested in isolation, so its position is free. Pick it so
`z` in the G4BL input [mm] / 1000 equals path length `s` in the OPALX input [m],
exactly. Nothing to fit, no offset in the notebook.

For WSX: the map's own z runs −1500…+1500 mm, so `place MAG z=1800` puts the field at
z = 300…3300 mm. OPALX places a map as `z_lab = ELEMEDGE + z_map`, so **`ELEMEDGE` is
the lab position of the map centre, not of the field start** — hence `ELEMEDGE = 1.800`.
Recording planes at 150 and 3450 mm, both in field-free drift.

**Same field.** The OPALX input points `FMAPFN` straight at the G4BL `.g4blmap`. No
conversion, no derived map files. OPALX reads the G4BL `cylinder` format natively
(`G4BL2DMagnetoStatic`) and does not normalise it, so `KS = 1.0` is G4BL's `current=1.`.
`ZREVERSE = TRUE` on the solenoid reproduces G4BL's `place ... rotation=Y180` (see
below); without it the on-axis peak lands at s = 2.150 m with the wrong sign instead of
at s = 1.450 m.

**Nothing else can differ.** Vacuum world, no iron, no decay, no stochastic processes,
no space charge, no scraping.

## What the notebooks report

Both open with a **setup diagram**: where the solenoid field sits, where the recording
planes are, the map centre (`ELEMEDGE`), the creation point and `ZSTOP`, over a plot of
the on-axis Bz. It renders before either code has run — the field panel just says so
until `g4bl_axis.txt` exists.

### wsx_solenoid/compare.ipynb

- per-particle exit `x, y, x', y', ζ, δ` from both codes and the difference
- the same at the entrance plane — pure plumbing, and it should be zero
- the 6×6 transfer matrix from each code, side by side, plus the element-wise
  difference and the symplectic residual
- the same three matrices as annotated heatmaps — G4BL and OPALX on a shared colour
  scale, the difference on its own — and a second panel showing the difference
  *relative to each element's own size*, so a small element that is badly wrong isn't
  hidden by a large element that is slightly off
- the Larmor rotation angle, **with sign**, three ways: from each matrix, from the x+
  particle's trajectory, and from ∫Bz ds / (2Bρ)
- OPALX `Bx_ref, By_ref, Bz_ref` from the `.stat` file against G4BL's on-axis
  `fieldntuple` dump, plotted per component. No tracking is involved, so a mismatch
  there is a field-import problem and everything below it is meaningless until it
  passes. `Bx` and `By` vanish on the axis of an axisymmetric magnet, so those two are
  an alignment check — nonzero means the reference orbit is off the magnet axis, or
  the map isn't axisymmetric about it.

`x'` and `y'` are the primary comparison: they're momentum ratios, so they don't depend
on either code's frame convention.

## Adding an element

Copy `wsx_solenoid/`, and change three things: the G4BL element (`fieldmap`/`genericquad`/…
plus its `place`), the OPALX element line, and the plane positions if the element has a
different length. The notebook is element-agnostic below the first two cells.

## Notes for whoever hits these next

- **`rotation=Y180`** in G4BL is `R = diag(−1, 1, −1)`: it mirrors the map in z *and*
  flips the sign of Bz (`Bz_lab(z) = −Bz_map(z_place − z)`, Br unchanged). The WSX map
  is strongly asymmetric in z — peak at map-local +350 mm — so this is not a
  symmetry you can ignore. Confirmed here by G4BL's own field dump: with
  `place MAG z=1800` the on-axis peak comes out at **z = 1450 mm** (= 1800 − 350) and
  `∫Bz dz = −0.259104 T·m`, i.e. mirrored and negative. OPALX reproduces it with
  `ZREVERSE = TRUE` on the solenoid, which mirrors the map in z and negates Bz at load
  time. `KS = −1` is *not* the same thing: it would negate Br too. With `ZREVERSE` the
  OPALX `Bz_ref` trace gives peak −0.260442 T at s = 1.4499 m and `∫Bz dz = −0.259104
  T·m`, matching G4BL to the `%.6g` output floor.
  `opalx-runs/data/wsx_solenoid_1D.map`, produced by
  `runs/mue4_fieldmaps/tools/g4blmap_solenoid_to_astra.py`, applies neither: its peak
  sits 700 mm the other side of centre and `∫Bz ds` comes out **positive**. The
  Larmor-angle sign is the sharpest test of it.
- **The Larmor angle has to be read off the matrix carefully.** `M[0:2,0:2] = cosθ·A`
  and `M[2:4,0:2] = −sinθ·A` for a common 2×2 `A`, so `atan2(−M[2,0], M[0,0])` gives θ
  only modulo π — and it lands on the wrong branch here, because `A` has a negative
  leading element. The notebook picks the branch by requiring `A[0,1] > 0` (a
  net-drifting channel). A much stronger solenoid, with Larmor phase advance past π,
  would need a different discriminant.
- **The symplectic residual sits around 1e-3 in both codes.** `ζ` and `δ` as defined
  here are conjugate only up to a factor, and the finite-difference steps are 1e-3, so
  that is the method, not the code. Compare the two codes' residuals to each other.
- **The pointwise `Bz` difference is a sampling offset, not a field error.** OPALX
  reports `Bz_ref` at the reference particle but tags it with that step's path length
  `s`; a sub-0.1 mm mismatch between the two is enough to show up as a few 1e-4 T
  wherever the field is steep. The tell is that peak `Bz` and `∫Bz ds` agree to 1e-6
  while the pointwise max is 1e-3 of peak. The notebook converts the difference into
  the longitudinal offset that would explain it — median 0.2 µm, 25 µm at the 99th
  percentile — so it can be read for what it is.
- **`CHARGE = 1` on the OPALX BEAM is mandatory.** `ParticleProperties` maps `MUON` to
  charge −1, and `Beam::execute` only consults that table when the attribute is absent.
  A charge sign error is indistinguishable from a map orientation error in every
  downstream number, so `run.log` is checked for `CHARGE      +e * 1`.
- **G4BL's tracking parameter is `maxStep`, capital S.** `mue4_WsxOn.g4bl` writes
  `param maxstep=10`, which is a different, unused parameter — that model is running on
  the 100 mm default.
- **`fieldmap` places no physical volume**, so nothing expands the G4BL world. Without
  the vacuum `WORLDPIPE` every track leaves immediately and the ntuples come out empty.
- **`start` must come before every `place`** in a G4BL input, or it aborts with
  `Invalid start` — it refuses to run once a centerline segment exists.
- **G4BL ASCII output is `%.6g`** — 6 significant digits. That's a ~5e-7 relative floor
  on every comparison; don't read anything into differences below it.
- **OPALX monitor `.h5` files have several `Step#N` groups**, one per tracking step
  during which particles crossed the plane. All of them have to be concatenated, then
  sorted by `id`. `time` is in seconds, `x` in m, `px` in βγ.
- **A monitor plane inside the field, or past `ZSTOP`, records nothing** and does not
  complain. The notebook checks the particle count at every plane.

## Reference numbers for the WSX on-axis map

```
peak Bz  0.260443 T at map-local z = +350 mm   (the map is NOT z-symmetric)
int Bz dz  0.259104 T*m  ->  Leff = 0.995 m
Brho       0.093398 T*m  ->  |Larmor angle| = 1.3871 rad = 79.5 deg
```
Sign of the Larmor angle is negative for the `rotation=Y180` placement.

## Where wsx_solenoid stands

Passing, with OPALX reading the `.g4blmap` natively.

| | OPALX | G4BL |
|---|---|---|
| peak Bz on axis | −0.260442 T at s = 1.4499 m | −0.260443 T at s = 1.4500 m |
| ∫Bz dz | −0.259104 T·m | −0.259104 T·m |

Particle by particle, all 13, ids matched:

| | entrance (s = 0.150 m) | exit (s = 3.450 m) |
|---|---|---|
| max \|Δx, Δy\| | 1.4e-15 mm | 7.0e-06 mm |
| max \|Δpx, Δpy\| | 1.6e-15 MeV/c | 4.2e-08 MeV/c |
| max \|Δpz\| | 4.5e-12 MeV/c | 1.5e-05 MeV/c |

Larmor rotation of the position vector: −10.5325° (OPALX) vs −10.5327° (G4BL),
worst particle 0.0005° apart. The entrance plane is pure plumbing and comes out at
round-off, as it should.

## The Gaussian case

`wsx_solenoid_gauss/` — same solenoid, same geometry, same frame; only the beam
differs. 15000 particles, sigma_x = 10 mm, sigma_y = 3 mm, sigma_x' = sigma_y' = 2 mrad,
monoenergetic at 28 MeV/c. Eight recording planes at z = 150, 600, 1050, 1500, 1950,
2400, 2850, 3450 mm — the first and last field-free, the middle six inside the field so
the rotation can be watched developing.

`compare_gauss.ipynb` writes one figure per plane with the same three panels as
`runs/processing/plot_monitors.py` (x-x', y-y', x-y), OPALX on the top row and G4BL on
the bottom. It also plots the measured x-y ellipse tilt at every plane against the
field-only prediction.

Panels are direct particle scatters with `alpha = 0.18` — lower than the house 0.35,
which at 15000 points saturates the core to a solid block. Axis limits are shared down
each column so the two rows can be read against each other.

**The sample mean is subtracted exactly** when the beam is generated. OPALX under
`FROMFILE` takes its reference orbit from the bunch mean, and a raw draw of 5000 is off
by ~sigma/sqrt(15000) (0.08 mm in x) — enough to tilt the OPALX reference away from
G4BL's fixed `reference particle` and offset every comparison. The same reason the
13-particle set has to stay symmetric.

Two things to know when reading the rotation plot:

- **Sign.** With the block convention `M = [[cos*A, sin*A], [-sin*A, cos*A]]`, the lab
  beam turns by `-theta_L`, not `+theta_L` — the (x,y) pair is acted on by `R(-theta)`.
  The notebook compares the tilt against `-int Bz ds / 2*Brho`.
- **The tilt is degenerate at a round waist.** At z = 2400 mm the beam is nearly
  circular (rms 2.03 x 2.13 mm) and the ellipse orientation carries little information.
  It happens to land close to the prediction at 15000 particles (+75.02 vs +75.34 deg)
  but wandered 14 deg off it at 5000 — it is sampling noise either way, not a result.
  Those planes are ringed in the plot.

Result: the two codes agree on the tilt to **<= 5e-4 deg** at all eight planes and on
rms x and y to the printed precision, with all 15000 particles surviving to every plane.
