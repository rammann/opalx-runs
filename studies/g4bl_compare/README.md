# g4bl_compare — checking OPALX elements against G4beamline

One folder per element. Each holds a G4BL input, an OPALX input describing the same
element, and a notebook that pushes the **same particles** through both and compares
the outputs particle by particle.

The point is attribution. Comparing whole beamlines only says "roughly agrees";
comparing one element in isolation says which element, and which term of its map, is
wrong.

```
g4bl_compare/
  wsx_solenoid/
    wsx_solenoid.g4bl    G4BL: the muE4 WSX solenoid alone
    wsx_solenoid.in      OPALX: the same solenoid alone
    compare.ipynb        writes the particle files, reads both outputs, compares
```

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
conversion, no derived map files. **This needs the `.g4blmap` reader** — without it
OPALX stops at `Couldn't determine type of fieldmap`. A commented fallback line in
`wsx_solenoid.in` points at the already-converted `opalx-runs/data/wsx_solenoid_1D.map`,
which is enough to exercise the comparison chain: that map carries *absolute* lab z, so
it needs `ELEMEDGE = 0.515` to cover the same 0.300–3.300 m span. It will not agree
physically — it has no Y180 applied (see below).

**Nothing else can differ.** Vacuum world, no iron, no decay, no stochastic processes,
no space charge, no scraping.

## What the notebook reports

- per-particle exit `x, y, x', y', ζ, δ` from both codes and the difference
- the same at the entrance plane — pure plumbing, and it should be zero
- the 6×6 transfer matrix from each code, side by side, plus the element-wise
  difference and the symplectic residual
- the Larmor rotation angle, **with sign**, three ways: from each matrix, from the x+
  particle's trajectory, and from ∫Bz ds / (2Bρ)
- OPALX `Bz_ref` from the `.stat` file against G4BL's on-axis `fieldntuple` dump —
  no tracking is involved, so a mismatch there is a field-import problem and
  everything below it is meaningless until it passes

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
  `∫Bz dz = −0.259104 T·m`, i.e. mirrored and negative.
  `opalx-runs/data/wsx_solenoid_1D.map`, produced by
  `runs/mue4_fieldmaps/tools/g4blmap_solenoid_to_astra.py`, applies neither: its peak
  sits 700 mm the other side of centre and `∫Bz ds` comes out **positive**. Running the
  fallback through this notebook reproduces exactly that — same |Larmor angle| (1.3883
  vs 1.3871 rad from the field), opposite sign. The Larmor-angle sign is the sharpest
  test of it.
- **The Larmor angle has to be read off the matrix carefully.** `M[0:2,0:2] = cosθ·A`
  and `M[2:4,0:2] = −sinθ·A` for a common 2×2 `A`, so `atan2(−M[2,0], M[0,0])` gives θ
  only modulo π — and it lands on the wrong branch here, because `A` has a negative
  leading element. The notebook picks the branch by requiring `A[0,1] > 0` (a
  net-drifting channel). A much stronger solenoid, with Larmor phase advance past π,
  would need a different discriminant.
- **The symplectic residual sits around 1e-3 in both codes.** `ζ` and `δ` as defined
  here are conjugate only up to a factor, and the finite-difference steps are 1e-3, so
  that is the method, not the code. Compare the two codes' residuals to each other.
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
