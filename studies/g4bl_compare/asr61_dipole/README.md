# asr61_dipole — muE4 ASR61 dipole, OPALX vs G4beamline

The first case in `g4bl_compare/` that tests a **bend** and the first that uses a
G4beamline **3D `grid`** map. `asr61_300d_track.g4blmap` is 60 x 25 x 553 points on a
10 x 10 x 5 mm grid, x = -200..390 mm, y = -120..120 mm, z = -1380..1380 mm, absolute
Tesla, 56 MB. OPALX reads it natively through `G4BL3DMagnetoStatic`; `FMAPSCALE` on the
bend is G4BL's deck-side `current=`, because neither reader normalises.

## Frame

`z` in the G4BL input [mm] / 1000 == path length `s` in the OPALX input [m].

OPALX places a map as `z_lab = ELEMEDGE + z_map`, so **`ELEMEDGE` is the lab position
of the map's own z = 0** — for this map its centre, not the field start.
`ELEMEDGE = 1.500` puts the field at s = 0.120..2.880 m, matching `place MAG z=1500`.

`ANGLE = 0` on the SBEND, and the G4BL centerline is straight. With `FMAPFN` the map
makes the field and `ANGLE` only places the elements downstream, so a straight geometry
keeps the element's local frame identical to the lab frame and the comparison is one to
one. The real `mue4_WsxOn.g4bl` instead wraps this magnet in two `cornerarc` commands
so the centerline follows the beam; reproducing that placement inside one OPALX bend
would need the map to sit on the bend's bisector, which the current convention does not
do.

The orbit bends towards +x and leaves the map sideways near x = 390 mm, where the map's
own box ends — in both codes. That is why the box reaches +390 mm on one side and only
-200 mm on the other.

## Running

```bash
python3 make_parts.py                             # parts.txt + beam.txt, one source

source /Users/rammann/Code/G4BL/G4beamline-3.08.app/Contents/root/bin/thisroot.sh
g4bl asr61_dipole.g4bl > g4bl.log 2>&1

mpirun -n 1 /Users/rammann/Code/OPALX/build/src/opalx asr61_dipole.in --info 1 > run.log 2>&1

python3 compare.py
```

One rank only — particle-to-row identity across ranks is not guaranteed. Both codes
must run from inside this folder; the map and particle paths are relative to the
working directory.

## What it showed

- On-axis field, 541 points from g4bl's `fieldntuple` against the same points from
  OPALX's reader: agreement to **5e-8 T**, which is g4bl's six-digit ASCII output
  precision. Peak `By = -0.0860619 T` at the map centre in both. The field import is
  exact; anything below is a tracking difference.
- Entrance plane at z = 50 mm: **zero** to the printed precision, i.e. the two codes
  are fed the same particles in the same frame.
- Exit plane at z = 2950 mm, after a 41.4 degree bend through 2.76 m of map:
  **|dx| <= 50 um out of 1.24 m** and **|dP| <= 0.6 keV/c out of 28 MeV/c**, both
  about 4e-5 relative. That is above g4bl's 10 um output quantisation, so it is a
  real difference — the two codes interpolate and integrate the same grid
  differently — but it is small enough that the map is clearly being read and applied
  the same way.
- The time step is converged: 2e-12 s reproduces 1e-11 s to five digits, so the
  residual is not OPALX's integration error.
- OPALX's monitors record 11 of the 13 particles, at both planes; g4bl records all 13.
  The two it drops are the lowest-momentum ones. It happens at the entrance plane too,
  before any field, so it is a monitor issue and not a field one.

`compare.py` matches on the `id` dataset. OPALX writes monitor rows in rank order, not
id order — comparing row to row pairs different particles and invents a difference of
tens of mm.
