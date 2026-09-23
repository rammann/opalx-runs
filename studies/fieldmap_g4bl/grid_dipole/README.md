# grid_dipole

Uniform dipole in the 3D cartesian grid format, on a posed FIELDMAP.

## Field

Uniform By = 0.058218 T over 1 m. A particle entering on the axis follows a circular arc of radius 5.7588 m and leaves at 10 degrees, offset 0.087489 m. A constant field is reproduced by trilinear interpolation exactly.

## Setup (`grid_dipole.in`)

- Map `../maps/dipole.g4blmap`, grid format, carried by a `FIELDMAP`.
- The map's own z = 0 sits at lab z = 1.0 m, so the field runs from
  0.5 m to 1.5 m along the line.
- Multiplier on the tabulated field: 1.0.
- 0.1 GeV electron, no space charge, time step 1e-13 s, tracking stops at
  2.0 m.
- 13 particles: a reference plus a +/- step in each coordinate, so the transfer
  matrix comes from centred differences.

## How the output is checked (`fmlib.py` -> `run_tests.py`)

Bend angle and transverse offset of the reference orbit against the arc, the field window against the placement, and |p| against round-off.
