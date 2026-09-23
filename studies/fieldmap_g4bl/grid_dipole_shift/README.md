# grid_dipole_shift

The same dipole map posed 0.25 m further downstream.

## Field

Uniform By = 0.058218 T over 1 m. A particle entering on the axis follows a circular arc of radius 5.7588 m and leaves at 10 degrees, offset 0.087489 m. A constant field is reproduced by trilinear interpolation exactly.

## Setup (`grid_dipole_shift.in`)

- Map `../maps/dipole.g4blmap`, grid format, carried by a `FIELDMAP`.
- The map's own z = 0 sits at lab z = 1.25 m, so the field runs from
  0.75 m to 1.75 m along the line.
- Multiplier on the tabulated field: 1.0.
- 0.1 GeV electron, no space charge, time step 1e-13 s, tracking stops at
  2.25 m.
- 13 particles: a reference plus a +/- step in each coordinate, so the transfer
  matrix comes from centred differences.

## How the output is checked (`fmlib.py` -> `run_tests.py`)

The field window must move by exactly 0.25 m and the bend angle must not change, which is what says the pose positions the map's own origin.
