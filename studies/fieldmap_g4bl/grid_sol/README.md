# grid_sol

Uniform Bz in the 3D cartesian grid format, on a posed FIELDMAP.

## Field

Uniform Bz = 0.1 T over 1 m, which turns the transverse momentum through k L = -0.298272 rad. There is no radial field at the ends, so this is not a focusing solenoid and not a divergence-free field; it is used because both formats carry it exactly and the motion has a closed form.

## Setup (`grid_sol.in`)

- Map `../maps/sol.g4blmap`, grid format, carried by a `FIELDMAP`.
- The map's own z = 0 sits at lab z = 1.0 m, so the field runs from
  0.5 m to 1.5 m along the line.
- Multiplier on the tabulated field: 1.0.
- 0.1 GeV electron, no space charge, time step 1e-13 s, tracking stops at
  2.0 m.
- 13 particles: a reference plus a +/- step in each coordinate, so the transfer
  matrix comes from centred differences.

## How the output is checked (`fmlib.py` -> `run_tests.py`)

The 4x4 transverse matrix against the closed form, and against the cylinder map holding the same field.
