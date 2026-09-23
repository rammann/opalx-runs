# cyl_sol_solenoid

The same cylinder map on a SOLENOID placed by ELEMEDGE, the other placement convention.

## Field

Uniform Bz = 0.1 T over 1 m, which turns the transverse momentum through k L = -0.298272 rad. There is no radial field at the ends, so this is not a focusing solenoid and not a divergence-free field; it is used because both formats carry it exactly and the motion has a closed form.

## Setup (`cyl_sol_solenoid.in`)

- Map `../maps/sol_cyl.g4blmap`, cylinder format, carried by a `SOLENOID`.
- The map's own z = 0 sits at lab z = 1.0 m, so the field runs from
  0.5 m to 1.5 m along the line.
- Multiplier on the tabulated field: 1.0.
- 0.1 GeV electron, no space charge, time step 1e-13 s, tracking stops at
  2.0 m.
- 13 particles: a reference plus a +/- step in each coordinate, so the transfer
  matrix comes from centred differences.

## How the output is checked (`fmlib.py` -> `run_tests.py`)

Must reproduce cyl_sol to round-off, which is what says ELEMEDGE and the pose put the map's own z = 0 at the same place, and that KS = 1 reproduces the map as written.
