# cyl_ramp

The asymmetric solenoid profile, read the way it is written.

## Field

A solenoid profile that is deliberately not symmetric in z, with the radial field from the usual first-order expansion. Turning it round is a different field, which is what makes it a fair test of ZREVERSE.

## Setup (`cyl_ramp.in`)

- Map `../maps/ramp_cyl.g4blmap`, cylinder format, carried by a `FIELDMAP`.
- The map's own z = 0 sits at lab z = 1.0 m, so the field runs from
  0.5 m to 1.5 m along the line.
- Multiplier on the tabulated field: 1.0.
- 0.1 GeV electron, no space charge, time step 1e-13 s, tracking stops at
  2.0 m.
- 13 particles: a reference plus a +/- step in each coordinate, so the transfer
  matrix comes from centred differences.

## How the output is checked (`fmlib.py` -> `run_tests.py`)

Provides the field profile the reversed cases are compared with; its own field window and |p| are checked.
