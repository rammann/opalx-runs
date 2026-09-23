# grid_efield_escale2

The same electric map with ESCALE = 2 on the element.

## Field

Uniform Ez = 1 MV/m over 1 m, written in the last three columns of the nine-column form. A particle crossing it gains q E dz, which for a field with no z dependence is -1.0000 MeV for this electron, whatever path it takes. Constant, so trilinear interpolation is exact. ESCALE = 2 is applied on top.

## Setup (`grid_efield_escale2.in`)

- Map `../maps/efield.g4blmap`, grid format, carried by a `FIELDMAP`.
- The map's own z = 0 sits at lab z = 1.0 m, so the field runs from
  0.5 m to 1.5 m along the line.
- Multiplier on the tabulated field: 1.0.
- 0.1 GeV electron, no space charge, time step 1e-13 s, tracking stops at
  2.0 m.
- 13 particles: a reference plus a +/- step in each coordinate, so the transfer
  matrix comes from centred differences.

## How the output is checked (`fmlib.py` -> `run_tests.py`)

The energy change must be exactly twice grid_efield's, which is what says ESCALE is a plain multiplier on the tabulated field.
