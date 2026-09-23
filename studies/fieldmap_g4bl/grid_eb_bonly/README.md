# grid_eb_bonly

The same combined map with ESCALE = 0, so only the magnetic half acts.

## Field

Uniform Ez = 1 MV/m over 1 m, written in the last three columns of the nine-column form. A particle crossing it gains q E dz, which for a field with no z dependence is -1.0000 MeV for this electron, whatever path it takes. Constant, so trilinear interpolation is exact. Here ESCALE = 0 switches the electric field off entirely.

## Setup (`grid_eb_bonly.in`)

- Map `../maps/eb.g4blmap`, grid format, carried by a `FIELDMAP`.
- The map's own z = 0 sits at lab z = 1.0 m, so the field runs from
  0.5 m to 1.5 m along the line.
- Multiplier on the tabulated field: 1.0.
- 0.1 GeV electron, no space charge, time step 1e-13 s, tracking stops at
  2.0 m.
- 13 particles: a reference plus a +/- step in each coordinate, so the transfer
  matrix comes from centred differences.

## How the output is checked (`fmlib.py` -> `run_tests.py`)

Must reproduce grid_dipole to round-off. That is the sharp test that the two scales are independent: the same file, one field switched off, has to track like the magnet-only map.
