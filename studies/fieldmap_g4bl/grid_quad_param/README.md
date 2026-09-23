# grid_quad_param

The plain quadrupole field with normB = 4 and current = 2 on the param line.

## Field

Quadrupole, By = g x and Bx = g y with g = 0.083816 T/m, so k1 = -0.2500 m^-2 (negative: the electron's charge makes a positive gradient defocus in x). Linear in one coordinate each, so trilinear interpolation is exact. This field is both divergence free and curl free, so its transfer matrix is symplectic. The param line asks for a factor normB/current = 2.

## Setup (`grid_quad_param.in`)

- Map `../maps/quad_param.g4blmap`, grid format, carried by a `FIELDMAP`.
- The map's own z = 0 sits at lab z = 1.0 m, so the field runs from
  0.5 m to 1.5 m along the line.
- Multiplier on the tabulated field: 1.0.
- 0.1 GeV electron, no space charge, time step 1e-13 s, tracking stops at
  2.0 m.
- 13 particles: a reference plus a +/- step in each coordinate, so the transfer
  matrix comes from centred differences.

## How the output is checked (`fmlib.py` -> `run_tests.py`)

Must reproduce grid_quad_2x to round-off. If either key were ignored the factor would come out 4 or 1/2 instead of 2.
