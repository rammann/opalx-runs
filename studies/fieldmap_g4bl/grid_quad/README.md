# grid_quad

Quadrupole in the 3D cartesian grid format, on a posed FIELDMAP.

## Field

Quadrupole, By = g x and Bx = g y with g = 0.083816 T/m, so k1 = -0.2500 m^-2 (negative: the electron's charge makes a positive gradient defocus in x). Linear in one coordinate each, so trilinear interpolation is exact. This field is both divergence free and curl free, so its transfer matrix is symplectic.

## Setup (`grid_quad.in`)

- Map `../maps/quad.g4blmap`, grid format, carried by a `FIELDMAP`.
- The map's own z = 0 sits at lab z = 1.0 m, so the field runs from
  0.5 m to 1.5 m along the line.
- Multiplier on the tabulated field: 1.0.
- 0.1 GeV electron, no space charge, time step 1e-13 s, tracking stops at
  2.0 m.
- 13 particles: a reference plus a +/- step in each coordinate, so the transfer
  matrix comes from centred differences.

## How the output is checked (`fmlib.py` -> `run_tests.py`)

The 4x4 transverse matrix against the closed form, and its symplecticity.
