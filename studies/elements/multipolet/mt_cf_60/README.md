# mt_cf_60

**MULTIPOLET** sector bend, **60 deg**, tanh fringe (LFRINGE = RFRINGE = 0.02 m). 0.1 GeV electron, no space charge.

## Physics

The design particle rides an arc of radius rho = L/theta = 0.9549 m turning through theta = 60 deg. Rigidity Brho = 0.3353 T m, so the dipole field is `TP[0] = Brho/rho = -0.3511 T`, `TP[1] = +0.5 T/m`. A magnetic field does no work, so |p| is conserved (test 1).

The body is a sector bend with the centre of curvature on the -x side, the same convention as SBEND, so the analytic map is the plain sector matrix.

Combined function: gradient TP[1] = 0.5 T/m, so k1 = q/|q| B1/Brho = -1.4914 1/m^2 (a positive gradient defocuses a negative charge in x, the same rule as MULTIPOLE K1). Expected: kx = 1/rho^2 + k1 = -0.3947, ky = -k1 = +1.4914; the map is the combined-function sector matrix (cos/sin for k > 0, cosh/sinh for k < 0) with dispersion (1 - cos(sqrt(kx) L))/(rho kx) (test 11).

## Setup (`mt_cf_60.in`)

- `FIELDSOLVER, TYPE = NONE`.
- `DISTRIBUTION, TYPE = FROMFILE` reading `parts.txt`: 13 particles, a reference plus a +/- step in each phase-space coordinate (eps x,y,z = 0.0001 m, eps x',y' = 0.0001 rad, eps delta = 0.001). The field reference momentum is the `REAL P0` deck variable.
- Lattice: `DRIFT(L=1) -> MULTIPOLET(L=1, ANGLE=60 deg) -> DRIFT`. Design faces at s = 1.000 .. 2.000 m.

## How the output is processed (`bendlib.py` -> `run_tests.py`)

- `Case.read_plane` reads the co-moving phase space from `mt_cf_60.h5` (positions [m] relative to the reference particle, momenta in beta*gamma rotated into the reference direction): x' = px/pz, y' = py/pz, delta = |p|/p0 - 1.
- `Case.transfer_map` takes a field-free plane on each side of the bend (field extent from `By_ref` in `mt_cf_60.stat`), projects every particle onto the reference transverse plane, strips the drift to the design face, and finite-differences the +/- pairs.
- The reference orbit (`ref_x/ref_z/ref_px/ref_pz`) gives the bend angle (test 2).
