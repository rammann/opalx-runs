# mt_fringe_30

**MULTIPOLET** sector bend, **30 deg**, tanh fringe (LFRINGE = RFRINGE = 0.02 m). 0.1 GeV electron, no space charge.

## Physics

The design particle rides an arc of radius rho = L/theta = 1.9099 m turning through theta = 30 deg. Rigidity Brho = 0.3353 T m, so the dipole field is `TP[0] = Brho/rho = -0.1755 T`. A magnetic field does no work, so |p| is conserved (test 1).

The body is a sector bend with the centre of curvature on the -x side, the same convention as SBEND, so the analytic map is the plain sector matrix.

First-order motion is the 6x6 map in (x, x', y, y', z, delta). The bending plane of a sector bend is

```
[ cos      rho sin ]   [ +0.8660  +0.9549 ]
[ -sin/rho cos     ] = [ -0.2618  +0.8660 ]
```

the vertical plane is a drift of the arc length L = 1.0000 m, and the dispersion is |x| = rho(1 - cos) delta = 0.2559 delta, |x'| = sin delta = 0.5000 delta (tests 3, 4). No vertical focusing to first order (test 6). A symplectic map satisfies M^T J M = J (test 7).

## Setup (`mt_fringe_30.in`)

- `FIELDSOLVER, TYPE = NONE`.
- `DISTRIBUTION, TYPE = FROMFILE` reading `parts.txt`: 13 particles, a reference plus a +/- step in each phase-space coordinate (eps x,y,z = 0.0001 m, eps x',y' = 0.0001 rad, eps delta = 0.001). The field reference momentum is the `REAL P0` deck variable.
- Lattice: `DRIFT(L=1) -> MULTIPOLET(L=1, ANGLE=30 deg) -> DRIFT`. Design faces at s = 1.000 .. 2.000 m.

## How the output is processed (`bendlib.py` -> `run_tests.py`)

- `Case.read_plane` reads the co-moving phase space from `mt_fringe_30.h5` (positions [m] relative to the reference particle, momenta in beta*gamma rotated into the reference direction): x' = px/pz, y' = py/pz, delta = |p|/p0 - 1.
- `Case.transfer_map` takes a field-free plane on each side of the bend (field extent from `By_ref` in `mt_fringe_30.stat`), projects every particle onto the reference transverse plane, strips the drift to the design face, and finite-differences the +/- pairs.
- The reference orbit (`ref_x/ref_z/ref_px/ref_pz`) gives the bend angle (test 2).
