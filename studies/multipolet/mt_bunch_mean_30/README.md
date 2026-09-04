# mt_bunch_mean_30

**MULTIPOLET** sector bend, **30 deg**, tanh fringe 0.02 m. Gaussian-bunch case for test 10 -- centroid shift from mean energy. No space charge; linear.

## Physics

The design particle rides an arc of radius rho = L/theta = 1.9099 m turning through theta = 30 deg. Rigidity Brho = 0.3353 T m, so the dipole field is `TP[0] = Brho/rho = -0.1755 T`. A magnetic field does no work, so |p| is conserved (test 1).

The body is a sector bend with the centre of curvature on the -x side, the same convention as SBEND, so the analytic map is the plain sector matrix.

The covariance transforms as Sigma_out = M Sigma_in M^T with the map M measured in `mt_fringe_30`. The bunch mean momentum is P0 (1 + <delta>) while the field reference stays P0, so the bunch is uniformly off-energy and its centroid disperses: <x>_out = D <delta> (test 10). OPALX's reference particle rides with its own bunch, so the shift is read from the reference orbit against the on-energy `mt_bunch_emit_30` orbit.

## Setup (`mt_bunch_mean_30.in`)

- `FIELDSOLVER, TYPE = NONE`.
- `DISTRIBUTION, TYPE = GAUSS`, N = 8000. sigma x,y = 0.001 m, sigma x',y' = 0.0002 rad, sigma z = 0.0005 m, no energy spread.
- Mean energy offset <delta> = 0.002 via `PC = P0 (1 + <delta>)`.
- Same lattice as `mt_fringe_30`.
