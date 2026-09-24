# mt_bunch_emit_60

**MULTIPOLET** sector bend, **60 deg**, tanh fringe 0.02 m. Gaussian-bunch case for test 8 -- emittance invariance. No space charge; linear.

## Physics

The design particle rides an arc of radius rho = L/theta = 0.9549 m turning through theta = 60 deg. Rigidity Brho = 0.3353 T m, so the dipole field is `TP[0] = Brho/rho = -0.3511 T`. A magnetic field does no work, so |p| is conserved (test 1).

The body is a sector bend with the centre of curvature on the -x side, the same convention as SBEND, so the analytic map is the plain sector matrix.

The covariance transforms as Sigma_out = M Sigma_in M^T with the map M measured in `mt_fringe_60`. With zero energy spread the projected RMS emittances eps_x and eps_y and the 6D emittance sqrt(det Sigma) are invariant under the symplectic map (test 8).

## Setup (`mt_bunch_emit_60.in`)

- `FIELDSOLVER, TYPE = NONE`.
- `DISTRIBUTION, TYPE = GAUSS`, N = 8000. sigma x,y = 0.001 m, sigma x',y' = 0.0002 rad, sigma z = 0.0005 m, no energy spread.
- Same lattice as `mt_fringe_60`.
