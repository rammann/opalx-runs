# mt_bunch_disp_30

**MULTIPOLET** sector bend, **30 deg**, tanh fringe 0.02 m. Gaussian-bunch case for test 9 -- RMS x growth from energy spread. No space charge; linear.

## Physics

The design particle rides an arc of radius rho = L/theta = 1.9099 m turning through theta = 30 deg. Rigidity Brho = 0.3353 T m, so the dipole field is `TP[0] = Brho/rho = -0.1755 T`. A magnetic field does no work, so |p| is conserved (test 1).

The body is a sector bend with the centre of curvature on the -x side, the same convention as SBEND, so the analytic map is the plain sector matrix.

The covariance transforms as Sigma_out = M Sigma_in M^T with the map M measured in `mt_fringe_30`. With an energy spread sigma_delta, dispersion adds to the horizontal size: sigma_x_out^2 = betatron part + (D sigma_delta)^2 (test 9), the (1,1) element of Sigma_out = M Sigma_in M^T.

## Setup (`mt_bunch_disp_30.in`)

- `FIELDSOLVER, TYPE = NONE`.
- `DISTRIBUTION, TYPE = GAUSS`, N = 8000. sigma x,y = 0.001 m, sigma x',y' = 0.0002 rad, sigma z = 0.0005 m, energy spread 0.001.
- Same lattice as `mt_fringe_30`.
