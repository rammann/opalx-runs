# MULTIPOLET as a bend

Does `MULTIPOLET` with `ANGLE != 0` act as a sector bend in opal-t tracking? Same method as
the SBEND/RBEND suite in `runs/bendtest/validation`: a 13-particle finite-difference transfer
map and Gaussian-bunch covariance tests, checked against analytic sector-bend optics. `SBEND`
with the Enge fringe is tracked next to it as the reference bend.

0.1 GeV kinetic electron, arc length 1 m, 30, 60 and 90 deg, no space charge.

## Short answer

Yes, after the element was rewritten to the structure of the other elements
(branch `537-update-multipolet-and-variablerfcavity-to-conform-with-opalx-architecture`).
**225 of 225 checks pass**, MULTIPOLET as close to the analytic optics as SBEND is.

The element was rewritten because the original did not work as a bend at all:

- It did not give the reference particle any field (no `applyToReferenceParticle`), so the
  design orbit went straight while the tracked particles bent away.
- Its field support was the bare body `[0, L)`, so the entrance half of the fringe was never
  applied: the field started at 0.53 of the flat top right at the face, the exit angle was
  0.6 % short, and the extracted map was not symplectic (residual 0.26).
- It added the local Frenet field to the particle without rotating it into the entrance
  frame, which gave a spurious vertical focusing `R43 = -sin(theta)/rho` (-0.26 at 30 deg,
  -0.89 at 60 deg).
- Its transverse coordinate ran the other way from its own geometry, so a gradient focused in
  both planes at once (impossible for a source-free field), the horizontal gradient sign was
  the opposite of every other magnet, and anything placed after a curved MULTIPOLET sat on the
  straight extension of the entrance line.
- Its longitudinal window was measured along the entrance tangent instead of the arc, so at
  90 deg the element never left the active set and the orbit threader looped forever.

All of those are gone. The rewritten element keeps its field math in
`opalx/src/AbsBeamline/MultipoleTFieldModel.h` (a header-only namespace like
`BendFieldModel.h`), its geometry in `MultipoleTRep`, and it overrides the same `ElementBase`
methods `SBend` does, so the tracker, the placement, the aperture scraping and the element
position output treat it as the bend it is.

## Run it

```bash
./run_all.sh                        # write the inputs, track all cases, run the tests
./run_all.sh --test-only            # re-run only the analysis on existing output
```

Tracking takes about 2 min with a Release build; with a Debug build each 8000-particle
case alone takes about 11 min. `run_all.sh` takes the binary from `$OPALX`, or the
workspace build `/Users/rammann/Code/OPALX/build/src/opalx` when that is not set, and
the conda `python` (`/opt/homebrew/Caskroom/miniconda/base/bin/python`, has
h5py/numpy/scipy/matplotlib; override with `PY=...`). The tracking output goes to
`output/elements/multipolet/<case>/`. The table goes to stdout and `results.txt`; figures
to `output/elements/multipolet/plots/30deg/` and `.../60deg/`, one per test.

## Files

| file | what it does |
|------|--------------|
| `make_decks.py` | writes every `<case>/<case>.in` (+ `parts.txt`, per-case `README.md`) and the `cases.json` manifest |
| `bendlib.py` | readers, finite-difference map builder, analytic (combined-function) sector matrix, symplecticity, covariance transport |
| `run_tests.py` | the 11 tests, prints the table, writes `results.txt`, calls `plot_tests.py` |
| `plot_tests.py` | one figure per test per angle |
| `run_all.sh` | one-command driver |
| `<case>/` | input, `parts.txt` (map cases), `README.md` (physics, setup, how the output is checked); the tracking output is in `output/elements/multipolet/<case>/` |

## Cases

| case | element | fringe | what for |
|------|---------|--------|----------|
| `mt_fringe_{30,60,90}` | MULTIPOLET | tanh, LFRINGE = RFRINGE = 0.02 m | primary map cases |
| `mt_sharp_{30,60,90}` | MULTIPOLET | tanh 0.005 m | sharper edge, same tests |
| `sbend_enge_{30,60,90}` | SBEND | Enge, HGAP = 0.01 m | reference bend |
| `mt_cf_{30,60}` | MULTIPOLET | tanh 0.02 m, TP[1] = 0.5 T/m | combined function (test 11) |
| `mt_bunch_{emit,disp,pure,mean}_{30,60}` | MULTIPOLET | tanh 0.02 m | 8000-particle bunches for tests 5, 8, 9, 10 |

All bends: arc length L = 1 m, rho = L/ANGLE, dipole field |B0| = Brho/rho (0.1755 T at
30 deg, 0.3511 T at 60 deg, 0.5266 T at 90 deg), Brho = 0.3353 T m. The MULTIPOLET decks set
`TP = {B0, B1}` with **B0 negative for the electron**, `HAPERT = 0.4`, `VAPERT = 0.2`,
`MAXFORDER = 3`.

## Conventions

`TP` is the physical mid-plane field profile in tesla, sign included: `By = TP[0] + TP[1] x`
in the flat top. `ANGLE` and `L` only set the arc geometry, with the centre of curvature on the
**-x** side, the same as SBEND, so a positive `ANGLE` turns the orbit towards -x. The deck has
to keep `TP[0]` consistent with that: for a negative charge and a positive `ANGLE`,
`TP[0] = -Brho/rho`, which is the same physical field SBEND builds internally from `ANGLE`,
`P0` and the charge.

The gradient follows the ordinary rule as well: `k1 = sign(q) TP[1] / Brho`, so a positive
`TP[1]` defocuses an electron horizontally, exactly as a positive `K1` does on a `MULTIPOLE`
or `QUADRUPOLE` (see the multipole study).

`HAPERT` and `VAPERT` are the **full** width and height of a rectangular aperture and are
installed as the element's aperture, so particles outside them are scraped like anywhere
else. They are alternatives to the generic `APERTURE` string, not additions to it.
`LFRINGE = RFRINGE = 0` is a hard edge, like `HGAP = 0` on an SBEND.

Phase space is handled as in the SBEND suite: the `.h5` holds co-moving positions and
beta*gamma momenta, we form (x, x' = px/pz, y, y' = py/pz, z, delta), read a field-free plane
on each side of the bend, drift-project every particle onto the reference plane, strip the
drift to the design faces, and centred-difference the +/- pairs. Both elements now share the
co-moving x convention, so the analytic map is the plain sector matrix.

## The 11 tests

| # | check | how |
|---|-------|-----|
| 1 | energy conservation | max per-particle change of \|p\|, first to last dump |
| 2 | bend angle | theta = int By_ref ds / Brho vs the orbit exit angle atan2(px, pz), and the exit angle vs the deck ANGLE |
| 3 | transfer matrix | 6x6 map by finite differences: bending block, dispersion R16 R26, vertical block, vs the sector matrix |
| 4 | dispersion | R16 = rho (1 - cos theta), R26 = sin theta |
| 5 | RMS z from R56 | energy-spread-only bunch: sigma_z = \|R56\| sigma_delta, sigma_x = \|D\| sigma_delta |
| 6 | vertical plane | sector bend: R43 within the fringe residual and R34 = L, MULTIPOLET next to SBEND |
| 7 | symplecticity | max\|M^T J M - J\| and the three 2x2 block determinants |
| 8 | emittance invariance | eps_x, eps_y and the 6D emittance in vs out |
| 9 | RMS x from energy spread | sigma_x,out from M Sigma M^T and from betatron + (D sigma_delta)^2 |
| 10 | centroid from mean energy | off-energy bunch orbit shift = D <delta> |
| 11 | combined function | full map of `mt_cf_*` vs the combined-function sector matrix (kx = 1/rho^2 + k1, ky = -k1) |

## Results (results.txt)

225 of 225 checks pass. The map cases, per angle and fringe:

| quantity | mt_fringe_30 | mt_sharp_30 | sbend_enge_30 | mt_fringe_60 | mt_fringe_90 | sbend_enge_90 |
|---|---|---|---|---|---|---|
| exit angle - ANGLE [rad] | -2.2e-5 | -1.4e-6 | -1.2e-5 | -1.5e-4 | -4.0e-4 | -2.3e-4 |
| R11 (analytic cos) | 0.8678 (0.8660) | 0.8665 | 0.8675 | 0.5063 (0.5) | 0.0109 (0) | — |
| R21 (-sin/rho) | -0.2585 (-0.2618) | -0.2610 | -0.2592 | -0.8989 (-0.9069) | -1.5693 (-1.5708) | — |
| R16 (rho(1-cos)) | 0.25587 (0.25587) | 0.25587 | 0.25587 | 0.47747 (0.47747) | 0.63662 (0.63662) | — |
| R43 (sector: 0) | 0.0017 | 0.0004 | 0.0038 | 0.0069 | 0.0158 | 0.0347 |
| symplectic residual | 3.1e-5 | 3.0e-6 | 5.4e-7 | 5.4e-5 | 1.3e-4 | 1.7e-6 |

Bunch tests (MULTIPOLET, tanh 0.02 m): emittances in/out agree to 1e-4, sigma_x from
M Sigma M^T to 1e-5, sigma_z = |R56| sigma_delta to 6e-4, centroid shift D <delta> to 3e-3.
Combined function (`mt_cf_*`): the whole 6x6 map matches the analytic combined-function sector
matrix to a few 1e-3, vertical block included — `R33 = 0.3431` against `0.3425` for
`ky = -k1`, where a mirrored transverse field would have given `1.843`.

### About R43

A thin-edge sector bend has no vertical focusing at all, but any real extended fringe leaves a
small residual that grows with the bend angle, so test 6 bounds it relative to the
rectangular-bend edge focusing `tan(theta/2)/rho` of the same geometry rather than against a
fixed number. Measured: about 1 % of that scale for the tanh fringe at every angle, about
2.5 % for the Enge fringe (the SBEND reference, whose fringe carries the FINT edge term the
tanh model does not have). The old absolute bound of 2e-2 happened to pass at 30 and 60 deg
and failed the SBEND reference at 90 deg.

## Tolerances

Set in `run_tests.py` (`TOL`): energy 1e-6, field-integral angle 2e-3 rad, exit angle vs
ANGLE 1e-3 rad, matrix elements 3e-2 relative with a 2e-2 absolute floor, dispersion 2e-2,
R43 5 % of tan(theta/2)/rho, symplectic residual 1e-4, emittance 5e-3 (6D 1e-2), RMS growth
2e-2, centroid 1e-2.
