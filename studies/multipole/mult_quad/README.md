# mult_quad

**MULTIPOLE**, L = 0.3 m. MULTIPOLE KN = {0, 2}: must reproduce `quad_defoc_x` exactly.
0.1 GeV kinetic electron (Bρ = 0.3353 T·m), no space charge.

## Physics

The normalised strength `K1` [m⁻²] is turned into a physical gradient `G = K1·Bρ` [T/m] (`OpalQuadrupole::update`, `OpalMultipole::update`), giving `By = G·x`, `Bx = G·y` for the normal component and `Bx = −Gs·x`, `By = +Gs·y` for the skew one.

The equations of motion are `x″ = −(q/p)·By`, `y″ = +(q/p)·Bx`, so the effective focusing strength carries the **sign of the charge**:

```
k = (q/p)·G = sign(q)·K1        x″ = −k·x,  y″ = +k·y
```

This beam is an **electron** (`CHARGE = −1`), so a positive `K1` **defocuses** in x and focuses in y, opposite to the MAD-X convention for a positive charge. Here k = -1·K1 = -2 m⁻².

With √|k| = 1.4142 m⁻¹ and L = 0.3 m the phase advance is √|k|·L = 0.4243 rad. The focusing plane (**x**) is

```
[  cos(√k L)      sin(√k L)/√k ]   [ +0.911342  +0.291081 ]
[ -√k sin(√k L)   cos(√k L)    ] = [ -0.582161  +0.911342 ]
```

and the defocusing plane (**y**) is the same with cos→cosh, sin→sinh and the R21 sign flipped:

```
[  cosh(√k L)     sinh(√k L)/√k ]   [ +1.091358  +0.309081 ]
[  √k sinh(√k L)  cosh(√k L)    ] = [ +0.618163  +1.091358 ]
```

A quadrupole has no dispersion (R16 = R26 = 0) and the longitudinal row is the drift one, R56 = L/γ² = 7.754e-06 m.

A physically correct (symplectic) map satisfies **MᵀJM = J** and, with the planes uncoupled, each 2×2 diagonal block has determinant 1.

## Setup (`mult_quad.in`)

- `FIELDSOLVER, TYPE = NONE` (no space charge).
- `DISTRIBUTION, TYPE = FROMFILE` reading `parts.txt`: **13 particles** — a reference at the design momentum plus a symmetric **± step in each of the six phase-space coordinates** (εx=εy=εz=0.0001 m, εx′=εy′=0.0001 rad, εδ=0.001). No `PC` on the BEAM (FROMFILE reads the momenta from the file); the field reference momentum comes from the `REAL P0` deck variable.
- Lattice: `DRIFT(L=0.5) → MULTIPOLE(L=0.3) → DRIFT(L=0.5)`. Faces at s = 0.500 .. 0.800 m. `MULTIPOLE` has **no fringe** — the field is exactly zero outside the faces.
- `DT = 1e-13` s. The hard edge is crossed between time steps, so the effective length is uncertain by ~βcΔt = 30 µm; at this DT every matrix element is converged to <1e-5 relative.

## How the output is processed (`multlib.py` → `run_tests.py`)

- `Case.read_plane` reads the co-moving phase space from `mult_quad.h5`: positions in the beam frame [m], momenta in β·γ rotated into the reference direction → x′ = px/pz, y′ = py/pz, δ = |p|/p₀ − 1.
- `Case.transfer_map` builds the 6×6 map: take a dump in the drift on each side of the magnet, project every particle onto the reference transverse plane (drift by −z), strip the exact drift from each plane to the face, then centered-finite-difference the ± pairs.
- **Tests fed:** 1 (energy), 2 (transfer matrix), 3 (symplecticity), 4 (element equivalence).
