# mult_skew

**MULTIPOLE**, L = 0.3 m. MULTIPOLE KS = {0, 2}: must reproduce `quad_skew` exactly.
0.1 GeV kinetic electron (Bρ = 0.3353 T·m), no space charge.

## Physics

The normalised strength `K1` [m⁻²] is turned into a physical gradient `G = K1·Bρ` [T/m] (`OpalQuadrupole::update`, `OpalMultipole::update`), giving `By = G·x`, `Bx = G·y` for the normal component and `Bx = −Gs·x`, `By = +Gs·y` for the skew one.

The equations of motion are `x″ = −(q/p)·By`, `y″ = +(q/p)·Bx`, so the effective focusing strength carries the **sign of the charge**:

```
k = (q/p)·G = sign(q)·K1        x″ = −k·x,  y″ = +k·y
```

This beam is an **electron** (`CHARGE = −1`), so a positive `K1` **defocuses** in x and focuses in y, opposite to the MAD-X convention for a positive charge. Here ks = -1·K1S = -2 m⁻².

A **skew** quadrupole couples the planes: `x″ = −ks·y`, `y″ = −ks·x`. In the 45°-rotated coordinates u = (x+y)/√2, v = (x−y)/√2 it decouples into `u″ = −ks·u`, `v″ = +ks·v` — i.e. an ordinary normal quadrupole of the same strength. So the analytic map is **Rᵀ·Q(ks)·R** with R the 45° rotation (`multlib.skew_quad_map`).

A physically correct (symplectic) map satisfies **MᵀJM = J**. The planes are coupled here, so the 2×2 diagonal blocks are *not* unit-determinant; the 4×4 transverse block is the one with determinant 1.

## Setup (`mult_skew.in`)

- `FIELDSOLVER, TYPE = NONE` (no space charge).
- `DISTRIBUTION, TYPE = FROMFILE` reading `parts.txt`: **13 particles** — a reference at the design momentum plus a symmetric **± step in each of the six phase-space coordinates** (εx=εy=εz=0.0001 m, εx′=εy′=0.0001 rad, εδ=0.001). No `PC` on the BEAM (FROMFILE reads the momenta from the file); the field reference momentum comes from the `REAL P0` deck variable.
- Lattice: `DRIFT(L=0.5) → MULTIPOLE(L=0.3) → DRIFT(L=0.5)`. Faces at s = 0.500 .. 0.800 m. `MULTIPOLE` has **no fringe** — the field is exactly zero outside the faces.
- `DT = 1e-13` s. The hard edge is crossed between time steps, so the effective length is uncertain by ~βcΔt = 30 µm; at this DT every matrix element is converged to <1e-5 relative.

## How the output is processed (`multlib.py` → `run_tests.py`)

- `Case.read_plane` reads the co-moving phase space from `mult_skew.h5`: positions in the beam frame [m], momenta in β·γ rotated into the reference direction → x′ = px/pz, y′ = py/pz, δ = |p|/p₀ − 1.
- `Case.transfer_map` builds the 6×6 map: take a dump in the drift on each side of the magnet, project every particle onto the reference transverse plane (drift by −z), strip the exact drift from each plane to the face, then centered-finite-difference the ± pairs.
- **Tests fed:** 1 (energy), 2 (transfer matrix), 3 (symplecticity), 4 (element equivalence).
