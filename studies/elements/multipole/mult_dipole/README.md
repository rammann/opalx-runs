# mult_dipole

**MULTIPOLE**, L = 0.3 m. MULTIPOLE KN = {0.05}: pure dipole, checks the steering angle.
0.1 GeV kinetic electron (Bρ = 0.3353 T·m), no space charge.

## Physics

A `KN = {K0}` entry is the **dipole** component: a uniform `By = K0·Bρ = 0.0168 T` that steers the reference orbit by **θ = K0·L = 0.0150 rad**. This case exists to check that single number against the tracked reference orbit, so it does not feed the matrix test.

A physically correct (symplectic) map satisfies **MᵀJM = J** and, with the planes uncoupled, each 2×2 diagonal block has determinant 1.

## Setup (`mult_dipole.in`)

- `FIELDSOLVER, TYPE = NONE` (no space charge).
- `DISTRIBUTION, TYPE = FROMFILE` reading `parts.txt`: **13 particles** — a reference at the design momentum plus a symmetric **± step in each of the six phase-space coordinates** (εx=εy=εz=0.0001 m, εx′=εy′=0.0001 rad, εδ=0.001). No `PC` on the BEAM (FROMFILE reads the momenta from the file); the field reference momentum comes from the `REAL P0` deck variable.
- Lattice: `DRIFT(L=0.5) → MULTIPOLE(L=0.3) → DRIFT(L=0.5)`. Faces at s = 0.500 .. 0.800 m. `MULTIPOLE` has **no fringe** — the field is exactly zero outside the faces.
- `DT = 1e-13` s. The hard edge is crossed between time steps, so the effective length is uncertain by ~βcΔt = 30 µm; at this DT every matrix element is converged to <1e-5 relative.

## How the output is processed (`multlib.py` → `run_tests.py`)

- `Case.read_plane` reads the co-moving phase space from `mult_dipole.h5`: positions in the beam frame [m], momenta in β·γ rotated into the reference direction → x′ = px/pz, y′ = py/pz, δ = |p|/p₀ − 1.
- The reference orbit (`ref_px`/`ref_pz` in `mult_dipole.stat`) gives the steering angle.
- **Tests fed:** 1 (energy), 5 (dipole steering angle).
