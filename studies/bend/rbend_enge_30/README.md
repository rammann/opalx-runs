# rbend_enge_30

**RBEND** bend at **30°**, Enge fringe (HGAP=0.01 m). 0.1 GeV electron, no space charge.

## Physics

The design particle rides an arc of radius **ρ = 1.9319 m** (ρ = L/(2 sin(θ/2))) turning through **θ = 30°**. The uniform dipole field is
**B₀ = (P₀/c)·k₀/q = Bρ/ρ = 0.1735 T**  (magnetic rigidity Bρ = 0.3353 T·m, k₀ = 1/ρ). A magnetic field does no work, so |p| is conserved (test 1).

First-order motion is the 6×6 map in (x, x′, y, y′, z, δ). A **rectangular bend** meets each pole face at the edge angle θ/2 = 15°. The two edges defocus horizontally by just enough to cancel the body focusing, so the **bending plane reduces to a drift** (R11=R22=1, R21=0, R12=ρ sinθ = 0.9659 m), while the same edges **focus vertically** (1/f = tan(θ/2)/ρ per face, test 6). Dispersion x = ρ(1−cosθ)·δ = 0.2588·δ (test 4).

A physically correct (symplectic) map satisfies **MᵀJM = J** and each 2×2 block has determinant 1 (test 7).

## Setup (`rbend_enge_30.in`)

- `FIELDSOLVER, TYPE = NONE` (no space charge).
- `DISTRIBUTION, TYPE = FROMFILE` reading `parts.txt`: **13 particles** — a reference at the design momentum plus a symmetric **± step in each of the six phase-space coordinates** (εx=εy=εz=0.0001 m, εx′=εy′=0.0001 rad, εδ=0.001). No `PC` on the BEAM (FROMFILE reads the momenta from the file); the field reference momentum comes from the `REAL P0` deck variable.
- Lattice: `DRIFT(L=1) → RBEND(L=1, ANGLE=30°, HGAP=0.01) → DRIFT`. Design faces at s = 1.000 .. 2.012 m; the Enge field spills a little beyond each face. Dumps are dense (`PSDUMPFREQ` small) so read planes land near the faces.

## How the output is processed (`bendlib.py` → `run_tests.py`)

- `Case.read_plane` reads the co-moving phase space from `rbend_enge_30.h5`: positions are in the beam frame [m], momenta in β·γ rotated into the reference direction → x′ = px/pz, y′ = py/pz, δ = |p|/p₀ − 1.
- `Case.transfer_map` builds the 6×6 map: pick a field-free plane on each side of the bend (the field extent comes from `By_ref` in `rbend_enge_30.stat`), project every particle onto the reference transverse plane (drift by −z), strip the exact drift from each plane to the design face, then centered-finite-difference the ± pairs.
- The reference orbit (`ref_x/ref_z/ref_px/ref_pz` in the stat) gives the bend angle (test 2) and the energy check (test 1).
- **Tests fed:** 1 (energy), 2 (bend angle from ∫B ds), 3 (transfer matrix), 4 (dispersion), 7 (symplecticity), and 6 (rbend-vs-sbend vertical focusing, the 30° pair) (and supplies R56/R16 for the bunch tests 5, 8-10).
