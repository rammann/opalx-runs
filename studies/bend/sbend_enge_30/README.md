# sbend_enge_30

**SBEND** bend at **30°**, Enge fringe (HGAP=0.01 m). 0.1 GeV electron, no space charge.

## Physics

The design particle rides an arc of radius **ρ = 1.9099 m** (ρ = L/θ) turning through **θ = 30°**. The uniform dipole field is
**B₀ = (P₀/c)·k₀/q = Bρ/ρ = 0.1755 T**  (magnetic rigidity Bρ = 0.3353 T·m, k₀ = 1/ρ). A magnetic field does no work, so |p| is conserved (test 1).

First-order motion is the 6×6 map in (x, x′, y, y′, z, δ). For a **sector bend** the bending plane is

```
[ cosθ      ρ sinθ ]   [ +0.8660  +0.9549 ]
[ -sinθ/ρ   cosθ   ] = [ -0.2618  +0.8660 ]
```

the non-bending (vertical) plane is a **drift** of the arc length L = 1.0000 m, and the dispersion is x = ρ(1−cosθ)·δ = 0.2559·δ, x′ = sinθ·δ = 0.5000·δ (tests 3, 4). A sector bend gives no vertical focusing to first order (test 6).

A physically correct (symplectic) map satisfies **MᵀJM = J** and each 2×2 block has determinant 1 (test 7).

## Setup (`sbend_enge_30.in`)

- `FIELDSOLVER, TYPE = NONE` (no space charge).
- `DISTRIBUTION, TYPE = FROMFILE` reading `parts.txt`: **13 particles** — a reference at the design momentum plus a symmetric **± step in each of the six phase-space coordinates** (εx=εy=εz=0.0001 m, εx′=εy′=0.0001 rad, εδ=0.001). No `PC` on the BEAM (FROMFILE reads the momenta from the file); the field reference momentum comes from the `REAL P0` deck variable.
- Lattice: `DRIFT(L=1) → SBEND(L=1, ANGLE=30°, HGAP=0.01) → DRIFT`. Design faces at s = 1.000 .. 2.000 m; the Enge field spills a little beyond each face. Dumps are dense (`PSDUMPFREQ` small) so read planes land near the faces.

## How the output is processed (`bendlib.py` → `run_tests.py`)

- `Case.read_plane` reads the co-moving phase space from `sbend_enge_30.h5`: positions are in the beam frame [m], momenta in β·γ rotated into the reference direction → x′ = px/pz, y′ = py/pz, δ = |p|/p₀ − 1.
- `Case.transfer_map` builds the 6×6 map: pick a field-free plane on each side of the bend (the field extent comes from `By_ref` in `sbend_enge_30.stat`), project every particle onto the reference transverse plane (drift by −z), strip the exact drift from each plane to the design face, then centered-finite-difference the ± pairs.
- The reference orbit (`ref_x/ref_z/ref_px/ref_pz` in the stat) gives the bend angle (test 2) and the energy check (test 1).
- **Tests fed:** 1 (energy), 2 (bend angle from ∫B ds), 3 (transfer matrix), 4 (dispersion), 7 (symplecticity), and 6 (rbend-vs-sbend vertical focusing, the 30° pair) (and supplies R56/R16 for the bunch tests 5, 8-10).
