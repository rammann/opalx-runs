# rbend_bunch_disp_30

**RBEND** bend at **30°**, Enge fringe. Gaussian-bunch case for **test 9 -- RMS x growth from energy spread**. No space charge; linear, no radiation.

## Physics

The design particle rides an arc of radius **ρ = 1.9319 m** (ρ = L/(2 sin(θ/2))) turning through **θ = 30°**. The uniform dipole field is
**B₀ = (P₀/c)·k₀/q = Bρ/ρ = 0.1735 T**  (magnetic rigidity Bρ = 0.3353 T·m, k₀ = 1/ρ). A magnetic field does no work, so |p| is conserved (test 1).

The bunch covariance transforms as **Σ_out = M·Σ_in·Mᵀ** with the same map M measured in the map case `rbend_enge_30`. This case is the special projection:

With an energy spread σδ, dispersion adds to the horizontal size: **σx²_out = (betatron part) + (D·σδ)²** (test 9), where D = R16 is the dispersion from the map. The full statement is Σ_out = M·Σ_in·Mᵀ.

## Setup (`rbend_bunch_disp_30.in`)

- `FIELDSOLVER, TYPE = NONE` (no space charge).
- `DISTRIBUTION, TYPE = GAUSS`, N = 8000 particles. Spreads: σx=σy=0.001 m, σx′=σy′=0.0002 rad, σz=0.0005 m, energy spread σδ=0.001.
- Same lattice as `rbend_enge_30`: `DRIFT → RBEND 30° (HGAP=0.01) → DRIFT`.

## How the output is processed (`bendlib.py` → `run_tests.py`)

- `Case.sigma_at_faces` gives Σ_in and Σ_out.
- Check σx_out against (M·Σ_in·Mᵀ)₁₁ (exact) and against the decomposition √(betatron + (D·σδ)²) using D from the matching map case.
