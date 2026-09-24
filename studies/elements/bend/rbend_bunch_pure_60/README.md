# rbend_bunch_pure_60

**RBEND** bend at **60°**, Enge fringe. Gaussian-bunch case for **test 5 -- RMS z growth from R56**. No space charge; linear, no radiation.

## Physics

The design particle rides an arc of radius **ρ = 1.0000 m** (ρ = L/(2 sin(θ/2))) turning through **θ = 60°**. The uniform dipole field is
**B₀ = (P₀/c)·k₀/q = Bρ/ρ = 0.3353 T**  (magnetic rigidity Bρ = 0.3353 T·m, k₀ = 1/ρ). A magnetic field does no work, so |p| is conserved (test 1).

The bunch covariance transforms as **Σ_out = M·Σ_in·Mᵀ** with the same map M measured in the map case `rbend_enge_60`. This case is the special projection:

With **energy spread only** (negligible transverse size) the growth is clean: **σx → |D|·σδ** and **σz → |R56|·σδ** (test 5). R56 is the path-length-vs-energy element of the map.

## Setup (`rbend_bunch_pure_60.in`)

- `FIELDSOLVER, TYPE = NONE` (no space charge).
- `DISTRIBUTION, TYPE = GAUSS`, N = 8000 particles. Spreads: σx=σy=σz≈1e-6 m (negligible transverse), energy spread σδ=0.001.
- Same lattice as `rbend_enge_60`: `DRIFT → RBEND 60° (HGAP=0.01) → DRIFT`.

## How the output is processed (`bendlib.py` → `run_tests.py`)

- `Case.sigma_at_faces` gives Σ_in, Σ_out; σδ = √Σ_in[δ,δ].
- Take R56 and D from the matching map case and check σz_out = |R56|·σδ and σx_out = |D|·σδ.
