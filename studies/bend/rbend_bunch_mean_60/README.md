# rbend_bunch_mean_60

**RBEND** bend at **60°**, Enge fringe. Gaussian-bunch case for **test 10 -- centroid shift from mean energy**. No space charge; linear, no radiation.

## Physics

The design particle rides an arc of radius **ρ = 1.0000 m** (ρ = L/(2 sin(θ/2))) turning through **θ = 60°**. The uniform dipole field is
**B₀ = (P₀/c)·k₀/q = Bρ/ρ = 0.3353 T**  (magnetic rigidity Bρ = 0.3353 T·m, k₀ = 1/ρ). A magnetic field does no work, so |p| is conserved (test 1).

The bunch covariance transforms as **Σ_out = M·Σ_in·Mᵀ** with the same map M measured in the map case `rbend_enge_60`. This case is the special projection:

The whole bunch is given a **mean energy offset ⟨δ⟩** (BEAM `PC = P0·(1+⟨δ⟩)` while the field reference stays P0), so it is uniformly off-energy against the on-energy design orbit and its centroid disperses: **⟨x⟩_out = D·⟨δ⟩** (test 10).

## Setup (`rbend_bunch_mean_60.in`)

- `FIELDSOLVER, TYPE = NONE` (no space charge).
- `DISTRIBUTION, TYPE = GAUSS`, N = 8000 particles. Spreads: σx=σy=0.001 m, σx′=σy′=0.0002 rad, σz=0.0005 m, zero energy spread.
- Mean energy offset ⟨δ⟩=0.002 applied via `PC = P0·(1+⟨δ⟩)` (field reference stays P0).
- Same lattice as `rbend_enge_60`: `DRIFT → RBEND 60° (HGAP=0.01) → DRIFT`.

## How the output is processed (`bendlib.py` → `run_tests.py`)

- OPALX's reference particle rides *with* its own (off-energy) bunch, so the shift shows up in the **reference orbit**, not in mean_x. `centroid_orbit_shift` compares this case's reference orbit to the on-design `rbend_bunch_emit_60` orbit and projects the separation onto the design perpendicular.
- Check that separation at the exit face equals D·⟨δ⟩.
