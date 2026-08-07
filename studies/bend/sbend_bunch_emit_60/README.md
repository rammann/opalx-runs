# sbend_bunch_emit_60

**SBEND** bend at **60°**, Enge fringe. Gaussian-bunch case for **test 8 -- emittance invariance**. No space charge; linear, no radiation.

## Physics

The design particle rides an arc of radius **ρ = 0.9549 m** (ρ = L/θ) turning through **θ = 60°**. The uniform dipole field is
**B₀ = (P₀/c)·k₀/q = Bρ/ρ = 0.3511 T**  (magnetic rigidity Bρ = 0.3353 T·m, k₀ = 1/ρ). A magnetic field does no work, so |p| is conserved (test 1).

The bunch covariance transforms as **Σ_out = M·Σ_in·Mᵀ** with the same map M measured in the map case `sbend_enge_60`. This case is the special projection:

With **zero energy spread** the transverse motion is decoupled from energy, so the projected RMS emittances εx = √(⟨x²⟩⟨x′²⟩−⟨xx′⟩²) and εy are **invariant** under the symplectic map (test 8). The 6D emittance √(det Σ) is invariant in all planes. Zero energy spread is used so dispersion does not inflate the projected emittance.

## Setup (`sbend_bunch_emit_60.in`)

- `FIELDSOLVER, TYPE = NONE` (no space charge).
- `DISTRIBUTION, TYPE = GAUSS`, N = 8000 particles. Spreads: σx=σy=0.001 m, σx′=σy′=0.0002 rad, σz=0.0005 m, zero energy spread.
- Same lattice as `sbend_enge_60`: `DRIFT → SBEND 60° (HGAP=0.01) → DRIFT`.

## How the output is processed (`bendlib.py` → `run_tests.py`)

- `Case.sigma_at_faces` reads Σ at the design faces (field-free planes, each particle drift-projected onto the face).
- Compute εx, εy and the 6D emittance at entrance vs exit; they must match (in/out ratio → 1).
