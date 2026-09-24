# collimator — COLLIMATOR removes exactly the particles outside its aperture

A grid of particles (`parts.txt`) goes through one COLLIMATOR; the particles outside the
aperture must be deleted and the rest must pass untouched.

| case | aperture | expected |
|---|---|---|
| `circle` | `CIRCLE(0.02)`, 57 particles | 32 deleted, 25 survive |
| `circle_np2` | the same on 2 MPI ranks (`--np 2`) | the same |
| `ellipse` | `ELLIPSE(0.03, 0.02)`, 49 particles | 24 deleted, 25 survive |
| `rectangle` | `RECTANGLE(0.02, 0.04)` | 20 deleted, 29 survive |
| `square` | `SQUARE(0.02)` | 24 deleted, 25 survive |
| `noop_flag` | `CIRCLE(0.02)` with `DELETEONTRANSVERSEEXIT = FALSE` | nothing deleted |

The cases were written by `make_decks.py` of the old collimator validation suite and
checked by its `run_tests.py`; neither is in this repo. Each case's `README.md` has the
expected result.

Run one case from the repo root; the output goes to the same path under `output/`:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
export OPALX=/Users/rammann/Code/OPALX/opalx/build_serial/src/opalx
$PY -m opalxruns.run studies/elements/collimator/circle/circle.in
$PY -m opalxruns.process_run studies/elements/collimator/circle      # plots and ParaView files for that run
```
