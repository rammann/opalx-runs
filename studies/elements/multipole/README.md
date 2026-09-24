# multipole — MULTIPOLE and QUADRUPOLE against the analytic transfer matrix

13 particles (`parts.txt`) through one element; the tracked transfer matrix is compared
with the analytic one. The MULTIPOLE cases must equal the matching QUADRUPOLE case.

| case | element |
|---|---|
| `quad_defoc_x`, `quad_foc_x`, `quad_strong` | QUADRUPOLE, `K1 = +2`, `-2`, `+8` |
| `quad_skew` | QUADRUPOLE, `K1S = +2` |
| `mult_quad` | MULTIPOLE `KN = {0, 2}`, must equal `quad_defoc_x` |
| `mult_skew` | MULTIPOLE `KS = {0, 2}`, must equal `quad_skew` |
| `mult_dipole` | MULTIPOLE `KN = {0.05}`, dipole steering |

Two things to know: the sign of `K1` follows the particle's charge, and a MULTIPOLE dipole
term has half the strength one might expect. The cases were written by `make_decks.py` of
the old multipole validation suite (with `multlib.py` and `run_tests.py`); none of those is
in this repo. Each case's `README.md` has the expected result.

Run one case from the repo root; the output goes to the same path under `output/`:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
$PY -m opalxruns.run studies/elements/multipole/quad_foc_x/quad_foc_x.in
$PY -m opalxruns.process_run studies/elements/multipole/quad_foc_x      # plots and ParaView files for that run
```
