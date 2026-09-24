# bend — SBEND and RBEND against the analytic result

| cases | what they check |
|---|---|
| `{sbend,rbend}_{30,90,180}deg` | one bend, the reference orbit against the analytic arc (hand-written) |
| `{sbend,rbend}_{enge,he}_{30,60,90}` | 13 particles, the transfer matrix against the analytic one; `enge` has fringe fields (`HGAP 0.01`), `he` hard edges (`HGAP 0`) |
| `{sbend,rbend}_bunch_{pure,emit,disp,mean}_{30,60}` | an 8000-particle Gaussian bunch: size, emittance, dispersion, centroid |
| `chicane` | four SBENDs (+, -, -, +), 200 electrons |

The generated cases were written by `make_decks.py` of the old bend validation suite and
checked by its `run_tests.py` (their headers and READMEs still name them). Those scripts are
not in this repo; `../multipolet/` holds their descendants. Here the expected numbers are in
each case's `README.md`.

Known issue: these inputs set `DESIGNENERGY` on the bends, which `OpalSBend::update` rejects
since OPALX PR #520; remove it to run them on a newer build.

Run one case from the repo root; the output goes to the same path under `output/`:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
$PY -m opalxruns.run studies/elements/bend/sbend_he_30/sbend_he_30.in
$PY -m opalxruns.process_run studies/elements/bend/sbend_he_30      # plots and ParaView files for that run
```
