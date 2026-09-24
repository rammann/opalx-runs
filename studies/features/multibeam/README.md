# multibeam — two beams in one run

`fodo`: a 10-cell FODO lattice built from MULTIPOLE quadrupoles, tracking an electron beam
and a proton beam together (`BEAMS = {BEAM1, BEAM2}`).

Run one case from the repo root; the output goes to the same path under `output/`:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
export OPALX=/Users/rammann/Code/OPALX/opalx/build_serial/src/opalx
$PY -m opalxruns.run studies/features/multibeam/fodo/fodo.in
$PY -m opalxruns.process_run studies/features/multibeam/fodo      # plots and ParaView files for that run
```
