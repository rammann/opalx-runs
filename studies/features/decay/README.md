# decay — muon and pion decay, and the spin they hand on

| case | what it is |
|---|---|
| `muon_decay` | 10^6 muons with `DECAY` on, along a long drift: the count against the exponential decay law |
| `muon_decay_with_electrons` | muon decay with the daughter electrons tracked as a second beam (`DAUGHTERBEAM`) |
| `pion_muon_electron_decay` | the pi -> mu -> e chain in three containers, with the polarization output |
| `polarized_muons_no_decay` | 10000 muons with polarization (0, 0, 1) in zero field, no decay: the spin baseline |

The notebooks that plotted these (PlotDecay, PlotDecayDaughters) were never committed to
this repo. Some numbers in the comments of `muon_decay_with_electrons.in` are out of date
(drift length, gamma, step count); the input's own values are the ones that run.

Run one case from the repo root; the output goes to the same path under `output/`:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
export OPALX=/Users/rammann/Code/OPALX/opalx/build_serial/src/opalx
$PY -m opalxruns.run studies/features/decay/muon_decay/muon_decay.in
$PY -m opalxruns.process_run studies/features/decay/muon_decay      # plots and ParaView files for that run
```
