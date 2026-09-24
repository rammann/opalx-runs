# spin — the Thomas-BMT spin integrator against analytic rotations

One muon each; all three read the field map `shared/BF_Uniform_Large.T7`.

| case | field and start | what is checked |
|---|---|---|
| `spin_bench_case1` | uniform longitudinal B = 1 T, spin transverse | the spin rotation against the analytic one |
| `spin_bench_case2` | uniform Bz, transverse launch | the angle between spin and momentum |
| `spin_bench_case3` | crossed E and B (SOLENOID plus CONSTANTEFIELDCAVITY) at the drift velocity | the spin in crossed fields |

Each input's header gives the numbers needed to rebuild the precession frequency. The
analysis notebook they mention (SpinBenchmark) was never committed to this repo. The spin
comparisons against G4beamline are in `../../g4bl/spin/`.

Run one case from the repo root; the output goes to the same path under `output/`:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
$PY -m opalxruns.run studies/features/spin/spin_bench_case1/spin_bench_case1.in
$PY -m opalxruns.process_run studies/features/spin/spin_bench_case1      # plots and ParaView files for that run
```
