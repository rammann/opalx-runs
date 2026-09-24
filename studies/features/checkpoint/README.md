# checkpoint — stop a run, restart it from its checkpoint file

| case | what it is |
|---|---|
| `mue4_ckpt` | the muE4 line (`beamlines/mue4/mue4_WsxOn` with 5000 muons from `shared/beam_WsxOn_5k.dat`), `CHECKPOINTFREQ = 2000`: kill and restart end to end |
| `repro_drift` | a 20 m drift with a two-segment `TRACK`: the smallest input that shows a restart resuming in the wrong step-size segment |

What was found (`repro_drift.in`'s header describes the segment case): a one-rank, one-segment restart
reproduces the uninterrupted run exactly; on several ranks the restart aborts in
`LossDataSink`; monitor data written before the checkpoint is lost; with several
`DT`/`ZSTOP` segments the restart resumes in the wrong one.

Run and restart from the repo root; the checkpoint file stays in the run folder:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
export OPALX=/Users/rammann/Code/OPALX/opalx/build_serial/src/opalx
$PY -m opalxruns.run studies/features/checkpoint/repro_drift/repro_drift.in
$PY -m opalxruns.run studies/features/checkpoint/repro_drift/repro_drift.in \
    --keep -- --restart repro_drift_checkpoint.h5
```
