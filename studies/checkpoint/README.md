# checkpoint — muE4 kill-and-restart test

End-to-end test of `OPTION, CHECKPOINTFREQ` + `opalx --restart <file>`, added on branch
`523-add-checkpointrestart-capability` (commit `aa518cbf`). The feature ships with two
file-format round-trip unit tests; nothing exercises an actual tracking run being killed and
resumed. This does.

## The contract

A restart is correct if the run it produces is indistinguishable from a run that was never
interrupted. Concretely: same number of `.stat` rows, a seam with no duplicated or skipped step,
the same final beam, and monitor files holding each particle crossing exactly once.

## The deck

[mue4_ckpt/mue4_ckpt.in](mue4_ckpt/mue4_ckpt.in) — the muE4 beamline from
[studies/mue4/mue4_WsxOn](../mue4/mue4_WsxOn/mue4_WsxOn.in), physics unchanged, with three edits:
5000 muons from `beam_WsxOn_5k.dat` instead of 14462, `OPTION, CHECKPOINTFREQ = 2000`, and
`STATDUMPFREQ = 1` kept so the `.stat` has one row per step and the seam is checkable exactly.

It was picked because it stresses what a unit test cannot: aperture scraping (5000 → 512, so the
checkpoint must restore a bunch that shrinks during tracking), a solenoid + 3 SBENDs + 12 quads,
and 20 `MONITOR` planes writing their own HDF5 files. `FROMFILE` input means no random numbers
anywhere, so an uninterrupted run and a killed-and-restarted run are directly comparable.

Run length is 25693 steps to `ZSTOP = 19.25 m`; about 25 s on one rank with a Release SERIAL
build. The killed runs are stopped after their third checkpoint, i.e. at global step 6000,
`s ≈ 4.49 m` — inside quad triplet 1, past 7 of the 20 monitors.

## Running it

```bash
./run_checkpoint_test.sh                  # all six scenarios, ~3 min
./run_checkpoint_test.sh ref r1           # or a subset
/opt/homebrew/Caskroom/miniconda/base/bin/python compare.py
```

`OPALX=/path/to/opalx` overrides the binary; `NCHECKPOINTS=N` changes how far the killed runs get.
Each scenario gets its own `work_*/` directory (gitignored) — a restart appends to the killed
run's own `.stat`/`.h5` after rewinding them, so it has to stay where it was killed.

| scenario | what it does |
|---|---|
| `work_ref` | uninterrupted reference, 1 rank |
| `work_r1` | SIGKILL mid-track, restart on 1 rank |
| `work_r2` | SIGKILL mid-track on 1 rank, restart on 2 ranks |
| `work_ctl2` | control: fresh 2-rank run, no checkpoint — is the *deck* fine on 2 ranks? |
| `work_r2b` | control: killed on 2 ranks, restarted on 2 ranks — is it the rank *change* or just multi-rank? |
| `work_nomon` | control: `r2b` with the MONITORs stripped — is the abort in the monitor path? |

`compare.py` checks resume point, `.stat` row count and monotonicity, every `.stat` column against
the reference, `.h5` dump count, final particle count, final phase space, and per-monitor crossing
counts. Differences are scaled by the column's own magnitude rather than element-wise, so a 1e-11
absolute difference on a coordinate that happens to sit near zero is not reported as catastrophic.
Particle ids are handed out per rank, so when the id sets disagree the final phase space is
compared as sorted distributions instead.

## Results

Release SERIAL build, Open MPI 5.0.9, macOS arm64.

**Tracking state restores exactly.** In all four restart scenarios the run resumes at global step
6000, produces all 25693 `.stat` rows with a strictly increasing `s` (no duplicated or skipped
step at the seam), 257 phase-space dumps, and 512 surviving particles. Every checked `.stat`
column matches the reference, and so does the final phase space:

| | `.stat` worst column | final phase space |
|---|---|---|
| `work_r1` (1 → 1 rank) | 3.7e-13 | 1.9e-11, matched by id |
| `work_r2` (1 → 2 ranks) | 2.2e-13 | 1.9e-11, matched by id |
| `work_r2b` (2 → 2 ranks) | 1.3e-08 | 5.4e-11, sorted distributions |
| `work_ctl2` (control, no restart) | 1.3e-08 | 5.4e-11, sorted distributions |

The 1e-8 figures are MPI reduction-order noise: `work_ctl2` never restarts at all and shows the
same number, so it is the cost of running on 2 ranks, not of restarting.

Two defects fall out.

### 1. Restart on more than one rank aborts at the end of the run

```
* Dump phase space of last step
*** An error occurred in MPI_Allreduce
*** MPI_ERR_TRUNCATE: message truncated
```

Exit code 15. It happens for both `1 → 2` (`work_r2`) and `2 → 2` (`work_r2b`), so it is
**multi-rank restart, not the rank change**. `work_ctl2` — the same deck, fresh, on 2 ranks —
completes normally, so the deck itself is fine on 2 ranks.

The abort lands after all tracking and after the final phase-space dump, so `.stat` and `.h5` are
complete and correct (see the table above). What is lost is every `MON_*.h5`: monitors flush at
`goOffline()`, which the abort pre-empts.

`work_nomon` repeats `work_r2b` with the 20 MONITORs stripped and exits 0, which puts the abort in
the monitor / `LossDataSink` path.

Reading that path, the mechanism is a rank-local early return guarding a collective.
`LossDataSink::splitSets` (`src/Structure/LossDataSink.cpp:1097`):

```cpp
if (numSets <= 1 || particles_m.empty()) {
    return;
}
```

`particles_m.empty()` is per-rank, and the body it skips contains
`ippl::Comm->allreduce(data.data(), 2 * numSets, ...)` at line 1135. A rank that recorded no
crossings for a monitor returns early while a rank that recorded some enters the allreduce, and
the collectives desynchronize. `LossDataSink::hasNoParticlesToDump` (line 705) does a global
reduction for exactly this reason and its comment says so; `splitSets` does not.

That also explains why a fresh run is fine and a restart is not: fresh, every monitor collects
thousands of crossings spread over both ranks, so no rank is ever empty. After a restart the
monitors near the resume point collect very few (see below — 16 crossings at `MON_06`), and those
can easily all land on one rank.

### 2. Monitor data recorded before the checkpoint is lost

Visible in `work_r1`, the clean single-rank case that otherwise passes everything:

| | reference | after restart |
|---|---|---|
| `MON_01_SOL_IN` … `MON_05_FS61_IN` | 2913–5000 crossings | file does not exist |
| `MON_06_FS61_OUT` | 2913 | 16 |
| `MON_07_Q1_IN` | 2913 | 961 |
| `MON_08_Q1_OUT` … `MON_20_Q4_OUT` | — | matches |

The checkpoint stores particle, tracking, reference and cavity state, but not `LossDataSink`
contents. Monitor crossings accumulate in memory and are only written at `goOffline()`, at the end
of the run — the SIGKILLed run left no `MON_*.h5` at all (`work_r1/killed_state/` holds only the
checkpoint). So a restart can only record planes the beam has yet to reach: everything upstream of
`s ≈ 4.49 m` is gone, and `MON_06`/`MON_07`, which straddle the resume point, keep only the
stragglers that had not yet crossed.

`DataSink::rewindToCheckpoint` handles `.stat`, `.lbal` and the phase-space `.h5` correctly; the
monitor sinks are simply not part of it.

## Notes

- `opalx --restart ck.h5 deck.in --info 2` fails with "Unknown argument": the positional input
  file is only accepted at `argv[1]` or last (`src/Main.cpp:300`). Put the deck first —
  `opalx deck.in --restart ck.h5 --info 2`.
- Only one checkpoint file exists at a time; each write atomically replaces
  `<basename>_checkpoint.h5`. Copy it aside to keep a history.
- Restart is refused outright for emitting distributions, `T0 > 0`, `EMITTEDFROMFILE`, and any
  `GLOBALPROCESSES` such as DECAY (`src/Track/TrackRun.cpp:365`). muE4 clears all of them.
- No checkpoint is written at end of run — only when `globalTrackStep % CHECKPOINTFREQ == 0`.
