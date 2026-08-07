# opalx-runs

OPALX input decks. One directory per case, grouped by topic under `studies/`.

## Running a case

OPALX resolves `FMAPFN` / `FNAME` relative to the working directory, so always `cd` into the
case first:

```bash
OPALX=/path/to/executable

cd studies/collimator/circle
mpirun -n 1 $OPALX circle.in --info 2
```

Output (`.h5`, `.stat`, `run.log`, `timing.dat`, `data/`) is written next to the deck and is
gitignored.