# muE4 strong scaling on daint

The whole muE4 line from `../full/` (22 field maps, nothing scrapes) with 16e6 muons, run
on 1, 2, 4 and 8 GH200 GPUs on daint, one rank per GPU. Every run uses the same input, so
the run times show how well one problem gets faster with more GPUs.

## What is different from `full/`

| | `full/` | here |
|---|---|---|
| particles | 10 000 | 16 000 000: `../full/parts.txt` repeated 1600 times |
| `DT` | 1e-12 s, about 250 000 steps | 1e-11 s, about 25 000 steps |
| MONITOR planes | 22 | none (each would write 16e6 particles to an `.h5` file) |
| phase space dumps | at the end | none (`PSDUMPFREQ = 0`) |
| statistics | every 10 steps | every 1000 steps |
| map paths | absolute, on the laptop | `maps/<file>`, relative to the run folder |

The particles repeat, but there is no space charge, so the copies do not interact and
each one costs the same work as in `full/`. Repeating keeps the bunch mean, which sets
OPALX's reference particle, exactly as in `full/`. 16e6 particles are 2e6 per GPU at 8
GPUs. With much fewer particles per GPU, the fixed cost of each step dominates and
adding GPUs stops helping. To change the count, set `N_PART` in `make_scaling.py`.

## Files

| file | what | where it runs |
|---|---|---|
| `make_scaling.py` | writes `scaling.in` and `lattice_scaling.in` from `../full/lattice_full.in` | laptop |
| `scaling.in`, `lattice_scaling.in` | generated, the OPALX input | |
| `copy_maps_to_daint.sh` | copies the 6 field maps (about 165 MB, not in git) into `maps/` of the clone on daint | laptop |
| `make_beam.sh` | writes `parts_16M.txt` (1.5 GB) from `../full/parts.txt` | daint |
| `submit.sh` | one Slurm job per GPU count, each in its own run folder | daint |
| `job.sh` | the job itself: `srun opalx scaling.in --info 2` | daint (compute node) |
| `collect.sh` | reads `timing.dat` of each run, prints a table with speedup | daint or laptop |

## Running

**1. daint, once: build OPALX for the GPUs.** The branch `general-fieldmap-element` has
the `FIELDMAP` element.

```bash
uenv start --view=default /capstor/store/cscs/cscs/public/uenvs/opal-x-gh200-mpich-gcc-2026-10-01.squashfs
cd $SCRATCH
git clone -b general-fieldmap-element git@github.com:OPALX-project/OPALX.git
mkdir build-cuda && cd build-cuda
cmake ../OPALX -DCMAKE_BUILD_TYPE=Release -DPLATFORMS=CUDA -DARCH=HOPPER90
make -j 32 opalx_exe         # not `make opalx`, which builds only the library
exit                         # leave the uenv
```

`submit.sh` expects the executable at `$SCRATCH/build-cuda/src/opalx`; set `OPALX` for
another one.

**2. daint, once: clone this repo.** Commit and push this folder first. Cloning from
daint needs a GitHub SSH key on daint.

```bash
cd $SCRATCH
git clone git@github.com:rammann/opalx-runs.git
```

Later changes come with `git pull` in `$SCRATCH/opalx-runs`.

**3. Laptop, once: copy the field maps.** They are not in git.

```bash
cd studies/g4bl/mue4/daint_scaling
./copy_maps_to_daint.sh      # DAINT=<ssh host>, default daint
```

**4. daint, once: write the particle file.**

```bash
cd $SCRATCH/opalx-runs/studies/g4bl/mue4/daint_scaling
./make_beam.sh
```

To change the particle count: set `N_PART` in `make_scaling.py`, run `$PY make_scaling.py`
on the laptop, commit, push, then `git pull` and `./make_beam.sh` on daint.

**5. daint: a short test on 1 GPU** (the first metre of the line, a few minutes):

```bash
GPUS=1 ZSTOP=1.0 ./submit.sh
```

Check in `slurm_<id>.out` that the run ends without an error and prints
`Rank 0: 16000000 local particles`, and that `timing.dat` was written. Multiply the
`steps` time from `./collect.sh` by 19.4 to estimate the 1-GPU run of the whole line;
if that is near the 30 min limit, set `TIME=01:00:00` for step 6.

**6. daint: the scaling runs.**

```bash
./submit.sh                  # 1, 2, 4, 8 GPUs; ACCOUNT=c41, TIME=00:30:00 by default
squeue --me
./collect.sh                 # once all jobs are done
```

Each call of `submit.sh` makes a new folder
`$SCRATCH/opalx-runs/output/g4bl/mue4/daint_scaling/runs_<date>_<time>/gpus_<n>/`, so
calling it twice gives a second set of runs to compare with the first.

## Reading the table

`collect.sh` prints the `Wall max` of `timing.dat` (the slowest rank), in seconds:

| column | timers |
|---|---|
| `main` | `mainTimer`: the whole run |
| `steps` | everything that runs every time step: `TIntegration1` + `TIntegration2` + `External field eval` + `computeMoments` + `updateParticle` |
| `field` | `External field eval` |
| `push` | `TIntegration1` + `TIntegration2` |
| `moments` | `computeMoments` |
| `speedup`, `eff` | $S(n) = T_\mathrm{steps}(1) / T_\mathrm{steps}(n)$, $E(n) = S(n)/n$ |

Use `steps` for the scaling, not `main`. `main` also counts reading the maps and the
particle file, and every rank reads the whole particle file (`FromFile::readFile`), so
that part does not get faster with more GPUs.

`computeMoments` runs every step even though the field solver is `TYPE = NONE`:
`ParallelTracker::computeSpaceChargeFields()` calls `calcBeamParameters()` before it asks
the solver for anything. In the local check below it was three quarters of the step time.

## Local check

The whole chain (`make_beam.sh`, `submit.sh`, `job.sh`, `collect.sh`) was run on the
laptop on 2026-10-02 with stand-ins for `sbatch` and `srun` that start `mpirun`, with
20 000 particles, `ZSTOP = 3.0`, SERIAL build, 1 and 2 ranks:

```
 GPUs      main     steps     field      push   moments  speedup    eff
    1      16.9      10.7       1.7       0.8       8.2     1.00   1.00
    2      10.7       6.2       0.9       0.5       4.5     1.72   0.86
```

These are CPU numbers for a tiny bunch. They only show that the scripts work and that the
input runs; the GPU numbers will be different.
