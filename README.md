# opalx-runs

OPALX input files and the checks built on them, grouped by what each check is compared
against. Most were written while adding or fixing an OPALX feature on a branch.

## Layout

```
opalx-runs/
  studies/
    elements/    one element checked against an analytic result
    features/    code features: checkpoint, decay, spin, multibeam, ring
    g4bl/        OPALX against G4beamline
    beamlines/   whole machines
  shared/        input files used by several topics (beams, field maps)
  output/        run output, not tracked; the same folder tree as studies/
  opalxruns/     shared Python: running cases, reading output, plots (see its README)
  tests/         unit tests for opalxruns
```

## Setup, once

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python   # has numpy, h5py, matplotlib, scipy, vtk
$PY -m pip install -e . --no-deps
```

| variable | what | default |
|---|---|---|
| `OPALX` | the opalx executable | none: required to run OPALX |
| `G4BL_APP` | the G4beamline app | `/Users/rammann/Code/G4BL/G4beamline-3.08.app` |
| `G4BL_FILES` | the muE4 G4beamline input and field maps (889 MB, not in this repo) | `../g4bl-files` |
| `PY` | the Python the `run_all.sh` scripts use | the miniconda one above |

## Running

One case, from the repo root:

```bash
export OPALX=/Users/rammann/Code/OPALX/opalx/build_serial/src/opalx
$PY -m opalxruns.run studies/elements/collimator/circle/circle.in          # --np 2 for 2 ranks
$PY -m opalxruns.process_run studies/elements/collimator/circle           # plots, ParaView files
```

The run happens in `output/elements/collimator/circle/`. The runner clears that folder's
own output (folders of other runs inside it stay), links in the input and every file the
input names, and starts OPALX there, so all of the output lands there. `--keep -- --restart <file>` resumes from a checkpoint. G4beamline
inputs (`.g4bl`) run the same way, and
`$PY -m opalxruns.run <case>.g4bl <case>.in` runs both codes, one after the other, in one folder.

Topics with scripts have a `run_all.sh`:

```bash
cd studies/elements/fieldmap
./run_all.sh               # write the cases, run them, check -> results.txt
./run_all.sh --test-only   # check the existing output again
```

`$PY -m opalxruns.refs` lists files that inputs name but that do not exist. Twelve are
expected: the `g4bl/elements/wsx_solenoid*` notebooks write their particle files
themselves, and the `g4bl/spin/uniform_bz/scan/` inputs read theirs from the folder above.
`$PY -m unittest discover -s tests` runs the unit tests.

## Topics

| topic | what it checks | how to run | result |
|---|---|---|---|
| `elements/bend` | SBEND, RBEND: orbit, transfer matrix, bunch moments against analytic results | one case | in each case's README |
| `elements/multipole` | MULTIPOLE and QUADRUPOLE transfer matrices against analytic ones | one case | in each case's README |
| `elements/multipolet` | MULTIPOLET as a bend, 11 tests against analytic results | `run_all.sh` | `results.txt`, 225/225 |
| `elements/aperture` | per-element `APERTURE` and bend `HGAP`/`HAPERT`: survivors against the Gaussian integral | one case | README table |
| `elements/collimator` | COLLIMATOR deletes exactly the particles outside its aperture | one case | README table |
| `elements/fieldmap` | G4beamline `grid` and `cylinder` map formats against analytic fields | `run_all.sh` | `results.txt`, 68/68 |
| `features/checkpoint` | stopping and restarting from a checkpoint | one case, then a restart | README |
| `features/decay` | muon and pion decay, daughter beams, polarization | one case | none written |
| `features/spin` | the Thomas-BMT spin integrator against analytic rotations | one case | none written |
| `features/multibeam` | two beams in one run | one case | none written |
| `features/ring` | several turns through the same elements (square ring, ISIS) | one case | README |
| `g4bl/elements` | one field map at a time, OPALX against G4beamline | `run_all.sh` | `results.txt`, 379/383 |
| `g4bl/mue4` | the whole muE4 line from G4beamline's maps, against G4beamline | README steps, `full/run_full.sh` | `full/full_results.txt` |
| `g4bl/spin` | spin precession in both codes against the closed form | `run_all.sh` | `results.txt`, `bz_results.txt` |
| `g4bl/report` | one HTML page from `g4bl/elements`, `g4bl/mue4` and `g4bl/spin` | `report/build_report.py` | `output/g4bl/report/index.html` |
| `beamlines/mue4` | the muE4 line in OPALX: hand translation and analytic twin | one case | none written |

## Adding a study

1. Put it in the group that says what it is compared against: an analytic result goes to
   `elements/`, a code feature to `features/`, G4beamline to `g4bl/`, a whole machine to
   `beamlines/`.
2. One folder per case, `studies/<group>/<topic>/<case>/<case>.in`, names in lower case with
   underscores.
3. Inputs name other files only through `FMAPFN`, `FNAME`, `CALL` (OPALX) or `file=`,
   `filename=` (G4beamline). Those files live in the case folder, its topic folder (`../`),
   `shared/`, or `g4bl-files/` (absolute paths written by a generator).
   `python -m opalxruns.refs` checks them.
4. Run through `opalxruns.run`, so output never lands in `studies/`.
5. A topic with scripts has `make_*.py` (writes the case folders and `cases.json`),
   `run_all.sh [--test-only]`, `run_tests.py` (writes `results.txt`), `plot_tests.py` and a
   `README.md`. Code that more than one topic needs goes into `opalxruns/`.
6. Track input files, scripts, `results.txt` and `*_data.json`. Do not track run output,
   except G4beamline results that take long to redo, kept in a `g4bl_reference/` folder.
