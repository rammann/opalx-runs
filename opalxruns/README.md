# opalxruns

Shared Python for the studies in this repo: plots and ParaView geometry from an
OPALX run directory, readers for OPALX and G4beamline files, and the helpers that
more than one study uses.

Install once, into the miniconda Python — it has `numpy`, `h5py`, `matplotlib`,
`pandas`, `scipy` and `vtk`; the system `python3` has none of them:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
$PY -m pip install -e /Users/rammann/Code/OPALX/opalx-runs --no-deps
```

## Plots and ParaView geometry for one run

One entry point, `process_run`, plus the individual scripts it calls. Elements
are drawn as hollow pipes with the element's real OPALX aperture as the bore
(see below).

```bash
$PY -m opalxruns.process_run studies/ring/square_ring             # everything
$PY -m opalxruns.process_run studies/ring/square_ring --dry-run   # what would run
```

## What a run directory needs

```
<run>/
  <base>.in                             deck
  <base>.stat                           SDDS statistics       -> stat plots
  <base>.h5                             bunch dumps           -> ParaView particles
  MON_*.h5                              monitor planes        -> monitor plots
  timing.dat                            timer table           -> timing plot
  data/<base>_DesignPath.dat            reference orbit       -> lab frame, ref orbit
      /<base>_ElementPositions.txt      element lab positions -> ParaView elements
```

Nothing is required. Each step whose inputs are missing is skipped with a note,
so a run with no monitors still gets everything else. Outputs go to
`<run>/plots/` and `<run>/paraview/`.

If a directory holds the output of **more than one deck** (as `mue4_analytical`
does), the basename of `data/` decides which run is processed — that is the one
OPALX wrote last. The others are named in the inventory and reachable with
`--base`.

## Scripts

| file | what it does |
|---|---|
| `process_run.py` | the run-all. Reports what it found, runs each step, keeps going if one fails |
| `plot_stat.py` | `.stat` overview: transmission, energy, RMS size, emittance, centroid, reference orbit, dispersion. `--columns a,b,c` plots any named columns |
| `plot_timing.py` | `timing.dat` as sorted bars with each timer's share of `mainTimer` |
| `plot_monitors.py` | per-monitor x–x', y–y', x–y phase space, plus `envelope.png` and `stats.csv` over all planes |
| `particles_to_vtk.py` | bunch dumps to `.pvd` + `.vtp`, in lab and/or co-moving coordinates |
| `elements_to_vtk.py` | element bodies and the reference orbit as `.vtp`, in lab coordinates. An element with an aperture becomes a hollow pipe (bore = the real aperture per OPALX rules: `APERTURE` string, or `HGAP`/`HAPERT` rectangle for bends, `--wall` thick); one without (hard-edge bend, bare element) a solid tube of `--default-aperture` radius. Cell arrays `element_type` and `element_id` (per-element index, names printed at run time) |
| `opalx_run.py` | shared: the `Run` file-discovery class, `timing.dat` / `DesignPath` / `ElementPositions` readers, co-moving→lab frame math |
| `opalx_diagnostics.py` | `.stat` and `.h5` readers, plus the notebook dropdown widgets. Imported by `fmlib.py`, `bendlib.py`, `cmplib.py` and `plot_layout.py` — keep its names stable |
| `paths.py` | where things are: the repo, `G4BL_FILES`, `G4BL_APP`, and `opalx()` from `$OPALX` |
| `h5.py` | `read_monitor`: every `Step#` group of an OPALX monitor file, sorted by id |
| `g4bl.py` | G4beamline files: writers for the `grid` and `cylinder` field map formats, `read_map_header`, `read_track_file` and `track_on_plane` for #BLTrackFile output |
| `case.py` | `Case`: one case folder of a scripted study — `.stat` table, reference orbit, dumps with their path lengths, `read_plane`, `trajectory`. `fmlib.py` and `bendlib.py` subclass it |
| `matrices.py` | `drift_matrix`, `symplectic_residual`, and `transfer_matrix` from tracked +/- particle pairs |
| `particles.py` | the 13-particle set (`map_particles`) and the FROMFILE particle file (`write_parts`) |
| `results.py` | `Results`, the PASS/FAIL table every scripted study writes to `results.txt` |
| `plotstyle.py` | the OPALX / G4BL colours, the shared matplotlib settings (`use`), `save_with_footer` |
| `mue4.py` | where every element of the muE4 G4beamline input ends up (the centreline walk) |
| `mue4_beam.py` | the muE4 reference muon (28 MeV/c) and the functions that assume it: G4beamline track files in beta*gamma, the six comparison coordinates, matching the two codes by particle id |

Every script also runs on its own, taking a run directory:

```bash
$PY -m opalxruns.plot_stat        studies/mue4/mue4_analytical --columns rms_x,rms_y,Dx
$PY -m opalxruns.plot_monitors    studies/mue4/mue4_analytical
$PY -m opalxruns.particles_to_vtk studies/mue4/mue4_analytical --frame lab --stride 5
$PY -m opalxruns.elements_to_vtk  studies/mue4/mue4_analytical --default-aperture 0.05
```

## Two things worth knowing

**Particles are dumped in the co-moving frame.** OPALX tracks with the origin on
the reference particle and local +z along the reference momentum, so the `x,y,z`
in the `.h5` are not lab coordinates. `particles_to_vtk.py` writes both:
`bunch_comoving.pvd` is the raw data, `bunch_lab.pvd` applies

```
R_lab = RefPartR(s) + Q(s) . (x,y,z)_local
P_lab =               Q(s) . (px,py,pz)_local
```

`RefPartR` is stored per step; the rotation `Q(s)` is not — `TaitBryantAngles` is
written but stubbed to zero. It is rebuilt by parallel transport of the
reference-momentum direction from `DesignPath.dat`, which reproduces OPALX's own
`toLabTrafo` including the accumulated frame roll. Without a `DesignPath` the
script warns and falls back to a single shortest-arc rotation, which is only
right for narrow, near-planar beams. Polarization, `E` and `B` are vectors in the
same frame and are rotated with the momentum.

**A monitor file holds several `Step#N` groups for one pass.**
`LossDataSink::splitSets` partitions the recorded batch by crossing time, so the
sets are consecutive slices of the same crossing, not repeats. Reading only
`Step#0` halves the answer — in `mue4_analytical` the first plane's two sets hold
100907 and 95972 particles with disjoint ids, and only their union (196881 of
200000 launched) is the real count. `plot_monitors.py` concatenates every set and
reports duplicate ids rather than hiding them.

## ParaView

1. **File → Open** `bunch_lab.pvd`, `<base>_elements.vtp`, `<base>_reforbit.vtp`;
   **Apply** each.
2. Elements: colour by `element_type` (`0=DRIFT 1=SOLENOID 2=DIPOLE
   3=QUADRUPOLE 4=MONITOR 5=OTHER`), opacity ≈ 0.3 to see the particles inside.
3. Particles: colour by `Ekin_MeV` or `Pmag`; a **Glyph** filter on the vector
   `P` gives momentum arrows. Spin runs also carry `Pol`, `E` and `B`.
4. Press play — the `.pvd` animates over the dump times.

`--stride N` keeps the output down on long runs. The binary `.vtp` costs about
87 bytes per particle per dump, so `mue4_analytical` (251 dumps, 200k particles
falling to 9k) writes ~635 MB per frame, ~1.3 GB for both.

## Notebooks

Interactive and left as they are. They import `opalx_diagnostics`, so run Jupyter
from this directory.

| notebook | what |
|---|---|
| `PlotPolarization.ipynb` | polarization columns from `.stat` / `.h5`, via the `StatSelector` / `H5Selector` dropdowns |
| `SpinBenchmark.ipynb` | the spin pusher against analytic Thomas–BMT precession: longitudinal **B**, transverse initial spin |
| `PlotDecay.ipynb` | `numParticles` against the analytic exponential muon decay |
| `PlotDecayDaughters.ipynb` | multi-container viewer: `numParticles(t)` per `*_cN.stat`, plus their sum as a conservation check |
| `PlotFM.ipynb` | field maps from `.T7` files |

All five still point at an **older directory layout** — `../output`,
`inputfiles/`, `opalx_runs/fieldmaps` — none of which exist under `runs/` any
more. Each has its input directory in a constant near the top (`STAT_DIR`,
`GLOB_PATTERN`, or an `output_root` argument); repoint that at the run directory
you want before running one.
