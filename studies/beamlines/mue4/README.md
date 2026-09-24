# mue4 — the muE4 beamline at PSI in OPALX

| case | what it is |
|---|---|
| `mue4_WsxOn` | hand translation of the G4beamline muE4 line with the WSX solenoid on: solenoid from a 1D map, 3 SBEND, 12 QUADRUPOLE, scraping collimators, 20 monitors, 14462 muons from `shared/beam_WsxOn.dat` |
| `mue4_analytical` | the analytic twin: 200000 muons from a Gaussian, analytic magnets, scraping apertures |
| `mue4_analytical_open` | the same with every aperture opened to 20 m, 20000 muons |
| `mue4_WsxOff` | an unfinished WSX-off translation built from MULTIPOLE elements; its header says it does not run |

The line built from G4beamline's own 3D field maps, placed by absolute pose and compared
particle by particle with G4beamline, is in `../../g4bl/mue4/`. The checkpoint test on this
line is `../../features/checkpoint/mue4_ckpt`.

Known issue: `mue4_WsxOn` and `mue4_analytical_open` set `DESIGNENERGY` on the bends, which
`OpalSBend::update` rejects since OPALX PR #520; remove it to run them on a newer build.

Run one case from the repo root; the output goes to the same path under `output/`:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
export OPALX=/Users/rammann/Code/OPALX/opalx/build_serial/src/opalx
$PY -m opalxruns.run studies/beamlines/mue4/mue4_analytical/mue4_analytical.in
$PY -m opalxruns.process_run studies/beamlines/mue4/mue4_analytical      # plots and ParaView files for that run
```
