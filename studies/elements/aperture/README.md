# aperture — per-element apertures cut the beam where they should

Each case sends a Gaussian beam of 10000 electrons through one element with an `APERTURE`
(or, for the bends, `HGAP`/`HAPERT`) and counts what survives. The expected counts come from
the Gaussian integral over the aperture shape; `GAUSS` cuts each axis at 3 sigma.

| case | element | expected survivors of 10000 |
|---|---|---|
| `drift_ellipse` | DRIFT, `ELLIPSE(0.02, 0.01)` | 8694 ± 34 |
| `drift_ellipse_np2` | the same on 2 MPI ranks (`--np 2`) | 8694 ± 34 |
| `drift_rectangle` | DRIFT, `RECTANGLE(0.02, 0.01)` | 9160 ± 28 |
| `drift_noaperture` | DRIFT, no aperture (control) | 10000 |
| `collimator_ellipse` | COLLIMATOR, same aperture as `drift_ellipse` | same as `drift_ellipse` |
| `rbend_gap`, `sbend_gap` | bend with `ANGLE = 0`, `HGAP` bounds y, `HAPERT` bounds x | 9160 ± 28 |

`aperture_analysis.ipynb` browses the `.stat` files, draws phase-space heat maps and writes
ParaView files for the runs under `output/elements/aperture/`.

Run one case from the repo root; the output goes to the same path under `output/`:

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
export OPALX=/Users/rammann/Code/OPALX/opalx/build_serial/src/opalx
$PY -m opalxruns.run studies/elements/aperture/drift_ellipse/drift_ellipse.in
$PY -m opalxruns.process_run studies/elements/aperture/drift_ellipse      # plots and ParaView files for that run
```
