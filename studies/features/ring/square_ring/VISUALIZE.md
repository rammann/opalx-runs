# square_ring ParaView visualization

Uses `opalxruns` (see the repo README); run it with the Python that has
`h5py`/`numpy`/`vtk` (miniconda base), from the repo root, after the run.

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
RUN=output/features/ring/square_ring

# particles (electron mass, both frames)
$PY -m opalxruns.particles_to_vtk $RUN --frame both --mass 0.00051099895

# element bodies + reference orbit
$PY -m opalxruns.elements_to_vtk $RUN
```

Outputs land in `output/features/ring/square_ring/paraview/`:

| file | what |
|---|---|
| `bunch_lab.pvd` | time series of the bunch in lab coordinates (open this, not the single `.vtp`s) |
| `bunch_comoving.pvd` | the same dumps in OPALX's co-moving frame |
| `square_ring_elements.vtp` | element bodies: 4 drift pipes with the real 0.1 m bore, 4 solid bend arcs; colored by `element_type` or `element_id` |
| `square_ring_reforbit.vtp` | reference orbit polyline (3 turns) |

## In ParaView

1. **File → Open** `bunch_lab.pvd`, `square_ring_elements.vtp`,
   `square_ring_reforbit.vtp`; **Apply** each.
2. Elements: color by `element_type` (`0=DRIFT 2=DIPOLE`), **Opacity ≈ 0.3**.
3. Particles: color by `Ekin_MeV` or `Pmag`; press play — the bunch goes around
   the ring three times.

## Notes

- Verified for this run: lab-frame bunch centroids match the design path to
  < 1 mm at every dump across all 3 turns (the frame reconstruction from
  `DesignPath.dat` handles the full 2π rotation per turn).
- Elements with an aperture are hollow pipes: the drift bore is the deck's
  `APERTURE = "ELLIPSE(0.1, 0.1)"` (half-width 0.05 m), the wall is `--wall`
  (0.01 m) thick, and the ends are open annuli — particles fly through the
  opening. The hard-edge bends (HGAP = 0) have no aperture in OPALX, so they
  render as solid arcs of `--default-aperture` radius.
- Dumps are every 100 steps (`PSDUMPFREQ`), 800 frames for 3 turns.
