# Ring studies — multiple passes through the same elements

## square_ring

Closed ring: 4 x 90° SBENDs (arc length 1 m, hard edge) alternating with 4 x 1 m
drifts, circumference 8 m. 0.1 GeV electron pencil bunch, no space charge.
`ZSTOP = 1000` is far beyond the tracked path, so `MAXSTEPS = 80000` alone sets the
number of turns (~3 at ~0.3 mm per step).

Run, from the repo root (output in `output/features/ring/square_ring/`):
```bash
PY=/opt/homebrew/Caskroom/miniconda/base/bin/python
$PY -m opalxruns.run studies/features/ring/square_ring/square_ring.in
$PY studies/features/ring/square_ring/plot_trajectory.py   # x-z path + closure offset per turn
```

ParaView (particles + element tubes, via `opalxruns`):
see [square_ring/VISUALIZE.md](square_ring/VISUALIZE.md).

### Verified results (2026-09-02, serial Release build)

- Run completes; all 200 particles survive 80000 steps, final s = 23.98 m (3.00 turns).
- IndexMap printout: 24 sections, the element sequence D1..B4 repeated three times at
  increasing path length; one log line
  `IndexMap: ring detected (element re-entered at s = 8.000158 m)` when the element
  position file is written.
- `data/square_ring_ElementPositions.sdds` ends at s = 8.00016 m — the first turn only
  (by design: the element-to-range map keeps the first crossing per element).
- `plot_trajectory.py`: the lab-frame path traces the square three times.
  Closure offset to the start: 3.4 mm after turn 1, 3.1 mm after turn 2.
  Known feature: hard-edge bends integrated with the Boris pusher do not close the
  orbit exactly; closed-orbit finding is future work.

### Why the drifts carry an explicit APERTURE

The default transverse bound of an element is 1e6 m. In a ring, the drift on the
opposite side then also claims the particle position (its local z window matches; the
transverse offset of ~2.3 m is far below 1e6 m). Harmless for tracking — drifts
contribute no field — but the element sets contain far-side drifts, so the ring
detection fires half a turn early and the element position file is cut there.
`APERTURE = "ELLIPSE(0.1, 0.1)"` on the drifts removes the overlap. Any ring input
should bound its drifts this way.

## isis_ring

The ISIS synchrotron: 10 superperiods, circumference 163.363 m, elements placed by
`ELEMEDGE`, a coasting 70 MeV proton pencil bunch, no space charge, about 3 turns
(`MAXSTEPS = 45000` at `DT = 1e-10`). Adapted from a legacy OPAL-cycl input, which is
kept as `isis_ring.in.orig`; the header of `isis_ring.in` lists every substitution
(main dipoles as SBEND, RF cavities as drifts, scaling models removed, drift apertures).

```bash
$PY -m opalxruns.run studies/features/ring/isis_ring/isis_ring.in   # output in output/features/ring/isis_ring/
```

## Follow-up work (not part of the multi-pass change)

- Closed-orbit finding / closure correction for rings.
- `data/*_ElementPositions.py` (3D element mesh from `MeshGenerator::write`) is not
  valid Python: each element's triangle list is followed by a spurious extra `], `,
  so the script fails to parse (`SyntaxError: unmatched ']'`). Pre-existing and not
  ring-specific — every study's mesh file has it (checked mue4, bendtest). The mesh
  placement itself is per element in lab coordinates and would handle a ring once the
  writer is fixed.
