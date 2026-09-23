# muE4 in OPALX from the G4beamline field maps, placed by absolute lab pose

## Short answer

The whole muE4 channel can be rebuilt in OPALX from G4beamline's own field maps, placed where
G4beamline puts them, and **the two codes agree.** Tracking the same 11 muons through the same
19.25 m channel in both, the transverse positions agree to **24 µm horizontally and 2.6 µm
vertically** at the end of the line, and the momenta to **0.75 keV/c out of 28 MeV/c**. Through
the solenoid, at 1.7 m, the agreement is sub-micrometre.

The geometry is generated from the `.g4bl` deck rather than typed. A 42-element, all-posed
lattice tracks end to end in about 12 seconds and conserves momentum exactly.

## Why absolute placement

G4beamline bends the *centreline coordinate system* with `cornerarc`, which is not a magnet, and
then places a tabulated map at the joint between two of them, so the map's frame bisects the total
bend. OPALX's `SBEND` derives its field from `ANGLE`, and `RBEND` derives the element's own
orientation from half that angle — self-consistent for an ideal magnet, wrong for a measured one,
because with a map the field decides the real bend and `ANGLE` is only a claim about it.

The `FIELDMAP` element added alongside this study takes its field and its box from the map and is
placed by an absolute 6D lab pose, deriving nothing. That is also how G4beamline works, so the
comparison becomes one to one.

## Files

| path | what it is |
|---|---|
| `tools/mue4lib.py` | parses a G4beamline deck and replays its centreline walk into lab poses |
| `tools/make_lattice.py` | emits `poses.json` and `lattice.in` from `mue4_WsxOn.g4bl` |
| `tools/test_geometry.py` | gate 0: checks the walk against G4beamline's own `g4bl.out` |
| `lattice.in` | generated — 42 posed elements plus the `MUE4` line |
| `poses.json` | generated — every pose, map, scale and aperture, with its source line |
| `probe/probe.in` | the structural probe: 5 muons through the whole line, OPALX only |
| `cmp/make_case.py` | writes the shared particle set and a runnable G4beamline twin |
| `cmp/compare.py` | 11 hand-picked muons, per particle |
| `gauss/make_gauss.py` | samples one Gaussian beam, writes it in both codes' formats |
| `gauss/compare.py` | 2000-muon Gaussian: tracking and transmission, separately |

Regenerate with `python3 tools/make_lattice.py`; check with `python3 tools/test_geometry.py`.

## Two things that make this exact rather than approximate

**The map frame is exactly the bisector.** There is no z-based tie-break in G4beamline. `place`
resolves against `segmentCLVector[currentCL]`, and `currentCL` moves only when a `cornerarc` line
is *executed* (`BLCoordinates.cc:358`). The muE4 dipole maps sit between the two `cornerarc`
commands in the deck, so the frame is rotated by exactly 20°, not by 20° plus the second arc's
first kink.

**A `cornerarc` ends exactly on the true arc.** It is three kinks, not an arc, but the two
straight pieces have length `R·θ/2` each and run at `θ/2 ∓ β`, so their sum is
`2·R·(θ/2)·sin(θ/2)·cos β` — and `β` is chosen so `cos β = sin(θ/2)/(θ/2)`. That collapses to
`R(1−cos θ)` across and `R sin θ` along, the chord of the true arc, for any angle. Only the middle
vertex is off the arc, and nothing in muE4 is placed there.

## Conventions

- **Everything is posed.** OPALX rejects a beamline that mixes `ELEMEDGE` and 6D poses, and
  `FIELDMAP` rejects `ELEMEDGE`, so the whole line uses `X, Y, Z, THETA, PHI, PSI`.
- **Names carry a zero-padded ordinal.** An all-posed lattice is sorted by *name*, not position —
  `fieldStart()` returns 0.0 for every posed element — so `E01_…`, `E02_…` keeps every dump in
  beam order.
- **`SCALE` is the placement's `current=`, verbatim.** The G4beamline readers store absolute
  Tesla and do not normalise, so `SCALE = 1` reproduces a bare `fieldmap` placement.
- **Rotations become the pose.** `rotation=Y180` folds into `THETA`; `Y180,Z180` is a 180° flip
  about x; `Z180` is one about z. `make_lattice.py` composes `R_frame · R_element` as a matrix and
  decomposes it into OPALX's `Ry(THETA)·Rx(PHI)·Rz(PSI)`, then **reconstructs the matrix and
  checks it**, so the branch choice is proven rather than argued.
- **`ZSTOP` is accumulated path length, not lab z.** On the design orbit G4beamline's centreline z
  *is* arc length, so `ZSTOP` and the deck's centreline coordinate are directly comparable.
- **Slits are centred on the centreline.** G4beamline builds each one from two iron blocks offset
  to either side; the opening is `|offset| − half the block's own size` along that axis. The
  deck's own arithmetic does not give this (`0.5*(102+800)` uses the block's *width* where its
  *height* matters), so the generator computes it from the geometry — which recovers the intended
  190 / 105 / 190 / 192 / 109 mm half-gaps.

## What G4beamline objects become

| G4beamline | OPALX | note |
|---|---|---|
| `fieldmap` | `FIELDMAP` | `SCALE` = the placement's `current=` |
| `virtualdetector` | `MONITOR` | `DELETEONTRANSVERSEEXIT = FALSE`, no `L` |
| jaw pair (`box` ×2) | one `COLLIMATOR` | slit as a `rectangle` aperture, centred |
| beam pipe (`tubs`) | `COLLIMATOR` | circular aperture at the inner radius |
| `ASR61_BOX1` | dropped | no material in the deck; marks the field region only |

## Established

0. **The two codes agree, particle by particle.** `cmp/` runs the same 11 muons through both,
   comparing at each detector in centreline coordinates:

   | plane | n | max \|dx\| | max \|dy\| | max \|dp\| |
   |---|---|---|---|---|
   | LEMSdet (0.54 m) | 11 | 0.0 µm | 0.0 µm | 0.43 eV/c |
   | MiddleOfSolenoid (1.29 m) | 11 | 0.0 µm | 0.0 µm | 0.31 eV/c |
   | EndOfSolenoid (1.70 m) | 11 | 0.6 µm | 0.0 µm | 0.21 eV/c |
   | **virtDet (19.25 m)** | 10 | **24.1 µm** | **2.6 µm** | **0.75 keV/c** |

   **Two corrections are needed or this measurement is meaningless**, and both were got wrong
   before they were got right:

   - *Compare in centreline coordinates, not global.* G4beamline writes 6 significant figures.
     In global coordinates the last detector is at z = 16549 mm, so the printed resolution is
     0.1 mm — which projected into the tilted detector plane is 36 µm of pure rounding noise,
     larger than the difference being measured. A global-frame comparison cannot resolve this at
     all; it returns the rounding floor and looks like a result. In centreline coordinates the
     transverse numbers are millimetres and the floor is ~0.01 µm.
   - *Drift the G4beamline record onto the plane.* A virtualdetector is a 1 mm thick volume that
     records on entry, so its z is 0.5 mm upstream of the plane an OPALX monitor interpolates
     to. At a 40 mrad angle that is 20 µm — again comparable to the difference. Each record is
     drifted forward along its own momentum first.

0b. **A 2000-muon Gaussian beam agrees too.** `gauss/` samples one distribution with a fixed
   seed and writes it in both formats, so neither code samples anything of its own. The deck's
   own beam parameters: point source, 130 mrad divergence in each plane, 28 MeV/c.

   | plane | in both | p50 | p95 | p99 | max | d\|p\| |
   |---|---|---|---|---|---|---|
   | LEMSdet (0.54 m) | 1429 | 0.1 µm | 0.2 µm | 0.4 µm | 0.5 µm | 0.05 keV/c |
   | MiddleOfSolenoid (1.29 m) | 1598 | 0.1 µm | 0.5 µm | 0.6 µm | 0.7 µm | 0.05 keV/c |
   | EndOfSolenoid (1.70 m) | 1598 | 0.3 µm | 4.0 µm | 6.9 µm | 25.6 µm | 0.05 keV/c |
   | **virtDet (19.25 m)** | 865 | **3.4 µm** | **22.1 µm** | 93.1 µm | 5412 µm | 0.05 keV/c |

   **`|p|` agrees to 50 eV/c out of 28 MeV/c for every particle at every plane** — two parts per
   million. That is the sharpest statement that both codes see the same field, and it is clean
   because `|p|` is conserved exactly in a magnetostatic field and so does not depend on where
   along the trajectory the record was taken. (The momentum *vector* does: a G4beamline
   virtualdetector records 0.5 mm upstream of the plane, and inside a field region that half
   millimetre of kick shows up as ~15 keV/c. Comparing vectors there measures the offset, not
   the physics.)

   One particle of 865 (0.12 %) diverges by 5.4 mm. It is not a material interaction: every
   G4beamline survivor still has exactly 28 MeV/c, and that particle's `|p|` matches OPALX's to
   0.05 keV/c. It is a large-amplitude trajectory sampling the nonlinear fringe of the maps,
   where a micrometre amplifies over the remaining metres. The criterion is therefore
   statistical — p95, p99 and an outlier fraction — with outliers named by id rather than hidden
   inside a wider tolerance.

   **Transmission agrees exactly through the solenoid and not at the end.** Two effects, both
   geometry rather than physics:

   | plane | g4bl | OPALX | in both | g4bl only |
   |---|---|---|---|---|
   | LEMSdet | 1429 | 1429 | 1429 | 0 |
   | MiddleOfSolenoid | 1598 | 1771 | 1598 | 0 |
   | EndOfSolenoid | 1598 | 1605 | 1598 | 0 |
   | virtDet | 1259 | 870 | 865 | 394 |

   A G4beamline virtualdetector is a finite disc that records only what hits it; an OPALX
   MONITOR is an unbounded plane. Applying LEMSdet's 112 mm radius to the OPALX records takes it
   from 199 to 145 against G4beamline's 145 — an exact match, and the predicted 72 % of a
   70 mm-sigma Gaussian.

   The 394 at the end are the collimator model. A G4beamline jaw is a finite iron block —
   FS61H spans |x| ∈ [190, 800] mm — so a wide-angle particle at |x| > 800 mm misses it. An
   OPALX aperture is one inside/outside test and kills everything beyond 190 mm. Measured by
   bracketing: OPALX keeps 79 of 200 with the apertures as modelled and 152 with them opened,
   while G4beamline keeps 122 — G4beamline sits between the two, exactly as an over-aggressive
   aperture model predicts. **This is a real limit of what OPALX apertures can express, and it
   means channel transmission is not yet comparable at the few-percent level.**

1. **Geometry, against G4beamline's own output.** All 47 placements reproduce the positions
   `g4bl.out` prints, to within half of its last printed digit. This includes the transverse
   offsets: the two FS61H jaws are 636 mm apart in global z purely from their ±495 mm x offset,
   and both come out right. This is a static check of *placement*, using an archived log; the
   tracking agreement in (0) is what confirms it dynamically.
2. **All six maps load.** `wsx_total` is `cylinder`; the rest are `grid`. None uses an `extend*`
   line, and every grid map is complete.
3. **`asr62shim_280_sm_track.g4blmap` is identically zero** — all 52,725 rows. Four of the six
   shim placements therefore contribute nothing. Only the ASR61 shim map carries real field, and
   it is substantial: max |By| = 0.31 T. `test_geometry.py` asserts this rather than assuming it.
4. **The quadrupole map is a nine-column file** (`x y z Bx By Bz Ex Ey Ez`) where the other five
   are six-column. This was found here and fixed in the reader — it now accepts nine columns only
   when the E values are zero, and says so plainly otherwise, instead of discarding them silently.
5. **The full line tracks.** 42 posed elements, 25,261 steps, ~12 s, exit clean:

   | quantity | result | expected |
   |---|---|---|
   | path length reached | 19.3997 m | `ZSTOP` 19.40, so `ZSTOP` fired |
   | final direction about y | +40.179° | +40.000° nominal |
   | momentum, start → end | 0.265005 → 0.265005 βγ | conserved exactly |
   | final position | (8.150, 0.000, 16.673) m | (8.143, 0, 16.663) from the generated geometry |
   | out-of-plane excursion | −0.1 mm | 0 |
   | transmission | 5 / 5 at all four monitors | — |

   The tracked orbit landing on the generated geometry is the dynamic confirmation of the static
   check in (1): the maps really do steer the beam along the line the walk predicted.

## Two problems this turned up

**The OPALX monitor drops the lowest-momentum particle.** At virtDet, OPALX recorded 10 of 11
while G4beamline recorded all 11. The missing one is id 8 at 27.72 MeV/c — the slowest in the
set. All 11 are alive in OPALX's own phase-space dump, so the particle is not lost, only
unrecorded. This matters for the production run: a momentum-dependent drop biases every ensemble
statistic, so `compare.py` names the dropped ids and their momenta rather than letting the
matched set quietly shrink.

**A comparison can be precision-limited and still look like a result.** The first version of
this comparison reported 35 µm, which was entirely G4beamline's ASCII rounding; the second
reported 18.5 µm vertically, which was entirely the uncorrected 0.5 mm detector-thickness drift.
Both passed their tolerance. Neither measured anything. Any tolerance here has to be checked
against the measurement floor before it means anything.

**The committed G4beamline deck cannot run on this machine.** Its `fieldmap` lines point at
`/Users/Andreas/Dropbox/...`, a path from the machine it was written on. `make_case.py` rewrites
them to the local copies; that is also why `g4bl.out` is an archived log rather than something
reproducible in place.

## Status

Done: the generator, gate 0 (geometry), the structural probe, and the per-particle track
comparison in `cmp/`.

Not done: the intermediate staged gates (solenoid alone, one dipole, triplet) — the full-line
comparison passing makes them a way to localise a future regression rather than a prerequisite;
the field gate (`DUMPEMFIELDS` against `fieldntuple`); and the production run on the real phase
space under the shared momentum window. The +40.179° final direction against a nominal +40° is
the measured fields' own deviation, now corroborated by both codes agreeing on it.

## Running

```bash
export OPALX=/Users/rammann/Code/OPALX/opalx/build_serial/src/opalx   # NOT build/src/opalx
python3 tools/test_geometry.py          # gate 0, no tracking
python3 tools/make_lattice.py           # regenerate lattice.in and poses.json
cd probe && mpirun -n 1 $OPALX probe.in --info 3 > run.log 2>&1

# the comparison
cd ../cmp && python3 make_case.py
export PATH="/Users/rammann/Code/G4BL/G4beamline-3.08.app/Contents/MacOS:$PATH"
g4bl mue4_cmp.g4bl > g4bl.log 2>&1
mpirun -n 1 $OPALX cmp.in --info 3 > run.log 2>&1
/opt/homebrew/Caskroom/miniconda/base/bin/python compare.py
```

`mue4_cmp.g4bl` keeps every `place` and `cornerarc` line byte-identical to `mue4_WsxOn.g4bl` --
the generator checks this -- and changes only the map paths, the beam, the physics list (decay
and secondaries off, since OPALX has no material physics) and `maxStep` (1 mm, converged; the
deck's own `param maxstep=10` is a lower-case typo that never took effect, so the archived run
used 100 mm).

Only `build_serial/src/opalx` has the `FIELDMAP` element; the older binary rejects `SCALE` at
parse time. Both codes resolve relative paths from the working directory, so `cd` into the case
first. One rank — particle-to-row identity across ranks is not guaranteed, and comparisons match
on the `id` dataset, never on row order.
