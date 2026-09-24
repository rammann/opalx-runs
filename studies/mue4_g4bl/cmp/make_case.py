#!/usr/bin/env python
"""make_case.py -- build the paired OPALX / G4beamline comparison case.

Writes three things from one source:

  parts.txt      the particle set, OPALX FROMFILE form
  beam.txt       the same particles, G4beamline #BLTrackFile form
  mue4_cmp.g4bl  a runnable twin of mue4_WsxOn.g4bl

The twin keeps every `place` and `cornerarc` line byte-identical -- that is the whole point,
since the geometry is what is under test -- and changes only what has to change:

  * map paths. The original deck points at /Users/Andreas/Dropbox/..., a path from the machine
    it was written on, so the deck as committed cannot run here at all. Rewritten to the local
    copies under g4bl-files/muE4/maps/.
  * the beam. The 1e7-muon phase-space file is replaced by the shared beam.txt, and the
    `rotation=Y270` that went with it is dropped because beam.txt is already in beam
    coordinates.
  * physics. Decay, secondaries and stochastic processes off, so the only thing that can move
    a particle is the field. OPALX has no material physics, so leaving them on would compare
    two different problems.
  * maxStep. The deck says `param maxstep=10`, but G4beamline's parameter is `maxStep` with a
    capital S -- the lower-case spelling never took effect and the archived run used the 100 mm
    default. Set to 1 mm here, converged, so any disagreement is attributable to the fields
    rather than to G4beamline's own integration error.
  * `collective` dropped. With no `spacecharge` command it never enabled collective mode
    anyway, so this changes nothing physical and removes a per-step monitor.
  * virtualdetectors write ascii in CENTERLINE coordinates, where the transverse numbers are
    millimetres and G4beamline's 6-significant-figure output resolves ~0.01 um.

Usage:  python3 make_case.py
"""

from __future__ import annotations

import re
from pathlib import Path

from opalxruns.paths import G4BL_FILES

HERE = Path(__file__).resolve().parent
G4BL_DIR = G4BL_FILES / "muE4"
SRC = G4BL_DIR / "g4bl_input" / "mue4_WsxOn.g4bl"
MAPS = G4BL_DIR / "maps"

MUON_MASS = 0.1056583755  # GeV
P0 = 0.028  # GeV/c, the channel's design momentum

# (x [mm], y [mm], dp/p). Kept small: these probe the optics, not the collimator edges, so
# nothing is scraped in either code and the comparison is not confused by jaw geometry.
OFFSETS = [
    (0.0, 0.0, 0.0),
    (5.0, 0.0, 0.0),
    (-5.0, 0.0, 0.0),
    (0.0, 5.0, 0.0),
    (0.0, -5.0, 0.0),
    (5.0, 5.0, 0.0),
    (-5.0, -5.0, 0.0),
    (0.0, 0.0, 0.01),
    (0.0, 0.0, -0.01),
    (5.0, -5.0, 0.005),
    (-5.0, 5.0, -0.005),
]


def write_beam() -> None:
    with open(HERE / "parts.txt", "w") as f:
        f.write(f"{len(OFFSETS)}\n")
        f.write("x y z px py pz\n")
        for x, y, dpp in OFFSETS:
            pz = P0 * (1.0 + dpp) / MUON_MASS  # beta*gamma
            f.write(f"{x * 1e-3:.9e} {y * 1e-3:.9e} 0.0 0.0 0.0 {pz:.9e}\n")

    with open(HERE / "beam.txt", "w") as f:
        f.write("#BLTrackFile written by make_case.py\n")
        f.write("#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\n")
        f.write("#mm mm mm MeV/c MeV/c MeV/c ns - - - - -\n")
        for i, (x, y, dpp) in enumerate(OFFSETS):
            pz = P0 * (1.0 + dpp) * 1000.0  # MeV/c
            f.write(f"{x:.9g} {y:.9g} 0 0 0 {pz:.9g} 0 -13 {i + 1} 1 0 1\n")

    print(f"wrote parts.txt and beam.txt ({len(OFFSETS)} muons)")


def write_twin() -> None:
    out = []
    for raw in SRC.read_text().splitlines():
        line = raw

        # Map paths: the committed deck points at a Dropbox path on another machine.
        if line.strip().startswith("fieldmap "):
            line = re.sub(r"file=\S+", lambda m: f"file={MAPS / Path(m.group(0)[5:]).name}", line)

        elif line.strip().startswith("physics "):
            line = "physics QGSP_BERT disable=Decay spinTracking=0 doStochastics=0"

        elif line.strip() == "collective":
            line = "# collective  -- dropped: no spacecharge command, so it never did anything"

        elif line.strip().startswith("param maxstep"):
            line = "param maxStep=1.0 deltaChord=0.01   # note the capital S; the deck's "
            line += "`maxstep` never took effect"

        elif line.strip().startswith("beam "):
            line = "beam ascii filename=beam.txt"

        elif line.strip().startswith("trackcuts"):
            line = "trackcuts keep=mu+ killSecondaries=1"

        elif line.strip().startswith("virtualdetector "):
            # Replace any format= the deck already set (virtDet asks for rootExtended) rather
            # than appending a second one.
            line = re.sub(r"\s*format=\S+", "", line.rstrip())
            line += " format=ascii coordinates=centerline"

        elif line.strip().startswith(("profile ", "zntuple ")):
            line = "# " + line

        out.append(line)

    (HERE / "mue4_cmp.g4bl").write_text("\n".join(out) + "\n")
    print("wrote mue4_cmp.g4bl")


if __name__ == "__main__":
    write_beam()
    write_twin()
