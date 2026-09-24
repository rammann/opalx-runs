#!/usr/bin/env python
"""make_lattice.py -- turn the G4beamline muE4 deck into OPALX elements.

Reads mue4_WsxOn.g4bl, replays its centreline walk (opalxruns/mue4.py), and writes:

  poses.json          every element's lab pose, map, scale and aperture
  lattice.in          the OPALX element definitions and LINE, for inclusion in a case deck

Everything is generated rather than typed, because the centreline walk is the error-prone
part and hand-copying 47 placements is how a wrong pose gets in unnoticed.

How G4beamline objects map onto OPALX ones:

  fieldmap          -> FIELDMAP, SCALE = the placement's own current=
  virtualdetector   -> MONITOR
  box pair (jaws)   -> one COLLIMATOR carrying the slit as a rectangle aperture
  tubs (beam pipe)  -> COLLIMATOR with a circular aperture at the inner radius
  box ASR61_BOX1    -> dropped: it has no material, it only marks the field region

The jaw pairs are the fiddly case. G4beamline builds a slit from two iron blocks offset to
either side, so the opening is |offset| - half the block's own size along that axis. The
numbers in the deck's arithmetic do not give this directly -- `0.5*(102+800)` uses the block's
width where its height is what matters -- so the gap is computed from the geometry instead,
which reproduces the intended 190/105/190/192/109 mm half-gaps.

Usage:  python3 tools/make_lattice.py [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from opalxruns import mue4 as m4
from opalxruns.paths import G4BL_FILES

G4BL_DIR = G4BL_FILES / "muE4"
DECK = G4BL_DIR / "g4bl_input" / "mue4_WsxOn.g4bl"
MAPS = G4BL_DIR / "maps"

MM = 1.0e-3

# Placed but not modelled, with the reason. Kept explicit so nothing is dropped by accident.
SKIP = {
    "ASR61_BOX11": "no material in the deck; marks the field region only",
    "ASR61_BOX12": "no material in the deck; marks the field region only",
}

# An aperture limit on the axis a slit does not constrain. Wide enough never to scrape,
# finite so the bounding box stays sane.
OPEN_M = 2.0


def jaw_pairs(placements):
    """Group the collimator jaws into slits.

    Two placements of the same solid, offset to either side of the same axis at the same
    centreline z, are one slit. Returns (base_name, axis, half_gap_mm, cl_z, solid).
    """
    slits = []
    seen = set()
    for i, a in enumerate(placements):
        if a.solid.kind != "box" or a.index in seen or a.name in SKIP:
            continue
        for b in placements[i + 1 :]:
            if b.solid is not a.solid or b.index in seen:
                continue
            if abs(b.cl_z - a.cl_z) > 1e-6:
                continue
            for axis, off_a, off_b, extent_key in (
                ("x", a.off_x, b.off_x, "width"),
                ("y", a.off_y, b.off_y, "height"),
            ):
                if abs(off_a) < 1e-9 or abs(off_a + off_b) > 1e-6:
                    continue
                half_block = 0.5 * float(a.solid.attrs[extent_key])
                half_gap = abs(off_a) - half_block
                slits.append((a.name.rsplit("_", 1)[0], axis, half_gap, a.cl_z, a.solid, a))
                seen.update({a.index, b.index})
                break
            if a.index in seen:
                break
    return slits


def build(placements, solids):
    """One record per OPALX element, in beam order."""
    slits = jaw_pairs(placements)
    slit_by_index = {rep.index: s for *s, rep in ((*s[:5], s[5]) for s in slits)}

    records = []
    for p in placements:
        if p.name in SKIP:
            continue

        common = {
            "g4bl_name": p.name,
            "g4bl_line": p.line,
            "cl_z_mm": p.cl_z,
            "X": p.lab[0] * MM,
            "Y": p.lab[1] * MM,
            "Z": p.lab[2] * MM,
            "THETA": p.pose[0],
            "PHI": p.pose[1],
            "PSI": p.pose[2],
            "frame_deg": p.frame_deg,
            "rotation": p.rotation,
        }

        if p.solid.kind == "fieldmap":
            fname = p.solid.attrs["file"]
            records.append(
                {
                    **common,
                    "type": "FIELDMAP",
                    "fmapfn": str(MAPS / Path(fname).name),
                    "scale": p.current if p.current is not None else 1.0,
                }
            )
        elif p.solid.kind == "virtualdetector":
            records.append({**common, "type": "MONITOR"})
        elif p.solid.kind == "tubs":
            inner = float(p.solid.attrs["innerRadius"]) * MM
            records.append(
                {
                    **common,
                    "type": "COLLIMATOR",
                    "L": float(p.solid.attrs["length"]) * MM,
                    "aperture": f"circle({2 * inner:.6f})",
                }
            )
        elif p.index in slit_by_index:
            base, axis, half_gap_mm, _, solid = slit_by_index[p.index]
            gap = 2.0 * half_gap_mm * MM
            # A slit is centred on the centreline, not on either jaw. `common` carries the
            # position of the jaw this record was built from, which is offset by half the
            # opening plus half the block; undo that so the aperture sits where the opening is.
            across = (math.cos(math.radians(p.frame_deg)), -math.sin(math.radians(p.frame_deg)))
            common = {
                **common,
                "X": (p.lab[0] - p.off_x * across[0]) * MM,
                "Y": 0.0,
                "Z": (p.lab[2] - p.off_x * across[1]) * MM,
            }
            aperture = (
                f"rectangle({gap:.6f},{2 * OPEN_M:.6f})"
                if axis == "x"
                else f"rectangle({2 * OPEN_M:.6f},{gap:.6f})"
            )
            records.append(
                {
                    **common,
                    "g4bl_name": base,
                    "type": "COLLIMATOR",
                    "L": float(solid.attrs["length"]) * MM,
                    "aperture": aperture,
                    "half_gap_mm": half_gap_mm,
                    "axis": axis,
                }
            )
        # the partner jaw of a slit, and anything unrecognised, falls through

    records.sort(key=lambda r: r["cl_z_mm"])
    for i, r in enumerate(records, start=1):
        r["name"] = f"E{i:02d}_{r['g4bl_name']}"
    return records


def emit_deck(records) -> str:
    """The OPALX element definitions and LINE."""
    out = [
        "// lattice.in -- GENERATED by tools/make_lattice.py, do not edit.",
        "//",
        "// Every element is placed by an absolute lab pose (X, Y, Z, THETA, PHI, PSI). OPALX",
        "// requires one placement convention per beamline and FIELDMAP rejects ELEMEDGE, so the",
        "// whole line is posed. Names carry a zero-padded ordinal because an all-posed lattice",
        "// is sorted by name, not by position (OpalBeamline.cpp fieldStart returns 0 for every",
        "// posed element), and this keeps every dump in beam order.",
        "",
    ]
    for r in records:
        pose = (
            f"X = {r['X']:.9f}, Y = {r['Y']:.9f}, Z = {r['Z']:.9f}, "
            f"THETA = {r['THETA']:.9f}, PHI = {r['PHI']:.9f}, PSI = {r['PSI']:.9f}"
        )
        note = f"  // g4bl {r['g4bl_name']} line {r['g4bl_line']}, centreline z = {r['cl_z_mm']:.4f} mm"
        if r["type"] == "FIELDMAP":
            out.append(f"{r['name']}: FIELDMAP, {pose},{note}")
            out.append(f'    FMAPFN = "{r["fmapfn"]}", SCALE = {r["scale"]:.6e};')
        elif r["type"] == "MONITOR":
            out.append(f"{r['name']}: MONITOR, {pose},{note}")
            out.append(f'    DELETEONTRANSVERSEEXIT = FALSE, OUTFN = "{r["name"]}";')
        elif r["type"] == "COLLIMATOR":
            out.append(f"{r['name']}: COLLIMATOR, {pose},{note}")
            out.append(f'    L = {r["L"]:.6f}, APERTURE = "{r["aperture"]}";')
        out.append("")

    out.append(f"MUE4: LINE = ({', '.join(r['name'] for r in records)});")
    out.append("")
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent.parent))
    args = ap.parse_args()
    out_dir = Path(args.out)

    placements, solids = m4.load(DECK)
    records = build(placements, solids)

    (out_dir / "poses.json").write_text(json.dumps(records, indent=2) + "\n")
    (out_dir / "lattice.in").write_text(emit_deck(records))

    kinds: dict[str, int] = {}
    for r in records:
        kinds[r["type"]] = kinds.get(r["type"], 0) + 1
    print(f"wrote {out_dir/'poses.json'} and {out_dir/'lattice.in'}")
    print(f"{len(records)} elements: " + ", ".join(f"{v} {k}" for k, v in sorted(kinds.items())))
    print()
    print(f"{'name':26s} {'type':11s} {'X [m]':>10s} {'Z [m]':>10s} {'THETA':>9s}")
    for r in records:
        print(
            f"{r['name']:26s} {r['type']:11s} {r['X']:10.4f} {r['Z']:10.4f} "
            f"{math.degrees(r['THETA']):9.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
