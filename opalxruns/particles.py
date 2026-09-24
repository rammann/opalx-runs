"""The 13-particle set the transfer-matrix studies track, and the FROMFILE
particle file OPALX reads it from."""

from __future__ import annotations

import math
from pathlib import Path


def map_particles(bg0: float, eps: dict) -> list[tuple[str, list[float]]]:
    """13 rows: a reference particle plus a symmetric +/- step in each of the six
    coordinates. Each row is (label, [x, px, y, py, z, pz]) with the momenta in
    beta*gamma. The steps in a pair cancel in the centred difference, so the
    transfer matrix is taken about the axis.

    ``bg0`` is the reference beta*gamma, ``eps`` the step per coordinate, keyed
    x, xp, y, yp, z, delta.
    """
    rows: list[tuple[str, list[float]]] = [("ref", [0, 0, 0, 0, 0, bg0])]

    def add(label, x=0.0, xp=0.0, y=0.0, yp=0.0, z=0.0, delta=0.0):
        pz = bg0 * (1.0 + delta)
        px, py = pz * math.tan(xp), pz * math.tan(yp)
        norm = math.sqrt(px * px + py * py + pz * pz)
        scale = bg0 * (1.0 + delta) / norm
        rows.append((label, [x, px * scale, y, py * scale, z, pz * scale]))

    add("x+", x=+eps["x"]);                 add("x-", x=-eps["x"])
    add("xp+", xp=+eps["xp"]);              add("xp-", xp=-eps["xp"])
    add("y+", y=+eps["y"]);                 add("y-", y=-eps["y"])
    add("yp+", yp=+eps["yp"]);              add("yp-", yp=-eps["yp"])
    add("z+", z=+eps["z"]);                 add("z-", z=-eps["z"])
    add("delta+", delta=+eps["delta"]);     add("delta-", delta=-eps["delta"])
    return rows


def write_parts(path: Path, rows) -> None:
    """The particle file for DISTRIBUTION TYPE=FROMFILE: the count, a header,
    then one x px y py z pz row per particle."""
    lines = [str(len(rows)), "x px y py z pz"]
    for _label, v in rows:
        lines.append(" ".join(f"{c:.12e}" for c in v))
    path.write_text("\n".join(lines) + "\n")
