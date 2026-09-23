#!/usr/bin/env python3
"""A plan view of the muE4 line: where the magnets and the recording planes are.

Seen from above. The path is the reference orbit OPALX actually tracked, read
from the .stat file, so it is the real bent geometry rather than a sketch. Each
field map is drawn at its own position and turned by its own angle; each
recording plane is drawn as a short bar across the orbit.
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "g4bl_compare"))
import cmplib

INK, MUTED = "#222222", "#777777"
# One fixed colour per kind of magnet, from the palette checked earlier in this work.
KIND = {"solenoid": "#2a78d6", "dipole": "#eb6834", "quadrupole": "#1baf7a"}
PLANE = "#4a3aa7"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 10, "axes.labelsize": 8,
                     "axes.edgecolor": MUTED, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": "#e4e4e4",
                     "grid.linewidth": 0.5, "axes.axisbelow": True,
                     "legend.frameon": False, "figure.facecolor": "white"})


def kind_of(fmapfn):
    n = Path(fmapfn).name
    if "wsx" in n:
        return "solenoid"
    if "qsm" in n:
        return "quadrupole"
    return "dipole"


def elements():
    maps, planes = [], []
    txt = (HERE / "lattice_full.in").read_text().splitlines()
    for i, ln in enumerate(txt):
        m = re.match(r"(E\d+)_(\w+): (FIELDMAP|MONITOR), X = ([-\d.]+), Y = ([-\d.]+), "
                     r"Z = ([-\d.]+), THETA = ([-\d.]+)", ln)
        if not m:
            continue
        x, z, th = float(m.group(4)), float(m.group(6)), float(m.group(7))
        clz = re.search(r"centreline z = ([\d.]+)", ln)
        clz = float(clz.group(1)) if clz else None
        if m.group(3) == "FIELDMAP":
            fn = re.search(r'FMAPFN = "([^"]+)"', txt[i + 1])
            maps.append(dict(x=x, z=z, th=th, kind=kind_of(fn.group(1)) if fn else "dipole",
                             cl=clz))
        else:
            planes.append(dict(x=x, z=z, th=th, cl=clz))
    return maps, planes


def main():
    maps, planes = elements()
    df = cmplib.read_stat(HERE / "full.stat")
    rz, rx = np.asarray(df["ref_z"]), np.asarray(df["ref_x"])

    fig, ax = plt.subplots(figsize=(10.6, 6.4))
    ax.plot(rz, rx, "-", color=INK, lw=1.3, zorder=2, label="the beam's path")

    seen = set()
    for m in maps:
        c = KIND[m["kind"]]
        # a short bar across the orbit, turned the way the magnet is turned
        L, W = 0.30, 0.16
        dz, dx = np.cos(m["th"]), np.sin(m["th"])
        ax.plot([m["z"] - W * dx, m["z"] + W * dx], [m["x"] + W * dz, m["x"] - W * dz],
                "-", color=c, lw=5, solid_capstyle="butt", alpha=0.85, zorder=3,
                label=m["kind"] if m["kind"] not in seen else None)
        seen.add(m["kind"])

    for i, p in enumerate(planes):
        dz, dx = np.cos(p["th"]), np.sin(p["th"])
        W = 0.26
        ax.plot([p["z"] - W * dx, p["z"] + W * dx], [p["x"] + W * dz, p["x"] - W * dz],
                "-", color=PLANE, lw=1.4, zorder=4,
                label="recording plane" if i == 0 else None)
        # four alternating offsets, because the first six planes sit inside two
        # metres of a seventeen metre line and one offset puts them on top of
        # each other
        off = [(0, 11), (0, -16), (0, 21), (0, -26)][i % 4]
        ax.annotate(f"{p['cl']:.0f}", (p["z"], p["x"]), textcoords="offset points",
                    xytext=off, ha="center", fontsize=5.4, color=PLANE,
                    arrowprops=dict(arrowstyle="-", color=PLANE, lw=0.4,
                                    shrinkA=0, shrinkB=2, alpha=0.5))

    ax.annotate("beam enters here", (rz[0], rx[0]), textcoords="offset points",
                xytext=(34, -30), fontsize=8, color=INK, ha="left",
                arrowprops=dict(arrowstyle="->", color=INK, lw=1))
    ax.set(xlabel="lab z [m]", ylabel="lab x [m]",
           title="muE4 seen from above: 22 field maps and 22 recording planes\n"
                 "numbers are the distance along the line in mm; "
                 "108 degrees of bending in three magnets")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(float(rz.min()) - 1.0, float(rz.max()) + 1.2)
    ax.set_ylim(float(rx.min()) - 1.4, float(rx.max()) + 1.4)
    ax.legend(loc="upper left", fontsize=7.5, ncol=2)
    fig.tight_layout()
    (HERE / "plots").mkdir(exist_ok=True)
    fig.savefig(HERE / "plots" / "layout.png", dpi=140)
    print(f"wrote plots/layout.png  "
          f"({len(maps)} field maps, {len(planes)} recording planes)")
    kinds = {}
    for m in maps:
        kinds[m["kind"]] = kinds.get(m["kind"], 0) + 1
    print("  ", kinds)


if __name__ == "__main__":
    main()
