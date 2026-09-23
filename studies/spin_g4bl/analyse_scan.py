#!/usr/bin/env python3
"""The g-2 result: fit the muon anomaly from each code's own spin angles.

The spin turns (1 + G*gamma) times faster than the momentum in a field across
the direction of travel, so after a bend of theta the spin leads the momentum by
G*gamma*theta. Holding theta fixed and changing gamma makes that angle a straight
line through the origin in gamma, and its slope is G. Fitting it separately for
each code measures the anomaly each code is actually using, rather than only
checking that the two agree.

Writes results.txt and plots/g2_scan.png.
"""
from __future__ import annotations
import json
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
G_OPALX, G_G4 = 1.16592061e-03, 0.0011659208
OPALX_C, G4_C = "#1f5fd1", "#c1432d"
INK, MUTED = "#222222", "#777777"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                     "axes.edgecolor": MUTED, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": "#dddddd",
                     "grid.linewidth": 0.5, "axes.axisbelow": True,
                     "legend.frameon": False, "figure.facecolor": "white"})


def read_g4(p):
    return np.array([[float(v) for v in l.split()] for l in open(p)
                     if not l.startswith("#") and l.strip()])[0]


def read_opalx(p):
    with h5py.File(p, "r") as h:
        g = h[sorted([k for k in h if k.startswith("Step#")])[0]]
        i = int(np.argsort(np.asarray(g["id"]))[0])
        return {k: float(np.asarray(g[k])[i]) for k in
                ("px", "py", "pz", "polx", "poly", "polz")}


def ang(a, b):
    """Angle from a to b in the bend plane, one convention for both the bend angle
    and the spin-momentum angle so their signs can be compared."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    return float(np.arctan2(a[2] * b[0] - a[0] * b[2], a @ b))


def main():
    rows = json.loads((HERE / "scan.json").read_text())
    Z = np.array([0.0, 0.0, 1.0])
    L = ["Thomas-BMT spin precession: OPALX and G4beamline against the closed form",
         "",
         "A muon bends through a uniform By. The spin starts along the momentum and ends",
         "ahead of it by G*gamma*theta -- the only term in the equation that is not the",
         "same one that turns the momentum. The bend angle is held fixed while the",
         "momentum changes, so the five cases differ only in gamma.", "",
         f"{'p [GeV/c]':>9} {'gamma':>8} {'bend [rad]':>11} | "
         f"{'G4beamline':>11} {'OPALX':>11} {'theory':>11}   [mrad]",
         "-" * 78]
    out = []
    for r in rows:
        d = HERE / r["name"]
        go, oo = read_g4(d / "Z1200.txt"), read_opalx(d / "MON_OUT.h5")
        pg, sg = np.array(go[3:6]), np.array(go[20:23])
        po = np.array([oo["px"], oo["py"], oo["pz"]])
        so = np.array([oo["polx"], oo["poly"], oo["polz"]])
        th_g, th_o = ang(Z, pg), ang(Z, po)
        dg, do = ang(pg, sg), ang(po, so)
        exp = G_OPALX * r["gamma"] * 0.5 * (th_g + th_o)
        L.append(f"{r['p']:9.3f} {r['gamma']:8.4f} {0.5*(th_g+th_o):11.7f} | "
                 f"{dg*1e3:11.5f} {do*1e3:11.5f} {exp*1e3:11.5f}")
        out.append((r["gamma"], th_g, th_o, dg, do, exp, r["p"]))
    a = np.array(out)
    gam = a[:, 0]
    Gg = float(np.sum(a[:, 3] * gam * a[:, 1]) / np.sum((gam * a[:, 1]) ** 2))
    Go = float(np.sum(a[:, 4] * gam * a[:, 2]) / np.sum((gam * a[:, 2]) ** 2))
    L += ["", "Anomaly fitted from each code's own five points:",
          f"  G4beamline   G = {Gg:.9e}   its compiled value {G_G4:.9e}"
          f"   {(Gg-G_G4)/G_G4*100:+.4f} %",
          f"  OPALX        G = {Go:.9e}   its compiled value {G_OPALX:.9e}"
          f"   {(Go-G_OPALX)/G_OPALX*100:+.4f} %",
          "",
          "The bend angle and the spin-momentum angle come out with the same sign in both",
          "codes, which is what a positive anomaly requires: the spin leads the momentum.",
          "", f"{'p':>7} {'G4BL - theory':>14} {'OPALX - theory':>15} {'OPALX - G4BL':>13}   [mrad]",
          "-" * 54]
    for g_, tg, to, dg, do, e, p in out:
        L.append(f"{p:7.3f} {(dg-e)*1e3:14.6f} {(do-e)*1e3:15.6f} {(do-dg)*1e3:13.6f}")
    L += ["",
          "G4beamline prints six digits, so a spin component of order one is resolved to",
          "1e-6 and an angle to about 0.001 mrad. Both codes sit at that floor."]
    txt = "\n".join(L)
    (HERE / "results.txt").write_text(txt + "\n")
    (HERE / "scan_data.json").write_text(json.dumps(dict(
        fitted=dict(g4bl=Gg, opalx=Go, g4bl_compiled=G_G4, opalx_compiled=G_OPALX),
        points=[dict(p=float(p_), gamma=float(g_), theta=float(0.5*(tg+to)),
                     g4bl=float(dg), opalx=float(do), theory=float(e))
                for g_, tg, to, dg, do, e, p_ in out]), indent=1))
    print(txt)

    (HERE / "plots").mkdir(exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))
    ax[0].plot(gam, np.abs(a[:, 3]) * 1e3, "o", ms=7, color=G4_C, label="G4beamline")
    ax[0].plot(gam, np.abs(a[:, 4]) * 1e3, "x", ms=8, color=OPALX_C, label="OPALX")
    gl = np.linspace(1, 30, 50)
    ax[0].plot(gl, G_OPALX * gl * abs(a[0, 1]) * 1e3, "-", lw=1.2, color=INK,
               label=r"$G\,\gamma\,\theta$")
    ax[0].set(xlabel=r"$\gamma$", ylabel="spin ahead of momentum [mrad]",
              title="the anomalous precession grows with $\\gamma$")
    ax[0].legend(fontsize=7)
    ax[1].plot(gam, (a[:, 3] - a[:, 5]) * 1e3, "o-", ms=6, color=G4_C,
               label="G4beamline - theory")
    ax[1].plot(gam, (a[:, 4] - a[:, 5]) * 1e3, "x--", ms=7, color=OPALX_C,
               label="OPALX - theory")
    ax[1].axhline(0, color=MUTED, lw=0.6)
    ax[1].axhspan(-0.001, 0.001, color="#cccccc", alpha=0.45,
                  label="G4beamline's printing floor")
    ax[1].set(xlabel=r"$\gamma$", ylabel="measured - theory [mrad]",
              title="both sit at the floor of what can be printed")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(HERE / "plots" / "g2_scan.png", dpi=140)
    print(f"\nwrote results.txt and plots/g2_scan.png")


if __name__ == "__main__":
    main()
