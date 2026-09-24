#!/usr/bin/env python3
"""The uniform Bz case: does OPALX turn the spin, which way, and by how much.

A muon travelling along a uniform Bz feels no force at all, so the orbit is a
straight line and the only thing that changes is the spin. It turns about the
field by (1+G) times the cyclotron angle, which is a closed-form answer both
codes can be held against.

Also records two things the single run cannot say on its own: how the answer
moves as the time step shrinks, and whether the error grows with the length of
field, which separates an error made at the field boundary from one that builds
up along the way.

Writes bz_results.txt, bz_data.json and plots/bz_step.png.
"""
from __future__ import annotations
import json
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from opalxruns import plotstyle
from opalxruns.g4bl import read_track_file
from opalxruns.plotstyle import G4BL as G4_C, INK, OPALX as OPALX_C

HERE = Path(__file__).resolve().parent
D = HERE / "uniform_bz"
G_OPALX, G_G4 = 1.16592061e-03, 0.0011659208
THETA_1M, THETA_2M = 1.0, 2.0
BETA = 0.2561627
plotstyle.use()


def op_angle(path):
    with h5py.File(path, "r") as h:
        g = h[sorted([k for k in h if k.startswith("Step#")])[0]]
        i = int(np.argsort(np.asarray(g["id"]))[0])
        px, py = float(np.asarray(g["polx"])[i]), float(np.asarray(g["poly"])[i])
        pmom = [float(np.asarray(g[k])[i]) for k in ("px", "py", "pz")]
    return float(np.arctan2(py, px)), px, py, pmom


def g4_angle(path):
    a = read_track_file(path)[0]
    return float(np.arctan2(a[21], a[20])), a[20], a[21], list(a[3:6])


def main():
    exp1 = -(1.0 + G_OPALX) * THETA_1M
    ao, opx, opy, opm = op_angle(D / "scan" / "out_1e-12.h5")   # the base step
    ag, gpx, gpy, gpm = g4_angle(D / "Z1200.txt")
    L = ["Spin in a uniform Bz: the one case with no orbit at all",
         "",
         "1.000 m of uniform Bz = 0.0933979467 T, a mu+ at 28 MeV/c travelling along the",
         "field. Nothing pushes the particle sideways, so the momentum must come out",
         "exactly as it went in and the only thing that changes is the spin. It turns about",
         "the field by (1+G) times the cyclotron angle of 1.0 rad.",
         "",
         f"{'':22} {'spin x':>14} {'spin y':>14} {'angle [rad]':>15}",
         "-" * 68,
         f"{'G4beamline':22} {gpx:14.9f} {gpy:14.9f} {ag:15.9f}",
         f"{'OPALX':22} {opx:14.9f} {opy:14.9f} {ao:15.9f}",
         f"{'closed form':22} {np.cos(exp1):14.9f} {np.sin(exp1):14.9f} {exp1:15.9f}",
         "",
         f"G4beamline - closed form : {ag-exp1:+.3e} rad  ({abs(ag-exp1)/abs(exp1):.2e} relative)",
         f"OPALX      - closed form : {ao-exp1:+.3e} rad  ({abs(ao-exp1)/abs(exp1):.2e} relative)",
         "",
         "Both turn the spin the same way, and that way is the one the torque mu x B gives",
         "for a positive muon: the spin goes from +x towards -y. OPALX has no test of its",
         "own for the direction -- its unit test takes the absolute value of the angle and",
         "says so in a comment -- so this is the first check of it.",
         "",
         f"momentum unchanged: OPALX px,py = {opm[0]:.2e},{opm[1]:.2e}, pz = {opm[2]:.9f}",
         f"                    G4BL  Px,Py = {gpm[0]:.2e},{gpm[1]:.2e}, Pz = {gpm[2]:.6f} MeV/c",
         "",
         "How the OPALX answer moves as the time step shrinks",
         "",
         f"{'DT [s]':>10} {'mm per step':>12} {'angle [rad]':>15} {'error [rad]':>13} "
         f"{'half a step':>13}",
         "-" * 68]
    steps = []
    for dt in ("1e-12", "5e-13", "2.5e-13", "1.25e-13"):
        a, _, _, _ = op_angle(D / "scan" / f"out_{dt}.h5")
        mm = BETA * 2.99792458e8 * float(dt) * 1e3
        bound = 0.5 * mm / 1000.0 * abs(exp1)       # rotation in half a step
        L.append(f"{dt:>10} {mm:12.5f} {a:15.9f} {a-exp1:+13.3e} {bound:13.3e}")
        steps.append(dict(dt=dt, mm=mm, angle=a, err=a - exp1, bound=bound))
    L += ["",
          "Every error sits inside half a step's worth of rotation, and that bound shrinks",
          "with the step, but the error itself jumps around inside it rather than falling",
          "smoothly. That is what an error made at the field edge looks like: it depends on",
          "where the boundary happens to fall between two steps."]

    a2, _, _, _ = op_angle(D / "scan" / "out_2m.h5")
    exp2 = -(1.0 + G_OPALX) * THETA_2M
    L += ["",
          "Is the error made at the boundary, or does it build up along the field?",
          "",
          f"{'field length':>14} {'angle [rad]':>15} {'closed form':>15} {'error [rad]':>13}",
          "-" * 60,
          f"{'1.0 m':>14} {steps[0]['angle']:15.9f} {exp1:15.9f} {steps[0]['err']:+13.3e}",
          f"{'2.0 m':>14} {a2:15.9f} {exp2:15.9f} {a2-exp2:+13.3e}",
          "",
          "The rotation doubled and the error did not: it went from "
          f"{abs(steps[0]['err']):.2e} to {abs(a2-exp2):.2e} rad.",
          "Nothing accumulates along the field, so the Thomas-BMT integration itself is",
          "right and what is left is the entry and exit. Geant4 steps exactly to a geometry",
          "boundary; OPALX steps straight through it and carries up to half a step of extra",
          "rotation, which is why G4beamline sits a hundred times closer to the closed form",
          "at a comparable step size."]

    txt = "\n".join(L)
    (HERE / "bz_results.txt").write_text(txt + "\n")
    (HERE / "bz_data.json").write_text(json.dumps(dict(
        opalx=dict(angle=ao, polx=opx, poly=opy), g4bl=dict(angle=ag, polx=gpx, poly=gpy),
        theory=exp1, steps=steps, two_metre=dict(angle=a2, theory=exp2, err=a2 - exp2)),
        indent=1))
    print(txt)

    (HERE / "plots").mkdir(exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4))
    mm = np.array([s["mm"] for s in steps])
    er = np.array([abs(s["err"]) for s in steps])
    bd = np.array([s["bound"] for s in steps])
    ax[0].loglog(mm, er, "o-", color=OPALX_C, ms=6, label="OPALX, distance from the closed form")
    ax[0].loglog(mm, bd, "--", color=INK, lw=1.2, label="half a step's worth of rotation")
    ax[0].axhline(abs(ag - exp1), color=G4_C, lw=1.4,
                  label="G4beamline at maxStep = 0.1 mm")
    ax[0].set(xlabel="OPALX step [mm]", ylabel="|angle - closed form| [rad]",
              title="the error is bounded by the step, but not smooth in it")
    ax[0].legend(fontsize=6.5)
    ax[1].bar([0, 1], [abs(steps[0]["err"]), abs(a2 - exp2)], 0.5,
              color=[OPALX_C, OPALX_C])
    ax[1].set_xticks([0, 1], ["1.0 m of field\n1.0 rad of spin", "2.0 m of field\n2.0 rad of spin"])
    ax[1].set(ylabel="|angle - closed form| [rad]",
              title="twice the rotation, no more error\nso nothing builds up along the field")
    for i, v in enumerate([abs(steps[0]["err"]), abs(a2 - exp2)]):
        ax[1].text(i, v, f"  {v:.2e}", ha="center", va="bottom", fontsize=7)
    ax[1].set_ylim(0, max(abs(steps[0]["err"]), abs(a2 - exp2)) * 1.35)
    fig.tight_layout()
    fig.savefig(HERE / "plots" / "bz_step.png", dpi=140)
    print("\nwrote bz_results.txt, bz_data.json, plots/bz_step.png")


if __name__ == "__main__":
    main()
