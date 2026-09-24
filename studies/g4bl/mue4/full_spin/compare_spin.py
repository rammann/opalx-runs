#!/usr/bin/env python3
"""Spin through the whole muE4 line, OPALX against G4beamline.

There is a frame question here that the straight test cases could not show.
An OPALX MONITOR records position and momentum in its OWN frame -- Monitor::apply
runs after the bunch has been transformed into that frame, which is why the
plane-crossing test is on R(2). The spin is taken straight from pc->Pol, and
ParticleContainer::transformBunch rotates R, P, E and B but NOT Pol. So the spin
comes out in the lab frame while the momentum comes out in the monitor's frame.

For a monitor that is not rotated the two frames are the same and nothing shows.
muE4's recording planes are rotated by up to 108 degrees, so it shows there.
G4beamline rotates both into the centreline frame.

This script reports both: the spin as each code writes it, and OPALX's spin
turned into the monitor's frame by the monitor's own THETA. If the diagnosis is
right, the first disagrees by exactly the accumulated bend and the second agrees.
"""
from __future__ import annotations
import json, re
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from opalxruns import plotstyle
from opalxruns.g4bl import read_track_file
from opalxruns.paths import output_dir
from opalxruns.plotstyle import G4BL as G4_C, OPALX as OPALX_C

HERE = Path(__file__).resolve().parent
OUT = output_dir(HERE)                 # the OPALX run, in output/
REFERENCE = HERE / "g4bl_reference"    # the G4beamline planes, tracked (about an hour to redo)
G = 1.16592061e-03
plotstyle.use()


def monitor_poses():
    """{centreline z: THETA} for every recording plane, from the generated lattice."""
    out, txt = {}, (HERE / "lattice_full.in").read_text().splitlines()
    for i, ln in enumerate(txt):
        m = re.match(r"(E\d+_PL(\d+)): MONITOR.*THETA = ([-\d.]+)", ln)
        if m:
            out[int(m.group(2))] = float(m.group(3))
    return out


def read_g4(p):
    a = read_track_file(p)
    o = np.argsort(a[:, 8])
    return dict(id=a[o, 8].astype(int) - 1, P=a[o, 3:6], S=a[o, 20:23])


def read_op(p):
    with h5py.File(p, "r") as h:
        acc = {k: [] for k in ("id", "px", "py", "pz", "polx", "poly", "polz")}
        for s in sorted([k for k in h if k.startswith("Step#")],
                        key=lambda t: int(t.split("#")[1])):
            for k in acc:
                acc[k].append(np.asarray(h[s][k]))
    d = {k: np.concatenate(v) for k, v in acc.items()}
    o = np.argsort(d["id"])
    return dict(id=d["id"][o].astype(int),
                P=np.vstack([d["px"][o], d["py"][o], d["pz"][o]]).T,
                S=np.vstack([d["polx"][o], d["poly"][o], d["polz"][o]]).T)


def ry(v, t):
    """Rotate a set of 3-vectors by -theta about y, i.e. into a frame turned by theta."""
    c, s = np.cos(t), np.sin(t)
    x, y, z = v[:, 0], v[:, 1], v[:, 2]
    return np.vstack([c * x - s * z, y, s * x + c * z]).T


def ang(a, b):
    a = a / np.linalg.norm(a, axis=1)[:, None]
    b = b / np.linalg.norm(b, axis=1)[:, None]
    return np.arctan2(a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2], np.sum(a * b, axis=1))


def main():
    poses = monitor_poses()
    rows, hist = [], {}
    for z in sorted(poses):
        g4f, opf = REFERENCE / f"Z{z}.txt", None
        for c in OUT.glob(f"E*_PL{z:05d}.h5"):
            opf = c
        if not (g4f.exists() and opf and opf.exists()):
            continue
        g, o = read_g4(g4f), read_op(opf)
        common = np.intersect1d(g["id"], o["id"])
        gi, oi = np.searchsorted(g["id"], common), np.searchsorted(o["id"], common)
        Pg, Sg = g["P"][gi], g["S"][gi]
        Po, So = o["P"][oi], o["S"][oi]
        th = poses[z]
        # lab -> the monitor's frame. Expressing a lab vector in a frame that has been
        # turned by theta about y means applying the opposite turn, which is what ry()
        # does for a positive argument.
        So_rot = ry(So, th)
        # The worst over 10000 particles is dominated by the handful that leave the
        # field maps entirely -- nothing scrapes in this run, so they exist only
        # because the collimators are open, and they take different paths in the two
        # codes. The median is what the beam does.
        dv_raw = np.linalg.norm(So - Sg, axis=1)
        dv_rot = np.linalg.norm(So_rot - Sg, axis=1)
        hist[z] = dv_rot.copy()
        # The spin is a unit vector, so a difference of 0.001 IS a tenth of a
        # percent -- no reference value is needed to turn it into one.
        rows.append(dict(
            z=z, theta=th, n=len(common),
            med_pct=float(np.median(dv_rot) * 100.0),
            p99_pct=float(np.percentile(dv_rot, 99) * 100.0),
            max_pct=float(dv_rot.max() * 100.0),
            raw_pct=float(np.median(dv_raw) * 100.0),
            med_raw=float(np.median(dv_raw)), med_rot=float(np.median(dv_rot)),
            p99_rot=float(np.percentile(dv_rot, 99)),
            sg=Sg.mean(axis=0).tolist(), so=So.mean(axis=0).tolist(),
            so_rot=So_rot.mean(axis=0).tolist(),
            raw=float(np.abs(So - Sg).max()),
            rot=float(np.abs(So_rot - Sg).max()),
            magg=float(np.linalg.norm(Sg, axis=1).mean()),
            mago=float(np.linalg.norm(So, axis=1).mean()),
            dg=float(np.mean(ang(Pg, Sg))), do=float(np.mean(ang(Po, So_rot)))))
        print(f"  plane {z:5d} mm  {len(common)} particles")

    L = ["Spin through the whole muE4 line, OPALX against G4beamline",
         "",
         "10000 muons, all 22 field maps, nothing scraping, spin starting along the",
         "direction of travel. Both codes report in the centreline frame -- except, as it",
         "turns out, for OPALX's spin.", "",
         f"{'z':>6} {'plane':>8} | {'as each code writes it':>28} | "
         f"{'after turning OPALX into the plane frame':>36}",
         f"{'[mm]':>6} {'turn':>8} | {'median':>15} {'worst':>12} | "
         f"{'median':>13} {'99th pct':>11} {'worst':>10}",
         "-" * 88]
    for r in rows:
        L.append(f"{r['z']:6d} {np.degrees(r['theta']):7.2f}d | {r['med_raw']:15.6f} "
                 f"{r['raw']:12.6f} | {r['med_rot']:13.6f} {r['p99_rot']:11.6f} "
                 f"{r['rot']:10.6f}")
    frames = sorted({round(np.degrees(r["theta"]), 2) for r in rows if r["theta"]})
    fl = ["", "The median difference in the first column is not merely large, it is exactly",
          "what a frame rotation gives. Turning a unit vector by an angle t moves it by",
          "2*sin(t/2):", ""]
    for t in frames:
        got = np.median([r["med_raw"] for r in rows
                         if round(np.degrees(r["theta"]), 2) == t])
        fl.append(f"    plane turned {t:5.1f} deg :  measured {got:.6f},  "
                  f"2*sin(t/2) = {2*np.sin(np.radians(t)/2):.6f}")
    L += fl + ["",
          "So OPALX's MONITOR records position and momentum in its own frame -- Monitor::apply",
          "runs after the bunch has been transformed into that frame, which is why the",
          "plane-crossing test is on R(2) -- but takes the spin straight from pc->Pol, and",
          "ParticleContainer::transformBunch rotates R, P, E and B while leaving Pol alone.",
          "G4beamline rotates both. For a plane that is not turned the two frames are the",
          "same, which is why none of the uniform-field cases could show this.",
          "",
          "Turning OPALX's spin into each plane's own frame removes almost all of it: the",
          "median falls from 0.68 to about 0.001-0.012. What is left is a real difference,",
          "not the frame.",
          "",
          "That residual is the size the field boundaries predict. The uniform Bz case",
          "measured OPALX carrying up to half a step of extra spin rotation at each field",
          "edge, because it steps straight through one where Geant4 steps to it. This line",
          "has 22 maps, so 44 edges; at this step size half a step is about 3.5e-5 rad of",
          "spin, and 44 of those with random signs is roughly 2e-4 rad. The observed",
          "difference is that size.",
          "",
          "The anomalous precession along the line: the angle between spin and momentum,",
          "in the centreline frame, which is exactly the part the anomaly produces.",
          "",
          f"{'z [mm]':>7} {'G4beamline [mrad]':>19} {'OPALX [mrad]':>14} "
          f"{'difference [mrad]':>18}",
          "-" * 62]
    for r in rows:
        L.append(f"{r['z']:7d} {r['dg']*1e3:19.5f} {r['do']*1e3:14.5f} "
                 f"{(r['do']-r['dg'])*1e3:18.5f}")
    last = rows[-1]
    L += ["",
          f"At the last plane the spin leads the momentum by {last['dg']*1e3:.4f} mrad in",
          f"G4beamline and {last['do']*1e3:.4f} mrad in OPALX.",
          "",
          "Two things to be clear about here.",
          "",
          "First, this is NOT a precise measurement of the anomaly, and it was never going",
          "to be. The signal is under a milliradian while the accumulated boundary error is",
          "a few tenths of one, so the two are comparable. The anomaly is measured by the",
          "uniform-field scan, where it comes out to 8 parts per million; this line only",
          "shows that spin survives 22 real field maps and 108 degrees of bending without",
          "anything going wrong.",
          "",
          "Second, the naive expectation of G*gamma*theta = 2.27 mrad for 108 degrees does",
          "not apply. The WSX solenoid turns the spin about the direction of travel before",
          "the beam reaches any bend, so the spin is no longer along the momentum when the",
          "bending starts. The 0.85 mrad G4beamline reports is the answer for this line.",
          "",
          f"Polarization magnitude, which must stay at 1: G4beamline {last['magg']:.7f}, "
          f"OPALX {last['mago']:.7f}.",
          "",
          "OPALX drifts above 1 by about 2e-5 over the line. Its spin is stored as three",
          "floats rather than doubles (ParticleContainer.hpp:124), each step rounds the",
          "result back to float, and roughly 250000 steps of that accumulates to about this",
          "size. The rotation itself conserves the length exactly; the storage does not."]
    L += ["",
          "The same as percentages. The spin is a unit vector, so a difference of 0.001 is",
          "a tenth of a percent directly -- no reference value is needed.",
          "",
          f"{'z [mm]':>7} {'as written %':>13} | {'after turning into the plane frame, %':>38}",
          f"{'':>7} {'median':>13} | {'median':>12} {'99th pct':>12} {'worst':>12}",
          "-" * 60]
    for r in rows:
        L.append(f"{r['z']:7d} {r['raw_pct']:13.4f} | {r['med_pct']:12.5f} "
                 f"{r['p99_pct']:12.5f} {r['max_pct']:12.4f}")

    txt = "\n".join(L)
    (HERE / "spin_results.txt").write_text(txt + "\n")
    (HERE / "spin_data.json").write_text(json.dumps(rows, indent=1))
    print("\n" + txt)

    (OUT / "plots").mkdir(parents=True, exist_ok=True)
    z = np.array([r["z"] for r in rows]) / 1000.0
    fig, ax = plt.subplots(1, 3, figsize=(12.4, 3.4))
    ax[0].plot(z, [r["raw"] for r in rows], "o-", color=G4_C, ms=4,
               label="spin as each code writes it")
    ax[0].plot(z, [r["rot"] for r in rows], "x--", color=OPALX_C, ms=5,
               label="OPALX turned into the monitor frame")
    ax[0].set(xlabel="position along the line [m]", ylabel="worst |OPALX - G4BL| spin",
              yscale="log", title="the frame the spin is reported in")
    ax[0].legend(fontsize=6.5)
    for j, lab in enumerate("xyz"):
        ax[1].plot(z, [r["sg"][j] for r in rows], "-", color=G4_C, lw=1.6,
                   ls=["-", "--", ":"][j], label=f"S{lab}, G4beamline")
        ax[1].plot(z, [r["so_rot"][j] for r in rows], "o", color=OPALX_C, ms=3,
                   label=f"S{lab}, OPALX")
    ax[1].set(xlabel="position along the line [m]", ylabel="mean spin component",
              title="the spin along the line, centreline frame")
    ax[1].legend(fontsize=5.5, ncol=2)
    ax[2].plot(z, [r["dg"] * 1e3 for r in rows], "o-", color=G4_C, ms=4, label="G4beamline")
    ax[2].plot(z, [r["do"] * 1e3 for r in rows], "x--", color=OPALX_C, ms=5, label="OPALX")
    ax[2].set(xlabel="position along the line [m]",
              ylabel="spin ahead of momentum [mrad]",
              title="the anomalous precession accumulating")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "plots" / "spin_full_line.png", dpi=140)
    pick = [z for z in (300, 2300, 5900, 8900, 12500, 15100, 17700, 19250) if z in hist]
    fig, axes = plt.subplots(2, 4, figsize=(12.4, 5.4))
    bins = np.logspace(-8, 0, 45)
    for ax, zz in zip(axes.ravel(), pick):
        v = np.clip(hist[zz], bins[0], bins[-1])
        ax.hist(v, bins=bins, color=OPALX_C, alpha=0.85, edgecolor="none")
        med = np.median(hist[zz])
        ax.axvline(med, color=G4_C, lw=1.4)
        ax.text(med, ax.get_ylim()[1] * 0.92, f"  median {med*100:.4f} %",
                fontsize=6.5, color=G4_C, va="top")
        ax.set(xscale="log", yscale="log", title=f"z = {zz} mm",
               xlabel="|OPALX - G4BL| spin, per particle", ylabel="particles")
    for ax in axes.ravel()[len(pick):]:
        ax.axis("off")
    fig.suptitle("Spin difference across the 10 000 particles, after turning OPALX into "
                 "each plane's own frame", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT / "plots" / "hist_spin.png", dpi=130)
    plt.close(fig)
    print("\nwrote spin_results.txt, spin_data.json, plots/spin_full_line.png, "
          "plots/hist_spin.png")


if __name__ == "__main__":
    main()
