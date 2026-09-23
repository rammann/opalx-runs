#!/usr/bin/env python3
"""Compare the two codes over the whole muE4 line, plane by plane.

Writes:
  full_results.txt   one row per plane: both codes' numbers and their difference
  plots/plane_<z>.png   phase space at that plane, OPALX above, G4beamline below
  plots/along_line.png  beam size, beam centre and emittance against position
  plots/growth.png      how the per-particle difference grows along the line
"""
from __future__ import annotations
import json, re, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "g4bl_compare"))
import cmplib

PLOTS = HERE / "plots"
PLOTS.mkdir(exist_ok=True)
OPALX, G4BL = "#1f5fd1", "#c1432d"
INK, MUTED, GRIDC = "#222222", "#777777", "#dddddd"
plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRIDC, "grid.linewidth": 0.5, "axes.axisbelow": True,
    "legend.frameon": False, "figure.facecolor": "white"})


def planes():
    """{centreline z: (g4bl file, opalx file)} for every plane both codes wrote."""
    out = {}
    for h5 in HERE.glob("E*_PL*.h5"):
        m = re.search(r"PL(\d+)", h5.name)
        if not m:
            continue
        z = int(m.group(1))
        txt = HERE / f"Z{z}.txt"
        if txt.exists():
            out[z] = (txt, h5)
    return dict(sorted(out.items()))


def emit(X, a, b):
    """rms emittance in the (a, b) plane, in m rad."""
    u, v = X[a] - X[a].mean(), X[b] - X[b].mean()
    return float(np.sqrt(max((u @ u) * (v @ v) - (u @ v) ** 2, 0.0)) / len(u))


def phase_space(z, Xg, Xo, n):
    fig, ax = plt.subplots(2, 3, figsize=(10.2, 6.2))
    cols = [(0, 1, "x [mm]", "x' [mrad]"), (2, 3, "y [mm]", "y' [mrad]"),
            (0, 2, "x [mm]", "y [mm]")]
    for r, (X, lbl, col) in enumerate(((Xo, "OPALX", OPALX), (Xg, "G4beamline", G4BL))):
        for c, (a, b, xl, yl) in enumerate(cols):
            A = X[a] * 1e3
            B = X[b] * (1e3 if b in (0, 2) else 1e3)
            ax[r, c].scatter(A, B, s=2.5, alpha=0.16, color=col, edgecolors="none")
            ax[r, c].set(xlabel=xl, ylabel=yl)
            if c == 0:
                ax[r, c].set_ylabel(f"{lbl}\n{yl}")
            if r == 0:
                ax[r, c].set_title(["x against x'", "y against y'",
                                    "x against y"][c])
    # Same limits down each column so the two rows can be read against each other, and
    # taken from percentiles rather than the extremes. Nothing scrapes in this run, so a
    # handful of particles out of 10000 reach a metre off axis; scaling to them squashes
    # the beam itself into a dot. How many fall outside the view is printed on the panel.
    clipped = 0
    for c, (a, b, _, _) in enumerate(cols):
        xs = np.concatenate([Xg[a], Xo[a]]) * 1e3
        ys = np.concatenate([Xg[b], Xo[b]]) * 1e3
        xl = np.percentile(xs, [0.2, 99.8])
        yl = np.percentile(ys, [0.2, 99.8])
        mx, my = 0.08 * max(np.ptp(xl), 1e-9), 0.08 * max(np.ptp(yl), 1e-9)
        out = int(((xs < xl[0] - mx) | (xs > xl[1] + mx)
                   | (ys < yl[0] - my) | (ys > yl[1] + my)).sum())
        clipped = max(clipped, out)
        for r in range(2):
            ax[r, c].set_xlim(xl[0] - mx, xl[1] + mx)
            ax[r, c].set_ylim(yl[0] - my, yl[1] + my)
        if c == 2:
            for r in range(2):
                ax[r, c].set_aspect("equal", adjustable="box")
    if clipped:
        ax[0, 0].text(0.02, 0.97, f"{clipped} of {2*n} points outside the view",
                      transform=ax[0, 0].transAxes, fontsize=6.5, color=MUTED,
                      va="top")
    fig.suptitle(f"muE4 at centreline z = {z} mm   —   {n} particles, matched one to one",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(PLOTS / f"plane_{z:05d}.png", dpi=110)
    plt.close(fig)


def main():
    P = planes()
    if not P:
        print("no planes found -- has the run finished?")
        return 1
    rows, z_list = [], []
    for z, (gt, oh) in P.items():
        dg = cmplib.read_bltrack(gt)
        do = cmplib.read_monitor(oh)
        Xg, Xo, ids, miss = cmplib.matched(dg, do)
        if len(ids) < 10:
            continue
        vp = np.abs(Xo[[0, 2]] - Xg[[0, 2]]).max(axis=0)
        va = np.abs(Xo[[1, 3]] - Xg[[1, 3]]).max(axis=0)
        # The real line has collimators; this run does not, so a few particles fly far
        # outside any physical aperture and outside every field map. They dominate the
        # rms while being exactly the particles muE4 would have removed. Both sets of
        # numbers are reported: everything, and only what fits through the tightest
        # collimator muE4 has (190 mm across, 105 mm high).
        keep = ((np.abs(Xg[0]) < 0.190) & (np.abs(Xg[2]) < 0.105)
                & (np.abs(Xo[0]) < 0.190) & (np.abs(Xo[2]) < 0.105))
        rows.append(dict(
            z=z, ng=len(dg["id"]), no=len(do["id"]), n=len(ids), missing=len(miss),
            kept=int(keep.sum()),
            cut_rms_x=(float(Xg[0][keep].std()), float(Xo[0][keep].std())),
            cut_rms_y=(float(Xg[2][keep].std()), float(Xo[2][keep].std())),
            cut_mean_x=(float(Xg[0][keep].mean()), float(Xo[0][keep].mean())),
            cut_mean_y=(float(Xg[2][keep].mean()), float(Xo[2][keep].mean())),
            cut_max_pos=float(vp[keep].max()) if keep.any() else float("nan"),
            mean_x=(Xg[0].mean(), Xo[0].mean()), mean_y=(Xg[2].mean(), Xo[2].mean()),
            rms_x=(Xg[0].std(), Xo[0].std()), rms_y=(Xg[2].std(), Xo[2].std()),
            rms_xp=(Xg[1].std(), Xo[1].std()), rms_yp=(Xg[3].std(), Xo[3].std()),
            emit_x=(emit(Xg, 0, 1), emit(Xo, 0, 1)),
            emit_y=(emit(Xg, 2, 3), emit(Xo, 2, 3)),
            med_pos=float(np.median(vp)), p99_pos=float(np.percentile(vp, 99)),
            max_pos=float(vp.max()),
            med_ang=float(np.median(va)), p99_ang=float(np.percentile(va, 99)),
            max_ang=float(va.max())))
        z_list.append(z)
        phase_space(z, Xg, Xo, len(ids))
        print(f"  plane {z:5d} mm  {len(ids)} matched")

    # ---- numbers ------------------------------------------------------------
    L = []
    L.append("muE4 whole line, field maps, nothing scraping")
    L.append("OPALX against G4beamline, same particles, matched one to one")
    L.append("")
    L.append(f"{'z':>6} {'G4BL':>6} {'OPALX':>6} {'miss':>5} | "
             f"{'mean x [mm]':>22} | {'rms x [mm]':>22} | {'emit x [mm mrad]':>22}")
    L.append(f"{'[mm]':>6} {'n':>6} {'n':>6} {'':>5} | "
             f"{'G4BL':>10} {'OPALX':>10} | {'G4BL':>10} {'OPALX':>10} | "
             f"{'G4BL':>10} {'OPALX':>10}")
    L.append("-" * 118)
    for r in rows:
        L.append(f"{r['z']:6d} {r['ng']:6d} {r['no']:6d} {r['missing']:5d} | "
                 f"{r['mean_x'][0]*1e3:10.4f} {r['mean_x'][1]*1e3:10.4f} | "
                 f"{r['rms_x'][0]*1e3:10.4f} {r['rms_x'][1]*1e3:10.4f} | "
                 f"{r['emit_x'][0]*1e6:10.4f} {r['emit_x'][1]*1e6:10.4f}")
    L.append("")
    L.append(f"{'z':>6} | {'difference in beam centre and size [um]':>42} | "
             f"{'per particle, position [um]':>34}")
    L.append(f"{'[mm]':>6} | {'mean x':>10} {'mean y':>10} {'rms x':>10} {'rms y':>10} | "
             f"{'median':>10} {'99th pct':>10} {'worst':>10}")
    L.append("-" * 96)
    for r in rows:
        L.append(f"{r['z']:6d} | "
                 f"{(r['mean_x'][1]-r['mean_x'][0])*1e6:10.2f} "
                 f"{(r['mean_y'][1]-r['mean_y'][0])*1e6:10.2f} "
                 f"{(r['rms_x'][1]-r['rms_x'][0])*1e6:10.2f} "
                 f"{(r['rms_y'][1]-r['rms_y'][0])*1e6:10.2f} | "
                 f"{r['med_pos']*1e6:10.2f} {r['p99_pos']*1e6:10.2f} {r['max_pos']*1e6:10.2f}")
    L.append("")
    L.append("The same, counting only what fits through the tightest collimator muE4 has")
    L.append("(190 mm across, 105 mm high). Nothing scrapes in this run, so the rows above")
    L.append("include particles the real line would have removed -- and those are exactly")
    L.append("the particles that leave the field maps, where the two codes part company.")
    L.append("")
    L.append(f"{'z':>6} {'kept':>6} | {'mean x [mm]':>22} | {'rms x [mm]':>22} | "
             f"{'rms y [mm]':>22} | {'worst':>9}")
    L.append(f"{'[mm]':>6} {'of n':>6} | {'G4BL':>10} {'OPALX':>10} | {'G4BL':>10} "
             f"{'OPALX':>10} | {'G4BL':>10} {'OPALX':>10} | {'[um]':>9}")
    L.append("-" * 104)
    for r in rows:
        L.append(f"{r['z']:6d} {r['kept']:6d} | "
                 f"{r['cut_mean_x'][0]*1e3:10.4f} {r['cut_mean_x'][1]*1e3:10.4f} | "
                 f"{r['cut_rms_x'][0]*1e3:10.4f} {r['cut_rms_x'][1]*1e3:10.4f} | "
                 f"{r['cut_rms_y'][0]*1e3:10.4f} {r['cut_rms_y'][1]*1e3:10.4f} | "
                 f"{r['cut_max_pos']*1e6:9.1f}")
    L.append("")
    last = rows[-1]
    L.append(f"At the last plane, z = {last['z']} mm:")
    L.append(f"  beam centre x  {last['mean_x'][0]*1e3:.4f} mm (G4BL) vs "
             f"{last['mean_x'][1]*1e3:.4f} mm (OPALX), difference "
             f"{(last['mean_x'][1]-last['mean_x'][0])*1e6:.1f} um")
    L.append(f"  beam size  x   {last['rms_x'][0]*1e3:.4f} mm vs {last['rms_x'][1]*1e3:.4f} mm, "
             f"{abs(last['rms_x'][1]-last['rms_x'][0])/last['rms_x'][0]*100:.3f} % apart")
    L.append(f"  beam size  y   {last['rms_y'][0]*1e3:.4f} mm vs {last['rms_y'][1]*1e3:.4f} mm, "
             f"{abs(last['rms_y'][1]-last['rms_y'][0])/last['rms_y'][0]*100:.3f} % apart")
    L.append(f"  emittance  x   {last['emit_x'][0]*1e6:.4f} vs {last['emit_x'][1]*1e6:.4f} mm mrad")
    L.append(f"  per particle   median {last['med_pos']*1e6:.1f} um, "
             f"worst {last['max_pos']*1e6:.1f} um")
    L.append("")
    L.append(f"Counting only the {last['kept']} particles inside the real aperture:")
    L.append(f"  beam size  x   {last['cut_rms_x'][0]*1e3:.4f} mm vs "
             f"{last['cut_rms_x'][1]*1e3:.4f} mm, "
             f"{abs(last['cut_rms_x'][1]-last['cut_rms_x'][0])/last['cut_rms_x'][0]*100:.3f} % apart")
    L.append(f"  beam size  y   {last['cut_rms_y'][0]*1e3:.4f} mm vs "
             f"{last['cut_rms_y'][1]*1e3:.4f} mm, "
             f"{abs(last['cut_rms_y'][1]-last['cut_rms_y'][0])/last['cut_rms_y'][0]*100:.3f} % apart")
    L.append(f"  worst particle {last['cut_max_pos']*1e6:.1f} um")
    txt = "\n".join(L)
    (HERE / "full_results.txt").write_text(txt + "\n")
    (HERE / "full_data.json").write_text(json.dumps(rows, indent=1))
    print("\n" + txt)

    # ---- plots along the line ----------------------------------------------
    z = np.array([r["z"] for r in rows]) / 1000.0
    fig, ax = plt.subplots(1, 3, figsize=(12.2, 3.3))
    for j, (key, lab) in enumerate((("rms_x", "rms x"), ("rms_y", "rms y"))):
        ax[0].plot(z, [r[key][0] * 1e3 for r in rows], "-", color=G4BL, lw=1.8,
                   ls=["-", "--"][j], label=f"{lab}, G4beamline")
        ax[0].plot(z, [r[key][1] * 1e3 for r in rows], "o", color=OPALX, ms=3,
                   label=f"{lab}, OPALX")
    ax[0].set(xlabel="centreline position [m]", ylabel="beam size [mm]", title="beam size")
    ax[0].legend(fontsize=6, ncol=2)
    for j, (key, lab) in enumerate((("mean_x", "centre x"), ("mean_y", "centre y"))):
        ax[1].plot(z, [r[key][0] * 1e3 for r in rows], color=G4BL, lw=1.8,
                   ls=["-", "--"][j], label=f"{lab}, G4beamline")
        ax[1].plot(z, [r[key][1] * 1e3 for r in rows], "o", color=OPALX, ms=3,
                   label=f"{lab}, OPALX")
    ax[1].set(xlabel="centreline position [m]", ylabel="beam centre [mm]",
              title="beam centre, in centreline coordinates")
    ax[1].legend(fontsize=6, ncol=2)
    for j, (key, lab) in enumerate((("emit_x", "x"), ("emit_y", "y"))):
        ax[2].plot(z, [r[key][0] * 1e6 for r in rows], color=G4BL, lw=1.8,
                   ls=["-", "--"][j], label=f"emittance {lab}, G4beamline")
        ax[2].plot(z, [r[key][1] * 1e6 for r in rows], "o", color=OPALX, ms=3,
                   label=f"emittance {lab}, OPALX")
    ax[2].set(xlabel="centreline position [m]", ylabel="rms emittance [mm mrad]",
              title="emittance")
    ax[2].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(PLOTS / "along_line.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.4))
    ax[0].plot(z, [r["med_pos"] * 1e6 for r in rows], "o-", color=OPALX, ms=3,
               label="median")
    ax[0].plot(z, [r["p99_pos"] * 1e6 for r in rows], "s--", color=G4BL, ms=3,
               label="99th percentile")
    ax[0].plot(z, [r["max_pos"] * 1e6 for r in rows], "^:", color=INK, ms=3, label="worst")
    ax[0].set(xlabel="centreline position [m]", ylabel="|OPALX - G4BL| [um]",
              yscale="log", title="position difference, per particle")
    ax[0].legend(fontsize=6)
    ax[1].plot(z, [r["med_ang"] * 1e6 for r in rows], "o-", color=OPALX, ms=3,
               label="median")
    ax[1].plot(z, [r["p99_ang"] * 1e6 for r in rows], "s--", color=G4BL, ms=3,
               label="99th percentile")
    ax[1].plot(z, [r["max_ang"] * 1e6 for r in rows], "^:", color=INK, ms=3, label="worst")
    ax[1].set(xlabel="centreline position [m]", ylabel="|OPALX - G4BL| [urad]",
              yscale="log", title="angle difference, per particle")
    ax[1].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(PLOTS / "growth.png", dpi=130)
    plt.close(fig)
    print(f"\nplots in {PLOTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
