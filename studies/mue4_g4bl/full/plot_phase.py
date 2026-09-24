#!/usr/bin/env python3
"""Phase space at each recording plane: density maps and projections, both codes.

Two figures per plane.

  phase_heat_<z>.png   x-x', y-y' and x-y as density, G4beamline on top, OPALX
                       below on the SAME colour scale, and their difference
                       underneath on a scale centred at zero.
  phase_proj_<z>.png   the one-dimensional projections, both codes drawn over
                       each other, with the difference underneath.

Both codes are binned on identical edges, otherwise subtracting one from the
other measures the binning rather than the physics. The edges come from
percentiles rather than the extremes: nothing scrapes in this run, so a handful
of particles out of 10000 reach a metre off axis and scaling to them would put
the whole beam in one bin.
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, TwoSlopeNorm

from opalxruns import mue4_beam
from opalxruns.h5 import read_monitor

HERE = Path(__file__).resolve().parent

PLOTS = HERE / "plots"
PLOTS.mkdir(exist_ok=True)
OPALX, G4BL, INK, MUTED = "#1f5fd1", "#c1432d", "#222222", "#777777"
NBIN = 48
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                     "axes.edgecolor": MUTED, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": False,
                     "legend.frameon": False, "figure.facecolor": "white"})

# x, x', y, y' in the six coordinates mue4_beam.canonical returns
PAIRS = [(0, 1, "x [mm]", "x' [mrad]"), (2, 3, "y [mm]", "y' [mrad]"),
         (0, 2, "x [mm]", "y [mm]")]
PROJ = [(0, "x [mm]"), (1, "x' [mrad]"), (2, "y [mm]"), (3, "y' [mrad]")]


def planes():
    out = {}
    for h5 in HERE.glob("E*_PL*.h5"):
        m = re.search(r"PL(\d+)", h5.name)
        z = int(m.group(1))
        t = HERE / f"Z{z}.txt"
        if t.exists():
            out[z] = (t, h5)
    return dict(sorted(out.items()))


def edges(a, b, n=NBIN):
    """Bin edges covering both codes, from percentiles so the halo cannot set them."""
    v = np.concatenate([a, b])
    # 1 to 99 rather than the extremes: nothing scrapes in this run, so a few
    # particles reach a metre off axis and would otherwise set the whole range
    lo, hi = np.percentile(v, [1.0, 99.0])
    if hi <= lo:
        lo, hi = v.min(), v.max() + 1e-12
    pad = 0.06 * (hi - lo)
    return np.linspace(lo - pad, hi + pad, n + 1)


def heat(z, Xg, Xo, n):
    fig, ax = plt.subplots(3, 3, figsize=(10.4, 9.0))
    for c, (a, b, xl, yl) in enumerate(PAIRS):
        sa = 1e3
        sb = 1e3
        xe = edges(Xg[a] * sa, Xo[a] * sa)
        ye = edges(Xg[b] * sb, Xo[b] * sb)
        Hg, _, _ = np.histogram2d(Xg[a] * sa, Xg[b] * sb, bins=[xe, ye])
        Ho, _, _ = np.histogram2d(Xo[a] * sa, Xo[b] * sb, bins=[xe, ye])
        ext = [xe[0], xe[-1], ye[0], ye[-1]]
        # a log scale, because one bin in the core holds a thousand particles and a
        # linear scale then paints everything else white
        vmax = max(Hg.max(), Ho.max())
        norm = LogNorm(vmin=1, vmax=max(vmax, 2))
        for r, (H, lab) in enumerate(((Hg, "G4beamline"), (Ho, "OPALX"))):
            im = ax[r, c].imshow(np.ma.masked_equal(H.T, 0), origin="lower",
                                 aspect="auto", extent=ext, cmap="Blues", norm=norm)
            ax[r, c].set(xlabel=xl, ylabel=(f"{lab}\n{yl}" if c == 0 else yl))
            if r == 0:
                ax[r, c].set_title(f"{xl.split(' ')[0]} against {yl.split(' ')[0]}")
            fig.colorbar(im, ax=ax[r, c], label="particles per bin" if c == 2 else None)
        D = Ho - Hg
        lim = max(abs(D).max(), 1.0)
        im = ax[2, c].imshow(np.ma.masked_equal(D.T, 0), origin="lower", aspect="auto",
                             extent=ext, cmap="RdBu_r",
                             norm=TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim))
        ax[2, c].set(xlabel=xl, ylabel=("OPALX - G4beamline\n" + yl) if c == 0 else yl,
                     title=f"difference, worst {int(abs(D).max())} particles per bin")
        fig.colorbar(im, ax=ax[2, c], label="difference in count" if c == 2 else None)
    fig.suptitle(f"muE4 phase space at z = {z} mm   —   {n:,} particles, "
                 f"the same ones in both codes, binned on identical edges", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(PLOTS / f"phase_heat_{z:05d}.png", dpi=110)
    plt.close(fig)


def proj(z, Xg, Xo, n):
    fig, ax = plt.subplots(2, 4, figsize=(13.0, 5.0),
                           gridspec_kw=dict(height_ratios=[2.4, 1]))
    for c, (j, lab) in enumerate(PROJ):
        e = edges(Xg[j] * 1e3, Xo[j] * 1e3, 60)
        hg, _ = np.histogram(Xg[j] * 1e3, bins=e)
        ho, _ = np.histogram(Xo[j] * 1e3, bins=e)
        mid = 0.5 * (e[1:] + e[:-1])
        ax[0, c].step(mid, hg, where="mid", color=G4BL, lw=1.6, label="G4beamline")
        ax[0, c].step(mid, ho, where="mid", color=OPALX, lw=1.0, ls="--", label="OPALX")
        ax[0, c].set(ylabel="particles" if c == 0 else None, title=lab)
        ax[0, c].grid(True, color="#e6e6e6", lw=0.5)
        if c == 0:
            ax[0, c].legend(fontsize=7)
        d = ho - hg
        ax[1, c].step(mid, d, where="mid", color=INK, lw=1.0)
        ax[1, c].axhline(0, color=MUTED, lw=0.6)
        ax[1, c].set(xlabel=lab, ylabel="OPALX - G4BL" if c == 0 else None)
        ax[1, c].grid(True, color="#e6e6e6", lw=0.5)
        worst = int(abs(d).max())
        ax[1, c].text(0.98, 0.9, f"worst {worst} of {n:,}", transform=ax[1, c].transAxes,
                      ha="right", va="top", fontsize=6.5, color=MUTED)
    fig.suptitle(f"muE4 at z = {z} mm: the phase space projected onto each coordinate",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(PLOTS / f"phase_proj_{z:05d}.png", dpi=110)
    plt.close(fig)


def main():
    P = planes()
    only = [int(a) for a in sys.argv[1:]] or list(P)
    for z in only:
        if z not in P:
            continue
        gt, oh = P[z]
        Xg, Xo, ids, _ = mue4_beam.matched(mue4_beam.read_bltrack(gt), read_monitor(oh))
        heat(z, Xg, Xo, len(ids))
        proj(z, Xg, Xo, len(ids))
        print(f"  z = {z:5d} mm  {len(ids):,} particles")
    print(f"\nwrote {2*len(only)} figures into {PLOTS}")


if __name__ == "__main__":
    main()
