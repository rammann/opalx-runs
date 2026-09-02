#!/usr/bin/env python
"""plot_stat.py -- overview plots of an OPALX ``.stat`` file.

The ``.stat`` file is SDDS-like: a header of ``&parameter`` / ``&column``
declarations, then one row per tracking dump. Columns are documented in
``Structure/StatWriter.cpp``; the ones plotted here are

    numParticles          transmission
    energy, dE            mean bunch energy and its spread
    rms_x/y/s             RMS beam size
    emit_x/y/s            normalised emittance
    mean_x/y              bunch centroid, in the co-moving frame
    ref_x/ref_z           the reference particle in the lab frame
    Dx, Dy                dispersion

Panels whose columns are missing, empty or all-NaN are dropped rather than drawn
blank, so a run that does not compute dispersion just gets a smaller figure.

A run with several particle containers writes ``<base>_cN.stat``; each gets its
own figure.

Usage
-----
    python plot_stat.py <run_dir>
    python plot_stat.py <run_dir> --columns rms_x,rms_y,Dx --vs s
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from opalx_run import Run, parse_opal_stat  # noqa: E402

# (title, y-label, [(column, label, scale)], kind)
#   kind "xy"    -> plot the two columns against each other, equal aspect
#   kind "vs_x"  -> plot each column against the x axis
PANELS = [
    ("transmission", "N macro particles",
     [("numParticles", "N", 1.0)], "vs_x"),
    ("mean energy", "E [MeV]",
     [("energy", "energy", 1.0), ("dE", "dE", 1.0)], "vs_x"),
    ("RMS beam size", "rms [mm]",
     [("rms_x", "x", 1e3), ("rms_y", "y", 1e3), ("rms_s", "s", 1e3)], "vs_x"),
    ("normalised emittance", "emit [mm mrad]",
     [("emit_x", "x", 1e6), ("emit_y", "y", 1e6), ("emit_s", "s", 1e6)], "vs_x"),
    ("bunch centroid (co-moving)", "mean [mm]",
     [("mean_x", "x", 1e3), ("mean_y", "y", 1e3)], "vs_x"),
    ("reference orbit (lab)", "ref_x [m]",
     [("ref_z", "ref_x", 1.0)], "xy"),
    ("dispersion", "D [m]",
     [("Dx", "Dx", 1.0), ("Dy", "Dy", 1.0)], "vs_x"),
]

COLORS = ["#1f5fd1", "#c1432d", "#1e8a4c", "#8a5cc1", "#c98a1e"]


def _usable(df, name):
    """True if the column exists and holds at least one finite value."""
    if name not in df.columns:
        return False
    v = df[name].to_numpy(dtype=float)
    return v.size > 0 and np.isfinite(v).any()


def _draw_panel(ax, df, panel, xcol):
    title, ylabel, series, kind = panel
    if kind == "xy":
        zc, _, _ = series[0]
        ax.plot(df[zc].to_numpy(), df["ref_x"].to_numpy(), "-", color=COLORS[0], lw=1.4)
        ax.set_xlabel("ref_z [m]")
        ax.set_ylabel(ylabel)
        ax.set_aspect("equal", adjustable="datalim")
    else:
        x = df[xcol].to_numpy()
        for i, (col, label, scale) in enumerate(series):
            if not _usable(df, col):
                continue
            ax.plot(x, df[col].to_numpy() * scale, "-", lw=1.3,
                    color=COLORS[i % len(COLORS)], label=label)
        ax.set_xlabel(f"{xcol} [{'m' if xcol == 's' else 'ns'}]")
        ax.set_ylabel(ylabel)
        if len(series) > 1:
            ax.legend(fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.grid(True, alpha=0.25)


def plot_overview(stat_path, outdir, xcol="s"):
    """Write ``<outdir>/<statbase>_overview.png``. Returns the path."""
    meta, df = parse_opal_stat(stat_path)
    if df.empty:
        raise ValueError(f"{Path(stat_path).name} has no data rows")
    if xcol not in df.columns:
        xcol = "s" if "s" in df.columns else df.columns[0]

    panels = []
    for panel in PANELS:
        _, _, series, kind = panel
        if kind == "xy":
            if _usable(df, "ref_x") and _usable(df, "ref_z"):
                panels.append(panel)
        elif any(_usable(df, c) for c, _, _ in series):
            panels.append(panel)
    if not panels:
        raise ValueError(f"{Path(stat_path).name} has no plottable columns")

    ncol = 2
    nrow = (len(panels) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.1 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, panel in zip(axes, panels):
        _draw_panel(ax, df, panel, xcol)
    for ax in axes[len(panels):]:
        ax.axis("off")

    base = Path(stat_path).stem
    subtitle = "  ·  ".join(
        f"{k} = {v}" for k, v in meta.items()
        if k in ("species", "processors", "revision")
    )
    fig.suptitle(f"{base}   —   {len(df)} dumps, s = {df[xcol].iloc[0]:.3f} … "
                 f"{df[xcol].iloc[-1]:.3f}", fontsize=12)
    if subtitle:
        fig.text(0.5, 0.005, subtitle, ha="center", va="bottom", fontsize=7,
                 family="monospace")
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])

    out = Path(outdir) / f"{base}_overview.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def plot_columns(stat_path, outdir, columns, xcol="s"):
    """One figure with the named columns against ``xcol``. Returns the path."""
    _, df = parse_opal_stat(stat_path)
    have = [c for c in columns if _usable(df, c)]
    missing = [c for c in columns if c not in have]
    if missing:
        print(f"  note: no usable data for {', '.join(missing)}")
    if not have:
        raise ValueError("none of the requested columns hold data")
    if xcol not in df.columns:
        raise ValueError(f"no '{xcol}' column in {Path(stat_path).name}")

    x = df[xcol].to_numpy()
    fig, axes = plt.subplots(len(have), 1, figsize=(9, 2.4 * len(have)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, col in zip(axes, have):
        ax.plot(x, df[col].to_numpy(), "-", lw=1.3, color=COLORS[0])
        ax.set_ylabel(col)
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel(f"{xcol} [{'m' if xcol == 's' else 'ns'}]")
    base = Path(stat_path).stem
    fig.suptitle(f"{base} — {', '.join(have)} vs {xcol}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out = Path(outdir) / f"{base}_columns.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?", default=".",
                    help="OPALX run directory (default: cwd)")
    ap.add_argument("--base", default=None,
                    help="run basename, if the directory holds more than one deck")
    ap.add_argument("--stat", default=None,
                    help="a specific .stat file (default: every one of this run)")
    ap.add_argument("--columns", default=None,
                    help="comma-separated extra columns for a second figure")
    ap.add_argument("--vs", dest="xcol", default="s", choices=("s", "t"),
                    help="x axis: path length or time (default: s)")
    ap.add_argument("--outdir", default=None,
                    help="output directory (default: <run>/plots)")
    args = ap.parse_args(argv)

    run = Run(args.run_dir, base=args.base)
    stats = [Path(args.stat)] if args.stat else run.stat_files
    if not stats:
        print(f"plot_stat: no .stat file in {run.dir}", file=sys.stderr)
        return 1

    outdir = Path(args.outdir) if args.outdir else run.plots_dir
    outdir.mkdir(parents=True, exist_ok=True)

    rc = 0
    for stat in stats:
        try:
            out = plot_overview(stat, outdir, args.xcol)
            print(f"  wrote {out}")
        except (ValueError, KeyError) as exc:
            print(f"  {stat.name}: {exc}", file=sys.stderr)
            rc = 1
            continue
        if args.columns:
            cols = [c.strip() for c in args.columns.split(",") if c.strip()]
            try:
                out = plot_columns(stat, outdir, cols, args.xcol)
                print(f"  wrote {out}")
            except (ValueError, KeyError) as exc:
                print(f"  {stat.name}: {exc}", file=sys.stderr)
                rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
