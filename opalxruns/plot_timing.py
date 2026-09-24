#!/usr/bin/env python
"""plot_timing.py -- where an OPALX run spent its wall time.

``timing.dat`` holds two tables: total wall time per timer, and, below it, the
same timers with max / min / avg across ranks. Both are read; on a multi-rank run
the max-min spread is drawn as an error bar, which is what shows load imbalance.

``mainTimer`` is the whole run and is usually several times the longest
sub-timer, so plotting it as a bar would squash everything else into the axis
origin. It goes in the title instead, and each sub-timer bar is annotated with
its seconds and its share of it. Bars are sorted longest first.

Note that the timers overlap and do not partition the run -- ``computeMoments``
and ``External field eval`` are both inside the tracking loop, for instance -- so
the shares can sum to more or less than 100%. The "unaccounted" figure in the
subtitle is ``mainTimer`` minus the sum of the rest and is only meaningful when
the timers happen to be disjoint.

Usage
-----
    python plot_timing.py <run_dir>
    python plot_timing.py <run_dir> --timing path/to/timing.dat --outdir /tmp
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from opalxruns.opalx_run import Run, parse_timing

MAIN = "mainTimer"
BAR = "#1f5fd1"
ZERO_BAR = "#c8ccd4"


def plot_timing(timing_path, outdir, title=None):
    """Write ``<outdir>/timing.png``. Returns the path."""
    rows = parse_timing(timing_path)
    if not rows:
        raise ValueError(f"no timer rows in {Path(timing_path).name}")

    main = next((r for r in rows if r["name"] == MAIN), None)
    total = main["wall_tot"] if main and main["wall_tot"] else None
    others = [r for r in rows if r["name"] != MAIN]
    if not others:
        raise ValueError(f"{Path(timing_path).name} only has {MAIN}")

    ranks = max(r["ranks"] or 1 for r in rows)

    # On a multi-rank run the "Wall tot" table usually lists mainTimer alone, so
    # most timers only have max/min/avg. Plot the average then, not the max --
    # the max would overstate every bar. Which one is in use goes in the footer.
    use_avg = any(r["wall_tot"] is None for r in others)

    def value(r):
        v = r["wall_tot"] if not use_avg else r["wall_avg"]
        if v is None:
            v = r["wall_avg"] if r["wall_avg"] is not None else r["wall_max"]
        return v or 0.0

    others.sort(key=value)                     # ascending: barh draws bottom-up
    names = [r["name"] for r in others]
    vals = np.array([value(r) for r in others])

    # min/max across ranks, as an asymmetric error bar around the plotted value
    err = None
    if ranks > 1 and all(r["wall_max"] is not None and r["wall_min"] is not None
                         for r in others):
        lo = vals - np.array([r["wall_min"] for r in others])
        hi = np.array([r["wall_max"] for r in others]) - vals
        err = np.vstack([np.clip(lo, 0, None), np.clip(hi, 0, None)])

    y = np.arange(len(others))
    fig, ax = plt.subplots(figsize=(9.5, 0.42 * len(others) + 2.4))
    ax.barh(y, vals, color=[ZERO_BAR if v == 0 else BAR for v in vals],
            xerr=err, error_kw=dict(ecolor="0.35", capsize=2, lw=0.9))

    span = vals.max() if vals.max() > 0 else 1.0
    for i, v in enumerate(vals):
        label = f"{v:.3g} s" + (f"  ({100 * v / total:.1f}%)" if total else "")
        ax.text(v + span * 0.015, i, label, va="center", fontsize=8)

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("wall time [s]")
    ax.set_xlim(0, span * 1.28)
    ax.grid(True, axis="x", alpha=0.25)

    head = title or f"OPALX timers — {Path(timing_path).parent.name}"
    if total:
        head += f"   ({MAIN} = {total:.4g} s)"
    ax.set_title(head, fontsize=12)

    bits = [f"{ranks} rank(s)"]
    bits.append("bars: wall avg across ranks" if use_avg else "bars: wall total")
    if total:
        bits.append(f"unaccounted = {total - vals.sum():.4g} s")
    if err is not None:
        bits.append("error bars: min–max across ranks")
    bits.append("timers overlap; shares need not sum to 100%")
    footer = textwrap.fill("  ·  ".join(bits), width=96,
                           break_long_words=False, break_on_hyphens=False)
    fig.text(0.5, 0.005, footer, ha="center", va="bottom",
             fontsize=7.5, family="monospace")

    fig.tight_layout(rect=[0, 0.02 + 0.012 * footer.count("\n"), 1, 1])
    out = Path(outdir) / "timing.png"
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
    ap.add_argument("--timing", default=None,
                    help="a specific timing.dat (default: the run's)")
    ap.add_argument("--outdir", default=None,
                    help="output directory (default: <run>/plots)")
    args = ap.parse_args(argv)

    run = Run(args.run_dir, base=args.base)
    timing = Path(args.timing) if args.timing else run.timing
    if timing is None:
        print(f"plot_timing: no timing.dat in {run.dir}", file=sys.stderr)
        return 1

    outdir = Path(args.outdir) if args.outdir else run.plots_dir
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        out = plot_timing(timing, outdir, title=f"OPALX timers — {run.base}")
    except ValueError as exc:
        print(f"plot_timing: {exc}", file=sys.stderr)
        return 1
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
