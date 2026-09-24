#!/usr/bin/env python
"""plot_monitors.py -- phase-space plots at each MONITOR plane of an OPALX run.

Every MONITOR in the deck writes one HDF5 file ``<monitor_name>.h5`` (H5Part
``Step#N`` layout) holding the 6D phase space of the particles that crossed that
plane, in the co-moving / centerline frame -- x and y transverse to the
reference orbit, the same frame as g4bl ``coordinates=centerline``.

For each monitor this writes an x-x', y-y', x-y panel (x' = px/pz,
y' = py/pz), plus an envelope plot and ``stats.csv``
(s, N, rms_x, rms_y, emit_x, emit_y) over all planes.

Every ``Step#N`` group is read, not just the first
--------------------------------------------------
``Monitor::goOffline`` calls ``save(numPassages_m)`` and ``LossDataSink::splitSets``
**partitions** the recorded batch by crossing time, so one pass through a plane
can land in several ``Step#N`` groups. They are consecutive slices of the same
crossing, not repeats: in ``mue4_analytical`` the two sets of
``MON_01_SOL_IN.h5`` hold 100907 and 95972 particles with disjoint ids and
adjoining time ranges, and only their union (196879 of 200000 launched) is the
real count.

Reading only ``Step#0`` therefore halves the answer. This script concatenates
every set. Ids must be unique across them; a repeat means a genuine second
crossing, so it is reported rather than silently deduplicated. Use ``--step`` to
restrict the plot to one set.

Usage
-----
    python plot_monitors.py <run_dir>
    python plot_monitors.py <run_dir> --step 0 --outdir /tmp/mon
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from opalxruns.opalx_run import Run, h5_steps

NEEDED = ("x", "y", "z", "px", "py", "pz")


def read_plane(path, step=None):
    """Read one monitor file, concatenating its ``Step#N`` sets.

    Parameters
    ----------
    path : Path
        The monitor ``.h5``.
    step : int, optional
        Read only this set (index into the sorted sets) instead of all of them.

    Returns
    -------
    dict or None
        ``x, y, z, px, py, pz`` arrays plus ``spos`` (mean over the sets read),
        ``nsets``, ``nread`` and any ``warnings``; ``None`` when the file holds
        no usable particles.
    """
    steps = h5_steps(path)
    if not steps:
        return None
    if step is not None:
        try:
            steps = [steps[step]]
        except IndexError:
            return None

    parts = {n: [] for n in NEEDED}
    ids, sposes = [], []
    nsets = 0
    with h5py.File(path, "r") as f:
        nsets = len([k for k in f if k.startswith("Step#")])
        for key in steps:
            g = f[key]
            if not all(n in g for n in NEEDED) or len(g["x"]) == 0:
                continue
            for n in NEEDED:
                parts[n].append(np.asarray(g[n][()], dtype=float).reshape(-1))
            if "id" in g:
                ids.append(np.asarray(g["id"][()], dtype=np.int64).reshape(-1))
            if "SPOS" in g.attrs:
                sposes.append(float(np.ravel(g.attrs["SPOS"])[0]))

    if not parts["x"]:
        return None

    d = {n: np.concatenate(v) for n, v in parts.items()}
    d["spos"] = float(np.mean(sposes)) if sposes else float("nan")
    d["nsets"] = nsets
    d["nread"] = len(steps)
    d["warnings"] = []

    # The sets partition one pass, so an id must not appear twice. A repeat is a
    # real second crossing (or a stale re-run appended to the file) -- say so,
    # never quietly drop it.
    if ids:
        all_ids = np.concatenate(ids)
        uniq, counts = np.unique(all_ids, return_counts=True)
        if len(uniq) != len(all_ids):
            d["warnings"].append(
                f"{len(all_ids) - len(uniq)} duplicate id(s) across "
                f"{len(steps)} Step group(s) (e.g. id {int(uniq[counts > 1][0])})")
    return d


def geom_emittance(u, up):
    """RMS geometric emittance ``sqrt(<u^2><up^2> - <u up>^2)``."""
    u = u - u.mean()
    up = up - up.mean()
    return float(np.sqrt(max(0.0, np.mean(u * u) * np.mean(up * up)
                             - np.mean(u * up) ** 2)))


def plot_plane(name, d, outdir):
    """Write ``<outdir>/<name>.png`` and return the plane's summary row."""
    x = d["x"] * 1e3          # m -> mm
    y = d["y"] * 1e3
    pz = d["pz"]
    safe = np.where(np.abs(pz) > 1e-12, pz, np.nan)
    xp = d["px"] / safe * 1e3  # rad -> mrad
    yp = d["py"] / safe * 1e3
    n = len(x)

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    for a, (u, v, xl, yl, tt) in zip(ax, [
            (x, xp, "x [mm]", "x' [mrad]", "x - x'"),
            (y, yp, "y [mm]", "y' [mrad]", "y - y'"),
            (x, y,  "x [mm]", "y [mm]",    "x - y")]):
        a.scatter(u, v, s=2, alpha=0.35, edgecolors="none", color="#1f5fd1")
        a.set_xlabel(xl)
        a.set_ylabel(yl)
        a.set_title(tt)
        a.grid(True, alpha=0.25)
    sp = d["spos"]
    head = f"{name}   (s = {sp:.3f} m, N = {n})" if np.isfinite(sp) \
        else f"{name}   (N = {n})"
    if d["nsets"] > 1:
        head += f"   [{d['nread']} of {d['nsets']} time-partitioned sets]"
    fig.suptitle(head, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(outdir, name + ".png"), dpi=110)
    plt.close(fig)

    return {
        "monitor": name,
        "s_m": sp,
        "N": n,
        "rms_x_mm": float(np.std(x)),
        "rms_y_mm": float(np.std(y)),
        "emit_x_mm_mrad": geom_emittance(x, xp),
        "emit_y_mm_mrad": geom_emittance(y, yp),
    }


def plot_envelope(rows, outdir, title):
    """RMS size and surviving count against s, over all planes."""
    s = [r["s_m"] for r in rows]
    fig, ax = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    ax[0].plot(s, [r["rms_x_mm"] for r in rows], "o-", label="rms x", color="#1f5fd1")
    ax[0].plot(s, [r["rms_y_mm"] for r in rows], "s-", label="rms y", color="#c1432d")
    ax[0].set_ylabel("rms size [mm]")
    ax[0].legend()
    ax[0].grid(True, alpha=0.25)
    ax[0].set_title(title)
    ax[1].plot(s, [r["N"] for r in rows], "o-", color="#1e8a4c")
    ax[1].set_ylabel("N surviving")
    ax[1].set_xlabel("s [m]")
    ax[1].grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "envelope.png"), dpi=120)
    plt.close(fig)


def write_stats_csv(rows, outdir):
    cols = ["monitor", "s_m", "N", "rms_x_mm", "rms_y_mm",
            "emit_x_mm_mrad", "emit_y_mm_mrad"]
    with open(os.path.join(outdir, "stats.csv"), "w") as o:
        o.write(",".join(cols) + "\n")
        for r in rows:
            vals = []
            for c in cols:
                v = r[c]
                if c == "monitor":
                    vals.append(v)
                elif isinstance(v, float) and np.isnan(v):
                    vals.append("")
                else:
                    vals.append("%.6g" % v)
            o.write(",".join(vals) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?", default=".",
                    help="OPALX run directory (default: cwd)")
    ap.add_argument("--base", default=None,
                    help="run basename, if the directory holds more than one deck")
    ap.add_argument("--step", type=int, default=None,
                    help="plot only this Step#N set (default: concatenate all, "
                         "which is what one pass through the plane amounts to)")
    ap.add_argument("--outdir", default=None,
                    help="output directory (default: <run>/plots/monitors)")
    args = ap.parse_args(argv)

    run = Run(args.run_dir, base=args.base)
    files = run.monitor_h5
    if not files:
        print(f"plot_monitors: no monitor .h5 in {run.dir}", file=sys.stderr)
        return 1

    outdir = Path(args.outdir) if args.outdir else run.plots_dir / "monitors"
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for path in files:
        name = path.stem
        d = read_plane(path, args.step)
        if d is None or len(d["x"]) == 0:
            print(f"  skip {name:<20s} (no particles)")
            continue
        for w in d["warnings"]:
            print(f"  WARNING {name}: {w}")
        row = plot_plane(name, d, outdir)
        rows.append(row)
        print("  %-20s s=%7.3f  N=%6d  rms_x=%6.1f  rms_y=%6.1f mm"
              % (name, row["s_m"], row["N"], row["rms_x_mm"], row["rms_y_mm"]))

    if not rows:
        print("plot_monitors: no monitor file had particles", file=sys.stderr)
        return 1

    rows.sort(key=lambda r: (np.inf if not np.isfinite(r["s_m"]) else r["s_m"]))
    write_stats_csv(rows, outdir)
    plot_envelope(rows, outdir,
                  f"{run.base} — beam envelope & transmission at monitor planes")

    print(f"  wrote {len(rows)} plane plots + envelope.png + stats.csv to {outdir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
