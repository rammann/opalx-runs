#!/usr/bin/env python
"""Compare checkpoint/restart runs against the uninterrupted reference.

A restart is correct if the run it produces is indistinguishable from a run that was
never interrupted: same number of .stat rows, a seam with no duplicated or skipped
step, the same final beam, and monitor files holding each particle crossing once.

Usage:
    python compare.py [work_ref] [work_r1 work_r2 ...]

Needs h5py. On this machine:
    /opt/homebrew/Caskroom/miniconda/base/bin/python compare.py
"""

import glob
import os
import re
import sys

import h5py
import numpy as np

BASE = "mue4_ckpt"

# 1 rank replays the identical arithmetic in the identical order, so it should agree
# to round-off. More ranks change the MPI reduction order, which is a real but tiny
# difference -- loosen for those rather than pretend it is exact.
TOL_SAME_RANKS = 1e-10
TOL_DIFF_RANKS = 1e-6

STAT_COLS_CHECKED = ["s", "numParticles", "energy", "rms_x", "rms_y", "rms_s"]


def worst_diff(a, b, tol):
    """Largest |a-b| relative to the scale of the column, not to each element.

    Dividing by |a| element-wise makes a 1e-11 absolute difference on a coordinate
    that happens to sit near zero look like a catastrophic relative error. Scale by
    the column's own magnitude instead, which is what "these two runs agree" means.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    scale = max(float(np.max(np.abs(a))), 1e-300)
    worst = float(np.max(np.abs(a - b))) / scale
    return worst, worst <= tol


def read_stat(path):
    """Read an OPALX SDDS .stat into {column name: array}."""
    names, rows = [], []
    with open(path) as fh:
        in_data = False
        skip = 0
        for line in fh:
            if not in_data:
                m = re.match(r"\s*name=(\w+),", line)
                if m:
                    names.append(m.group(1))
                elif line.startswith("&data"):
                    in_data = True
                    # &data block header, then the parameter values (processors,
                    # revision, species) that precede the first data row.
                    skip = 3 + 3
                continue
            if skip:
                skip -= 1
                continue
            parts = line.split()
            if parts:
                rows.append(parts)

    # The parameters share the namespace with the columns; columns are whatever the
    # data rows actually carry, counted from the end of the name list.
    ncol = len(rows[0])
    cols = names[-ncol:]
    arr = np.array(rows, dtype=float)
    return {c: arr[:, i] for i, c in enumerate(cols)}, arr.shape[0]


def last_h5_step(path):
    """Particles of the final Step# in an OPALX phase-space file, sorted by id."""
    with h5py.File(path, "r") as f:
        steps = sorted(f.keys(), key=lambda k: int(k.split("#")[1]))
        g = f[steps[-1]]
        ids = np.asarray(g["id"])
        order = np.argsort(ids)
        cols = {k: np.asarray(g[k])[order] for k in ("x", "y", "z", "px", "py", "pz")}
        return cols, ids[order], len(steps)


def monitor_counts(d):
    """Number of recorded crossings in each MON_*.h5 in a run directory."""
    out = {}
    for path in sorted(glob.glob(os.path.join(d, "MON_*.h5"))):
        try:
            with h5py.File(path, "r") as f:
                n = sum(len(f[s]["id"]) for s in f.keys() if "id" in f[s])
        except Exception as exc:  # a truncated file is itself the finding
            n = f"UNREADABLE ({exc.__class__.__name__})"
        out[os.path.basename(path)] = n
    return out


def restart_step(d):
    """(step the checkpoint was written at, step the restart resumed from)."""
    wrote = resumed = None
    run_log = os.path.join(d, "run.log")
    if os.path.exists(run_log):
        for line in open(run_log, errors="replace"):
            m = re.search(r"Wrote checkpoint.*after global step (\d+)", line)
            if m:
                wrote = int(m.group(1))
    restart_log = os.path.join(d, "restart.log")
    if os.path.exists(restart_log):
        for line in open(restart_log, errors="replace"):
            m = re.search(r"Restored checkpoint.*at global step (\d+)", line)
            if m:
                resumed = int(m.group(1))
    return wrote, resumed


def restart_failed(d):
    """Any fatal message in the restart log -- MPI abort, exception, refusal."""
    path = os.path.join(d, "restart.log")
    if not os.path.exists(path):
        return None
    pat = re.compile(
        r"does not yet support|MPI_ERR|error occurred in|ERRORS_ARE_FATAL"
        r"|calling \"abort\"|OpalException|terminate called"
    )
    for line in open(path, errors="replace"):
        if pat.search(line):
            return line.strip()
    return None


class Report:
    def __init__(self):
        self.rows = []
        self.failed = False

    def check(self, name, ok, detail=""):
        self.rows.append(("PASS" if ok else "FAIL", name, detail))
        if not ok:
            self.failed = True

    def note(self, name, detail):
        self.rows.append(("--  ", name, detail))

    def show(self, title):
        print(f"\n{title}")
        print("-" * len(title))
        for status, name, detail in self.rows:
            print(f"  [{status}] {name:<34} {detail}")


def compare(ref_dir, run_dir, tol):
    r = Report()

    # A control run was never interrupted, so the restart-only checks do not apply.
    is_restart = os.path.exists(os.path.join(run_dir, "restart.log"))

    if is_restart:
        fatal = restart_failed(run_dir)
        r.check("0 restart ran without aborting", fatal is None, fatal or "")

        wrote, resumed = restart_step(run_dir)
        r.check(
            "1 resumed at the checkpoint step",
            wrote is not None and wrote == resumed,
            f"checkpoint written at {wrote}, restart resumed at {resumed}",
        )
    else:
        r.note("0 not a restart", "control run, restart-only checks skipped")

    ref_stat, ref_n = read_stat(os.path.join(ref_dir, f"{BASE}.stat"))
    run_stat, run_n = read_stat(os.path.join(run_dir, f"{BASE}.stat"))

    # A duplicated seam shows up as s going backwards or repeating; a skipped seam as
    # a jump. Compare the whole s column against the reference instead of guessing a
    # threshold for "a jump".
    s_ref, s_run = ref_stat["s"], run_stat["s"]
    r.check("2 .stat row count matches reference", ref_n == run_n, f"{run_n} vs {ref_n}")
    r.check(
        "3 .stat s column is strictly increasing",
        bool(np.all(np.diff(s_run) > 0)),
        f"{int(np.sum(np.diff(s_run) <= 0))} non-increasing steps",
    )
    if ref_n == run_n:
        for col in STAT_COLS_CHECKED:
            worst, ok = worst_diff(ref_stat[col], run_stat[col], tol)
            r.check(f"4 .stat {col} matches everywhere", ok, f"max diff/scale {worst:.3e}")

    ref_last, ref_ids, ref_steps = last_h5_step(os.path.join(ref_dir, f"{BASE}.h5"))
    run_last, run_ids, run_steps = last_h5_step(os.path.join(run_dir, f"{BASE}.h5"))
    r.check("5 .h5 dump count matches", ref_steps == run_steps, f"{run_steps} vs {ref_steps}")
    same_n = len(ref_last["x"]) == len(run_last["x"])
    r.check(
        "6 final particle count matches",
        same_n,
        f"{len(run_last['x'])} vs {len(ref_last['x'])}",
    )
    if same_n:
        # Particle ids are handed out per rank, so the same physical muon carries a
        # different id in a 1-rank and a 2-rank run. When the id sets disagree, compare
        # the sorted coordinate distributions instead -- that still catches a wrong
        # final beam without pretending the labels line up.
        by_id = np.array_equal(ref_ids, run_ids)
        worst = 0.0
        for k in ref_last:
            a, b = (ref_last[k], run_last[k]) if by_id else (
                np.sort(ref_last[k]), np.sort(run_last[k]))
            worst = max(worst, worst_diff(a, b, tol)[0])
        how = "matched by id" if by_id else "id sets differ, compared as sorted distributions"
        r.check("7 final phase space matches", worst <= tol, f"max diff/scale {worst:.3e} ({how})")

    ref_mon, run_mon = monitor_counts(ref_dir), monitor_counts(run_dir)
    missing = sorted(set(ref_mon) - set(run_mon))
    extra = sorted(set(run_mon) - set(ref_mon))
    r.check(
        "8 same set of monitor files",
        not missing and not extra,
        f"missing {missing} extra {extra}" if (missing or extra) else "",
    )
    bad = {k: (run_mon[k], ref_mon[k]) for k in ref_mon if k in run_mon and run_mon[k] != ref_mon[k]}
    r.check(
        "9 monitor crossing counts match",
        not bad,
        "; ".join(f"{k}: {v[0]} vs {v[1]}" for k, v in list(bad.items())[:6]) if bad else "",
    )
    return r


def main():
    args = sys.argv[1:]
    here = os.path.dirname(os.path.abspath(__file__))
    ref = args[0] if args else os.path.join(here, "work_ref")
    # work_nomon runs a deliberately different deck (monitors stripped), so it is a
    # standalone diagnostic -- its exit code is the result, not a diff against ref.
    runs = args[1:] or sorted(
        d
        for d in glob.glob(os.path.join(here, "work_*"))
        if os.path.basename(d) not in ("work_ref", "work_nomon")
    )

    failed = False
    for run in runs:
        if not os.path.exists(os.path.join(run, f"{BASE}.stat")):
            print(f"\nskipping {os.path.basename(run)}: no {BASE}.stat")
            continue
        # ctl2 is a fresh 2-rank run, not a restart; r2/r2b restart on 2 ranks.
        tol = TOL_SAME_RANKS if os.path.basename(run) == "work_r1" else TOL_DIFF_RANKS
        rep = compare(ref, run, tol)
        rep.show(f"{os.path.basename(run)}  vs  {os.path.basename(ref)}   (tol {tol:g})")
        failed |= rep.failed

    print()
    print("RESULT:", "FAIL" if failed else "PASS")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
