#!/usr/bin/env python
"""process_run.py -- run every post-processing step on one OPALX run directory.

Works out what the directory contains, then runs the steps whose inputs are
present. A step with no inputs is skipped with a note; a step that fails is
reported and the rest still run, so one broken output does not cost you the
others.

Steps
-----
=========  ===================================  ================================
step       needs                                writes
=========  ===================================  ================================
stat       ``<base>.stat``                      ``plots/<base>_overview.png``
timing     ``timing.dat``                       ``plots/timing.png``
monitors   ``MON_*.h5``                         ``plots/monitors/``
elements   deck + ``data/*_ElementPositions``   ``paraview/<base>_elements.vtp``
                                                ``paraview/<base>_reforbit.vtp``
particles  ``<base>.h5``                        ``paraview/bunch_{lab,comoving}.pvd``
=========  ===================================  ================================

Usage
-----
    python process_run.py <run_dir>
    python process_run.py <run_dir> --only stat,timing
    python process_run.py <run_dir> --skip particles
    python process_run.py <run_dir> --frame lab --stride 5
    python process_run.py <run_dir> --dry-run

Run with a Python that has numpy/h5py/matplotlib/vtk (miniconda base), not the
system ``python3``.
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from opalx_run import Run  # noqa: E402

STEPS = ("stat", "timing", "monitors", "elements", "particles")


# --------------------------------------------------------------------------- #
# steps -- each returns None, or a reason string when its inputs are missing   #
# --------------------------------------------------------------------------- #
def _why_skip(run, step):
    if step == "stat" and not run.stat_files:
        return "no .stat file"
    if step == "timing" and run.timing is None:
        return "no timing.dat"
    if step == "monitors" and not run.monitor_h5:
        return "no monitor .h5"
    if step == "elements":
        if run.element_positions is None:
            return "no data/*_ElementPositions.txt"
        if run.deck is None:
            return "no input deck"
    if step == "particles" and run.bunch_h5 is None:
        return "no phase-space .h5"
    return None


def _run_stat(run, args):
    import plot_stat

    argv = [str(run.dir), "--vs", args.vs]
    if args.base:
        argv += ["--base", args.base]
    if args.columns:
        argv += ["--columns", args.columns]
    return plot_stat.main(argv)


def _run_timing(run, args):
    import plot_timing

    argv = [str(run.dir)]
    if args.base:
        argv += ["--base", args.base]
    return plot_timing.main(argv)


def _run_monitors(run, args):
    import plot_monitors

    argv = [str(run.dir), "--step", str(args.step)]
    if args.base:
        argv += ["--base", args.base]
    return plot_monitors.main(argv)


def _run_elements(run, args):
    import elements_to_vtk

    argv = [str(run.dir), "--default-aperture", str(args.default_aperture)]
    if args.base:
        argv += ["--base", args.base]
    return elements_to_vtk.main(argv)


def _run_particles(run, args):
    import particles_to_vtk

    argv = [str(run.dir), "--frame", args.frame, "--stride", str(args.stride),
            "--mass", str(args.mass), "--quiet"]
    if args.base:
        argv += ["--base", args.base]
    return particles_to_vtk.main(argv)


RUNNERS = {
    "stat": _run_stat,
    "timing": _run_timing,
    "monitors": _run_monitors,
    "elements": _run_elements,
    "particles": _run_particles,
}


# --------------------------------------------------------------------------- #
def _parse_step_list(value, flag):
    names = [s.strip() for s in value.split(",") if s.strip()]
    bad = [n for n in names if n not in STEPS]
    if bad:
        raise SystemExit(f"{flag}: unknown step(s) {', '.join(bad)} "
                         f"(pick from {', '.join(STEPS)})")
    return names


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", nargs="?", default=".",
                    help="OPALX run directory (default: cwd)")
    ap.add_argument("--base", default=None,
                    help="run basename, if the directory holds more than one deck")
    ap.add_argument("--only", default=None,
                    help=f"comma-separated subset of: {', '.join(STEPS)}")
    ap.add_argument("--skip", default=None,
                    help="comma-separated steps to leave out")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what was found and what would run, then stop")

    g = ap.add_argument_group("step options")
    g.add_argument("--vs", default="s", choices=("s", "t"),
                   help="stat: x axis (default: s)")
    g.add_argument("--columns", default=None,
                   help="stat: extra columns for a second figure")
    g.add_argument("--step", type=int, default=-1,
                   help="monitors: which set to plot, -1 = newest (default: -1)")
    g.add_argument("--default-aperture", type=float, default=0.05,
                   help="elements: tube radius [m] with no explicit aperture "
                        "(default: 0.05)")
    g.add_argument("--frame", choices=("lab", "comoving", "both"), default="both",
                   help="particles: frame(s) to write (default: both)")
    g.add_argument("--stride", type=int, default=1,
                   help="particles: write every Nth dump (default: 1)")
    g.add_argument("--mass", type=float, default=0.1056583755,
                   help="particles: species rest mass [GeV] for Ekin_MeV "
                        "(default: muon)")
    args = ap.parse_args(argv)

    run = Run(args.run_dir, base=args.base)

    wanted = _parse_step_list(args.only, "--only") if args.only else list(STEPS)
    if args.skip:
        skipped = _parse_step_list(args.skip, "--skip")
        wanted = [s for s in wanted if s not in skipped]

    print(f"run: {run.dir}")
    for label, value in run.inventory():
        print(f"  {label:<12s} {value}")
    print()

    plan = []
    for step in STEPS:
        if step not in wanted:
            plan.append((step, "not requested"))
            continue
        plan.append((step, _why_skip(run, step)))

    print("steps:")
    for step, reason in plan:
        print(f"  {step:<10s} {'skip — ' + reason if reason else 'run'}")
    print()

    if args.dry_run:
        return 0

    failed, ran = [], []
    for step, reason in plan:
        if reason:
            continue
        print(f"--- {step} ---")
        try:
            rc = RUNNERS[step](run, args)
        except Exception:  # noqa: BLE001 - one broken step must not stop the rest
            traceback.print_exc()
            rc = 1
        if rc:
            failed.append(step)
        else:
            ran.append(step)
        print()

    print("summary")
    print(f"  ran:     {', '.join(ran) if ran else '-'}")
    print(f"  skipped: {', '.join(s for s, r in plan if r) or '-'}")
    if failed:
        print(f"  FAILED:  {', '.join(failed)}")
    outputs = [d for d in (run.dir / "plots", run.dir / "paraview") if d.is_dir()]
    for d in outputs:
        print(f"  output:  {d}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
