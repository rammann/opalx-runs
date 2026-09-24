"""Run OPALX or G4beamline so that the output lands in output/, never in studies/.

A run folder under output/ (same path as the case folder under studies/) gets links
to the input file and to every file the input names, and the code runs with the run
folder as its working directory. OPALX names .h5, .stat and checkpoint files after
the input path it is given (OpalData::getInputBasename) and writes data/ and
timing.dat into the working directory, so all of it lands in the run folder.
G4beamline writes into the working directory too.

    python -m opalxruns.run studies/elements/collimator/circle/circle.in [--np 2]
    python -m opalxruns.run <input> --keep -- --restart <stem>_checkpoint.h5
    python -m opalxruns.run <case>.g4bl <case>.in       # both codes, one run folder

A run started this way clears its run folder first (its files, data/, plots/ and
paraview/; folders of other runs inside it stay), because a stale .h5 left by a
failed run would otherwise read as a fresh result. --keep leaves it as it is.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from opalxruns import paths, refs

# Folders a run writes itself (OPALX's data/, process_run's plots/ and paraview/).
# Emptying a run folder clears these; any other folder inside it belongs to another
# run (a g4bl/elements stage, a spin dt_* run) and is left alone.
RUN_FOLDERS = {"data", "plots", "paraview", "monitor_plots"}


def _inside(path: Path, root: Path) -> bool:
    return Path(os.path.normpath(path)).is_relative_to(Path(os.path.normpath(root)))


def _link(link: Path, target: Path) -> None:
    if link.is_symlink() and link.resolve() == target.resolve():
        return
    if link.exists() or link.is_symlink():
        raise FileExistsError(f"{link} exists and is not a link to {target}")
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(target)


def prepare(run_dir, inputs, ref_dir=None, keep=False, root=None) -> Path:
    """Make ``run_dir`` ready for a run of ``inputs`` and return it.

    The folder is emptied first unless ``keep``: its files and the folders in
    RUN_FOLDERS go, other folders (other runs) stay. Each input is linked into it under
    its file name, and so is every file the inputs name (see opalxruns.refs), at the
    same relative path, looked up from ``ref_dir`` (default: the input's own folder).
    A path that already reaches the file from the run folder is not linked; since
    output/ has the same depth as studies/, that is the case for shared/.

    Nothing is ever emptied or written outside ``root`` (default: output/), and a
    file the inputs name but that does not exist is an error.
    """
    root = Path(root or paths.OUTPUT)
    run_dir = Path(os.path.abspath(run_dir))
    if not _inside(run_dir, root) or Path(os.path.normpath(run_dir)) == Path(os.path.normpath(root)):
        raise ValueError(f"run folder {run_dir} is outside {root}")

    links, missing = [], []
    for inp in (Path(os.path.abspath(i)) for i in inputs):
        work = Path(os.path.abspath(ref_dir)) if ref_dir else inp.parent
        links.append((run_dir / inp.name, inp))
        for ref in refs.references(inp, work):
            target = Path(os.path.normpath(work / ref))
            if not target.exists():
                missing.append(f"{ref} (named by {inp.name}, looked up in {work})")
                continue
            if Path(ref).is_absolute():
                continue
            mirror = Path(os.path.normpath(run_dir / ref))
            if mirror == target:
                continue
            if not _inside(mirror, root):
                raise ValueError(f"{ref} named by {inp.name} would need a link outside {root}")
            links.append((mirror, target))
    if missing:
        raise FileNotFoundError("missing: " + "; ".join(missing))

    if run_dir.exists() and not keep:
        for entry in run_dir.iterdir():
            if entry.is_dir() and not entry.is_symlink():
                if entry.name in RUN_FOLDERS:
                    shutil.rmtree(entry)
            else:
                entry.unlink()
    run_dir.mkdir(parents=True, exist_ok=True)
    for link, target in links:
        _link(link, target)
    return run_dir


def run_opalx(run_dir, name: str, np: int = 1, args=(), info: int = 1, log: str = "run.log") -> int:
    """Run OPALX on the input ``name`` (a file in the run folder) with the run folder
    as working directory; the output of the run goes to ``log`` there."""
    cmd = ["mpirun", "-n", str(np), str(paths.opalx()), name, "--info", str(info), *args]
    with open(Path(run_dir) / log, "w") as f:
        return subprocess.run(cmd, cwd=run_dir, stdout=f, stderr=subprocess.STDOUT).returncode


def run_g4bl(run_dir, name: str, log: str = "g4bl.log") -> int:
    """Run G4beamline on ``name`` (a file in the run folder) with the run folder as
    working directory; its output goes to ``log`` there."""
    with open(Path(run_dir) / log, "w") as f:
        return subprocess.run([str(paths.g4bl()), name], cwd=run_dir,
                              stdout=f, stderr=subprocess.STDOUT).returncode


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Run OPALX (.in) or G4beamline (.g4bl) inputs in one run folder under "
                    "output/, in the order given. Arguments after -- go to OPALX.")
    ap.add_argument("inputs", type=Path, nargs="+", metavar="input",
                    help="OPALX .in or G4beamline .g4bl files under studies/")
    ap.add_argument("--np", type=int, default=1, help="MPI ranks for OPALX (default 1)")
    ap.add_argument("--keep", action="store_true", help="do not empty the run folder first")
    ap.add_argument("--run-dir", type=Path,
                    help="run folder (default: output/ + the first input's folder)")
    ap.add_argument("--ref-dir", type=Path, help="folder the inputs' files are looked up from")
    if argv is None:
        argv = sys.argv[1:]
    extra = []
    if "--" in argv:
        i = argv.index("--")
        argv, extra = argv[:i], argv[i + 1:]
    a = ap.parse_args(argv)

    inputs = [p.resolve() for p in a.inputs]
    # Find both codes before anything is emptied: a missing binary must not cost the
    # last good output, or an hour of G4beamline before OPALX is found missing.
    if any(i.suffix != ".g4bl" for i in inputs):
        paths.opalx()
    if any(i.suffix == ".g4bl" for i in inputs) and not paths.g4bl().is_file():
        raise RuntimeError(f"no G4beamline at {paths.g4bl()}; set G4BL_APP")
    run_dir = a.run_dir or paths.output_dir(inputs[0].parent)
    prepare(run_dir, inputs, ref_dir=a.ref_dir, keep=a.keep)
    print(f"output: {run_dir}")
    for inp in inputs:
        start = time.time()
        if inp.suffix == ".g4bl":
            rc, log = run_g4bl(run_dir, inp.name), "g4bl.log"
        else:
            rc, log = run_opalx(run_dir, inp.name, np=a.np, args=extra), "run.log"
        print(f"  {inp.name:24s} exit {rc}  {time.time() - start:7.1f} s  log: {log}")
        if rc != 0:
            return rc
    return 0


if __name__ == "__main__":
    sys.exit(main())
