"""Where things are. The only place in the Python code that knows machine paths.

Environment variables override the defaults:
  OPALX       the opalx executable (required to run OPALX; there is no default)
  G4BL_APP    the G4beamline app (default: the one in /Users/rammann/Code/G4BL)
  G4BL_FILES  folder with the muE4 G4beamline input and field maps
              (default: g4bl-files/ next to this repo)
"""

from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
STUDIES = REPO / "studies"
SHARED = REPO / "shared"
# Run output, in the same folder tree as studies/. It must be a real folder, not a
# link: inputs reach shared/ through ../, which only works at the same depth.
OUTPUT = REPO / "output"

# Generators write this path into tracked input files, so resolve() keeps the
# spelling they have always used.
G4BL_FILES = Path(os.environ.get("G4BL_FILES", REPO.parent / "g4bl-files")).expanduser().resolve()

G4BL_APP = Path(os.environ.get("G4BL_APP", "/Users/rammann/Code/G4BL/G4beamline-3.08.app")).expanduser()


def output_dir(folder) -> Path:
    """The folder under output/ that belongs to a folder under studies/."""
    return OUTPUT / Path(folder).resolve().relative_to(STUDIES)


def g4bl() -> Path:
    """The g4bl executable inside the G4beamline app."""
    return G4BL_APP / "Contents" / "MacOS" / "g4bl"


def main(argv=None) -> int:
    """python -m opalxruns.paths <folder> ...: print the output/ folder of each."""
    import sys

    for folder in (argv if argv is not None else sys.argv[1:]):
        print(output_dir(folder))
    return 0


def opalx() -> Path:
    """The opalx executable named by $OPALX.

    There is no default on purpose: build trees move, and a stale binary that
    still runs leaves old output on disk that reads as a passing result.
    """
    value = os.environ.get("OPALX", "")
    if not value:
        raise RuntimeError("set OPALX to the opalx executable, "
                           "e.g. export OPALX=/Users/rammann/Code/OPALX/opalx/build_serial/src/opalx")
    path = Path(value).expanduser()
    if not (path.is_file() and os.access(path, os.X_OK)):
        raise RuntimeError(f"OPALX={value} is not an executable file")
    return path


if __name__ == "__main__":
    raise SystemExit(main())
