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

# Generators write this path into tracked input files, so resolve() keeps the
# spelling they have always used.
G4BL_FILES = Path(os.environ.get("G4BL_FILES", REPO.parent / "g4bl-files")).expanduser().resolve()

G4BL_APP = Path(os.environ.get("G4BL_APP", "/Users/rammann/Code/G4BL/G4beamline-3.08.app")).expanduser()


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
