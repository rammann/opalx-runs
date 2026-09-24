"""Where things are. The only place in the Python code that knows machine paths.

Environment variables override the defaults:
  OPALX       the opalx executable (default: build/src/opalx in the workspace, the
              folder that holds this repo)
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

# The workspace build: /Users/rammann/Code/OPALX/build, configured from opalx/.
OPALX_DEFAULT = REPO.parent / "build" / "src" / "opalx"


def output_dir(folder) -> Path:
    """The folder under output/ that belongs to a folder under studies/."""
    return OUTPUT / Path(folder).resolve().relative_to(STUDIES)


def g4bl() -> Path:
    """The g4bl executable inside the G4beamline app."""
    return G4BL_APP / "Contents" / "MacOS" / "g4bl"


def main(argv=None) -> int:
    """python -m opalxruns.paths <folder> ...   print the output/ folder of each
    python -m opalxruns.paths --opalx          print the opalx executable runs use"""
    import sys

    argv = argv if argv is not None else sys.argv[1:]
    if argv == ["--opalx"]:
        try:
            print(opalx())
        except RuntimeError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1
        return 0
    for folder in argv:
        print(output_dir(folder))
    return 0


def opalx() -> Path:
    """The opalx executable: $OPALX if it is set, else the workspace build.

    Rebuild that with `cd /Users/rammann/Code/OPALX/build && make -j8 opalx_exe` --
    `make opalx` builds only the library and leaves a stale executable in place.
    """
    value = os.environ.get("OPALX", "")
    if not value:
        if not (OPALX_DEFAULT.is_file() and os.access(OPALX_DEFAULT, os.X_OK)):
            raise RuntimeError(f"OPALX is not set and the workspace build {OPALX_DEFAULT} "
                               "is missing; build it or set OPALX to an opalx executable")
        return OPALX_DEFAULT
    path = Path(value).expanduser()
    if not (path.is_file() and os.access(path, os.X_OK)):
        raise RuntimeError(f"OPALX={value} is not an executable file")
    return path


if __name__ == "__main__":
    raise SystemExit(main())
