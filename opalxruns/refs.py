"""Which files an input file names, and which of them are missing.

OPALX input (.in): the values of FMAPFN, FNAME and CALL's FILE. A CALLed file is
read too. OUTFN and FILE_NAME name files OPALX writes, so they do not count.
G4beamline input (.g4bl): file= and filename= on fieldmap and beam lines;
fieldntuple writes its file.

Comments are removed first (// and /* */ in OPALX, # in G4beamline), so a path
in a comment does not count. Both codes open these paths from the working
folder, not from the folder of the file that names them, so every path is
checked against one working folder (by default the input's own folder).

    python -m opalxruns.refs [folder ...]    # missing files of every input under studies/
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

OPALX_READS = re.compile(r'\b(?:FMAPFN|FNAME|FILE)\s*=\s*"([^"]*)"', re.IGNORECASE)
G4BL_READS = re.compile(r'\b(?:file|filename)=(?:"([^"]*)"|(\S+))')
G4BL_READING_COMMANDS = {"fieldmap", "beam"}


def _opalx_text(path: Path) -> str:
    text = re.sub(r"/\*.*?\*/", " ", path.read_text(errors="replace"), flags=re.DOTALL)
    return re.sub(r"//[^\n]*", "", text)


def _opalx_refs(path: Path, work_dir: Path, seen: set[Path]) -> list[str]:
    out = []
    for value in OPALX_READS.findall(_opalx_text(path)):
        out.append(value)
        called = work_dir / value
        if value.lower().endswith(".in") and called.is_file() and called.resolve() not in seen:
            seen.add(called.resolve())
            out += _opalx_refs(called, work_dir, seen)
    return out


def _g4bl_refs(path: Path) -> list[str]:
    out = []
    for line in path.read_text(errors="replace").splitlines():
        line = re.sub(r"(^|\s)#.*", "", line)
        words = line.split()
        if words and words[0] in G4BL_READING_COMMANDS:
            out += [quoted or bare for quoted, bare in G4BL_READS.findall(line)]
    return out


def references(input_file: Path, work_dir: Path | None = None) -> list[str]:
    """Every path the input names, as written, in order, without repeats."""
    input_file = Path(input_file)
    work_dir = Path(work_dir) if work_dir else input_file.parent
    if input_file.suffix == ".g4bl":
        found = _g4bl_refs(input_file)
    else:
        found = _opalx_refs(input_file, work_dir, {input_file.resolve()})
    return list(dict.fromkeys(found))


def missing(input_file: Path, work_dir: Path | None = None) -> list[str]:
    """The references that do not exist, looked up from the working folder."""
    input_file = Path(input_file)
    work_dir = Path(work_dir) if work_dir else input_file.parent
    return [r for r in references(input_file, work_dir) if not (work_dir / r).exists()]


def main(argv=None) -> int:
    from opalxruns.paths import REPO

    folders = [Path(a) for a in (argv if argv is not None else sys.argv[1:])] or [REPO / "studies"]
    inputs = sorted(p for f in folders for p in f.rglob("*") if p.suffix in (".in", ".g4bl"))
    n = 0
    for path in inputs:
        for ref in missing(path):
            print(f"{path.relative_to(REPO) if path.is_relative_to(REPO) else path}: {ref}")
            n += 1
    print(f"{n} missing file(s) named by {len(inputs)} input file(s)")
    return 1 if n else 0


if __name__ == "__main__":
    raise SystemExit(main())
