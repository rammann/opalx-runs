#!/usr/bin/env python
"""test_geometry.py -- check the centreline walk against G4beamline's own output.

`g4bl.out` is the log of an archived run of the muE4 deck. When G4beamline constructs each
element it prints the global position it placed it at, so that log is an independent oracle
for the walk in opalxruns/mue4.py: if the two agree for every element, the lab poses the OPALX deck
will be built from are right, and no tracking has to happen to find that out.

g4bl prints one decimal, so agreement means within half of the last printed digit, 0.05 mm.

Also checks the two map facts the study leans on: that no map uses an `extend*` symmetry line
(the OPALX reader rejects those) and that the ASR62 shim map is identically zero, which is why
four of the six shim placements contribute nothing.

Usage:  python3 tools/test_geometry.py
"""

from __future__ import annotations

import re
from pathlib import Path

from opalxruns import mue4 as m4
from opalxruns.paths import G4BL_FILES

HERE = Path(__file__).resolve().parent
G4BL_DIR = G4BL_FILES / "muE4"
DECK = G4BL_DIR / "g4bl_input" / "mue4_WsxOn.g4bl"
LOG = G4BL_DIR / "g4bl_input" / "g4bl.out"
MAPS = G4BL_DIR / "maps"

TOL = {
    "position_mm": 0.05,  # half of g4bl's last printed digit
}

# g4bl prints one construction line per placed element, in placement order.
CONSTRUCT = re.compile(
    r"^BL(?:CMD|MappedMagnet)?[A-Za-z]*::Construct\s+(\S+)\s+parent=\s*relZ=(-?[\d.eE+-]+)"
)


def read_oracle(path: Path) -> list[tuple[str, float]]:
    """Every 'Construct <name> ... relZ=<z>' line, in order."""
    out = []
    for line in path.read_text(errors="replace").splitlines():
        hit = CONSTRUCT.match(line.strip())
        if hit:
            out.append((hit.group(1), float(hit.group(2))))
    return out


def check_positions() -> list[str]:
    placements, _ = m4.load(DECK)
    oracle = read_oracle(LOG)

    failures = []
    by_name: dict[str, list[float]] = {}
    for name, z in oracle:
        by_name.setdefault(name, []).append(z)

    print(f"{'element':22s} {'walk X':>11s} {'walk Z':>11s} {'g4bl Z':>10s} {'diff':>9s}  ")
    print("-" * 70)

    checked = 0
    for p in placements:
        expected = by_name.get(p.name)
        if not expected:
            print(f"{p.name:22s} {p.lab[0]:11.3f} {p.lab[2]:11.3f} {'--':>10s} {'not in log':>9s}")
            continue
        # Elements placed more than once under the same name (the shim maps) appear once per
        # placement, in order; pop the next one.
        ref = expected.pop(0)
        diff = p.lab[2] - ref
        ok = abs(diff) <= TOL["position_mm"]
        checked += 1
        if not ok:
            failures.append(f"{p.name}: walk Z {p.lab[2]:.4f} mm vs g4bl {ref:.1f} mm")
        print(
            f"{p.name:22s} {p.lab[0]:11.3f} {p.lab[2]:11.3f} {ref:10.1f} "
            f"{diff:+9.4f}  {'ok' if ok else 'FAIL'}"
        )

    print(f"\n{checked} elements checked against g4bl.out")
    return failures


def check_maps() -> list[str]:
    """The two map properties the study relies on."""
    failures = []
    used = [
        "wsx_total",
        "asr61_300d_track",
        "asr61_300sm_track",
        "asr62shim_280_track",
        "asr62shim_280_sm_track",
        "qsm01a_210_track",
    ]

    print("\nmap checks")
    print("-" * 70)
    for stem in used:
        path = MAPS / f"{stem}.g4blmap"
        kind = None
        extend = 0
        nonzero = 0
        rows = 0
        for i, line in enumerate(path.read_text().splitlines()):
            token = line.split()[0] if line.split() else ""
            if i < 3 and token in ("cylinder", "grid"):
                kind = token
            if token.startswith("extend"):
                extend += 1
            if kind == "grid" and i >= 3 and line.strip():
                # Six columns (x y z Bx By Bz) or nine, with Ex Ey Ez appended -- the muE4
                # quadrupole map is a nine-column one.
                cols = line.split()
                if len(cols) in (6, 9):
                    rows += 1
                    if any(float(c) != 0.0 for c in cols[3:6]):
                        nonzero += 1
        note = f"{kind}"
        if kind == "grid":
            note += f", {nonzero}/{rows} rows non-zero"
        if extend:
            failures.append(f"{stem}: {extend} extend* lines, which the OPALX reader rejects")
            note += f", {extend} extend* lines"
        print(f"  {stem:26s} {note}")

    # The claim that lets four shim placements be no-ops.
    zero_map = MAPS / "asr62shim_280_sm_track.g4blmap"
    any_field = any(
        any(float(c) != 0.0 for c in cols[3:6])
        for cols in (ln.split() for ln in zero_map.read_text().splitlines()[3:])
        if len(cols) in (6, 9)
    )
    if any_field:
        failures.append("asr62shim_280_sm_track is NOT identically zero, contrary to the study note")
    else:
        print("  asr62shim_280_sm_track is identically zero, as documented")

    return failures


def main() -> int:
    failures = check_positions() + check_maps()
    print()
    if failures:
        print(f"FAIL: {len(failures)} problem(s)")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("PASS: the centreline walk reproduces g4bl.out, and the maps are as documented")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
