#!/usr/bin/env bash
#
# run_all.sh -- write the G4beamline field maps and inputs, track them all
#               through OPALX, then run the comparison tests.
#
#     ./run_all.sh              # generate + track + test
#     ./run_all.sh --test-only  # skip tracking, just re-run the analysis
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"          # /Users/rammann/Code/OPALX

# Find the opalx binary. Build trees move, and a stale hardcoded path is the
# worst failure mode here: the run aborts, the old output stays on disk and
# reads as a passing result. Override with OPALX_BIN=... ./run_all.sh
if [[ -z "${OPALX_BIN:-}" ]]; then
    for cand in "$ROOT/opalx/build_serial/src/opalx" "$ROOT/build/src/opalx" \
                "$ROOT/opalx/build/src/opalx"; do
        [[ -x "$cand" ]] && { OPALX_BIN="$cand"; break; }
    done
fi
if [[ -z "${OPALX_BIN:-}" || ! -x "$OPALX_BIN" ]]; then
    echo "ERROR: no opalx binary found. Set OPALX_BIN=/path/to/opalx" >&2
    exit 2
fi

# This study needs the FIELDMAP element and the 3D grid reader, which are newer
# than some build trees lying around. Checking here turns a confusing parse
# error in 14 run.log files into one message. Note the executable target is
# opalx_exe -- `make opalx` builds only the static library and silently leaves
# a stale executable in place.
# grep -c rather than grep -q: with `set -o pipefail`, grep -q exits as soon as
# it matches, strings takes SIGPIPE, and the pipeline reports failure on the
# success path.
have_fieldmap="$(strings "$OPALX_BIN" | grep -c "a FIELDMAP element takes no L" || true)"
if [[ "$have_fieldmap" -eq 0 ]]; then
    echo "ERROR: $OPALX_BIN has no FIELDMAP element." >&2
    echo "       Rebuild with: cd <build dir> && make -j8 opalx_exe" >&2
    exit 2
fi

# numpy/h5py/matplotlib live in the conda base env, not the system python3.
PY="${PY:-/opt/homebrew/Caskroom/miniconda/base/bin/python}"

TEST_ONLY=0
[[ "${1:-}" == "--test-only" ]] && TEST_ONLY=1
cd "$HERE"

if [[ $TEST_ONLY -eq 0 ]]; then
    echo "== writing the field maps and inputs =="
    "$PY" make_inputs.py

    echo "== tracking every case through OPALX (one rank, no space charge) =="
    for d in */; do
        name="$(basename "$d")"
        [[ -f "$d/$name.in" ]] || continue
        # Clear the previous output first: a stale .h5 left by a failed run
        # reads as a fresh result and passes.
        rm -f "$d/$name.h5" "$d/$name.stat" "$d"/*.h5 "$d"/*.stat
        printf "  %-20s ... " "$name"
        if ( cd "$d" && mpirun -n 1 "$OPALX_BIN" "$name.in" --info 1 > run.log 2>&1 ); then
            echo "done"
        else
            echo "FAILED (see $d/run.log)"
        fi
    done
fi

echo "== running the comparison tests =="
"$PY" run_tests.py
