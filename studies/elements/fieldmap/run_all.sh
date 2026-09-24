#!/usr/bin/env bash
#
# run_all.sh -- write the G4beamline field maps and inputs, track them all
#               through OPALX, then run the comparison tests.
#
#     ./run_all.sh                        # generate + track + test
#     ./run_all.sh --test-only            # skip tracking, just re-run the analysis
#     (the opalx binary is $OPALX, or /Users/rammann/Code/OPALX/build/src/opalx if unset)
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# numpy/h5py/matplotlib live in the conda base env, not the system python3.
PY="${PY:-/opt/homebrew/Caskroom/miniconda/base/bin/python}"

TEST_ONLY=0
[[ "${1:-}" == "--test-only" ]] && TEST_ONLY=1
cd "$HERE"

if [[ $TEST_ONLY -eq 0 ]]; then
    # The binary: $OPALX if it is set, else the workspace build,
    # /Users/rammann/Code/OPALX/build/src/opalx. opalxruns.paths decides, so this
    # script and the runner use the same file.
    OPALX_BIN="$("$PY" -m opalxruns.paths --opalx)" || exit 2
    echo "opalx: $OPALX_BIN"

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
        echo "       Rebuild with: cd $(dirname "$(dirname "$OPALX_BIN")") && make -j8 opalx_exe" >&2
        exit 2
    fi

    echo "== writing the field maps and inputs =="
    "$PY" make_inputs.py

    echo "== tracking every case through OPALX (one rank, no space charge) =="
    OUT="$("$PY" -m opalxruns.paths "$HERE")"
    for d in */; do
        name="$(basename "$d")"
        [[ -f "$d/$name.in" ]] || continue
        printf "  %-20s ... " "$name"
        # Output goes to output/, same folders as here. The runner empties the
        # case's run folder first: a stale .h5 left by a failed run would read as
        # a fresh result and pass.
        if "$PY" -m opalxruns.run "$d/$name.in" > /dev/null; then
            echo "done"
        else
            echo "FAILED (see $OUT/$name/run.log)"
        fi
    done
fi

echo "== running the comparison tests =="
"$PY" run_tests.py
