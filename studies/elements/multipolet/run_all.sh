#!/usr/bin/env bash
#
# run_all.sh -- generate the MULTIPOLET bend-study decks, track them all through
#               OPALX, then run the analytic-comparison tests.
#     ./run_all.sh                        # generate + track + test
#     ./run_all.sh --test-only            # skip tracking, just re-run the analysis
#     (the opalx binary is $OPALX, or /Users/rammann/Code/OPALX/build/src/opalx if unset)
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# numpy/h5py/scipy/matplotlib live in the conda base env, not the system python3.
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

    echo "== generating decks =="
    "$PY" make_decks.py

    echo
    echo "== tracking all cases through OPALX (one rank, no space charge) =="
    SKIP="$("$PY" -c 'import json;print(" ".join(m["name"] for m in json.load(open("cases.json")) if m.get("skip")))')"
    OUT="$("$PY" -m opalxruns.paths "$HERE")"
    for d in */; do
        name="$(basename "$d")"
        [[ -f "$d/$name.in" ]] || continue
        if [[ " $SKIP " == *" $name "* ]]; then
            printf "  %-20s ... skipped (hangs OPALX, see its README)\n" "$name"
            continue
        fi
        printf "  %-20s ... " "$name"
        # output goes to output/, same folders as here
        "$PY" -m opalxruns.run "$d/$name.in" > /dev/null \
            && echo "done" || echo "FAILED (see $OUT/$name/run.log)"
    done
fi

echo
echo "== running the validation tests =="
"$PY" run_tests.py
