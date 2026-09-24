#!/usr/bin/env bash
#
# Generate every case, run both codes on it, test, and plot.
#
#   ./run_all.sh                 everything
#   ./run_all.sh --test-only     skip both codes, just re-run the analysis
#   ./run_all.sh --pair-only     skip the 20000-particle stage
#   ./run_all.sh asr61_dipole    one case (repeatable)
#
# Both codes resolve relative paths from the working directory, so every run
# happens inside the case's own folder.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-/opt/homebrew/Caskroom/miniconda/base/bin/python}"
G4BL_APP="${G4BL_APP:-/Users/rammann/Code/G4BL/G4beamline-3.08.app}"

TEST_ONLY=0
PAIR_ONLY=0
ONLY_CASES=()
for arg in "$@"; do
    case "$arg" in
        --test-only) TEST_ONLY=1 ;;
        --pair-only) PAIR_ONLY=1 ;;
        --*) echo "unknown option $arg" >&2; exit 2 ;;
        *) ONLY_CASES+=("$arg") ;;
    esac
done

# --- OPALX and G4beamline --------------------------------------------------
# Only needed to run the codes; --test-only works without either.
if [[ "$TEST_ONLY" -eq 0 ]]; then
    # The binary comes from $OPALX and nowhere else. The CMake target is opalx_exe.
    # `make opalx` builds only the static library and leaves whatever executable
    # was there before, which is how a run silently uses a binary from weeks ago.
    OPALX_BIN="${OPALX:?set OPALX to the opalx executable}"

    # Every deck here places elements with a FIELDMAP, which an older binary rejects
    # at parse time with a message about SCALE. Refuse to run rather than produce a
    # log full of parse errors.
    # grep -c, not grep -q: with `set -o pipefail`, grep -q exits on the first match,
    # strings then dies of SIGPIPE, and the pipeline reports failure on a binary that
    # is perfectly fine. grep -c reads to the end.
    if [[ "$(strings "$OPALX_BIN" | grep -c "a FIELDMAP element takes no L")" -eq 0 ]]; then
        echo "ERROR: $OPALX_BIN has no FIELDMAP element." >&2
        echo "       Rebuild with: cd <build dir> && make -j8 opalx_exe" >&2
        exit 1
    fi
    echo "opalx:  $OPALX_BIN"

    export PATH="$G4BL_APP/Contents/MacOS:$PATH"
    command -v g4bl >/dev/null || { echo "g4bl not on PATH" >&2; exit 1; }
fi
echo "python: $PY"

# --- cases -----------------------------------------------------------------
"$PY" "$HERE/make_cases.py"

# macOS ships bash 3.2, which has no mapfile and errors on expanding an empty
# array under `set -u`, so both are done the long way.
CASES=()
if [[ "${#ONLY_CASES[@]}" -gt 0 ]]; then
    CASES=("${ONLY_CASES[@]}")
else
    while IFS= read -r line; do
        CASES+=("$line")
    done < <("$PY" -c "
import json
print('\n'.join(c['name'] for c in json.load(open('$HERE/cases.json'))))")
fi

OUT="$("$PY" -m opalxruns.paths "$HERE")"

run_one() {   # <case> <stem> <what>: both codes, in output/.../<case>/<what>/
    local c="$1" stem="$2" what="$3"
    echo "  $what"
    # Each stage has its own run folder, because the stages write the same plane
    # and field file names. The runner empties it first.
    "$PY" -m opalxruns.run "$HERE/$c/$stem.g4bl" "$HERE/$c/$stem.in" --run-dir "$OUT/$c/$what" \
        > /dev/null || { echo "    FAILED, see the logs in $OUT/$c/$what" >&2; return 1; }
    # OPALX writes its field dumps into data/; the tests read them next to the planes.
    if [[ "$what" == "pair" ]]; then
        cp -f "$OUT/$c/pair"/data/opalx_field_*.dat "$OUT/$c/pair"/ 2>/dev/null || true
    fi
}

if [[ "$TEST_ONLY" -eq 0 ]]; then
    for c in "${CASES[@]}"; do
        echo "== $c"
        run_one "$c" "$c" pair
        # Cases sensitive enough to resolve below their own step error run the pair
        # stage again at half the step, so run_tests can extrapolate rather than
        # assume the step is fine.
        [[ -f "$HERE/$c/${c}_fine.in" ]] && run_one "$c" "${c}_fine" fine
        [[ "$PAIR_ONLY" -eq 1 ]] || run_one "$c" "${c}_gauss" gauss
    done
fi

"$PY" "$HERE/run_tests.py" "${CASES[@]}"
