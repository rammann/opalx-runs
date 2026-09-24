#!/usr/bin/env bash
#
# run_all.sh -- write the spin cases, run G4beamline and OPALX on each, then
#               compare both with the closed form.
#
#     ./run_all.sh                        # write + run + analyse
#     ./run_all.sh --test-only            # just re-run the analysis
#     (the opalx binary is $OPALX, or /Users/rammann/Code/OPALX/build/src/opalx if unset)
#
# Output goes to output/g4bl/spin/, same folders as here. uniform_bz/ gets one run
# folder per input: uniform_bz (both codes), dt_<dt> (OPALX, one per scan/dt_*.in)
# and bz2m (OPALX over 2 m). The analysis reads the monitor file MON_OUT.h5 and the
# G4beamline plane Z1200.txt from those folders.
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
    OUT="$("$PY" -m opalxruns.paths "$HERE")"

    echo "== writing the cases =="
    "$PY" make_scan.py

    echo "== G4beamline and OPALX on every case =="
    for d in by_p*/; do
        name="$(basename "$d")"
        printf "  %-14s ... " "$name"
        "$PY" -m opalxruns.run "$d/$name.g4bl" "$d/$name.in" > /dev/null \
            && echo "done" || echo "FAILED (see $OUT/$name)"
    done

    # The scan inputs sit in scan/ but read the field map and particles of
    # uniform_bz/, so every input here has its files looked up there.
    U=uniform_bz
    run_case() {   # <run folder name> <input> ...
        local dir="$1"; shift
        printf "  %-14s ... " "$dir"
        "$PY" -m opalxruns.run "$@" --ref-dir "$U" --run-dir "$OUT/$U/$dir" > /dev/null \
            && echo "done" || echo "FAILED (see $OUT/$U/$dir)"
    }
    run_case uniform_bz "$U/uniform_bz.g4bl" "$U/uniform_bz.in"
    for f in "$U"/scan/dt_*.in; do
        run_case "$(basename "$f" .in)" "$f"
    done
    run_case bz2m "$U/bz2m.in"
fi

echo "== comparing with the closed form =="
"$PY" analyse_scan.py
"$PY" analyse_bz.py
