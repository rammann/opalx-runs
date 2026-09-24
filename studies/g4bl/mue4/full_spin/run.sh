#!/usr/bin/env bash
# Run the whole muE4 line in both codes, with spin (G4beamline integrates 12 components). About an hour, G4beamline being the slow one.
# Output goes to output/g4bl/mue4/full_spin/. The comparison reads the G4beamline planes
# from g4bl_reference/, which is tracked because this run is slow:
#   ./run.sh                      run both codes
#   ./run.sh --update-reference   also copy the new G4beamline planes there
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-/opt/homebrew/Caskroom/miniconda/base/bin/python}"
# opalx: $OPALX, else the workspace build /Users/rammann/Code/OPALX/build/src/opalx
OPALX_BIN="$("$PY" -m opalxruns.paths --opalx)" || exit 2
echo "opalx: $OPALX_BIN"
OUT="$("$PY" -m opalxruns.paths "$HERE")"
"$PY" -m opalxruns.run "$HERE/full_spin.g4bl" "$HERE/full_spin.in"
echo "   $(ls "$OUT"/Z*.txt 2>/dev/null | wc -l | tr -d ' ') G4beamline planes, $(ls "$OUT"/E*_PL*.h5 2>/dev/null | wc -l | tr -d ' ') OPALX planes"
if [[ "${1:-}" == "--update-reference" ]]; then
    rm -f "$HERE"/g4bl_reference/Z*.txt
    cp "$OUT"/Z*.txt "$HERE"/g4bl_reference/
    echo "   copied the G4beamline planes into g4bl_reference/"
fi
