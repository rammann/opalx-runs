#!/usr/bin/env bash
# Run the whole muE4 line in both codes. About an hour, G4beamline being the slow one.
# Output goes to output/g4bl/mue4/full/. The comparison reads the G4beamline planes
# from g4bl_reference/, which is tracked because this run is slow:
#   ./run_full.sh                      run both codes
#   ./run_full.sh --update-reference   also copy the new G4beamline planes there
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-/opt/homebrew/Caskroom/miniconda/base/bin/python}"
: "${OPALX:?set OPALX to the opalx executable}"
OUT="$("$PY" -m opalxruns.paths "$HERE")"
"$PY" -m opalxruns.run "$HERE/full.g4bl" "$HERE/full.in"
echo "   $(ls "$OUT"/Z*.txt 2>/dev/null | wc -l | tr -d ' ') G4beamline planes, $(ls "$OUT"/E*_PL*.h5 2>/dev/null | wc -l | tr -d ' ') OPALX planes"
if [[ "${1:-}" == "--update-reference" ]]; then
    rm -f "$HERE"/g4bl_reference/Z*.txt
    cp "$OUT"/Z*.txt "$HERE"/g4bl_reference/
    echo "   copied the G4beamline planes into g4bl_reference/"
fi
