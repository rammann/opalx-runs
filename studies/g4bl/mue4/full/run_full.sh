#!/usr/bin/env bash
# Run the whole muE4 line in both codes. About an hour, G4beamline being the slow one.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
G4BL_APP="${G4BL_APP:-/Users/rammann/Code/G4BL/G4beamline-3.08.app}"
export PATH="$G4BL_APP/Contents/MacOS:$PATH"
: "${OPALX:?set OPALX to the opalx executable}"
rm -f Z*.txt E*.h5 full.h5 full.stat full.lbal timing.dat
rm -rf data
echo "== G4beamline"
S=$(date +%s); g4bl full.g4bl > g4bl.log 2>&1
echo "   $(( $(date +%s) - S )) s, $(ls Z*.txt 2>/dev/null | wc -l | tr -d ' ') planes"
echo "== OPALX"
S=$(date +%s); mpirun -n 1 "$OPALX" full.in --info 1 > run.log 2>&1
echo "   $(( $(date +%s) - S )) s, $(ls E*_PL*.h5 2>/dev/null | wc -l | tr -d ' ') planes"
