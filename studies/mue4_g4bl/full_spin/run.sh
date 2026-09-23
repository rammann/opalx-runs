#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export PATH="/Users/rammann/Code/G4BL/G4beamline-3.08.app/Contents/MacOS:$PATH"
rm -f Z*.txt E*.h5 full_spin.h5 full_spin.stat full_spin.lbal timing.dat
rm -rf data
echo "== G4beamline (spin tracking, 12-component integration)"
S=$(date +%s); g4bl full_spin.g4bl > g4bl.log 2>&1
echo "   $(( $(date +%s) - S )) s, $(ls Z*.txt 2>/dev/null | wc -l | tr -d ' ') planes"
echo "== OPALX"
S=$(date +%s); mpirun -n 1 /Users/rammann/Code/OPALX/build/src/opalx full_spin.in --info 1 > run.log 2>&1
echo "   $(( $(date +%s) - S )) s, $(ls E*_PL*.h5 2>/dev/null | wc -l | tr -d ' ') planes"
