#!/usr/bin/env bash
#
# copy_maps_to_daint.sh -- copy the six field maps the lattice names into maps/ of
#                          this folder in the opalx-runs clone on daint. The maps
#                          are not in git; everything else comes with the clone.
#                          Run on the laptop.
#
#     ./copy_maps_to_daint.sh
#     DAINT=myalias DAINT_ROOT=/capstor/scratch/cscs/me/opalx-runs ./copy_maps_to_daint.sh
#
#   DAINT       ssh host name of daint (default: daint)
#   DAINT_ROOT  the opalx-runs clone on daint (default: $SCRATCH/opalx-runs there)
#   G4BL_FILES  local folder with muE4/maps (default: g4bl-files next to this repo)
#
# About 165 MB.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
REL="${HERE#"$REPO"/}"                                  # studies/g4bl/mue4/daint_scaling
MAPS="${G4BL_FILES:-$REPO/../g4bl-files}/muE4/maps"
DAINT="${DAINT:-daint}"
if [[ -z "${DAINT_ROOT:-}" ]]; then
    # A login shell (bash -l), because a plain ssh command may not set SCRATCH.
    remote_scratch="$(ssh "$DAINT" 'bash -lc "echo \$SCRATCH"' | tail -1)"
    [[ -n "$remote_scratch" ]] || { echo "ERROR: SCRATCH is empty on $DAINT; set DAINT_ROOT" >&2; exit 1; }
    DAINT_ROOT="$remote_scratch/opalx-runs"
fi

DEST="$DAINT_ROOT/$REL/maps"
echo "to $DAINT:$DEST"

# The maps the lattice names, from FMAPFN = "maps/<file>".
maps=()
while IFS= read -r f; do
    [[ -f "$MAPS/$f" ]] || { echo "ERROR: $MAPS/$f not found" >&2; exit 1; }
    maps+=("$MAPS/$f")
done < <(sed -n 's/.*FMAPFN = "maps\/\([^"]*\)".*/\1/p' "$HERE/lattice_scaling.in" | sort -u)

ssh "$DAINT" "mkdir -p '$DEST'"
rsync -av "${maps[@]}" "$DAINT:$DEST/"
