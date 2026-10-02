#!/usr/bin/env bash
#
# make_beam.sh -- write the particle file scaling.in reads: ../full/parts.txt
#                 repeated until it holds n_particles particles.
#
#     ./make_beam.sh          # run once on daint, in the opalx-runs clone
#
# The count and the file name are read from scaling.in, so change N_PART in
# make_scaling.py, not here. 16e6 particles make a file of about 1.5 GB.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$HERE/../full/parts.txt"

n_part="$(sed -n 's/^REAL n_particles = \([0-9]*\);.*/\1/p' "$HERE/scaling.in")"
fname="$(sed -n 's/.*FNAME = "\([^"]*\)".*/\1/p' "$HERE/scaling.in")"
n_base="$(head -1 "$BASE" | tr -d '[:space:]')"

if [[ -z "$n_part" || -z "$fname" || -z "$n_base" ]]; then
    echo "ERROR: could not read n_particles / FNAME from scaling.in or the count from $BASE" >&2
    exit 1
fi
if (( n_part % n_base != 0 )); then
    echo "ERROR: n_particles = $n_part is not a multiple of $n_base" >&2
    exit 1
fi
copies=$(( n_part / n_base ))

out="$HERE/$fname"
if [[ -f "$out" && "$(head -1 "$out" | tr -d '[:space:]')" == "$n_part" ]]; then
    echo "$fname already holds $n_part particles"
    exit 0
fi

echo "writing $fname: $n_base particles x $copies = $n_part"
# Header: the count, then the column names. Then the particle rows of parts.txt
# (everything after its two header lines), $copies times.
{
    echo "$n_part"
    echo "x y z px py pz"
    awk -v c="$copies" 'NR > 2 { row[NR] = $0 }
                        END { for (i = 0; i < c; i++) for (j = 3; j <= NR; j++) print row[j] }' "$BASE"
} > "$out.tmp"
mv "$out.tmp" "$out"
echo "done: $(du -h "$out" | cut -f1)"
