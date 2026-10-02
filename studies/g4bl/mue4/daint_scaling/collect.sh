#!/usr/bin/env bash
#
# collect.sh -- print the timings of one submit.sh call as a table.
#
#     ./collect.sh                                    # the newest runs_* folder
#     ./collect.sh ../../../../output/g4bl/mue4/daint_scaling/runs_20261002_120000
#
# Times are the "Wall max" column of timing.dat (the slowest rank), in seconds.
#   main     mainTimer: the whole run, including reading the maps and the particle
#            file, which every rank does in full and which does not get faster
#   steps    the timers that run every time step, added up:
#            TIntegration1 + TIntegration2 + External field eval + computeMoments
#            + updateParticle
#   field    External field eval
#   push     TIntegration1 + TIntegration2
#   moments  computeMoments. Runs every step even with FIELDSOLVER TYPE = NONE,
#            because computeSpaceChargeFields() calls calcBeamParameters() first.
#   speedup  steps(fewest GPUs) / steps(n) * (fewest GPUs)
#   eff      speedup / n
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
REL="${HERE#"$REPO"/}"

if [[ $# -ge 1 ]]; then
    RUN_BASE="$1"
else
    runs=("$REPO/output/${REL#studies/}"/runs_*)
    [[ -d "${runs[${#runs[@]}-1]}" ]] || { echo "ERROR: no runs_* folder found" >&2; exit 1; }
    RUN_BASE="${runs[${#runs[@]}-1]}"                      # the names sort by date and time
fi
echo "$RUN_BASE"

# Wall max of one timer: the last line starting with its name (timing.dat lists
# "Wall tot" first, then "Wall max / min / avg"); Wall max is third from the end.
# 0 if the timer is not there.
wall_max() {
    grep "^$2" "$1" | tail -1 | awk '{ v = $(NF - 2) } END { print (v == "" ? 0 : v) }'
}

for d in "$RUN_BASE"/gpus_*; do
    g="${d##*_}"
    t="$d/timing.dat"
    if [[ ! -f "$t" ]]; then
        echo "$g"
        continue
    fi
    echo "$g $(wall_max "$t" mainTimer) $(wall_max "$t" TIntegration1)" \
         "$(wall_max "$t" TIntegration2) $(wall_max "$t" "External field eval")" \
         "$(wall_max "$t" computeMoments) $(wall_max "$t" updateParticle)"
done | sort -n | awk '
    BEGIN {
        printf "%5s %9s %9s %9s %9s %9s %8s %6s\n",
               "GPUs", "main", "steps", "field", "push", "moments", "speedup", "eff"
    }
    NF < 7 { printf "%5s  no timing.dat (still running or failed)\n", $1; next }
    {
        push = $3 + $4; field = $5; moments = $6
        steps = push + field + moments + $7
        if (!g0) { g0 = $1; t0 = steps }
        s = t0 / steps * g0
        printf "%5d %9.1f %9.1f %9.1f %9.1f %9.1f %8.2f %6.2f\n",
               $1, $2, steps, field, push, moments, s, s / $1
    }'
