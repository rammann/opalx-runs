#!/usr/bin/env bash
#
# End-to-end checkpoint/restart test for OPALX, driven from the muE4 deck.
#
# Three runs, each in its own scratch directory next to this script (so the deck's
# ../../../data/... paths still resolve):
#
#   work_ref  uninterrupted reference run, 1 rank
#   work_r1   run killed mid-track, then restarted from the checkpoint on 1 rank
#   work_r2   run killed mid-track on 1 rank, then restarted on 2 ranks
#
# Two controls, to attribute any failure correctly:
#
#   work_ctl2  fresh uninterrupted 2-rank run, no checkpoint involved -- separates
#              "restart is broken on 2 ranks" from "this deck is broken on 2 ranks"
#   work_r2b   killed on 2 ranks and restarted on 2 ranks -- separates "restart on
#              multiple ranks" from "restart with a *changed* rank count"
#   work_nomon same as r2b but with the 20 MONITORs stripped out -- attributes the
#              multi-rank restart abort to the monitor / LossDataSink path
#
# A restart appends to the killed run's own .stat/.h5 after rewinding them to the
# checkpoint, so each scenario must stay in the directory where it was killed.
#
# The reference and the to-be-killed runs are launched as a bare binary (singleton
# MPI, 1 rank) rather than through mpirun, so the kill is a plain SIGKILL to a known
# pid instead of a fight with the mpirun process tree. Only the 2-rank restart needs
# mpirun. One rank is one rank either way.
#
# Usage:  ./run_checkpoint_test.sh [ref|r1|r2 ...]     (default: all three)
#   OPALX=/path/to/opalx  ./run_checkpoint_test.sh
#   NCHECKPOINTS=3        how many checkpoints to let the killed runs write first

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OPALX="${OPALX:-/Users/rammann/Code/OPALX/build/src/opalx}"
DECK="$HERE/mue4_ckpt/mue4_ckpt.in"
BASE="mue4_ckpt"
CKPT="${BASE}_checkpoint.h5"
NCHECKPOINTS="${NCHECKPOINTS:-3}"

[[ -x "$OPALX" ]] || { echo "opalx not executable: $OPALX" >&2; exit 1; }
[[ -f "$DECK"  ]] || { echo "deck not found: $DECK"        >&2; exit 1; }

# Fresh work dir with the deck copied in. Removing it is what makes a scenario start
# from nothing -- a stale .stat would otherwise be appended to.
prepare() {
    local dir="$HERE/work_$1"
    rm -rf "$dir"
    mkdir -p "$dir"
    cp "$DECK" "$dir/$BASE.in"
    printf '%s' "$dir"
}

# Start a run in the background, return once it has written N checkpoints, then
# SIGKILL it. macOS has no timeout(1), and a fixed sleep lands on a different step on
# every machine; polling the log is deterministic.
run_until_checkpoints_then_kill() {
    local dir="$1" n="$2" nranks="${3:-1}"
    if [[ "$nranks" == "1" ]]; then
        ( cd "$dir" && exec "$OPALX" "$BASE.in" --info 2 > run.log 2>&1 ) &
    else
        ( cd "$dir" && exec mpirun -n "$nranks" "$OPALX" "$BASE.in" --info 2 > run.log 2>&1 ) &
    fi
    local pid=$!

    local waited=0
    while kill -0 "$pid" 2>/dev/null; do
        local got
        got=$(grep -c "Wrote checkpoint" "$dir/run.log" 2>/dev/null || true)
        if [[ "${got:-0}" -ge "$n" ]]; then
            # Under mpirun the pid is the launcher, so kill the ranks too.
            [[ "$nranks" == "1" ]] || pkill -9 -f "$OPALX $BASE.in" 2>/dev/null
            kill -9 "$pid" 2>/dev/null
            wait "$pid" 2>/dev/null
            echo "  SIGKILL after $n checkpoints ($(echo "$waited * 0.2" | bc)s of tracking)"
            return 0
        fi
        # A Release build reaches the third checkpoint in seconds, so poll finely --
        # at 1 s granularity the run can finish before the kill lands.
        sleep 0.2
        waited=$((waited + 1))
        if [[ $waited -gt 3000 ]]; then
            echo "  TIMEOUT waiting for $n checkpoints" >&2
            kill -9 "$pid" 2>/dev/null
            return 1
        fi
    done

    # Finished before writing n checkpoints: CHECKPOINTFREQ is too coarse, or the beam
    # died early. Either way there is no mid-track state to restart from.
    echo "  ERROR: run ended before $n checkpoints were written" >&2
    return 1
}

# Freeze the killed state so the comparison can separate what the restart inherited
# from the crash from what it rebuilt.
snapshot_kill_state() {
    local dir="$1"
    grep "Wrote checkpoint" "$dir/run.log" | tail -1 \
        | sed -E 's/.*after global step ([0-9]+).*/\1/' > "$dir/kill_step.txt"
    wc -l < "$dir/$BASE.stat" | tr -d ' ' > "$dir/kill_stat_lines.txt"
    mkdir -p "$dir/killed_state"
    cp "$dir/$CKPT" "$dir/killed_state/" 2>/dev/null
    for f in "$dir"/MON_*.h5; do
        [[ -e "$f" ]] && cp "$f" "$dir/killed_state/"
    done
    echo "  last checkpoint: global step $(cat "$dir/kill_step.txt"); .stat at kill: $(cat "$dir/kill_stat_lines.txt") lines"
}

scenario_ref() {
    echo "== work_ref: uninterrupted reference, 1 rank =="
    local dir rc
    dir=$(prepare ref)
    ( cd "$dir" && exec "$OPALX" "$BASE.in" --info 2 > run.log 2>&1 )
    rc=$?
    echo "  exit $rc; $(grep -c '^' "$dir/$BASE.stat" 2>/dev/null || echo 0) .stat lines; \
$(grep -c 'Wrote checkpoint' "$dir/run.log" 2>/dev/null || echo 0) checkpoints written"
}

scenario_ctl2() {
    echo "== work_ctl2: fresh uninterrupted run on 2 ranks, no checkpoint restart =="
    local dir rc
    dir=$(prepare ctl2)
    ( cd "$dir" && exec mpirun -n 2 "$OPALX" "$BASE.in" --info 2 > run.log 2>&1 )
    rc=$?
    echo "  exit $rc; $(grep -c '^' "$dir/$BASE.stat" 2>/dev/null || echo 0) .stat lines"
}

# Same deck with every MONITOR definition and every reference to one in the LINE
# removed. Everything else -- solenoid, bends, quads, scraping apertures -- stays.
strip_monitors() {
    python3 - "$1" <<'PY'
import re, sys
p = sys.argv[1]
s = open(p).read()
s = "\n".join(l for l in s.split("\n") if not re.match(r"\s*mon_\d+_\w+\s*:\s*MONITOR", l))
s = re.sub(r"\bmon_\d+_\w+\s*,\s*", "", s)          # LINE entries with a following comma
s = re.sub(r",(\s*)\bmon_\d+_\w+\b", r"\1", s)      # a trailing LINE entry
s = re.sub(r",(\s*)\)", r"\1)", s)                  # any comma left dangling before ')'
open(p, "w").write(s)
assert "MONITOR" not in s.split("*/", 1)[1], "monitors survived the strip"
PY
}

scenario_nomon() {
    echo "== work_nomon: r2b repeated with the MONITORs removed =="
    local dir rc
    dir=$(prepare nomon)
    strip_monitors "$dir/$BASE.in" || { echo "  could not strip monitors" >&2; return 1; }

    run_until_checkpoints_then_kill "$dir" "$NCHECKPOINTS" 2 || return 1
    snapshot_kill_state "$dir"
    ( cd "$dir" && exec mpirun -n 2 "$OPALX" "$BASE.in" \
        --restart "$CKPT" --info 2 > restart.log 2>&1 )
    rc=$?
    echo "  restart exit $rc  (0 here means the abort in work_r2b comes from the monitors)"
}

scenario_restart() {
    local tag="$1" nranks="$2" killranks="${3:-1}"
    echo "== work_$tag: kill mid-track on $killranks rank(s), restart on $nranks rank(s) =="
    local dir rc
    dir=$(prepare "$tag")

    run_until_checkpoints_then_kill "$dir" "$NCHECKPOINTS" "$killranks" || return 1
    snapshot_kill_state "$dir"

    # Deck first: Main.cpp only accepts the positional input file at argv[1] or last.
    if [[ "$nranks" == "1" ]]; then
        ( cd "$dir" && exec "$OPALX" "$BASE.in" --restart "$CKPT" --info 2 > restart.log 2>&1 )
    else
        ( cd "$dir" && exec mpirun -n "$nranks" "$OPALX" "$BASE.in" \
            --restart "$CKPT" --info 2 > restart.log 2>&1 )
    fi
    rc=$?
    echo "  restart exit $rc"
    grep -E "Restored checkpoint|does not yet support|Exception|error message" \
        "$dir/restart.log" | head -5
}

targets=("$@")
[[ ${#targets[@]} -eq 0 ]] && targets=(ref r1 r2 ctl2 r2b nomon)

for t in "${targets[@]}"; do
    case "$t" in
        ref)   scenario_ref ;;
        r1)    scenario_restart r1  1 1 ;;
        r2)    scenario_restart r2  2 1 ;;
        r2b)   scenario_restart r2b 2 2 ;;
        ctl2)  scenario_ctl2 ;;
        nomon) scenario_nomon ;;
        *)     echo "unknown scenario: $t" >&2; exit 1 ;;
    esac
done

echo
echo "Now run:  /opt/homebrew/Caskroom/miniconda/base/bin/python compare.py"
