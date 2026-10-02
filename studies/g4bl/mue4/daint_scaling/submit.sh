#!/usr/bin/env bash
#
# submit.sh -- submit the muE4 strong scaling runs on daint: the same input on
#              1, 2, 4 and 8 GPUs, one rank per GPU. Run on daint.
#
#     ./submit.sh                       # 1, 2, 4, 8 GPUs, the whole line
#     GPUS=1 ZSTOP=1.0 ./submit.sh      # quick test: 1 GPU, first metre only
#
#   OPALX    the CUDA opalx executable (default: $SCRATCH/build-cuda/src/opalx)
#   UENV     the uenv image the runs use (default: the opal-x GH200 image below)
#   ACCOUNT  Slurm account (default: c41)
#   TIME     time limit per job (default: 00:30:00)
#   GPUS     GPU counts to run (default: "1 2 4 8")
#   ZSTOP    if set, replaces ZSTOP in the copied input (metres)
#
# Every call writes a new folder output/g4bl/mue4/daint_scaling/runs_<date>_<time>/
# with one subfolder gpus_<n> per run. collect.sh reads the timings from there.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
REL="${HERE#"$REPO"/}"

OPALX="${OPALX:-$SCRATCH/build-cuda/src/opalx}"
UENV="${UENV:-/capstor/store/cscs/cscs/public/uenvs/opal-x-gh200-mpich-gcc-2026-10-01.squashfs}"
ACCOUNT="${ACCOUNT:-c41}"
TIME="${TIME:-00:30:00}"
GPUS="${GPUS:-1 2 4 8}"
ZSTOP="${ZSTOP:-}"
GPUS_PER_NODE=4

beam="$(sed -n 's/.*FNAME = "\([^"]*\)".*/\1/p' "$HERE/scaling.in")"

# Check everything before submitting anything.
[[ -x "$OPALX" ]] || { echo "ERROR: OPALX=$OPALX is not an executable" >&2; exit 1; }
[[ -f "$HERE/$beam" ]] || { echo "ERROR: $beam missing; run ./make_beam.sh first" >&2; exit 1; }
[[ -d "$REPO/maps" ]] || { echo "ERROR: $REPO/maps missing" >&2; exit 1; }
# grep -c, not grep -q: with pipefail, grep -q can make strings exit with SIGPIPE.
if command -v strings > /dev/null \
        && [[ "$(strings "$OPALX" | grep -cx BSCALE || true)" -eq 0 ]]; then
    echo "ERROR: $OPALX has no FIELDMAP element with BSCALE; build general-fieldmap-element" >&2
    exit 1
fi

RUN_BASE="$REPO/output/${REL#studies/}/runs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUN_BASE"
echo "runs in $RUN_BASE"

export OPALX
for g in $GPUS; do
    nodes=$(( (g + GPUS_PER_NODE - 1) / GPUS_PER_NODE ))
    per_node=$(( g < GPUS_PER_NODE ? g : GPUS_PER_NODE ))
    if (( nodes * per_node != g )); then
        echo "ERROR: $g GPUs do not fill whole nodes of $GPUS_PER_NODE" >&2
        exit 1
    fi

    run_dir="$RUN_BASE/gpus_$g"
    mkdir -p "$run_dir"
    cp "$HERE/scaling.in" "$HERE/lattice_scaling.in" "$run_dir/"
    ln -s "$REPO/maps" "$run_dir/maps"
    ln -s "$HERE/$beam" "$run_dir/$beam"
    if [[ -n "$ZSTOP" ]]; then
        # -i.bak works with both GNU and BSD sed.
        sed -i.bak "s/ZSTOP = {[^}]*}/ZSTOP = {$ZSTOP}/" "$run_dir/scaling.in"
        rm "$run_dir/scaling.in.bak"
    fi

    job_id="$(sbatch --parsable \
        --job-name="mue4_${g}gpu" \
        --output="$run_dir/slurm_%j.out" \
        --error="$run_dir/slurm_%j.err" \
        --time="$TIME" \
        --partition=normal \
        --account="$ACCOUNT" \
        --nodes="$nodes" \
        --ntasks-per-node="$per_node" \
        --gpus-per-task=1 \
        --cpus-per-task=8 \
        --exclusive \
        --uenv="$UENV" \
        --view=default \
        --chdir="$run_dir" \
        --export=ALL \
        "$HERE/job.sh")"
    echo "  $g GPU(s) on $nodes node(s): job $job_id"
done
