#!/bin/bash -l
#
# job.sh -- one OPALX run, started by submit.sh through sbatch. The run folder,
# node count and GPU count come from the sbatch options; OPALX from the environment.
#
set -euo pipefail

export MPICH_GPU_SUPPORT_ENABLED=1

echo "============================================================"
echo "muE4 strong scaling run"
echo "Job ID:            ${SLURM_JOB_ID}"
echo "Nodes:             ${SLURM_JOB_NUM_NODES}"
echo "Ranks (GPUs):      ${SLURM_NTASKS}"
echo "Working directory: $(pwd)"
echo "Node list:         ${SLURM_JOB_NODELIST}"
echo "opalx:             ${OPALX}"
echo "============================================================"

srun "${OPALX}" scaling.in --info 2
