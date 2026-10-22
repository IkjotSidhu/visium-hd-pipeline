#!/bin/bash
# =============================================================================
# LSF job body for ONE pipeline step at ONE resolution (Minerva / Rocky9).
# Normally submitted by run_pipeline.sh, which overrides the resources below
# per step and passes STEP / CONFIG / RES through the environment:
#
#   bsub -env "all, STEP=03_cluster, CONFIG=/abs/config.yaml, RES=8" < submit_step.sh
# =============================================================================
#BSUB -P acc_naiklab
#BSUB -q premium
#BSUB -n 1
#BSUB -W 04:00
#BSUB -R "rusage[mem=32000]"
#BSUB -R "span[hosts=1]"
#BSUB -J visiumhd_step

# No `set -e`: `module` is a shell function whose nonzero returns would abort the job
set -uo pipefail

PIPELINE_DIR="${PIPELINE_DIR:?PIPELINE_DIR not set}"
echo "Job ${LSB_JOBID:-local} | ${STEP} | res ${RES} | host $(hostname) | $(date)"

# R/4.4.1 Rocky9 module only - the project library is ABI-compiled against it
module purge
module load R/4.4.1

# One BLAS/OpenMP thread per process; parallelism comes from future workers.
# Without this, OpenBLAS spawns a thread per physical core and jobs segfault.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLAS_NUM_THREADS=1

Rscript "${PIPELINE_DIR}/steps/${STEP}.R" "${CONFIG}" "${RES}"
EXIT_CODE=$?
echo "Exit code: ${EXIT_CODE} | $(date)"
exit ${EXIT_CODE}
