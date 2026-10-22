#!/bin/bash
# =============================================================================
# Submit the pipeline as chained LSF jobs: step 01 once, then steps 02-07 for
# every resolution in bin_sizes (each resolution is its own chain, so 8 um,
# 16 um and cells run in parallel). Each job waits for the previous step.
#
# Usage:
#   bash run_pipeline.sh config.yaml                 # all steps, all resolutions
#   bash run_pipeline.sh config.yaml -s 03-05        # only steps 03 to 05
#   bash run_pipeline.sh config.yaml -r 8,cells      # only these resolutions
#   bash run_pipeline.sh config.yaml -n              # dry run: print bsub commands
#   bash run_pipeline.sh config.yaml -t 0.1          # test mode: 10% of each replicate -> output_test/
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$(realpath "${1:?Usage: run_pipeline.sh config.yaml [-s 02-05] [-r 8,16,cells] [-t 0.1] [-n]}")"
shift
STEPS="01-07"; RES_FILTER=""; DRY=0; TEST=""
while getopts "s:r:t:n" opt; do
  case $opt in
    s) STEPS="$OPTARG" ;;
    r) RES_FILTER="$OPTARG" ;;
    t) TEST="$OPTARG" ;;
    n) DRY=1 ;;
    *) exit 1 ;;
  esac
done
FIRST="${STEPS%-*}"; LAST="${STEPS#*-}"
export VHD_TEST_FRACTION="$TEST"   # read by R/utils.R (empty = full data)

# Plan from the config (validated in R): one line per job = step res cores mem walltime queue project logdir
module purge >/dev/null 2>&1; module load R/4.4.1 >/dev/null 2>&1
PLAN="$(Rscript -e "source('${PIPELINE_DIR}/R/utils.R'); cat(plan_jobs(load_config('${CONFIG}'), '${FIRST}', '${LAST}', '${RES_FILTER}'), sep = '\n')")"

declare -A LAST_JOB   # last submitted job id per resolution chain
STEP01_JOB=""
while read -r STEP RES CORES MEM WALL QUEUE PROJECT LOGDIR; do
  [[ -z "$STEP" ]] && continue
  mkdir -p "$LOGDIR"

  # Dependency: previous step of this resolution, or step 01 for the first job in a chain
  DEP="${LAST_JOB[$RES]:-$STEP01_JOB}"
  WAIT=(); [[ -n "$DEP" ]] && WAIT=(-w "done(${DEP})")

  CMD=(bsub -P "$PROJECT" -q "$QUEUE" -n "$CORES" -W "$WALL" -R "rusage[mem=${MEM}]" -R "span[hosts=1]"
       -J "vhd_${STEP}_${RES}" -oo "${LOGDIR}/${STEP}_${RES}.out" -eo "${LOGDIR}/${STEP}_${RES}.err" "${WAIT[@]}"
       -env "all, STEP=${STEP}, CONFIG=${CONFIG}, RES=${RES}, PIPELINE_DIR=${PIPELINE_DIR}${TEST:+, VHD_TEST_FRACTION=${TEST}}")
  if [[ $DRY -eq 1 ]]; then
    echo "${CMD[*]} < ${PIPELINE_DIR}/submit_step.sh"; JOB="DRY_${STEP}_${RES}"
  else
    OUT="$("${CMD[@]}" < "${PIPELINE_DIR}/submit_step.sh" 2>&1)" || true
    JOB="$(grep -oE 'Job <[0-9]+>' <<< "$OUT" | grep -oE '[0-9]+' || true)"
    if [[ -z "$JOB" ]]; then echo "ERROR: bsub failed for ${STEP} (${RES}):" >&2; echo "$OUT" >&2; exit 1; fi
    echo "Submitted ${STEP} (${RES}) as job ${JOB}${DEP:+ after ${DEP}}"
  fi

  if [[ "$STEP" == 01_* ]]; then STEP01_JOB="$JOB"; else LAST_JOB[$RES]="$JOB"; fi
done <<< "$PLAN"
