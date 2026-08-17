#!/usr/bin/env bash
#SBATCH --job-name=chm13_wg20_col_timing
#SBATCH --output=/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/crisprapido/benchmark/performance/chm13v2_whole_genome_20_guides/slurm_logs/%x_%A_%a.out
#SBATCH --error=/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/crisprapido/benchmark/performance/chm13v2_whole_genome_20_guides/slurm_logs/%x_%A_%a.err
#SBATCH --array=1-10%1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
PACKAGE_DIR="${REPO_ROOT}/benchmark/performance/chm13v2_whole_genome_20_guides"
[[ -d "${PACKAGE_DIR}" ]] || exit 1
LOG_DIR="${PACKAGE_DIR}/slurm_logs"; RAW_DIR="${PACKAGE_DIR}/raw"
case "${PACKAGE_DIR}:${LOG_DIR}:${RAW_DIR}" in *:/var/spool/*|/var/spool/*|*:/tmp/*|/tmp/*) exit 1;; esac
mkdir -p "${LOG_DIR}"
if [[ "${PATH_DEBUG:-0}" == "1" ]]; then
  printf 'REPO_ROOT=%s\nPACKAGE_DIR=%s\nLOG_DIR=%s\nRAW_DIR=%s\nRUN_ROOT=%s\n' \
    "${REPO_ROOT}" "${PACKAGE_DIR}" "${LOG_DIR}" "${RAW_DIR}" "${RUN_ROOT:-${RAW_DIR}/${RUN_ID:-debug}}"
  printf 'ITERATION=%s\nTIMING_DESIGN=%s\nCACHE_STATE_LABEL=%s\n' \
    "${ITERATION:-NA}" "${TIMING_DESIGN:-NA}" "${CACHE_STATE_LABEL:-NA}"
  exit 0
fi
export BATCH_ID="batch_$(printf '%03d' "${SLURM_ARRAY_TASK_ID}")"
export MODE=columba PHASE=timing ITERATION="${ITERATION:-measured_1}" PACKAGE_DIR REPO_ROOT
export RUN_ROOT="${RUN_ROOT:-${RAW_DIR}/${RUN_ID:?RUN_ID is required}}"
exec "${PACKAGE_DIR}/run_batch.sh"
