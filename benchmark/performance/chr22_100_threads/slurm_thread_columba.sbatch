#!/usr/bin/env bash
#SBATCH --job-name=chr22_100_threads_columba
#SBATCH --output=/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/crisprapido/benchmark/performance/chr22_100_threads/slurm_logs/%x_%A_%a.out
#SBATCH --error=/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/crisprapido/benchmark/performance/chr22_100_threads/slurm_logs/%x_%A_%a.err
#SBATCH --array=1-10%3
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
PACKAGE_DIR="${REPO_ROOT}/benchmark/performance/chr22_100_threads"
[[ -d "${PACKAGE_DIR}" ]] || { echo "Missing benchmark package: ${PACKAGE_DIR}" >&2; exit 1; }
LOG_DIR="${PACKAGE_DIR}/slurm_logs"
RAW_DIR="${PACKAGE_DIR}/raw"
case "${PACKAGE_DIR}:${LOG_DIR}:${RAW_DIR}" in
  *:/var/spool/*|/var/spool/*|*:/tmp/*|/tmp/*)
    echo "Unsafe Slurm runtime path resolved under /var/spool or /tmp" >&2
    echo "REPO_ROOT=${REPO_ROOT}" >&2
    echo "PACKAGE_DIR=${PACKAGE_DIR}" >&2
    echo "LOG_DIR=${LOG_DIR}" >&2
    echo "RAW_DIR=${RAW_DIR}" >&2
    exit 1
    ;;
esac
mkdir -p "${LOG_DIR}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is not set}}"
if [[ "${THREADS}" != "${SLURM_CPUS_PER_TASK}" ]]; then
  echo "THREADS (${THREADS}) must match SLURM_CPUS_PER_TASK (${SLURM_CPUS_PER_TASK})" >&2
  exit 2
fi
if [[ "${PATH_DEBUG:-0}" == "1" ]]; then
  printf 'REPO_ROOT=%s\n' "${REPO_ROOT}"
  printf 'PACKAGE_DIR=%s\n' "${PACKAGE_DIR}"
  printf 'LOG_DIR=%s\n' "${LOG_DIR}"
  printf 'RAW_DIR=%s\n' "${RAW_DIR}"
  printf 'RUN_ROOT=%s\n' "${RUN_ROOT:-${RAW_DIR}/${RUN_ID:-slurm_debug}}"
  printf 'THREADS=%s\n' "${THREADS}"
  printf 'SLURM_CPUS_PER_TASK=%s\n' "${SLURM_CPUS_PER_TASK}"
  exit 0
fi
BATCH_ID="batch_$(printf '%03d' "${SLURM_ARRAY_TASK_ID}")"
export BATCH_ID MODE="columba" PHASE="timing" ITERATION="${ITERATION:-measured_1}" THREADS PACKAGE_DIR REPO_ROOT
export RUN_ROOT="${RUN_ROOT:-${RAW_DIR}/${RUN_ID:-slurm_threads_${THREADS}_$(date -u +%Y%m%dT%H%M%SZ)}}"
exec "${PACKAGE_DIR}/run_batch.sh"
