#!/usr/bin/env bash
set -euo pipefail
SOURCE_PATH="${BASH_SOURCE[0]}"
if [[ "${SOURCE_PATH}" != /* ]]; then SOURCE_PATH="${PWD}/${SOURCE_PATH}"; fi
PACKAGE_DIR="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)"
DRY_RUN="${DRY_RUN:-1}"
RUN_ID_PREFIX="${RUN_ID_PREFIX:-chr22_100_threads_tux05_$(date -u +%Y%m%dT%H%M%SZ)}"
THREAD_SET="${THREAD_SET:-1 2 4 8}"
SBATCH_PARTITION="${SBATCH_PARTITION:-tux}"
SBATCH_NODE="${SBATCH_NODE:-tux05}"
SBATCH_COMMON=(-p "${SBATCH_PARTITION}" -w "${SBATCH_NODE}")
[[ -n "${SBATCH_ACCOUNT:-}" ]] && SBATCH_COMMON+=(--account "${SBATCH_ACCOUNT}")
[[ -n "${SBATCH_QOS:-}" ]] && SBATCH_COMMON+=(--qos "${SBATCH_QOS}")
[[ -n "${SBATCH_EXTRA:-}" ]] && read -r -a EXTRA_ARGS <<< "${SBATCH_EXTRA}" && SBATCH_COMMON+=("${EXTRA_ARGS[@]}")
concurrency_for_threads() {
  case "$1" in
    1|2) echo 3 ;;
    4) echo 2 ;;
    *) echo 1 ;;
  esac
}
baseline_mem_for_threads() {
  case "$1" in
    1) echo 1G ;;
    2) echo 2G ;;
    4) echo 4G ;;
    *) echo 8G ;;
  esac
}
columba_mem_for_threads() {
  case "$1" in
    1) echo 3G ;;
    2) echo 4G ;;
    4) echo 6G ;;
    *) echo 10G ;;
  esac
}
for threads in ${THREAD_SET}; do
  cap="$(concurrency_for_threads "${threads}")"
  base_mem="$(baseline_mem_for_threads "${threads}")"
  col_mem="$(columba_mem_for_threads "${threads}")"
  run_id="${RUN_ID_PREFIX}_t${threads}"
  commands=(
    "RUN_ID=${run_id}_smoke_baseline THREADS=${threads} ITERATION=smoke sbatch -c ${threads} --mem=${base_mem} --array=1-1 ${PACKAGE_DIR}/slurm_thread_baseline.sbatch"
    "RUN_ID=${run_id}_smoke_columba THREADS=${threads} ITERATION=smoke sbatch -c ${threads} --mem=${col_mem} --array=1-1 ${PACKAGE_DIR}/slurm_thread_columba.sbatch"
    "RUN_ID=${run_id}_timing THREADS=${threads} ITERATION=measured_1 sbatch -c ${threads} --mem=${base_mem} --array=1-10%${cap} ${PACKAGE_DIR}/slurm_thread_baseline.sbatch"
    "RUN_ID=${run_id}_timing THREADS=${threads} ITERATION=measured_1 sbatch -c ${threads} --mem=${col_mem} --array=1-10%${cap} ${PACKAGE_DIR}/slurm_thread_columba.sbatch"
    "RUN_ID=${run_id}_timing THREADS=${threads} ITERATION=measured_2 sbatch -c ${threads} --mem=${col_mem} --array=1-10%${cap} ${PACKAGE_DIR}/slurm_thread_columba.sbatch"
    "RUN_ID=${run_id}_timing THREADS=${threads} ITERATION=measured_3 sbatch -c ${threads} --mem=${col_mem} --array=1-10%${cap} ${PACKAGE_DIR}/slurm_thread_columba.sbatch"
  )
  for cmd in "${commands[@]}"; do
    rendered="${cmd/sbatch/sbatch $(printf '%q ' "${SBATCH_COMMON[@]}")}"
    if [[ "${DRY_RUN}" == "1" ]]; then echo "${rendered}"; else eval "${rendered}"; fi
  done
done
