#!/usr/bin/env bash
set -euo pipefail

SOURCE_PATH="${BASH_SOURCE[0]}"
if [[ "${SOURCE_PATH}" != /* ]]; then SOURCE_PATH="${PWD}/${SOURCE_PATH}"; fi
PACKAGE_DIR="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)"
REPO_ROOT="$(cd "${PACKAGE_DIR}/../../.." && pwd)"
RAW_DIR="${PACKAGE_DIR}/raw"

RUN_ID="${RUN_ID:-chm13_wg20_controlled_timing_$(date -u +%Y%m%dT%H%M%SZ)}"
PARTITION="${SBATCH_PARTITION:-tux}"
NODE="${SBATCH_NODE:-tux05}"
SUBMIT="${SUBMIT:-0}"
TIMING_DESIGN="${TIMING_DESIGN:-sequential_batches_and_replicates}"

BASELINE_SCRIPT="${PACKAGE_DIR}/slurm_measured_baseline.sbatch"
COLUMBA_SCRIPT="${PACKAGE_DIR}/slurm_measured_columba.sbatch"
RUN_ROOT="${RAW_DIR}/${RUN_ID}"

[[ "${SUBMIT}" == "0" || "${SUBMIT}" == "1" ]] || {
  echo "SUBMIT must be 0 or 1" >&2
  exit 2
}
[[ -d "${REPO_ROOT}/.git" ]] || {
  echo "Repository root not found: ${REPO_ROOT}" >&2
  exit 2
}
[[ -x "${BASELINE_SCRIPT}" && -x "${COLUMBA_SCRIPT}" ]] || {
  echo "Measured Slurm scripts are missing or not executable" >&2
  exit 2
}
grep -qx '#SBATCH --array=1-10%1' "${BASELINE_SCRIPT}"
grep -qx '#SBATCH --array=1-10%1' "${COLUMBA_SCRIPT}"

base_export="ALL,RUN_ID=${RUN_ID},ITERATION=measured_1,TIMING_DESIGN=${TIMING_DESIGN},CACHE_STATE_LABEL=baseline_reference_scan"
col1_export="ALL,RUN_ID=${RUN_ID},ITERATION=measured_1,TIMING_DESIGN=${TIMING_DESIGN},CACHE_STATE_LABEL=cold_start_candidate"
col2_export="ALL,RUN_ID=${RUN_ID},ITERATION=measured_2,TIMING_DESIGN=${TIMING_DESIGN},CACHE_STATE_LABEL=warm_cache_candidate"
col3_export="ALL,RUN_ID=${RUN_ID},ITERATION=measured_3,TIMING_DESIGN=${TIMING_DESIGN},CACHE_STATE_LABEL=warm_cache_candidate"

print_plan() {
  printf 'cd %q\n' "${REPO_ROOT}"
  printf 'RUN_ID=%q\n' "${RUN_ID}"
  printf 'BASELINE_JOB_ID=$(sbatch --parsable -p %q -w %q --array=1-10%%1 --export=%q %q)\n' \
    "${PARTITION}" "${NODE}" "${base_export}" "${BASELINE_SCRIPT}"
  printf 'BASELINE_JOB_ID=${BASELINE_JOB_ID%%%%;*}\n'
  printf 'COLUMBA_1_JOB_ID=$(sbatch --parsable -p %q -w %q --array=1-10%%1 --dependency=afterok:${BASELINE_JOB_ID} --export=%q %q)\n' \
    "${PARTITION}" "${NODE}" "${col1_export}" "${COLUMBA_SCRIPT}"
  printf 'COLUMBA_1_JOB_ID=${COLUMBA_1_JOB_ID%%%%;*}\n'
  printf 'COLUMBA_2_JOB_ID=$(sbatch --parsable -p %q -w %q --array=1-10%%1 --dependency=afterok:${COLUMBA_1_JOB_ID} --export=%q %q)\n' \
    "${PARTITION}" "${NODE}" "${col2_export}" "${COLUMBA_SCRIPT}"
  printf 'COLUMBA_2_JOB_ID=${COLUMBA_2_JOB_ID%%%%;*}\n'
  printf 'COLUMBA_3_JOB_ID=$(sbatch --parsable -p %q -w %q --array=1-10%%1 --dependency=afterok:${COLUMBA_2_JOB_ID} --export=%q %q)\n' \
    "${PARTITION}" "${NODE}" "${col3_export}" "${COLUMBA_SCRIPT}"
  printf 'COLUMBA_3_JOB_ID=${COLUMBA_3_JOB_ID%%%%;*}\n'
  printf 'printf "baseline=%%s columba_1=%%s columba_2=%%s columba_3=%%s\\n" "${BASELINE_JOB_ID}" "${COLUMBA_1_JOB_ID}" "${COLUMBA_2_JOB_ID}" "${COLUMBA_3_JOB_ID}"\n'
}

if [[ "${SUBMIT}" != "1" ]]; then
  echo "# Dry run only. Set SUBMIT=1 to submit this dependency chain."
  echo "# The first Columba panel is a cold-start candidate; caches are not forcibly cleared."
  print_plan
  exit 0
fi

if [[ -e "${RUN_ROOT}" ]]; then
  echo "Refusing to reuse existing run directory: ${RUN_ROOT}" >&2
  exit 1
fi

cd "${REPO_ROOT}"
baseline_raw="$(sbatch --parsable -p "${PARTITION}" -w "${NODE}" --array=1-10%1 \
  --export="${base_export}" "${BASELINE_SCRIPT}")"
baseline_job="${baseline_raw%%;*}"
columba_1_raw="$(sbatch --parsable -p "${PARTITION}" -w "${NODE}" --array=1-10%1 \
  --dependency="afterok:${baseline_job}" --export="${col1_export}" "${COLUMBA_SCRIPT}")"
columba_1_job="${columba_1_raw%%;*}"
columba_2_raw="$(sbatch --parsable -p "${PARTITION}" -w "${NODE}" --array=1-10%1 \
  --dependency="afterok:${columba_1_job}" --export="${col2_export}" "${COLUMBA_SCRIPT}")"
columba_2_job="${columba_2_raw%%;*}"
columba_3_raw="$(sbatch --parsable -p "${PARTITION}" -w "${NODE}" --array=1-10%1 \
  --dependency="afterok:${columba_2_job}" --export="${col3_export}" "${COLUMBA_SCRIPT}")"
columba_3_job="${columba_3_raw%%;*}"

printf 'RUN_ID=%s\n' "${RUN_ID}"
printf 'baseline=%s\n' "${baseline_job}"
printf 'columba_measured_1=%s dependency=afterok:%s cache_state=cold_start_candidate\n' \
  "${columba_1_job}" "${baseline_job}"
printf 'columba_measured_2=%s dependency=afterok:%s cache_state=warm_cache_candidate\n' \
  "${columba_2_job}" "${columba_1_job}"
printf 'columba_measured_3=%s dependency=afterok:%s cache_state=warm_cache_candidate\n' \
  "${columba_3_job}" "${columba_2_job}"
