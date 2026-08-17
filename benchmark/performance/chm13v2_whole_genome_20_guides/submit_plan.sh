#!/usr/bin/env bash
set -euo pipefail

SOURCE_PATH="${BASH_SOURCE[0]}"
if [[ "${SOURCE_PATH}" != /* ]]; then SOURCE_PATH="${PWD}/${SOURCE_PATH}"; fi
PACKAGE_DIR="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)"
RUN_ID_PREFIX="${RUN_ID_PREFIX:-chm13_wg20_$(date -u +%Y%m%dT%H%M%SZ)}"
PARTITION="${SBATCH_PARTITION:-tux}"
NODE="${SBATCH_NODE:-tux05}"

printf 'Smoke baseline:\nRUN_ID=%s_smoke sbatch -p %q -w %q --array=1-1%%1 %q\n' "${RUN_ID_PREFIX}" "${PARTITION}" "${NODE}" "${PACKAGE_DIR}/slurm_correctness_baseline.sbatch"
printf 'Smoke Columba:\nRUN_ID=%s_smoke sbatch -p %q -w %q --array=1-1%%1 %q\n' "${RUN_ID_PREFIX}" "${PARTITION}" "${NODE}" "${PACKAGE_DIR}/slurm_correctness_columba.sbatch"
printf 'Full correctness baseline:\nRUN_ID=%s_correctness sbatch -p %q -w %q %q\n' "${RUN_ID_PREFIX}" "${PARTITION}" "${NODE}" "${PACKAGE_DIR}/slurm_correctness_baseline.sbatch"
printf 'Full correctness Columba:\nRUN_ID=%s_correctness sbatch -p %q -w %q %q\n' "${RUN_ID_PREFIX}" "${PARTITION}" "${NODE}" "${PACKAGE_DIR}/slurm_correctness_columba.sbatch"
printf 'Controlled sequential timing plan (dry run):\n'
RUN_ID="${RUN_ID_PREFIX}_controlled_timing" SBATCH_PARTITION="${PARTITION}" SBATCH_NODE="${NODE}" \
  "${PACKAGE_DIR}/submit_controlled_timing.sh"
