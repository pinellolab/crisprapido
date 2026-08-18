#!/usr/bin/env bash
set -euo pipefail
SOURCE_PATH="${BASH_SOURCE[0]}"
if [[ "${SOURCE_PATH}" != /* ]]; then
  SOURCE_PATH="${PWD}/${SOURCE_PATH}"
fi
PACKAGE_DIR="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)"
DRY_RUN="${DRY_RUN:-1}"
RUN_ID_PREFIX="${RUN_ID_PREFIX:-$(basename "${PACKAGE_DIR}")_tux05_$(date -u +%Y%m%dT%H%M%SZ)}"
SBATCH_PARTITION="${SBATCH_PARTITION:-tux}"
SBATCH_NODE="${SBATCH_NODE:-tux05}"
SBATCH_ARGS=(-p "${SBATCH_PARTITION}" -w "${SBATCH_NODE}")
[[ -n "${SBATCH_ACCOUNT:-}" ]] && SBATCH_ARGS+=(--account "${SBATCH_ACCOUNT}")
[[ -n "${SBATCH_QOS:-}" ]] && SBATCH_ARGS+=(--qos "${SBATCH_QOS}")
[[ -n "${SBATCH_EXTRA:-}" ]] && read -r -a EXTRA_ARGS <<< "${SBATCH_EXTRA}" && SBATCH_ARGS+=("${EXTRA_ARGS[@]}")
commands=(
  "RUN_ID=${RUN_ID_PREFIX}_smoke_baseline ITERATION=smoke sbatch --array=1-1 ${PACKAGE_DIR}/slurm_measured_baseline.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_smoke_columba ITERATION=smoke sbatch --array=1-1 ${PACKAGE_DIR}/slurm_measured_columba.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_timing ITERATION=measured_1 sbatch ${PACKAGE_DIR}/slurm_measured_baseline.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_timing ITERATION=measured_1 sbatch ${PACKAGE_DIR}/slurm_measured_columba.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_timing ITERATION=measured_2 sbatch ${PACKAGE_DIR}/slurm_measured_columba.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_timing ITERATION=measured_3 sbatch ${PACKAGE_DIR}/slurm_measured_columba.sbatch"
)
for cmd in "${commands[@]}"; do
  rendered="${cmd/sbatch/sbatch $(printf '%q ' "${SBATCH_ARGS[@]}")}"
  if [[ "${DRY_RUN}" == "1" ]]; then
    echo "${rendered}"
  else
    eval "${rendered}"
  fi
done
