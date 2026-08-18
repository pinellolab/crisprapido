#!/usr/bin/env bash
set -euo pipefail
SOURCE_PATH="${BASH_SOURCE[0]}"
if [[ "${SOURCE_PATH}" != /* ]]; then
  SOURCE_PATH="${PWD}/${SOURCE_PATH}"
fi
PACKAGE_DIR="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)"
SCRIPT_DIR="${PACKAGE_DIR}"
DRY_RUN="${DRY_RUN:-1}"
RUN_ID_PREFIX="${RUN_ID_PREFIX:-$(basename "${SCRIPT_DIR}")_$(date -u +%Y%m%dT%H%M%SZ)}"
SBATCH_ARGS=()
[[ -n "${SBATCH_PARTITION:-}" ]] && SBATCH_ARGS+=(--partition "${SBATCH_PARTITION}")
[[ -n "${SBATCH_ACCOUNT:-}" ]] && SBATCH_ARGS+=(--account "${SBATCH_ACCOUNT}")
[[ -n "${SBATCH_QOS:-}" ]] && SBATCH_ARGS+=(--qos "${SBATCH_QOS}")
[[ -n "${SBATCH_EXTRA:-}" ]] && read -r -a EXTRA_ARGS <<< "${SBATCH_EXTRA}" && SBATCH_ARGS+=("${EXTRA_ARGS[@]}")
commands=(
  "RUN_ID=${RUN_ID_PREFIX}_correctness sbatch ${SCRIPT_DIR}/slurm_correctness_baseline.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_correctness sbatch ${SCRIPT_DIR}/slurm_correctness_columba.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_timing ITERATION=measured_1 sbatch ${SCRIPT_DIR}/slurm_measured_baseline.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_timing ITERATION=measured_1 sbatch ${SCRIPT_DIR}/slurm_measured_columba.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_timing ITERATION=measured_2 sbatch ${SCRIPT_DIR}/slurm_measured_columba.sbatch"
  "RUN_ID=${RUN_ID_PREFIX}_timing ITERATION=measured_3 sbatch ${SCRIPT_DIR}/slurm_measured_columba.sbatch"
)
for cmd in "${commands[@]}"; do
  if [[ ${#SBATCH_ARGS[@]} -gt 0 ]]; then
    rendered="${cmd/sbatch/sbatch $(printf '%q ' "${SBATCH_ARGS[@]}")}"
  else
    rendered="${cmd}"
  fi
  if [[ "${DRY_RUN}" == "1" ]]; then
    echo "${rendered}"
  else
    eval "${rendered}"
  fi
done
