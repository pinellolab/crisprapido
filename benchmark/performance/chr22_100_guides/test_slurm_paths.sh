#!/usr/bin/env bash
set -euo pipefail
SOURCE_PATH="${BASH_SOURCE[0]}"
if [[ "${SOURCE_PATH}" != /* ]]; then
  SOURCE_PATH="${PWD}/${SOURCE_PATH}"
fi
PACKAGE_DIR="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)"
REPO_ROOT="$(cd "${PACKAGE_DIR}/../../.." && pwd)"
OUT="${OUT:-${PACKAGE_DIR}/raw/path_test_$(date -u +%Y%m%dT%H%M%SZ)}"
rm -rf "${OUT}"
mkdir -p "${OUT}"
(
  cd /tmp
  DRY_RUN=1 RUN_ID_PREFIX=path_test "${PACKAGE_DIR}/submit_plan.sh" > "${OUT}/submit_plan.txt"
)
if grep -E '(/var/spool|^slurm_logs|[[:space:]]slurm_logs)' "${OUT}/submit_plan.txt"; then
  echo "submission plan contains unsafe relative/spool path" >&2
  exit 1
fi
for script in "${PACKAGE_DIR}"/*.sbatch; do
  if grep -E '^#SBATCH --(output|error)=[^/]' "${script}"; then
    echo "relative SBATCH log path in ${script}" >&2
    exit 1
  fi
  if grep -E 'SOURCE_PATH=.*BASH_SOURCE|dirname.*BASH_SOURCE' "${script}"; then
    echo "sbatch runtime path still depends on BASH_SOURCE in ${script}" >&2
    exit 1
  fi
  (
    cd /tmp
    SLURM_SUBMIT_DIR="${REPO_ROOT}" SLURM_ARRAY_TASK_ID=1 PATH_DEBUG=1 bash "${script}" > "${OUT}/$(basename "${script}").paths"
  )
  cat "${OUT}/$(basename "${script}").paths"
  if grep -E '^.*=/var/spool|^.*=/tmp' "${OUT}/$(basename "${script}").paths"; then
    echo "resolved path points to /var/spool or /tmp in ${script}" >&2
    exit 1
  fi
  grep -q "^REPO_ROOT=${REPO_ROOT}$" "${OUT}/$(basename "${script}").paths"
  grep -q "^PACKAGE_DIR=${PACKAGE_DIR}$" "${OUT}/$(basename "${script}").paths"
done
echo "path dry-run passed for ${PACKAGE_DIR}"
