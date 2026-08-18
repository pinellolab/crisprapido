#!/usr/bin/env bash
set -euo pipefail

SOURCE_PATH="${BASH_SOURCE[0]}"
if [[ "${SOURCE_PATH}" != /* ]]; then SOURCE_PATH="${PWD}/${SOURCE_PATH}"; fi
PACKAGE_DIR="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)"
REPO_ROOT="$(cd "${PACKAGE_DIR}/../../.." && pwd)"
OUT="$(mktemp -d /tmp/chm13_wg20_paths.XXXXXX)"
trap 'rm -rf "${OUT}"' EXIT

for script in "${PACKAGE_DIR}"/*.sbatch; do
  if grep -E '^#SBATCH --(output|error)=[^/]' "${script}"; then
    echo "relative Slurm output path: ${script}" >&2
    exit 1
  fi
  if grep -q 'BASH_SOURCE' "${script}"; then
    echo "sbatch runtime path depends on BASH_SOURCE: ${script}" >&2
    exit 1
  fi
  case "$(basename "${script}")" in
    slurm_measured_*.sbatch)
      grep -qx '#SBATCH --array=1-10%1' "${script}" || {
        echo "measured array is not sequential: ${script}" >&2
        exit 1
      }
      ;;
  esac
  (
    cd /tmp
    SLURM_SUBMIT_DIR="${REPO_ROOT}" SLURM_ARRAY_TASK_ID=1 PATH_DEBUG=1 \
      bash "${script}" > "${OUT}/$(basename "${script}").paths"
  )
  grep -q "^REPO_ROOT=${REPO_ROOT}$" "${OUT}/$(basename "${script}").paths"
  grep -q "^PACKAGE_DIR=${PACKAGE_DIR}$" "${OUT}/$(basename "${script}").paths"
  if grep -E '^.*=(/var/spool|/tmp)' "${OUT}/$(basename "${script}").paths"; then
    echo "unsafe runtime path: ${script}" >&2
    exit 1
  fi
done

(
  cd /tmp
  RUN_ID_PREFIX=path_test "${PACKAGE_DIR}/submit_plan.sh" > "${OUT}/submit_plan.txt"
)
if grep -E '(/var/spool|[[:space:]]slurm_logs)' "${OUT}/submit_plan.txt"; then
  echo "unsafe submission plan path" >&2
  exit 1
fi

(
  cd /tmp
  RUN_ID=path_test_controlled "${PACKAGE_DIR}/submit_controlled_timing.sh" \
    > "${OUT}/controlled_timing.txt"
)
grep -q -- '--array=1-10%1' "${OUT}/controlled_timing.txt"
grep -q 'dependency=afterok:${BASELINE_JOB_ID}' "${OUT}/controlled_timing.txt"
grep -q 'dependency=afterok:${COLUMBA_1_JOB_ID}' "${OUT}/controlled_timing.txt"
grep -q 'dependency=afterok:${COLUMBA_2_JOB_ID}' "${OUT}/controlled_timing.txt"
grep -q 'CACHE_STATE_LABEL=cold_start_candidate' "${OUT}/controlled_timing.txt"
if grep -E '(/var/spool|[[:space:]]slurm_logs)' "${OUT}/controlled_timing.txt"; then
  echo "unsafe controlled timing path" >&2
  exit 1
fi
echo "Slurm path test passed"
