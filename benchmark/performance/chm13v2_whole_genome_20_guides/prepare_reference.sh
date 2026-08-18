#!/usr/bin/env bash
set -euo pipefail

SOURCE_PATH="${BASH_SOURCE[0]}"
if [[ "${SOURCE_PATH}" != /* ]]; then
  SOURCE_PATH="${PWD}/${SOURCE_PATH}"
fi
PACKAGE_DIR="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)"
REPO_ROOT="$(cd "${PACKAGE_DIR}/../../.." && pwd)"

CHM13_FASTA="${CHM13_FASTA:-${REPO_ROOT}/../data/real_reference/chm13v2.fa}"
COLUMBA_INDEX="${COLUMBA_INDEX:-${REPO_ROOT}/../results/chm13_whole_genome_index/chm13v2}"
COLUMBA_BUILD="${COLUMBA_BUILD:-${REPO_ROOT}/../columba/build_Vanilla/columba_build}"
BUILD_INDEX="${BUILD_INDEX:-0}"

if [[ ! -r "${CHM13_FASTA}" ]]; then
  echo "Missing full CHM13v2 FASTA: ${CHM13_FASTA}" >&2
  echo "Supply CHM13_FASTA or place the complete assembly at the default path." >&2
  echo "The available chr2 and chr22 files are not a whole-genome substitute." >&2
  exit 2
fi

python3 "${PACKAGE_DIR}/preflight.py" --reference "${CHM13_FASTA}" --index-prefix "${COLUMBA_INDEX}"

if [[ "${BUILD_INDEX}" != "1" ]]; then
  echo "Index audit complete. Set BUILD_INDEX=1 to construct a missing index." >&2
  exit 0
fi
if [[ ! -x "${COLUMBA_BUILD}" ]]; then
  echo "Missing Columba index builder: ${COLUMBA_BUILD}" >&2
  exit 2
fi

INDEX_DIR="$(dirname "${COLUMBA_INDEX}")"
mkdir -p "${INDEX_DIR}"
if python3 "${PACKAGE_DIR}/preflight.py" --reference "${CHM13_FASTA}" --index-prefix "${COLUMBA_INDEX}" --require-ready >/dev/null 2>&1; then
  echo "Complete index already exists: ${COLUMBA_INDEX}" >&2
  exit 0
fi

LOG_PREFIX="${COLUMBA_INDEX}_index_build"
start_epoch="$(date +%s)"
set +e
"${COLUMBA_BUILD}" -f "${CHM13_FASTA}" -r "${COLUMBA_INDEX}" \
  > "${LOG_PREFIX}.stdout.txt" 2> "${LOG_PREFIX}.stderr.txt"
status=$?
set -e
end_epoch="$(date +%s)"
{
  printf 'command\t%q -f %q -r %q\n' "${COLUMBA_BUILD}" "${CHM13_FASTA}" "${COLUMBA_INDEX}"
  printf 'exit_status\t%s\n' "${status}"
  printf 'elapsed_seconds\t%s\n' "$((end_epoch - start_epoch))"
} > "${LOG_PREFIX}.status.tsv"
if [[ "${status}" != "0" ]]; then
  echo "Columba index construction failed with exit ${status}" >&2
  exit "${status}"
fi
python3 "${PACKAGE_DIR}/preflight.py" --reference "${CHM13_FASTA}" --index-prefix "${COLUMBA_INDEX}" --require-ready

