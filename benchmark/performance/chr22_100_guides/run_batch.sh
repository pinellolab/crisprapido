#!/usr/bin/env bash
set -euo pipefail

if [[ -n "${PACKAGE_DIR:-}" ]]; then
  PACKAGE_DIR="$(cd "${PACKAGE_DIR}" && pwd)"
  REPO_ROOT="$(cd "${REPO_ROOT:?REPO_ROOT is not set}" && pwd)"
else
  SOURCE_PATH="${BASH_SOURCE[0]}"
  if [[ "${SOURCE_PATH}" != /* ]]; then
    SOURCE_PATH="${PWD}/${SOURCE_PATH}"
  fi
  PACKAGE_DIR="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)"
  REPO_ROOT="$(cd "${PACKAGE_DIR}/../../.." && pwd)"
fi
SCRIPT_DIR="${PACKAGE_DIR}"
RAW_DIR="${PACKAGE_DIR}/raw"
case "${PACKAGE_DIR}:${REPO_ROOT}:${RAW_DIR}" in
  *:/var/spool/*|/var/spool/*|*:/tmp/*|/tmp/*)
    echo "Unsafe runtime path resolved under /var/spool or /tmp" >&2
    echo "PACKAGE_DIR=${PACKAGE_DIR}" >&2
    echo "REPO_ROOT=${REPO_ROOT}" >&2
    echo "RAW_DIR=${RAW_DIR}" >&2
    exit 1
    ;;
esac
RUN_ROOT="${RUN_ROOT:-${RAW_DIR}/${RUN_ID:-manual_$(date -u +%Y%m%dT%H%M%SZ)}}"
BATCH_ID="${BATCH_ID:-${1:-}}"
MODE="${MODE:-${2:-baseline}}"
PHASE="${PHASE:-${3:-correctness}}"
ITERATION="${ITERATION:-${4:-pilot}}"
FORCE="${FORCE:-0}"

if [[ -z "${BATCH_ID}" ]]; then
  echo "Usage: BATCH_ID=batch_001 MODE=baseline PHASE=correctness ITERATION=pilot $0" >&2
  exit 2
fi

CRISPRAPIDO_BIN="${CRISPRAPIDO_BIN:-${REPO_ROOT}/target/release/crisprapido}"
COLUMBA_BIN="${COLUMBA_BIN:-${REPO_ROOT}/../columba/build_Vanilla/columba}"
REFERENCE="${REFERENCE:-${REPO_ROOT}/../data/real_reference/chm13v2_chr22.fa}"
COLUMBA_INDEX="${COLUMBA_INDEX:-${REPO_ROOT}/../results/chm13_chr22_index/chm13v2_chr22}"
GUIDE_TSV="${GUIDE_TSV:-${PACKAGE_DIR}/chr22_guides.tsv}"
BATCH_TSV="${BATCH_TSV:-${PACKAGE_DIR}/batches.tsv}"
PAM="${PAM:-GG}"
MAX_MISMATCHES="${MAX_MISMATCHES:-0}"
MAX_BULGES="${MAX_BULGES:-1}"
MAX_BULGE_SIZE="${MAX_BULGE_SIZE:-2}"
MIN_MATCH_FRACTION="${MIN_MATCH_FRACTION:-0.75}"
THREADS="${THREADS:-1}"
SAMPLE_INTERVAL="${SAMPLE_INTERVAL:-0.05}"
PYTHON="${PYTHON:-python3}"

mkdir -p "${RUN_ROOT}/batches/${PHASE}/${MODE}/${ITERATION}" "${RUN_ROOT}/metadata"
BATCH_DIR="${RUN_ROOT}/batches/${PHASE}/${MODE}/${ITERATION}/${BATCH_ID}"
DONE="${BATCH_DIR}/SUCCESS"
if [[ -e "${DONE}" && "${FORCE}" != "1" ]]; then
  echo "Skipping completed successful batch: ${BATCH_DIR}" >&2
  exit 0
fi
if [[ -e "${BATCH_DIR}" && "${FORCE}" != "1" ]]; then
  echo "Refusing to overwrite incomplete batch directory: ${BATCH_DIR}" >&2
  exit 1
fi
rm -rf "${BATCH_DIR}"
mkdir -p "${BATCH_DIR}/guides"

if [[ ! -x "${CRISPRAPIDO_BIN}" ]]; then echo "Missing executable: ${CRISPRAPIDO_BIN}" >&2; exit 1; fi
if [[ "${MODE}" == "columba" && ! -x "${COLUMBA_BIN}" ]]; then echo "Missing executable: ${COLUMBA_BIN}" >&2; exit 1; fi
if [[ ! -r "${REFERENCE}" ]]; then echo "Missing reference: ${REFERENCE}" >&2; exit 1; fi
if [[ ! -r "${GUIDE_TSV}" || ! -r "${BATCH_TSV}" ]]; then echo "Missing guide or batch TSV" >&2; exit 1; fi

awk -F'\t' -v id="${BATCH_ID}" 'NR>1 && $1==id {print $5}' "${BATCH_TSV}" > "${BATCH_DIR}/guide_ids.txt"
if [[ ! -s "${BATCH_DIR}/guide_ids.txt" ]]; then echo "Unknown batch ID: ${BATCH_ID}" >&2; exit 1; fi
tr ',' '\n' < "${BATCH_DIR}/guide_ids.txt" > "${BATCH_DIR}/guide_order.txt"
sha256sum "${BATCH_DIR}/guide_order.txt" > "${BATCH_DIR}/guide_order.sha256"

{
  printf 'date_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'git_commit\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  printf 'git_branch\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse --abbrev-ref HEAD)"
  printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-NA}"
  printf 'slurm_array_task_id\t%s\n' "${SLURM_ARRAY_TASK_ID:-NA}"
  printf 'batch_id\t%s\n' "${BATCH_ID}"
  printf 'mode\t%s\n' "${MODE}"
  printf 'phase\t%s\n' "${PHASE}"
  printf 'iteration\t%s\n' "${ITERATION}"
  printf 'reference\t%s\n' "${REFERENCE}"
  printf 'columba_index\t%s\n' "${COLUMBA_INDEX}"
  printf 'pam\t%s\n' "${PAM}"
  printf 'max_mismatches\t%s\n' "${MAX_MISMATCHES}"
  printf 'max_bulges\t%s\n' "${MAX_BULGES}"
  printf 'max_bulge_size\t%s\n' "${MAX_BULGE_SIZE}"
  printf 'min_match_fraction\t%s\n' "${MIN_MATCH_FRACTION}"
  printf 'threads\t%s\n' "${THREADS}"
} > "${BATCH_DIR}/environment.tsv"

printf 'guide_id\tmode\tphase\titeration\texit_status\twall_seconds\tuser_seconds\tsystem_seconds\tpeak_rss_kib\tpaf_records\tstdout_sha256\tstderr_sha256\n' > "${BATCH_DIR}/guide_results.tsv"

while IFS= read -r guide_id; do
  guide_seq="$(awk -F'\t' -v id="${guide_id}" 'NR>1 && $1==id {print $2}' "${GUIDE_TSV}")"
  if [[ -z "${guide_seq}" ]]; then echo "Missing guide sequence: ${guide_id}" >&2; exit 1; fi
  out_dir="${BATCH_DIR}/guides/${guide_id}"
  mkdir -p "${out_dir}"
  paf="${out_dir}/output.paf"
  stderr="${out_dir}/stderr.txt"
  metrics="${out_dir}/metrics.tsv"
  cmd=("${CRISPRAPIDO_BIN}" -r "${REFERENCE}" -g "${guide_seq}" -p "${PAM}" -m "${MAX_MISMATCHES}" -b "${MAX_BULGES}" -z "${MAX_BULGE_SIZE}" -f "${MIN_MATCH_FRACTION}" -t "${THREADS}")
  if [[ "${MODE}" == "columba" ]]; then
    cmd+=(--columba-bin "${COLUMBA_BIN}" --columba-index "${COLUMBA_INDEX}")
  fi
  printf '%q ' "${cmd[@]}" > "${out_dir}/command.txt"
  set +e
  env -u RUSTFLAGS -u LIBRARY_PATH -u LD_LIBRARY_PATH "${PYTHON}" "${PACKAGE_DIR}/collect_peak_rss.py" --stdout "${paf}" --stderr "${stderr}" --metrics "${metrics}" --sample-interval "${SAMPLE_INTERVAL}" -- "${cmd[@]}"
  status=$?
  set -e
  wall="$(awk -F'\t' '$1=="wall_seconds" {print $2}' "${metrics}" 2>/dev/null || echo NA)"
  user="$(awk -F'\t' '$1=="user_seconds" {print $2}' "${metrics}" 2>/dev/null || echo NA)"
  sys="$(awk -F'\t' '$1=="system_seconds" {print $2}' "${metrics}" 2>/dev/null || echo NA)"
  rss="$(awk -F'\t' '$1=="peak_rss_kib" {print $2}' "${metrics}" 2>/dev/null || echo NA)"
  out_sha="$(awk -F'\t' '$1=="stdout_sha256" {print $2}' "${metrics}" 2>/dev/null || sha256sum "${paf}" | awk '{print $1}')"
  err_sha="$(awk -F'\t' '$1=="stderr_sha256" {print $2}' "${metrics}" 2>/dev/null || sha256sum "${stderr}" | awk '{print $1}')"
  records="$(wc -l < "${paf}")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${guide_id}" "${MODE}" "${PHASE}" "${ITERATION}" "${status}" "${wall}" "${user}" "${sys}" "${rss}" "${records}" "${out_sha}" "${err_sha}" >> "${BATCH_DIR}/guide_results.tsv"
  if [[ "${status}" != "0" ]]; then
    echo "Guide failed: ${guide_id} exit ${status}" >&2
    exit "${status}"
  fi
done < "${BATCH_DIR}/guide_order.txt"

touch "${DONE}"
