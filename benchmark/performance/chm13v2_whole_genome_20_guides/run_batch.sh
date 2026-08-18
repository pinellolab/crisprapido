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
RAW_DIR="${PACKAGE_DIR}/raw"
case "${PACKAGE_DIR}:${REPO_ROOT}:${RAW_DIR}" in
  *:/var/spool/*|/var/spool/*|*:/tmp/*|/tmp/*)
    echo "Unsafe runtime path resolved under /var/spool or /tmp" >&2
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
  echo "BATCH_ID is required" >&2
  exit 2
fi

CRISPRAPIDO_BIN="${CRISPRAPIDO_BIN:-${REPO_ROOT}/target/release/crisprapido}"
COLUMBA_BIN="${COLUMBA_BIN:-${REPO_ROOT}/../columba/build_Vanilla/columba}"
REFERENCE="${CHM13_FASTA:-${REPO_ROOT}/../data/real_reference/chm13v2.fa}"
COLUMBA_INDEX="${COLUMBA_INDEX:-${REPO_ROOT}/../results/chm13_whole_genome_index/chm13v2}"
GUIDE_TSV="${GUIDE_TSV:-${PACKAGE_DIR}/guides.tsv}"
BATCH_TSV="${BATCH_TSV:-${PACKAGE_DIR}/batches.tsv}"
VALIDATION_TSV="${VALIDATION_TSV:-${PACKAGE_DIR}/validation_summary.tsv}"
PAM="${PAM:-GG}"
MAX_MISMATCHES="${MAX_MISMATCHES:-0}"
MAX_BULGES="${MAX_BULGES:-1}"
MAX_BULGE_SIZE="${MAX_BULGE_SIZE:-2}"
MIN_MATCH_FRACTION="${MIN_MATCH_FRACTION:-0.75}"
THREADS="${THREADS:-1}"
SAMPLE_INTERVAL="${SAMPLE_INTERVAL:-0.05}"
PYTHON="${PYTHON:-python3}"

if [[ ! -x "${CRISPRAPIDO_BIN}" ]]; then echo "Missing executable: ${CRISPRAPIDO_BIN}" >&2; exit 2; fi
if [[ "${MODE}" == "columba" && ! -x "${COLUMBA_BIN}" ]]; then echo "Missing executable: ${COLUMBA_BIN}" >&2; exit 2; fi
if [[ ! -r "${REFERENCE}" ]]; then echo "Missing whole-genome reference: ${REFERENCE}" >&2; exit 2; fi
if [[ ! -r "${GUIDE_TSV}" || ! -r "${BATCH_TSV}" ]]; then echo "Panel has not been prepared" >&2; exit 2; fi
if [[ ! -r "${VALIDATION_TSV}" ]]; then echo "Exact validation has not been completed" >&2; exit 2; fi
if [[ "$(awk -F'\t' 'NR>1 && $10=="pass" {n++} END {print n+0}' "${VALIDATION_TSV}")" != "20" ]]; then
  echo "Exact validation does not contain 20 passing guides" >&2
  exit 2
fi
for suffix in .brt .bwt .cct .fsid .headerSN.bin .meta .pos .rev.brt .sa.4 .sa.bv.4 .sna .txt.bin; do
  if [[ ! -r "${COLUMBA_INDEX}${suffix}" ]]; then
    echo "Incomplete Columba index, missing ${COLUMBA_INDEX}${suffix}" >&2
    exit 2
  fi
done

mkdir -p "${RUN_ROOT}/batches/${PHASE}/${MODE}/${ITERATION}" "${RUN_ROOT}/metadata"
BATCH_DIR="${RUN_ROOT}/batches/${PHASE}/${MODE}/${ITERATION}/${BATCH_ID}"
DONE="${BATCH_DIR}/SUCCESS"
if [[ -e "${DONE}" && "${FORCE}" != "1" ]]; then
  echo "Skipping completed successful batch: ${BATCH_DIR}" >&2
  exit 0
fi
if [[ -e "${BATCH_DIR}" && "${FORCE}" != "1" ]]; then
  echo "Refusing to overwrite incomplete batch: ${BATCH_DIR}" >&2
  exit 1
fi
rm -rf "${BATCH_DIR}"
mkdir -p "${BATCH_DIR}/guides"

awk -F'\t' -v id="${BATCH_ID}" 'NR>1 && $1==id {print $5}' "${BATCH_TSV}" > "${BATCH_DIR}/guide_ids.csv"
if [[ ! -s "${BATCH_DIR}/guide_ids.csv" ]]; then echo "Unknown batch: ${BATCH_ID}" >&2; exit 2; fi
tr ',' '\n' < "${BATCH_DIR}/guide_ids.csv" > "${BATCH_DIR}/guide_order.txt"
sha256sum "${BATCH_DIR}/guide_order.txt" > "${BATCH_DIR}/guide_order.sha256"

{
  printf 'date_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'git_commit\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  printf 'git_branch\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse --abbrev-ref HEAD)"
  printf 'git_status_short\t%s\n' "$(git -C "${REPO_ROOT}" status --short | tr '\n' ';')"
  printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-NA}"
  printf 'slurm_array_job_id\t%s\n' "${SLURM_ARRAY_JOB_ID:-NA}"
  printf 'slurm_array_task_id\t%s\n' "${SLURM_ARRAY_TASK_ID:-NA}"
  printf 'slurm_job_dependency\t%s\n' "${SLURM_JOB_DEPENDENCY:-NA}"
  printf 'slurm_partition\t%s\n' "${SLURM_JOB_PARTITION:-NA}"
  printf 'slurm_node_list\t%s\n' "${SLURM_JOB_NODELIST:-NA}"
  printf 'slurm_cpus_per_task\t%s\n' "${SLURM_CPUS_PER_TASK:-NA}"
  printf 'run_id\t%s\n' "${RUN_ID:-NA}"
  printf 'timing_design\t%s\n' "${TIMING_DESIGN:-NA}"
  printf 'cache_state_label\t%s\n' "${CACHE_STATE_LABEL:-NA}"
  printf 'batch_id\t%s\n' "${BATCH_ID}"
  printf 'mode\t%s\n' "${MODE}"
  printf 'phase\t%s\n' "${PHASE}"
  printf 'iteration\t%s\n' "${ITERATION}"
  printf 'reference\t%s\n' "${REFERENCE}"
  printf 'columba_index\t%s\n' "${COLUMBA_INDEX}"
  printf 'crisprapido_sha256\t%s\n' "$(sha256sum "${CRISPRAPIDO_BIN}" | awk '{print $1}')"
  printf 'columba_sha256\t%s\n' "$(sha256sum "${COLUMBA_BIN}" | awk '{print $1}')"
  printf 'guide_tsv_sha256\t%s\n' "$(sha256sum "${GUIDE_TSV}" | awk '{print $1}')"
  printf 'pam\t%s\n' "${PAM}"
  printf 'max_mismatches\t%s\n' "${MAX_MISMATCHES}"
  printf 'max_bulges\t%s\n' "${MAX_BULGES}"
  printf 'max_bulge_size\t%s\n' "${MAX_BULGE_SIZE}"
  printf 'min_match_fraction\t%s\n' "${MIN_MATCH_FRACTION}"
  printf 'threads\t%s\n' "${THREADS}"
  printf 'candidate_e\t%s\n' "$((MAX_MISMATCHES + MAX_BULGES * MAX_BULGE_SIZE))"
} > "${BATCH_DIR}/environment.tsv"

printf 'guide_id\tmode\tphase\titeration\texit_status\twall_seconds\tuser_seconds\tsystem_seconds\tpeak_rss_kib\tpaf_records\tstdout_sha256\tstderr_sha256\n' > "${BATCH_DIR}/guide_results.tsv"
while IFS= read -r guide_id; do
  guide_seq="$(awk -F'\t' -v id="${guide_id}" 'NR>1 && $1==id {print $2}' "${GUIDE_TSV}")"
  if [[ -z "${guide_seq}" ]]; then echo "Missing guide: ${guide_id}" >&2; exit 2; fi
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
  env -u RUSTFLAGS -u LIBRARY_PATH -u LD_LIBRARY_PATH "${PYTHON}" "${PACKAGE_DIR}/collect_peak_rss.py" \
    --stdout "${paf}" --stderr "${stderr}" --metrics "${metrics}" --sample-interval "${SAMPLE_INTERVAL}" -- "${cmd[@]}"
  status=$?
  set -e
  wall="$(awk -F'\t' '$1=="wall_seconds" {print $2}' "${metrics}")"
  user="$(awk -F'\t' '$1=="user_seconds" {print $2}' "${metrics}")"
  system="$(awk -F'\t' '$1=="system_seconds" {print $2}' "${metrics}")"
  rss="$(awk -F'\t' '$1=="peak_rss_kib" {print $2}' "${metrics}")"
  out_sha="$(awk -F'\t' '$1=="stdout_sha256" {print $2}' "${metrics}")"
  err_sha="$(awk -F'\t' '$1=="stderr_sha256" {print $2}' "${metrics}")"
  records="$(wc -l < "${paf}")"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${guide_id}" "${MODE}" "${PHASE}" "${ITERATION}" "${status}" "${wall}" "${user}" "${system}" "${rss}" "${records}" "${out_sha}" "${err_sha}" \
    >> "${BATCH_DIR}/guide_results.tsv"
  if [[ "${status}" != "0" ]]; then exit "${status}"; fi
done < "${BATCH_DIR}/guide_order.txt"
touch "${DONE}"

