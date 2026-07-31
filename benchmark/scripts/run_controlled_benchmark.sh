#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${BENCHMARK_DIR}/.." && pwd)"

OUT_DIR="${1:-${BENCHMARK_DIR}/controlled}"

COLUMBA_BIN="${COLUMBA_BIN:-${REPO_ROOT}/../columba/build_Vanilla/columba}"
COLUMBA_INDEX="${COLUMBA_INDEX:-${REPO_ROOT}/../results/controlled_index/controlled_reference}"
CONTROLLED_REFERENCE="${CONTROLLED_REFERENCE:-${REPO_ROOT}/../data/controlled_reference.fa}"
CONTROLLED_SAM_DIR="${CONTROLLED_SAM_DIR:-${REPO_ROOT}/../results}"

GUIDE="GAGTCCGAGCAGAAGAAGAA"
PAM="GG"
MIN_MATCH_FRACTION="0.75"
THREADS="1"
CRISPRAPIDO_BIN="${REPO_ROOT}/target/release/crisprapido"

CONFIGS=(
  "A	0	0	0"
  "B	1	0	0"
  "C	0	1	1"
  "D	0	1	2"
  "E	1	1	2"
  "F	2	1	2"
)

if [[ ! -x "${CRISPRAPIDO_BIN}" ]]; then
  echo "Missing executable: ${CRISPRAPIDO_BIN}" >&2
  echo "Build it first with: cargo +stable build --release" >&2
  exit 1
fi

if [[ ! -x "${COLUMBA_BIN}" ]]; then
  echo "Missing or non-executable Columba binary: ${COLUMBA_BIN}" >&2
  exit 1
fi

if [[ ! -f "${CONTROLLED_REFERENCE}" ]]; then
  echo "Missing controlled reference: ${CONTROLLED_REFERENCE}" >&2
  exit 1
fi

has_index_file=false
for suffix in meta bwt brt rev.brt sa.4 sa.bv.4 cct fsid sna pos headerSN.bin txt.bin; do
  if [[ -e "${COLUMBA_INDEX}.${suffix}" ]]; then
    has_index_file=true
    break
  fi
done
if [[ "${has_index_file}" != true ]]; then
  echo "Missing Columba index for prefix: ${COLUMBA_INDEX}" >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"
SUMMARY="${OUT_DIR}/summary.tsv"
printf 'config	m	b	z	f	candidate_e	manual_candidates	automatic_candidates	manual_records	automatic_records	manual_exit	automatic_exit	paf_byte_identical	stderr_diff
' > "${SUMMARY}"

failed=false

count_mapped_candidates() {
  local sam="$1"
  perl -F'\t' -lane 'next if /^@/ || /^\s*$/; next if @F < 11; $flag=$F[1]+0; $c++ if (($flag & 4) == 0); END { print $c+0 }' "${sam}"
}

run_crisprapido() {
  local stdout_path="$1"
  local stderr_path="$2"
  shift 2

  set +e
  env \
    -u RUSTFLAGS \
    -u LIBRARY_PATH \
    -u LD_LIBRARY_PATH \
    "${CRISPRAPIDO_BIN}" "$@" > "${stdout_path}" 2> "${stderr_path}"
  local status=$?
  set -e
  return "${status}"
}

run_columba_manual_generation() {
  local query_fasta="$1"
  local candidate_e="$2"
  local sam_path="$3"
  local stdout_path="$4"
  local stderr_path="$5"

  set +e
  env \
    -u RUSTFLAGS \
    -u LIBRARY_PATH \
    -u LD_LIBRARY_PATH \
    "${COLUMBA_BIN}" \
      -r "${COLUMBA_INDEX}" \
      -f "${query_fasta}" \
      -a all \
      -m edit \
      -e "${candidate_e}" \
      -t "${THREADS}" \
      -o "${sam_path}" > "${stdout_path}" 2> "${stderr_path}"
  local status=$?
  set -e
  return "${status}"
}

for config in "${CONFIGS[@]}"; do
  IFS=$'\t' read -r label max_mismatches max_bulges max_bulge_size <<< "${config}"
  candidate_e=$((max_mismatches + max_bulges * max_bulge_size))
  config_id="${label}_m${max_mismatches}_b${max_bulges}_z${max_bulge_size}_e${candidate_e}"
  config_dir="${OUT_DIR}/${config_id}"
  mkdir -p "${config_dir}"

  query_fasta="${config_dir}/guide.fa"
  manual_sam="${config_dir}/manual.sam"
  columba_stdout="${config_dir}/manual_columba.stdout.txt"
  columba_stderr="${config_dir}/manual_columba.stderr.txt"
  manual_paf="${config_dir}/manual.paf"
  manual_stderr="${config_dir}/manual.stderr.txt"
  automatic_paf="${config_dir}/automatic.paf"
  automatic_stderr="${config_dir}/automatic.stderr.txt"

  printf '>guide\n%s\n' "${GUIDE}" > "${query_fasta}"

  if run_columba_manual_generation \
    "${query_fasta}" \
    "${candidate_e}" \
    "${manual_sam}" \
    "${columba_stdout}" \
    "${columba_stderr}"
  then
    columba_exit=0
  else
    columba_exit=$?
    failed=true
  fi

  if [[ "${columba_exit}" -eq 0 ]]; then
    manual_candidates="$(count_mapped_candidates "${manual_sam}")"
  else
    manual_candidates="not_observed"
  fi

  if run_crisprapido \
    "${manual_paf}" \
    "${manual_stderr}" \
    -r "${CONTROLLED_REFERENCE}" \
    -g "${GUIDE}" \
    -p "${PAM}" \
    -m "${max_mismatches}" \
    -b "${max_bulges}" \
    -z "${max_bulge_size}" \
    -f "${MIN_MATCH_FRACTION}" \
    -t "${THREADS}" \
    --columba-sam "${manual_sam}"
  then
    manual_exit=0
  else
    manual_exit=$?
  fi

  if run_crisprapido \
    "${automatic_paf}" \
    "${automatic_stderr}" \
    -r "${CONTROLLED_REFERENCE}" \
    -g "${GUIDE}" \
    -p "${PAM}" \
    -m "${max_mismatches}" \
    -b "${max_bulges}" \
    -z "${max_bulge_size}" \
    -f "${MIN_MATCH_FRACTION}" \
    -t "${THREADS}" \
    --columba-bin "${COLUMBA_BIN}" \
    --columba-index "${COLUMBA_INDEX}"
  then
    automatic_exit=0
  else
    automatic_exit=$?
  fi

  manual_records="$(wc -l < "${manual_paf}")"
  automatic_records="$(wc -l < "${automatic_paf}")"
  automatic_candidates="not_observed"

  if cmp -s "${manual_paf}" "${automatic_paf}"; then
    paf_byte_identical="yes"
  else
    paf_byte_identical="no"
    failed=true
  fi

  if cmp -s "${manual_stderr}" "${automatic_stderr}"; then
    stderr_diff="none"
  else
    stderr_diff="qname_or_other"
  fi

  if [[ "${manual_exit}" -ne 0 || "${automatic_exit}" -ne 0 ]]; then
    failed=true
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${label}" \
    "${max_mismatches}" \
    "${max_bulges}" \
    "${max_bulge_size}" \
    "${MIN_MATCH_FRACTION}" \
    "${candidate_e}" \
    "${manual_candidates}" \
    "${automatic_candidates}" \
    "${manual_records}" \
    "${automatic_records}" \
    "${manual_exit}" \
    "${automatic_exit}" \
    "${paf_byte_identical}" \
    "${stderr_diff}" >> "${SUMMARY}"
done

cat "${SUMMARY}"

if [[ "${failed}" == true ]]; then
  echo "Controlled benchmark failed" >&2
  exit 1
fi
