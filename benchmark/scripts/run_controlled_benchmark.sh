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
CRISPRAPIDO_BIN="${REPO_ROOT}/target/release/crisprapido"

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
printf 'k\tmanual_records\tautomatic_records\tpaf_byte_identical\tmanual_exit\tautomatic_exit\tmanual_candidates\tautomatic_candidates\tstderr_diff\n' > "${SUMMARY}"

failed=false

count_manual_candidates() {
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

for k in 0 1 2 3 4; do
  k_dir="${OUT_DIR}/k${k}"
  mkdir -p "${k_dir}"

  manual_paf="${k_dir}/manual.paf"
  manual_stderr="${k_dir}/manual.stderr.txt"
  automatic_paf="${k_dir}/automatic.paf"
  automatic_stderr="${k_dir}/automatic.stderr.txt"
  manual_sam="${CONTROLLED_SAM_DIR}/controlled_k${k}.sam"

  if [[ ! -f "${manual_sam}" ]]; then
    echo "Missing controlled SAM for k=${k}: ${manual_sam}" >&2
    exit 1
  fi

  if run_crisprapido \
    "${manual_paf}" \
    "${manual_stderr}" \
    -r "${CONTROLLED_REFERENCE}" \
    -g "${GUIDE}" \
    -p "${PAM}" \
    -m "${k}" \
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
    -m "${k}" \
    --columba-bin "${COLUMBA_BIN}" \
    --columba-index "${COLUMBA_INDEX}"
  then
    automatic_exit=0
  else
    automatic_exit=$?
  fi

  manual_records="$(wc -l < "${manual_paf}")"
  automatic_records="$(wc -l < "${automatic_paf}")"
  manual_candidates="$(count_manual_candidates "${manual_sam}")"
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

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${k}" \
    "${manual_records}" \
    "${automatic_records}" \
    "${paf_byte_identical}" \
    "${manual_exit}" \
    "${automatic_exit}" \
    "${manual_candidates}" \
    "${automatic_candidates}" \
    "${stderr_diff}" >> "${SUMMARY}"
done

cat "${SUMMARY}"

if [[ "${failed}" == true ]]; then
  echo "Controlled benchmark failed" >&2
  exit 1
fi
