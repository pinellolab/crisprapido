#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

OUT_DIR="${1:-${SCRIPT_DIR}/raw/local_$(date -u +%Y%m%dT%H%M%SZ)}"
CRISPRAPIDO_BIN="${CRISPRAPIDO_BIN:-${REPO_ROOT}/target/release/crisprapido}"
COLUMBA_BIN="${COLUMBA_BIN:-${REPO_ROOT}/../columba/build_Vanilla/columba}"
COLUMBA_BUILD="${COLUMBA_BUILD:-${REPO_ROOT}/../columba/build_Vanilla/columba_build}"
REFERENCE="${SCRIPT_DIR}/reference.fa"
EXPECTED="${SCRIPT_DIR}/expected_hits.tsv"
GUIDE="GAGTCCGAGCAGAAGAAGAA"
PAM="GG"
THREADS="1"

if [[ -e "${OUT_DIR}" ]]; then
  echo "Refusing to overwrite existing output directory: ${OUT_DIR}" >&2
  exit 1
fi
for executable in "${CRISPRAPIDO_BIN}" "${COLUMBA_BIN}" "${COLUMBA_BUILD}"; do
  [[ -x "${executable}" ]] || { echo "Missing executable: ${executable}" >&2; exit 1; }
done

python3 "${SCRIPT_DIR}/prepare_fixture.py" --output-dir "${SCRIPT_DIR}" >/dev/null
python3 "${SCRIPT_DIR}/correctness_oracle.py" \
  --reference "${REFERENCE}" \
  --expected "${EXPECTED}" \
  --validate-fixture

mkdir -p "${OUT_DIR}/index"
INDEX_PREFIX="${OUT_DIR}/index/mismatch_reference"
"${COLUMBA_BUILD}" -f "${REFERENCE}" -r "${INDEX_PREFIX}" \
  >"${OUT_DIR}/index_build.stdout.txt" \
  2>"${OUT_DIR}/index_build.stderr.txt"

printf 'config\tmax_mismatches\tmax_bulges\tmax_bulge_size\texpected_loci\tbaseline_records\tbaseline_valid\tbaseline_recovered\tbaseline_missing\tbaseline_invalid\tcolumba_records\tcolumba_valid\tcolumba_recovered\tcolumba_missing\tcolumba_invalid\tlocus_sets_identical\tbaseline_deterministic\tcolumba_deterministic\n' \
  >"${OUT_DIR}/summary.tsv"
failed=false

while IFS=$'\t' read -r config m b z f candidate_e expected_loci; do
  config_dir="${OUT_DIR}/${config}"
  mkdir -p "${config_dir}"

  common=(
    -r "${REFERENCE}"
    -g "${GUIDE}"
    -p "${PAM}"
    -m "${m}"
    -b "${b}"
    -z "${z}"
    -f "${f}"
    -t "${THREADS}"
  )

  for replicate in 1 2; do
    env -u RUSTFLAGS -u LIBRARY_PATH -u LD_LIBRARY_PATH \
      "${CRISPRAPIDO_BIN}" "${common[@]}" \
      >"${config_dir}/baseline_${replicate}.paf" \
      2>"${config_dir}/baseline_${replicate}.stderr.txt"
    env -u RUSTFLAGS -u LIBRARY_PATH -u LD_LIBRARY_PATH \
      "${CRISPRAPIDO_BIN}" "${common[@]}" \
      --columba-bin "${COLUMBA_BIN}" \
      --columba-index "${INDEX_PREFIX}" \
      >"${config_dir}/columba_${replicate}.paf" \
      2>"${config_dir}/columba_${replicate}.stderr.txt"
  done

  cmp "${config_dir}/baseline_1.paf" "${config_dir}/baseline_2.paf"
  cmp "${config_dir}/columba_1.paf" "${config_dir}/columba_2.paf"

  if ! python3 "${SCRIPT_DIR}/correctness_oracle.py" \
    --reference "${REFERENCE}" \
    --expected "${EXPECTED}" \
    --baseline "${config_dir}/baseline_1.paf" \
    --columba "${config_dir}/columba_1.paf" \
    --max-mismatches "${m}" \
    --summary-out "${config_dir}/summary.tsv" \
    --loci-out "${config_dir}/loci.tsv"
  then
    failed=true
  fi

  awk -F '\t' -v OFS='\t' -v config="${config}" \
    'NR == 2 {print config,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,"yes","yes"}' \
    "${config_dir}/summary.tsv" >>"${OUT_DIR}/summary.tsv"
done < <(tail -n +2 "${SCRIPT_DIR}/configs.tsv")

cat "${OUT_DIR}/summary.tsv"

if [[ "${failed}" == true ]]; then
  echo "Mismatch validation failed" >&2
  exit 1
fi
