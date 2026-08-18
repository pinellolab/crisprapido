#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

OUT_DIR="${1:-${SCRIPT_DIR}/raw/local_$(date -u +%Y%m%dT%H%M%SZ)}"
CRISPRAPIDO_BIN="${CRISPRAPIDO_BIN:-${REPO_ROOT}/target/release/crisprapido}"
COLUMBA_BIN="${COLUMBA_BIN:-${REPO_ROOT}/../columba/build_Vanilla/columba}"
REFERENCE="${CHR22_REFERENCE:-${REPO_ROOT}/../data/real_reference/chm13v2_chr22.fa}"
COLUMBA_INDEX="${COLUMBA_INDEX:-${REPO_ROOT}/../results/chm13_chr22_index/chm13v2_chr22}"
GUIDES="${SCRIPT_DIR}/guides.tsv"
CONFIGS="${SCRIPT_DIR}/configs.tsv"
PAM="GG"
THREADS="1"

if [[ -e "${OUT_DIR}" ]]; then
  echo "Refusing to overwrite existing output directory: ${OUT_DIR}" >&2
  exit 1
fi

for path in "${CRISPRAPIDO_BIN}" "${COLUMBA_BIN}" "${REFERENCE}" "${GUIDES}" "${CONFIGS}"; do
  [[ -e "${path}" ]] || { echo "Missing required path: ${path}" >&2; exit 1; }
done
[[ -x "${CRISPRAPIDO_BIN}" ]] || { echo "Not executable: ${CRISPRAPIDO_BIN}" >&2; exit 1; }
[[ -x "${COLUMBA_BIN}" ]] || { echo "Not executable: ${COLUMBA_BIN}" >&2; exit 1; }

index_found=false
for suffix in meta bwt brt rev.brt sa.4 sa.bv.4 cct fsid sna pos headerSN.bin txt.bin; do
  if [[ -e "${COLUMBA_INDEX}.${suffix}" ]]; then
    index_found=true
    break
  fi
done
[[ "${index_found}" == true ]] || { echo "Missing Columba index: ${COLUMBA_INDEX}" >&2; exit 1; }

mkdir -p "${OUT_DIR}"
failures=0

{
  printf 'field\tvalue\n'
  printf 'date_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'hostname\t%s\n' "$(hostname)"
  printf 'git_commit\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  printf 'git_status_short\t%s\n' "$(git -C "${REPO_ROOT}" status --short | tr '\n' ';')"
  printf 'crisprapido_bin\t%s\n' "${CRISPRAPIDO_BIN}"
  printf 'crisprapido_sha256\t%s\n' "$(sha256sum "${CRISPRAPIDO_BIN}" | awk '{print $1}')"
  printf 'columba_bin\t%s\n' "${COLUMBA_BIN}"
  printf 'columba_sha256\t%s\n' "$(sha256sum "${COLUMBA_BIN}" | awk '{print $1}')"
  printf 'reference\t%s\n' "${REFERENCE}"
  printf 'reference_sha256\t%s\n' "$(sha256sum "${REFERENCE}" | awk '{print $1}')"
  printf 'columba_index\t%s\n' "${COLUMBA_INDEX}"
  printf 'guides_sha256\t%s\n' "$(sha256sum "${GUIDES}" | awk '{print $1}')"
  printf 'configs_sha256\t%s\n' "$(sha256sum "${CONFIGS}" | awk '{print $1}')"
  printf 'pam\t%s\n' "${PAM}"
  printf 'threads\t%s\n' "${THREADS}"
} >"${OUT_DIR}/environment.tsv"

printf 'config\tguide_id\tmode\treplicate\texit_status\tpaf_records\tstdout_sha256\tstderr_sha256\n' \
  >"${OUT_DIR}/command_results.tsv"
printf 'config\tguide_id\tmode\tbyte_identical\n' >"${OUT_DIR}/determinism.tsv"

run_one() {
  local config="$1" guide_id="$2" guide_sequence="$3" mode="$4" replicate="$5"
  local max_mismatches="$6" max_bulges="$7" max_bulge_size="$8" min_match_fraction="$9"
  local guide_dir="${OUT_DIR}/${config}/${guide_id}"
  local paf="${guide_dir}/${mode}_${replicate}.paf"
  local stderr="${guide_dir}/${mode}_${replicate}.stderr.txt"
  local status_file="${guide_dir}/${mode}_${replicate}.exit_status.txt"
  local command_file="${guide_dir}/${mode}_${replicate}.command.txt"
  local command=(
    "${CRISPRAPIDO_BIN}"
    -r "${REFERENCE}" -g "${guide_sequence}" -p "${PAM}"
    -m "${max_mismatches}" -b "${max_bulges}" -z "${max_bulge_size}"
    -f "${min_match_fraction}" -t "${THREADS}"
  )
  if [[ "${mode}" == columba ]]; then
    command+=(--columba-bin "${COLUMBA_BIN}" --columba-index "${COLUMBA_INDEX}")
  fi

  mkdir -p "${guide_dir}"
  printf '%q ' "${command[@]}" >"${command_file}"
  printf '\n' >>"${command_file}"
  set +e
  env -u RUSTFLAGS -u LIBRARY_PATH -u LD_LIBRARY_PATH "${command[@]}" \
    >"${paf}" 2>"${stderr}"
  local status=$?
  set -e
  printf '%s\n' "${status}" >"${status_file}"
  if [[ "${status}" -ne 0 ]]; then
    failures=$((failures + 1))
  fi

  local paf_records=0
  if [[ -s "${paf}" ]]; then
    paf_records="$(wc -l <"${paf}")"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${config}" "${guide_id}" "${mode}" "${replicate}" "${status}" "${paf_records}" \
    "$(sha256sum "${paf}" | awk '{print $1}')" \
    "$(sha256sum "${stderr}" | awk '{print $1}')" \
    >>"${OUT_DIR}/command_results.tsv"
}

while IFS=$'\t' read -r config max_mismatches max_bulges max_bulge_size min_match_fraction candidate_e; do
  calculated_e=$((max_mismatches + max_bulges * max_bulge_size))
  if [[ "${candidate_e}" -ne "${calculated_e}" ]]; then
    echo "Candidate threshold mismatch for ${config}: expected ${calculated_e}, manifest has ${candidate_e}" >&2
    exit 1
  fi

  while IFS=$'\t' read -r guide_id guide_sequence chromosome start strand pam exact_hits copy_class source_panel; do
    for replicate in 1 2; do
      run_one "${config}" "${guide_id}" "${guide_sequence}" baseline "${replicate}" \
        "${max_mismatches}" "${max_bulges}" "${max_bulge_size}" "${min_match_fraction}"
      run_one "${config}" "${guide_id}" "${guide_sequence}" columba "${replicate}" \
        "${max_mismatches}" "${max_bulges}" "${max_bulge_size}" "${min_match_fraction}"
    done

    for mode in baseline columba; do
      first="${OUT_DIR}/${config}/${guide_id}/${mode}_1.paf"
      second="${OUT_DIR}/${config}/${guide_id}/${mode}_2.paf"
      if cmp -s "${first}" "${second}"; then
        identical=yes
      else
        identical=no
        failures=$((failures + 1))
      fi
      printf '%s\t%s\t%s\t%s\n' "${config}" "${guide_id}" "${mode}" "${identical}" \
        >>"${OUT_DIR}/determinism.tsv"
    done
  done < <(tail -n +2 "${GUIDES}")
done < <(tail -n +2 "${CONFIGS}")

set +e
python3 "${SCRIPT_DIR}/chr22_mismatch_oracle.py" \
  --reference "${REFERENCE}" --guides "${GUIDES}" --configs "${CONFIGS}" \
  --run-root "${OUT_DIR}" --summary-out "${OUT_DIR}/summary.tsv" \
  --by-mismatch-out "${OUT_DIR}/by_mismatch.tsv" \
  --by-strand-out "${OUT_DIR}/by_strand.tsv" \
  --loci-out "${OUT_DIR}/loci.tsv" --records-out "${OUT_DIR}/records.tsv"
oracle_status=$?
set -e

cat "${OUT_DIR}/summary.tsv"
cat "${OUT_DIR}/by_mismatch.tsv"
cat "${OUT_DIR}/by_strand.tsv"

if [[ "${failures}" -ne 0 || "${oracle_status}" -ne 0 ]]; then
  echo "chr22 mismatch validation failed: command/determinism failures=${failures}, oracle_status=${oracle_status}" >&2
  exit 1
fi
