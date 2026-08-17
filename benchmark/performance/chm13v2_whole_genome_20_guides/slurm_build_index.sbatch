#!/usr/bin/env bash
#SBATCH --job-name=chm13_wg_index
#SBATCH --output=/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/crisprapido/benchmark/performance/chm13v2_whole_genome_20_guides/slurm_logs/%x_%j.out
#SBATCH --error=/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/crisprapido/benchmark/performance/chm13v2_whole_genome_20_guides/slurm_logs/%x_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=04:00:00

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}"
REPO_ROOT="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
PACKAGE_DIR="${REPO_ROOT}/benchmark/performance/chm13v2_whole_genome_20_guides"
[[ -d "${PACKAGE_DIR}" ]] || { echo "Missing package: ${PACKAGE_DIR}" >&2; exit 1; }
mkdir -p "${PACKAGE_DIR}/slurm_logs"
if [[ "${PATH_DEBUG:-0}" == "1" ]]; then
  printf 'REPO_ROOT=%s\nPACKAGE_DIR=%s\nLOG_DIR=%s\nRAW_DIR=%s\n' "${REPO_ROOT}" "${PACKAGE_DIR}" "${PACKAGE_DIR}/slurm_logs" "${PACKAGE_DIR}/raw"
  exit 0
fi
export BUILD_INDEX=1
exec "${PACKAGE_DIR}/prepare_reference.sh"

