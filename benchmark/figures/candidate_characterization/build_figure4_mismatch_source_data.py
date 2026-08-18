#!/usr/bin/env python3
"""Build compact Figure 4 mismatch source data from finalized oracle summaries."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


FIGURE_DIR = Path(__file__).resolve().parent
REPO_ROOT = FIGURE_DIR.parents[2]
MISMATCH_PACKAGE = (
    REPO_ROOT / "benchmark" / "performance" / "mismatch_validation" / "chr22_10_guides"
)
DEFAULT_OUTPUT = FIGURE_DIR / "figure4_mismatch_source_data.tsv"
SOURCE_FILES = ("summary.tsv", "by_mismatch.tsv", "by_strand.tsv")
FIELDNAMES = [
    "record_type",
    "config",
    "max_mismatches",
    "max_bulges",
    "max_bulge_size",
    "min_match_fraction",
    "candidate_e",
    "guide_count",
    "mismatch_count",
    "strand",
    "baseline_valid_loci",
    "columba_recovered_baseline_valid_loci",
    "baseline_missing_loci",
    "columba_valid_loci",
    "columba_only_valid_loci",
    "columba_invalid_loci",
    "baseline_command_failures",
    "columba_command_failures",
    "baseline_deterministic",
    "columba_deterministic",
    "source_file",
    "source_sha256",
    "config_manifest_sha256",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No rows in {path}")
    return rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def one_by_config(rows: list[dict[str, str]], config: str, source: str) -> dict[str, str]:
    matches = [row for row in rows if row["config"] == config]
    if len(matches) != 1:
        raise ValueError(f"{source}: expected one row for {config}, found {len(matches)}")
    return matches[0]


def make_row(config: dict[str, str], source_file: str, source_hash: str) -> dict[str, str]:
    row = {field: "" for field in FIELDNAMES}
    row.update(
        {
            "config": config["config"],
            "max_mismatches": config["max_mismatches"],
            "max_bulges": config["max_bulges"],
            "max_bulge_size": config["max_bulge_size"],
            "min_match_fraction": config["min_match_fraction"],
            "candidate_e": config["candidate_e"],
            "source_file": source_file,
            "source_sha256": source_hash,
            "config_manifest_sha256": config["config_manifest_sha256"],
        }
    )
    return row


def build_rows(run_root: Path) -> list[dict[str, str]]:
    configs_path = MISMATCH_PACKAGE / "configs.tsv"
    configs = read_tsv(configs_path)
    config_hash = sha256(configs_path)
    for config in configs:
        config["config_manifest_sha256"] = config_hash
        calculated_e = int(config["max_mismatches"]) + int(config["max_bulges"]) * int(
            config["max_bulge_size"]
        )
        if calculated_e != int(config["candidate_e"]):
            raise ValueError(f"Incorrect candidate_e for {config['config']}")

    source_paths = {name: run_root / name for name in SOURCE_FILES}
    for path in source_paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    source_hashes = {name: sha256(path) for name, path in source_paths.items()}
    summary_rows = read_tsv(source_paths["summary.tsv"])
    mismatch_rows = read_tsv(source_paths["by_mismatch.tsv"])
    strand_rows = read_tsv(source_paths["by_strand.tsv"])

    output: list[dict[str, str]] = []
    for config in configs:
        name = config["config"]
        summary = one_by_config(summary_rows, name, "summary.tsv")
        summary_out = make_row(config, "summary.tsv", source_hashes["summary.tsv"])
        summary_out.update(
            {
                "record_type": "summary",
                "guide_count": summary["guide_count"],
                "baseline_valid_loci": summary["baseline_valid_loci"],
                "columba_recovered_baseline_valid_loci": summary[
                    "columba_recovered_baseline_valid_loci"
                ],
                "baseline_missing_loci": summary["baseline_missing_loci"],
                "columba_valid_loci": summary["columba_valid_loci"],
                "columba_only_valid_loci": summary["columba_only_valid_loci"],
                "columba_invalid_loci": summary["columba_invalid_alignment_loci"],
                "baseline_command_failures": summary["baseline_command_failures"],
                "columba_command_failures": summary["columba_command_failures"],
                "baseline_deterministic": summary["baseline_deterministic"],
                "columba_deterministic": summary["columba_deterministic"],
            }
        )
        output.append(summary_out)

        selected_mismatch = [row for row in mismatch_rows if row["config"] == name]
        selected_strand = [row for row in strand_rows if row["config"] == name]
        for source_name, record_type, selected, grouping_field in (
            ("by_mismatch.tsv", "mismatch_count", selected_mismatch, "mismatch_count"),
            ("by_strand.tsv", "strand", selected_strand, "strand"),
        ):
            for source_row in selected:
                row = make_row(config, source_name, source_hashes[source_name])
                row.update(
                    {
                        "record_type": record_type,
                        grouping_field: source_row[grouping_field],
                        "baseline_valid_loci": source_row["baseline_valid_loci"],
                        "columba_recovered_baseline_valid_loci": source_row[
                            "columba_recovered_baseline_valid_loci"
                        ],
                        "baseline_missing_loci": source_row["baseline_missing_loci"],
                        "columba_valid_loci": source_row["columba_valid_loci"],
                        "columba_only_valid_loci": source_row["columba_only_valid_loci"],
                    }
                )
                output.append(row)

        for selected, label in ((selected_mismatch, "mismatch"), (selected_strand, "strand")):
            for field in (
                "baseline_valid_loci",
                "columba_recovered_baseline_valid_loci",
                "baseline_missing_loci",
                "columba_valid_loci",
                "columba_only_valid_loci",
            ):
                grouped_total = sum(int(row[field]) for row in selected)
                if grouped_total != int(summary[field]):
                    raise ValueError(
                        f"{name}: {label} total for {field} is {grouped_total}, "
                        f"expected {summary[field]}"
                    )
    return output


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        type=Path,
        required=True,
        help="Finalized chr22 mismatch-validation output directory",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = build_rows(args.run_root.resolve())
    write_tsv(args.output, rows)
    print(f"wrote {args.output} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
