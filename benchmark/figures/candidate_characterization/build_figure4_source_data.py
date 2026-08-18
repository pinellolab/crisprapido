#!/usr/bin/env python3
"""Build compact Figure 4 per-locus data from finalized correctness outputs."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
from collections import defaultdict
from pathlib import Path


FIGURE_DIR = Path(__file__).resolve().parent
REPO_ROOT = FIGURE_DIR.parents[2]
BENCHMARK_DIR = (
    REPO_ROOT
    / "benchmark"
    / "performance"
    / "chm13v2_whole_genome_20_guides"
)
DEFAULT_RUN_ROOT = (
    BENCHMARK_DIR / "raw" / "chm13_wg20_correctness_20260814T164635Z"
)
DEFAULT_REFERENCE = REPO_ROOT.parent / "data" / "real_reference" / "chm13v2.fa"
DEFAULT_OUTPUT = FIGURE_DIR / "figure4_source_data.tsv"
SHARED_ORACLE_PATH = (
    REPO_ROOT
    / "benchmark"
    / "performance"
    / "chr22_500_guides"
    / "correctness_oracle.py"
)

OUTPUT_FIELDS = [
    "record_id",
    "guide_id",
    "guide_sequence",
    "reference_name",
    "reference_start_0based",
    "reference_end_0based",
    "strand",
    "comparison_class",
    "shared_with_baseline",
    "requested_pam",
    "observed_pam",
    "requested_pam_match",
    "oracle_valid",
    "intended_valid",
    "oriented_target_sequence",
    "canonical_cigar",
    "canonical_edit_distance",
    "canonical_event_class",
    "canonical_mismatches",
    "canonical_gap_groups",
    "canonical_max_gap_size",
    "canonical_guide_insertion_bases",
    "canonical_reference_deletion_bases",
    "reported_cigar",
    "reported_edit_distance",
    "reported_mismatches",
    "reported_gap_groups",
    "reported_max_gap_size",
    "reported_alignment_score",
    "cfd_score",
    "matched_baseline_start_0based",
    "matched_baseline_end_0based",
    "matched_baseline_cigar",
    "source_paf",
    "source_line",
]


def load_shared_oracle():
    spec = importlib.util.spec_from_file_location("figure4_shared_oracle", SHARED_ORACLE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load oracle: {SHARED_ORACLE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SHARED = load_shared_oracle()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No data rows in {path}")
    return rows


def fasta_records(path: Path):
    name = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def cigar_edit_bases(cigar: str) -> int:
    return sum(
        int(size)
        for size, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if op in "IDX"
    )


def cigar_op_bases(cigar: str, requested_op: str) -> int:
    return sum(
        int(size)
        for size, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if op == requested_op
    )


def parse_paf(path: Path, guide_id: str, guide: str, mode: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF row has fewer than 12 columns: {path}:{line_number}")
            tags = parse_tags(fields)
            missing_tags = sorted({"as", "nm", "ng", "bs", "cg", "cf"} - set(tags))
            if missing_tags:
                raise ValueError(
                    f"Missing PAF tags {','.join(missing_tags)}: {path}:{line_number}"
                )
            rows.append(
                {
                    "guide_id": guide_id,
                    "guide": guide,
                    "mode": mode,
                    "rname": fields[5],
                    "strand": fields[4],
                    "start": int(fields[7]),
                    "end": int(fields[8]),
                    "cg": tags["cg"],
                    "reported_as": int(tags["as"]),
                    "reported_nm": int(tags["nm"]),
                    "reported_ng": int(tags["ng"]),
                    "reported_bs": int(tags["bs"]),
                    "cfd": float(tags["cf"]),
                    "source_paf": path.relative_to(REPO_ROOT).as_posix(),
                    "source_line": line_number,
                }
            )
    return rows


def evaluate_row(row: dict[str, object], reference: str, requested_pam: str) -> None:
    start = int(row["start"])
    end = int(row["end"])
    raw_target = reference[start:end]
    target = SHARED.reverse_complement(raw_target) if row["strand"] == "-" else raw_target
    if row["strand"] == "+":
        observed_pam = reference[end : end + len(requested_pam)]
    elif start >= len(requested_pam):
        observed_pam = SHARED.reverse_complement(
            reference[start - len(requested_pam) : start]
        )
    else:
        observed_pam = ""

    oracle_valid, oracle_cigar, oracle_class = SHARED.exists_valid_alignment(
        row["guide"], target
    )
    canonical_stats = SHARED.cigar_stats(oracle_cigar) if oracle_valid else {
        "mismatches": 0,
        "gap_groups": 0,
        "max_gap": 0,
    }
    guide_insertions = cigar_op_bases(oracle_cigar, "I") if oracle_valid else 0
    reference_deletions = cigar_op_bases(oracle_cigar, "D") if oracle_valid else 0
    if oracle_class == "exact":
        event_class = "exact"
    elif oracle_class == "guide_insertion":
        event_class = f"guide_insertion_{guide_insertions}nt"
    elif oracle_class == "reference_deletion":
        event_class = f"reference_deletion_{reference_deletions}nt"
    else:
        event_class = "invalid"

    row.update(
        {
            **SHARED.cigar_stats(row["cg"]),
            "target": target,
            "pam": observed_pam,
            "observed_pam": observed_pam,
            "oracle_valid": oracle_valid,
            "oracle_cigar": oracle_cigar,
            "oracle_class": oracle_class,
            "canonical_edit_distance": cigar_edit_bases(oracle_cigar) if oracle_valid else "",
            "canonical_event_class": event_class,
            "canonical_mismatches": canonical_stats["mismatches"],
            "canonical_gap_groups": canonical_stats["gap_groups"],
            "canonical_max_gap": canonical_stats["max_gap"],
            "canonical_guide_insertions": guide_insertions,
            "canonical_reference_deletions": reference_deletions,
            "pam_valid": observed_pam == requested_pam,
            "intended_valid": oracle_valid and observed_pam == requested_pam,
        }
    )


def evaluate_rows(reference_path: Path, rows: list[dict[str, object]], requested_pam: str) -> None:
    by_reference: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_reference[str(row["rname"])].append(row)
    found: set[str] = set()
    for name, sequence in fasta_records(reference_path):
        if name not in by_reference:
            continue
        found.add(name)
        for row in by_reference[name]:
            evaluate_row(row, sequence, requested_pam)
    missing = sorted(set(by_reference) - found)
    if missing:
        raise RuntimeError(f"PAF references absent from FASTA: {','.join(missing)}")


def batch_map() -> dict[str, str]:
    result: dict[str, str] = {}
    for row in read_tsv(BENCHMARK_DIR / "batches.tsv"):
        for guide_id in row["guide_ids"].split(","):
            result[guide_id] = row["batch_id"]
    return result


def load_rows(run_root: Path) -> tuple[dict[str, list[dict[str, object]]], dict[str, list[dict[str, object]]]]:
    batches = batch_map()
    baseline_by_guide: dict[str, list[dict[str, object]]] = {}
    columba_by_guide: dict[str, list[dict[str, object]]] = {}
    for guide_row in read_tsv(BENCHMARK_DIR / "guides.tsv"):
        guide_id = guide_row["guide_id"]
        guide = guide_row["guide_sequence"].upper()
        batch_id = batches[guide_id]
        mode_rows = {}
        for mode in ("baseline", "columba"):
            path = (
                run_root
                / "batches"
                / "correctness"
                / mode
                / "pilot"
                / batch_id
                / "guides"
                / guide_id
                / "output.paf"
            )
            if not path.is_file():
                raise FileNotFoundError(path)
            mode_rows[mode] = parse_paf(path, guide_id, guide, mode)
        baseline_by_guide[guide_id] = mode_rows["baseline"]
        columba_by_guide[guide_id] = mode_rows["columba"]
    return baseline_by_guide, columba_by_guide


def classify_loci(
    baseline_by_guide: dict[str, list[dict[str, object]]],
    columba_by_guide: dict[str, list[dict[str, object]]],
    requested_pam: str,
    max_bulge_size: int,
) -> None:
    equivalence = SHARED.Oracle("", requested_pam, max_bulge_size)
    for guide_id, columba_rows in columba_by_guide.items():
        baseline_rows = baseline_by_guide[guide_id]
        valid_baseline = [row for row in baseline_rows if row["intended_valid"]]
        valid_columba_indices = [
            index for index, row in enumerate(columba_rows) if row["intended_valid"]
        ]
        matched_columba: dict[int, dict[str, object]] = {}
        for baseline in valid_baseline:
            match = next(
                (
                    index
                    for index in valid_columba_indices
                    if index not in matched_columba
                    and equivalence.equivalent(baseline, columba_rows[index])
                ),
                None,
            )
            if match is not None:
                matched_columba[match] = baseline

        for index, row in enumerate(columba_rows):
            baseline = matched_columba.get(index)
            if not row["oracle_valid"]:
                comparison_class = "oracle_invalid"
                shared = "not_applicable"
            elif not row["pam_valid"]:
                comparison_class = "non_requested_pam"
                shared = "not_applicable"
            elif baseline is not None:
                comparison_class = "shared_baseline"
                shared = "yes"
            else:
                comparison_class = "columba_only"
                shared = "no"
            row["comparison_class"] = comparison_class
            row["shared_with_baseline"] = shared
            row["matched_baseline"] = baseline


def yes_no(value: object) -> str:
    return "yes" if value else "no"


def output_row(row: dict[str, object], ordinal: int, requested_pam: str) -> dict[str, object]:
    baseline = row.get("matched_baseline")
    return {
        "record_id": f"columba_{ordinal:06d}",
        "guide_id": row["guide_id"],
        "guide_sequence": row["guide"],
        "reference_name": row["rname"],
        "reference_start_0based": row["start"],
        "reference_end_0based": row["end"],
        "strand": row["strand"],
        "comparison_class": row["comparison_class"],
        "shared_with_baseline": row["shared_with_baseline"],
        "requested_pam": requested_pam,
        "observed_pam": row["observed_pam"],
        "requested_pam_match": yes_no(row["pam_valid"]),
        "oracle_valid": yes_no(row["oracle_valid"]),
        "intended_valid": yes_no(row["intended_valid"]),
        "oriented_target_sequence": row["target"],
        "canonical_cigar": row["oracle_cigar"],
        "canonical_edit_distance": row["canonical_edit_distance"],
        "canonical_event_class": row["canonical_event_class"],
        "canonical_mismatches": row["canonical_mismatches"],
        "canonical_gap_groups": row["canonical_gap_groups"],
        "canonical_max_gap_size": row["canonical_max_gap"],
        "canonical_guide_insertion_bases": row["canonical_guide_insertions"],
        "canonical_reference_deletion_bases": row["canonical_reference_deletions"],
        "reported_cigar": row["cg"],
        "reported_edit_distance": cigar_edit_bases(str(row["cg"])),
        "reported_mismatches": row["reported_nm"],
        "reported_gap_groups": row["reported_ng"],
        "reported_max_gap_size": row["reported_bs"],
        "reported_alignment_score": row["reported_as"],
        "cfd_score": f"{float(row['cfd']):.4f}",
        "matched_baseline_start_0based": baseline["start"] if baseline else "",
        "matched_baseline_end_0based": baseline["end"] if baseline else "",
        "matched_baseline_cigar": baseline["cg"] if baseline else "",
        "source_paf": row["source_paf"],
        "source_line": row["source_line"],
    }


def validate_counts(
    baseline_by_guide: dict[str, list[dict[str, object]]],
    output_rows: list[dict[str, object]],
) -> None:
    summary = read_tsv(BENCHMARK_DIR / "correctness_summary.tsv")[0]
    counts = {
        "baseline_raw_records": sum(len(rows) for rows in baseline_by_guide.values()),
        "columba_raw_records": len(output_rows),
        "baseline_valid_loci": sum(
            row["intended_valid"]
            for rows in baseline_by_guide.values()
            for row in rows
        ),
        "columba_valid_loci": sum(row["intended_valid"] == "yes" for row in output_rows),
        "shared_baseline_loci": sum(row["comparison_class"] == "shared_baseline" for row in output_rows),
        "columba_only_valid_loci": sum(row["comparison_class"] == "columba_only" for row in output_rows),
        "columba_invalid_records": sum(row["oracle_valid"] == "no" for row in output_rows),
        "columba_non_gg_pam_records": sum(
            row["oracle_valid"] == "yes" and row["requested_pam_match"] == "no"
            for row in output_rows
        ),
    }
    counts["baseline_missing_from_columba"] = (
        counts["baseline_valid_loci"] - counts["shared_baseline_loci"]
    )
    for field, observed in counts.items():
        expected = int(summary[field])
        if observed != expected:
            raise ValueError(f"{field}: extracted {observed}, expected {expected}")
    if any(int(row["canonical_mismatches"]) != 0 for row in output_rows if row["oracle_valid"] == "yes"):
        raise ValueError("An oracle-valid row has a canonical mismatch under m=0")


def build(run_root: Path, reference_path: Path, output: Path) -> list[dict[str, object]]:
    requested_pam = "GG"
    baseline_by_guide, columba_by_guide = load_rows(run_root)
    all_rows = [
        row
        for rows_by_guide in (baseline_by_guide, columba_by_guide)
        for rows in rows_by_guide.values()
        for row in rows
    ]
    evaluate_rows(reference_path, all_rows, requested_pam)
    classify_loci(baseline_by_guide, columba_by_guide, requested_pam, max_bulge_size=2)
    output_rows = [
        output_row(row, ordinal, requested_pam)
        for ordinal, row in enumerate(
            (row for rows in columba_by_guide.values() for row in rows), start=1
        )
    ]
    validate_counts(baseline_by_guide, output_rows)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    return output_rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Figure 4 per-locus data from retained correctness PAFs."
    )
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    rows = build(args.run_root.resolve(), args.reference.resolve(), args.output.resolve())
    print(f"wrote {args.output} ({len(rows)} Columba records)")


if __name__ == "__main__":
    main()
