#!/usr/bin/env python3
import argparse
import csv
import hashlib
import importlib.util
import shutil
from collections import Counter, defaultdict
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
SHARED_PATH = SCRIPT_DIR.parent / "chr22_500_guides" / "correctness_oracle.py"
SPEC = importlib.util.spec_from_file_location("chr22_oracle", SHARED_PATH)
SHARED = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SHARED)


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_records(path):
    name = None
    parts = []
    with path.open() as handle:
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


def parse_tags(fields):
    result = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            result[parts[0]] = parts[2]
    return result


def parse_paf(path, guide_id, guide):
    rows = []
    if not path.exists():
        return rows
    with path.open() as handle:
        for raw in handle:
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            tags = parse_tags(fields)
            rows.append(
                {
                    "guide_id": guide_id,
                    "guide": guide,
                    "rname": fields[5],
                    "strand": fields[4],
                    "start": int(fields[7]),
                    "end": int(fields[8]),
                    "cg": tags.get("cg", ""),
                    "line": raw.rstrip("\n"),
                }
            )
    return rows


def evaluate(row, reference, pam):
    raw_target = reference[row["start"] : row["end"]]
    target = SHARED.reverse_complement(raw_target) if row["strand"] == "-" else raw_target
    if row["strand"] == "+":
        observed_pam = reference[row["end"] : row["end"] + len(pam)]
    elif row["start"] >= len(pam):
        observed_pam = SHARED.reverse_complement(reference[row["start"] - len(pam) : row["start"]])
    else:
        observed_pam = ""
    oracle_valid, oracle_cigar, oracle_class = SHARED.exists_valid_alignment(row["guide"], target)
    return {
        **row,
        **SHARED.cigar_stats(row["cg"]),
        "target": target,
        "pam": observed_pam,
        "oracle_valid": oracle_valid,
        "oracle_cigar": oracle_cigar,
        "oracle_class": oracle_class,
        "pam_valid": observed_pam == pam,
        "intended_valid": oracle_valid and observed_pam == pam,
    }


def evaluate_rows(reference_path, rows, pam):
    by_reference = defaultdict(list)
    for row in rows:
        by_reference[row["rname"]].append(row)
    found = set()
    for name, sequence in fasta_records(reference_path):
        if name not in by_reference:
            continue
        found.add(name)
        for row in by_reference[name]:
            row.update(evaluate(row, sequence, pam))
    missing = sorted(set(by_reference) - found)
    if missing:
        raise RuntimeError(f"PAF references absent from FASTA: {','.join(missing)}")


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(run_root, output):
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "size_bytes", "sha256"])
        for path in sorted(run_root.rglob("*")):
            if path.is_file():
                writer.writerow([str(path.relative_to(run_root)), path.stat().st_size, file_sha256(path)])


def batch_map():
    result = {}
    for row in read_tsv(SCRIPT_DIR / "batches.tsv"):
        for guide_id in row["guide_ids"].split(","):
            result[guide_id] = row["batch_id"]
    return result


def metric_exit(path):
    with path.open() as handle:
        values = dict(line.rstrip("\n").split("\t", 1) for line in handle if "\t" in line)
    return values.get("exit_status", "missing")


def match_valid_rows(baseline_rows, columba_rows, equivalence):
    used = set()
    shared = 0
    for baseline in baseline_rows:
        match = next(
            (
                index
                for index, columba in enumerate(columba_rows)
                if index not in used and equivalence.equivalent(baseline, columba)
            ),
            None,
        )
        if match is not None:
            used.add(match)
            shared += 1
    return shared, used


def copy_aggregate_outputs(run_root):
    summary = run_root / "summary"
    for name in ["aggregate_status.tsv", "batch_summary.tsv", "guide_results.tsv"]:
        shutil.copyfile(summary / name, SCRIPT_DIR / name)


def run_oracle(run_root, reference_path, pam, max_bulge_size):
    run_root = run_root.resolve()
    copy_aggregate_outputs(run_root)
    guide_rows = read_tsv(SCRIPT_DIR / "guides.tsv")
    guide_to_batch = batch_map()
    per_guide_raw = {}
    all_rows = []
    for guide_row in guide_rows:
        guide_id = guide_row["guide_id"]
        guide = guide_row["guide_sequence"].upper()
        batch_id = guide_to_batch[guide_id]
        mode_rows = {}
        for mode in ["baseline", "columba"]:
            guide_dir = run_root / "batches" / "correctness" / mode / "pilot" / batch_id / "guides" / guide_id
            rows = parse_paf(guide_dir / "output.paf", guide_id, guide)
            all_rows.extend(rows)
            mode_rows[mode] = (guide_dir, rows)
        per_guide_raw[guide_id] = mode_rows
    evaluate_rows(reference_path, all_rows, pam)

    equivalence = SHARED.Oracle("", pam, max_bulge_size)
    per_guide = []
    all_columba = []
    for guide_row in guide_rows:
        guide_id = guide_row["guide_id"]
        baseline_dir, baseline_rows = per_guide_raw[guide_id]["baseline"]
        columba_dir, columba_rows = per_guide_raw[guide_id]["columba"]
        all_columba.extend(columba_rows)
        valid_baseline = [row for row in baseline_rows if row["intended_valid"]]
        valid_columba = [row for row in columba_rows if row["intended_valid"]]
        shared, used = match_valid_rows(valid_baseline, valid_columba, equivalence)
        per_guide.append(
            {
                "guide_id": guide_id,
                "canonical_chromosome": guide_row["canonical_chromosome"],
                "target_copy_class": guide_row["target_copy_class"],
                "baseline_exit": metric_exit(baseline_dir / "metrics.tsv"),
                "columba_exit": metric_exit(columba_dir / "metrics.tsv"),
                "baseline_raw_records": len(baseline_rows),
                "columba_raw_records": len(columba_rows),
                "baseline_valid_loci": len(valid_baseline),
                "columba_valid_loci": len(valid_columba),
                "shared_baseline_loci": shared,
                "baseline_missing_from_columba": len(valid_baseline) - shared,
                "columba_only_valid_loci": len(valid_columba) - len(used),
                "columba_invalid_records": sum(not row["oracle_valid"] for row in columba_rows),
                "columba_non_gg_pam_records": sum(row["oracle_valid"] and not row["pam_valid"] for row in columba_rows),
                "columba_classes": dict(Counter(row["oracle_class"] for row in columba_rows)),
            }
        )

    aggregate_keys = [
        "baseline_raw_records", "columba_raw_records", "baseline_valid_loci",
        "columba_valid_loci", "shared_baseline_loci", "baseline_missing_from_columba",
        "columba_only_valid_loci", "columba_invalid_records", "columba_non_gg_pam_records",
    ]
    aggregate = {key: sum(int(row[key]) for row in per_guide) for key in aggregate_keys}
    exit_ok = all(row["baseline_exit"] == "0" and row["columba_exit"] == "0" for row in per_guide)
    timing_eligible = (
        len(per_guide) == 20
        and exit_ok
        and aggregate["baseline_missing_from_columba"] == 0
        and aggregate["columba_invalid_records"] == 0
    )
    summary = run_root / "summary"
    per_fields = list(per_guide[0])
    for path in [SCRIPT_DIR / "per_guide_correctness.tsv", summary / "per_guide_correctness.tsv"]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, per_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(per_guide)
    summary_row = {**aggregate, "timing_eligible": "yes" if timing_eligible else "no"}
    for path in [SCRIPT_DIR / "correctness_summary.tsv", summary / "correctness_summary.tsv"]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, list(summary_row), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(summary_row)

    aggregate_results = read_tsv(summary / "guide_results.tsv")
    run_rows = []
    for mode in ["baseline", "columba"]:
        rows = [row for row in aggregate_results if row["mode"] == mode]
        run_rows.append(
            {
                "mode": mode,
                "guide_count": len(rows),
                "exit_failures": sum(row["exit_status"] != "0" for row in rows),
                "total_paf_records": sum(int(row["paf_records"]) for row in rows),
                "sum_wall_seconds": f"{sum(float(row['wall_seconds']) for row in rows):.6f}",
                "sum_user_seconds": f"{sum(float(row['user_seconds']) for row in rows):.6f}",
                "sum_system_seconds": f"{sum(float(row['system_seconds']) for row in rows):.6f}",
                "max_peak_rss_kib": max(int(row["peak_rss_kib"]) for row in rows),
            }
        )
    for path in [SCRIPT_DIR / "run_summary.tsv", summary / "run_summary.tsv"]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, list(run_rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(run_rows)
    write_manifest(run_root, SCRIPT_DIR / "raw_correctness_manifest.tsv")
    return summary_row


def run_self_test():
    SHARED.run_self_test()
    row = {"guide_id": "g", "guide": "A" * 20, "rname": "1", "strand": "+", "start": 2, "end": 22, "cg": "20="}
    result = evaluate(row, "CC" + "A" * 20 + "GG", "GG")
    assert result["intended_valid"]
    assert result["target"] == "A" * 20 and result["pam"] == "GG"
    print("multi-contig self-test ok")


def main():
    parser = argparse.ArgumentParser(description="Independent whole-genome correctness oracle.")
    parser.add_argument("--run-root")
    parser.add_argument("--reference", default=REPO_ROOT / "../data/real_reference/chm13v2.fa")
    parser.add_argument("--pam", default="GG")
    parser.add_argument("--max-bulge-size", type=int, default=2)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    if not args.run_root:
        raise SystemExit("--run-root is required")
    print(run_oracle(Path(args.run_root), Path(args.reference), args.pam, args.max_bulge_size))


if __name__ == "__main__":
    main()

