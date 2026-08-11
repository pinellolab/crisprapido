#!/usr/bin/env python3
import argparse
import csv
import hashlib
import re
import shutil
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
DEFAULT_RUN_ID = "chr22_500_guides_20260810T040146_correctness_final"

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    return seq.translate(COMP)[::-1].upper()


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    name = None
    parts = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            else:
                parts.append(line)
    return "".join(parts).upper()


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tags(fields):
    parsed = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            parsed[parts[0]] = parts[2]
    return parsed


def cigar_stats(cigar):
    matches = 0
    mismatches = 0
    query_consumed = 0
    ref_consumed = 0
    gap_groups = 0
    max_gap = 0
    in_gap = False
    current_gap = 0
    ops = set()

    for size, op in ((int(size), op) for size, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)):
        ops.add(op)
        if op in "M=":
            matches += size
            query_consumed += size
            ref_consumed += size
            if in_gap:
                max_gap = max(max_gap, current_gap)
                in_gap = False
                current_gap = 0
        elif op == "X":
            mismatches += size
            query_consumed += size
            ref_consumed += size
            if in_gap:
                max_gap = max(max_gap, current_gap)
                in_gap = False
                current_gap = 0
        elif op == "I":
            query_consumed += size
            if not in_gap:
                gap_groups += 1
                in_gap = True
                current_gap = 0
            current_gap += size
        elif op == "D":
            ref_consumed += size
            if not in_gap:
                gap_groups += 1
                in_gap = True
                current_gap = 0
            current_gap += size
    if in_gap:
        max_gap = max(max_gap, current_gap)

    return {
        "matches": matches,
        "mismatches": mismatches,
        "gap_groups": gap_groups,
        "max_gap": max_gap,
        "query_consumed": query_consumed,
        "ref_consumed": ref_consumed,
        "ops": ops,
    }


def exists_valid_alignment(guide, target):
    guide_len = len(guide)
    target_len = len(target)
    if target_len == guide_len and guide == target:
        return True, "20=", "exact"

    if target_len in (guide_len - 1, guide_len - 2):
        gap = guide_len - target_len
        for pos in range(0, guide_len - gap + 1):
            if guide[:pos] + guide[pos + gap :] == target:
                return (
                    True,
                    (f"{pos}=" if pos else "")
                    + f"{gap}I"
                    + (f"{guide_len - pos - gap}=" if guide_len - pos - gap else ""),
                    "guide_insertion",
                )

    if target_len in (guide_len + 1, guide_len + 2):
        gap = target_len - guide_len
        for pos in range(0, guide_len + 1):
            if target[:pos] == guide[:pos] and target[pos + gap :] == guide[pos:]:
                return (
                    True,
                    (f"{pos}=" if pos else "") + f"{gap}D" + (f"{guide_len - pos}=" if guide_len - pos else ""),
                    "reference_deletion",
                )

    return False, "", "invalid"


class Oracle:
    def __init__(self, reference, pam, max_bulge_size):
        self.reference = reference
        self.pam = pam.upper()
        self.max_bulge_size = max_bulge_size

    def parse_paf(self, path, guide_id, guide):
        rows = []
        if not path.exists():
            return rows
        with path.open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                parsed_tags = tags(fields)
                rows.append(
                    self.evaluate(
                        {
                            "guide_id": guide_id,
                            "guide": guide,
                            "rname": fields[5],
                            "strand": fields[4],
                            "start": int(fields[7]),
                            "end": int(fields[8]),
                            "cg": parsed_tags.get("cg", ""),
                            "line": line.rstrip("\n"),
                        }
                    )
                )
        return rows

    def evaluate(self, row):
        raw = self.reference[row["start"] : row["end"]]
        target = reverse_complement(raw) if row["strand"] == "-" else raw
        if row["strand"] == "+":
            pam = self.reference[row["end"] : row["end"] + len(self.pam)]
        else:
            pam = reverse_complement(self.reference[row["start"] - len(self.pam) : row["start"]]) if row["start"] >= len(self.pam) else ""
        oracle_valid, oracle_cigar, oracle_class = exists_valid_alignment(row["guide"], target)
        return {
            **row,
            **cigar_stats(row["cg"]),
            "target": target,
            "pam": pam,
            "oracle_valid": oracle_valid,
            "oracle_cigar": oracle_cigar,
            "oracle_class": oracle_class,
            "pam_valid": pam == self.pam,
            "intended_valid": oracle_valid and pam == self.pam,
        }

    def ordinary_equivalent(self, baseline, columba):
        if baseline["guide_id"] != columba["guide_id"]:
            return False
        if baseline["rname"] != columba["rname"] or baseline["strand"] != columba["strand"]:
            return False
        if baseline["pam"] != columba["pam"]:
            return False
        if not baseline["intended_valid"] or not columba["intended_valid"]:
            return False
        if max(abs(baseline["start"] - columba["start"]), abs(baseline["end"] - columba["end"])) > self.max_bulge_size:
            return False
        return not (baseline["end"] <= columba["start"] or columba["end"] <= baseline["start"])

    def conservative_indel_equivalent(self, baseline, columba):
        if baseline["guide_id"] != columba["guide_id"]:
            return False
        if baseline["rname"] != columba["rname"] or baseline["strand"] != columba["strand"]:
            return False
        if baseline["pam"] != self.pam or columba["pam"] != self.pam:
            return False
        if not baseline["intended_valid"] or not columba["intended_valid"]:
            return False
        if baseline["start"] != columba["start"]:
            return False
        if baseline["end"] <= columba["start"] or columba["end"] <= baseline["start"]:
            return False
        if baseline["mismatches"] != 0 or columba["mismatches"] != 0:
            return False
        if baseline["gap_groups"] > 1 or columba["gap_groups"] > 1:
            return False
        if baseline["max_gap"] > self.max_bulge_size or columba["max_gap"] > self.max_bulge_size:
            return False
        if not (("D" in baseline["ops"] and "I" in columba["ops"]) or ("I" in baseline["ops"] and "D" in columba["ops"])):
            return False
        if abs(baseline["end"] - columba["end"]) > baseline["max_gap"] + columba["max_gap"]:
            return False
        baseline_target = baseline["target"]
        columba_target = columba["target"]
        return (
            baseline_target.startswith(columba_target)
            or columba_target.startswith(baseline_target)
            or baseline_target.endswith(columba_target)
            or columba_target.endswith(baseline_target)
        )

    def equivalent(self, baseline, columba):
        return self.ordinary_equivalent(baseline, columba) or self.conservative_indel_equivalent(baseline, columba)


def batch_for_guide(guide_id):
    ordinal = int(guide_id.rsplit("_", 1)[1])
    return f"batch_{((ordinal - 1) // 25) + 1:03d}"


def metric_exit(path):
    text = path.read_text()
    return text.split("exit_status\t")[1].split("\n")[0]


def write_manifest(run_root, output_path):
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "size_bytes", "sha256"])
        for path in sorted(run_root.rglob("*")):
            if path.is_file():
                writer.writerow([str(path.relative_to(SCRIPT_DIR)), path.stat().st_size, file_sha256(path)])


def copy_batch_summaries(run_root):
    summary_dir = run_root / "summary"
    for name in ["aggregate_status.tsv", "batch_summary.tsv", "guide_results.tsv"]:
        shutil.copyfile(summary_dir / name, SCRIPT_DIR / name)


def run_oracle(run_root, reference_path, pam, max_bulge_size):
    run_root = run_root.resolve()
    summary_dir = run_root / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    copy_batch_summaries(run_root)

    reference = read_fasta(reference_path)
    oracle = Oracle(reference, pam, max_bulge_size)
    guides = [(row["guide_id"], row["guide_sequence"].upper()) for row in read_tsv(SCRIPT_DIR / "guides.tsv")]

    per_guide = []
    all_columba = []
    for guide_id, guide in guides:
        batch_id = batch_for_guide(guide_id)
        baseline_dir = run_root / "batches" / "correctness" / "baseline" / "pilot" / batch_id / "guides" / guide_id
        columba_dir = run_root / "batches" / "correctness" / "columba" / "pilot" / batch_id / "guides" / guide_id
        baseline_rows = oracle.parse_paf(baseline_dir / "output.paf", guide_id, guide)
        columba_rows = oracle.parse_paf(columba_dir / "output.paf", guide_id, guide)
        all_columba.extend(columba_rows)

        valid_baseline = [row for row in baseline_rows if row["intended_valid"]]
        valid_columba = [row for row in columba_rows if row["intended_valid"]]
        used = set()
        shared = 0
        for baseline_row in valid_baseline:
            hit = None
            for index, columba_row in enumerate(valid_columba):
                if index not in used and oracle.equivalent(baseline_row, columba_row):
                    hit = index
                    break
            if hit is not None:
                used.add(hit)
                shared += 1

        per_guide.append(
            {
                "guide_id": guide_id,
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
        "baseline_raw_records",
        "columba_raw_records",
        "baseline_valid_loci",
        "columba_valid_loci",
        "shared_baseline_loci",
        "baseline_missing_from_columba",
        "columba_only_valid_loci",
        "columba_invalid_records",
        "columba_non_gg_pam_records",
    ]
    aggregate = {key: sum(int(row[key]) for row in per_guide) for key in aggregate_keys}
    exit_ok = all(row["baseline_exit"] == "0" and row["columba_exit"] == "0" for row in per_guide)
    timing_eligible = exit_ok and aggregate["baseline_missing_from_columba"] == 0 and aggregate["columba_invalid_records"] == 0

    per_fields = [
        "guide_id",
        "baseline_exit",
        "columba_exit",
        *aggregate_keys,
        "columba_classes",
    ]
    for path in [SCRIPT_DIR / "per_guide_correctness.tsv", summary_dir / "per_guide_correctness.tsv"]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, per_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(per_guide)

    summary_fields = [*aggregate_keys, "timing_eligible"]
    summary_row = {**aggregate, "timing_eligible": "yes" if timing_eligible else "no"}
    for path in [SCRIPT_DIR / "correctness_summary.tsv", summary_dir / "correctness_summary.tsv"]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, summary_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerow(summary_row)

    guide_results = read_tsv(summary_dir / "guide_results.tsv")
    run_rows = []
    for mode in ["baseline", "columba"]:
        rows = [row for row in guide_results if row["mode"] == mode]
        run_rows.append(
            {
                "mode": mode,
                "guide_count": len(rows),
                "exit_failures": sum(1 for row in rows if row["exit_status"] != "0"),
                "total_paf_records": sum(int(row["paf_records"]) for row in rows),
                "sum_wall_seconds": f"{sum(float(row['wall_seconds']) for row in rows):.6f}",
                "sum_user_seconds": f"{sum(float(row['user_seconds']) for row in rows):.6f}",
                "sum_system_seconds": f"{sum(float(row['system_seconds']) for row in rows):.6f}",
                "max_peak_rss_kib": max(int(row["peak_rss_kib"]) for row in rows),
            }
        )
    run_fields = ["mode", "guide_count", "exit_failures", "total_paf_records", "sum_wall_seconds", "sum_user_seconds", "sum_system_seconds", "max_peak_rss_kib"]
    for path in [SCRIPT_DIR / "run_summary.tsv", summary_dir / "run_summary.tsv"]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, run_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(run_rows)

    oracle_rows = [
        {
            "guide_id": row["guide_id"],
            "strand": row["strand"],
            "start": row["start"],
            "end": row["end"],
            "cg": row["cg"],
            "oracle_cigar": row["oracle_cigar"],
            "oracle_class": row["oracle_class"],
            "pam": row["pam"],
            "pam_valid": row["pam_valid"],
            "oracle_valid": row["oracle_valid"],
            "intended_valid": row["intended_valid"],
            "target": row["target"],
        }
        for row in all_columba
    ]
    with (summary_dir / "oracle_columba_rows.tsv").open("w", newline="") as handle:
        fields = ["guide_id", "strand", "start", "end", "cg", "oracle_cigar", "oracle_class", "pam", "pam_valid", "oracle_valid", "intended_valid", "target"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(oracle_rows)

    write_manifest(run_root, SCRIPT_DIR / "raw_correctness_manifest.tsv")
    return summary_row


def make_row(guide_id, strand, start, end, cigar, target, pam, guide="A" * 20):
    row = {
        "guide_id": guide_id,
        "guide": guide,
        "rname": "22",
        "strand": strand,
        "start": start,
        "end": end,
        "cg": cigar,
        "target": target,
        "pam": pam,
        "oracle_valid": True,
        "pam_valid": pam == "GG",
        "intended_valid": pam == "GG",
    }
    row.update(cigar_stats(cigar))
    return row


def run_self_test():
    oracle = Oracle("A" * 100, "GG", 2)

    known_pairs = [
        ("18=1D2=", "2I18=", "AGCTGATGTAACAATTGTCTA", "TGATGTAACAATTGTCTA"),
        ("17=1D3=", "1=2I17=", "GGAGTGGAATGGAGTGGAATG", "GTGGAATGGAGTGGAATG"),
        ("18=1D2=", "2=2I16=", "CGTCGGCCTCCCAAAGTGCTG", "CGGCCTCCCAAAGTGCTG"),
    ]
    for baseline_cigar, columba_cigar, baseline_target, columba_target in known_pairs:
        baseline = make_row("guide", "-", 100, 121, baseline_cigar, baseline_target, "GG")
        columba = make_row("guide", "-", 100, 118, columba_cigar, columba_target, "GG")
        assert oracle.conservative_indel_equivalent(baseline, columba)

    baseline = make_row("guide", "-", 100, 121, "18=1D2=", "AGCTGATGTAACAATTGTCTA", "GG")
    assert not oracle.conservative_indel_equivalent(baseline, make_row("guide", "-", 100, 118, "2I18=", "TGATGTAACAATTGTCTA", "AG"))
    assert not oracle.conservative_indel_equivalent(baseline, make_row("guide", "+", 100, 118, "2I18=", "TGATGTAACAATTGTCTA", "GG"))
    assert not oracle.conservative_indel_equivalent(baseline, make_row("guide", "-", 100, 116, "2I18=", "TGATGTAACAATTGTCTA", "GG"))
    assert not oracle.conservative_indel_equivalent(baseline, make_row("guide", "-", 105, 123, "2I18=", "TGATGTAACAATTGTCTA", "GG"))
    assert not oracle.conservative_indel_equivalent(baseline, make_row("guide", "-", 100, 118, "1X1I18=", "TGATGTAACAATTGTCTA", "GG"))
    assert not oracle.conservative_indel_equivalent(baseline, make_row("guide", "-", 100, 118, "1I1=1I17=", "TGATGTAACAATTGTCTA", "GG"))
    print("self-test ok")


def parse_args():
    parser = argparse.ArgumentParser(description="Run chr22_500 biological-locus correctness oracle.")
    parser.add_argument("--run-root", default=str(SCRIPT_DIR / "raw" / DEFAULT_RUN_ID))
    parser.add_argument("--reference", default=str(REPO_ROOT / "../data/real_reference/chm13v2_chr22.fa"))
    parser.add_argument("--pam", default="GG")
    parser.add_argument("--max-bulge-size", type=int, default=2)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.self_test:
        run_self_test()
        return
    summary = run_oracle(Path(args.run_root), Path(args.reference), args.pam, args.max_bulge_size)
    print(summary)


if __name__ == "__main__":
    main()
