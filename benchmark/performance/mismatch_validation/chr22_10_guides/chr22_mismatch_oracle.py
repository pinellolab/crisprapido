#!/usr/bin/env python3
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
SUMMARY_FIELDS = [
    "config",
    "max_mismatches",
    "guide_count",
    "baseline_records",
    "baseline_valid_loci",
    "columba_records",
    "columba_valid_loci",
    "columba_recovered_baseline_valid_loci",
    "baseline_missing_loci",
    "columba_only_valid_loci",
    "baseline_invalid_alignment_loci",
    "columba_invalid_alignment_loci",
    "baseline_nonrequested_pam_loci",
    "columba_nonrequested_pam_loci",
    "baseline_command_failures",
    "columba_command_failures",
    "baseline_deterministic",
    "columba_deterministic",
]


def reverse_complement(sequence):
    return sequence.translate(COMP)[::-1].upper()


def read_tsv(path):
    with Path(path).open(encoding="ascii") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):
    records = {}
    name = None
    parts = []
    with Path(path).open(encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper()
    return records


def parse_tags(fields):
    parsed = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            parsed[parts[0]] = parts[2]
    return parsed


def expand_cigar(cigar):
    operations = []
    digits = ""
    for character in cigar:
        if character.isdigit():
            digits += character
            continue
        if character not in "MIDNSHP=X":
            return None
        length = int(digits) if digits else 1
        if length <= 0:
            return None
        operations.extend(character for _ in range(length))
        digits = ""
    if digits:
        return None
    return "".join(operations)


def compact_cigar(operations):
    if not operations:
        return ""
    chunks = []
    current = operations[0]
    length = 1
    for operation in operations[1:]:
        if operation == current:
            length += 1
        else:
            chunks.append("{}{}".format(length, current))
            current = operation
            length = 1
    chunks.append("{}{}".format(length, current))
    return "".join(chunks)


def safe_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def locus_key(guide_id, reference, strand, start, end):
    return (guide_id, reference, strand, start, end)


def validate_paf(path, mode, config, guide_row, references, max_mismatches):
    guide_id = guide_row["guide_id"]
    guide = guide_row["guide_sequence"].upper()
    details = []
    valid_loci = {}
    invalid_alignment_loci = set()
    nonrequested_pam_loci = set()
    record_count = 0

    with Path(path).open(encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            record_count += 1
            fields = line.split("\t")
            reasons = []
            if len(fields) < 12:
                invalid_alignment_loci.add((guide_id, "malformed", line_number))
                details.append(
                    {
                        "config": config,
                        "mode": mode,
                        "guide_id": guide_id,
                        "reference": "unknown",
                        "start": "unknown",
                        "end": "unknown",
                        "strand": "unknown",
                        "mismatch_count": "unknown",
                        "pam": "unknown",
                        "reported_cigar": "missing",
                        "canonical_cigar": "unknown",
                        "cigar_equivalence": "none",
                        "alignment_valid": "no",
                        "requested_pam": "no",
                        "valid_locus": "no",
                        "reason": "malformed_paf",
                    }
                )
                continue

            reference = fields[5]
            strand = fields[4]
            start = safe_int(fields[7])
            end = safe_int(fields[8])
            tags = parse_tags(fields)
            target = None
            pam = None
            mismatch_count = None
            canonical_operations = None

            if reference not in references:
                reasons.append("missing_reference")
            if strand not in ("+", "-"):
                reasons.append("invalid_strand")
            if start is None or end is None or start < 0 or end < start:
                reasons.append("invalid_coordinates")
            elif reference in references and end > len(references[reference]):
                reasons.append("coordinates_out_of_bounds")

            if not reasons:
                raw_target = references[reference][start:end]
                target = raw_target if strand == "+" else reverse_complement(raw_target)
                if len(target) != len(guide):
                    reasons.append("reference_span_not_guide_length")
                else:
                    canonical_operations = "".join(
                        "=" if guide_base == target_base else "X"
                        for guide_base, target_base in zip(guide, target)
                    )
                    mismatch_count = canonical_operations.count("X")
                    if mismatch_count > max_mismatches:
                        reasons.append("mismatch_limit")

                if strand == "+":
                    pam_end = end + 2
                    if pam_end <= len(references[reference]):
                        pam = references[reference][end:pam_end]
                else:
                    if start >= 2:
                        pam = reverse_complement(references[reference][start - 2 : start])

            reported_nm = safe_int(tags.get("nm"))
            reported_ng = safe_int(tags.get("ng"))
            reported_bs = safe_int(tags.get("bs"))
            reported_cigar = tags.get("cg", "")
            expanded_cigar = expand_cigar(reported_cigar)

            if mismatch_count is not None and reported_nm != mismatch_count:
                reasons.append("nm_disagrees_with_reference")
            if reported_ng != 0:
                reasons.append("reported_gap_group")
            if reported_bs != 0:
                reasons.append("reported_gap_size")
            if expanded_cigar is None:
                reasons.append("invalid_cigar")
            elif canonical_operations is not None:
                allowed_cigars = {canonical_operations}
                if strand == "-":
                    allowed_cigars.add(canonical_operations[::-1])
                if expanded_cigar not in allowed_cigars:
                    reasons.append("cigar_disagrees_with_reference")
                if any(operation not in "=X" for operation in expanded_cigar):
                    reasons.append("non_mismatch_cigar_operation")

            reported_matches = safe_int(fields[9])
            if mismatch_count is not None and reported_matches != len(guide) - mismatch_count:
                reasons.append("paf_match_count")

            key = locus_key(guide_id, reference, strand, start, end)
            alignment_valid = not reasons
            requested_pam = pam == "GG"
            valid_locus = alignment_valid and requested_pam
            if alignment_valid and not requested_pam:
                nonrequested_pam_loci.add(key)
            elif not alignment_valid:
                invalid_alignment_loci.add(key)

            cigar_equivalence = "none"
            if expanded_cigar is not None and canonical_operations is not None:
                if expanded_cigar == canonical_operations:
                    cigar_equivalence = "canonical"
                elif strand == "-" and expanded_cigar == canonical_operations[::-1]:
                    cigar_equivalence = "reverse_strand_operation_order"

            detail = {
                "config": config,
                "mode": mode,
                "guide_id": guide_id,
                "reference": reference,
                "start": start,
                "end": end,
                "strand": strand,
                "mismatch_count": mismatch_count if mismatch_count is not None else "unknown",
                "pam": pam if pam is not None else "missing",
                "reported_cigar": reported_cigar if reported_cigar else "missing",
                "canonical_cigar": compact_cigar(canonical_operations)
                if canonical_operations is not None
                else "unknown",
                "cigar_equivalence": cigar_equivalence,
                "alignment_valid": "yes" if alignment_valid else "no",
                "requested_pam": "yes" if requested_pam else "no",
                "valid_locus": "yes" if valid_locus else "no",
                "reason": "accepted"
                if valid_locus
                else ("nonrequested_pam" if alignment_valid else ",".join(sorted(set(reasons)))),
            }
            details.append(detail)
            if valid_locus:
                valid_loci[key] = detail

    return {
        "record_count": record_count,
        "valid_loci": valid_loci,
        "invalid_alignment_loci": invalid_alignment_loci,
        "nonrequested_pam_loci": nonrequested_pam_loci,
        "details": details,
    }


def panel_is_valid(guides, references):
    errors = []
    for row in guides:
        reference = row["chromosome"]
        start = int(row["zero_based_protospacer_start"])
        guide = row["guide_sequence"].upper()
        if row["copy_class"] != "unique" or row["exact_k0_hits"] != "1":
            errors.append("{} is not an exact unique guide".format(row["guide_id"]))
            continue
        if reference not in references:
            errors.append("{} references missing contig {}".format(row["guide_id"], reference))
            continue
        sequence = references[reference]
        if sequence[start : start + len(guide)] != guide:
            errors.append("{} guide sequence does not match source coordinate".format(row["guide_id"]))
        if sequence[start + len(guide) : start + len(guide) + 2] != "GG":
            errors.append("{} source coordinate lacks GG PAM".format(row["guide_id"]))
    return errors


def command_failures(config_dir, guides, mode):
    failures = 0
    for guide in guides:
        for replicate in (1, 2):
            path = config_dir / guide["guide_id"] / "{}_{}.exit_status.txt".format(
                mode, replicate
            )
            if not path.exists() or path.read_text(encoding="ascii").strip() != "0":
                failures += 1
    return failures


def deterministic(config_dir, guides, mode):
    for guide in guides:
        first = config_dir / guide["guide_id"] / "{}_1.paf".format(mode)
        second = config_dir / guide["guide_id"] / "{}_2.paf".format(mode)
        if not first.exists() or not second.exists() or first.read_bytes() != second.read_bytes():
            return False
    return True


def aggregate(args):
    references = read_fasta(args.reference)
    guides = read_tsv(args.guides)
    configs = read_tsv(args.configs)
    errors = panel_is_valid(guides, references)
    if errors:
        raise ValueError("; ".join(errors))

    run_root = Path(args.run_root)
    summary_rows = []
    mismatch_rows = []
    strand_rows = []
    locus_rows = []
    record_rows = []
    failed = False

    for config in configs:
        config_name = config["config"]
        max_mismatches = int(config["max_mismatches"])
        config_dir = run_root / config_name
        mode_results = {}

        for mode in ("baseline", "columba"):
            combined = {
                "record_count": 0,
                "valid_loci": {},
                "invalid_alignment_loci": set(),
                "nonrequested_pam_loci": set(),
                "details": [],
            }
            for guide in guides:
                paf = config_dir / guide["guide_id"] / "{}_1.paf".format(mode)
                result = validate_paf(
                    paf, mode, config_name, guide, references, max_mismatches
                )
                combined["record_count"] += result["record_count"]
                combined["valid_loci"].update(result["valid_loci"])
                combined["invalid_alignment_loci"].update(result["invalid_alignment_loci"])
                combined["nonrequested_pam_loci"].update(result["nonrequested_pam_loci"])
                combined["details"].extend(result["details"])
            mode_results[mode] = combined
            record_rows.extend(combined["details"])

        baseline = mode_results["baseline"]
        columba = mode_results["columba"]
        baseline_keys = set(baseline["valid_loci"])
        columba_keys = set(columba["valid_loci"])
        shared = baseline_keys & columba_keys
        missing = baseline_keys - columba_keys
        additional = columba_keys - baseline_keys
        baseline_failures = command_failures(config_dir, guides, "baseline")
        columba_failures = command_failures(config_dir, guides, "columba")
        baseline_deterministic = deterministic(config_dir, guides, "baseline")
        columba_deterministic = deterministic(config_dir, guides, "columba")

        summary_rows.append(
            {
                "config": config_name,
                "max_mismatches": max_mismatches,
                "guide_count": len(guides),
                "baseline_records": baseline["record_count"],
                "baseline_valid_loci": len(baseline_keys),
                "columba_records": columba["record_count"],
                "columba_valid_loci": len(columba_keys),
                "columba_recovered_baseline_valid_loci": len(shared),
                "baseline_missing_loci": len(missing),
                "columba_only_valid_loci": len(additional),
                "baseline_invalid_alignment_loci": len(baseline["invalid_alignment_loci"]),
                "columba_invalid_alignment_loci": len(columba["invalid_alignment_loci"]),
                "baseline_nonrequested_pam_loci": len(baseline["nonrequested_pam_loci"]),
                "columba_nonrequested_pam_loci": len(columba["nonrequested_pam_loci"]),
                "baseline_command_failures": baseline_failures,
                "columba_command_failures": columba_failures,
                "baseline_deterministic": "yes" if baseline_deterministic else "no",
                "columba_deterministic": "yes" if columba_deterministic else "no",
            }
        )

        for mismatch_count in range(max_mismatches + 1):
            baseline_at_distance = {
                key
                for key in baseline_keys
                if baseline["valid_loci"][key]["mismatch_count"] == mismatch_count
            }
            columba_at_distance = {
                key
                for key in columba_keys
                if columba["valid_loci"][key]["mismatch_count"] == mismatch_count
            }
            mismatch_rows.append(
                {
                    "config": config_name,
                    "mismatch_count": mismatch_count,
                    "baseline_valid_loci": len(baseline_at_distance),
                    "columba_recovered_baseline_valid_loci": len(
                        baseline_at_distance & columba_keys
                    ),
                    "baseline_missing_loci": len(baseline_at_distance - columba_keys),
                    "columba_valid_loci": len(columba_at_distance),
                    "columba_only_valid_loci": len(columba_at_distance - baseline_keys),
                }
            )

        for strand in ("+", "-"):
            baseline_on_strand = {key for key in baseline_keys if key[2] == strand}
            columba_on_strand = {key for key in columba_keys if key[2] == strand}
            strand_rows.append(
                {
                    "config": config_name,
                    "strand": strand,
                    "baseline_valid_loci": len(baseline_on_strand),
                    "columba_recovered_baseline_valid_loci": len(
                        baseline_on_strand & columba_keys
                    ),
                    "baseline_missing_loci": len(baseline_on_strand - columba_keys),
                    "columba_valid_loci": len(columba_on_strand),
                    "columba_only_valid_loci": len(columba_on_strand - baseline_keys),
                }
            )

        for key in sorted(baseline_keys | columba_keys):
            detail = baseline["valid_loci"].get(key) or columba["valid_loci"][key]
            locus_rows.append(
                {
                    "config": config_name,
                    "guide_id": key[0],
                    "reference": key[1],
                    "strand": key[2],
                    "start": key[3],
                    "end": key[4],
                    "mismatch_count": detail["mismatch_count"],
                    "pam": detail["pam"],
                    "baseline": "yes" if key in baseline_keys else "no",
                    "columba": "yes" if key in columba_keys else "no",
                    "classification": "shared"
                    if key in shared
                    else ("baseline_only" if key in missing else "columba_only"),
                }
            )

        if missing or columba["invalid_alignment_loci"] or baseline_failures or columba_failures:
            failed = True
        if not baseline_deterministic or not columba_deterministic:
            failed = True

    write_tsv(args.summary_out, summary_rows, SUMMARY_FIELDS)
    write_tsv(
        args.by_mismatch_out,
        mismatch_rows,
        [
            "config",
            "mismatch_count",
            "baseline_valid_loci",
            "columba_recovered_baseline_valid_loci",
            "baseline_missing_loci",
            "columba_valid_loci",
            "columba_only_valid_loci",
        ],
    )
    write_tsv(
        args.by_strand_out,
        strand_rows,
        [
            "config",
            "strand",
            "baseline_valid_loci",
            "columba_recovered_baseline_valid_loci",
            "baseline_missing_loci",
            "columba_valid_loci",
            "columba_only_valid_loci",
        ],
    )
    write_tsv(
        args.loci_out,
        locus_rows,
        [
            "config",
            "guide_id",
            "reference",
            "strand",
            "start",
            "end",
            "mismatch_count",
            "pam",
            "baseline",
            "columba",
            "classification",
        ],
    )
    write_tsv(
        args.records_out,
        record_rows,
        [
            "config",
            "mode",
            "guide_id",
            "reference",
            "start",
            "end",
            "strand",
            "mismatch_count",
            "pam",
            "reported_cigar",
            "canonical_cigar",
            "cigar_equivalence",
            "alignment_valid",
            "requested_pam",
            "valid_locus",
            "reason",
        ],
    )
    return 1 if failed else 0


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--guides", required=True)
    parser.add_argument("--configs", required=True)
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--by-mismatch-out", required=True)
    parser.add_argument("--by-strand-out", required=True)
    parser.add_argument("--loci-out", required=True)
    parser.add_argument("--records-out", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(aggregate(parse_args()))
