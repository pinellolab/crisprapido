#!/usr/bin/env python3
"""Independent mismatch-only oracle for baseline and Columba PAF output."""

import argparse
import csv
import re
from pathlib import Path


CIGAR_PATTERN = re.compile(r"(\d+)([=XIDM])")


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_fasta(path):
    records = {}
    name = None
    pieces = []
    with path.open(encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(pieces).upper()
                name = line[1:].split()[0]
                pieces = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence before first header")
                pieces.append(line)
    if name is not None:
        records[name] = "".join(pieces).upper()
    return records


def read_expected(path):
    with path.open(encoding="ascii", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compact_cigar(guide, target):
    operations = ["=" if left == right else "X" for left, right in zip(guide, target)]
    result = []
    current = operations[0]
    length = 0
    for operation in operations:
        if operation == current:
            length += 1
        else:
            result.append("{}{}".format(length, current))
            current = operation
            length = 1
    result.append("{}{}".format(length, current))
    return "".join(result)


def reverse_cigar_operation_order(cigar):
    operations = CIGAR_PATTERN.findall(cigar)
    if "".join("{}{}".format(length, operation) for length, operation in operations) != cigar:
        raise ValueError("Invalid compact CIGAR {}".format(cigar))
    return "".join("{}{}".format(length, operation) for length, operation in reversed(operations))


def parse_tags(fields):
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_paf(path):
    rows = []
    with path.open(encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise ValueError("{}:{} has fewer than 12 PAF fields".format(path, line_number))
            rows.append(
                {
                    "source_line": line_number,
                    "reference_name": fields[5],
                    "start": int(fields[7]),
                    "end": int(fields[8]),
                    "strand": fields[4],
                    "tags": parse_tags(fields[12:]),
                }
            )
    return rows


def extract_target_and_pam(reference, start, end, strand):
    if start < 0 or end > len(reference) or start >= end:
        raise ValueError("Invalid reference interval {}..{}".format(start, end))
    target = reference[start:end]
    if strand == "+":
        pam = reference[end : end + 2]
    elif strand == "-":
        if start < 2:
            pam = ""
        else:
            pam = reverse_complement(reference[start - 2 : start])
        target = reverse_complement(target)
    else:
        raise ValueError("Invalid strand {}".format(strand))
    return target, pam


def validate_fixture(references, expected_rows):
    errors = []
    guide_sequences = {row["guide_sequence"] for row in expected_rows}
    if len(guide_sequences) != 1:
        errors.append("fixture must contain exactly one guide sequence")
        return errors
    guide = next(iter(guide_sequences))
    expected_keys = set()
    for row in expected_rows:
        reference_name = row["reference_name"]
        start = int(row["reference_start_0based"])
        end = int(row["reference_end_0based"])
        strand = row["strand"]
        key = (reference_name, start, end, strand)
        if key in expected_keys:
            errors.append("duplicate expected locus {}".format(key))
        expected_keys.add(key)
        reference = references.get(reference_name)
        if reference is None:
            errors.append("missing reference {}".format(reference_name))
            continue
        target, pam = extract_target_and_pam(reference, start, end, strand)
        mismatches = sum(left != right for left, right in zip(guide, target))
        if len(target) != len(guide):
            errors.append("{} target length is not {}".format(reference_name, len(guide)))
        if target != row["oriented_target_sequence"]:
            errors.append("{} target sequence differs from expected".format(reference_name))
        if pam != row["expected_pam"]:
            errors.append("{} PAM {} != {}".format(reference_name, pam, row["expected_pam"]))
        if mismatches != int(row["expected_mismatches"]):
            errors.append("{} mismatch count differs".format(reference_name))
        if compact_cigar(guide, target) != row["expected_cigar"]:
            errors.append("{} canonical CIGAR differs".format(reference_name))

    for reference_name, reference in references.items():
        for start in range(0, len(reference) - len(guide) + 1):
            end = start + len(guide)
            genomic = reference[start:end]
            for strand, target in (("+", genomic), ("-", reverse_complement(genomic))):
                mismatches = sum(left != right for left, right in zip(guide, target))
                if mismatches <= 3 and (reference_name, start, end, strand) not in expected_keys:
                    errors.append(
                        "unexpected <=3-mismatch locus {}:{}-{}:{}".format(
                            reference_name, start, end, strand
                        )
                    )
    return errors


def validate_paf(mode, path, references, expected_rows, max_mismatches):
    expected = {
        (
            row["reference_name"],
            int(row["reference_start_0based"]),
            int(row["reference_end_0based"]),
            row["strand"],
        ): row
        for row in expected_rows
        if int(row["expected_mismatches"]) <= max_mismatches
    }
    guide = expected_rows[0]["guide_sequence"]
    observed = {}
    details = []
    for row in parse_paf(path):
        key = (row["reference_name"], row["start"], row["end"], row["strand"])
        reasons = []
        reference = references.get(row["reference_name"])
        if reference is None:
            reasons.append("unknown_reference")
            target = ""
            pam = ""
        else:
            target, pam = extract_target_and_pam(
                reference, row["start"], row["end"], row["strand"]
            )
        if len(target) != len(guide):
            reasons.append("reference_span_not_guide_length")
            mismatches = -1
            canonical_cigar = "NA"
        else:
            mismatches = sum(left != right for left, right in zip(guide, target))
            canonical_cigar = compact_cigar(guide, target)
        if mismatches > max_mismatches:
            reasons.append("mismatch_threshold")
        if pam != "GG":
            reasons.append("pam")
        tags = row["tags"]
        if tags.get("ng") != "0" or tags.get("bs") != "0":
            reasons.append("reported_gap")
        if tags.get("nm") != str(mismatches):
            reasons.append("reported_nm")
        reported_cigar = tags.get("cg", "missing")
        equivalent_cigars = {canonical_cigar}
        if row["strand"] == "-" and canonical_cigar != "NA":
            equivalent_cigars.add(reverse_cigar_operation_order(canonical_cigar))
        if reported_cigar not in equivalent_cigars:
            reasons.append("reported_cigar")
        if key not in expected:
            reasons.append("unexpected_locus")
        if key in observed:
            reasons.append("duplicate_locus")
        observed[key] = row
        details.append(
            {
                "mode": mode,
                "reference_name": row["reference_name"],
                "reference_start_0based": row["start"],
                "reference_end_0based": row["end"],
                "strand": row["strand"],
                "oracle_mismatches": mismatches,
                "oracle_cigar": canonical_cigar,
                "oracle_pam": pam,
                "reported_nm": tags.get("nm", "missing"),
                "reported_ng": tags.get("ng", "missing"),
                "reported_bs": tags.get("bs", "missing"),
                "reported_cigar": reported_cigar,
                "cigar_equivalence": "canonical"
                if reported_cigar == canonical_cigar
                else (
                    "reverse_strand_operation_order"
                    if reported_cigar in equivalent_cigars
                    else "none"
                ),
                "valid": "yes" if not reasons else "no",
                "reason": ",".join(reasons) if reasons else "accepted",
            }
        )
    missing = sorted(set(expected) - set(observed))
    invalid = [row for row in details if row["valid"] == "no"]
    return {
        "records": len(details),
        "valid": len(details) - len(invalid),
        "expected": len(expected),
        "recovered": len(expected) - len(missing),
        "missing": len(missing),
        "invalid": len(invalid),
        "details": details,
    }


def write_tsv(path, rows, fieldnames):
    with path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--validate-fixture", action="store_true")
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--columba", type=Path)
    parser.add_argument("--max-mismatches", type=int)
    parser.add_argument("--summary-out", type=Path)
    parser.add_argument("--loci-out", type=Path)
    args = parser.parse_args()

    references = read_fasta(args.reference)
    expected_rows = read_expected(args.expected)
    fixture_errors = validate_fixture(references, expected_rows)
    if fixture_errors:
        for error in fixture_errors:
            print("ERROR\t{}".format(error))
        raise SystemExit(1)
    print("fixture_valid=yes")

    if args.validate_fixture and not args.baseline and not args.columba:
        return
    required = (args.baseline, args.columba, args.max_mismatches, args.summary_out, args.loci_out)
    if any(value is None for value in required):
        parser.error("PAF validation requires baseline, columba, max-mismatches, and output paths")

    baseline = validate_paf(
        "baseline", args.baseline, references, expected_rows, args.max_mismatches
    )
    columba = validate_paf(
        "columba", args.columba, references, expected_rows, args.max_mismatches
    )
    summary = [
        {
            "max_mismatches": args.max_mismatches,
            "max_bulges": 0,
            "max_bulge_size": 0,
            "expected_loci": baseline["expected"],
            "baseline_records": baseline["records"],
            "baseline_valid": baseline["valid"],
            "baseline_recovered": baseline["recovered"],
            "baseline_missing": baseline["missing"],
            "baseline_invalid": baseline["invalid"],
            "columba_records": columba["records"],
            "columba_valid": columba["valid"],
            "columba_recovered": columba["recovered"],
            "columba_missing": columba["missing"],
            "columba_invalid": columba["invalid"],
            "locus_sets_identical": "yes"
            if baseline["missing"] == 0
            and columba["missing"] == 0
            and baseline["invalid"] == 0
            and columba["invalid"] == 0
            else "no",
        }
    ]
    write_tsv(args.summary_out, summary, list(summary[0]))
    details = baseline["details"] + columba["details"]
    write_tsv(args.loci_out, details, list(details[0]))

    failed = any(
        result[key] != 0
        for result in (baseline, columba)
        for key in ("missing", "invalid")
    )
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
