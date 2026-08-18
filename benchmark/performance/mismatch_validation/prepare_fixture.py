#!/usr/bin/env python3
"""Generate the deterministic mismatch-only validation fixture."""

import argparse
import csv
import hashlib
from pathlib import Path


GUIDE_ID = "mismatch_guide"
GUIDE = "GAGTCCGAGCAGAAGAAGAA"
PAM = "GG"
FLANK_LENGTH = 16
MAX_TESTED_MISMATCHES = 3

CASES = (
    ("exact_forward", "+", ()),
    ("exact_reverse", "-", ()),
    ("mismatch1_forward_pos3", "+", (3,)),
    ("mismatch1_forward_pos10", "+", (10,)),
    ("mismatch1_reverse_pos17", "-", (17,)),
    ("mismatch2_forward_pos5_15", "+", (5, 15)),
    ("mismatch2_forward_pos1_20", "+", (1, 20)),
    ("mismatch2_reverse_pos4_17", "-", (4, 17)),
    ("mismatch3_forward_pos2_10_19", "+", (2, 10, 19)),
    ("mismatch3_reverse_pos6_12_18", "-", (6, 12, 18)),
    ("mismatch4_negative_forward", "+", (2, 7, 13, 19)),
    ("mismatch4_negative_reverse", "-", (3, 8, 14, 20)),
)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mutate(sequence, positions):
    substitutions = {"A": "C", "C": "G", "G": "T", "T": "A"}
    bases = list(sequence)
    for position in positions:
        bases[position - 1] = substitutions[bases[position - 1]]
    return "".join(bases)


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


def hamming(left, right):
    return sum(a != b for a, b in zip(left, right))


def hash_bases(label, salt):
    digest = hashlib.sha256("{}:{}".format(label, salt).encode("ascii")).digest()
    return "".join("ACGT"[byte & 3] for byte in digest[:FLANK_LENGTH])


def nearby_hits(sequence):
    hits = []
    for start in range(0, len(sequence) - len(GUIDE) + 1):
        genomic = sequence[start : start + len(GUIDE)]
        for strand, target in (("+", genomic), ("-", reverse_complement(genomic))):
            distance = hamming(GUIDE, target)
            if distance <= MAX_TESTED_MISMATCHES:
                hits.append((start, start + len(GUIDE), strand, distance))
    return hits


def build_case(reference_name, strand, positions):
    target = mutate(GUIDE, positions)
    for salt in range(100000):
        left = hash_bases(reference_name + ":left", salt)
        right = hash_bases(reference_name + ":right", salt)
        if strand == "+":
            start = len(left)
            sequence = left + target + PAM + right
        else:
            start = len(left) + len(PAM)
            sequence = left + reverse_complement(PAM) + reverse_complement(target) + right
        end = start + len(GUIDE)
        expected = [(start, end, strand, len(positions))] if len(positions) <= 3 else []
        if nearby_hits(sequence) == expected:
            return sequence, start, end, target, salt
    raise RuntimeError("Could not create an isolated locus for {}".format(reference_name))


def write_fixture(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    reference_path = output_dir / "reference.fa"
    guides_path = output_dir / "guides.fa"
    expected_path = output_dir / "expected_hits.tsv"

    rows = []
    records = []
    for reference_name, strand, positions in CASES:
        sequence, start, end, target, salt = build_case(reference_name, strand, positions)
        mismatch_count = len(positions)
        records.append((reference_name, sequence))
        rows.append(
            {
                "guide_id": GUIDE_ID,
                "guide_sequence": GUIDE,
                "reference_name": reference_name,
                "reference_start_0based": start,
                "reference_end_0based": end,
                "strand": strand,
                "expected_mismatches": mismatch_count,
                "expected_cigar": compact_cigar(GUIDE, target),
                "expected_alignment_class": "exact" if mismatch_count == 0 else "substitution",
                "expected_pam": PAM,
                "mutation_positions_1based": ",".join(str(value) for value in positions) or "none",
                "oriented_target_sequence": target,
                "accepted_m1": "yes" if mismatch_count <= 1 else "no",
                "accepted_m2": "yes" if mismatch_count <= 2 else "no",
                "accepted_m3": "yes" if mismatch_count <= 3 else "no",
                "control_role": "expected" if mismatch_count <= 3 else "negative_outside_m3",
                "flank_salt": salt,
            }
        )

    guides_path.write_text(">{}\n{}\n".format(GUIDE_ID, GUIDE), encoding="ascii")
    with reference_path.open("w", encoding="ascii") as handle:
        for name, sequence in records:
            handle.write(">{}\n{}\n".format(name, sequence))

    fieldnames = list(rows[0])
    with expected_path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    return reference_path, guides_path, expected_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for reference.fa, guides.fa, and expected_hits.tsv",
    )
    args = parser.parse_args()
    paths = write_fixture(args.output_dir)
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
