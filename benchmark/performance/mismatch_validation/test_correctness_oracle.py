#!/usr/bin/env python3
"""Unit tests for the standalone mismatch-only oracle."""

import tempfile
import unittest
from pathlib import Path

import correctness_oracle as oracle


PACKAGE_DIR = Path(__file__).resolve().parent


class CorrectnessOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.references = oracle.read_fasta(PACKAGE_DIR / "reference.fa")
        cls.expected = oracle.read_expected(PACKAGE_DIR / "expected_hits.tsv")

    def test_fixture_has_no_unexpected_three_mismatch_locus(self):
        self.assertEqual(oracle.validate_fixture(self.references, self.expected), [])

    def test_reverse_target_and_pam_are_oriented(self):
        row = next(item for item in self.expected if item["reference_name"] == "exact_reverse")
        target, pam = oracle.extract_target_and_pam(
            self.references[row["reference_name"]],
            int(row["reference_start_0based"]),
            int(row["reference_end_0based"]),
            row["strand"],
        )
        self.assertEqual(target, row["guide_sequence"])
        self.assertEqual(pam, "GG")

    def test_one_mismatch_paf_is_valid(self):
        row = next(
            item for item in self.expected if item["reference_name"] == "mismatch1_forward_pos3"
        )
        paf = (
            "Guide\t20\t0\t20\t{strand}\t{reference}\t54\t{start}\t{end}\t19\t20\t255\t"
            "as:i:3\tnm:i:1\tng:i:0\tbs:i:0\tcg:Z:{cigar}\n"
        ).format(
            strand=row["strand"],
            reference=row["reference_name"],
            start=row["reference_start_0based"],
            end=row["reference_end_0based"],
            cigar=row["expected_cigar"],
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "one.paf"
            path.write_text(paf, encoding="ascii")
            result = oracle.validate_paf(
                "test", path, self.references, [row], max_mismatches=1
            )
        self.assertEqual(result["missing"], 0)
        self.assertEqual(result["invalid"], 0)

    def test_reported_gap_is_rejected(self):
        row = next(item for item in self.expected if item["reference_name"] == "exact_forward")
        paf = (
            "Guide\t20\t0\t20\t+\t{reference}\t54\t{start}\t{end}\t20\t20\t255\t"
            "as:i:0\tnm:i:0\tng:i:1\tbs:i:1\tcg:Z:20=\n"
        ).format(
            reference=row["reference_name"],
            start=row["reference_start_0based"],
            end=row["reference_end_0based"],
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gap.paf"
            path.write_text(paf, encoding="ascii")
            result = oracle.validate_paf(
                "test", path, self.references, [row], max_mismatches=1
            )
        self.assertEqual(result["invalid"], 1)
        self.assertIn("reported_gap", result["details"][0]["reason"])

    def test_reverse_strand_operation_order_is_equivalent(self):
        row = next(
            item for item in self.expected if item["reference_name"] == "mismatch1_reverse_pos17"
        )
        reversed_cigar = oracle.reverse_cigar_operation_order(row["expected_cigar"])
        paf = (
            "Guide\t20\t0\t20\t-\t{reference}\t54\t{start}\t{end}\t19\t20\t255\t"
            "as:i:3\tnm:i:1\tng:i:0\tbs:i:0\tcg:Z:{cigar}\n"
        ).format(
            reference=row["reference_name"],
            start=row["reference_start_0based"],
            end=row["reference_end_0based"],
            cigar=reversed_cigar,
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reverse.paf"
            path.write_text(paf, encoding="ascii")
            result = oracle.validate_paf(
                "test", path, self.references, [row], max_mismatches=1
            )
        self.assertEqual(result["invalid"], 0)
        self.assertEqual(
            result["details"][0]["cigar_equivalence"], "reverse_strand_operation_order"
        )


if __name__ == "__main__":
    unittest.main()
