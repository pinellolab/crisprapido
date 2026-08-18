#!/usr/bin/env python3
import tempfile
import unittest
from pathlib import Path

import chr22_mismatch_oracle as oracle


class OracleTests(unittest.TestCase):
    def setUp(self):
        self.guide = {
            "guide_id": "g1",
            "guide_sequence": "ACGTACGTACGTACGTACGT",
        }
        self.references = {
            "22": "TT" + "ACGTACGTACGTACGTACGA" + "GG" + "TT",
            "non_gg": "TT" + "ACGTACGTACGTACGTACGA" + "TG" + "TT",
        }

    def validate(self, row, max_mismatches=1):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.paf"
            path.write_text(row + "\n", encoding="ascii")
            return oracle.validate_paf(
                path, "test", "M1", self.guide, self.references, max_mismatches
            )

    def test_forward_mismatch_is_valid(self):
        result = self.validate(
            "Guide\t20\t0\t20\t+\t22\t26\t2\t22\t19\t20\t255\t"
            "nm:i:1\tng:i:0\tbs:i:0\tcg:Z:19=1X"
        )
        self.assertEqual(len(result["valid_loci"]), 1)
        self.assertEqual(len(result["invalid_alignment_loci"]), 0)

    def test_reverse_operation_order_is_narrowly_equivalent(self):
        reference = "CC" + oracle.reverse_complement("TCGTACGTACGTACGTACGT")
        reference += "TT"
        self.references["reverse"] = reference
        row = (
            "Guide\t20\t0\t20\t-\treverse\t24\t2\t22\t19\t20\t255\t"
            "nm:i:1\tng:i:0\tbs:i:0\tcg:Z:19=1X"
        )
        result = self.validate(row)
        self.assertEqual(len(result["valid_loci"]), 1)

    def test_nonrequested_pam_is_not_valid(self):
        result = self.validate(
            "Guide\t20\t0\t20\t+\tnon_gg\t26\t2\t22\t19\t20\t255\t"
            "nm:i:1\tng:i:0\tbs:i:0\tcg:Z:19=1X"
        )
        self.assertEqual(len(result["valid_loci"]), 0)
        self.assertEqual(len(result["nonrequested_pam_loci"]), 1)

    def test_gap_is_invalid_when_bulges_are_disabled(self):
        result = self.validate(
            "Guide\t20\t0\t20\t+\t22\t26\t2\t22\t19\t20\t255\t"
            "nm:i:1\tng:i:1\tbs:i:1\tcg:Z:1I19="
        )
        self.assertEqual(len(result["valid_loci"]), 0)
        self.assertEqual(len(result["invalid_alignment_loci"]), 1)

    def test_mismatch_above_threshold_is_invalid(self):
        self.references["two"] = "TT" + "TCGTACGTACGTACGTACGA" + "GG"
        result = self.validate(
            "Guide\t20\t0\t20\t+\ttwo\t24\t2\t22\t18\t20\t255\t"
            "nm:i:2\tng:i:0\tbs:i:0\tcg:Z:1X18=1X",
            max_mismatches=1,
        )
        self.assertEqual(len(result["valid_loci"]), 0)
        self.assertEqual(len(result["invalid_alignment_loci"]), 1)


if __name__ == "__main__":
    unittest.main()
