#!/usr/bin/env python3
"""Unit tests + optional BAM integration checks for tools/bismark_x_tags.py.

Usage:
  python tests/test_bismark_x_tags.py
  python tests/test_bismark_x_tags.py --se-bam sample_SE.bam --pe-bam sample_PE.bam
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
import unittest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from tools.bismark_x_tags import (
    strand_id_from_conversions,
    validate_paired_group,
    validate_single_end_record,
    validate_bam,
)


class FakeRec:
    def __init__(self, qname, seq, tags, is_reverse=False, is_read1=False, is_read2=False):
        self.query_name = qname
        self.query_sequence = seq
        self._tags = tags
        self.is_reverse = is_reverse
        self.is_read1 = is_read1
        self.is_read2 = is_read2
        self.is_unmapped = False
        self.is_secondary = False
        self.is_supplementary = False

    def has_tag(self, tag):
        return tag in self._tags

    def get_tag(self, tag):
        return self._tags[tag]


class TestTagRules(unittest.TestCase):
    def test_strand_id_mapping(self):
        self.assertEqual(strand_id_from_conversions("CT", "CT"), "OT")
        self.assertEqual(strand_id_from_conversions("GA", "GA"), "CTOB")
        self.assertEqual(strand_id_from_conversions("GA", "CT"), "CTOT")
        self.assertEqual(strand_id_from_conversions("CT", "GA"), "OB")

    def test_single_end_valid(self):
        rec = FakeRec("r1", "ACGT", {"XM": "....", "XR": "CT", "XG": "CT"}, is_reverse=False)
        self.assertEqual(validate_single_end_record(rec), [])

    def test_single_end_invalid_combo(self):
        rec = FakeRec("r2", "ACGT", {"XM": "....", "XR": "CT", "XG": "GA"}, is_reverse=False)
        self.assertTrue(validate_single_end_record(rec))

    def test_paired_valid(self):
        r1 = FakeRec("p1", "ACGT", {"XM": "....", "XR": "CT", "XG": "CT"}, is_read1=True)
        r2 = FakeRec("p1", "ACGT", {"XM": "....", "XR": "GA", "XG": "CT"}, is_read2=True, is_reverse=True)
        self.assertEqual(validate_paired_group([r1, r2]), [])

    def test_paired_invalid(self):
        r1 = FakeRec("p2", "ACGT", {"XM": "....", "XR": "CT", "XG": "CT"}, is_read1=True)
        r2 = FakeRec("p2", "ACGT", {"XM": "....", "XR": "CT", "XG": "CT"}, is_read2=True, is_reverse=True)
        self.assertTrue(validate_paired_group([r1, r2]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--se-bam")
    parser.add_argument("--pe-bam")
    args, remaining = parser.parse_known_args()

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestTagRules)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    rc = 0 if result.wasSuccessful() else 1

    if args.se_bam:
        checked, nerr, errs = validate_bam(args.se_bam, mode="single")
        print(f"[SE BAM] checked={checked} errors={nerr}")
        for e in errs[:10]:
            print(f"  - {e.qname}: {e.reason}")
        rc = 1 if nerr else rc

    if args.pe_bam:
        checked, nerr, errs = validate_bam(args.pe_bam, mode="paired")
        print(f"[PE BAM] checked={checked} errors={nerr}")
        for e in errs[:10]:
            print(f"  - {e.qname}: {e.reason}")
        rc = 1 if nerr else rc

    raise SystemExit(rc)
