#!/usr/bin/env python3
"""Unit tests + optional BAM integration checks for tools/bismark_x_tags.py."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
import unittest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from tools.bismark_x_tags import (
    generate_xm_from_alignment,
    qpos_to_rpos_from_cigar,
    strand_id_from_conversions,
    infer_output_minus_strand,
)


class FakeRec:
    def __init__(self, seq, cigartuples, start=0):
        self.query_sequence = seq
        self.cigartuples = cigartuples
        self.reference_start = start




class FakeTaggedRec:
    def __init__(self, tags, is_reverse=False):
        self._tags = tags
        self.is_reverse = is_reverse

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

    def test_cigar_mapping_with_insertion_and_softclip(self):
        # 2M1I1M1S against ref from 10
        rec = FakeRec(seq="ACGTT", cigartuples=[(0, 2), (1, 1), (0, 1), (4, 1)], start=10)
        self.assertEqual(qpos_to_rpos_from_cigar(rec), [10, 11, None, 12, None])

    def test_generate_xm_ct_cpg(self):
        # ref: A C G T ; read: A C G T => C at pos1 in CpG context, methylated => Z
        qseq = "ACGT"
        q2r = [0, 1, 2, 3]
        ref = "ACGT"
        xm = generate_xm_from_alignment(
            query_seq=qseq,
            xr="CT",
            qpos_to_rpos=q2r,
            is_reverse=False,
            ref_base_fetcher=lambda p: ref[p] if 0 <= p < len(ref) else None,
        )
        self.assertEqual(xm, ".Z..")

    def test_generate_xm_ct_unmethylated_chh(self):
        # ref at target C followed by A/A => CHH ; read has T at C pos => h
        qseq = "ATAA"
        q2r = [0, 1, 2, 3]
        ref = "ACAA"
        xm = generate_xm_from_alignment(
            query_seq=qseq,
            xr="CT",
            qpos_to_rpos=q2r,
            is_reverse=False,
            ref_base_fetcher=lambda p: ref[p] if 0 <= p < len(ref) else None,
        )
        self.assertEqual(xm, ".h..")


    def test_infer_output_minus_strand_with_ys(self):
        rec = FakeTaggedRec({"YS": "OB", "XR": "CT", "XG": "GA"}, is_reverse=False)
        self.assertTrue(infer_output_minus_strand(rec))

        rec2 = FakeTaggedRec({"YS": "OT", "XR": "CT", "XG": "CT"}, is_reverse=True)
        self.assertFalse(infer_output_minus_strand(rec2))
    def test_generate_xm_ga_unmethylated_cpg(self):
        # GA mode: G->A at read pos2, upstream in oriented ref is C => z
        qseq = "AAAA"
        q2r = [0, 1, 2, 3]
        ref = "CCGA"
        xm = generate_xm_from_alignment(
            query_seq=qseq,
            xr="GA",
            qpos_to_rpos=q2r,
            is_reverse=False,
            ref_base_fetcher=lambda p: ref[p] if 0 <= p < len(ref) else None,
        )
        self.assertEqual(xm, "..z.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--se-bam")
    parser.add_argument("--pe-bam")
    parser.add_argument("--reference")
    args, _ = parser.parse_known_args()

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestTagRules)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    rc = 0 if result.wasSuccessful() else 1

    if args.reference and (args.se_bam or args.pe_bam):
        from tools.bismark_x_tags import validate_bam_with_reference

        if args.se_bam:
            checked, nerr, errs = validate_bam_with_reference(args.se_bam, args.reference)
            print(f"[SE BAM] checked={checked} errors={nerr}")
            for e in errs[:10]:
                print(f"  - {e.qname}: {e.reason}")
            rc = 1 if nerr else rc

        if args.pe_bam:
            checked, nerr, errs = validate_bam_with_reference(args.pe_bam, args.reference)
            print(f"[PE BAM] checked={checked} errors={nerr}")
            for e in errs[:10]:
                print(f"  - {e.qname}: {e.reason}")
            rc = 1 if nerr else rc

    raise SystemExit(rc)
