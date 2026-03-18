#!/usr/bin/env python3
"""Validate Bismark XM/XR/XG tags for single-end and paired-end BAM records.

This script implements tag conventions used in Bismark's Perl source:
- single-end: `single_end_SAM_output` and `extract_corresponding_genomic_sequence_single_end`
- paired-end: `paired_end_SAM_output` and `extract_corresponding_genomic_sequence_paired_end`
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable



def _require_pysam():
    try:
        import pysam  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("pysam is required for BAM I/O. Install with: pip install pysam") from exc
    return pysam

VALID_XM_CHARS = set(".zZxXhHuU")
VALID_CONVERSIONS = {"CT", "GA"}


@dataclass
class ValidationError:
    qname: str
    reason: str


def strand_id_from_conversions(xr: str, xg: str) -> str:
    mapping = {
        ("CT", "CT"): "OT",
        ("GA", "GA"): "CTOB",
        ("GA", "CT"): "CTOT",
        ("CT", "GA"): "OB",
    }
    return mapping.get((xr, xg), "INVALID")


def validate_common(record) -> list[str]:
    errs: list[str] = []
    if not record.has_tag("XM"):
        errs.append("missing XM tag")
    else:
        xm = record.get_tag("XM")
        if set(xm) - VALID_XM_CHARS:
            errs.append(f"XM contains invalid chars: {''.join(sorted(set(xm) - VALID_XM_CHARS))}")
        seq = record.query_sequence or ""
        if len(xm) != len(seq):
            errs.append(f"XM length ({len(xm)}) != read length ({len(seq)})")

    if not record.has_tag("XR"):
        errs.append("missing XR tag")
    elif record.get_tag("XR") not in VALID_CONVERSIONS:
        errs.append(f"XR must be CT/GA, got {record.get_tag('XR')}")

    if not record.has_tag("XG"):
        errs.append("missing XG tag")
    elif record.get_tag("XG") not in VALID_CONVERSIONS:
        errs.append(f"XG must be CT/GA, got {record.get_tag('XG')}")
    return errs


def validate_single_end_record(record) -> list[str]:
    errs = validate_common(record)
    if errs:
        return errs

    xr = record.get_tag("XR")
    xg = record.get_tag("XG")
    is_reverse = record.is_reverse

    if not is_reverse and (xr, xg) not in {("CT", "CT"), ("GA", "GA")}:
        errs.append("single-end forward record must be (XR,XG) in {(CT,CT),(GA,GA)}")

    if is_reverse and (xr, xg) not in {("CT", "GA"), ("GA", "CT")}:
        errs.append("single-end reverse record must be (XR,XG) in {(CT,GA),(GA,CT)}")

    return errs


def validate_paired_group(records: Iterable) -> list[ValidationError]:
    grouped = {1: [], 2: []}
    errors: list[ValidationError] = []

    for rec in records:
        mate = 1 if rec.is_read1 else 2 if rec.is_read2 else 0
        if mate in (1, 2):
            grouped[mate].append(rec)

    if not grouped[1] or not grouped[2]:
        qname = next(iter(records)).query_name
        return [ValidationError(qname, "paired group missing read1/read2")]

    r1 = grouped[1][0]
    r2 = grouped[2][0]

    for rec in (r1, r2):
        for msg in validate_common(rec):
            errors.append(ValidationError(rec.query_name, msg))

    if errors:
        return errors

    xr1, xg1 = r1.get_tag("XR"), r1.get_tag("XG")
    xr2, xg2 = r2.get_tag("XR"), r2.get_tag("XG")

    if xg1 != xg2:
        errors.append(ValidationError(r1.query_name, f"read1/read2 XG mismatch: {xg1} vs {xg2}"))
        return errors

    if (xr1, xr2, xg1) not in {
        ("CT", "GA", "CT"),
        ("GA", "CT", "CT"),
        ("GA", "CT", "GA"),
        ("CT", "GA", "GA"),
    }:
        errors.append(ValidationError(r1.query_name, f"invalid pair XR/XG combination: R1={xr1}, R2={xr2}, XG={xg1}"))

    return errors


def validate_bam(path: str, mode: str = "auto", max_errors: int = 50) -> tuple[int, int, list[ValidationError]]:
    checked = 0
    errors: list[ValidationError] = []

    pysam = _require_pysam()
    with pysam.AlignmentFile(path, "rb") as bam:
        if mode == "single":
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                checked += 1
                for msg in validate_single_end_record(rec):
                    errors.append(ValidationError(rec.query_name, msg))
                    if len(errors) >= max_errors:
                        return checked, len(errors), errors
        else:
            by_name: dict[str, list] = defaultdict(list)
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                if mode == "auto" and not rec.is_paired:
                    checked += 1
                    for msg in validate_single_end_record(rec):
                        errors.append(ValidationError(rec.query_name, msg))
                        if len(errors) >= max_errors:
                            return checked, len(errors), errors
                    continue

                by_name[rec.query_name].append(rec)

            for qname, recs in by_name.items():
                checked += len(recs)
                for err in validate_paired_group(recs):
                    errors.append(err)
                    if len(errors) >= max_errors:
                        return checked, len(errors), errors

    return checked, len(errors), errors


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Bismark XM/XR/XG tag definitions against BAM records")
    parser.add_argument("bam", help="Input BAM (Bismark output)")
    parser.add_argument("--mode", choices=["auto", "single", "paired"], default="auto")
    parser.add_argument("--max-errors", type=int, default=50)
    args = parser.parse_args()

    checked, n_errors, errors = validate_bam(args.bam, mode=args.mode, max_errors=args.max_errors)
    print(f"Checked records: {checked}")

    if n_errors == 0:
        print("PASS: all checked records satisfy XM/XR/XG rules")
        return 0

    print(f"FAIL: found {n_errors} issue(s), showing up to {args.max_errors}")
    for err in errors:
        print(f"- {err.qname}: {err.reason}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
