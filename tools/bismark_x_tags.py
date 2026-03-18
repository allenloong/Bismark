#!/usr/bin/env python3
"""Generate and validate Bismark XM/XR/XG tags from BAM + reference FASTA.

Goal: reproduce Bismark XM exactly from alignment, XR/XG and reference sequence.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Callable, Iterable

VALID_XM_CHARS = set(".zZxXhHuU")
VALID_CONVERSIONS = {"CT", "GA"}


def _require_pysam():
    try:
        import pysam  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("pysam is required for BAM/FASTA I/O. Install with: pip install pysam") from exc
    return pysam


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


def _ctx_symbol_ct(methylated: bool, d1: str | None, d2: str | None) -> str:
    if d1 == "G":
        return "Z" if methylated else "z"
    if d1 in (None, "N", "X"):
        return "U" if methylated else "u"
    if d2 == "G":
        return "X" if methylated else "x"
    if d2 in (None, "N", "X"):
        return "U" if methylated else "u"
    return "H" if methylated else "h"


def _ctx_symbol_ga(methylated: bool, u1: str | None, u2: str | None) -> str:
    if u1 == "C":
        return "Z" if methylated else "z"
    if u1 in (None, "N", "X"):
        return "U" if methylated else "u"
    if u2 == "C":
        return "X" if methylated else "x"
    if u2 in (None, "N", "X"):
        return "U" if methylated else "u"
    return "H" if methylated else "h"


def generate_xm_from_alignment(
    query_seq: str,
    xr: str,
    qpos_to_rpos: list[int | None],
    is_reverse: bool,
    ref_base_fetcher: Callable[[int], str | None],
) -> str:
    """Generate XM in BAM SEQ orientation.

    `qpos_to_rpos` must be length == query length and contain reference positions
    for each query base (None for insertion/soft-clipped/etc.).
    """
    if xr not in VALID_CONVERSIONS:
        raise ValueError(f"XR must be CT/GA, got {xr}")
    if len(qpos_to_rpos) != len(query_seq):
        raise ValueError("qpos_to_rpos length must equal query length")

    query_seq = query_seq.upper()
    out: list[str] = []

    def oriented_base(rpos: int, offset: int) -> str | None:
        pos = rpos - offset if is_reverse else rpos + offset
        b = ref_base_fetcher(pos)
        return b.upper() if b else None

    for qpos, qbase in enumerate(query_seq):
        rpos = qpos_to_rpos[qpos]
        if rpos is None:
            out.append(".")
            continue

        rbase = oriented_base(rpos, 0)
        if rbase is None:
            out.append(".")
            continue

        if xr == "CT":
            if qbase == rbase and rbase == "C":
                out.append(_ctx_symbol_ct(True, oriented_base(rpos, +1), oriented_base(rpos, +2)))
            elif rbase == "C" and qbase == "T":
                out.append(_ctx_symbol_ct(False, oriented_base(rpos, +1), oriented_base(rpos, +2)))
            else:
                out.append(".")
        else:  # xr == 'GA'
            if qbase == rbase and rbase == "G":
                out.append(_ctx_symbol_ga(True, oriented_base(rpos, -1), oriented_base(rpos, -2)))
            elif rbase == "G" and qbase == "A":
                out.append(_ctx_symbol_ga(False, oriented_base(rpos, -1), oriented_base(rpos, -2)))
            else:
                out.append(".")

    return "".join(out)


def qpos_to_rpos_from_cigar(record) -> list[int | None]:
    """Create query-position -> reference-position mapping from CIGAR.

    Includes soft-clips and insertions as None positions (as Bismark effectively treats
    these as padded X and emits '.' in XM).
    """
    if record.cigartuples is None:
        return [None] * len(record.query_sequence)

    mapping: list[int | None] = []
    qpos = 0
    rpos = record.reference_start

    for op, length in record.cigartuples:
        if op in (0, 7, 8):  # M, =, X
            for _ in range(length):
                mapping.append(rpos)
                qpos += 1
                rpos += 1
        elif op in (1, 4):  # I, S
            mapping.extend([None] * length)
            qpos += 1 * length
        elif op in (2, 3):  # D, N
            rpos += length
        elif op in (5, 6):  # H, P
            continue
        else:
            raise ValueError(f"Unsupported CIGAR op: {op}")

    # keep parity with query sequence length
    qlen = len(record.query_sequence or "")
    if len(mapping) < qlen:
        mapping.extend([None] * (qlen - len(mapping)))
    elif len(mapping) > qlen:
        mapping = mapping[:qlen]

    return mapping


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


def compare_xm_for_record(record, ref_fetcher: Callable[[str, int], str | None]) -> tuple[str, str]:
    xr = record.get_tag("XR")
    chrom = record.reference_name
    mapping = qpos_to_rpos_from_cigar(record)

    def fetch(pos: int) -> str | None:
        if pos < 0:
            return None
        return ref_fetcher(chrom, pos)

    expected = generate_xm_from_alignment(
        query_seq=record.query_sequence or "",
        xr=xr,
        qpos_to_rpos=mapping,
        is_reverse=record.is_reverse,
        ref_base_fetcher=fetch,
    )
    observed = record.get_tag("XM")
    return expected, observed


def validate_bam_with_reference(path: str, fasta_path: str, max_errors: int = 50) -> tuple[int, int, list[ValidationError]]:
    pysam = _require_pysam()
    checked = 0
    errors: list[ValidationError] = []

    with pysam.AlignmentFile(path, "rb") as bam, pysam.FastaFile(fasta_path) as fa:

        def fetch(chrom: str, pos: int) -> str | None:
            try:
                return fa.fetch(chrom, pos, pos + 1)
            except Exception:
                return None

        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue

            checked += 1
            for msg in validate_common(rec):
                errors.append(ValidationError(rec.query_name, msg))
                if len(errors) >= max_errors:
                    return checked, len(errors), errors

            if not rec.has_tag("XR") or not rec.has_tag("XM"):
                continue

            exp, obs = compare_xm_for_record(rec, fetch)
            if exp != obs:
                errors.append(ValidationError(rec.query_name, f"XM mismatch: expected={exp} observed={obs}"))
                if len(errors) >= max_errors:
                    return checked, len(errors), errors

    return checked, len(errors), errors


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate/validate Bismark XM from BAM + XR/XG + reference")
    parser.add_argument("bam", help="Input BAM (Bismark output)")
    parser.add_argument("--reference", required=True, help="Reference FASTA used for Bismark alignment")
    parser.add_argument("--max-errors", type=int, default=50)
    args = parser.parse_args()

    checked, n_errors, errors = validate_bam_with_reference(args.bam, args.reference, max_errors=args.max_errors)
    print(f"Checked records: {checked}")

    if n_errors == 0:
        print("PASS: generated XM is identical to BAM XM for all checked records")
        return 0

    print(f"FAIL: found {n_errors} issue(s), showing up to {args.max_errors}")
    for err in errors:
        print(f"- {err.qname}: {err.reason}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
