#!/usr/bin/env python3
"""Generate and validate Bismark XM/XR/XG tags from BAM + reference FASTA."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Callable

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


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def strand_id_from_conversions(xr: str, xg: str) -> str:
    mapping = {
        ("CT", "CT"): "OT",
        ("GA", "GA"): "CTOB",
        ("GA", "CT"): "CTOT",
        ("CT", "GA"): "OB",
    }
    return mapping.get((xr, xg), "INVALID")


def qpos_to_rpos_from_cigar(record) -> list[int | None]:
    if record.cigartuples is None:
        return [None] * len(record.query_sequence or "")

    mapping: list[int | None] = []
    rpos = record.reference_start

    for op, length in record.cigartuples:
        if op in (0, 7, 8):  # M, =, X
            for _ in range(length):
                mapping.append(rpos)
                rpos += 1
        elif op in (1, 4):  # I, S
            mapping.extend([None] * length)
        elif op in (2, 3):  # D, N
            rpos += length
        elif op in (5, 6):  # H, P
            continue
        else:
            raise ValueError(f"Unsupported CIGAR op: {op}")

    qlen = len(record.query_sequence or "")
    if len(mapping) < qlen:
        mapping.extend([None] * (qlen - len(mapping)))
    elif len(mapping) > qlen:
        mapping = mapping[:qlen]
    return mapping


def _methylation_call_core(seq: str, genomic: str, xr: str) -> str:
    seq_bases = list(seq.upper())
    g = list(genomic.upper())
    out: list[str] = []

    if xr == "CT":
        for i, b in enumerate(seq_bases):
            gb = g[i]
            if b == gb:
                if gb == "C":
                    d1 = g[i + 1]
                    if d1 == "G":
                        out.append("Z")
                    elif d1 in ("N", "X"):
                        out.append("U")
                    else:
                        d2 = g[i + 2]
                        if d2 == "G":
                            out.append("X")
                        elif d2 in ("N", "X"):
                            out.append("U")
                        else:
                            out.append("H")
                else:
                    out.append(".")
            elif gb == "C" and b == "T":
                d1 = g[i + 1]
                if d1 == "G":
                    out.append("z")
                elif d1 in ("N", "X"):
                    out.append("u")
                else:
                    d2 = g[i + 2]
                    if d2 == "G":
                        out.append("x")
                    elif d2 in ("N", "X"):
                        out.append("u")
                    else:
                        out.append("h")
            else:
                out.append(".")
    elif xr == "GA":
        for i, b in enumerate(seq_bases):
            gb = g[i + 2]
            if b == gb:
                if gb == "G":
                    u1 = g[i + 1]
                    if u1 == "C":
                        out.append("Z")
                    elif u1 in ("N", "X"):
                        out.append("U")
                    else:
                        u2 = g[i]
                        if u2 == "C":
                            out.append("X")
                        elif u2 in ("N", "X"):
                            out.append("U")
                        else:
                            out.append("H")
                else:
                    out.append(".")
            elif gb == "G" and b == "A":
                u1 = g[i + 1]
                if u1 == "C":
                    out.append("z")
                elif u1 in ("N", "X"):
                    out.append("u")
                else:
                    u2 = g[i]
                    if u2 == "C":
                        out.append("x")
                    elif u2 in ("N", "X"):
                        out.append("u")
                    else:
                        out.append("h")
            else:
                out.append(".")
    else:
        raise ValueError(f"XR must be CT/GA, got {xr}")

    return "".join(out)


def generate_xm_from_alignment(
    query_seq: str,
    xr: str,
    qpos_to_rpos: list[int | None],
    is_reverse: bool,
    ref_base_fetcher: Callable[[int], str | None],
) -> str:
    """Reconstruct XM following Bismark methylation_call and SAM-output orientation."""
    if xr not in VALID_CONVERSIONS:
        raise ValueError(f"XR must be CT/GA, got {xr}")

    seq_bam = (query_seq or "").upper()
    if len(qpos_to_rpos) != len(seq_bam):
        raise ValueError("qpos_to_rpos length must equal query length")

    # Mirror Bismark SAM output path: methylation call is generated on internal
    # read orientation, while reverse-strand output is represented by reversed XM.
    if is_reverse:
        seq_work = revcomp(seq_bam)
        map_work = list(reversed(qpos_to_rpos))
    else:
        seq_work = seq_bam
        map_work = qpos_to_rpos

    # Build genomic sequence aligned per query position, using X for I/S (None mapping).
    aligned_g: list[str] = []
    mapped_positions: list[int] = []
    for r in map_work:
        if r is None:
            aligned_g.append("X")
        else:
            b = ref_base_fetcher(r)
            aligned_g.append((b or "N").upper())
            mapped_positions.append(r)

    if not mapped_positions:
        methcall = "." * len(seq_work)
    else:
        first = mapped_positions[0]
        last = mapped_positions[-1]

        # In work orientation (+), CT needs +2 at 3' end; GA needs +2 at 5' end.
        extra_5 = []
        extra_3 = []
        if xr == "CT":
            extra_3 = [ref_base_fetcher(last + 1) or "N", ref_base_fetcher(last + 2) or "N"]
            genomic = "".join(aligned_g + [b.upper() for b in extra_3])
        else:  # GA
            extra_5 = [ref_base_fetcher(first - 2) or "N", ref_base_fetcher(first - 1) or "N"]
            genomic = "".join([b.upper() for b in extra_5] + aligned_g)

        methcall = _methylation_call_core(seq_work, genomic, xr)

    return methcall[::-1] if is_reverse else methcall


def validate_common(record) -> list[str]:
    errs: list[str] = []
    if not record.has_tag("XM"):
        errs.append("missing XM tag")
    else:
        xm = record.get_tag("XM")
        if set(xm) - VALID_XM_CHARS:
            errs.append(f"XM contains invalid chars: {''.join(sorted(set(xm) - VALID_XM_CHARS))}")
        if len(xm) != len(record.query_sequence or ""):
            errs.append("XM length != query length")

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
    return expected, record.get_tag("XM")


def validate_bam_with_reference(path: str, fasta_path: str, max_errors: int = 50) -> tuple[int, int, list[ValidationError]]:
    pysam = _require_pysam()
    checked = 0
    errors: list[ValidationError] = []

    with pysam.AlignmentFile(path, "rb") as bam, pysam.FastaFile(fasta_path) as fa:

        fasta_refs = set(fa.references)
        chrom_cache: dict[str, str | None] = {}

        def resolve_chrom(chrom: str) -> str | None:
            if chrom in chrom_cache:
                return chrom_cache[chrom]
            candidates = [chrom]
            if chrom.startswith("chr"):
                candidates.append(chrom[3:])
            else:
                candidates.append(f"chr{chrom}")
            for c in candidates:
                if c in fasta_refs:
                    chrom_cache[chrom] = c
                    return c
            chrom_cache[chrom] = None
            return None

        def fetch(chrom: str, pos: int) -> str | None:
            resolved = resolve_chrom(chrom)
            if resolved is None:
                return None
            try:
                return fa.fetch(resolved, pos, pos + 1)
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
    parser = argparse.ArgumentParser(description="Reconstruct Bismark XM from BAM + XR/XG + reference")
    parser.add_argument("bam", help="Input BAM")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--max-errors", type=int, default=50)
    args = parser.parse_args()

    checked, n_errors, errors = validate_bam_with_reference(args.bam, args.reference, args.max_errors)
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
