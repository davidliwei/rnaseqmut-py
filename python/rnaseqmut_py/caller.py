from __future__ import annotations

from typing import Iterable

from .models import Variant

try:
    import pysam
except ImportError as e:
    raise RuntimeError("pysam is required for rnaseqmut-py") from e


def has_indel(read: pysam.AlignedSegment) -> bool:
    return any(op in (1, 2) for op, _ in (read.cigartuples or []))


def parse_md_mismatches(read: pysam.AlignedSegment) -> list[tuple[int, str, str, int]]:
    if read.is_unmapped or read.query_sequence is None or not read.has_tag("MD"):
        return []
    md = read.get_tag("MD")
    seq = read.query_sequence
    aligned_pairs = read.get_aligned_pairs(matches_only=False)
    ref_to_q = {r + 1: q for q, r in aligned_pairs if q is not None and r is not None}

    out: list[tuple[int, str, str, int]] = []
    ref_pos = read.reference_start + 1
    i = 0
    while i < len(md):
        if md[i].isdigit():
            j = i
            while j < len(md) and md[j].isdigit():
                j += 1
            ref_pos += int(md[i:j])
            i = j
            continue
        if md[i] == '^':
            j = i + 1
            while j < len(md) and md[j].isalpha():
                j += 1
            ref_pos += j - (i + 1)
            i = j
            continue
        if md[i].isalpha():
            ref = md[i]
            qpos = ref_to_q.get(ref_pos)
            if qpos is not None and 0 <= qpos < len(seq):
                alt = seq[qpos]
                out.append((ref_pos, ref, alt, qpos))
            ref_pos += 1
            i += 1
            continue
        i += 1
    return out


def discover_variants(bam: pysam.AlignmentFile, args) -> list[Variant]:
    found: set[Variant] = set()
    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        if read.is_paired and not read.is_proper_pair:
            continue
        if (not args.with_indel_read) and has_indel(read):
            continue
        try:
            if read.get_tag("NM") > args.max_mismatch:
                continue
        except KeyError:
            pass
        for pos1, ref, alt, qpos in parse_md_mismatches(read):
            if qpos < args.mut_span or (read.query_length - qpos) < args.mut_span:
                continue
            if (not args.output_n) and alt == "N":
                continue
            if (not args.with_indel) and (len(ref) != 1 or len(alt) != 1):
                continue
            chrom = bam.get_reference_name(read.reference_id)
            found.add(Variant(chrom, pos1, ref, alt))
    return sorted(found, key=lambda v: (v.chrom, v.pos, v.ref, v.alt))


def load_given_variants(path: str) -> list[Variant]:
    variants: list[Variant] = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c, p, r, a = line.strip().split()[:4]
            variants.append(Variant(c, int(p), r, a))
    return variants


def count_variant(bam: pysam.AlignmentFile, var: Variant, args) -> tuple[int, int, int, int]:
    ref_f = ref_r = alt_f = alt_r = 0
    for col in bam.pileup(var.chrom, var.pos - 1, var.pos, truncate=True, stepper="all"):
        if col.reference_pos != var.pos - 1:
            continue
        for pr in col.pileups:
            read = pr.alignment
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_paired and not read.is_proper_pair:
                continue
            if (not args.with_indel_read) and has_indel(read):
                continue
            try:
                if read.get_tag("NM") > args.max_mismatch:
                    continue
            except KeyError:
                pass
            qpos = pr.query_position
            if qpos < args.mut_span or (read.query_length - qpos) < args.mut_span:
                continue
            base = read.query_sequence[qpos]
            fwd = not read.is_reverse
            if base == var.ref:
                if fwd:
                    ref_f += 1
                else:
                    ref_r += 1
            elif base == var.alt:
                if fwd:
                    alt_f += 1
                else:
                    alt_r += 1
    return ref_f, ref_r, alt_f, alt_r
