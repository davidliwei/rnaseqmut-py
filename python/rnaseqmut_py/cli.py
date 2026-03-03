from __future__ import annotations

import argparse
import sys

import pysam

from . import __version__
from .caller import count_variant, discover_variants, load_given_variants
from .io import sanitize_sample_name, write_vcf_header


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="rnaseqmut")
    p.add_argument("bam_file")
    p.add_argument("-t", "--use_mdtag", action="store_true")
    p.add_argument("-r", "--ref_fasta", default="")
    p.add_argument("-l", "--mutation_list", default="")
    p.add_argument("-k", "--with_indel_read", action="store_true")
    p.add_argument("-d", "--with_indel", action="store_true")
    p.add_argument("-s", "--max_mismatch", type=int, default=1)
    p.add_argument("-n", "--output_n", action="store_true")
    p.add_argument("-i", "--min_read", type=int, default=1)
    p.add_argument("-m", "--mut_span", type=int, default=4)
    p.add_argument("--vcf_output", default="")
    p.add_argument("--version", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.version:
        print(__version__)
        return 0

    bam = pysam.AlignmentFile(args.bam_file, "rb")
    if not bam.has_index():
        print("Error: BAM index (.bai) is required", file=sys.stderr)
        return 1

    variants = load_given_variants(args.mutation_list) if args.mutation_list else discover_variants(bam, args)

    vcf_path = args.vcf_output or f"{args.bam_file}.vcf"
    with open(vcf_path, "w") as vcf:
        sample = sanitize_sample_name(args.bam_file)
        write_vcf_header(vcf, sample, source=f"rnaseqmut-py-{__version__}")

        for var in variants:
            ref_f, ref_r, alt_f, alt_r = count_variant(bam, var, args)
            alt_count = alt_f + alt_r
            if (not args.mutation_list) and alt_count < args.min_read:
                continue
            print(f"{var.chrom}\t{var.pos}\t{var.ref}\t{var.alt}\t{ref_f}\t{ref_r}\t{alt_f}\t{alt_r}")
            dp = ref_f + ref_r + alt_f + alt_r
            vcf.write(
                f"{var.chrom}\t{var.pos}\t.\t{var.ref}\t{var.alt}\t.\tPASS\tDP={dp}\t"
                f"GT:DP:AD:ADF:ADR\t0/1:{dp}:{ref_f+ref_r},{alt_count}:{ref_f},{alt_f}:{ref_r},{alt_r}\n"
            )

    bam.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
