from __future__ import annotations

import re
from pathlib import Path


def sanitize_sample_name(path: str) -> str:
    name = Path(path).stem
    name = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
    return name or "SAMPLE"


def write_vcf_header(out, sample: str, source: str) -> None:
    out.write("##fileformat=VCFv4.2\n")
    out.write(f"##source={source}\n")
    out.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at this locus (REF+ALT)\">\n")
    out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    out.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for this sample\">\n")
    out.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths for REF and ALT\">\n")
    out.write("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Forward-strand allele depths for REF and ALT\">\n")
    out.write("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Reverse-strand allele depths for REF and ALT\">\n")
    out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
