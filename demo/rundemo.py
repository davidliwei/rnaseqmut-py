#!/usr/bin/env python3
from __future__ import annotations
import subprocess
from pathlib import Path
import sys


def run(cmd: str, cwd: Path):
    print(f"#### COMMAND LINE: {cmd}")
    subprocess.run(cmd, shell=True, check=True, cwd=str(cwd))


def main() -> int:
    demo_dir = Path(__file__).resolve().parent
    base = demo_dir.parent
    data_dir = demo_dir / "data"
    results = demo_dir / "results"
    results.mkdir(exist_ok=True)
    for p in results.glob("*"):
        if p.is_file() and (p.name.endswith(".txt") or p.name.endswith(".vcf")):
            p.unlink()

    labels = "NORMAL1,NORMAL2,TUMOR1,TUMOR2"
    bam_files = sorted(data_dir.glob("*.bam"))

    print("\n####### Step 1, de-novo mutation calling ##########")
    for bam in bam_files:
        out = results / f"{bam.name}.1st.txt"
        run(f"python3 {base/'bin'/'rnaseqmut'} {bam} > {out}", demo_dir)

    print("\n####### Step 2, merging mutations in Step 1 into a candidate mutation list ##########")
    run(f"{base/'script'/'merge1stfile'} {results}/*.1st.txt > {results/'ALLMUTLIST.txt'}", demo_dir)

    print("\n####### Step 3, calling mutations again using the given list in Step 2 ##########")
    for bam in bam_files:
        out = results / f"{bam.name}.2nd.txt"
        run(f"python3 {base/'bin'/'rnaseqmut'} -l {results/'ALLMUTLIST.txt'} {bam} > {out}", demo_dir)

    print("\n####### Step 4, merging mutations in Step 3 into a big table ##########")
    run(f"python3 {base/'script'/'merge2ndvcf.py'} -l {labels} {results}/*.2nd.txt > {results/'ALLMUT.txt'}", demo_dir)

    print("\n####### Step 5, custom filtering based on mutations in Step 4 ##########")
    run(
        f"python3 {base/'script'/'filtermut.py'} -d 10 -f 0.0 -b 0 -c 0,1 -l {labels} < {results/'ALLMUT.txt'} > {results/'ALLMUT_FILTERED.vcf'}",
        demo_dir,
    )
    print("\n####### DEMO completed successfully. Check results/ALLMUT_FILTERED.vcf ##########")
    return 0


if __name__ == "__main__":
    sys.exit(main())
