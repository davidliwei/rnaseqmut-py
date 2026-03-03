from __future__ import annotations

import subprocess
from pathlib import Path


def run(cmd: str, cwd: Path):
    print(cmd)
    subprocess.run(cmd, shell=True, check=True, cwd=str(cwd))


def main() -> int:
    project = Path(__file__).resolve().parents[2]
    demo = project / "demo"
    data = demo / "data"
    results = demo / "results"
    results.mkdir(exist_ok=True)

    for p in results.glob("*.txt"):
        p.unlink()
    for p in results.glob("*.vcf"):
        p.unlink()

    labels = "NORMAL1,NORMAL2,TUMOR1,TUMOR2"
    bams = sorted(data.glob("*.bam"))

    runner = f"PYTHONPATH={project/'python'} python3 -m rnaseqmut_py.cli"

    for bam in bams:
        out = results / f"{bam.name}.1st.txt"
        run(f"{runner} {bam} > {out}", demo)

    run(f"{project/'script'/'merge1stfile'} {results}/*.1st.txt > {results/'ALLMUTLIST.txt'}", demo)

    for bam in bams:
        out = results / f"{bam.name}.2nd.txt"
        run(f"{runner} -l {results/'ALLMUTLIST.txt'} {bam} > {out}", demo)

    run(f"python3 {project/'script'/'merge2ndvcf.py'} -l {labels} {results}/*.2nd.txt > {results/'ALLMUT.txt'}", demo)

    run(
        f"python3 {project/'script'/'filtermut.py'} -d 10 -f 0.0 -b 0 -c 0,1 -l {labels} < {results/'ALLMUT.txt'} > {results/'ALLMUT_FILTERED.vcf'}",
        demo,
    )

    assert (results / "ALLMUT_FILTERED.vcf").exists(), "Demo output not found"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
