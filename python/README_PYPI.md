# rnaseqmut-py

Python-only port of the `rnaseqmut` mutation calling workflow.

- No C++ compiler needed
- Input: sorted/indexed BAM (`.bam` + `.bai`)
- Output: tab-delimited calls on STDOUT and VCF file

## Install

```bash
pip install rnaseqmut-py
```

## Usage

```bash
rnaseqmut sample.bam > sample.txt
```

```bash
rnaseqmut -l mutation_list.txt sample.bam > sample.2nd.txt
```
