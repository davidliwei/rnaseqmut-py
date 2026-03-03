from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int  # 1-based
    ref: str
    alt: str
