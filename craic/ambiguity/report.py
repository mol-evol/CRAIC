"""Find and summarise the parts of an alignment a human should look at.

Turns the per-column reliability track into (1) a ranked list of ambiguous
*hotspots* to jump to, and (2) a one-glance confidence summary for the whole
alignment — the QC you'd attach to a submission.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np


@dataclass
class Hotspot:
    col_start: int          # inclusive (0-based column)
    col_stop: int           # exclusive
    mean_reliability: float

    @property
    def width(self) -> int:
        return self.col_stop - self.col_start

    @property
    def severity(self) -> float:
        return self.width * (1.0 - self.mean_reliability)


def hotspots(col_reliability: np.ndarray, threshold: float = 0.5,
             min_len: int = 2) -> List[Hotspot]:
    """Contiguous runs of columns below ``threshold``, ranked worst-first."""
    rel = np.asarray(col_reliability, dtype=float)
    out: List[Hotspot] = []
    c = 0
    n = len(rel)
    while c < n:
        v = rel[c]
        if v == v and v < threshold:
            start = c
            while c < n and rel[c] == rel[c] and rel[c] < threshold:
                c += 1
            if c - start >= min_len:
                out.append(Hotspot(start, c, float(np.nanmean(rel[start:c]))))
        else:
            c += 1
    out.sort(key=lambda h: -h.severity)
    return out


@dataclass
class Summary:
    n_columns: int
    scored_columns: int
    mean_reliability: float
    pct_high: float                  # fraction of scored columns >= threshold
    threshold: float
    hotspots: List[Hotspot] = field(default_factory=list)

    def as_text(self) -> str:
        lines = [
            f"Columns: {self.n_columns}  (scored: {self.scored_columns})",
            f"Mean reliability: {self.mean_reliability:.2f}",
            f"High-confidence columns (≥ {self.threshold:.2f}): {self.pct_high * 100:.0f}%",
            f"Ambiguous hotspots: {len(self.hotspots)}",
        ]
        for h in self.hotspots[:12]:
            lines.append(f"  columns {h.col_start + 1}–{h.col_stop}  "
                         f"(width {h.width}, mean {h.mean_reliability:.2f})")
        return "\n".join(lines)


def summary(col_reliability: np.ndarray, threshold: float = 0.5) -> Summary:
    rel = np.asarray(col_reliability, dtype=float)
    scored = rel[~np.isnan(rel)]
    pct = float((scored >= threshold).mean()) if scored.size else 0.0
    return Summary(
        n_columns=len(rel),
        scored_columns=int(scored.size),
        mean_reliability=float(scored.mean()) if scored.size else 0.0,
        pct_high=pct,
        threshold=threshold,
        hotspots=hotspots(rel, threshold),
    )
