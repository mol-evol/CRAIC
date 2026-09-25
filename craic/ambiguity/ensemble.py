"""Alignment as an ensemble, not a single answer.

Re-align a region many ways (the sandbox's alternative engines / gap regimes *are*
an ensemble), then measure how often each pair of residues ends up aligned. The
co-occurrence matrix complements the pair-HMM posterior: the posterior is what one
model believes; co-occurrence is what a crowd of real aligners actually did.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np

from ..domain import Alignment, Level
from ..engines import AlignerEngine
from . import sandbox


def _res_at(row: str) -> List[int]:
    out, r = [], 0
    for ch in row:
        if ch == "-":
            out.append(-1)
        else:
            out.append(r)
            r += 1
    return out


def cooccurrence(
    aln: Alignment,
    seq_i: int,
    seq_j: int,
    col_start: int,
    col_stop: int,
    engines: Sequence[AlignerEngine],
    codon_aware: bool = False,
    table: int = 1,
    level: Level = Level.NT,
    params_by_key: Optional[dict] = None,
) -> Tuple[np.ndarray, str, str, int]:
    """Fraction of alternative alignments in which residue a of seq_i is aligned to
    residue b of seq_j. Returns (matrix [ni x nj], chars_i, chars_j, n_alternatives)."""
    alts = sandbox.alternatives(aln, col_start, col_stop, engines,
                                codon_aware=codon_aware, table=table, level=level,
                                params_by_key=params_by_key)
    id_i, id_j = aln.ids[seq_i], aln.ids[seq_j]
    sub_i = aln.rows[seq_i][col_start:col_stop].replace("-", "")
    sub_j = aln.rows[seq_j][col_start:col_stop].replace("-", "")
    ni, nj = len(sub_i), len(sub_j)
    co = np.zeros((ni, nj))
    used = 0
    for alt in alts:
        block = alt.block
        if id_i not in block.ids or id_j not in block.ids:
            continue
        ri = block.rows[block.ids.index(id_i)]
        rj = block.rows[block.ids.index(id_j)]
        if len(ri) != len(rj):
            continue
        ai, aj = _res_at(ri), _res_at(rj)
        for c in range(len(ri)):
            a, b = ai[c], aj[c]
            if 0 <= a < ni and 0 <= b < nj:
                co[a, b] += 1.0
        used += 1
    if used:
        co /= used
    return co, sub_i, sub_j, used
