"""Posterior explorer: the pair-HMM posterior matrix as alternative homologies.

Instead of one hard column assignment, show P(residue i of A ~ residue j of B)
for a chosen pair of sequences over a chosen column window. Bright off-diagonal
mass = the alignment has plausible alternative homology hypotheses there.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .. import accel, domain
from ..domain import Alignment, Alphabet


def model_for(aln: Alignment, identity: float = 0.9) -> accel.EmissionModel:
    kind = "protein" if aln.alphabet == Alphabet.PROTEIN else "dna"
    return accel.emission_model(kind, identity=identity)


# Probing one residue, or profiling one sequence, averages a full pair-HMM
# posterior against every *other* sequence — each costs about li*lj cells. On a
# big alignment that is thousands of full matrices, which can freeze the app and
# thrash memory. Cap the total work by averaging over an evenly-spaced sample of
# partners instead; for small alignments this is every other sequence (unchanged).
# The Rust core is ~100x faster than the NumPy fallback, so the budget tracks the
# backend that's actually installed.
MAX_PROBE_CELLS = 30_000_000 if accel.HAVE_RUST else 1_000_000


def _partner_indices(n_seqs: int, seq_i: int, li: int) -> List[int]:
    partners = [j for j in range(n_seqs) if j != seq_i]
    budget = max(1, MAX_PROBE_CELLS // max(1, li * li))
    if len(partners) <= budget:
        return partners
    idx = np.linspace(0, len(partners) - 1, budget).round().astype(int)
    return [partners[k] for k in dict.fromkeys(idx.tolist())]


@dataclass
class PosteriorView:
    seq_i: int
    seq_j: int
    matrix: np.ndarray            # (len_i, len_j) restricted to the window
    res_i: List[int]             # original residue indices for rows
    res_j: List[int]             # original residue indices for cols
    chars_i: str
    chars_j: str
    col_start: int
    col_stop: int


def region_posterior(
    aln: Alignment,
    seq_i: int,
    seq_j: int,
    col_start: int,
    col_stop: int,
    model: Optional[accel.EmissionModel] = None,
    delta: float = 0.02,
    epsilon: float = 0.5,
    mode: str = "posterior",
) -> PosteriorView:
    """Posterior submatrix for the two sequences over columns [start, stop).

    mode="posterior" shows the raw posterior; mode="residual" zeroes the cells
    the current alignment already commits to (residues sharing a column), so the
    confident diagonal disappears and only the *alternative* homology mass — the
    genuinely ambiguous part — remains visible.
    """
    model = model or model_for(aln)
    row_i = aln.rows[seq_i]
    row_j = aln.rows[seq_j]
    map_i = domain.residue_index(row_i)
    map_j = domain.residue_index(row_j)
    ri = [map_i[c] for c in range(col_start, col_stop) if map_i[c] >= 0]
    rj = [map_j[c] for c in range(col_start, col_stop) if map_j[c] >= 0]
    if not ri or not rj:
        return PosteriorView(seq_i, seq_j, np.zeros((len(ri), len(rj))), ri, rj,
                             "", "", col_start, col_stop)
    # Window the pair-HMM to the region's residues (plus a context margin) rather
    # than computing the full sequence-length matrix — keeps a drag-selection on a
    # long alignment fast and bounded, independent of total sequence length.
    seq_i_u = row_i.replace("-", "")
    seq_j_u = row_j.replace("-", "")
    margin = 60
    lo_i, hi_i = max(0, ri[0] - margin), min(len(seq_i_u), ri[-1] + 1 + margin)
    lo_j, hi_j = max(0, rj[0] - margin), min(len(seq_j_u), rj[-1] + 1 + margin)
    win = accel.posterior_matrix(seq_i_u[lo_i:hi_i], seq_j_u[lo_j:hi_j],
                                 model, delta=delta, epsilon=epsilon)
    sub = win[np.ix_([r - lo_i for r in ri], [c - lo_j for c in rj])].copy()
    if mode == "residual":
        rc_i = domain.column_index(row_i)
        rc_j = domain.column_index(row_j)
        for a, ra in enumerate(ri):
            ca = rc_i[ra]
            for b, rb in enumerate(rj):
                if rc_j[rb] == ca:
                    sub[a, b] = 0.0
    chars_i = "".join(row_i[c] for c in range(col_start, col_stop) if map_i[c] >= 0)
    chars_j = "".join(row_j[c] for c in range(col_start, col_stop) if map_j[c] >= 0)
    return PosteriorView(seq_i, seq_j, sub, ri, rj, chars_i, chars_j, col_start, col_stop)


def region_profile(
    aln: Alignment,
    seq_i: int,
    col_start: int,
    col_stop: int,
    model: Optional[accel.EmissionModel] = None,
) -> PosteriorView:
    """One sequence vs the rest: rows are residues of seq_i, columns are
    alignment columns, value = mean posterior that the residue is homologous to
    whatever sits in that column across all other sequences.

    Unlike the pairwise view this reflects the whole MSA, so a residue's row
    shows its committed column *and* any other columns it could plausibly belong
    to — the alignment's own uncertainty, projected onto its columns.
    """
    model = model or model_for(aln)
    row_i = aln.rows[seq_i]
    map_i = domain.residue_index(row_i)
    ri = [map_i[c] for c in range(col_start, col_stop) if map_i[c] >= 0]
    if not ri:
        return PosteriorView(seq_i, -1, np.zeros((0, 0)), ri, [], "", "",
                             col_start, col_stop)
    L = aln.length
    li = sum(1 for ch in row_i if ch != "-") or 1
    ri_arr = np.array(ri)
    acc = np.zeros((len(ri), L))
    cnt = np.zeros(L)
    for j in _partner_indices(aln.n_seqs, seq_i, li):
        P = accel.posterior_matrix(row_i, aln.rows[j], model)
        if P.size == 0:
            continue
        res_at = np.array(domain.residue_index(aln.rows[j]))
        cols_valid = np.where(res_at >= 0)[0]
        if cols_valid.size == 0:
            continue
        acc[:, cols_valid] += P[np.ix_(ri_arr, res_at[cols_valid])]
        cnt[cols_valid] += 1
    reg = slice(col_start, col_stop)
    denom = cnt[reg].copy()
    denom[denom == 0] = 1.0
    sub = acc[:, reg] / denom
    sub[:, cnt[reg] == 0] = 0.0
    chars_i = "".join(row_i[c] for c in range(col_start, col_stop) if map_i[c] >= 0)
    res_cols = list(range(col_start, col_stop))
    return PosteriorView(seq_i, -1, sub, ri, res_cols, chars_i, "", col_start, col_stop)


def probe_residue(
    aln: Alignment,
    seq_i: int,
    res_index: int,
    model: Optional[accel.EmissionModel] = None,
    window: int = 120,
) -> np.ndarray:
    """'Where else could this residue go?' For one residue of seq_i, the mean
    posterior (over a sample of other sequences) that it is homologous to each
    alignment column. Returns an array of length aln.length (nan where no data).

    A residue can only realistically shift a short distance, so the pair-HMM is
    restricted to a ``window`` of columns either side of the residue's own column
    (with a context margin). That keeps each comparison O(window^2) instead of
    O(sequence_length^2), so a click stays responsive on long alignments instead
    of trying to score homology against the entire sequence.
    """
    model = model or model_for(aln)
    row_i = aln.rows[seq_i]
    L = aln.length
    cols_of_res = domain.column_index(row_i)
    if res_index not in cols_of_res:
        raise IndexError(
            f"sequence {seq_i} has {len(cols_of_res)} residues; "
            f"there is no residue {res_index} to probe"
        )
    ci = cols_of_res[res_index]                        # the residue's own column
    c0, c1 = max(0, ci - window), min(L, ci + window + 1)
    seq_i_u = row_i.replace("-", "")
    margin = 60
    lo_i = max(0, res_index - window - margin)
    hi_i = min(len(seq_i_u), res_index + window + margin + 1)
    sub_i = seq_i_u[lo_i:hi_i]
    local_i = res_index - lo_i                          # clicked residue within sub_i
    acc = np.zeros(L)
    cnt = np.zeros(L)
    for j in _partner_indices(aln.n_seqs, seq_i, hi_i - lo_i):
        row_j = aln.rows[j]
        map_j = domain.residue_index(row_j)
        cols_j = [c for c in range(c0, c1) if map_j[c] >= 0]
        if not cols_j:
            continue
        rj = [map_j[c] for c in cols_j]
        seq_j_u = row_j.replace("-", "")
        lo_j = max(0, rj[0] - margin)
        hi_j = min(len(seq_j_u), rj[-1] + 1 + margin)
        P = accel.posterior_matrix(sub_i, seq_j_u[lo_j:hi_j], model)
        if local_i >= P.shape[0]:
            continue
        prow = P[local_i]
        for c, r in zip(cols_j, rj):
            jl = r - lo_j
            if 0 <= jl < prow.shape[0]:
                acc[c] += prow[jl]
                cnt[c] += 1
    out = np.full(L, np.nan)
    nz = cnt > 0
    out[nz] = acc[nz] / cnt[nz]
    return out


def ambiguity_profile(
    aln: Alignment,
    seq_i: int,
    seq_j: int,
    model: Optional[accel.EmissionModel] = None,
) -> np.ndarray:
    """Per-residue ambiguity of seq_i vs seq_j = 1 - max posterior over partners.

    1.0 means residue i has no confident homolog in j (highly ambiguous).
    """
    model = model or model_for(aln)
    n_res_i = sum(1 for ch in aln.rows[seq_i] if ch != "-")
    P = accel.posterior_matrix(aln.rows[seq_i], aln.rows[seq_j], model)
    if P.size == 0:
        # An empty partner means no evidence about any residue of i, not "no
        # residues" — return one fully-ambiguous entry per residue so the shape
        # contract holds for callers that index by residue.
        return np.ones(n_res_i)
    return 1.0 - P.max(axis=1)
