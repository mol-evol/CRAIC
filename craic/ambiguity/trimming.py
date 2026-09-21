"""Model-free column trimming and sequence-outlier detection.

Faithful reimplementations of the algorithms behind widely-used tools, offered as
complements to CRAIC's model-based posterior reliability (which is derived from
the pair-HMM). These use only gaps and residue conservation -- no probabilistic
model -- so they are cheap and always available:

* ``gap_score`` / ``gap_threshold_mask``  -- MSA_trimmer-style gappiness.
* ``similarity_score``                     -- per-column conservation.
* ``gappyout_mask`` / ``trimal_strict_mask`` -- trimAl-style automatic trimming.
* ``gblocks_mask``                         -- Gblocks-style conserved blocks.
* ``sequence_scores`` / ``outlier_sequences`` -- EvalMSA / OD-seq-style outliers.

These are reimplementations of the published *algorithms*, not the original code,
and are not guaranteed identical to the reference tools; for authoritative
filtering, run the originals (trimAl, Gblocks, EvalMSA).
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from ..domain import Alignment, Alphabet


# --------------------------------------------------------------------------- #
# Encoding + similarity
# --------------------------------------------------------------------------- #

def _char_matrix(aln: Alignment) -> np.ndarray:
    return np.array([list(r) for r in aln.rows])


def _similarity(alphabet) -> Tuple[dict, np.ndarray]:
    """(index, normalised [0,1] pairwise similarity matrix).

    BLOSUM62 deliberately, even though CRAIC's own emission model uses BLOSUM45:
    this is the trimAl similarity score, reimplemented to match what trimAl
    computes, and its matrix is part of that definition rather than a CRAIC
    parameter. Changing it here would make the trimming tracks disagree with the
    published rules they are named after.
    """
    if alphabet == Alphabet.PROTEIN:
        from Bio.Align import substitution_matrices
        bl = substitution_matrices.load("BLOSUM62")
        syms = list(bl.alphabet)
        arr = np.array(bl, dtype=float)
        lo, hi = arr.min(), arr.max()
        norm = (arr - lo) / (hi - lo) if hi > lo else np.eye(len(syms))
        return {c: i for i, c in enumerate(syms)}, norm
    syms = "ACGT"
    return {c: i for i, c in enumerate(syms)}, np.eye(4)


def _encode(aln: Alignment):
    """(n, L) int array of residue indices, -1 for gap / unknown."""
    idx, sim = _similarity(aln.alphabet)
    nuc = aln.alphabet != Alphabet.PROTEIN
    n, L = aln.n_seqs, aln.length
    E = np.full((n, L), -1, dtype=int)
    for i, row in enumerate(aln.rows):
        for c, ch in enumerate(row):
            if ch == "-":
                continue
            u = ch.upper()
            if nuc and u == "U":
                u = "T"
            j = idx.get(u)
            if j is not None:
                E[i, c] = j
    return E, sim


# --------------------------------------------------------------------------- #
# Per-column scores
# --------------------------------------------------------------------------- #

def gap_score(aln: Alignment) -> np.ndarray:
    """Fraction of non-gap residues per column (trimAl gap score). 1 = no gaps."""
    M = _char_matrix(aln)
    return (M != "-").mean(axis=0)


def similarity_score(aln: Alignment) -> np.ndarray:
    """Mean pairwise residue similarity per column in [0, 1] (trimAl-style).
    Gaps excluded; a column with 0/1 residues scores 0/1."""
    E, sim = _encode(aln)
    n, L = E.shape
    out = np.zeros(L)
    for c in range(L):
        col = E[:, c]
        present = col[col >= 0]
        if present.size < 2:
            out[c] = 1.0 if present.size == 1 else 0.0
            continue
        sub = sim[np.ix_(present, present)]
        iu = np.triu_indices(present.size, 1)
        out[c] = float(sub[iu].mean())
    return out


# --------------------------------------------------------------------------- #
# Column masks (keep = True)
# --------------------------------------------------------------------------- #

def gap_threshold_mask(aln: Alignment, max_gap: float = 0.5) -> np.ndarray:
    """MSA_trimmer-style: keep columns whose gap fraction <= ``max_gap``."""
    return (1.0 - gap_score(aln)) <= max_gap + 1e-9


def _min_block(keep: np.ndarray, min_len: int) -> np.ndarray:
    """Drop kept runs shorter than ``min_len``."""
    out = keep.copy()
    L = len(keep)
    i = 0
    while i < L:
        if keep[i]:
            j = i
            while j < L and keep[j]:
                j += 1
            if j - i < min_len:
                out[i:j] = False
            i = j
        else:
            i += 1
    return out


def gappyout_mask(aln: Alignment, max_removed: float = 0.5) -> np.ndarray:
    """trimAl gappyout: keep the low-gap mode of the (bimodal) gap-score
    distribution, cutting at its largest jump.

    ``max_removed`` bounds the search for that jump. Without it, a distribution
    whose largest jump sits at the *top* — e.g. a hundred columns at gap-score
    0.99 and one at 1.0 — puts the threshold above nearly every column and trims
    the alignment away. gappyout is meant to separate a gappy mode from a clean
    one, so a cut that discards most of the alignment is a sign that there is no
    such separation, not a valid cut. trimAl's own implementation reads the
    slope of the whole distribution; bounding the cut is the simpler guard with
    the same effect on this failure case.
    """
    g = gap_score(aln)
    order = np.sort(g)
    L = order.size
    if L < 3:
        return g >= np.median(g) - 1e-9
    diffs = np.diff(order)
    # A cut after sorted position i keeps L - (i + 1) columns; restrict i so that
    # at least (1 - max_removed) of the columns survive.
    last = max(0, int(np.floor(L * max_removed)) - 1)
    cut = int(np.argmax(diffs[:last + 1]))
    thr = (order[cut] + order[cut + 1]) / 2.0
    return g >= thr


def trimal_strict_mask(aln: Alignment, gap_thr: float = 0.5,
                       sim_thr: float = 0.35, min_block: int = 3) -> np.ndarray:
    """trimAl strict-style: keep columns passing both a gap and a similarity
    threshold, then remove blocks shorter than ``min_block``."""
    keep = ((1.0 - gap_score(aln)) <= gap_thr) & (similarity_score(aln) >= sim_thr)
    return _min_block(keep, min_block)


def gblocks_mask(aln: Alignment, b1: Optional[int] = None, b2: Optional[int] = None,
                 b3: int = 8, b4: int = 10, allow_gaps: str = "none") -> np.ndarray:
    """Gblocks-style conserved blocks (Castresana 2000).

    b1 = min sequences for a conserved column (default >50%), b2 = for a flanking
    column (default 85%), b3 = max run of non-conserved columns, b4 = min block
    length, allow_gaps in {none, half, all}.
    """
    n, L = aln.n_seqs, aln.length
    if b1 is None:
        b1 = n // 2 + 1
    if b2 is None:
        b2 = int(np.ceil(0.85 * n))
    max_gap = {"none": 0.0, "half": 0.5, "all": 1.0}.get(allow_gaps, 0.0)
    M = _char_matrix(aln)
    nongap_ok = (1.0 - gap_score(aln)) <= max_gap + 1e-9

    conserved = np.zeros(L, bool)
    flank = np.zeros(L, bool)
    for c in range(L):
        col = [ch.upper() for ch in M[:, c] if ch != "-"]
        if not col:
            continue
        mx = int(np.unique(col, return_counts=True)[1].max())
        conserved[c] = mx >= b1
        flank[c] = mx >= b2

    good = conserved & nongap_ok
    flank = flank & nongap_ok

    # reject runs of >b3 contiguous non-good columns; shorter runs stay inside blocks
    eligible = good.copy()
    i = 0
    while i < L:
        if not good[i]:
            j = i
            while j < L and not good[j]:
                j += 1
            eligible[i:j] = (j - i) <= b3
            i = j
        else:
            i += 1

    # contiguous eligible segments, trimmed to flank positions, length >= b4
    keep = np.zeros(L, bool)
    i = 0
    while i < L:
        if eligible[i]:
            j = i
            while j < L and eligible[j]:
                j += 1
            s, e = i, j - 1
            while s < j and not flank[s]:
                s += 1
            while e >= s and not flank[e]:
                e -= 1
            if e >= s and (e - s + 1) >= b4:
                keep[s:e + 1] = True
            i = j
        else:
            i += 1
    return keep


# --------------------------------------------------------------------------- #
# Sequence-level outliers (EvalMSA / OD-seq-style)
# --------------------------------------------------------------------------- #

def sequence_scores(aln: Alignment) -> np.ndarray:
    """Per-sequence agreement: mean pairwise column similarity of each sequence
    against the rest, in [0, 1]. Low = outlier."""
    E, sim = _encode(aln)
    n, L = E.shape
    scores = np.zeros(n)
    for i in range(n):
        tot, cnt = 0.0, 0
        ei = E[i]
        for j in range(n):
            if i == j:
                continue
            ej = E[j]
            both = (ei >= 0) & (ej >= 0)
            if both.any():
                tot += float(sim[ei[both], ej[both]].mean())
                cnt += 1
        scores[i] = tot / cnt if cnt else 0.0
    return scores


def outlier_sequences(aln: Alignment, z: float = 3.0) -> Tuple[np.ndarray, np.ndarray]:
    """Return (per-sequence score, boolean outlier flag). A sequence is an outlier
    if its score is more than ``z`` robust (MAD) standard deviations below the median."""
    s = sequence_scores(aln)
    med = float(np.median(s))
    mad = float(np.median(np.abs(s - med)))
    if mad <= 0.0:
        # No robust spread at all: the sequences agree with each other equally
        # well (near-identical input is the common case). Substituting a tiny
        # epsilon here would put the threshold immediately below the median and
        # flag every below-median sequence as an outlier. With no scale there is
        # no outlier to find, so report none.
        return s, np.zeros(s.shape, dtype=bool)
    flags = s < med - z * 1.4826 * mad
    return s, flags
