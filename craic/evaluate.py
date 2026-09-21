"""Accuracy of an alignment against a known-true reference.

Sum-of-pairs and total-column scores, per-column correctness, and the ROC AUC of
a score predicting that correctness. These used to live in ``benchmarks/`` where
only the benchmark harness could reach them; they are here so that the GUI's
truth mode, the command line and the benchmark all compute accuracy the same
way. A number shown to a student and a number printed in a paper being produced
by two different implementations is exactly the drift worth designing out.

A "reference" here is a trusted alignment of the *same sequences* — a structural
reference such as BAliBASE, or the known-true alignment of a simulated dataset.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .domain import (Alignment, IdMismatch, column_index, match_ids,
                     residue_index, ungap)


# --------------------------------------------------------------------------- #
# Alignment accuracy (SP / TC)
# --------------------------------------------------------------------------- #

def _pairs(rows: Sequence[str], core: Optional[np.ndarray] = None) -> set:
    """Set of homologous residue pairs ((seq_i, res_i), (seq_j, res_j)).

    ``core``, if given, is a boolean mask over the columns of ``rows``; only
    pairs asserted by core columns are counted.
    """
    ridx = [residue_index(r) for r in rows]
    pairs = set()
    for c in range(len(rows[0])):
        if core is not None and not core[c]:
            continue
        present = [(si, ridx[si][c]) for si in range(len(rows)) if rows[si][c] != "-"]
        for a in range(len(present)):
            for b in range(a + 1, len(present)):
                pairs.add((present[a], present[b]))
    return pairs


def _columns(rows: Sequence[str], core: Optional[np.ndarray] = None) -> set:
    """Each fully-specified column as a frozenset of (seq, residue-index)."""
    ridx = [residue_index(r) for r in rows]
    cols = set()
    for c in range(len(rows[0])):
        if core is not None and not core[c]:
            continue
        present = frozenset((si, ridx[si][c]) for si in range(len(rows)) if rows[si][c] != "-")
        if len(present) >= 2:
            cols.add(present)
    return cols


def sp_score(true_rows: Sequence[str], inf_rows: Sequence[str],
             core: Optional[np.ndarray] = None) -> Tuple[float, float]:
    """(recall, precision) of homologous residue pairs. Recall is the SP score.

    ``core`` is a boolean mask over the *reference* columns. Published BAliBASE
    SP figures are restricted to the reference's core blocks, so a score
    computed over all columns is not comparable with them — pass the core mask
    when one is available, and label the result accordingly when it is not.
    """
    T, I = _pairs(true_rows, core), _pairs(inf_rows)
    inter = len(T & I)
    recall = inter / len(T) if T else 1.0
    prec = inter / len(I) if I else 1.0
    return recall, prec


def tc_score(true_rows: Sequence[str], inf_rows: Sequence[str],
             core: Optional[np.ndarray] = None) -> float:
    """Total-column score: fraction of reference columns reproduced exactly.

    See :func:`sp_score` on ``core``.
    """
    T, I = _columns(true_rows, core), _columns(inf_rows)
    return len(T & I) / len(T) if T else 1.0


def col_correct_fraction(inf_rows: Sequence[str], true_rows: Sequence[str]) -> np.ndarray:
    """Per inferred column, the fraction of its asserted pairs that are true.

    ``nan`` for a column that asserts no pair at all (a singleton column): it
    makes no claim, so it can be neither right nor wrong. Those columns are not
    counted as errors anywhere downstream.
    """
    true_pairs = _pairs(true_rows)
    ridx = [residue_index(r) for r in inf_rows]
    out = []
    for c in range(len(inf_rows[0])):
        present = [(si, ridx[si][c]) for si in range(len(inf_rows)) if inf_rows[si][c] != "-"]
        pairs = [(present[a], present[b])
                 for a in range(len(present)) for b in range(a + 1, len(present))]
        out.append(np.mean([p in true_pairs for p in pairs]) if pairs else np.nan)
    return np.array(out)


# --------------------------------------------------------------------------- #
# Does a score predict correctness?
# --------------------------------------------------------------------------- #

def _midranks(values: np.ndarray) -> np.ndarray:
    """Ranks of ``values``, ties sharing their average rank.

    Reliability scores tie heavily — conserved columns pile up at 1.0 — and
    assigning tied values arbitrary distinct ranks makes the AUC depend on input
    order. The Mann-Whitney statistic is defined with midranks, so that is what
    is used.
    """
    order = values.argsort(kind="mergesort")
    ranks = np.empty(values.size, dtype=float)
    ranks[order] = np.arange(1, values.size + 1, dtype=float)
    srt = values[order]
    i = 0
    while i < srt.size:
        j = i + 1
        while j < srt.size and srt[j] == srt[i]:
            j += 1
        if j - i > 1:
            ranks[order[i:j]] = (i + 1 + j) / 2.0
        i = j
    return ranks


def auc(score: np.ndarray, positive: np.ndarray) -> float:
    """ROC AUC: probability a random positive outranks a random negative.

    Ties count as half, via midranks, so a score that cannot discriminate at all
    returns 0.5 rather than an order-dependent number.
    """
    score = np.asarray(score, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    n_pos, n_neg = int(positive.sum()), int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = _midranks(score)
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


# --------------------------------------------------------------------------- #
# One alignment against one reference
# --------------------------------------------------------------------------- #

@dataclass
class Accuracy:
    """How an inferred alignment compares with a trusted reference."""

    sp: float                      #: recall of the reference's homologous pairs
    precision: float               #: fraction of asserted pairs that are true
    tc: float                      #: fraction of reference columns reproduced exactly
    col_correct: np.ndarray        #: per inferred column, fraction of its pairs that are true
    n_scorable: int                #: columns that assert at least one pair

    @property
    def correct_mask(self) -> np.ndarray:
        """Columns that are entirely correct (nan columns count as not-correct,
        since they assert nothing to be correct about)."""
        return np.nan_to_num(self.col_correct, nan=0.0) >= 0.999

    def reliability_auc(self, scores: np.ndarray) -> float:
        """How well a per-column score predicts which columns are actually right."""
        scores = np.asarray(scores, dtype=float)
        ok = ~np.isnan(scores) & ~np.isnan(self.col_correct)
        if ok.sum() < 4:
            return float("nan")
        return auc(scores[ok], self.col_correct[ok] >= 0.999)


def matched_rows(aln: Alignment, reference: Alignment) -> List[str]:
    """The reference's rows, reordered onto ``aln``'s own sequence order.

    The two are matched on sequence id (tolerating the whitespace truncation some
    aligners apply), so that a reference whose sequences are listed differently —
    which is usual — is still usable. If the two do not describe the same
    sequences, or describe them with different residues, that is an error rather
    than a very low score: a reference for the wrong data would otherwise read as
    an alignment that is entirely wrong.
    """
    mapping = match_ids(aln, reference)                     # raises IdMismatch
    by_id = dict(zip(reference.ids, reference.rows))
    true_rows = [by_id[mapping[rid]] for rid in aln.ids]

    for rid, inferred, truth in zip(aln.ids, aln.rows, true_rows):
        if ungap(inferred).upper() != ungap(truth).upper():
            raise IdMismatch(
                f"sequence {rid!r} differs between the alignment and the reference "
                f"({len(ungap(inferred))} vs {len(ungap(truth))} residues); "
                "the reference must be an alignment of the same sequences"
            )
    return true_rows


def compare_to_reference(aln: Alignment, reference: Alignment,
                         core: Optional[np.ndarray] = None) -> Accuracy:
    """Score ``aln`` against a trusted ``reference`` alignment of the same sequences.

    See :func:`matched_rows` on how the two are paired up.
    """
    true_rows = matched_rows(aln, reference)
    inf_rows = list(aln.rows)
    sp, prec = sp_score(true_rows, inf_rows, core)
    tc = tc_score(true_rows, inf_rows, core)
    frac = col_correct_fraction(inf_rows, true_rows)
    return Accuracy(sp=sp, precision=prec, tc=tc, col_correct=frac,
                    n_scorable=int((~np.isnan(frac)).sum()))


# --------------------------------------------------------------------------- #
# Where every residue *should* have been
# --------------------------------------------------------------------------- #
#
# The column scores above say how wrong a column is, which is enough to find
# trouble but not enough to learn from it: a column at 0.8 names neither the
# residue that broke it nor the answer it got wrong. The inferred and the true
# alignment hold exactly the same residues, though, so every residue has a true
# column, and that one lookup is what everything below is built from.
#
# The two alignments cannot be shown side by side column-for-column — they have
# different numbers of columns — so the true column travels with the residue
# instead of with the screen position.


def true_column_grid(inf_rows: Sequence[str], true_rows: Sequence[str]) -> np.ndarray:
    """``(n_seqs, n_cols)`` of true-alignment columns, ``-1`` where a cell is a gap.

    Entry ``[s, c]`` is the column the *reference* puts the residue that the
    inferred alignment placed at row ``s``, column ``c``. Two residues are
    homologous exactly when their entries are equal, which is the whole of the
    comparison — every function below is a different way of reading this grid.
    """
    n, L = len(inf_rows), len(inf_rows[0])
    grid = np.full((n, L), -1, dtype=int)
    for s, (inferred, truth) in enumerate(zip(inf_rows, true_rows)):
        col_of = column_index(truth)                   # residue index -> true column
        for c, r in enumerate(residue_index(inferred)):
            if r >= 0:
                grid[s, c] = col_of[r]
    return grid


def cell_correct_fraction(inf_rows: Sequence[str], true_rows: Sequence[str],
                          grid: Optional[np.ndarray] = None) -> np.ndarray:
    """Per residue, the fraction of its partners in that column that are homologous.

    The per-cell counterpart of :func:`col_correct_fraction`, and exactly
    consistent with it: a column's score is the mean of its cells' scores. One
    misplaced residue in an otherwise sound
    column scores 0 while its neighbours stay near 1, which is what makes the
    culprit visible rather than merely the damage.

    ``nan`` for a gap, and for a residue with no partners at all — a residue
    alone in its column asserts nothing, so it can be neither right nor wrong.
    """
    grid = true_column_grid(inf_rows, true_rows) if grid is None else grid
    out = np.full(grid.shape, np.nan, dtype=float)
    for c in range(grid.shape[1]):
        col = grid[:, c]
        present = np.flatnonzero(col >= 0)
        if present.size < 2:
            continue
        counts = Counter(int(v) for v in col[present])
        for s in present:
            out[s, c] = (counts[int(col[s])] - 1) / (present.size - 1)
    return out


@dataclass
class Placement:
    """One residue, where it sits, and where the reference says it belongs."""

    seq: int                       #: row index in the inferred alignment
    seq_id: str                    #: that row's identifier
    residue: int                   #: index in the ungapped sequence, 0-based
    char: str                      #: the residue itself
    true_col: int                  #: the column the reference puts it in
    col: int                       #: the column it currently occupies
    target_col: Optional[int] = None
    """The inferred column holding most of its true partners — where it would
    have to move to join them. ``None`` if it has no true partners to join."""

    @property
    def offset(self) -> Optional[int]:
        """Columns to the right (negative: to the left) it would have to move."""
        return None if self.target_col is None else self.target_col - self.col


@dataclass
class ColumnTruth:
    """What the reference says about one inferred column.

    ``agree`` are the residues the reference groups together and the alignment
    got right; ``intruders`` are residues placed here that belong elsewhere; and
    ``missing`` are residues that belong here but were placed elsewhere. A
    correct column has neither of the latter two.
    """

    col: int
    true_col: Optional[int]        #: the reference column this one mostly represents
    agree: List[Placement] = field(default_factory=list)
    intruders: List[Placement] = field(default_factory=list)
    missing: List[Placement] = field(default_factory=list)

    @property
    def correct(self) -> bool:
        return not self.intruders and not self.missing

    def describe(self) -> str:
        """The same thing in words, for the command line and for tooltips."""
        if self.true_col is None:
            return f"Column {self.col + 1} is empty."
        def names(ps):
            return ", ".join(f"{p.seq_id}/{p.char}{p.residue + 1}" for p in ps)
        lines = [f"Column {self.col + 1} — reference column {self.true_col + 1}"]
        if self.agree:
            lines.append(f"Correctly grouped: {names(self.agree)}")
        if self.intruders:
            for p in self.intruders:
                where = ("no true partners in this alignment" if p.offset is None
                         else f"belongs {abs(p.offset)} column"
                              f"{'s' if abs(p.offset) != 1 else ''} to the "
                              f"{'right' if p.offset > 0 else 'left'}, "
                              f"in column {p.target_col + 1}")
                lines.append(f"Does not belong here: {p.seq_id}/{p.char}"
                             f"{p.residue + 1} — {where}")
        if self.missing:
            for p in self.missing:
                lines.append(f"Should be here but is not: {p.seq_id}/{p.char}"
                             f"{p.residue + 1} — currently in column {p.col + 1}")
        if self.correct:
            lines.append("This column is exactly right.")
        return "\n".join(lines)


def explain_column(inf_rows: Sequence[str], true_rows: Sequence[str], col: int,
                   ids: Optional[Sequence[str]] = None,
                   grid: Optional[np.ndarray] = None) -> ColumnTruth:
    """Account for one inferred column against the reference.

    Pass ``grid`` from :func:`true_column_grid` to explain many columns without
    rebuilding the lookup each time.
    """
    grid = true_column_grid(inf_rows, true_rows) if grid is None else grid
    n = len(inf_rows)
    ids = list(ids) if ids is not None else [f"seq_{i + 1}" for i in range(n)]
    ridx = [residue_index(r) for r in inf_rows]

    def placement(s: int, c: int, true_col: int, target: Optional[int] = None) -> Placement:
        return Placement(seq=s, seq_id=ids[s], residue=ridx[s][c], char=inf_rows[s][c],
                         true_col=true_col, col=c, target_col=target)

    here = grid[:, col]
    present = [int(s) for s in np.flatnonzero(here >= 0)]
    if not present:
        return ColumnTruth(col=col, true_col=None)

    # The reference column this one mostly represents. Ties go to the leftmost,
    # so the answer never depends on dictionary order.
    counts = Counter(int(here[s]) for s in present)
    best = max(counts.values())
    true_col = min(t for t, k in counts.items() if k == best)

    def column_of(s: int, t: int) -> Optional[int]:
        """Which inferred column holds sequence ``s``'s residue from true column ``t``."""
        found = np.flatnonzero(grid[s] == t)
        return int(found[0]) if found.size else None

    def target_for(s: int, t: int) -> Optional[int]:
        """Where this residue's true partners currently sit, by majority."""
        elsewhere = Counter()
        for other in range(n):
            if other == s:
                continue
            c = column_of(other, t)
            if c is not None:
                elsewhere[c] += 1
        if not elsewhere:
            return None
        most = max(elsewhere.values())
        return min(c for c, k in elsewhere.items() if k == most)

    out = ColumnTruth(col=col, true_col=true_col)
    for s in present:
        t = int(here[s])
        if t == true_col:
            out.agree.append(placement(s, col, t))
        else:
            out.intruders.append(placement(s, col, t, target_for(s, t)))

    for s in range(n):
        if here[s] == true_col:
            continue
        c = column_of(s, true_col)
        if c is not None:
            out.missing.append(placement(s, c, true_col))
    return out
