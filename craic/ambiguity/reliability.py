"""Per-column / per-residue reliability + live masking.

Two complementary, well-grounded signals:

* **consistency** (TCS / heads-style): for each alignment column, average the
  pair-HMM posterior P(res_i ~ res_j) over the residue pairs that column
  asserts are homologous. High = the unaligned-sequence evidence supports the
  column. Uses the CRAIC core directly.
* **perturbation** (guide-tree / gap-regime sensitivity, in the spirit of
  GUIDANCE but not equivalent to it): re-align the same sequences under an
  ensemble of perturbed guide trees across a few gap regimes, and measure what
  fraction of each residue's asserted homologies survive. See
  :func:`perturbation` for what the ensemble does and does not sample.

Both live in [0, 1]; ``combined`` averages them. A threshold turns either into a
column mask that the viewer previews live before you commit.

Two properties this module is careful about, because getting either wrong
produces a plausible-looking number that is wrong:

* Both scores are computed under parameters **estimated from the alignment
  being scored** (:func:`craic.progressive.estimate_params`), not under a fixed
  default. Scoring a 40%-identity family under a 90%-identity model makes every
  column look unreliable and biases masking toward destroying exactly the
  divergent data that masking decisions are about.
* A failed ensemble is an error, not a score of zero. If no replicate
  completes, :func:`perturbation` raises :class:`EnsembleFailure` rather than
  returning all-nan, which downstream would be indistinguishable from "every
  column is maximally unreliable" and would mask the whole alignment away.
"""

from __future__ import annotations

import warnings
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from .. import accel, domain, progressive
from ..domain import Alignment, Alphabet

ResKey = Tuple[int, int]  # (sequence index, residue index)


class EnsembleFailure(RuntimeError):
    """Every replicate in the perturbation ensemble failed.

    Raised rather than returning an all-nan score, which a threshold would read
    as "mask everything".
    """


def _col_nanmean(mat: np.ndarray) -> np.ndarray:
    """Column-wise nanmean that returns nan (not a warning) for empty columns."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(mat, axis=0)


# --------------------------------------------------------------------------- #
# Consistency (uses the pair-HMM posterior core)
# --------------------------------------------------------------------------- #

def consistency(
    aln: Alignment,
    model: Optional[accel.EmissionModel] = None,
    max_pairs: int = 300,
    seed: int = 0,
    delta: Optional[float] = None,
    epsilon: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (col_score[length], cell_score[n, length] with nan at gaps).

    ``model``, ``delta`` and ``epsilon`` default to the values estimated from
    ``aln`` itself rather than to library defaults; pass them explicitly only
    to score under a model of your own choosing.
    """
    est_model, est_delta, est_epsilon = progressive.estimate_params(aln.rows, aln.alphabet)
    if model is None:
        model = est_model
    if delta is None:
        delta = est_delta
    if epsilon is None:
        epsilon = est_epsilon

    n = aln.n_seqs
    L = aln.length
    maps = [domain.residue_index(r) for r in aln.rows]

    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if len(pairs) > max_pairs:
        rng = np.random.default_rng(seed)
        pairs = [pairs[k] for k in rng.choice(len(pairs), max_pairs, replace=False)]

    # accumulate per-cell sums and counts
    cell_sum = np.zeros((n, L))
    cell_cnt = np.zeros((n, L))
    col_sum = np.zeros(L)
    col_cnt = np.zeros(L)

    for i, j in pairs:
        P = accel.posterior_matrix(aln.rows[i], aln.rows[j], model, delta, epsilon)
        if P.size == 0:
            continue
        mi, mj = maps[i], maps[j]
        for c in range(L):
            ri, rj = mi[c], mj[c]
            if ri < 0 or rj < 0:
                continue
            p = float(P[ri, rj])
            col_sum[c] += p
            col_cnt[c] += 1
            cell_sum[i, c] += p
            cell_cnt[i, c] += 1
            cell_sum[j, c] += p
            cell_cnt[j, c] += 1

    col = np.divide(col_sum, col_cnt, out=np.full(L, np.nan), where=col_cnt > 0)
    cell = np.divide(cell_sum, cell_cnt, out=np.full((n, L), np.nan), where=cell_cnt > 0)
    return col, cell


# --------------------------------------------------------------------------- #
# Perturbation (guide-tree / gap-regime sensitivity)
# --------------------------------------------------------------------------- #

# Pair-HMM gap-open probabilities cycled across the replicates, spanning roughly
# an order of magnitude either side of the usual estimate so that the ensemble
# samples gap-model uncertainty as well as guide-tree uncertainty. These are
# tool defaults, not calibrated quantities: they set the scale of the score, so
# they are exposed as an argument and documented rather than buried.
_PERTURB_DELTAS = (0.01, 0.03, 0.08)


def perturbation(
    aln: Alignment,
    alphabet: Optional[Alphabet] = None,
    n_replicates: int = 16,
    seed: int = 0,
    deltas: Tuple[float, ...] = _PERTURB_DELTAS,
    stats: Optional[Dict[str, int]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (col_score[length], cell_score[n, length] with nan at gaps).

    Re-align the same sequences under ``n_replicates`` perturbed guide trees
    (cycling through the gap regimes in ``deltas``) and, for each residue, score
    the fraction of its reference-asserted homologies that survive across the
    ensemble. A residue with no asserted partners (a singleton column) has no
    residue-pair evidence and is left nan rather than scored a free 1.0.

    **What this is.** A sensitivity analysis of one engine: how stable are the
    homologies this alignment asserts, when the guide tree is perturbed and the
    gap model is varied? Low scores mark residues whose placement is an artefact
    of a particular tree or gap cost.

    **What this is not.** It is not GUIDANCE. GUIDANCE bootstraps alignment
    columns to build perturbed guide trees and re-aligns with the *same* aligner
    that produced the reference, over ~100 replicates. Here the replicates come
    from CRAIC's built-in engine under plain posterior decoding (no consistency
    pass), which is a weaker aligner than the one that may have produced the
    reference and a *different* aligner if the reference came from MAFFT or
    PRANK. In that case the score conflates genuine alignment uncertainty with
    systematic between-method difference. Read it as a stability probe, and read
    the disagreement map for between-method evidence.

    Raises :class:`EnsembleFailure` if no replicate completes. If some fail, a
    warning is issued and the surviving replicates are used; ``analyse`` records
    the count on the report.
    """
    alphabet = alphabet or aln.alphabet
    seqs = [domain.ungap(r) for r in aln.rows]
    records = [(str(i), s) for i, s in enumerate(seqs)]

    model, _est_delta, _est_epsilon = progressive.estimate_params(aln.rows, alphabet)

    wanted = max(1, n_replicates)
    alts: List[Alignment] = []
    first_error: Optional[BaseException] = None
    for i in range(wanted):
        try:
            alts.append(progressive.align(
                records, alphabet, model=model,
                delta=deltas[i % len(deltas)], epsilon=0.5,
                estimate=False, consistency_iters=0, guide_seed=seed + i + 1))
        except Exception as exc:          # noqa: BLE001 - recorded, then re-raised if total
            if first_error is None:
                first_error = exc
    n_failed = wanted - len(alts)
    if not alts:
        raise EnsembleFailure(
            f"all {wanted} perturbation replicates failed; "
            f"first error: {first_error!r}"
        ) from first_error
    if n_failed:
        warnings.warn(
            f"{n_failed} of {wanted} perturbation replicates failed; "
            "scores are based on the remainder",
            RuntimeWarning,
            stacklevel=2,
        )
    if stats is not None:
        stats["n_ok"], stats["n_failed"] = len(alts), n_failed

    n = aln.n_seqs
    L = aln.length

    # reference membership uses *positional* sequence index; alts use str(index)
    def mem_by_index(a: Alignment) -> Dict[ResKey, int]:
        mem: Dict[ResKey, int] = {}
        for row_id, row in zip(a.ids, a.rows):
            si = int(row_id)
            r = 0
            for c, ch in enumerate(row):
                if ch != "-":
                    mem[(si, r)] = c
                    r += 1
        return mem

    ref_mem = {(si, r): c for (si, r), c in domain.membership(aln).items()}
    alt_mems = [mem_by_index(a) for a in alts]

    # partners asserted by the reference, per residue
    partners: Dict[ResKey, List[ResKey]] = defaultdict(list)
    maps = [domain.residue_index(r) for r in aln.rows]
    for c in range(L):
        present = [(si, maps[si][c]) for si in range(n) if maps[si][c] >= 0]
        for a_idx in range(len(present)):
            for b_idx in range(len(present)):
                if a_idx != b_idx:
                    partners[present[a_idx]].append(present[b_idx])

    res_score: Dict[ResKey, float] = {}
    for x, plist in partners.items():
        if not alt_mems or not plist:
            continue  # no ensemble, or a singleton column: no residue-pair evidence -> nan
        tot = ok = 0
        for y in plist:
            for mem in alt_mems:
                tot += 1
                cx, cy = mem.get(x), mem.get(y)
                if cx is not None and cx == cy:
                    ok += 1
        if tot:
            res_score[x] = ok / tot

    cell = np.full((n, L), np.nan)
    for (si, r), score in res_score.items():
        c = ref_mem.get((si, r))
        if c is not None:
            cell[si, c] = score
    col = _col_nanmean(cell)
    return col, cell


# --------------------------------------------------------------------------- #
# Combined report + masking
# --------------------------------------------------------------------------- #

def keep_mask(
    scores: np.ndarray,
    threshold: float,
    unscored: str = "drop",
) -> np.ndarray:
    """Columns to keep: ``scores >= threshold``.

    This is the *one* definition of a reliability mask in CRAIC. It previously
    existed twice with opposite handling of unscoreable (nan) columns — the GUI
    kept them, the benchmark dropped them — so the mask a user exported was not
    the mask the benchmark evaluated.

    A nan column is one with no residue-pair evidence at all (a singleton
    column, or one whose ensemble evidence is missing). ``unscored="drop"`` is
    the default and the conservative choice for downstream phylogenetics: a
    column that could not be assessed is not evidence. ``unscored="keep"``
    retains them.

    Raises ``ValueError`` if *every* column is unscoreable, which means the
    reliability analysis failed rather than that the alignment is worthless.
    """
    if unscored not in ("drop", "keep"):
        raise ValueError(f"unscored must be 'drop' or 'keep', not {unscored!r}")
    scores = np.asarray(scores, dtype=float)
    if scores.size and bool(np.all(np.isnan(scores))):
        raise ValueError(
            "every column is unscored, so no threshold is meaningful; "
            "the reliability analysis did not produce usable scores"
        )
    fill = 0.0 if unscored == "drop" else 1.0
    return np.nan_to_num(scores, nan=fill) >= threshold


@dataclass
class Reliability:
    length: int
    col_consistency: np.ndarray
    col_perturbation: np.ndarray
    col_combined: np.ndarray
    cell_consistency: np.ndarray
    cell_perturbation: np.ndarray
    #: replicates that completed / failed in the perturbation ensemble
    n_replicates_ok: int = 0
    n_replicates_failed: int = 0

    @property
    def cell_combined(self) -> np.ndarray:
        """Per-residue confidence (n x length), nan at gaps — for colouring the
        alignment by reliability."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return np.nanmean(
                np.stack([self.cell_consistency, self.cell_perturbation]), axis=0
            )

    def column_scores(self, which: str = "combined") -> np.ndarray:
        return {
            "combined": self.col_combined,
            "consistency": self.col_consistency,
            "perturbation": self.col_perturbation,
        }[which]

    def keep_mask(self, threshold: float, which: str = "combined",
                  unscored: str = "drop") -> np.ndarray:
        return keep_mask(self.column_scores(which), threshold, unscored)


def analyse(
    aln: Alignment,
    do_perturbation: bool = True,
    model: Optional[accel.EmissionModel] = None,
    n_replicates: int = 16,
) -> Reliability:
    col_c, cell_c = consistency(aln, model=model)
    stats: Dict[str, int] = {}
    if do_perturbation:
        col_p, cell_p = perturbation(aln, n_replicates=n_replicates, stats=stats)
    else:
        col_p = np.full(aln.length, np.nan)
        cell_p = np.full((aln.n_seqs, aln.length), np.nan)
    stack = np.vstack([col_c, col_p])
    combined = _col_nanmean(stack)
    return Reliability(
        aln.length, col_c, col_p, combined, cell_c, cell_p,
        n_replicates_ok=stats.get("n_ok", 0),
        n_replicates_failed=stats.get("n_failed", 0),
    )


def apply_mask(aln: Alignment, keep_mask: np.ndarray) -> Alignment:
    keep = [c for c in range(aln.length) if keep_mask[c]]
    rows = ["".join(r[c] for c in keep) for r in aln.rows]
    return Alignment(list(aln.ids), rows, aln.alphabet, coding=None, meta=dict(aln.meta))
