"""Benchmark-only metrics.

The accuracy primitives (SP, TC, per-column correctness, AUC) are in
``craic.evaluate`` and re-exported here; what remains in this module is the part
only the benchmark needs:

* Masking gain against a random-mask null and rule-based controls.
* Tree error: Robinson-Foulds distance between the NJ tree of an alignment and a
  reference topology.
* The masking-method comparison (CRAIC reliability vs trimAl / Gblocks / gap).
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

# The accuracy primitives live in the installed package, not here, so that the
# GUI's truth mode, the command line and this harness all score an alignment the
# same way. Re-exported for the benchmark scripts that import them from here.
from craic.evaluate import (  # noqa: F401
    _columns, _midranks, _pairs, auc, col_correct_fraction, sp_score, tc_score,
)


# --------------------------------------------------------------------------- #
# Reliability vs truth
# --------------------------------------------------------------------------- #

def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def masking_gain(
    rel: np.ndarray,
    correct_frac: np.ndarray,
    threshold: float = 0.5,
    controls: Optional[Dict[str, np.ndarray]] = None,
    n_null: int = 200,
    seed: int = 0,
) -> dict:
    """Does masking on ``rel`` raise mean column accuracy *more than chance*?

    The naive comparison — mean accuracy over all columns versus over the
    retained ones — is not a result. Any score correlated with column difficulty
    raises retained accuracy simply by removing hard columns, so
    ``acc_retained > acc_all`` holds for column gap fraction, for column entropy,
    and indeed for any score with a trace of signal. Reporting it as evidence
    that reliability masking works is circular.

    So the retained accuracy is compared against two controls, each removing
    *exactly the same number of columns*:

    * ``acc_retained_random`` — the mean over ``n_null`` random masks. This is
      the floor: what removing that many arbitrary columns buys.
    * ``acc_retained_<name>`` for each score in ``controls`` — e.g. gap fraction.
      This is the bar that matters: a model-based reliability score has to beat
      the trivial rule-based ones, not just beat nothing.

    ``z`` places the observed retained accuracy on the random-mask null
    distribution. It is the number to report.
    """
    rel = np.asarray(rel, dtype=float)
    correct_frac = np.asarray(correct_frac, dtype=float)
    ok = ~np.isnan(rel) & ~np.isnan(correct_frac)
    rel_o, frac_o = rel[ok], correct_frac[ok]
    n = frac_o.size
    retained = rel_o >= threshold
    n_masked = int((~retained).sum())

    out = {
        "acc_all": float(frac_o.mean()) if n else float("nan"),
        "acc_retained": float(frac_o[retained].mean()) if retained.any() else float("nan"),
        "n_masked": n_masked,
        "n_cols": int(n),
        "acc_retained_random": float("nan"),
        "acc_retained_random_sd": float("nan"),
        "z": float("nan"),
    }
    if n == 0 or n_masked == 0 or n_masked >= n:
        # Nothing was masked, or everything was: no contrast to control for.
        return out

    rng = np.random.default_rng(seed)
    null = np.empty(n_null)
    for k in range(n_null):
        keep = np.ones(n, dtype=bool)
        keep[rng.choice(n, n_masked, replace=False)] = False
        null[k] = frac_o[keep].mean()
    sd = float(null.std(ddof=1))
    out["acc_retained_random"] = float(null.mean())
    out["acc_retained_random_sd"] = sd
    out["z"] = float((out["acc_retained"] - null.mean()) / sd) if sd > 0 else float("nan")

    for name, score in (controls or {}).items():
        s = np.asarray(score, dtype=float)[ok]
        # Remove the same number of columns, taking the lowest-scoring ones.
        cut = np.argsort(np.nan_to_num(s, nan=-np.inf), kind="mergesort")[n_masked:]
        keep = np.zeros(n, dtype=bool)
        keep[cut] = True
        out[f"acc_retained_{name}"] = float(frac_o[keep].mean()) if keep.any() else float("nan")
    return out


# --------------------------------------------------------------------------- #
# Tree error (NJ + Robinson-Foulds)
# --------------------------------------------------------------------------- #

def _p_distance_matrix(names: List[str], rows: List[str]):
    from Bio.Phylo.TreeConstruction import DistanceMatrix
    n = len(rows)
    lower = []
    for i in range(n):
        row = []
        for j in range(i):
            m = [(a, b) for a, b in zip(rows[i], rows[j]) if a != "-" and b != "-"]
            d = (sum(a != b for a, b in m) / len(m)) if m else 1.0
            row.append(d)
        row.append(0.0)
        lower.append(row)
    return DistanceMatrix(names=list(names), matrix=lower)


def nj_tree_bipartitions(names: List[str], rows: List[str]) -> set:
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    dm = _p_distance_matrix(names, rows)
    tree = DistanceTreeConstructor().nj(dm)
    all_leaves = frozenset(names)
    ref = min(all_leaves)
    bips = set()
    for clade in tree.get_nonterminals():
        side = frozenset(t.name for t in clade.get_terminals())
        if 2 <= len(side) <= len(all_leaves) - 2:
            canon = side if ref not in side else frozenset(all_leaves - side)
            bips.add(canon)
    return bips


def true_tree_bipartitions(tree: dict, names: List[str]) -> set:
    all_leaves = frozenset(names)
    ref = min(all_leaves)
    bips = set()

    def descend(node) -> frozenset:
        if not node["children"]:
            return frozenset([node["name"]])
        leaves: set = set()
        for child, _ in node["children"]:
            leaves |= descend(child)
        side = frozenset(leaves)
        if 2 <= len(side) <= len(all_leaves) - 2:
            canon = side if ref not in side else frozenset(all_leaves - side)
            bips.add(canon)
        return side

    descend(tree)
    return bips


def rf_distance(true_bips: set, inf_bips: set, normalized: bool = True) -> float:
    rf = len(true_bips ^ inf_bips)
    if not normalized:
        return float(rf)
    denom = len(true_bips) + len(inf_bips)
    return float(rf / denom) if denom else 0.0


# --------------------------------------------------------------------------- #
# Masking-method comparison: CRAIC reliability vs trimAl / Gblocks / gap
# --------------------------------------------------------------------------- #

def _mask_stats(keep, frac, correct, a, n_null: int = 200, seed: int = 0):
    """Score one keep-mask against the truth, with a matched random baseline.

    ``acc_kept_random`` is the mean accuracy of a random mask retaining the same
    number of columns. Without it, ``acc_kept`` alone says nothing: removing any
    columns at all raises it.
    """
    scored = ~np.isnan(frac)
    ks = keep & scored
    corr = scored & correct
    n, n_keep = int(scored.sum()), int(ks.sum())
    stats = {
        "kept": float(n_keep / n) if n else float("nan"),
        "acc_all": float(np.nanmean(frac[scored])) if n else float("nan"),
        "acc_kept": float(np.nanmean(frac[ks])) if ks.any() else float("nan"),
        "precision": float(correct[ks].mean()) if ks.any() else float("nan"),   # kept cols that are correct
        "recall": float(keep[corr].mean()) if corr.any() else float("nan"),      # correct cols that are kept
        "auc": a,
        "acc_kept_random": float("nan"),
    }
    if 0 < n_keep < n:
        rng = np.random.default_rng(seed)
        frac_o = frac[scored]
        null = np.empty(n_null)
        for k in range(n_null):
            pick = rng.choice(n, n_keep, replace=False)
            null[k] = frac_o[pick].mean()
        stats["acc_kept_random"] = float(null.mean())
    return stats


def compare_masking(aln, true_rows, do_perturbation: bool = False):
    """For one inferred alignment (with known-true reference), score every masking
    method by how well it removes the truly mis-aligned columns.

    Two things this deliberately does not do. It does not apply one threshold to
    several unrelated scales — a 0.5 cut means something different for a
    posterior, a gap fraction and a BLOSUM-normalised similarity, so the
    resulting retention figures were never comparable. And it does not leave the
    rule-based methods without an AUC: a binary mask is a (coarse) classifier
    and has a perfectly well-defined AUC, so all six methods are now scored on
    the same threshold-free footing and the comparison can actually be lost.
    """
    from craic.ambiguity import reliability, trimming
    frac = col_correct_fraction(list(aln.rows), true_rows)
    correct = frac >= 0.999
    rep = reliability.analyse(aln, do_perturbation=do_perturbation)
    rel = rep.col_combined if do_perturbation else rep.col_consistency

    def scored_auc(score):
        ok = (~np.isnan(frac)) & (~np.isnan(np.asarray(score, dtype=float)))
        return auc(np.asarray(score, float)[ok], correct[ok]) if ok.sum() > 3 else float("nan")

    out = {}
    for name, score in [("CRAIC reliability", rel),
                        ("gap (MSA_trimmer)", trimming.gap_score(aln)),
                        ("similarity (trimAl)", trimming.similarity_score(aln))]:
        # Each continuous score is cut at its own median, so every method here
        # retains a comparable share of columns; the AUC is threshold-free and
        # is the number to compare.
        s = np.asarray(score, dtype=float)
        thr = float(np.nanmedian(s))
        out[name] = _mask_stats(np.nan_to_num(s, nan=-np.inf) >= thr,
                                frac, correct, scored_auc(s))
    for name, keep in [("trimAl gappyout", trimming.gappyout_mask(aln)),
                       ("trimAl strict", trimming.trimal_strict_mask(aln)),
                       ("Gblocks", trimming.gblocks_mask(aln))]:
        keep = np.asarray(keep, dtype=bool)
        out[name] = _mask_stats(keep, frac, correct, scored_auc(keep.astype(float)))
    return out
