"""Tests for the benchmark layer.

This is the code that produces every number in the paper, and in v0.1 it had no
tests at all. The metrics are small pure functions over known inputs, so they are
cheap to pin down exactly — which matters more here than anywhere else in the
codebase, because a silently wrong metric does not crash, it publishes.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "benchmarks"))

import metrics as M          # noqa: E402
from craic import simulate as S              # noqa: E402

from craic.domain import Alphabet     # noqa: E402


# --------------------------------------------------------------------------- #
# SP / TC against hand-computed answers
# --------------------------------------------------------------------------- #

def test_sp_and_tc_are_perfect_on_an_identical_alignment():
    rows = ["ACGT", "A-GT", "ACGT"]
    assert M.sp_score(rows, rows) == (1.0, 1.0)
    assert M.tc_score(rows, rows) == 1.0


def test_sp_counts_the_right_pairs_by_hand():
    #        col0 col1
    true = ["AC", "AC"]          # pairs: (s0r0,s1r0), (s0r1,s1r1) -> 2
    inf = ["AC-", "-AC"]         # s0 residues at cols 0,1; s1 at cols 1,2
    #                              shared column 1 pairs s0r1 with s1r0 -> not true
    recall, prec = M.sp_score(true, inf)
    assert recall == 0.0
    assert prec == 0.0


def test_tc_requires_the_whole_column():
    true = ["AC", "AC", "AC"]
    # third sequence misplaced: no reference column is reproduced exactly
    inf = ["AC-", "AC-", "-AC"]
    assert M.tc_score(true, inf) == 0.0


def test_core_mask_restricts_sp_and_tc_to_core_columns():
    true = ["AC", "AD"]
    inf = ["AC", "AD"]
    core = np.array([True, False])
    # with the second column excluded there is one reference pair, still matched
    assert M.sp_score(true, inf, core)[0] == 1.0
    # and a mismatch outside the core no longer counts against the score
    inf_bad = ["AC-", "A-D"]
    assert M.sp_score(true, inf_bad, core)[0] == 1.0
    assert M.sp_score(true, inf_bad)[0] == 0.5


def test_col_correct_fraction_marks_the_wrong_column():
    true = ["AC", "AC"]
    inf = ["AC", "AC"]
    frac = M.col_correct_fraction(inf, true)
    assert list(frac) == [1.0, 1.0]
    shifted = M.col_correct_fraction(["AC-", "-AC"], true)
    assert np.isnan(shifted[0]) and np.isnan(shifted[2])   # singleton columns
    assert shifted[1] == 0.0                                # wrong pairing


# --------------------------------------------------------------------------- #
# AUC — including the tie handling that reliability scores actually exercise
# --------------------------------------------------------------------------- #

def test_auc_perfect_and_inverted():
    score = np.array([0.1, 0.2, 0.8, 0.9])
    pos = np.array([False, False, True, True])
    assert M.auc(score, pos) == 1.0
    assert M.auc(-score, pos) == 0.0


def test_auc_of_a_constant_score_is_one_half_not_order_dependent():
    """The case that matters: reliability ties heavily at 1.0 in conserved
    columns. Without midranks this returns an arbitrary number that depends on
    the order of the input."""
    score = np.ones(6)
    pos = np.array([True, True, True, False, False, False])
    assert M.auc(score, pos) == 0.5
    assert M.auc(score, pos[::-1]) == 0.5


def test_auc_with_partial_ties_matches_the_hand_computed_value():
    # scores 1,1,0 for positives and 1,0 for negatives.
    # Mann-Whitney U = sum over pairs of [pos>neg] + 0.5*[pos==neg]
    #   (1,1): 0.5   (1,0): 1     (1,1): 0.5   (1,0): 1     (0,1): 0   (0,0): 0.5
    #   total = 3.5 out of 6 pairs
    score = np.array([1.0, 1.0, 0.0, 1.0, 0.0])
    pos = np.array([True, True, True, False, False])
    assert M.auc(score, pos) == pytest.approx(3.5 / 6.0)


def test_auc_is_nan_when_one_class_is_empty():
    assert np.isnan(M.auc(np.array([1.0, 2.0]), np.array([True, True])))


# --------------------------------------------------------------------------- #
# The masking control — the point of the v0.2 rewrite
# --------------------------------------------------------------------------- #

def test_masking_gain_reports_a_random_null_at_the_same_column_count():
    rng = np.random.default_rng(0)
    frac = rng.random(200)
    rel = frac.copy()                    # a perfect predictor
    out = M.masking_gain(rel, frac, threshold=0.5, n_null=50, seed=1)
    assert out["n_masked"] > 0
    # a perfect predictor must beat the random mask that removes as many columns
    assert out["acc_retained"] > out["acc_retained_random"]
    assert out["z"] > 3


def test_a_useless_score_does_not_beat_the_random_null():
    """The regression guard for the tautology. A reliability score with no
    information still raises acc_retained above acc_all, because it removes
    columns; it must not beat the matched random mask."""
    rng = np.random.default_rng(3)
    frac = rng.random(300)
    rel = rng.random(300)                # independent of the truth
    out = M.masking_gain(rel, frac, threshold=0.5, n_null=200, seed=2)
    assert out["acc_retained"] > out["acc_all"] - 0.05   # the misleading comparison
    assert abs(out["z"]) < 3                             # the honest one


def test_masking_gain_handles_nothing_masked():
    frac = np.full(10, 0.9)
    out = M.masking_gain(np.ones(10), frac, threshold=0.5)
    assert out["n_masked"] == 0
    assert np.isnan(out["z"])


def test_masking_gain_control_uses_the_same_number_of_columns():
    frac = np.linspace(0, 1, 100)
    rel = frac.copy()
    gap = np.random.default_rng(0).random(100)
    out = M.masking_gain(rel, frac, 0.5, controls={"gap": gap}, n_null=20)
    assert "acc_retained_gap" in out
    assert not np.isnan(out["acc_retained_gap"])


# --------------------------------------------------------------------------- #
# Tree metrics
# --------------------------------------------------------------------------- #

def test_rf_distance_of_identical_bipartition_sets_is_zero():
    bips = {frozenset({"a", "b"})}
    assert M.rf_distance(bips, bips) == 0.0
    assert M.rf_distance(bips, set()) == 1.0


def test_nj_tree_recovers_an_obvious_grouping():
    names = ["a", "b", "c", "d"]
    rows = ["AAAACCCC", "AAAACCCC", "GGGGTTTT", "GGGGTTTT"]
    bips = M.nj_tree_bipartitions(names, rows)
    assert frozenset({"c", "d"}) in bips or frozenset({"a", "b"}) in bips


# --------------------------------------------------------------------------- #
# The simulator's own ground truth must be internally consistent
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("kwargs", [
    {},
    {"rate_alpha": 0.5},
    {"indel_zipf": 1.7},
    {"rate_alpha": 0.8, "indel_zipf": 2.0},
])
def test_simulated_truth_is_self_consistent(kwargs):
    d = S.simulate(taxa=6, root_len=80, seed=4, **kwargs)
    rows, seqs = d["true_rows"], d["seqs"]
    assert len({len(r) for r in rows}) == 1              # rectangular
    assert len(rows) == len(seqs) == len(d["names"])
    for (name, ungapped), row in zip(seqs, rows):
        # the stated ungapped sequence is exactly the row without gaps
        assert row.replace("-", "") == ungapped
    # every column holds at least one residue
    for c in range(len(rows[0])):
        assert any(r[c] != "-" for r in rows)


def test_simulation_is_reproducible_from_the_seed():
    a = S.simulate(taxa=5, root_len=60, seed=11)
    b = S.simulate(taxa=5, root_len=60, seed=11)
    assert a["true_rows"] == b["true_rows"]


def test_rate_variation_changes_the_outcome():
    plain = S.simulate(taxa=5, root_len=120, seed=7)
    varied = S.simulate(taxa=5, root_len=120, seed=7, rate_alpha=0.3)
    assert plain["true_rows"] != varied["true_rows"]


# --------------------------------------------------------------------------- #
# End to end: the masking comparison runs and scores every method
# --------------------------------------------------------------------------- #

def test_compare_masking_scores_every_method_including_the_rule_based_ones():
    from craic.engines import BuiltinProgressive
    # A genuinely inferred alignment, so that some columns are wrong and the AUC
    # is defined. (Scoring the true alignment against itself has no negatives.)
    d = S.simulate(taxa=5, root_len=90, seed=2, indel_rate=3.0, bmax=0.4)
    aln = BuiltinProgressive().align(d["seqs"], Alphabet.DNA)
    res = M.compare_masking(aln, d["true_rows"])
    assert set(res) == {"CRAIC reliability", "gap (MSA_trimmer)", "similarity (trimAl)",
                        "trimAl gappyout", "trimAl strict", "Gblocks"}
    for name, stats in res.items():
        # every method gets an AUC now, so the comparison can be lost
        assert not np.isnan(stats["auc"]), f"{name} has no AUC"


# --------------------------------------------------------------------------- #
# Residue masks in the tree comparison
# --------------------------------------------------------------------------- #

def test_matched_and_random_residue_masks_remove_exactly_k_residues():
    from craic.domain import Alignment

    aln = Alignment(["a", "b", "c"], ["ACGT-", "AC-TT", "ACGTT"], Alphabet.DNA)
    cells = np.array([[0.9, 0.1, 0.8, 0.2, np.nan],
                      [0.9, 0.3, np.nan, 0.7, 0.6],
                      [0.9, 0.4, 0.5, 0.6, np.nan]])
    low = M.lowest_residues(aln, cells, 3)
    assert low == {"a": [1, 3], "b": [1]}                # 0.1, 0.2, 0.3; nan never picked
    rnd = M.random_residues(aln, 4, np.random.default_rng(0))
    assert sum(len(v) for v in rnd.values()) == 4


def test_nj_distance_skips_masked_residues_like_gaps():
    names = ["a", "b", "c", "d"]
    rows = ["AAAAAAAA", "AAAAAAAT", "TTTTAAAA", "TTTTAAAT"]
    masked = ["AAAAAAAA", "AAAAAAAN", "TTTTAAAA", "TTTTAAAN"]
    # masking b's and d's differing site as missing makes a~b and c~d identical
    dm = M._p_distance_matrix(names, masked, missing="-N")
    assert dm["a", "b"] == 0.0 and dm["c", "d"] == 0.0
    assert M.nj_tree_bipartitions(names, rows) == M.nj_tree_bipartitions(names, masked, "-N")


def test_resuming_an_old_results_file_rewrites_it_under_the_current_header(tmp_path, monkeypatch):
    import csv
    import run_benchmark as RB

    out = tmp_path / "old.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["mode", "condition", "item", "aligner", "sp"])
        w.writerow(["sim", "c", "rep0", "X", "0.5"])
    from craic import engines

    # main() sets these globals; let monkeypatch put them back for later tests
    for obj, name in ((RB, "_ONLY"), (RB, "_TREE"), (engines.ExternalAligner, "timeout"),
                      (engines.BuiltinProgressive, "benchmark_effort"),
                      (engines.BuiltinProgressive, "consistency_mem_gb")):
        monkeypatch.setattr(obj, name, getattr(obj, name))
    monkeypatch.setattr(sys, "argv", ["run_benchmark.py", "--sim", "--reps", "0", "--engines",
                                      "builtin", "--out", str(out)])
    RB.main()
    rows = list(csv.DictReader(open(out)))
    assert rows[0]["sp"] == "0.5" and "rf_resid_matched" in rows[0]
