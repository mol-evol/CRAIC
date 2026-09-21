"""Regression tests for the v0.2 correctness fixes.

Each test here corresponds to a specific way v0.1 could produce a plausible but
wrong answer. They are kept together, and named for the failure rather than the
function, so that the reason each one exists survives.
"""

import warnings

import numpy as np
import pytest

from craic import accel, progressive
from craic.ambiguity import disagreement, posterior, reliability, sandbox, trimming
from craic.domain import Alignment, Alphabet, CodingSpec


def _aln(rows, ids=None, alphabet=Alphabet.DNA, **kw):
    ids = ids or [f"s{i}" for i in range(len(rows))]
    return Alignment(list(ids), list(rows), alphabet, **kw)


# --------------------------------------------------------------------------- #
# Gap characters
# --------------------------------------------------------------------------- #

def test_dot_and_question_gaps_are_normalised_on_construction():
    """v0.1 counted '.' and '?' as residues in every residue-index walk, which
    shifts every index downstream and corrupts every score derived from it."""
    aln = _aln(["AC.GT", "AC?GT"])
    assert aln.rows == ["AC-GT", "AC-GT"]
    assert aln.gap_fraction(2) == 1.0


def test_normalised_gaps_give_the_same_reliability_as_dashes():
    a = reliability.consistency(_aln(["ACGTACGT", "ACG.ACGT", "ACGTACGT"]))[0]
    b = reliability.consistency(_aln(["ACGTACGT", "ACG-ACGT", "ACGTACGT"]))[0]
    np.testing.assert_allclose(a, b, equal_nan=True)


# --------------------------------------------------------------------------- #
# Pair-HMM parameter validation
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("delta", [0.5, 0.7, 0.0, -0.1])
def test_invalid_delta_is_rejected_rather_than_producing_nan(delta):
    """log(1 - 2*delta) is undefined at delta >= 0.5. The Rust core turned that
    into NaN that propagated silently through every posterior."""
    model = accel.emission_model("dna")
    with pytest.raises(ValueError):
        accel.posterior_matrix("ACGT", "ACGT", model, delta, 0.5)


@pytest.mark.parametrize("epsilon", [1.0, 1.5, -0.01])
def test_invalid_epsilon_is_rejected(epsilon):
    model = accel.emission_model("dna")
    with pytest.raises(ValueError):
        accel.posterior_matrix("ACGT", "ACGT", model, 0.02, epsilon)


def test_valid_parameters_still_produce_finite_posteriors():
    model = accel.emission_model("dna")
    P = accel.posterior_matrix("ACGTACGT", "ACGTACGT", model, 0.02, 0.5)
    assert np.isfinite(P).all()
    assert (P >= 0).all() and (P <= 1).all()


# --------------------------------------------------------------------------- #
# Reliability: estimated parameters, loud failure, one mask
# --------------------------------------------------------------------------- #

def test_consistency_uses_parameters_estimated_from_the_alignment():
    """v0.1 scored every alignment under a fixed 90%-identity model, which makes
    divergent families look unreliable and biases masking against them."""
    rows = ["ACGTACGTAA", "TGCATGCATT", "ACGAACGTAG"]
    aln = _aln(rows)
    est_model, est_delta, est_eps = progressive.estimate_params(rows, Alphabet.DNA)
    from_default = reliability.consistency(
        aln, model=accel.emission_model("dna"), delta=0.02, epsilon=0.5)[0]
    from_estimate = reliability.consistency(aln)[0]
    estimated_explicitly = reliability.consistency(
        aln, model=est_model, delta=est_delta, epsilon=est_eps)[0]
    # the default is what it now agrees with
    np.testing.assert_allclose(from_estimate, estimated_explicitly, equal_nan=True)
    # and that is a different answer from the old fixed model
    assert not np.allclose(np.nan_to_num(from_estimate),
                           np.nan_to_num(from_default))


def test_total_ensemble_failure_raises_instead_of_masking_everything(monkeypatch):
    """The single most dangerous line in v0.1: a bare except in the replicate
    loop meant a fully failed ensemble returned all-nan, which keep_mask read as
    'every column is maximally unreliable' and trimming obeyed."""
    aln = _aln(["ACGTACGT", "ACGTACGT", "ACGAACGT"])

    def boom(*a, **kw):
        raise RuntimeError("engine exploded")

    monkeypatch.setattr(progressive, "align", boom)
    with pytest.raises(reliability.EnsembleFailure):
        reliability.perturbation(aln, n_replicates=4)


def test_partial_ensemble_failure_warns_and_reports_the_count():
    aln = _aln(["ACGTACGT", "ACGTACGT", "ACGAACGT"])
    real = progressive.align
    calls = {"n": 0}

    def flaky(*a, **kw):
        calls["n"] += 1
        if calls["n"] % 2 == 0:
            raise RuntimeError("intermittent")
        return real(*a, **kw)

    stats = {}
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(progressive, "align", flaky)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            reliability.perturbation(aln, n_replicates=4, stats=stats)
    assert stats["n_failed"] == 2 and stats["n_ok"] == 2
    assert any("replicates failed" in str(w.message) for w in caught)


def test_keep_mask_refuses_an_all_unscored_score():
    with pytest.raises(ValueError):
        reliability.keep_mask(np.full(5, np.nan), 0.5)


def test_keep_mask_nan_policy_is_explicit_and_shared():
    """v0.1 had two definitions with opposite NaN handling — the GUI kept
    unscored columns, the benchmark dropped them — so the exported mask was not
    the evaluated mask."""
    scores = np.array([0.9, np.nan, 0.1])
    np.testing.assert_array_equal(
        reliability.keep_mask(scores, 0.5), [True, False, False])
    np.testing.assert_array_equal(
        reliability.keep_mask(scores, 0.5, unscored="keep"), [True, True, False])
    with pytest.raises(ValueError):
        reliability.keep_mask(scores, 0.5, unscored="sometimes")


def test_analyse_records_the_replicate_counts():
    aln = _aln(["ACGTACGTAC", "ACGTACGTAC", "ACGAACGTAC"])
    rep = reliability.analyse(aln, do_perturbation=True, n_replicates=3)
    assert rep.n_replicates_ok == 3
    assert rep.n_replicates_failed == 0


# --------------------------------------------------------------------------- #
# Trimming
# --------------------------------------------------------------------------- #

def test_zero_mad_flags_no_outliers():
    """With identical sequences the MAD is 0; v0.1 substituted 1e-9 and then
    flagged every sequence scoring below the median."""
    aln = _aln(["ACGTACGT"] * 5)
    scores, flags = trimming.outlier_sequences(aln)
    assert not flags.any()


def test_gappyout_does_not_trim_almost_everything_on_a_top_heavy_jump():
    """g = [0.99]*100 + [1.0] put v0.1's cut at 0.995 and removed 100 of 101
    columns: the largest jump sat at the top of the distribution."""
    n_cols = 100
    # 100 columns with one gap out of 100 sequences, then one with none
    rows = []
    for i in range(100):
        row = ["A"] * n_cols + ["A"]
        if i < 1:
            for c in range(n_cols):
                row[c] = "-"
        rows.append("".join(row))
    keep = trimming.gappyout_mask(_aln(rows))
    assert keep.sum() >= 0.5 * len(keep)


def test_gappyout_still_separates_a_genuinely_bimodal_distribution():
    clean = ["ACGTACGT" + "A" * 8] * 6
    gappy = [r[:8] + "-" * 8 for r in clean]
    rows = clean[:3] + gappy[:3]
    keep = trimming.gappyout_mask(_aln(rows))
    assert keep[:8].all()          # the clean half survives


# --------------------------------------------------------------------------- #
# Disagreement map
# --------------------------------------------------------------------------- #

def test_truncated_ids_do_not_read_as_total_disagreement():
    """MAFFT truncates FASTA ids at whitespace. In v0.1 a miss incremented the
    pair denominator but not the numerator, so every column scored 0 agreement —
    a wrong result indistinguishable from a real one."""
    ref = _aln(["ACGT", "ACGT"], ids=["seq1 human haemoglobin", "seq2 mouse"])
    other = _aln(["ACGT", "ACGT"], ids=["seq1", "seq2"])
    agree = disagreement.compare(ref, other)
    assert np.nanmin(agree) == 1.0


def test_genuinely_different_sequences_still_raise():
    ref = _aln(["ACGT", "ACGT"], ids=["a", "b"])
    other = _aln(["ACGT", "ACGT"], ids=["x", "y"])
    with pytest.raises(disagreement.IdMismatch):
        disagreement.compare(ref, other)


# --------------------------------------------------------------------------- #
# Edge cases
# --------------------------------------------------------------------------- #

def test_empty_input_to_align_is_an_error_not_a_valueerror_from_max():
    with pytest.raises(ValueError, match="no sequences"):
        progressive.align([], Alphabet.DNA)


def test_duplicate_ids_are_rejected():
    with pytest.raises(ValueError, match="duplicate sequence id"):
        progressive.align([("a", "ACGT"), ("a", "ACGT")], Alphabet.DNA)


def test_probe_residue_out_of_range_raises_indexerror():
    aln = _aln(["ACGT", "ACGT"])
    with pytest.raises(IndexError):
        posterior.probe_residue(aln, 0, 99)


def test_splice_drops_the_coding_annotation_when_it_would_shift_the_frame():
    """A block one column wider re-phases every codon downstream while
    CodingSpec.frame stays put: a wrong amino-acid view of correct nucleotides."""
    aln = _aln(["ATGAAACCC", "ATGAAACCC"], coding=CodingSpec(frame=0))
    in_frame = _aln(["GGG", "GGG"], coding=CodingSpec(frame=0))
    out_of_frame = _aln(["GGGG", "GGGG"], coding=CodingSpec(frame=0))
    assert sandbox.splice(aln, 3, 6, in_frame).coding is not None
    assert sandbox.splice(aln, 3, 6, out_of_frame).coding is None


def test_ambiguity_profile_returns_one_entry_per_residue():
    aln = _aln(["ACGT", "----"])
    prof = posterior.ambiguity_profile(aln, 0, 1)
    assert prof.shape == (4,)


# --------------------------------------------------------------------------- #
# GUI: the posterior explorer must not offer a sequence against itself
# --------------------------------------------------------------------------- #

def test_sequence_pickers_exclude_each_other():
    """A sequence versus itself is the trivial diagonal. v0.2 withdraws the
    choice instead of accepting it and then refusing to plot it."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication, QComboBox
    from craic.gui.panels import exclusive_pair

    app = QApplication.instance() or QApplication([])          # noqa: F841
    a, b = QComboBox(), QComboBox()
    ids = ["s0", "s1", "s2"]
    a.addItems(ids)
    b.addItems(ids)
    b.setCurrentIndex(1)
    exclusive_pair(a, b)

    # each combo's entry for the other's selection is not selectable
    assert not b.model().item(a.currentIndex()).isEnabled()
    assert not a.model().item(b.currentIndex()).isEnabled()
    # …and every other entry still is
    assert a.model().item(0).isEnabled()

    # a collision moves the second picker rather than being rejected later
    a.setCurrentIndex(1)
    exclusive_pair(a, b)
    assert a.currentIndex() != b.currentIndex()


def test_exclusive_pair_is_a_noop_for_a_single_sequence():
    """With one sequence there is no valid pair; disabling the only entry would
    leave an unusable control, so the pair logic stands aside."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication, QComboBox
    from craic.gui.panels import exclusive_pair

    app = QApplication.instance() or QApplication([])          # noqa: F841
    a, b = QComboBox(), QComboBox()
    a.addItems(["only"])
    b.addItems(["only"])
    exclusive_pair(a, b)
    assert a.model().item(0).isEnabled()
    assert b.model().item(0).isEnabled()


# --------------------------------------------------------------------------- #
# GUI: controls that would do nothing are not offered
# --------------------------------------------------------------------------- #

def test_colour_schemes_are_amino_acid_only():
    """All four schemes colour amino acids; residue_color ignores the scheme for
    nucleotides, so the control must report itself inapplicable there."""
    from craic.domain import Alphabet as A, Level as L
    from craic.gui import colors

    assert not colors.scheme_applies(L.NT, A.DNA)
    assert not colors.scheme_applies(L.NT, A.RNA)
    assert colors.scheme_applies(L.CODON, A.DNA)     # codons colour by translation
    assert colors.scheme_applies(L.AA, A.DNA)
    assert colors.scheme_applies(L.AA, A.PROTEIN)


def test_scheme_applies_matches_what_residue_color_actually_does():
    """Guard against the two drifting apart: where scheme_applies says False,
    changing the scheme must not change any colour."""
    from craic.domain import Alphabet as A, Level as L
    from craic.gui import colors

    for level, alphabet, cells in ((L.NT, A.DNA, "ACGT"), (L.AA, A.PROTEIN, "AKLW")):
        differs = any(
            colors.residue_color(c, level, alphabet, scheme="clustal")
            != colors.residue_color(c, level, alphabet, scheme="taylor")
            for c in cells
        )
        assert differs == colors.scheme_applies(level, alphabet)


def test_score_legend_paints_the_same_ramp_the_canvas_uses():
    """The key is only useful if it is the same ramp. Sample the rendered bar and
    check its ends and middle against colors.score_color."""
    pytest.importorskip("PySide6")
    from PySide6.QtCore import QPoint
    from PySide6.QtGui import QImage
    from PySide6.QtWidgets import QApplication

    from craic.gui import colors
    from craic.gui.panels import ScoreLegend

    app = QApplication.instance() or QApplication([])          # noqa: F841
    w = ScoreLegend()
    img = QImage(w.width(), w.height(), QImage.Format_ARGB32)
    img.fill(0xFF000000)
    w.render(img, QPoint(0, 0))

    y = w.height() // 2
    x0 = ScoreLegend._PAD
    for frac in (0.0, 0.5, 1.0):
        x = x0 + round(frac * (w._bar_w - 1))
        got = img.pixelColor(x, y)
        want = colors.score_color(frac)
        # within one unit per channel: the bar is drawn column by column
        assert abs(got.red() - want.red()) <= 1
        assert abs(got.green() - want.green()) <= 1
        assert abs(got.blue() - want.blue()) <= 1
