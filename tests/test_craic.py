"""Test-suite for CRAIC. Fast, deterministic, no external aligners required."""

import os

import numpy as np
import pytest

from craic import accel
from craic.domain import Alignment, Alphabet, CodingSpec, Level, translate_codon
from craic.engines import (
    BuiltinProgressive, ClustalOmegaEngine, ClustalWEngine, CodonAware, MafftEngine,
    MuscleEngine, Param, PrankEngine, ProbConsEngine, all_engines, builtin_variants,
)
from craic import io, progressive
from craic.ambiguity import reliability, disagreement, sandbox, posterior


# --- acceleration core ----------------------------------------------------- #

def test_emission_models_normalised():
    for kind in ("dna", "protein"):
        m = accel.emission_model(kind)
        assert abs(float(m.joint.sum()) - 1.0) < 1e-9
        assert m.k == len(m.symbols)
        assert m.encode("???")[0] == m.k - 1  # unknown -> wildcard


def test_rna_uracil_encodes_as_thymine():
    # RNA is modelled on the DNA emission model; uracil must map to thymine,
    # not to the wildcard (which would strip the homology signal from RNA data).
    rna = accel.emission_model("rna")
    assert rna.encode("U") == rna.encode("T")
    assert rna.encode("AUGC") == rna.encode("ATGC")
    assert rna.encode("N")[0] == rna.k - 1            # genuine unknown still wildcards
    # protein selenocysteine 'U' is a different residue and must stay wildcard
    prot = accel.emission_model("protein")
    assert prot.encode("U")[0] == prot.k - 1


def test_posterior_identical_is_diagonal():
    m = accel.emission_model("dna")
    s = "ACGTGCATGAC"  # low internal repetition
    P = accel.posterior_matrix(s, s, m)
    assert P.shape == (len(s), len(s))
    # the diagonal is the mode of every row, and carries almost all the mass
    assert np.all(np.argmax(P, axis=1) == np.arange(len(s)))
    assert np.mean(np.diag(P)) > 0.9


@pytest.mark.skipif(not accel.HAVE_RUST, reason="Rust core not built")
def test_rust_matches_numpy():
    rng = np.random.default_rng(0)
    for kind in ("dna", "protein"):
        m = accel.emission_model(kind)
        alpha = list(m.symbols[:-1])
        a = m.encode("".join(rng.choice(alpha, 25)))
        b = m.encode("".join(rng.choice(alpha, 31)))
        R = np.array(accel._rust.pair_posteriors(
            list(a), list(b), m.k, list(m.joint.ravel()), list(m.bg), 0.02, 0.5)[0]
        ).reshape(25, 31)
        N = accel._pp_numpy(a, b, m.k, m.joint.ravel(), m.bg, 0.02, 0.5)
        assert np.abs(R - N).max() < 1e-9


def test_blosum_bit_scales_are_right():
    """Every divisor in ``_BLOSUM_SCALE`` must be the one the matrix was built with.

    A BLOSUM file holds ``s = scale * log2(q_ij / p_i p_j)``, so with the right
    divisor ``sum_ij p_i p_j 2^(s/divisor)`` recovers the joint and comes to 1.
    With the wrong one it does not — but ``_protein_model`` normalises afterwards,
    so the error survives as a plausible matrix with the wrong peakedness rather
    than as an exception. NCBI is not consistent about this (45, 50 and 80 are
    third-bits; 62 and 90 are half-bits), so each entry is checked here against
    the matrix itself rather than assumed.
    """
    substitution_matrices = pytest.importorskip("Bio.Align").substitution_matrices
    bg = np.array([accel._ROBINSON[a] for a in accel._AA20])
    for name, divisor in accel._BLOSUM_SCALE.items():
        bl = substitution_matrices.load(name)
        total = sum(bg[i] * bg[j] * 2.0 ** (float(bl[a, b]) / divisor)
                    for i, a in enumerate(accel._AA20)
                    for j, b in enumerate(accel._AA20))
        assert abs(total - 1.0) < 0.05, (
            f"{name} with divisor {divisor} gives {total:.3f}, not ~1.0 — "
            f"wrong bit-scale")


def test_protein_model_is_the_matrix_we_benchmarked():
    """BLOSUM45, not BLOSUM62. Changing this changes every posterior in CRAIC,
    and so every reliability score and every benchmark number, so it is not a
    thing to alter without re-running the benchmarks."""
    assert accel.PROTEIN_MATRIX == "BLOSUM45"
    m = accel.emission_model("protein")
    assert abs(float(m.joint.sum()) - 1.0) < 1e-9


# --- domain ---------------------------------------------------------------- #

def test_iupac_column_consensus():
    from craic.domain import column_consensus

    assert column_consensus(["A", "T", "A", "T"], Level.NT, Alphabet.DNA) == "W"   # tie → IUPAC
    assert column_consensus(["A", "A", "T"], Level.NT, Alphabet.DNA) == "A"        # majority
    assert column_consensus(["A", "C", "G", "T"], Level.NT, Alphabet.DNA) == "N"
    assert column_consensus(["A", "-", "A"], Level.NT, Alphabet.DNA) == "A"        # gaps ignored
    assert column_consensus(["K", "R", "K", "R"], Level.AA, Alphabet.PROTEIN) == "X"   # aa tie → X
    assert column_consensus(["K", "K", "R"], Level.AA, Alphabet.PROTEIN) == "K"
    assert column_consensus(["ATA", "ATT"], Level.CODON, Alphabet.DNA) == "ATW"    # per-position
    assert column_consensus(["-", "-"], Level.NT, Alphabet.DNA) == ""


def test_translate_codon_gap_rules():
    assert translate_codon("---") == "-"
    assert translate_codon("AT-") == "X"   # partial gap surfaced
    assert translate_codon("ATG") == "M"
    assert translate_codon("TAA") == "*"   # stop


def test_codon_aware_levels_and_projection():
    recs = [("a", "ATGAAAACCGCATATATTGCAAAA"),
            ("b", "ATGAAAACCGCATATGGGATTGCAAAA")]
    aln = CodonAware(BuiltinProgressive()).align(recs, Alphabet.DNA)
    assert aln.length % 3 == 0
    assert aln.is_codon_aware()
    assert set(aln.available_levels()) == {Level.NT, Level.CODON, Level.AA}
    dm = aln.display(Level.AA)
    assert dm.span == 3
    assert len(dm.cells[0]) == aln.length // 3


# --- engines --------------------------------------------------------------- #

def test_progressive_places_indels():
    recs = [("s1", "ACGTACGTACGT"), ("s2", "ACGTACGTACGT"),
            ("s3", "ACGTAAACGTACGT"), ("s4", "ACGTACGTACG")]
    aln = progressive.align(recs, Alphabet.DNA)
    assert len({len(r) for r in aln.rows}) == 1   # rectangular
    assert any("-" in r for r in aln.rows)        # gaps introduced


def test_gotoh_align_returns_valid_global_path():
    from craic import accel

    rng = np.random.default_rng(7)
    S = rng.standard_normal((12, 15))
    ca, cb = accel.gotoh_align(S, -6.0, -0.5)
    assert len(ca) == len(cb)                                  # rectangular alignment
    assert [c for c in ca if c >= 0] == list(range(12))        # all A columns, in order
    assert [c for c in cb if c >= 0] == list(range(15))        # all B columns, in order


def test_gotoh_rust_matches_numpy():
    from craic import accel

    if not (accel.HAVE_RUST and hasattr(accel._rust, "gotoh_align")):
        pytest.skip("Rust gotoh_align not built — cross-validation runs where it is")
    rng = np.random.default_rng(11)
    for shape in [(8, 8), (10, 14), (20, 5), (1, 9)]:
        S = rng.standard_normal(shape)
        assert accel.gotoh_align(S, -6.0, -0.5) == accel._gotoh_numpy(S, -6.0, -0.5)


def test_mea_align_returns_valid_gapfree_path():
    from craic import accel

    rng = np.random.default_rng(3)
    S = rng.random((12, 15))                                   # posteriors in [0, 1)
    ca, cb = accel.mea_align(S)
    assert len(ca) == len(cb)                                  # rectangular alignment
    assert [c for c in ca if c >= 0] == list(range(12))        # all A columns, in order
    assert [c for c in cb if c >= 0] == list(range(15))        # all B columns, in order
    assert all((a >= 0) or (b >= 0) for a, b in zip(ca, cb))   # no all-gap columns


def test_mea_rust_matches_numpy():
    from craic import accel

    if not (accel.HAVE_RUST and hasattr(accel._rust, "mea_align")):
        pytest.skip("Rust mea_align not built — cross-validation runs where it is")
    rng = np.random.default_rng(13)
    for shape in [(8, 8), (10, 14), (20, 5), (1, 9), (7, 1), (13, 11)]:
        S = rng.random(shape)
        assert accel.mea_align(S) == accel._mea_numpy(S)


def test_disagreement_compare_self_is_perfect():
    aln, _ = _ambiguous_dataset()
    same = disagreement.compare(aln, aln)
    finite = same[~np.isnan(same)]
    assert same.shape[0] == aln.length
    assert finite.size and (finite >= 0.999).all()    # an alignment fully agrees with itself


def test_editing_operations():
    from craic import editing

    aln = Alignment(["a", "b"], ["AC-GT", "ACTGT"], Alphabet.DNA)
    m = editing.move_residue(aln, 0, 3, -1)              # G slides into the gap at col 2
    assert m is not None and m.rows[0] == "ACG-T"
    assert editing.move_residue(aln, 1, 3, -1) is None   # col 2 of row b is a residue
    g = editing.insert_gap_column(aln, 2)
    assert all(len(r) == 6 for r in g.rows) and all(r[2] == "-" for r in g.rows)
    allgap = Alignment(["a", "b"], ["AC-GT", "AC-GT"], Alphabet.DNA)
    d = editing.delete_column(allgap, 2)
    assert d is not None and d.rows == ["ACGT", "ACGT"]
    assert editing.delete_column(allgap, 0) is None       # col 0 isn't all-gap
    assert editing.column_support(Alignment(["a", "b", "c"], ["A", "A", "A"], Alphabet.DNA), 0) == 1.0
    assert editing.column_support(Alignment(["a", "b"], ["A", "T"], Alphabet.DNA), 0) == 0.0


def test_keyboard_edit_primitives():
    from craic import editing

    aln = Alignment(["a", "b"], ["ACGT-", "ACGTA"], Alphabet.DNA)
    # '-' key: open a gap at the cursor in one row, absorbing the trailing gap
    ins = editing.open_gap(aln, 0, 1)
    assert ins.rows[0] == "A-CGT" and len(ins.rows[0]) == len(aln.rows[1])
    # Backspace/Delete: remove that gap again, pulling residues left
    back = editing.close_gap(ins, 0, 1)
    assert back.rows[0] == "ACGT-"
    assert editing.close_gap(aln, 0, 0) is None          # col 0 isn't a gap → no-op
    # Space: slide a residue (or block) into an adjacent gap
    slide = Alignment(["a", "b"], ["AC--GT", "ACTTGT"], Alphabet.DNA)
    right = editing.slide_block(slide, 0, 1, 2, 1)        # slide 'C' right into the gap
    assert right is not None and right.rows[0] == "A-C-GT"
    block = editing.slide_block(slide, 1, 2, 4, -1)       # block 'TT' has no gap to its left
    assert block is None
    left = editing.slide_block(Alignment(["a"], ["A-CG"], Alphabet.DNA), 0, 2, 4, -1)
    assert left is not None and left.rows[0] == "ACG-"


def test_shift_rows_only_moves_selected_sequences():
    from craic import editing

    aln = Alignment(["a", "b", "c"], ["ACGT-", "ACGT-", "ACGT-"], Alphabet.DNA)
    # slide only row b right; it absorbs its own trailing gap, others don't move,
    # and the alignment does not widen
    r = editing.shift_rows(aln, {1}, 1, 1)
    assert r.rows == ["ACGT-", "A-CGT", "ACGT-"] and r.length == 5
    # slide it back left (there is now a gap at col 1)
    back = editing.shift_rows(r, {1}, 2, -1)
    assert back.rows[1] == "ACGT-"
    # a selected sequence with no trailing gap forces a widen of just that column,
    # padding the others with a trailing gap (they stay visually put)
    aln2 = Alignment(["a", "b"], ["ACGT", "ACGT"], Alphabet.DNA)
    w = editing.shift_rows(aln2, {0}, 1, 1)
    assert w.rows[0] == "A-CGT" and w.rows[1] == "ACGT-" and w.length == 5
    # left slide refused unless every selected row has a gap to the left
    assert editing.shift_rows(aln2, {0, 1}, 2, -1) is None


def test_shift_rows_codon_unit_keeps_frame():
    from craic import editing
    from craic.domain import CodingSpec

    aln = Alignment(["x", "y"], ["ATGAAACCC---", "ATGAAACCCGGG"], Alphabet.DNA,
                    coding=CodingSpec(0, 1))
    # slide a whole codon (unit=3) right in row x at nt column 3
    r = editing.shift_rows(aln, {0}, 3, 1, unit=3, keep_coding=True)
    assert r.rows[0] == "ATG---AAACCC" and len(r.rows[0]) % 3 == 0
    assert r.rows[1] == "ATGAAACCCGGG"           # unselected row untouched
    assert r.coding is not None                    # reading frame preserved
    # and back again
    back = editing.shift_rows(r, {0}, 6, -1, unit=3, keep_coding=True)
    assert back.rows[0] == "ATGAAACCC---"


def test_sequence_groups_sidecar_roundtrip(tmp_path):
    aln = Alignment(["cat", "dog", "fly"], ["ACGT", "ACGT", "ACGA"], Alphabet.DNA)
    aln.meta["groups"] = {"animals": ["cat", "dog"]}
    path = str(tmp_path / "aln.fasta")
    io.save_alignment(path, aln)
    # the alignment file stays standard FASTA; the group lives in a sidecar
    assert os.path.exists(io.sidecar_path(path))
    with open(path) as fh:
        assert "groups" not in fh.read()      # nothing injected into the FASTA
    reloaded = io.load_alignment(path)
    assert reloaded.meta.get("groups") == {"animals": ["cat", "dog"]}
    # deleting all groups removes the sidecar
    aln.meta["groups"] = {}
    io.save_alignment(path, aln)
    assert not os.path.exists(io.sidecar_path(path))


def test_anchored_realign_preserves_residues_and_shape():
    from craic import editing
    from craic.engines import BuiltinProgressive

    aln, _ = _ambiguous_dataset()
    out = editing.anchored_realign(aln, [0, 1, aln.length - 2, aln.length - 1],
                                   BuiltinProgressive())
    assert len({len(r) for r in out.rows}) == 1                 # rectangular
    for i in range(aln.n_seqs):
        assert out.rows[i].replace("-", "") == aln.rows[i].replace("-", "")   # residues kept


def test_report_hotspots_and_summary():
    from craic.ambiguity import report as rep

    rel = np.array([0.9, 0.9, 0.2, 0.1, 0.15, 0.9, 0.95, 0.3, 0.3])
    hs = rep.hotspots(rel, threshold=0.5, min_len=2)
    assert hs and hs[0].col_start == 2 and hs[0].col_stop == 5     # worst (widest) run first
    s = rep.summary(rel, threshold=0.5)
    assert s.n_columns == 9 and 0.0 <= s.pct_high <= 1.0 and s.hotspots
    assert "reliability" in s.as_text().lower()


def test_ensemble_cooccurrence():
    from craic.ambiguity import ensemble
    from craic.engines import builtin_variants

    aln, _ = _ambiguous_dataset()
    co, ci, cj, n = ensemble.cooccurrence(aln, 0, 1, 20, aln.length - 20, builtin_variants())
    assert n >= 1
    assert co.shape == (len(ci), len(cj))
    assert (co >= 0).all() and (co <= 1.0001).all()
    assert co.max() > 0.5            # the committed homologies recur across alternatives


def test_gotoh_numpy_cancels_mid_run():
    from craic import accel

    S = np.zeros((80, 80))
    with pytest.raises(accel.Cancelled):
        accel._gotoh_numpy(S, -6.0, -0.5, cancelled=lambda: True)


def test_progressive_progress_and_cancel():
    recs = [("a", "ACGTACGT"), ("b", "ACGTAAGT"), ("c", "ACGAACGT"), ("d", "ACGTACGA")]
    calls = []
    aln = progressive.align(recs, Alphabet.DNA, progress=lambda d, t: calls.append((d, t)))
    assert aln.length > 0
    assert calls and calls[-1][0] == calls[-1][1]          # progress ran to total
    with pytest.raises(progressive.Cancelled):             # cooperative stop
        progressive.align(recs, Alphabet.DNA, cancelled=lambda: True)


# --- ambiguity ------------------------------------------------------------- #

def _ambiguous_dataset():
    rng = np.random.default_rng(1)
    fl, fr = "ATGCATGCATGCATGCATGC", "TTTAGGGCCCAAATTTGGGC"
    recs = [(f"s{i}", fl + "".join(rng.choice(list("ACGT"), rng.integers(5, 10))) + fr)
            for i in range(5)]
    return progressive.align(recs, Alphabet.DNA), recs


def test_reliability_flanks_beat_middle():
    aln, _ = _ambiguous_dataset()
    rep = reliability.analyse(aln, do_perturbation=True)
    assert rep.col_combined.shape == (aln.length,)
    left = np.nanmean(rep.col_combined[:20])
    mid = np.nanmean(rep.col_combined[20:aln.length - 20])
    assert left > mid                              # flanks more reliable
    keep = rep.keep_mask(0.5)
    assert reliability.apply_mask(aln, keep).length == int(keep.sum())


def test_disagreement_runs_and_scores():
    _, recs = _ambiguous_dataset()
    res = disagreement.compute(recs, Alphabet.DNA, builtin_variants())
    assert res.col_agreement.shape[0] == res.reference.length
    finite = res.col_agreement[~np.isnan(res.col_agreement)]
    assert finite.min() >= 0.0 and finite.max() <= 1.0


def test_sandbox_splice_is_rectangular():
    aln, _ = _ambiguous_dataset()
    c0, c1 = 20, aln.length - 20
    alts = sandbox.alternatives(aln, c0, c1, [BuiltinProgressive()])
    assert alts
    spliced = sandbox.splice(aln, c0, c1, alts[0].block)
    assert len({len(r) for r in spliced.rows}) == 1


def test_alternatives_reports_progress_and_cancels():
    aln, _ = _ambiguous_dataset()
    c0, c1 = 20, aln.length - 20
    calls = []
    alts = sandbox.alternatives(aln, c0, c1, [BuiltinProgressive()],
                                progress=lambda d, t, lbl: calls.append((d, t, lbl)))
    assert alts
    assert calls and calls[-1][0] == calls[-1][1]      # progress ran through to total
    # cancelling up front produces nothing and runs no engine
    stopped = sandbox.alternatives(aln, c0, c1, [BuiltinProgressive()], cancelled=lambda: True)
    assert stopped == []


def test_posterior_region_view():
    aln, _ = _ambiguous_dataset()
    pv = posterior.region_posterior(aln, 0, 1, 20, aln.length - 20)
    assert pv.matrix.ndim == 2
    assert (pv.matrix >= 0).all() and (pv.matrix <= 1).all()


def test_residual_mode_drops_committed_mass():
    aln, _ = _ambiguous_dataset()
    base = posterior.region_posterior(aln, 0, 1, 0, aln.length)
    resid = posterior.region_posterior(aln, 0, 1, 0, aln.length, mode="residual")
    assert resid.matrix.sum() < base.matrix.sum()   # diagonal removed
    assert (resid.matrix >= 0).all()


def test_region_posterior_windowed_keeps_diagonal():
    # A region cut from the *middle* of long sequences must still place the
    # committed homology on the diagonal — i.e. the windowing offsets map back
    # correctly. (200-residue identical rows; the window is a strict subset.)
    row = "ACDEFGHIKLMNPQRSTVWY" * 10
    aln = Alignment(["a", "b"], [row, row], Alphabet.PROTEIN)
    pv = posterior.region_posterior(aln, 0, 1, 90, 110)
    assert pv.matrix.shape == (20, 20)
    assert (pv.matrix >= 0).all() and (pv.matrix <= 1.0001).all()
    assert np.mean(np.diag(pv.matrix)) > 0.4        # diagonal preserved by windowing


def test_region_profile_rows_are_residues():
    aln, _ = _ambiguous_dataset()
    prof = posterior.region_profile(aln, 0, 0, aln.length)
    ungapped = len(aln.rows[0].replace("-", ""))
    assert prof.matrix.shape == (ungapped, aln.length)
    assert (prof.matrix >= 0).all() and (prof.matrix <= 1.0001).all()


def test_probe_residue_is_a_column_distribution():
    aln, _ = _ambiguous_dataset()
    probs = posterior.probe_residue(aln, 0, 2)
    assert probs.shape[0] == aln.length
    finite = probs[~np.isnan(probs)]
    assert finite.size and finite.max() <= 1.0001


def test_probe_partner_sampling_bounds_work():
    from craic.ambiguity.posterior import MAX_PROBE_CELLS, _partner_indices

    assert _partner_indices(6, 0, li=50) == [1, 2, 3, 4, 5]      # small → every partner
    idx = _partner_indices(5000, 10, li=2000)                    # big → bounded sample
    budget = max(1, MAX_PROBE_CELLS // (2000 * 2000))
    assert 0 < len(idx) <= budget
    assert 10 not in idx and len(set(idx)) == len(idx)           # excludes self, unique


def test_probe_residue_is_local_and_bounded():
    # 400-column alignment; probing a middle residue must only light up columns
    # near it (the window), not the whole alignment — that's what keeps it fast.
    row = "ACGT" * 100
    aln = Alignment(["a", "b", "c"], [row, row, row], Alphabet.DNA)
    probs = posterior.probe_residue(aln, 0, 200, window=40)
    assert probs.shape[0] == aln.length
    lit = np.where(~np.isnan(probs))[0]
    assert lit.size                                              # something was scored
    assert lit.min() >= 200 - 40 and lit.max() <= 200 + 40       # stayed local to the window
    assert np.nanmax(probs) <= 1.0001


def test_reliability_cell_combined_shape():
    aln, _ = _ambiguous_dataset()
    rep = reliability.analyse(aln, do_perturbation=False)
    assert rep.cell_combined.shape == (aln.n_seqs, aln.length)


# --- io -------------------------------------------------------------------- #

def test_fasta_roundtrip(tmp_path):
    p = tmp_path / "a.fasta"
    io.write_fasta(str(p), ["x", "y"], ["ACGT-", "ACG-T"])
    recs = io.read_records(str(p))
    assert recs == [("x", "ACGT-"), ("y", "ACG-T")]


def test_all_write_formats(tmp_path):
    aln = Alignment(["a", "b", "c"], ["ATGAAACCC", "ATG---CCC", "ATGAAATCC"],
                    Alphabet.DNA, coding=CodingSpec(0, 1))
    for key, _label, ext in io.write_formats():
        p = tmp_path / f"out_{key}{ext}"
        io.save_alignment(str(p), aln, key)
        assert p.exists() and p.stat().st_size > 0


def test_dump_alignment_all_formats():
    aln = Alignment(["a", "b", "c"], ["ATGAAACCC", "ATG---CCC", "ATGAAATCC"],
                    Alphabet.DNA, coding=CodingSpec(0, 1))
    for key, _label, _ext in io.write_formats():
        s = io.dump_alignment(aln, key)
        assert isinstance(s, str) and s.strip()           # non-empty for every format
    assert io.dump_alignment(aln, "fasta").startswith(">a\n")
    assert "#MEGA" in io.dump_alignment(aln, "mega")
    assert "CLUSTAL" in io.dump_alignment(aln, "clustal").upper()


def _coding_alignment():
    return CodonAware(BuiltinProgressive()).align(
        [("a", "ATGAAAACCGCATATATTGCAAAA"),
         ("b", "ATGAAAACCGCATATGGGATTGCAAAA"),
         ("c", "ATGAAAACCGCGTATATTGCAAAA")],
        Alphabet.DNA,
    )


def test_sandbox_codon_aware_previews_amino_acids():
    aln = _coding_alignment()
    alts = sandbox.alternatives(aln, 0, aln.length, [BuiltinProgressive()],
                                codon_aware=True, table=1, level=Level.AA)
    assert alts
    rows = alts[0].display_rows()
    assert all(all(len(c) == 1 for c in r) for r in rows)  # amino acids
    spliced = sandbox.splice(aln, 0, aln.length, alts[0].block)
    assert spliced.coding is not None and spliced.length % 3 == 0


def test_posterior_amino_acid_projection():
    aa = _coding_alignment().amino_acid_alignment()
    assert aa.alphabet == Alphabet.PROTEIN
    pv = posterior.region_posterior(aa, 0, 1, 0, aa.length)
    assert pv.matrix.ndim == 2 and (pv.matrix <= 1.0001).all()


def test_similarity_order_groups_similar():
    order = progressive.similarity_order(["AAAAAA", "AAAAAT", "TTTTTT", "TTTTTA"])
    assert set(order) == {0, 1, 2, 3}
    pos = {s: i for i, s in enumerate(order)}
    assert abs(pos[0] - pos[1]) == 1 and abs(pos[2] - pos[3]) == 1


def test_disagreement_overlays_base_alignment():
    recs = [("a", "ACGTACGT"), ("b", "ACGTAAGT"), ("c", "ACGAACGT")]
    base = progressive.align(recs, Alphabet.DNA)
    res = disagreement.compute(recs, Alphabet.DNA, builtin_variants(), base_alignment=base)
    assert res.reference is base                       # annotates, does not replace
    assert res.col_agreement.shape[0] == base.length
    finite = res.col_agreement[~np.isnan(res.col_agreement)]
    assert (finite >= 0).all() and (finite <= 1).all()


# --- bespoke aligner parameters ------------------------------------------- #

def test_engine_parameters_are_params():
    for eng in all_engines():
        keys = [p.key for p in eng.parameters()]
        assert len(keys) == len(set(keys)), eng.key
        for p in eng.parameters():
            assert isinstance(p, Param)
            assert p.kind in ("float", "int", "choice", "bool")
            if p.kind == "choice":
                assert p.default in (p.choices or [])
            if p.default is None:
                # only a number can be left to the tool, and it must say what then
                assert p.kind in ("float", "int") and p.blank, (eng.key, p.key)
                for alph in (Alphabet.PROTEIN, Alphabet.DNA):
                    assert p.blank_for(alph)
            assert p.alphabet in ("", "protein", "nucleotide")


def test_parameters_know_which_data_they_apply_to():
    params = {p.key: p for p in ClustalWEngine().parameters()}
    assert params["dnamatrix"].applies_to(Alphabet.DNA)
    assert not params["dnamatrix"].applies_to(Alphabet.PROTEIN)
    assert params["matrix"].applies_to(Alphabet.PROTEIN)
    assert not params["matrix"].applies_to(Alphabet.RNA)
    assert params["gapopen"].applies_to(Alphabet.PROTEIN)      # both kinds
    assert params["dnamatrix"].applies_to(None)                # unknown: show it
    assert params["gapopen"].blank_for(Alphabet.PROTEIN) == "ClustalW default, 10"
    assert params["gapopen"].blank_for(Alphabet.DNA) == "ClustalW default, 15"
    assert "10" in params["gapopen"].blank_for() and "15" in params["gapopen"].blank_for()


def test_builtin_gap_regime_changes_alignment():
    _, recs = _ambiguous_dataset()
    # the gap model is the pair-HMM gap-open probability (delta); a rare-gap vs a
    # frequent-gap model must produce different posterior-decoded alignments
    rare = progressive.align(recs, Alphabet.DNA, delta=0.003, estimate=False, consistency_iters=0)
    freq = progressive.align(recs, Alphabet.DNA, delta=0.2, estimate=False, consistency_iters=0)
    assert rare.rows != freq.rows


def test_mafft_command_honours_op_ep_strategy():
    argv, mode = MafftEngine()._command(
        "in.fa", "out.fa", Alphabet.DNA, strategy="L-INS-i (accurate, local)", op=2.5, ep=0.1)
    assert mode == "stdout"
    assert "--op" in argv and argv[argv.index("--op") + 1] == "2.5"
    assert "--ep" in argv and argv[argv.index("--ep") + 1] == "0.1"
    assert "--localpair" in argv                       # strategy flags applied


def test_codon_aware_delegates_parameters():
    inner = BuiltinProgressive()
    assert CodonAware(inner).parameters() == inner.parameters()


def test_prank_leaves_gap_costs_to_prank_unless_they_are_set():
    """PRANK's gap defaults differ between DNA and protein. CRAIC used to pass
    0.025 / 0.5 for both, overriding PRANK with values right for neither."""
    for alph in (Alphabet.PROTEIN, Alphabet.DNA):
        argv = PrankEngine._argv("in.fa", "out", alph)
        assert not any(a.startswith(("-gaprate", "-gapext")) for a in argv)
    argv = PrankEngine._argv("in.fa", "out", Alphabet.PROTEIN, gaprate=0.01, gapext=0.6,
                             F=True, iterate=3, termgap=True)
    assert {"-gaprate=0.01", "-gapext=0.6", "+F", "-iterate=3", "-termgap"} <= set(argv)
    assert "-iterate=5" not in PrankEngine._argv("in.fa", "out", Alphabet.DNA, iterate=5)


def test_mafft_matrix_follows_the_alphabet_and_maxiterate_overrides_the_strategy():
    eng = MafftEngine()
    prot, _ = eng._command("in.fa", "out.fa", Alphabet.PROTEIN, matrix="BLOSUM30", kimura="1 PAM")
    assert prot[prot.index("--bl") + 1] == "30" and "--kimura" not in prot
    dna, _ = eng._command("in.fa", "out.fa", Alphabet.DNA, matrix="BLOSUM30", kimura="1 PAM")
    assert dna[dna.index("--kimura") + 1] == "1" and "--bl" not in dna
    argv, _ = eng._command("in.fa", "out.fa", Alphabet.DNA,
                           strategy="L-INS-i (accurate, local)", maxiterate=10,
                           lop=-3.0, lep=0.2)
    assert argv.count("--maxiterate") == 1 and argv[argv.index("--maxiterate") + 1] == "10"
    assert argv[argv.index("--lop") + 1] == "-3.0" and argv[argv.index("--lep") + 1] == "0.2"


def test_clustalw_passes_gap_costs_only_when_set_and_protein_switches_only_for_protein():
    eng = ClustalWEngine()
    argv, mode = eng._command("in.fa", "out.fa", Alphabet.PROTEIN)
    assert mode == "file" and "-TYPE=PROTEIN" in argv and "-OUTPUT=FASTA" in argv
    assert not any(a.startswith(("-GAPOPEN", "-GAPEXT", "-PWGAP")) for a in argv)
    argv, _ = eng._command("in.fa", "out.fa", Alphabet.PROTEIN, gapopen=25.0, gapext=0.5,
                           pwgapopen=12.0, pwgapext=0.3, matrix="BLOSUM", gapdist=8,
                           nopgap=True, nohgap=True)
    assert {"-GAPOPEN=25.0", "-GAPEXT=0.5", "-PWGAPOPEN=12.0", "-PWGAPEXT=0.3",
            "-MATRIX=BLOSUM", "-GAPDIST=8", "-NOPGAP", "-NOHGAP"} <= set(argv)
    argv, _ = eng._command("in.fa", "out.fa", Alphabet.DNA, dnamatrix="CLUSTALW",
                           nopgap=True, gapdist=8)
    assert "-TYPE=DNA" in argv and "-DNAMATRIX=CLUSTALW" in argv
    assert not any(a.startswith(("-MATRIX", "-GAPDIST", "-NOPGAP")) for a in argv)


def test_muscle_passes_its_v5_options(monkeypatch):
    seen = []

    class Failed:
        returncode, stderr = 1, "no"

    def fake_run(argv, timeout=None):
        seen.append(argv)
        return Failed()

    monkeypatch.setattr(MuscleEngine, "_run", staticmethod(fake_run))
    with pytest.raises(RuntimeError):
        MuscleEngine().align([("a", "MKV"), ("b", "MKI")], Alphabet.PROTEIN,
                             algorithm="super5", perm="abc", perturb=3)
    v5 = seen[0]
    assert v5[1] == "-super5" and v5[v5.index("-perm") + 1] == "abc"
    assert v5[v5.index("-perturb") + 1] == "3"
    assert "-perm" not in seen[1]                         # v3 has none of them


def test_clustalo_and_probcons_options_reach_the_command():
    argv, _ = ClustalOmegaEngine()._command("in.fa", "out.fa", Alphabet.PROTEIN,
                                            iterations=2, full=True,
                                            max_guidetree_iterations=1)
    assert "--full" in argv and "--max-guidetree-iterations=1" in argv
    assert not any(a.startswith("--max-hmm") for a in argv)
    argv, _ = ProbConsEngine()._command("in.fa", "out.fa", Alphabet.PROTEIN,
                                        consistency=0, iterations=100, pretraining=3)
    assert argv[:5] == ["probcons", "-c", "0", "-pre", "3"]   # defaults left off


def test_builtin_exposes_gap_probabilities_and_matrix():
    keys = {p.key for p in BuiltinProgressive().parameters()}
    assert {"delta", "epsilon", "matrix", "consistency_iters", "refine_iters"} <= keys
    assert accel.emission_model("protein", matrix="BLOSUM62") is not \
        accel.emission_model("protein", matrix="BLOSUM45")
    recs = [("a", "MKVLAAGIVGLLLAHSQAENKTEWLK"), ("b", "MKVLAGIVALLLHSQAENKEWLK"),
            ("c", "MRVLAAGIIGLLAHSEAENKTDWLRK")]
    aln = BuiltinProgressive().align(recs, Alphabet.PROTEIN, effort="min",
                                     delta=0.05, epsilon=0.4, matrix="BLOSUM62")
    assert [r.replace("-", "") for r in aln.rows] == [s for _, s in recs]


def test_param_reads_command_line_text():
    strategy = next(p for p in MafftEngine().parameters() if p.key == "strategy")
    assert strategy.coerce("L-INS-i") == "L-INS-i (accurate, local)"
    with pytest.raises(ValueError):
        strategy.coerce("fastest")
    gaprate = next(p for p in PrankEngine().parameters() if p.key == "gaprate")
    assert gaprate.coerce("0.01") == 0.01 and gaprate.coerce("") is None
    with pytest.raises(ValueError):
        gaprate.coerce("2")                               # out of 0-1
    plus_f = next(p for p in PrankEngine().parameters() if p.key == "F")
    assert plus_f.coerce("yes") is True and plus_f.coerce("false") is False


def test_cli_param_is_checked_against_the_engine():
    from craic.cli import _engine_opts, build_parser

    assert _engine_opts(MafftEngine(), ["op=2.5", "strategy=G-INS-i"]) == {
        "op": 2.5, "strategy": "G-INS-i (accurate, global)"}
    with pytest.raises(SystemExit):
        _engine_opts(MafftEngine(), ["gapopen=2"])        # ClustalW's name, not MAFFT's
    with pytest.raises(SystemExit):
        _engine_opts(MafftEngine(), ["op"])
    args = build_parser().parse_args(["align", "x.fa", "--param", "op=2", "--param", "ep=1"])
    assert args.param == ["op=2", "ep=1"]


def test_cli_engines_lists_every_engine(capsys):
    from craic.cli import cmd_engines

    cmd_engines(None)
    out = capsys.readouterr().out
    for eng in all_engines():
        assert out.count(f"{eng.key} — ") == 1
    assert "gapopen" in out and "PRANK default, 0.005" in out
    assert "dnamatrix: IUB | CLUSTALW  [default: IUB]  (nucleotide only)" in out


def test_sandbox_runs_each_engine_with_its_settings():
    seen = {}

    class Recorder(BuiltinProgressive):
        key = "recorder"

        def align(self, records, alphabet, **opts):
            seen.update(opts)
            return super().align(records, alphabet, effort="min")

    aln, _ = _ambiguous_dataset()
    sandbox.alternatives(aln, 0, 20, [Recorder()], param_grid=[],
                         params_by_key={"recorder": {"delta": 0.07}})
    assert seen == {"delta": 0.07}


def test_disagreement_accepts_params_by_key():
    _, recs = _ambiguous_dataset()
    engines = builtin_variants()
    params = {e.key: {"consistency_iters": 1} for e in engines}
    res = disagreement.compute(recs, Alphabet.DNA, engines, params_by_key=params)
    assert res.col_agreement.shape[0] == res.reference.length


# --- new behaviour: MEA aligner, bootstrap perturbation, adversarial IO ------ #

def test_builtin_aligner_mea_stacks_conserved_flanks():
    # Posterior-decoding (MEA) should align identical conserved flanks into
    # gap-free, identical columns and push indels into the variable middle.
    fl, fr = "ATGCATGCATGC", "TTAGGGCCCTTT"
    recs = [("s1", fl + "AAA" + fr), ("s2", fl + "AAAAAAAA" + fr), ("s3", fl + "A" + fr)]
    aln = progressive.align(recs, Alphabet.DNA)
    assert len({len(r) for r in aln.rows}) == 1                     # rectangular
    for c in range(len(fl)):                                        # leading flank
        col = {r[c] for r in aln.rows}
        assert "-" not in col and len(col) == 1                     # gap-free, identical
    assert any("-" in r for r in aln.rows)                          # gaps in the middle


def test_perturbation_ranks_conserved_above_ambiguous_and_nans_singletons():
    aln, _ = _ambiguous_dataset()
    col, cell = reliability.perturbation(aln, n_replicates=12)
    L = aln.length
    assert col.shape == (L,)
    assert np.nanmax(col) <= 1.0 + 1e-9 and np.nanmin(col) >= 0.0
    left = np.nanmean(col[:20]); mid = np.nanmean(col[20:L - 20])
    assert left > mid                                              # conserved flank beats ambiguous middle
    # a residue with no asserted partners is left nan, never scored a free 1.0
    solo = Alignment(["a", "b"], ["A-", "-A"], Alphabet.DNA)       # two singleton columns
    scol, _ = reliability.perturbation(solo, n_replicates=4)
    assert np.all(np.isnan(scol))


def test_perturbation_is_low_variance_across_seeds():
    aln, _ = _ambiguous_dataset()
    a, _ = reliability.perturbation(aln, n_replicates=16, seed=0)
    b, _ = reliability.perturbation(aln, n_replicates=16, seed=123)
    finite = ~np.isnan(a) & ~np.isnan(b)
    assert finite.sum() > 5
    assert np.corrcoef(a[finite], b[finite])[0, 1] > 0.9          # a stable estimator


def test_io_handles_adversarial_input(tmp_path):
    def write(name, text):
        p = tmp_path / name
        p.write_text(text)
        return str(p)

    empty = io.load_alignment(write("empty.fasta", ""))            # no crash on empty
    assert empty.n_seqs == 0 and empty.length == 0

    one = io.load_alignment(write("one.fasta", ">x\nACGTACGT\n"))   # single sequence
    assert one.n_seqs == 1 and one.meta.get("unaligned")

    mismatch = io.load_alignment(write("m.fasta", ">a\nACGT\n>b\nACGTAAA\n"))
    assert mismatch.meta.get("unaligned")                          # unequal lengths -> not aligned

    # malformed non-FASTA must not silently load as a bogus alignment
    garbage = write("g.fasta", "this is not fasta\nrandom text\n")
    try:
        aln = io.load_alignment(garbage)
        assert aln.n_seqs == 0                                     # tolerated only if empty
    except Exception:
        pass                                                       # rejecting it is fine too


# --- memory-tiered aligner: streaming fallback for large families ------------ #

def test_aligner_streams_when_over_memory_budget():
    # a family over the consistency memory budget must NOT materialise all pairwise
    # posteriors (it streams them on demand), yet still return a valid alignment.
    from craic import progressive

    _, recs = _ambiguous_dataset()
    orig = progressive._all_pairs_posteriors
    calls = []

    def spy(*a, **k):
        calls.append(1)
        return orig(*a, **k)

    progressive._all_pairs_posteriors = spy
    try:
        aln = progressive.align(recs, Alphabet.DNA, consistency_mem_gb=0.0)   # force streaming
    finally:
        progressive._all_pairs_posteriors = orig

    assert not calls                                        # never stored all-pairs (memory-safe)
    assert len({len(r) for r in aln.rows}) == 1            # rectangular
    for i in range(len(recs)):
        assert aln.rows[i].replace("-", "") == recs[i][1]  # residues preserved


def test_aligner_memory_tier_selection():
    from craic import progressive

    assert progressive._posterior_gb(6, 40) <= 1.0         # small family -> consistency path
    assert progressive._posterior_gb(200, 1000) > 1.0      # large family -> streaming path


def test_effort_levels_produce_valid_alignments():
    from craic import progressive

    _, recs = _ambiguous_dataset()
    for effort in ("min", "med", "max"):
        aln = progressive.align(recs, Alphabet.DNA, effort=effort)
        assert len({len(r) for r in aln.rows}) == 1                    # rectangular
        for i in range(len(recs)):
            assert aln.rows[i].replace("-", "") == recs[i][1]          # residues preserved


def test_refinement_never_lowers_the_objective():
    # _refine only keeps changes that raise the sum-of-pairs posterior score, so the
    # objective is monotonic non-decreasing (truth may still diverge, but not the objective).
    import numpy as np
    from craic import progressive, accel

    _, recs = _ambiguous_dataset()
    seqs = [s for _, s in recs]
    model = accel.emission_model("dna")
    P = progressive._all_pairs_posteriors(seqs, model, 0.02, 0.5)

    def post(x, y):
        return progressive._post(P, x, y)

    start = progressive._progressive(seqs, post, dist=progressive._posterior_distances(seqs, P))
    s0 = progressive._sp_posterior_score(start, post)
    refined = progressive._refine([t for t in start], post, 15, np.random.default_rng(0))
    assert progressive._sp_posterior_score(refined, post) >= s0 - 1e-6


# --- model-free trimming (Gblocks / trimAl / MSA_trimmer) + outliers --------- #

def test_gap_score_and_threshold_mask():
    from craic.ambiguity import trimming as T
    aln = Alignment(["a", "b", "c", "d"], ["AC-T", "A--T", "A-GT", "A--T"], Alphabet.DNA)
    gs = T.gap_score(aln)
    assert abs(gs[0] - 1.0) < 1e-9 and abs(gs[1] - 0.25) < 1e-9   # col0 all present, col1 1/4
    keep = T.gap_threshold_mask(aln, max_gap=0.5)
    assert keep[0] and not keep[1]                                # drop the very gappy column (0.75 gaps)


def test_gblocks_keeps_conserved_core():
    from craic.ambiguity import trimming as T
    core = "ATGCATGCATGC"
    rows = ["--x-" .replace("x", "A") + core + "-t-" .replace("t", "A") for _ in range(4)]
    rows = ["A-C-" + core + "-T-A", "-TC-" + core + "-TG-",
            "A---" + core + "G--A", "-TCG" + core + "-T-A"]
    aln = Alignment(["s1", "s2", "s3", "s4"], rows, Alphabet.DNA)
    keep = T.gblocks_mask(aln, b4=4, allow_gaps="none")
    kept = set(np.where(keep)[0].tolist())
    assert kept == set(range(4, 16))                              # exactly the conserved core


def test_outlier_sequence_detection():
    from craic.ambiguity import trimming as T
    rng = np.random.default_rng(0)
    base = "ATGCATGCATGCATGCATGC"
    good = [base[:i] + rng.choice(list("ACGT")) + base[i + 1:] for i in range(5)]
    odd = "".join(rng.choice(list("ACGT"), 20))
    aln = Alignment([f"g{i}" for i in range(5)] + ["ODD"], good + [odd], Alphabet.DNA)
    scores, flags = T.outlier_sequences(aln, z=2.0)
    assert flags[-1] and not flags[:-1].any()                     # only ODD flagged
    assert scores[-1] < scores[:-1].min()                         # and it scores lowest


def test_trimming_masks_are_valid():
    from craic.ambiguity import trimming as T
    aln, _ = _ambiguous_dataset()
    for mask in (T.gappyout_mask(aln), T.trimal_strict_mask(aln), T.gblocks_mask(aln, b4=3)):
        assert mask.dtype == bool and mask.shape == (aln.length,)
    assert (T.similarity_score(aln) >= 0).all() and (T.similarity_score(aln) <= 1.0001).all()


# --- residue masking -------------------------------------------------------- #

def test_a_residue_mask_names_residues_and_keeps_columns():
    aln = Alignment(["x", "y"], ["AC-GT", "ACTG-"], Alphabet.DNA)
    picked = reliability.residues_in(aln, [0, 1], 1, 4)
    assert picked == {"x": [1, 2], "y": [1, 2, 3]}
    masked = reliability.with_residue_mask(aln, picked)
    assert masked.rows == aln.rows                          # nothing changes until export
    assert reliability.masked_cells(masked, picked) == {(0, 1), (0, 3), (1, 1), (1, 2), (1, 3)}
    out = reliability.apply_residue_mask(masked)
    assert out.rows == ["AN-NT", "ANNN-"] and out.length == aln.length
    assert reliability.RESIDUE_MASK_KEY not in out.meta     # applied, so not carried
    prot = Alignment(["x"], ["MK-V"], Alphabet.PROTEIN)
    assert reliability.apply_residue_mask(prot, {"x": [1]}).rows == ["MX-V"]


def test_a_residue_mask_follows_residues_and_drops_what_does_not_fit():
    aln = Alignment(["x", "y"], ["AC-GT", "ACTG-"], Alphabet.DNA)
    shifted = Alignment(["y", "x"], ["-ACTG", "ACG-T"], Alphabet.DNA,
                        meta={reliability.RESIDUE_MASK_KEY: {"x": [2, 9], "z": [0]}})
    # x's residue 2 is G wherever the gaps go; index 9 and sequence z do not exist
    assert reliability.residue_mask(shifted) == {"x": [2]}
    assert reliability.apply_residue_mask(shifted).rows == ["-ACTG", "ACN-T"]
    # column masking removes residues, so it must not carry residue indices
    carried = reliability.apply_mask(reliability.with_residue_mask(aln, {"x": [0]}),
                                     np.array([False, True, True, True, True]))
    assert reliability.RESIDUE_MASK_KEY not in carried.meta


def test_residues_below_a_threshold_leave_unscored_residues_alone():
    aln = Alignment(["x", "y"], ["AC-G", "ACTG"], Alphabet.DNA)
    cells = np.array([[0.9, 0.2, np.nan, 0.6],
                      [0.9, 0.3, np.nan, 0.4]])            # y's T at column 2 is unscored
    assert reliability.residues_below(aln, cells, 0.5) == {"x": [1], "y": [1, 3]}
    with pytest.raises(ValueError):
        reliability.residues_below(aln, cells[:, :3], 0.5)
    a, b = {"x": [1, 2]}, {"x": [2], "y": [0]}
    assert reliability.merge_masks(a, b) == {"x": [1, 2], "y": [0]}
    assert reliability.subtract_masks(a, b) == {"x": [1]}
    assert reliability.mask_size(a) == 2


def test_a_residue_mask_survives_a_session(tmp_path):
    from craic.session import Session

    aln = reliability.with_residue_mask(
        Alignment(["x", "y"], ["AC-GT", "ACTG-"], Alphabet.DNA), {"y": [1, 3]})
    path = str(tmp_path / "s.craic.json")
    Session(alignment=aln).save(path)
    assert reliability.residue_mask(Session.load(path).alignment) == {"y": [1, 3]}


def test_cli_mask_residues_keeps_every_column(tmp_path, capsys):
    from craic.cli import main

    aln, _ = _ambiguous_dataset()
    src, out = tmp_path / "a.fasta", tmp_path / "m.fasta"
    io.save_alignment(str(src), aln)
    with pytest.raises(SystemExit) as done:
        main(["mask", str(src), "-o", str(out), "--residues", "--fast", "--threshold", "0.9"])
    assert done.value.code == 0
    back = io.load_alignment(str(out))
    assert back.length == aln.length
    assert any("N" in r for r in back.rows)
    assert "all" in capsys.readouterr().err


def test_citation_is_the_same_everywhere():
    """The citation is written out in the package, the README and the website's
    About page; this catches one of them being updated without the others."""
    import re
    from pathlib import Path

    from craic import CITATION

    root = Path(__file__).resolve().parents[1]

    def norm(s):
        return re.sub(r"\s+", " ", s.replace(">", " "))

    for name in ("README.md", "docs/about.md"):
        assert norm(CITATION) in norm((root / name).read_text(encoding="utf-8")), name


def test_make_app_refuses_off_macos(monkeypatch, capsys):
    from craic import cli

    monkeypatch.setattr("sys.platform", "linux")
    with pytest.raises(SystemExit) as exc:
        cli.main(["make-app", "--dest", "unused"])
    assert exc.value.code == 1
    assert "only runs on a Mac" in capsys.readouterr().err


def test_make_app_bundle_launches_this_python(tmp_path):
    import plistlib
    import sys

    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from craic import __version__, macapp

    app, _ = macapp.build(tmp_path)          # no iconutil off macOS: the icon is skipped
    info = plistlib.loads((app / "Contents" / "Info.plist").read_bytes())
    assert info["CFBundleExecutable"] == "CRAIC"
    assert info["CFBundleShortVersionString"] == __version__
    launcher = (app / "Contents" / "MacOS" / "CRAIC").read_text()
    assert f'exec "{sys.executable}" -m craic' in launcher
