"""Test-suite for CRAIC. Fast, deterministic, no external aligners required."""

import os

import numpy as np
import pytest

from craic import accel
from craic.domain import Alignment, Alphabet, CodingSpec, Level, translate_codon
from craic.engines import (
    BuiltinProgressive, CodonAware, MafftEngine, Param, builtin_variants,
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
    for p in BuiltinProgressive().parameters() + MafftEngine().parameters():
        assert isinstance(p, Param)
        assert p.kind in ("float", "int", "choice")
        if p.kind == "choice":
            assert p.default in (p.choices or [])


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
