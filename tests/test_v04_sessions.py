"""Tests for v0.4: the session document, unsaved-change tracking and the
crash-recovery autosave.

The thing being protected here is a curation session — twenty minutes of hand
judgement that nothing else can reconstruct. Each test below corresponds to a
way that work could be lost.
"""

import json
import os

import pytest

from craic import session as session_mod
from craic.domain import Alignment, Alphabet, CodingSpec, IdMismatch
from craic.session import Session, SessionError


def _aln(rows=("ACG-T", "ACGTT"), ids=("a", "b"), **kw):
    return Alignment(list(ids), list(rows), Alphabet.DNA, **kw)


# --------------------------------------------------------------------------- #
# The session document
# --------------------------------------------------------------------------- #

def test_a_session_round_trips_everything_that_matters(tmp_path):
    aln = _aln(coding=CodingSpec(frame=0, table=1),
               meta={"history": ["opened x", "aligned"]})
    sess = Session(alignment=aln,
                   reference=_aln(rows=("ACGT-", "ACGTT")),
                   reference_label="simulated truth",
                   annotations=[{"name": "loop", "start": 1, "stop": 3}],
                   groups={"clade A": ["a"]},
                   source_path="/data/x.fasta",
                   view={"track": 11, "threshold": 55})
    path = str(tmp_path / ("x.fasta" + session_mod.SUFFIX))
    sess.save(path)
    back = Session.load(path)

    assert back.alignment.rows == aln.rows
    assert back.alignment.ids == aln.ids
    assert back.alignment.coding == CodingSpec(frame=0, table=1)
    assert back.alignment.meta["history"] == ["opened x", "aligned"]
    assert back.reference.rows == ["ACGT-", "ACGTT"]
    assert back.reference_label == "simulated truth"
    assert back.annotations == [{"name": "loop", "start": 1, "stop": 3}]
    assert back.groups == {"clade A": ["a"]}
    assert back.source_path == "/data/x.fasta"
    assert back.view["track"] == 11 and back.view["threshold"] == 55


def test_unserialisable_metadata_is_dropped_not_fatal(tmp_path):
    """Derived state can always be recomputed; refusing to save because of it
    would turn a cosmetic problem into lost work."""
    import numpy as np

    aln = _aln(meta={"history": ["a"], "scores": np.arange(3)})
    path = str(tmp_path / "s.craic.json")
    Session(alignment=aln).save(path)
    back = Session.load(path)
    assert back.alignment.meta["history"] == ["a"]
    assert "scores" not in back.alignment.meta


def test_a_session_is_readable_json(tmp_path):
    path = str(tmp_path / "s.craic.json")
    Session(alignment=_aln()).save(path)
    data = json.loads(open(path).read())
    assert data["format"] == "craic-session"
    assert data["format_version"] == session_mod.FORMAT_VERSION
    assert data["alignment"]["rows"] == ["ACG-T", "ACGTT"]


def test_a_foreign_json_file_is_refused(tmp_path):
    path = tmp_path / "other.json"
    path.write_text('{"hello": 1}')
    with pytest.raises(SessionError, match="not a CRAIC session"):
        Session.load(str(path))


def test_a_newer_format_is_refused_rather_than_misread(tmp_path):
    path = tmp_path / "future.craic.json"
    path.write_text(json.dumps({"format": "craic-session",
                                "format_version": session_mod.FORMAT_VERSION + 5}))
    with pytest.raises(SessionError, match="newer CRAIC"):
        Session.load(str(path))


def test_saving_is_atomic_and_leaves_no_partial_file(tmp_path):
    """The autosave rewrites this file constantly, and a crash mid-write is
    exactly when it matters, so a half-written recovery file must be impossible."""
    path = str(tmp_path / "s.craic.json")
    Session(alignment=_aln()).save(path)
    Session(alignment=_aln(rows=("ACGTT", "ACGTT"))).save(path)
    assert sorted(os.listdir(tmp_path)) == ["s.craic.json"]
    assert Session.load(path).alignment.rows == ["ACGTT", "ACGTT"]


def test_is_session_file_recognises_one_without_parsing_it_all(tmp_path):
    good = tmp_path / ("a.fasta" + session_mod.SUFFIX)
    Session(alignment=_aln()).save(str(good))
    bad = tmp_path / "b.json"
    bad.write_text('{"not": "ours"}')
    plain = tmp_path / "c.fasta"
    plain.write_text(">a\nACGT\n")
    assert session_mod.is_session_file(str(good))
    assert not session_mod.is_session_file(str(bad))
    assert not session_mod.is_session_file(str(plain))


# --------------------------------------------------------------------------- #
# Unsaved changes in the GUI
# --------------------------------------------------------------------------- #

def _win():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.app import CraicWindow

    QApplication.instance() or QApplication([])
    win = CraicWindow()
    win._autosave_timer.stop()          # tests drive the autosave explicitly
    return win


def test_editing_marks_the_document_dirty_and_the_title_shows_it():
    from craic import editing

    win = _win()
    win._set_alignment(_aln(rows=("ACGTAC", "ACG-AC", "ACGTAC"),
                            ids=("a", "b", "c")), dirty=False)
    assert win._dirty is False
    assert not win.windowTitle().startswith("•")

    new = editing.insert_gap_column(win.aln, 2)
    win._apply_edit(new, "insert gap column")
    assert win._dirty is True
    assert win.windowTitle().startswith("•")


def test_closing_with_unsaved_work_asks_first(monkeypatch, tmp_path):
    from PySide6.QtGui import QCloseEvent
    from PySide6.QtWidgets import QMessageBox

    from craic.gui import app as app_mod

    win = _win()
    win._set_alignment(_aln(), dirty=True)

    asked = {"n": 0}

    def fake_warning(*_a, **_kw):
        asked["n"] += 1
        return QMessageBox.Cancel

    monkeypatch.setattr(app_mod.QMessageBox, "warning", staticmethod(fake_warning))
    ev = QCloseEvent()
    win.closeEvent(ev)
    assert asked["n"] == 1
    assert not ev.isAccepted()              # Cancel means stay open


def test_discarding_closes_and_saving_writes_a_session(monkeypatch, tmp_path):
    from PySide6.QtGui import QCloseEvent
    from PySide6.QtWidgets import QMessageBox

    from craic.gui import app as app_mod

    # Discard
    win = _win()
    win._set_alignment(_aln(), dirty=True)
    monkeypatch.setattr(app_mod.QMessageBox, "warning",
                        staticmethod(lambda *a, **k: QMessageBox.Discard))
    ev = QCloseEvent()
    win.closeEvent(ev)
    assert ev.isAccepted()

    # Save
    win = _win()
    win._set_alignment(_aln(), dirty=True)
    target = str(tmp_path / "saved.craic.json")
    win._session_path = target
    monkeypatch.setattr(app_mod.QMessageBox, "warning",
                        staticmethod(lambda *a, **k: QMessageBox.Save))
    ev = QCloseEvent()
    win.closeEvent(ev)
    assert ev.isAccepted()
    assert Session.load(target).alignment.rows == ["ACG-T", "ACGTT"]


def test_a_clean_document_closes_without_a_prompt(monkeypatch):
    from PySide6.QtGui import QCloseEvent

    from craic.gui import app as app_mod

    win = _win()
    win._set_alignment(_aln(), dirty=False)
    monkeypatch.setattr(app_mod.QMessageBox, "warning",
                        staticmethod(lambda *a, **k: pytest.fail("should not ask")))
    ev = QCloseEvent()
    win.closeEvent(ev)
    assert ev.isAccepted()


def test_saving_a_session_makes_the_document_clean(tmp_path):
    win = _win()
    win._set_alignment(_aln(), dirty=True)
    win._session_path = str(tmp_path / "s.craic.json")
    assert win._save_session() is True
    assert win._dirty is False
    assert not win.windowTitle().startswith("•")


# --------------------------------------------------------------------------- #
# Crash recovery
# --------------------------------------------------------------------------- #

def _redirect_autosave(win, tmp_path):
    path = str(tmp_path / "recovery.craic.json")
    win._autosave_path = lambda: path
    return path


def test_autosave_writes_the_session_and_a_clean_close_removes_it(monkeypatch, tmp_path):
    from PySide6.QtGui import QCloseEvent

    win = _win()
    path = _redirect_autosave(win, tmp_path)
    win._set_alignment(_aln(), dirty=True)
    win._write_autosave()
    assert os.path.exists(path)
    assert Session.load(path).autosaved is True

    # A clean exit removes it, so the presence of one always means a crash.
    win._set_dirty(False)
    win.closeEvent(QCloseEvent())
    assert not os.path.exists(path)


def test_autosave_does_not_write_beside_the_users_data(tmp_path):
    """Opening someone's reference alignment and looking at it must not leave
    files in their directory."""
    win = _win()
    data_dir = tmp_path / "their_data"
    data_dir.mkdir()
    (data_dir / "ref.fasta").write_text(">a\nACGT\n>b\nACGT\n")
    win._path = str(data_dir / "ref.fasta")
    win._set_alignment(_aln(), dirty=True)
    win._write_autosave()
    assert sorted(os.listdir(data_dir)) == ["ref.fasta"]
    assert "craic" in win._autosave_path()


def test_recovery_offers_the_interrupted_session_and_restores_it(monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox

    from craic.gui import app as app_mod

    # a "crashed" session on disk
    path = str(tmp_path / "recovery.craic.json")
    Session(alignment=_aln(rows=("AAAA", "ACGT"), ids=("x", "y")),
            reference=_aln(rows=("AAAA", "ACGT"), ids=("x", "y")),
            reference_label="truth", source_path="/data/x.fasta",
            autosaved=True).save(path)

    win = _win()
    win._autosave_path = lambda: path
    monkeypatch.setattr(app_mod.QMessageBox, "question",
                        staticmethod(lambda *a, **k: QMessageBox.Yes))
    assert win.maybe_recover() is True
    assert win.aln.rows == ["AAAA", "ACGT"]
    assert win.reference is not None and win._reference_label == "truth"
    # recovered work is not yet saved anywhere, so it stays dirty
    assert win._dirty is True


def test_declining_recovery_removes_the_file(monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox

    from craic.gui import app as app_mod

    path = str(tmp_path / "recovery.craic.json")
    Session(alignment=_aln(), autosaved=True).save(path)
    win = _win()
    win._autosave_path = lambda: path
    monkeypatch.setattr(app_mod.QMessageBox, "question",
                        staticmethod(lambda *a, **k: QMessageBox.No))
    assert win.maybe_recover() is False
    assert not os.path.exists(path)


def test_no_recovery_file_means_no_prompt(tmp_path):
    win = _win()
    win._autosave_path = lambda: str(tmp_path / "nothing.craic.json")
    assert win.maybe_recover() is False


def test_a_corrupt_recovery_file_is_discarded_quietly(tmp_path):
    path = tmp_path / "recovery.craic.json"
    path.write_text("{ this is not json")
    win = _win()
    win._autosave_path = lambda: str(path)
    assert win.maybe_recover() is False
    assert not path.exists()


# --------------------------------------------------------------------------- #
# Opening and round-tripping through the window
# --------------------------------------------------------------------------- #

def test_open_path_recognises_a_session(tmp_path):
    win = _win()
    path = str(tmp_path / ("x.fasta" + session_mod.SUFFIX))
    Session(alignment=_aln(rows=("GGGG", "GGGA"), ids=("p", "q")),
            source_path="/data/x.fasta").save(path)
    win.open_path(path)
    assert win.aln.ids == ["p", "q"]
    assert win._session_path == path
    assert win._dirty is False


def test_a_full_window_round_trip_keeps_the_working_state(tmp_path):
    from craic.gui.app import _TRUTH_TRACK

    win = _win()
    aln = _aln(rows=("ACGTAC", "ACG-AC"), ids=("a", "b"))
    win._set_alignment(aln, dirty=False)
    win._set_reference(_aln(rows=("ACGTAC", "ACG-AC"), ids=("a", "b")), "truth")
    win._groups = {"mine": ["a"]}
    win._annotations = [{"name": "site", "start": 0, "stop": 2}]
    win.thr.setValue(72)

    path = str(tmp_path / "round.craic.json")
    win._session_path = path
    assert win._save_session() is True

    other = _win()
    other.open_session_path(path)
    assert other.aln.rows == aln.rows
    assert other.reference is not None
    assert other._reference_label == "truth"
    assert other._groups == {"mine": ["a"]}
    assert other._annotations == [{"name": "site", "start": 0, "stop": 2}]
    assert other.thr.value() == 72
    assert other.track_combo.currentIndex() == _TRUTH_TRACK
    assert other._dirty is False


def test_an_expensive_track_is_not_resumed_automatically(tmp_path):
    """Restoring a session should not silently start a minutes-long analysis."""
    win = _win()
    win._set_alignment(_aln(), dirty=False)
    path = str(tmp_path / "s.craic.json")
    sess = win._current_session()
    sess.view["track"] = 2                       # Reliability
    sess.save(path)

    other = _win()
    other.open_session_path(path)
    assert other.track_combo.currentIndex() == 0


# --------------------------------------------------------------------------- #
# Quitting is not the only way to lose a session
# --------------------------------------------------------------------------- #

def test_generating_a_dataset_asks_before_discarding_unsaved_work(monkeypatch):
    """Generate replaces the document as completely as quitting does, and used
    to do it silently — which is how an aligned dataset disappeared."""
    from PySide6.QtWidgets import QMessageBox

    from craic.gui import app as app_mod

    win = _win()
    win._set_alignment(_aln(rows=("ACGTAC", "ACG-AC"), ids=("a", "b")), dirty=True)
    before = win.aln.rows

    asked = {"n": 0}

    def cancel(*_a, **_k):
        asked["n"] += 1
        return QMessageBox.Cancel

    monkeypatch.setattr(app_mod.QMessageBox, "warning", staticmethod(cancel))
    monkeypatch.setattr(app_mod, "SimulateDialog",
                        lambda *a, **k: pytest.fail("dialog opened despite Cancel"))
    win._generate_dataset()
    assert asked["n"] == 1
    assert win.aln.rows == before          # the work is still there


def test_opening_another_file_asks_first(monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox

    from craic.gui import app as app_mod

    win = _win()
    win._set_alignment(_aln(), dirty=True)
    monkeypatch.setattr(app_mod.QMessageBox, "warning",
                        staticmethod(lambda *a, **k: QMessageBox.Cancel))
    monkeypatch.setattr(app_mod.QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: pytest.fail("asked for a file")))
    win._open()
    win._open_session()


def test_a_clean_document_is_replaced_without_a_prompt(monkeypatch):
    from craic.gui import app as app_mod

    win = _win()
    win._set_alignment(_aln(), dirty=False)
    monkeypatch.setattr(app_mod.QMessageBox, "warning",
                        staticmethod(lambda *a, **k: pytest.fail("should not ask")))
    assert win._confirm_discard("Doing something") is True


def test_realigning_can_be_undone():
    """Realigning throws away every hand edit since the last alignment. That is
    as destructive as any single edit, and was the one with no way back."""
    win = _win()
    first = _aln(rows=("ACGTAC", "ACG-AC"), ids=("a", "b"))
    win._set_alignment(first, dirty=False)
    win._commit(_aln(rows=("ACGTAC", "A-CGAC"), ids=("a", "b")), "aligned · MAFFT")
    assert win.aln.rows == ["ACGTAC", "A-CGAC"]
    win._undo()
    assert win.aln.rows == first.rows


# --------------------------------------------------------------------------- #
# The document model
# --------------------------------------------------------------------------- #

def test_the_document_owns_the_state_and_bumps_a_generation():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.document import Document

    QApplication.instance() or QApplication([])
    doc = Document()
    assert doc.is_empty and doc.generation == 0 and doc.dirty is False

    doc.set_alignment(_aln())
    assert doc.generation == 1 and doc.dirty is True
    doc.set_alignment(_aln(rows=("ACGTT", "ACGTT")), dirty=False)
    assert doc.generation == 2 and doc.dirty is False


def test_the_document_rescores_when_the_alignment_changes():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.document import Document

    QApplication.instance() or QApplication([])
    doc = Document()
    aln = _aln(rows=("ACGTAC", "ACG-AC"), ids=("a", "b"))
    doc.set_alignment(aln, dirty=False)
    doc.set_reference(_aln(rows=("ACGTAC", "ACG-AC"), ids=("a", "b")), "truth")
    assert doc.truth is not None and doc.truth.sp == 1.0

    doc.set_alignment(_aln(rows=("ACGTAC-", "ACG-AC-"), ids=("a", "b")))
    assert doc.truth is not None                 # the reference survives a change


def test_the_document_drops_a_reference_for_other_sequences():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.document import Document

    QApplication.instance() or QApplication([])
    doc = Document()
    doc.set_alignment(_aln(ids=("a", "b")), dirty=False)
    doc.set_reference(_aln(ids=("a", "b")), "truth")
    doc.set_alignment(_aln(ids=("x", "y")))      # different sequences entirely
    assert doc.reference is None and doc.truth is None


def test_a_reference_for_other_sequences_is_refused_up_front():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.document import Document

    QApplication.instance() or QApplication([])
    doc = Document()
    doc.set_alignment(_aln(ids=("a", "b")), dirty=False)
    with pytest.raises(IdMismatch):
        doc.set_reference(_aln(ids=("x", "y")), "wrong")


def test_the_window_exposes_the_document_under_the_old_names():
    """The window's aln / reference / truth / _dirty are the document's, so the
    refactor did not require every caller to change at once."""
    win = _win()
    win._set_alignment(_aln(), dirty=False)
    assert win.aln is win.doc.alignment
    assert win._dirty is win.doc.dirty
    assert win._doc_generation == win.doc.generation


# --------------------------------------------------------------------------- #
# The track registry
# --------------------------------------------------------------------------- #

def test_tracks_are_described_once_and_the_cheap_set_follows():
    from craic.ambiguity import trimming
    from craic.gui import tracks as track_mod

    defs = track_mod.build(lambda aln: None, trimming)
    labels = [t.label for t in defs]
    assert len(set(t.key for t in defs)) == len(defs)      # keys are unique
    assert labels[0] == "(none)"
    # the expensive ones are exactly the background analyses
    expensive = [t.key for t in defs if not t.cheap and t.kind != "none"]
    assert set(expensive) == {"reliability", "consistency", "perturbation", "agreement"}
    assert track_mod.index_of(defs, "truth") == len(defs) - 1


def test_adding_a_track_does_not_need_three_edits():
    """The regression guard for the bug this replaced: the combo box, the
    cheap-track rule and the dispatch all read from one list."""
    from craic.gui import app as app_mod

    assert app_mod._TRACKS == [t.label for t in app_mod._TRACK_DEFS]
    assert app_mod._TRUTH_TRACK == app_mod._TRACK_DEFS.index(
        next(t for t in app_mod._TRACK_DEFS if t.key == "truth"))
    assert app_mod._CHEAP_TRACKS == frozenset(
        i for i, t in enumerate(app_mod._TRACK_DEFS) if t.cheap and i != 0)
