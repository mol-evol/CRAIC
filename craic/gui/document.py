"""The document: what the user is working on, separate from how it is shown.

The window used to own this state directly, mixed in with a hundred and eighty
widget attributes. That is not merely untidy — two bugs came out of it:

* Background analyses had no way to tell that the alignment had moved on, so a
  reliability report computed for one alignment could be painted onto another.
  The fix is a generation counter, and a counter belongs to the thing that
  changes, not to the window that draws it.
* "Everything that replaces the document" was not a concept anywhere in the
  code, so the guard against discarding unsaved work was added to one of the
  four places that needed it and missed the other three.

So the document is an object. It owns the alignment, the reference it is being
scored against and the resulting accuracy, the groups and annotations, where it
came from, whether it has unsaved changes, and a generation that increases
whenever the alignment is replaced. It emits ``changed`` when any of that moves,
and it knows how to become a :class:`craic.session.Session` and back.

It deliberately knows nothing about tracks, colours, dialogs or Qt widgets.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from PySide6.QtCore import QObject, Signal

from .. import evaluate, session as session_mod
from ..domain import Alignment, IdMismatch


class Document(QObject):
    """The alignment being worked on, and everything that travels with it."""

    #: the alignment was replaced (a realignment, an edit, a new file)
    changed = Signal()
    #: the unsaved-changes state moved
    dirty_changed = Signal(bool)
    #: the reference or the accuracy against it moved
    truth_changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._alignment: Optional[Alignment] = None
        self._reference: Optional[Alignment] = None
        self.reference_label: str = ""
        self.truth: Optional[evaluate.Accuracy] = None
        #: The reference's rows on *this* alignment's sequence order, the grid of
        #: true columns behind every cell, and each residue's own correctness.
        #: Recomputed with :attr:`truth` so that the per-residue view and the
        #: column inspector move with an edit exactly as the track does.
        self.truth_rows: Optional[List[str]] = None
        self.truth_grid = None
        self.truth_cells = None
        self.groups: Dict[str, List[str]] = {}
        self.annotations: List[dict] = []
        self.path: str = ""                 # the alignment file, for Save
        self.session_path: str = ""         # the session document
        self._dirty = False
        #: Bumped on every alignment change. Background work captures it when it
        #: starts and its result is discarded if it has moved on since.
        self.generation = 0

    # -- the alignment ------------------------------------------------------ #
    @property
    def alignment(self) -> Optional[Alignment]:
        return self._alignment

    def set_alignment(self, aln: Alignment, dirty: bool = True) -> None:
        """Adopt a new alignment.

        ``dirty`` is False only when the document is being *replaced* rather than
        changed — opening a file, restoring a session — since those are not
        unsaved work.
        """
        self._alignment = aln
        self.generation += 1
        self.rescore()
        self.changed.emit()
        self.set_dirty(dirty)

    @property
    def is_empty(self) -> bool:
        return self._alignment is None

    # -- unsaved changes ---------------------------------------------------- #
    @property
    def dirty(self) -> bool:
        return self._dirty

    def set_dirty(self, dirty: bool = True) -> None:
        dirty = bool(dirty)
        self._dirty = dirty
        self.dirty_changed.emit(dirty)

    # -- the reference, and accuracy against it ----------------------------- #
    @property
    def reference(self) -> Optional[Alignment]:
        return self._reference

    def set_reference(self, reference: Optional[Alignment], label: str = "") -> None:
        """Attach a trusted alignment. Raises :class:`IdMismatch` if it does not
        describe the same sequences, rather than scoring it as a bad alignment."""
        if reference is not None and self._alignment is not None:
            evaluate.compare_to_reference(self._alignment, reference)   # may raise
        self._reference = reference
        self.reference_label = label
        self.rescore()

    def rescore(self) -> None:
        """Recompute accuracy against the reference, if there is one.

        Called on every alignment change, so that realigning or editing updates
        what truth mode shows rather than leaving a stale answer on screen. The
        reference survives a realignment on purpose; it stops applying only when
        the document becomes a different set of sequences, and then it is dropped
        quietly, because that means the user opened something else.
        """
        if self._reference is None or self._alignment is None:
            self._forget_truth()
        else:
            try:
                self.truth = evaluate.compare_to_reference(self._alignment, self._reference)
                self.truth_rows = evaluate.matched_rows(self._alignment, self._reference)
                self.truth_grid = evaluate.true_column_grid(self._alignment.rows,
                                                            self.truth_rows)
                self.truth_cells = evaluate.cell_correct_fraction(
                    self._alignment.rows, self.truth_rows, grid=self.truth_grid)
            except IdMismatch:
                self._reference = None
                self.reference_label = ""
                self._forget_truth()
        self.truth_changed.emit()

    def _forget_truth(self) -> None:
        self.truth = self.truth_rows = self.truth_grid = self.truth_cells = None

    def explain_column(self, col: int) -> Optional[evaluate.ColumnTruth]:
        """What the reference says about one column, or ``None`` without a reference."""
        if self.truth_rows is None or self._alignment is None:
            return None
        if not 0 <= col < self._alignment.length:
            return None
        return evaluate.explain_column(self._alignment.rows, self.truth_rows, col,
                                       ids=self._alignment.ids, grid=self.truth_grid)

    def clear_reference(self) -> None:
        self._reference = None
        self.reference_label = ""
        self._forget_truth()
        self.truth_changed.emit()

    # -- sessions ----------------------------------------------------------- #
    def to_session(self, view: Optional[dict] = None,
                   autosaved: bool = False) -> session_mod.Session:
        return session_mod.Session(
            alignment=self._alignment,
            reference=self._reference,
            reference_label=self.reference_label,
            annotations=list(self.annotations),
            groups={k: list(v) for k, v in self.groups.items()},
            source_path=self.path,
            view=dict(view or {}),
            autosaved=autosaved,
        )

    def restore_context(self, sess: session_mod.Session, session_path: str = "") -> None:
        """Restore everything in a session *except* the alignment.

        Split out because the window adopts the alignment through its own path
        (which also refreshes widgets); this puts the rest in place first, so
        that when the alignment lands the reference is already there to score
        against.
        """
        self.groups = {k: list(v) for k, v in sess.groups.items()}
        self.annotations = list(sess.annotations)
        self._reference = sess.reference
        self.reference_label = sess.reference_label
        self.truth = None
        self.path = sess.source_path
        self.session_path = session_path

    def load_session(self, sess: session_mod.Session, session_path: str = "") -> None:
        """Restore a whole session, alignment included."""
        self.restore_context(sess, session_path)
        self.set_alignment(sess.alignment, dirty=False)

    # -- description -------------------------------------------------------- #
    def display_name(self) -> str:
        import os

        for candidate in (self.path, self.session_path):
            if candidate:
                return os.path.basename(candidate)
        return ""
