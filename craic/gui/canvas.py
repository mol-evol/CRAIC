"""The alignment canvas: a custom-painted, virtualized alignment viewer.

Only the visible cells are drawn (viewport culling), so it stays responsive on
large alignments. A display column's pixel width scales with its span, so a codon
column is three nucleotide columns wide and stays readable (and the horizontal
scale is identical across nt / codon / aa views). A fixed name gutter, a
per-column score track, a ruler, and an optional consensus row sit around the
grid. The same drawing routine renders both the viewport and a full image export.
"""

from __future__ import annotations

from typing import List, Optional, Set, Tuple

import numpy as np
from PySide6.QtCore import QRect, Qt, QTimer, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QFontMetrics, QImage, QPainter
from PySide6.QtWidgets import QAbstractScrollArea

from ..domain import Alignment, Level, column_consensus
from . import colors


class AlignmentCanvas(QAbstractScrollArea):
    selectionChanged = Signal(int, int)      # nucleotide coords of the selection
    residueClicked = Signal(int, int)        # (sequence index, nucleotide column)
    editRequested = Signal(str)              # keyboard edit op (see keyPressEvent)
    cursorMoved = Signal(int, int)           # edit cursor (row, display column)
    rowsSelected = Signal(int)               # how many sequences are selected for editing
    selectionCleared = Signal()              # column + row selection was cleared (Esc)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.aln: Optional[Alignment] = None
        self.level = Level.NT
        self._dm = None

        self.cell_w = 15           # width of one *nucleotide* column
        self.cell_h = 17
        self.name_w = 150
        self.track_h = 22
        self.ruler_h = 16

        self._nt_scores: Optional[np.ndarray] = None
        self._score_label = ""
        self._nt_keep: Optional[np.ndarray] = None
        #: Nucleotide columns the user has pinned as trusted (see set_pins).
        self._pins: set = set()
        self._masked: Set[Tuple[int, int]] = set()      # masked residues, (row, nt column)
        self._sel_cols: Optional[Tuple[int, int]] = None
        self._anchor: Optional[int] = None

        self.color_mode = "residue"          # or "confidence", or "truth"
        self.scheme = "clustal"
        self._cell_scores: Optional[np.ndarray] = None
        self._probe: Optional[np.ndarray] = None
        self._probe_cell: Optional[Tuple[int, int]] = None
        self._cursor: Optional[Tuple[int, int]] = None   # edit cursor (row, display col)
        self._sel_rows: Set[int] = set()                 # selected sequences (for editing)
        self._row_anchor: Optional[int] = None           # for Shift-click range select
        self._press: Optional[Tuple[Optional[int], Optional[int]]] = None
        self._press_x = 0
        self._press_mod = Qt.NoModifier
        self._moved = False
        #: While true, nothing in the canvas can start an edit — the keyboard
        #: slide keys and the drag-to-slide are both refused at source. Set when
        #: what is on screen is not the document being worked on (the true
        #: alignment, shown for reference), where an edit would otherwise be
        #: applied to the working alignment's coordinates behind the user's back.
        self.read_only = False
        self._slide = False          # mouse-drag sliding the selected sequences
        self._slide_col = 0          # last column the slide-drag has reached
        # Auto-scroll while drag-selecting past the viewport edge.
        self._drag_x = 0
        self._autoscroll_dir = 0
        self._autoscroll = QTimer(self)
        self._autoscroll.setInterval(30)
        self._autoscroll.timeout.connect(self._autoscroll_tick)

        self.show_consensus = False
        self._consensus: Optional[List[str]] = None
        self._annotations: List[dict] = []   # column features: {name, start, stop(nt), color}
        self._find: Set[Tuple[int, int]] = set()
        self._find_row: Optional[int] = None
        self._outlier_rows: Set[int] = set()             # EvalMSA-style flagged sequences

        f = QFont("Menlo")
        f.setStyleHint(QFont.TypeWriter)
        f.setPixelSize(11)
        self.setFont(f)
        self._fm = QFontMetrics(f)
        self._letter = QFont(f)          # residue font; rescaled with the cell size
        self._rescale_font()

        self.viewport().setMouseTracking(True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.setFocusPolicy(Qt.StrongFocus)      # so the canvas receives edit keys

    # -- data -------------------------------------------------------------- #
    # -- public read-only accessors -------------------------------------- #
    # Controllers read these instead of reaching into the canvas's privates,
    # keeping the view/controller boundary explicit.
    @property
    def dm(self):
        """Current DisplayMatrix (cells + nt-coordinate map) being rendered, or None."""
        return self._dm

    @property
    def nt_scores(self) -> Optional[np.ndarray]:
        """The active track's per-nucleotide scores, or None."""
        return self._nt_scores

    @property
    def score_label(self) -> str:
        """Human-readable label of the active track (for status text)."""
        return self._score_label

    @property
    def keep_mask(self) -> Optional[np.ndarray]:
        """Current per-column keep mask (True = retained), or None."""
        return self._nt_keep

    def set_alignment(self, aln: Alignment, level: Optional[Level] = None) -> None:
        self.aln = aln
        levels = aln.available_levels()
        # Prefer the requested level, else keep the current one (so an edit while
        # viewing protein/codon doesn't snap the view back to nucleotide), else default.
        if level in levels:
            self.level = level
        elif self.level not in levels:
            self.level = levels[0]
        self._dm = aln.display(self.level)
        self._nt_scores = self._nt_keep = self._sel_cols = None
        self._pins = set()
        self._masked = set()
        self._outlier_rows = set()
        self._cell_scores = self._probe = self._probe_cell = None
        self._cursor = None
        self._sel_rows = set()
        self._row_anchor = None
        self._consensus = None
        self._find = set()
        self._find_row = None
        if self.show_consensus:
            self._compute_consensus()
        self._update_scrollbars()
        self.viewport().update()

    def set_level(self, level: Level) -> None:
        if not self.aln or level not in self.aln.available_levels():
            return
        self.level = level
        self._dm = self.aln.display(level)
        self._consensus = None
        self._find = set()
        self._find_row = None
        if self.show_consensus:
            self._compute_consensus()
        self._update_scrollbars()
        self.viewport().update()

    def set_column_scores(self, nt_scores: Optional[np.ndarray], label: str = "") -> None:
        self._nt_scores = nt_scores
        self._score_label = label
        self.viewport().update()

    def set_keep_mask(self, nt_keep: Optional[np.ndarray]) -> None:
        self._nt_keep = nt_keep
        self.viewport().update()

    def set_pins(self, nt_cols) -> None:
        """Mark nucleotide columns the user has pinned as trusted.

        Drawn as a solid band above the grid rather than a tint over it, so that
        a pin stays visible whatever colour mode the residues are in — pins are
        the user's own assertion and must not be confused with a computed score.
        """
        self._pins = set(int(c) for c in nt_cols)
        self.viewport().update()

    @property
    def pins(self) -> set:
        return set(self._pins)

    def set_residue_mask(self, cells) -> None:
        """Residues the user has masked, as ``(row, nucleotide column)`` cells.

        Drawn greyed and hatched: still visible, since the column is kept, but
        plainly not part of what will be exported.
        """
        self._masked = set(cells)
        self.viewport().update()

    def _cell_masked(self, r: int, c: int) -> bool:
        if not self._masked:
            return False
        s = self._dm.nt_start[c]
        return any((r, s + k) in self._masked for k in range(self._dm.span))

    def _col_pinned(self, c: int) -> bool:
        if not self._pins or self._dm is None:
            return False
        s = self._dm.nt_start[c]
        return any((s + k) in self._pins for k in range(self._dm.span))

    def set_outlier_rows(self, rows) -> None:
        """Highlight whole sequences flagged as poorly-aligned outliers."""
        self._outlier_rows = set(rows)
        self.viewport().update()

    def set_read_only(self, on: bool) -> None:
        self.read_only = bool(on)
        self._slide = False

    def set_color_mode(self, mode: str) -> None:
        self.color_mode = mode
        self.viewport().update()

    def set_scheme(self, scheme: str) -> None:
        self.scheme = scheme
        self.viewport().update()

    def set_cell_scores(self, cell_scores: Optional[np.ndarray]) -> None:
        """Per-residue scores to colour by, shaped (n_seqs, n_columns).

        The shape is checked here rather than trusted: paintEvent indexes this
        array directly, so an array left over from a previous alignment shows up
        as wrong colours or an IndexError deep in the paint loop. A mismatch is
        a bug upstream, so it is refused loudly at the boundary.
        """
        if cell_scores is not None and getattr(self, "aln", None) is not None:
            expected = (self.aln.n_seqs, self.aln.length)
            actual = tuple(np.shape(cell_scores))
            if actual != expected:
                raise ValueError(
                    f"cell scores are shaped {actual}, but the alignment is {expected}"
                )
        self._cell_scores = cell_scores
        self.viewport().update()

    def set_probe(self, nt_probs: Optional[np.ndarray], seq: int = -1, nt_col: int = -1) -> None:
        self._probe = nt_probs
        self._probe_cell = (seq, nt_col) if nt_probs is not None else None
        self.viewport().update()

    def set_consensus(self, show: bool) -> None:
        self.show_consensus = bool(show)
        if self.show_consensus and self._consensus is None and self._dm:
            self._compute_consensus()
        self._update_scrollbars()
        self.viewport().update()

    def scroll_to_nt(self, nt_col: int) -> None:
        """Put the edit cursor on the display column holding ``nt_col``.

        Callers that work in alignment coordinates — the column inspector, which
        reports columns of the canonical model — should not have to know that a
        codon or amino-acid view packs three of them into one display column.
        """
        rng = self._nt_to_cols(nt_col, nt_col + 1)
        if rng is None:
            return
        cur = self.cursor()
        self.set_cursor(cur[0] if cur else 0, rng[0])

    def scroll_to_col(self, c: int) -> None:
        ox, _ = self._grid_origin()
        target = c * self._col_w() - max(0, (self.viewport().width() - ox)) // 2
        self.horizontalScrollBar().setValue(max(0, int(target)))
        self.viewport().update()

    def selected_records(self) -> List[Tuple[str, str]]:
        """The visible cells (current level), sliced to the column selection if
        any, as ``(id, sequence)`` records — the basis for Copy in any format."""
        if not self.aln or not self._dm:
            return []
        c0, c1 = self._sel_cols if self._sel_cols else (0, self._dm.n_cols)
        return [
            (self.aln.ids[r], "".join(self._dm.cells[r][c] for c in range(c0, c1)))
            for r in range(self.aln.n_seqs)
        ]

    def selection_nt(self) -> Optional[Tuple[int, int]]:
        if not self._sel_cols or not self._dm:
            return None
        c0, c1 = self._sel_cols
        return self._dm.nt_start[c0], self._dm.nt_start[c1 - 1] + self._dm.span

    def search(self, query: str) -> bool:
        """Find a sequence by name, or a residue/codon motif. Highlights and
        scrolls to the first hit; returns True if anything matched."""
        self._find = set()
        self._find_row = None
        q = (query or "").upper().strip()
        if not q or not self.aln or not self._dm:
            self.viewport().update()
            return False
        for i, name in enumerate(self.aln.ids):
            if q in name.upper():
                self._find_row = i
                self._ensure_row_visible(i)
                self.viewport().update()
                return True
        first = None
        for r in range(self.aln.n_seqs):
            s = ""
            char_to_col: List[int] = []
            for c, cell in enumerate(self._dm.cells[r]):
                if cell.strip("-"):
                    up = cell.upper()
                    s += up
                    char_to_col += [c] * len(up)
            idx = s.find(q)
            while idx != -1:
                for k in range(idx, idx + len(q)):
                    self._find.add((r, char_to_col[k]))
                    if first is None:
                        first = (r, char_to_col[k])
                idx = s.find(q, idx + 1)
        if first:
            self._ensure_row_visible(first[0])
            self.scroll_to_col(first[1])
        self.viewport().update()
        return bool(self._find)

    def render_full(self) -> QImage:
        """Render the entire alignment at the current cell size to an image."""
        if not self.aln or not self._dm:
            return QImage()
        ox, oy = self._grid_origin()
        w = min(20000, ox + self._dm.n_cols * self._col_w())
        h = min(20000, oy + self.aln.n_seqs * self.cell_h
                + self._consensus_h() + self._annot_h())
        img = QImage(int(w), int(h), QImage.Format_RGB32)
        img.fill(colors.BG)
        p = QPainter(img)
        p.setFont(self.font())
        self._paint(p, int(w), int(h), 0, 0)
        p.end()
        return img

    def _rescale_font(self) -> None:
        """Size the residue font to the current cell, so zooming changes the text."""
        size = max(6, min(self.cell_h - 3, int(self.cell_w * 1.45)))
        self._letter.setPixelSize(size)

    def set_cell_size(self, w: int, h: int) -> None:
        """Set the per-nucleotide cell size (and rescale the font) — drives zoom."""
        self.cell_w = int(min(40, max(3, w)))
        self.cell_h = int(min(44, max(4, h)))
        self._rescale_font()
        self._update_scrollbars()
        self.viewport().update()

    # -- geometry ---------------------------------------------------------- #
    def _grid_origin(self) -> Tuple[int, int]:
        return self.name_w, self.track_h + self.ruler_h

    def _consensus_h(self) -> int:
        return self.cell_h if (self.show_consensus and self.aln) else 0

    def _annot_h(self) -> int:
        return max(16, self.cell_h) if (self._annotations and self.aln) else 0

    def set_annotations(self, anns: List[dict]) -> None:
        """Column features to draw in the bottom band: ``{name, start, stop, color}``
        in nucleotide coordinates (fixed by position)."""
        self._annotations = list(anns or [])
        self._update_scrollbars()
        self.viewport().update()

    def _nt_to_cols(self, nt_start: int, nt_stop: int) -> Optional[Tuple[int, int]]:
        """Display-column range (inclusive) covering nucleotide [nt_start, nt_stop)."""
        if not self._dm or nt_stop <= nt_start:
            return None
        span = self._dm.span
        cols = [c for c in range(self._dm.n_cols)
                if self._dm.nt_start[c] < nt_stop and self._dm.nt_start[c] + span > nt_start]
        return (cols[0], cols[-1]) if cols else None

    def _col_w(self) -> int:
        """Pixel width of one display column (3 nucleotides wide for codons/aa)."""
        return self.cell_w * (self._dm.span if self._dm else 1)

    def _update_scrollbars(self) -> None:
        cols = self._dm.n_cols if self._dm else 0
        rows = self.aln.n_seqs if self.aln else 0
        col_w = self._col_w()
        vw, vh = self.viewport().width(), self.viewport().height()
        avail_w = max(0, vw - self.name_w)
        avail_h = max(0, vh - self.track_h - self.ruler_h
                      - self._consensus_h() - self._annot_h())
        hbar, vbar = self.horizontalScrollBar(), self.verticalScrollBar()
        hbar.setRange(0, max(0, cols * col_w - avail_w))
        hbar.setPageStep(avail_w)
        hbar.setSingleStep(col_w)
        vbar.setRange(0, max(0, rows * self.cell_h - avail_h))
        vbar.setPageStep(avail_h)
        vbar.setSingleStep(self.cell_h)

    def resizeEvent(self, e):
        self._update_scrollbars()
        super().resizeEvent(e)

    def scrollContentsBy(self, dx, dy):
        self.viewport().update()

    def _ensure_row_visible(self, r: int) -> None:
        _, oy = self._grid_origin()
        y = r * self.cell_h
        vbar = self.verticalScrollBar()
        avail = max(0, self.viewport().height() - oy - self._consensus_h() - self._annot_h())
        if y < vbar.value():
            vbar.setValue(y)
        elif y + self.cell_h > vbar.value() + avail:
            vbar.setValue(max(0, y - avail + self.cell_h))

    def _ensure_col_visible(self, c: int) -> None:
        ox, _ = self._grid_origin()
        cw = self._col_w()
        x = c * cw
        hbar = self.horizontalScrollBar()
        avail = max(0, self.viewport().width() - ox)
        if x < hbar.value():
            hbar.setValue(x)
        elif x + cw > hbar.value() + avail:
            hbar.setValue(max(0, x - avail + cw))

    def _compute_consensus(self) -> None:
        cells = self._dm.cells
        n = len(cells)
        cons = []
        for c in range(self._dm.n_cols):
            cc = column_consensus([cells[r][c] for r in range(n)], self.level, self.aln.alphabet)
            cons.append(cc if cc else "-" * self._dm.span)
        self._consensus = cons

    def _agg_vec(self, arr: np.ndarray, c: int) -> float:
        s = self._dm.nt_start[c]
        seg = arr[s : s + self._dm.span]
        seg = seg[~np.isnan(seg)] if len(seg) else seg
        return float(np.mean(seg)) if len(seg) else float("nan")

    def _agg_score(self, c: int) -> float:
        return self._agg_vec(self._nt_scores, c) if self._nt_scores is not None else float("nan")

    def _agg_cell(self, r: int, c: int) -> float:
        return self._agg_vec(self._cell_scores[r], c) if self._cell_scores is not None else float("nan")

    def _agg_probe(self, c: int) -> float:
        return self._agg_vec(self._probe, c) if self._probe is not None else float("nan")

    def _col_kept(self, c: int) -> bool:
        if self._nt_keep is None or self._dm is None:
            return True
        s = self._dm.nt_start[c]
        seg = self._nt_keep[s : s + self._dm.span]
        return bool(np.mean(seg) >= 0.5) if len(seg) else True

    # -- painting ---------------------------------------------------------- #
    def paintEvent(self, e):
        p = QPainter(self.viewport())
        p.fillRect(self.viewport().rect(), colors.BG)
        if not self.aln or not self._dm:
            p.setPen(colors.MUTED)
            p.drawText(self.viewport().rect(), Qt.AlignCenter,
                       "Open an alignment or sequences to begin")
            return
        self._paint(p, self.viewport().width(), self.viewport().height(),
                    self.horizontalScrollBar().value(), self.verticalScrollBar().value())

    def _paint(self, p: QPainter, vw: int, vh: int, hoff: int, voff: int) -> None:
        ox, oy = self._grid_origin()
        cons_h = self._consensus_h()
        annot_h = self._annot_h()
        rows_bottom = vh - cons_h - annot_h
        ncols, nrows = self._dm.n_cols, self.aln.n_seqs
        # Column-wide overlays (selection, mask) should stop at the last sequence,
        # not run to the bottom of the viewport when the alignment is short.
        grid_bottom = min(rows_bottom, oy + nrows * self.cell_h - voff)
        table = self.aln.coding.table if self.aln.coding else 1
        scheme = self.scheme
        col_w = self._col_w()
        span = self._dm.span

        first_col = max(0, hoff // col_w)
        last_col = min(ncols, (hoff + (vw - ox)) // col_w + 1)
        first_row = max(0, voff // self.cell_h)
        last_row = min(nrows, (voff + (rows_bottom - oy)) // self.cell_h + 1)
        draw_letters = col_w >= 12 and self.cell_h >= 9

        # ---- cells (the residue-probe tint is drawn *under* the letters) ----
        p.save()
        p.setClipRect(QRect(ox, oy, vw - ox, rows_bottom - oy))
        p.setFont(self._letter)
        probe_on = self._probe is not None
        for c in range(first_col, last_col):
            x = ox + c * col_w - hoff
            kept = self._col_kept(c)
            pv = self._agg_probe(c) if probe_on else float("nan")
            tint = pv == pv and pv > 0.03
            for r in range(first_row, last_row):
                y = oy + r * self.cell_h - voff
                cell = self._dm.cells[r][c]
                is_gap = cell.strip("-") == ""
                if is_gap:
                    col = colors.GAP
                elif self.color_mode in ("confidence", "truth") and self._cell_scores is not None:
                    col = colors.score_color(self._agg_cell(r, c))
                else:
                    col = colors.residue_color(cell, self.level, self.aln.alphabet, table, scheme)
                masked = not is_gap and self._cell_masked(r, c)
                if masked:
                    col = colors.MASKED
                if not kept:
                    col = col.darker(220)
                p.fillRect(x, y, col_w - 1, self.cell_h - 1, col)
                if masked:
                    p.fillRect(x, y, col_w - 1, self.cell_h - 1,
                               QBrush(colors.MUTED, Qt.BDiagPattern))
                if tint and not is_gap:
                    # Gamma-boost so faint alternative columns are visible, not
                    # just the residue's own (always-strongest) column.
                    a = int(55 + 195 * (min(1.0, pv) ** 0.5))
                    p.fillRect(x, y, col_w - 1, self.cell_h - 1, QColor(216, 94, 240, a))
                if draw_letters and not is_gap:
                    p.setPen(colors.MUTED if masked else QColor(20, 24, 28))
                    p.drawText(QRect(x, y, col_w - 1, self.cell_h - 1), Qt.AlignCenter, cell)
                if (r, c) in self._find:
                    p.setPen(QColor(255, 214, 0))
                    p.drawRect(x, y, col_w - 2, self.cell_h - 2)
            if not kept:
                # Mask: clearly dim excluded columns, but only over the sequences.
                p.fillRect(x, oy, col_w - 1, grid_bottom - oy, QColor(8, 10, 14, 150))
        # outline the probed residue (its candidate columns are the tinted ones)
        if probe_on and self._probe_cell and self._probe_cell[0] >= 0:
            pseq, pntcol = self._probe_cell
            pc = (pntcol - self._dm.nt_start[0]) // max(1, span)
            px = ox + pc * col_w - hoff
            py = oy + pseq * self.cell_h - voff
            p.setPen(QColor(255, 255, 255))
            p.drawRect(px, py, col_w - 1, self.cell_h - 1)
        # edit cursor (keyboard editing): a bold amber box on the focused cell
        if self._cursor is not None:
            cr, cc = self._cursor
            if first_row <= cr < last_row and first_col <= cc < last_col:
                cx = ox + cc * col_w - hoff
                cy = oy + cr * self.cell_h - voff
                pen = p.pen()
                p.setPen(QColor(255, 170, 0))
                p.drawRect(cx, cy, col_w - 1, self.cell_h - 1)
                p.drawRect(cx + 1, cy + 1, col_w - 3, self.cell_h - 3)
                p.setPen(pen)
        # selected sequences (edit target): a faint amber band across each row
        for r in self._sel_rows:
            if first_row <= r < last_row:
                ry = oy + r * self.cell_h - voff
                p.fillRect(ox, ry, vw - ox, self.cell_h - 1, QColor(255, 170, 0, 28))
        # pinned columns: a solid band at the top of the grid, always visible
        if self._pins:
            for c in range(first_col, last_col):
                if self._col_pinned(c):
                    px = ox + c * col_w - hoff
                    p.fillRect(px, oy, col_w, grid_bottom - oy, QColor(255, 196, 60, 22))
                    p.fillRect(px, oy, col_w, 3, QColor(255, 196, 60, 230))
        # selection overlay
        if self._sel_cols:
            c0, c1 = self._sel_cols
            sx = ox + c0 * col_w - hoff
            sw = (c1 - c0) * col_w
            p.fillRect(sx, oy, sw, grid_bottom - oy, QColor(120, 170, 255, 45))
            p.setPen(QColor(120, 170, 255, 200))
            p.drawRect(sx, oy, sw - 1, grid_bottom - oy - 1)
        p.restore()

        # ---- consensus band ----
        if cons_h:
            cy = rows_bottom
            p.save()
            p.setClipRect(QRect(ox, cy, vw - ox, cons_h))
            p.setFont(self._letter)
            p.fillRect(ox, cy, vw - ox, cons_h, colors.BG)
            cons = self._consensus or []
            for c in range(first_col, last_col):
                if c >= len(cons):
                    break
                x = ox + c * col_w - hoff
                cell = cons[c]
                if cell.strip("-") == "":
                    col = colors.GAP
                else:
                    col = colors.residue_color(cell, self.level, self.aln.alphabet, table, scheme)
                p.fillRect(x, cy + 1, col_w - 1, cons_h - 2, col)
                if draw_letters and cell.strip("-"):
                    p.setPen(QColor(20, 24, 28))
                    p.drawText(QRect(x, cy, col_w - 1, cons_h), Qt.AlignCenter, cell)
            p.restore()

        # ---- annotation band (column features) ----
        if annot_h:
            ay = rows_bottom + cons_h
            p.save()
            p.setClipRect(QRect(ox, ay, vw - ox, annot_h))
            p.setFont(self.font())
            p.fillRect(ox, ay, vw - ox, annot_h, colors.PANEL)
            for ann in self._annotations:
                cols = self._nt_to_cols(int(ann.get("start", 0)), int(ann.get("stop", 0)))
                if not cols:
                    continue
                c0, c1 = cols
                x0 = ox + c0 * col_w - hoff
                x1 = ox + (c1 + 1) * col_w - hoff
                if x1 < ox or x0 > vw:
                    continue
                col = QColor(ann.get("color") or "#5b8def")
                p.fillRect(int(x0), ay + 2, int(x1 - x0) - 1, annot_h - 4, col)
                name = str(ann.get("name", ""))
                if name:
                    p.setPen(QColor(20, 24, 28))
                    label = self._fm.elidedText(name, Qt.ElideRight, int(x1 - x0) - 6)
                    p.drawText(QRect(int(x0) + 3, ay, int(x1 - x0) - 6, annot_h),
                               Qt.AlignVCenter, label)
            p.restore()

        # ---- ruler ----
        p.save()
        p.setClipRect(QRect(ox, self.track_h, vw - ox, self.ruler_h))
        p.fillRect(ox, self.track_h, vw - ox, self.ruler_h, colors.PANEL)
        p.setPen(colors.MUTED)
        step = max(1, round(55 / col_w)) * (5 if span > 1 else 10)
        for c in range(first_col, last_col):
            if c % step == 0:
                x = ox + c * col_w - hoff
                p.drawLine(x, self.track_h + self.ruler_h - 4, x, self.track_h + self.ruler_h)
                p.drawText(x + 2, self.track_h + self.ruler_h - 4, str(c + 1))
        p.restore()

        # ---- score track ----
        p.save()
        p.setClipRect(QRect(ox, 0, vw - ox, self.track_h))
        p.fillRect(ox, 0, vw - ox, self.track_h, colors.PANEL)
        if self._nt_scores is not None:
            for c in range(first_col, last_col):
                x = ox + c * col_w - hoff
                p.fillRect(x, 2, col_w - 1, self.track_h - 4, colors.score_color(self._agg_score(c)))
        p.restore()

        # ---- name gutter ----
        p.save()
        p.setClipRect(QRect(0, oy, self.name_w, rows_bottom - oy))
        p.fillRect(0, oy, self.name_w, rows_bottom - oy, colors.PANEL)
        for r in range(first_row, last_row):
            y = oy + r * self.cell_h - voff
            if r in self._sel_rows:
                p.fillRect(0, y, self.name_w, self.cell_h, QColor(255, 170, 0, 80))
                p.fillRect(0, y, 3, self.cell_h, QColor(255, 170, 0))
            if r == self._find_row:
                p.fillRect(0, y, self.name_w, self.cell_h, QColor(255, 214, 0, 70))
            if r in self._outlier_rows:
                p.fillRect(0, y, 5, self.cell_h, QColor(220, 40, 40))
            elided = self._fm.elidedText(self.aln.ids[r], Qt.ElideRight, self.name_w - 10)
            p.setPen(colors.TEXT)
            p.drawText(QRect(6, y, self.name_w - 10, self.cell_h), Qt.AlignVCenter, elided)
        p.restore()

        # ---- consensus gutter label ----
        if cons_h:
            p.fillRect(0, rows_bottom, self.name_w, cons_h, colors.GRID)
            p.setPen(colors.TEXT)
            p.drawText(QRect(6, rows_bottom, self.name_w - 8, cons_h), Qt.AlignVCenter, "consensus")

        # ---- annotation gutter label ----
        if annot_h:
            ay = rows_bottom + cons_h
            p.fillRect(0, ay, self.name_w, annot_h, colors.GRID)
            p.setPen(colors.TEXT)
            p.drawText(QRect(6, ay, self.name_w - 8, annot_h), Qt.AlignVCenter, "annotations")

        # ---- corner box ----
        p.fillRect(0, 0, self.name_w, oy, colors.GRID)
        p.setPen(colors.TEXT)
        p.drawText(QRect(6, 0, self.name_w - 8, self.track_h), Qt.AlignVCenter,
                   self._score_label or "no track")
        p.setPen(colors.MUTED)
        p.drawText(QRect(6, self.track_h, self.name_w - 8, self.ruler_h),
                   Qt.AlignVCenter, f"position ({self.level.value})")

    # -- interaction ------------------------------------------------------- #
    def _col_at(self, x: int) -> Optional[int]:
        ox, _ = self._grid_origin()
        if x < ox or not self._dm:
            return None
        c = (x - ox + self.horizontalScrollBar().value()) // self._col_w()
        return int(min(max(0, c), self._dm.n_cols - 1))

    def _row_at(self, y: int) -> Optional[int]:
        _, oy = self._grid_origin()
        if y < oy or not self.aln:
            return None
        r = (y - oy + self.verticalScrollBar().value()) // self.cell_h
        if r < 0 or r >= self.aln.n_seqs:
            return None
        return int(r)

    def select_all(self) -> None:
        """Select every column of the alignment."""
        if not self._dm:
            return
        self._sel_cols = (0, self._dm.n_cols)
        self.viewport().update()
        nt = self.selection_nt()
        if nt:
            self.selectionChanged.emit(nt[0], nt[1])

    def select_nt_range(self, nt0: int, nt1: int) -> None:
        """Select (and scroll to) the display columns covering nucleotide [nt0, nt1)."""
        if not self._dm:
            return
        cols = [c for c in range(self._dm.n_cols) if nt0 <= self._dm.nt_start[c] < nt1]
        if not cols:
            return
        self._sel_cols = (cols[0], cols[-1] + 1)
        self._ensure_col_visible(cols[0])
        self.viewport().update()
        nt = self.selection_nt()
        if nt:
            self.selectionChanged.emit(nt[0], nt[1])

    # -- edit cursor & row selection --------------------------------------- #
    def cursor(self) -> Optional[Tuple[int, int]]:
        """The edit cursor as (row, display column), or None."""
        return self._cursor

    def selected_rows(self) -> Set[int]:
        """The sequences chosen for editing (falls back to the cursor's row)."""
        if self._sel_rows:
            return set(self._sel_rows)
        return {self._cursor[0]} if self._cursor else set()

    def set_selected_rows(self, rows: Set[int]) -> None:
        if self.aln:
            rows = {r for r in rows if 0 <= r < self.aln.n_seqs}
        self._sel_rows = set(rows)
        self.viewport().update()
        self.rowsSelected.emit(len(self._sel_rows))

    def clear_selection(self) -> None:
        """Drop both the column range and the sequence selection (the Esc key)."""
        had = bool(self._sel_cols) or bool(self._sel_rows)
        self._sel_cols = None
        self._sel_rows = set()
        self._row_anchor = None
        self.viewport().update()
        if had:
            self.rowsSelected.emit(0)
            self.selectionCleared.emit()

    def _select_row(self, row: int, toggle: bool) -> None:
        """Plain click selects just this row; Cmd/Ctrl-click toggles it in/out."""
        if toggle:
            sel = set(self._sel_rows)
            sel.discard(row) if row in sel else sel.add(row)
        else:
            sel = {row}
        self.set_selected_rows(sel)

    def _click_select_row(self, row: int, toggle: bool, shift: bool) -> None:
        """Plain = just this row; Cmd/Ctrl = toggle; Shift = range from the anchor.

        A plain click *inside* an existing multi-selection keeps it (so you can
        position the cursor within a group without losing the group); a plain click
        on a sequence outside the selection starts fresh with just that one.
        """
        if shift and self._row_anchor is not None:
            lo, hi = sorted((self._row_anchor, row))
            self.set_selected_rows(set(range(lo, hi + 1)))   # anchor kept, so you can extend
        elif toggle:
            self._select_row(row, True)
            self._row_anchor = row
        elif row in self._sel_rows and len(self._sel_rows) > 1:
            self._row_anchor = row                           # keep the group, just move on
        else:
            self._select_row(row, False)
            self._row_anchor = row

    def set_cursor(self, row: int, col: int) -> None:
        if not self._dm:
            return
        row = int(min(max(0, row), self.aln.n_seqs - 1))
        col = int(min(max(0, col), self._dm.n_cols - 1))
        self._cursor = (row, col)
        self._ensure_row_visible(row)
        self._ensure_col_visible(col)
        self.viewport().update()
        self.cursorMoved.emit(row, col)

    def keyPressEvent(self, e):
        if not self.aln or not self._dm:
            super().keyPressEvent(e)
            return
        key, mod = e.key(), e.modifiers()
        zoom_mod = bool(mod & (Qt.ControlModifier | Qt.MetaModifier))
        # Ctrl/Cmd + '+' / '-' / '0' change the residue font size (zoom).
        if zoom_mod and key in (Qt.Key_Plus, Qt.Key_Equal):
            self.zoom(1); return
        if zoom_mod and key in (Qt.Key_Minus, Qt.Key_Underscore):
            self.zoom(-1); return
        if zoom_mod and key == Qt.Key_0:
            self.reset_zoom(); return
        if key == Qt.Key_Escape:
            self.clear_selection(); return
        if self._cursor is None:
            self.set_cursor(0, 0)
            return
        r, c = self._cursor
        if key == Qt.Key_Left:
            self.set_cursor(r, c - 1); return
        if key == Qt.Key_Right:
            self.set_cursor(r, c + 1); return
        if key == Qt.Key_Up:
            self.set_cursor(r - 1, c); return
        if key == Qt.Key_Down:
            self.set_cursor(r + 1, c); return
        # Editing keys ignore Ctrl/Cmd (those are zoom, above). '-' / Space slide
        # the selected sequences right (open a gap); Backspace / Delete / Shift+Space
        # slide them left (close a gap to the left).
        if not zoom_mod and not self.read_only:
            if key in (Qt.Key_Minus, Qt.Key_Space) and not (mod & Qt.ShiftModifier):
                self.editRequested.emit("slide_right"); return
            if key in (Qt.Key_Backspace, Qt.Key_Delete) or (
                    key == Qt.Key_Space and mod & Qt.ShiftModifier):
                self.editRequested.emit("slide_left"); return
        super().keyPressEvent(e)

    def _extend_to_x(self, x: int) -> None:
        """Grow the column selection to the cursor, clamped to the visible grid."""
        if self._anchor is None or not self._dm:
            return
        ox, _ = self._grid_origin()
        c = self._col_at(min(max(x, ox), self.viewport().width() - 1))
        if c is None:
            return
        if c != self._anchor:
            self._moved = True
        lo, hi = min(self._anchor, c), max(self._anchor, c)
        self._sel_cols = (lo, hi + 1)
        self.viewport().update()

    def _update_autoscroll(self, x: int) -> None:
        ox, _ = self._grid_origin()
        edge = 24
        if x > self.viewport().width() - edge:
            self._autoscroll_dir = 1
        elif x < ox + edge:
            self._autoscroll_dir = -1
        else:
            self._autoscroll_dir = 0
        if self._autoscroll_dir and not self._autoscroll.isActive():
            self._autoscroll.start()
        elif not self._autoscroll_dir and self._autoscroll.isActive():
            self._autoscroll.stop()

    def _autoscroll_tick(self) -> None:
        if self._anchor is None or not self._autoscroll_dir:
            self._autoscroll.stop()
            return
        bar = self.horizontalScrollBar()
        bar.setValue(bar.value() + self._autoscroll_dir * self._col_w())
        self._extend_to_x(self._drag_x)

    def mousePressEvent(self, e):
        x = int(e.position().x())
        c = self._col_at(x)
        r = self._row_at(int(e.position().y()))
        self._press = (c, r)
        self._press_x = x
        self._press_mod = e.modifiers()
        self._moved = False
        # Drag to slide: grab a selected (amber) sequence, or Alt-drag any sequence.
        alt = bool(e.modifiers() & Qt.AltModifier)
        grab = (not self.read_only and r is not None and c is not None
                and (r in self._sel_rows or alt))
        if grab:
            if r not in self._sel_rows:
                self._select_row(r, toggle=False)
            self.set_cursor(r, c)
            self._slide = True
            self._slide_col = c
            self._anchor = None
        else:
            self._slide = False
            self._anchor = c

    def mouseMoveEvent(self, e):
        x = int(e.position().x())
        if self._slide:
            self._moved = True
            ox, _ = self._grid_origin()
            c = self._col_at(min(max(x, ox), self.viewport().width() - 1))
            if c is None or self._cursor is None:
                return
            # Emit one slide per column boundary crossed; the app applies it and
            # moves the cursor, so we read the cursor back to detect a refused slide.
            while c > self._slide_col:
                before = self._cursor[1]
                self.editRequested.emit("slide_right")
                if self._cursor[1] == before:
                    break
                self._slide_col = self._cursor[1]
            while c < self._slide_col:
                before = self._cursor[1]
                self.editRequested.emit("slide_left")
                if self._cursor[1] == before:
                    break
                self._slide_col = self._cursor[1]
            return
        if self._anchor is None:
            return
        self._drag_x = x
        self._extend_to_x(x)
        self._update_autoscroll(x)

    def mouseReleaseEvent(self, e):
        self._autoscroll.stop()
        self._autoscroll_dir = 0
        sliding = self._slide
        self._slide = False
        self._anchor = None
        if self._moved:
            if not sliding:
                nt = self.selection_nt()
                if nt:
                    self.selectionChanged.emit(nt[0], nt[1])
            return
        if self._press:
            c, r = self._press
            toggle = bool(self._press_mod & (Qt.MetaModifier | Qt.ControlModifier))
            shift = bool(self._press_mod & Qt.ShiftModifier)
            if r is not None and c is None and self._press_x < self.name_w:
                # Clicking a sequence name selects it (Cmd toggles, Shift = range).
                self._click_select_row(r, toggle, shift)
                col = self._cursor[1] if self._cursor else 0
                self.set_cursor(r, col)
            elif c is not None and r is not None and self._dm is not None:
                # A plain click places the edit cursor, selects that one sequence,
                # and probes the residue (its homology cloud). Cmd/Ctrl-click adds a
                # sequence; Shift-click selects a contiguous range of sequences.
                self.set_cursor(r, c)
                self._click_select_row(r, toggle, shift)
                if not toggle and not shift:
                    self.residueClicked.emit(r, self._dm.nt_start[c])

    def wheelEvent(self, e):
        if e.modifiers() & (Qt.ControlModifier | Qt.MetaModifier):
            self.zoom(1 if e.angleDelta().y() > 0 else -1)
        else:
            super().wheelEvent(e)

    def zoom(self, delta: int) -> None:
        """Grow/shrink the cells (and the residue font) by ``delta`` steps."""
        self.set_cell_size(self.cell_w + delta, self.cell_h + delta)

    def reset_zoom(self) -> None:
        self.set_cell_size(15, 17)
