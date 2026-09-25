"""Dockable panels: the posterior explorer, the realignment sandbox, and the
column inspector."""

from __future__ import annotations

import html
from typing import List, Optional, Sequence, Tuple

from PySide6.QtCore import QRect, Qt, Signal
from PySide6.QtGui import QPainter
from PySide6.QtWidgets import (
    QComboBox, QHBoxLayout, QLabel, QListWidget, QListWidgetItem,
    QPlainTextEdit, QProgressDialog, QPushButton, QTextBrowser, QVBoxLayout,
    QWidget,
)

from ..ambiguity import posterior as post_mod
from ..ambiguity import sandbox as sandbox_mod
from ..domain import Alignment, Level
from ..engines import AlignerEngine
from ..evaluate import ColumnTruth, Placement
from . import colors
from .workers import run_cancellable

# The posterior explorer inspects *local* ambiguous regions. Building and painting
# the pair-HMM matrix for a huge selection (e.g. Select all on a long alignment)
# would allocate millions of cells on the GUI thread and freeze the app, so we cap
# the region and ask the user to narrow it instead.
MAX_POSTERIOR_COLS = 200


def exclusive_pair(combo_a: QComboBox, combo_b: QComboBox) -> None:
    """Stop two sequence pickers from selecting the same sequence.

    A sequence does have a posterior against itself, but it is the trivial bright
    diagonal: it answers no question, and every consumer of the pair has to
    special-case it. Rather than let the user choose it and then refuse the
    result, the choice is withdrawn — each combo's entry for the other's current
    selection is greyed out, and a collision is resolved by moving the second
    picker on by one.

    Call it whenever either selection changes, and once after populating.
    """
    n = combo_a.count()
    if n < 2 or combo_b.count() != n:
        return                       # one sequence: there is no valid pair to offer
    if combo_a.currentIndex() == combo_b.currentIndex():
        combo_b.blockSignals(True)
        combo_b.setCurrentIndex((combo_b.currentIndex() + 1) % n)
        combo_b.blockSignals(False)
    taken_a, taken_b = combo_a.currentIndex(), combo_b.currentIndex()
    for combo, taken in ((combo_a, taken_b), (combo_b, taken_a)):
        model = combo.model()
        for k in range(combo.count()):
            item = model.item(k) if hasattr(model, "item") else None
            if item is not None:
                item.setEnabled(k != taken)


# --------------------------------------------------------------------------- #
# Posterior explorer
# --------------------------------------------------------------------------- #

class ScoreLegend(QWidget):
    """Colour key for the reliability ramp.

    The confidence colouring and the column score track both use
    ``colors.score_color``, a continuous red-amber-green ramp. Without a key,
    a residue painted some shade of orange tells the user only that it is not at
    an extreme — not whether it is better or worse than the orange beside it. So
    the ramp is drawn here with its end points labelled, plus the separate grey
    used for columns that could not be scored at all, which is not a low score
    and must not be read as one.
    """

    _PAD = 12          # room for the "0" and "1" end labels, outside the ramp

    #: What the ramp means, per colouring mode. The bar is identical in both —
    #: what changes is the claim it makes, and reading a residue as "unreliable"
    #: when it is in fact "wrongly placed" is exactly the confusion to avoid.
    _MEANING = {
        "confidence":
            "Confidence: red = the alignment is guessing here, green = well "
            "supported. Grey means the column could not be scored (no "
            "residue-pair evidence), which is not the same as a low score.",
        "truth":
            "Correct placement: green = this residue sits with the partners the "
            "known-true alignment gives it, red = it does not. Grey means the "
            "residue has no partners in its column, so it claims nothing.",
    }

    def __init__(self, parent=None, width: int = 120, height: int = 20):
        super().__init__(parent)
        self._bar_w = width
        self.setFixedSize(width + 2 * self._PAD + 62, height)
        self.set_meaning("confidence")

    def set_meaning(self, mode: str) -> None:
        """Say what the ramp is measuring right now."""
        self.setToolTip(self._MEANING.get(mode, self._MEANING["confidence"]))

    def paintEvent(self, _ev):
        from PySide6.QtGui import QFont

        p = QPainter(self)
        p.setRenderHint(QPainter.TextAntialiasing, True)
        h = self.height()
        bar_h = max(8, h - 8)
        y = (h - bar_h) // 2
        # (y >= 1 so the outline drawn one pixel above the ramp stays in view)

        f = QFont()
        f.setPixelSize(9)
        p.setFont(f)

        # End labels sit *outside* the ramp: printed on top of it they land on
        # saturated red and green and are barely readable.
        x0 = self._PAD
        p.setPen(colors.MUTED)
        p.drawText(QRect(0, y, x0 - 2, bar_h), Qt.AlignVCenter | Qt.AlignRight, "0")
        p.drawText(QRect(x0 + self._bar_w + 2, y, x0 - 2, bar_h),
                   Qt.AlignVCenter | Qt.AlignLeft, "1")

        # The ramp, sampled a pixel at a time so it matches score_color exactly
        # rather than approximating it with a QLinearGradient's own interpolation.
        for x in range(self._bar_w):
            p.fillRect(QRect(x0 + x, y, 1, bar_h),
                       colors.score_color(x / max(1, self._bar_w - 1)))
        # Outline the ramp without painting over its end colours, which are the
        # two the user most needs to recognise.
        p.setPen(colors.GRID)
        p.drawRect(QRect(x0 - 1, y - 1, self._bar_w + 1, bar_h + 1))

        # the unscored swatch, set apart from the ramp
        sw_x = x0 + self._bar_w + 2 * self._PAD
        p.fillRect(QRect(sw_x, y, bar_h, bar_h), colors.score_color(float("nan")))
        p.setPen(colors.GRID)
        p.drawRect(QRect(sw_x, y, bar_h - 1, bar_h - 1))
        p.setPen(colors.MUTED)
        p.drawText(QRect(sw_x + bar_h + 4, y, 44, bar_h),
                   Qt.AlignVCenter | Qt.AlignLeft, "n/a")
        p.end()


class PosteriorHeatmap(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.view: Optional[post_mod.PosteriorView] = None
        self.cmap = "Viridis"
        self.setMinimumHeight(220)

    def set_view(self, view):
        self.view = view
        self.update()

    def set_cmap(self, name: str):
        self.cmap = name
        self.update()

    def paintEvent(self, e):
        p = QPainter(self)
        p.fillRect(self.rect(), colors.BG)
        v = self.view
        if v is None or v.matrix.size == 0:
            p.setPen(colors.MUTED)
            p.drawText(self.rect(), Qt.AlignCenter,
                       "Select a column range to see\nalternative homology hypotheses")
            return
        M = v.matrix
        ri, rj = M.shape
        margin = 16
        legend_h = 18
        gw = self.width() - 2 * margin
        gh = self.height() - 2 * margin - legend_h
        cw = max(2, gw // max(1, rj))
        ch = max(2, gh // max(1, ri))
        cw = ch = min(cw, ch, 22)
        x0 = margin
        y0 = margin
        for i in range(ri):
            for j in range(rj):
                p.fillRect(x0 + j * cw, y0 + i * ch, cw - 1, ch - 1,
                           colors.colormap_color(float(M[i, j]), self.cmap))
        if cw >= 11:
            p.setPen(colors.MUTED)
            for j, ch_j in enumerate(v.chars_j[:rj]):
                p.drawText(QRect(x0 + j * cw, 0, cw, margin), Qt.AlignCenter, ch_j)
            for i, ch_i in enumerate(v.chars_i[:ri]):
                p.drawText(QRect(0, y0 + i * ch, margin, ch), Qt.AlignCenter, ch_i)
        self._draw_legend(p, margin, legend_h)

    def _draw_legend(self, p, margin: int, legend_h: int):
        """A horizontal colorbar (0 → 1) so the colours are readable as values."""
        w = self.width() - 2 * margin
        if w <= 20:
            return
        y = self.height() - legend_h
        bar_h = 8
        for k in range(w):
            p.fillRect(margin + k, y, 1, bar_h, colors.colormap_color(k / (w - 1), self.cmap))
        p.setPen(colors.MUTED)
        p.drawText(QRect(margin, y + bar_h, 24, legend_h - bar_h),
                   Qt.AlignLeft | Qt.AlignVCenter, "0")
        p.drawText(QRect(margin + w - 24, y + bar_h, 24, legend_h - bar_h),
                   Qt.AlignRight | Qt.AlignVCenter, "1")


class PosteriorPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.aln: Optional[Alignment] = None
        self.region: Optional[Tuple[int, int]] = None
        self.aa_mode = False
        self._aa_aln: Optional[Alignment] = None

        self.seq_i = QComboBox()
        self.seq_j = QComboBox()
        self.seq_i.currentIndexChanged.connect(self._on_pair_changed)
        self.seq_j.currentIndexChanged.connect(self._on_pair_changed)
        self.mode = QComboBox()
        self.mode.addItems(["Posterior", "Residual (minus current)", "Profile (i vs rest)"])
        self.mode.currentIndexChanged.connect(self._refresh)
        self.info = QLabel("pair-HMM posterior P(res i ~ res j)")
        # Theme default colour (readable on light *or* dark); italic keeps it a caption.
        self.info.setStyleSheet("font-style:italic;")
        self.info.setWordWrap(True)
        self.heatmap = PosteriorHeatmap()
        self.cmap = QComboBox()
        self.cmap.addItems(colors.COLORMAP_NAMES)
        self.cmap.setToolTip("Heatmap colour scheme")
        self.cmap.currentTextChanged.connect(self.heatmap.set_cmap)

        top = QHBoxLayout()
        top.addWidget(QLabel("seq"))
        top.addWidget(self.seq_i, 1)
        top.addWidget(QLabel("vs"))
        top.addWidget(self.seq_j, 1)
        moderow = QHBoxLayout()
        moderow.addWidget(QLabel("mode"))
        moderow.addWidget(self.mode, 1)
        moderow.addWidget(QLabel("colour"))
        moderow.addWidget(self.cmap, 1)
        lay = QVBoxLayout(self)
        lay.addLayout(top)
        lay.addLayout(moderow)
        lay.addWidget(self.info)
        lay.addWidget(self.heatmap, 1)

    def set_alignment(self, aln: Alignment):
        self.aln = aln
        self._aa_aln = None
        self.region = None                # forget the previous file's selection
        for combo in (self.seq_i, self.seq_j):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(aln.ids)
            combo.blockSignals(False)
        if aln.n_seqs > 1:
            self.seq_j.setCurrentIndex(1)
        exclusive_pair(self.seq_i, self.seq_j)
        self.heatmap.set_view(None)
        self.info.setText("select a column range")

    def _on_pair_changed(self, *args):
        exclusive_pair(self.seq_i, self.seq_j)
        self._refresh()

    def set_region(self, nt0: int, nt1: int):
        self.region = (nt0, nt1)
        self._refresh()

    def clear_region(self):
        self.region = None
        self.heatmap.set_view(None)
        self.info.setText("select a column range")

    def set_view(self, aa_mode: bool):
        """When True (viewing codons/aa of a coding alignment), compute the
        posterior in amino-acid space so it matches what's on screen."""
        self.aa_mode = aa_mode
        self._refresh()

    def _refresh(self, *args):
        if not self.aln or self.region is None:
            return
        c0, c1 = self.region
        if self.aa_mode and self.aln.coding is not None:
            if self._aa_aln is None:
                self._aa_aln = self.aln.amino_acid_alignment()
            work = self._aa_aln
            a0, a1 = self.aln.codon_col(c0), self.aln.codon_col(c1)
            unit = "amino acids"
        else:
            work = self.aln
            a0, a1 = c0, c1
            unit = "residues"
        if a1 - a0 > MAX_POSTERIOR_COLS:
            self.heatmap.set_view(None)
            self.info.setText(
                f"Region too large for the posterior view ({a1 - a0} columns). "
                f"Select up to {MAX_POSTERIOR_COLS} columns to explore alternatives.")
            return
        i = self.seq_i.currentIndex()
        j = self.seq_j.currentIndex()
        mode_idx = self.mode.currentIndex()
        self.seq_j.setEnabled(mode_idx != 2)
        try:
            if mode_idx == 2:  # profile: i vs the rest
                view = post_mod.region_profile(work, i, a0, a1)
            else:
                if i < 0 or j < 0 or i == j:
                    self.heatmap.set_view(None)
                    self.info.setText("pick two different sequences")
                    return
                m = "residual" if mode_idx == 1 else "posterior"
                view = post_mod.region_posterior(work, i, j, a0, a1, mode=m)
        except Exception:
            self.heatmap.set_view(None)
            return
        self.heatmap.set_view(view)
        if not view.matrix.size:
            self.info.setText("no residues in this region")
            return
        if mode_idx == 2:
            self.info.setText(
                f"profile · {work.ids[i]} vs the rest ({unit})  |  "
                f"{view.matrix.shape[0]} × {view.matrix.shape[1]}  |  "
                "each row shows where that residue could sit"
            )
        elif mode_idx == 1:
            self.info.setText(
                f"residual ({unit})  |  committed alignment removed  |  "
                f"strongest remaining alternative = {float(view.matrix.max()):.2f}"
            )
        else:
            off = view.matrix.copy()
            for d in range(min(off.shape)):
                off[d, d] = 0
            self.info.setText(
                f"posterior ({unit})  |  {view.matrix.shape[0]}x{view.matrix.shape[1]}  |  "
                f"max off-diagonal homology prob = {off.max():.2f}"
            )


# --------------------------------------------------------------------------- #
# Realignment sandbox
# --------------------------------------------------------------------------- #

class SandboxPanel(QWidget):
    applied = Signal(object)  # new full Alignment after splice

    def __init__(self, parent=None):
        super().__init__(parent)
        self.aln: Optional[Alignment] = None
        self.engines: Sequence[AlignerEngine] = []
        #: engine -> its parameter values. The window sets this to its own
        #: settings lookup, so alternatives use what the ⚙ dialog holds.
        self.params_for = lambda engine: {}
        self.region: Optional[Tuple[int, int]] = None
        self._alts: List[sandbox_mod.Alternative] = []
        self.codon_aware = False
        self.table = 1
        self.view_level = Level.NT
        self._progress: Optional[QProgressDialog] = None
        self._cancel_token = None
        self._cancelled = False

        self.header = QLabel("Select a column range, then generate alternatives")
        # Use the theme's default text colour (readable on light *or* dark panels);
        # bold keeps it prominent as a status line without hard-coding a colour.
        self.header.setStyleSheet("font-weight:600;")
        self.gen_btn = QPushButton("Generate alternatives")
        self.gen_btn.clicked.connect(self._generate)
        self.gen_btn.setEnabled(False)
        self.alt_list = QListWidget()
        self.alt_list.currentRowChanged.connect(self._preview)
        self.preview = QPlainTextEdit()
        self.preview.setReadOnly(True)
        self.preview.setLineWrapMode(QPlainTextEdit.NoWrap)   # scroll long rows, don't wrap
        self.preview.setStyleSheet("font-family:Menlo,monospace; background:#14181d; color:#dfe6ee;")
        self.apply_btn = QPushButton("Apply selected → splice into alignment")
        self.apply_btn.clicked.connect(self._apply)
        self.apply_btn.setEnabled(False)

        lay = QVBoxLayout(self)
        lay.addWidget(self.header)
        lay.addWidget(self.gen_btn)
        lay.addWidget(self.alt_list, 1)
        lay.addWidget(self.preview, 1)
        lay.addWidget(self.apply_btn)

    def set_context(self, aln: Alignment, engines: Sequence[AlignerEngine]):
        # a new alignment fully resets the sandbox — no memory of the old file
        self.aln = aln
        self.engines = engines
        self.region = None
        self._alts = []
        self.alt_list.clear()
        self.preview.clear()
        self.header.setText("Select a column range, then generate alternatives")
        self.gen_btn.setEnabled(False)
        self.apply_btn.setEnabled(False)

    def set_region(self, nt0: int, nt1: int):
        self.region = (nt0, nt1)
        self.header.setText(f"Region: columns {nt0+1}–{nt1}  ({nt1-nt0} nt)")
        self.gen_btn.setEnabled(self.aln is not None and nt1 - nt0 >= 2)

    def clear_region(self):
        self.region = None
        self.header.setText("Select a column range, then generate alternatives")
        self.gen_btn.setEnabled(False)

    def set_view(self, codon_aware: bool, table: int, level: Level):
        """Tell the sandbox what level the user is viewing, so alternatives are
        aligned codon-aware and previewed at that level."""
        self.codon_aware = codon_aware
        self.table = table
        self.view_level = level

    def _generate(self):
        if not self.aln or self.region is None:
            return
        c0, c1 = self.region
        self.gen_btn.setEnabled(False)
        self._cancelled = False
        engines = list(self.engines)
        params_by_key = {e.key: self.params_for(e) for e in engines}
        aln = self.aln
        codon_aware, table, level = self.codon_aware, self.table, self.view_level

        dlg = QProgressDialog("Generating alternatives…", "Cancel", 0, 0, self)
        dlg.setWindowTitle("Realignment sandbox")
        dlg.setWindowModality(Qt.WindowModal)
        dlg.setMinimumDuration(0)        # show at once — these runs can be slow
        dlg.setAutoClose(False)
        dlg.setAutoReset(False)
        dlg.canceled.connect(self._cancel_generate)
        self._progress = dlg

        def work(report, cancelled):
            return sandbox_mod.alternatives(aln, c0, c1, engines, codon_aware=codon_aware,
                                            table=table, level=level,
                                            params_by_key=params_by_key,
                                            progress=report, cancelled=cancelled)

        self._cancel_token = run_cancellable(work, self._on_alts, self._on_err,
                                             self._on_progress)
        dlg.show()

    def _on_progress(self, done: int, total: int, label: str):
        dlg = self._progress
        if dlg is None:
            return
        if total > 0 and dlg.maximum() != total:
            dlg.setRange(0, total)
        if not self._cancelled:
            dlg.setLabelText(label)
        # A modal QProgressDialog re-enters the event loop inside setValue(), which
        # can deliver the 'done' signal and close the dialog — so do it last and
        # touch nothing on the dialog afterwards.
        dlg.setValue(done)

    def _cancel_generate(self):
        self._cancelled = True
        if self._cancel_token is not None:
            self._cancel_token.cancel()
        if self._progress is not None:
            self._progress.setLabelText("Cancelling — finishing the current engine…")

    def _close_progress(self):
        if self._progress is not None:
            # closing emits canceled() (Qt treats it as a cancel) — drop the
            # connection first so a normal finish isn't mistaken for a cancel.
            try:
                self._progress.canceled.disconnect(self._cancel_generate)
            except (RuntimeError, TypeError):
                pass
            self._progress.close()
            self._progress = None
        self._cancel_token = None

    def _on_alts(self, alts):
        self._close_progress()
        self.gen_btn.setEnabled(True)
        if self._cancelled:
            self._cancelled = False
            return
        self._alts = sorted(alts, key=lambda a: -a.score)
        self.alt_list.clear()
        for a in self._alts:
            QListWidgetItem(f"{a.label}   —   width {a.width}, SP-id {a.score:.2f}", self.alt_list)
        if self._alts:
            self.alt_list.setCurrentRow(0)

    def _on_err(self, msg):
        self._close_progress()
        self.gen_btn.setEnabled(True)
        if self._cancelled:
            self._cancelled = False
            return
        self.preview.setPlainText("Error:\n" + msg)

    def _preview(self, row: int):
        self.apply_btn.setEnabled(0 <= row < len(self._alts))
        if not (0 <= row < len(self._alts)):
            return
        alt = self._alts[row]
        rows = alt.display_rows()
        w = max((len(i) for i in alt.block.ids), default=0)
        text = "\n".join(f"{i:<{w}}  {r}" for i, r in zip(alt.block.ids, rows))
        self.preview.setPlainText(text)

    def _apply(self):
        row = self.alt_list.currentRow()
        if not (self.aln and self.region and 0 <= row < len(self._alts)):
            return
        c0, c1 = self.region
        alt = self._alts[row]
        new_aln = sandbox_mod.splice(self.aln, c0, c1, alt.block)
        hist = list(self.aln.meta.get("history", []))
        hist.append(f"realigned columns {c0+1}–{c1} · {alt.label}")
        new_aln.meta["history"] = hist
        self.applied.emit(new_aln)


class ColumnInspector(QWidget):
    """The known-true answer for one column, in words.

    Colour tells you *that* a column is wrong; this says what the right answer
    was. For the column under the cursor it names the residues the reference
    groups together, which of them the alignment got into place, which residues
    do not belong there and how far they would have to move, and which belong
    there but were put somewhere else. Every column it mentions is a link, so
    following a displaced residue to where it actually sits is one click.

    Note that a column can score a perfect 1.0 and still be listed here as
    incomplete: the score asks whether the pairs a column asserts are true, and a
    column that quietly left a residue out asserts nothing false. Seeing both at
    once is the point.
    """

    #: a column the user asked to be taken to
    jumpRequested = Signal(int)

    _GREEN, _RED, _AMBER, _GREY = "#1a7f37", "#b42318", "#b54708", "#6b7280"

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        self._browser = QTextBrowser()
        self._browser.setOpenLinks(False)          # links are columns, not URLs
        self._browser.anchorClicked.connect(self._on_anchor)
        lay.addWidget(self._browser)
        self.clear("Open or generate an alignment whose true answer is known.")

    # -- rendering ---------------------------------------------------------- #
    def clear(self, why: str) -> None:
        self._browser.setHtml(
            f"<p style='color:{self._GREY}'>{html.escape(why)}</p>")

    def set_column(self, truth: Optional[ColumnTruth]) -> None:
        if truth is None:
            self.clear("No known answer for this alignment. Generate a dataset, "
                       "or load a trusted alignment of the same sequences, from "
                       "the Teach menu.")
            return
        if truth.true_col is None:
            self.clear(f"Column {truth.col + 1} is empty.")
            return

        verdict = ("is exactly right" if truth.correct else "does not match the reference")
        colour = self._GREEN if truth.correct else self._RED
        parts = [
            f"<p style='margin:0 0 8px 0'><b>Column {truth.col + 1}</b> {verdict}."
            f"<br><span style='color:{self._GREY}'>Reference column "
            f"{truth.true_col + 1}.</span></p>"
        ]
        if truth.agree:
            parts.append(self._section(
                "Correctly grouped here", self._GREEN,
                [self._residue(p) for p in truth.agree]))
        if truth.intruders:
            parts.append(self._section(
                "Does not belong here", self._RED,
                [f"{self._residue(p)} — {self._displacement(p)}"
                 for p in truth.intruders]))
        if truth.missing:
            parts.append(self._section(
                "Belongs here but is elsewhere", self._AMBER,
                [f"{self._residue(p)} — currently in {self._link(p.col)}"
                 for p in truth.missing]))
        if truth.correct:
            parts.append(f"<p style='color:{colour};margin:8px 0 0 0'>"
                         "Every residue the reference puts in this column is here, "
                         "and nothing else is.</p>")
        # Residues are named by their position in the *sequence*, columns by
        # their position in the *alignment*, and both are plain numbers — so say
        # once which is which rather than let the reader guess.
        parts.append(f"<p style='color:{self._GREY};margin:12px 0 0 0;"
                     "font-size:small'>A5 means the 5th residue of that sequence. "
                     "Column numbers refer to the alignment.</p>")
        self._browser.setHtml("".join(parts))

    def _section(self, title: str, colour: str, items: Sequence[str]) -> str:
        rows = "".join(f"<li style='margin-bottom:2px'>{it}</li>" for it in items)
        return (f"<p style='color:{colour};margin:8px 0 2px 0'><b>{title}</b></p>"
                f"<ul style='margin:0 0 0 -20px'>{rows}</ul>")

    def _residue(self, p: Placement) -> str:
        return (f"<code>{html.escape(p.seq_id)}</code> "
                f"{html.escape(p.char)}{p.residue + 1}")

    def _displacement(self, p: Placement) -> str:
        if p.offset is None:
            return "its true partners are not in this alignment"
        n = abs(p.offset)
        way = "right" if p.offset > 0 else "left"
        return (f"belongs {n} column{'' if n == 1 else 's'} to the {way}, "
                f"in {self._link(p.target_col)}")

    def _link(self, col: int) -> str:
        return f"<a href='col:{col}'>column {col + 1}</a>"

    # -- links -------------------------------------------------------------- #
    def _on_anchor(self, url) -> None:
        text = url.toString()
        if text.startswith("col:"):
            self.jumpRequested.emit(int(text[4:]))
