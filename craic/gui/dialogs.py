"""The dialogs.

Small, self-contained forms: aligner parameters, the figure exporter, and the
teaching-dataset generator. They were in ``app.py``, which is already the largest
file in the project and does not need to also hold three unrelated widgets that
nothing else in the window touches.
"""

from __future__ import annotations

from PySide6.QtGui import QDoubleValidator, QIntValidator
from PySide6.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QLabel,
    QLineEdit, QSpinBox,
)

from .panels import exclusive_pair


def _fmt_params(params: dict) -> str:
    return ", ".join(f"{k}={v}" for k, v in params.items())


class SimulateDialog(QDialog):
    """Generate a dataset whose true alignment is known exactly.

    The point of this in a workbench is teaching, not benchmarking: with the
    truth in hand a student can align the sequences, turn on truth mode, and see
    for themselves which columns the aligner got wrong and whether the
    reliability score knew. Nothing else in the tool can show that, because with
    real data nobody knows the answer.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Generate a dataset with a known answer")
        self.setMinimumWidth(460)
        form = QFormLayout(self)

        self.taxa = QSpinBox(); self.taxa.setRange(3, 200); self.taxa.setValue(8)
        form.addRow("Sequences", self.taxa)
        self.length = QSpinBox(); self.length.setRange(30, 5000); self.length.setValue(200)
        self.length.setSingleStep(50)
        form.addRow("Root length (nt)", self.length)

        self.divergence = QDoubleSpinBox()
        self.divergence.setRange(0.02, 1.5); self.divergence.setSingleStep(0.05)
        self.divergence.setValue(0.28); self.divergence.setDecimals(2)
        self.divergence.setToolTip(
            "Upper bound on branch lengths. Higher is more divergent and harder "
            "to align; the reliability signal is strongest in the middle of this "
            "range, where there is real but recoverable error.")
        form.addRow("Divergence", self.divergence)

        self.indel = QDoubleSpinBox()
        self.indel.setRange(0.0, 20.0); self.indel.setSingleStep(0.5); self.indel.setValue(2.0)
        self.indel.setToolTip("Insertion/deletion rate. Indels are what makes "
                              "alignment ambiguous, so this is the dial that "
                              "creates hard regions.")
        form.addRow("Indel rate", self.indel)

        self.rate_alpha = QDoubleSpinBox()
        self.rate_alpha.setRange(0.0, 5.0); self.rate_alpha.setSingleStep(0.1)
        self.rate_alpha.setValue(0.0); self.rate_alpha.setDecimals(2)
        self.rate_alpha.setToolTip(
            "Gamma shape for among-site rate variation; 0 means every site "
            "evolves at the same rate. Uniform rates make the data closer to "
            "what CRAIC's own model assumes, so a non-zero value here is the "
            "harder and more realistic test.")
        form.addRow("Rate variation (alpha)", self.rate_alpha)

        self.seed = QSpinBox(); self.seed.setRange(0, 10**6); self.seed.setValue(0)
        self.seed.setToolTip("Same seed, same dataset — so a class can all work "
                             "on identical data, or you can return to one.")
        form.addRow("Seed", self.seed)

        note = QLabel("The unaligned sequences are loaded ready to align. The true "
                      "alignment is kept as the reference, so truth mode works "
                      "straight away.")
        note.setWordWrap(True)
        note.setStyleSheet("color: gray;")
        form.addRow("", note)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.button(QDialogButtonBox.Ok).setText("Generate")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

    def values(self) -> dict:
        return {"taxa": self.taxa.value(), "root_len": self.length.value(),
                "seed": self.seed.value(), "indel_rate": self.indel.value(),
                "bmax": self.divergence.value(), "rate_alpha": self.rate_alpha.value()}


def _fmt_bound(v) -> str:
    """A bound as a user would write it: 1000, not 1000.0."""
    if v is None:
        return "…"
    return str(int(v)) if float(v).is_integer() else str(v)


def _clamp(value, p):
    """Hold a value inside the parameter's bounds.

    The validator stops most bad input at the keystroke, but paste and
    programmatic edits can still get past it, and an out-of-range value costs the
    user a failed alignment rather than a warning.
    """
    if p.lo is not None:
        value = max(value, type(value)(p.lo))
    if p.hi is not None:
        value = min(value, type(value)(p.hi))
    return value


class EngineParamsDialog(QDialog):
    """A small form to edit one aligner's parameters."""

    def __init__(self, engine, values, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"{engine.label} — parameters")
        self._widgets = {}
        form = QFormLayout(self)
        params = engine.parameters()
        if not params:
            form.addRow(QLabel("This aligner has no tunable parameters in CRAIC."))
        for p in params:
            v = values.get(p.key, p.default)
            if p.kind == "choice":
                w = QComboBox()
                w.addItems(p.choices or [])
                if v in (p.choices or []):
                    w.setCurrentText(v)
            else:
                # Plain text field (with a numeric validator) rather than a spin
                # box — lets you type any value directly, on any platform.
                w = QLineEdit(str(v))
                # A parameter that names bounds gets them: the aligner will
                # reject an out-of-range value anyway, and it is better to refuse
                # the keystroke than to fail after the run has started.
                if p.kind == "int":
                    lo = int(p.lo) if p.lo is not None else -1000000
                    hi = int(p.hi) if p.hi is not None else 1000000
                    w.setValidator(QIntValidator(lo, hi, w))
                else:
                    lo = float(p.lo) if p.lo is not None else -1e6
                    hi = float(p.hi) if p.hi is not None else 1e6
                    val = QDoubleValidator(lo, hi, 6, w)
                    val.setNotation(QDoubleValidator.StandardNotation)
                    w.setValidator(val)
                if p.lo is not None or p.hi is not None:
                    w.setPlaceholderText(f"{_fmt_bound(p.lo)} – {_fmt_bound(p.hi)}")
            if p.help:
                w.setToolTip(p.help)
            self._widgets[p.key] = (p, w)
            form.addRow(p.label, w)
        buttons = QDialogButtonBox(
            QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self._restore)
        form.addRow(buttons)

    def _restore(self):
        for _key, (p, w) in self._widgets.items():
            if p.kind == "choice":
                w.setCurrentText(str(p.default))
            else:
                w.setText(str(p.default))

    def values(self) -> dict:
        out = {}
        for key, (p, w) in self._widgets.items():
            if p.kind == "choice":
                out[key] = w.currentText()
            elif p.kind == "int":
                try:
                    out[key] = _clamp(int(w.text()), p)
                except ValueError:
                    out[key] = int(p.default)
            else:
                try:
                    out[key] = _clamp(float(w.text()), p)
                except ValueError:
                    out[key] = float(p.default)
        return out


_FIG_KINDS = [
    ("Alignment figure", "alignment"),
    ("Confidence report", "report"),
    ("Uncertainty logo", "logo"),
    ("Homology arcs (two sequences)", "arcs"),
    ("Ensemble co-occurrence (two sequences)", "cooc"),
    ("Homology card (one residue)", "card"),
    ("Interactive HTML", "html"),
]

_FIG_HELP = {
    "alignment": "The whole alignment, coloured and (if computed) shaded by confidence.",
    "report": "Per-column reliability profile with ambiguous hotspots and a QC summary.",
    "logo": "A sequence logo whose columns shrink/fade where the alignment is uncertain.",
    "arcs": "Two sequences as tracks, with arcs linking residues that could be homologous.",
    "cooc": "Heatmap of how often two sequences' residues co-align across alternative alignments.",
    "card": "The 'cloud' of columns one residue could align to, with confidence.",
    "html": "A self-contained, shareable HTML view with hover tooltips.",
}


class FigureExportDialog(QDialog):
    """Ask the user exactly what figure to export, and over which sequences /
    region — so paired figures no longer silently use the first two sequences."""

    def __init__(self, aln, seq_i, seq_j, probe, has_selection, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Export figure")
        self.setMinimumWidth(420)
        self._aln = aln
        form = QFormLayout(self)

        self.kind = QComboBox()
        for label, key in _FIG_KINDS:
            self.kind.addItem(label, key)
        form.addRow("Figure", self.kind)

        self.help = QLabel()
        self.help.setWordWrap(True)
        self.help.setStyleSheet("color: gray;")
        form.addRow("", self.help)

        ids = list(aln.ids)
        self.seq_i = QComboBox(); self.seq_i.addItems(ids)
        self.seq_j = QComboBox(); self.seq_j.addItems(ids)
        self.seq_i.setCurrentIndex(min(seq_i, len(ids) - 1))
        self.seq_j.setCurrentIndex(min(seq_j, len(ids) - 1))
        # A figure of a sequence against itself is the trivial diagonal, so the
        # two pickers exclude each other rather than producing one on request.
        exclusive_pair(self.seq_i, self.seq_j)
        self.seq_i.currentIndexChanged.connect(
            lambda *_: exclusive_pair(self.seq_i, self.seq_j))
        self.seq_j.currentIndexChanged.connect(
            lambda *_: exclusive_pair(self.seq_j, self.seq_i))
        form.addRow("Sequence A", self.seq_i)
        form.addRow("Sequence B", self.seq_j)

        self.residue = QSpinBox(); self.residue.setMinimum(1)
        self.seq_card = probe[0] if probe else 0
        self._card_res = (probe[1] + 1) if probe else 1
        form.addRow("Residue #", self.residue)

        self.region = QComboBox()
        self.region.addItem("Whole alignment", False)
        if has_selection:
            self.region.addItem("Current selection", True)
            self.region.setCurrentIndex(1)
        form.addRow("Region", self.region)

        self.theme = QComboBox(); self.theme.addItems(["Light", "Dark"])
        form.addRow("Theme", self.theme)

        self.fmt = QComboBox(); self.fmt.addItems(["pdf", "svg", "png"])
        form.addRow("Format", self.fmt)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

        self.kind.currentIndexChanged.connect(self._sync)
        self._sync()

    def _sync(self):
        kind = self.kind.currentData()
        self.help.setText(_FIG_HELP.get(kind, ""))
        paired = kind in ("arcs", "cooc")
        card = kind == "card"
        self.seq_i.setEnabled(paired or card)
        self.seq_j.setEnabled(paired)
        self.residue.setEnabled(card)
        if card:
            self.seq_i.setCurrentIndex(min(self.seq_card, self.seq_i.count() - 1))
            row = self._aln.rows[self.seq_i.currentIndex()]
            self.residue.setMaximum(max(1, len(row.replace("-", ""))))
            self.residue.setValue(min(self._card_res, self.residue.maximum()))
        self.region.setEnabled(kind in ("logo", "arcs", "cooc"))
        self.fmt.clear()
        self.fmt.addItems(["html"] if kind == "html" else ["pdf", "svg", "png"])

    def values(self) -> dict:
        return {
            "kind": self.kind.currentData(),
            "seq_i": self.seq_i.currentIndex(),
            "seq_j": self.seq_j.currentIndex(),
            "residue": self.residue.value() - 1,
            "selection": bool(self.region.currentData()),
            "theme": self.theme.currentText(),
            "fmt": self.fmt.currentText(),
        }
