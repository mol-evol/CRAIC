#!/usr/bin/env python3
"""Regenerate every screenshot in the documentation.

Screenshots rot. A control is renamed, a panel gains a tab, the colour ramp
changes, and the images in the docs quietly start describing a version that no
longer exists — usually noticed by a reader rather than by us. So the images are
not taken by hand: this script drives the real window over fixed data and writes
them all, and it is re-run whenever the interface changes.

    python scripts/make_screenshots.py            # -> docs/assets/

It runs offscreen (``QT_QPA_PLATFORM=offscreen``), so it needs no display and
works in CI. The data is seeded, so re-running it produces the same pictures and
a diff means something really moved.
"""

from __future__ import annotations

import os
import sys
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from PySide6.QtWidgets import QApplication                         # noqa: E402

import numpy as np                                                 # noqa: E402

from craic import simulate                                         # noqa: E402
from craic.domain import Alignment, Alphabet, Level                # noqa: E402
from craic.engines import BuiltinProgressive                       # noqa: E402
from craic.gui.app import (CraicWindow, _COLOUR_CONFIDENCE,        # noqa: E402
                           _COLOUR_TRUTH, _TRACKS, _TRUTH_TRACK)

OUT = os.path.join(REPO, "docs", "assets")
WIDTH = 1400
#: Bigger than the default so the residue letters survive being scaled into a page.
CELL = (13, 17)
#: One dataset for every teaching shot, so the pictures agree with each other and
#: a reader can follow one example across pages.
SEED = 18


def _app():
    return QApplication.instance() or QApplication([])


def silence_dialogs():
    """Answer every modal dialog before it opens.

    Some of what this script drives is allowed to ask a question — marking a
    nucleotide alignment as coding warns about broken reading frames, for
    instance. With no display there is nobody to click OK, so the dialog blocks
    the script for ever. Their answers do not change the pictures, so they are
    stubbed out here rather than tiptoed around in the window code.
    """
    from PySide6.QtWidgets import QMessageBox
    for name, answer in (("information", QMessageBox.Ok), ("warning", QMessageBox.Ok),
                         ("critical", QMessageBox.Ok), ("question", QMessageBox.Yes)):
        setattr(QMessageBox, name, staticmethod(lambda *a, _r=answer, **k: _r))


def _settle(app, n: int = 4):
    """Let Qt lay out and paint before grabbing, or the shot is half-drawn."""
    for _ in range(n):
        app.processEvents()
        time.sleep(0.02)


def wait_until(app, predicate, timeout: float = 180.0, what: str = "a background task"):
    """Pump the event loop until ``predicate`` holds.

    The analyses run on worker threads and deliver their results through
    signals, so a fixed number of ``processEvents`` calls is a race: it passes on
    a fast machine and silently captures an empty track on a slow one. Waiting on
    the condition instead makes the script's output the same either way.
    """
    deadline = time.time() + timeout
    while time.time() < deadline:
        app.processEvents()
        if predicate():
            _settle(app)
            return
        time.sleep(0.05)
    raise SystemExit(f"timed out after {timeout:.0f}s waiting for {what}")


def fit(win, rows: int, height: int = 0) -> None:
    """Size the window to the data, so a shot of eight sequences is not mostly
    empty canvas."""
    win.resize(WIDTH, height or min(820, 190 + rows * win.canvas.cell_h))


def shot(app, win, name: str):
    _settle(app)
    path = os.path.join(OUT, name)
    win.grab().save(path)
    print(f"  {name}", flush=True)


def track_index(key_label: str) -> int:
    for i, label in enumerate(_TRACKS):
        if label.lower().startswith(key_label.lower()):
            return i
    raise SystemExit(f"no track called {key_label!r}; have {_TRACKS}")


def simulated_pair(seed: int = SEED):
    """An unaligned dataset and its known-true alignment."""
    # Divergence and indel load chosen so the built-in aligner gets most of it
    # right and a clear minority wrong (SP ~0.76, TC ~0.51). A dataset it aligns
    # perfectly makes an all-green picture that teaches nothing, and one it
    # destroys is just as useless.
    d = simulate.simulate(taxa=7, root_len=150, bmax=0.9, indel_rate=0.1,
                          rate_alpha=0.8, seed=seed)
    truth = Alignment(list(d["names"]), list(d["true_rows"]), Alphabet.DNA)
    return d, truth


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    app = _app()
    silence_dialogs()
    from craic import io

    # ---- the plain viewer, on the shipped coding example ------------------- #
    # This one is deliberately easy: it is showing what the viewer looks like,
    # not what the analyses find.
    print("viewer", flush=True)
    win = CraicWindow()
    win.canvas.set_cell_size(*CELL)
    win.show()
    example = os.path.join(REPO, "examples", "coding_genes.fasta")
    aln = BuiltinProgressive().align(io.read_records(example), Alphabet.DNA)
    win._set_alignment(aln)
    win.canvas.set_cell_size(*CELL)
    fit(win, aln.n_seqs)
    win.track_combo.setCurrentIndex(track_index("conservation"))
    shot(app, win, "viewer-nucleotide.png")

    win.coding_chk.setChecked(True)
    levels = [win.level_combo.itemData(i) for i in range(win.level_combo.count())]
    if Level.AA in levels:
        win.level_combo.setCurrentIndex(levels.index(Level.AA))
        shot(app, win, "viewer-amino-acid.png")
    win.close()

    # ---- everything else, on a dataset with real ambiguity in it ----------- #
    # An easy alignment makes an all-green reliability track and a clean
    # posterior diagonal: pictures that show the controls but none of the work.
    print("analysis", flush=True)
    d, truth = simulated_pair()
    win = CraicWindow()
    win.canvas.set_cell_size(*CELL)
    win.show()
    aligned = BuiltinProgressive().align(d["seqs"], Alphabet.DNA)
    win._set_alignment(aligned)
    win.canvas.set_cell_size(*CELL)
    fit(win, aligned.n_seqs)

    win.track_combo.setCurrentIndex(track_index("reliability"))
    wait_until(app, lambda: win.canvas.nt_scores is not None,
               what="the reliability analysis")
    shot(app, win, "reliability-track.png")
    win.color_combo.setCurrentIndex(_COLOUR_CONFIDENCE)
    wait_until(app, lambda: win.canvas.color_mode == "confidence",
               what="confidence colouring")
    shot(app, win, "confidence-colouring.png")

    print("masking", flush=True)
    win.thr.setValue(60)
    shot(app, win, "masking-preview.png")
    win.thr.setValue(0)
    win.color_combo.setCurrentIndex(0)

    print("ambiguity tools", flush=True)
    win.dock.show()
    fit(win, aligned.n_seqs, height=480)
    win.canvas.select_nt_range(40, 110)
    win.tabs.setCurrentWidget(win.sandbox)
    # Actually generate the alternatives: a shot of an empty list shows the panel
    # but not the thing the panel is for.
    win.sandbox._generate()
    wait_until(app, lambda: win.sandbox.alt_list.count() > 0,
               timeout=240, what="the realignment alternatives")
    win.sandbox._close_progress()
    _settle(app, 6)
    shot(app, win, "realignment-sandbox.png")
    win.tabs.setCurrentWidget(win.posterior)
    wait_until(app, lambda: win.posterior.heatmap.view is not None,
               timeout=120, what="the posterior matrix")
    shot(app, win, "posterior-explorer.png")
    win.dock.hide()
    win.canvas.select_nt_range(0, 0)
    fit(win, aligned.n_seqs)

    # ---- teaching: the same dataset, with its answer attached -------------- #
    print("teaching", flush=True)
    win._set_reference(truth, "simulated truth")
    win.track_combo.setCurrentIndex(_TRUTH_TRACK)
    shot(app, win, "truth-track.png")

    win.color_combo.setCurrentIndex(_COLOUR_TRUTH)
    shot(app, win, "truth-colouring.png")

    wrong = [c for c in range(win.doc.alignment.length)
             if np.nan_to_num(win.doc.truth.col_correct[c], nan=1.0) < 0.999]

    def instructive(col) -> bool:
        """One residue in the wrong place, with somewhere definite to go, and a
        gap where it should have been — the clearest case to read."""
        e = win.doc.explain_column(col)
        return (len(e.intruders) == 1 and e.intruders[0].target_col is not None
                and len(e.missing) >= 1 and len(e.agree) >= 3)

    pick = next((c for c in wrong if instructive(c)),
                next((c for c in wrong
                      if any(p.target_col is not None
                             for p in win.doc.explain_column(c).intruders)), None))
    if pick is not None:
        win.dock.show()
        fit(win, aligned.n_seqs, height=480)
        win.tabs.setCurrentWidget(win.inspector)
        win._jump_to_column(pick)
        shot(app, win, "column-inspector.png")
        win.dock.hide()
        fit(win, aligned.n_seqs)

    win.show_truth_act.setChecked(True)
    shot(app, win, "true-alignment.png")
    win.show_truth_act.setChecked(False)

    print(f"\nwrote to {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
