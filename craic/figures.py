"""Publication-quality figures for CRAIC.

One vector rendering engine: every figure is painted with a ``QPainter`` onto a
``QPaintDevice``, so the *same* drawing code produces a PNG (raster), an SVG
(vector), or a PDF (vector, paginated) with no extra dependencies — Qt's
``QSvgGenerator`` and ``QPdfWriter`` are both ``QPaintDevice``s.

The headline figure is the *confidence-aware alignment*: a traditional wrapped
block alignment where each residue's colour is faded toward the page when the
column is unreliable, so the trustworthy core stands out and the ambiguous parts
recede. The homology card / arcs / fuzzy logo (see the ``card``/``arcs``/``logo``
functions) reuse the same engine.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence

import numpy as np
from PySide6.QtCore import QRectF, QSize, Qt
from PySide6.QtGui import (
    QColor, QFont, QFontMetricsF, QImage, QPainter, QPainterPath, QPdfWriter, QPen,
)
from PySide6.QtSvg import QSvgGenerator

from .domain import Alignment, Level, column_consensus
from .gui import colors


# --------------------------------------------------------------------------- #
# Themes
# --------------------------------------------------------------------------- #

@dataclass
class Theme:
    bg: QColor
    ink: QColor          # primary text
    muted: QColor        # secondary text / rules
    name: str = "Light"


LIGHT = Theme(QColor("#ffffff"), QColor("#1b2128"), QColor("#8a97a6"), "Light")
DARK = Theme(colors.BG, QColor("#dfe6ee"), colors.MUTED, "Dark")
THEMES = {"Light": LIGHT, "Dark": DARK}


@dataclass
class FigureOptions:
    level: Level = Level.NT
    scheme: str = "clustal"
    theme: str = "Light"
    columns_per_block: int = 60
    cell_w: float = 13.0
    cell_h: float = 15.0
    show_consensus: bool = True
    show_track: bool = True
    title: str = ""
    # per-nucleotide-column reliability in [0, 1] (or None); low → faded cells.
    reliability: Optional[np.ndarray] = None
    table: int = 1


# --------------------------------------------------------------------------- #
# Colour helpers
# --------------------------------------------------------------------------- #

def _blend(a: QColor, b: QColor, t: float) -> QColor:
    t = max(0.0, min(1.0, t))
    return QColor(round(a.red() * (1 - t) + b.red() * t),
                  round(a.green() * (1 - t) + b.green() * t),
                  round(a.blue() * (1 - t) + b.blue() * t))


def _fade(color: QColor, rel: float, bg: QColor) -> QColor:
    """Fade a residue colour toward the page background as reliability drops."""
    if rel != rel:  # nan → treat as fully shown
        return color
    return _blend(bg, color, 0.18 + 0.82 * max(0.0, min(1.0, rel)))


# --------------------------------------------------------------------------- #
# Alignment figure
# --------------------------------------------------------------------------- #

def _agg_reliability(rel: Optional[np.ndarray], nt_start: List[int], span: int,
                     length: int) -> Optional[np.ndarray]:
    if rel is None:
        return None
    out = np.full(len(nt_start), np.nan)
    for c, s in enumerate(nt_start):
        seg = rel[s:min(length, s + span)]
        seg = seg[~np.isnan(seg)] if seg.size else seg
        if seg.size:
            out[c] = float(seg.mean())
    return out


def _consensus(dm, n_seqs: int, level: Level, alphabet) -> List[str]:
    return [column_consensus([dm.cells[r][c] for r in range(n_seqs)], level, alphabet)
            for c in range(dm.n_cols)]


def _draw_alignment(p: QPainter, opt: FigureOptions, aln: Alignment, width: float,
                    new_page: Optional[Callable] = None, page_h: float = 0.0) -> float:
    """Paint the alignment; return the total height used (single-surface devices)."""
    th = THEMES.get(opt.theme, LIGHT)
    dm = aln.display(opt.level)
    rel = _agg_reliability(opt.reliability, dm.nt_start, dm.span, aln.length)
    consensus = _consensus(dm, aln.n_seqs, opt.level, aln.alphabet) if opt.show_consensus else None

    f = QFont("Menlo")
    f.setStyleHint(QFont.TypeWriter)
    f.setPixelSize(int(opt.cell_h * 0.62))
    p.setFont(f)
    fm = QFontMetricsF(f)
    name_w = max((fm.horizontalAdvance(i) for i in aln.ids), default=40.0) + 12.0

    cw, ch = opt.cell_w * dm.span, opt.cell_h
    margin = 24.0
    grid_w = width - 2 * margin - name_w
    per_block = max(1, min(opt.columns_per_block, int(grid_w / cw))) if cw else opt.columns_per_block
    ruler_h = ch
    track_h = ch * 0.6 if opt.show_track else 0.0
    cons_h = ch if opt.show_consensus else 0.0
    block_h = ruler_h + track_h + aln.n_seqs * ch + cons_h + ch  # + gap between blocks

    y = margin
    if opt.title:
        p.setPen(th.ink)
        big = QFont(f); big.setPixelSize(int(opt.cell_h * 0.95)); big.setBold(True)
        p.setFont(big)
        p.drawText(QRectF(margin, y, width - 2 * margin, opt.cell_h * 1.2),
                   Qt.AlignLeft | Qt.AlignVCenter, opt.title)
        p.setFont(f)
        y += opt.cell_h * 1.6

    n_blocks = (dm.n_cols + per_block - 1) // per_block
    for b in range(n_blocks):
        c0 = b * per_block
        c1 = min(dm.n_cols, c0 + per_block)
        if new_page is not None and y + block_h > page_h - margin and b > 0:
            new_page()
            y = margin
        x0 = margin + name_w
        # ruler
        p.setPen(th.muted)
        for c in range(c0, c1):
            ntpos = dm.nt_start[c] + 1
            if ntpos % 10 == 0 or c == c0:
                p.drawText(QRectF(x0 + (c - c0) * cw, y, cw * 4, ruler_h),
                           Qt.AlignLeft | Qt.AlignVCenter, str(ntpos))
        y += ruler_h
        # reliability/conservation track
        if opt.show_track and rel is not None:
            for c in range(c0, c1):
                v = rel[c]
                col = colors.score_color(v) if v == v else QColor("#d9dee4")
                p.fillRect(QRectF(x0 + (c - c0) * cw, y, cw - 0.5, track_h - 1), col)
            y += track_h
        elif opt.show_track:
            y += track_h
        # sequence rows
        for r in range(aln.n_seqs):
            p.setPen(th.ink)
            p.drawText(QRectF(margin, y, name_w - 6, ch),
                       Qt.AlignLeft | Qt.AlignVCenter, aln.ids[r])
            for c in range(c0, c1):
                cell = dm.cells[r][c]
                gap = cell.strip("-") == ""
                rx = x0 + (c - c0) * cw
                if gap:
                    base = QColor("#eef1f4") if th is LIGHT else colors.GAP
                else:
                    base = colors.residue_color(cell, opt.level, aln.alphabet,
                                                opt.table, opt.scheme)
                    if rel is not None and rel[c] == rel[c]:
                        base = _fade(base, rel[c], th.bg)
                p.fillRect(QRectF(rx, y, cw - 0.5, ch - 1), base)
                if not gap and cw >= 9:
                    p.setPen(QColor("#20242a"))
                    p.drawText(QRectF(rx, y, cw - 0.5, ch - 1), Qt.AlignCenter, cell)
            y += ch
        # consensus
        if consensus is not None:
            p.setPen(th.muted)
            p.drawText(QRectF(margin, y, name_w - 6, ch),
                       Qt.AlignLeft | Qt.AlignVCenter, "consensus")
            for c in range(c0, c1):
                ch_c = consensus[c]
                if ch_c:
                    p.setPen(th.ink)
                    p.drawText(QRectF(x0 + (c - c0) * cw, y, cw - 0.5, ch),
                               Qt.AlignCenter, ch_c)
            y += ch
        y += ch  # gap between blocks

    # legend
    p.setPen(th.muted)
    note = (f"{aln.n_seqs} sequences × {aln.length} columns · level: {opt.level.value}"
            + ("  ·  faded = low reliability" if rel is not None else ""))
    p.drawText(QRectF(margin, y, width - 2 * margin, ch), Qt.AlignLeft | Qt.AlignVCenter, note)
    y += ch + margin
    return y


def _figure_height(opt: FigureOptions, aln: Alignment, width: float) -> float:
    dm = aln.display(opt.level)
    cw = opt.cell_w * dm.span
    name_w = 120.0
    margin = 24.0
    grid_w = width - 2 * margin - name_w
    per_block = max(1, min(opt.columns_per_block, int(grid_w / cw))) if cw else opt.columns_per_block
    n_blocks = (dm.n_cols + per_block - 1) // per_block
    ch = opt.cell_h
    block_h = ch + (ch * 0.6 if opt.show_track else 0) + aln.n_seqs * ch \
        + (ch if opt.show_consensus else 0) + ch
    title_h = opt.cell_h * 1.6 if opt.title else 0
    return margin * 2 + title_h + n_blocks * block_h + ch


def _column_consensus(aln: Alignment, level: Level) -> List[str]:
    dm = aln.display(level)
    return [column_consensus([dm.cells[r][c] for r in range(aln.n_seqs)], level, aln.alphabet)
            for c in range(dm.n_cols)]


def homology_card(aln: Alignment, seq_i: int, res_index: int, path: str,
                  fmt: Optional[str] = None, theme: str = "Light",
                  window: int = 120, level: Level = Level.NT) -> str:
    """The 'cloud of potential homologs' for one residue, as a printable card:
    where else (which columns) this residue could align, with what confidence."""
    from .ambiguity import posterior as post
    th = THEMES.get(theme, LIGHT)
    fmt = (fmt or os.path.splitext(path)[1].lstrip(".")).lower()

    probs = post.probe_residue(aln, seq_i, res_index, window=window)
    row_i = aln.rows[seq_i]
    ci = [c for c, ch in enumerate(row_i) if ch != "-"][res_index]   # residue's column
    letter = row_i[ci]
    cand = [c for c in range(len(probs)) if probs[c] == probs[c] and probs[c] > 0.01]
    if ci not in cand:
        cand.append(ci)
    cand.sort()
    cons = _column_consensus(aln, Level.NT)

    cellw, top, barh = 26.0, 70.0, 150.0
    width = max(560.0, 2 * 24 + len(cand) * cellw)
    height = top + barh + 90.0

    def draw(p: QPainter, w: float, h: float):
        p.setPen(th.ink)
        big = QFont("Helvetica"); big.setPixelSize(17); big.setBold(True)
        p.setFont(big)
        p.drawText(QRectF(24, 18, w - 48, 26), Qt.AlignLeft | Qt.AlignVCenter,
                   f"Homology cloud — residue {res_index + 1} ({letter}) of {aln.ids[seq_i]}")
        small = QFont("Helvetica"); small.setPixelSize(11)
        p.setFont(small); p.setPen(th.muted)
        p.drawText(QRectF(24, 44, w - 48, 18), Qt.AlignLeft | Qt.AlignVCenter,
                   "bar height = mean posterior it aligns to that column; current column outlined")
        base = top + barh
        x = 24.0
        for c in cand:
            pv = float(probs[c]) if probs[c] == probs[c] else 0.0
            bh = barh * max(0.0, min(1.0, pv))
            p.fillRect(QRectF(x, base - bh, cellw - 4, bh),
                       colors.colormap_color(pv, "Viridis"))
            if c == ci:
                p.setPen(QColor("#d23b3b")); p.drawRect(QRectF(x - 1, base - barh, cellw - 2, barh))
            p.setPen(th.ink); p.setFont(small)
            p.drawText(QRectF(x, base + 2, cellw - 4, 16), Qt.AlignHCenter, cons[c] or "-")
            p.setPen(th.muted)
            p.drawText(QRectF(x - 6, base + 20, cellw + 8, 14), Qt.AlignHCenter, str(c + 1))
            if pv >= 0.08:
                p.setPen(th.ink)
                p.drawText(QRectF(x - 4, base - bh - 15, cellw + 4, 13), Qt.AlignHCenter,
                           f"{pv:.2f}")
            x += cellw

    return _render_single(path, fmt, width, height, th, draw)


def homology_arcs(aln: Alignment, seq_i: int, seq_j: int, col_start: int, col_stop: int,
                  path: str, fmt: Optional[str] = None, theme: str = "Light",
                  min_post: float = 0.04) -> str:
    """Two sequences as tracks, with arcs linking residues that could be homologous,
    opacity/width ∝ posterior. A confident region shows tight near-vertical arcs;
    an ambiguous one fans out into a 'braid' — alignment uncertainty made visual."""
    from .ambiguity import posterior as post
    th = THEMES.get(theme, LIGHT)
    fmt = (fmt or os.path.splitext(path)[1].lstrip(".")).lower()
    pv = post.region_posterior(aln, seq_i, seq_j, col_start, col_stop)
    M = pv.matrix
    ni, nj = M.shape
    ci, cj = pv.chars_i, pv.chars_j

    cw, cellh, arch = 24.0, 18.0, 150.0
    # left gutter sized to the sequence names so labels never sit under the tracks
    gutter = max(70.0, 8.0 + 7.4 * max(len(aln.ids[seq_i]), len(aln.ids[seq_j])))
    right = 24.0
    top_pad = 48.0
    band = max(ni, nj) * cw                       # the wider of the two tracks
    width = gutter + band + right
    height = top_pad + cellh * 2 + arch + 22.0

    def cell(p, x, y, ch):
        col = colors.residue_color(ch, Level.NT, aln.alphabet, 1, "clustal")
        p.fillRect(QRectF(x, y, cw - 2, cellh), col)
        p.setPen(QColor("#20242a"))
        p.drawText(QRectF(x, y, cw - 2, cellh), Qt.AlignCenter, ch)

    def draw(p: QPainter, w: float, h: float):
        p.setRenderHint(QPainter.Antialiasing, True)
        p.setPen(th.ink)
        title = QFont("Helvetica"); title.setPixelSize(15); title.setBold(True)
        p.setFont(title)
        p.drawText(QRectF(12, 10, w - 24, 22), Qt.AlignLeft,
                   f"Homology arcs — {aln.ids[seq_i]} vs {aln.ids[seq_j]}  "
                   f"(cols {col_start + 1}–{col_stop})")
        topy = top_pad
        boty = topy + cellh + arch
        xi0 = gutter + (band - ni * cw) / 2       # centre each track within the band,
        xj0 = gutter + (band - nj * cw) / 2       # to the right of the name gutter
        for a in range(ni):
            for b in range(nj):
                v = float(M[a, b])
                if v < min_post:
                    continue
                xa = xi0 + a * cw + cw / 2
                xb = xj0 + b * cw + cw / 2
                midy = (topy + cellh + boty) / 2
                pth = QPainterPath()
                pth.moveTo(xa, topy + cellh)
                pth.cubicTo(xa, midy, xb, midy, xb, boty)
                col = colors.colormap_color(v, "Viridis")
                col.setAlphaF(min(1.0, 0.12 + 0.88 * v))
                pen = QPen(col); pen.setWidthF(0.7 + 3.0 * v)
                p.setPen(pen); p.setBrush(Qt.NoBrush); p.drawPath(pth)
        p.setFont(QFont("Menlo", 9))
        for a, ch in enumerate(ci):
            cell(p, xi0 + a * cw, topy, ch)
        for b, ch in enumerate(cj):
            cell(p, xj0 + b * cw, boty, ch)
        p.setPen(th.muted); p.setFont(QFont("Helvetica", 9))
        p.drawText(QRectF(8, topy, gutter - 12, cellh), Qt.AlignLeft | Qt.AlignVCenter, aln.ids[seq_i])
        p.drawText(QRectF(8, boty, gutter - 12, cellh), Qt.AlignLeft | Qt.AlignVCenter, aln.ids[seq_j])

    return _render_single(path, fmt, width, height, th, draw)


def uncertainty_logo(aln: Alignment, col_start: int, col_stop: int, path: str,
                     reliability: Optional[np.ndarray] = None, fmt: Optional[str] = None,
                     theme: str = "Light", level: Level = Level.NT) -> str:
    """A sequence logo whose per-column stack height is scaled by reliability, so
    confidently-aligned columns stand tall and ambiguous ones shrink/fade."""
    from collections import Counter
    th = THEMES.get(theme, LIGHT)
    fmt = (fmt or os.path.splitext(path)[1].lstrip(".")).lower()
    dm = aln.display(level)
    # col_start/col_stop are alignment (nt) coordinates; pick the display columns
    # that fall in that range so callers can pass the on-screen selection directly.
    cols = [c for c in range(dm.n_cols) if col_start <= dm.nt_start[c] < col_stop]
    if not cols:
        cols = list(range(dm.n_cols))
    rel = _agg_reliability(reliability, dm.nt_start, dm.span, aln.length)

    cw, maxh = 22.0, 150.0
    width = max(360.0, 2 * 24 + len(cols) * cw)
    height = maxh + 70.0

    def draw(p: QPainter, w: float, h: float):
        base = 30.0 + maxh
        p.setPen(th.ink)
        title = QFont("Helvetica"); title.setPixelSize(15); title.setBold(True)
        p.setFont(title)
        p.drawText(QRectF(24, 8, w - 48, 20), Qt.AlignLeft, "Uncertainty-weighted logo")
        x = 24.0
        glyph = QFont("Helvetica"); glyph.setBold(True)
        for c in cols:
            freqs = Counter(dm.cells[r][c] for r in range(aln.n_seqs)
                            if dm.cells[r][c].strip("-"))
            tot = sum(freqs.values())
            scale = (rel[c] if (rel is not None and rel[c] == rel[c]) else 1.0)
            y = base
            if tot:
                for sym, n in sorted(freqs.items(), key=lambda kv: kv[1]):
                    hh = maxh * (n / tot) * scale
                    if hh < 1.2:
                        continue
                    col = colors.residue_color(sym, level, aln.alphabet, 1, "clustal")
                    p.save()
                    p.translate(x, y - hh)
                    p.scale((cw - 3) / 10.0, hh / 12.0)   # stretch glyph into its slot
                    p.setFont(glyph); p.setPen(col)
                    glyph.setPixelSize(12); p.setFont(glyph)
                    p.drawText(QRectF(0, 0, 10, 12), Qt.AlignCenter, sym)
                    p.restore()
                    y -= hh
            p.setPen(th.muted); p.setFont(QFont("Helvetica", 7))
            p.drawText(QRectF(x - 4, base + 4, cw + 4, 14), Qt.AlignHCenter, str(dm.nt_start[c] + 1))
            x += cw

    return _render_single(path, fmt, width, height, th, draw)


def export_html(aln: Alignment, path: str, reliability: Optional[np.ndarray] = None,
                level: Level = Level.NT, scheme: str = "clustal",
                title: str = "CRAIC alignment", table: int = 1) -> str:
    """A self-contained, shareable HTML view: coloured residues, a reliability
    strip, sticky names, and hover tooltips (column position + reliability)."""
    dm = aln.display(level)
    rel = _agg_reliability(reliability, dm.nt_start, dm.span, aln.length)

    def hx(qc: QColor) -> str:
        return f"#{qc.red():02x}{qc.green():02x}{qc.blue():02x}"

    css = (
        "body{font-family:-apple-system,Helvetica,Arial,sans-serif;background:#fff;"
        "color:#1b2128;margin:18px;}"
        ".aln{font-family:Menlo,monospace;font-size:13px;white-space:nowrap;overflow:auto;"
        "border:1px solid #e3e7ec;padding:4px;}"
        ".r{display:block;height:16px;}"
        ".nm{display:inline-block;width:120px;color:#555;position:sticky;left:0;background:#fff;}"
        ".c{display:inline-block;width:14px;text-align:center;}"
        ".trk{display:inline-block;width:14px;height:9px;}"
        "h1{font-size:16px;margin:0 0 4px;} .lg{color:#8a97a6;font-size:12px;margin:6px 0 10px;}"
    )
    out = [f"<!doctype html><html><head><meta charset='utf-8'><title>{title}</title>",
           f"<style>{css}</style></head><body><h1>{title}</h1>",
           f"<div class='lg'>{aln.n_seqs} sequences × {aln.length} columns · "
           "hover a cell for its column &amp; reliability</div><div class='aln'>"]
    if rel is not None:
        cells = ["<span class='nm'>reliability</span>"]
        for c in range(dm.n_cols):
            v = rel[c]
            col = colors.score_color(v) if v == v else QColor("#d9dee4")
            tip = f"col {dm.nt_start[c] + 1}: {v:.2f}" if v == v else f"col {dm.nt_start[c] + 1}"
            cells.append(f"<span class='trk' style='background:{hx(col)}' title='{tip}'></span>")
        out.append("<span class='r'>" + "".join(cells) + "</span>")
    for r in range(aln.n_seqs):
        cells = [f"<span class='nm'>{aln.ids[r]}</span>"]
        for c in range(dm.n_cols):
            cell = dm.cells[r][c]
            gap = cell.strip("-") == ""
            bg = "#eef1f4" if gap else hx(colors.residue_color(cell, level, aln.alphabet, table, scheme))
            tip = f"col {dm.nt_start[c] + 1}"
            if rel is not None and rel[c] == rel[c]:
                tip += f" · rel {rel[c]:.2f}"
            cells.append(f"<span class='c' style='background:{bg}' title='{tip}'>"
                         f"{'-' if gap else cell}</span>")
        out.append("<span class='r'>" + "".join(cells) + "</span>")
    out.append("</div></body></html>")
    with open(path, "w") as fh:
        fh.write("".join(out))
    return path


def confidence_report(aln: Alignment, col_reliability: np.ndarray, path: str,
                      fmt: Optional[str] = None, theme: str = "Light",
                      threshold: float = 0.5, title: str = "Alignment confidence report") -> str:
    """A per-column reliability profile with ambiguous hotspots shaded, plus a
    text QC summary — the report you'd attach to a submission."""
    from .ambiguity import report as rep_mod
    th = THEMES.get(theme, LIGHT)
    fmt = (fmt or os.path.splitext(path)[1].lstrip(".")).lower()
    rel = np.asarray(col_reliability, dtype=float)
    s = rep_mod.summary(rel, threshold)
    n = max(1, len(rel))

    margin = 28.0
    plot_w = 1040.0
    bar_area = plot_w - 2 * margin
    bw = bar_area / n
    maxh = 150.0
    text_lines = s.as_text().splitlines()
    text_h = 18.0 * (len(text_lines) + 1)
    width = plot_w
    height = margin + 26 + maxh + 28 + text_h + margin

    def draw(p: QPainter, w: float, h: float):
        p.setPen(th.ink)
        t = QFont("Helvetica"); t.setPixelSize(16); t.setBold(True); p.setFont(t)
        p.drawText(QRectF(margin, 10, w - 2 * margin, 22), Qt.AlignLeft, title)
        base = 38 + maxh
        # hotspot bands
        for hs in s.hotspots:
            x = margin + hs.col_start * bw
            p.fillRect(QRectF(x, 38, max(1.0, hs.width * bw), maxh), QColor(214, 59, 59, 38))
        # threshold line
        p.setPen(QColor(150, 150, 160))
        ty = base - threshold * maxh
        p.drawLine(int(margin), int(ty), int(margin + bar_area), int(ty))
        # reliability bars
        for c in range(n):
            v = rel[c]
            if v != v:
                continue
            bh = max(0.0, min(1.0, v)) * maxh
            p.fillRect(QRectF(margin + c * bw, base - bh, max(0.6, bw - 0.3), bh),
                       colors.score_color(v))
        p.setPen(th.muted); sf = QFont("Helvetica"); sf.setPixelSize(10); p.setFont(sf)
        p.drawText(QRectF(margin, base + 4, 200, 14), Qt.AlignLeft, "column → ")
        p.drawText(QRectF(w - margin - 120, base + 4, 120, 14), Qt.AlignRight, "red = hotspot")
        # summary text
        p.setPen(th.ink); mono = QFont("Menlo"); mono.setPixelSize(11); p.setFont(mono)
        y = base + 28
        for line in text_lines:
            p.drawText(QRectF(margin, y, w - 2 * margin, 16), Qt.AlignLeft | Qt.AlignVCenter, line)
            y += 16

    return _render_single(path, fmt, width, height, th, draw)


def matrix_heatmap(M: np.ndarray, row_labels: Sequence[str], col_labels: Sequence[str],
                   path: str, fmt: Optional[str] = None, theme: str = "Light",
                   title: str = "", cmap: str = "Viridis", subtitle: str = "") -> str:
    """A labelled heatmap (used for the ensemble co-occurrence map)."""
    th = THEMES.get(theme, LIGHT)
    fmt = (fmt or os.path.splitext(path)[1].lstrip(".")).lower()
    M = np.asarray(M, dtype=float)
    ni, nj = (M.shape + (0, 0))[:2] if M.ndim == 2 else (0, 0)
    cell, left, top = 18.0, 40.0, 56.0
    width = max(360.0, left + nj * cell + 30)
    height = top + ni * cell + 54

    def draw(p: QPainter, w: float, h: float):
        p.setPen(th.ink)
        t = QFont("Helvetica"); t.setPixelSize(15); t.setBold(True); p.setFont(t)
        p.drawText(QRectF(20, 8, w - 40, 20), Qt.AlignLeft, title or "Heatmap")
        if subtitle:
            p.setPen(th.muted); s = QFont("Helvetica"); s.setPixelSize(11); p.setFont(s)
            p.drawText(QRectF(20, 28, w - 40, 16), Qt.AlignLeft, subtitle)
        gl = QFont("Menlo"); gl.setPixelSize(10)
        for a in range(ni):
            for b in range(nj):
                p.fillRect(QRectF(left + b * cell, top + a * cell, cell - 0.5, cell - 0.5),
                           colors.colormap_color(float(M[a, b]), cmap))
            p.setPen(th.ink); p.setFont(gl)
            if a < len(row_labels):
                p.drawText(QRectF(2, top + a * cell, left - 6, cell), Qt.AlignRight | Qt.AlignVCenter,
                           row_labels[a])
        p.setPen(th.ink); p.setFont(gl)
        for b in range(nj):
            if b < len(col_labels):
                p.drawText(QRectF(left + b * cell, top - 14, cell, 12), Qt.AlignHCenter, col_labels[b])

    return _render_single(path, fmt, width, height, th, draw)


def _render_single(path: str, fmt: str, width: float, height: float, theme: Theme,
                   draw: Callable[[QPainter, float, float], None]) -> str:
    """Paint a single-surface figure (no pagination) to png/svg/pdf."""
    if fmt == "svg":
        gen = QSvgGenerator()
        gen.setFileName(path)
        gen.setSize(QSize(int(width), int(height)))
        gen.setViewBox(QRectF(0, 0, width, height))
        p = QPainter(gen)
        p.fillRect(QRectF(0, 0, width, height), theme.bg)
        draw(p, width, height)
        p.end()
    elif fmt == "pdf":
        from PySide6.QtCore import QMarginsF, QSizeF
        from PySide6.QtGui import QPageSize
        w = QPdfWriter(path)
        w.setResolution(150)
        w.setPageSize(QPageSize(QSizeF(width / 150.0, height / 150.0), QPageSize.Inch))
        w.setPageMargins(QMarginsF(0, 0, 0, 0))
        p = QPainter(w)
        pw, ph = float(w.width()), float(w.height())
        p.fillRect(QRectF(0, 0, pw, ph), theme.bg)
        p.scale(pw / width, ph / height)
        draw(p, width, height)
        p.end()
    else:
        img = QImage(int(width), int(height), QImage.Format_ARGB32)
        img.fill(theme.bg)
        p = QPainter(img)
        draw(p, width, height)
        p.end()
        img.save(path, "PNG")
    return path


def export_alignment(aln: Alignment, path: str, fmt: Optional[str] = None,
                     opt: Optional[FigureOptions] = None, width: float = 1100.0) -> str:
    """Render the alignment figure to ``path`` as svg / pdf / png (by extension)."""
    opt = opt or FigureOptions()
    fmt = (fmt or os.path.splitext(path)[1].lstrip(".")).lower()
    th = THEMES.get(opt.theme, LIGHT)
    height = _figure_height(opt, aln, width)

    if fmt == "pdf":   # paginated across A4 pages
        from PySide6.QtCore import QMarginsF
        from PySide6.QtGui import QPageSize
        w = QPdfWriter(path)
        w.setResolution(150)
        w.setPageSize(QPageSize(QPageSize.A4))
        w.setPageMargins(QMarginsF(0, 0, 0, 0))
        p = QPainter(w)
        pw, ph = float(w.width()), float(w.height())
        p.fillRect(QRectF(0, 0, pw, ph), th.bg)
        _draw_alignment(p, opt, aln, pw, new_page=w.newPage, page_h=ph)
        p.end()
        return path
    _render_single(path, fmt, width, height, th,
                   lambda p, w, h: _draw_alignment(p, opt, aln, w))
    return path
