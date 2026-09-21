"""Manual alignment editing — pure, rectangular-preserving operations.

Every function returns a *new* ``Alignment`` (or ``None`` if the edit is not
applicable), which makes undo/redo trivial: the GUI keeps a stack of alignments.
Edits that change column structure drop the coding annotation (a hand-made gap can
break the reading frame); residue nudges keep it, since they don't move columns.
"""

from __future__ import annotations

from typing import List, Optional, Sequence

from .domain import Alignment


def _aln(aln: Alignment, rows: List[str], keep_coding: bool = False) -> Alignment:
    return Alignment(list(aln.ids), rows, aln.alphabet,
                     coding=aln.coding if keep_coding else None)


def move_residue(aln: Alignment, row: int, col: int, direction: int) -> Optional[Alignment]:
    """Slide the residue at (row, col) into an adjacent gap (direction -1/+1)."""
    if not (0 <= row < aln.n_seqs):
        return None
    s = aln.rows[row]
    j = col + direction
    if not (0 <= col < len(s) and 0 <= j < len(s)):
        return None
    if s[col] == "-" or s[j] != "-":     # need a residue to move and a gap to move into
        return None
    chars = list(s)
    chars[col], chars[j] = chars[j], chars[col]
    rows = list(aln.rows)
    rows[row] = "".join(chars)
    return _aln(aln, rows, keep_coding=True)   # columns unchanged → frame intact


def insert_gap_column(aln: Alignment, col: int) -> Alignment:
    """Open a gap in every sequence at ``col`` (alignment grows by one column)."""
    col = max(0, min(col, aln.length))
    return _aln(aln, [r[:col] + "-" + r[col:] for r in aln.rows])


def delete_column(aln: Alignment, col: int) -> Optional[Alignment]:
    """Remove column ``col`` only if it is all gaps."""
    if not (0 <= col < aln.length):
        return None
    if any(r[col] != "-" for r in aln.rows):
        return None
    return _aln(aln, [r[:col] + r[col + 1:] for r in aln.rows])


def open_gap(aln: Alignment, row: int, col: int) -> Optional[Alignment]:
    """Open a gap at (row, col) in one sequence, absorbing a trailing gap to keep
    the alignment rectangular; if the row has no trailing gap to give, fall back to
    opening a gap column for every sequence."""
    if not (0 <= row < aln.n_seqs):
        return None
    s = aln.rows[row]
    new = s[:col] + "-" + s[col:]
    if new.endswith("-"):
        rows = list(aln.rows)
        rows[row] = new[:-1]
        return _aln(aln, rows, keep_coding=True)
    return insert_gap_column(aln, col)


def _drop_trailing_gap_cols(rows: List[str], limit: int) -> List[str]:
    """Strip up to ``limit`` shared trailing columns that are all gaps."""
    for _ in range(limit):
        if rows and all(len(r) > 1 and r[-1] == "-" for r in rows):
            rows = [r[:-1] for r in rows]
        else:
            break
    return rows


def shift_rows(aln: Alignment, rows: Sequence[int], col: int, direction: int,
               unit: int = 1, keep_coding: bool = False) -> Optional[Alignment]:
    """Slide *only the selected sequences* by one step at nucleotide column ``col``.

    ``unit`` is the step size in nucleotide columns: 1 at the nucleotide level, 3
    at the codon / amino-acid level (so whole codons move and the frame is kept).

    direction +1 (right): open ``unit`` gaps at ``col`` in each selected row,
    pushing its residues right; other rows are padded with trailing gaps to stay
    rectangular, and wholly-gap trailing columns are then trimmed (so a sequence
    with a trailing gap absorbs it instead of widening the whole alignment).

    direction -1 (left): delete ``unit`` gap columns ending at ``col`` in each
    selected row, pulling its residues left. Only allowed if *every* selected row
    is all gaps across those columns; returns None otherwise.
    """
    sel = set(r for r in rows if 0 <= r < aln.n_seqs)
    if not sel or direction == 0 or unit < 1:
        return None
    pad = "-" * unit
    out = list(aln.rows)
    if direction > 0:
        col = max(0, min(col, aln.length))
        for r in range(aln.n_seqs):
            s = aln.rows[r]
            out[r] = (s[:col] + pad + s[col:]) if r in sel else (s + pad)
        out = _drop_trailing_gap_cols(out, unit)   # absorb trailing gaps, don't widen
    else:
        j = col - unit
        if j < 0:
            return None
        if any(aln.rows[r][j:col] != pad for r in sel):
            return None
        for r in range(aln.n_seqs):                # width stays fixed on a left slide
            s = aln.rows[r]
            out[r] = (s[:j] + s[col:] + pad) if r in sel else s
    return _aln(aln, out, keep_coding=keep_coding)


def close_gap(aln: Alignment, row: int, col: int) -> Optional[Alignment]:
    """Delete the gap at (row, col) in one sequence, pulling its residues left and
    appending a trailing gap to keep the alignment rectangular. Inverse of
    ``open_gap``; a no-op (None) if the cell isn't a gap."""
    if not (0 <= row < aln.n_seqs):
        return None
    s = aln.rows[row]
    if not (0 <= col < len(s)) or s[col] != "-":
        return None
    rows = list(aln.rows)
    rows[row] = s[:col] + s[col + 1:] + "-"
    return _aln(aln, rows, keep_coding=True)


def slide_block(aln: Alignment, row: int, col_lo: int, col_hi: int,
                direction: int) -> Optional[Alignment]:
    """Slide one row's residues in columns [col_lo, col_hi) by one column into the
    adjacent gap (direction -1/+1). Needs a gap in the slot the block moves into;
    that gap ends up filling the slot the block vacates. Returns None if blocked."""
    if not (0 <= row < aln.n_seqs):
        return None
    s = aln.rows[row]
    lo, hi = max(0, col_lo), min(len(s), col_hi)
    if lo >= hi:
        return None
    if direction > 0:
        if hi >= len(s) or s[hi] != "-":     # need a gap just right of the block
            return None
        region = s[lo:hi + 1]
        region = region[-1] + region[:-1]    # gap rotates to the front
        new = s[:lo] + region + s[hi + 1:]
    elif direction < 0:
        if lo == 0 or s[lo - 1] != "-":      # need a gap just left of the block
            return None
        region = s[lo - 1:hi]
        region = region[1:] + region[0]      # gap rotates to the back
        new = s[:lo - 1] + region + s[hi:]
    else:
        return None
    rows = list(aln.rows)
    rows[row] = new
    return _aln(aln, rows, keep_coding=True)


def column_support(aln: Alignment, col: int) -> float:
    """Sum-of-pairs identity of one column (1.0 = all identical / <2 residues)."""
    if not (0 <= col < aln.length):
        return 0.0
    present = [r[col].upper() for r in aln.rows if r[col] != "-"]
    if len(present) < 2:
        return 1.0
    tot = match = 0
    for a in range(len(present)):
        for b in range(a + 1, len(present)):
            tot += 1
            if present[a] == present[b]:
                match += 1
    return match / tot if tot else 1.0


def anchored_realign(aln: Alignment, anchors: Sequence[int], engine, **opts) -> Alignment:
    """Realign the spans *between* anchor columns, keeping the anchors fixed.

    Constraint-based curation: pin the columns you trust (anchors) and let the
    engine re-solve only the uncertain stretches in between.
    """
    ids = list(aln.ids)
    anchor_set = sorted({a for a in anchors if 0 <= a < aln.length})
    pieces = []           # ('anchor', col) or ('block', start, stop)
    prev = 0
    for a in anchor_set:
        if a > prev:
            pieces.append(("block", prev, a))
        pieces.append(("anchor", a))
        prev = a + 1
    if prev < aln.length:
        pieces.append(("block", prev, aln.length))

    out = {i: "" for i in ids}
    # Realigning the blocks between anchors changes their widths, so the anchors
    # land at new column indices. Record them, so a caller holding these columns
    # as pins can follow them instead of keeping indices that now point at
    # whatever the realignment moved into that position.
    moved: List[int] = []
    width = 0
    for piece in pieces:
        if piece[0] == "anchor":
            c = piece[1]
            moved.append(width)
            width += 1
            for i, r in zip(ids, aln.rows):
                out[i] += r[c]
        else:
            _, s, e = piece
            recs = [(i, r[s:e].replace("-", "")) for i, r in zip(ids, aln.rows)]
            nonempty = [(i, sub) for i, sub in recs if sub]
            if len(nonempty) >= 2:
                block = engine.align(nonempty, aln.alphabet, **opts)
                bmap = dict(zip(block.ids, block.rows))
                bw = block.length
            elif len(nonempty) == 1:
                bmap = {nonempty[0][0]: nonempty[0][1]}
                bw = len(nonempty[0][1])
            else:
                bmap, bw = {}, 0
            for i in ids:
                out[i] += bmap.get(i, "-" * bw)
            width += bw
    new = _aln(aln, [out[i] for i in ids])
    if moved:
        new.meta["anchors"] = moved
    return new
