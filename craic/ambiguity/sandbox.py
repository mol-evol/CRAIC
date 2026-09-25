"""Local realignment sandbox: the human-in-the-loop core.

Select a hard block of columns, re-align just that block under alternative
engines / parameters, compare the alternatives, and splice the chosen one back
in without disturbing the rest of the alignment.

The block is re-aligned at the level you are *viewing*: when you look at a coding
alignment as codons or amino acids, the alternatives are aligned codon-aware and
shown at that level, so you never compare a nucleotide re-alignment while reading
amino acids.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

from ..domain import Alignment, Level
from ..engines import AlignerEngine, BuiltinProgressive, CodonAware

Record = Tuple[str, str]


def extract_block(aln: Alignment, col_start: int, col_stop: int) -> List[Record]:
    """Ungapped subsequences within [col_start, col_stop). Empty rows dropped."""
    out: List[Record] = []
    for name, row in zip(aln.ids, aln.rows):
        sub = row[col_start:col_stop].replace("-", "")
        if sub:
            out.append((name, sub))
    return out


def _sp_identity(rows: Sequence[str]) -> float:
    """Cheap quality proxy: mean identity over aligned non-gap pairs."""
    if len(rows) < 2:
        return 0.0
    width = len(rows[0])
    total = match = 0
    for c in range(width):
        present = [r[c] for r in rows if r[c] != "-"]
        for a in range(len(present)):
            for b in range(a + 1, len(present)):
                total += 1
                if present[a].upper() == present[b].upper():
                    match += 1
    return match / total if total else 0.0


@dataclass
class Alternative:
    label: str
    block: Alignment                 # the re-aligned block (nucleotide if coding)
    display_level: Level = Level.NT  # the level it should be shown at

    def display_rows(self) -> List[str]:
        if self.display_level in (Level.CODON, Level.AA) and self.block.coding:
            dm = self.block.display(self.display_level)
            return ["".join(c) for c in dm.cells]
        return list(self.block.rows)

    @property
    def width(self) -> int:
        rows = self.display_rows()
        return len(rows[0]) if rows else 0

    @property
    def score(self) -> float:
        return _sp_identity(self.display_rows())


def realign_block(
    aln: Alignment,
    col_start: int,
    col_stop: int,
    engine: AlignerEngine,
    **opts,
) -> Optional[Alignment]:
    records = extract_block(aln, col_start, col_stop)
    if len(records) < 1:
        return None
    return engine.align(records, aln.alphabet, **opts)


def alternatives(
    aln: Alignment,
    col_start: int,
    col_stop: int,
    engines: Sequence[AlignerEngine],
    codon_aware: bool = False,
    table: int = 1,
    level: Level = Level.NT,
    param_grid: Optional[Sequence[float]] = None,
    progress: Optional[Callable[[int, int, str], None]] = None,
    cancelled: Optional[Callable[[], bool]] = None,
    params_by_key: Optional[dict] = None,
) -> List[Alternative]:
    """Labelled block alignments from several engines / gap settings.

    With ``codon_aware`` the block is re-aligned in amino-acid space (translate →
    align → back-translate), so the result keeps the reading frame and can be
    shown as codons or amino acids.

    ``params_by_key`` maps an engine's key to the parameters it is run with (the
    settings dialog's values); an engine not in it runs with its defaults. The
    built-in gap-probability grid is a fixed contrast and is not affected.

    ``progress(done, total, label)`` is called before each engine runs, and
    ``cancelled()`` is polled between engines so a long run can be stopped early
    (the in-flight engine finishes first — cancellation is cooperative).
    """
    records = extract_block(aln, col_start, col_stop)
    out: List[Alternative] = []
    if len(records) < 2:
        return out

    def wrap(e: AlignerEngine) -> AlignerEngine:
        return CodonAware(e, table) if codon_aware else e

    avail = [e for e in engines if e.available()]
    grid = list(param_grid or [0.01, 0.05])   # pair-HMM gap-open probabilities
    total = len(avail) + len(grid)
    done = 0

    for eng in avail:
        if cancelled is not None and cancelled():
            return out
        if progress is not None:
            progress(done, total, f"Aligning with {eng.label}…")
        try:
            params = (params_by_key or {}).get(eng.key, {})
            block = wrap(eng).align(records, aln.alphabet, **params)
            out.append(Alternative(eng.label, block, level))
        except Exception:
            pass
        done += 1

    builtin = BuiltinProgressive()
    for d in grid:
        if cancelled is not None and cancelled():
            return out
        if progress is not None:
            progress(done, total, f"built-in (gap prob {d})…")
        try:
            block = wrap(builtin).align(records, aln.alphabet, delta=d, epsilon=0.5,
                                        estimate=False, consistency_iters=1)
            out.append(Alternative(f"built-in (gap prob {d})", block, level))
        except Exception:
            pass
        done += 1

    if progress is not None:
        progress(total, total, "Done")
    return out


def splice(aln: Alignment, col_start: int, col_stop: int, block: Alignment) -> Alignment:
    """Replace columns [col_start, col_stop) with `block`, re-padding all rows.

    Sequences absent from the block receive gaps for the new block width, so the
    result stays rectangular.

    The coding annotation is preserved only if both the alignment and the block
    are codon-aware *and* the splice does not change the width modulo 3. A block
    that is wider or narrower by a non-multiple of three re-phases every codon
    downstream while ``CodingSpec.frame`` stays where it was, which silently
    produces a wrong amino-acid view of correct nucleotides — so in that case the
    annotation is dropped and the user re-sets the frame deliberately.
    """
    width = block.length
    in_frame = (width - (col_stop - col_start)) % 3 == 0
    block_by_id = dict(zip(block.ids, block.rows))
    if len(block_by_id) != len(block.ids):
        raise ValueError("the replacement block has duplicate sequence ids")
    new_rows: List[str] = []
    for name, row in zip(aln.ids, aln.rows):
        mid = block_by_id.get(name, "-" * width)
        if len(mid) != width:  # safety pad
            mid = mid.ljust(width, "-")
        new_rows.append(row[:col_start] + mid + row[col_stop:])
    keep_coding = aln.coding if (aln.coding and block.coding and in_frame) else None
    return Alignment(list(aln.ids), new_rows, aln.alphabet,
                     coding=keep_coding, meta=dict(aln.meta))
