"""Core domain model.

A single canonical representation is stored at the sequences' native level
(nucleotide rows for DNA/RNA, amino-acid rows for protein). Codon and
amino-acid *views* are derived projections of a coding nucleotide alignment, so
the GUI can switch levels live without ever holding three copies of the data.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Sequence, Tuple


class Alphabet(Enum):
    DNA = "dna"
    RNA = "rna"
    PROTEIN = "protein"

    @property
    def is_nucleotide(self) -> bool:
        return self in (Alphabet.DNA, Alphabet.RNA)


class Level(Enum):
    NT = "nt"
    CODON = "codon"
    AA = "aa"


_NUC = set("ACGTUN")

#: Characters that mean "no residue here". ``-`` is the canonical form; ``.``
#: (Stockholm/Pfam insert columns) and ``?`` (NEXUS missing data) are accepted
#: on input and normalised to ``-`` by :class:`Alignment`, so that the rest of
#: the codebase can test ``ch != "-"`` and be right. Normalising in one place
#: rather than teaching every residue-index walk about three characters is
#: deliberate: a single missed site silently shifts every residue index and
#: corrupts every reliability score computed from it.
GAP_CHARS = frozenset("-.?")
_GAP = GAP_CHARS  # backwards-compatible alias

_GAP_NORMALISE = str.maketrans({c: "-" for c in GAP_CHARS if c != "-"})


def normalise_gaps(row: str) -> str:
    """Map every accepted gap character onto the canonical ``-``."""
    return row.translate(_GAP_NORMALISE)


def ungap(row: str) -> str:
    """The residues of a row, with every gap character removed."""
    return normalise_gaps(row).replace("-", "")


# --------------------------------------------------------------------------- #
# Residue indexing
# --------------------------------------------------------------------------- #
#
# Every part of CRAIC that reasons about homology has to answer two questions:
# which residue sits in this column, and which column holds this residue. Those
# walks used to be reimplemented privately in six modules. They are trivial, and
# that was the problem: trivial code gets copied, copies drift, and a drifted
# copy shifts every residue index downstream and silently corrupts every score
# computed from it — which is exactly the shape of the v0.2 gap-character bug.
# One implementation, here, next to the definition of what a gap is.

def residue_index(row: str) -> List[int]:
    """Per column, the index of its residue in the ungapped sequence, or -1."""
    out: List[int] = []
    r = 0
    for ch in row:
        if ch == "-":
            out.append(-1)
        else:
            out.append(r)
            r += 1
    return out


def residue_index_array(row: str):
    """:func:`residue_index` as a NumPy array, for vectorised callers."""
    import numpy as np

    return np.asarray(residue_index(row), dtype=int)


def column_index(row: str) -> Dict[int, int]:
    """The inverse of :func:`residue_index`: residue index -> column."""
    return {r: c for c, r in enumerate(residue_index(row)) if r >= 0}


def membership(aln: "Alignment", by: str = "index") -> Dict[Tuple[Any, int], int]:
    """``(sequence, residue) -> column`` for every non-gap cell.

    ``by="index"`` keys on the sequence's position, ``by="id"`` on its
    identifier. Both are needed: scores computed per row want the position,
    while anything comparing two alignments has to match on the name.
    """
    if by not in ("index", "id"):
        raise ValueError(f"by must be 'index' or 'id', not {by!r}")
    keys = range(aln.n_seqs) if by == "index" else aln.ids
    mem: Dict[Tuple[Any, int], int] = {}
    for key, row in zip(keys, aln.rows):
        for c, r in enumerate(residue_index(row)):
            if r >= 0:
                mem[(key, r)] = c
    return mem


# IUPAC nucleotide ambiguity: the code for a *set* of bases.
_IUPAC = {
    frozenset("A"): "A", frozenset("C"): "C", frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M",
    frozenset("CGT"): "B", frozenset("AGT"): "D", frozenset("ACT"): "H", frozenset("ACG"): "V",
    frozenset("ACGT"): "N",
}


def _iupac_nt(bases: Sequence[str]) -> str:
    """IUPAC code for the most-common base(s); ties give an ambiguity code."""
    from collections import Counter
    norm = [b.upper().replace("U", "T") for b in bases if b not in _GAP]
    if not norm:
        return "-"
    cnt = Counter(norm)
    top = max(cnt.values())
    winners = frozenset(b for b, n in cnt.items() if n == top)
    return _IUPAC.get(winners, "N")


def _majority_or_x(vals: Sequence[str], unknown: str = "X") -> str:
    """The unique most-common symbol, or ``unknown`` when there's a tie."""
    from collections import Counter
    cnt = Counter(v.upper() for v in vals)
    top = max(cnt.values())
    winners = [v for v, n in cnt.items() if n == top]
    return winners[0] if len(winners) == 1 else unknown


def column_consensus(cells: Sequence[str], level: Level, alphabet: Alphabet) -> str:
    """Consensus for one display column's cells, respecting the viewing level:

    * nucleotide → IUPAC code (a tie of A/T gives W, etc.)
    * codon      → per-position IUPAC across the codons (e.g. ``ATW``)
    * amino acid → the majority residue, or ``X`` when there is no majority
    """
    vals = [c for c in cells if c.strip("-")]
    if not vals:
        return ""
    if level == Level.CODON:
        width = max(len(v) for v in vals)
        out = []
        for p in range(width):
            bases = [v[p] for v in vals if p < len(v) and v[p] not in _GAP]
            out.append(_iupac_nt(bases) if bases else "-")
        return "".join(out)
    if alphabet == Alphabet.PROTEIN or level == Level.AA:
        return _majority_or_x(vals)
    return _iupac_nt(vals)


def detect_alphabet(seqs: Sequence[str]) -> Alphabet:
    """Heuristic alphabet detection from residue composition."""
    counts: Dict[str, int] = {}
    total = 0
    for s in seqs:
        for c in s.upper():
            if c in _GAP:
                continue
            counts[c] = counts.get(c, 0) + 1
            total += 1
    if total == 0:
        return Alphabet.PROTEIN
    nuc = sum(v for c, v in counts.items() if c in _NUC)
    if nuc / total >= 0.9:
        has_u = counts.get("U", 0) > 0
        has_t = counts.get("T", 0) > 0
        return Alphabet.RNA if (has_u and not has_t) else Alphabet.DNA
    return Alphabet.PROTEIN


@dataclass(frozen=True)
class CodingSpec:
    """How a nucleotide alignment maps to codons."""

    frame: int = 0   # leading alignment columns to skip before the first codon
    table: int = 1   # NCBI genetic-code table id


def _build_table(table_id: int) -> Dict[str, str]:
    from Bio.Data import CodonTable

    t = CodonTable.unambiguous_dna_by_id[table_id]
    d = dict(t.forward_table)
    for stop in t.stop_codons:
        d[stop] = "*"
    return d


_TABLE_CACHE: Dict[int, Dict[str, str]] = {}


def genetic_codes() -> List[Tuple[int, str]]:
    """Every NCBI translation table Biopython knows, as ``(id, name)``.

    The list is not hard-coded here because it is not ours to curate: NCBI adds
    tables, and anything shipped as a fixed list goes stale silently. Table 1 is
    forced to the front so the standard code is the first thing offered.
    """
    from Bio.Data import CodonTable

    codes = []
    for tid, tbl in sorted(CodonTable.unambiguous_dna_by_id.items()):
        name = next((n for n in tbl.names if n and not n.startswith("SGC")), f"table {tid}")
        codes.append((tid, name))
    return sorted(codes, key=lambda c: (c[0] != 1, c[0]))


def code_name(table_id: int) -> str:
    """The NCBI name for a table id, or a bare label if it is unknown."""
    for tid, name in genetic_codes():
        if tid == table_id:
            return name
    return f"table {table_id}"


def translate_codon(codon: str, table_id: int = 1) -> str:
    """Translate one (possibly gapped) codon to a single residue.

    '---' -> '-' (clean gap), any partial gap -> 'X' (frameshift / indel
    ambiguity, surfaced rather than hidden), unknown -> 'X', stop -> '*'.
    """
    codon = codon.upper().replace("U", "T")
    if codon == "---":
        return "-"
    if "-" in codon:
        return "X"
    if table_id not in _TABLE_CACHE:
        _TABLE_CACHE[table_id] = _build_table(table_id)
    return _TABLE_CACHE[table_id].get(codon, "X")


@dataclass
class Alignment:
    ids: List[str]
    rows: List[str]
    alphabet: Alphabet
    coding: Optional[CodingSpec] = None
    meta: Dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.rows:
            return
        # Canonicalise gap characters once, here, so that every residue-index
        # walk downstream can test ``ch != "-"``. See GAP_CHARS.
        self.rows = [normalise_gaps(r) for r in self.rows]
        width = len(self.rows[0])
        if any(len(r) != width for r in self.rows):
            raise ValueError("alignment rows are not equal length")
        if len(self.ids) != len(self.rows):
            raise ValueError("ids and rows count mismatch")

    # -- basic shape ------------------------------------------------------- #
    @property
    def n_seqs(self) -> int:
        return len(self.rows)

    @property
    def length(self) -> int:
        return len(self.rows[0]) if self.rows else 0

    def column(self, j: int) -> str:
        return "".join(r[j] for r in self.rows)

    def gap_fraction(self, j: int) -> float:
        col = self.column(j)
        return sum(c in _GAP for c in col) / len(col)

    # -- levels ------------------------------------------------------------ #
    def available_levels(self) -> List[Level]:
        if self.alphabet == Alphabet.PROTEIN:
            return [Level.AA]
        if self.coding is not None:
            return [Level.NT, Level.CODON, Level.AA]
        return [Level.NT]

    def n_codons(self) -> int:
        if not self.coding:
            return 0
        return max(0, (self.length - self.coding.frame) // 3)

    def _codon_cols(self, c: int) -> Tuple[int, int, int]:
        start = self.coding.frame + 3 * c
        return start, start + 1, start + 2

    def is_codon_aware(self) -> bool:
        """True if no codon contains a partial gap (frame is intact)."""
        if not self.coding:
            return False
        for r in self.rows:
            for c in range(self.n_codons()):
                a, b, d = self._codon_cols(c)
                cod = r[a] + r[b] + r[d]
                gaps = cod.count("-")
                if gaps not in (0, 3):
                    return False
        return True

    def frame_break_columns(self) -> List[int]:
        """Codon columns where at least one row has a partial-gap codon."""
        if not self.coding:
            return []
        out = []
        for c in range(self.n_codons()):
            a, b, d = self._codon_cols(c)
            broken = any(
                0 < (r[a] + r[b] + r[d]).count("-") < 3 for r in self.rows
            )
            if broken:
                out.append(c)
        return out

    def internal_stops(self) -> List[Tuple[int, int]]:
        """``(row, codon)`` for every stop codon that is not the last codon of
        its sequence.

        A stop in the middle of a coding sequence usually means one of three
        things, and the user needs to be told which is being assumed: the wrong
        genetic code (mitochondrial and Mollicute genes read TGA as tryptophan,
        the standard code reads it as a stop), the wrong reading frame, or a
        pseudogene. One check catches all three, which is why it is worth making
        even though CRAIC cannot tell them apart.

        Trailing gaps are ignored, so a short sequence padded out to the
        alignment width is not reported as stopping early.
        """
        if not self.coding:
            return []
        out: List[Tuple[int, int]] = []
        for ri, row in enumerate(self.rows):
            last = -1
            for c in range(self.n_codons()):
                a, b, d = self._codon_cols(c)
                if (row[a] + row[b] + row[d]) != "---":
                    last = c
            for c in range(self.n_codons()):
                if c >= last:
                    break                      # the final codon may legitimately stop
                a, b, d = self._codon_cols(c)
                if translate_codon(row[a] + row[b] + row[d], self.coding.table) == "*":
                    out.append((ri, c))
        return out

    def display(self, level: Level) -> "DisplayMatrix":
        """Project to a renderable matrix of string cells.

        Returns ``cells[row][col]`` plus, for each display column, the starting
        nucleotide index so the GUI can map selections back to nt coordinates.
        """
        if level == Level.NT or self.alphabet == Alphabet.PROTEIN:
            cells = [list(r) for r in self.rows]
            nt_start = list(range(self.length))
            return DisplayMatrix(cells, nt_start, span=1)

        if not self.coding:
            raise ValueError("codon/aa view requires a CodingSpec")

        nc = self.n_codons()
        nt_start = [self.coding.frame + 3 * c for c in range(nc)]
        cells: List[List[str]] = []
        for r in self.rows:
            row_cells = []
            for c in range(nc):
                a, b, d = self._codon_cols(c)
                cod = r[a] + r[b] + r[d]
                if level == Level.CODON:
                    row_cells.append(cod)
                else:  # AA
                    row_cells.append(translate_codon(cod, self.coding.table))
            cells.append(row_cells)
        return DisplayMatrix(cells, nt_start, span=3)

    def translated_rows(self) -> List[str]:
        dm = self.display(Level.AA)
        return ["".join(c) for c in dm.cells]

    def amino_acid_alignment(self) -> "Alignment":
        """A protein Alignment of the translated rows, for amino-acid-level tools.

        For protein data this is the alignment itself; for a coding nucleotide
        alignment each amino-acid column corresponds to one codon column.
        """
        if self.alphabet == Alphabet.PROTEIN:
            return self
        if not self.coding:
            raise ValueError("amino-acid view requires a CodingSpec")
        return Alignment(list(self.ids), self.translated_rows(), Alphabet.PROTEIN)

    def codon_col(self, nt_col: int) -> int:
        """Map a nucleotide column to its codon / amino-acid column index."""
        if not self.coding:
            return nt_col
        return (nt_col - self.coding.frame) // 3

    # -- io helpers -------------------------------------------------------- #
    def slice_columns(self, start: int, stop: int) -> "Alignment":
        return Alignment(
            list(self.ids),
            [r[start:stop] for r in self.rows],
            self.alphabet,
            coding=None,
            meta=dict(self.meta),
        )

    @classmethod
    def from_records(
        cls,
        records: Sequence[Tuple[str, str]],
        alphabet: Optional[Alphabet] = None,
        coding: Optional[CodingSpec] = None,
    ) -> "Alignment":
        ids = [r[0] for r in records]
        rows = [r[1] for r in records]
        alph = alphabet or detect_alphabet(rows)
        return cls(ids, rows, alph, coding=coding)


@dataclass
class DisplayMatrix:
    cells: List[List[str]]   # cells[row][col]
    nt_start: List[int]      # nt index where each display column begins
    span: int                # nt columns per display column (1 or 3)

    @property
    def n_cols(self) -> int:
        return len(self.nt_start)


# --------------------------------------------------------------------------- #
# Matching two alignments of the same sequences
# --------------------------------------------------------------------------- #

class IdMismatch(ValueError):
    """Two alignments of supposedly the same sequences do not share sequence ids."""


def _short(rid: str) -> str:
    """The identifier as a whitespace-truncating tool would have written it."""
    return rid.split()[0] if rid.split() else rid


def match_ids(reference: Alignment, other: Alignment) -> Dict[str, str]:
    """Map each reference id onto the corresponding id in ``other``.

    Several aligners (MAFFT among them) truncate FASTA identifiers at the first
    whitespace. Matching on the raw string then misses every sequence, and
    because a miss increments the pair denominator without incrementing the
    numerator, *every column silently scores zero agreement* — a wrong result
    that looks exactly like a real one. So: try the ids as given, then try them
    truncated, and if neither matches, say so.
    """
    other_ids = set(other.ids)
    if set(reference.ids) <= other_ids:
        return {rid: rid for rid in reference.ids}

    by_short: Dict[str, List[str]] = {}
    for oid in other.ids:
        by_short.setdefault(_short(oid), []).append(oid)
    mapped = {}
    for rid in reference.ids:
        cands = by_short.get(_short(rid), [])
        if len(cands) == 1:
            mapped[rid] = cands[0]
    if len(mapped) == len(reference.ids):
        return mapped

    missing = [rid for rid in reference.ids if rid not in mapped]
    raise IdMismatch(
        f"{len(missing)} of {len(reference.ids)} sequence ids in the reference "
        f"have no unique counterpart in the other alignment "
        f"(first unmatched: {missing[0]!r}); the two alignments cannot be compared"
    )
