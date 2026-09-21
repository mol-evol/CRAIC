"""Reading and writing alignments / unaligned record sets."""

from __future__ import annotations

import os
from typing import List, Optional, Sequence, Tuple

from .domain import Alignment, Alphabet, CodingSpec, detect_alphabet

Record = Tuple[str, str]

_EXT_FMT = {
    ".fa": "fasta", ".fasta": "fasta", ".fna": "fasta", ".faa": "fasta",
    ".phy": "phylip-relaxed", ".phylip": "phylip-relaxed",
    ".aln": "clustal", ".sto": "stockholm", ".stk": "stockholm",
    ".nex": "nexus", ".nexus": "nexus",
}


def guess_format(path: str) -> str:
    return _EXT_FMT.get(os.path.splitext(path)[1].lower(), "fasta")


def read_records(path: str, fmt: Optional[str] = None) -> List[Record]:
    from Bio import SeqIO

    fmt = fmt or guess_format(path)
    return [(r.id, str(r.seq)) for r in SeqIO.parse(path, fmt)]


def looks_aligned(records: Sequence[Record]) -> bool:
    if not records:
        return False
    lengths = {len(s) for _, s in records}
    has_gap = any("-" in s for _, s in records)
    return len(lengths) == 1 and (has_gap or len(records) > 1)


# Named sequence groups and column annotations are stored in a small companion
# file next to the alignment, so the alignment file itself stays a standard format.
def sidecar_path(path: str) -> str:
    return path + ".craic.json"


def _read_sidecar(path: str) -> dict:
    import json

    side = sidecar_path(path)
    if not os.path.exists(side):
        return {}
    try:
        with open(side) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def load_groups(path: str) -> dict:
    """Read named sequence groups: ``{name: [sequence_id, ...]}`` (or empty)."""
    groups = _read_sidecar(path).get("groups", {})
    return {str(k): [str(i) for i in v] for k, v in groups.items() if v}


def load_annotations(path: str) -> list:
    """Read column annotations: ``[{name, start, stop, color}, ...]`` (or empty)."""
    anns = _read_sidecar(path).get("annotations", [])
    out = []
    for a in anns:
        try:
            out.append({"name": str(a.get("name", "")),
                        "start": int(a["start"]), "stop": int(a["stop"]),
                        "color": str(a.get("color", "#5b8def"))})
        except (KeyError, ValueError, TypeError):
            continue
    return out


def write_sidecar(path: str, groups: dict, annotations: list) -> None:
    """Write groups + annotations to the sidecar (or remove it if both are empty)."""
    import json

    side = sidecar_path(path)
    data = {}
    if groups:
        data["groups"] = {k: list(v) for k, v in groups.items()}
    if annotations:
        data["annotations"] = list(annotations)
    if not data:
        if os.path.exists(side):
            try:
                os.remove(side)
            except OSError:
                pass
        return
    with open(side, "w") as fh:
        json.dump(data, fh, indent=2)


def load_alignment(
    path: str,
    fmt: Optional[str] = None,
    alphabet: Optional[Alphabet] = None,
    coding: Optional[CodingSpec] = None,
) -> Alignment:
    recs = read_records(path, fmt)
    if not recs:
        return Alignment([], [], alphabet or Alphabet.PROTEIN, coding=coding)
    alph = alphabet or detect_alphabet([s for _, s in recs])
    groups = load_groups(path)
    annotations = load_annotations(path)
    if looks_aligned(recs):
        aln = Alignment.from_records(recs, alphabet=alph, coding=coding)
        if groups:
            aln.meta["groups"] = groups
        if annotations:
            aln.meta["annotations"] = annotations
        return aln
    # Not aligned: keep as a degenerate single-column-less alignment is wrong;
    # callers that load unaligned data should run an engine. We still return a
    # padded alignment so the viewer can show the raw sequences.
    width = max(len(s) for _, s in recs)
    padded = [(i, s.ljust(width, "-")) for i, s in recs]
    aln = Alignment.from_records(padded, alphabet=alph, coding=coding)
    aln.meta["unaligned"] = True
    if groups:
        aln.meta["groups"] = groups
    if annotations:
        aln.meta["annotations"] = annotations
    return aln


def write_fasta(path: str, ids: Sequence[str], rows: Sequence[str], wrap: int = 60) -> None:
    with open(path, "w") as fh:
        for name, seq in zip(ids, rows):
            fh.write(f">{name}\n")
            if wrap and wrap > 0:
                for k in range(0, len(seq), wrap):
                    fh.write(seq[k : k + wrap] + "\n")
            else:
                fh.write(seq + "\n")


# Write formats offered in the Save dialog: (key, human label, default extension).
WRITE_FORMATS: List[Tuple[str, str, str]] = [
    ("fasta", "FASTA", ".fasta"),
    ("clustal", "Clustal", ".aln"),
    ("phylip", "PHYLIP (interleaved)", ".phy"),
    ("phylip-relaxed", "PHYLIP (relaxed)", ".phy"),
    ("phylip-sequential", "PHYLIP (sequential)", ".phy"),
    ("stockholm", "Stockholm", ".sto"),
    ("nexus", "NEXUS", ".nex"),
    ("mega", "MEGA", ".meg"),
]

_WRITE_EXT = {
    ".fa": "fasta", ".fasta": "fasta", ".aln": "clustal",
    ".phy": "phylip-relaxed", ".phylip": "phylip-relaxed",
    ".sto": "stockholm", ".stk": "stockholm",
    ".nex": "nexus", ".nexus": "nexus", ".meg": "mega", ".mega": "mega",
}


def write_formats() -> List[Tuple[str, str, str]]:
    return list(WRITE_FORMATS)


def guess_write_format(path: str) -> str:
    return _WRITE_EXT.get(os.path.splitext(path)[1].lower(), "fasta")


def _molecule_type(aln: Alignment) -> str:
    if aln.alphabet == Alphabet.PROTEIN:
        return "protein"
    if aln.alphabet == Alphabet.RNA:
        return "RNA"
    return "DNA"


def _to_msa(aln: Alignment):
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    mol = _molecule_type(aln)
    records = []
    for name, row in zip(aln.ids, aln.rows):
        rec = SeqRecord(Seq(row), id=str(name), name=str(name), description="")
        rec.annotations["molecule_type"] = mol
        records.append(rec)
    msa = MultipleSeqAlignment(records)
    msa.annotations["molecule_type"] = mol  # NEXUS writer reads this
    return msa


def dump_alignment(aln: Alignment, fmt: str = "fasta", wrap: int = 0) -> str:
    """Serialise `aln` to a string. `fmt` is a key from WRITE_FORMATS.

    Single source of truth for serialisation — both Save (to a file) and Copy (to
    the clipboard) go through here.
    """
    fmt = (fmt or "fasta").lower()
    if fmt == "fasta":
        out: List[str] = []
        for name, seq in zip(aln.ids, aln.rows):
            out.append(f">{name}")
            if wrap and wrap > 0:
                out += [seq[k : k + wrap] for k in range(0, len(seq), wrap)]
            else:
                out.append(seq)
        return "\n".join(out) + "\n"
    if fmt == "mega":
        dt = "Protein" if aln.alphabet == Alphabet.PROTEIN else "DNA"
        out = ["#MEGA", "!Title CRAIC alignment;", f"!Format DataType={dt} indel=-;", ""]
        for name, row in zip(aln.ids, aln.rows):
            out.append(f"#{str(name).replace(' ', '_')}")
            out.append(row)
        return "\n".join(out) + "\n"
    from io import StringIO

    from Bio import AlignIO

    buf = StringIO()
    AlignIO.write(_to_msa(aln), buf, fmt)
    return buf.getvalue()


def save_alignment(path: str, aln: Alignment, fmt: Optional[str] = None) -> None:
    """Write `aln` to `path`. `fmt` is a key from WRITE_FORMATS; if omitted it is
    inferred from the file extension."""
    fmt = (fmt or guess_write_format(path)).lower()
    wrap = 60 if fmt == "fasta" else 0
    with open(path, "w") as fh:
        fh.write(dump_alignment(aln, fmt, wrap=wrap))
    # Named groups + column annotations travel in a companion sidecar, leaving the
    # alignment file itself a clean, standard format.
    write_sidecar(path, aln.meta.get("groups", {}), aln.meta.get("annotations", []))
