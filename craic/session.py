"""The working session as a document.

A curation session is more than its alignment. It is also the reference
alignment truth mode is scoring against, the column annotations, the named
sequence groups, the provenance log of how the alignment was arrived at, and the
view state that makes any of it legible. Exporting FASTA or NEXUS keeps the
residues and silently discards all of that — fine as an *export*, useless as a
way to put work down and pick it up again.

So a session is its own small document: readable JSON, no dependencies, and
versioned so that an old file can be recognised rather than misread. The
alignment can still be saved to any standard format at any time; the two are
separate acts, because they answer different questions ("give this to another
program" versus "let me carry on tomorrow").

The module is deliberately free of Qt, so that the command line, the tests and
the crash-recovery autosave can all use it without a display.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from . import __version__
from .domain import Alignment, Alphabet, CodingSpec

#: Bumped when the on-disk shape changes incompatibly. A reader that does not
#: recognise the major version refuses the file rather than guessing at it.
FORMAT_VERSION = 1

#: Suffix for an explicit session document saved beside the alignment.
SUFFIX = ".craic.json"


class SessionError(ValueError):
    """A session file could not be read as one."""


def session_path_for(alignment_path: str) -> str:
    """The session document that belongs beside ``alignment_path``."""
    return alignment_path + SUFFIX


def _alignment_to_dict(aln: Alignment) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "ids": list(aln.ids),
        "rows": list(aln.rows),
        "alphabet": aln.alphabet.value,
    }
    if aln.coding is not None:
        out["coding"] = {"frame": aln.coding.frame, "table": aln.coding.table}
    # Only JSON-safe metadata travels; anything else is derived and can be
    # recomputed, so storing it would be a way to go stale rather than a saving.
    meta = {}
    for key, value in (aln.meta or {}).items():
        try:
            json.dumps(value)
        except (TypeError, ValueError):
            continue
        meta[key] = value
    if meta:
        out["meta"] = meta
    return out


def _alignment_from_dict(d: Optional[Dict[str, Any]]) -> Optional[Alignment]:
    if not d:
        return None
    try:
        alphabet = Alphabet(d["alphabet"])
        coding = CodingSpec(**d["coding"]) if d.get("coding") else None
        return Alignment(list(d["ids"]), list(d["rows"]), alphabet,
                         coding=coding, meta=dict(d.get("meta", {})))
    except (KeyError, TypeError, ValueError) as exc:
        raise SessionError(f"the alignment in this session is unreadable: {exc}") from exc


@dataclass
class Session:
    """Everything needed to resume a curation session."""

    alignment: Optional[Alignment] = None
    #: the trusted alignment truth mode scores against, and what to call it
    reference: Optional[Alignment] = None
    reference_label: str = ""
    annotations: List[dict] = field(default_factory=list)
    groups: Dict[str, List[str]] = field(default_factory=dict)
    #: where the alignment itself was last read from or written to
    source_path: str = ""
    #: view state — which overlay, the mask threshold, the viewing level
    view: Dict[str, Any] = field(default_factory=dict)
    #: set by the autosave; an explicit save clears it
    autosaved: bool = False
    saved_at: float = 0.0
    craic_version: str = ""

    # -- serialisation ----------------------------------------------------- #
    def to_dict(self) -> Dict[str, Any]:
        return {
            "format": "craic-session",
            "format_version": FORMAT_VERSION,
            "craic_version": self.craic_version or __version__,
            "saved_at": self.saved_at or time.time(),
            "autosaved": bool(self.autosaved),
            "source_path": self.source_path,
            "alignment": _alignment_to_dict(self.alignment) if self.alignment else None,
            "reference": _alignment_to_dict(self.reference) if self.reference else None,
            "reference_label": self.reference_label,
            "annotations": list(self.annotations),
            "groups": {k: list(v) for k, v in self.groups.items()},
            "view": dict(self.view),
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Session":
        if not isinstance(d, dict) or d.get("format") != "craic-session":
            raise SessionError("this file is not a CRAIC session")
        version = d.get("format_version")
        if not isinstance(version, int) or version > FORMAT_VERSION:
            raise SessionError(
                f"this session was written by a newer CRAIC (format {version}); "
                f"this one understands up to {FORMAT_VERSION}")
        return cls(
            alignment=_alignment_from_dict(d.get("alignment")),
            reference=_alignment_from_dict(d.get("reference")),
            reference_label=d.get("reference_label", "") or "",
            annotations=list(d.get("annotations") or []),
            groups={k: list(v) for k, v in (d.get("groups") or {}).items()},
            source_path=d.get("source_path", "") or "",
            view=dict(d.get("view") or {}),
            autosaved=bool(d.get("autosaved")),
            saved_at=float(d.get("saved_at") or 0.0),
            craic_version=d.get("craic_version", "") or "",
        )

    # -- files -------------------------------------------------------------- #
    def save(self, path: str) -> None:
        """Write the session, atomically.

        Via a temporary file and a rename, because the autosave writes this
        repeatedly and a crash during the write is exactly the moment the file
        matters most — a half-written recovery file is worse than none.
        """
        self.saved_at = time.time()
        self.craic_version = __version__
        tmp = path + ".part"
        directory = os.path.dirname(os.path.abspath(path))
        os.makedirs(directory, exist_ok=True)
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=1)
        os.replace(tmp, path)

    @classmethod
    def load(cls, path: str) -> "Session":
        try:
            with open(path, encoding="utf-8") as fh:
                data = json.load(fh)
        except OSError as exc:
            raise SessionError(f"cannot read {path}: {exc}") from exc
        except ValueError as exc:
            raise SessionError(f"{os.path.basename(path)} is not valid JSON: {exc}") from exc
        return cls.from_dict(data)

    # -- description -------------------------------------------------------- #
    def describe(self) -> str:
        """One-line summary, for a recovery prompt or the command line."""
        if self.alignment is None:
            return "an empty session"
        name = os.path.basename(self.source_path) if self.source_path else "unsaved sequences"
        bits = [f"{name}: {self.alignment.n_seqs} sequences x "
                f"{self.alignment.length} columns"]
        if self.reference is not None:
            bits.append(f"reference: {self.reference_label or 'loaded'}")
        history = (self.alignment.meta or {}).get("history") or []
        if history:
            bits.append(f"{len(history)} steps of history")
        return "; ".join(bits)


def is_session_file(path: str) -> bool:
    """Whether ``path`` looks like a session document, cheaply."""
    if not path.endswith(SUFFIX) and not path.endswith(".json"):
        return False
    try:
        with open(path, encoding="utf-8") as fh:
            return '"craic-session"' in fh.read(400)
    except OSError:
        return False
