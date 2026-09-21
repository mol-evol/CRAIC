"""The column overlays, described once.

A "track" is a per-column statistic painted above the alignment. There are a
dozen, and they used to be three parallel encodings of the integers 0 to 11: a
list of labels for the combo box, a set of which ones were cheap enough to
recompute after an edit, and an if-chain dispatching on the index. Inserting a
track in the middle silently broke the other two, and the author of this module
did exactly that at least once.

So each track is one object that knows its own label, what it costs, and how it
is produced. The combo box, the dispatch and the keep-through-an-edit rule all
read from the same list.

Tracks differ in *how* they are produced, not just in what they compute, which is
what ``kind`` records:

``score``      a per-column array, computed immediately from the alignment
``mask``       a boolean keep-mask chosen by the method itself (the rule-based
               trimmers), shown over conservation for context
``reliability``  needs the reliability report, which runs in the background
``agreement``    needs the multi-aligner comparison, also background
``truth``        needs a reference alignment to have been loaded
``none``         no overlay
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional


@dataclass(frozen=True)
class Track:
    """One column overlay."""

    key: str
    label: str
    kind: str
    #: For ``score`` and ``mask`` tracks: how to compute it from an Alignment.
    compute: Optional[Callable] = None
    #: Which field of the reliability report a ``reliability`` track shows.
    field: str = ""

    @property
    def cheap(self) -> bool:
        """Can this be recomputed after every edit without making editing
        unusable? Score and mask tracks take a fraction of a second; the
        reliability family and the aligner-agreement map take seconds to minutes
        and run in the background, so they are dropped on an edit and recomputed
        when asked for again."""
        return self.kind in ("score", "mask", "truth")


def build(conservation, trimming) -> List[Track]:
    """The track list. ``conservation`` is the per-column conservation function
    and ``trimming`` the rule-based trimming module; both are passed in so this
    module stays free of imports from the window it serves."""
    return [
        Track("none", "(none)", "none"),
        Track("conservation", "Conservation", "score", conservation),
        Track("reliability", "Reliability", "reliability", field="col_combined"),
        Track("consistency", "Consistency", "reliability", field="col_consistency"),
        Track("perturbation", "Perturbation", "reliability", field="col_perturbation"),
        Track("agreement", "Aligner agreement", "agreement"),
        Track("gap", "Gap fraction (MSA_trimmer)", "score", trimming.gap_score),
        Track("similarity", "Similarity (trimAl)", "score", trimming.similarity_score),
        Track("gappyout", "trimAl gappyout", "mask", trimming.gappyout_mask),
        Track("strict", "trimAl strict", "mask", trimming.trimal_strict_mask),
        Track("gblocks", "Gblocks", "mask", trimming.gblocks_mask),
        Track("truth", "Reference correctness (truth mode)", "truth"),
    ]


def index_of(tracks: List[Track], key: str) -> int:
    """Position of a track by key, so nothing has to name an integer."""
    for i, t in enumerate(tracks):
        if t.key == key:
            return i
    raise KeyError(f"no such track: {key!r}")
