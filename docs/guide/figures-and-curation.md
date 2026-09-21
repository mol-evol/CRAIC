# Figures, reports & curation

Everything CRAIC knows about alignment uncertainty can now be turned into a
publication-quality figure, a QC report, or a hands-on edit. All figures are drawn
through one vector engine, so each exports identically to **PNG, SVG, or PDF** (and,
for the alignment view, an interactive **HTML** page) — no extra software needed.

## Exporting figures

**File ▾ → Export figure ▸** offers:

| Figure | What it shows |
| --- | --- |
| **Alignment figure** | The wrapped block alignment with residues *faded where the column is unreliable* — the trustworthy core stands out and ambiguous stretches recede. Includes a reliability strip and consensus. |
| **Confidence report** | A per-column reliability profile with the ambiguous **hotspots** shaded, plus a text QC summary (mean reliability, % high-confidence columns, the hotspot coordinates). |
| **Sequence logo (selection)** | An *uncertainty-weighted* logo: confident columns stand tall, ambiguous ones shrink into multi-letter stacks. |
| **Homology arcs (explorer pair)** | The two posterior-explorer sequences as tracks, with arcs linking residues that could be homologous — opacity ∝ posterior. A confident region is a tight braid; an ambiguous one fans out. |
| **Co-occurrence map (selection)** | How often each pair of residues ends up aligned across the *ensemble* of alternative alignments — the crowd's view of the region, complementing the single-model posterior. |
| **Homology card (last-clicked residue)** | The "cloud" for one residue: which columns it could occupy and with what confidence. Click a residue first to choose it. |
| **Interactive HTML** | A self-contained, shareable page: coloured residues, a reliability strip, sticky names, and hover tooltips. |

Figures use the level you're viewing (nt / codon / amino acid) and, where relevant,
the current column selection.

## Navigating ambiguity

**Edit ▾ → Go to next ambiguous region** (Ctrl/⌘+G) jumps the selection to each
reliability hotspot in turn, worst first — so you go straight to the columns that
need a human instead of scrolling.

## Comparing two alignments

**File ▾ → Compare with alignment…** loads another alignment of the same sequences
and shows, as the column track, how often the two place residues together — red
where they disagree. A natural extension of the multi-aligner agreement track to an
*external* alignment (e.g. yours vs a collaborator's, or before/after editing).

## Manual editing

The fastest way to curate is the keyboard, and it acts on the **selected
sequences only** — so you can realign one (or a few) sequences without disturbing
the rest. Click a cell to place the edit cursor (a bold amber box) and select that
sequence; Cmd/Ctrl-click more cells or names to add them, or Shift-click to grab a
contiguous range. Selected sequences slide together and stay fixed relative to one
another. If you reuse the same set often, save it as a named **sequence group**
(Edit ▸ Sequence groups) and re-select it in one click. Then:

| Key | Action |
| --- | --- |
| ← → ↑ ↓ | Move the cursor |
| `-` (minus) or Space | Slide the selected sequences right (open a gap; a trailing gap is absorbed, so the alignment doesn't widen) |
| Backspace / Delete / Shift+Space | Slide the selected sequences left (consume the gap to the left, if every selected sequence has one) |

Holding `-` or Backspace auto-repeats, so you can slide continuously. You can also
**drag** a selected (amber) sequence with the mouse to slide it — one column per
column dragged — or Alt/Option-drag any sequence to grab it directly.

The unselected sequences stay exactly where they are. Editing works at **every**
view level: at the nucleotide level a slide moves one base, and at the codon or
amino-acid level it moves a whole codon (three bases), keeping the reading frame
intact. The view you're in stays put after each edit.

The same operations are also on **Edit ▾ → Edit alignment** (handy for discovery),
all recorded on the provenance log and fully undoable:

| Action | |
| --- | --- |
| **Undo / Redo** | Ctrl/⌘+Z, Ctrl/⌘+Shift+Z — covers every edit, splice, sort, and degap. |
| **Nudge residue left / right** | Ctrl/⌘+, and Ctrl/⌘+. — slide the last-clicked residue into an adjacent gap; the status bar shows the column's support score change so you see whether the edit helped. |
| **Insert gap column** | Open a gap in every sequence at the selected column. |
| **Delete empty column** | Remove an all-gap column. |
| **Pin selected columns** (`Ctrl/⌘ P`) | Mark the selected columns as trusted. Pinned columns are drawn with a gold bar above the grid. |
| **Unpin selected columns** (`Ctrl/⌘ ⇧ P`) / **Clear all pins** | Remove pins. |
| **Realign around pinned columns** | Keep the pinned columns fixed and let the engine re-solve only the uncertain stretches between them. The pins move with their columns, so the blocks you trust stay pinned afterwards. |

Manual edits to gap structure drop the coding annotation (a hand-made gap can break
the reading frame); residue nudges keep it, since they don't move columns.

## Column annotations

Label regions of the alignment — α-helices, active sites, domains, anything — as
named, coloured **column features**. Select a column range (drag), then **Edit ▾ →
Annotations → Add annotation from selection…** and give it a name. Annotations show
in a labelled colour band along the bottom of the alignment.

* **Fixed by position.** An annotation marks a column range; if you later edit the
  alignment, it stays at those columns (it does not follow the residues).
* **Select to highlight.** **Edit ▾ → Annotations → Select “…”** highlights that
  feature's columns and scrolls to it.
* **Export the regions only.** **Export annotated columns only…** writes a new
  alignment containing just the columns covered by your annotations (the union).
* **Saved with the alignment.** Annotations live in the same sidecar file as
  sequence groups (`<file>.craic.json`), so the alignment file stays standard and
  the annotations come back when you reopen it.

Annotations and the conservation **Mask** slider are independent channels and can
be used together: the slider dims low-scoring columns and exports by score, while
annotations colour-band specific regions and export by region — pick whichever
export you want from its own menu, so the two never conflict.
