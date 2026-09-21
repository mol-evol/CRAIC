# Keyboard & mouse

## Keyboard

| Shortcut | Action |
| --- | --- |
| Ctrl/⌘ + O | Open sequences or an alignment |
| Ctrl/⌘ + S | Save the current alignment |
| Ctrl/⌘ + C | Copy the selection (or whole alignment) as FASTA |
| Ctrl/⌘ + P | Pin the selected columns as trusted |
| Ctrl/⌘ + Shift + P | Unpin the selected columns |
| Ctrl/⌘ + Z | Undo the last edit |
| Ctrl/⌘ + Shift + Z | Redo |
| Ctrl/⌘ + T | Show the true alignment, read-only (needs a reference); again to go back |
| Ctrl/⌘ + `+` | Bigger text (zoom in) |
| Ctrl/⌘ + `-` | Smaller text (zoom out) |
| Ctrl/⌘ + `0` | Reset text size |
| Esc | Clear the column and sequence selection |

The residue font scales with the cells, so zooming makes the letters bigger or
smaller. You can also use **View ▾ → Bigger / Smaller / Reset text size**, or
Ctrl/⌘ + scroll over the alignment.

## Editing with the keyboard

Editing acts on the **selected sequences only**. Click a cell to place the edit
cursor (a bold amber box) and select that one sequence; Cmd/Ctrl-click more cells
(or sequence names) to add them to the selection. Selected sequences are tinted
amber. Then slide them with the keyboard. Sliding works at every view level: one
base at the nucleotide level, or a whole codon (keeping the reading frame) at the
codon / amino-acid level. Every edit is undoable.

| Gesture / Key | Action |
| --- | --- |
| Click a cell or name | Place the cursor and select that one sequence |
| Cmd/Ctrl-click a cell or name | Add / remove a sequence from the selection |
| Shift-click a cell or name | Select a contiguous range of sequences from the last click |
| **Drag a selected (amber) sequence** | Slide it (and any other selected sequences) left/right with the mouse — one column per column you drag |
| **Alt/Option-drag any sequence** | Grab and slide that sequence even if it wasn't selected |
| ← → ↑ ↓ | Move the edit cursor |
| `-` (minus) or Space | Slide the selected sequences one column **right** (opens a gap; absorbs a trailing gap rather than widening the whole alignment) |
| Backspace / Delete / Shift+Space | Slide the selected sequences one column **left** (consumes the gap to the left; refused if any selected sequence has no gap there) |
| Hold `-` / Backspace (etc.) | The key auto-repeats, so holding it slides continuously |

!!! tip
    To realign one sequence against the rest, click it and either drag it into
    place or tap / hold `-` and Backspace. To move several at once, Cmd-click each
    one first (or Shift-click a range) — they slide together, keeping their
    alignment relative to one another, while the unselected sequences stay put.

### Sequence groups

If you keep moving the same set of sequences together (say, all the animals),
save them as a named group: select them, then **Edit ▸ Sequence groups ▸ Save
selection as group…**. Later, **Edit ▸ Sequence groups ▸ Select “name”** re-selects
all its members in one click, ready to slide together. Groups are stored in a small
companion file next to the alignment (`<file>.craic.json`) — the alignment file
itself stays a standard FASTA / Stockholm / etc. — and are restored when you reopen
the alignment.

## Mouse

| Gesture | Action |
| --- | --- |
| Scroll | Pan the alignment vertically |
| Ctrl/⌘ + scroll | Zoom the cells and text in and out |
| Click-drag across columns | Select a column range (drives the sandbox and posterior explorer) |
| Single-click a residue | Probe it (paint the columns it could also align to) and place the edit cursor there |
| Drag a dock title bar | Move or float the **Ambiguity tools** panel |

!!! tip
    A single click and a click-drag do different things: a *click* probes one residue without
    changing your column selection; a *drag* sets the selection that both ambiguity panels use.
