# Sessions and unsaved work

## Why a session is not an alignment file

What you have at the end of an hour's curation is more than a set of aligned
residues. It is also:

- the reference alignment truth mode is scoring against
- the column annotations you marked up
- the named sequence groups
- the **provenance log** — the record of which engine produced the alignment and
  every edit made since, which is what makes a hand-built alignment reproducible
- the view state: which overlay, what threshold, which level you were reading at

FASTA holds the first of those. NEXUS holds the residues and a charset. Neither
holds the rest, which is exactly why CRAIC does not auto-save to them on the way
out: it would preserve the cheap part and silently drop the part that took the
judgement.

## Saving a session

**File → Save session** (`Ctrl+S`) writes a `.craic.json` document holding all of
it. **Open session…** restores it — and a session also opens from the ordinary
Open dialog, since CRAIC recognises one by its content.

`Ctrl+S` deliberately saves the *session* rather than exporting the alignment,
because it is the operation that loses nothing, and so is the right thing for a
reflex to reach for. Exporting to a standard format is **Save alignment as…**
(`Ctrl+Shift+S`) and is unchanged: that is what you hand to RAxML, IQ-TREE or a
collaborator.

The file is plain JSON. You can read it, diff it, and put it in version control
alongside the data.

## Unsaved changes

A `•` appears in the title bar as soon as there are changes that are not on disk.
Closing with unsaved work asks whether to save, discard, or stay — and says which
of the two kinds of save keeps what.

Exporting the alignment to a standard format also clears the mark: the residues,
groups and annotations are on disk (the last two in the sidecar file beside it).
If a truth-mode reference is loaded, CRAIC says that the reference is not part of
a standard format, rather than letting you discover it tomorrow.

Opening another file, opening another session and generating a teaching dataset
all replace the document just as completely as quitting does, so each asks the
same question first. Realigning does not — but it is undoable, which it needs to
be, since realigning discards every hand edit made since the last alignment.

## If CRAIC stops unexpectedly

The session is autosaved a few seconds after each change into CRAIC's own
application-data directory — **never** beside your data, so opening someone's
reference alignment and looking at it leaves nothing behind in their folder.

On the next launch, if CRAIC did not exit cleanly, it offers the interrupted
session back, describing what it holds and when it was last touched. Recovered
work is still marked unsaved, because it is: save it somewhere deliberate.

A clean exit removes the recovery file, so being offered one always means
something actually went wrong.

## From the command line

```bash
craic session work.fasta.craic.json
```

prints what the session holds, including every step of the provenance log; `-o`
exports its alignment to any supported format.
