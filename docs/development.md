# Architecture & development

## Layout

```
craic/
  accel.py          dispatch to the Rust core or the NumPy mirror (identical results)
  domain.py         canonical model + nt / codon / aa projection, genetic codes
  io.py             FASTA / PHYLIP / Clustal / Stockholm / NEXUS
  progressive.py    the built-in progressive aligner
  engines.py        AlignerEngine interface + external wrappers + codon-aware decorator
  ambiguity/        reliability · disagreement · sandbox · posterior
  gui/              PySide6 front-end (canvas, tracks, panels, workers)
src/lib.rs          the pair-HMM core (PyO3)
tests/              pytest suite (core + GUI smoke tests)
```

The GUI and the ambiguity tools depend only on small interfaces (`AlignerEngine`, the
`accel` dispatch), so pieces are swappable. The Rust core and its NumPy fallback satisfy the
same contract — callers can't tell which answered.

## Build and test

```bash
./build_rust.sh                      # venv + deps + Rust core
.venv/bin/python -m pytest tests     # core tests; GUI smoke tests run if PySide6 imports
```

The Rust core builds via `maturin` into `craic._accel`. Without Rust, everything still runs on
the NumPy mirror.

## How the pieces fit

Everything in CRAIC reduces to one primitive and three layers on top of it.

```
              craic/accel.py  ──  dispatch
                 │                    │
        src/lib.rs (Rust)      NumPy mirror        identical results,
        via PyO3                                   cross-validated to ~1e-9
                 └────────┬───────────┘
                          │
                pair-HMM posteriors  P[i,j]
                          │
     ┌────────────────────┼────────────────────┐
     │                    │                    │
  progressive.py     ambiguity/           evaluate.py
  (MEA aligner)      reliability,         (accuracy against
                     posterior,            a known answer)
                     sandbox,
                     disagreement,
                     trimming
     │                    │                    │
     └────────────────────┼────────────────────┘
                          │
                  domain.py  ──  the canonical model
                          │
            ┌─────────────┴─────────────┐
         gui/                          cli.py
      (PySide6)                     (never imports Qt)
```

### The primitive

`accel.pair_posteriors` returns the matrix of match probabilities for two
sequences under a three-state pair-HMM. Every uncertainty measure in the program
is a different reduction of that matrix, which is what keeps the numbers mutually
consistent: the score the posterior explorer draws, the score the reliability
track shows and the score the masking slider thresholds are all the same
quantity, not three implementations that happen to agree.

The Rust core releases the GIL for the duration of the kernel (`py.allow_threads`),
so the GUI stays responsive while an analysis runs — measured at 178 time slices
delivered to another Python thread during one 566 ms kernel call.

### The canonical model

`domain.py` holds one representation — gapped rows over an alphabet, plus an
optional coding annotation — and *projects* it to whatever level is being viewed:
nucleotides, codons, or the encoded amino acids. Nothing downstream keeps a
second copy, so switching the view cannot desynchronise anything, and an edit made
in the amino-acid view is an edit to the nucleotides underneath.

It also owns the residue-index API that everything comparing two alignments needs:

| function | meaning |
|---|---|
| `residue_index(row)` | column → index in the ungapped sequence, or `-1` |
| `column_index(row)` | the inverse: residue → column |
| `membership(aln)` | `(sequence, residue) → column` for the whole alignment |

These were reimplemented in five places before 0.4.1. They are one place now, and
anything that walks the correspondence between an alignment's columns and its
residues should call them rather than write the loop again.

### Accuracy lives in one module

`evaluate.py` computes SP, TC, per-column correctness, per-residue correctness and
the reliability AUC, and it is used by the GUI's truth mode, the command line and
the benchmark harness alike. This is deliberate: a number shown to a student and a
number printed in a paper being produced by two different implementations is
exactly the drift worth designing out. `benchmarks/metrics.py` re-exports these
rather than defining its own.

`true_column_grid(inf_rows, true_rows)` is the object the whole teaching layer
reads: entry `[s, c]` is the column the reference puts the residue that the
inferred alignment placed at row `s`, column `c`. Two residues are homologous
exactly when their entries are equal.

### Two front ends, one of which must stay headless

`cli.py` never imports PySide6, and there is a test that asserts it — a headless
subcommand that pulls in Qt fails on a cluster node with no display, which is
precisely where the command line is wanted. Anything shared between the GUI and
the CLI therefore belongs below `gui/`, not inside it.

## Extending it

**Add an alignment engine.** Subclass `AlignerEngine` in `craic/engines.py`, implement
`available()` and `align(records, alphabet)`, and add it to `all_engines()`. It immediately
appears in the engine menu, the sandbox, and the disagreement comparison.

**Add an ambiguity view.** Add a module under `craic/ambiguity/` that consumes
`accel.posterior_matrix(...)` and returns plain arrays; wire a panel under `craic/gui/`.

## Optional: API reference from docstrings

The code is documented with module/function docstrings, so you can auto-generate an API
reference. Install the plugin and enable it:

```bash
.venv/bin/python -m pip install "mkdocstrings[python]"
```

Then uncomment the `mkdocstrings` block in `mkdocs.yml`, add a page such as `docs/reference.md`
containing:

```markdown
# API reference
::: craic.ambiguity.posterior
```

and add it to the `nav`. Run docs from the project venv so the package is importable.


## The GUI, structurally

`craic/gui/` separates what the user is working on from how it is drawn:

- **`document.py`** — the `Document`: alignment, the reference truth mode scores
  against and the resulting accuracy, groups, annotations, file paths, the
  unsaved-changes flag, and a generation counter that increases whenever the
  alignment is replaced. It emits signals; it knows nothing about widgets.
- **`tracks.py`** — every column overlay described once: its label, its kind
  (score, mask, reliability, agreement, truth) and how it is computed. The combo
  box, the dispatch and the rule about which overlays survive an edit all read
  from this one list.
- **`app.py`** — the window. Owns widgets, menus and the wiring between them.
  `win.aln`, `win.reference`, `win.truth` and `win._dirty` are thin properties
  onto the document.
- **`canvas.py`**, **`panels.py`**, **`dialogs.py`**, **`colors.py`** — the
  painted alignment view, the dockable ambiguity tools, the small forms, and the
  colour schemes.

Two rules worth keeping:

1. **Background work captures `doc.generation` when it starts** and its result is
   discarded if the generation has moved on. A reliability report computed for
   one alignment must never be painted onto another.
2. **Anything that replaces the document calls `_confirm_discard` first.** There
   are four such paths — quit, open file, open session, generate dataset — and
   they are the same question, so they share one prompt.
