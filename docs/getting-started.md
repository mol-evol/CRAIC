# Getting started

CRAIC runs on macOS, Linux, and Windows. It needs **Python 3.9+**. The Rust acceleration
core is optional — without it CRAIC falls back to an identical NumPy implementation, just
slower.

## Install and launch

=== "One command (recommended)"

    ```bash
    cd path/to/craic
    ./build_rust.sh      # makes a .venv, installs deps, builds the Rust core if Rust is present
    ./craic.command      # or:  .venv/bin/python -m craic
    ```

=== "By hand (NumPy only)"

    ```bash
    python3 -m venv .venv
    .venv/bin/python -m pip install numpy biopython PySide6
    .venv/bin/python -m craic examples/coding_genes.fasta
    ```

!!! note "macOS, Homebrew, and conda"
    If your `python3` is Homebrew's, a system-wide `pip install` fails with
    *externally-managed-environment* (PEP 668). The virtualenv above avoids that, and
    `build_rust.sh` also unsets `CONDA_PREFIX` for its own run so `maturin` won't refuse to
    build inside a conda shell.

## A double-clickable app (macOS)

Launched with `python -m craic`, macOS shows the Python executable (e.g. *python3.13*) in the
Dock and app switcher, because a bare script isn't an application bundle. To get a proper
**CRAIC** name and icon everywhere — and a Finder/Dock launcher — build a small `.app`:

```bash
python scripts/make_app.py        # writes CRAIC.app in the project root
```

It wires the bundle to the Python you ran it with, so use the interpreter that has CRAIC
installed (e.g. `.venv/bin/python scripts/make_app.py`). Re-run it if you move the project or
change environments. Drag the resulting `CRAIC.app` to your Dock or `/Applications`. (The
bundle is git-ignored; the script is the source of truth.)

The status bar shows **`core: rust`** once the native core is built, otherwise **`core: numpy`**.

## Your first alignment

1. **Open** a file (FASTA, PHYLIP, Clustal, Stockholm, or NEXUS) — or load the bundled
   `examples/coding_genes.fasta`.
2. Choose an **engine** and press **Align**. Tick **align as protein** for protein-coding
   nucleotides — CRAIC translates, aligns the amino acids, and threads the nucleotides
   back through, keeping the reading frame intact and switching to the protein view.
3. Toggle **View** between *nucleotide*, *codon*, and *amino acid*.
4. Use the **Track** dropdown to overlay a per-column statistic (Conservation, Reliability,
   Aligner agreement, …); it computes whatever it needs on demand the first time you pick it.
5. **Click-drag** a column range to drive the realignment sandbox and the posterior explorer;
   **single-click** a residue to probe where else it could align.

## Then try it with the answer known

The fastest way to understand what the reliability score is claiming is to check it against
a case where the right answer exists. **Teach → Generate dataset with known answer…** evolves
sequences down a random tree and loads them *unaligned*. Align them, then compare the
**Reliability** track — a prediction — against the **Reference correctness** track, which is
the outcome. Colour by **Correct placement** to see which individual residues went wrong, and
open the **Column inspector** to read what the right answer was.

See [Teaching with a known answer](guide/teaching.md).

## Optional: external aligners

CRAIC auto-detects **MAFFT**, **MUSCLE**, **Clustal Omega**, and **PRANK** on your `PATH`
and offers them in the engine menu and the disagreement comparison. None are required — the
built-in progressive aligner always works.
