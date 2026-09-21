# Contributing to CRAIC

Thanks for your interest in improving CRAIC. Bug reports, feature ideas, and pull
requests are all welcome.

## Reporting bugs / requesting features

Open an issue at https://github.com/mol-evol/craic/issues. For bugs, please include:

- your OS and Python version,
- how you installed CRAIC (PyPI wheel, or from source),
- what you did, what you expected, and what happened,
- the contents of the status bar (it shows `core: rust` or `core: numpy`), and a
  small example alignment if you can share one.

## Development setup

CRAIC is a Python package with a small Rust acceleration core (built with
[maturin](https://www.maturin.rs/)). The Rust core is optional at runtime — there is
a NumPy fallback — but you need a Rust toolchain to build it.

```bash
git clone https://github.com/mol-evol/craic
cd craic
python3 -m venv .venv && source .venv/bin/activate
pip install maturin
maturin develop --release          # builds the Rust core into the venv
pip install -e ".[test,docs]"
python -m craic examples/coding_genes.fasta
```

No Rust toolchain? Install the Python dependencies only and CRAIC will run on the
NumPy fallback:

```bash
pip install numpy biopython PySide6
python -m craic examples/coding_genes.fasta
```

## Running the tests

```bash
pytest tests
```

The GUI tests run headlessly (`QT_QPA_PLATFORM=offscreen`) and are skipped
automatically if PySide6 is unavailable. Please keep the suite green and add a test
for any behaviour you change.

## Building the docs

```bash
pip install mkdocs-material
mkdocs serve        # live preview at http://127.0.0.1:8000
mkdocs build --strict
```

## Pull requests

- Keep changes focused and the test suite passing.
- Match the existing style: clear, simple code (KISS), no speculative features.
- Update the docs and add tests alongside code changes.

By contributing, you agree that your contributions are licensed under the MIT License.
