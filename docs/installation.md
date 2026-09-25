# Installation

CRAIC runs on macOS, Windows, and Linux. It needs **Python 3.9 or newer**.

The fast path is a prebuilt wheel from PyPI — no compiler or Rust toolchain required.
The Rust acceleration core is bundled in the wheel; if you install from source without
Rust, CRAIC still runs on a built-in NumPy fallback (the status bar shows `core: numpy`
instead of `core: rust`).

Prefer a **double-clickable app** with no Python at all? See
[Download a ready-to-run app](#download-a-ready-to-run-app) below.

## Download a ready-to-run app

Each [release](https://github.com/mol-evol/craic/releases) includes standalone bundles
that need **no Python and no install** — just download, unpack, and double-click:

| Platform | Download | Run |
| --- | --- | --- |
| macOS | `CRAIC-macos.dmg` | Open the `.dmg`, drag **CRAIC** to Applications |
| Windows | `CRAIC-windows.zip` | Unzip, open the folder, double-click **CRAIC.exe** |
| Linux | `CRAIC-linux.tar.gz` | Extract, then run `./CRAIC/CRAIC` |

These apps are **not code-signed** (that needs paid developer certificates), so the
operating system shows a one-time warning the first time you open one. This is expected —
here's how to get past it:

!!! warning "macOS — first launch"
    The app isn't signed by an Apple-registered developer, so macOS blocks it the
    first time ("Apple could not verify 'CRAIC' is free of malware…"). To allow it,
    once:

    1. Drag **CRAIC** to Applications, double-click it, and click **Done** on the warning.
    2. Open **System Settings → Privacy & Security** and scroll down to **Security**,
       where it says "CRAIC" was blocked. Click **Open Anyway**, enter your password,
       then click **Open Anyway** again.

    After that it opens normally. (On macOS 14 and earlier, right-click the app and
    choose **Open** instead; macOS 15 removed that shortcut.)

    **No warning at all:** if you have Python, install from PyPI (below) and run
    `craic make-app`. It builds the app on your own Mac, so macOS never marks it as
    downloaded.

    If macOS says the app is "damaged" or won't open, clear the quarantine flag in
    Terminal (this is safe — it only removes the download-quarantine attribute):

    ```bash
    xattr -dr com.apple.quarantine /Applications/CRAIC.app
    ```

!!! warning "Windows — first launch"
    Windows SmartScreen may say "Windows protected your PC". Click **More info →
    Run anyway**.

!!! note "Linux"
    If the binary isn't executable after extracting, run `chmod +x CRAIC/CRAIC`.

If you'd rather avoid the warnings entirely, install from PyPI instead (below) — the
`pip` route is never blocked by the OS.

## Install from PyPI (recommended)

```bash
pip install craic-msa
craic                      # launch the GUI
craic path/to/alignment.fasta   # …or open a file straight away
```

We recommend installing into a virtual environment so CRAIC's dependencies (NumPy,
Biopython, PySide6) don't clash with other tools:

```bash
python3 -m venv craic-env
# macOS / Linux:
source craic-env/bin/activate
# Windows (PowerShell):
craic-env\Scripts\Activate.ps1

pip install craic-msa
craic
```

### Platform notes

- **macOS** — works on both Apple Silicon and Intel. If `craic` isn't found after
  install, make sure your virtual environment is activated, or launch with
  `python -m craic`.
- **Windows** — use the same commands in PowerShell or Command Prompt. If `craic`
  isn't on your `PATH`, use `python -m craic`.
- **Linux** — PySide6 needs a few system libraries for Qt. On Debian/Ubuntu:

  ```bash
  sudo apt-get install -y libegl1 libgl1 libxkbcommon0 libdbus-1-3
  ```

## Install from source

You'll need [Rust](https://rustup.rs) to build the acceleration core (or skip it and
use the NumPy fallback).

```bash
git clone https://github.com/mol-evol/craic
cd craic
python3 -m venv .venv && source .venv/bin/activate
pip install maturin
maturin develop --release        # compiles the Rust core into the venv
pip install -e .
craic examples/coding_genes.fasta
```

**No Rust?** Install just the Python dependencies and run on the NumPy fallback:

```bash
pip install numpy biopython PySide6
python -m craic examples/coding_genes.fasta
```

## Optional: external aligners

CRAIC auto-detects these on your `PATH` and offers them in the engine menu and the
multi-aligner disagreement map: **MAFFT**, **MUSCLE**, **Clustal Omega**, **ProbCons**,
**PRANK** and **ClustalW**. None are required — the built-in aligner always works — but
installing one or more gives you stronger alignments and a more informative disagreement
comparison. All six are on Bioconda:

```bash
conda install -c bioconda mafft muscle clustalo probcons prank clustalw
```

## Verifying the install

```bash
craic --version
```

Launch CRAIC and open `examples/coding_genes.fasta`. The status bar at the bottom
shows the active core (`core: rust` or `core: numpy`) and confirms everything loaded.
