#!/usr/bin/env bash
# Set up CRAIC in a self-contained virtualenv and build the optional Rust core.
#
# Using a venv sidesteps macOS / Homebrew "externally-managed-environment"
# (PEP 668) errors and keeps CRAIC isolated from your system Python. Re-running
# is safe — it reuses the existing .venv.
#
#   ./build_rust.sh
#
# When it finishes it prints the exact command to launch CRAIC.
set -e
cd "$(dirname "$0")"
PROJ="$(pwd)"

# Keep the virtualenv OFF cloud-synced folders (Dropbox / OneDrive / iCloud /
# Google Drive). Cloud sync dehydrates the venv's thousands of files, and then
# every Python startup hangs re-downloading them — which looks exactly like the
# script freezing. The code can live in the cloud; the environment must be local.
case "$PROJ" in
  *CloudStorage*|*Dropbox*|*OneDrive*|*"Google Drive"*|*"Mobile Documents"*)
    VENV="$HOME/.craic/venv-v0.5"
    echo "Project is in a cloud-synced folder — putting the venv in $VENV instead." ;;
  *)
    VENV="$PROJ/.venv" ;;
esac
mkdir -p "$(dirname "$VENV")"

# 1. virtualenv (activation-free: we drive it via VIRTUAL_ENV + PATH so this
#    works under sh, bash, or zsh)
if [ ! -d "$VENV" ]; then
  echo "Creating virtualenv in .venv ..."
  python3 -m venv "$VENV"
fi
export VIRTUAL_ENV="$VENV"
export PATH="$VENV/bin:$PATH"
# maturin refuses to run when both a venv and a conda env look active. We have
# deliberately chosen the venv, so neutralise conda for this script's run only
# (this does not touch your interactive shell).
unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER CONDA_SHLVL 2>/dev/null || true

# 2. runtime dependencies — install only what's missing. A blanket reinstall of
#    PySide6 is a large (~100 MB) download and is the slow part of re-running this
#    script, so skip it when the imports already work.
if python -c "import numpy, Bio, PySide6" 2>/dev/null; then
  echo "Dependencies already present in .venv — skipping install."
else
  echo "Installing dependencies into .venv ..."
  echo "(PySide6 is a large download — the first run can take a few minutes; let it finish.)"
  python -m pip install --upgrade pip
  python -m pip install numpy biopython PySide6
fi

# 3. optional Rust acceleration core
if command -v cargo >/dev/null 2>&1; then
  echo "Rust found — building the acceleration core (maturin) ..."
  # Build artifacts (thousands of files) also stay off any cloud-synced folder.
  export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$HOME/.craic/rust-target}"
  python -m pip install maturin
  maturin develop --release
else
  echo "Rust not found — skipping the native core (CRAIC will use the NumPy"
  echo "fallback). To get the ~50x speedup later: install Rust from"
  echo "https://rustup.rs and re-run this script."
fi

# 4. verify and report
python - <<'PY'
from craic import accel
print("CRAIC core backend:", accel.backend())
PY

echo
echo "Setup complete. Launch CRAIC with:"
echo "    \"$VENV/bin/python\" -m craic \"$PROJ/examples/coding_genes.fasta\""
echo "or just:"
echo "    \"$VENV/bin/python\" -m craic"
