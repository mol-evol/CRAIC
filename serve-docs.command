#!/usr/bin/env bash
# Live-preview the CRAIC documentation. Uses the project's .venv and installs
# mkdocs-material into it on first run. Double-click in Finder or run ./serve-docs.command
DIR="$(cd "$(dirname "$0")" && pwd)"
# Match build_rust.sh: the venv lives outside cloud-synced folders.
case "$DIR" in
  *CloudStorage*|*Dropbox*|*OneDrive*|*"Google Drive"*|*"Mobile Documents"*)
    PY="$HOME/.craic/venv-v0.5/bin/python" ;;
  *)
    PY="$DIR/.venv/bin/python" ;;
esac

if [ ! -x "$PY" ]; then
  echo "CRAIC isn't set up yet — run ./build_rust.sh first."
  read -r -p "Press Return to close… " _
  exit 1
fi

if ! "$PY" -m mkdocs --version >/dev/null 2>&1; then
  echo "Installing mkdocs-material ..."
  "$PY" -m pip install mkdocs-material
fi

# The API reference page is generated from the docstrings, which needs
# mkdocstrings; without it mkdocs refuses the config outright.
if ! "$PY" -c "import mkdocstrings" >/dev/null 2>&1; then
  echo "Installing mkdocstrings ..."
  "$PY" -m pip install "mkdocstrings[python]"
fi

cd "$DIR"
echo "Serving docs at http://127.0.0.1:8000  (Ctrl-C to stop)"
exec "$PY" -m mkdocs serve
