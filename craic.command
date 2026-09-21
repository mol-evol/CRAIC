#!/usr/bin/env bash
# CRAIC launcher — always runs with the project's bundled virtualenv, so it
# works no matter which Python (conda, Homebrew, system) is active in your shell.
# Double-click in Finder, or run ./craic.command [optional_file] in a terminal.
DIR="$(cd "$(dirname "$0")" && pwd)"
# Match build_rust.sh: the venv lives outside cloud-synced folders.
case "$DIR" in
  *CloudStorage*|*Dropbox*|*OneDrive*|*"Google Drive"*|*"Mobile Documents"*)
    PY="$HOME/.craic/venv-v0.5/bin/python" ;;
  *)
    PY="$DIR/.venv/bin/python" ;;
esac

if [ ! -x "$PY" ]; then
  echo "CRAIC isn't set up yet: no .venv found in"
  echo "  $DIR"
  echo "Run  ./build_rust.sh  first, then try again."
  read -r -p "Press Return to close… " _
  exit 1
fi

cd "$DIR"
exec "$PY" -m craic "$@"
