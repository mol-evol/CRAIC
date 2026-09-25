"""Build a double-clickable CRAIC.app on this Mac (``craic make-app``).

The bundle is a small launcher that runs this same Python (``python -m craic``),
so it needs CRAIC installed here, e.g. with ``pip install craic-msa``. It is made
on the user's own machine, so macOS never marks it as downloaded and Gatekeeper
does not block it, unlike the ready-made app on the Releases page. It also gives
CRAIC its own name and icon in the Dock and app switcher, where a bare
``python -m craic`` shows the Python executable.

Re-run it after moving or deleting the Python environment it points at.
"""
from __future__ import annotations

import os
import plistlib
import subprocess
import sys
import tempfile
from pathlib import Path

APP_NAME = "CRAIC"
BUNDLE_ID = "io.github.mol-evol.craic"
DEFAULT_DEST = Path.home() / "Applications"


def _render_icns(dest_icns: Path) -> bool:
    """Render the app icon into a proper .icns via macOS ``iconutil``."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    try:
        from PySide6.QtCore import Qt
        from PySide6.QtWidgets import QApplication

        from .gui.app import _make_app_icon
    except Exception as exc:
        print(f"  ! could not import Qt to draw the icon: {exc}")
        return False

    # A QApplication, not a QGuiApplication: a later window in the same process
    # needs the widgets flavour, and Qt allows only one application object.
    QApplication.instance() or QApplication([])
    master = _make_app_icon(1024)
    with tempfile.TemporaryDirectory() as td:
        iconset = Path(td) / "craic.iconset"
        iconset.mkdir()
        for base in (16, 32, 128, 256, 512):
            for scale, suffix in ((1, ""), (2, "@2x")):
                px = base * scale
                img = master.scaled(px, px, Qt.KeepAspectRatio, Qt.SmoothTransformation)
                img.save(str(iconset / f"icon_{base}x{base}{suffix}.png"), "PNG")
        try:
            subprocess.run(["iconutil", "-c", "icns", str(iconset), "-o", str(dest_icns)],
                           check=True)
        except (subprocess.CalledProcessError, FileNotFoundError) as exc:
            print(f"  ! iconutil failed ({exc}); the app will use a default icon")
            return False
    return dest_icns.exists() and dest_icns.stat().st_size > 0


def build(dest: Path) -> tuple:
    """Write ``CRAIC.app`` into ``dest``; returns (app path, whether the icon went in)."""
    from . import __version__

    app = dest / f"{APP_NAME}.app"
    contents = app / "Contents"
    macos = contents / "MacOS"
    resources = contents / "Resources"
    for d in (macos, resources):
        d.mkdir(parents=True, exist_ok=True)

    launcher = macos / APP_NAME
    launcher.write_text(
        "#!/bin/bash\n"
        # A Finder-launched app has no terminal, so keep a log of why it failed.
        'exec >>"$HOME/Library/Logs/CRAIC.log" 2>&1\n'
        'echo "--- $(date) launching CRAIC"\n'
        # Run from $HOME: if CRAIC was started from a folder in Dropbox / iCloud /
        # Documents, macOS privacy protection would deny a Finder-launched app
        # access to it and the app would die silently.
        'cd "$HOME"\n'
        f'exec "{sys.executable}" -m craic "$@"\n')
    launcher.chmod(0o755)

    have_icon = _render_icns(resources / "craic.icns")
    info = {
        "CFBundleName": APP_NAME,
        "CFBundleDisplayName": APP_NAME,
        "CFBundleExecutable": APP_NAME,
        "CFBundleIdentifier": BUNDLE_ID,
        "CFBundlePackageType": "APPL",
        "CFBundleInfoDictionaryVersion": "6.0",
        "CFBundleShortVersionString": __version__,
        "CFBundleVersion": __version__,
        "LSMinimumSystemVersion": "10.14",
        "NSHighResolutionCapable": True,
    }
    if have_icon:
        info["CFBundleIconFile"] = "craic"
    with open(contents / "Info.plist", "wb") as fh:
        plistlib.dump(info, fh)
    (contents / "PkgInfo").write_text("APPL????")
    return app, have_icon


def _warn_if_editable() -> None:
    """``pip install -e`` leaves CRAIC in its project folder. That works from a
    terminal, but a Finder-launched app gets no access to Dropbox / iCloud /
    Documents under macOS privacy protection, so it would exit with no visible
    error."""
    import craic

    src = Path(craic.__file__).resolve()
    if Path(sys.prefix).resolve() not in src.parents:
        print(f"  ! CRAIC is installed editable, from {src.parent}")
        print("    The app may not be able to read that folder when launched from Finder.")
        print(f'    To fix:  "{sys.executable}" -m pip install --no-deps --force-reinstall <project>')


def _report_engines() -> None:
    """Which external aligners the app will find."""
    from .engines import available_engines

    names = [e.label for e in available_engines()]
    print(f"  aligners found: {', '.join(names)}")


def _codesign(app: Path) -> None:
    """Ad-hoc sign the bundle. Without a signature macOS cannot attribute the app
    to anything, so privacy protection silently denies it access to Dropbox /
    iCloud / Documents / Desktop; an ad-hoc signature gives it an identity, so
    macOS asks the user instead."""
    try:
        subprocess.run(["codesign", "--force", "--deep", "--sign", "-", str(app)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        err = getattr(exc, "stderr", b"") or b""
        print(f"  ! ad-hoc signing failed ({err.decode(errors='replace').strip() or exc})")
        print("    macOS may deny the app access to Dropbox/Documents without asking.")


def _refresh_finder(app: Path) -> None:
    """Register the bundle with LaunchServices and bump its mtime so Finder shows
    the new icon instead of a cached generic one."""
    lsregister = ("/System/Library/Frameworks/CoreServices.framework/Frameworks/"
                  "LaunchServices.framework/Support/lsregister")
    try:
        subprocess.run([lsregister, "-f", str(app)], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.utime(app, None)
    except Exception:
        pass


def main(dest: Path = DEFAULT_DEST) -> int:
    if sys.platform != "darwin":
        print("craic make-app builds a macOS app, so it only runs on a Mac. On Windows "
              "and Linux, start CRAIC with the 'craic' command.", file=sys.stderr)
        return 1
    dest.mkdir(parents=True, exist_ok=True)
    _warn_if_editable()
    app, have_icon = build(dest)
    _report_engines()
    _codesign(app)
    _refresh_finder(app)
    print(f"\nBuilt {app}")
    if have_icon:
        print("  If Finder shows a generic icon it is caching the old one; run: killall Finder")
    else:
        print("  The icon could not be made (see above); the app still runs.")
    print("Open it from there, or drag it onto your Dock. It runs this Python:")
    print(f"  {sys.executable}")
    print("so run 'craic make-app' again if you move or delete that environment.")
    return 0
