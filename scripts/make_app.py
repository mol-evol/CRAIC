#!/usr/bin/env python3
"""Build a double-clickable CRAIC.app bundle on macOS.

A ``python -m craic`` launch shows the Python executable ("python3.13") in the
Dock and Cmd-Tab because it is not an app bundle. This wraps CRAIC in a minimal
.app so it reads "CRAIC" everywhere and can be launched from Finder / the Dock.

Run this with the interpreter CRAIC is installed in; the bundle launches that
same interpreter. CRAIC must be installed into it normally (not editable) --
see the warning this script prints if it is not.

    python scripts/make_app.py [--dest DIR]
"""
from __future__ import annotations

import argparse
import plistlib
import subprocess
import sys
import tempfile
from pathlib import Path

APP_NAME = "CRAIC"
BUNDLE_ID = "io.github.mol-evol.craic"
ROOT = Path(__file__).resolve().parent.parent


def _version() -> str:
    try:
        from craic import __version__
        return __version__
    except Exception:
        return "0.1.0"


def _render_icns(dest_icns: Path) -> bool:
    """Render the runtime app icon into a proper .icns via macOS ``iconutil``."""
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    try:
        from PySide6.QtCore import Qt
        from PySide6.QtGui import QGuiApplication

        from craic.gui.app import _make_app_icon
    except Exception as exc:
        print(f"  ! could not import Qt/CRAIC to render the icon: {exc}")
        return False

    QGuiApplication.instance() or QGuiApplication([])
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
    return True


def build(dest: Path) -> Path:
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
        # Run from $HOME and import CRAIC from the interpreter's own
        # site-packages. Nothing in the project folder is touched: if the project
        # sits in Dropbox / iCloud / Documents, macOS privacy protection denies a
        # Finder-launched app access to it and the app would die silently.
        'cd "$HOME"\n'
        f'exec "{sys.executable}" -m craic "$@"\n')
    launcher.chmod(0o755)

    icns = resources / "craic.icns"
    have_icon = _render_icns(icns)
    if have_icon and not (icns.exists() and icns.stat().st_size > 0):
        have_icon = False

    info = {
        "CFBundleName": APP_NAME,
        "CFBundleDisplayName": APP_NAME,
        "CFBundleExecutable": APP_NAME,
        "CFBundleIdentifier": BUNDLE_ID,
        "CFBundlePackageType": "APPL",
        "CFBundleInfoDictionaryVersion": "6.0",
        "CFBundleShortVersionString": _version(),
        "CFBundleVersion": _version(),
        "LSMinimumSystemVersion": "10.14",
        "NSHighResolutionCapable": True,
    }
    if have_icon:
        info["CFBundleIconFile"] = "craic"
    with open(contents / "Info.plist", "wb") as fh:
        plistlib.dump(info, fh)
    (contents / "PkgInfo").write_text("APPL????")
    return app, have_icon, (icns.stat().st_size if have_icon else 0)


def _check_not_editable() -> None:
    """Warn if CRAIC is installed editable, i.e. still living in the project folder.

    ``maturin develop`` and ``pip install -e`` leave the package where it was
    built. That is fine from a terminal, but a Finder-launched app gets no access
    to Dropbox / iCloud / Documents under macOS privacy protection, so it would
    exit immediately with no visible error. Call this BEFORE adding the project
    to sys.path, or the check answers itself.
    """
    try:
        import craic
        src = Path(craic.__file__).resolve()
    except Exception:
        print("  ! CRAIC is not importable from this interpreter.")
        print(f'    Run:  "{sys.executable}" -m pip install --no-deps --force-reinstall .')
        return
    if Path(sys.prefix).resolve() not in src.parents:
        print(f"  ! CRAIC is installed editable, from {src.parent}")
        print("    The app will not be able to read that folder when launched from Finder.")
        print(f'    Run:  "{sys.executable}" -m pip install --no-deps --force-reinstall .')


def _report_engines() -> None:
    """List which external aligners the bundle will be able to see."""
    from shutil import which
    found = {t: which(t) for t in ("mafft", "muscle", "probcons", "prank", "clustalo")}
    found["clustalw"] = which("clustalw2") or which("clustalw")
    ok = [t for t, p in found.items() if p]
    missing = [t for t, p in found.items() if not p]
    print(f"  engines on PATH: {', '.join(ok) if ok else 'none'}")
    if missing:
        print(f"  not on PATH:     {', '.join(missing)}")
        print("    These will be greyed out in the app. Re-run this script from a shell")
        print("    where they are on PATH (e.g. with the right conda env activated).")


def _codesign(app: Path) -> None:
    """Ad-hoc sign the bundle.

    Without a code signature macOS cannot attribute the app to anything, so TCC
    (privacy) denies it access to protected folders -- Dropbox / iCloud /
    Documents / Desktop -- *silently*, with no prompt. The app then appears to do
    nothing at all. An ad-hoc signature gives it a stable identity, so macOS asks
    the user instead of refusing.
    """
    try:
        subprocess.run(["codesign", "--force", "--deep", "--sign", "-", str(app)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        print("  signed: ad-hoc ok")
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        err = getattr(exc, "stderr", b"") or b""
        print(f"  ! ad-hoc signing failed ({err.decode(errors='replace').strip() or exc})")
        print("    macOS may deny the app access to Dropbox/Documents without a prompt.")


def _refresh_finder(app: Path) -> None:
    """Register the bundle with LaunchServices and bump its mtime so Finder re-reads
    its icon instead of showing a cached generic one."""
    import os
    lsregister = ("/System/Library/Frameworks/CoreServices.framework/Frameworks/"
                  "LaunchServices.framework/Support/lsregister")
    try:
        subprocess.run([lsregister, "-f", str(app)], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass
    try:
        os.utime(app, None)
    except Exception:
        pass


def main() -> None:
    if sys.platform != "darwin":
        sys.exit("This builds a macOS .app bundle; run it on macOS.")
    ap = argparse.ArgumentParser(description="Build a CRAIC.app bundle")
    ap.add_argument("--dest", default=str(ROOT),
                    help="directory to write CRAIC.app into (default: project root)")
    args = ap.parse_args()
    _check_not_editable()                       # before ROOT joins sys.path
    sys.path.insert(0, str(ROOT))               # so `craic` imports even if not pip-installed
    app, have_icon, size = build(Path(args.dest))
    _report_engines()
    _codesign(app)
    _refresh_finder(app)
    print(f"\nBuilt {app}")
    if have_icon:
        print(f"  icon:  embedded ok  (Contents/Resources/craic.icns, {size // 1024} KB)")
        print("  If Finder still shows a generic icon it is caching it \u2014 run:  killall Finder")
    else:
        print("  icon:  NOT embedded \u2014 see the warning above (PySide6 render or iconutil).")
        print("         The app still runs; it just falls back to the default icon.")
    print("Double-click it in Finder, or drag it onto your Dock / into /Applications.")


if __name__ == "__main__":
    main()
