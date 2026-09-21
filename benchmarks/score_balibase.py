#!/usr/bin/env python3
"""Score BAliBASE test alignments with the *official* ``bali_score``.

``run_benchmark.py`` computes its own SP and TC over all columns. Those numbers
are internally consistent and fine for comparing conditions, but they are not
the numbers published BAliBASE tables report: the official score is restricted
to the reference's core blocks, which are the regions the structural alignment
is actually confident about. Quoting an all-column figure beside a published
core-block one is a silent apples-to-oranges comparison, so anything destined
for a paper goes through this script.

It reads the ``bali_manifest.tsv`` that ``run_benchmark.py --write-alignments``
leaves behind, converts each test alignment to the GCG/MSF layout the official
scorer insists on, runs the scorer, and collects the results.

    python benchmarks/score_balibase.py aligned_out/bali_manifest.tsv \\
        --bali-score /path/to/bali_score -o balibase_official.csv

Getting the scorer: BAliBASE 3 ships its source in ``bb3_release/bali_score_src``.
The bundled ``libexpat.a`` is a 2008 Alpha binary and will not link, so build
against the system expat instead::

    cd bb3_release/bali_score_src
    cc -O -U_FORTIFY_SOURCE -o bali_score readxml.c init.c util.c bali_score.c -lm -lexpat

(``-U_FORTIFY_SOURCE`` because the 2008 sources trip modern glibc's
``sprintf`` bounds checking and abort otherwise.) Verify it with a reference
against itself, which must give 1.000 / 1.000::

    ./bali_score ../RV11/BB11001.xml ../RV11/BB11001.msf
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile

#: The official reader truncates names at 30 characters and splits on the first
#: run of spaces, so a name must be short and must not contain a space.
MAX_NAME = 30
#: Residues per chunk in the written MSF. The reader takes one line at a time up
#: to its own 10k line limit, so wrapping is about staying well inside that.
WRAP = 50


def read_fasta(path: str):
    ids, rows, cur = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if cur:
                    rows.append("".join(cur))
                    cur = []
                ids.append(line[1:].split()[0] if len(line) > 1 else "")
            elif line:
                cur.append(line.strip())
    if cur:
        rows.append("".join(cur))
    return ids, rows


def write_msf(path: str, ids, rows) -> None:
    """The minimum the official reader accepts: a header, ``//``, then one
    block of ``name sequence`` lines per chunk, blocks separated by blank lines.
    """
    width = max((len(r) for r in rows), default=0)
    names = [re.sub(r"\s+", "_", i)[:MAX_NAME] for i in ids]
    pad = max((len(n) for n in names), default=1) + 2
    with open(path, "w") as fh:
        fh.write("PileUp\n\n")
        fh.write(f"   MSF: {width:>4}  Type: P    Check:  0   ..\n\n")
        for n, r in zip(names, rows):
            fh.write(f" Name: {n} oo  Len: {len(r):>4}  Check:  0  Weight:  1.00\n")
        fh.write("\n//\n\n")
        for start in range(0, width, WRAP):
            for n, r in zip(names, rows):
                fh.write(f"{n:<{pad}}{r[start:start + WRAP]}\n")
            fh.write("\n")



def read_manifest(path: str):
    with open(path) as fh:
        return [tuple(r) for r in list(csv.reader(fh, delimiter="\t"))[1:] if len(r) >= 3]


def from_directory(aligned_dir: str, reference_dir: str):
    """Pair up every written alignment with its reference, without a manifest.

    ``run_benchmark.py`` writes the manifest when it finishes, so a run that is
    still going — or that was interrupted — leaves an incomplete one. The file
    names carry everything needed (``FAMILY__engine.fasta``), so the pairing can
    be rebuilt from the directory instead.
    """
    if not reference_dir:
        raise SystemExit("scoring a directory needs --reference-dir")
    refs = {}
    for root, _dirs, names in os.walk(reference_dir):
        for n in names:
            if n.endswith(".xml"):
                refs[n[:-4]] = os.path.join(root, n)
    out = []
    for name in sorted(os.listdir(aligned_dir)):
        if not name.endswith(".fasta") or "__" not in name:
            continue
        family, engine = name[:-len(".fasta")].split("__", 1)
        if family in refs:
            out.append((os.path.join(aligned_dir, name), refs[family], engine))
    return out


_SP = re.compile(r"SP score\s*=\s*([0-9.]+)")
_TC = re.compile(r"TC score\s*=\s*([0-9.]+)")


def score_one(scorer: str, reference_xml: str, test_fasta: str):
    """(SP, TC) from the official scorer, or (None, None) if it refused."""
    ids, rows = read_fasta(test_fasta)
    if not ids:
        return None, None, "empty test alignment"
    with tempfile.NamedTemporaryFile("w", suffix=".msf", delete=False) as tmp:
        msf = tmp.name
    try:
        write_msf(msf, ids, rows)
        proc = subprocess.run([scorer, reference_xml, msf],
                              capture_output=True, text=True, timeout=300)
    except (subprocess.TimeoutExpired, OSError) as exc:
        return None, None, str(exc)
    finally:
        os.unlink(msf)
    out = proc.stdout + proc.stderr
    sp, tc = _SP.search(out), _TC.search(out)
    if not sp or not tc:
        return None, None, out.strip().splitlines()[-1] if out.strip() else "no score"
    return float(sp.group(1)), float(tc.group(1)), ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest",
                    help="bali_manifest.tsv from run_benchmark.py, or the "
                         "--write-alignments directory itself (with --reference-dir)")
    ap.add_argument("--reference-dir",
                    help="where the BAliBASE .xml references live; required when "
                         "scoring a directory rather than a manifest")
    ap.add_argument("--bali-score", default="bali_score",
                    help="path to the compiled official scorer")
    ap.add_argument("-o", "--out", default="balibase_official.csv")
    args = ap.parse_args(argv)

    rows = (from_directory(args.manifest, args.reference_dir)
            if os.path.isdir(args.manifest) else read_manifest(args.manifest))
    if not rows:
        print("nothing to score", file=sys.stderr)
        return 1
    print(f"scoring {len(rows)} alignments", flush=True)

    results, failures = [], []
    for n, (test, reference, aligner) in enumerate(rows, 1):
        family = os.path.basename(reference).rsplit(".", 1)[0]
        sp, tc, why = score_one(args.bali_score, reference, test)
        if sp is None:
            failures.append((family, aligner, why))
        else:
            results.append({"family": family, "set": family[:4], "aligner": aligner,
                            "sp": f"{sp:.4f}", "tc": f"{tc:.4f}"})
        if n % 100 == 0:
            print(f"  {n}/{len(rows)} scored", flush=True)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["family", "set", "aligner", "sp", "tc"])
        w.writeheader()
        w.writerows(results)
    print(f"wrote {len(results)} scores to {args.out}")

    if failures:
        print(f"\n{len(failures)} could not be scored:", file=sys.stderr)
        for family, aligner, why in failures[:20]:
            print(f"  {family} / {aligner}: {why}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
