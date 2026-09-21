"""Command line.

``craic`` and ``craic FILE`` launch the workbench, exactly as they always have.
Everything else is a subcommand that runs headless:

    craic align seqs.fasta -o aln.fasta
    craic score aln.fasta -o scores.tsv
    craic mask  aln.fasta -o masked.fasta --threshold 0.5
    craic trim  aln.fasta -o trimmed.fasta --method gappyout
    craic simulate -o truth.fasta --taxa 12 --divergence 0.35
    craic session work.fasta.craic.json            # what is in a saved session

The headless path never imports PySide6. That is the point of it: a reliability
score you can only obtain by opening a window is a score that cannot go in a
pipeline, on a cluster node, or in a reproducible workflow — which is most of
where alignment curation actually happens.
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional, Sequence

_LAUNCH_HELP = """\
Run with no arguments to open the workbench, or with a single alignment file to
open it with that file loaded. The subcommands below all run headless.
"""


# --------------------------------------------------------------------------- #
# Shared helpers
# --------------------------------------------------------------------------- #

def _load(path: str, alphabet: Optional[str] = None):
    """Read an alignment, or exit with a message a user can act on."""
    from .domain import Alphabet
    from . import io

    alph = {"dna": Alphabet.DNA, "rna": Alphabet.RNA,
            "protein": Alphabet.PROTEIN}.get((alphabet or "").lower())
    try:
        aln = io.load_alignment(path, alphabet=alph)
    except Exception as exc:
        raise SystemExit(f"craic: cannot read {path}: {exc}")
    if aln.n_seqs == 0:
        raise SystemExit(f"craic: {path} contains no sequences")
    return aln


def _write(aln, out: Optional[str], fmt: Optional[str] = None) -> None:
    from . import io

    if out in (None, "-"):
        sys.stdout.write(io.dump_alignment(aln, fmt or "fasta", wrap=60))
        return
    io.save_alignment(out, aln, fmt)
    print(f"wrote {out}", file=sys.stderr)


def _table(rows: Sequence[Sequence[object]], header: Sequence[str],
           out: Optional[str]) -> None:
    """Write a tab-separated table to a file or stdout."""
    lines = ["\t".join(header)]
    lines += ["\t".join("" if v is None else str(v) for v in r) for r in rows]
    text = "\n".join(lines) + "\n"
    if out in (None, "-"):
        sys.stdout.write(text)
    else:
        with open(out, "w", encoding="utf-8") as fh:
            fh.write(text)
        print(f"wrote {out}", file=sys.stderr)


def _fmt(x, nd: int = 4) -> str:
    try:
        f = float(x)
    except (TypeError, ValueError):
        return ""
    return "" if f != f else f"{f:.{nd}f}"


def _engine(key: str):
    from . import engines

    for eng in engines.available_engines():
        if eng.key == key:
            return eng
    have = ", ".join(e.key for e in engines.available_engines())
    raise SystemExit(f"craic: engine {key!r} is not available (have: {have})")


# --------------------------------------------------------------------------- #
# Subcommands
# --------------------------------------------------------------------------- #

def cmd_align(args) -> int:
    from .domain import Alphabet, detect_alphabet
    from . import io

    records = io.read_records(args.input)
    if not records:
        raise SystemExit(f"craic: {args.input} contains no sequences")
    alph = {"dna": Alphabet.DNA, "rna": Alphabet.RNA,
            "protein": Alphabet.PROTEIN}.get((args.alphabet or "").lower())
    if alph is None:
        alph = detect_alphabet([s for _, s in records])
    eng = _engine(args.engine)
    if args.codon:
        from .engines import CodonAware
        eng = CodonAware(eng, table=args.code)
    opts = {} if args.engine != "builtin" else {"effort": args.effort}
    aln = eng.align(records, alph, **opts)
    _warn_internal_stops(aln)
    _write(aln, args.out, args.format)
    return 0


def _warn_internal_stops(aln) -> None:
    """Say so when a coding alignment stops in the middle of a sequence.

    Almost always one of three things: the wrong genetic code, the wrong reading
    frame, or a pseudogene. CRAIC cannot tell which, but silence is the worst of
    the options — mitochondrial and Mollicute genes read TGA as tryptophan, and
    under the standard code every one of those becomes a stop that carries no
    homology signal at all.
    """
    from .domain import code_name

    if not getattr(aln, "coding", None):
        return
    stops = aln.internal_stops()
    if not stops:
        return
    rows = sorted({r for r, _ in stops})
    table = aln.coding.table
    print(f"craic: {len(stops)} internal stop codon(s) in {len(rows)} of "
          f"{aln.n_seqs} sequences, reading the {code_name(table)} code "
          f"(table {table}).", file=sys.stderr)
    if table == 1:
        print("craic: if these are mitochondrial or Mollicute sequences, try "
              "--code 2 / --code 4 (`craic codes` lists them); otherwise check "
              "the reading frame.", file=sys.stderr)
    return


def cmd_codes(args) -> int:
    """The genetic-code tables, so that --code N does not need a web search."""
    from .domain import genetic_codes

    for tid, name in genetic_codes():
        print(f"{tid:>3}  {name}")
    return 0


def cmd_score(args) -> int:
    """Per-column reliability, as a table."""
    import numpy as np

    from .ambiguity import reliability as rel_mod
    from .ambiguity import trimming

    aln = _load(args.input, args.alphabet)
    rep = rel_mod.analyse(aln, do_perturbation=not args.fast,
                          n_replicates=args.replicates)
    gap = trimming.gap_score(aln)
    sim = trimming.similarity_score(aln)

    reference = None
    if args.reference:
        from .evaluate import compare_to_reference
        from .domain import IdMismatch
        try:
            reference = compare_to_reference(aln, _load(args.reference, args.alphabet))
        except IdMismatch as exc:
            raise SystemExit(f"craic: {exc}")

    header = ["column", "consistency", "perturbation", "combined",
              "gap_score", "similarity"]
    if reference is not None:
        header.append("correct_fraction")
    rows = []
    for c in range(aln.length):
        row = [c + 1, _fmt(rep.col_consistency[c]), _fmt(rep.col_perturbation[c]),
               _fmt(rep.col_combined[c]), _fmt(gap[c]), _fmt(sim[c])]
        if reference is not None:
            row.append(_fmt(reference.col_correct[c]))
        rows.append(row)
    _table(rows, header, args.out)

    if rep.n_replicates_failed:
        print(f"craic: {rep.n_replicates_failed} perturbation replicates failed",
              file=sys.stderr)
    if reference is not None:
        a = reference.reliability_auc(rep.col_combined)
        print(f"SP={reference.sp:.4f} TC={reference.tc:.4f} "
              f"precision={reference.precision:.4f} reliability_AUC={_fmt(a)}",
              file=sys.stderr)
    else:
        scored = int(np.sum(~np.isnan(rep.col_combined)))
        print(f"{scored}/{aln.length} columns scored", file=sys.stderr)
    return 0


def cmd_mask(args) -> int:
    from .ambiguity import reliability as rel_mod

    aln = _load(args.input, args.alphabet)
    rep = rel_mod.analyse(aln, do_perturbation=not args.fast,
                          n_replicates=args.replicates)
    try:
        keep = rep.keep_mask(args.threshold, which=args.which, unscored=args.unscored)
    except ValueError as exc:
        raise SystemExit(f"craic: {exc}")
    masked = rel_mod.apply_mask(aln, keep)
    if masked.length == 0:
        raise SystemExit(
            f"craic: a threshold of {args.threshold} removes every column; "
            "nothing was written")
    _write(masked, args.out, args.format)
    print(f"kept {int(keep.sum())}/{aln.length} columns "
          f"({args.which} >= {args.threshold})", file=sys.stderr)
    return 0


_TRIMMERS = {
    "gappyout": "trimAl gappyout: keep the low-gap mode of the gap distribution",
    "strict": "trimAl strict: gap and similarity thresholds, then drop short blocks",
    "gblocks": "Gblocks-style conserved blocks",
}


def cmd_trim(args) -> int:
    import numpy as np

    from .ambiguity import reliability as rel_mod
    from .ambiguity import trimming

    aln = _load(args.input, args.alphabet)
    fn = {"gappyout": trimming.gappyout_mask,
          "strict": trimming.trimal_strict_mask,
          "gblocks": trimming.gblocks_mask}[args.method]
    keep = np.asarray(fn(aln), dtype=bool)
    trimmed = rel_mod.apply_mask(aln, keep)
    if trimmed.length == 0:
        raise SystemExit(f"craic: {args.method} removed every column; nothing was written")
    _write(trimmed, args.out, args.format)
    print(f"kept {int(keep.sum())}/{aln.length} columns ({args.method})", file=sys.stderr)
    return 0


def cmd_simulate(args) -> int:
    from .domain import Alignment, Alphabet
    from . import io, simulate

    d = simulate.simulate(taxa=args.taxa, root_len=args.length, seed=args.seed,
                          indel_rate=args.indel, bmax=args.divergence,
                          rate_alpha=args.rate_alpha, indel_zipf=args.indel_zipf)
    truth = Alignment(list(d["names"]), list(d["true_rows"]), Alphabet.DNA)
    _write(truth, args.out, args.format)
    if args.unaligned:
        io.write_fasta(args.unaligned, [n for n, _ in d["seqs"]],
                       [s for _, s in d["seqs"]], wrap=60)
        print(f"wrote {args.unaligned}", file=sys.stderr)
    print(f"{args.taxa} sequences, true alignment {truth.length} columns",
          file=sys.stderr)
    return 0


def cmd_session(args) -> int:
    """Describe a saved session, and optionally export its alignment.

    A session is readable JSON by design, but "which alignment is in here and
    how did it get that way" is the question actually being asked, so it gets an
    answer without a text editor.
    """
    from .session import Session, SessionError

    try:
        sess = Session.load(args.input)
    except SessionError as exc:
        raise SystemExit(f"craic: {exc}")

    print(sess.describe())
    if sess.source_path:
        print(f"alignment from: {sess.source_path}")
    if sess.craic_version:
        import time as _time
        when = _time.strftime("%Y-%m-%d %H:%M", _time.localtime(sess.saved_at))
        print(f"saved by CRAIC {sess.craic_version} at {when}"
              + (" (autosave)" if sess.autosaved else ""))
    if sess.groups:
        print(f"groups: {', '.join(sorted(sess.groups))}")
    if sess.annotations:
        print(f"annotations: {len(sess.annotations)}")
    history = (sess.alignment.meta or {}).get("history") if sess.alignment else None
    for i, step in enumerate(history or [], 1):
        print(f"  {i}. {step}")
    if args.out and sess.alignment is not None:
        _write(sess.alignment, args.out, args.format)
    return 0


# --------------------------------------------------------------------------- #
# Argument parsing
# --------------------------------------------------------------------------- #

_SUBCOMMANDS = ("align", "score", "mask", "trim", "simulate", "session", "codes")


def _add_common(p, alphabet: bool = True) -> None:
    p.add_argument("-o", "--out", help="output file ('-' or omitted: stdout)")
    if alphabet:
        p.add_argument("--alphabet", choices=["dna", "rna", "protein"],
                       help="force the alphabet instead of detecting it from the data")


def _add_analysis(p) -> None:
    """Options that control how the reliability analysis is run."""
    p.add_argument("--fast", action="store_true",
                   help="consistency only; skip the perturbation ensemble")
    p.add_argument("--replicates", type=int, default=16,
                   help="perturbation replicates (default: 16)")


def _add_masking(p) -> None:
    """Options that control which columns a threshold keeps."""
    _add_analysis(p)
    p.add_argument("--threshold", type=float, default=0.5,
                   help="keep columns scoring at least this (default: 0.5)")
    p.add_argument("--which", choices=["combined", "consistency", "perturbation"],
                   default="combined", help="which reliability signal to threshold")
    p.add_argument("--unscored", choices=["drop", "keep"], default="drop",
                   help="what to do with columns that could not be scored at all "
                        "(default: drop — a column that could not be assessed is "
                        "not evidence)")


def build_parser() -> argparse.ArgumentParser:
    from . import __version__

    ap = argparse.ArgumentParser(
        prog="craic",
        description="CRAIC — a multiple sequence alignment workbench built for "
                    "ambiguous regions.",
        epilog=_LAUNCH_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"CRAIC {__version__}")
    ap.add_argument("--core", action="store_true",
                    help="report which acceleration core is in use, and exit")
    sub = ap.add_subparsers(dest="command")

    p = sub.add_parser("align", help="align sequences with any available engine")
    p.add_argument("input", help="unaligned sequences")
    _add_common(p)
    p.add_argument("--engine", default="builtin",
                   help="builtin, mafft, muscle, clustalo or prank (default: builtin)")
    p.add_argument("--effort", choices=["min", "med", "max"], default="med",
                   help="built-in engine compute level (default: med)")
    p.add_argument("--codon", action="store_true",
                   help="align coding nucleotides in amino-acid space and back-translate")
    p.add_argument("--code", type=int, default=1, metavar="N",
                   help="NCBI genetic-code table for --codon (default: 1, the "
                        "standard code; 2 = vertebrate mitochondrial, 4 = mould "
                        "mitochondrial / Mycoplasma / Spiroplasma, 11 = bacterial). "
                        "`craic codes` lists them all")
    p.add_argument("--format", help="output format (default: from the extension)")
    p.set_defaults(func=cmd_align)

    p = sub.add_parser("codes", help="list the NCBI genetic-code tables")
    p.set_defaults(func=cmd_codes)

    p = sub.add_parser("score", help="per-column reliability as a table")
    p.add_argument("input", help="an alignment")
    _add_common(p)
    _add_analysis(p)
    p.add_argument("--reference",
                   help="a trusted alignment of the same sequences; adds per-column "
                        "correctness and reports SP, TC and the reliability AUC")
    p.set_defaults(func=cmd_score)

    p = sub.add_parser("mask", help="drop columns below a reliability threshold")
    p.add_argument("input", help="an alignment")
    _add_common(p)
    _add_masking(p)
    p.add_argument("--format", help="output format (default: from the extension)")
    p.set_defaults(func=cmd_mask)

    p = sub.add_parser("trim", help="rule-based trimming (trimAl / Gblocks style)")
    p.add_argument("input", help="an alignment")
    _add_common(p)
    p.add_argument("--method", choices=list(_TRIMMERS), default="gappyout",
                   help="; ".join(f"{k}: {v}" for k, v in _TRIMMERS.items()))
    p.add_argument("--format", help="output format (default: from the extension)")
    p.set_defaults(func=cmd_trim)

    p = sub.add_parser("simulate", help="generate sequences with a known true alignment")
    _add_common(p, alphabet=False)
    p.add_argument("--taxa", type=int, default=8)
    p.add_argument("--length", type=int, default=200, help="root sequence length")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--divergence", type=float, default=0.28,
                   help="branch-length upper bound; higher is more divergent")
    p.add_argument("--indel", type=float, default=2.0, help="indel rate")
    p.add_argument("--rate-alpha", dest="rate_alpha", type=float, default=0.0,
                   help="gamma shape for among-site rate variation; 0 = uniform rates")
    p.add_argument("--indel-zipf", dest="indel_zipf", type=float, default=0.0,
                   help="Zipf exponent for indel lengths; 0 = geometric")
    p.add_argument("--unaligned", help="also write the unaligned sequences here")
    p.add_argument("--format", help="output format (default: from the extension)")
    p.set_defaults(func=cmd_simulate)

    p = sub.add_parser("session", help="describe a saved session, or export its alignment")
    p.add_argument("input", help="a .craic.json session document")
    p.add_argument("-o", "--out", help="also write the session's alignment here")
    p.add_argument("--format", help="output format (default: from the extension)")
    p.set_defaults(func=cmd_session)
    return ap


def main(argv: Optional[List[str]] = None) -> None:
    """Dispatch: a subcommand (or a help/version flag) runs headless; anything
    else — no arguments at all, or a file to open — launches the workbench, which
    is how ``craic`` and ``craic my.fasta`` have always behaved."""
    argv = list(sys.argv[1:] if argv is None else argv)

    if argv and argv[0] == "--core":
        from . import accel
        print(accel.backend())
        return

    headless = bool(argv) and (argv[0] in _SUBCOMMANDS
                               or argv[0] in ("-h", "--help", "--version"))
    if not headless:
        _launch_gui(argv)
        return

    args = build_parser().parse_args(argv)
    raise SystemExit(args.func(args))


def _launch_gui(argv: List[str]) -> None:
    try:
        from .gui.app import main as gui_main
    except ImportError as exc:                 # pragma: no cover - needs PySide6 absent
        raise SystemExit(
            f"craic: the workbench needs PySide6, which is not installed ({exc}).\n"
            "Install it with 'pip install PySide6', or use the headless "
            "subcommands — run 'craic --help' to see them.")
    sys.argv = [sys.argv[0]] + argv
    gui_main()


if __name__ == "__main__":
    main()
