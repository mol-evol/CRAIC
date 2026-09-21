#!/usr/bin/env python3
"""CRAIC alignment / reliability benchmark.

Two ground-truth sources, the same metrics for each:

* ``--sim``   a simulation sweep (divergence x indel rate x #taxa x replicates)
              with an exactly known true alignment and tree.
* ``--reference GLOB``  real reference alignments (BAliBASE, HOMSTRAD, PREFAB,
              ...); any format Biopython can read. Sequences are re-aligned and
              scored against the reference.

For every (condition, aligner) it reports:
  * alignment accuracy  : SP (sum-of-pairs) and TC (total-column) vs truth;
  * reliability quality : ROC AUC of CRAIC's per-column reliability predicting
                          which columns are actually mis-aligned;
  * masking gain        : mean column accuracy of the kept columns, against a
                          random mask removing the same number of columns and
                          against a gap-score mask of the same size. The bare
                          all-columns vs kept-columns contrast is not reported
                          as a result: it rises for any score correlated with
                          column difficulty and so proves nothing;
  * tree error          : Robinson-Foulds of the NJ tree to the true topology,
                          full alignment vs reliability-masked (sim mode only).

Aligners are auto-detected: the CRAIC built-in always runs; MAFFT / MUSCLE /
Clustal Omega / PRANK are included automatically if on the PATH.

Runs on the NumPy fallback, but is much faster with the compiled Rust core.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
import time
from typing import Dict, List, Optional

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                  # benchmarks/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root -> craic
import metrics as M            # noqa: E402
from craic import simulate as S                                    # noqa: E402

from craic.domain import Alphabet, Alignment                       # noqa: E402
from craic.ambiguity import reliability as rel_mod                 # noqa: E402
from craic.ambiguity import trimming                               # noqa: E402
from craic import engines, io                                      # noqa: E402
from craic import __version__ as craic_version                     # noqa: E402


# Families skipped, and (condition, item, aligner, error) for every engine failure.
# Both are reported at the end: a systematically failing aligner used to vanish
# from the results, leaving the means averaged over unequal, unreported subsets.
_SKIPPED: List[tuple] = []
_FAILURES: List[tuple] = []


def detect_alphabet(rows: List[str]) -> Alphabet:
    """Protein or nucleotide, read off the reference alignment itself."""
    from craic.domain import detect_alphabet as _detect
    return _detect([r.replace("-", "") for r in rows])


def _r(x, nd: int = 4):
    """Round for the CSV, writing an empty cell rather than 'nan'."""
    try:
        return round(float(x), nd) if float(x) == float(x) else ""
    except (TypeError, ValueError):
        return ""


def aligners():
    """CRAIC built-in plus any external engine found on the PATH."""
    return [e for e in engines.available_engines()]


def aligner_versions(engs) -> Dict[str, str]:
    """Version string for each engine, so a result set names the binaries that
    produced it. A benchmark whose competitor versions are unrecorded is not
    reproducible."""
    out = {}
    for eng in engs:
        binary = getattr(eng, "binary", "")
        if not binary:
            out[eng.label] = f"craic {craic_version}"
            continue
        for flag in ("--version", "-version", "-v"):
            try:
                proc = eng._run([binary, flag], timeout=20)
            except Exception:
                continue
            text = (proc.stdout or proc.stderr or "").strip().splitlines()
            if text:
                out[eng.label] = text[0][:120]
                break
        else:
            out[eng.label] = "unknown"
    return out


def _core_mask(raw_rows: List[str]) -> Optional[np.ndarray]:
    """Core-block mask from the case convention, or None if the file has none.

    Several reference-alignment distributions mark unreliable (non-core) regions
    by lower-casing them. Published BAliBASE SP/TC figures are computed over core
    blocks only, so uppercasing the data on the way in — as this script used to —
    silently converts the score into something that cannot be compared with any
    published table. Where the convention is present, we honour it; where it is
    not, we return None and the caller labels the score as full-column rather
    than pretending otherwise. Either way the official ``bali_score`` manifest
    remains the citable path.
    """
    has_lower = any(any(c.islower() for c in r) for r in raw_rows)
    has_upper = any(any(c.isupper() for c in r) for r in raw_rows)
    if not (has_lower and has_upper):
        return None
    width = len(raw_rows[0])
    # A column is core if every residue in it is upper case.
    return np.array([
        all((not r[c].isalpha()) or r[c].isupper() for r in raw_rows)
        for c in range(width)
    ], dtype=bool)


def _read_balibase_xml(path: str):
    """BAliBASE .xml -> (names, gapped rows, core mask or None).

    The reference alignment lives in the <seq-data> of each <sequence>
    (BAliBASE ships no .msf in recent releases).
    """
    import xml.etree.ElementTree as ET
    try:
        root = ET.parse(path).getroot()
    except Exception:
        return None
    names, raw = [], []
    for sq in root.findall(".//sequence"):
        nm = (sq.findtext("seq-name") or "").strip()
        data = "".join((sq.findtext("seq-data") or "").split()).replace(".", "-")
        if nm and data:
            names.append(nm)
            raw.append(data)
    if len(raw) >= 3 and len(set(len(r) for r in raw)) == 1:
        return names, [r.upper() for r in raw], _core_mask(raw)
    return None


def _read_reference(path: str):
    """Reference MSA -> (names, gapped rows, core mask or None).

    BAliBASE .xml or any Biopython-readable format.
    """
    if path.lower().endswith(".xml"):
        return _read_balibase_xml(path)
    from Bio import AlignIO
    fmt = io.guess_format(path)
    tried = []
    for f in (fmt, "fasta", "clustal", "msf", "stockholm", "phylip-relaxed"):
        try:
            aln = AlignIO.read(path, f)
        except Exception as exc:                       # wrong format guess, keep trying
            tried.append(f"{f}: {exc}")
            continue
        names = [r.id for r in aln]
        raw = [str(r.seq).replace(".", "-") for r in aln]
        if raw and len(set(len(x) for x in raw)) == 1:
            return names, [r.upper() for r in raw], _core_mask(raw)
        tried.append(f"{f}: parsed but rows are ragged")
    print(f"  ! could not read {path}; tried {'; '.join(tried[:3])}", flush=True)
    return None


def evaluate(true_rows, names, true_tree, aln: Alignment, threshold: float,
             do_perturbation: bool, core: Optional[np.ndarray] = None) -> dict:
    inf_rows = list(aln.rows)
    sp, _ = M.sp_score(true_rows, inf_rows, core)
    tc = M.tc_score(true_rows, inf_rows, core)
    rep = rel_mod.analyse(aln, do_perturbation=do_perturbation)
    rel = rep.col_combined if do_perturbation else rep.col_consistency
    frac = M.col_correct_fraction(inf_rows, true_rows)
    ok = ~np.isnan(rel) & ~np.isnan(frac)
    positive = frac[ok] >= 0.999
    a = M.auc(rel[ok], positive) if ok.sum() > 3 else float("nan")
    # The random-mask null and a gap-score control at matched column count: the
    # bare acc_all vs acc_retained contrast is not evidence on its own.
    mg = M.masking_gain(rel, frac, threshold,
                        controls={"gap": trimming.gap_score(aln)})
    row = {"sp": round(sp, 4), "tc": round(tc, 4),
           "score_scope": "core" if core is not None else "all-columns",
           "rel_auc": round(a, 4) if a == a else "",
           "pos_rate": round(float(positive.mean()), 4) if positive.size else "",
           "acc_all": round(mg["acc_all"], 4), "acc_retained": round(mg["acc_retained"], 4),
           "acc_retained_random": _r(mg["acc_retained_random"]),
           "acc_retained_gap": _r(mg.get("acc_retained_gap", float("nan"))),
           "mask_z": _r(mg["z"]),
           "n_masked": mg["n_masked"], "n_cols": mg["n_cols"]}
    if true_tree is not None and aln.n_seqs >= 4:
        tb = M.true_tree_bipartitions(true_tree, names)
        row["rf_full"] = round(M.rf_distance(tb, M.nj_tree_bipartitions(names, inf_rows)), 4)
        keep = rep.keep_mask(threshold)
        masked = rel_mod.apply_mask(aln, keep)
        if masked.length >= aln.n_seqs:
            row["rf_masked"] = round(M.rf_distance(tb, M.nj_tree_bipartitions(names, list(masked.rows))), 4)
        else:
            row["rf_masked"] = ""
    else:
        row["rf_full"] = row["rf_masked"] = ""
    return row


def run_sim(args, writer, summary, done):
    for bmax in args.divergence:
        for ir in args.indel:
            for taxa in args.taxa:
                for rep in range(args.reps):
                    cond, item = f"div{bmax}_indel{ir}_taxa{taxa}", f"rep{rep}"
                    pending = [e for e in aligners() if (cond, item, e.label) not in done]
                    if not pending:
                        continue
                    d = S.simulate(taxa=taxa, root_len=args.length, seed=rep,
                                   indel_rate=ir, bmax=bmax,
                                   rate_alpha=args.rate_alpha, indel_zipf=args.indel_zipf)
                    for eng in pending:
                        # Time the alignment and nothing else. Previously the
                        # clock ran across evaluate() too — reliability analysis,
                        # NJ tree, RF distance — which is not aligner runtime and
                        # scales with alignment length, so it differed
                        # systematically between aligners.
                        t0 = time.perf_counter()
                        try:
                            aln = eng.align(d["seqs"], Alphabet.DNA)
                        except Exception as exc:
                            print(f"  ! {eng.label} failed: {exc}", flush=True)
                            _FAILURES.append((cond, item, eng.label, str(exc)[:200]))
                            continue
                        secs = round(time.perf_counter() - t0, 2)
                        row = evaluate(d["true_rows"], d["names"], d["tree"], aln,
                                       args.threshold, args.perturbation)
                        row.update(mode="sim", condition=cond, item=item, aligner=eng.label,
                                   n_seqs=aln.n_seqs, secs=secs)
                        writer.writerow(row); summary.append(row)
                        print(f"  [{cond}/{item}] {eng.label}: "
                              f"SP={row['sp']} TC={row['tc']} relAUC={row['rel_auc']} "
                              f"RF {row['rf_full']}->{row['rf_masked']} ({row['secs']}s)", flush=True)


def run_reference(args, writer, summary, done):
    files = sorted(f for pat in args.reference for f in glob.glob(pat, recursive=True))
    if not files:
        print(f"no reference files matched: {args.reference}"); return
    if args.write_alignments:
        os.makedirs(args.write_alignments, exist_ok=True)
    for path in files:
        base = os.path.splitext(os.path.basename(path))[0]
        cond = os.path.basename(path)
        pending = [e for e in aligners() if (cond, "ref", e.label) not in done]
        if not pending:
            continue
        ref = _read_reference(path)
        if not ref:
            _SKIPPED.append((cond, "unreadable")); continue
        names, true_rows, core = ref
        seqs = [(n, r.replace("-", "")) for n, r in zip(names, true_rows)]
        if len(seqs) < 3 or min(len(s) for _, s in seqs) < 5:
            _SKIPPED.append((cond, "fewer than 3 sequences, or a sequence under 5 residues"))
            continue
        # Detect the alphabet from the data. Defaulting to DNA meant that running
        # a protein benchmark without --protein scored every family under a DNA
        # emission model, silently producing meaningless numbers.
        alph = detect_alphabet(true_rows) if args.protein is None else (
            Alphabet.PROTEIN if args.protein else Alphabet.DNA)
        maxlen = max(len(s) for _, s in seqs)
        gb = (len(seqs) * (len(seqs) - 1) // 2) * (maxlen ** 2) * 8 / 1e9
        if gb > args.max_mem_gb:
            # NB this bound comes from CRAIC's dense-posterior cost; the external
            # aligners have no such cost. Skipping the family for every engine
            # keeps the means comparable, but the count must be reported.
            print(f"  ~ skip {cond}: {len(seqs)} seqs x {maxlen} res (~{gb:.0f} GB posteriors) "
                  f"> --max-mem-gb {args.max_mem_gb}", flush=True)
            _SKIPPED.append((cond, f"~{gb:.0f} GB posteriors > --max-mem-gb {args.max_mem_gb}"))
            continue
        for eng in pending:
            t0 = time.perf_counter()
            try:
                aln = eng.align(seqs, alph)
            except Exception as exc:
                print(f"  ! {eng.label} failed on {cond}: {exc}", flush=True)
                _FAILURES.append((cond, "ref", eng.label, str(exc)[:200]))
                continue
            secs = round(time.perf_counter() - t0, 2)
            row = evaluate(true_rows, names, None, aln, args.threshold, args.perturbation, core)
            row.update(mode="reference", condition=cond, item="ref",
                       aligner=eng.label, n_seqs=aln.n_seqs, secs=secs)
            writer.writerow(row); summary.append(row)
            print(f"  [{cond}] {eng.label}: SP={row['sp']} TC={row['tc']} "
                  f"relAUC={row['rel_auc']}", flush=True)
            if args.write_alignments:
                key = re.sub(r"[^A-Za-z0-9]+", "_", eng.key)
                outp = os.path.join(args.write_alignments, base + "__" + key + ".fasta")
                io.write_fasta(outp, list(aln.ids), list(aln.rows), wrap=0)
    if args.write_alignments:
        _rebuild_manifest(args.write_alignments, files)


def _rebuild_manifest(outdir, files):
    """(Re)build the bali_score manifest from every alignment present on disk, so
    it is complete even across resumed runs."""
    man = []
    engs = aligners()
    for path in files:
        base = os.path.splitext(os.path.basename(path))[0]
        for eng in engs:
            key = re.sub(r"[^A-Za-z0-9]+", "_", eng.key)
            fp = os.path.join(outdir, base + "__" + key + ".fasta")
            if os.path.exists(fp):
                man.append((fp, os.path.abspath(path), eng.label))
    if man:
        _write_manifest(outdir, man)


def _write_manifest(outdir, manifest):
    """Write a tab-separated manifest of (test alignment, reference xml, aligner)
    for official scoring with bali-score. Paths may contain spaces, so the
    manifest is tab-delimited and best consumed by the Python runner in the
    benchmarks README."""
    tab, nl = chr(9), chr(10)
    mpath = os.path.join(outdir, "bali_manifest.tsv")
    with open(mpath, "w") as fh:
        fh.write(tab.join(["test_fasta", "reference_xml", "aligner"]) + nl)
        for t, r, a in manifest:
            fh.write(tab.join([t, r, a]) + nl)
    print(nl + "wrote " + str(len(manifest)) + " test alignments + bali_manifest.tsv to " + outdir + "/")
    print("score with the official bali-score - see benchmarks/README.md (Official BAliBASE scoring)")


def _report_exclusions():
    """Say what is missing from the means, and why."""
    if _SKIPPED:
        print(f"\n{len(_SKIPPED)} families skipped (excluded from every aligner's mean):")
        reasons = {}
        for _, why in _SKIPPED:
            reasons[why] = reasons.get(why, 0) + 1
        for why, k in sorted(reasons.items(), key=lambda kv: -kv[1]):
            print(f"    {k:4d}  {why}")
    if _FAILURES:
        print(f"\n{len(_FAILURES)} engine failures (these rows are absent, so the "
              f"means below are over unequal subsets):")
        per = {}
        for _, _, label, err in _FAILURES:
            per.setdefault(label, []).append(err)
        for label, errs in per.items():
            print(f"    {label:32s} {len(errs):4d}  e.g. {errs[0]}")


def print_summary(summary):
    if not summary:
        return
    print("\n=== summary (means by aligner) ===")
    by = {}
    for r in summary:
        by.setdefault(r["aligner"], []).append(r)
    for aligner, rows in by.items():
        def m(k):
            vals = [r[k] for r in rows if isinstance(r.get(k), (int, float))]
            return f"{np.mean(vals):.3f}" if vals else "-"
        print(f"  {aligner:32s} SP={m('sp')} TC={m('tc')} relAUC={m('rel_auc')} "
              f"kept={m('acc_retained')} vs random={m('acc_retained_random')} "
              f"vs gap={m('acc_retained_gap')} z={m('mask_z')} "
              f"RF_full={m('rf_full')} RF_masked={m('rf_masked')}  (n={len(rows)})")
    scopes = {r.get("score_scope") for r in summary}
    if scopes:
        print(f"  SP/TC scope: {', '.join(sorted(str(x) for x in scopes))}"
              " (all-columns figures are NOT comparable with published"
              " core-block BAliBASE tables; use the bali_score manifest for those)")
    print("""
Reading this:
  relAUC > 0.5            the reliability score carries information about which
                          columns are mis-aligned. This is the headline number.
  kept vs random          masking is only doing work if the columns it keeps are
                          more accurate than the same NUMBER of columns chosen at
                          random. 'kept > acc_all' on its own is not evidence:
                          removing hard columns raises the mean whatever rule
                          picks them.
  kept vs gap             and it has to beat the trivial gap-fraction rule to be
                          worth the compute.
  z                       the same comparison in units of the random-mask spread.
  RF_masked <= RF_full    masking improves the downstream tree (sim mode only).""")


def _engine_by_key(key: str):
    """The engine to build alignments with in --masking mode."""
    for eng in aligners():
        if eng.key == key:
            return eng
    have = ", ".join(e.key for e in aligners())
    raise SystemExit(f"--mask-engine {key!r} is not available; have: {have}")


def run_masking(args):
    """Compare masking methods on each reference/simulated family: which best
    removes the truly mis-aligned columns."""
    order = ["CRAIC reliability", "gap (MSA_trimmer)", "similarity (trimAl)",
             "trimAl gappyout", "trimAl strict", "Gblocks"]
    # Which engine produces the alignments being masked. Scoring CRAIC's
    # reliability only on CRAIC's own alignments gives it home-field advantage:
    # the score derives from the same posterior model that built the alignment,
    # while trimAl and Gblocks see a foreign one. --mask-engine lets the whole
    # comparison be repeated on a MAFFT or MUSCLE alignment, which is the
    # version a referee will ask for.
    engine = _engine_by_key(args.mask_engine)
    fields = ["family", "engine", "method", "auc", "kept", "acc_all", "acc_kept",
              "acc_kept_random", "precision", "recall"]
    done = set()
    mode = "w"
    if os.path.exists(args.out) and os.path.getsize(args.out) > 0 and not args.fresh:
        with open(args.out) as fh:
            for r in csv.DictReader(fh):
                done.add(r["family"])
        mode = "a"
        print(f"resuming: {len(done)} families already in {args.out}")

    def families():
        if args.reference:
            for path in sorted(f for pat in args.reference for f in glob.glob(pat, recursive=True)):
                ref = _read_reference(path)
                if not ref:
                    continue
                names, true_rows, _core = ref
                yield os.path.basename(path), [(n, r.replace("-", "")) for n, r in zip(names, true_rows)], true_rows
        else:
            for seed in range(args.reps):
                d = S.simulate(taxa=args.taxa[0], root_len=args.length, seed=seed,
                               indel_rate=args.indel[0], bmax=args.divergence[0],
                               rate_alpha=args.rate_alpha, indel_zipf=args.indel_zipf)
                yield f"sim{seed}", d["seqs"], d["true_rows"]

    agg = {m: [] for m in order}
    with open(args.out, mode, newline="", buffering=1) as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if mode == "w":
            w.writeheader()
        for fam, seqs, true_rows in families():
            if fam in done or len(seqs) < 3 or min(len(s) for _, s in seqs) < 5:
                continue
            maxlen = max(len(s) for _, s in seqs)
            if (len(seqs) * (len(seqs) - 1) // 2) * (maxlen ** 2) * 8 / 1e9 > args.max_mem_gb:
                print(f"  ~ skip {fam} (too big)", flush=True); continue
            alph = detect_alphabet(true_rows) if args.protein is None else (
                Alphabet.PROTEIN if args.protein else Alphabet.DNA)
            try:
                aln = engine.align(seqs, alph)
            except Exception as exc:
                print(f"  ! {fam} align failed: {exc}", flush=True); continue
            res = M.compare_masking(aln, true_rows)
            for m in order:
                v = res[m]
                w.writerow({"family": fam, "engine": engine.label, "method": m,
                            **{k: _r(v[k]) for k in
                               ("auc", "kept", "acc_all", "acc_kept",
                                "acc_kept_random", "precision", "recall")}})
                agg[m].append(v)
            best = max(order, key=lambda m: res[m]["auc"] if res[m]["auc"] == res[m]["auc"] else -1)
            print(f"  [{fam}] best column-AUC: {best} ({res[best]['auc']:.3f})", flush=True)

    print("\n=== masking methods (means; AUC = how well it ranks correct columns) ===")
    print(f"{'method':22s} {'AUC':>6} {'kept':>6} {'acc_kept':>9} {'vs random':>10} "
          f"{'precision':>10} {'recall':>7}")
    for m in order:
        vs = agg[m]
        def mn(f):
            x = [v[f] for v in vs if v[f] == v[f]]
            return float(np.mean(x)) if x else float("nan")
        au = mn("auc")
        print(f"{m:22s} {'   -  ' if au != au else f'{au:6.3f}'} {mn('kept'):6.2f} "
              f"{mn('acc_kept'):9.3f} {mn('acc_kept_random'):10.3f} "
              f"{mn('precision'):10.3f} {mn('recall'):7.3f}")
    print("\nHigher AUC = better at telling correct columns from wrong ones. Every method,"
          "\nincluding the rule-based ones, is scored on the same threshold-free footing,"
          "\nand 'vs random' is what the same number of randomly-kept columns would score."
          f"\nAlignments were built with: {engine.label}.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sim", action="store_true", help="run the simulation sweep")
    ap.add_argument("--reference", nargs="+", metavar="GLOB", help="reference alignment file(s)")
    ap.add_argument("--taxa", type=int, nargs="+", default=[8])
    ap.add_argument("--length", type=int, default=200, help="root sequence length (sim)")
    ap.add_argument("--reps", type=int, default=5, help="replicates per condition (sim)")
    ap.add_argument("--divergence", type=float, nargs="+", default=[0.22, 0.35],
                    help="branch-length upper bounds (higher = more divergent)")
    ap.add_argument("--indel", type=float, nargs="+", default=[1.5, 3.0], help="indel rates (sim)")
    ap.add_argument("--rate-alpha", type=float, default=0.0,
                    help="gamma shape for among-site rate variation (sim); 0 = uniform rates. "
                         "The uniform default makes CRAIC's pair-HMM nearly well-specified on "
                         "this data, so report a run with rate variation too.")
    ap.add_argument("--indel-zipf", type=float, default=0.0,
                    help="Zipf exponent for indel lengths (sim); 0 = geometric. >1 gives the "
                         "heavy tail real indel distributions have.")
    ap.add_argument("--threshold", type=float, default=0.5, help="reliability masking threshold")
    ap.add_argument("--perturbation", action="store_true",
                    help="use combined consistency+perturbation reliability (slower)")
    ap.add_argument("--protein", dest="protein", action="store_true", default=None,
                    help="force protein; default is to detect the alphabet from the data")
    ap.add_argument("--dna", dest="protein", action="store_false",
                    help="force nucleotide instead of detecting it")
    ap.add_argument("--write-alignments", metavar="DIR",
                    help="write each engine alignment (FASTA, reference names) for official bali_score")
    ap.add_argument("--max-mem-gb", type=float, default=1.0,
                    help="skip a family if its dense pairwise posteriors would exceed this many GB "
                         "(guards the built-in aligner against out-of-memory on huge alignments)")
    ap.add_argument("--engine-timeout", type=float, default=600.0,
                    help="per external-aligner timeout in seconds (skips hangs); 0 = none")
    ap.add_argument("--effort", default="max", choices=["min", "med", "max"],
                    help="compute level for CRAIC's built-in engine. Default 'max', which "
                         "includes iterative refinement: benchmarking at 'med' measures the "
                         "engine with refinement switched off, which understates the method "
                         "it implements rather than just this implementation of it.")
    ap.add_argument("--fresh", action="store_true", help="overwrite --out instead of resuming")
    ap.add_argument("--masking", action="store_true",
                    help="compare column-masking methods (CRAIC reliability vs trimAl/Gblocks/gap) "
                         "instead of aligners")
    ap.add_argument("--mask-engine", default="builtin",
                    help="engine that builds the alignments to be masked in --masking mode "
                         "(builtin, mafft, muscle, clustalo, prank). Repeat the comparison on a "
                         "foreign alignment: scoring CRAIC's reliability only on CRAIC's own "
                         "alignments flatters it.")
    ap.add_argument("--out", default="benchmark_results.csv")
    ap.add_argument("--quick", action="store_true", help="tiny preset for a smoke test")
    args = ap.parse_args()

    if args.quick:
        args.taxa, args.length, args.reps = [6], 150, 2
        args.divergence, args.indel = [0.25], [2.0]
    if not args.sim and not args.reference:
        args.sim = True

    engines.ExternalAligner.timeout = args.engine_timeout or None
    engines.BuiltinProgressive.consistency_mem_gb = args.max_mem_gb
    engines.BuiltinProgressive.benchmark_effort = args.effort

    if args.masking:
        run_masking(args)
        return

    done = set()
    mode = "w"
    if os.path.exists(args.out) and os.path.getsize(args.out) > 0 and not args.fresh:
        with open(args.out) as fh:
            for r in csv.DictReader(fh):
                done.add((r.get("condition"), r.get("item"), r.get("aligner")))
        mode = "a"
        print(f"resuming: {len(done)} (item, aligner) results already in {args.out}")

    engs = aligners()
    print("aligners:", ", ".join(e.label for e in engs))
    for label, ver in aligner_versions(engs).items():
        print(f"    {label:32s} {ver}")
    fields = ["mode", "condition", "item", "aligner", "n_seqs", "n_cols",
              "sp", "tc", "score_scope", "rel_auc", "pos_rate",
              "acc_all", "acc_retained", "acc_retained_random", "acc_retained_gap",
              "mask_z", "n_masked", "rf_full", "rf_masked", "secs"]
    summary: List[dict] = []
    with open(args.out, mode, newline="", buffering=1) as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if mode == "w":
            writer.writeheader()
        if args.sim:
            print("\n--- simulation sweep ---")
            run_sim(args, writer, summary, done)
        if args.reference:
            print("\n--- reference benchmark ---")
            run_reference(args, writer, summary, done)
    print_summary(summary)
    _report_exclusions()
    print(f"\nwrote {len(summary)} new rows to {args.out}")


if __name__ == "__main__":
    main()
