"""Alignment engines behind one interface (Dependency-Inversion boundary).

The GUI and the ambiguity tools depend only on ``AlignerEngine``; concrete
engines (the built-in progressive aligner, plus subprocess wrappers for MAFFT /
MUSCLE / Clustal Omega / PRANK) are interchangeable. A ``CodonAware`` decorator
turns any engine into a reading-frame-preserving codon aligner via
translate -> align amino acids -> back-translate.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from . import progressive
from .domain import Alignment, Alphabet, CodingSpec, translate_codon

Record = Tuple[str, str]


@dataclass
class Param:
    """A tunable parameter an engine exposes to the GUI settings dialog."""

    key: str
    label: str
    kind: str                       # "float" | "int" | "choice"
    default: object
    choices: Optional[List[str]] = None
    help: str = ""
    #: Inclusive bounds for a numeric parameter, where the tool it is passed to
    #: has them. The dialog refuses values outside these rather than letting the
    #: aligner reject them after the run has started.
    lo: Optional[float] = None
    hi: Optional[float] = None


class AlignerEngine(ABC):
    key: str = "engine"
    label: str = "engine"

    @abstractmethod
    def available(self) -> bool: ...

    @abstractmethod
    def align(self, records: Sequence[Record], alphabet: Alphabet, **opts) -> Alignment: ...

    def parameters(self) -> List[Param]:
        """Tunable parameters exposed in the GUI settings dialog. Default: none."""
        return []

    def supports(self, alphabet: Alphabet) -> bool:
        """Can this engine align that alphabet? Most can align anything.

        Asked *before* a run is started, so an engine that cannot take the data
        is refused with an explanation rather than by a traceback from inside a
        worker thread.
        """
        return True

    def __repr__(self) -> str:  # pragma: no cover
        return f"<{self.__class__.__name__} {self.key}>"


# --------------------------------------------------------------------------- #
# Built-in
# --------------------------------------------------------------------------- #

class BuiltinProgressive(AlignerEngine):
    key = "builtin"
    label = "CRAIC built-in (progressive)"
    consistency_mem_gb = 1.0   # RAM budget (GB proxy) for the accurate consistency path;
                               # families larger than this stream. Real peak RAM ~2-2.5x this.
    #: Effort used when the caller names none. The GUI and CLI pass one explicitly;
    #: the benchmark harness sets this so that a result set records, and does not
    #: silently choose, the compute level it measured.
    benchmark_effort = "med"

    def available(self) -> bool:
        return True

    def align(self, records, alphabet, effort: Optional[str] = None, consistency_iters=None,
              delta=None, epsilon=None, estimate: bool = True, progress=None,
              cancelled=None, **opts):
        effort = self.benchmark_effort if effort is None else effort
        return progressive.align(records, alphabet, effort=effort, delta=delta,
                                 epsilon=epsilon, estimate=estimate,
                                 consistency_iters=consistency_iters,
                                 consistency_mem_gb=self.consistency_mem_gb,
                                 progress=progress, cancelled=cancelled)

    def parameters(self):
        return [
            Param("effort", "Compute level", "choice", "med", ["min", "med", "max"],
                  help="min = fast streaming posterior decoding (scales to large "
                       "alignments); med = consistency transformation + posterior "
                       "guide tree (default, most accurate per unit time); "
                       "max = also runs iterative refinement (slower)."),
        ]


# --------------------------------------------------------------------------- #
# External wrappers
# --------------------------------------------------------------------------- #

_SEARCH_PATH: Optional[str] = None


def _search_path() -> str:
    """PATH to search for external aligners. A GUI / .app launch on macOS inherits
    only a minimal PATH (``/usr/bin:/bin:...``), so tools installed via conda,
    Homebrew, MacPorts or pip go undetected. Augment it with the usual locations."""
    global _SEARCH_PATH
    if _SEARCH_PATH is None:
        extra = [
            os.path.join(sys.prefix, "bin"),                       # active venv / conda env
            os.path.join(os.environ.get("CONDA_PREFIX", ""), "bin"),
            "/opt/homebrew/bin", "/usr/local/bin",                 # Homebrew (arm / intel)
            os.path.expanduser("~/.local/bin"),
            "/opt/local/bin",                                      # MacPorts
        ]
        seen, dirs = set(), []
        for d in os.environ.get("PATH", "").split(os.pathsep) + extra:
            if d and d not in seen and os.path.isdir(d):
                seen.add(d)
                dirs.append(d)
        _SEARCH_PATH = os.pathsep.join(dirs)
    return _SEARCH_PATH


def _which(binary: str) -> Optional[str]:
    """Locate an external tool on the augmented search path (absolute path or None)."""
    return shutil.which(binary, path=_search_path())


class ExternalAligner(AlignerEngine):
    binary: str = ""
    timeout = None   # per-call subprocess timeout (s); None = no limit (set by benchmarks)

    def available(self) -> bool:
        return _which(self.binary) is not None

    def _command(self, infile: str, outfile: str, alphabet: Alphabet, **opts) -> Tuple[List[str], str]:
        """Return (argv, output_mode) where mode is 'stdout' or 'file'."""
        raise NotImplementedError

    @classmethod
    def _run(cls, argv: List[str], timeout=None):
        """Run an external aligner on the same augmented PATH that found it.

        ``available()`` resolves the binary against ``_search_path()``, so a
        subprocess launched against the *inherited* PATH can fail to find a tool
        the engine has just reported as available — which is exactly what
        happens in a macOS ``.app`` launch, the case ``_search_path`` exists to
        fix. Every external invocation goes through here so the two cannot drift
        apart again.
        """
        exe = _which(argv[0]) or argv[0]
        return subprocess.run(argv, capture_output=True, text=True, timeout=timeout,
                              executable=exe, env={**os.environ, "PATH": _search_path()})

    def align(self, records, alphabet, **opts):
        from .io import read_records, write_fasta

        opts.pop("progress", None)            # external aligners run as one subprocess;
        opts.pop("cancelled", None)           # they don't report progress or cancel mid-run
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            outfile = os.path.join(d, "out.fasta")
            safe = [(f"s{i}", s.replace("-", "")) for i, (_, s) in enumerate(records)]
            write_fasta(infile, [x[0] for x in safe], [x[1] for x in safe], wrap=0)
            argv, mode = self._command(infile, outfile, alphabet, **opts)
            proc = self._run(argv, self.timeout)
            if proc.returncode != 0:
                raise RuntimeError(f"{self.label} failed: {proc.stderr[:400]}")
            if mode == "stdout":
                with open(outfile, "w") as fh:
                    fh.write(proc.stdout)
            parsed = read_records(outfile, "fasta")
            # restore original ids by position
            idmap = {f"s{i}": rec[0] for i, rec in enumerate(records)}
            parsed = [(idmap.get(pid, pid), seq.upper()) for pid, seq in parsed]
            order = {rec[0]: i for i, rec in enumerate(records)}
            parsed.sort(key=lambda p: order.get(p[0], 0))
            return Alignment.from_records(parsed, alphabet=alphabet)


class MafftEngine(ExternalAligner):
    key = "mafft"
    label = "MAFFT"
    binary = "mafft"

    _STRATEGY = {
        "auto": ["--auto"],
        "L-INS-i (accurate, local)": ["--localpair", "--maxiterate", "1000"],
        "G-INS-i (accurate, global)": ["--globalpair", "--maxiterate", "1000"],
        "E-INS-i (long gaps)": ["--genafpair", "--maxiterate", "1000"],
        "FFT-NS-2 (fast)": ["--retree", "2"],
    }

    def parameters(self):
        return [
            Param("strategy", "Strategy", "choice", "auto", list(self._STRATEGY),
                  help="MAFFT preset / iterative-refinement strategy."),
            Param("op", "Gap open (--op)", "float", 1.53),
            Param("ep", "Gap extend (--ep)", "float", 0.0),
        ]

    def _command(self, infile, outfile, alphabet, **opts):
        flags = self._STRATEGY.get(opts.get("strategy", "auto"), ["--auto"])
        op = float(opts.get("op", 1.53))
        ep = float(opts.get("ep", 0.0))
        argv = ["mafft", *flags, "--op", str(op), "--ep", str(ep), "--quiet", infile]
        return argv, "stdout"


class MuscleEngine(ExternalAligner):
    key = "muscle"
    label = "MUSCLE"
    binary = "muscle"

    def parameters(self):
        return []   # MUSCLE v5 doesn't expose gap penalties on the command line

    def align(self, records, alphabet, **opts):
        from .io import read_records, write_fasta

        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            outfile = os.path.join(d, "out.fasta")
            safe = [(f"s{i}", s.replace("-", "")) for i, (_, s) in enumerate(records)]
            write_fasta(infile, [x[0] for x in safe], [x[1] for x in safe], wrap=0)
            # try v5 syntax, then v3. A v5 attempt that fails to even start
            # (wrong binary, timeout) must not abort the v3 attempt, so the
            # failure is recorded and the loop continues.
            last_error = "no attempt ran"
            for argv in (
                ["muscle", "-align", infile, "-output", outfile],
                ["muscle", "-in", infile, "-out", outfile],
            ):
                try:
                    proc = self._run(argv, self.timeout)
                except (OSError, subprocess.SubprocessError) as exc:
                    last_error = str(exc)
                    continue
                if proc.returncode == 0 and os.path.exists(outfile):
                    break
                last_error = proc.stderr[:400]
            else:
                raise RuntimeError(f"MUSCLE failed: {last_error}")
            parsed = read_records(outfile, "fasta")
            idmap = {f"s{i}": rec[0] for i, rec in enumerate(records)}
            parsed = [(idmap.get(pid, pid), seq.upper()) for pid, seq in parsed]
            order = {rec[0]: i for i, rec in enumerate(records)}
            parsed.sort(key=lambda p: order.get(p[0], 0))
            return Alignment.from_records(parsed, alphabet=alphabet)


class ClustalOmegaEngine(ExternalAligner):
    key = "clustalo"
    label = "Clustal Omega"
    binary = "clustalo"

    def parameters(self):
        return [Param("iterations", "Combined iterations", "int", 0, lo=0, hi=100,
                      help="Guide-tree/HMM refinement iterations; 0 = default.")]

    def _command(self, infile, outfile, alphabet, **opts):
        argv = ["clustalo", "-i", infile, "-o", outfile, "--force", "--outfmt=fasta"]
        it = int(opts.get("iterations", 0))
        if it > 0:
            argv += ["--iterations", str(it)]
        return argv, "file"


class ProbConsEngine(ExternalAligner):
    """ProbCons (Do et al., 2005), the method CRAIC's built-in engine follows.

    Present so that the built-in engine can be compared with the thing it is a
    reimplementation of, rather than only with aligners of a different lineage.
    Without it, a benchmark table reporting the built-in engine below MAFFT
    invites the reader to conclude that posterior-consistency alignment is the
    weaker approach, when the honest question is how much of the gap is the
    method and how much is this implementation of it.

    Protein only: ProbCons has no nucleotide model.
    """

    key = "probcons"
    label = "ProbCons"
    binary = "probcons"

    def supports(self, alphabet: Alphabet) -> bool:
        return not alphabet.is_nucleotide

    def align(self, records, alphabet, **opts):
        if not self.supports(alphabet):
            raise RuntimeError(
                "ProbCons aligns protein sequences only. For coding nucleotides, "
                "tick 'align as protein' to translate, align and back-translate; "
                "otherwise choose another engine.")
        return super().align(records, alphabet, **opts)

    def parameters(self):
        return [Param("iterations", "Refinement iterations", "int", 100,
                      lo=0, hi=1000,
                      help="Rounds of random-bipartition iterative refinement "
                           "(ProbCons -ir). ProbCons accepts 0 to 1000; the "
                           "published default is 100.")]

    def _command(self, infile, outfile, alphabet, **opts):
        argv = ["probcons"]
        it = int(opts.get("iterations", 100))
        if it != 100:
            argv += ["-ir", str(it)]
        return argv + [infile], "stdout"


class PrankEngine(ExternalAligner):
    key = "prank"
    label = "PRANK (phylogeny-aware)"
    binary = "prank"

    def parameters(self):
        return [
            Param("gaprate", "Gap opening rate", "float", 0.025, lo=0.0, hi=1.0,
                  help="PRANK -gaprate. A probability, so 0 to 1."),
            Param("gapext", "Gap extension prob.", "float", 0.5, lo=0.0, hi=1.0,
                  help="PRANK -gapext. A probability, so 0 to 1."),
        ]

    def align(self, records, alphabet, **opts):
        from .io import read_records, write_fasta

        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            prefix = os.path.join(d, "out")
            safe = [(f"s{i}", s.replace("-", "")) for i, (_, s) in enumerate(records)]
            write_fasta(infile, [x[0] for x in safe], [x[1] for x in safe], wrap=0)
            argv = ["prank", f"-d={infile}", f"-o={prefix}",
                    f"-gaprate={float(opts.get('gaprate', 0.025))}",
                    f"-gapext={float(opts.get('gapext', 0.5))}"]
            if alphabet.is_nucleotide:
                argv.append("-DNA")
            proc = self._run(argv, self.timeout)
            best = prefix + ".best.fas"
            if not os.path.exists(best):
                best = prefix + ".2.fas"
            if not os.path.exists(best):
                raise RuntimeError(f"PRANK failed: {proc.stderr[:400]}")
            parsed = read_records(best, "fasta")
            idmap = {f"s{i}": rec[0] for i, rec in enumerate(records)}
            parsed = [(idmap.get(pid, pid), seq.upper()) for pid, seq in parsed]
            order = {rec[0]: i for i, rec in enumerate(records)}
            parsed.sort(key=lambda p: order.get(p[0], 0))
            return Alignment.from_records(parsed, alphabet=alphabet)


# --------------------------------------------------------------------------- #
# Codon-aware decorator (translate -> align aa -> back-translate)
# --------------------------------------------------------------------------- #

class CodonAware(AlignerEngine):
    """Wrap any engine to align coding nucleotides in amino-acid space.

    Assumes reading frame 0 and intact frames (no internal frameshifts). The
    result is a codon-aware nucleotide alignment carrying a CodingSpec, so the
    nt/codon/aa views are all immediately valid.
    """

    def __init__(self, inner: AlignerEngine, table: int = 1):
        self.inner = inner
        self.table = table
        self.key = f"codon:{inner.key}"
        self.label = f"{inner.label} + codon-aware"

    def available(self) -> bool:
        return self.inner.available()

    def parameters(self):
        return self.inner.parameters()

    def align(self, records, alphabet, **opts):
        aa_records: List[Record] = []
        codon_lists = {}
        for name, seq in records:
            s = seq.replace("-", "").upper().replace("U", "T")
            n = len(s) // 3
            codons = [s[3 * i : 3 * i + 3] for i in range(n)]
            aa = "".join(translate_codon(c, self.table) for c in codons)
            aa_records.append((name, aa))
            codon_lists[name] = codons

        aa_aln = self.inner.align(aa_records, Alphabet.PROTEIN, **opts)

        nt_rows = []
        for name, aa_row in zip(aa_aln.ids, aa_aln.rows):
            codons = codon_lists[name]
            ci = 0
            buf = []
            for ch in aa_row:
                if ch == "-":
                    buf.append("---")
                else:
                    buf.append(codons[ci] if ci < len(codons) else "---")
                    ci += 1
            nt_rows.append("".join(buf))
        aln = Alignment(list(aa_aln.ids), nt_rows, Alphabet.DNA,
                        coding=CodingSpec(frame=0, table=self.table))
        order = {rec[0]: i for i, rec in enumerate(records)}
        pairs = sorted(zip(aln.ids, aln.rows), key=lambda p: order.get(p[0], 0))
        return Alignment([p[0] for p in pairs], [p[1] for p in pairs], Alphabet.DNA,
                         coding=CodingSpec(frame=0, table=self.table))


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #

class _BuiltinVariant(BuiltinProgressive):
    """Built-in aligner pinned to a particular pair-HMM gap-open probability.

    Used so the disagreement map is informative even when no external aligner is
    installed: different gap models give genuinely different MEA alignments that
    disagree exactly in the hard regions.
    """

    def __init__(self, delta: float, label: str):
        self._delta = delta
        self.label = label
        self.key = f"builtin:d{delta}"

    def align(self, records, alphabet, **opts):
        return progressive.align(records, alphabet, delta=self._delta, epsilon=0.5,
                                 estimate=False, consistency_iters=1)


def builtin_variants() -> List[AlignerEngine]:
    # Vary the pair-HMM gap-open probability (delta): a rarer-gap vs freer-gap
    # model produces genuinely different posterior-decoded alignments, so the
    # disagreement map stays informative with no external aligner installed.
    return [
        _BuiltinVariant(0.005, "built-in · rare gaps"),
        _BuiltinVariant(0.03, "built-in · default gaps"),
        _BuiltinVariant(0.1, "built-in · frequent gaps"),
    ]


def all_engines() -> List[AlignerEngine]:
    return [BuiltinProgressive(), MafftEngine(), MuscleEngine(),
            ClustalOmegaEngine(), ProbConsEngine(), PrankEngine()]


def available_engines() -> List[AlignerEngine]:
    return [e for e in all_engines() if e.available()]


def get_engine(key: str) -> Optional[AlignerEngine]:
    for e in all_engines():
        if e.key == key:
            return e
    return None
