# Concepts & methods

The short version. For the complete specification — the pair-HMM's parameters,
how they are estimated, the decoding criterion, the ensemble design and the
complexity of each step — see **[Methods in full](methods.md)**.

## The pair-HMM core

Everything CRAIC knows about alignment uncertainty comes from one primitive: the **posterior
match probabilities** of a three-state pair-HMM (match / insert / delete, affine gaps),
computed by forward–backward in log space — the ProbCons idea. For two sequences it yields a
matrix `P[i, j] = P(residue i of A is homologous to residue j of B)`.

That single matrix powers:

- the **posterior explorer** (it renders `P` directly),
- the **consistency** reliability score (the average `P` over the residue pairs a column
  asserts),
- the built-in **maximum-expected-accuracy** aligner (it decodes `P` for the alignment that
  maximises the expected number of correctly-aligned residue pairs).

The core is implemented in **Rust** (`src/lib.rs`, via PyO3) and is biology-agnostic: it takes
integer-encoded sequences plus an emission model and returns the matrix. A **NumPy** mirror in
`craic/accel.py` reproduces it line for line, and the test-suite cross-validates the two to
~1e-9. If the Rust core isn't compiled, CRAIC simply uses the mirror.

### Emission models

The emission model (built in Python, not Rust) supplies the joint match probabilities and
background frequencies. Protein uses **BLOSUM45** converted to a proper joint via Robinson
background frequencies; nucleotides use an identity-driven model. An unknown residue maps to a
neutral wildcard, contributing no spurious homology signal.

## Reliability: two signals

- **Consistency** (in the spirit of T-Coffee's TCS and ZORRO): for each column, the mean
  pair-HMM posterior over the residue pairs it groups together. High = the unaligned-sequence
  evidence supports the column.
- **Perturbation** (a stability probe in the spirit of GUIDANCE, but not equivalent to it —
  16 replicates of CRAIC's own engine, not ~100 of the aligner that built the alignment):
  re-align the same sequences across an ensemble
  of bootstrapped guide trees — guide-tree uncertainty is a dominant source of alignment error —
  cycling through a few gap regimes, and measure what fraction of each residue's asserted
  homologies survive. A residue with no asserted partners is left undefined, not scored a free 1.

Both live in `[0, 1]`; the combined score averages them, and the mask is just a threshold.

## Disagreement

Run several aligners, project them onto a reference's coordinates, and for each column measure
how often its homology assertions are reproduced by the other alignments. Columns every method
agrees on are trustworthy; the rest are flagged.

## The built-in progressive aligner

A deliberately modest, dependency-free aligner that is *always* available, built entirely on the
pair-HMM posterior. It **estimates** the emission divergence and gap parameters from the data
(no assumption that sequences are ~90% identical), computes every pairwise posterior, applies the
**ProbCons consistency transformation** (re-estimating each posterior using every third
sequence), builds its **guide tree from the posteriors themselves** (a far better signal than
k-mer overlap), and aligns progressively by **pure maximum-expected-accuracy** decoding — so gap
placement falls out of the model rather than from tuned penalties. A **compute level** trades speed
for accuracy: *min* streams plain decoding (scales to large alignments), *med* (default) adds the
consistency transformation and posterior guide tree, and *max* also runs iterative refinement. On simulated data with a
known-true alignment, estimating parameters and applying the consistency transformation raises
mean sum-of-pairs accuracy from ~0.62 (plain decoding) to ~0.82. For large families, where
holding every pairwise posterior in memory would be prohibitive, it automatically **streams** --
computing each posterior on demand and discarding it, so peak memory is a single matrix -- trading
the consistency transformation for scale. It has no sequence weighting and is meant to always be available, not to beat MAFFT; for production-grade alignment,
wire in an external engine.

## Other reliability signals

Alongside its model-based reliability, CRAIC also implements the common model-free
column trimmers as complementary signals: gap-fraction (MSA_trimmer-style),
trimAl-style gap and similarity scores with `gappyout` automatic trimming, and
Gblocks-style conserved blocks -- plus EvalMSA / OD-seq-style per-sequence outlier
detection. These are faithful reimplementations of the published algorithms
(`craic.ambiguity.trimming`); on benchmarks the pair-HMM posterior reliability
distinguishes mis-aligned columns better than the model-free methods, but the
latter are fast, familiar, and always available.

## Codon-aware alignment

For protein-coding nucleotides, CRAIC translates → aligns the amino acids → back-translates, so
gaps land on codon boundaries and the reading frame is preserved (assumes frame 0 and intact
frames; partial-gap codons are surfaced as `X`).

## Method lineage

The ideas build on ProbCons (posterior decoding), T-Coffee/TCS and ZORRO (consistency),
GUIDANCE (perturbation via guide-tree bootstrapping), PRANK (phylogeny-aware gaps), and MACSE (codon-aware
alignment). CRAIC's contribution is to put these uncertainty signals *in front of the user,
interactively*, in one viewer.
