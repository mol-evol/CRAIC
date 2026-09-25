# Glossary

Terms as CRAIC uses them. Where a word is used differently elsewhere in the
literature, that is noted, because most of the confusion around alignment
reliability is vocabulary rather than statistics.

**Alignment** — a rectangular arrangement of sequences in which each column
asserts that its residues descend from the same ancestral position. Every column
is a hypothesis; see *homology*.

**Ambiguous region** — a stretch where the evidence does not strongly prefer one
arrangement over its neighbours. Usually created by insertions and deletions:
substitutions alone rarely make a column doubtful.

**Codon-aware alignment** — aligning protein-coding nucleotides by translating,
aligning the amino acids, and threading the nucleotides back through the result,
so gaps land on codon boundaries and the reading frame survives.

**Combined score** — CRAIC's headline per-column reliability: the mean of the
*consistency* and *perturbation* scores. In `[0, 1]`.

**Conservation** — how similar the residues in a column are. Not reliability: a
column can be perfectly conserved and completely misplaced. CRAIC shows it as a
separate track for exactly this reason.

**Consistency score** — per column, the mean pair-HMM posterior over the residue
pairs the column groups together. High means the unaligned-sequence evidence
supports the column. In the spirit of T-Coffee's TCS and of ZORRO, but computed
over CRAIC's own posteriors and so not numerically interchangeable with either.

**Consistency transformation** — the ProbCons idea of re-estimating each pairwise
posterior using every third sequence, so that evidence from the whole family
informs each pair. Costs O(N³) in sequences and buys a large accuracy
improvement.

**Core block** — in a structural reference set such as BAliBASE, the region the
reference is actually confident about. Published BAliBASE SP and TC figures are
restricted to core blocks, so an all-column figure is not comparable with them.

**Correct placement** — CRAIC's per-residue view of accuracy against a known
answer: the fraction of a residue's partners in its column that are genuinely
homologous to it. The per-residue counterpart of *reference correctness*, and
exactly consistent with it — a column's score is the mean of its residues'.

**Disagreement / aligner agreement** — how often the other installed aligners
reproduce each column of the current alignment. A model-free proxy for
uncertainty: where methods differ, something is being guessed.

**Gappyout / strict** — trimAl's automatic column-selection rules, reimplemented
in `craic.ambiguity.trimming`. They choose their own threshold from the
distribution of gap scores rather than taking one from the user.

**Genetic code** — the NCBI translation table used to turn codons into amino
acids. CRAIC carries one per alignment and defaults to table 1 (standard).
Mitochondrial and Mollicute genomes differ, most consequentially in reading TGA
as tryptophan rather than as a stop.

**Guide tree** — the order in which a progressive aligner merges sequences.
Errors in it are a dominant source of alignment error, which is why perturbing it
is informative.

**Homology** — here, positional homology: two residues are homologous if they
descend from the same ancestral position. A column asserts homology among all its
residues; SP and TC both count homologous *pairs* and *columns* respectively.

**Internal stop codon** — a stop before a sequence's last codon. Usually means
the wrong genetic code, the wrong reading frame, or a pseudogene; CRAIC reports
them but cannot say which.

**Intruder** — in the column inspector, a residue placed in a column that the
reference puts somewhere else.

**Masking** — removing columns judged unreliable before a downstream analysis,
or, with a **residue mask**, only the residues judged unreliable (written as
missing data, `N` or `X`, with the column kept).
CRAIC previews it live (dimmed columns) and exports the result. Whether it helps
is an empirical question, and a threshold that removes hard columns will raise
mean accuracy whatever rule picks them — hence the random and gap-matched
controls in the benchmark.

**MEA (maximum expected accuracy)** — decoding an alignment that maximises the
expected number of correctly aligned residue pairs under the posterior, rather
than one that maximises a sum of substitution scores minus gap penalties. Gap
placement falls out of the model instead of being tuned.

**Missing residue** — in the column inspector, a residue the reference puts in
the column being inspected but which the alignment placed elsewhere. A column can
score a perfect 1.0 and still have one: the score asks whether the homologies a
column *asserts* are true, and a column that left a residue out asserts nothing
false.

**Pair-HMM** — the three-state (match / insert / delete) hidden Markov model with
affine gaps whose forward–backward posteriors are the single primitive
everything in CRAIC is built on.

**Perturbation score** — per column, the fraction of its asserted homologies that
survive re-alignment across an ensemble of bootstrapped guide trees and gap
regimes. Inspired by GUIDANCE but not equivalent to it: 16 replicates of CRAIC's
own engine, not ~100 of the aligner that built the alignment.

**Posterior (match probability)** — `P[i, j]`, the probability that residue *i*
of one sequence is homologous to residue *j* of another, under the pair-HMM. The
posterior explorer draws this matrix directly.

**Precision** — of an alignment against a reference, the fraction of the
homologous pairs it *asserts* that are true. The complement of SP, which is
recall. Reported because a gappy alignment can have high precision and terrible
recall.

**Reference alignment** — a trusted alignment of the same sequences: a structural
reference such as BAliBASE, or the known-true alignment of a simulated dataset.

**Reference correctness** — per column, the fraction of the homologous pairs it
asserts that are true, given a reference. `nan` for a column that asserts no pair
at all, which can be neither right nor wrong.

**Reliability** — how well supported a column is, as distinct from how conserved
it is. CRAIC's reliability is the combined score.

**Reliability AUC** — the ROC area under the curve of a reliability score
predicting which columns are actually correct, given a known answer. 0.5 is
chance. This is the number that says whether the reliability claim is worth
anything on a given dataset.

**Session** — CRAIC's own document (`*.craic.json`): the alignment plus the
reference, annotations, sequence groups, provenance log and view state. A FASTA
export keeps the residues and drops all of that.

**SP (sum-of-pairs) score** — the fraction of the reference's homologous residue
pairs that the alignment recovers. Recall, in other words.

**TC (total column) score** — the fraction of the reference's columns that the
alignment reproduces exactly. A much harsher measure than SP: one misplaced
residue loses the whole column.

**Truth mode** — working with a reference attached, so that accuracy is scored
and displayed live. CRAIC's teaching centre of gravity.
