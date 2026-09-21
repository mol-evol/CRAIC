# Background & intent

## The problem

A multiple sequence alignment is an inference, not a measurement. Every column
is a hypothesis that the residues in it descend from the same ancestral
position, and some of those hypotheses are very much better supported than
others. In conserved regions the answer is obvious and every method agrees. In
regions with indels, in fast-evolving loops, near the ends of sequences, the
answer is a guess — and the guess is made silently, by software, and then handed
downstream as though it were data.

What happens downstream is the actual problem. Phylogenetic inference, selection
tests, ancestral reconstruction and structural mapping all take the alignment as
given. The error does not announce itself; it propagates. Bootstrap values on a
tree quantify the sampling variance of the tree *given the alignment* — they say
nothing about whether the alignment underneath was right, and they will happily
give strong support to a clade that exists only because a hundred residues were
placed in the wrong columns.

This is well known, and it has been well known for a long time. The field's
response has been a set of good methods for measuring alignment reliability —
transitive consistency (T-Coffee's TCS), guide-tree bootstrapping (GUIDANCE),
posterior-probability masking (ZORRO, Divvier), model-free trimming (trimAl,
Gblocks, ClipKIT), and alignment ensembles (MUSCLE 5). What has not happened is
that these measures became part of how people *look* at alignments.

## Why the measures don't get used

The reliability methods are, almost without exception, command-line programs or
web services that consume an alignment and emit a file of numbers. Using one
means leaving whatever you were doing, running a second program, and then
finding a way to look at the result. The commonest outcome is that the step is
skipped: an alignment is built, eyeballed in a viewer, and passed on.

The viewers, meanwhile, are excellent at everything except this. Jalview,
AliView, SeaView, MEGA and UGENE are all actively maintained and all let you
read and edit an alignment fluently. What they show you by default is
*conservation* — how similar the residues in a column are — which is a different
thing from reliability and, in the cases that matter most, close to the opposite
of it. A column of identical residues in a well-anchored block is both conserved
and reliable. A column of identical residues that three of the aligners would
have placed somewhere else entirely is conserved and unreliable, and nothing on
screen distinguishes the two.

Jalview is the important exception and it is worth being precise about what it
does. Since version 2.8 (2012) it has been able to import a T-Coffee
`score_ascii` file and colour the alignment by it, and GUIDANCE has shipped
Jalview-compatible output since 2011. So reliability *can* be looked at in a
GUI, and has been able to be for over a decade. What Jalview does not do is
compute it: the score has to be produced elsewhere, by a command-line round trip,
and what arrives is a static colouring of one particular alignment. Edit the
alignment and the colouring is stale. Realign it and the file is meaningless.
There is nothing to interrogate — no way to ask *why* a column scored badly, or
what the alternative was.

## What CRAIC is for

CRAIC's premise is that alignment uncertainty should be a live object in the
workbench rather than an artefact imported into it. Concretely, that means four
things that follow from computing the uncertainty in the same process that holds
the alignment:

**It is always available.** The reliability analysis is a menu item, not a
second program. There is no format conversion, no round trip, and no reason to
skip the step.

**It survives editing.** Move a block by hand and the scores are recomputed.
That turns curation from a blind activity into a closed loop: you can see
whether what you just did made the alignment better supported or worse.

**It can be interrogated.** A low score is a starting point, not a verdict. The
posterior explorer shows the actual matrix of homology probabilities for a pair
of sequences, so you can see the alternative the aligner rejected. The
realignment sandbox re-aligns just the ambiguous block under several engines and
gap regimes, so you can see what the competing hypotheses look like. The
disagreement map shows where the installed aligners actually part company.

**It can be checked against the truth.** CRAIC will generate sequences down a
random tree, so the true alignment is known exactly, and then show you per column
and per residue what the aligner got wrong — beside the reliability score that
tried to predict it. This is the part that makes the reliability claim testable
by the user rather than taken on trust, and it is also, in our experience, the
fastest way to teach what an alignment error is.

The last of these has a partial precedent worth naming: **SuiteMSA** (Anderson
*et al.*, 2011) wrapped indel-Seq-Gen in a Java GUI and compared an inferred
alignment against the simulated truth position by position. It is the closest
prior art to CRAIC's teaching mode, and the distinction is that SuiteMSA compared
alignments to each other, whereas CRAIC puts the comparison next to a
*statistical* uncertainty estimate — the point being to ask whether the estimate
predicted the error. SuiteMSA's published download URL has been returning 403 for
some time, so it also appears no longer to be obtainable.

## What CRAIC is not

It is not a production aligner. The built-in progressive engine exists so that
the program works out of the box and so that the perturbation analysis has
something to re-align with; it is a ProbCons-style maximum-expected-accuracy
implementation with a consistency transform, it is honest, and on BAliBASE it is
beaten by MAFFT. Anyone aligning real data for publication should point CRAIC at
MAFFT, MUSCLE, Clustal Omega or PRANK, all of which it will drive if they are on
your `PATH`. See [Validation](validation.md) for the numbers.

It is not a phylogenetics package. It does not infer trees except as a
diagnostic inside the benchmark, and it has no intention of competing with
IQ-TREE, RAxML or MrBayes. It ends where the alignment ends.

It is not a claim to have invented alignment reliability. Every score it
computes is an implementation of somebody else's published idea, and the
[methods page](concepts.md) says whose. The contribution is putting them where a
person can see them, argue with them, and check them against a known answer.

## Who it is for

Two audiences, which turn out to want the same thing.

**Researchers** curating an alignment before a downstream analysis, who need to
know which regions they should not trust and would like to do something about it
without leaving the program.

**Students**, and the people teaching them. Alignment is usually taught as a
solved preliminary — press Align, accept the output, move on to the tree — partly
because with real data there is no answer to compare against, so there is nothing
to discover. Given a dataset whose answer is known, a student can watch alignment
accuracy collapse as divergence rises, find the point where it does, test whether
masking the flagged columns actually helps, and repair a wrong column by hand.
That is a different lesson from being told that alignments contain errors.
