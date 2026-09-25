# CRAIC

**Interactive interrogation of multiple sequence alignment uncertainty.**

CRAIC (*Conserved-Region Alignment by Iterative Convergence*) is a free desktop
workbench for building, viewing and, above all, *interrogating* multiple sequence
alignments. It treats an alignment as a set of hypotheses to stress-test, not a settled
answer, and it is built around the ambiguously aligned regions where downstream
phylogenetics so often goes quietly wrong.

[Download the app](installation.md#download-a-ready-to-run-app){ .md-button .md-button--primary }
[Install with pip](installation.md){ .md-button }
[Get started](getting-started.md){ .md-button }
[How to cite](about.md#how-to-cite-craic){ .md-button }

![The CRAIC alignment viewer, nucleotide view with a reliability track](assets/reliability-track.png){ width="100%" }

## What it does

<div class="grid cards" markdown>

-   **See where the alignment can be trusted**

    ---

    A reliability score for every column and residue, from one pair-HMM, recomputed as
    you edit. Preview a mask with a slider, or mask single residues by hand.

    [Ambiguity tools](guide/ambiguity-tools.md)

-   **Compare aligners**

    ---

    Run MAFFT, MUSCLE, Clustal Omega, ProbCons, PRANK, ClustalW or the built-in
    aligner, with their settings exposed, and colour each column by how often they
    agree.

    [Engines & settings](reference/engines-and-formats.md)

-   **Realign one region**

    ---

    Select a block, realign it several ways, compare the alternatives side by side, and
    splice the best back in without disturbing the rest.

    [Realignment sandbox](guide/ambiguity-tools.md)

-   **Look inside the model**

    ---

    The posterior explorer shows the alternative ways a residue could align, not just
    the one the aligner chose.

    [Posterior explorer](guide/ambiguity-tools.md)

-   **Teach with a known answer**

    ---

    Simulate sequences whose true alignment is known, align them, and see exactly which
    residues each aligner, or each student, got wrong.

    [Teaching](guide/teaching.md)

-   **Script it**

    ---

    `craic score`, `mask`, `align`, `trim` and `simulate` run headless, on a laptop, a
    cluster node or in a pipeline.

    [Command line](reference/command-line.md)

</div>

Nucleotides, codons and amino acids come from one representation, so you can align at
one level and read the result at another.

## Does the reliability score work?

On BAliBASE 3 and on simulations with a known answer, the score picks out misaligned
columns well. Deleting the columns it flags still makes trees worse, and so does masking
the flagged residues, so CRAIC is built for inspecting and curating an alignment rather
than filtering it automatically. The evidence is in **[Validation](validation.md)**, and
what the score cannot do is in **[Limitations](limitations.md)**.

Wondering why this exists when you already have an alignment viewer? Read
**[Background & intent](background.md)**. Wondering how it works? Read
**[Concepts & methods](concepts.md)**.

## Citing CRAIC

!!! quote "If you use CRAIC, please cite"
    McInerney J. 2026. CRAIC: interactive interrogation of multiple sequence alignment
    uncertainty. bioRxiv preprint (link to follow).

The paper is under review; the **[About page](about.md)** has the current citation,
BibTeX, and the methods to cite alongside it.

---

CRAIC is written by **[James McInerney](https://mol-evol.github.io/)**, University of
Liverpool. It is free and open source under the MIT licence; the code is on
[GitHub](https://github.com/mol-evol/craic).
