# Command line

`craic` with no arguments, or with a file, opens the workbench as it always has:

```bash
craic
craic examples/coding_genes.fasta
```

Everything else is a subcommand that runs **headless**. The headless path never
imports PySide6, so it works on a cluster node, in a container, or inside a
pipeline with no display and no Qt installed.

```bash
craic --help            # the subcommands
craic --version
craic --core            # 'rust' or 'numpy': which acceleration core is active
```

Output goes to standard output when `-o` is omitted, and progress and summary
lines go to standard error, so the subcommands pipe cleanly.

---

## `craic align`

Align sequences with any available engine.

```bash
craic align seqs.fasta -o aln.fasta
craic align seqs.fasta -o aln.fasta --engine mafft
craic align cds.fasta  -o aln.fasta --codon        # align as protein, back-translate
craic align mito.fasta -o aln.fasta --codon --code 2   # vertebrate mitochondrial
```

| Option | Meaning |
|---|---|
| `--engine` | `builtin` (default), `mafft`, `muscle`, `clustalo`, `prank` — whichever are on your `PATH` |
| `--effort` | built-in engine compute level: `min`, `med` (default), `max` |
| `--codon` | translate → align amino acids → back-translate, so gaps respect the reading frame |
| `--code N` | NCBI genetic-code table for `--codon` (default 1, the standard code) |
| `--alphabet` | force `dna`, `rna` or `protein` instead of detecting it |
| `--format` | output format; inferred from the extension otherwise |

!!! warning "Set `--code` for organelle and Mollicute data"
    Vertebrate mitochondria are table 2, *Mycoplasma* and *Spiroplasma* table 4,
    bacteria and plastids table 11. All of these read **TGA as tryptophan**; the
    standard code reads it as a stop, and a stop is not one of the twenty
    residues, so it contributes no homology signal at all. CRAIC warns when a
    coding alignment contains internal stop codons, which is the usual symptom.

## `craic codes`

List the NCBI genetic-code tables, so `--code N` does not need a web search.

```bash
craic codes
```

## `craic score`

Per-column reliability as a tab-separated table: `column`, `consistency`,
`perturbation`, `combined`, `gap_score`, `similarity`.

```bash
craic score aln.fasta -o scores.tsv
craic score aln.fasta --fast -o scores.tsv                    # consistency only
craic score aln.fasta --reference truth.fasta -o scores.tsv   # add truth
```

With `--reference`, a `correct_fraction` column is added and SP, TC, precision
and the reliability AUC are printed to stderr. The reference must be an alignment
of the same sequences; one that is not is refused with an explanation rather than
scored as a very bad alignment.

| Option | Meaning |
|---|---|
| `--fast` | skip the perturbation ensemble; consistency only, much quicker |
| `--replicates` | perturbation replicates (default 16) |
| `--reference` | a trusted alignment of the same sequences |

## `craic mask`

Drop columns below a reliability threshold.

```bash
craic mask aln.fasta -o masked.fasta --threshold 0.5
craic mask aln.fasta -o masked.fasta --which consistency --unscored keep
```

| Option | Meaning |
|---|---|
| `--threshold` | keep columns scoring at least this (default 0.5) |
| `--which` | `combined` (default), `consistency` or `perturbation` |
| `--unscored` | `drop` (default) or `keep` — what to do with columns that had no residue-pair evidence at all |
| `--fast`, `--replicates` | as for `score` |

`--unscored drop` is the default and the conservative choice: a column that could
not be assessed is not evidence. It is the same definition the GUI's mask slider
and the benchmark both use — there is one implementation, so the mask you export
here is the mask that was evaluated.

A threshold that would remove every column is an error, not an empty file.

## `craic trim`

Rule-based trimming, as a comparison for the model-based mask.

```bash
craic trim aln.fasta -o trimmed.fasta --method gappyout
```

`--method` is `gappyout` (default), `strict` or `gblocks`. These are
dependency-free reimplementations of the published algorithms; for authoritative
numbers run the original trimAl or Gblocks.

## `craic simulate`

Generate sequences whose true alignment is known exactly.

```bash
craic simulate --taxa 12 --length 300 --divergence 0.4 --indel 3 \
    -o truth.fasta --unaligned seqs.fasta
```

| Option | Meaning |
|---|---|
| `--taxa`, `--length`, `--seed` | how many sequences, root length, and the seed (same seed, same data) |
| `--divergence` | branch-length upper bound; higher is more divergent |
| `--indel` | indel rate — the dial that creates alignment ambiguity |
| `--rate-alpha` | gamma shape for among-site rate variation; 0 = uniform rates |
| `--indel-zipf` | Zipf exponent for indel lengths; 0 = geometric |
| `--unaligned` | also write the unaligned sequences, ready to align |

`-o` writes the *true* alignment. Use it as the `--reference` for `craic score`.

!!! note
    The default settings — uniform rates, geometric indel lengths — are close to
    what CRAIC's own pair-HMM assumes, which flatters CRAIC relative to aligners
    with empirically-tuned parameters. `--rate-alpha` and `--indel-zipf` break
    that correspondence and are the harder test.

## `craic session`

Describe a saved session — what alignment it holds, whether a truth-mode
reference is attached, and the full provenance log of how the alignment was
arrived at.

```bash
craic session work.fasta.craic.json
craic session work.fasta.craic.json -o alignment.nex   # export its alignment
```

A session is readable JSON by design, but "what is in here and how did it get
that way" is the question actually being asked, so it gets an answer without a
text editor.

---

## A worked pipeline

Simulate, align with two engines, and ask which columns each got wrong:

```bash
craic simulate --taxa 10 --divergence 0.35 -o truth.fasta --unaligned seqs.fasta

for engine in builtin mafft; do
    craic align seqs.fasta --engine "$engine" -o "aln_$engine.fasta"
    craic score "aln_$engine.fasta" --reference truth.fasta \
        --fast -o "scores_$engine.tsv"
done
```

Each `score` run prints its SP, TC and reliability AUC to stderr, and the tables
hold reliability and correctness column by column.
