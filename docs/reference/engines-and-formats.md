# Engines & file formats

## Alignment engines

CRAIC discovers engines at startup: the built-in aligner is always present, and each external
tool is offered if its binary is found on your `PATH`.

| Engine | Requires | How CRAIC runs it |
| --- | --- | --- |
| **CRAIC built-in (progressive)** | nothing | ProbCons-style pair-HMM posteriors → consistency → progressive alignment by posterior decoding. See [Concepts](../concepts.md). Always available. |
| **MAFFT** | `mafft` | `mafft <strategy> --op --ep --lop --lep <matrix> --quiet` |
| **MUSCLE** | `muscle` | `muscle -align IN -output OUT` (v5; `-super5`, `-perm`, `-perturb` when set), falling back to `muscle -in IN -out OUT` (v3) |
| **Clustal Omega** | `clustalo` | `clustalo -i IN -o OUT --force --outfmt=fasta` |
| **ProbCons** | `probcons` | `probcons IN` — protein only |
| **PRANK (phylogeny-aware)** | `prank` | `prank -d=IN -o=OUT` (`-DNA` added for nucleotides) |
| **ClustalW** | `clustalw2` or `clustalw` | `clustalw2 -INFILE=IN -OUTFILE=OUT -OUTPUT=FASTA -TYPE=…` |

### Parameters

The **⚙** button next to the engine dropdown opens a settings dialog for the selected engine,
with a **Restore Defaults** button; hover over a field for what it does and its range. The
dialog shows only the settings that apply to the sequences the engine will be handed: with
nucleotides it hides protein matrices, and with **align as protein** ticked it hides the DNA
settings, since the engine then sees amino acids. Hidden settings keep their values. Blank
fields say what the tool will use for that kind of data. On the
command line the same settings are `--param KEY=VALUE` (see
[`craic align`](command-line.md#craic-align)), and `craic engines` lists every key.

| Engine | Settings (key: default) |
| --- | --- |
| Built-in | compute level `effort`: med · gap open probability `delta`: *estimated from the data* · gap extension probability `epsilon`: *estimated from the data* · protein `matrix`: BLOSUM45 (also 50, 62, 80, 90) · `consistency_iters`, `refine_iters`: *set by the compute level* |
| MAFFT | `strategy`: auto (L-INS-i, G-INS-i, E-INS-i, FFT-NS-2) · `op`: 1.53 · `ep`: 0.0 · local-pair `lop`: −2.0 and `lep`: 0.1 (used by L-INS-i and E-INS-i) · `maxiterate`: *set by the strategy* · protein `matrix`: BLOSUM62 (30, 45, 80, JTT, transmembrane) · DNA `kimura`: 200 PAM (20, 1) |
| MUSCLE 5 | `algorithm`: align (or super5) · guide-tree `perm`: none (abc, acb, bca) · `perturb` seed: 0 = none. MUSCLE 5 has no gap penalties on its command line. |
| Clustal Omega | combined `iterations`: 0 · `full` distance matrix: off · `max_guidetree_iterations`, `max_hmm_iterations`: *same as iterations*. Clustal Omega aligns profile HMMs and has no gap penalties to set. |
| ProbCons | consistency passes `consistency`: 2 (0–5) · refinement `iterations`: 100 (0–1000) · `pretraining`: 0 (0–20). Its gap probabilities live in a parameter file and are not offered. |
| PRANK | `gaprate`, `gapext`: *PRANK's default* (DNA 0.025 / 0.75, protein 0.005 / 0.5) · `F` (+F, keep insertions as insertions): off · `iterate`: 5 · `termgap`: off |
| ClustalW | multiple-alignment `gapopen`, `gapext` and pairwise `pwgapopen`, `pwgapext`: *ClustalW's default* (protein 10 / 0.2 and 10 / 0.1, DNA 15 / 6.66) · protein `matrix` series: GONNET (BLOSUM, PAM, ID) · `dnamatrix`: IUB (CLUSTALW) · `gapdist`: 4 · `nopgap`, `nohgap` (turn off residue-specific and hydrophilic gap adjustment): off |

A setting shown in *italics* is left blank by default, and CRAIC then passes nothing for it,
so the tool uses its own value. That matters where the tool's default depends on the data:
PRANK and ClustalW use different gap costs for DNA and protein, and one number filled in for
both would be wrong for one of them.

ClustalW is there for teaching. Its successors model gaps inside profile HMMs, so there is no
single gap penalty to change; ClustalW keeps the textbook opening and extension costs, separately
for the pairwise stage that builds the guide tree and the progressive stage, which makes it the
engine to use to show what gap costs do.

Settings are remembered per engine for the session and apply everywhere that engine runs: the
**Align** button, **Realign around pinned columns**, the realignment sandbox, the ensemble
co-occurrence figure and the **Aligner agreement** track. The alignment's **History** records the
settings each run used, leaving out the ones left to the tool.

Any engine can **align as protein** (the toolbar checkbox): CRAIC translates each
sequence, aligns the amino acids with that engine, and back-translates onto codon
boundaries — keeping the original nucleotides in memory and switching the view to amino acid.
This assumes reading frame 0 and intact frames; partial-gap codons are surfaced as `X`. The
genetic code defaults to the NCBI **standard** table (id 1). Codon-aware wrapping inherits the
inner engine's parameters.

!!! note "Aligner agreement with no external tools"
    If fewer than two engines are installed, the **Aligner agreement** track falls back to three
    gap regimes of the built-in aligner — gap-open probability 0.005 (*rare gaps*), 0.03
    (*default*) and 0.1 (*frequent gaps*) — which disagree precisely in the hard regions.

## File formats

### Reading

Format is inferred from the extension:

| Extensions | Format |
| --- | --- |
| `.fa` `.fasta` `.fna` `.faa` | FASTA |
| `.phy` `.phylip` | PHYLIP (relaxed) |
| `.aln` | Clustal |
| `.sto` `.stk` | Stockholm |
| `.nex` `.nexus` | NEXUS |

Anything else is tried as FASTA. If the records aren't all the same length (i.e. not yet
aligned), CRAIC pads them for display and prompts you to **Align**.

### Writing

**Save…** offers a range of formats — the format comes from the type you choose in the dialog:

| Format | Extension |
| --- | --- |
| FASTA | `.fasta` |
| Clustal | `.aln` |
| PHYLIP — interleaved / relaxed / sequential | `.phy` |
| Stockholm | `.sto` |
| NEXUS | `.nex` |
| MEGA | `.meg` |

**Export masked…** writes FASTA. **Copy** (Ctrl/⌘+C) puts the current selection (or the whole
alignment) on the clipboard as FASTA at the level you're viewing, and **Export image…** saves
the current view as a PNG. Writing uses Biopython, so adding another format is a small change in
`craic/io.py`.

## Command line

```bash
python -m craic [FILE]      # launch the GUI, optionally opening FILE
python -m craic --version   # print the version and exit
```

Two double-clickable launchers live in the project root and always use the bundled `.venv`:

| Launcher | Purpose |
| --- | --- |
| `craic.command` | Launch the workbench. |
| `serve-docs.command` | Serve this documentation locally. |
