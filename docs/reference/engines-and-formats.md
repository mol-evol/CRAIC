# Engines & file formats

## Alignment engines

CRAIC discovers engines at startup: the built-in aligner is always present, and each external
tool is offered if its binary is found on your `PATH`.

| Engine | Requires | How CRAIC runs it | Tunable parameters (⚙) |
| --- | --- | --- | --- |
| **CRAIC built-in (progressive)** | nothing | k-mer distances → UPGMA order → Gotoh profile alignment. See [Concepts](../concepts.md). Always available. | gap open (−6), gap extend (−0.5) |
| **MAFFT** | `mafft` | `mafft <strategy> --op <op> --ep <ep> --quiet` | strategy (auto / L-INS-i / G-INS-i / E-INS-i / FFT-NS-2), `--op` (1.53), `--ep` (0.0) |
| **MUSCLE** | `muscle` | tries `muscle -align IN -output OUT` (v5), then `muscle -in IN -out OUT` (v3) | — (v5 exposes none on the CLI) |
| **Clustal Omega** | `clustalo` | `clustalo -i IN -o OUT --force --outfmt=fasta` | combined iterations (0 = default) |
| **PRANK (phylogeny-aware)** | `prank` | `prank -d=IN -o=OUT` (`-DNA` added for nucleotides) | gap rate (0.025), gap extension (0.5) |

### Parameters

The **⚙** button next to the engine dropdown opens a per-engine settings dialog (gap costs,
MAFFT strategy, and so on) with a **Restore Defaults** button. Choices are remembered per engine
for the session. They apply both when you press **Align** and when the **Aligner agreement**
track re-runs every installed engine — so a comparison reflects the exact settings you've dialled
in, while leaving you free to tweak one engine without touching the others. The alignment's
**History** records the parameters each run used.

Any engine can **align as protein** (the toolbar checkbox): CRAIC translates each
sequence, aligns the amino acids with that engine, and back-translates onto codon
boundaries — keeping the original nucleotides in memory and switching the view to amino acid.
This assumes reading frame 0 and intact frames; partial-gap codons are surfaced as `X`. The
genetic code defaults to the NCBI **standard** table (id 1). Codon-aware wrapping inherits the
inner engine's parameters.

!!! note "Aligner agreement with no external tools"
    If fewer than two engines are installed, the **Aligner agreement** track falls back to three
    gap regimes of the built-in aligner — *gentle* (open −4, extend −0.3), *default* (−6, −0.5),
    and *strict* (−10, −1.0) — which disagree precisely in the hard regions.

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
