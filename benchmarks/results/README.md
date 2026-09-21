# Benchmark results shipped with 0.5.0

The raw output behind [docs/validation.md](../../docs/validation.md), so the
tables there can be checked without re-running anything. The commands that
produced them are given at the foot of that page.

| file | what it is |
|---|---|
| `balibase_coreblock.csv` | official `bali_score` SP and TC over the reference core blocks — the figures comparable with published BAliBASE tables |
| `balibase_allcolumn.csv` | the harness's own per-family output: all-column SP and TC, the per-column reliability AUC, and the masking controls (kept / matched-random / gap-fraction, with a z score) |
| `simulation.csv` | the 360-dataset simulation sweep: same metrics plus Robinson–Foulds distance of the neighbour-joining tree, full alignment versus masked |

Environment: MAFFT 7.505, MUSCLE 5.1, Clustal Omega 1.2.4, BAliBASE 3.0, CRAIC
0.5.0 with the compiled Rust core.

Two scope notes that matter for reading these numbers:

* The `rel_auc` column is the **consistency** score alone. The runs were made
  without `--perturbation`, which is roughly an order of magnitude more
  expensive; the combined score the interface shows by default adds the
  perturbation ensemble on top.
* BAliBASE coverage is RV11 and RV12 in full plus 17 families of RV20. Larger
  families exceed the compute budget — the consistency analysis is cubic in the
  number of sequences — and families that any aligner failed to complete are
  excluded from **every** aligner's mean, so the comparison stays like for like.
