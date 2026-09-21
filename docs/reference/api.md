# Python API

CRAIC is a library as well as a workbench. Everything the GUI does is available
from Python, which is the supported way to script it or to build on it.

```python
from craic import io, progressive, evaluate
from craic.ambiguity import reliability
from craic.domain import Alphabet

records = io.read_records("seqs.fasta")
aln = progressive.align(records, Alphabet.PROTEIN)

report = reliability.analyse(aln, do_perturbation=True)
print(report.col_combined)                  # per-column reliability, 0..1

truth = io.load_alignment("truth.fasta", alphabet=Alphabet.PROTEIN)
acc = evaluate.compare_to_reference(aln, truth)
print(acc.sp, acc.tc, acc.reliability_auc(report.col_combined))
```

Nothing below imports Qt, so all of it works headless.

---

## The canonical model

::: craic.domain
    options:
      members:
        - Alphabet
        - Level
        - Alignment
        - CodingSpec
        - genetic_codes
        - code_name
        - translate_codon
        - residue_index
        - column_index
        - membership
        - normalise_gaps
        - ungap
        - match_ids

---

## Reading and writing

::: craic.io
    options:
      members:
        - read_records
        - load_alignment
        - write_fasta

---

## The acceleration layer

The single primitive everything else is built on. The Rust core and the NumPy
mirror satisfy the same contract and are cross-validated to ~1e-9; callers cannot
tell which answered.

::: craic.accel
    options:
      members:
        - backend
        - pair_posteriors

---

## The built-in aligner

::: craic.progressive
    options:
      members:
        - align
        - estimate_params
        - consistency_transform
        - sparse_available

---

## Reliability

::: craic.ambiguity.reliability
    options:
      members:
        - analyse
        - consistency
        - perturbation
        - keep_mask
        - EnsembleFailure

---

## Accuracy against a known answer

::: craic.evaluate
    options:
      members:
        - compare_to_reference
        - Accuracy
        - matched_rows
        - sp_score
        - tc_score
        - col_correct_fraction
        - cell_correct_fraction
        - true_column_grid
        - explain_column
        - ColumnTruth
        - Placement
        - auc

---

## Simulation

::: craic.simulate
    options:
      members:
        - simulate

---

## Model-free trimming

Faithful reimplementations of the published column-selection rules, offered as
complementary signals to the model-based reliability score.

::: craic.ambiguity.trimming
    options:
      members:
        - gap_score
        - similarity_score
        - gappyout
        - strict
        - gblocks

---

## Sessions

::: craic.session
    options:
      members:
        - Session
        - is_session_file
        - SessionError
