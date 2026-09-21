"""Generate small example datasets with deliberately ambiguous regions.

Run:  python examples/make_examples.py
Produces unaligned FASTA files in the examples/ directory.
"""

import os
import random

from Bio.Data import CodonTable

HERE = os.path.dirname(__file__)
rng = random.Random(42)

# one codon per amino acid (first synonymous) from the standard table
_T = CodonTable.unambiguous_dna_by_id[1]
_CODON = {}
for codon, aa in _T.forward_table.items():
    _CODON.setdefault(aa, codon)


def backtranslate(protein: str) -> str:
    return "".join(_CODON.get(a, "NNN") for a in protein)


def mutate_protein(seq: str, rate: float) -> str:
    aas = "ARNDCQEGHILKMFPSTWYV"
    out = []
    for a in seq:
        out.append(rng.choice(aas) if rng.random() < rate else a)
    return "".join(out)


def coding_genes(n=8):
    """Coding DNA: conserved cores flanking a variable-length, variable loop."""
    core_a = "MKAILVTLLAGVSPDEFGHIKL"
    core_b = "MNPQRSTVWYACDEFGHIKLMN"
    recs = []
    for i in range(n):
        loop_len = rng.randint(2, 9)
        loop = "".join(rng.choice("GSAPTNQ") for _ in range(loop_len))
        prot = mutate_protein(core_a, 0.06) + loop + mutate_protein(core_b, 0.06)
        recs.append((f"gene_{i+1}", backtranslate(prot)))
    return recs


def rrna_loops(n=8):
    """Non-coding RNA-like: conserved stems with variable hypervariable loops."""
    stem5 = "GGCUACGUAGCUAGCUGAUCG".replace("U", "T")
    stem3 = "CGAUCAGCUAGCUACGUAGCC".replace("U", "T")
    spacer = "AAUUGGCCAAUUGGCC".replace("U", "T")
    recs = []
    for i in range(n):
        v1 = "".join(rng.choice("ACGT") for _ in range(rng.randint(4, 14)))
        v2 = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 10)))
        seq = stem5 + v1 + spacer + v2 + stem3
        recs.append((f"rRNA_{i+1}", seq))
    return recs


def write(path, recs):
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n{seq}\n")
    print("wrote", path, f"({len(recs)} seqs)")


if __name__ == "__main__":
    write(os.path.join(HERE, "coding_genes.fasta"), coding_genes())
    write(os.path.join(HERE, "rrna_loops.fasta"), rrna_loops())
