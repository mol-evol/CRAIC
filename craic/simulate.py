"""Ground-truth sequence simulator.

Used by the benchmark harness and by the workbench itself: a dataset whose true
alignment is known exactly is the only way to show, rather than assert, what a
reliability score is doing — which makes it as much a teaching device as a
benchmarking one.

Evolves DNA down a random binary tree with JC69 substitutions and insertions /
deletions, tracking the TRUE alignment (which residues are homologous) and the
TRUE topology. Divergence is controlled by the branch-length range and indel
load by ``indel_rate``. Because the homology of every residue is known exactly,
this is the clean way to evaluate a per-column reliability score.

**A caveat that belongs in any write-up using these numbers.** The default
generative model here — JC69, uniform base frequencies, one rate for every site,
geometric indel lengths — is close to the model CRAIC's pair-HMM assumes, so
CRAIC's inference model is nearly well-specified on this data while MAFFT and
MUSCLE carry empirical matrices tuned for real sequences. That is a structural
advantage to CRAIC, not a bug, but it must be disclosed. ``rate_alpha`` turns on
gamma-distributed among-site rate variation and ``indel_zipf`` gives indel
lengths a heavy tail, both of which break the correspondence; run the sweep with
and without them and report both.
"""

from __future__ import annotations

from typing import Dict

import numpy as np

BASES = "ACGT"


def random_tree(taxa: int, rng, bmin: float = 0.04, bmax: float = 0.28) -> dict:
    """Random rooted binary tree; node = {name, children: [(child, brlen), ...]}."""
    nodes = [{"name": f"t{i}", "children": []} for i in range(taxa)]
    forest = list(nodes)
    while len(forest) > 1:
        i, j = rng.choice(len(forest), 2, replace=False)
        a, b = forest[i], forest[j]
        parent = {"name": None, "children": [(a, float(rng.uniform(bmin, bmax))),
                                             (b, float(rng.uniform(bmin, bmax)))]}
        forest = [forest[k] for k in range(len(forest)) if k not in (i, j)] + [parent]
    return forest[0]


def _subst(base: str, t: float, rng) -> str:
    p = 0.75 * (1.0 - np.exp(-4.0 / 3.0 * t))      # JC69 probability of a change
    if rng.random() < p:
        return rng.choice([c for c in BASES if c != base])
    return base


def _indel_length(rng, zipf: float) -> int:
    """Indel length: geometric by default, Zipf (power-law) when ``zipf`` > 1.

    Real indel length distributions are heavy-tailed; a geometric distribution
    with p=0.5 almost never produces the long indels that make alignment hard.
    """
    if zipf and zipf > 1.0:
        return int(min(50, rng.zipf(zipf)))
    return 1 + int(rng.geometric(0.5))


class _Ctr:
    def __init__(self):
        self.id = 0
        self.pos: Dict[int, float] = {}

    def new(self, pos: float) -> int:
        i = self.id
        self.pos[i] = pos
        self.id += 1
        return i


def _evolve(sites, t, rng, ctr, indel_rate, rates=None, indel_zipf=0.0):
    # ``rates`` gives each site id its own relative substitution rate, so that
    # among-site rate variation is inherited consistently down the tree.
    s = [[i, _subst(b, t * (rates.get(i, 1.0) if rates is not None else 1.0), rng)]
         for i, b in sites]
    n_events = rng.poisson(indel_rate * t * max(1, len(s) / 50))
    for _ in range(int(n_events)):
        if rng.random() < 0.5 and len(s) > 2:                     # deletion
            L = min(len(s) - 1, _indel_length(rng, indel_zipf))
            k = rng.integers(0, len(s) - L + 1)
            del s[k:k + L]
        else:                                                     # insertion
            L = _indel_length(rng, indel_zipf)
            k = int(rng.integers(0, len(s) + 1))
            pl = ctr.pos[s[k - 1][0]] if k > 0 else (ctr.pos[s[0][0]] - 1 if s else 0.0)
            pr = ctr.pos[s[k][0]] if k < len(s) else (pl + 1.0)
            ins = [[ctr.new(pl + (m + 1) / (L + 1) * (pr - pl)), rng.choice(list(BASES))]
                   for m in range(L)]
            s[k:k] = ins
    return s


def simulate(taxa: int = 8, root_len: int = 200, seed: int = 0,
             indel_rate: float = 2.0, bmin: float = 0.04, bmax: float = 0.28,
             rate_alpha: float = 0.0, indel_zipf: float = 0.0) -> dict:
    """Return dict(names, seqs=[(name, ungapped)], true_rows=[gapped], tree).

    ``rate_alpha`` > 0 draws each site's relative substitution rate from a
    gamma(alpha, 1/alpha), i.e. among-site rate variation with mean 1; smaller
    alpha means more heterogeneity. ``indel_zipf`` > 1 replaces the geometric
    indel-length distribution with a Zipf one. Both default to off, which
    reproduces earlier results exactly; both should be exercised before claiming
    a result generalises — see the module docstring.
    """
    rng = np.random.default_rng(seed)
    tree = random_tree(taxa, rng, bmin, bmax)
    ctr = _Ctr()
    root = [[ctr.new(float(i)), rng.choice(list(BASES))] for i in range(root_len)]
    rates = None
    if rate_alpha and rate_alpha > 0:
        rates = {i: float(rng.gamma(rate_alpha, 1.0 / rate_alpha)) for i, _ in root}
    leaves: Dict[str, list] = {}

    def rec(node, sites):
        if not node["children"]:
            leaves[node["name"]] = sites
            return
        for child, bl in node["children"]:
            kids = _evolve(sites, bl, rng, ctr, indel_rate, rates, indel_zipf)
            if rates is not None:
                for i, _ in kids:                  # newly inserted sites need a rate
                    rates.setdefault(i, float(rng.gamma(rate_alpha, 1.0 / rate_alpha)))
            rec(child, kids)

    rec(tree, root)
    names = sorted(leaves, key=lambda n: int(n[1:]))
    allids = sorted({i for s in leaves.values() for i, _ in s}, key=lambda i: (ctr.pos[i], i))
    colidx = {i: c for c, i in enumerate(allids)}
    rows = []
    for nm in names:
        row = ["-"] * len(allids)
        for i, b in leaves[nm]:
            row[colidx[i]] = b
        rows.append("".join(row))
    seqs = [(nm, "".join(b for _, b in leaves[nm])) for nm in names]
    return dict(names=names, seqs=seqs, true_rows=rows, tree=tree)
