"""Dinucleotide-preserving shuffle (Altschul & Erickson 1985; as in uShuffle, Jiang et
al. 2008) and decoy small RNAs.

The sequence is a walk on a graph whose vertices are nucleotides and whose edges are
its dinucleotides. A random Eulerian walk with the same first base gives a sequence with
exactly the same dinucleotide counts (hence base counts and last base). The last edge
leaving each vertex is drawn so that these edges form a tree pointing to the final base,
which guarantees the walk uses every edge.
"""

from __future__ import annotations

import random
from collections import defaultdict


def dinucleotide_shuffle(seq: str, rng: random.Random) -> str:
    if len(seq) < 3:
        return seq
    edges: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    last = seq[-1]
    vertices = set(edges) - {last}
    while True:
        last_edge = {v: rng.choice(edges[v]) for v in vertices}
        if all(_reaches(v, last, last_edge) for v in vertices):
            break
    walk_edges: dict[str, list[str]] = {}
    for v, succ in edges.items():
        succ = list(succ)
        if v in last_edge:
            succ.remove(last_edge[v])
            rng.shuffle(succ)
            succ.append(last_edge[v])
        else:
            rng.shuffle(succ)
        walk_edges[v] = succ[::-1]  # pop from the end
    out = [seq[0]]
    v = seq[0]
    while walk_edges.get(v):
        v = walk_edges[v].pop()
        out.append(v)
    return "".join(out)


def _reaches(v: str, root: str, last_edge: dict[str, str]) -> bool:
    seen = set()
    while v != root:
        if v in seen:
            return False
        seen.add(v)
        v = last_edge[v]
    return True


def seed_of(seq: str, seed: tuple[int, int]) -> str:
    return seq[seed[0] - 1:seed[1]]


def make_decoys(candidates: list[tuple[str, str]], n: int, seed: tuple[int, int] = (2, 8),
                random_seed: int = 1, max_tries: int = 100) -> list[tuple[str, str, str]]:
    """Return (decoy_id, candidate_id, sequence) with ``n`` decoys per candidate.

    A decoy is redrawn when its seed equals the seed of any real candidate, so decoys
    never share the targets of a real small RNA by construction.
    """
    rng = random.Random(random_seed)
    real_seeds = {seed_of(s, seed) for _, s in candidates}
    out = []
    for cid, seq in candidates:
        for i in range(n):
            decoy = dinucleotide_shuffle(seq, rng)
            for _ in range(max_tries):
                if seed_of(decoy, seed) not in real_seeds:
                    break
                decoy = dinucleotide_shuffle(seq, rng)
            out.append((f"{cid}__decoy{i + 1:03d}", cid, decoy))
    return out
