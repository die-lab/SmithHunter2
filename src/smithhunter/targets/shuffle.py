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


def shuffle_outside_seed(seq: str, seed: tuple[int, int], rng: random.Random) -> str:
    """Keep guide positions ``seed`` (and the bases 5' of it) and shuffle the rest, which
    keeps the dinucleotides of the 3' part. Real and decoy then share the same seed sites
    and differ only in what pairs beyond the seed."""
    head, tail = seq[:seed[1]], seq[seed[1]:]
    return head + dinucleotide_shuffle(tail, rng) if len(tail) > 2 else seq


DECOY_MODES = ("shuffle", "keep_seed", "match_seed_gc")


def gc_count(seq: str) -> int:
    return sum(1 for c in seq if c in "GC")


def make_decoys(candidates: list[tuple[str, str]], n: int, seed: tuple[int, int] = (2, 8),
                random_seed: int = 1, max_tries: int = 100,
                keep_seed: bool = False, mode: str = "shuffle") -> list[tuple[str, str, str]]:
    """Return (decoy_id, candidate_id, sequence) with ``n`` decoys per candidate.

    ``mode``:

    ``shuffle`` (default): dinucleotide shuffles of the whole sequence. A decoy is redrawn
    when its seed equals the seed of any real candidate, so decoys never share the targets of
    a real small RNA by construction. The null is "a random sequence of the same
    composition"; its seed differs from the real one, so seed GC content is a confounder.

    ``match_seed_gc``: as ``shuffle``, but a decoy is also redrawn until its seed has the
    same number of G+C as the real seed (the closest one after ``max_tries``), which removes
    the largest seed-composition difference between real and decoy sites.

    ``keep_seed``: only the part 3' of the seed is shuffled (the seed and its 5' side stay),
    so every decoy has exactly the real seed sites and the comparison isolates the pairing
    beyond the seed. Decoys identical to the real sequence, or to an earlier decoy of the
    same candidate, are redrawn.

    ``keep_seed=True`` is the older spelling of ``mode="keep_seed"``.
    """
    if keep_seed:
        mode = "keep_seed"
    if mode not in DECOY_MODES:
        raise ValueError(f"decoy mode must be one of {DECOY_MODES}, got {mode!r}")
    rng = random.Random(random_seed)
    real_seeds = {seed_of(s, seed) for _, s in candidates}
    out = []
    for cid, seq in candidates:
        seen = {seq}
        want_gc = gc_count(seed_of(seq, seed))
        for i in range(n):
            best, best_gap, last = None, None, None
            for _ in range(max_tries):
                if mode == "keep_seed":
                    decoy = last = shuffle_outside_seed(seq, seed, rng)
                    if decoy not in seen:
                        best = decoy
                        break
                    continue
                decoy = last = dinucleotide_shuffle(seq, rng)
                if seed_of(decoy, seed) in real_seeds:
                    continue
                if mode == "shuffle":
                    best = decoy
                    break
                gap = abs(gc_count(seed_of(decoy, seed)) - want_gc)
                if gap == 0:
                    best = decoy
                    break
                if best_gap is None or gap < best_gap:
                    best, best_gap = decoy, gap
            decoy = best if best is not None else (last if last is not None else seq)
            seen.add(decoy)
            out.append((f"{cid}__decoy{i + 1:03d}", cid, decoy))
    return out
