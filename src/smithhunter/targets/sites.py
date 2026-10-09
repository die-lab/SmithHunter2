"""Seed site scan with canonical site classes.

Guide positions are numbered from the 5' end of the small RNA; target position t1 is
opposite guide position 1, so it lies just 3' of the seed match on the target. With
the seed at guide 2-8 (Bartel 2018; McGeary et al. 2019):

    8mer         match to guide 2-8 + A at t1
    7mer-m8      match to guide 2-8
    7mer-A1      match to guide 2-7 + A at t1
    6mer         match to guide 2-7
    offset-6mer  match to guide 3-8

Every site gets only its best class. With another seed range the class is
``seed<start>-<end>`` (a perfect match to that range). Context features follow Grimson
et al. 2007: distance from the region ends, AU content of the flanks, 3' supplementary
pairing of guide 13-16.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from .regions import revcomp

CANONICAL = ("8mer", "7mer-m8", "7mer-A1", "6mer", "offset-6mer")
FLANK_AU = 30
NEAR_STOP = 15


@dataclass
class Site:
    start: int  # 0-based, in the region
    end: int
    t1: int
    site_class: str


def occurrences(text: str, pattern: str):
    i = text.find(pattern)
    while i != -1:
        yield i
        i = text.find(pattern, i + 1)


def site_classes(seed: tuple[int, int]) -> tuple[str, ...]:
    return CANONICAL if tuple(seed) == (2, 8) else (f"seed{seed[0]}-{seed[1]}",)


def scan(guide: str, target: str, seed: tuple[int, int] = (2, 8),
         classes: tuple[str, ...] | None = None) -> list[Site]:
    """Sites of ``guide`` (DNA alphabet, 5'->3') in ``target`` (sense, 5'->3')."""
    guide = guide.upper().replace("U", "T")
    a, b = seed
    if tuple(seed) != (2, 8):
        pattern = revcomp(guide[a - 1:b])
        name = f"seed{a}-{b}"
        if len(pattern) != b - a + 1:
            return []
        return [Site(i, i + len(pattern), i + len(pattern) + a - 1, name)
                for i in occurrences(target, pattern)]
    if len(guide) < 8:
        return []
    allowed = set(classes or CANONICAL)
    h28, h27, h38 = revcomp(guide[1:8]), revcomp(guide[1:7]), revcomp(guide[2:8])
    anchors = {i + 6 for i in occurrences(target, h27)} | \
              {i + 7 for i in occurrences(target, h38)}
    sites = []
    n = len(target)
    for t1 in sorted(anchors):
        a1 = t1 < n and target[t1] == "A"
        m8 = t1 >= 7 and target[t1 - 7:t1] == h28
        m27 = target[t1 - 6:t1] == h27
        m38 = t1 >= 7 and target[t1 - 7:t1 - 1] == h38
        for cls, ok, span in (("8mer", m8 and a1, (t1 - 7, t1 + 1)),
                              ("7mer-m8", m8, (t1 - 7, t1)),
                              ("7mer-A1", m27 and a1, (t1 - 6, t1 + 1)),
                              ("6mer", m27, (t1 - 6, t1)),
                              ("offset-6mer", m38, (t1 - 7, t1 - 1))):
            if ok:
                if cls in allowed:
                    sites.append(Site(span[0], span[1], t1, cls))
                break
    return sites


def supplementary_pairing(guide: str, target: str, t1: int, max_offset: int = 4) -> int:
    """Longest run of Watson-Crick pairs between guide 13-16 and the target upstream of
    the seed, allowing the target to loop out up to ``max_offset`` bases."""
    guide = guide.upper().replace("U", "T")
    if len(guide) < 16:
        return 0
    pair = {"A": "T", "C": "G", "G": "C", "T": "A"}
    best = 0
    for offset in range(max_offset + 1):
        run = 0
        for k in range(13, 17):
            p = t1 - (k - 1) - offset
            ok = 0 <= p < len(target) and pair.get(guide[k - 1]) == target[p]
            run = run + 1 if ok else 0
            best = max(best, run)
    return best


def context(guide: str, target: str, site: Site, region_type: str) -> dict:
    lo, hi = max(0, site.start - FLANK_AU), min(len(target), site.end + FLANK_AU)
    flanks = target[lo:site.start] + target[site.end:hi]
    au = (flanks.count("A") + flanks.count("T")) / len(flanks) if flanks else 0.0
    dist5, dist3 = site.start, len(target) - site.end
    return {
        "dist_5p": dist5, "dist_3p": dist3, "au_flanks": round(au, 3),
        "supp_pairing": supplementary_pairing(guide, target, site.t1),
        "near_stop": int(region_type == "three_prime_UTR" and dist5 < NEAR_STOP),
    }


def count_sites(guide: str, targets, seed=(2, 8), classes=None) -> Counter:
    """Sites per class over an iterable of target sequences."""
    counts: Counter = Counter()
    for seq in targets:
        for s in scan(guide, seq, seed, classes):
            counts[s.site_class] += 1
    return counts
