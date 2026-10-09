"""Assign every unique read to its genome of origin.

bowtie is run with ``--best --strata -k K -m K``, so the alignments reported for a
read are all equally good (same number of mismatches). The read is then:

- ``assigned``   all best hits fall in one focus genome
- ``excluded``   all best hits fall in exclude genomes
- ``ambiguous``  best hits fall in more than one genome, at least one of them focus
                 (e.g. a NUMT, or rRNA conserved between bacteria)
- ``too_many_hits``  more than K best hits (reported by bowtie through ``--max``)
- ``unmapped``   no alignment

Ambiguous reads follow the ``ambiguous`` policy: ``discard`` drops them, ``report``
keeps them out of loci but lists them, ``split`` shares the read among all its hits.
Within one genome, multi-mapping reads are shared among their hits (``fractional``)
or dropped (``unique``).
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from itertools import groupby
from typing import Iterable, Iterator

from .seqio import open_text

POLICIES = ("discard", "report", "split")
MULTIMAP = ("fractional", "unique")
CLASSES = ("assigned", "excluded", "ambiguous", "too_many_hits", "unmapped", "multimap_dropped")

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(frozen=True)
class Hit:
    ref: str
    start: int  # 0-based, inclusive
    end: int  # 0-based, exclusive
    strand: str
    mismatches: int


@dataclass
class Contig:
    ref: str
    genome: str
    contig: str
    length: int
    circular: bool


def reference_length(cigar: str) -> int:
    return sum(int(n) for n, op in _CIGAR.findall(cigar) if op in "MDN=X")


def parse_sam(lines: Iterable[str]) -> Iterator[tuple[str, Hit]]:
    for line in lines:
        if line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        flag = int(fields[1])
        if flag & 4:
            continue
        start = int(fields[3]) - 1
        mismatches = 0
        for tag in fields[11:]:
            if tag.startswith("NM:i:"):
                mismatches = int(tag[5:])
                break
        yield fields[0], Hit(fields[2], start, start + reference_length(fields[5]),
                             "-" if flag & 16 else "+", mismatches)


def fold_circular(hit: Hit, contig: Contig) -> Hit:
    """Map an alignment that lies entirely in the appended overhang back to the origin."""
    if contig.circular and hit.start >= contig.length:
        return Hit(hit.ref, hit.start - contig.length, hit.end - contig.length,
                   hit.strand, hit.mismatches)
    return hit


def group_hits(sam_lines: Iterable[str]) -> Iterator[tuple[str, list[Hit]]]:
    """Group consecutive records by read; bowtie must run with ``--reorder``."""
    seen = set()
    for read_id, items in groupby(parse_sam(sam_lines), key=lambda x: x[0]):
        if read_id in seen:
            raise ValueError(f"read {read_id} is not contiguous in the SAM; run bowtie with --reorder")
        seen.add(read_id)
        yield read_id, [hit for _, hit in items]


def classify(hits: list[Hit], contigs: dict[str, Contig], roles: dict[str, str],
             policy: str, multimap: str):
    """Return (class, genomes, [(genome, hit, weight)]) for one read."""
    unique_hits = sorted({fold_circular(h, contigs[h.ref]) for h in hits},
                         key=lambda h: (h.ref, h.start, h.strand))
    by_genome: dict[str, list[Hit]] = defaultdict(list)
    for h in unique_hits:
        by_genome[contigs[h.ref].genome].append(h)
    genomes = sorted(by_genome)
    focus = [g for g in genomes if roles[g] == "focus"]

    if not focus:
        return "excluded", genomes, []
    if len(genomes) == 1:
        own = by_genome[focus[0]]
        if multimap == "unique" and len(own) > 1:
            return "multimap_dropped", genomes, []
        return "assigned", genomes, [(focus[0], h, 1.0 / len(own)) for h in own]
    # ambiguous between genomes
    if policy == "split":
        weight = 1.0 / len(unique_hits)
        placed = [(g, h, weight) for g in focus for h in by_genome[g]]
        return "ambiguous", genomes, placed
    return "ambiguous", genomes, []


def assign(sam_path: str, excess_ids: set[str], all_ids: list[str], contigs: dict[str, Contig],
           roles: dict[str, str], policy: str, multimap: str):
    """Yield (read_id, class, genomes, placed hits) for every unique read."""
    if policy not in POLICIES:
        raise ValueError(f"ambiguous policy must be one of {POLICIES}")
    if multimap not in MULTIMAP:
        raise ValueError(f"multimap must be one of {MULTIMAP}")
    mapped = set()
    with open_text(sam_path) as fh:
        for read_id, hits in group_hits(fh):
            mapped.add(read_id)
            cls, genomes, placed = classify(hits, contigs, roles, policy, multimap)
            yield read_id, cls, genomes, placed
    for read_id in all_ids:
        if read_id not in mapped:
            yield read_id, "too_many_hits" if read_id in excess_ids else "unmapped", [], []
