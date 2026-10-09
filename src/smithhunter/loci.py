"""Define small RNA loci on focus genomes, quantify them and score their ends.

Loci are built in two steps.

1. Blocks: runs of overlapping alignments on one strand of one contig (reads closer
   than ``merge_gap`` bases are joined). On circular contigs a block may cross the origin.
2. Peaks: inside a block, small RNAs show up as positions where many reads start (their
   5' end is precise), while degradation fragments start anywhere. The 5' positions are
   taken from the most to the least abundant; each one, with the reads starting within
   ``peak_window`` bases of it, becomes a ``peak`` locus when that read count is
   significantly above the local background (Poisson test, ``peak_pvalue``). The
   background is the mean number of reads starting at each position within
   ``background_flank`` bases, leaving out the positions of stronger peaks already called. Reads claimed by no peak are chained again into
   ``background`` loci, which are reported but never pass the filters. A weak peak close
   to a much stronger one (``satellite_distance``, ``satellite_fraction``) is taken as a
   5' variant of it and joins its locus. The ``isolation`` of a peak locus is the share
   of the reads starting within ``satellite_distance`` bases of the peak that belong to
   it: low values mean a region with scattered 5' ends, or a minor neighbour of a peak.

A sample supports a peak locus when its reads are significantly above that sample's own
local background and, after subtracting the background, reach ``min_count`` and ``min_rpm``.

Without step 2, degradation fragments chain every small RNA of an expressed transcript
into one locus as long as the transcript. This replaces the sequence-identity
clustering of SmithHunter v0, which could join reads from different positions.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field

from .annotate import Annotation, describe
from .ends import end_score
from .seqio import open_text, read_tsv


@dataclass
class Member:
    read_id: str
    start: int
    end: int
    weight: float
    count: float = 0.0  # weight x reads in all samples
    counts: tuple = ()  # weight x reads per sample


@dataclass
class Locus:
    genome: str
    contig: str
    strand: str
    start: int
    end: int
    members: list[Member] = field(default_factory=list)
    wraps: bool = False
    kind: str = "peak"
    background: list[float] = field(default_factory=list)  # expected reads per sample
    isolation: float = 1.0  # share of the reads starting near the peak that are in the locus


def merge_members(members: list[Member], gap: int) -> list[tuple[int, int, list[Member]]]:
    blocks: list[list] = []
    for m in sorted(members, key=lambda m: (m.start, m.end)):
        if blocks and m.start <= blocks[-1][1] + gap:
            blocks[-1][1] = max(blocks[-1][1], m.end)
            blocks[-1][2].append(m)
        else:
            blocks.append([m.start, m.end, [m]])
    return [tuple(b) for b in blocks]


def build_loci(genome, contig, strand, members, gap, length, circular) -> list[Locus]:
    blocks = merge_members(members, gap)
    loci = [Locus(genome, contig, strand, s, e, ms) for s, e, ms in blocks]
    if circular and len(loci) > 1 and loci[-1].end + gap - length >= loci[0].start:
        first, last = loci.pop(0), loci[-1]
        shifted = [Member(m.read_id, m.start + length, m.end + length, m.weight, m.count, m.counts)
                   for m in first.members]
        last.members.extend(shifted)
        last.end = max(last.end, first.end + length)
    for locus in loci:
        locus.wraps = circular and locus.end > length
    return loci


def five_prime(m: Member, strand: str) -> int:
    return m.start if strand == "+" else m.end - 1


def poisson_significant(k: float, background: float, pvalue: float) -> bool:
    """True when P(X >= k) < pvalue for X ~ Poisson(background)."""
    k = math.floor(k + 1e-9)
    if k <= 0:
        return False
    if background <= 0:
        return True
    if k > background + 10 * math.sqrt(background) + 10:
        return True
    term = cdf = math.exp(-background)
    for i in range(1, k):
        term *= background / i
        cdf += term
    return 1.0 - cdf < pvalue


def split_peaks(block: Locus, window: int = 2, pvalue: float = 0.001, flank: int = 50,
                gap: int = 0, satellite_distance: int = 10,
                satellite_fraction: float = 0.1) -> list[Locus]:
    """Split one block of overlapping reads into peak loci and background loci.

    A significant peak within ``satellite_distance`` bases of a stronger peak, with less
    than ``satellite_fraction`` of its reads, is a satellite: its reads join that peak.
    """
    profile: dict[int, float] = defaultdict(float)
    per_sample: dict[int, list[float]] = {}
    n = max((len(m.counts) for m in block.members), default=0)
    for m in block.members:
        pos = five_prime(m, block.strand)
        profile[pos] += m.count
        acc = per_sample.setdefault(pos, [0.0] * n)
        for i, c in enumerate(m.counts):
            acc[i] += c
    claimed: dict[int, int] = {}
    peaks: list[int] = []
    peak_reads: list[float] = []
    backgrounds: list[list[float]] = []
    for pos in sorted(profile, key=lambda p: (-profile[p], p)):
        if pos in claimed:
            continue
        near = [q for q in range(pos - window, pos + window + 1)
                if q in profile and q not in claimed]
        reads = sum(profile[q] for q in near)
        around = [q for q in range(pos - flank, pos + flank + 1)
                  if abs(q - pos) > window and q not in claimed]
        scale = (2 * window + 1) / len(around) if around else 0.0
        background = sum(profile.get(q, 0.0) for q in around) * scale
        if not poisson_significant(reads, background, pvalue):
            continue
        sample_bg = [sum(per_sample[q][i] for q in around if q in per_sample) * scale
                     for i in range(n)]
        parent = next((i for i, p in enumerate(peaks) if abs(p - pos) <= satellite_distance
                       and reads < satellite_fraction * peak_reads[i]), None)
        if parent is None:
            parent = len(peaks)
            peaks.append(pos)
            peak_reads.append(0.0)
            backgrounds.append([0.0] * n)
        peak_reads[parent] += reads
        backgrounds[parent] = [a + b for a, b in zip(backgrounds[parent], sample_bg)]
        for q in near:
            claimed[q] = parent

    groups: list[list[Member]] = [[] for _ in peaks]
    rest: list[Member] = []
    for m in block.members:
        idx = claimed.get(five_prime(m, block.strand))
        (rest if idx is None else groups[idx]).append(m)
    loci = []
    for pos, ms, bg in zip(peaks, groups, backgrounds):
        nearby = sum(profile.get(q, 0.0) for q in range(pos - satellite_distance,
                                                       pos + satellite_distance + 1))
        own = sum(m.count for m in ms)
        loci.append(Locus(block.genome, block.contig, block.strand, min(m.start for m in ms),
                          max(m.end for m in ms), ms, background=bg,
                          isolation=own / max(own, nearby) if own else 0.0))
    loci += [Locus(block.genome, block.contig, block.strand, s, e, ms, kind="background")
             for s, e, ms in merge_members(rest, gap)]
    return loci


def fold_origin(locus: Locus, length: int, circular: bool) -> Locus:
    """Bring a locus that lies past the end of a circular contig back to the origin."""
    if circular and locus.start >= length:
        locus.members = [Member(m.read_id, m.start - length, m.end - length, m.weight, m.count,
                                m.counts)
                         for m in locus.members]
        locus.start -= length
        locus.end -= length
    locus.wraps = circular and locus.end > length
    return locus


def load_counts(path: str):
    counts: dict[str, tuple[str, list[float]]] = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[3:]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            counts[f[0]] = (f[1], [float(x) for x in f[3:]])
    return samples, counts


def discover(alignments_path, counts_path, libraries_path, genomes_path, contigs_path,
             merge_gap=0, min_rpm=5.0, min_count=1, min_samples=None, n_thre=0.5, penalty=0.1,
             min_five=0.5, min_three=0.0, peaks=True, peak_window=2, peak_pvalue=0.001,
             background_flank=50, satellite_distance=10, satellite_fraction=0.1,
             min_isolation=0.5):
    samples, counts = load_counts(counts_path)
    library = {r["sample"]: float(r["reads_kept"]) for r in read_tsv(libraries_path)}
    genomes = {g["name"]: g for g in read_tsv(genomes_path)}
    contigs = {(c["genome"], c["contig"]): c for c in read_tsv(contigs_path)}
    if min_samples is None:
        min_samples = max(1, len(samples) - 1)

    grouped: dict[tuple, list[Member]] = defaultdict(list)
    genome_totals: dict[str, list[float]] = defaultdict(lambda: [0.0] * len(samples))
    with open_text(alignments_path) as fh:
        fh.readline()
        for line in fh:
            rid, genome, contig, start, end, strand, weight, _ = line.rstrip("\n").split("\t")
            g = genomes[genome]
            seq, row = counts[rid]
            if not int(g["min_length"]) <= len(seq) <= int(g["max_length"]):
                continue
            w = float(weight)
            grouped[(genome, contig, strand)].append(
                Member(rid, int(start), int(end), w, w * sum(row), tuple(w * c for c in row)))
            totals = genome_totals[genome]
            for i, c in enumerate(row):
                totals[i] += w * c

    annotations = {name: Annotation.from_gff(g["annotation"])
                   for name, g in genomes.items() if g["annotation"]}

    loci: list[Locus] = []
    for (genome, contig, strand), members in sorted(grouped.items()):
        info = contigs[(genome, contig)]
        length, circular = int(info["length"]), info["circular"] == "1"
        for block in build_loci(genome, contig, strand, members, merge_gap, length, circular):
            if not peaks:
                loci.append(block)
                continue
            for locus in split_peaks(block, peak_window, peak_pvalue, background_flank, merge_gap,
                                     satellite_distance, satellite_fraction):
                loci.append(fold_origin(locus, length, circular))
    loci.sort(key=lambda l: (l.genome, l.contig, l.start, l.strand))

    rows, reads = [], []
    serial: dict[str, int] = defaultdict(int)
    for locus in loci:
        serial[locus.genome] += 1
        locus_id = f"{locus.genome}_{serial[locus.genome]:05d}"
        length = int(contigs[(locus.genome, locus.contig)]["length"])
        per_sample = [0.0] * len(samples)
        per_read: dict[str, float] = defaultdict(float)
        five: dict[int, float] = defaultdict(float)
        three: dict[int, float] = defaultdict(float)
        for m in locus.members:
            seq, row = counts[m.read_id]
            for i, c in enumerate(row):
                per_sample[i] += m.weight * c
            w = m.weight * sum(row)
            per_read[m.read_id] += w
            five_pos, three_pos = (m.start, m.end - 1) if locus.strand == "+" else (m.end - 1, m.start)
            five[five_pos] += w
            three[three_pos] += w
            reads.append({"locus_id": locus_id, "read_id": m.read_id, "sequence": seq,
                          "start": m.start + 1, "end": m.end, "weight": round(m.weight, 6),
                          "count": round(w, 3)})
        rep_id = max(per_read, key=lambda r: (per_read[r], r))
        rep = max((m for m in locus.members if m.read_id == rep_id), key=lambda m: m.weight)
        rpm = [c / library[s] * 1e6 if library.get(s) else 0.0 for c, s in zip(per_sample, samples)]
        totals = genome_totals[locus.genome]
        rpm_genome = [c / t * 1e6 if t else 0.0 for c, t in zip(per_sample, totals)]
        five_score, five_dom = end_score(five, n_thre, penalty)
        three_score, three_dom = end_score(three, n_thre, penalty)
        # A sample counts when the locus is expressed above the local background in it.
        background = locus.background or [0.0] * len(samples)
        n_pass = 0
        for c, b, s in zip(per_sample, background, samples):
            net = c - b
            net_rpm = net / library[s] * 1e6 if library.get(s) else 0.0
            if net >= min_count and net_rpm >= min_rpm and (
                    b <= 0 or poisson_significant(c, b, peak_pvalue)):
                n_pass += 1

        ann = annotations.get(locus.genome)
        if ann:
            feats = ann.overlaps(locus.contig, locus.start, min(locus.end, length))
            if locus.wraps:
                feats += ann.overlaps(locus.contig, 0, locus.end - length)
            annotation = describe(feats, locus.strand)
        else:
            annotation = {"annotation_class": "NA", "annotation_types": "", "annotation_names": "",
                          "annotation_orientation": ""}

        expr_pass = n_pass >= min_samples
        ends_pass = (five_score > min_five and three_score >= min_three
                     and locus.isolation >= min_isolation)
        candidate = locus.kind == "peak"
        row = {
            "locus_id": locus_id, "genome": locus.genome, "contig": locus.contig,
            "start": locus.start + 1, "end": locus.end - length if locus.wraps else locus.end,
            "strand": locus.strand, "span": locus.end - locus.start, "wraps_origin": int(locus.wraps),
            "locus_type": locus.kind,
            "rep_id": rep_id, "rep_sequence": counts[rep_id][0], "rep_length": len(counts[rep_id][0]),
            "rep_start": rep.start % length + 1, "rep_end": (rep.end - 1) % length + 1,
            "unique_reads": len(per_read), "total_count": round(sum(per_sample), 3),
            "background_count": round(sum(background), 3),
            "samples_passing_rpm": n_pass,
            "five_score": round(five_score, 4), "five_dominant_fraction": round(five_dom, 4),
            "three_score": round(three_score, 4), "three_dominant_fraction": round(three_dom, 4),
            "isolation": round(locus.isolation, 4),
            **annotation,
            "expression_pass": int(expr_pass), "ends_pass": int(ends_pass),
            "pass": int(candidate and expr_pass and ends_pass),
        }
        for s, c, r, rg in zip(samples, per_sample, rpm, rpm_genome):
            row[f"count_{s}"] = round(c, 3)
            row[f"rpm_{s}"] = round(r, 3)
            row[f"rpm_genome_{s}"] = round(rg, 3)
        rows.append(row)
    return samples, rows, reads


LOCUS_COLUMNS = [
    "locus_id", "genome", "contig", "start", "end", "strand", "span", "wraps_origin",
    "locus_type", "rep_id", "rep_sequence", "rep_length", "rep_start", "rep_end", "unique_reads",
    "total_count", "background_count", "samples_passing_rpm", "five_score", "five_dominant_fraction",
    "three_score", "three_dominant_fraction", "isolation", "annotation_class", "annotation_types",
    "annotation_names", "annotation_orientation", "expression_pass", "ends_pass", "pass",
]
READ_COLUMNS = ["locus_id", "read_id", "sequence", "start", "end", "weight", "count"]


def locus_columns(samples: list[str]) -> list[str]:
    return (LOCUS_COLUMNS + [f"count_{s}" for s in samples] + [f"rpm_{s}" for s in samples]
            + [f"rpm_genome_{s}" for s in samples])


def write_fasta(path: str, rows: list[dict]) -> None:
    with open_text(path, "wt") as fh:
        for r in rows:
            if r["pass"]:
                fh.write(f">{r['locus_id']} {r['genome']}|{r['contig']}:{r['rep_start']}-"
                         f"{r['rep_end']}({r['strand']}) locus={r['start']}-{r['end']} "
                         f"count={r['total_count']} class={r['annotation_class']}\n"
                         f"{r['rep_sequence']}\n")
