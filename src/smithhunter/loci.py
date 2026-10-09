"""Define small RNA loci on focus genomes, quantify them and score their ends.

A locus is a run of overlapping alignments on one strand of one contig (reads
closer than ``merge_gap`` bases are joined). On circular contigs a locus may cross
the origin. This replaces the sequence-identity clustering of SmithHunter v0,
which could join reads from different positions and ignored terminal gaps.
"""

from __future__ import annotations

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


@dataclass
class Locus:
    genome: str
    contig: str
    strand: str
    start: int
    end: int
    members: list[Member] = field(default_factory=list)
    wraps: bool = False


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
        shifted = [Member(m.read_id, m.start + length, m.end + length, m.weight)
                   for m in first.members]
        last.members.extend(shifted)
        last.end = max(last.end, first.end + length)
    for locus in loci:
        locus.wraps = circular and locus.end > length
    return loci


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
             merge_gap=0, min_rpm=5.0, min_samples=None, n_thre=0.5, penalty=0.1,
             min_five=0.5, min_three=0.0):
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
            grouped[(genome, contig, strand)].append(Member(rid, int(start), int(end), w))
            totals = genome_totals[genome]
            for i, c in enumerate(row):
                totals[i] += w * c

    annotations = {name: Annotation.from_gff(g["annotation"])
                   for name, g in genomes.items() if g["annotation"]}

    loci: list[Locus] = []
    for (genome, contig, strand), members in sorted(grouped.items()):
        info = contigs[(genome, contig)]
        loci.extend(build_loci(genome, contig, strand, members, merge_gap,
                               int(info["length"]), info["circular"] == "1"))
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
        n_pass = sum(r >= min_rpm for r in rpm)

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
        ends_pass = five_score > min_five and three_score >= min_three
        row = {
            "locus_id": locus_id, "genome": locus.genome, "contig": locus.contig,
            "start": locus.start + 1, "end": locus.end - length if locus.wraps else locus.end,
            "strand": locus.strand, "span": locus.end - locus.start, "wraps_origin": int(locus.wraps),
            "rep_id": rep_id, "rep_sequence": counts[rep_id][0], "rep_length": len(counts[rep_id][0]),
            "rep_start": rep.start % length + 1, "rep_end": (rep.end - 1) % length + 1,
            "unique_reads": len(per_read), "total_count": round(sum(per_sample), 3),
            "samples_passing_rpm": n_pass,
            "five_score": round(five_score, 4), "five_dominant_fraction": round(five_dom, 4),
            "three_score": round(three_score, 4), "three_dominant_fraction": round(three_dom, 4),
            **annotation,
            "expression_pass": int(expr_pass), "ends_pass": int(ends_pass),
            "pass": int(expr_pass and ends_pass),
        }
        for s, c, r, rg in zip(samples, per_sample, rpm, rpm_genome):
            row[f"count_{s}"] = round(c, 3)
            row[f"rpm_{s}"] = round(r, 3)
            row[f"rpm_genome_{s}"] = round(rg, 3)
        rows.append(row)
    return samples, rows, reads


LOCUS_COLUMNS = [
    "locus_id", "genome", "contig", "start", "end", "strand", "span", "wraps_origin",
    "rep_id", "rep_sequence", "rep_length", "rep_start", "rep_end", "unique_reads",
    "total_count", "samples_passing_rpm", "five_score", "five_dominant_fraction",
    "three_score", "three_dominant_fraction", "annotation_class", "annotation_types",
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
