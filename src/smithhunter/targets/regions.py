"""Target regions (3'UTR, 5'UTR, CDS, whole transcript) from a genome + GFF3 or from
a transcripts FASTA.

Genome + GFF3: the exons of each mRNA are joined in transcript order (reverse
complemented on the minus strand). When an mRNA has no exon features, its CDS and UTR
features are used as exons. Regions come from the CDS span in transcript coordinates:
5'UTR before it, 3'UTR after it, so annotations without UTR features work too. Every
region keeps its exon blocks, to convert positions back to the genome.

Transcripts FASTA: an optional table (transcript, region, start, end; 1-based,
inclusive) gives the regions; without it, each transcript is one ``transcript`` region.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from urllib.parse import unquote

from ..seqio import open_text, read_fasta, read_tsv, write_tsv

REGION_TYPES = ("three_prime_UTR", "five_prime_UTR", "CDS", "transcript")
TRANSCRIPT_TYPES = {"mRNA", "transcript"}
COMPLEMENT = str.maketrans("ACGTUN", "TGCAAN")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


@dataclass
class Transcript:
    id: str
    gene: str
    contig: str = ""
    strand: str = "+"
    exons: list[tuple[int, int]] = field(default_factory=list)  # genomic, 0-based half-open
    cds: list[tuple[int, int]] = field(default_factory=list)
    utr: list[tuple[int, int]] = field(default_factory=list)


@dataclass
class Region:
    id: str
    transcript: str
    gene: str
    region: str
    sequence: str
    contig: str = ""
    strand: str = ""
    # genomic blocks (0-based half-open) in transcript order, covering the region
    blocks: list[tuple[int, int]] = field(default_factory=list)

    def genomic(self, pos: int) -> int | None:
        """Genomic 0-based coordinate of region position ``pos`` (0-based)."""
        offset = pos
        for s, e in self.blocks:
            n = e - s
            if offset < n:
                return s + offset if self.strand == "+" else e - 1 - offset
            offset -= n
        return None


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def parse_gff(path: str) -> list[Transcript]:
    transcripts: dict[str, Transcript] = {}
    parents: dict[str, str] = {}  # mRNA id -> gene name
    genes: dict[str, str] = {}
    children: list[tuple[str, str, int, int]] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = {k: unquote(v) for k, v in
                     (kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)}
            ftype, start, end = f[2], int(f[3]) - 1, int(f[4])
            if ftype == "gene":
                genes[attrs.get("ID", "")] = attrs.get("Name") or attrs.get("ID", "")
            elif ftype in TRANSCRIPT_TYPES:
                tid = attrs.get("ID") or f"{f[0]}:{start}-{end}"
                transcripts[tid] = Transcript(tid, "", f[0], f[6])
                parents[tid] = attrs.get("Parent", "")
                transcripts[tid].gene = attrs.get("Name", "")
            elif ftype in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
                for parent in attrs.get("Parent", "").split(","):
                    children.append((parent, ftype, start, end))
    for parent, ftype, start, end in children:
        t = transcripts.get(parent)
        if t is None:
            continue
        {"exon": t.exons, "CDS": t.cds}.get(ftype, t.utr).append((start, end))
    for tid, t in transcripts.items():
        gene_id = parents.get(tid, "")
        t.gene = genes.get(gene_id) or gene_id or t.gene or tid
        t.exons = merge(t.exons or (t.cds + t.utr))
    return [t for t in transcripts.values() if t.exons]


def to_transcript(t: Transcript, pos: int) -> int:
    """Transcript coordinate (0-based) of genomic position ``pos`` lying in an exon."""
    offset = 0
    ordered = t.exons if t.strand == "+" else t.exons[::-1]
    for s, e in ordered:
        if s <= pos < e:
            return offset + (pos - s if t.strand == "+" else e - 1 - pos)
        offset += e - s
    raise ValueError(f"{t.id}: position {pos} is not exonic")


def blocks_for(t: Transcript, start: int, end: int) -> list[tuple[int, int]]:
    """Genomic blocks, in transcript order, of transcript interval [start, end)."""
    out, offset = [], 0
    ordered = t.exons if t.strand == "+" else t.exons[::-1]
    for s, e in ordered:
        n = e - s
        a, b = max(start, offset), min(end, offset + n)
        if a < b:
            if t.strand == "+":
                out.append((s + a - offset, s + b - offset))
            else:
                out.append((e - (b - offset), e - (a - offset)))
        offset += n
    return out


def regions_from_genome(genome_fasta: str, gff: str, wanted: list[str],
                        whole_if_missing: bool = False) -> list[Region]:
    contigs = dict(read_fasta(genome_fasta))
    regions = []
    for t in parse_gff(gff):
        if t.contig not in contigs:
            continue
        seq = "".join(contigs[t.contig][s:e] for s, e in t.exons).upper()
        if t.strand == "-":
            seq = revcomp(seq)
        length = len(seq)
        spans: dict[str, tuple[int, int]] = {"transcript": (0, length)}
        cds = merge(c for c in t.cds if any(s <= c[0] < e for s, e in t.exons))
        if cds:
            ends = [to_transcript(t, cds[0][0]), to_transcript(t, cds[-1][1] - 1)]
            cs, ce = min(ends), max(ends) + 1
            spans["CDS"] = (cs, ce)
            spans["five_prime_UTR"] = (0, cs)
            spans["three_prime_UTR"] = (ce, length)
        found = False
        for name in wanted:
            a, b = spans.get(name, (0, 0))
            if b > a:
                found = True
                regions.append(Region(f"{t.id}|{name}", t.id, t.gene, name, seq[a:b],
                                      t.contig, t.strand, blocks_for(t, a, b)))
        if not found and whole_if_missing and "transcript" not in wanted:
            regions.append(Region(f"{t.id}|transcript", t.id, t.gene, "transcript", seq,
                                  t.contig, t.strand, blocks_for(t, 0, length)))
    return regions


def regions_from_transcripts(fasta: str, table: str = "", wanted=None,
                             whole_if_missing: bool = False) -> list[Region]:
    seqs = {name: seq.upper().replace("U", "T") for name, seq in read_fasta(fasta)}
    if not table:
        return [Region(f"{n}|transcript", n, n, "transcript", s) for n, s in seqs.items()]
    by_transcript = defaultdict(list)
    for row in read_tsv(table):
        by_transcript[row["transcript"]].append(row)
    regions = []
    for name, seq in seqs.items():
        found = False
        for row in by_transcript.get(name, []):
            if wanted and row["region"] not in wanted:
                continue
            a, b = int(row["start"]) - 1, int(row["end"])
            if not 0 <= a < b <= len(seq):
                raise ValueError(f"{table}: region {row} outside transcript {name}")
            regions.append(Region(f"{name}|{row['region']}", name, row.get("gene") or name,
                                  row["region"], seq[a:b]))
            found = True
        if not found and whole_if_missing:
            regions.append(Region(f"{name}|transcript", name, name, "transcript", seq))
    return regions


REGION_COLUMNS = ["region_id", "transcript", "gene", "region", "length", "contig", "strand",
                  "blocks", "sequence"]


def write_regions(path: str, regions: list[Region]) -> None:
    write_tsv(path, [{"region_id": r.id, "transcript": r.transcript, "gene": r.gene,
                      "region": r.region, "length": len(r.sequence), "contig": r.contig,
                      "strand": r.strand,
                      "blocks": ",".join(f"{s + 1}-{e}" for s, e in r.blocks),
                      "sequence": r.sequence} for r in regions], REGION_COLUMNS)


def read_regions(path: str) -> list[Region]:
    out = []
    for row in read_tsv(path):
        blocks = []
        for b in filter(None, row["blocks"].split(",")):
            s, e = b.split("-")
            blocks.append((int(s) - 1, int(e)))
        out.append(Region(row["region_id"], row["transcript"], row["gene"], row["region"],
                          row["sequence"], row["contig"], row["strand"], blocks))
    return out
