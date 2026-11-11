"""Target regions (3'UTR, 5'UTR, CDS, whole transcript) from a genome + GFF3 or from
a transcripts FASTA.

Genome + GFF3: the exons of each mRNA are joined in transcript order (reverse
complemented on the minus strand). When an mRNA has no exon features, its CDS and UTR
features are used as exons. Regions come from the CDS span in transcript coordinates:
5'UTR before it, 3'UTR after it, so annotations without UTR features work too. Every
region keeps its exon blocks, to convert positions back to the genome.

Transcripts FASTA: an optional table (transcript, region, start, end; 1-based,
inclusive) gives the regions; without it, each transcript is one ``transcript`` region.

Bacterial mode (``bacterial_windows``): for every annotated CDS a window ``[-U, +D]`` around
the start codon (region ``start_window``), on the CDS strand, clipped at the ends of linear
contigs or wrapped around circular ones. Coordinates relative to the start codon have no
position 0: the A of the start codon is +1 and the base before it -1. Such regions carry
``anchor``, the 0-based index of the start codon's first base in the window sequence.
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
    # bacterial start-codon windows only (None elsewhere): 0-based index of the first base of
    # the start codon in ``sequence``; nt of the part upstream of it that lie in another CDS of
    # the same strand (overlapping operon genes); nt of the window that belong to this CDS
    anchor: int | None = None
    upstream_cds: int | None = None
    cds_nt: int | None = None

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


ANCHOR_COLUMNS = ["anchor", "upstream_cds", "cds_nt"]


def write_regions(path: str, regions: list[Region]) -> None:
    """regions.tsv. The anchor columns (bacterial windows) are added only when used, so the
    table of the other modes keeps its layout. ``anchor`` is written 1-based."""
    with_anchor = any(r.anchor is not None for r in regions)
    rows = []
    for r in regions:
        row = {"region_id": r.id, "transcript": r.transcript, "gene": r.gene,
               "region": r.region, "length": len(r.sequence), "contig": r.contig,
               "strand": r.strand, "blocks": ",".join(f"{s + 1}-{e}" for s, e in r.blocks),
               "sequence": r.sequence}
        if with_anchor:
            row["anchor"] = "" if r.anchor is None else r.anchor + 1
            row["upstream_cds"] = "" if r.upstream_cds is None else r.upstream_cds
            row["cds_nt"] = "" if r.cds_nt is None else r.cds_nt
        rows.append(row)
    write_tsv(path, rows, REGION_COLUMNS + (ANCHOR_COLUMNS if with_anchor else []))


def read_regions(path: str) -> list[Region]:
    out = []
    for row in read_tsv(path):
        blocks = []
        for b in filter(None, row["blocks"].split(",")):
            s, e = b.split("-")
            blocks.append((int(s) - 1, int(e)))
        def opt(key, shift=0):
            v = row.get(key)
            return int(v) + shift if v not in (None, "") else None

        out.append(Region(row["region_id"], row["transcript"], row["gene"], row["region"],
                          row["sequence"], row["contig"], row["strand"], blocks,
                          opt("anchor", -1), opt("upstream_cds"), opt("cds_nt")))
    return out


# ---- bacterial mode: windows around start codons ---------------------------------------

@dataclass
class Cds:
    id: str                      # transcript column: locus_tag when present, else the CDS ID
    gene: str                    # gene symbol when annotated, else the locus tag
    contig: str
    strand: str
    segments: list[tuple[int, int]] = field(default_factory=list)  # 0-based half-open
    partial_start: bool = False  # start codon not annotated (5' partial)
    pseudo: bool = False

    @property
    def start(self) -> int:
        """Genomic 0-based position of the first base of the start codon."""
        if self.strand == "+":
            return min(s for s, _ in self.segments)
        return max(e for _, e in self.segments) - 1

    @property
    def length(self) -> int:
        return sum(e - s for s, e in self.segments)


def _attrs(text: str) -> dict[str, str]:
    return {k: unquote(v) for k, v in (kv.split("=", 1) for kv in text.split(";") if "=" in kv)}


def parse_cds(path: str) -> list[Cds]:
    """CDS features of a GFF3 (NCBI, Prokka, Bakta or mRNA/CDS style). CDS lines with the same
    ID form one CDS (e.g. a ribosomal frameshift); its start codon is at the 5' end of the first
    segment in strand order. The gene symbol comes from the CDS ``gene`` attribute, else from the
    parent gene's ``gene``/``Name``, else the locus tag. CDS without a start codon
    (``start_range``, or phase 1/2 on the first segment) and ``pseudo=true`` are flagged."""
    names: dict[str, str] = {}       # feature ID -> symbol of the gene / mRNA features
    parents: dict[str, str] = {}
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = _attrs(f[8])
            if f[2] == "CDS":
                rows.append((f, a))
            elif a.get("ID"):
                names[a["ID"]] = a.get("gene") or a.get("Name") or a.get("locus_tag") or ""
                parents[a["ID"]] = a.get("Parent", "")
    cds: dict[tuple, Cds] = {}
    first_phase: dict[tuple, tuple[int, str]] = {}
    for n, (f, a) in enumerate(rows):
        start, end, strand = int(f[3]) - 1, int(f[4]), f[6]
        if strand not in ("+", "-"):
            continue
        key = (f[0], strand, a["ID"]) if a.get("ID") else (f[0], strand, f"line{n}")
        sym = a.get("gene", "")
        parent, hops = a.get("Parent", "").split(",")[0], 0
        while not sym and parent and hops < 3:     # CDS -> mRNA -> gene
            sym = names.get(parent, "")
            parent, hops = parents.get(parent, ""), hops + 1
        locus = a.get("locus_tag", "")
        ident = locus or a.get("ID") or f"{f[0]}:{start + 1}-{end}({strand})"
        c = cds.get(key)
        if c is None:
            c = cds[key] = Cds(ident, sym or locus or a.get("Name") or ident, f[0], strand)
        c.segments.append((start, end))
        c.pseudo |= a.get("pseudo", "").lower() == "true"
        c.partial_start |= "start_range" in a
        edge = start if strand == "+" else -end   # smallest = first segment in strand order
        if key not in first_phase or edge < first_phase[key][0]:
            first_phase[key] = (edge, f[7] if f[7] in ("0", "1", "2") else "0")
    for key, c in cds.items():
        c.segments.sort()
        if first_phase[key][1] != "0":
            c.partial_start = True
    return list(cds.values())


def segments_of(lo: int, hi: int, length: int, circular: bool) -> list[tuple[int, int]]:
    """Plus-strand blocks (0-based half-open, in increasing order along the window) covering the
    possibly out-of-range interval [lo, hi) of a contig of ``length`` bases: cut at the contig
    ends when linear, continued at the other end when circular (more than once if the interval
    is longer than the contig)."""
    if not circular:
        lo, hi = max(lo, 0), min(hi, length)
        return [(lo, hi)] if lo < hi else []
    out, pos = [], lo
    while pos < hi:
        k = pos % length
        n = min(length - k, hi - pos)
        out.append((k, k + n))
        pos += n
    return out


def _upstream_coding(blocks, entry, own) -> int:
    """Number of bases of ``blocks`` covered by the union of the CDS in ``entry`` (their
    intervals sorted by start with the start list and the longest length) other than CDS
    number ``own``."""
    from bisect import bisect_left, bisect_right
    intervals, starts, max_len = entry
    covered = 0
    for a, b in blocks:
        i, j = bisect_left(starts, a - max_len), bisect_right(starts, b)
        clipped = [(max(s, a), min(e, b)) for s, e, k in intervals[i:j]
                   if k != own and e > a and s < b]
        covered += sum(e - s for s, e in merge(clipped))
    return covered


def bacterial_windows(genome_fasta: str, gff: str, window: tuple[int, int] = (-200, 100),
                      circular=False, min_length: int = 30, stats: dict | None = None
                      ) -> list[Region]:
    """One region ``start_window`` per CDS: ``window[0]`` .. ``window[1]`` relative to the start
    codon (-200..+100 = 200 nt upstream plus the first 100 nt counted from the A of the start
    codon), on the strand of the CDS.

    ``circular``: True (every contig), False, or a collection of contig names that are circular.
    A window reaching beyond the end of a circular contig continues at its other end; on a linear
    contig it is cut there (the start codon keeps its place, see ``anchor``), and dropped if
    shorter than ``min_length``. CDS without an annotated start codon and pseudogenes are
    skipped; ``stats`` (a dict) receives the counts.

    Operons: windows of neighbouring genes overlap, and the upstream part of a window may lie
    inside the preceding CDS on the same strand. Nothing is removed (sRNAs also pair inside
    upstream coding sequences, and the annotation of the 5'UTR is not known); ``upstream_cds``
    reports how many upstream nt are coding sequence of another gene."""
    up, down = -window[0], window[1]
    if window[0] >= 0 or down < 1:
        raise ValueError(f"window must be [-U, +D] with U >= 1 and D >= 1, got {list(window)}")
    contigs = {n: s.upper() for n, s in read_fasta(genome_fasta)}
    if circular is True:
        circ = set(contigs)
    elif not circular:
        circ = set()
    else:
        circ = set(circular)
    stats = stats if stats is not None else {}
    stats.update(cds=0, windows=0, partial=0, pseudo=0, no_contig=0, too_short=0, clipped=0,
                 wrapped=0)
    features = parse_cds(gff)
    coding: dict[tuple, list] = defaultdict(list)    # (contig, strand) -> (start, end, CDS)
    for i, c in enumerate(features):
        coding[(c.contig, c.strand)].extend((s, e, i) for s, e in c.segments)
    index = {}
    for k, v in coding.items():
        v.sort()
        index[k] = (v, [s for s, _, _ in v], max(e - s for s, e, _ in v))
    used: dict[str, int] = {}
    regions = []
    for ci, c in sorted(enumerate(features), key=lambda t: (t[1].contig, t[1].start)):
        stats["cds"] += 1
        if c.contig not in contigs:
            stats["no_contig"] += 1
            continue
        if c.pseudo:
            stats["pseudo"] += 1
            continue
        if c.partial_start:
            stats["partial"] += 1
            continue
        seq, length, is_circ = contigs[c.contig], len(contigs[c.contig]), c.contig in circ
        g = c.start
        if c.strand == "+":
            lo, hi, ulo, uhi = g - up, g + down, g - up, g
        else:
            lo, hi, ulo, uhi = g - down + 1, g + up + 1, g + 1, g + up + 1
        blocks = segments_of(lo, hi, length, is_circ)
        up_blocks = segments_of(ulo, uhi, length, is_circ)
        window_seq = "".join(seq[a:b] for a, b in blocks)
        if c.strand == "-":
            window_seq, blocks = revcomp(window_seq), blocks[::-1]
        n_up = sum(b - a for a, b in up_blocks)
        if len(window_seq) < min_length:
            stats["too_short"] += 1
            continue
        stats["clipped"] += len(window_seq) < up + down
        stats["wrapped"] += len(blocks) > 1
        n = used[c.id] = used.get(c.id, 0) + 1
        tid = c.id if n == 1 else f"{c.id}#{n}"
        regions.append(Region(f"{tid}|start_window", tid, c.gene, "start_window", window_seq,
                              c.contig, c.strand, blocks, anchor=n_up,
                              upstream_cds=_upstream_coding(up_blocks, index[(c.contig, c.strand)], ci),
                              cds_nt=min(c.length, len(window_seq) - n_up)))
        stats["windows"] += 1
    return regions
