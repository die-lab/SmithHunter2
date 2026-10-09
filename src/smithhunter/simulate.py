"""Simulate a small RNA-seq experiment with known truth, to benchmark module A.

The simulator writes synthetic genomes, a GFF3 for the mitochondrion, FASTQ files,
a ready-to-run configuration and the truth tables used by ``smithhunter evaluate``.

Genomes
- ``mito`` (focus, circular): an AT-rich genome tiled with protein genes, 22 tRNAs,
  two rRNAs and a control region, on both strands.
- ``nuclear`` (exclude): random sequence holding NUMTs (mitochondrial segments, one
  identical and one differing by a single base inside a smallRNA), a repeat with more
  copies than ``max_hits``, piRNA-like clusters and nuclear smallRNAs.
- ``bacterium`` (exclude): a contaminant with a multi-copy rRNA operon.

Library content, with defaults taken from the small RNA-seq literature
- smallRNA isoforms: 5' ends are precise (about 90% canonical, shifts of 1-2 nt),
  3' ends vary more (about 55% canonical, -3..+3 nt); about 10% of reads carry 1-2
  non-templated 3' nucleotides, mostly A and U (Morin et al. 2008 Genome Res;
  Neilsen et al. 2012 Trends Genet; Wyman et al. 2011 Genome Res; Burroughs et al. 2010).
- Replicate counts are negative binomial (dispersion 0.1), expression is log-uniform
  over two orders of magnitude.
- Sequencing errors are substitutions, about 0.2% per base rising towards the end of
  the read (Schirmer et al. 2016 BMC Bioinformatics).
- Degradation fragments of rRNAs and mRNAs, sense strand, 16-45 nt.
- tRNA-derived fragments: 5' fragments and halves, and 3' fragments carrying the
  non-templated CCA tail (Kumar et al. 2014 BMC Biol; Selitsky & Sethupathy 2015).
- piRNA-like reads (24-30 nt, 1U bias), reads from unrelated organisms (unmapped)
  and adapter dimers or very short inserts.
- Reads are built as insert + 3' adapter + random bases, cut to the read length;
  with ``--paired`` R2 is the reverse complement of the insert + the R2 adapter.

Scenarios placed on the mitochondrion, each with the outcome module A should give
(``expected`` column of ``truth/features.tsv``): ordinary smallRNAs, one inside the
most expressed rRNA, one across the origin, a sense/antisense pair, two smallRNAs
10 nt apart, two overlapping by 6 nt, one copied
identically into a NUMT (ambiguous reads: missed with ``ambiguous: report``), one in a
NUMT with one difference (assigned to mito by ``--best --strata``), one with a fixed
SNP against the reference, one with heteroplasmy in one sample, one too weak to pass
and one expressed in a single sample.
"""

from __future__ import annotations

import gzip
import math
import random
from dataclasses import dataclass, field
from pathlib import Path

from .seqio import write_tsv

ADAPTER_R1 = "TGGAATTCTCGGGTGCCAAGG"
ADAPTER_R2 = "GATCGTCGGACTGTAGAACTCTGAAC"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# Offsets of the 5' start (positive: later, shorter) and of the 3' end (positive: longer).
FIVE_SHIFT = {0: 0.90, -1: 0.035, 1: 0.035, -2: 0.015, 2: 0.015}
THREE_SHIFT = {-3: 0.03, -2: 0.07, -1: 0.15, 0: 0.55, 1: 0.13, 2: 0.05, 3: 0.02}
TAIL_RATE = 0.10
TAIL_BASES = {"A": 0.45, "T": 0.40, "C": 0.08, "G": 0.07}
DISPERSION = 0.1

# Fraction of each library taken by each background source; mitochondrial smallRNAs
# get their own expression values (RPM) on top of this.
COMPOSITION = {
    "mito_trf": 0.04, "nuclear_srna": 0.25, "nuclear_pirna": 0.10,
    "nuclear_degradation": 0.08, "nuclear_repeat": 0.02, "bacterium": 0.05,
    "foreign": 0.08, "short": 0.04,
}


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


@dataclass
class Params:
    outdir: str = "sim"
    seed: int = 1
    samples: int = 3
    reads: int = 200_000
    read_length: int = 75
    paired: bool = False
    mito_length: int = 16_000
    nuclear_length: int = 2_000_000
    bacterium_length: int = 300_000
    srnas: int = 30
    degradation: float = 0.08
    error_rate: float = 0.002
    max_hits: int = 50
    ambiguous: str = "report"


@dataclass
class Source:
    """Something that produces reads: one smallRNA, one degraded transcript, a pool."""
    id: str
    category: str
    group: str
    genome: str = ""
    contig: str = ""
    start: int = 0  # 0-based; end may exceed the contig length on circular contigs
    end: int = 0
    strand: str = "+"
    expected: str = ""
    means: list[float] = field(default_factory=list)
    template: str = ""  # oriented sequence with FLANK bases on each side (smallRNAs)
    alt: dict[int, str] = field(default_factory=dict)  # template position -> base
    alt_fraction: list[float] = field(default_factory=list)
    counts: list[int] = field(default_factory=list)
    sequence: str = ""


FLANK = 5


class Genome:
    def __init__(self, name: str, contig: str, seq: str, circular: bool):
        self.name, self.contig, self.seq, self.circular = name, contig, seq, circular

    def region(self, start: int, end: int, strand: str = "+") -> str:
        n = len(self.seq)
        if self.circular:
            s = "".join(self.seq[i % n] for i in range(start, end))
        else:
            s = self.seq[max(0, start):min(n, end)]
        return s if strand == "+" else revcomp(s)


class Simulator:
    def __init__(self, p: Params):
        self.p = p
        self.rng = random.Random(p.seed)
        self.samples = [f"sim{i + 1}" for i in range(p.samples)]
        self.sources: list[Source] = []
        

    # sequence helpers -----------------------------------------------------------------

    def random_seq(self, n: int, at: float) -> str:
        weights = [at / 2, (1 - at) / 2, (1 - at) / 2, at / 2]
        return "".join(self.rng.choices("ACGT", weights, k=n))

    def weighted(self, table: dict):
        return self.rng.choices(list(table), list(table.values()))[0]

    def poisson(self, lam: float) -> int:
        if lam <= 0:
            return 0
        if lam > 30:
            return max(0, round(self.rng.gauss(lam, math.sqrt(lam))))
        limit, k, prod = math.exp(-lam), 0, self.rng.random()
        while prod > limit:
            k += 1
            prod *= self.rng.random()
        return k

    def neg_binomial(self, mean: float) -> int:
        if mean <= 0:
            return 0
        return self.poisson(self.rng.gammavariate(1 / DISPERSION, mean * DISPERSION))

    def per_sample(self, rpm: float) -> list[float]:
        return [rpm * self.p.reads / 1e6] * self.p.samples

    # genomes ----------------------------------------------------------------------------

    def build_mito(self) -> Genome:
        p = self.p
        seq = self.random_seq(p.mito_length, at=0.77)
        layout = [("tRNA", f"trn{i + 1}", self.rng.randint(66, 72)) for i in range(22)]
        layout += [("rRNA", "rrnL", 1300), ("rRNA", "rrnS", 800)]
        cds_names = [f"cds{i + 1}" for i in range(13)]
        raw = [self.rng.uniform(0.3, 1.6) for _ in cds_names]
        spacer = 8
        free = p.mito_length - 1060 - sum(l for _, _, l in layout) - spacer * (len(layout) + 13)
        if free < 13 * 150:
            raise ValueError("mito_length is too small for the simulated gene layout")
        layout += [("CDS", n, int(free * r / sum(raw))) for n, r in zip(cds_names, raw)]
        self.rng.shuffle(layout)
        pos = 60  # keep the origin free of genes, for the smallRNA across it
        self.mito_features = []
        for ftype, name, length in layout:
            strand = "+" if self.rng.random() < 0.7 else "-"
            self.mito_features.append((ftype, name, pos, pos + length, strand))
            pos += length + spacer
        self.mito_features.append(("D_loop", "control_region", pos, p.mito_length, "+"))
        return Genome("mito", "chrM", seq, circular=True)

    def free_position(self, mito: Genome, length: int, taken: list[tuple], strand_of=None,
                      inside=None) -> int:
        """Pick a start away from tRNAs, the origin and other smallRNAs."""
        n = len(mito.seq)
        trnas = [(s, e) for t, _, s, e, _ in self.mito_features if t == "tRNA"]
        lo, hi = inside if inside else (40, n - length - 40)
        for _ in range(10_000):
            s = self.rng.randint(lo, hi - length)
            e = s + length
            if any(s < te + 20 and e > ts - 20 for ts, te in trnas):
                continue
            if any(s < te + 40 and e > ts - 40 for ts, te in taken):
                continue
            taken.append((s, e))
            return s
        raise ValueError("could not place all smallRNAs; lower --srnas or raise --mito-length")

    def srna(self, genome: Genome, sid: str, category: str, start: int, length: int,
             strand: str, means: list[float], expected: str, group: str = "mito_srna") -> Source:
        src = Source(sid, category, group, genome.name, genome.contig, start, start + length,
                     strand, expected, means)
        src.template = genome.region(start - FLANK, start + length + FLANK, strand)
        src.sequence = src.template[FLANK:FLANK + length]
        self.sources.append(src)
        return src

    def add_mito_sources(self, mito: Genome, nuclear: list[str]):
        p, rng = self.p, self.rng
        taken: list[tuple] = []
        n = len(mito.seq)

        def rpm():
            return 10 ** rng.uniform(math.log10(300), math.log10(20_000))

        for i in range(p.srnas):
            length = rng.randint(19, 26)
            s = self.free_position(mito, length, taken)
            self.srna(mito, f"srna_{i + 1:03d}", "srna", s, length, rng.choice("+-"),
                      self.per_sample(rpm()), "pass")

        rrnl = next(f for f in self.mito_features if f[1] == "rrnL")
        s = self.free_position(mito, 22, taken, inside=(rrnl[2] + 50, rrnl[3] - 50))
        self.srna(mito, "srna_in_rrna", "srna_in_rrna", s, 22, rrnl[4],
                  self.per_sample(rpm()), "pass")

        self.srna(mito, "srna_origin", "srna_origin", n - 9, 22, "+", self.per_sample(3000), "pass")

        s = self.free_position(mito, 22, taken)
        self.srna(mito, "srna_antisense_a", "srna_antisense_pair", s, 22, "+",
                  self.per_sample(2000), "pass")
        self.srna(mito, "srna_antisense_b", "srna_antisense_pair", s, 22, "-",
                  self.per_sample(1500), "pass")

        s = self.free_position(mito, 54, taken)
        self.srna(mito, "srna_close_a", "srna_close_pair", s, 22, "+", self.per_sample(2000), "pass")
        self.srna(mito, "srna_close_b", "srna_close_pair", s + 32, 22, "+",
                  self.per_sample(1000), "pass")

        s = self.free_position(mito, 38, taken)
        self.srna(mito, "srna_overlap_a", "srna_overlap_pair", s, 22, "+",
                  self.per_sample(2000), "pass")
        self.srna(mito, "srna_overlap_b", "srna_overlap_pair", s + 16, 22, "+",
                  self.per_sample(1000), "pass")

        # NUMTs: the nuclear copies (smallRNA +- 150 nt) are written by build_nuclear; the
        # whole copied segment is kept free of other smallRNAs.
        s = self.free_position(mito, 322, taken) + 150
        expected = "pass" if p.ambiguous == "split" else "missed"
        self.srna(mito, "srna_numt_identical", "srna_numt_identical", s, 22, "+",
                  self.per_sample(3000), expected)
        nuclear.append(("identical", mito.region(s - 150, s + 172), None))
        s = self.free_position(mito, 322, taken) + 150
        self.srna(mito, "srna_numt_1mm", "srna_numt_1mm", s, 22, "+", self.per_sample(3000), "pass")
        nuclear.append(("1mm", mito.region(s - 150, s + 172), 150 + 11))

        # Fixed SNP: the reads carry a base that differs from the reference.
        s = self.free_position(mito, 22, taken)
        src = self.srna(mito, "srna_snp", "srna_snp", s, 22, "+", self.per_sample(3000), "pass")
        pos = FLANK + 8
        src.alt = {pos: rng.choice([b for b in "ACGT" if b != src.template[pos]])}
        src.alt_fraction = [1.0] * p.samples
        s = self.free_position(mito, 22, taken)
        src = self.srna(mito, "srna_heteroplasmy", "srna_heteroplasmy", s, 22, "+",
                        self.per_sample(3000), "pass")
        pos = FLANK + 12
        src.alt = {pos: rng.choice([b for b in "ACGT" if b != src.template[pos]])}
        src.alt_fraction = [0.5] + [0.0] * (p.samples - 1)

        s = self.free_position(mito, 22, taken)
        low = 1.5 / (p.reads / 1e6)
        self.srna(mito, "srna_low", "srna_low", s, 22, "+", self.per_sample(low), "fail")
        if p.samples >= 2:
            s = self.free_position(mito, 22, taken)
            means = [3000 * p.reads / 1e6] + [0.0] * (p.samples - 1)
            self.srna(mito, "srna_one_sample", "srna_one_sample", s, 22, "-", means, "fail")

        # tRNA fragments: 5' fragments map; 3' fragments with CCA mostly do not (up to 3
        # mismatches, fewer when the bases after the tRNA happen to match).
        trnas = [f for f in self.mito_features if f[0] == "tRNA"]
        weights = [rng.lognormvariate(0, 1.2) for _ in trnas]
        total = COMPOSITION["mito_trf"] * p.reads
        for (ftype, name, s, e, strand), w in zip(trnas, weights):
            mean = total * w / sum(weights)
            for kind, share in (("trf5", 0.6), ("trf3_cca", 0.4)):
                # weakly expressed fragments may or may not pass the filters
                expected = "pass" if mean * share >= 20 else "any"
                if kind == "trf3_cca":
                    after = mito.region(e, e + 3) if strand == "+" else mito.region(s - 3, s, "-")
                    # reads still align when the genome after the tRNA almost spells CCA
                    expected = "fail" if sum(a != b for a, b in zip(after, "CCA")) > 1 else "any"
                src = Source(f"{kind}_{name}", kind, "mito_trf", mito.name, mito.contig, s, e,
                             strand, expected, [mean * share] * p.samples)
                src.template = mito.region(s, e, strand)
                src.sequence = src.template
                # truth interval: the genomic part of the fragments
                if kind == "trf5":
                    src.start, src.end = (s, s + 35) if strand == "+" else (e - 35, e)
                else:
                    src.start, src.end = (e - 19, e) if strand == "+" else (s, s + 19)
                self.sources.append(src)

        # Degradation of rRNAs (half of the reads) and mRNAs, sense strand.
        rrnas = [f for f in self.mito_features if f[0] == "rRNA"]
        cds = [f for f in self.mito_features if f[0] == "CDS"]
        total = p.degradation * p.reads
        for group, feats in (("rRNA", rrnas), ("CDS", cds)):
            ws = [rng.lognormvariate(0, 0.8) for _ in feats]
            for (ftype, name, s, e, strand), w in zip(feats, ws):
                mean = total / 2 * w / sum(ws)
                self.sources.append(Source(f"degradation_{name}", "degradation", "mito_degradation",
                                           mito.name, mito.contig, s, e, strand, "background",
                                           [mean] * p.samples, mito.region(s, e, strand)))

    def build_nuclear(self, numts: list) -> Genome:
        p, rng = self.p, self.rng
        seq = list(self.random_seq(p.nuclear_length, at=0.6))
        used: list[tuple] = []

        def place(length: int) -> int:
            for _ in range(10_000):
                s = rng.randint(1000, p.nuclear_length - length - 1000)
                if all(s + length < us or s > ue for us, ue in used):
                    used.append((s - 100, s + length + 100))
                    return s
            raise ValueError("nuclear genome too small for the simulated features")

        for kind, segment, mm_pos in numts:
            segment = list(segment)
            if mm_pos is not None:
                segment[mm_pos] = rng.choice([b for b in "ACGT" if b != segment[mm_pos]])
            s = place(len(segment))
            seq[s:s + len(segment)] = segment
        repeat = self.random_seq(300, at=0.6)
        repeat_pos = []
        for _ in range(p.max_hits + 30):
            s = place(len(repeat))
            seq[s:s + len(repeat)] = repeat
            repeat_pos.append(s)
        nuc = Genome("nuclear", "chrN", "".join(seq), circular=False)

        weights = [rng.lognormvariate(0, 1.5) for _ in range(40)]
        for i, w in enumerate(weights):
            length = rng.randint(20, 24)
            s = place(length)
            mean = COMPOSITION["nuclear_srna"] * p.reads * w / sum(weights)
            self.srna(nuc, f"nuclear_srna_{i + 1:03d}", "nuclear_srna", s, length,
                      rng.choice("+-"), [mean] * p.samples, "excluded", group="nuclear_srna")
        for i in range(2):
            s = place(10_000)
            mean = COMPOSITION["nuclear_pirna"] * p.reads / 2
            self.sources.append(Source(f"pirna_cluster_{i + 1}", "pirna", "nuclear_pirna",
                                       nuc.name, nuc.contig, s, s + 10_000, "+", "excluded",
                                       [mean] * p.samples, nuc.region(s, s + 10_000)))
        ws = [rng.lognormvariate(0, 1) for _ in range(50)]
        for i, w in enumerate(ws):
            length = rng.randint(1000, 3000)
            s = place(length)
            strand = rng.choice("+-")
            mean = COMPOSITION["nuclear_degradation"] * p.reads * w / sum(ws)
            self.sources.append(Source(f"nuclear_gene_{i + 1:03d}", "degradation",
                                       "nuclear_degradation", nuc.name, nuc.contig, s, s + length,
                                       strand, "excluded", [mean] * p.samples,
                                       nuc.region(s, s + length, strand)))
        self.sources.append(Source("nuclear_repeat", "repeat", "nuclear_repeat", nuc.name,
                                   nuc.contig, repeat_pos[0], repeat_pos[0] + 300, "+",
                                   "too_many_hits",
                                   [COMPOSITION["nuclear_repeat"] * p.reads] * p.samples, repeat))
        return nuc

    def build_bacterium(self) -> Genome:
        p, rng = self.p, self.rng
        seq = list(self.random_seq(p.bacterium_length, at=0.5))
        operon = self.random_seq(1500, at=0.45)
        starts = [int(p.bacterium_length * f) for f in (0.1, 0.4, 0.7)]
        for s in starts:
            seq[s:s + 1500] = operon
        mean = COMPOSITION["bacterium"] * p.reads
        self.sources.append(Source("bacterium_rrna", "degradation", "bacterium", "bacterium",
                                   "chrB", starts[0], starts[0] + 1500, "+", "excluded",
                                   [mean] * p.samples, operon))
        return Genome("bacterium", "chrB", "".join(seq), circular=True)

    def add_pools(self):
        p = self.p
        self.sources.append(Source("foreign", "foreign", "foreign", expected="unmapped",
                                   means=[COMPOSITION["foreign"] * p.reads] * p.samples))
        self.foreign = [self.random_seq(self.rng.randint(18, 32), at=0.5) for _ in range(2000)]
        self.foreign_w = [self.rng.lognormvariate(0, 1.5) for _ in self.foreign]
        self.sources.append(Source("short", "short", "short", expected="filtered",
                                   means=[COMPOSITION["short"] * p.reads] * p.samples))

    # reads ------------------------------------------------------------------------------

    def insert(self, src: Source, sample: int) -> str:
        rng = self.rng
        cat = src.category
        if src.template and (cat.startswith("srna") or cat == "nuclear_srna"):
            tmpl = src.template
            if src.alt and rng.random() < src.alt_fraction[sample]:
                tmpl = list(tmpl)
                for pos, base in src.alt.items():
                    tmpl[pos] = base
                tmpl = "".join(tmpl)
            length = src.end - src.start
            a = FLANK + self.weighted(FIVE_SHIFT)
            b = FLANK + length + self.weighted(THREE_SHIFT)
            seq = tmpl[a:b]
            if rng.random() < TAIL_RATE:
                seq += "".join(self.weighted(TAIL_BASES) for _ in range(rng.choice((1, 1, 1, 2))))
            return seq
        if cat == "trf5":
            return src.template[:rng.choice((rng.randint(18, 22), rng.randint(30, 35)))]
        if cat == "trf3_cca":
            return src.template[-rng.randint(15, 19):] + "CCA"
        if cat == "pirna":
            length = rng.randint(24, 30)
            s = rng.randint(0, len(src.template) - length)
            seq = src.template[s:s + length]
            seq = "T" + seq[1:] if rng.random() < 0.8 else seq
            return seq if rng.random() < 0.5 else revcomp(seq)
        if cat in ("degradation", "repeat"):
            length = rng.randint(16, 45)
            s = rng.randint(0, len(src.template) - length)
            return src.template[s:s + length]
        if cat == "foreign":
            return rng.choices(self.foreign, self.foreign_w)[0]
        if cat == "short":
            return self.random_seq(rng.randint(0, 15), at=0.5)
        raise ValueError(f"no read model for {cat}")

    def errors(self, read: str) -> str:
        rate, n = self.p.error_rate, len(read)
        out = list(read)
        for i in range(n):
            if self.rng.random() < rate * (1 + 2 * i / n):
                out[i] = self.rng.choice([b for b in "ACGT" if b != out[i]])
        return "".join(out)

    def build_read(self, insert: str, adapter: str) -> str:
        tail = self.random_seq(self.p.read_length, at=0.5)
        return self.errors((insert + adapter + tail)[:self.p.read_length])

    # main -------------------------------------------------------------------------------

    def run(self):
        p = self.p
        if p.ambiguous not in ("report", "discard", "split"):
            raise ValueError("ambiguous must be report, discard or split")
        out = Path(p.outdir)
        for sub in ("data", "config", "truth"):
            (out / sub).mkdir(parents=True, exist_ok=True)

        mito = self.build_mito()
        numts: list = []
        self.add_mito_sources(mito, numts)
        nuc = self.build_nuclear(numts)
        bact = self.build_bacterium()
        self.add_pools()
        for g in (mito, nuc, bact):
            with open(out / "data" / f"{g.name}.fasta", "w") as fh:
                fh.write(f">{g.contig}\n")
                for i in range(0, len(g.seq), 80):
                    fh.write(g.seq[i:i + 80] + "\n")
        with open(out / "data" / "mito.gff3", "w") as fh:
            fh.write("##gff-version 3\n")
            for ftype, name, s, e, strand in self.mito_features:
                fh.write(f"{mito.contig}\tsimulate\t{ftype}\t{s + 1}\t{e}\t.\t{strand}\t.\t"
                         f"ID={name};Name={name}\n")

        for src in self.sources:
            src.counts = [self.neg_binomial(m) for m in src.means]
        quality = "F" * p.read_length
        for i, sample in enumerate(self.samples):
            reads = []
            for src in self.sources:
                reads.extend(self.insert(src, i) for _ in range(src.counts[i]))
            self.rng.shuffle(reads)
            r1 = gzip.open(out / "data" / f"{sample}_1.fastq.gz", "wt", compresslevel=1)
            r2 = (gzip.open(out / "data" / f"{sample}_2.fastq.gz", "wt", compresslevel=1)
                  if p.paired else None)
            for j, ins in enumerate(reads):
                r1.write(f"@{sample}_{j}\n{self.build_read(ins, ADAPTER_R1)}\n+\n{quality}\n")
                if r2:
                    r2.write(f"@{sample}_{j}\n{self.build_read(revcomp(ins), ADAPTER_R2)}\n"
                             f"+\n{quality}\n")
            r1.close()
            if r2:
                r2.close()

        self.write_truth(out, {g.name: g for g in (mito, nuc, bact)})
        self.write_config(out)

    def write_truth(self, out: Path, genomes: dict[str, Genome]):
        rows = []
        for src in self.sources:
            g = genomes.get(src.genome)
            n = len(g.seq) if g else 0
            wraps = bool(g and g.circular and src.end > n)
            rows.append({
                "feature_id": src.id, "category": src.category, "group": src.group,
                "genome": src.genome, "contig": src.contig,
                "start": src.start + 1 if g else "", "end": (src.end - n if wraps else src.end) if g else "",
                "strand": src.strand if g else "", "wraps_origin": int(wraps),
                "length": src.end - src.start if g else "",
                "sequence": src.sequence if src.category.startswith("srna") else "",
                "expected": src.expected, "total_reads": sum(src.counts),
                **{f"reads_{s}": c for s, c in zip(self.samples, src.counts)},
            })
        columns = ["feature_id", "category", "group", "genome", "contig", "start", "end", "strand",
                   "wraps_origin", "length", "sequence", "expected", "total_reads",
                   *[f"reads_{s}" for s in self.samples]]
        write_tsv(str(out / "truth" / "features.tsv"), rows, columns)

        groups: dict[str, list[int]] = {}
        for src in self.sources:
            acc = groups.setdefault(src.group, [0] * self.p.samples)
            for i, c in enumerate(src.counts):
                acc[i] += c
        write_tsv(str(out / "truth" / "composition.tsv"),
                  [{"group": g, **dict(zip(self.samples, v))} for g, v in sorted(groups.items())],
                  ["group", *self.samples])

    def write_config(self, out: Path):
        p = self.p
        with open(out / "config" / "samples.tsv", "w") as fh:
            fh.write("sample\tfq1\tfq2\n")
            for s in self.samples:
                fq2 = f"data/{s}_2.fastq.gz" if p.paired else ""
                fh.write(f"{s}\tdata/{s}_1.fastq.gz\t{fq2}\n")
        (out / "config" / "config.yaml").write_text(f"""\
# Written by `smithhunter simulate` (seed {p.seed}). Run from this directory:
#   snakemake -s <SmithHunter2>/workflow/Snakefile --cores 2
#   smithhunter evaluate --truth truth --results results
samples: config/samples.tsv
outdir: results
threads: 2

genomes:
  mito:
    fasta: data/mito.fasta
    role: focus
    topology: circular
    annotation: data/mito.gff3
    targets: ago
  nuclear:
    fasta: data/nuclear.fasta
    role: exclude
  bacterium:
    fasta: data/bacterium.fasta
    role: exclude
    topology: circular

ambiguous: {p.ambiguous}

trimming:
  adapter_r1: {ADAPTER_R1}
  adapter_r2: {ADAPTER_R2}
  min_length: 18
  max_length: 35

mapping:
  mismatches: 1
  max_hits: {p.max_hits}

loci:
  merge_gap: 0
  peaks: true
  peak_window: 2
  peak_pvalue: 0.001
  background_flank: 50
  multimap: fractional
  min_rpm: 5
  min_count: 5
  min_samples: null

ends:
  n_thre: 0.5
  penalty: 0.1
  min_five_score: 0.5
  min_three_score: 0.0
""")


def simulate(p: Params) -> Simulator:
    sim = Simulator(p)
    sim.run()
    return sim
