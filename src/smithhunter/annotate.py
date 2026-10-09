"""Overlap loci with GFF3 features of their own genome."""

from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from urllib.parse import unquote

from .seqio import open_text

SKIP_TYPES = {"region", "chromosome", "source", "databank_entry", "contig", "scaffold"}
# Most specific first: the first type found among the overlaps becomes the class.
CLASS_PRIORITY = ("tRNA", "rRNA", "snoRNA", "snRNA", "miRNA", "pre_miRNA", "piRNA", "ncRNA",
                  "lnc_RNA", "rep_origin", "D_loop", "repeat_region", "five_prime_UTR",
                  "three_prime_UTR", "CDS", "exon", "mRNA", "transcript", "gene")


@dataclass
class Feature:
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive
    strand: str
    type: str
    name: str


class Annotation:
    def __init__(self, features: dict[str, list[Feature]]):
        self.features = {c: sorted(fs, key=lambda f: f.start) for c, fs in features.items()}
        self.starts = {c: [f.start for f in fs] for c, fs in self.features.items()}
        self.max_len = {c: max((f.end - f.start for f in fs), default=0)
                        for c, fs in self.features.items()}

    @classmethod
    def from_gff(cls, path: str) -> "Annotation":
        features: dict[str, list[Feature]] = defaultdict(list)
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] in SKIP_TYPES:
                    continue
                attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                name = attrs.get("Name") or attrs.get("gene") or attrs.get("ID") or "."
                features[f[0]].append(Feature(int(f[3]) - 1, int(f[4]), f[6], f[2], unquote(name)))
        return cls(features)

    def overlaps(self, contig: str, start: int, end: int) -> list[Feature]:
        fs = self.features.get(contig)
        if not fs:
            return []
        lo = bisect_left(self.starts[contig], start - self.max_len[contig])
        hi = bisect_left(self.starts[contig], end)
        return [f for f in fs[lo:hi] if f.start < end and f.end > start]


def describe(features: list[Feature], strand: str) -> dict[str, str]:
    if not features:
        return {"annotation_class": "intergenic", "annotation_types": "", "annotation_names": "",
                "annotation_orientation": ""}
    types = sorted({f.type for f in features})
    cls = next((t for t in CLASS_PRIORITY if t in types), types[0])
    best = [f for f in features if f.type == cls]
    if any(f.strand == strand for f in best):
        orientation = "sense"
    elif all(f.strand not in "+-" for f in best):
        orientation = "unstranded"
    else:
        orientation = "antisense"
    return {
        "annotation_class": cls,
        "annotation_types": ",".join(types),
        "annotation_names": ",".join(sorted({f.name for f in features})),
        "annotation_orientation": orientation,
    }
