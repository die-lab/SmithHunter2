"""Overlap loci with GFF3 features of their own genome."""

from __future__ import annotations

import re
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from urllib.parse import unquote

from .seqio import open_text

SKIP_TYPES = {"region", "chromosome", "source", "databank_entry", "contig", "scaffold"}
# Most specific first: the first type found among the overlaps becomes the class.
# tmRNA, RNase_P_RNA, SRP_RNA, antisense_RNA are the NCBI types of prokaryotic housekeeping /
# regulatory RNAs: without them such loci fell through to their child "exon".
CLASS_PRIORITY = ("tRNA", "rRNA", "tmRNA", "RNase_P_RNA", "RNase_MRP_RNA", "SRP_RNA",
                  "snoRNA", "snRNA", "miRNA", "pre_miRNA", "piRNA", "antisense_RNA", "ncRNA",
                  "lnc_RNA", "rep_origin", "D_loop", "repeat_region", "five_prime_UTR",
                  "three_prime_UTR", "CDS", "exon", "mRNA", "transcript", "gene")

# Abundant structural (housekeeping) RNAs: their degradation fragments are frequent and pass
# expression and end filters (wAlbB benchmark: 25% of the candidates, 58% of their reads).
STRUCTURAL_TYPES = {"tRNA", "rRNA", "tmRNA", "RNase_P_RNA", "RNase_MRP_RNA", "SRP_RNA",
                    "snoRNA", "snRNA", "scRNA", "telomerase_RNA", "Y_RNA", "vault_RNA"}
# ncRNA features recognised by gene name: tmRNA (ssrA), 6S RNA (ssrS), RNase P (rnpB),
# SRP / 4.5S (ffs), 7SK / 7SL; rRNA names when the type is generic.
STRUCTURAL_NAMES = re.compile(
    r"^(ssrA|ssrS|rnpB|ffs|rrs|rrl|rrf|RN7SK|RN7SL\d*|RNY\d*|5S|5\.8S|16S|18S|23S|28S|6S|4\.5S)$"
    r"|tmRNA|transfer-messenger|RNase P|signal recognition particle|\b6S RNA|ribosomal RNA",
    re.IGNORECASE)


@dataclass
class Feature:
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive
    strand: str
    type: str
    name: str
    product: str = ""


# Types whose name / product may identify a structural RNA (never protein-coding features:
# "ribosomal RNA methyltransferase" is a CDS product).
_NAMED_RNA_TYPES = {"ncRNA", "misc_RNA", "exon", "gene", "transcript", "primary_transcript"}


def is_structural(f: Feature) -> bool:
    if f.type in STRUCTURAL_TYPES:
        return True
    return f.type in _NAMED_RNA_TYPES and bool(
        STRUCTURAL_NAMES.search(f.name) or (f.product and STRUCTURAL_NAMES.search(f.product)))


def structural(features: list[Feature], strand: str) -> str:
    """'sense' / 'antisense' when the locus overlaps a structural RNA, '' otherwise.

    Sense wins when structural RNAs lie on both strands; unstranded features count as sense.
    """
    hits = [f for f in features if is_structural(f)]
    if not hits:
        return ""
    if any(f.strand == strand or f.strand not in "+-" for f in hits):
        return "sense"
    return "antisense"


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
                features[f[0]].append(Feature(int(f[3]) - 1, int(f[4]), f[6], f[2], unquote(name),
                                              unquote(attrs.get("product", ""))))
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
