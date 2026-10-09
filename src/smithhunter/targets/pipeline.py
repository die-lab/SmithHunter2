"""Seed site scan of candidates and decoys over one target set, with site-count FDR.

Identical region sequences (isoforms sharing a UTR) are scanned once. Counts for the
FDR use unique sequences, so isoforms do not inflate them; the site table lists every
region that carries the site.
"""

from __future__ import annotations

from collections import Counter, defaultdict

from ..seqio import read_fasta, write_tsv
from .fdr import site_count_fdr
from .regions import Region
from .sites import context, scan, site_classes


def read_candidates(path: str) -> list[tuple[str, str, str]]:
    """(id, genome, sequence) from module A's candidates.fasta or any FASTA.

    The genome is taken from a ``genome|contig:...`` second header field; it is empty
    when the header has none.
    """
    out = []
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:].split()
            elif header and line:
                genome = header[1].split("|")[0] if len(header) > 1 and "|" in header[1] else ""
                out.append((header[0], genome, line.upper().replace("U", "T")))
                header = None
    return out


SITE_COLUMNS = ["candidate", "region_id", "transcript", "gene", "region", "site_class",
                "site_start", "site_end", "genomic", "dist_5p", "dist_3p", "au_flanks",
                "supp_pairing", "near_stop", "class_fdr"]
TARGET_COLUMNS_HEAD = ["candidate", "transcript", "gene", "sites", "best_class", "best_class_fdr"]
FDR_COLUMNS = ["candidate", "site_class", "real_sites", "decoy_mean", "decoy_sd",
               "enrichment", "fdr"]


def run_sites(candidates: list[tuple[str, str]], decoys: list[tuple[str, str, str]],
              regions: list[Region], seed=(2, 8), classes=None):
    classes = tuple(classes or site_classes(seed))
    by_seq: dict[str, list[Region]] = defaultdict(list)
    for r in regions:
        by_seq[r.sequence].append(r)
    unique = list(by_seq)

    decoys_of: dict[str, list[str]] = defaultdict(list)
    for _, cid, seq in decoys:
        decoys_of[cid].append(seq)

    site_rows, fdr_rows = [], []
    per_target: dict[tuple, Counter] = defaultdict(Counter)
    pooled_real: Counter = Counter()
    pooled_decoy: list[Counter] = []
    for cid, guide in candidates:
        real: Counter = Counter()
        found = []
        for seq in unique:
            for site in scan(guide, seq, seed, classes):
                real[site.site_class] += 1
                found.append((seq, site))
        decoy_counts = []
        for dseq in decoys_of.get(cid, []):
            c: Counter = Counter()
            for seq in unique:
                for site in scan(dseq, seq, seed, classes):
                    c[site.site_class] += 1
            decoy_counts.append(c)
        rows = site_count_fdr(real, decoy_counts, classes)
        class_fdr = {r["site_class"]: r["fdr"] for r in rows}
        fdr_rows += [{"candidate": cid, **r} for r in rows]
        pooled_real.update(real)
        for i, c in enumerate(decoy_counts):
            if i < len(pooled_decoy):
                pooled_decoy[i].update(c)
            else:
                pooled_decoy.append(Counter(c))

        for seq, site in found:
            for r in by_seq[seq]:
                g = r.genomic(site.start)
                g_end = r.genomic(site.end - 1)
                genomic = (f"{r.contig}:{min(g, g_end) + 1}-{max(g, g_end) + 1}({r.strand})"
                           if g is not None and g_end is not None else "")
                site_rows.append({
                    "candidate": cid, "region_id": r.id, "transcript": r.transcript,
                    "gene": r.gene, "region": r.region, "site_class": site.site_class,
                    "site_start": site.start + 1, "site_end": site.end, "genomic": genomic,
                    **context(guide, seq, site, r.region),
                    "class_fdr": class_fdr.get(site.site_class, ""),
                })
                per_target[(cid, r.transcript, r.gene)][site.site_class] += 1

    fdr_rows += [{"candidate": "all", **r}
                 for r in site_count_fdr(pooled_real, pooled_decoy, classes)]

    target_rows = []
    fdr_of = {(r["candidate"], r["site_class"]): r["fdr"] for r in fdr_rows}
    for (cid, transcript, gene), counts in per_target.items():
        best = next(c for c in classes if counts.get(c))
        target_rows.append({
            "candidate": cid, "transcript": transcript, "gene": gene,
            "sites": sum(counts.values()), "best_class": best,
            "best_class_fdr": fdr_of.get((cid, best), ""),
            **{f"n_{c}": counts.get(c, 0) for c in classes},
        })
    order = {c: i for i, c in enumerate(classes)}
    target_rows.sort(key=lambda r: (r["candidate"], order[r["best_class"]], -r["sites"]))
    return classes, site_rows, target_rows, fdr_rows


def write_outputs(classes, site_rows, target_rows, fdr_rows, sites_path, targets_path,
                  fdr_path) -> None:
    write_tsv(sites_path, site_rows, SITE_COLUMNS)
    write_tsv(targets_path, target_rows, TARGET_COLUMNS_HEAD + [f"n_{c}" for c in classes])
    write_tsv(fdr_path, fdr_rows, FDR_COLUMNS)


def read_decoys(path: str) -> list[tuple[str, str, str]]:
    out = []
    for name, seq in read_fasta(path):
        out.append((name, name.split("__decoy")[0], seq))
    return out
