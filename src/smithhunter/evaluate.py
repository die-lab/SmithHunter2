"""Compare module A results with the truth written by ``smithhunter simulate``.

Every simulated smallRNA (and tRNA fragment) of a focus genome is matched to the
loci on the same contig and strand. A locus belongs to the smallRNA when its
representative (most abundant) read starts within ``OWN_DISTANCE`` bases of the
smallRNA's 5' end. Observed outcomes:

- ``pass`` / ``fail``  its own locus, holding no other simulated smallRNA, passes / fails
- ``merged``      its locus also contains another simulated smallRNA, or it has no own
                  locus but lies inside a locus that belongs to something else
- ``background``  its reads are only in a background locus (no 5' peak called)
- ``split``       more than one locus of its own, the second with at least
                  ``SPLIT_FRACTION`` of the reads of the first
- ``partial``     overlapping loci exist but none is its own or contains it
- ``missed``      no overlapping locus

``as_expected`` compares it with the outcome the simulator planned. For loci found
the table reports how far they extend beyond the smallRNA (``extra_span``, large when
degradation fragments chain into the locus), whether the representative sequence is
the canonical one, and the 5' offset of the representative read. Passing loci that
contain no simulated smallRNA are listed as false positives, with the source that
explains them (degradation, tRNA fragments).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from .seqio import read_tsv, write_tsv

TARGET_CATEGORIES = ("srna", "trf5", "trf3_cca")
OWN_DISTANCE = 3  # a locus belongs to a smallRNA if its representative starts this close
SPLIT_FRACTION = 0.1  # a second own locus splits the smallRNA if it has this share of reads


def interval(row: dict, length: int, circular: bool) -> tuple[int, int]:
    """0-based half-open interval; on circular contigs the end may exceed the length."""
    start, end = int(row["start"]) - 1, int(row["end"])
    if circular and int(row.get("wraps_origin") or 0):
        end += length
    return start, end


def overlaps(a: tuple[int, int], b: tuple[int, int], length: int, circular: bool) -> bool:
    shifts = (0, length, -length) if circular else (0,)
    return any(a[0] < b[1] + s and a[1] > b[0] + s for s in shifts)


def overlap_length(a, b, length, circular) -> int:
    shifts = (0, length, -length) if circular else (0,)
    return max(max(0, min(a[1], b[1] + s) - max(a[0], b[0] + s)) for s in shifts)


def covers(outer, inner, length, circular) -> bool:
    shifts = (0, length, -length) if circular else (0,)
    return any(outer[0] <= inner[0] + s and outer[1] >= inner[1] + s for s in shifts)


def five_offset(feature: dict, locus: dict, length: int) -> int:
    """Distance of the representative read's 5' end from the smallRNA's 5' end."""
    if feature["strand"] == "+":
        d = int(locus["rep_start"]) - int(feature["start"])
    else:
        d = int(feature["end"]) - int(locus["rep_end"])
    return (d + length // 2) % length - length // 2


def is_target(row: dict) -> bool:
    return row["category"].startswith(TARGET_CATEGORIES)


def evaluate(truth_dir: str, results_dir: str):
    features = read_tsv(f"{truth_dir}/features.tsv")
    loci = read_tsv(f"{results_dir}/discovery/loci.tsv")
    contigs = {(c["genome"], c["contig"]): (int(c["length"]), c["circular"] == "1")
               for c in read_tsv(f"{results_dir}/reference/contigs.tsv")}
    focus = {g["name"] for g in read_tsv(f"{results_dir}/reference/genomes.tsv")
             if g["role"] == "focus"}

    by_key: dict[tuple, list[dict]] = defaultdict(list)
    for locus in loci:
        length, circular = contigs[(locus["genome"], locus["contig"])]
        locus["_iv"] = interval(locus, length, circular)
        by_key[(locus["genome"], locus["contig"], locus["strand"])].append(locus)

    targets = [f for f in features if f["genome"] in focus and is_target(f)]
    srna_like = [f for f in targets if f["category"].startswith("srna")]
    for f in features:
        if (f["genome"], f["contig"]) in contigs:
            length, circular = contigs[(f["genome"], f["contig"])]
            f["_iv"] = interval(f, length, circular)

    rows = []
    for f in targets:
        length, circular = contigs[(f["genome"], f["contig"])]
        hits = [l for l in by_key[(f["genome"], f["contig"], f["strand"])]
                if overlaps(f["_iv"], l["_iv"], length, circular)]
        row = {k: f[k] for k in ("feature_id", "category", "start", "end", "strand", "length",
                                 "expected", "total_reads")}
        own = sorted((l for l in hits if abs(five_offset(f, l, length)) <= OWN_DISTANCE),
                     key=lambda l: -float(l["total_count"]))
        # weak loci next to the main one (e.g. background starting nearby) do not split it
        own = [l for l in own if float(l["total_count"]) >= SPLIT_FRACTION * float(own[0]["total_count"])]
        others = [o for o in srna_like if o is not f and o["strand"] == f["strand"]
                  and o["contig"] == f["contig"]]
        if len(own) > 1:
            observed = "split"
        elif own:
            locus = own[0]
            if locus.get("locus_type") == "background":
                observed = "background"
            elif any(covers(locus["_iv"], o["_iv"], length, circular) for o in others):
                observed = "merged"
            else:
                observed = "pass" if locus["pass"] == "1" else "fail"
        elif any(covers(l["_iv"], f["_iv"], length, circular) for l in hits):
            covering = [l for l in hits if covers(l["_iv"], f["_iv"], length, circular)]
            observed = "background" if all(l.get("locus_type") == "background"
                                           for l in covering) else "merged"
        else:
            observed = "partial" if hits else "missed"
        best = own[0] if own else max(hits, key=lambda l: float(l["total_count"]), default=None)
        if best:
            row.update({
                "locus_id": ",".join(l["locus_id"] for l in hits),
                "locus_start": best["start"], "locus_end": best["end"],
                "extra_span": int(best["span"]) - int(f["length"]),
                "rep_sequence": best["rep_sequence"],
                "rep_is_canonical": int(best["rep_sequence"] == f["sequence"]) if f["sequence"] else "",
                "rep_five_offset": five_offset(f, best, length),
                "locus_type": best.get("locus_type", ""),
                "locus_count": best["total_count"],
                "five_score": best["five_score"], "three_score": best["three_score"],
                "annotation_class": best["annotation_class"],
            })
        row["observed"] = observed
        if f["expected"] == "fail":  # none of its own loci passes
            row["as_expected"] = int(not any(l["pass"] == "1" for l in own))
        else:
            row["as_expected"] = int(observed == f["expected"] or f["expected"] == "any")
        rows.append(row)

    false_pos = []
    for locus in loci:
        if locus["pass"] != "1":
            continue
        length, circular = contigs[(locus["genome"], locus["contig"])]
        same = [f for f in features if "_iv" in f and f["genome"] == locus["genome"]
                and f["contig"] == locus["contig"] and f["strand"] == locus["strand"]
                and overlaps(f["_iv"], locus["_iv"], length, circular)]
        if any(f["category"].startswith("srna") for f in same):
            continue
        sources = sorted(same, key=lambda f: -overlap_length(f["_iv"], locus["_iv"], length,
                                                             circular))
        false_pos.append({
            "locus_id": locus["locus_id"], "start": locus["start"], "end": locus["end"],
            "strand": locus["strand"], "span": locus["span"], "total_count": locus["total_count"],
            "rep_sequence": locus["rep_sequence"], "five_score": locus["five_score"],
            "annotation_class": locus["annotation_class"],
            "explained_by": ",".join(f["feature_id"] for f in sources) or "none",
        })
    return rows, false_pos, features


ROW_COLUMNS = ["feature_id", "category", "start", "end", "strand", "length", "expected",
               "observed", "as_expected", "total_reads", "locus_id", "locus_start", "locus_end",
               "locus_type", "extra_span", "rep_sequence", "rep_is_canonical", "rep_five_offset",
               "locus_count", "five_score", "three_score", "annotation_class"]
FP_COLUMNS = ["locus_id", "start", "end", "strand", "span", "total_count", "rep_sequence",
              "five_score", "annotation_class", "explained_by"]


def report(rows, false_pos, truth_dir, results_dir) -> str:
    lines = []
    srna = [r for r in rows if r["category"].startswith("srna")]
    ok = sum(r["as_expected"] for r in srna)
    lines.append(f"Simulated smallRNAs on focus genomes: {len(srna)}, "
                 f"outcome as expected: {ok} ({ok / max(1, len(srna)):.0%})")
    lines.append("")
    lines.append(f"{'category':<22}{'n':>4}  {'expected':<9} observed")
    by_cat: dict[tuple, Counter] = defaultdict(Counter)
    for r in rows:
        by_cat[(r["category"], r["expected"])][r["observed"]] += 1
    for (cat, exp), obs in sorted(by_cat.items()):
        detail = ", ".join(f"{k} {v}" for k, v in obs.most_common())
        lines.append(f"{cat:<22}{sum(obs.values()):>4}  {exp:<9} {detail}")

    found = [r for r in srna if r["observed"] in ("pass", "fail", "merged")]
    if found:
        canon = [r for r in found if r.get("rep_is_canonical") != ""]
        exact = sum(int(r["rep_is_canonical"]) for r in canon)
        embedded = sum(1 for r in found if int(r["extra_span"]) > 10)
        lines.append("")
        lines.append(f"Loci found: representative = canonical sequence in {exact}/{len(canon)}; "
                     f"5' offset of representative = 0 in "
                     f"{sum(1 for r in found if int(r['rep_five_offset']) == 0)}/{len(found)}; "
                     f"locus longer than smallRNA + 10 nt in {embedded}/{len(found)}")

    unexpected = [r for r in rows if not r["as_expected"]]
    if unexpected:
        lines.append("")
        lines.append("Not as expected:")
        for r in unexpected:
            extra = f" locus {r.get('locus_id', '')} span +{r.get('extra_span', '')}" \
                if r.get("locus_id") else ""
            lines.append(f"  {r['feature_id']:<22} expected {r['expected']:<7} "
                         f"observed {r['observed']:<7}{extra}")

    lines.append("")
    lines.append(f"Passing loci without a simulated smallRNA: {len(false_pos)}")
    for fp in false_pos:
        lines.append(f"  {fp['locus_id']} {fp['start']}-{fp['end']}({fp['strand']}) "
                     f"span {fp['span']} count {fp['total_count']} "
                     f"{fp['annotation_class']} <- {fp['explained_by']}")

    lines.append("")
    lines.append("Simulated composition (reads before trimming):")
    lines.extend("  " + "\t".join(r.values()) for r in read_tsv(f"{truth_dir}/composition.tsv"))
    lines.append("Observed origin classes:")
    lines.extend("  " + "\t".join(r.values())
                 for r in read_tsv(f"{results_dir}/origin/origin_summary.tsv"))
    return "\n".join(lines) + "\n"


def parse_genomic(text: str):
    """'chrN:10-17(+)' -> (contig, start, end), 1-based inclusive."""
    contig, rest = text.rsplit(":", 1)
    a, b = rest.split("(")[0].split("-")
    return contig, int(a), int(b)


def evaluate_targets(truth_dir: str, results_dir: str, rows: list[dict]):
    """Planted seed sites (truth/targets.tsv) found in module B's sites.tsv, per target set."""
    planted_path = Path(truth_dir) / "targets.tsv"
    site_files = sorted(Path(results_dir).glob("targets/*/sites.tsv"))
    if not planted_path.exists() or not site_files:
        return [], ""
    from .targets.pipeline import read_candidates
    from .targets.sites import CANONICAL
    by_rep = {seq: cid for cid, _, seq in
              read_candidates(f"{results_dir}/discovery/candidates.fasta")}
    candidate_of = {r["feature_id"]: by_rep.get(r.get("rep_sequence", "")) for r in rows
                    if r["observed"] == "pass"}
    sites = [s for path in site_files for s in read_tsv(str(path))]
    order = {c: i for i, c in enumerate(CANONICAL)}
    out = []
    for p in read_tsv(str(planted_path)):
        cid = candidate_of.get(p["feature_id"])
        start, end = int(p["genomic_start"]), int(p["genomic_end"])
        best = ""
        for s in sites:
            if s["candidate"] != cid or s["transcript"] != p["transcript"] or not s["genomic"]:
                continue
            _, a, b = parse_genomic(s["genomic"])
            if a <= end and b >= start and (not best or order[s["site_class"]] < order[best]):
                best = s["site_class"]
        found = bool(best) and order[best] <= order[p["site_class"]]
        out.append({**p, "candidate": cid or "", "found_class": best, "found": int(found)})
    n = sum(r["found"] for r in out)
    text = f"Planted target sites found with their class (or better): {n}/{len(out)}\n"
    text += "".join(f"  missing: {r['feature_id']} {r['site_class']} in {r['transcript']} "
                    f"(candidate {r['candidate'] or 'none'})\n" for r in out if not r["found"])
    return out, text


def run(truth_dir: str, results_dir: str, outdir: str) -> str:
    Path(outdir).mkdir(parents=True, exist_ok=True)
    rows, false_pos, _ = evaluate(truth_dir, results_dir)
    write_tsv(f"{outdir}/smallrnas.tsv", rows, ROW_COLUMNS)
    write_tsv(f"{outdir}/false_positives.tsv", false_pos, FP_COLUMNS)
    text = report(rows, false_pos, truth_dir, results_dir)
    planted, target_text = evaluate_targets(truth_dir, results_dir, rows)
    if planted:
        write_tsv(f"{outdir}/target_sites.tsv", planted, list(planted[0]))
        text += "\n" + target_text
    Path(f"{outdir}/summary.txt").write_text(text)
    return text
