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

``as_expected`` compares it with the outcome the simulator planned for the genome-choice
mode of the run (``assignment.mode``: ``best``, ``severe`` or ``evil``). The mode, the
``evil`` tolerance and the ``ambiguous`` policy are read from the configuration the run
used (``results/report/provenance.yaml``, else ``config/config.yaml`` next to the results)
and select the ``expected_<mode>`` column of ``truth/features.tsv``. The reads of every
simulated smallRNA are also followed through the origin classes (``origin_by_feature.tsv``),
so that reads and loci removed by ``evil`` (class ``evil_rejected``) are reported. For loci found
the table reports how far they extend beyond the smallRNA (``extra_span``, large when
degradation fragments chain into the locus), whether the representative sequence is
the canonical one, and the 5' offset of the representative read. Passing loci that
contain no simulated smallRNA are listed as false positives, with the source that
explains them (degradation, tRNA fragments).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from .seqio import open_text, read_tsv, write_tsv

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


def read_yaml_scalars(path: str) -> dict[str, str]:
    """Scalars of a block-style YAML file as {'a.b.c': 'text'} (enough for the configs
    written by ``simulate`` and for ``provenance.yaml``; no third-party parser)."""
    out: dict[str, str] = {}
    stack: list[tuple[int, str]] = []
    with open(path) as fh:
        for line in fh:
            text = line.split(" #")[0].rstrip()
            body = text.lstrip()
            if not body or body.startswith(("#", "- ")) or ":" not in body:
                continue
            indent = len(text) - len(body)
            key, _, value = body.partition(":")
            while stack and stack[-1][0] >= indent:
                stack.pop()
            stack.append((indent, key.strip().strip('"')))
            value = value.strip()
            if value:
                out[".".join(k for _, k in stack)] = value.strip('"')
    return out


def run_settings(results_dir: str) -> dict:
    """Genome-choice settings the run used: mode, evil_max_diff, ambiguous, source."""
    settings = {"mode": "best", "evil_max_diff": 0, "ambiguous": "report",
                "source": "defaults (no configuration found)"}
    candidates = [(Path(results_dir) / "report/provenance.yaml", "config."),
                  (Path(results_dir).resolve().parent / "config/config.yaml", "")]
    for path, prefix in candidates:
        if not path.exists():
            continue
        values = read_yaml_scalars(str(path))
        settings["mode"] = values.get(f"{prefix}assignment.mode", "best")
        settings["evil_max_diff"] = int(values.get(f"{prefix}assignment.evil_max_diff", 0))
        settings["ambiguous"] = values.get(f"{prefix}ambiguous", "report")
        settings["source"] = str(path)
        break
    return settings


def expected_column(settings: dict) -> str:
    """Column of truth/features.tsv holding the outcome expected under the run's mode."""
    mode = settings["mode"]
    if mode == "evil" and settings["evil_max_diff"] == 0:
        mode = "severe"  # no tolerance: evil is severe
    return f"expected_{mode}"


def evaluate(truth_dir: str, results_dir: str, settings: dict | None = None):
    features = read_tsv(f"{truth_dir}/features.tsv")
    column = expected_column(settings or run_settings(results_dir))
    for f in features:  # truth written before the modes existed has only `expected`
        f["expected"] = f.get(column) or f["expected"]
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


ORIGIN_CLASSES = ("assigned", "ambiguous", "evil_rejected", "excluded", "multimap_dropped",
                  "too_many_hits", "unmapped")
SEED = 8  # reads are matched to a smallRNA through a shared 8-mer, then compared base by base
MIN_OVERLAP = 15
MAX_DIFF = 2


def origin_by_feature(features: list[dict], results_dir: str) -> list[dict]:
    """Reads of each simulated smallRNA of a focus genome per origin class.

    The classes of the reads carry no position (an ``evil_rejected`` read is placed
    nowhere), so reads are matched to smallRNAs by sequence: at least ``MIN_OVERLAP``
    bases aligned without gaps to the canonical sequence with at most ``MAX_DIFF``
    differences. Reads cut by fastp below the length range or with more errors are not
    seen, so the classes add up to a little less than the simulated reads; degradation
    fragments or tRNA fragments covering the same bases are counted with the smallRNA.
    """
    origin_path = Path(results_dir) / "origin/read_origin.tsv.gz"
    counts_path = Path(results_dir) / "collapsed/counts.tsv.gz"
    targets = [f for f in features if f["category"].startswith("srna") and f.get("sequence")]
    if not targets or not origin_path.exists() or not counts_path.exists():
        return []
    seeds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for i, f in enumerate(targets):
        seq = f["sequence"]
        for j in range(len(seq) - SEED + 1):
            seeds[seq[j:j + SEED]].append((i, j))

    def owner(read: str):
        best = None
        tried = set()
        for k in range(len(read) - SEED + 1):
            for i, j in seeds.get(read[k:k + SEED], ()):
                shift = j - k  # read position + shift = position in the smallRNA
                if (i, shift) in tried:
                    continue
                tried.add((i, shift))
                seq = targets[i]["sequence"]
                lo, hi = max(0, -shift), min(len(read), len(seq) - shift)
                diff = sum(read[x] != seq[x + shift] for x in range(lo, hi))
                if hi - lo >= MIN_OVERLAP and diff <= MAX_DIFF:
                    key = (diff - (hi - lo), i)
                    best = min(best, key) if best else key
        return best[1] if best else None

    with open_text(str(origin_path)) as fh:
        fh.readline()
        classes = dict(line.rstrip("\n").split("\t")[:2] for line in fh)
    table = [Counter() for _ in targets]
    with open_text(str(counts_path)) as fh:
        fh.readline()
        for line in fh:
            rid, seq, _, *per_sample = line.rstrip("\n").split("\t")
            i = owner(seq)
            if i is not None:
                table[i][classes.get(rid, "unmapped")] += sum(int(x) for x in per_sample)
    rows = []
    for f, c in zip(targets, table):
        rows.append({"feature_id": f["feature_id"], "category": f["category"],
                     "genome": f["genome"], "simulated_reads": f["total_reads"],
                     "matched_reads": sum(c.values()),
                     **{cls: c.get(cls, 0) for cls in ORIGIN_CLASSES},
                     "other": sum(v for k, v in c.items() if k not in ORIGIN_CLASSES)})
    return rows


ORIGIN_COLUMNS = ["feature_id", "category", "genome", "simulated_reads", "matched_reads",
                  *ORIGIN_CLASSES, "other"]


def mode_report(rows: list[dict], origin_rows: list[dict], settings: dict, truth_dir: str,
                results_dir: str) -> str:
    """Genome-choice mode of the run, what it did to the NUMT scenarios, evil_rejected."""
    mode = settings["mode"]
    detail = f" (evil_max_diff {settings['evil_max_diff']})" if mode == "evil" else ""
    lines = [f"Genome-choice mode of the run: {mode}{detail}, ambiguous policy "
             f"{settings['ambiguous']} [{settings['source']}]",
             f"Expected outcomes taken from column {expected_column(settings)}"]
    assumed_path = Path(truth_dir) / "assignment.tsv"
    if assumed_path.exists():
        assumed = {r["setting"]: r["value"] for r in read_tsv(str(assumed_path))}
        if assumed.get("ambiguous", settings["ambiguous"]) != settings["ambiguous"]:
            lines.append(f"WARNING: the truth assumes ambiguous policy {assumed['ambiguous']}, "
                         f"the run used {settings['ambiguous']}: outcomes of ties between "
                         "genomes in mode best may differ from the expected ones")

    observed = {r["feature_id"]: r for r in rows}
    shown = [o for o in origin_rows if "numt" in o["category"] or o["evil_rejected"]]
    if shown:
        lines.append("")
        lines.append("Reads of smallRNAs with a copy in an exclude genome, by origin class "
                     "(all samples; every smallRNA is in origin_by_feature.tsv):")
        lines.append(f"  {'feature':<21}{'expected':<9}{'observed':<9}{'reads':>6}"
                     f"{'assigned':>9}{'ambiguous':>10}{'evil_rejected':>14}{'excluded':>9}"
                     f"{'other':>6}")
        for o in shown:
            r = observed.get(o["feature_id"], {})
            other = o["matched_reads"] - sum(o[c] for c in ("assigned", "ambiguous",
                                                            "evil_rejected", "excluded"))
            lines.append(f"  {o['feature_id']:<21}{r.get('expected', ''):<9}"
                         f"{r.get('observed', ''):<9}{o['matched_reads']:>6}{o['assigned']:>9}"
                         f"{o['ambiguous']:>10}{o['evil_rejected']:>14}{o['excluded']:>9}"
                         f"{other:>6}")

    summary_path = Path(results_dir) / "origin/origin_summary.tsv"
    total = 0
    if summary_path.exists():
        for r in read_tsv(str(summary_path)):
            if r["class"] == "evil_rejected":
                total += sum(int(v) for k, v in r.items() if k not in ("class", "genomes"))
    by_feature = sorted(((o["evil_rejected"], o["feature_id"]) for o in origin_rows
                         if o["evil_rejected"]), reverse=True)
    explained = sum(n for n, _ in by_feature)
    lines.append("")
    if mode != "evil" and not total:
        lines.append("evil_rejected reads: 0 (the class exists only in mode evil)")
    else:
        lines.append(f"evil_rejected reads: {total} in the run, {explained} from simulated "
                     f"smallRNAs, {total - explained} from other sources"
                     + ("" if mode == "evil" else f"  (unexpected in mode {mode})"))
        for n, fid in by_feature:
            r = observed.get(fid, {})
            lines.append(f"  {fid:<21}{n:>6} reads  expected {r.get('expected', '?')}, "
                         f"observed {r.get('observed', '?')}")
        # lost to evil: most of its reads were rejected (a stray rejected read is not a cause)
        mostly = {o["feature_id"] for o in origin_rows
                  if 2 * o["evil_rejected"] > o["matched_reads"]}
        lost = [r["feature_id"] for r in rows if r["feature_id"] in mostly
                and r["observed"] != "pass"]
        lines.append("smallRNAs without a passing locus after evil rejection: "
                     + (", ".join(lost) if lost else "none"))
    return "\n".join(lines) + "\n"


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


def candidate_ids(results_dir: str, rows: list[dict]) -> dict[str, str]:
    """Simulated smallRNA -> module A candidate id (through the representative read)."""
    from .targets.pipeline import read_candidates
    path = Path(results_dir) / "discovery" / "candidates.fasta"
    if not path.exists():
        return {}
    by_rep = {seq: cid for cid, _, seq in read_candidates(str(path))}
    return {r["feature_id"]: by_rep.get(r.get("rep_sequence", "")) for r in rows
            if r["observed"] == "pass"}


def matching_site(sites: list[dict], cid, planted: dict, order: dict) -> dict | None:
    """The site of ``cid`` with the best class overlapping the planted site."""
    start, end = int(planted["genomic_start"]), int(planted["genomic_end"])
    best = None
    for s in sites:
        if s["candidate"] != cid or s["transcript"] != planted["transcript"] or not s["genomic"]:
            continue
        _, a, b = parse_genomic(s["genomic"])
        if a <= end and b >= start and (best is None or
                                        order[s["site_class"]] < order[best["site_class"]]):
            best = s
    return best


def evaluate_targets(truth_dir: str, results_dir: str, rows: list[dict]):
    """Planted target sites (truth/targets.tsv) found in module B's sites.tsv, and, when
    the sites were scored (scored_sites.tsv), their support and the share of calls
    (support >= 1) that are not planted sites."""
    planted_path = Path(truth_dir) / "targets.tsv"
    site_files = sorted(Path(results_dir).glob("targets/*/sites.tsv"))
    if not planted_path.exists() or not site_files:
        return [], ""
    from .targets.sites import CANONICAL
    candidate_of = candidate_ids(results_dir, rows)
    sites = [s for path in site_files for s in read_tsv(str(path))]
    scored_files = sorted(Path(results_dir).glob("targets/*/scored_sites.tsv"))
    scored = [s for path in scored_files for s in read_tsv(str(path))]
    order = {c: i for i, c in enumerate(CANONICAL)}
    out = []
    for p in read_tsv(str(planted_path)):
        cid = candidate_of.get(p["feature_id"])
        hit = matching_site(sites, cid, p, order)
        best = hit["site_class"] if hit else ""
        found = bool(best) and order[best] <= order[p["site_class"]]
        row = {**p, "pairing": p.get("pairing", "seed"), "candidate": cid or "",
               "found_class": best, "found": int(found)}
        if scored:
            s = matching_site(scored, cid, p, order)
            for col in ("intarna_E", "intarna_q", "rnahybrid_p", "rnahybrid_q", "support"):
                row[col] = s.get(col, "") if s else ""
        out.append(row)
    n = sum(r["found"] for r in out)
    text = f"Planted target sites found with their class (or better): {n}/{len(out)}\n"
    text += "".join(f"  missing: {r['feature_id']} {r['site_class']} in {r['transcript']} "
                    f"(candidate {r['candidate'] or 'none'})\n" for r in out if not r["found"])
    if scored:
        text += scoring_report(out, scored)
    return out, text


def scoring_report(planted: list[dict], scored: list[dict]) -> str:
    """Support of the planted sites by kind of pairing, and how many calls are false.

    Seed-only sites are not expected to stand out: a seed match is as frequent for the
    decoys. Sites with pairing beyond the seed should. Every call (support >= 1) that is
    not a planted site is false, as the simulated transcripts are random sequence."""
    text = "Planted target sites by pairing, support >= 1 / 2 (q <= fdr):\n"
    for kind in ("extensive", "supplementary", "seed"):
        group = [r for r in planted if r["pairing"] == kind]
        if group:
            s1 = sum(1 for r in group if r["support"] not in ("", None) and int(r["support"]) >= 1)
            s2 = sum(1 for r in group if r["support"] not in ("", None) and int(r["support"]) >= 2)
            text += f"  {kind:<14} {s1}/{len(group)} {s2}/{len(group)}\n"
    planted_keys = {(r["candidate"], r["transcript"]) for r in planted}
    calls = [s for s in scored if s.get("support", "") not in ("", None) and int(s["support"]) >= 1]
    true = [s for s in calls if (s["candidate"], s["transcript"]) in planted_keys
            and any(r["candidate"] == s["candidate"] and r["transcript"] == s["transcript"]
                    and r.get("support") not in ("", None) and int(r["support"]) >= 1
                    and _same_site(r, s) for r in planted)]
    n_false = len(calls) - len(true)
    share = f"{n_false / len(calls):.2f}" if calls else "NA"
    text += (f"Scored sites called (support >= 1): {len(calls)}, of which not planted: "
             f"{n_false} (false discovery proportion {share})\n")
    return text


def _same_site(planted: dict, site: dict) -> bool:
    if not site.get("genomic"):
        return False
    _, a, b = parse_genomic(site["genomic"])
    return a <= int(planted["genomic_end"]) and b >= int(planted["genomic_start"])


def evaluate_precursors(truth_dir: str, results_dir: str, rows: list[dict]):
    """Planted hairpins (truth/hairpins.tsv) against precursors.tsv: the best window of the
    smallRNA's candidate should be a significant single hairpin with the smallRNA in the
    planted arm; the other candidates give the background rate of p <= 0.05."""
    truth = Path(truth_dir) / "hairpins.tsv"
    table = Path(results_dir) / "precursors" / "precursors.tsv"
    if not truth.exists() or not table.exists():
        return [], ""
    windows = read_tsv(str(table))
    best = {w["candidate"]: w for w in windows if w.get("best") == "1"}
    candidate_of = candidate_ids(results_dir, rows)
    out = []
    for h in read_tsv(str(truth)):
        cid = candidate_of.get(h["feature_id"])
        w = best.get(cid, {})
        p = float(w["pvalue"]) if w.get("pvalue") not in (None, "", "NA") else None
        ok = (p is not None and p <= 0.05 and w.get("single_hairpin") == "1"
              and w.get("arm") == h["arm"])
        out.append({**h, "candidate": cid or "", "window": w.get("window", ""),
                    "pvalue": w.get("pvalue", ""), "mfei": w.get("mfei", ""),
                    "single_hairpin": w.get("single_hairpin", ""), "arm_found": w.get("arm", ""),
                    "found": int(ok)})
    planted = {r["candidate"] for r in out}
    others = [w for c, w in best.items() if c not in planted and w.get("pvalue") not in ("", "NA")]
    low = sum(1 for w in others if float(w["pvalue"]) <= 0.05)
    n = sum(r["found"] for r in out)
    text = (f"Planted hairpins found (best window p <= 0.05, single hairpin, smallRNA in the "
            f"planted arm): {n}/{len(out)}\n")
    text += "".join(f"  missing: {r['feature_id']} (candidate {r['candidate'] or 'none'}, "
                    f"window {r['window']}, p {r['pvalue']}, arm {r['arm_found']})\n"
                    for r in out if not r["found"])
    if others:
        text += (f"Other candidates with best-window p <= 0.05: {low}/{len(others)} "
                 f"(about 0.1 expected: the best of two windows)\n")
    return out, text


def run(truth_dir: str, results_dir: str, outdir: str) -> str:
    Path(outdir).mkdir(parents=True, exist_ok=True)
    settings = run_settings(results_dir)
    rows, false_pos, features = evaluate(truth_dir, results_dir, settings)
    write_tsv(f"{outdir}/smallrnas.tsv", rows, ROW_COLUMNS)
    write_tsv(f"{outdir}/false_positives.tsv", false_pos, FP_COLUMNS)
    text = report(rows, false_pos, truth_dir, results_dir)
    focus = {g["name"] for g in read_tsv(f"{results_dir}/reference/genomes.tsv")
             if g["role"] == "focus"}
    origin_rows = origin_by_feature([f for f in features if f["genome"] in focus], results_dir)
    if origin_rows:
        write_tsv(f"{outdir}/origin_by_feature.tsv", origin_rows, ORIGIN_COLUMNS)
    text += "\n" + mode_report(rows, origin_rows, settings, truth_dir, results_dir)
    planted, target_text = evaluate_targets(truth_dir, results_dir, rows)
    if planted:
        write_tsv(f"{outdir}/target_sites.tsv", planted, list(planted[0]))
        text += "\n" + target_text
    hairpins, hairpin_text = evaluate_precursors(truth_dir, results_dir, rows)
    if hairpins:
        write_tsv(f"{outdir}/hairpins.tsv", hairpins, list(hairpins[0]))
        text += "\n" + hairpin_text
    Path(f"{outdir}/summary.txt").write_text(text)
    return text
