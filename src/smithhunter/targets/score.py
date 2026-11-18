"""Module B step 2: energy scores for seed sites, decoy-based q-values and support.

Candidates and their decoys go through the identical scan and the identical predictors
(IntaRNA, RNAhybrid). Each unique (small RNA, region sequence, site) is scored once;
isoforms sharing a region sequence share the score and count once for the FDR.

For each predictor and site class, real and decoy standardised scores are pooled and each
real site gets a q-value: the lowest FDR among the thresholds it passes, with
FDR = ((decoy sites passing + 1) / decoys per candidate) / real sites passing (the + 1
keeps a q-value from being 0 when no decoy passes, which 20 decoys cannot show). A second
q-value uses only the candidate's own sites and decoys. ``support`` is the number of
predictors with pooled q <= ``fdr``. Nothing is filtered: every score stays in the table.

The score ranked by each predictor is its energy (IntaRNA total energy E, RNAhybrid MFE),
standardised against the candidate's own decoy sites (``<predictor>_z``) before pooling,
so that the composition of a small RNA does not decide its q-values. RNAhybrid's p-value
(calibrated per candidate on the target regions) is reported, not ranked.

The FDR measures how specific the procedure is, not the biology: a site that beats the
decoys is a real complementarity signal, not proof of regulation.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from bisect import bisect_right
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

from ..seqio import read_tsv, write_tsv
from . import intarna, rnahybrid
from .fdr import qvalues
from .regions import Region
from .sites import Site, scan, site_classes

WEAK = ("6mer", "offset-6mer")
PREDICTORS = ("intarna", "rnahybrid")
RANKED = {"intarna": "intarna_E", "rnahybrid": "rnahybrid_mfe"}
SCORE_COLUMNS = {
    "intarna": ["intarna_E", "intarna_z", "intarna_E_hybrid", "intarna_ED_target",
                "intarna_ED_query", "intarna_target", "intarna_guide", "intarna_duplex"],
    "rnahybrid": ["rnahybrid_mfe", "rnahybrid_z", "rnahybrid_p", "rnahybrid_target",
                  "rnahybrid_pairs", "rnahybrid_on_site", "rnahybrid_duplex"],
}


@dataclass(frozen=True)
class Job:
    """One site of one small RNA (real or decoy) on one unique region sequence."""
    guide_id: str
    candidate: str
    decoy: bool
    guide: str
    target: int  # index of the unique region sequence
    start: int
    end: int
    t1: int
    site_class: str

    @property
    def site(self) -> Site:
        return Site(self.start, self.end, self.t1, self.site_class)


def scored_classes(classes, score_weak: bool) -> tuple[str, ...]:
    return tuple(c for c in classes if score_weak or c not in WEAK)


def collect_jobs(candidates: list[tuple[str, str]], decoys: list[tuple[str, str, str]],
                 unique: list[str], seed=(2, 8), classes=None) -> list[Job]:
    """Sites of every candidate and decoy on the unique region sequences."""
    ids = {cid for cid, _ in candidates}
    guides = [(cid, cid, False, seq) for cid, seq in candidates]
    guides += [(did, cid, True, seq) for did, cid, seq in decoys if cid in ids]
    jobs = []
    for gid, cid, is_decoy, guide in guides:
        for i, seq in enumerate(unique):
            for s in scan(guide, seq, seed, classes):
                jobs.append(Job(gid, cid, is_decoy, guide, i, s.start, s.end, s.t1,
                                s.site_class))
    return jobs


def run_tool(cmd: list[str]) -> str:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed ({proc.returncode}): {proc.stderr.strip()[:500]}\n"
                           f"command: {' '.join(cmd)}")
    return proc.stdout


def require(binary: str) -> None:
    if not shutil.which(binary):
        sys.exit(f"{binary} not found on PATH. Create the environment from environment.yml "
                 "(or run snakemake with --sdm conda), or drop the predictor from "
                 "targets.predictors.")


def intarna_predictor(unique: list[str], seed, temperature):
    """Jobs -> scores. The accessibility of each region is computed once and sliced for
    each of its sites (see ``intarna``)."""
    require("IntaRNA")

    def predict_all(jobs: list[Job], threads: int = 1) -> list[dict]:
        by_target: dict[int, list[int]] = defaultdict(list)
        for k, job in enumerate(jobs):
            by_target[job.target].append(k)
        results: list[dict] = [{} for _ in jobs]
        with tempfile.TemporaryDirectory() as tmp, \
                ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            def table(t: int) -> list[str]:
                path = os.path.join(tmp, f"acc{t}")
                run_tool(intarna.accessibility_command(jobs[by_target[t][0]].guide, unique[t],
                                                       path, temperature))
                with open(path) as fh:
                    lines = fh.readlines()
                os.remove(path)
                return lines

            def one(k: int, lines: list[str]) -> dict:
                job = jobs[k]
                lo, hi = intarna.window(unique[job.target], job.site, len(job.guide), seed)
                path = os.path.join(tmp, f"job{k}")
                with open(path, "w") as fh:
                    fh.write(intarna.slice_accessibility(lines, lo, hi))
                cmd, offset = intarna.command(job.guide, unique[job.target], job.site, seed,
                                              temperature, accessibility=path)
                try:
                    return intarna.parse_csv(run_tool(cmd), offset)
                finally:
                    os.remove(path)

            # regions in batches, so the tables in memory stay bounded
            targets = list(by_target)
            for i in range(0, len(targets), 4 * max(1, threads)):
                batch = targets[i:i + 4 * max(1, threads)]
                tables = dict(zip(batch, pool.map(table, batch)))
                ks = [k for t in batch for k in by_target[t]]
                for k, res in zip(ks, pool.map(lambda k: one(k, tables[jobs[k].target]), ks)):
                    results[k] = res
        return results
    return predict_all


def calibrate(guides: dict[str, str], unique: list[str], n_windows: int,
              random_seed: int, threads: int) -> dict[str, tuple[float, float]]:
    """Extreme value parameters (xi, theta) of each small RNA, from RNAcalibrate on
    windows sampled from the target regions."""
    require("RNAcalibrate")
    if not guides:
        return {}
    longest = max(len(g) for g in guides.values())
    length = longest + 2 * rnahybrid.FLANK + 9
    windows = rnahybrid.calibration_windows(unique, length, n_windows, random_seed)
    if not windows or not guides:
        return {}
    names = sorted(guides)
    chunks = [names[i::threads] for i in range(threads) if names[i::threads]]
    with tempfile.TemporaryDirectory() as tmp:
        tfa = os.path.join(tmp, "windows.fa")
        with open(tfa, "w") as fh:
            fh.writelines(f">w{i}\n{w}\n" for i, w in enumerate(windows))

        def one(k: int) -> dict:
            qfa = os.path.join(tmp, f"q{k}.fa")
            with open(qfa, "w") as fh:
                fh.writelines(f">{n}\n{guides[n]}\n" for n in chunks[k])
            return rnahybrid.parse_calibrate(run_tool(rnahybrid.calibrate_command(
                tfa, qfa, max(len(w) for w in windows), longest)))
        with ThreadPoolExecutor(max_workers=threads) as pool:
            parts = list(pool.map(one, range(len(chunks))))
    return {k: v for part in parts for k, v in part.items()}


def rnahybrid_predictor(unique: list[str], seed, params: dict[str, tuple[float, float]] | None):
    """``params``: (xi, theta) per candidate, used for the candidate and for its decoys:
    a decoy has the length and dinucleotide composition of its candidate, and sharing the
    parameters maps the energies of both to p-values with the same function. None: no
    calibration, p-values are not reported."""
    require("RNAhybrid")

    def predict(job: Job) -> dict:
        win, offset = rnahybrid.window(unique[job.target], job.site, len(job.guide), seed)
        par = params.get(job.candidate) if params else None
        out = run_tool(rnahybrid.command(job.guide, win, job.site, seed, par))
        return rnahybrid.parse_compact(out, job.site, offset, seed, with_pvalue=par is not None)

    def predict_all(jobs: list[Job], threads: int = 1) -> list[dict]:
        out: list[dict] = []
        with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            for i in range(0, len(jobs), CHUNK):   # bounded number of pending calls
                out.extend(pool.map(predict, jobs[i:i + CHUNK]))
        return out
    return predict_all


def score_jobs(jobs: list[Job], predictors: dict, threads: int = 1) -> dict[Job, dict]:
    """Run every predictor (a function (jobs, threads) -> one dict of scores per job) on
    every job; the scores of all predictors are merged per job."""
    scores: dict[Job, dict] = {j: {} for j in jobs}
    for predict_all in predictors.values():
        for job, result in zip(jobs, predict_all(jobs, threads)):
            scores[job].update(result)
    return scores


MIN_DECOY_SCORES = 20
CHUNK = 10000


def standardise(jobs: list[Job], scores: dict[Job, dict], predictor: str, column: str) -> None:
    """Add ``<predictor>_z``: the score in standard deviations from the mean of the
    candidate's decoy sites (all classes). Lower is better, as for energies.

    Decoys have the length and dinucleotide composition of their candidate, so this takes
    out what composition does to the energies (a GC-rich small RNA pairs more strongly
    with anything). Without it, pooling raw energies over candidates lets the decoys of
    GC-rich candidates hide the true sites of AT-rich ones.

    A real site is measured against all of its candidate's decoys; a decoy site against
    the other decoys only (leave one decoy out). Otherwise the decoys would set their own
    scale, their z-scores could not reach the tail that a real site of no interest
    reaches, and the q-values would be too low. A candidate with fewer than
    ``MIN_DECOY_SCORES`` decoy scores once any one decoy is left out gets no z (neither its
    real nor its decoy sites), nor does a site whose reference scores are all equal.
    """
    # sums per candidate and per decoy, so that leaving one decoy out is a subtraction
    total: dict[str, list[float]] = defaultdict(lambda: [0, 0.0, 0.0])
    own: dict[str, list[float]] = defaultdict(lambda: [0, 0.0, 0.0])
    for j in jobs:
        value = scores[j].get(column, "")
        if j.decoy and value != "":
            for acc in (total[j.candidate], own[j.guide_id]):
                acc[0] += 1
                acc[1] += value
                acc[2] += value * value
    # a candidate is standardised only if every one of its decoys can be, so that its
    # real sites never face a pool without its decoys
    largest: dict[str, float] = defaultdict(float)
    for j in jobs:
        if j.decoy and j.guide_id in own:
            largest[j.candidate] = max(largest[j.candidate], own[j.guide_id][0])
    eligible = {c for c, (n, _, _) in total.items() if n - largest[c] >= MIN_DECOY_SCORES}
    for j in jobs:
        value = scores[j].get(column, "")
        if value == "" or j.candidate not in eligible:
            continue
        n, s1, s2 = total[j.candidate]
        if j.decoy:
            n, s1, s2 = n - own[j.guide_id][0], s1 - own[j.guide_id][1], s2 - own[j.guide_id][2]
        if n < MIN_DECOY_SCORES:
            continue
        m = s1 / n
        var = s2 / n - m * m
        if var <= 1e-12:
            continue
        scores[j][f"{predictor}_z"] = round((value - m) / var ** 0.5, 4)


def add_qvalues(jobs: list[Job], scores: dict[Job, dict], ranked: dict[str, str],
                n_decoys: int, fdr: float) -> list[dict]:
    """Add ``<predictor>_z``, ``<predictor>_q``, ``<predictor>_q_candidate`` and
    ``support`` to the scores of the real sites; return the FDR curves (real and decoy
    counts at each threshold).

    ``_q`` pools candidates on the standardised score (``standardise``), per site class.
    ``_q_candidate`` uses the candidate's own decoys only, on the raw score."""
    curves = []
    by_class: dict[str, list[Job]] = defaultdict(list)
    for j in jobs:
        by_class[j.site_class].append(j)
    for predictor, column in ranked.items():
        z = f"{predictor}_z"
        standardise(jobs, scores, predictor, column)
        for cls, group in by_class.items():
            have = [j for j in group if scores[j].get(z, "") != ""]
            real = [j for j in have if not j.decoy]
            decoy = [scores[j][z] for j in have if j.decoy]
            q = qvalues([scores[j][z] for j in real], decoy, n_decoys, pseudo=1)
            for j, value in zip(real, q):
                scores[j][f"{predictor}_q"] = value
            per_candidate: dict[str, list[Job]] = defaultdict(list)
            for j in group:
                if scores[j].get(column, "") != "":
                    per_candidate[j.candidate].append(j)
            for members in per_candidate.values():
                mine = [j for j in members if not j.decoy]
                theirs = [scores[j][column] for j in members if j.decoy]
                for j, value in zip(mine, qvalues([scores[j][column] for j in mine], theirs,
                                                  n_decoys, pseudo=1)):
                    scores[j][f"{predictor}_q_candidate"] = value
            real_scores = sorted(scores[j][z] for j in real)
            decoy.sort()
            for threshold in sorted(set(real_scores)):
                n_real = bisect_right(real_scores, threshold)
                n_decoy = bisect_right(decoy, threshold)
                expected = (n_decoy + 1) / max(1, n_decoys)
                curves.append({"predictor": predictor, "score": z, "site_class": cls,
                               "threshold": threshold, "real_sites": n_real,
                               "decoy_sites": n_decoy, "expected_false": round(expected, 3),
                               "fdr": round(min(1.0, expected / n_real), 4)})
    for j in jobs:
        if not j.decoy and scores[j]:
            scores[j]["support"] = sum(
                1 for p in ranked if scores[j].get(f"{p}_q", "") != ""
                and scores[j][f"{p}_q"] <= fdr)
    return curves


def run_scoring(candidates: list[tuple[str, str]], decoys: list[tuple[str, str, str]],
                regions: list[Region], site_rows: list[dict], predictors: dict,
                ranked: dict[str, str], seed=(2, 8), classes=None, score_weak: bool = False,
                fdr: float = 0.05, threads: int = 1):
    """Score the sites of step 1 (``site_rows``, from sites.tsv) and their decoy sites.

    ``predictors`` maps a name to a function that, given the unique region sequences,
    returns the predictor (jobs, threads) -> scores; it is a parameter so that tests can run without the external tools.
    Returns (scored site rows, target rows, decoy rows, FDR curves, columns).
    """
    classes = tuple(classes or site_classes(seed))
    wanted = scored_classes(classes, score_weak)
    unique = list(dict.fromkeys(r.sequence for r in regions))
    index = {seq: i for i, seq in enumerate(unique)}
    jobs = [j for j in collect_jobs(candidates, decoys, unique, seed, classes)
            if j.site_class in wanted]
    n_decoys = max((sum(1 for _, c, _ in decoys if c == cid) for cid, _ in candidates),
                   default=0)
    scores = score_jobs(jobs, {name: make(unique) for name, make in predictors.items()},
                        threads)
    curves = add_qvalues(jobs, scores, ranked, n_decoys, fdr)

    columns = [c for p in predictors for c in SCORE_COLUMNS.get(p, [])]
    q_columns = [c for p in ranked for c in (f"{p}_q", f"{p}_q_candidate")] + ["support"]
    region_seq = {r.id: r.sequence for r in regions}
    real = {(j.guide_id, j.target, j.start, j.site_class): scores[j]
            for j in jobs if not j.decoy}
    order = {c: i for i, c in enumerate(classes)}
    out = []
    for row in site_rows:
        seq = region_seq.get(row["region_id"])
        key = (row["candidate"], index.get(seq), int(row["site_start"]) - 1, row["site_class"])
        s = real.get(key, {})
        out.append({**row, **{c: s.get(c, "") for c in columns + q_columns}})

    def rank(r):
        e = r.get("intarna_E", "")
        return (r["candidate"], order.get(r["site_class"], len(order)),
                -(r["support"] if r["support"] != "" else -1), e if e != "" else 0.0)
    out.sort(key=rank)

    decoy_rows = [{"decoy": j.guide_id, "candidate": j.candidate, "region_sequence": j.target,
                   "site_class": j.site_class, "site_start": j.start + 1,
                   **{c: scores[j].get(c, "") for c in columns}}
                  for j in jobs if j.decoy]
    return out, summarise_targets(out, classes, ranked), decoy_rows, curves, columns + q_columns


def summarise_targets(rows: list[dict], classes, ranked: dict[str, str]) -> list[dict]:
    """One row per candidate x transcript: sites per class, best scores, best q, support."""
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        groups[(r["candidate"], r["transcript"], r["gene"])].append(r)
    order = {c: i for i, c in enumerate(classes)}
    out = []
    for (cid, transcript, gene), members in groups.items():
        def best(column):
            values = [m[column] for m in members if m.get(column, "") != ""]
            return min(values) if values else ""
        support = [m["support"] for m in members if m.get("support", "") != ""]
        row = {"candidate": cid, "transcript": transcript, "gene": gene, "sites": len(members),
               "best_class": min((m["site_class"] for m in members),
                                 key=lambda c: order.get(c, len(order))),
               **{f"n_{c}": sum(1 for m in members if m["site_class"] == c) for c in classes}}
        for p, column in ranked.items():
            row[f"best_{column}"] = best(column)
            row[f"best_{p}_q"] = best(f"{p}_q")
        row["max_support"] = max(support) if support else ""
        row["sites_supported"] = sum(1 for s in support if ranked and s == len(ranked))
        out.append(row)
    out.sort(key=lambda r: (r["candidate"], -(r["max_support"] if r["max_support"] != "" else -1),
                            order.get(r["best_class"], len(order))))
    return out


CURVE_COLUMNS = ["predictor", "score", "site_class", "threshold", "real_sites", "decoy_sites",
                 "expected_false", "fdr"]


def write_scored(site_rows, target_rows, decoy_rows, curves, score_columns, classes, ranked,
                 sites_path, targets_path, decoys_path, curves_path) -> None:
    from .pipeline import SITE_COLUMNS
    write_tsv(sites_path, site_rows, SITE_COLUMNS + score_columns)
    target_columns = (["candidate", "transcript", "gene", "sites", "best_class"]
                      + [f"n_{c}" for c in classes]
                      + [c for p, col in ranked.items() for c in (f"best_{col}", f"best_{p}_q")]
                      + ["max_support", "sites_supported"])
    write_tsv(targets_path, target_rows, target_columns)
    decoy_columns = ["decoy", "candidate", "region_sequence", "site_class", "site_start"] + \
        [c for c in score_columns if not c.endswith(("_q", "_q_candidate")) and c != "support"]
    write_tsv(decoys_path, decoy_rows, decoy_columns)
    write_tsv(curves_path, curves, CURVE_COLUMNS)


SECONDS_PER_JOB = 0.045   # IntaRNA + RNAhybrid per site, measured on the simulation
SECONDS_PER_KB = 0.7      # IntaRNA accessibility table, per kb of region with sites


def estimate(candidates, decoys, regions, seed=(2, 8), classes=None,
             score_weak=False) -> str:
    """Sites to score per class (real and decoy) and a rough single-core run time."""
    classes = tuple(classes or site_classes(seed))
    wanted = scored_classes(classes, score_weak)
    unique = list(dict.fromkeys(r.sequence for r in regions))
    jobs = collect_jobs(candidates, decoys, unique, seed, classes)
    lines = ["class\treal_sites\tdecoy_sites\tscored"]
    for c in classes:
        real = sum(1 for j in jobs if j.site_class == c and not j.decoy)
        decoy = sum(1 for j in jobs if j.site_class == c and j.decoy)
        lines.append(f"{c}\t{real}\t{decoy}\t{'yes' if c in wanted else 'no'}")
    scored = [j for j in jobs if j.site_class in wanted]
    kb = sum(len(unique[t]) for t in {j.target for j in scored}) / 1000
    hours = (len(scored) * SECONDS_PER_JOB + kb * SECONDS_PER_KB) / 3600
    lines.append(f"{len(scored)} sites to score on {kb:.0f} kb of regions: about "
                 f"{hours:.1f} core-hours with both predictors (rough)")
    return "\n".join(lines)


def run(candidates, decoys, regions, sites_path: str, wanted: list[str], seed=(2, 8),
        classes=None, score_weak=False, fdr=0.05, threads=1, temperature=None,
        calibration_windows=2000, random_seed=1, log=print):
    """Entry point of ``score-sites``: build the predictors and score everything."""
    unknown = set(wanted) - set(PREDICTORS)
    if unknown:
        sys.exit(f"unknown predictors {sorted(unknown)}; available: {', '.join(PREDICTORS)}")
    classes = tuple(classes or site_classes(seed))
    ranked = {p: RANKED[p] for p in PREDICTORS if p in wanted}
    predictors = {}
    if "intarna" in wanted:
        predictors["intarna"] = lambda unique: intarna_predictor(unique, seed, temperature)
    if "rnahybrid" in wanted:
        params = None
        if calibration_windows:
            unique = list(dict.fromkeys(r.sequence for r in regions))
            guides = dict(candidates)
            params = calibrate(guides, unique, calibration_windows, random_seed,
                               max(1, threads))
            failed = sorted(set(guides) - set(params))
            log(f"RNAcalibrate: parameters for {len(params)} of {len(guides)} candidates "
                f"({calibration_windows} windows of the target regions)")
            if failed:
                log(f"  no fit (no p-value; MFE and q-values are kept) for: "
                    f"{', '.join(failed[:20])}"
                    + (" ..." if len(failed) > 20 else ""))
        else:
            log("RNAhybrid calibration is off: no p-values")
        predictors["rnahybrid"] = lambda unique, p=params: rnahybrid_predictor(unique, seed, p)
    result = run_scoring(candidates, decoys, regions, read_tsv(sites_path), predictors, ranked,
                         seed, classes, score_weak, fdr, threads)
    return (*result, classes, ranked)

