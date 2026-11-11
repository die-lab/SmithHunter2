"""Module B step 2: score seed sites with IntaRNA and RNAhybrid, q-values from decoys.

Every site of every real candidate and of its decoys is scored identically (the decoys
are dinucleotide-shuffled copies, see ``shuffle.py``). Scores are kept as they come out of
the tools; trust comes from the decoys: for each predictor and site class the q-value of
a real site is the lowest FDR over the thresholds it passes, with

    FDR(s) = (decoy sites scoring <= s / decoys per candidate) / (real sites scoring <= s)

computed pooled over candidates and, alongside, for the candidate alone (``*_q_cand``).
``support`` is the number of counted predictors whose pooled q-value is at most ``fdr``; it
is a count of agreeing predictors, not a confidence. Calls are the UNION: a site is called
when ANY counted predictor passes (``called_by`` lists them, ``best_q`` is the smallest q), so
a predictor that does not suit the data never removes a call. Which predictors are computed
and which count is policy, see ``predictors.py`` (``predictors``, ``support_from``).

Two-stage mode (``intarna_prefilter`` < 1, both predictors): RNAhybrid, ~15x cheaper, scores
every site first; IntaRNA then scores the sites whose MFE is in the best share of each class,
the cutoff being one pooled threshold per class over real and decoy sites together. Real and
decoy sites meet the same rule, so the q-values stay valid; they are conditional on passing it,
and sites below the cutoff have no IntaRNA values (``prefilter_keep_classes`` are always
scored by IntaRNA in full; the sites IntaRNA never saw are counted in ``prefilter_summary.tsv``).

Identical region sequences (isoforms sharing a UTR) are scored once and reported for
every region that carries them. Only the best class of each site is scored, and by
default only the canonical strong classes (8mer, 7mer-m8, 7mer-A1): the 6mers are far
more frequent and cost most of the time (``score_weak_sites``).
"""

from __future__ import annotations

import hashlib
import os
import pickle
import shutil
import statistics
import time
from array import array
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field, replace

from ..seqio import write_tsv
from . import intarna, rnahybrid
from . import predictors as pred
from .fdr import qvalues
from .regions import Region
from .sites import Site, context, scan, site_classes

CALIBRATION_LENGTH = (500, 300)  # mean, sd of the random sequences of RNAcalibrate
STRONG = ("8mer", "7mer-m8", "7mer-A1")
CLASS_RANK = {c: i for i, c in enumerate(("8mer", "7mer-m8", "7mer-A1", "6mer", "offset-6mer"))}
PREDICTORS = {"intarna": "intarna_E", "rnahybrid": "rnahybrid_mfe"}


@dataclass
class ScoreConfig:
    predictors: tuple[str, ...] = ("intarna", "rnahybrid")
    score_weak_sites: bool = False
    temperature: float = 37.0
    intarna_flank: int = 100
    intarna_prefilter: float = 1.0  # <1: IntaRNA only for this best share of sites by RNAhybrid MFE
    rnahybrid_flank: int = 15
    calibrate_samples: int = 1000   # RNAcalibrate sample size; 0 = use the 3utr_human preset
    workers: int = 1
    fdr: float = 0.05
    max_error_share: float = 0.05   # stop the run when more sites than this fail
    timeout: float = 300            # seconds per tool call; a hang counts as an error of the site
    resume: bool = True             # reuse checkpoints in the work directory when they match
    fdr_pseudocount: bool = True    # (D + 1) / T: one pseudo decoy site in every FDR estimate
    rnahybrid_require_site: bool = True  # RNAhybrid duplex must cover the seed of this site
    intarna_bin: str = "IntaRNA"
    rnahybrid_bin: str = "RNAhybrid"
    rnacalibrate_bin: str = "RNAcalibrate"
    intarna_extra: list[str] = field(default_factory=list)
    rnahybrid_batch: int = 2000     # sites per RNAhybrid process (0: one process per site)
    intarna_batch: int = 50         # windows per IntaRNA process (0: one process per site);
    # small RNAs that share a window are scored in one process too. Both give the same
    # numbers as one process per site, so neither is part of the run signature.
    # predictor policy (predictors.py)
    mode: str = "ago"                       # ago | bacterial | cross-kingdom (warnings, defaults)
    support_from: tuple[str, ...] | None = None  # predictors that count for calls; None = all scored
    prefilter_by: str = "rnahybrid"         # the cheap predictor that decides what IntaRNA sees
    prefilter_keep_classes: tuple[str, ...] = ()  # classes IntaRNA always scores in full
    report_dir: str = ""                    # where predictor_report.tsv / prefilter_summary.tsv go


def scored_classes(classes, cfg: ScoreConfig):
    return tuple(c for c in classes if cfg.score_weak_sites or c in STRONG or c.startswith("seed"))


def check_tools(cfg: ScoreConfig) -> None:
    """Fail early, before any long step, when a needed program is not on PATH."""
    need = []
    if "intarna" in cfg.predictors:
        need.append(cfg.intarna_bin)
    if "rnahybrid" in cfg.predictors:
        need.append(cfg.rnahybrid_bin)
        if cfg.calibrate_samples > 0:
            need.append(cfg.rnacalibrate_bin)
    missing = [b for b in need if shutil.which(b) is None]
    if missing:
        raise SystemExit(f"score-sites: not on PATH: {', '.join(missing)} "
                         "(activate the environment with IntaRNA and RNAhybrid)")


def summarise_errors(scores, log, max_share: float = 0.05) -> None:
    errors = [(job, res["_error"]) for job, res in scores.items() if "_error" in res]
    if not errors:
        return
    log(f"WARNING: {len(errors)} of {len(scores)} sites failed and are left unscored; "
        f"first: {errors[0][1]}")
    if len(errors) > max_share * len(scores):
        raise RuntimeError(f"{len(errors)}/{len(scores)} sites failed (> {max_share:.0%}); "
                           f"first error: {errors[0][1]}")
    for job, _ in errors:
        scores[job].pop("_error")  # keep what the other predictor produced for this site


def run_signature(stage: str, guides, unique, seed, classes, cfg: ScoreConfig) -> str:
    """Identifies what a checkpoint holds: the small RNAs and decoys, the target sequences and
    every setting that changes a score. A checkpoint with another signature is ignored."""
    h = hashlib.sha1()
    h.update(repr((stage, seed, classes, sorted(guides.items()))).encode())
    h.update(repr((len(unique), sum(map(len, unique)), unique[:3], unique[-3:])).encode())
    settings = [cfg.rnahybrid_flank, cfg.rnahybrid_require_site, cfg.calibrate_samples > 0]
    if stage == "final":
        settings += [cfg.predictors, cfg.intarna_flank, cfg.intarna_prefilter, cfg.temperature,
                     cfg.intarna_extra]
        if cfg.intarna_prefilter < 1.0 and cfg.prefilter_keep_classes:  # changes which sites IntaRNA sees
            settings.append(sorted(cfg.prefilter_keep_classes))
    h.update(repr(settings).encode())
    return h.hexdigest()


def load_checkpoint(path: str, signature: str):
    """Scores saved by an earlier run with the same signature, or None. The files are pickles
    written by this program into the work directory; do not point it at files from elsewhere."""
    try:
        with open(path, "rb") as fh:
            data = pickle.load(fh)
    except (OSError, EOFError, pickle.UnpicklingError, AttributeError, ImportError):
        return None
    return data["scores"] if isinstance(data, dict) and data.get("signature") == signature else None


def save_checkpoint(path: str, signature: str, scores) -> None:
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        pickle.dump({"signature": signature, "scores": scores}, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)


def score_one(guide: str, target: str, site, xi_theta, cfg: ScoreConfig) -> dict:
    out: dict = {}
    if "intarna" in cfg.predictors:
        out.update(intarna.score_site(guide, target, site, cfg.intarna_flank, cfg.temperature,
                                      cfg.intarna_bin, cfg.intarna_extra, cfg.timeout) or {})
    if "rnahybrid" in cfg.predictors:
        out.update(rnahybrid.score_site(guide, target, site, xi_theta, cfg.rnahybrid_flank,
                                        cfg.rnahybrid_bin, cfg.rnahybrid_require_site,
                                        cfg.timeout) or {})
    return out


_PER_SITE = score_one  # a replaced score_one (tests, experiments) keeps the per-site path

# state shared with worker processes (inherited on fork)
_STATE: dict = {}


def _init(state):
    _STATE.update(state)


def _error(exc) -> dict:
    return {"_error": f"{type(exc).__name__}: {exc}"[:200]}


def _run_job(job):
    """Result of one site (or ``{"_error": ...}``): one process per predictor and site."""
    gid, si, start, end, t1, cls = job
    site = Site(start, end, t1, cls)
    s = _STATE
    try:
        return score_one(s["guides"][gid], s["unique"][si], site, s["xi"].get(gid), s["cfg"])
    except Exception as exc:  # one bad window must not kill a run of hours
        return _error(exc)


def _work(job):
    return job, _run_job(job)


def plan_intarna_calls(jobs, idxs, state, cfg: ScoreConfig, size: int) -> list[list[int]]:
    """Groups of job indices that one IntaRNA process can score (``intarna.score_pairs``):
    - jobs on the SAME WINDOW (same region, window bounds and layout: the same t1, guide
      length and seed range) from different small RNAs: one target, several queries, so the
      target accessibility, which is most of the cost, is computed once;
    - the rest, grouped by small RNA and layout (``--tRegion``/``--seedTRange`` are the same
      for every target of a call) in calls of at most ``size`` windows. Windows of another
      layout (sites near a region end) are never mixed in and never padded."""
    unique, guides, flank = state["unique"], state["guides"], cfg.intarna_flank
    by_window: dict[tuple, list[int]] = {}
    for i in idxs:
        gid, si, start, end, t1, cls = jobs[i]
        lo, hi, treg, seedt = intarna.window_coords(len(unique[si]), Site(start, end, t1, cls),
                                                    len(guides[gid]), flank)
        by_window.setdefault((si, lo, hi, treg, seedt, intarna.seed_range(cls)), []).append(i)
    calls: list[list[int]] = []
    singles: dict[tuple, list[int]] = {}
    for key, ii in by_window.items():
        if len(ii) > 1:
            calls += [ii[k:k + 100] for k in range(0, len(ii), 100)]
        else:
            singles.setdefault((jobs[ii[0]][0],) + key[3:], []).append(ii[0])
    for ii in singles.values():
        calls += [ii[k:k + size] for k in range(0, len(ii), size)]
    return calls


def _intarna_call(jobs) -> list:
    """Results (dict, or an Exception for a site that could not be scored) of jobs that share
    a window or a small RNA and a layout, in ONE IntaRNA process."""
    s = _STATE
    cfg, guides, unique = s["cfg"], s["guides"], s["unique"]
    geo = [intarna.window_coords(len(unique[j[1]]), Site(j[2], j[3], j[4], j[5]),
                                 len(guides[j[0]]), cfg.intarna_flank) for j in jobs]
    layouts = {(g[2], g[3], intarna.seed_range(j[5])) for g, j in zip(geo, jobs)}
    if len(layouts) != 1:  # cannot happen when the calls come from plan_intarna_calls
        out = []
        for j in jobs:
            try:
                out.append(intarna.score_site(guides[j[0]], unique[j[1]], Site(*j[2:5], j[5]),
                                              cfg.intarna_flank, cfg.temperature, cfg.intarna_bin,
                                              cfg.intarna_extra, cfg.timeout) or {})
            except Exception as exc:  # noqa: BLE001
                out.append(exc)
        return out
    windex: dict[tuple, int] = {}
    qindex: dict[str, int] = {}
    for j, g in zip(jobs, geo):
        windex.setdefault((j[1], g[0], g[1]), len(windex))
        qindex.setdefault(j[0], len(qindex))
    windows = [unique[si][lo:hi] for si, lo, hi in windex]
    pairs = intarna.score_pairs([guides[g] for g in qindex], windows, jobs[0][5], geo[0][2],
                                geo[0][3], cfg.temperature, cfg.intarna_bin, cfg.intarna_extra,
                                cfg.timeout, s.get("tmpdir"))
    out = []
    for j, g in zip(jobs, geo):
        r = pairs[(windex[(j[1], g[0], g[1])], qindex[j[0]])]
        if r is None:
            r = {}
        elif not isinstance(r, Exception):  # region coordinates (1-based), as score_site
            r["intarna_t_start"] += g[0]
            r["intarna_t_end"] += g[0]
        out.append(r)
    return out


def _work_intarna(task):
    """A task is a list of calls (lists of jobs); the result is a list of lists of result
    dicts, aligned."""
    return [[_error(r) if isinstance(r, Exception) else r for r in _intarna_call(call)]
            for call in task]


def _work_batch(jobs):
    """Results, aligned with ``jobs``, of a chunk of sites that share the small RNA and the
    seed range. IntaRNA is still one process per site; RNAhybrid scores the whole chunk in one
    process (``rnahybrid.score_sites``). A site that fails in either predictor gets
    ``{"_error": ...}`` and nothing else, as in the per-site path."""
    s = _STATE
    cfg = s["cfg"]
    gid = jobs[0][0]
    guide = s["guides"][gid]
    sites = [Site(j[2], j[3], j[4], j[5]) for j in jobs]
    targets = [s["unique"][j[1]] for j in jobs]
    out: list = [{} for _ in jobs]
    failed: list = [None] * len(jobs)
    if "intarna" in cfg.predictors and cfg.intarna_batch > 0:
        for call in plan_intarna_calls(jobs, range(len(jobs)), s, cfg, cfg.intarna_batch):
            for i, r in zip(call, _intarna_call([jobs[i] for i in call])):
                if isinstance(r, Exception):
                    failed[i] = r
                else:
                    out[i].update(r)
    elif "intarna" in cfg.predictors:
        for i, (target, site) in enumerate(zip(targets, sites)):
            try:
                out[i].update(intarna.score_site(guide, target, site, cfg.intarna_flank,
                                                 cfg.temperature, cfg.intarna_bin,
                                                 cfg.intarna_extra, cfg.timeout) or {})
            except Exception as exc:  # noqa: BLE001
                failed[i] = exc
    if "rnahybrid" in cfg.predictors:
        live = [i for i in range(len(jobs)) if failed[i] is None]
        try:
            res = rnahybrid.score_sites(
                guide, [(targets[i], sites[i]) for i in live], s["xi"].get(gid),
                cfg.rnahybrid_flank, cfg.rnahybrid_bin, cfg.rnahybrid_require_site,
                cfg.timeout, tmpdir=s.get("tmpdir"))
        except Exception as exc:  # noqa: BLE001  (e.g. no space for the temporary FASTA)
            res = [exc] * len(live)
        for i, r in zip(live, res):
            if isinstance(r, Exception):
                failed[i] = r
            else:
                out[i].update(r or {})
    return [_error(failed[i]) if failed[i] is not None else out[i] for i in range(len(jobs))]


def _calibrate_work(gid):
    s = _STATE
    # RNAcalibrate on random sequences with the dinucleotide frequencies of the target set, at
    # its default lengths (500 +- 300 nt) and with the same helix constraint. Calibrating at the
    # real window length (~55 nt) was tried and fails: below ~100 nt the extreme-value fit gives
    # a negative xi (the method was designed for long targets); RNAhybrid normalises the MFE by
    # ln(target length * query length) and so applies these parameters to the short window.
    guide = s["guides"][gid]
    return gid, rnahybrid.calibrate(guide, s["targets_fasta"], s["seed"], CALIBRATION_LENGTH[0],
                                    CALIBRATION_LENGTH[1], s["cfg"].calibrate_samples,
                                    s["cfg"].rnacalibrate_bin)


def prefilter_cutoffs(scores, share: float) -> dict[str, float]:
    """Per site class, the RNAhybrid MFE at the ``share`` quantile (lower is better) over
    all scored sites, real and decoy."""
    by_class: dict[str, list[float]] = defaultdict(list)
    for job, res in scores.items():
        if res.get("rnahybrid_mfe") is not None:
            by_class[job[5]].append(res["rnahybrid_mfe"])
    out = {}
    for cls, vals in by_class.items():
        vals.sort()
        out[cls] = vals[max(0, min(len(vals) - 1, int(share * len(vals)) - 1))]
    return out


def parallel_map(fn, items, cfg, state, chunksize=None, log=None, label="", sizes=None):
    """``map`` over worker processes, results in the order of ``items``; logs progress and an
    ETA about every 5% when ``log`` is given and there is a lot to do. ``sizes`` (the number of
    sites in each item, for items that are chunks of sites) makes the progress count sites.
    ``chunksize`` is how many items travel to a worker together (default: 1/20 of a worker's
    share, at most 64: big enough to cut the overhead, small enough to balance the load)."""
    items = items if isinstance(items, list) else list(items)
    n = sum(sizes) if sizes is not None else len(items)
    step = max(1, n // 20) if log and n >= 2000 else 0
    t0 = time.time()
    out = []
    progress = {"done": 0, "next": step}

    def tick(k):
        progress["done"] += k
        done = progress["done"]
        if step and done >= progress["next"] and done < n:
            rate = done / (time.time() - t0)
            log(f"{label}{done}/{n} ({done / n:.0%}), ETA {(n - done) / rate / 60:.0f} min")
            progress["next"] = done - done % step + step

    if cfg.workers <= 1:
        _init(state)
        for i, item in enumerate(items):
            out.append(fn(item))
            tick(sizes[i] if sizes is not None else 1)
        return out
    if chunksize is None:
        chunksize = max(1, min(64, len(items) // (cfg.workers * 20)))
    with ProcessPoolExecutor(cfg.workers, initializer=_init, initargs=(state,)) as ex:
        for i, res in enumerate(ex.map(fn, items, chunksize=chunksize)):
            out.append(res)
            tick(sizes[i] if sizes is not None else 1)
    return out


def _plan_work(gid):
    """Sites of one small RNA in every sequence, flat: (sequence index, start, end, t1, class
    index) per site, as an array (a compact message instead of a list of tuples)."""
    s = _STATE
    guide, classes = s["guides"][gid], s["classes"]
    cls_index = {c: k for k, c in enumerate(classes)}
    flat = array("l")
    for si, seq in enumerate(s["unique"]):
        for site in scan(guide, seq, s["seed"], classes):
            flat.extend((si, site.start, site.end, site.t1, cls_index[site.site_class]))
    return flat


def plan_jobs(guides: dict[str, str], unique: list[str], seed, classes, cfg=None):
    """All scorable sites: (guide id, sequence index, start, end, t1, class). The scan runs in
    ``cfg.workers`` processes (one task per small RNA; with 800+ small RNAs and 40,000
    sequences it is minutes in a single process). The tuples share their integer objects
    (about 8 M sites are held in memory, and again as the keys of the scores)."""
    classes = tuple(classes)
    state = {"guides": guides, "unique": unique, "seed": seed, "classes": classes}
    workers = min(cfg.workers, len(guides)) if cfg else 1
    pool_cfg = replace(cfg, workers=workers) if cfg else ScoreConfig(workers=1)
    gids = list(guides)
    ints = list(range(max(max(map(len, unique), default=0), len(unique)) + 2))
    jobs = []
    for gid, flat in zip(gids, parallel_map(_plan_work, gids, pool_cfg, state, chunksize=1)):
        for k in range(0, len(flat), 5):
            si, start, end, t1, c = flat[k:k + 5]
            jobs.append((gid, ints[si], ints[start], ints[end], ints[t1], classes[c]))
    return jobs


def _chunks(jobs, size: int):
    """Indices into ``jobs`` in chunks of at most ``size`` sites that share the small RNA and
    the seed range (the helix constraint), in order of first appearance."""
    groups: dict[tuple, array] = {}
    seed_of: dict[str, tuple] = {}
    for i, job in enumerate(jobs):
        cls = job[5]
        if cls not in seed_of:
            seed_of[cls] = intarna.seed_range(cls)
        groups.setdefault((job[0], seed_of[cls]), array("l")).append(i)
    return [idx[k:k + size] for idx in groups.values() for k in range(0, len(idx), size)]


def score_jobs(jobs, cfg: ScoreConfig, state, log=None, label="") -> list:
    """Result dict of every job, in the order of ``jobs``. Batch paths (same numbers as one
    process per site, much cheaper):
    - RNAhybrid: chunks of ``cfg.rnahybrid_batch`` sites of one small RNA and seed range, one
      process per chunk;
    - IntaRNA alone (the second stage of the two-stage mode): calls planned by
      ``plan_intarna_calls`` (shared windows, windows of one layout), ``cfg.intarna_batch``
      windows per process;
    - both in one pass: RNAhybrid chunks of at most 50 sites, IntaRNA batched inside the chunk.
    A replaced ``score_one`` (tests) or a batch size of 0 keeps the per-site path."""
    ia, rh = "intarna" in cfg.predictors, "rnahybrid" in cfg.predictors
    if score_one is not _PER_SITE or (ia and cfg.intarna_batch <= 0) or (
            rh and cfg.rnahybrid_batch <= 0):
        return parallel_map(_run_job, jobs, cfg, state, log=log, label=label)
    results: list = [None] * len(jobs)
    if not rh:
        calls = plan_intarna_calls(jobs, range(len(jobs)), state, cfg, cfg.intarna_batch)
        tasks, cur, n = [], [], 0
        for call in calls:   # ~40 sites per message
            cur.append(call)
            n += len(call)
            if n >= 40:
                tasks.append(cur)
                cur, n = [], 0
        if cur:
            tasks.append(cur)
        done = parallel_map(_work_intarna, [[[jobs[i] for i in call] for call in t] for t in tasks],
                            cfg, state, chunksize=1, log=log, label=label,
                            sizes=[sum(map(len, t)) for t in tasks])
        for task, res in zip(tasks, done):
            for call, rr in zip(task, res):
                for i, r in zip(call, rr):
                    results[i] = r
        return results
    # a chunk with IntaRNA inside is dominated by its IntaRNA calls: keep it small so the
    # progress and the end of the run stay balanced
    size = min(cfg.rnahybrid_batch, 50) if "intarna" in cfg.predictors else cfg.rnahybrid_batch
    chunks = _chunks(jobs, size)
    tasks = [[jobs[i] for i in idx] for idx in chunks]
    done = parallel_map(_work_batch, tasks, cfg, state, chunksize=1, log=log, label=label,
                        sizes=[len(idx) for idx in chunks])
    for idx, res in zip(chunks, done):
        for i, r in zip(idx, res):
            results[i] = r
    return results


def dry_run_summary(jobs, cfg: ScoreConfig) -> str:
    per_class = Counter(j[5] for j in jobs)
    n = len(jobs) * len(cfg.predictors)
    lines = [f"{len(jobs)} sites to score ({n} predictor runs):"]
    lines += [f"  {c}: {k}" for c, k in sorted(per_class.items())]
    # measured on the shared server: RNAhybrid 0.01 s per process, ~0.0005 s per site in a batch;
    # IntaRNA ~0.1 s per window (target accessibility), ~0.09 s when batched, much less when
    # small RNAs share windows (keep_seed decoys); only the best share goes to IntaRNA in
    # two-stage mode
    per_rh = 0.0005 if cfg.rnahybrid_batch > 0 else 0.01
    per_ia = 0.09 if cfg.intarna_batch > 0 else 0.1
    two_stage = cfg.intarna_prefilter < 1.0 and len(cfg.predictors) == 2
    est = 0.0
    if "rnahybrid" in cfg.predictors:
        est += len(jobs) * per_rh
    if "intarna" in cfg.predictors:
        est += len(jobs) * (cfg.intarna_prefilter if two_stage else 1.0) * per_ia
    est /= max(1, cfg.workers)
    lines.append(f"rough time on {cfg.workers} worker(s): {est / 3600:.1f} h (RNAhybrid "
                 f"{per_rh} s/site, IntaRNA {per_ia} s/site; less if small RNAs share windows)")
    return "\n".join(lines)


def curve_rows(predictor, cls, real, decoy, weights):
    """FDR along the score for the report: one row per distinct real score.
    ``weights`` are the decoy weights (1 / decoys of the candidate), so the column is
    the expected number of false sites per candidate."""
    from bisect import bisect_right
    from itertools import accumulate
    rows = []
    r = sorted(real)
    pairs = sorted(zip(decoy, weights))
    d = [x for x, _ in pairs]
    cum = list(accumulate(w for _, w in pairs))
    for x in sorted(set(r)):
        pr, k = bisect_right(r, x), bisect_right(d, x)
        false_sites = cum[k - 1] if k else 0.0
        rows.append({"predictor": predictor, "site_class": cls, "threshold": x,
                     "real_sites": pr, "decoy_sites_per_candidate": round(false_sites, 4),
                     "fdr": round(min(1.0, false_sites / pr), 4)})
    return rows


def run_scoring(candidates, decoys, regions: list[Region], seed=(2, 8), classes=None,
                cfg: ScoreConfig | None = None, workdir: str = ".", dry_run: bool = False,
                log=print):
    """Score all sites. ``candidates`` are (id, sequence); ``decoys`` (id, candidate, seq)."""
    cfg = cfg or ScoreConfig()
    pred.check_config(cfg, log)
    classes = scored_classes(tuple(classes or site_classes(seed)), cfg)
    by_seq: dict[str, list[Region]] = defaultdict(list)
    for r in regions:
        by_seq[r.sequence].append(r)
    unique = list(by_seq)

    guides = {cid: seq for cid, seq in candidates}
    decoy_of = {}
    for did, cid, seq in decoys:
        if cid in guides:
            guides[did] = seq
            decoy_of[did] = cid
    n_decoys = Counter(decoy_of.values())
    if not n_decoys:
        raise ValueError("score-sites needs decoys (decoys.fasta has none for these candidates): "
                         "without them there is no false discovery rate")
    without = [cid for cid, _ in candidates if n_decoys[cid] == 0]
    if without:
        log(f"WARNING: {len(without)} candidate(s) have no decoys and get no q-values "
            f"(e.g. {without[0]}); rerun the decoys step after editing the candidate list")

    jobs = plan_jobs(guides, unique, seed, classes, cfg)
    log(dry_run_summary(jobs, cfg))
    if dry_run:
        return None
    check_tools(cfg)
    os.makedirs(workdir, exist_ok=True)

    state = {"guides": guides, "unique": unique, "xi": {}, "cfg": cfg, "seed": seed,
             "tmpdir": workdir}
    if "rnahybrid" in cfg.predictors and cfg.calibrate_samples > 0:
        os.makedirs(workdir, exist_ok=True)
        lens = [len(s) for s in unique] or [100]
        state["targets_fasta"] = os.path.join(workdir, "calibration_targets.fa")
        with open(state["targets_fasta"], "w") as fh:
            for i, s in enumerate(unique):
                fh.write(f">r{i}\n{intarna.rna(s)}\n")
        state["mean"] = int(statistics.mean(lens))
        state["sd"] = max(1, int(statistics.pstdev(lens)))
        # only the real small RNAs: decoys are compared by MFE / energy, not by p-value,
        # so they use the preset and skip the (costly) calibration
        needed = sorted({j[0] for j in jobs if j[0] not in decoy_of})
        t0 = time.time()
        for gid, xt in parallel_map(_calibrate_work, needed, cfg, state, 1):
            if xt:
                state["xi"][gid] = xt
        log(f"RNAcalibrate: {len(state['xi'])}/{len(needed)} real small RNAs in {time.time() - t0:.0f} s")

    os.makedirs(workdir, exist_ok=True)
    ck_final = os.path.join(workdir, "checkpoint_scores.pkl")
    ck_stage1 = os.path.join(workdir, "checkpoint_stage1.pkl")
    sig_final = run_signature("final", guides, unique, seed, classes, cfg)
    sig_stage1 = run_signature("stage1", guides, unique, seed, classes, cfg)
    scores = load_checkpoint(ck_final, sig_final) if cfg.resume else None
    if scores is not None:
        log(f"resumed all {len(scores)} site scores from {ck_final}")
    else:
        t0 = time.time()
        two_stage = (cfg.intarna_prefilter < 1.0 and "intarna" in cfg.predictors
                     and "rnahybrid" in cfg.predictors)
        if two_stage:
            s1 = replace(cfg, predictors=("rnahybrid",))
            scores = load_checkpoint(ck_stage1, sig_stage1) if cfg.resume else None
            if scores is not None:
                log(f"resumed the RNAhybrid scores of {len(scores)} sites from {ck_stage1}")
            else:
                scores = dict(zip(jobs, score_jobs(jobs, s1, {**state, "cfg": s1}, log,
                                                   "RNAhybrid ")))
                save_checkpoint(ck_stage1, sig_stage1, scores)
            cutoff = prefilter_cutoffs(scores, cfg.intarna_prefilter)
            keep = pred.prefilter_keep(jobs, scores, cutoff, cfg.prefilter_keep_classes)
            log(f"RNAhybrid scored {len(jobs)} sites in {time.time() - t0:.0f} s; "
                f"{len(keep)} pass the prefilter (best {cfg.intarna_prefilter:.0%} per class, "
                f"cutoffs {', '.join(f'{c} {v:.1f}' for c, v in sorted(cutoff.items()))})")
            s2 = replace(cfg, predictors=("intarna",))
            for job, res in zip(keep, score_jobs(keep, s2, {**state, "cfg": s2}, log,
                                                 "IntaRNA ")):
                scores[job].update(res)
        else:
            scores = dict(zip(jobs, score_jobs(jobs, cfg, state, log, "sites ")))
        log(f"scored {len(jobs)} sites in {time.time() - t0:.0f} s")
        summarise_errors(scores, log, cfg.max_error_share)
        save_checkpoint(ck_final, sig_final, scores)

    # q-values per predictor and class, pooled and per candidate
    qpool: dict[tuple, float] = {}
    qcand: dict[tuple, float] = {}
    curves = []
    for pname, col in PREDICTORS.items():
        if pname not in cfg.predictors:
            continue
        unscored = sum(1 for res in scores.values() if res.get(col) is None)
        log(f"{pname}: {unscored} of {len(scores)} sites without a score "
            f"({unscored / max(1, len(scores)):.1%}; no interaction, filtered or failed)")
        for cls in classes:
            real, real_keys, decoy, weights = [], [], [], []
            decoy_by_cand: dict[str, list[float]] = defaultdict(list)
            real_by_cand: dict[str, list[tuple]] = defaultdict(list)
            for job, res in scores.items():
                if job[5] != cls or res.get(col) is None:
                    continue
                parent = decoy_of.get(job[0])
                if parent is not None:
                    decoy.append(res[col])
                    weights.append(1.0 / n_decoys[parent])
                    decoy_by_cand[parent].append(res[col])
                elif n_decoys[job[0]] > 0:   # without decoys a site cannot get a q-value
                    real.append(res[col])
                    real_keys.append(job)
                    real_by_cand[job[0]].append(job)
            if not real:
                continue
            plus = (len(n_decoys) / sum(n_decoys.values())) if cfg.fdr_pseudocount else 0.0
            for key, q in zip(real_keys, qvalues(real, decoy, 0, decoy_weights=weights,
                                                 plus_one=plus)):
                qpool[(pname, key)] = q
            curves += curve_rows(pname, cls, real, decoy, weights)
            for cid, keys in real_by_cand.items():
                vals = [scores[k][col] for k in keys]
                plus_c = 1.0 / n_decoys[cid] if cfg.fdr_pseudocount else 0.0
                for k, q in zip(keys, qvalues(vals, decoy_by_cand.get(cid, []), n_decoys[cid],
                                              plus_one=plus_c)):
                    qcand[(pname, k)] = q

    counted = pred.counted_predictors(cfg.predictors, cfg.support_from)
    report_predictors(scores, decoy_of, n_decoys, classes, qpool, counted, cfg, log)

    site_rows = []
    for job, res in scores.items():
        gid, si, start, end, t1, cls = job
        if gid in decoy_of:
            continue
        site = Site(start, end, t1, cls)
        seq = unique[si]
        support, called_by, best_q = pred.calls_of(
            {p: qpool.get((p, job)) for p in counted}, counted, cfg.fdr)
        for r in by_seq[seq]:
            g, g_end = r.genomic(start), r.genomic(end - 1)
            genomic = (f"{r.contig}:{min(g, g_end) + 1}-{max(g, g_end) + 1}({r.strand})"
                       if g is not None and g_end is not None else "")
            row = {"candidate": gid, "region_id": r.id, "transcript": r.transcript,
                   "gene": r.gene, "region": r.region, "site_class": cls,
                   "site_start": start + 1, "site_end": end, "genomic": genomic,
                   **context(guides[gid], seq, site, r.region), **res, "support": support,
                   "called_by": called_by, "best_q": best_q}
            for p in cfg.predictors:
                row[f"{p}_q"] = qpool.get((p, job), "")
                row[f"{p}_q_cand"] = qcand.get((p, job), "")
            site_rows.append(row)
    order = {c: i for i, c in enumerate(classes)}
    site_rows.sort(key=lambda r: (r["candidate"], -r["support"], order[r["site_class"]],
                                  site_energy(r, cfg)))

    target_rows = summarise_targets(site_rows, cfg)
    return classes, site_rows, target_rows, curves


def site_energy(row, cfg: ScoreConfig) -> float:
    """Sort key within a class: the IntaRNA energy, or the RNAhybrid MFE when IntaRNA is not
    computed (never mixed: the two scales differ); 0 for a site without a score."""
    col = "intarna_E" if "intarna" in cfg.predictors else "rnahybrid_mfe"
    return row[col] if row.get(col) is not None else 0.0


def report_predictors(scores, decoy_of, n_decoys, classes, qpool, counted, cfg: ScoreConfig, log):
    """Suitability diagnostic and prefilter cost, logged and (with ``cfg.report_dir``) written to
    ``predictor_report.tsv`` and ``prefilter_summary.tsv``. Advice only: calls are not touched."""
    rows = pred.diagnose(scores, decoy_of, n_decoys, cfg.predictors, classes, qpool, cfg.fdr,
                         counted)
    if rows:
        pred.log_diagnosis(rows, cfg.fdr, log)
    if pred.two_stage(cfg):
        pf = pred.prefilter_summary(scores, decoy_of, prefilter_cutoffs(scores, cfg.intarna_prefilter),
                                    cfg.prefilter_keep_classes, classes)
        pred.log_prefilter(pf, cfg, log)
    else:
        pf = pred.prefilter_off_rows(scores, decoy_of, cfg)
    if cfg.report_dir:
        os.makedirs(cfg.report_dir, exist_ok=True)
        pred.write_report(rows, os.path.join(cfg.report_dir, "predictor_report.tsv"))
        pred.write_prefilter(pf, os.path.join(cfg.report_dir, "prefilter_summary.tsv"))


def summarise_targets(site_rows, cfg: ScoreConfig):
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in site_rows:
        groups[(r["candidate"], r["transcript"], r["gene"])].append(r)
    out = []
    for (cid, transcript, gene), rows in groups.items():
        def best(col, fn=min):
            vals = [r[col] for r in rows if r.get(col) not in (None, "")]
            return fn(vals) if vals else ""
        # highest support, then the strongest class, then the lowest IntaRNA energy
        top = min(rows, key=lambda r: (-r["support"], CLASS_RANK.get(r["site_class"], 99),
                                       r["intarna_E"] if r.get("intarna_E") is not None else 0.0))
        out.append({
            "candidate": cid, "transcript": transcript, "gene": gene,
            "scored_sites": len(rows), "support": top["support"],
            # union over the sites of the gene: one counted predictor at one site is a call
            "called_by": ",".join(p for p in cfg.predictors
                                  if any(p in r.get("called_by", "").split(",") for r in rows)),
            "best_q": best("best_q"),
            "best_site_class": top["site_class"],
            "best_intarna_E": best("intarna_E"), "best_rnahybrid_mfe": best("rnahybrid_mfe"),
            "best_intarna_q": best("intarna_q"), "best_rnahybrid_q": best("rnahybrid_q"),
            "sites_q_pass": sum(1 for r in rows if r["support"] > 0),
        })
    out.sort(key=lambda r: (r["candidate"], -r["support"], r["best_q"] == "",
                            r["best_q"] if r["best_q"] != "" else 1.0))
    return out


SCORE_COLUMNS = ["candidate", "region_id", "transcript", "gene", "region", "site_class",
                 "site_start", "site_end", "genomic", "dist_5p", "dist_3p", "au_flanks",
                 "supp_pairing", "near_stop",
                 "intarna_E", "intarna_ED_target", "intarna_ED_query", "intarna_Pu_target",
                 "intarna_Pu_query", "intarna_E_hybrid", "intarna_t_start", "intarna_t_end",
                 "intarna_q_start", "intarna_q_end", "intarna_duplex", "intarna_q",
                 "intarna_q_cand",
                 "rnahybrid_mfe", "rnahybrid_p", "rnahybrid_pos", "rnahybrid_at_site",
                 "rnahybrid_q", "rnahybrid_q_cand", "support", "called_by", "best_q"]
TARGET_SCORE_COLUMNS = ["candidate", "transcript", "gene", "scored_sites", "support",
                        "called_by", "best_q",
                        "best_site_class", "best_intarna_E", "best_intarna_q",
                        "best_rnahybrid_mfe", "best_rnahybrid_q", "sites_q_pass"]
CURVE_COLUMNS = ["predictor", "site_class", "threshold", "real_sites",
                 "decoy_sites_per_candidate", "fdr"]


def write_scores(site_rows, target_rows, curves, sites_path, targets_path, curves_path):
    write_tsv(sites_path, site_rows, SCORE_COLUMNS)
    write_tsv(targets_path, target_rows, TARGET_SCORE_COLUMNS)
    write_tsv(curves_path, curves, CURVE_COLUMNS)
