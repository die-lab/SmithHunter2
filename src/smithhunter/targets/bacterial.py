"""Module B, bacterial mode: small RNA -> mRNA start-codon windows with IntaRNA, q-values from
decoys.

Model (bacterial trans-encoded sRNAs, usually Hfq-dependent; no AGO, no guide strand): an sRNA
pairs with an mRNA in a stretch of 6-12+ nt that can lie anywhere in the translation-initiation
region, with G:U pairs; the targets of trusted sRNA-mRNA pairs fall in -150..+100 around the
start codon. So there is no seed *position* and no site scan: each candidate is scored against
the whole window of every CDS (``regions.bacterial_windows``) with IntaRNA, which finds the best
interaction in the window under its own seed rule (``seed_bp`` base pairs, anywhere, or none) and
local accessibility of both partners (sliding window ``acc_window``/``acc_length``). RNAhybrid
and the site classes are not used (RNAhybrid alone is far worse than IntaRNA for bacterial sRNAs,
so it is not a prefilter either).

The unit of the false discovery rate is **one best interaction per (candidate, window)**:
the E of the best IntaRNA interaction of a real small RNA with a window is compared with the same
number for its decoys (dinucleotide shuffles, ``shuffle.make_decoys``, scored identically). The
q-value of a real window is the lowest FDR over the thresholds it passes,

    FDR(E) = (decoy windows with best E <= E, each weighted 1 / decoys of its candidate
              + one pseudo decoy) / (real windows with best E <= E)

pooled over candidates (``intarna_q``) and, alongside, for each candidate against its own decoys
(``intarna_q_cand``). One interaction per window, not one per seed match, keeps windows with many
possible seeds from counting many times, and the null (a random sequence of the same
composition) is tested against the same number of opportunities. There is no ``support``: with a
single predictor ``intarna_q <= fdr`` is the call (``pass_fdr``).

Identical window sequences (repeated genes, IS elements) are scored once and reported for every
region that carries them. The run is resumable: results are saved in the work directory during the
run (pickles written by this program; do not point it at files from elsewhere) and reused when
the inputs and settings are identical. Results of decoys are reduced to their energy.

Coordinates: ``intarna_t_start``/``intarna_t_end`` are 1-based positions in the window;
``rel_start``/``rel_end`` are relative to the start codon (A of the start codon = +1, the base
before it = -1, no 0); ``hit_region`` is ``5UTR`` for an interaction entirely upstream of the start
codon (5'UTR/RBS, or the end of the preceding gene, see ``upstream_cds``), ``CDS`` for one entirely
from +1 on and ``spans_start`` for one that crosses it; ``hit_rbs`` is 1 when it overlaps the
ribosome-binding zone (default -20..-1).
"""

from __future__ import annotations

import hashlib
import os
import pickle
import shutil
import time
from bisect import bisect_right
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from itertools import accumulate

from ..seqio import write_tsv
from . import intarna
from .fdr import qvalues
from .regions import Region


@dataclass
class BacterialConfig:
    temperature: float = 37.0
    seed_bp: int = 7             # base pairs of the IntaRNA seed (any position); 6-7 is typical
    no_seed: bool = False        # IntaRNA --noSeed: no seed requirement at all
    acc_window: int = 150        # IntaRNA --tAccW (local accessibility, sliding window)
    acc_length: int = 100        # IntaRNA --tAccL (maximal base-pair span)
    rbs_zone: tuple[int, int] = (-20, -1)   # relative to the start codon, for ``hit_rbs``
    workers: int = 1
    fdr: float = 0.05
    max_error_share: float = 0.05   # stop the run when more windows than this fail
    timeout: float = 300            # seconds per IntaRNA call; a hang is an error of that window
    resume: bool = True
    checkpoint_seconds: float = 60  # save partial results at most this often
    fdr_pseudocount: bool = True    # (D + 1) / T
    intarna_bin: str = "IntaRNA"
    intarna_extra: list[str] = field(default_factory=list)


# ---- coordinates -----------------------------------------------------------------------

def to_rel(pos1: int, anchor: int) -> int:
    """Position relative to the start codon of the 1-based window position ``pos1``; ``anchor``
    is the 0-based index of the first base of the start codon. There is no position 0."""
    return pos1 - anchor if pos1 > anchor else pos1 - 1 - anchor


def hit_region(rel_start: int, rel_end: int) -> str:
    if rel_end < 0:
        return "5UTR"
    if rel_start > 0:
        return "CDS"
    return "spans_start"


def derive(res: dict, anchor: int, cfg: BacterialConfig) -> dict:
    """Interaction coordinates relative to the start codon and the part of the window hit."""
    if res.get("intarna_t_start") is None:
        return {}
    rs, re_ = to_rel(res["intarna_t_start"], anchor), to_rel(res["intarna_t_end"], anchor)
    z0, z1 = cfg.rbs_zone
    return {"rel_start": rs, "rel_end": re_, "hit_region": hit_region(rs, re_),
            "hit_rbs": int(rs <= z1 and re_ >= z0)}


# ---- IntaRNA on one window -------------------------------------------------------------

def build_command(query: str, window: str, cfg: BacterialConfig) -> list[str]:
    """IntaRNA over the whole window: best interaction, local accessibility, a seed of
    ``seed_bp`` base pairs anywhere (or none)."""
    acc_w = min(cfg.acc_window, len(window)) if cfg.acc_window else 0
    acc_l = min(cfg.acc_length, acc_w) if cfg.acc_length and acc_w else cfg.acc_length
    cmd = [cfg.intarna_bin, "-q", intarna.rna(query), "-t", intarna.rna(window),
           "--outMode", "C", "--outCsvCols", ",".join(intarna.CSV_COLUMNS),
           "--temperature", str(cfg.temperature), "--threads", "1",
           "--tAccW", str(acc_w), "--tAccL", str(acc_l)]
    cmd += ["--noSeed"] if cfg.no_seed else ["--seedBP", str(cfg.seed_bp)]
    return cmd + list(cfg.intarna_extra)


def score_window(guide: str, window: str, cfg: BacterialConfig) -> dict | None:
    """Best IntaRNA interaction of ``guide`` with ``window`` (None when there is none)."""
    proc = intarna.run_tool(build_command(guide, window, cfg), cfg.timeout)
    if proc.returncode != 0:
        raise RuntimeError(f"IntaRNA failed ({proc.returncode}): {proc.stderr.strip()[:300]}")
    return intarna.parse_csv(proc.stdout)


# ---- parallel map with checkpoints (own copy of the approach of score.py) --------------

_STATE: dict = {}


def _init(state):
    _STATE.update(state)


def _work(job):
    gid, wi = job
    s = _STATE
    try:
        res = score_window(s["guides"][gid], s["windows"][wi][0], s["cfg"]) or {}
    except Exception as exc:  # one bad window must not kill a run of hours
        return job, {"_error": f"{type(exc).__name__}: {exc}"[:200]}
    if gid in s["decoy_ids"]:      # only the energy of a decoy is used: keeps millions of results small
        res = {"intarna_E": res["intarna_E"]} if "intarna_E" in res else {}
    return job, res


def run_signature(guides: dict, windows: list, cfg: BacterialConfig) -> str:
    h = hashlib.sha1()
    h.update(repr(sorted(guides.items())).encode())
    h.update(repr(windows).encode())
    h.update(repr((cfg.temperature, cfg.seed_bp, cfg.no_seed, cfg.acc_window, cfg.acc_length,
                   cfg.intarna_extra)).encode())
    return h.hexdigest()


class Journal:
    """Checkpoint of a run as an append-only series of pickles in ``directory``: every save writes
    only the results obtained since the previous one, so saving stays cheap for millions of
    results. A run with another ``signature`` (other inputs or settings) ignores and replaces
    the files. The files are pickles written by this program; do not point it at files from
    elsewhere."""

    def __init__(self, directory: str, signature: str):
        self.dir, self.signature, self.n = directory, signature, 0
        self.sig_file = os.path.join(directory, "signature.txt")

    def _segments(self):
        return sorted(f for f in os.listdir(self.dir) if f.startswith("seg-") and f.endswith(".pkl"))

    def load(self, resume: bool = True) -> dict:
        """The saved results of a run with this signature ({} when ``resume`` is false or the
        files belong to another run, which are then replaced)."""
        os.makedirs(self.dir, exist_ok=True)
        try:
            with open(self.sig_file) as fh:
                same = resume and fh.read().strip() == self.signature
        except OSError:
            same = False
        out: dict = {}
        if same:
            for name in self._segments():
                try:
                    with open(os.path.join(self.dir, name), "rb") as fh:
                        out.update(pickle.load(fh))
                except (OSError, EOFError, pickle.UnpicklingError, AttributeError, ImportError):
                    continue     # a segment cut short by a kill: its jobs are scored again
        else:
            for name in self._segments():
                os.remove(os.path.join(self.dir, name))
            with open(self.sig_file, "w") as fh:
                fh.write(self.signature + chr(10))
        self.n = len(self._segments())
        return {j: r for j, r in out.items() if "_error" not in r}    # failed ones are retried

    def append(self, items: dict) -> None:
        if not items:
            return
        self.n += 1
        path = os.path.join(self.dir, f"seg-{self.n:06d}.pkl")
        with open(path + ".tmp", "wb") as fh:
            pickle.dump(items, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(path + ".tmp", path)


def run_jobs(jobs, state, cfg: BacterialConfig, done: dict, journal: Journal, log,
             label="windows "):
    """Score the jobs missing from ``done`` (filled in place) with ``cfg.workers`` processes and
    save new results to the journal at most every ``cfg.checkpoint_seconds`` and at the end."""
    todo = [j for j in jobs if j not in done]
    n, t0, last_save = len(todo), time.time(), time.time()
    step = max(50, cfg.workers * 20)
    pending: dict = {}
    pool = None
    try:
        if cfg.workers > 1:
            pool = ProcessPoolExecutor(cfg.workers, initializer=_init, initargs=(state,))
        else:
            _init(state)
        for lo in range(0, n, step):
            part = todo[lo:lo + step]
            results = pool.map(_work, part, chunksize=4) if pool else map(_work, part)
            for job, res in results:
                done[job] = pending[job] = res
            if time.time() - last_save >= cfg.checkpoint_seconds:
                journal.append(pending)
                pending, last_save = {}, time.time()
                rate = (lo + len(part)) / (last_save - t0)
                log(f"{label}{lo + len(part)}/{n} ({(lo + len(part)) / n:.0%}), "
                    f"ETA {(n - lo - len(part)) / rate / 60:.0f} min")
    finally:
        if pool:
            pool.shutdown(wait=True, cancel_futures=True)
    journal.append(pending)


def summarise_errors(scores: dict, log, max_share: float) -> None:
    errors = [(job, res["_error"]) for job, res in scores.items() if "_error" in res]
    if not errors:
        return
    log(f"WARNING: {len(errors)} of {len(scores)} windows failed and are left unscored; "
        f"first: {errors[0][1]}")
    if len(errors) > max_share * len(scores):
        raise RuntimeError(f"{len(errors)}/{len(scores)} windows failed (> {max_share:.0%}); "
                           f"first error: {errors[0][1]}")
    for job, _ in errors:
        scores[job] = {}   # unscored


def check_tools(cfg: BacterialConfig) -> None:
    if shutil.which(cfg.intarna_bin) is None:
        raise SystemExit(f"score-windows: not on PATH: {cfg.intarna_bin} "
                         "(activate the environment with IntaRNA)")


# ---- the run ---------------------------------------------------------------------------

def curve_rows(real, decoy, weights):
    """FDR along the score: one row per distinct real score. ``weights`` are the decoy weights
    (1 / decoys of the candidate), so the column is the expected number of false windows per
    candidate. Same layout as the fdr_curves.tsv of score-sites (site_class = ``window``)."""
    rows, r = [], sorted(real)
    pairs = sorted(zip(decoy, weights))
    d = [x for x, _ in pairs]
    cum = list(accumulate(w for _, w in pairs))
    for x in sorted(set(r)):
        pr, k = bisect_right(r, x), bisect_right(d, x)
        false = cum[k - 1] if k else 0.0
        rows.append({"predictor": "intarna", "site_class": "window", "threshold": x,
                     "real_sites": pr, "decoy_sites_per_candidate": round(false, 4),
                     "fdr": round(min(1.0, false / pr), 4)})
    return rows


def genomic_span(r: Region, t_start: int, t_end: int) -> str:
    """``contig:start-end(strand)`` (1-based) of window positions t_start..t_end (1-based), as
    plus-strand coordinates. Start > end means the interaction crosses the origin of a circular
    contig (it runs clockwise from start to end)."""
    g, g_end = r.genomic(t_start - 1), r.genomic(t_end - 1)
    if g is None or g_end is None:
        return ""
    a, b = (g, g_end) if r.strand == "+" else (g_end, g)    # on '-' the window runs backwards
    return f"{r.contig}:{a + 1}-{b + 1}({r.strand})"


def plan_jobs(guides: dict, n_windows: int):
    return [(gid, wi) for gid in guides for wi in range(n_windows)]


def dry_run_summary(guides, n_windows: int, cfg: BacterialConfig) -> str:
    n = len(guides) * n_windows
    return (f"{n_windows} unique windows x {len(guides)} small RNAs (real + decoys) = {n} IntaRNA "
            f"runs; rough time at 0.15 s per run on {cfg.workers} worker(s): "
            f"{n * 0.15 / max(1, cfg.workers) / 3600:.1f} h")


def run_window_scoring(candidates, decoys, regions: list[Region],
                       cfg: BacterialConfig | None = None, workdir: str = ".",
                       dry_run: bool = False, log=print):
    """Score every candidate and decoy against every window. ``candidates`` are (id, sequence),
    ``decoys`` (id, candidate, sequence), ``regions`` the start windows (with ``anchor``).
    Returns (window_rows, target_rows, curves), or None for a dry run."""
    cfg = cfg or BacterialConfig()
    if not cfg.no_seed and cfg.seed_bp < 2:
        raise ValueError("seed_bp must be at least 2 (use no_seed for no seed)")
    missing = [r.id for r in regions if r.anchor is None]
    if missing:
        raise ValueError(f"{missing[0]}: not a start-codon window (no anchor); make the "
                         "regions with `regions --mode bacterial`")
    by_seq: dict[tuple, list[Region]] = defaultdict(list)
    for r in regions:
        by_seq[(r.sequence, r.anchor)].append(r)
    windows = list(by_seq)

    guides = {cid: seq for cid, seq in candidates}
    decoy_of = {}
    for did, cid, seq in decoys:
        if cid in guides:
            guides[did] = seq
            decoy_of[did] = cid
    n_decoys = Counter(decoy_of.values())
    if not n_decoys:
        raise ValueError("score-windows needs decoys (decoys.fasta has none for these "
                         "candidates): without them there is no false discovery rate")
    without = [cid for cid, _ in candidates if n_decoys[cid] == 0]
    if without:
        log(f"WARNING: {len(without)} candidate(s) have no decoys and get no q-values "
            f"(e.g. {without[0]})")
    short = [g for g, s in guides.items() if not cfg.no_seed and len(s) < cfg.seed_bp]
    if short:
        raise ValueError(f"{short[0]} is shorter than seed_bp ({cfg.seed_bp})")

    log(dry_run_summary(guides, len(windows), cfg))
    if dry_run:
        return None
    check_tools(cfg)

    os.makedirs(workdir, exist_ok=True)
    journal = Journal(os.path.join(workdir, "windows_checkpoint"), run_signature(guides, windows, cfg))
    scores = journal.load(cfg.resume)
    jobs = plan_jobs(guides, len(windows))
    if scores:
        log(f"resumed {len(scores)} of {len(jobs)} window scores from {journal.dir}")
    t0 = time.time()
    state = {"guides": guides, "windows": windows, "cfg": cfg, "decoy_ids": set(decoy_of)}
    run_jobs(jobs, state, cfg, scores, journal, log)
    log(f"scored {len(jobs)} (small RNA, window) pairs in {time.time() - t0:.0f} s")
    summarise_errors(scores, log, cfg.max_error_share)   # the saved checkpoint keeps the errors

    # q-values: one best interaction per (small RNA, window)
    real, real_keys, decoy, weights = [], [], [], []
    real_by_cand: dict[str, list[tuple]] = defaultdict(list)
    decoy_by_cand: dict[str, list[float]] = defaultdict(list)
    for job, res in scores.items():
        e = res.get("intarna_E")
        if e is None:
            continue
        parent = decoy_of.get(job[0])
        if parent is not None:
            decoy.append(e)
            weights.append(1.0 / n_decoys[parent])
            decoy_by_cand[parent].append(e)
        elif n_decoys[job[0]] > 0:
            real.append(e)
            real_keys.append(job)
            real_by_cand[job[0]].append(job)
    log(f"intarna: {sum(1 for r in scores.values() if r.get('intarna_E') is None)} of "
        f"{len(scores)} pairs without an interaction")
    qpool, qcand, curves = {}, {}, []
    if real:
        plus = (len(n_decoys) / sum(n_decoys.values())) if cfg.fdr_pseudocount else 0.0
        for key, q in zip(real_keys, qvalues(real, decoy, 0, decoy_weights=weights,
                                             plus_one=plus)):
            qpool[key] = q
        curves = curve_rows(real, decoy, weights)
        for cid, keys in real_by_cand.items():
            plus_c = 1.0 / n_decoys[cid] if cfg.fdr_pseudocount else 0.0
            for k, q in zip(keys, qvalues([scores[k]["intarna_E"] for k in keys],
                                          decoy_by_cand.get(cid, []), n_decoys[cid],
                                          plus_one=plus_c)):
                qcand[k] = q

    rows = []
    for job, res in scores.items():
        gid, wi = job
        if gid in decoy_of:
            continue
        seq, anchor = windows[wi]
        extra = derive(res, anchor, cfg)
        q = qpool.get(job, "")
        for r in by_seq[windows[wi]]:
            row = {"candidate": gid, "region_id": r.id, "transcript": r.transcript,
                   "gene": r.gene, "region": r.region, "contig": r.contig, "strand": r.strand,
                   "window_length": len(seq), "upstream_cds": r.upstream_cds,
                   "cds_nt": r.cds_nt, **extra, **res,
                   "intarna_q": q, "intarna_q_cand": qcand.get(job, ""),
                   "pass_fdr": int(q != "" and q <= cfg.fdr)}
            if res.get("intarna_t_start") is not None:
                row["genomic"] = genomic_span(r, res["intarna_t_start"], res["intarna_t_end"])
            rows.append(row)
    rows.sort(key=lambda r: (r["candidate"], r["intarna_q"] == "",
                             r["intarna_q"] if r["intarna_q"] != "" else 1.0,
                             r.get("intarna_E") if r.get("intarna_E") is not None else 0.0))
    return rows, summarise_targets(rows), curves


def summarise_targets(rows):
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        groups[(r["candidate"], r["transcript"], r["gene"])].append(r)
    out = []
    for (cid, transcript, gene), rs in groups.items():
        scored = [r for r in rs if r.get("intarna_E") is not None]
        top = min(scored, key=lambda r: r["intarna_E"]) if scored else {}

        def best_q(col):
            vals = [r[col] for r in rs if r.get(col) not in (None, "")]
            return min(vals) if vals else ""

        out.append({"candidate": cid, "transcript": transcript, "gene": gene,
                    "scored_sites": len(scored), "best_intarna_E": top.get("intarna_E", ""),
                    "best_intarna_q": best_q("intarna_q"),
                    "best_intarna_q_cand": best_q("intarna_q_cand"),
                    "hit_region": top.get("hit_region", ""), "rel_start": top.get("rel_start", ""),
                    "rel_end": top.get("rel_end", ""), "hit_rbs": top.get("hit_rbs", ""),
                    "windows_q_pass": sum(r["pass_fdr"] for r in rs)})
    out.sort(key=lambda r: (r["candidate"], r["best_intarna_q"] == "",
                            r["best_intarna_q"] if r["best_intarna_q"] != "" else 1.0))
    return out


WINDOW_COLUMNS = ["candidate", "region_id", "transcript", "gene", "region", "contig", "strand",
                  "window_length", "upstream_cds", "cds_nt", "genomic", "rel_start", "rel_end",
                  "hit_region", "hit_rbs",
                  "intarna_E", "intarna_ED_target", "intarna_ED_query", "intarna_Pu_target",
                  "intarna_Pu_query", "intarna_E_hybrid", "intarna_t_start", "intarna_t_end",
                  "intarna_q_start", "intarna_q_end", "intarna_duplex", "intarna_q",
                  "intarna_q_cand", "pass_fdr"]
TARGET_COLUMNS = ["candidate", "transcript", "gene", "scored_sites", "best_intarna_E",
                  "best_intarna_q", "best_intarna_q_cand", "hit_region", "rel_start", "rel_end",
                  "hit_rbs", "windows_q_pass"]
CURVE_COLUMNS = ["predictor", "site_class", "threshold", "real_sites",
                 "decoy_sites_per_candidate", "fdr"]


def write_windows(window_rows, target_rows, curves, windows_path, targets_path, curves_path):
    write_tsv(windows_path, window_rows, WINDOW_COLUMNS)
    write_tsv(targets_path, target_rows, TARGET_COLUMNS)
    write_tsv(curves_path, curves, CURVE_COLUMNS)
