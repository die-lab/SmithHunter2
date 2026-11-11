"""Which predictors run, which of them count, and whether they suit the data (module B step 2).

IntaRNA and RNAhybrid are not equally good on every kind of data (bacterial sRNA targets:
IntaRNA is far better; miRNA-like seed sites: both carry some signal, see the development
log, section 6). This module keeps the policy in one place and never throws a call away
silently:

- ``resolve_predictors``: the predictors that are COMPUTED for a target set. Precedence:
  the target set, then the global ``targets.predictors``, then the default of its mode
  (``ago`` both, ``bacterial`` IntaRNA only, ``cross-kingdom`` both). A predictor that is not
  computed costs nothing and its columns stay empty (never 0).
- ``counted_predictors``: the predictors that COUNT when calls are formed (``support_from``).
  A predictor can be computed and reported but ignored. Calls are the UNION: a site or gene is
  called when ANY counted predictor has q <= fdr, so a disagreeing tool never removes a call.
  ``support`` is a count of agreeing predictors, not a confidence.
- ``site_calls``: ``support``, ``called_by`` and ``best_q`` of one site.
- the prefilter helpers (``prefilter_keep``, ``prefilter_summary``) and ``check_config``,
  which says out loud what the RNAhybrid prefilter costs in IntaRNA-only candidates.
- ``diagnose``: per predictor and site class, whether its scores separate the real sites from
  the decoy sites (``predictor_report.tsv``). It only ADVISES; nothing is dropped from calls
  on its basis: the user decides through ``support_from`` and ``predictors``.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from itertools import accumulate

from ..seqio import write_tsv

SCORE_COLUMN = {"intarna": "intarna_E", "rnahybrid": "rnahybrid_mfe"}
KNOWN = tuple(SCORE_COLUMN)
MODES = ("ago", "bacterial", "cross-kingdom")
# Evidence: development log, section 6 (bacterial: IntaRNA with local accessibility does far
# better than RNAhybrid alone; miRNA-like and cross-kingdom: both, then ask what the data say).
MODE_PREDICTORS = {"ago": ("intarna", "rnahybrid"),
                   "bacterial": ("intarna",),
                   "cross-kingdom": ("intarna", "rnahybrid")}

AUC_TOLERANCE = 0.02   # |AUC - 0.5| within this and no call: "no separation"
MIN_SCORED = 30        # fewer scored real sites than this: the AUC is not reliable
LOW_COVERAGE = 0.9     # share of real sites with a score below which the report says so

REPORT_COLUMNS = ["predictor", "site_class", "counted", "real_sites", "real_scored",
                  "scored_share", "decoy_scored", "median_real", "median_decoy", "auc",
                  "real_q_le_fdr", "verdict", "note"]
PREFILTER_COLUMNS = ["site_class", "group", "sites", "rnahybrid_scored", "no_rnahybrid_score",
                     "rnahybrid_cutoff", "passed_cutoff", "kept_by_class_rule",
                     "scored_by_intarna", "never_scored_by_intarna", "share_never_scored", "note"]

PREFILTER_WARNING = (
    "WARNING: the RNAhybrid prefilter decides which sites IntaRNA sees. A site with a poor "
    "RNAhybrid MFE (or no RNAhybrid duplex on the seed) is NEVER scored by IntaRNA, so a "
    "candidate that only IntaRNA would find is lost, and the decoys cannot measure that loss. "
    "Do not use the prefilter when IntaRNA is the better-suited tool for your data (e.g. "
    "bacterial targets); prefilter_keep_classes lists the classes that are always scored by "
    "IntaRNA in full.")


# ---- policy ---------------------------------------------------------------------------------

def as_tuple(value) -> tuple[str, ...] | None:
    """A config value (list, 'a,b' string or None) as a tuple of names; None stays None."""
    if value is None:
        return None
    if isinstance(value, str):
        value = value.split(",")
    return tuple(dict.fromkeys(str(v).strip() for v in value if str(v).strip()))


def resolve_predictors(mode: str = "ago", set_value=None, global_value=None) -> tuple[str, ...]:
    """Predictors computed for a target set: the set's own ``predictors``, else the global
    ``targets.predictors``, else the default of its mode."""
    if mode not in MODES:
        raise ValueError(f"unknown mode {mode!r} (one of {', '.join(MODES)})")
    for value in (set_value, global_value):
        chosen = as_tuple(value)
        if chosen is not None:
            unknown = [p for p in chosen if p not in KNOWN]
            if unknown:
                raise ValueError(f"unknown predictor(s): {', '.join(unknown)} "
                                 f"(known: {', '.join(KNOWN)})")
            return chosen
    return MODE_PREDICTORS[mode]


def counted_predictors(predictors, support_from=None) -> tuple[str, ...]:
    """The scored predictors that count when calls are formed. ``support_from`` None = all of
    them. Names that are not scored are ignored here (``check_config`` says so)."""
    predictors = tuple(predictors)
    wanted = as_tuple(support_from)
    if wanted is None:
        return predictors
    unknown = [p for p in wanted if p not in KNOWN]
    if unknown:
        raise ValueError(f"support_from: unknown predictor(s): {', '.join(unknown)}")
    counted = tuple(p for p in predictors if p in wanted)
    if predictors and not counted:
        raise ValueError(f"support_from {list(wanted)} contains none of the scored predictors "
                         f"{list(predictors)}: no site could ever be called")
    return counted


def two_stage(cfg) -> bool:
    return (cfg.intarna_prefilter < 1.0 and "intarna" in cfg.predictors
            and "rnahybrid" in cfg.predictors)


def check_config(cfg, log=print) -> None:
    """Validate the predictor settings (ValueError) and log what they imply."""
    unknown = [p for p in cfg.predictors if p not in KNOWN]
    if unknown:
        raise ValueError(f"unknown predictor(s): {', '.join(unknown)}")
    if cfg.mode not in MODES:
        raise ValueError(f"unknown mode {cfg.mode!r} (one of {', '.join(MODES)})")
    if cfg.prefilter_by != "rnahybrid":
        raise ValueError(f"prefilter_by {cfg.prefilter_by!r}: only 'rnahybrid' is implemented "
                         "(IntaRNA is the expensive predictor, so filtering with it saves nothing)")
    counted = counted_predictors(cfg.predictors, cfg.support_from)
    for p in cfg.predictors:
        if p not in counted:
            log(f"{p} is computed and reported but IGNORED when forming calls "
                f"(support_from = {list(counted)})")
    absent = [p for p in (as_tuple(cfg.support_from) or ()) if p not in cfg.predictors]
    if absent:
        log(f"support_from names {absent}, which {'is' if len(absent) == 1 else 'are'} not "
            f"scored here (predictors = {list(cfg.predictors)}); ignored")
    if cfg.mode == "bacterial" and "rnahybrid" in cfg.predictors:
        log("note: RNAhybrid alone is much worse than IntaRNA for bacterial targets; its score is "
            "reported, add support_from: [intarna] to keep it out of the calls")
    if cfg.intarna_prefilter < 1.0:
        if not two_stage(cfg):
            log(f"WARNING: intarna_prefilter {cfg.intarna_prefilter} is ignored: the prefilter "
                f"needs both predictors (predictors = {list(cfg.predictors)})")
        else:
            log(PREFILTER_WARNING)
            if cfg.mode == "bacterial":
                log("WARNING: the prefilter is not meant for bacterial targets (IntaRNA is the "
                    "better tool there); set intarna_prefilter: 1.0")
            if "rnahybrid" not in counted:
                log("WARNING: RNAhybrid is ignored for calls but still decides which sites "
                    "IntaRNA sees (prefilter)")


def calls_of(qs: dict, counted, fdr: float) -> tuple[int, str, float | str]:
    """``(support, called_by, best_q)`` from ``qs`` = {predictor: pooled q or None}.

    support   number of counted predictors with q <= fdr (a count of agreeing predictors, not
              a confidence);
    called_by those predictors, comma separated, in the order of ``counted``;
    best_q    the smallest q over the counted predictors that have one, '' if none has.
    The union rule: one counted predictor is enough to call a site."""
    called = [p for p in counted if qs.get(p) is not None and qs[p] <= fdr]
    present = [qs[p] for p in counted if qs.get(p) is not None]
    return len(called), ",".join(called), (min(present) if present else "")


# ---- prefilter ------------------------------------------------------------------------------

def passes_prefilter(job, res, cutoff, keep_classes) -> bool:
    """True when IntaRNA scores this site in two-stage mode: the class is always scored, or the
    RNAhybrid MFE is at or below the pooled cutoff of the class."""
    if job[5] in keep_classes:
        return True
    mfe = res.get("rnahybrid_mfe")
    return mfe is not None and mfe <= cutoff[job[5]]


def prefilter_keep(jobs, scores, cutoff, keep_classes=()) -> list:
    """The jobs IntaRNA scores in two-stage mode (see ``passes_prefilter``)."""
    keep_classes = frozenset(keep_classes)
    return [j for j in jobs if passes_prefilter(j, scores[j], cutoff, keep_classes)]


def prefilter_summary(scores, decoy_of, cutoff, keep_classes, classes) -> list[dict]:
    """Per class and group (real / decoy), how many sites IntaRNA never saw, recomputed from the
    final scores (the RNAhybrid values are all there, so this also works for a resumed run)."""
    keep_classes = frozenset(keep_classes)
    acc: dict[tuple, dict] = {}
    for job, res in scores.items():
        group = "decoy" if job[0] in decoy_of else "real"
        a = acc.setdefault((job[5], group), dict.fromkeys(
            ("sites", "scored", "passed"), 0))
        a["sites"] += 1
        if res.get("rnahybrid_mfe") is not None:
            a["scored"] += 1
            if res["rnahybrid_mfe"] <= cutoff.get(job[5], float("-inf")):
                a["passed"] += 1
    rows = []
    for cls in classes:
        for group in ("real", "decoy"):
            a = acc.get((cls, group))
            if not a:
                continue
            # IntaRNA ran on the sites that passed the cutoff and, for an always-scored class, on
            # all the others too (``rescued``: sites it would have missed without that rule)
            ran = a["sites"] if cls in keep_classes else a["passed"]
            by_rule = ran - a["passed"]
            never = a["sites"] - ran
            rows.append({
                "site_class": cls, "group": group, "sites": a["sites"],
                "rnahybrid_scored": a["scored"], "no_rnahybrid_score": a["sites"] - a["scored"],
                "rnahybrid_cutoff": cutoff.get(cls, ""), "passed_cutoff": a["passed"],
                "kept_by_class_rule": by_rule, "scored_by_intarna": ran,
                "never_scored_by_intarna": never,
                "share_never_scored": round(never / a["sites"], 4),
                "note": "class always scored by IntaRNA" if cls in keep_classes else ""})
    return rows


def prefilter_off_rows(scores, decoy_of, cfg) -> list[dict]:
    """The summary when no prefilter ran: every site went to every predictor in use."""
    n_real = sum(1 for job in scores if job[0] not in decoy_of)
    note = ("no prefilter: every site was sent to every computed predictor"
            if "intarna" in cfg.predictors else "IntaRNA was not computed")
    return [{"site_class": "all", "group": "real", "sites": n_real, "note": note},
            {"site_class": "all", "group": "decoy", "sites": len(scores) - n_real, "note": note}]


def log_prefilter(rows, cfg, log) -> None:
    real = [r for r in rows if r["group"] == "real" and r["site_class"] != "all"]
    decoy = [r for r in rows if r["group"] == "decoy" and r["site_class"] != "all"]
    n, never = sum(r["sites"] for r in real), sum(r["never_scored_by_intarna"] for r in real)
    nd, neverd = sum(r["sites"] for r in decoy), sum(r["never_scored_by_intarna"] for r in decoy)
    nors = sum(r["no_rnahybrid_score"] for r in real)
    log(f"prefilter cost: {never} of {n} real sites ({never / max(1, n):.1%}) and {neverd} of "
        f"{nd} decoy sites ({neverd / max(1, nd):.1%}) were never scored by IntaRNA "
        f"({nors} real sites had no RNAhybrid score); always scored: "
        f"{list(cfg.prefilter_keep_classes) or 'no class'}. IntaRNA q-values are conditional on "
        "passing the RNAhybrid cutoff; an IntaRNA-only candidate in those sites is lost.")


# ---- suitability diagnostic -------------------------------------------------------------------

def weighted_auc(real, decoy, weights, higher_is_better: bool = False):
    """Probability that a random real site scores better than a random decoy site (ties count
    half; decoy sites weighted as in the q-values), or None without data. 0.5 = no separation."""
    if not real or not decoy:
        return None
    sign = -1.0 if higher_is_better else 1.0
    pairs = sorted(zip((sign * x for x in decoy), weights))
    d = [x for x, _ in pairs]
    cum = [0.0] + list(accumulate(w for _, w in pairs))
    total = cum[-1]
    if total <= 0:
        return None
    acc = 0.0
    for x in real:
        x = sign * x
        lo, hi = bisect_left(d, x), bisect_right(d, x)
        acc += (total - cum[hi]) + 0.5 * (cum[hi] - cum[lo])   # decoys worse + half the ties
    return acc / (len(real) * total)


def weighted_median(values, weights=None):
    if not values:
        return None
    if weights is None:
        s = sorted(values)
        n = len(s)
        return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2
    pairs = sorted(zip(values, weights))
    cum = list(accumulate(w for _, w in pairs))
    return pairs[min(bisect_left(cum, cum[-1] / 2), len(pairs) - 1)][0]


def advise(n_real, n_decoy, auc, calls, share, fdr) -> tuple[str, str]:
    """(verdict code, plain-language note) for one predictor and class."""
    if n_real == 0:
        return "no_scores", "no real site has a score from this predictor: nothing to assess"
    if n_decoy == 0:
        return "no_scores", "no decoy site has a score: no comparison with decoys is possible"
    if n_real < MIN_SCORED:
        verdict, note = "too_few", (f"only {n_real} scored real sites: too few for a reliable "
                                    "comparison with decoys")
    elif auc < 0.5 - AUC_TOLERANCE:
        verdict, note = "worse_than_decoys", (
            f"real sites score WORSE than decoy sites (AUC {auc:.2f}): this predictor does not "
            "suit this data set, or the decoys are not a fair null here (check temperature and "
            "decoy_mode)")
    elif abs(auc - 0.5) <= AUC_TOLERANCE + 1e-12 and calls == 0:
        verdict, note = "no_separation", ("no separation from decoys: this predictor adds "
                                          "nothing for this data set")
    elif abs(auc - 0.5) <= AUC_TOLERANCE + 1e-12:
        verdict, note = "calls_without_separation", (
            f"no overall separation from decoys (AUC {auc:.2f}) but {calls} site(s) reach "
            f"q <= {fdr:g}: look at these calls before trusting them")
    elif calls == 0:
        verdict, note = "separates_no_calls", (
            f"scores are better than decoys on average (AUC {auc:.2f}) but no site reaches "
            f"q <= {fdr:g}: useful to rank genes, adds no calls at this FDR")
    else:
        verdict, note = "separates", (f"separates from decoys (AUC {auc:.2f}) and {calls} site(s) "
                                      f"reach q <= {fdr:g}")
    if share < LOW_COVERAGE:
        note += (f"; only {share:.0%} of the real sites have a score (prefilter or failures), "
                 "so this holds for those sites only")
    return verdict, note


def diagnose(scores, decoy_of, n_decoys, predictors, classes, qpool, fdr,
             counted=None) -> list[dict]:
    """Rows of ``predictor_report.tsv``: for each predictor and site class, the numbers of real
    and decoy sites with a score, the median score of each, the AUC of real vs decoy scores,
    the real sites with q <= ``fdr`` and a note. Uses the same sites, scores and decoy weights
    as the q-values. Sites are those of the unique target sequences (isoforms sharing a UTR
    count once). Advice only: nothing is removed from any call."""
    counted = tuple(predictors) if counted is None else tuple(counted)
    rows = []
    for p in predictors:
        col = SCORE_COLUMN[p]
        acc: dict[str, dict] = defaultdict(lambda: {"all": 0, "real": [], "decoy": [], "w": [],
                                                    "calls": 0})
        for job, res in scores.items():
            a = acc[job[5]]
            v = res.get(col)
            parent = decoy_of.get(job[0])
            if parent is not None:
                if v is not None:
                    a["decoy"].append(v)
                    a["w"].append(1.0 / n_decoys[parent])
            elif n_decoys[job[0]] > 0:       # without decoys a site gets no q-value either
                a["all"] += 1
                if v is not None:
                    a["real"].append(v)
                    q = qpool.get((p, job))
                    if q is not None and q <= fdr:
                        a["calls"] += 1
        for cls in classes:
            a = acc.get(cls)
            if not a or not a["all"]:
                continue
            real, decoy, w = a["real"], a["decoy"], a["w"]
            auc = weighted_auc(real, decoy, w)
            share = len(real) / a["all"]
            verdict, note = advise(len(real), len(decoy), auc, a["calls"], share, fdr)
            med_r, med_d = weighted_median(real), weighted_median(decoy, w)
            rows.append({
                "predictor": p, "site_class": cls, "counted": "yes" if p in counted else "no",
                "real_sites": a["all"], "real_scored": len(real), "scored_share": round(share, 4),
                "decoy_scored": len(decoy),
                "median_real": "" if med_r is None else round(med_r, 3),
                "median_decoy": "" if med_d is None else round(med_d, 3),
                "auc": "" if auc is None else round(auc, 4),
                "real_q_le_fdr": a["calls"], "verdict": verdict,
                "note": note + ("" if p in counted else
                                " [reported only: ignored when forming calls]")})
    return rows


def log_diagnosis(rows, fdr, log) -> None:
    """Three lines for two predictors: a header and one line per predictor."""
    log("predictor suitability (AUC of real vs decoy sites, 0.5 = no separation; advice only, "
        "calls are not changed: use predictors / support_from to act on it):")
    by_p: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_p[r["predictor"]].append(r)
    for p, rs in by_p.items():
        aucs = ", ".join(f"{r['site_class']} {r['auc']:.2f}" for r in rs if r["auc"] != "")
        calls = sum(r["real_q_le_fdr"] for r in rs)
        verdicts = {r["verdict"] for r in rs}
        if verdicts <= {"no_separation"}:
            text = "no separation from decoys: this predictor adds nothing for this data set"
        elif verdicts & {"separates", "separates_no_calls", "calls_without_separation"}:
            text = ("separates from decoys" if verdicts & {"separates", "separates_no_calls"}
                    else "calls without overall separation: check them")
        elif "worse_than_decoys" in verdicts:
            text = "real sites score worse than decoys: probably unsuited to this data"
        else:
            text = next(r["note"] for r in rs)
        log(f"  {p}: AUC {aucs or 'n/a'}; {calls} real site(s) at q <= {fdr:g}"
            f"{'' if rs[0]['counted'] == 'yes' else ' (ignored for calls)'}: {text}")


def write_report(rows, path: str) -> None:
    write_tsv(path, rows, REPORT_COLUMNS)


def write_prefilter(rows, path: str) -> None:
    write_tsv(path, rows, PREFILTER_COLUMNS)
