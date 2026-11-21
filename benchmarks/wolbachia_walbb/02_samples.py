#!/usr/bin/env python3
"""Split the downloaded runs into the two lines and write one run directory per line.

    BENCH=~/bench_walbb python3 02_samples.py [--threads 8] [--adapter auto|none|SEQ]

Reads $BENCH/runs.tsv (from 01_download.sh), decides for every run whether it is WB2
(wAlbB-infected) or WB2.tet (tetracycline-cleared) from its ENA titles, checks the 3' adapter on
the first reads of each fastq, and writes:

    $BENCH/infected/config/{config.yaml,samples.tsv}   9 WB2 libraries
    $BENCH/tet/config/{config.yaml,samples.tsv}        9 WB2.tet libraries
    $BENCH/groups.tsv                                  run -> line, age (edit and re-run if wrong)
    $BENCH/adapter_check.tsv                           adapter share and read lengths per run

If the titles do not tell the lines apart, the script stops: write groups.tsv by hand
(columns run, group with group = infected | tet) and run it again; an existing groups.tsv wins.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter
from pathlib import Path

# Common 3' adapters of small RNA libraries. The PRJNA671731 runs carry the Illumina universal
# adapter, not the TruSeq small RNA one: checking only the latter reported 0% and left the choice
# to fastp autodetection (which found the right one, 98.7% of the reads trimmed).
ADAPTERS = {
    "truseq_small_rna": "TGGAATTCTCGGGTGCCAAGG",
    "illumina_universal": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",   # also NEBNext small RNA
    "qiaseq_mirna": "AACTGTAGGCACCATCAAT",
    "nextflex_v3_solid": "TCGTATGCCGTCTTCTGCTTG",
}
PROBE_LEN = 12
TRUSEQ = ADAPTERS["truseq_small_rna"]
SMALL_RNA = {"mirna-seq", "ncrna-seq"}


def read_table(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def describe(run: dict) -> str:
    return " | ".join(run.get(k, "") for k in
                      ("sample_title", "sample_alias", "experiment_title", "library_name"))


def guess_group(run: dict) -> str | None:
    text = describe(run).lower()
    if re.search(r"tet|tetracycline|cured|cleared|uninfected|wolbachia[- ]free", text):
        return "tet"
    if re.search(r"wb2|walbb|infected|wolbachia", text):
        return "infected"
    return None


def guess_age(run: dict) -> str:
    m = re.search(r"(\d+)\s*-?\s*(?:d\b|day|days|dpe)", describe(run).lower())
    return m.group(1) if m else ""


def adapter_check(fastq: Path, n: int = 200_000) -> dict:
    """Share of the first n reads containing the first 12 nt of each known adapter."""
    lengths, hits, reads = Counter(), Counter(), 0
    probes = {name: seq[:PROBE_LEN] for name, seq in ADAPTERS.items()}
    with gzip.open(fastq, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            seq = line.strip()
            reads += 1
            lengths[len(seq)] += 1
            for name, probe in probes.items():
                if probe in seq:
                    hits[name] += 1
            if reads >= n:
                break
    best, best_n = (hits.most_common(1)[0] if hits else ("", 0))
    mode_len = lengths.most_common(1)[0][0] if lengths else 0
    return {"reads_checked": reads, "adapter": best,
            "adapter_fraction": round(best_n / max(reads, 1), 3),
            "modal_length": mode_len, "max_length": max(lengths, default=0),
            **{f"fraction_{name}": round(hits[name] / max(reads, 1), 3) for name in ADAPTERS}}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench", default=os.environ.get("BENCH", str(Path.home() / "bench_walbb")))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--adapter", default="auto", help="auto (default), none, or a sequence")
    args = ap.parse_args()
    bench = Path(args.bench).expanduser().resolve()
    here = Path(__file__).resolve().parent

    runs = read_table(bench / "runs.tsv")
    strategies = {r.get("library_strategy", "").lower() for r in runs}
    if strategies & SMALL_RNA:
        dropped = [r["run_accession"] for r in runs if r.get("library_strategy", "").lower() not in SMALL_RNA]
        runs = [r for r in runs if r.get("library_strategy", "").lower() in SMALL_RNA]
        if dropped:
            print(f"not small RNA (library_strategy), ignored: {', '.join(dropped)}")

    groups_path = bench / "groups.tsv"
    if groups_path.exists():
        groups = {r["run"]: r["group"].strip() for r in read_table(groups_path)}
        print(f"using {groups_path}")
    else:
        groups = {r["run_accession"]: guess_group(r) for r in runs}
        with open(groups_path, "w") as fh:
            fh.write("run\tgroup\tage_days\tdescription\n")
            for r in runs:
                fh.write(f"{r['run_accession']}\t{groups[r['run_accession']] or ''}\t{guess_age(r)}\t{describe(r)}\n")
        print(f"wrote {groups_path}")

    for r in runs:
        g = groups.get(r["run_accession"]) or "?"
        print(f"  {r['run_accession']}\t{g:8s}\t{describe(r)}")
    missing = [r["run_accession"] for r in runs if groups.get(r["run_accession"]) not in ("infected", "tet")]
    counts = Counter(groups.get(r["run_accession"]) for r in runs)
    # the paper cites SRR12893564/SRR12893577 as WB2.tet and SRR12893581/SRR12893570 as WB2
    known = {"SRR12893564": "tet", "SRR12893577": "tet", "SRR12893581": "infected", "SRR12893570": "infected"}
    wrong = [f"{k} should be {v}" for k, v in known.items() if k in groups and groups[k] != v]
    if missing or wrong or counts["infected"] != 9 or counts["tet"] != 9:
        print(f"\nSTOP: lines not resolved (infected={counts['infected']}, tet={counts['tet']}, "
              f"unassigned={missing}, inconsistent with the paper: {wrong}).\n"
              f"Fix the 'group' column of {groups_path} (infected | tet) and run again.", file=sys.stderr)
        sys.exit(1)

    rows, fractions = [], []
    for r in runs:
        fq = bench / "fastq" / f"{r['run_accession']}.fastq.gz"
        if not fq.exists():
            sys.exit(f"missing {fq}: run 01_download.sh first")
        res = adapter_check(fq)
        fractions.append(res["adapter_fraction"])
        rows.append({"run": r["run_accession"], "group": groups[r["run_accession"]], **res})
    with open(bench / "adapter_check.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    def median_of(values):
        return sorted(values)[len(values) // 2]

    # the adapter with the highest median share over all runs, and that share
    medians = {name: median_of([r[f"fraction_{name}"] for r in rows]) for name in ADAPTERS}
    name = max(medians, key=medians.get)
    median = medians[name]
    if args.adapter == "auto":
        adapter = ADAPTERS[name] if median >= 0.3 else ""
        why = (f"{name} adapter found in {median:.0%} of reads (median)" if adapter else
               f"no known adapter in more than {median:.0%} of reads: fastp autodetection; "
               "check adapter_check.tsv (already trimmed if modal_length is ~22)")
    else:
        adapter = "" if args.adapter == "none" else args.adapter.upper()
        why = "from --adapter"
    print(f"\nadapter_r1: {adapter or '(autodetect)'}  [{why}]")

    template = (here / "config_template.yaml").read_text()
    for group in ("infected", "tet"):
        run_dir = bench / group
        (run_dir / "config").mkdir(parents=True, exist_ok=True)
        for link in ("data", "fastq"):
            p = run_dir / link
            if not p.exists():
                p.symlink_to(Path("..") / link)
        text = (template.replace("@GROUP@", group).replace("@THREADS@", str(args.threads))
                .replace("@ADAPTER@", adapter))
        (run_dir / "config" / "config.yaml").write_text(text)
        with open(run_dir / "config" / "samples.tsv", "w") as fh:
            fh.write("sample\tfq1\tfq2\n")
            for r in runs:
                if groups[r["run_accession"]] == group:
                    fh.write(f"{r['run_accession']}\t{bench / 'fastq' / (r['run_accession'] + '.fastq.gz')}\t\n")
        print(f"wrote {run_dir}/config/config.yaml and samples.tsv")
    print(f"\nnext: bash {here / '03_run.sh'}")


if __name__ == "__main__":
    main()
