#!/usr/bin/env python3
"""Cross-study concordance of passing loci on the same genome (a proxy for precision).

    python3 07_concordance.py --run aag2=DIR --run jw18=DIR --run ovary=DIR [--genome wolbachia]
                              [--five-tolerance 2]

DIR is a results dir (DIR/discovery/loci.tsv) of the infected line, e.g. bench/infected/results
or bench/variants/<name>/infected. Studies must use the same genome sequence (same contig names
and coordinates). A passing locus is *concordant* when another study has a passing locus on the
same strand whose representative 5' end lies within --five-tolerance nt (same small RNA), and
*overlapping* when the loci merely overlap. Real processed small RNAs are expected to recur across
tissues and labs; noise, sample-specific degradation and isolated fragments are not, so a setting
that raises the concordant share without losing the known small RNAs is more precise. (Abundant
degradation fragments also recur: read it together with the structural share.)
"""
from __future__ import annotations

import argparse
import importlib.util
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("compare", HERE / "04_compare.py")
compare = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(compare)


def five(r: dict) -> int:
    return int(r["rep_start"]) if r["strand"] == "+" else int(r["rep_end"])


def load(path: Path, genome: str) -> list[dict]:
    return [r for r in compare.read_tsv(path / "discovery/loci.tsv")
            if r["genome"] == genome and r["pass"] == "1"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, help="label=results_dir (repeatable)")
    ap.add_argument("--genome", default="wolbachia")
    ap.add_argument("--five-tolerance", type=int, default=2)
    args = ap.parse_args()
    runs = {}
    for item in args.run:
        label, _, path = item.partition("=")
        runs[label] = load(Path(path).expanduser(), args.genome)
    index = {}   # label -> (contig, strand) -> sorted 5' ends
    for label, rows in runs.items():
        d = defaultdict(list)
        for r in rows:
            d[(r["contig"], r["strand"])].append(five(r))
        index[label] = {k: sorted(v) for k, v in d.items()}

    tol = args.five_tolerance
    out = []
    for label, rows in runs.items():
        others = [o for o in runs if o != label]
        conc = over = 0
        conc_reads = total_reads = 0.0
        for r in rows:
            key, f = (r["contig"], r["strand"]), five(r)
            n = sum(any(abs(f - g) <= tol for g in index[o].get(key, ())) for o in others)
            m = sum(any(compare.overlaps(r, q) for q in runs[o]) for o in others)
            reads = float(r["total_count"])
            total_reads += reads
            if n:
                conc += 1
                conc_reads += reads
            over += bool(m)
        k = len(rows)
        out.append([label, k, f"{conc} ({conc / k:.0%})" if k else "-",
                    f"{over} ({over / k:.0%})" if k else "-",
                    f"{conc_reads / total_reads:.0%}" if total_reads else "-"])
    print(compare.table(["study", "passing", f"concordant (5' within {tol} nt in another study)",
                         "overlapping another study", "reads in concordant loci"], out))


if __name__ == "__main__":
    main()
