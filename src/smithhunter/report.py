"""Run report: read funnel, coverage of focus genomes and provenance.

The tables feed the plots and the MultiQC report of rule ``report``
(``workflow/rules/report.smk``). The reporting layer follows the one of the
smithhunter24 port of v0 module A (read funnel, provenance.yaml, ggplot2 plots),
rewritten for loci on several genomes.

Funnel stages, per sample:

- ``raw``: reads in the input (R1 for paired-end data), from the fastp report;
- ``trimmed``: reads left by fastp; ``length_range``: reads in the length range;
- ``mapped``: reads aligned to any genome (every origin class except ``unmapped``);
- ``focus``: reads placed on focus genomes, i.e. in any locus (with ``ambiguous:
  split`` this includes the shares of ambiguous reads);
- ``in_peaks`` / ``in_passing``: reads in peak loci / in loci passing all filters.
"""

from __future__ import annotations

import json
import math
import re
from array import array
from collections import defaultdict

from .loci import load_counts
from .seqio import open_text, read_tsv

FUNNEL_STAGES = (
    ("raw", "Raw reads"),
    ("trimmed", "After trimming"),
    ("length_range", "In the length range"),
    ("mapped", "Aligned to any genome"),
    ("focus", "On focus genomes"),
    ("in_peaks", "In peak loci"),
    ("in_passing", "In passing loci"),
)


def fastp_raw_reads(path: str) -> int:
    with open(path) as fh:
        return int(json.load(fh)["read1_before_filtering"]["total_reads"])


def fastp_version(path: str) -> str:
    with open(path) as fh:
        return str(json.load(fh)["summary"].get("fastp_version", "unknown"))


def funnel(samples: list[str], fastp_jsons: list[str], libraries_path: str, summary_path: str,
           loci_path: str) -> list[dict]:
    """Return one row per funnel stage: stage, label, one count per sample, total."""
    values: dict[str, list[float]] = {stage: [0.0] * len(samples) for stage, _ in FUNNEL_STAGES}
    for i, path in enumerate(fastp_jsons):
        values["raw"][i] = fastp_raw_reads(path)
    libraries = {r["sample"]: r for r in read_tsv(libraries_path)}
    for i, s in enumerate(samples):
        values["trimmed"][i] = float(libraries[s]["reads_in"])
        values["length_range"][i] = float(libraries[s]["reads_kept"])
    for r in read_tsv(summary_path):
        if r["class"] != "unmapped":
            for i, s in enumerate(samples):
                values["mapped"][i] += float(r[s])
    for r in read_tsv(loci_path):
        for i, s in enumerate(samples):
            count = float(r[f"count_{s}"])
            values["focus"][i] += count
            if r["locus_type"] == "peak":
                values["in_peaks"][i] += count
            if r["pass"] == "1":
                values["in_passing"][i] += count
    rows = []
    for stage, label in FUNNEL_STAGES:
        row = {"stage": stage, "label": label}
        for s, v in zip(samples, values[stage]):
            row[s] = _number(v)
        row["total"] = _number(sum(values[stage]))
        rows.append(row)
    return rows


def _number(value: float):
    return int(value) if float(value).is_integer() else round(value, 1)


def write_mqc_funnel(path: str, samples: list[str], rows: list[dict]) -> None:
    """MultiQC custom content: one table row per sample, one column per stage."""
    with open(path, "w") as fh:
        fh.write("# id: 'smithhunter_funnel'\n# section_name: 'SmithHunter2 read funnel'\n"
                 "# description: 'Reads left after each step (fractional where reads are "
                 "shared among alignments).'\n# plot_type: 'table'\n")
        fh.write("Sample\t" + "\t".join(r["label"] for r in rows) + "\n")
        for s in samples:
            fh.write(s + "\t" + "\t".join(str(r[s]) for r in rows) + "\n")


def write_mqc_origin(path: str, samples: list[str], summary_path: str) -> None:
    """MultiQC custom content: origin classes as a stacked bar per sample."""
    groups: dict[str, list[float]] = defaultdict(lambda: [0.0] * len(samples))
    for r in read_tsv(summary_path):
        key = r["class"] if not r["genomes"] else f"{r['class']} ({r['genomes']})"
        for i, s in enumerate(samples):
            groups[key][i] += float(r[s])
    with open(path, "w") as fh:
        fh.write("# id: 'smithhunter_origin'\n# section_name: 'SmithHunter2 read origin'\n"
                 "# description: 'Genome of origin of the reads in the length range.'\n"
                 "# plot_type: 'bargraph'\n# pconfig:\n#     id: 'smithhunter_origin_plot'\n"
                 "#     title: 'Read origin'\n#     ylab: 'Reads'\n")
        fh.write("Sample\t" + "\t".join(groups) + "\n")
        for i, s in enumerate(samples):
            fh.write(s + "\t" + "\t".join(str(_number(v[i])) for v in groups.values()) + "\n")


def coverage(alignments_path: str, counts_path: str, genomes_path: str, contigs_path: str):
    """Depth and 5' end counts per position of focus genome contigs.

    Yields rows (genome, contig, position 1-based, strand, depth, five_ends, per-sample
    depth) for positions with depth > 0. Reads outside the length range of their genome
    are left out, as in loci. Alignments across the origin of a circular contig wrap.
    """
    samples, counts = load_counts(counts_path)
    genomes = {g["name"]: g for g in read_tsv(genomes_path)}
    focus = {n for n, g in genomes.items() if g["role"] == "focus"}
    lengths = {(c["genome"], c["contig"]): int(c["length"]) for c in read_tsv(contigs_path)
               if c["genome"] in focus}
    n = len(samples)
    depth: dict[tuple, list[array]] = {}
    five: dict[tuple, array] = {}

    def tracks(key, length):
        if key not in depth:
            # difference arrays: +w at start, -w after end
            depth[key] = [array("d", bytes(8 * (length + 1))) for _ in range(n)]
            five[key] = array("d", bytes(8 * length))
        return depth[key], five[key]

    with open_text(alignments_path) as fh:
        fh.readline()
        for line in fh:
            rid, genome, contig, start, end, strand, weight, _ = line.rstrip("\n").split("\t")
            if genome not in focus:
                continue
            g = genomes[genome]
            seq, row = counts[rid]
            if not int(g["min_length"]) <= len(seq) <= int(g["max_length"]):
                continue
            length = lengths[(genome, contig)]
            d, f = tracks((genome, contig, strand), length)
            s, e, w = int(start), int(end), float(weight)
            five_pos = s if strand == "+" else e - 1
            f[five_pos % length] += w * sum(row)
            pieces = [(s % length, e - (s - s % length))]
            if pieces[0][1] > length:  # wraps the origin
                pieces = [(pieces[0][0], length), (0, pieces[0][1] - length)]
            for a, b in pieces:
                for i, c in enumerate(row):
                    if c:
                        d[i][a] += w * c
                        d[i][b] -= w * c
    for (genome, contig, strand) in sorted(depth):
        d, f = depth[(genome, contig, strand)], five[(genome, contig, strand)]
        running = [0.0] * n
        for pos in range(lengths[(genome, contig)]):
            for i in range(n):
                running[i] += d[i][pos]
            total = sum(running)
            if total > 1e-9:
                yield ([genome, contig, pos + 1, strand, round(total, 3), round(f[pos], 3)]
                       + [round(x, 3) for x in running])


def coverage_columns(samples: list[str]) -> list[str]:
    return (["genome", "contig", "position", "strand", "depth", "five_ends"]
            + [f"depth_{s}" for s in samples])


def env_pins(paths: list[str]) -> dict[str, str]:
    """Package specs listed in conda env files, e.g. {'fastp': '=0.23.4'}."""
    pins = {}
    for path in paths:
        in_deps = False
        with open(path) as fh:
            for line in fh:
                text = line.split("#")[0].rstrip()
                if text.startswith("dependencies:"):
                    in_deps = True
                    continue
                m = re.match(r"\s*-\s*([A-Za-z0-9_.\-]+)\s*([<>=!].*)?$", text)
                if in_deps and m:
                    pins[m.group(1)] = (m.group(2) or "").strip()
                elif in_deps and text and not text.startswith(" "):
                    in_deps = False
    return pins


def loci_summary(loci_path: str) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = defaultdict(lambda: {"loci": 0, "peak_loci": 0, "passing": 0})
    for r in read_tsv(loci_path):
        g = out[r["genome"]]
        g["loci"] += 1
        g["peak_loci"] += r["locus_type"] == "peak"
        g["passing"] += r["pass"] == "1"
    return dict(out)


_PLAIN_KEY = re.compile(r"^[A-Za-z_][A-Za-z0-9_.\-]*$")


def _scalar(value) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return ".nan"
        if math.isinf(value):
            return ".inf" if value > 0 else "-.inf"
        text = repr(value)
        if "e" in text and "." not in text.split("e")[0]:
            mantissa, exp = text.split("e")
            text = f"{mantissa}.0e{exp}"  # YAML 1.1 reads 1e-05 as a string
        return text
    if isinstance(value, dict):
        return "{}"
    if isinstance(value, list):
        return "[]"
    return json.dumps(str(value))  # a JSON string is a valid YAML string


def _key(key) -> str:
    key = str(key)
    return key if _PLAIN_KEY.match(key) and key.lower() not in ("true", "false", "null", "yes",
                                                                 "no", "on", "off") \
        else json.dumps(key)


def yaml_lines(value, indent: int = 0) -> list[str]:
    """Block-style YAML for dicts, lists and scalars (no third-party dependency)."""
    pad = "  " * indent
    lines = []
    if isinstance(value, dict):
        for k, v in value.items():
            if isinstance(v, (dict, list)) and v:
                lines.append(f"{pad}{_key(k)}:")
                lines += yaml_lines(v, indent + 1)
            else:
                lines.append(f"{pad}{_key(k)}: {_scalar(v)}")
    elif isinstance(value, list):
        for v in value:
            if isinstance(v, (dict, list)) and v:
                nested = yaml_lines(v, indent + 1)
                lines.append(f"{pad}- {nested[0].lstrip()}")
                lines += nested[1:]
            else:
                lines.append(f"{pad}- {_scalar(v)}")
    else:
        lines.append(pad + _scalar(value))
    return lines


def provenance(config: dict, samples: list[str], env_files: list[str], fastp_jsons: list[str],
               loci_path: str, software: dict) -> dict:
    """Everything needed to trace a result back to how it was produced."""
    versions = dict(software)
    if fastp_jsons:
        versions["fastp (ran)"] = fastp_version(fastp_jsons[0])
    return {
        "software": versions,
        "conda_pins": env_pins(env_files),
        "samples": list(samples),
        "results": loci_summary(loci_path),
        "config": config,
    }


def write_yaml(path: str, data: dict, header: str = "") -> None:
    with open(path, "w") as fh:
        if header:
            fh.write("".join(f"# {line}\n" for line in header.splitlines()))
        fh.write("\n".join(yaml_lines(data)) + "\n")
