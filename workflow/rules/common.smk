import csv
import json
from pathlib import Path

OUT = Path(config.get("outdir", "results"))
SRC = Path(workflow.basedir).parent / "src"
SH = f"PYTHONPATH={SRC} python -m smithhunter"

TRIM = config.get("trimming", {})
MIN_LEN = int(TRIM.get("min_length", 18))
MAX_LEN = int(TRIM.get("max_length", 35))
MAPPING = config.get("mapping", {})
LOCI = config.get("loci", {})
ENDS = config.get("ends", {})
THREADS = int(config.get("threads", 8))
GENOMES = config["genomes"]

# Reads are kept if they fit the widest length range among genomes.
GLOBAL_MIN = min([MIN_LEN] + [int(g["length"][0]) for g in GENOMES.values() if g.get("length")])
GLOBAL_MAX = max([MAX_LEN] + [int(g["length"][1]) for g in GENOMES.values() if g.get("length")])


def load_samples(path):
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("sample")]
    samples = {}
    for r in rows:
        name = r["sample"].strip()
        if name in samples:
            raise ValueError(f"{path}: duplicated sample {name!r}")
        if not r.get("fq1"):
            raise ValueError(f"{path}: sample {name!r} has no fq1")
        samples[name] = {"fq1": r["fq1"].strip(), "fq2": (r.get("fq2") or "").strip()}
    if not samples:
        raise ValueError(f"{path}: no samples")
    return samples


SAMPLES = load_samples(config["samples"])


def fastq_inputs(wildcards):
    s = SAMPLES[wildcards.sample]
    return {"fq1": s["fq1"], "fq2": s["fq2"]} if s["fq2"] else {"fq1": s["fq1"]}


def fastp_args(wildcards, input, output):
    s = SAMPLES[wildcards.sample]
    args = [f"-i {input.fq1}"]
    if s["fq2"]:
        # Small RNA inserts are shorter than the reads: merge the pair and keep the merged read.
        args += [f"-I {input.fq2}", "--merge", f"--merged_out {output.reads}"]
        if TRIM.get("adapter_r2"):
            args.append(f"--adapter_sequence_r2 {TRIM['adapter_r2']}")
    else:
        args.append(f"-o {output.reads}")
    if TRIM.get("adapter_r1"):
        args.append(f"--adapter_sequence {TRIM['adapter_r1']}")
    args += [f"--length_required {max(10, GLOBAL_MIN - 2)}", "--trim_poly_g"]
    if TRIM.get("extra"):
        args.append(TRIM["extra"])
    return " ".join(args)
