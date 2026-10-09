import csv
import json
import sys
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
    args.append(f"-o {output.reads}")
    if s["fq2"]:
        # Keep R1 only; R2 serves overlap-based error correction. Merging (--merge) is not used:
        # fastp needs >= 30 bp of overlap, so most small RNA inserts would be lost.
        args += [f"-I {input.fq2}", "-O /dev/null", "--correction"]
        if TRIM.get("adapter_r2"):
            args.append(f"--adapter_sequence_r2 {TRIM['adapter_r2']}")
    if TRIM.get("adapter_r1"):
        args.append(f"--adapter_sequence {TRIM['adapter_r1']}")
    args += [f"--length_required {max(10, GLOBAL_MIN - 2)}", "--trim_poly_g"]
    if TRIM.get("extra"):
        args.append(TRIM["extra"])
    return " ".join(args)


# Module B. Every target set uses the general seed-based logic (mode "ago") unless it is
# declared "bacterial". A focus genome is tested against the sets it lists in
# `target_sets`; with no list, against every set (the old `targets: ago | bacterial`
# field selects the sets of that mode, `targets: none` none).
TARGET_SETS = config.get("target_sets") or {}
TARGETS_CFG = config.get("targets") or {}
CANDIDATES = TARGETS_CFG.get("candidates", str(OUT / "discovery/candidates.fasta"))
SEED = "-".join(str(x) for x in TARGETS_CFG.get("seed", [2, 8]))


def set_mode(name):
    mode = TARGET_SETS[name].get("mode", "ago")
    if mode not in ("ago", "bacterial"):
        raise ValueError(f"config: target set {name!r} mode must be ago or bacterial")
    return mode


def target_sets_of(genome):
    entry = GENOMES[genome]
    if entry.get("role", "focus") != "focus":
        return []
    if "target_sets" in entry:
        unknown = set(entry["target_sets"]) - set(TARGET_SETS)
        if unknown:
            raise ValueError(f"config: genome {genome!r} lists unknown target sets {unknown}")
        return list(entry["target_sets"])
    legacy = entry.get("targets")
    if legacy == "none":
        return []
    if legacy in ("ago", "bacterial"):
        return [s for s in TARGET_SETS if set_mode(s) == legacy]
    return list(TARGET_SETS)


def genomes_of_set(name):
    return [g for g in GENOMES if name in target_sets_of(g)]


AGO_SETS = [s for s in TARGET_SETS if set_mode(s) == "ago" and genomes_of_set(s)]
BACTERIAL_SETS = [s for s in TARGET_SETS if set_mode(s) == "bacterial"]
