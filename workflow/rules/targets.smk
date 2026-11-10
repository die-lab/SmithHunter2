wildcard_constraints:
    tset="[^/]+",


def regions_args(wildcards):
    ts = TARGET_SETS[wildcards.tset]
    wanted = ts.get("regions", ["three_prime_UTR"])
    args = [f"--regions {','.join(wanted)}"]
    if ts.get("genome"):
        args += [f"--genome {ts['genome']}", f"--annotation {ts['annotation']}"]
    elif ts.get("transcripts"):
        args.append(f"--transcripts {ts['transcripts']}")
        if ts.get("region_table"):
            args.append(f"--table {ts['region_table']}")
    else:
        raise ValueError(f"config: target set {wildcards.tset!r} needs genome + annotation "
                         "or transcripts")
    if ts.get("whole_transcript_if_no_utr"):
        args.append("--whole-if-missing")
    return " ".join(args)


def regions_inputs(wildcards):
    ts = TARGET_SETS[wildcards.tset]
    keys = ("genome", "annotation", "transcripts", "region_table")
    return [ts[k] for k in keys if ts.get(k)]


rule target_regions:
    input:
        regions_inputs,
    output:
        OUT / "targets/{tset}/regions.tsv",
    params:
        args=regions_args,
    log:
        OUT / "logs/targets/{tset}.regions.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} regions {params.args} --out {output} > {log} 2>&1"


rule decoys:
    input:
        CANDIDATES,
    output:
        OUT / "targets/decoys.fasta",
    params:
        n=TARGETS_CFG.get("decoys", 20),
        mode=TARGETS_CFG.get("decoy_mode", "shuffle"),
    log:
        OUT / "logs/targets/decoys.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} decoys --candidates {input} --n {params.n} --seed {SEED} --decoy-mode {params.mode} --out {output} "
        "> {log} 2>&1"


rule seed_sites:
    input:
        candidates=CANDIDATES,
        decoys=OUT / "targets/decoys.fasta",
        regions=OUT / "targets/{tset}/regions.tsv",
    output:
        sites=OUT / "targets/{tset}/sites.tsv",
        targets=OUT / "targets/{tset}/targets.tsv",
        fdr=OUT / "targets/{tset}/site_fdr.tsv",
    params:
        genomes=lambda wc: ",".join(genomes_of_set(wc.tset)),
        classes=",".join(TARGETS_CFG.get("site_classes", [])),
    log:
        OUT / "logs/targets/{tset}.sites.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} sites --candidates {input.candidates} --decoys {input.decoys} "
        "--regions {input.regions} --genomes {params.genomes} --seed {SEED} "
        "--classes '{params.classes}' --sites {output.sites} --targets {output.targets} "
        "--fdr {output.fdr} > {log} 2>&1"


PREDICTORS = TARGETS_CFG.get("predictors", ["intarna", "rnahybrid"])


rule score_sites:
    input:
        candidates=CANDIDATES,
        decoys=OUT / "targets/decoys.fasta",
        regions=OUT / "targets/{tset}/regions.tsv",
    output:
        sites=OUT / "targets/{tset}/sites_scored.tsv",
        targets=OUT / "targets/{tset}/targets_scored.tsv",
        curves=OUT / "targets/{tset}/fdr_curves.tsv",
    params:
        genomes=lambda wc: ",".join(genomes_of_set(wc.tset)),
        classes=",".join(TARGETS_CFG.get("site_classes", [])),
        predictors=",".join(PREDICTORS),
        weak="--score-weak-sites" if TARGETS_CFG.get("score_weak_sites", False) else "",
        temperature=("--temperature " + str(TARGETS_CFG["temperature"])
                     if "temperature" in TARGETS_CFG else ""),
        calibrate=TARGETS_CFG.get("calibrate_samples", 1000),
        fdr=TARGETS_CFG.get("fdr", 0.05),
        prefilter=TARGETS_CFG.get("intarna_prefilter", 1.0),
        workdir=lambda wc: OUT / f"targets/{wc.tset}/work",
    threads: 8
    log:
        OUT / "logs/targets/{tset}.score.log",
    conda:
        "../envs/rna.yaml"
    shell:
        "{SH} score-sites --candidates {input.candidates} --decoys {input.decoys} "
        "--regions {input.regions} --genomes {params.genomes} --seed {SEED} "
        "--classes '{params.classes}' --predictors {params.predictors} {params.weak} "
        "{params.temperature} --calibrate-samples {params.calibrate} "
        "--fdr {params.fdr} --intarna-prefilter {params.prefilter} --workers {threads} --workdir {params.workdir} "
        "--sites {output.sites} --targets {output.targets} --curves {output.curves} "
        "> {log} 2>&1"


if BACTERIAL_SETS:
    print(f"SmithHunter2: bacterial mode is not implemented yet; target sets "
          f"{', '.join(BACTERIAL_SETS)} are skipped.", file=sys.stderr)


rule targets:
    input:
        expand(OUT / "targets/{tset}/{f}.tsv", tset=AGO_SETS, f=["sites", "targets", "site_fdr"]),
        expand(OUT / "targets/{tset}/{f}.tsv", tset=AGO_SETS if PREDICTORS else [],
               f=["sites_scored", "targets_scored", "fdr_curves"]),
