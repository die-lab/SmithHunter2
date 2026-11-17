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
    log:
        OUT / "logs/targets/decoys.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} decoys --candidates {input} --n {params.n} --seed {SEED} --out {output} "
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


if BACTERIAL_SETS:
    print(f"SmithHunter2: bacterial mode is not implemented yet; target sets "
          f"{', '.join(BACTERIAL_SETS)} are skipped.", file=sys.stderr)


PREDICTORS = TARGETS_CFG.get("predictors", ["intarna", "rnahybrid"]) or []
unknown_predictors = set(PREDICTORS) - {"intarna", "rnahybrid"}
if unknown_predictors:
    raise ValueError(f"config: targets.predictors: unknown {sorted(unknown_predictors)} "
                     "(available: intarna, rnahybrid)")


def score_args(wildcards):
    args = []
    if TARGETS_CFG.get("score_weak_sites"):
        args.append("--score-weak-sites")
    if TARGETS_CFG.get("temperature") is not None:
        args.append(f"--temperature {TARGETS_CFG['temperature']}")
    for key, flag in (("fdr", "--fdr"), ("intarna_context", "--context"),
                      ("calibration_windows", "--calibration-windows")):
        if TARGETS_CFG.get(key) is not None:
            args.append(f"{flag} {TARGETS_CFG[key]}")
    return " ".join(args)


# Step 2: energy predictors on the seed sites of candidates and decoys, q-values, support.
rule score_sites:
    input:
        candidates=CANDIDATES,
        decoys=OUT / "targets/decoys.fasta",
        regions=OUT / "targets/{tset}/regions.tsv",
        sites=OUT / "targets/{tset}/sites.tsv",
    output:
        sites=OUT / "targets/{tset}/scored_sites.tsv",
        targets=OUT / "targets/{tset}/scored_targets.tsv",
        decoys=OUT / "targets/{tset}/decoy_scores.tsv",
        curves=OUT / "targets/{tset}/fdr_curves.tsv",
    params:
        genomes=lambda wc: ",".join(genomes_of_set(wc.tset)),
        classes=",".join(TARGETS_CFG.get("site_classes", [])),
        predictors=",".join(PREDICTORS),
        args=score_args,
    threads: config.get("threads", 4)
    log:
        OUT / "logs/targets/{tset}.score.log",
    conda:
        "../envs/rna.yaml"
    shell:
        "{SH} score-sites --candidates {input.candidates} --decoys {input.decoys} "
        "--regions {input.regions} --sites {input.sites} --genomes {params.genomes} "
        "--seed {SEED} --classes '{params.classes}' --predictors {params.predictors} "
        "{params.args} --threads {threads} --scored-sites {output.sites} "
        "--scored-targets {output.targets} --decoy-scores {output.decoys} "
        "--curves {output.curves} > {log} 2>&1"


rule targets:
    input:
        expand(OUT / "targets/{tset}/{f}.tsv", tset=AGO_SETS, f=["sites", "targets", "site_fdr"]),
        expand(OUT / "targets/{tset}/{f}.tsv", tset=AGO_SETS if PREDICTORS else [],
               f=["scored_sites", "scored_targets", "fdr_curves"]),
