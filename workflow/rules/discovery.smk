import shlex


def min_samples_arg():
    value = LOCI.get("min_samples")
    out = "" if value is None else f"--min-samples {int(value)}"
    scope = LOCI.get("min_samples_scope", "all")
    if scope not in ("all", "group"):
        raise ValueError("config: loci.min_samples_scope must be all or group")
    if scope == "group":
        if not SAMPLE_GROUPS or not all(SAMPLE_GROUPS.values()):
            raise ValueError("config: loci.min_samples_scope group needs a `group` value for "
                             "every sample in the samples table")
        out += f" --min-samples-scope group --samples-table {shlex.quote(str(config['samples']))}"
    return out


rule discover:
    input:
        alignments=OUT / "origin/alignments.tsv.gz",
        counts=OUT / "collapsed/counts.tsv.gz",
        libraries=OUT / "collapsed/libraries.tsv",
        genomes=OUT / "reference/genomes.tsv",
        contigs=OUT / "reference/contigs.tsv",
        annotations=[g["annotation"] for g in GENOMES.values() if g.get("annotation")],
        # the group column is read with min_samples_scope group: edits must re-run discover
        samples=[config["samples"]] if LOCI.get("min_samples_scope", "all") == "group" else [],
    output:
        loci=OUT / "discovery/loci.tsv",
        reads=OUT / "discovery/locus_reads.tsv.gz",
        fasta=OUT / "discovery/candidates.fasta",
    params:
        merge_gap=LOCI.get("merge_gap", 0),
        min_rpm=LOCI.get("min_rpm", 5),
        min_rpm_genome=LOCI.get("min_rpm_genome", 0),
        min_count=LOCI.get("min_count", 1),
        min_samples=min_samples_arg(),
        peaks="" if LOCI.get("peaks", True) else "--no-peaks",
        peak_window=LOCI.get("peak_window", 2),
        peak_pvalue=LOCI.get("peak_pvalue", 0.001),
        background_flank=LOCI.get("background_flank", 50),
        satellite_distance=LOCI.get("satellite_distance", 10),
        satellite_fraction=LOCI.get("satellite_fraction", 0.1),
        n_thre=ENDS.get("n_thre", 0.5),
        penalty=ENDS.get("penalty", 0.1),
        min_five=ENDS.get("min_five_score", 0.5),
        min_three=ENDS.get("min_three_score", 0.0),
        min_isolation=ENDS.get("min_isolation", 0.5),
        five_filter=ENDS.get("five_filter", "dominant"),
        min_five_fraction=ENDS.get("min_five_fraction", 0.5),
        five_jitter=ENDS.get("five_jitter", 1),
        exclude_structural=LOCI.get("exclude_structural", "none"),
        min_exact_fraction=LOCI.get("min_exact_fraction", 0),
    log:
        OUT / "logs/discover.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} discover --alignments {input.alignments} --counts {input.counts} "
        "--libraries {input.libraries} --genomes {input.genomes} --contigs {input.contigs} "
        "--merge-gap {params.merge_gap} --min-rpm {params.min_rpm} --min-count {params.min_count} "
        "--min-rpm-genome {params.min_rpm_genome} "
        "{params.min_samples} {params.peaks} --peak-window {params.peak_window} "
        "--peak-pvalue {params.peak_pvalue} --background-flank {params.background_flank} "
        "--satellite-distance {params.satellite_distance} "
        "--satellite-fraction {params.satellite_fraction} "
        "--n-thre {params.n_thre} --penalty {params.penalty} --min-five {params.min_five} "
        "--min-three {params.min_three} --min-isolation {params.min_isolation} "
        "--five-filter {params.five_filter} --min-five-fraction {params.min_five_fraction} "
        "--five-jitter {params.five_jitter} --exclude-structural {params.exclude_structural} "
        "--min-exact-fraction {params.min_exact_fraction} "
        "--loci {output.loci} --reads {output.reads} "
        "--fasta {output.fasta} > {log} 2>&1"
