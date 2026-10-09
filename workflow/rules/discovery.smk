def min_samples_arg():
    value = LOCI.get("min_samples")
    return "" if value is None else f"--min-samples {int(value)}"


rule discover:
    input:
        alignments=OUT / "origin/alignments.tsv.gz",
        counts=OUT / "collapsed/counts.tsv.gz",
        libraries=OUT / "collapsed/libraries.tsv",
        genomes=OUT / "reference/genomes.tsv",
        contigs=OUT / "reference/contigs.tsv",
        annotations=[g["annotation"] for g in GENOMES.values() if g.get("annotation")],
    output:
        loci=OUT / "discovery/loci.tsv",
        reads=OUT / "discovery/locus_reads.tsv.gz",
        fasta=OUT / "discovery/candidates.fasta",
    params:
        merge_gap=LOCI.get("merge_gap", 0),
        min_rpm=LOCI.get("min_rpm", 5),
        min_samples=min_samples_arg(),
        n_thre=ENDS.get("n_thre", 0.5),
        penalty=ENDS.get("penalty", 0.1),
        min_five=ENDS.get("min_five_score", 0.5),
        min_three=ENDS.get("min_three_score", 0.0),
    log:
        OUT / "logs/discover.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} discover --alignments {input.alignments} --counts {input.counts} "
        "--libraries {input.libraries} --genomes {input.genomes} --contigs {input.contigs} "
        "--merge-gap {params.merge_gap} --min-rpm {params.min_rpm} {params.min_samples} "
        "--n-thre {params.n_thre} --penalty {params.penalty} --min-five {params.min_five} "
        "--min-three {params.min_three} --loci {output.loci} --reads {output.reads} "
        "--fasta {output.fasta} > {log} 2>&1"
