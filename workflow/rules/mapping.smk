rule build_reference:
    input:
        [g["fasta"] for g in GENOMES.values()],
    output:
        fasta=OUT / "reference/combined.fa",
        genomes=OUT / "reference/genomes.tsv",
        contigs=OUT / "reference/contigs.tsv",
    params:
        spec=json.dumps(GENOMES),
    log:
        OUT / "logs/build_reference.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} build-ref --spec {params.spec:q} --min-length {MIN_LEN} --max-length {MAX_LEN} "
        "--fasta {output.fasta} --genomes {output.genomes} --contigs {output.contigs} > {log} 2>&1"


rule bowtie_index:
    input:
        OUT / "reference/combined.fa",
    output:
        directory(OUT / "reference/index"),
    log:
        OUT / "logs/bowtie_index.log",
    threads: THREADS
    conda:
        "../envs/bowtie.yaml"
    shell:
        "mkdir -p {output} && bowtie-build --threads {threads} {input} {output}/combined > {log} 2>&1"


# --best --strata -k K -m K: report every equally best alignment, and set aside reads
# with more than K of them. --reorder keeps the alignments of a read together.
rule align:
    input:
        reads=OUT / "collapsed/unique.fa",
        index=OUT / "reference/index",
    output:
        sam=OUT / "mapping/hits.sam.gz",
        excess=OUT / "mapping/too_many_hits.fa",
    params:
        mismatches=MAPPING.get("mismatches", 1),
        max_hits=MAPPING.get("max_hits", 50),
    log:
        OUT / "logs/align.log",
    threads: THREADS
    conda:
        "../envs/bowtie.yaml"
    shell:
        "bowtie -f -v {params.mismatches} -k {params.max_hits} -m {params.max_hits} "
        "--best --strata --reorder -p {threads} --sam --max {output.excess} "
        "{input.index}/combined {input.reads} 2> {log} | gzip -c > {output.sam}; "
        "touch {output.excess}"


rule assign_origin:
    input:
        sam=OUT / "mapping/hits.sam.gz",
        excess=OUT / "mapping/too_many_hits.fa",
        counts=OUT / "collapsed/counts.tsv.gz",
        genomes=OUT / "reference/genomes.tsv",
        contigs=OUT / "reference/contigs.tsv",
    output:
        alignments=OUT / "origin/alignments.tsv.gz",
        origin=OUT / "origin/read_origin.tsv.gz",
        summary=OUT / "origin/origin_summary.tsv",
        lengths=OUT / "origin/length_distribution.tsv",
    params:
        ambiguous=config.get("ambiguous", "report"),
        multimap=LOCI.get("multimap", "fractional"),
    log:
        OUT / "logs/assign_origin.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} assign --sam {input.sam} --excess {input.excess} --counts {input.counts} "
        "--genomes {input.genomes} --contigs {input.contigs} --ambiguous {params.ambiguous} "
        "--multimap {params.multimap} --alignments {output.alignments} --origin {output.origin} "
        "--summary {output.summary} --lengths {output.lengths} > {log} 2>&1"
