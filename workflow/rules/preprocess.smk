rule fastp:
    input:
        unpack(fastq_inputs),
    output:
        reads=OUT / "trimmed/{sample}.fastq.gz",
        json=OUT / "qc/{sample}.fastp.json",
        html=OUT / "qc/{sample}.fastp.html",
    params:
        args=fastp_args,
    log:
        OUT / "logs/fastp/{sample}.log",
    threads: min(4, THREADS)
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp {params.args} --thread {threads} --json {output.json} --html {output.html} "
        "> {log} 2>&1"


rule collapse:
    input:
        expand(OUT / "trimmed/{sample}.fastq.gz", sample=SAMPLES),
    output:
        fasta=OUT / "collapsed/unique.fa",
        counts=OUT / "collapsed/counts.tsv.gz",
        libraries=OUT / "collapsed/libraries.tsv",
    params:
        samples=lambda wc, input: " ".join(f"{s}={p}" for s, p in zip(SAMPLES, input)),
    log:
        OUT / "logs/collapse.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} collapse --sample {params.samples} --min-length {GLOBAL_MIN} "
        "--max-length {GLOBAL_MAX} --fasta {output.fasta} --counts {output.counts} "
        "--libraries {output.libraries} > {log} 2>&1"
