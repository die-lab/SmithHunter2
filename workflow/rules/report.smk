# Run report. The read funnel and provenance.yaml are plain Python and part of rule `all`;
# plots (R/ggplot2) and MultiQC need their own conda envs and run with rule `report`:
#   snakemake -s workflow/Snakefile --cores 2 --sdm conda report
# Adapted from the reporting layer of the smithhunter24 port of v0 module A.

REPORT = OUT / "report"
RSCRIPTS = Path(workflow.basedir) / "scripts"
ENV_FILES = sorted(str(p) for p in (Path(workflow.basedir) / "envs").glob("*.yaml"))


def git_state():
    import subprocess
    repo = str(Path(workflow.basedir).parent)
    try:
        commit = subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"], capture_output=True,
                                text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", repo, "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return commit + ("-dirty" if dirty else "")


rule read_funnel:
    input:
        fastp=expand(OUT / "qc/{sample}.fastp.json", sample=SAMPLES),
        libraries=OUT / "collapsed/libraries.tsv",
        summary=OUT / "origin/origin_summary.tsv",
        loci=OUT / "discovery/loci.tsv",
    output:
        funnel=REPORT / "read_funnel.tsv",
        mqc_funnel=REPORT / "multiqc_input/read_funnel_mqc.tsv",
        mqc_origin=REPORT / "multiqc_input/read_origin_mqc.tsv",
    log:
        OUT / "logs/report/read_funnel.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} funnel --fastp {input.fastp} --libraries {input.libraries} "
        "--summary {input.summary} --loci {input.loci} --out {output.funnel} "
        "--mqc-funnel {output.mqc_funnel} --mqc-origin {output.mqc_origin} > {log} 2>&1"


# Runs inside the Snakemake process (no conda env): it records the config and versions of
# this run, and the code it needs is stdlib only.
rule provenance:
    input:
        loci=OUT / "discovery/loci.tsv",
        fastp=expand(OUT / "qc/{sample}.fastp.json", sample=SAMPLES),
        envs=ENV_FILES,
    output:
        REPORT / "provenance.yaml",
    run:
        import datetime
        from snakemake import __version__ as snakemake_version
        sys.path.insert(0, str(SRC))
        from smithhunter import __version__ as smithhunter_version
        from smithhunter.report import provenance, write_yaml
        software = {"smithhunter": smithhunter_version, "git_commit": git_state(),
                    "snakemake": snakemake_version, "python": sys.version.split()[0]}
        data = provenance(dict(config), list(SAMPLES), list(input.envs), list(input.fastp),
                          input.loci, software)
        stamp = datetime.datetime.now().isoformat(timespec="seconds")
        write_yaml(output[0], data, header=f"SmithHunter2 run provenance, written {stamp}")


rule coverage_tracks:
    input:
        alignments=OUT / "origin/alignments.tsv.gz",
        counts=OUT / "collapsed/counts.tsv.gz",
        genomes=OUT / "reference/genomes.tsv",
        contigs=OUT / "reference/contigs.tsv",
    output:
        REPORT / "coverage.tsv.gz",
    log:
        OUT / "logs/report/coverage.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} coverage --alignments {input.alignments} --counts {input.counts} "
        "--genomes {input.genomes} --contigs {input.contigs} --out {output} > {log} 2>&1"


rule plot_reads:
    input:
        funnel=REPORT / "read_funnel.tsv",
        summary=OUT / "origin/origin_summary.tsv",
        lengths=OUT / "origin/length_distribution.tsv",
        script=RSCRIPTS / "plot_reads.R",
    output:
        REPORT / "reads.pdf",
    log:
        OUT / "logs/report/plot_reads.log",
    conda:
        "../envs/r.yaml"
    shell:
        "Rscript {input.script} {input.funnel} {input.summary} {input.lengths} {output} "
        "> {log} 2>&1"


rule plot_genomes:
    input:
        coverage=REPORT / "coverage.tsv.gz",
        loci=OUT / "discovery/loci.tsv",
        script=RSCRIPTS / "plot_genomes.R",
    output:
        REPORT / "genomes.pdf",
    log:
        OUT / "logs/report/plot_genomes.log",
    conda:
        "../envs/r.yaml"
    shell:
        "Rscript {input.script} {input.coverage} {input.loci} {output} > {log} 2>&1"


rule plot_loci:
    input:
        reads=OUT / "discovery/locus_reads.tsv.gz",
        loci=OUT / "discovery/loci.tsv",
        script=RSCRIPTS / "plot_loci.R",
    output:
        REPORT / "loci.pdf",
    log:
        OUT / "logs/report/plot_loci.log",
    conda:
        "../envs/r.yaml"
    shell:
        "Rscript {input.script} {input.reads} {input.loci} {output} > {log} 2>&1"


rule multiqc:
    input:
        fastp=expand(OUT / "qc/{sample}.fastp.json", sample=SAMPLES),
        align_log=OUT / "logs/align.log",
        mqc=[REPORT / "multiqc_input/read_funnel_mqc.tsv",
             REPORT / "multiqc_input/read_origin_mqc.tsv"],
    output:
        REPORT / "multiqc_report.html",
    log:
        OUT / "logs/report/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force --no-data-dir -o {REPORT} -n multiqc_report.html "
        "{input.fastp} {input.align_log} {input.mqc} > {log} 2>&1"


rule report:
    input:
        REPORT / "read_funnel.tsv",
        REPORT / "provenance.yaml",
        REPORT / "reads.pdf",
        REPORT / "genomes.pdf",
        REPORT / "loci.pdf",
        REPORT / "multiqc_report.html",
