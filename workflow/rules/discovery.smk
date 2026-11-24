import shlex

sys.path.insert(0, str(SRC))
from smithhunter.settings import discover_args, needs_groups  # noqa: E402  (see settings.py)


def discover_settings():
    if needs_groups(config):
        if not SAMPLE_GROUPS or not all(SAMPLE_GROUPS.values()):
            raise ValueError("config: loci.min_samples_scope group and loci.control_group need "
                             "a `group` value for every sample in the samples table")
    return shlex.join(discover_args(config, config["samples"]))


rule discover:
    input:
        alignments=OUT / "origin/alignments.tsv.gz",
        counts=OUT / "collapsed/counts.tsv.gz",
        libraries=OUT / "collapsed/libraries.tsv",
        genomes=OUT / "reference/genomes.tsv",
        contigs=OUT / "reference/contigs.tsv",
        annotations=[g["annotation"] for g in GENOMES.values() if g.get("annotation")],
        # the group column is read (min_samples_scope group, control_group): edits re-run discover
        samples=[config["samples"]] if needs_groups(config) else [],
    output:
        loci=OUT / "discovery/loci.tsv",
        reads=OUT / "discovery/locus_reads.tsv.gz",
        fasta=OUT / "discovery/candidates.fasta",
    params:
        settings=discover_settings(),
    log:
        OUT / "logs/discover.log",
    conda:
        "../envs/python.yaml"
    shell:
        "{SH} discover --alignments {input.alignments} --counts {input.counts} "
        "--libraries {input.libraries} --genomes {input.genomes} --contigs {input.contigs} "
        "{params.settings} "
        "--loci {output.loci} --reads {output.reads} "
        "--fasta {output.fasta} > {log} 2>&1"
