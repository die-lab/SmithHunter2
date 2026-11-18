# Precursor folding (module B). Run separately, like `targets`:
#   snakemake -s workflow/Snakefile --cores 8 precursors
PRECURSORS = config.get("precursors") or {}


def precursor_args(wildcards):
    # Checked here, not when the workflow is parsed, so that module A and `targets`
    # run with no `precursors:` block in the config.
    if PRECURSORS.get("temperature") is None:
        raise ValueError(
            "config: precursors.temperature is required for the `precursors` rule and has no "
            "default (folding temperature in Celsius, e.g. 25 as in SmithHunter v0 for "
            "bivalves, 37 for mammals). Add to the config:  precursors: {temperature: 25}")
    flanks = PRECURSORS.get("flanks", [15, 50])
    if not isinstance(flanks, list) or len(flanks) != 2:
        raise ValueError("config: precursors.flanks must be [short, long], e.g. [15, 50]")
    return (f"--temperature {float(PRECURSORS['temperature'])} "
            f"--flanks {int(flanks[0])} {int(flanks[1])} "
            f"--shuffles {int(PRECURSORS.get('shuffles', 1000))}")


rule fold_precursors:
    input:
        loci=OUT / "discovery/loci.tsv",
        candidates=CANDIDATES,
        genomes=OUT / "reference/genomes.tsv",
        contigs=OUT / "reference/contigs.tsv",
        fasta=[g["fasta"] for g in GENOMES.values() if g.get("role", "focus") == "focus"],
    output:
        OUT / "precursors/precursors.tsv",
    params:
        args=precursor_args,
    threads: THREADS
    log:
        OUT / "logs/precursors/fold.log",
    conda:
        "../envs/rna.yaml"
    shell:
        "{SH} precursors --loci {input.loci} --candidates {input.candidates} "
        "--genomes {input.genomes} --contigs {input.contigs} {params.args} "
        "--threads {threads} --out {output} > {log} 2>&1"


rule precursors:
    input:
        OUT / "precursors/precursors.tsv",
