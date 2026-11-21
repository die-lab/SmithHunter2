#!/usr/bin/env python3
"""Download and configure one more infected vs uninfected Wolbachia small RNA study.

    python3 prepare.py --study jw18 --bench ~/bench_wolb/jw18 [--data ~/bench_wolb/dmel_data]

Same layout as benchmarks/wolbachia_walbb (lines `infected` and `tet`, each a module A run), so
../wolbachia_walbb/03_run.sh and 04_compare.py work unchanged (GENOME_FASTA=data/focus.fasta).
Genomes go to --data (shared by studies with the same host and strain, so the bowtie index can
be copied between them) and are linked as <bench>/data. Re-runnable: finished downloads are
skipped, partial ones resumed.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = "run_accession,sample_alias,sample_title,experiment_title,library_strategy,read_count,fastq_ftp"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype={fmt}&retmode=text"
GFF3 = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id={acc}&retmode=text"

DMEL = {"url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/"
               "GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz",
        "mito": "NC_024511.2"}
AAEG = {"url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/"
               "GCF_002204515.2_AaegL5.0_genomic.fna.gz",
        "mito": "NC_035159.1"}
TTRU = {"url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/476/895/GCA_028476895.1_ASM2847689v1/"
               "GCA_028476895.1_ASM2847689v1_genomic.fna.gz",   # 3 chromosomes, no mitochondrion
        "mito": "NC_024874.1", "mito_separate": True}
WMEL = {"accession": "NC_002978.6"}   # also the closest complete reference for wMelPop
WTTRU = {"accession": "NZ_CP099592.1"}   # GCF_024205405.1, one circular chromosome

# group -> (field, regex) on the ENA run table; every other run is ignored
STUDIES = {
    "jw18": {   # Rainey et al. 2016 PLoS Pathog; JW18 cells +/- wMel (tetracycline-cured), no SFV
        "accession": "PRJEB9710", "host": DMEL, "focus": WMEL,
        "groups": {"infected": ("sample_alias", r"^W0[1-5]$"), "tet": ("sample_alias", r"^T0[1-5]$")},
        "adapter": "",   # reads already trimmed (~21 nt)
    },
    "ovary": {  # He et al. 2019 BMC Genomics; ovaries +/- wMel (tetracycline-cured, >= 6 generations)
        "accession": "PRJNA439370", "host": DMEL, "focus": WMEL,
        "groups": {"infected": ("experiment_title", r"small RNA-seq of wMel-infected"),
                   "tet": ("experiment_title", r"small RNA-seq of uninfected")},
        "adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    },
    "aag2": {   # Mayoral et al. 2014 (PLoS One, PNAS: WsnRNA-46/59); Aag2 cells +/- wMelPop-CLA,
                # cytoplasm and nucleus, one library each (not replicates): min_samples = 1
        "accession": "PRJNA238976", "host": AAEG, "focus": WMEL,
        "groups": {"infected": ("sample_title", r"^Pop-"), "tet": ("sample_title", r"^Aag2-")},
        "adapter": "TGGAATTCTCGGGTGCCAAGG",
    },
    "mite": {   # Tetranychus truncatus adults (3 female + 3 male pools per line), wTtru;
                # iScience 2026 doi:10.1016/j.isci.2025.114419 (12 WsnRNAs, positive_wTtru.tsv).
                # The infected runs SRR3385484x are re-uploads of SRR3323042x (same biosamples,
                # identical first reads): only the first set is used.
        "accession": "PRJNA1253370", "host": TTRU, "focus": WTTRU,
        "groups": {"infected": ("run_accession", r"^SRR3323(0427|0428|0429|0430|0433|0434)$"),
                   "tet": ("run_accession", r"^SRR3323(0423|0424|0425|0426|0431|0432)$")},
        "adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",   # NEBNext small RNA
    },
}


def get(url: str, tries: int = 5) -> bytes:
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001 - network: retry anything
            print(f"  retry {i + 1}/{tries}: {e}", file=sys.stderr)
            time.sleep(30 * (i + 1))
    sys.exit(f"cannot fetch {url}")


def curl(url: str, out: Path) -> None:
    """Resumable download with long pauses (ENA answers 403 to fast retries)."""
    if out.exists() and out.stat().st_size > 0:
        print(f"have {out.name}")
        return
    part = out.with_name(out.name + ".part")
    for i in range(30):
        r = subprocess.run(["curl", "-fsSL", "-C", "-", "-o", str(part), url])
        if r.returncode == 0:
            part.rename(out)
            print(f"got  {out.name}")
            return
        if r.returncode == 33:   # server refuses ranges: restart
            part.unlink(missing_ok=True)
        wait = min(600, 60 * (i + 1))
        print(f"  curl {r.returncode} on {out.name}, retry in {wait}s", file=sys.stderr, flush=True)
        time.sleep(wait)
    sys.exit(f"download failed: {url}")


def fetch_checked(url: str, out: Path, prefix: bytes) -> None:
    """Download to a temporary name, keep only if it starts with prefix (error pages are HTML)."""
    if out.exists():
        return
    body = get(url)
    if not body.startswith(prefix):
        sys.exit(f"{url}: unexpected content {body[:80]!r} (not saved)")
    tmp = out.with_name(out.name + ".tmp")
    tmp.write_bytes(body)
    tmp.rename(out)


def genomes(spec: dict, data: Path) -> None:
    data.mkdir(parents=True, exist_ok=True)
    focus = spec["focus"]["accession"]
    fetch_checked(EUTILS.format(acc=focus, fmt="fasta"), data / "focus.fasta", b">")
    fetch_checked(GFF3.format(acc=focus), data / "focus.gff3", b"##gff-version")
    host = data / "host.fna.gz"
    curl(spec["host"]["url"], host)
    mt, nuc = data / "host_mt.fasta", data / "host_nuc.fasta"
    if spec["host"].get("mito_separate"):   # assembly without mitochondrion: fetch it apart
        fetch_checked(EUTILS.format(acc=spec["host"]["mito"], fmt="fasta"), mt, b">")
        if not nuc.exists():
            tmp = nuc.with_name(nuc.name + ".tmp")
            with gzip.open(host, "rt") as fh, open(tmp, "w") as fn:
                for line in fh:
                    fn.write(line)
            tmp.rename(nuc)
    elif not (mt.exists() and nuc.exists()):
        print("split host into nuclear and mitochondrial")
        mito = ">" + spec["host"]["mito"]
        with gzip.open(host, "rt") as fh, open(mt, "w") as fm, open(nuc, "w") as fn:
            out = fn
            for line in fh:
                if line.startswith(">"):
                    out = fm if line.split()[0] == mito else fn
                out.write(line)
        if mt.stat().st_size == 0:
            sys.exit(f"mito contig {mito} not found")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", required=True, choices=sorted(STUDIES))
    ap.add_argument("--bench", required=True, type=Path)
    ap.add_argument("--data", type=Path, help="shared genome folder (default <bench>/data)")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    spec = STUDIES[args.study]
    bench = args.bench.expanduser().resolve()
    data = (args.data or bench / "data").expanduser().resolve()
    (bench / "fastq").mkdir(parents=True, exist_ok=True)

    table = get(f"{ENA}?accession={spec['accession']}&result=read_run&fields={FIELDS}&format=tsv").decode()
    (bench / "runs.tsv").write_text(table)
    runs = list(csv.DictReader(io.StringIO(table), delimiter="\t"))
    chosen: dict[str, list[dict]] = {g: [] for g in spec["groups"]}
    for r in runs:
        hits = [g for g, (field, rx) in spec["groups"].items() if re.search(rx, r.get(field, ""))]
        if len(hits) > 1:
            sys.exit(f"{r['run_accession']} matches several groups: {hits}")
        if hits:
            chosen[hits[0]].append(r)
    with open(bench / "groups.tsv", "w") as fh:
        fh.write("run\tgroup\tsample_alias\tdescription\n")
        for g, rs in chosen.items():
            for r in sorted(rs, key=lambda r: r["run_accession"]):
                fh.write(f"{r['run_accession']}\t{g}\t{r['sample_alias']}\t{r['experiment_title']}\n")
                print(f"  {r['run_accession']}\t{g}\t{r['sample_alias']}\t{r['read_count']}")
    if any(not rs for rs in chosen.values()):
        sys.exit("a group has no runs: check STUDIES")
    no_fastq = [r["run_accession"] for rs in chosen.values() for r in rs if not r["fastq_ftp"].strip()]
    if no_fastq:
        sys.exit(f"no ENA fastq for: {', '.join(no_fastq)}")

    for rs in chosen.values():
        for r in rs:
            url = r["fastq_ftp"].split(";")[0]
            curl(f"https://{url}", bench / "fastq" / f"{r['run_accession']}.fastq.gz")
    genomes(spec, data)
    if not (bench / "data").exists():
        (bench / "data").symlink_to(data)

    template = (HERE.parent / "wolbachia_walbb/config_template.yaml").read_text()
    for old, new in (("data/wAlbB.fasta", "data/focus.fasta"), ("data/wAlbB.gff3", "data/focus.gff3"),
                     ("data/AaegL5_mt.fasta", "data/host_mt.fasta"),
                     ("data/AaegL5_nuc.fasta", "data/host_nuc.fasta")):
        if old not in template:
            sys.exit(f"config template changed: {old} not found")
        template = template.replace(old, new)
    for group, rs in chosen.items():
        d = bench / group
        (d / "config").mkdir(parents=True, exist_ok=True)
        for name in ("data", "fastq"):
            if not (d / name).exists():
                (d / name).symlink_to(Path("..") / name)
        text = (template.replace("@GROUP@", group).replace("@THREADS@", str(args.threads))
                .replace("@ADAPTER@", spec["adapter"]))
        (d / "config/config.yaml").write_text(text)
        with open(d / "config/samples.tsv", "w") as fh:
            fh.write("sample\tfq1\tfq2\n")
            for r in sorted(rs, key=lambda r: r["run_accession"]):
                alias = r["sample_alias"]
                name = (f"{alias}_{r['run_accession']}" if re.fullmatch(r"[A-Za-z0-9_.-]{1,20}", alias)
                        else r["run_accession"])
                fh.write(f"{name}\t{bench / 'fastq' / (r['run_accession'] + '.fastq.gz')}\t\n")
        print(f"wrote {d}/config")
    print("next: BENCH=<bench> GENOME_FASTA=data/focus.fasta CORES=8 bash ../wolbachia_walbb/03_run.sh")


if __name__ == "__main__":
    main()
