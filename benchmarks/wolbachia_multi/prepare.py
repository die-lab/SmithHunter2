#!/usr/bin/env python3
"""Download and configure one more infected vs uninfected symbiont small RNA study
(Wolbachia; Buchnera with a symbiont-free tissue as the negative).

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
MPER = {"url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/856/785/"
               "GCF_001856785.1_MPER_G0061.0/GCF_001856785.1_MPER_G0061.0_genomic.fna.gz",
        "mito": "NC_029727.1"}   # Myzus persicae clone G006 (one of the Feng 2018 lines)
WMEL = {"accession": "NC_002978.6"}   # also the closest complete reference for wMelPop
WTTRU = {"accession": "NZ_CP099592.1"}   # GCF_024205405.1, one circular chromosome
# Buchnera of M. persicae G002: chromosome NZ_CP002701.1 + plasmid pLeu NZ_CP002702.1
BUCH_G002 = {"assembly": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/521/565/"
                         "GCF_000521565.1_ASM52156v1/GCF_000521565.1_ASM52156v1"}
HSAP = {"url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/"
               "GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
        "mito": "NC_012920.1"}
BTAB_Q = {"url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/994/315/GCA_003994315.1_CAAS_Bemtab1.0/"
                 "GCA_003994315.1_CAAS_Bemtab1.0_genomic.fna.gz",   # isolate CAASIVF2012BtQ, no mito
          "mito": "KU663131.1", "mito_separate": True}
CTR_L2B = {"accession": "NC_010280.1,NC_010285.1"}   # C. trachomatis L2b/UCH-1 chromosome + pL2
MTB_H37RV = {"accession": "NC_000962.3"}
MMUS = {"url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/"
               "GCF_000001635.27_GRCm39_genomic.fna.gz",
        "mito": "NC_005089.1"}
LPN_PHILA = {"accession": "NC_002942.5"}   # L. pneumophila Philadelphia 1
# Cardinium of B. tabaci (China, CAAS): 3 scaffolds (993,648 / 2,448 / 713 bp): linear
CARD_BTAB = {"assembly": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/300/865/"
                         "GCF_004300865.1_CAAS_CanCar1.0/GCF_004300865.1_CAAS_CanCar1.0", "linear": True}
PORTIERA_MED = {"accession": "NZ_CP016304.1"}   # obligate symbiont of B. tabaci MED, in both lines
TRUSEQ_SMALL = "TGGAATTCTCGGGTGCCAAGG"
ILLUMINA_UNIVERSAL = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

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
        "sample_groups": ("sample_alias", {"male": r"M\d$", "female": r"F\d$"}),
    },
    "buchnera_myzus": {   # Feng et al. 2018 Mol Ecol (PMID 29271121): M. persicae lines G002,
                # G006, BTIRed; bacteriocytes (Buchnera) vs gut (no Buchnera) of the same lines.
                # `tet` holds the gut libraries: the negative control is a tissue, not a cured line.
                # G006 runs are paired 2 x 101 with the mates swapped between runs: the mate that
                # carries the TruSeq small RNA 3' adapter is the sense strand and the only one used.
                # G002 / BTIRed: single 52 nt, Illumina universal adapter. Both go to fastp.
        "accession": "PRJNA395356", "host": MPER, "focus": BUCH_G002,
        "groups": {"infected": ("experiment_title", r"Myzus persicae \S+: Bacteriocyte$"),
                   "tet": ("experiment_title", r"Myzus persicae \S+: Gut$")},
        "adapter": "", "adapters": [TRUSEQ_SMALL, ILLUMINA_UNIVERSAL], "sense_adapter": TRUSEQ_SMALL[:12],
        "focus_name": "buchnera",   # FOCUS=buchnera for 03_run.sh / --genome for 05_variants.py
    },
    "chlamydia": {   # Furuse et al. 2014 PLoS One (PMID 25184567): small RNA-seq of infected cells,
                # total and RISC (Ago-IP) libraries, reads already trimmed. `infected`: C. trachomatis
                # in HeLa; `tet`: M. tuberculosis in THP-1 (human cells without Chlamydia: the
                # negative, and reads of another bacterium, not in the reference, as off-target test).
                # Known sRNAs: positive_ctr.tsv (Albrecht et al. 2010, L2b coordinates).
        "accession": "PRJNA248348", "host": HSAP, "focus": CTR_L2B,
        "groups": {"infected": ("experiment_title", r"Chlamydia trachomatis in HeLa"),
                   "tet": ("experiment_title", r"Mycobacterium tuberculosis in THP-1")},
        "adapter": "", "focus_name": "chlamydia",
    },
    "mtb": {   # Furuse et al. 2014 (as chlamydia, roles swapped): `infected` M. tuberculosis in
                # THP-1 (total + RISC), `tet` C. trachomatis in HeLa (human cells without Mtb).
                # Known sRNAs: the 19 RefSeq ncRNAs of H37Rv (positive_mtb.tsv).
        "accession": "PRJNA248348", "host": HSAP, "focus": MTB_H37RV,
        "groups": {"infected": ("experiment_title", r"Mycobacterium tuberculosis in THP-1"),
                   "tet": ("experiment_title", r"Chlamydia trachomatis in HeLa")},
        "adapter": "", "focus_name": "mtb",
    },
    "legionella": {   # Furuse et al. 2014: `infected` L. pneumophila in mouse RAW264.7 macrophages
                # (total + RISC), `tet` RAW264.7 with M. marinum or M. smegmatis (no Legionella).
        "accession": "PRJNA248348", "host": MMUS, "focus": LPN_PHILA,
        "groups": {"infected": ("experiment_title", r"Legionella pneumophila in RAW"),
                   "tet": ("experiment_title", r"Mycobacterium (marinum|smegmatis) in RAW")},
        "adapter": "", "focus_name": "legionella",
    },
    "cardinium": {   # Li et al. 2018 Front Physiol (PMID 30483149): B. tabaci Q, Cardinium-infected
                # vs uninfected of the same genetic background, 2 + 2 libraries. No known Cardinium
                # small RNAs: negative control only. Portiera (in both lines) is an exclude genome.
        "accession": "PRJNA321074", "host": BTAB_Q, "focus": CARD_BTAB,
        "groups": {"infected": ("sample_alias", r"^CI\d$"), "tet": ("sample_alias", r"^CU\d$")},
        "adapter": ILLUMINA_UNIVERSAL, "focus_name": "cardinium",
        "extra_genomes": {"portiera": PORTIERA_MED},
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


def gunzip_checked(src: Path, out: Path, prefix: str) -> None:
    if out.exists():
        return
    tmp = out.with_name(out.name + ".tmp")
    with gzip.open(src, "rt") as fh, open(tmp, "w") as fo:
        first = fh.readline()
        if not first.startswith(prefix):
            sys.exit(f"{src}: unexpected content {first[:80]!r}")
        fo.write(first)
        for line in fh:
            fo.write(line)
    tmp.rename(out)


def bacterium(spec: dict, data: Path, name: str) -> None:
    """data/<name>.fasta and .gff3 from an assembly (all its sequences) or accession(s)."""
    if "assembly" in spec:   # whole assembly: chromosome and plasmids / scaffolds
        base = spec["assembly"]
        curl(base + "_genomic.fna.gz", data / f"{name}.fna.gz")
        curl(base + "_genomic.gff.gz", data / f"{name}.gff.gz")
        gunzip_checked(data / f"{name}.fna.gz", data / f"{name}.fasta", ">")
        gunzip_checked(data / f"{name}.gff.gz", data / f"{name}.gff3", "##gff-version")
    else:   # one accession or several, comma separated
        acc = spec["accession"]
        fetch_checked(EUTILS.format(acc=acc, fmt="fasta"), data / f"{name}.fasta", b">")
        fetch_checked(GFF3.format(acc=acc), data / f"{name}.gff3", b"##gff-version")
        want = set(acc.split(","))
        got = {line[1:].split()[0] for line in open(data / f"{name}.fasta") if line.startswith(">")}
        if got != want:
            sys.exit(f"{name}.fasta: sequences {sorted(got)}, expected {sorted(want)}")


def genomes(spec: dict, data: Path) -> None:
    data.mkdir(parents=True, exist_ok=True)
    bacterium(spec["focus"], data, "focus")
    for name, extra in spec.get("extra_genomes", {}).items():
        bacterium(extra, data, name)
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


def head_reads(url: str, n: int = 2000) -> list[str]:
    """First n read sequences of a gzipped fastq (streamed, no file written)."""
    req = urllib.request.Request(f"https://{url}")
    with urllib.request.urlopen(req, timeout=120) as r, gzip.GzipFile(fileobj=r) as gz:
        reads = []
        for i, line in enumerate(io.TextIOWrapper(gz)):
            if i % 4 == 1:
                reads.append(line.strip())
                if len(reads) >= n:
                    break
    return reads


def sense_mate(urls: list[str], adapter: str) -> str:
    """Single-end runs: the only file. Paired small RNA runs: the mate whose reads carry the
    3' adapter (the sense strand of the small RNA), the other one reads the 5' adapter."""
    singles = [u for u in urls if not re.search(r"_[12]\.fastq\.gz$", u)]
    mates = [u for u in urls if re.search(r"_[12]\.fastq\.gz$", u)]
    if not mates:
        return urls[0]
    if not adapter:
        return singles[0] if singles else mates[0]
    share = {}
    for u in mates:
        reads = head_reads(u)
        share[u] = sum(adapter in s for s in reads) / max(1, len(reads))
    best = max(share, key=share.get)
    print(f"  {Path(best).name}: sense mate ({', '.join(f'{Path(u).name} {v:.0%}' for u, v in share.items())} "
          f"with {adapter})")
    if share[best] < 0.5 or sorted(share.values())[0] > 0.1:
        sys.exit(f"cannot tell the sense mate of {urls}: {share}")
    return best


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
            urls = [u for u in r["fastq_ftp"].split(";") if u]
            url = sense_mate(urls, spec.get("sense_adapter", ""))
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
        if spec.get("focus_name"):
            if "\n  wolbachia:\n" not in text:
                sys.exit("config template changed: genome wolbachia not found")
            text = text.replace("\n  wolbachia:\n", f"\n  {spec['focus_name']}:\n")
        if spec["focus"].get("linear"):   # scaffolds: no wrap around their ends
            text = text.replace("    topology: circular", "    topology: linear", 1)
        for name, extra in spec.get("extra_genomes", {}).items():
            if "\n  host_mt:\n" not in text:
                sys.exit("config template changed: host_mt not found")
            block = (f"\n  {name}:\n    fasta: data/{name}.fasta\n    role: exclude\n"
                     f"    topology: {'linear' if extra.get('linear') else 'circular'}\n")
            text = text.replace("\n  host_mt:\n", block + "  host_mt:\n", 1)
        if spec.get("adapters"):   # libraries with different adapters in one line
            fa = bench / "adapters.fa"
            fa.write_text("".join(f">adapter{i + 1}\n{a}\n" for i, a in enumerate(spec["adapters"])))
            if '  extra: ""' not in text:
                sys.exit("config template changed: trimming extra not found")
            text = text.replace('  extra: ""', f'  extra: "--adapter_fasta {fa}"')
        (d / "config/config.yaml").write_text(text)
        field, patterns = spec.get("sample_groups", ("", {}))
        with open(d / "config/samples.tsv", "w") as fh:
            fh.write("sample\tfq1\tfq2\tgroup\n")   # group: for loci.min_samples_scope group
            for r in sorted(rs, key=lambda r: r["run_accession"]):
                alias = r["sample_alias"]
                name = (f"{alias}_{r['run_accession']}" if re.fullmatch(r"[A-Za-z0-9_.-]{1,20}", alias)
                        else r["run_accession"])
                grp = next((g for g, rx in patterns.items() if re.search(rx, r.get(field, ""))), "")
                fh.write(f"{name}\t{bench / 'fastq' / (r['run_accession'] + '.fastq.gz')}\t\t{grp}\n")
        print(f"wrote {d}/config")
    focus = f"FOCUS={spec['focus_name']} " if spec.get("focus_name") else ""
    print(f"next: BENCH=<bench> {focus}GENOME_FASTA=data/focus.fasta CORES=8 bash ../wolbachia_walbb/03_run.sh")


if __name__ == "__main__":
    main()
