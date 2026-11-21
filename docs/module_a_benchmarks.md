# Module A on a bacterial symbiont: real-data benchmarks and decisions

Started 2026-11-21 (branch `wolbachia-benchmark`). Question: on the maintainer's *Wolbachia* data
module A returned 300+ loci, most passing, where v0 found 105 clusters. Which part of the
pipeline is responsible, and which settings predict real symbiont small RNAs best?

## 1. Data

Every study compares an infected line with an uninfected one (tetracycline-cured, or the
naturally uninfected parental line). A *Wolbachia* locus passing in the uninfected line is a
false positive. Known *Wolbachia* small RNAs are the positive control where the strain has them.

| bench | study | host / strain | infected + uninfected | Wolbachia share of the infected reads | known small RNAs |
|---|---|---|---|---|---|
| `walbb` | PRJNA671731, Sci Rep 2022 | *Ae. aegypti* adults / wAlbB | 9 + 9 | 1.9% | none (WsnRNA-46 has no homolog in wAlbB: best 3 mm / 27 nt, 0 reads) |
| `aag2` | PRJNA238976, Mayoral 2014 | Aag2 cells / wMelPop (wMel reference) | 2 + 2 (cytoplasm, nucleus) | 15% | WsnRNA-46, WsnRNA-59 (Rfam RF02625 / RF02626 precursors) |
| `jw18` | PRJEB9710, Rainey 2016 | *Drosophila* JW18 cells / wMel | 5 + 5 | 20% | WsnRNA-46, WsnRNA-59 |
| `ovary` | PRJNA439370, He 2019 | *D. melanogaster* ovaries / wMel | 3 + 3 | 0.03% | WsnRNA-46, WsnRNA-59 |
| `mite` | PRJNA1253370, iScience 2026 | *Tetranychus truncatus* adults / wTtru | 6 + 6 (3 female + 3 male) | (running) | 12 published WsnRNAs (`positive_wTtru.tsv`) |

Scripts: `benchmarks/wolbachia_walbb/` (01-04 for wAlbB, 05 variants, 06 summary, 07 concordance),
`benchmarks/wolbachia_multi/prepare.py` (other studies). Search for datasets (2026-11-21): no
negative-control small RNA data were found for *Ae. albopictus*, *Culex* wPip or Aa23/RML12 cells;
other candidates: *Laodelphax* PRJNA1185640 (3 + 3), *Ostrinia* PRJDB14352 (cured and naturally
uninfected, no replicates), *Brugia* wBm (no uninfected control).

## 2. Findings (repository defaults unless stated)

1. **No false positives in any study** (0 passing loci in every uninfected line, every variant
   tried). The genome assignment keeps host reads off the symbiont; the "300 vs 105" gap is not
   an assignment problem.
2. **Structural RNA fragments are a large share of the candidates**: wAlbB 40 of 157 loci (58% of
   their reads), Aag2 67 of 570 (88%, one 16S rRNA fragment with 1.7M reads), JW18 104 of 1957.
3. **Classification bug**: tmRNA, RNase P, SRP RNAs were classified as `exon` (their child feature)
   because the NCBI prokaryotic types were missing from `CLASS_PRIORITY`. Fixed.
4. **Adapter check bug** in `02_samples.py` (only TruSeq small RNA probed; PRJNA671731 carries the
   Illumina universal adapter). fastp autodetection had found the right one (98.7% trimmed).
5. **Positive controls**: WsnRNA-59 passes in Aag2 (270 reads) and JW18 (2,001); WsnRNA-46 passes
   in JW18 (2,332) but not in Aag2 (~170 reads split in several 5' peaks of the hairpin arm:
   fails expression and isolation 0.45 < 0.5; recovered only with `min_rpm` 2 and
   `min_isolation` 0.3, at the cost of 1,590 candidates). In ovaries neither has more than 5 reads.
6. **The library-RPM threshold follows the symbiont load**: JW18 (20% of the reads) passes 1,957
   loci, wAlbB (1.9%) 157, ovaries (0.03%) 0. A per-library threshold is not a per-genome one.
7. Stricter `min_rpm` cuts candidates fast (wAlbB 125 / 61 / 16 at 5 / 10 / 20 with `sense`;
   Aag2 504 / 259 / 132) but **20 loses WsnRNA-59 in Aag2**; `min_count` has no effect at these
   depths (5 RPM is already 30-100 reads).

## 3. Decisions

| decision | evidence | status |
|---|---|---|
| `loci.exclude_structural: none / sense / any`, column `structural` always written; default `none` (results unchanged, simulated tRFs are targets), `sense` advised for symbionts | finding 2: `sense` removes 32-75 loci per study, loses no positive control and adds no false positive | done (755a6c2, review fixes 2dd94d4) |
| structural detection by type, gbkey, ncrna_class; products only on RNA features | code review: eukaryotic mRNA exons carry protein products ("signal recognition particle 54 kDa protein") | done (2dd94d4) |
| `loci.min_rpm_genome` (per million reads of the locus' genome, on top of `min_rpm`), default 0 | finding 6 | implemented (926c570); value under evaluation (section 4) |

## 4. Open: expression threshold

Grid on wAlbB, Aag2, JW18 (`variants_grid2`): `sense` x {`min_rpm` 10, `min_rpm_genome` 25, 50,
100, 200}; criteria: positives kept, false positives 0, candidates, cross-study concordance of
passing loci between Aag2 and JW18 (same wMel reference; `07_concordance.py`).
