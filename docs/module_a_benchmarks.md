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
| `loci.min_rpm_genome` (per million reads of the locus' genome, on top of `min_rpm`), default 0 | findings 6, 9 | done (926c570); 100 in the profiles |
| `loci.min_exact_fraction`, column `exact_fraction`, default 0 | finding 11 | done (104ae10); 0.5 in the profiles (Aag2 with the wMel reference for wMelPop loses 6 of 504 loci, keeps WsnRNA-59) |
| defaults unchanged; two profiles `config/profiles/symbiont_{precise,sensitive}.yaml` added after the main config | section 5 | done |

## 4. Expression threshold and off-target reads (grids 2-5)

8. **Spider mite (12 published WsnRNAs, ~9 distinct loci; Wolbachia 0.35% of the reads)**: at the
   defaults 22 candidates, 1 known WsnRNA passing (7542, 3,468 reads). Every other known WsnRNA
   exists as a locus whose representative read is exactly the published mature sequence, but fails
   the expression filter (5 RPM of the library = 50-62 reads per sample). WsnRNA-744 is
   male-biased (149-181 reads in males, 22-54 in females): with mixed-sex pools `min_samples`
   n - 1 rejects it.
9. **Library RPM is wrong both ways for a symbiont**: lenient at high load (JW18 20%: 1,957
   candidates), strict at low load (mite 0.35%: 22; ovary 0.03%: 0). `min_rpm_genome` follows the
   load: with `sense`, `min_rpm_genome` 100 gives JW18 706 (both positives kept), Aag2 177 (59 kept),
   wAlbB unchanged (125, low load); 200 loses the wMel positives (they sit at ~120-200 reads per
   million Wolbachia reads in one library).
10. **The library RPM is what guards the uninfected line**: lowering `min_rpm` to 1 for the mite's
    sake gives 2-3 passing loci in the JW18 cured line. They are the most abundant Wolbachia small
    RNAs of the infected line (320,706 reads for the top one) at 0.03% vs 20% of the library: cross-talk
    between multiplexed libraries or residual infection, not misassigned host reads.
11. **Off-target organisms**: the mite cured line has 2 passing loci at the defaults, both far more
    abundant in the cured than in the infected line (a tRNA: 7,371 vs 26 reads; dnaA antisense:
    653 vs 1), i.e. reads of another bacterium on conserved sequence. The dnaA locus has 0% reads
    without mismatches; the known WsnRNAs 85-97%. `min_exact_fraction` 0.5 removes all 5 false
    positives of the mite at `min_rpm` 1 / `min_rpm_genome` 200 / `min_samples` 3 (232 -> 205
    candidates) and keeps the 4 recovered WsnRNA loci (7542, 6108/3640/4274, 89, 744); in JW18 it
    does not touch the cross-talk loci (their reads are exact).
12a. Library threshold 2 or 3 RPM (with `min_rpm_genome` 100, `min_exact_fraction` 0.5, `sense`):
    mite 33 / 24 candidates, 0 false positives, WsnRNA-7542 and -744; JW18 694, 1 cross-talk locus
    in the cured line; wAlbB 334 / 211 (from 109 at 5 RPM: low-load studies gain many candidates
    that no control can judge). Only `min_rpm` 1 with `min_samples` 3 recovers the 6108 family and
    WsnRNA-89 in the mite (205 candidates, 0 false positives).
12. Cross-study concordance (Aag2 vs JW18, same wMel reference) does not rise with stricter
    thresholds (23% -> 18% of Aag2 loci): it is not a usable precision proxy here (different hosts
    and strains, wMelPop vs wMel). The top JW18 locus `TAGAGGCACTGGGATGACAAAAAA` (intergenic,
    229440-229489 -) carries the `GCACTGGGATGAC` motif of the WsnRNA-46 / mite 6108 family.

## 5. Profiles (grid 6; candidates / passing loci in the uninfected line / known small RNAs)

| study | defaults | `symbiont_precise` (sense, exact 0.5, rpm_genome 100, rpm 5) | `symbiont_sensitive` (rpm 2) |
|---|---|---|---|
| wAlbB | 157 / 0 / – | 109 / 0 / – | 334 / 0 / – |
| Aag2 | 570 / 0 / 59 | 172 / 0 / 59 | 172 / 0 / 59 |
| JW18 | 1,957 / 0 / 46, 59 | 694 / 0 / 46, 59 | 694 / 1 / 46, 59 |
| ovary | 0 / 0 / – | 0 / 0 / – | 1 / 0 / – |
| mite | 22 / 2 / 7542 | 14 / 0 / 7542 | 33 / 0 / 7542, 744 |

Decision: repository defaults unchanged (the simulation treats tRNA fragments as targets and the
mitochondrial use case has no evidence here); `symbiont_precise` is the advised profile for a
bacterial focus genome: zero false positives in five studies, every known small RNA found at the
defaults kept, 1.5-3x fewer candidates. For the maintainer's 300+ vs 105 question: run module A
with `symbiont_precise` and compare again; the remaining gap is expected from v0 clustering by
sequence (several v2 peaks can be one v0 cluster and vice versa).

## 6. Open

- Sensitivity at low load is limited by the library threshold; the 6108 family and WsnRNA-89 of
  the mite need `min_rpm` 1 + `min_samples` 3. A replicate-aware rule (e.g. reproducibility across
  pools instead of fixed RPM) could recover them without the cross-talk cost; not tried.
- WsnRNA-46 in Aag2 (~170 reads split in 5' isoforms of the hairpin arm) fails isolation: a
  precursor-aware grouping of arm isoforms would help hairpin-derived small RNAs.
- No v0 comparison was run (about a day per line; only with the maintainer's permission).
