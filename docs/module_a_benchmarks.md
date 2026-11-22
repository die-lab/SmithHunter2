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

### Condition-aware `min_samples` (bea877c)

`loci.min_samples_scope: group` with a `group` column in `samples.tsv`: a locus passes the
expression filter in group size − 1 samples of at least one group. Mite (sex as group), all with 0
passing loci in the cured line:

| variant | candidates | WsnRNA loci |
|---|---|---|
| `symbiont_precise` | 14 | 7542 |
| `symbiont_precise` + group | 37 | 7542, 744 |
| `symbiont_sensitive` + group | 112 | 7542, 744 |
| `symbiont_precise` + group + `min_rpm` 1 | 238 | 7542, 744, 6108 family, 89 |

Advice: with mixed conditions use `group` rather than lowering `min_samples` by hand.

## 6. Open

- Sensitivity at low load is limited by the library threshold; the 6108 family and WsnRNA-89 of
  the mite need `min_rpm` 1 (with group scope: 238 candidates, no false positive here, but in JW18
  `min_rpm` 1 lets cross-talk of the top symbiont small RNAs pass in the cured line).
- WsnRNA-46 in Aag2: checked, grouping of 5' isoforms **not implemented**. The precursor
  (1204493-1204642 +) carries five peak loci with 5' ends at 1204507, 511, 521, 525 and 532
  (21, 52, 50, 6 and 9 reads; ragged 5' ends over 25 nt, not two hairpin arms). Merged they
  would reach 7.2 / 6.9 RPM of the library (passes the default 5) but only 38 / 51 reads per
  million Wolbachia reads, far below `min_rpm_genome` 100 of the advised profiles: grouping
  would rescue it at the defaults only, while merging neighbouring peaks risks joining
  distinct small RNAs everywhere else. One control at very low load is not enough evidence.
- No v0 comparison was run (about a day per line; only with the maintainer's permission).

## 7. Tooling after the review (2026-11-22)

- **One list of discover settings** (`src/smithhunter/settings.py`): `discovery.smk`,
  `05_variants.py` and the regression tests build the `discover` arguments from it; tests check
  it against the command line parser and `config/config.yaml`.
- **Real Snakemake run** of `discover` with `symbiont_precise` + `min_samples_scope: group` on
  the spider mite infected line (upstream results linked, only `discover` planned and run):
  24,470 loci, every column identical to the `05_variants.py` result of the same settings,
  37 passing Wolbachia loci in both.
- **Regression fixtures** (`tests/data/{jw18,mite}_{infected,tet}`, ~400 kB, made by
  `benchmarks/wolbachia_multi/make_fixture.py`): reads within 400 nt of WsnRNA-46/59, the JW18
  cross-talk loci, WsnRNA-7542/744 and the two mite false positives; reads left out are summed
  into one filler read so `rpm_genome` stays exact. Checked against the full runs (59 / 10 / 21 /
  3 loci identical). `tests/test_benchmark_regression.py` asserts the results of sections 2-5.
- **Clusters** (`cluster_id`, `cluster_size`, `cluster_count`; report only): Aag2 at the
  defaults has 570 passing loci in 513 clusters (57 are further 5' isoforms of a passing locus);
  WsnRNA-46 shows as one cluster of 6 peaks, 157 reads, none passing.

## 8. Beyond Wolbachia, and profile v2 (2026-11-22)

Studies (`prepare.py --study ...`; same layout, `FOCUS=<genome>` for `03_run.sh`):

| study | design | finding |
|---|---|---|
| `chlamydia` (Furuse 2014, PRJNA248348) | *C. trachomatis* in HeLa (total + Ago-IP) vs human THP-1 with *M. tuberculosis*; host GRCh38; 11 sRNAs of Albrecht 2010 (`positive_ctr.tsv`, L2b coordinates) | defaults: 1,837 candidates, **14 false positives, all rRNA** (Mtb rRNA reads on the two Chlamydia rRNA operons, sense and antisense; exact fraction 0.01-0.78), 8 of 11 sRNAs |
| `buchnera_myzus` (Feng 2018, PRJNA395356) | *M. persicae* bacteriocytes vs gut of lines G002, G006, BTIRed; *Buchnera* G002 + pLeu. G006 runs paired with swapped mates: the mate with the 3' adapter is used | the gut holds 0.5% *Buchnera* reads (bacteriocytes 86%): 76 passing loci, 73 also in the bacteriocytes. Contamination, not misassignment: **not a clean negative** |
| `cardinium` (Li 2018, PRJNA321074) | *B. tabaci* Q, Cardinium-infected vs uninfected, 2 + 2; *Portiera* as exclude genome | the "uninfected" line has more *Cardinium* reads (0.74% vs 0.57%), same exact fraction (0.71 vs 0.74): it carries Cardinium (or labels are swapped). **Not a usable negative** |

Profile grid on the six clean studies (candidates / passing in the negative / known sRNAs):

| study | v1 precise (sense, rpm_genome 100) | v2 precise (any, rpm_genome 0) | v2 sensitive (rpm 2) |
|---|---|---|---|
| wAlbB | 109 / 0 / – | 104 / 0 / – | 332 / 0 / – |
| Aag2 | 172 / 0 / 59 | 497 / 0 / 59 | 1,289 / 0 / 59 |
| JW18 | 694 / 0 / 46, 59 | 1,808 / 0 / 46, 59 | 3,694 / 1 / 46, 59 |
| ovary | 0 / 0 / – | 0 / 0 / – | 1 / 0 / – |
| mite | 14 / 0 / 7542 | 14 / 0 / 7542 (group: 37 / 0 / +744) | 33 / 0 / 7542, 744 |
| Chlamydia | 270 / **2** / **4 of 11** | 1,628 / 0 / 8 of 11 | 3,566 / 0 / 11 of 11 |

Decisions:
- **`min_rpm_genome` off in both profiles.** In no study did it remove a false positive (they all
  go with the structural and exact-read filters); it only cut candidates at high load, and in
  Chlamydia (23% of the library) it dropped ctrR3/4/6/7, Northern-validated sRNAs. The candidate
  reduction it gave in JW18 / Aag2 was not shown to be a precision gain.
- **`exclude_structural: any` in both profiles.** `sense` left the 2 antisense rRNA loci of the
  Chlamydia negative; `any` removes them and loses no known small RNA in any study (none lies
  antisense to a structural RNA).
- Regression fixtures `tests/data/chlamydia_{infected,tet}` (positives and both rRNA operons) pin
  these results.

Open: a clean negative is the exception for symbionts (contaminated tissue, residual infection,
cross-talk): an explicit infected vs control comparison (enrichment per locus) is the next filter
to try, on Buchnera gut and JW18 cross-talk. The *A. pisum* Buchnera benchmark with the Hansen &
Degnan 2014 / Thairu & Hansen 2019 sRNAs waits for their supplementary tables (PMC refuses
automated downloads).
