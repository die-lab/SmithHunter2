import gzip

from smithhunter.evaluate import evaluate, interval, overlaps
from smithhunter.seqio import read_fasta, read_tsv, write_tsv
from smithhunter.simulate import ADAPTER_R1, FLANK, Params, Simulator, revcomp, simulate


def small_params(tmp_path, **kw):
    return Params(outdir=str(tmp_path / "sim"), reads=3000, nuclear_length=600_000,
                  bacterium_length=20_000, srnas=5, **kw)


def test_simulation_writes_consistent_truth(tmp_path):
    sim = simulate(small_params(tmp_path))
    out = tmp_path / "sim"
    mito = dict(read_fasta(str(out / "data/mito.fasta")))["chrM"]
    features = {f["feature_id"]: f for f in read_tsv(str(out / "truth/features.tsv"))}

    # canonical sequences are the genome (or its reverse complement) at the stated place
    for f in features.values():
        if f["category"] in ("srna", "srna_close_pair", "srna_numt_1mm"):
            s, e = int(f["start"]) - 1, int(f["end"])
            seq = mito[s:e] if f["strand"] == "+" else revcomp(mito[s:e])
            assert seq == f["sequence"]
    origin = features["srna_origin"]
    assert origin["wraps_origin"] == "1"
    assert origin["sequence"] == mito[int(origin["start"]) - 1:] + mito[:int(origin["end"])]

    # the identical NUMT is in the nuclear genome, the 1-mismatch one is not
    nuclear = dict(read_fasta(str(out / "data/nuclear.fasta")))["chrN"]
    assert features["srna_numt_identical"]["sequence"] in nuclear
    assert features["srna_numt_1mm"]["sequence"] not in nuclear

    # FASTQ: reads carry the adapter after the insert, counts match the truth
    with gzip.open(out / "data/sim1_1.fastq.gz", "rt") as fh:
        lines = fh.read().splitlines()
    reads = lines[1::4]
    assert len(reads) == sum(s.counts[0] for s in sim.sources)
    assert sum(ADAPTER_R1[:10] in r for r in reads) > 0.9 * len(reads)
    config = (out / "config/config.yaml").read_text()
    assert "data/mito.gff3" in config and "ambiguous: report" in config


def test_srna_isoforms_follow_the_model(tmp_path):
    sim = Simulator(small_params(tmp_path))
    mito = sim.build_mito()
    src = sim.srna(mito, "x", "srna", 1000, 22, "+", [100.0], "pass")
    inserts = [sim.insert(src, 0) for _ in range(5000)]
    canonical = sum(i == src.sequence for i in inserts) / len(inserts)
    five_exact = sum(i.startswith(src.sequence[:10]) for i in inserts) / len(inserts)
    assert 0.35 < canonical < 0.6  # 0.90 (5') x 0.55 (3') x 0.9 (no tail)
    assert five_exact > 0.85
    assert src.template[FLANK:FLANK + 22] == src.sequence


def write_results(tmp_path, loci):
    res = tmp_path / "results"
    for sub in ("discovery", "reference", "origin"):
        (res / sub).mkdir(parents=True, exist_ok=True)
    write_tsv(str(res / "reference/contigs.tsv"),
              [{"ref": "mito|chrM", "genome": "mito", "contig": "chrM", "length": 1000,
                "circular": 1}], ["ref", "genome", "contig", "length", "circular"])
    write_tsv(str(res / "reference/genomes.tsv"), [{"name": "mito", "role": "focus"}],
              ["name", "role"])
    cols = ["locus_id", "genome", "contig", "start", "end", "strand", "span", "wraps_origin",
            "rep_sequence", "rep_start", "rep_end", "total_count", "five_score", "three_score",
            "annotation_class", "pass"]
    write_tsv(str(res / "discovery/loci.tsv"),
              [dict(zip(cols, l)) for l in loci], cols)
    return str(res)


def test_evaluate_classifies_outcomes(tmp_path):
    truth = tmp_path / "truth"
    truth.mkdir()
    cols = ["feature_id", "category", "genome", "contig", "start", "end", "strand",
            "wraps_origin", "length", "sequence", "expected", "total_reads"]
    feats = [
        ("a", "srna", 101, 122, "+", 0, "AAAA", "pass"),
        ("b", "srna", 301, 322, "+", 0, "", "pass"),     # merged with c
        ("c", "srna", 330, 351, "+", 0, "", "pass"),
        ("d", "srna", 501, 522, "+", 0, "", "pass"),     # missed
        ("e", "srna", 990, 12, "+", 1, "", "pass"),      # across the origin
    ]
    write_tsv(str(truth / "features.tsv"),
              [dict(zip(cols, (f, c, "mito", "chrM", s, e, st, w, 22, q, x, 10)))
               for f, c, s, e, st, w, q, x in feats], cols)
    results = write_results(tmp_path, [
        ("m1", "mito", "chrM", 100, 125, "+", 26, 0, "AAAA", 101, 122, 50, 1, 1, "NA", 1),
        ("m2", "mito", "chrM", 290, 360, "+", 71, 0, "", 301, 322, 50, 1, 1, "NA", 1),
        ("m3", "mito", "chrM", 101, 122, "-", 22, 0, "", 101, 122, 50, 1, 1, "NA", 1),
        ("m4", "mito", "chrM", 991, 10, "+", 20, 1, "", 991, 10, 50, 1, 1, "NA", 0),
    ])
    rows, false_pos, _ = evaluate(str(truth), results)
    observed = {r["feature_id"]: r["observed"] for r in rows}
    assert observed == {"a": "pass", "b": "merged", "c": "merged", "d": "missed", "e": "fail"}
    a = next(r for r in rows if r["feature_id"] == "a")
    assert a["rep_is_canonical"] == 1 and a["rep_five_offset"] == 0 and a["extra_span"] == 4
    assert [fp["locus_id"] for fp in false_pos] == ["m3"]  # antisense, no truth on that strand


def test_circular_overlap():
    assert overlaps((995, 1010), (0, 5), 1000, True)
    assert not overlaps((995, 1000), (0, 5), 1000, False)
    assert interval({"start": "991", "end": "10", "wraps_origin": "1"}, 1000, True) == (990, 1010)
