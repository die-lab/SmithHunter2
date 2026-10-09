"""Collapse trimmed reads of all samples into one table of unique sequences.

Unique sequences are aligned once, which is much faster than aligning every read
and keeps per-sample counts available for every downstream step.
"""

from __future__ import annotations

from .seqio import open_text, read_fastq_sequences, write_tsv


def collapse(samples: dict[str, list[str]], min_length: int, max_length: int):
    """Return (sample names, [(sequence, per-sample counts)] sorted by abundance, stats)."""
    names = list(samples)
    counts: dict[str, list[int]] = {}
    stats = {name: {"sample": name, "reads_in": 0, "reads_kept": 0} for name in names}
    for idx, name in enumerate(names):
        for path in samples[name]:
            for seq in read_fastq_sequences(path):
                stats[name]["reads_in"] += 1
                seq = seq.upper().replace("U", "T")
                if not min_length <= len(seq) <= max_length or "N" in seq:
                    continue
                row = counts.get(seq)
                if row is None:
                    row = counts[seq] = [0] * len(names)
                row[idx] += 1
                stats[name]["reads_kept"] += 1
    ordered = sorted(counts.items(), key=lambda kv: (-sum(kv[1]), kv[0]))
    return names, ordered, [stats[n] for n in names]


def unique_id(index: int) -> str:
    return f"u{index:09d}"


def write_outputs(names, ordered, stats, fasta_path, counts_path, libraries_path) -> None:
    with open_text(fasta_path, "wt") as fa, open_text(counts_path, "wt") as tsv:
        tsv.write("\t".join(["read_id", "sequence", "length", *names]) + "\n")
        for i, (seq, row) in enumerate(ordered, start=1):
            rid = unique_id(i)
            fa.write(f">{rid}\n{seq}\n")
            tsv.write("\t".join([rid, seq, str(len(seq)), *map(str, row)]) + "\n")
    write_tsv(libraries_path, stats, ["sample", "reads_in", "reads_kept"])
