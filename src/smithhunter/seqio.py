"""Minimal FASTA/FASTQ/TSV helpers (plain or gzip-compressed)."""

from __future__ import annotations

import csv
import gzip
from typing import IO, Iterator


def open_text(path: str, mode: str = "rt") -> IO[str]:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)  # type: ignore[return-value]
    return open(path, mode)


def read_fasta(path: str) -> Iterator[tuple[str, str]]:
    """Yield (name, sequence); the name is the header up to the first whitespace."""
    name = None
    chunks: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def read_fasta_names(path: str) -> Iterator[str]:
    for name, _ in read_fasta(path):
        yield name


def read_fastq_sequences(path: str) -> Iterator[str]:
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip()


def read_tsv(path: str) -> list[dict[str, str]]:
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: str, rows: list[dict], columns: list[str]) -> None:
    with open_text(path, "wt") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
