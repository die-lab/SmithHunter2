"""Build one combined reference from all user-declared genomes.

Every contig is renamed ``<genome>|<contig>`` so that each alignment can be traced
back to its genome. Circular contigs get the first ``overhang`` bases appended, so
reads spanning the origin align end-to-end; coordinates are folded back later.
"""

from __future__ import annotations

from dataclasses import dataclass

from .seqio import open_text, read_fasta, write_tsv

SEP = "|"
ROLES = ("focus", "exclude")
TOPOLOGIES = ("linear", "circular")
INTERACTION_MODES = ("ago", "bacterial", "none")


@dataclass
class GenomeSpec:
    name: str
    fasta: str
    role: str = "focus"
    circular: bool = False
    annotation: str = ""
    min_length: int = 0
    max_length: int = 0
    targets: str = "none"


def parse_genomes(spec: dict, default_min: int, default_max: int) -> list[GenomeSpec]:
    """Validate the ``genomes:`` block of the configuration."""
    if not spec:
        raise ValueError("config: 'genomes' must declare at least one genome")
    genomes = []
    for name, entry in spec.items():
        if SEP in name or not name:
            raise ValueError(f"config: genome name {name!r} must be non-empty and contain no '{SEP}'")
        if "fasta" not in entry:
            raise ValueError(f"config: genome {name!r} has no 'fasta'")
        role = entry.get("role", "focus")
        topology = entry.get("topology", "linear")
        targets = entry.get("targets", "none") or "none"
        if role not in ROLES:
            raise ValueError(f"config: genome {name!r} role must be one of {ROLES}, got {role!r}")
        if topology not in TOPOLOGIES:
            raise ValueError(f"config: genome {name!r} topology must be one of {TOPOLOGIES}")
        if targets not in INTERACTION_MODES:
            raise ValueError(f"config: genome {name!r} targets must be one of {INTERACTION_MODES}")
        length = entry.get("length") or [default_min, default_max]
        if len(length) != 2 or not 0 < int(length[0]) <= int(length[1]):
            raise ValueError(f"config: genome {name!r} length must be [min, max]")
        genomes.append(GenomeSpec(
            name=name, fasta=entry["fasta"], role=role, circular=topology == "circular",
            annotation=entry.get("annotation") or "", min_length=int(length[0]),
            max_length=int(length[1]), targets=targets,
        ))
    if not any(g.role == "focus" for g in genomes):
        raise ValueError("config: at least one genome must have role 'focus'")
    return genomes


def build_reference(genomes: list[GenomeSpec], fasta_out: str, overhang: int):
    """Write the combined FASTA and return the contig table rows."""
    contigs = []
    seen = set()
    with open_text(fasta_out, "wt") as out:
        for g in genomes:
            for contig, seq in read_fasta(g.fasta):
                if SEP in contig:
                    raise ValueError(f"{g.fasta}: contig name {contig!r} contains '{SEP}'")
                ref = f"{g.name}{SEP}{contig}"
                if ref in seen:
                    raise ValueError(f"{g.fasta}: duplicated contig {contig!r}")
                seen.add(ref)
                seq = seq.upper()
                extended = seq + seq[:min(overhang, len(seq))] if g.circular else seq
                out.write(f">{ref}\n")
                for i in range(0, len(extended), 80):
                    out.write(extended[i:i + 80] + "\n")
                contigs.append({"ref": ref, "genome": g.name, "contig": contig,
                                "length": len(seq), "circular": int(g.circular)})
    return contigs


def write_tables(genomes: list[GenomeSpec], contigs: list[dict], genomes_out: str, contigs_out: str):
    write_tsv(genomes_out, [g.__dict__ | {"circular": int(g.circular)} for g in genomes],
              ["name", "fasta", "role", "circular", "annotation", "min_length", "max_length", "targets"])
    write_tsv(contigs_out, contigs, ["ref", "genome", "contig", "length", "circular"])
