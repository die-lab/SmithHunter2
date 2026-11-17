"""Precursor folding: could a candidate small RNA come from a hairpin, as a miRNA would?

For each candidate two windows are cut on its strand, so that the candidate can sit in
either arm of a hairpin (``flanks = [short, long]``, default [15, 50] as in v0):

- window ``3p``: ``long`` bases upstream (5') and ``short`` downstream (3') of the
  candidate, which would then be in the 3' arm;
- window ``5p``: ``short`` upstream and ``long`` downstream (candidate in the 5' arm).

Coordinates are handled as 0-based half-open internally and written 1-based inclusive,
as in ``loci.tsv``. Windows wrap the origin of circular contigs (written with
start > end on the forward strand, ``wraps_origin`` = 1) and are clipped at the ends of
linear ones (``clipped`` = 1). Minus-strand windows are reverse-complemented, so every
sequence, structure and in-window position reads 5' -> 3' on the RNA.

Per window:

- ``mfe``: RNAfold minimum free energy (kcal/mol, <= 0) at the configured temperature;
- ``amfe`` = -MFE / length x 100 and ``mfei`` = AMFE / GC%, both written as positive
  numbers (Zhang et al. 2006: pre-miRNAs are typically above MFEI ~0.85); MFEI is NA
  when the window has no G or C;
- ``pvalue`` = (k + 1) / (n + 1), with k the number of the n dinucleotide-shuffled
  windows whose MFE is <= the real one (as in randfold, Bonnet et al. 2004, with the
  +1 correction so that it is never 0; the smallest value is 1 / (n + 1));
- the hairpin check of ``hairpin_check``.

The shuffles of a window are drawn from a generator seeded with the random seed, the
candidate id and the window name, so they depend neither on the other candidates nor on
the number of threads. Each window is one RNAfold process (real sequence + shuffles on
stdin); ``threads`` such processes run at a time.

The best window of a candidate has the lowest p-value, then the lowest MFE, then is the
``3p`` one; both windows are kept.
"""

from __future__ import annotations

import random
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

from .seqio import read_fasta, read_tsv
from .targets.shuffle import dinucleotide_shuffle

WINDOWS = ("3p", "5p")  # arm the candidate would occupy
DOMINANT_FRACTION = 0.8  # share of the base pairs in one stem-loop for a "single hairpin"
_COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")

COLUMNS = ["candidate", "genome", "contig", "strand", "candidate_start", "candidate_end",
           "window", "window_start", "window_end", "wraps_origin", "clipped", "length",
           "candidate_from", "candidate_to", "sequence", "structure", "mfe", "amfe", "gc",
           "mfei", "shuffles", "pvalue", "hairpin_loops", "stem_pairs", "stem_fraction",
           "single_hairpin", "in_arm", "arm", "paired_fraction", "best"]


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------- windows

@dataclass
class Window:
    name: str          # "3p" or "5p"
    start: int         # 0-based start on the forward strand, 0 <= start < contig length
    end: int           # 0-based exclusive end; > contig length when the window wraps
    cand_from: int     # 0-based offset of the candidate in the window, on the RNA strand
    cand_to: int       # exclusive
    wraps: bool = False
    clipped: bool = False

    @property
    def length(self) -> int:
        return self.end - self.start


def to_interval(start: int, end: int, length: int) -> tuple[int, int]:
    """1-based inclusive (``loci.tsv``; end < start across the origin) -> 0-based half-open.

    The end of an interval across the origin is returned beyond the contig length.
    """
    start, end = int(start), int(end)
    if not (1 <= start <= length and 1 <= end <= length):
        raise ValueError(f"interval {start}-{end} outside a contig of {length} nt")
    return start - 1, end + length if end < start else end


def make_windows(start: int, end: int, strand: str, length: int, circular: bool,
                 flanks: tuple[int, int] = (15, 50)) -> list[Window]:
    """The two windows around the 0-based half-open candidate ``[start, end)``.

    ``end`` may exceed ``length`` for a candidate across the origin of a circular contig.
    """
    short, long = flanks
    if strand not in "+-" or len(strand) != 1:
        raise ValueError(f"strand must be + or -, got {strand!r}")
    if short < 0 or long < 0:
        raise ValueError("flanks must not be negative")
    if not 0 <= start < end or (end > length and not circular) or end - start > length:
        raise ValueError(f"candidate [{start}, {end}) does not fit a contig of {length} nt")
    out = []
    for name, up, down in (("3p", long, short), ("5p", short, long)):
        # upstream = 5' of the RNA: lower coordinates on +, higher on -
        left, right = (up, down) if strand == "+" else (down, up)
        s, e = start - left, end + right
        clipped = False
        if circular:
            excess = (e - s) - length  # a window cannot be longer than the contig
            if excess > 0:
                clipped = True
                cut = min(right, (excess + 1) // 2)
                e -= cut
                s += excess - cut
        else:
            clipped = s < 0 or e > length
            s, e = max(s, 0), min(e, length)
        offset = start - s if strand == "+" else e - end
        shift = (s // length) * length if circular else 0  # bring the start into the contig
        s, e = s - shift, e - shift
        out.append(Window(name, s, e, offset, offset + (end - start), e > length, clipped))
    return out


def window_sequence(contig_seq: str, window: Window, strand: str) -> str:
    """Sequence of the window on its strand, 5' -> 3'."""
    length = len(contig_seq)
    seq = contig_seq[window.start:min(window.end, length)]
    if window.end > length:
        seq += contig_seq[:window.end - length]
    seq = seq.upper()
    return seq if strand == "+" else reverse_complement(seq)


# ---------------------------------------------------------------------------- RNAfold

def rnafold_command(temperature: float, executable: str = "RNAfold") -> list[str]:
    """RNAfold reading FASTA from stdin, no PostScript files, at ``temperature`` (Celsius)."""
    return [executable, "--noPS", "-T", f"{float(temperature):g}"]


_STRUCTURE = re.compile(r"^([.()]+)\s+\(\s*(-?\d+(?:\.\d+)?)\)$")


def parse_rnafold(text: str) -> list[tuple[str, str, str, float]]:
    """(name, sequence, structure, MFE) for each record of RNAfold's output.

    The sequence is as RNAfold prints it (U instead of T); the name is empty for input
    without headers.
    """
    out = []
    name, seq = "", None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name, seq = line[1:].split()[0] if len(line) > 1 else "", None
            continue
        m = _STRUCTURE.match(line)
        if m and seq is not None:
            if len(m.group(1)) != len(seq):
                raise ValueError(f"RNAfold output: structure and sequence of {name!r} differ in length")
            out.append((name, seq, m.group(1), float(m.group(2))))
            name, seq = "", None
        elif seq is None:
            seq = line
        else:
            raise ValueError(f"RNAfold output: unexpected line {line!r}")
    if seq is not None:
        raise ValueError(f"RNAfold output: no structure for {name or seq!r}")
    return out


def run_rnafold(sequences: list[str], temperature: float,
                executable: str = "RNAfold") -> list[tuple[str, float]]:
    """(structure, MFE) of each sequence, in order, from one RNAfold process."""
    if not sequences:
        return []
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"{executable} not found on PATH: install ViennaRNA (conda: viennarna >= 2.6, see "
            "environment.yml) or run snakemake with --software-deployment-method conda")
    stdin = "".join(f">s{i}\n{s}\n" for i, s in enumerate(sequences))
    proc = subprocess.run(rnafold_command(temperature, path), input=stdin, capture_output=True,
                          text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"RNAfold failed (exit {proc.returncode}): {proc.stderr.strip()[-500:]}")
    records = parse_rnafold(proc.stdout)
    if [r[0] for r in records] != [f"s{i}" for i in range(len(sequences))]:
        raise RuntimeError(f"RNAfold returned {len(records)} records for {len(sequences)} sequences")
    return [(structure, mfe) for _, _, structure, mfe in records]


# ---------------------------------------------------------------------------- measures

def gc_percent(seq: str) -> float:
    return 100.0 * sum(seq.upper().count(b) for b in "GC") / len(seq) if seq else 0.0


def amfe(mfe: float, length: int) -> float:
    """Adjusted MFE: -MFE per 100 nt (positive for a folded sequence)."""
    return -mfe / length * 100 if length else 0.0


def mfei(mfe: float, seq: str) -> float | None:
    """MFE index = AMFE / GC%; None when the sequence has no G or C."""
    gc = gc_percent(seq)
    return amfe(mfe, len(seq)) / gc if gc > 0 else None


def shuffle_pvalue(real: float, shuffled: list[float]) -> float | None:
    """(k + 1) / (n + 1), k = shuffled MFEs <= the real one; None without shuffles."""
    if not shuffled:
        return None
    k = sum(1 for x in shuffled if x <= real + 1e-9)
    return (k + 1) / (len(shuffled) + 1)


def shuffled_windows(seq: str, n: int, random_seed: int, candidate: str, window: str) -> list[str]:
    """``n`` dinucleotide shuffles, reproducible for this candidate and window alone."""
    rng = random.Random(f"{random_seed}:{candidate}:{window}")
    return [dinucleotide_shuffle(seq, rng) for _ in range(n)]


# ---------------------------------------------------------------------------- hairpin

def pair_table(structure: str) -> dict[int, int]:
    """Partner of every paired position of a dot-bracket string (0-based, both ways)."""
    pairs, stack = {}, []
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if not stack:
                raise ValueError(f"unbalanced structure: {structure}")
            j = stack.pop()
            pairs[i], pairs[j] = j, i
        elif c != ".":
            raise ValueError(f"unexpected character {c!r} in structure")
    if stack:
        raise ValueError(f"unbalanced structure: {structure}")
    return pairs


def hairpin_check(structure: str, cand_from: int, cand_to: int,
                  dominant_fraction: float = DOMINANT_FRACTION) -> dict:
    """Describe the structure as a hairpin holding the candidate ``[cand_from, cand_to)``.

    A stem-loop is a hairpin loop with the base pairs that enclose it and no other
    hairpin loop (its stem, with bulges and internal loops, up to the first branching).
    The dominant stem-loop is the one with most base pairs.

    - ``hairpin_loops``: number of hairpin loops (1 = unbranched structure);
    - ``stem_pairs``, ``stem_fraction``: base pairs of the dominant stem-loop and their
      share of all base pairs;
    - ``single_hairpin``: 1 when ``stem_fraction`` >= ``dominant_fraction`` (0.8);
    - ``arm``: ``5p`` / ``3p`` when the candidate bases paired in the dominant stem-loop
      are all on that side, ``both`` when the candidate spans the loop, NA when none is
      paired there; ``in_arm`` = 1 for ``5p`` or ``3p``;
    - ``paired_fraction``: share of the candidate bases that are paired (anywhere).
    """
    pairs = pair_table(structure)
    opening = sorted(i for i, j in pairs.items() if i < j)
    cand = range(cand_from, cand_to)
    out = {"hairpin_loops": 0, "stem_pairs": 0, "stem_fraction": 0.0, "single_hairpin": 0,
           "in_arm": 0, "arm": "NA",
           "paired_fraction": sum(1 for i in cand if i in pairs) / len(cand) if len(cand) else 0.0}
    if not opening:
        return out
    # hairpin loops: pairs with no pair inside
    loops = [i for i in opening if not any(k in pairs for k in range(i + 1, pairs[i]))]
    stems = {}
    for loop in loops:
        others = [x for x in loops if x != loop]
        stems[loop] = [i for i in opening if i <= loop and pairs[i] >= pairs[loop]
                       and not any(i < x and pairs[x] < pairs[i] for x in others)]
    best = max(loops, key=lambda x: (len(stems[x]), -x))
    five = set(stems[best])
    three = {pairs[i] for i in five}
    in_five = any(i in five for i in cand)
    in_three = any(i in three for i in cand)
    arm = "both" if in_five and in_three else "5p" if in_five else "3p" if in_three else "NA"
    fraction = len(five) / len(opening)
    out.update(hairpin_loops=len(loops), stem_pairs=len(five), stem_fraction=fraction,
               single_hairpin=int(fraction >= dominant_fraction - 1e-9),
               in_arm=int(arm in ("5p", "3p")), arm=arm)
    return out


# ---------------------------------------------------------------------------- pipeline

def read_candidate_ids(path: str) -> list[str]:
    """Ids of a (possibly edited) candidates FASTA."""
    return [name for name, _ in read_fasta(path)]


def select_loci(loci_path: str, ids: list[str] | None = None,
                anchor: str = "representative") -> tuple[list[dict], list[str]]:
    """Candidates from ``loci.tsv``: the loci named in ``ids``, else those with pass = 1.

    ``anchor`` chooses the interval taken as the candidate: the most abundant read
    (``representative``, the sequence in candidates.fasta) or the whole ``locus``.
    Returns the candidates and the ids that are not in the table.
    """
    if anchor not in ("representative", "locus"):
        raise ValueError("anchor must be representative or locus")
    keys = ("rep_start", "rep_end") if anchor == "representative" else ("start", "end")
    rows = {r["locus_id"]: r for r in read_tsv(loci_path)}
    wanted = list(dict.fromkeys(ids)) if ids is not None else [k for k, r in rows.items()
                                                               if r.get("pass") == "1"]
    out = [{"candidate": k, "genome": rows[k]["genome"], "contig": rows[k]["contig"],
            "strand": rows[k]["strand"], "start": int(rows[k][keys[0]]),
            "end": int(rows[k][keys[1]])} for k in wanted if k in rows]
    return out, [k for k in wanted if k not in rows]


def load_contigs(genomes_path: str, contigs_path: str, needed: set[tuple[str, str]]) -> dict:
    """(genome, contig) -> (sequence, circular) for the contigs that carry candidates."""
    fasta = {g["name"]: g["fasta"] for g in read_tsv(genomes_path)}
    circular = {(c["genome"], c["contig"]): c["circular"] == "1" for c in read_tsv(contigs_path)}
    out = {}
    for genome in sorted({g for g, _ in needed}):
        if genome not in fasta:
            raise ValueError(f"genome {genome!r} of the candidates is not in {genomes_path}")
        for contig, seq in read_fasta(fasta[genome]):
            if (genome, contig) in needed:
                out[(genome, contig)] = (seq, circular.get((genome, contig), False))
    missing = needed - set(out)
    if missing:
        raise ValueError(f"contigs of the candidates not found in the genome FASTA: {sorted(missing)}")
    return out


def fold_candidates(candidates: list[dict], contigs: dict, temperature: float,
                    flanks: tuple[int, int] = (15, 50), shuffles: int = 1000,
                    random_seed: int = 1, threads: int = 1, executable: str = "RNAfold",
                    fold=run_rnafold) -> list[dict]:
    """Rows of ``precursors.tsv``: two per candidate, in the order of ``candidates``.

    ``fold(sequences, temperature, executable)`` returns (structure, MFE) per sequence;
    it is a parameter so the pipeline can be tested without RNAfold.
    """
    jobs = []
    for c in candidates:
        seq, circular = contigs[(c["genome"], c["contig"])]
        start, end = to_interval(c["start"], c["end"], len(seq))
        for w in make_windows(start, end, c["strand"], len(seq), circular, flanks):
            jobs.append((c, w, window_sequence(seq, w, c["strand"]), len(seq)))

    def work(job):
        c, w, wseq, _ = job
        return fold([wseq] + shuffled_windows(wseq, shuffles, random_seed, c["candidate"], w.name),
                    temperature, executable)

    with ThreadPoolExecutor(max_workers=max(1, int(threads))) as pool:
        folded = list(pool.map(work, jobs))

    rows = []
    for (c, w, wseq, length), result in zip(jobs, folded):
        structure, mfe = result[0]
        index = mfei(mfe, wseq)
        pvalue = shuffle_pvalue(mfe, [m for _, m in result[1:]])
        hp = hairpin_check(structure, w.cand_from, w.cand_to)
        rows.append({
            "candidate": c["candidate"], "genome": c["genome"], "contig": c["contig"],
            "strand": c["strand"], "candidate_start": c["start"], "candidate_end": c["end"],
            "window": w.name, "window_start": w.start + 1,
            "window_end": w.end - length if w.wraps else w.end,
            "wraps_origin": int(w.wraps), "clipped": int(w.clipped), "length": w.length,
            "candidate_from": w.cand_from + 1, "candidate_to": w.cand_to,
            "sequence": wseq, "structure": structure, "mfe": f"{mfe:.2f}",
            "amfe": f"{amfe(mfe, len(wseq)):.2f}", "gc": f"{gc_percent(wseq):.1f}",
            "mfei": "NA" if index is None else f"{index:.3f}",
            "shuffles": len(result) - 1, "pvalue": "NA" if pvalue is None else f"{pvalue:.4g}",
            **hp, "stem_fraction": f"{hp['stem_fraction']:.3f}",
            "paired_fraction": f"{hp['paired_fraction']:.3f}",
            "_rank": (pvalue if pvalue is not None else 1.0, mfe, WINDOWS.index(w.name)),
        })
    for i in range(0, len(rows), len(WINDOWS)):
        pair = rows[i:i + len(WINDOWS)]
        best = min(pair, key=lambda r: r["_rank"])
        for r in pair:
            r["best"] = int(r is best)
    return rows
