"""IntaRNA (>= 3) on one seed site: command builders and CSV parser.

The interaction is anchored to the site found by the scan: the seed is pinned on the
guide (``--seedQRange``) and on the target (``--seedTRange``), with as many base pairs
as the site class has and no G:U pairs. On the target the interaction may use the site
and the stretch upstream of it where the 3' part of the guide pairs (``window``: guide
length + ``FLANK`` nt upstream of the seed match, ``FLANK`` nt downstream).

Accessibility (ED) is most of IntaRNA's run time and depends only on the target, so it
is computed once per region on the whole region (``accessibility_command``) and reused
for every small RNA and site: each call gets the window only, with the matching slice of
the ED table (``slice_accessibility``, ``--tAcc E``). The energies are those of the whole
region (tested) at a cost per site that does not grow with region length.

Every site gets a score: the cut-offs that make IntaRNA report nothing for weak
interactions (``--outMaxE``, ``--seedMaxE``, both 0 by default) are lifted. Otherwise
weak sites go missing and the decoy distribution is truncated.

IntaRNA numbers sequence 1 = target, sequence 2 = query, positions from 1.
"""

from __future__ import annotations

from .sites import Site

CSV_COLUMNS = ("start1", "end1", "start2", "end2", "E", "E_hybrid", "ED1", "ED2",
               "seedStart1", "seedEnd1", "seedStart2", "seedEnd2", "subseqDP", "hybridDP")
FLANK = 15


def seed_pairing(site: Site, seed: tuple[int, int] = (2, 8)) -> tuple[int, int, int, int]:
    """Paired stretch of a site: (guide_from, guide_to, target_from, target_to).

    Guide positions are 1-based from the 5' end; target positions are 0-based inclusive
    in the region, ``target_to`` pairing with ``guide_from``.
    """
    t1 = site.t1
    if site.site_class in ("8mer", "7mer-m8"):
        return 2, 8, t1 - 7, t1 - 1
    if site.site_class in ("7mer-A1", "6mer"):
        return 2, 7, t1 - 6, t1 - 1
    if site.site_class == "offset-6mer":
        return 3, 8, t1 - 7, t1 - 2
    a, b = seed
    return a, b, site.start, site.end - 1


def window(target: str, site: Site, guide_length: int, seed: tuple[int, int] = (2, 8),
           flank: int = FLANK) -> tuple[int, int]:
    """Stretch of the target (0-based, end excluded) the interaction may use."""
    _, _, t_from, t_to = seed_pairing(site, seed)
    return max(0, t_from - guide_length - flank), min(len(target), t_to + 2 + flank)


def accessibility_command(guide: str, target: str, out_path: str,
                          temperature: float | None = None,
                          binary: str = "IntaRNA") -> list[str]:
    """IntaRNA call that writes the ED table of the whole target to ``out_path``. Any
    query does (the table depends on the target only); its prediction is discarded."""
    cmd = [binary, "-q", guide, "-t", target, "--out", f"tAcc:{out_path}",
           "--outMode", "C", "--threads", "1"]
    if temperature is not None:
        cmd += ["--temperature", str(temperature)]
    return cmd


def slice_accessibility(lines: list[str], lo: int, hi: int) -> str:
    """ED table of target[lo:hi] cut from the table of the whole target.

    A row holds, for the region ending at position i, the ED of its last l bases
    (l = 1, 2, ...). Rows are renumbered from 1 and lengths reaching before ``lo`` are
    set to NA, as IntaRNA writes them at the start of a sequence.
    """
    header = [ln for ln in lines if ln.lstrip().startswith("#")]
    rows = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    out = list(header)
    for k, row in enumerate(rows[lo:hi], start=1):
        values = row.rstrip("\n").rstrip("\t").split("\t")[1:]
        values = values[:k] + ["NA"] * (len(values) - k)
        out.append("\t".join([str(k)] + values) + "\t\n")   # IntaRNA ends rows with a tab
    return "".join(out)


def command(guide: str, target: str, site: Site, seed: tuple[int, int] = (2, 8),
            temperature: float | None = None, binary: str = "IntaRNA",
            accessibility: str | None = None, clip: bool = True) -> tuple[list[str], int]:
    """IntaRNA command for one site, and the offset to add to its target positions to
    get positions in the region (1-based).

    ``clip``: pass only the window (with ``accessibility``, the matching slice of the
    region's ED table; without it IntaRNA computes ED on the window alone). Not clipped:
    the whole region, with the window as ``--tRegion``.
    """
    g_from, g_to, t_from, t_to = seed_pairing(site, seed)
    lo, hi = window(target, site, len(guide), seed)
    if clip:
        sequence, offset, region = target[lo:hi], lo, []
    else:
        sequence, offset, region = target, 0, ["--tRegion", f"{lo + 1}-{hi}"]
    cmd = [binary, "-q", guide, "-t", sequence, *region,
           "--seedQRange", f"{g_from}-{g_to}",
           "--seedTRange", f"{t_from - offset + 1}-{t_to - offset + 1}",
           "--seedBP", str(g_to - g_from + 1), "--seedNoGU",
           "--seedMaxE", "999", "--outMaxE", "999",
           "--outMode", "C", "--outCsvCols", ",".join(CSV_COLUMNS), "--threads", "1"]
    if accessibility:
        cmd += ["--tAcc", "E", "--tAccFile", accessibility]
    if temperature is not None:
        cmd += ["--temperature", str(temperature)]
    return cmd, offset


def parse_csv(text: str, offset: int = 0) -> dict:
    """Scores of the best interaction in IntaRNA's CSV output; {} when there is none.

    ``intarna_target`` is the paired target stretch in region coordinates (1-based).
    """
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if len(lines) < 2:
        return {}
    row = dict(zip(lines[0].split(";"), lines[1].split(";")))
    return {
        "intarna_E": float(row["E"]),
        "intarna_E_hybrid": float(row["E_hybrid"]),
        "intarna_ED_target": float(row["ED1"]),
        "intarna_ED_query": float(row["ED2"]),
        "intarna_target": f"{int(row['start1']) + offset}-{int(row['end1']) + offset}",
        "intarna_guide": f"{row['start2']}-{row['end2']}",
        "intarna_duplex": row["hybridDP"],
    }
