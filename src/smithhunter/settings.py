"""The `discover` settings read from the config, in one place.

The workflow (workflow/rules/discovery.smk) and the benchmark tools (05_variants.py, the
regression tests) build the `smithhunter discover` arguments from this list, so a setting added
here reaches all of them; a setting missing here cannot be used from the config.
"""
from __future__ import annotations

# (command line flag, config section, config key, default); None = flag left out
DISCOVER_SETTINGS = [
    ("--merge-gap", "loci", "merge_gap", 0),
    ("--min-rpm", "loci", "min_rpm", 5),
    ("--min-rpm-genome", "loci", "min_rpm_genome", 0),
    ("--min-count", "loci", "min_count", 1),
    ("--min-samples", "loci", "min_samples", None),
    ("--peak-window", "loci", "peak_window", 2),
    ("--peak-pvalue", "loci", "peak_pvalue", 0.001),
    ("--background-flank", "loci", "background_flank", 50),
    ("--satellite-distance", "loci", "satellite_distance", 10),
    ("--satellite-fraction", "loci", "satellite_fraction", 0.1),
    ("--cluster-gap", "loci", "cluster_gap", 0),
    ("--transcript-gap", "loci", "transcript_gap", 0),
    ("--control-group", "loci", "control_group", None),
    ("--min-enrichment", "loci", "min_enrichment", 0),
    ("--exclude-structural", "loci", "exclude_structural", "none"),
    ("--min-exact-fraction", "loci", "min_exact_fraction", 0),
    ("--n-thre", "ends", "n_thre", 0.5),
    ("--penalty", "ends", "penalty", 0.1),
    ("--min-five", "ends", "min_five_score", 0.5),
    ("--min-three", "ends", "min_three_score", 0.0),
    ("--min-isolation", "ends", "min_isolation", 0.5),
    ("--five-filter", "ends", "five_filter", "dominant"),
    ("--min-five-fraction", "ends", "min_five_fraction", 0.5),
    ("--five-jitter", "ends", "five_jitter", 1),
]
# read here but not a plain flag/value pair
SPECIAL = {("loci", "peaks"), ("loci", "min_samples_scope")}
KNOWN = {(section, key) for _, section, key, _ in DISCOVER_SETTINGS} | SPECIAL
SCOPES = ("all", "group")


def value(config: dict, section: str, key: str, default=None):
    return (config.get(section) or {}).get(key, default)


def effective(config: dict) -> dict[str, object]:
    """section.key -> value used by discover (defaults filled in)."""
    out = {f"{s}.{k}": value(config, s, k, d) for _, s, k, d in DISCOVER_SETTINGS}
    out["loci.peaks"] = bool(value(config, "loci", "peaks", True))
    out["loci.min_samples_scope"] = value(config, "loci", "min_samples_scope", "all")
    return out


def discover_args(config: dict, samples_table: str | None = None) -> list[str]:
    """Arguments of `smithhunter discover` for these settings (inputs and outputs excluded).

    samples_table: the samples.tsv with a `group` column, needed with min_samples_scope group
    and with control_group.
    """
    used = effective(config)
    args: list[str] = []
    for flag, section, key, _ in DISCOVER_SETTINGS:
        v = used[f"{section}.{key}"]
        if v is None:
            continue
        if key == "min_samples":
            v = int(v)
        args += [flag, str(v)]
    if not used["loci.peaks"]:
        args.append("--no-peaks")
    scope = used["loci.min_samples_scope"]
    if scope not in SCOPES:
        raise ValueError("config: loci.min_samples_scope must be all or group")
    if scope == "group":
        args += ["--min-samples-scope", "group"]
    if scope == "group" or used["loci.control_group"]:
        if not samples_table:
            raise ValueError("config: loci.min_samples_scope group and loci.control_group need "
                             "the samples table")
        args += ["--samples-table", str(samples_table)]
    return args


def needs_groups(config: dict) -> bool:
    """True when discover reads the `group` column of the samples table."""
    used = effective(config)
    return used["loci.min_samples_scope"] == "group" or bool(used["loci.control_group"])
