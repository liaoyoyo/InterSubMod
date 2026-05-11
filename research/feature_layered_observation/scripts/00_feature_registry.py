"""Feature Registry helper for Thread F layered observation project.

Central catalog of ~137 features (117 ISM TSV + 11 VCF + 9 BAM QC).
Downstream Phase C agents import `get_feature_metadata(feature_id)` to
retrieve C++ source, computation formula, filter condition, etc.

Usage
-----
>>> from research.feature_layered_observation.scripts.feature_registry import (
...     load_registry, get_feature_metadata, list_features_by_group, get_group_summary
... )
>>> md = get_feature_metadata('HPFineNGroups')
>>> md['source_file']
'src/core/LabelTest.cpp'
>>> md['computation']
'Number of non-empty buckets among {HP1, HP1-1, HP2, HP2-1}; reads with HP0/HP3 excluded'
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

# Registry path (project-relative)
REGISTRY_PATH = (
    Path(__file__).resolve().parents[1] / "data" / "feature_registry.tsv"
)

# Group names (G1..G10) for pretty-print
GROUP_NAMES: dict[str, str] = {
    "G1": "Read Counting & Coverage",
    "G2": "LOH & Subclone Annotation",
    "G3": "Caller VCF Signals",
    "G4": "BAM Read-Level QC",
    "G5": "HP Merged (2-bucket)",
    "G6": "HP Fine (4-bucket LOH-phasing)",
    "G7": "Cluster & Global Methylation Structure",
    "G8": "Entropy, Epipolymorphism & PerCpgASM",
    "G9": "Allele-Specific Methylation (ASM)",
    "G10": "Quality Summary & Verification",
}


_REGISTRY_CACHE: pd.DataFrame | None = None


def load_registry(path: Path | str | None = None) -> pd.DataFrame:
    """Load the feature registry TSV into a DataFrame (cached).

    Parameters
    ----------
    path : optional override; defaults to ``data/feature_registry.tsv``.

    Returns
    -------
    DataFrame indexed by ``feature_id`` with columns: group, source_file,
    source_line, computation, upstream_deps, filter_condition, data_type,
    range, prior_conclusion.
    """
    global _REGISTRY_CACHE
    if _REGISTRY_CACHE is None or path is not None:
        p = Path(path) if path else REGISTRY_PATH
        if not p.exists():
            raise FileNotFoundError(f"Feature registry not found: {p}")
        df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        if df["feature_id"].duplicated().any():
            dups = df.loc[df["feature_id"].duplicated(), "feature_id"].tolist()
            raise ValueError(f"Duplicate feature_id entries in registry: {dups}")
        df = df.set_index("feature_id", drop=False)
        if path is None:
            _REGISTRY_CACHE = df
        return df
    return _REGISTRY_CACHE


def get_feature_metadata(feature_id: str) -> dict[str, Any]:
    """Return metadata for a single feature as a dict.

    Raises ``KeyError`` if the feature is not in the registry.
    """
    df = load_registry()
    if feature_id not in df.index:
        raise KeyError(
            f"Feature '{feature_id}' not in registry. "
            f"Available: {len(df)} features across {df['group'].nunique()} groups."
        )
    row = df.loc[feature_id].to_dict()
    # Add group long name
    row["group_name"] = GROUP_NAMES.get(row.get("group", ""), "UNKNOWN")
    return row


def list_features_by_group(group: str) -> list[str]:
    """List all feature_ids within a group (e.g., 'G1', 'G6')."""
    df = load_registry()
    mask = df["group"] == group
    return df.loc[mask, "feature_id"].tolist()


def get_group_summary() -> pd.DataFrame:
    """Return a DataFrame summarizing feature counts and data-type mix per group."""
    df = load_registry()
    rows = []
    for g, sub in df.groupby("group"):
        rows.append(
            {
                "group": g,
                "group_name": GROUP_NAMES.get(g, "UNKNOWN"),
                "n_features": len(sub),
                "n_continuous": (sub["data_type"] == "continuous").sum(),
                "n_categorical": (sub["data_type"] == "categorical").sum(),
                "n_binary": (sub["data_type"] == "binary").sum(),
                "n_ordinal": (sub["data_type"] == "ordinal").sum(),
                "n_inferred": sub["prior_conclusion"].str.contains(
                    "INFERRED_FROM_NAME", na=False
                ).sum()
                + sub["computation"].str.contains("INFERRED_FROM_NAME", na=False).sum(),
            }
        )
    return pd.DataFrame(rows).sort_values("group")


def features_needing_source_lookup() -> list[str]:
    """Return feature_ids flagged with INFERRED_FROM_NAME (TO DO: resolve)."""
    df = load_registry()
    mask = df["prior_conclusion"].str.contains(
        "INFERRED_FROM_NAME", na=False
    ) | df["computation"].str.contains("INFERRED_FROM_NAME", na=False)
    return df.loc[mask, "feature_id"].tolist()


def main() -> None:
    """CLI smoke test: print group summary + a sample lookup."""
    df = load_registry()
    print(f"Registry loaded: {len(df)} features")
    print("\nGroup summary:")
    print(get_group_summary().to_string(index=False))
    print("\nSample lookup: HPFineNGroups")
    md = get_feature_metadata("HPFineNGroups")
    for k, v in md.items():
        print(f"  {k}: {v}")
    inf = features_needing_source_lookup()
    print(f"\n{len(inf)} features flagged INFERRED_FROM_NAME:")
    for f in inf:
        print(f"  - {f}")


if __name__ == "__main__":
    main()
