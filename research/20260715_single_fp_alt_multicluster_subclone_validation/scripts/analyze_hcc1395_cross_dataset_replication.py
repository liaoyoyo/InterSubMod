#!/usr/bin/env python3
"""Audit exact-site replication between the two HCC1395 ONT datasets."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


ENDPOINTS = [
    "stable_null_multigroup",
    "residual_unexplained_multigroup",
    "phase_anchored_robust_epigenetic_candidate",
]


def active(row: dict[str, str], field: str) -> bool:
    return row.get(field, "").lower() == "true"


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--site-results",
        type=Path,
        default=root / "results/focal_alt_multicluster/latest_full_v1/latest_site_results_with_topology.tsv",
    )
    parser.add_argument(
        "--summary", type=Path, default=root / "results/hcc1395_cross_dataset_replication.json"
    )
    args = parser.parse_args()

    with args.site_results.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    key = lambda row: (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
    first = {key(row): row for row in rows if row["sample"] == "HCC1395"}
    second = {key(row): row for row in rows if row["sample"] == "HCC1395_DORADO"}
    common = set(first) & set(second)

    endpoints = {}
    for field in ENDPOINTS:
        first_yes = {site for site in common if active(first[site], field)}
        second_yes = {site for site in common if active(second[site], field)}
        union = first_yes | second_yes
        endpoints[field] = {
            "hcc1395_n": len(first_yes),
            "dorado_n": len(second_yes),
            "both_n": len(first_yes & second_yes),
            "either_n": len(union),
            "jaccard": len(first_yes & second_yes) / len(union) if union else None,
        }

    replicated_e4 = []
    for site in sorted(common):
        if active(first[site], ENDPOINTS[-1]) and active(second[site], ENDPOINTS[-1]):
            replicated_e4.append(
                {
                    "site": f"{site[0]}:{site[1]}:{site[2]}>{site[3]}",
                    "HCC1395": {
                        name: first[site].get(name)
                        for name in (
                            "n_alt_raw",
                            "cluster_sizes",
                            "alt_hp_counts",
                            "caller_af",
                            "caller_gq",
                            "longphase_filter_transition",
                            "layered_topology_context",
                        )
                    },
                    "HCC1395_DORADO": {
                        name: second[site].get(name)
                        for name in (
                            "n_alt_raw",
                            "cluster_sizes",
                            "alt_hp_counts",
                            "caller_af",
                            "caller_gq",
                            "longphase_filter_transition",
                            "layered_topology_context",
                        )
                    },
                }
            )

    summary = {
        "schema_name": "intersubmod.hcc1395_cross_dataset_exact_site_replication",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_hcc1395_sites": len(first),
        "n_hcc1395_dorado_sites": len(second),
        "n_common_sites": len(common),
        "n_common_both_evaluable": sum(
            first[site]["analysis_status"] == "evaluable"
            and second[site]["analysis_status"] == "evaluable"
            for site in common
        ),
        "endpoints": endpoints,
        "replicated_e4": replicated_e4,
        "guardrail": (
            "Exact-site replication supports a recurrent local epigenetic pattern, but both sites remain "
            "truth-labeled FP and do not independently establish a cellular subclone or lineage topology."
        ),
        "pass": len(common) == 204 and len(replicated_e4) == 2,
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
