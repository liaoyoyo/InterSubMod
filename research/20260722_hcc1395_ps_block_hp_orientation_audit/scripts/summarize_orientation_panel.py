#!/usr/bin/env python3
"""Summarize a fixed panel of PS-orientation probe JSON documents."""

from __future__ import annotations

import argparse
import csv
import glob
import json
from collections import Counter
from pathlib import Path
from typing import Any


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe-glob", required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()

    paths = [Path(path) for path in sorted(glob.glob(args.probe_glob))]
    if not paths:
        raise FileNotFoundError(args.probe_glob)
    rows: list[dict[str, Any]] = []
    changed_fields = Counter()
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            document = json.load(handle)
        configurations = document.get("configurations", [])
        comparisons = document.get("comparisons_to_native", [])
        if len(configurations) != 2 or len(comparisons) != 1:
            raise ValueError(f"panel expects exactly two relative orientations: {path}")
        native = configurations[0]["region"]
        alternative = configurations[1]["region"]
        changed = comparisons[0].get("changed_fields", {})
        changed_fields.update(changed.keys())
        phase_counts = document["read_evidence"]["phase_set_counts"]
        ordered_counts = sorted((int(value) for value in phase_counts.values()), reverse=True)
        rows.append(
            {
                "region": document["scope"]["region"],
                "k": int(document["scope"]["k"]),
                "ps_read_counts": ",".join(str(value) for value in ordered_counts),
                "reads_with_known_ps": int(document["read_evidence"]["reads_with_known_ps"]),
                "reads_without_ps": int(document["read_evidence"]["reads_without_ps"]),
                "canonical_reproduction": bool(document["canonical_reproduction"]["all_pass"]),
                "orientation_sensitive": bool(document["orientation_sensitive"]),
                "native_primary_families": native["n_primary_families"],
                "alternative_primary_families": alternative["n_primary_families"],
                "native_T": native["T"],
                "alternative_T": alternative["T"],
                "native_Topo": native["Topo"],
                "alternative_Topo": alternative["Topo"],
                "native_structural_class": native["structural_class"],
                "alternative_structural_class": alternative["structural_class"],
                "native_morphology_options": ",".join(native.get("morphology_options", [])),
                "alternative_morphology_options": ",".join(alternative.get("morphology_options", [])),
                "changed_fields": ",".join(sorted(changed)),
                "source_json": str(path.resolve()),
            }
        )

    checks = {
        "panel_size_is_12": len(rows) == 12,
        "all_native_results_reproduce_canonical": all(row["canonical_reproduction"] for row in rows),
        "all_explicit_flips_change_signature": all(row["orientation_sensitive"] for row in rows),
        "all_have_two_nonempty_ps": all(
            len(row["ps_read_counts"].split(",")) == 2
            and all(int(value) >= 8 for value in row["ps_read_counts"].split(","))
            for row in rows
        ),
        "all_k_le_4": all(row["k"] <= 4 for row in rows),
    }
    output = {
        "schema_name": "intersubmod.hcc1395_ps_orientation_enriched_panel",
        "schema_version": "1.0.0",
        "scope": {
            "sample": "HCC1395",
            "panel_type": "enriched sensitivity stress panel; not a prevalence estimate",
            "selection_rule": (
                "nPS=2, retained k<=4, each PS has >=8 tagged reads, and the two blocks have "
                "opposite numeric-HP majorities; top 12 by balance x majority-purity score"
            ),
            "orientation_test": "anchor first PS; explicitly swap HP1/HP2 on second PS only",
            "claim_boundary": (
                "12/12 sensitivity proves the failure mode is real in enriched HCC1395 cases; "
                "it must not be extrapolated to the 865 mixed-PS regions"
            ),
        },
        "counts": {
            "panel_regions": len(rows),
            "canonical_reproduction_pass": sum(row["canonical_reproduction"] for row in rows),
            "orientation_sensitive": sum(row["orientation_sensitive"] for row in rows),
            "changed_fields": dict(sorted(changed_fields.items())),
        },
        "regions": rows,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"output_json": str(args.output_json), "output_tsv": str(args.output_tsv), **output["counts"], "all_checks_pass": output["all_checks_pass"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
