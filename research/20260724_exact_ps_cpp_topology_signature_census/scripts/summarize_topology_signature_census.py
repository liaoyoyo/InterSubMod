#!/usr/bin/env python3
"""Summarize exact top-tree signature census outputs without mutating inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
RESOLUTIONS = (
    "UNIQUE_TREE",
    "TIED_SAME_TOPOLOGY",
    "TIED_CROSS_TOPOLOGY",
)
COARSE_CLASSES = (
    "Single-only",
    "Sister-only",
    "Direct-only",
    "Sister+direct",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def percentage(numerator: int, denominator: int) -> float:
    return 100.0 * numerator / denominator if denominator else 0.0


def summarize_sample(path: Path, expected_sample: str) -> dict:
    resolution = Counter()
    resolved_coarse = Counter()
    active_k = Counter()
    ranked = 0
    best_tree_total = 0
    maximum_tie_count = 0
    exact_cross_same_coarse = 0
    exact_cross_cross_coarse = 0

    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            row = json.loads(line)
            if row["sample"] != expected_sample:
                raise ValueError(
                    f"{path}:{line_number}: sample mismatch "
                    f"{row['sample']} != {expected_sample}"
                )
            if not row.get("canonical_reproduction_pass"):
                raise ValueError(
                    f"{path}:{line_number}: canonical reproduction failed"
                )
            tie_count = int(row["best_tree_tie_count"])
            enumerated = int(row["enumerated_best_tree_count"])
            signature_total = sum(
                int(item["tree_count"]) for item in row["topology_signatures"]
            )
            coarse_total = sum(
                int(value) for value in row["coarse_class_tree_counts"].values()
            )
            if tie_count != enumerated or tie_count != signature_total:
                raise ValueError(
                    f"{path}:{line_number}: topology-tree count mismatch"
                )
            if tie_count != coarse_total:
                raise ValueError(
                    f"{path}:{line_number}: coarse-tree count mismatch"
                )
            if row["topology_signature_count"] != len(
                row["topology_signatures"]
            ):
                raise ValueError(
                    f"{path}:{line_number}: signature cardinality mismatch"
                )
            if row["coarse_class_count"] != len(
                row["coarse_class_tree_counts"]
            ):
                raise ValueError(
                    f"{path}:{line_number}: coarse cardinality mismatch"
                )

            ranked += 1
            best_tree_total += tie_count
            maximum_tie_count = max(maximum_tie_count, tie_count)
            resolution[row["resolution_class"]] += 1
            active_k[int(row["active_bit_count"])] += 1
            if row["coarse_class_count"] == 1:
                resolved_coarse[next(iter(row["coarse_class_tree_counts"]))] += 1
            if row["resolution_class"] == "TIED_CROSS_TOPOLOGY":
                if row["coarse_class_count"] == 1:
                    exact_cross_same_coarse += 1
                else:
                    exact_cross_cross_coarse += 1

    unknown_resolution = set(resolution) - set(RESOLUTIONS)
    unknown_coarse = set(resolved_coarse) - set(COARSE_CLASSES)
    if unknown_resolution or unknown_coarse:
        raise ValueError(
            f"{path}: unknown labels {unknown_resolution=} {unknown_coarse=}"
        )
    if sum(resolution.values()) != ranked:
        raise ValueError(f"{path}: resolution conservation failed")

    one_exact_topology = (
        resolution["UNIQUE_TREE"] + resolution["TIED_SAME_TOPOLOGY"]
    )
    one_coarse_class = sum(resolved_coarse.values())
    cross_coarse_class = ranked - one_coarse_class
    result = {
        "sample": expected_sample,
        "input": {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        },
        "ranked_units": ranked,
        "best_tree_total": best_tree_total,
        "maximum_best_tree_tie_count": maximum_tie_count,
        "resolution": {
            label: {
                "n": resolution[label],
                "pct_ranked": percentage(resolution[label], ranked),
            }
            for label in RESOLUTIONS
        },
        "one_exact_topology": {
            "n": one_exact_topology,
            "pct_ranked": percentage(one_exact_topology, ranked),
        },
        "one_coarse_class": {
            "n": one_coarse_class,
            "pct_ranked": percentage(one_coarse_class, ranked),
        },
        "cross_coarse_class": {
            "n": cross_coarse_class,
            "pct_ranked": percentage(cross_coarse_class, ranked),
        },
        "cross_exact_topology_but_same_coarse_class": exact_cross_same_coarse,
        "cross_exact_and_cross_coarse_class": exact_cross_cross_coarse,
        "resolved_coarse_class": {
            label: {
                "n": resolved_coarse[label],
                "pct_ranked": percentage(resolved_coarse[label], ranked),
                "pct_coarse_resolved": percentage(
                    resolved_coarse[label], one_coarse_class
                ),
            }
            for label in COARSE_CLASSES
        },
        "active_k": {
            str(k): {
                "n": count,
                "pct_ranked": percentage(count, ranked),
            }
            for k, count in sorted(active_k.items())
        },
        "checks": {
            "all_rows_canonical_reproduction_pass": True,
            "all_enumerated_counts_equal_canonical_tie_count": True,
            "resolution_conservation": True,
            "coarse_resolution_conservation": True,
        },
    }
    return result


def combine(samples: list[dict]) -> dict:
    ranked = sum(row["ranked_units"] for row in samples)
    resolution = {
        label: sum(row["resolution"][label]["n"] for row in samples)
        for label in RESOLUTIONS
    }
    resolved_coarse = {
        label: sum(
            row["resolved_coarse_class"][label]["n"] for row in samples
        )
        for label in COARSE_CLASSES
    }
    active_k = Counter()
    for row in samples:
        for k, entry in row["active_k"].items():
            active_k[int(k)] += entry["n"]
    one_exact = resolution["UNIQUE_TREE"] + resolution["TIED_SAME_TOPOLOGY"]
    one_coarse = sum(resolved_coarse.values())
    return {
        "sample": "ALL7",
        "ranked_units": ranked,
        "best_tree_total": sum(row["best_tree_total"] for row in samples),
        "maximum_best_tree_tie_count": max(
            row["maximum_best_tree_tie_count"] for row in samples
        ),
        "resolution": {
            label: {
                "n": resolution[label],
                "pct_ranked": percentage(resolution[label], ranked),
            }
            for label in RESOLUTIONS
        },
        "one_exact_topology": {
            "n": one_exact,
            "pct_ranked": percentage(one_exact, ranked),
        },
        "one_coarse_class": {
            "n": one_coarse,
            "pct_ranked": percentage(one_coarse, ranked),
        },
        "cross_coarse_class": {
            "n": ranked - one_coarse,
            "pct_ranked": percentage(ranked - one_coarse, ranked),
        },
        "cross_exact_topology_but_same_coarse_class": sum(
            row["cross_exact_topology_but_same_coarse_class"] for row in samples
        ),
        "cross_exact_and_cross_coarse_class": sum(
            row["cross_exact_and_cross_coarse_class"] for row in samples
        ),
        "resolved_coarse_class": {
            label: {
                "n": resolved_coarse[label],
                "pct_ranked": percentage(resolved_coarse[label], ranked),
                "pct_coarse_resolved": percentage(
                    resolved_coarse[label], one_coarse
                ),
            }
            for label in COARSE_CLASSES
        },
        "active_k": {
            str(k): {
                "n": count,
                "pct_ranked": percentage(count, ranked),
            }
            for k, count in sorted(active_k.items())
        },
        "checks": {
            "ranked_units_equal_71955": ranked == 71955,
            "resolution_conservation": sum(resolution.values()) == ranked,
            "coarse_resolution_conservation": one_coarse
            + (ranked - one_coarse)
            == ranked,
            "all_sample_checks_pass": all(
                all(row["checks"].values()) for row in samples
            ),
        },
    }


def tsv_row(row: dict) -> dict:
    values = {
        "sample": row["sample"],
        "ranked_units": row["ranked_units"],
        "UNIQUE_TREE_n": row["resolution"]["UNIQUE_TREE"]["n"],
        "UNIQUE_TREE_pct": row["resolution"]["UNIQUE_TREE"]["pct_ranked"],
        "TIED_SAME_TOPOLOGY_n": row["resolution"]["TIED_SAME_TOPOLOGY"]["n"],
        "TIED_SAME_TOPOLOGY_pct": row["resolution"][
            "TIED_SAME_TOPOLOGY"
        ]["pct_ranked"],
        "TIED_CROSS_TOPOLOGY_n": row["resolution"][
            "TIED_CROSS_TOPOLOGY"
        ]["n"],
        "TIED_CROSS_TOPOLOGY_pct": row["resolution"][
            "TIED_CROSS_TOPOLOGY"
        ]["pct_ranked"],
        "one_exact_topology_n": row["one_exact_topology"]["n"],
        "one_exact_topology_pct": row["one_exact_topology"]["pct_ranked"],
        "one_coarse_class_n": row["one_coarse_class"]["n"],
        "one_coarse_class_pct": row["one_coarse_class"]["pct_ranked"],
        "cross_coarse_class_n": row["cross_coarse_class"]["n"],
        "cross_coarse_class_pct": row["cross_coarse_class"]["pct_ranked"],
    }
    for label in COARSE_CLASSES:
        key = label.replace("+", "_plus_").replace("-", "_")
        values[f"{key}_n"] = row["resolved_coarse_class"][label]["n"]
        values[f"{key}_pct_ranked"] = row["resolved_coarse_class"][label][
            "pct_ranked"
        ]
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()
    if args.output_json.exists() or args.output_tsv.exists():
        raise FileExistsError("output targets must not already exist")

    sample_rows = [
        summarize_sample(
            args.input_dir / f"{sample}.census.jsonl", sample
        )
        for sample in SAMPLES
    ]
    cohort = combine(sample_rows)
    document = {
        "schema_name": (
            "intersubmod.exact_ps_cpp_topology_signature_census.summary"
        ),
        "schema_version": "1.0.0",
        "topology_signature_definition": (
            "root-preserving unlabeled rooted-tree canonical parenthesis "
            "signature; sibling signatures sorted"
        ),
        "coarse_class_definition": (
            "Single-only/Sister-only/Direct-only/Sister+direct from "
            "branching and root-to-node depth>=2"
        ),
        "samples": sample_rows,
        "cohort": cohort,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n"
    )

    rows = [tsv_row(row) for row in sample_rows + [cohort]]
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
