#!/usr/bin/env python3
"""Scan one frozen M2 extraction for exact 10/01/11 edge-score eligibility."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import pathlib
from collections import Counter, defaultdict
from itertools import combinations


def read_gzip_rows(path: pathlib.Path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extraction-dir", required=True)
    parser.add_argument("--dataset", default="HCC1395_DORADO")
    parser.add_argument("--chrom", default="chr6")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    root = pathlib.Path(args.extraction_dir).resolve()
    prefix = f"{args.dataset}.{args.chrom}"
    components_path = root / f"{prefix}.components.tsv.gz"
    membership_path = root / f"{prefix}.site_component_membership.tsv.gz"
    calls_path = root / f"{prefix}.molecule_sparse_calls.tsv.gz"

    components: dict[str, dict] = {}
    by_threshold = Counter()
    for row in read_gzip_rows(components_path):
        if row["inference_role"] != "PRIMARY_PS_AWARE" or int(row["k"]) != 2:
            continue
        component_id = row["component_id"]
        family = row["linkage_basis"].removeprefix("PS_HP")
        components[component_id] = {
            "component_id": component_id,
            "threshold": int(row["threshold"]),
            "family": family,
            "phase_set": row["phase_set"],
            "start1": int(row["start1"]),
            "end1": int(row["end1"]),
            "sites": [],
        }
        by_threshold[int(row["threshold"])] += 1

    for row in read_gzip_rows(membership_path):
        component = components.get(row["component_id"])
        if component is not None:
            component["sites"].append(int(row["site_index"]))
    bad_membership = [
        component_id
        for component_id, component in components.items()
        if len(set(component["sites"])) != 2
    ]
    if bad_membership:
        raise RuntimeError(f"k=2 components with non-two-site membership: {bad_membership[:5]}")

    lookup: dict[tuple[str, str, tuple[int, int]], list[str]] = defaultdict(list)
    for component_id, component in components.items():
        pair = tuple(sorted(set(component["sites"])))
        component["sites"] = list(pair)
        lookup[(component["family"], component["phase_set"], pair)].append(component_id)

    state_counts: dict[str, Counter] = defaultdict(Counter)
    n_call_rows = 0
    n_fixed_pairs_checked = 0
    calls_header: list[str] = []
    with gzip.open(calls_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        calls_header = list(reader.fieldnames or [])
        for row in reader:
            n_call_rows += 1
            family = row["hp_family"]
            phase_set = row["phase_set"]
            indices = [int(value) for value in row["site_indices"].split(",") if value]
            codes = list(row["call_codes"])
            fixed = {
                index: code
                for index, code in zip(indices, codes)
                if code in {"R", "A"}
            }
            for left, right in combinations(sorted(fixed), 2):
                n_fixed_pairs_checked += 1
                component_ids = lookup.get((family, phase_set, (left, right)), ())
                if not component_ids:
                    continue
                state = fixed[left] + fixed[right]
                for component_id in component_ids:
                    state_counts[component_id][state] += 1

    threshold_rows = []
    for threshold in sorted(by_threshold):
        subset = [
            component
            for component in components.values()
            if component["threshold"] == threshold
        ]
        eligible_ge1 = []
        eligible_ge_threshold = []
        state_presence = Counter()
        max_min_count = 0
        for component in subset:
            counts = state_counts[component["component_id"]]
            values = [counts["AR"], counts["RA"], counts["AA"]]
            state_presence[sum(value > 0 for value in values)] += 1
            max_min_count = max(max_min_count, min(values))
            if min(values) >= 1:
                eligible_ge1.append(component["component_id"])
            if min(values) >= threshold:
                eligible_ge_threshold.append(component["component_id"])
        threshold_rows.append(
            {
                "threshold": threshold,
                "n_primary_k2_components": len(subset),
                "n_with_AR_RA_AA_each_ge_1": len(eligible_ge1),
                "n_with_AR_RA_AA_each_ge_threshold": len(eligible_ge_threshold),
                "max_over_components_of_min_AR_RA_AA_count": max_min_count,
                "n_components_by_number_of_present_target_states": {
                    str(key): value for key, value in sorted(state_presence.items())
                },
                "eligible_component_ids_ge_1": eligible_ge1,
            }
        )

    receipt = {
        "schema": "intersubmod.solver_methyl_edge_probe.real_feasibility.v1",
        "scope": "HCC1395_DORADO_CHR6_PARTIAL",
        "inputs": {
            "components": str(components_path),
            "membership": str(membership_path),
            "molecule_sparse_calls": str(calls_path),
        },
        "contract": (
            "One primary PS-aware k=2 component is edge-score state-eligible only "
            "if exact AR, RA and AA molecules all exist in the same HP family and PS."
        ),
        "calls_table_has_mm_ml_columns": (
            "MM" in calls_header or "ML" in calls_header
        ),
        "n_molecule_call_rows": n_call_rows,
        "n_fixed_site_pairs_checked": n_fixed_pairs_checked,
        "thresholds": threshold_rows,
        "n_any_threshold_components_AR_RA_AA_ge_1": sum(
            row["n_with_AR_RA_AA_each_ge_1"] for row in threshold_rows
        ),
        "status": (
            "ABSTAIN_NO_EVALUABLE_REAL_EDGE_UNIT"
            if not any(row["n_with_AR_RA_AA_each_ge_1"] for row in threshold_rows)
            else "STATE_ELIGIBLE_MM_ML_JOIN_STILL_REQUIRED"
        ),
        "claim_ceiling": (
            "Feasibility scan only; even a state-eligible component would still "
            "require raw MM/ML extraction, common-CpG, coverage and CN/LOH gates."
        ),
    }
    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": receipt["status"],
                "thresholds": threshold_rows,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
