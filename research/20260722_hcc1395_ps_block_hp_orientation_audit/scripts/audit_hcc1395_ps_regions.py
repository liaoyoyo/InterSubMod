#!/usr/bin/env python3
"""Audit HCC1395 retained regions for phase-set mixing.

This script is read-only with respect to canonical inputs.  It joins the
canonical MLHP region records to the current-v5 topology sidecar by the exact
``chrom:start-end`` key and writes a machine-readable census plus a region TSV.
It does not infer cross-PS orientation and therefore does not relabel HP tags.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


def percentage(numerator: int, denominator: int) -> float | None:
    return None if denominator == 0 else round(100.0 * numerator / denominator, 4)


def region_key(record: dict[str, Any]) -> str:
    return f"{record['chrom']}:{record['start']}-{record['end']}"


def hp_family(raw_hp: str) -> str:
    if raw_hp.startswith("1"):
        return "1"
    if raw_hp.startswith("2"):
        return "2"
    if raw_hp in {"3", "4"}:
        return raw_hp
    return "none"


def nested_counter() -> defaultdict[str, Counter[str]]:
    return defaultdict(Counter)


def counter_json(counter: Counter[Any]) -> dict[str, int]:
    return {str(key): int(value) for key, value in sorted(counter.items(), key=lambda item: str(item[0]))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mlhp-glob", required=True)
    parser.add_argument("--current-topology", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()

    mlhp_paths = [Path(path) for path in sorted(glob.glob(args.mlhp_glob))]
    if not mlhp_paths:
        raise FileNotFoundError(args.mlhp_glob)

    groups: list[dict[str, Any]] = []
    aggregate_phase_set_counts = Counter()
    aggregate_tag_counts = Counter()
    expected_group_count = 0
    for path in mlhp_paths:
        with path.open(encoding="utf-8") as handle:
            document = json.load(handle)
        if document.get("sample") != "HCC1395":
            raise ValueError(f"unexpected sample in {path}: {document.get('sample')}")
        groups.extend(document.get("groups", []))
        expected_group_count += int(document.get("n_groups_analyzed", 0))
        aggregate_phase_set_counts.update(
            document.get("read_tag_census", {}).get("phase_set_region_counts", {})
        )
        for field in (
            "alignment_group_exposures",
            "sidecar_exact_matches",
            "sidecar_missing",
            "sidecar_conflicts",
            "alignment_identity_allele_conflicts",
        ):
            aggregate_tag_counts[field] += int(
                document.get("read_tag_census", {}).get(field, 0)
            )

    with args.current_topology.open(encoding="utf-8") as handle:
        current = json.load(handle)
    if current.get("sample") != "HCC1395":
        raise ValueError(f"unexpected current topology sample: {current.get('sample')}")
    current_by_region = {record["region"]: record for record in current.get("regions", [])}
    if len(current_by_region) != len(current.get("regions", [])):
        raise ValueError("duplicate current-v5 region key")

    group_by_region: dict[str, dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []
    ps_distribution = Counter()
    structural_by_ps_class = nested_counter()
    morphology_by_ps_class = nested_counter()
    selection_by_ps_class = nested_counter()
    capped_by_ps_class = nested_counter()
    k_by_ps_class = nested_counter()
    hp_multiplicity_by_ps_class = nested_counter()
    mixed_thresholds = Counter()
    mixed_family_exposure = Counter()

    for group in groups:
        key = region_key(group)
        if key in group_by_region:
            raise ValueError(f"duplicate MLHP region key: {key}")
        group_by_region[key] = group
        topology = current_by_region.get(key)
        if topology is None:
            raise ValueError(f"MLHP region absent from current-v5: {key}")

        ps_counts = {str(k): int(v) for k, v in (group.get("phase_set_counts") or {}).items()}
        n_ps = len(ps_counts)
        if n_ps != int(group.get("n_unique_phase_sets", 0)):
            raise ValueError(f"n_unique_phase_sets mismatch: {key}")
        ps_class = "none" if n_ps == 0 else "one" if n_ps == 1 else "multiple"
        ps_distribution[n_ps] += 1
        structural_by_ps_class[ps_class][topology.get("structural_class", "missing")] += 1
        morphology_by_ps_class[ps_class][topology.get("morphology_class", "missing")] += 1
        selection_by_ps_class[ps_class][topology.get("selection_class", "missing")] += 1
        capped = int(group.get("pre_cap_n_sSNV", group.get("n_sSNV", 0))) > int(group.get("n_sSNV", 0))
        capped_by_ps_class[ps_class][str(capped).lower()] += 1
        k_by_ps_class[ps_class][str(group.get("n_sSNV"))] += 1
        hp_multiplicity = len(topology.get("primary_families", []))
        hp_multiplicity_by_ps_class[ps_class][str(hp_multiplicity)] += 1

        ordered_ps = sorted(ps_counts.items(), key=lambda item: (-item[1], item[0]))
        total_ps_reads = sum(ps_counts.values())
        second_ps_reads = ordered_ps[1][1] if len(ordered_ps) > 1 else 0
        second_ps_share = 0.0 if total_ps_reads == 0 else second_ps_reads / total_ps_reads

        phase_family_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
        for token, count in (group.get("phase_set_HP_counts") or {}).items():
            ps, raw_hp = str(token).split("|", 1)
            phase_family_counts[ps][hp_family(raw_hp)] += int(count)
        known_ps_with_hp12 = sum(
            (counts.get("1", 0) + counts.get("2", 0)) > 0
            for counts in phase_family_counts.values()
        )
        ps_with_hp1 = sum(counts.get("1", 0) > 0 for counts in phase_family_counts.values())
        ps_with_hp2 = sum(counts.get("2", 0) > 0 for counts in phase_family_counts.values())
        ps_with_both_hp12 = sum(
            counts.get("1", 0) > 0 and counts.get("2", 0) > 0
            for counts in phase_family_counts.values()
        )

        if n_ps > 1:
            mixed_thresholds["all"] += 1
            mixed_thresholds["second_ps_reads_ge_3"] += second_ps_reads >= 3
            mixed_thresholds["second_ps_reads_ge_10"] += second_ps_reads >= 10
            mixed_thresholds["second_ps_share_ge_0.10"] += second_ps_share >= 0.10
            mixed_thresholds["at_least_two_ps_with_hp12_reads"] += known_ps_with_hp12 >= 2
            mixed_thresholds["current_primary"] += bool(topology.get("primary_families"))
            mixed_thresholds["current_complete_primary"] += topology.get("structural_class") in {
                "exact_and_topology_unique",
                "topology_unique_exact_multiple",
                "topology_multiple_exact_multiple",
            }
            mixed_thresholds["current_incomplete"] += topology.get("structural_class") == "incomplete"
            mixed_thresholds["current_no_primary"] += topology.get("structural_class") == "no_primary_lineage"
            mixed_family_exposure["hp1_crosses_ps"] += ps_with_hp1 >= 2
            mixed_family_exposure["hp2_crosses_ps"] += ps_with_hp2 >= 2
            mixed_family_exposure["either_hp_family_crosses_ps"] += ps_with_hp1 >= 2 or ps_with_hp2 >= 2
            mixed_family_exposure["both_hp_families_cross_ps"] += ps_with_hp1 >= 2 and ps_with_hp2 >= 2
            mixed_family_exposure["at_least_two_ps_each_contain_both_hp_families"] += ps_with_both_hp12 >= 2

        rows.append(
            {
                "region": key,
                "chrom": group["chrom"],
                "start": int(group["start"]),
                "end": int(group["end"]),
                "span": int(group.get("span", int(group["end"]) - int(group["start"]))),
                "k": int(group.get("n_sSNV", 0)),
                "pre_cap_k": int(group.get("pre_cap_n_sSNV", group.get("n_sSNV", 0))),
                "capped": capped,
                "n_ps": n_ps,
                "ps_class": ps_class,
                "total_ps_reads": total_ps_reads,
                "second_ps_reads": second_ps_reads,
                "second_ps_share": round(second_ps_share, 6),
                "known_ps_with_hp12": known_ps_with_hp12,
                "ps_with_hp1": ps_with_hp1,
                "ps_with_hp2": ps_with_hp2,
                "ps_with_both_hp12": ps_with_both_hp12,
                "n_alignment_exposures": int(
                    group.get("read_tag_diagnostics", {}).get("alignment_group_exposures", 0)
                ),
                "n_primary_families": hp_multiplicity,
                "primary_families": ",".join(topology.get("primary_families", [])),
                "structural_class": topology.get("structural_class"),
                "selection_class": topology.get("selection_class"),
                "morphology_class": topology.get("morphology_class"),
                "T_region": topology.get("C"),
                "Topo_region": topology.get("Topo"),
                "phase_set_counts_json": json.dumps(ps_counts, sort_keys=True, separators=(",", ":")),
                "phase_set_hp_counts_json": json.dumps(
                    group.get("phase_set_HP_counts") or {}, sort_keys=True, separators=(",", ":")
                ),
            }
        )

    missing_in_mlhp = sorted(set(current_by_region) - set(group_by_region))
    rows.sort(key=lambda row: (int(row["chrom"][3:]), row["start"], row["end"]))
    mixed_rows = [row for row in rows if row["ps_class"] == "multiple"]
    representative_rows = sorted(
        (
            row
            for row in mixed_rows
            if not row["capped"]
            and row["structural_class"] != "incomplete"
            and row["n_primary_families"] > 0
            and row["second_ps_reads"] >= 3
        ),
        key=lambda row: (
            row["morphology_class"] != "direct_and_sister",
            -row["second_ps_reads"],
            -row["second_ps_share"],
            row["region"],
        ),
    )[:20]

    checks = {
        "mlhp_part_count_is_5": len(mlhp_paths) == 5,
        "mlhp_declared_group_sum_matches_loaded": expected_group_count == len(groups),
        "current_W_tree_matches_loaded": int(current["summary"]["W_tree"]) == len(groups),
        "region_key_sets_match": not missing_in_mlhp and len(current_by_region) == len(group_by_region),
        "aggregate_ps_class_matches_recount": {
            label: int(aggregate_phase_set_counts.get(label, 0))
            == sum(row["ps_class"] == label for row in rows)
            for label in ("none", "one", "multiple")
        },
        "mixed_partition_conserves": (
            mixed_thresholds["current_complete_primary"]
            + mixed_thresholds["current_incomplete"]
            + mixed_thresholds["current_no_primary"]
            == mixed_thresholds["all"]
        ),
        "exact_sidecar_join": (
            aggregate_tag_counts["alignment_group_exposures"]
            == aggregate_tag_counts["sidecar_exact_matches"]
            and aggregate_tag_counts["sidecar_missing"] == 0
            and aggregate_tag_counts["sidecar_conflicts"] == 0
            and aggregate_tag_counts["alignment_identity_allele_conflicts"] == 0
        ),
    }
    all_checks_pass = all(
        value if isinstance(value, bool) else all(value.values()) for value in checks.values()
    )

    output = {
        "schema_name": "intersubmod.hcc1395_ps_region_audit",
        "schema_version": "1.0.0",
        "scope": {
            "sample": "HCC1395",
            "chromosomes": "chr1-chr22",
            "unit": "retained MLHP region",
            "mixed_ps_definition": "n_unique_phase_sets > 1 among reads retained in the region",
            "interpretation_boundary": (
                "phase-set counts establish exposure to independent phase blocks; they do not establish "
                "the correct cross-block orientation or prove that topology changed"
            ),
        },
        "inputs": {
            "mlhp_parts": [str(path.resolve()) for path in mlhp_paths],
            "current_topology": str(args.current_topology.resolve()),
        },
        "counts": {
            "W_tree": len(rows),
            "W_primary": sum(bool(row["n_primary_families"]) for row in rows),
            "ps_class": counter_json(Counter(row["ps_class"] for row in rows)),
            "n_ps_distribution": counter_json(ps_distribution),
            "orientation_configuration_space": {
                "including_native": sum(
                    count * (2 ** (int(n_ps) - 1))
                    for n_ps, count in ps_distribution.items()
                    if int(n_ps) > 1
                ),
                "alternative_flips_only": sum(
                    count * ((2 ** (int(n_ps) - 1)) - 1)
                    for n_ps, count in ps_distribution.items()
                    if int(n_ps) > 1
                ),
                "global_swap_removed": True,
            },
            "mixed_ps_thresholds": counter_json(mixed_thresholds),
            "mixed_ps_family_exposure": counter_json(mixed_family_exposure),
            "mixed_ps_secondary_support": {
                "second_ps_reads_median": statistics.median(
                    [row["second_ps_reads"] for row in mixed_rows]
                ),
                "second_ps_share_median": round(
                    statistics.median([row["second_ps_share"] for row in mixed_rows]), 6
                ),
            },
            "read_tag_join": counter_json(aggregate_tag_counts),
        },
        "percentages": {
            "mixed_over_W_tree": percentage(len(mixed_rows), len(rows)),
            "mixed_primary_over_W_primary": percentage(
                mixed_thresholds["current_primary"],
                sum(bool(row["n_primary_families"]) for row in rows),
            ),
            "mixed_complete_primary_over_mixed": percentage(
                mixed_thresholds["current_complete_primary"], mixed_thresholds["all"]
            ),
            "mixed_complete_primary_over_mixed_primary": percentage(
                mixed_thresholds["current_complete_primary"], mixed_thresholds["current_primary"]
            ),
            "mixed_second_ps_ge3_over_mixed": percentage(
                mixed_thresholds["second_ps_reads_ge_3"], mixed_thresholds["all"]
            ),
            "mixed_second_share_ge10pct_over_mixed": percentage(
                mixed_thresholds["second_ps_share_ge_0.10"], mixed_thresholds["all"]
            ),
            "mixed_any_hp_family_crosses_ps_over_mixed": percentage(
                mixed_family_exposure["either_hp_family_crosses_ps"], mixed_thresholds["all"]
            ),
        },
        "cross_tabs": {
            "structural_by_ps_class": {
                key: counter_json(value) for key, value in structural_by_ps_class.items()
            },
            "morphology_by_ps_class": {
                key: counter_json(value) for key, value in morphology_by_ps_class.items()
            },
            "selection_by_ps_class": {
                key: counter_json(value) for key, value in selection_by_ps_class.items()
            },
            "capped_by_ps_class": {
                key: counter_json(value) for key, value in capped_by_ps_class.items()
            },
            "k_by_ps_class": {key: counter_json(value) for key, value in k_by_ps_class.items()},
            "hp_multiplicity_by_ps_class": {
                key: counter_json(value) for key, value in hp_multiplicity_by_ps_class.items()
            },
        },
        "representative_candidates": representative_rows,
        "checks": checks,
        "all_checks_pass": all_checks_pass,
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    fieldnames = list(rows[0]) if rows else []
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)

    print(json.dumps({"output_json": str(args.output_json), "output_tsv": str(args.output_tsv), **output["counts"], "all_checks_pass": all_checks_pass}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
