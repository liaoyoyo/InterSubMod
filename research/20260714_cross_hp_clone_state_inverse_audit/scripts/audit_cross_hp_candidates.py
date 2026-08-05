#!/usr/bin/env python3
"""Audit cross-HP direct+sister patterns and conditional clone-mixture readiness."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "1.1.0"
MINREAD_DEFAULT = 3
COUNT_FIELDS = (
    "regions_total",
    "regions_two_primary_hp",
    "regions_pair_complete",
    "regions_two_primary_hp_single_ps",
    "regions_two_primary_hp_two_site",
    "regions_two_primary_hp_any_cross_hp_same_site_alt",
    "regions_two_primary_hp_all_sites_cross_hp_same_site_alt",
    "regions_two_primary_hp_two_site_all_sites_cross_hp_same_site_alt",
    "regions_direct_sister_shape_invariant",
    "regions_direct_sister_shape_invariant_pair_complete",
    "regions_direct_sister_tree_unique",
    "regions_direct_sister_tree_unique_pair_complete",
    "regions_direct_sister_shape_invariant_any_collision",
    "regions_direct_sister_shape_invariant_pair_complete_any_collision",
    "regions_direct_sister_tree_unique_any_collision",
    "regions_direct_sister_tree_unique_pair_complete_any_collision",
    "regions_direct_sister_shape_invariant_all_sites_collision",
    "regions_direct_sister_tree_unique_all_sites_collision",
    "regions_observed_hp1_direct_hp2_sister_broad",
    "regions_observed_hp1_direct_only_hp2_sister_only",
    "regions_observed_either_direct_only_sister_only",
    "regions_observed_hp1_direct_only_hp2_sister_only_pair_complete_single_ps",
    "regions_observed_hp1_direct_only_hp2_sister_only_pair_complete_single_ps_cn_neutral_proxy",
    "regions_two_site_direct_sister_tree_unique",
    "regions_two_site_retained_catalog_exact",
    "regions_two_site_plugin_nonnegative",
    "regions_pattern_level_inverse_ready",
    "regions_two_primary_hp_cn_screened_neutral_proxy",
    "regions_cn_ps_screened_probe",
)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def state_of(node: str, k: int) -> str:
    if node == "ROOT":
        return "R" * k
    return node[2:] if node.startswith("H_") else node


def classify_tree(tree: dict[str, Any], k: int) -> dict[str, Any]:
    edges = [(str(parent), str(child)) for parent, child in tree.get("edges", [])]
    adjacency: dict[str, list[str]] = defaultdict(list)
    indegree: Counter[str] = Counter()
    nodes = {"ROOT"}
    for parent, child in edges:
        adjacency[parent].append(child)
        indegree[child] += 1
        nodes.update((parent, child))

    depths = {"ROOT": 0}
    queue: deque[str] = deque(["ROOT"])
    while queue:
        parent = queue.popleft()
        for child in adjacency.get(parent, []):
            if child not in depths or depths[child] > depths[parent] + 1:
                depths[child] = depths[parent] + 1
                queue.append(child)

    connected = len(depths) == len(nodes)
    max_depth = max(depths.values(), default=0)
    branching_nodes = sorted(node for node, children in adjacency.items() if len(children) >= 2)
    has_sister = bool(branching_nodes)
    has_direct = max_depth >= 2
    if not connected:
        shape_class = "disconnected_or_invalid"
    elif has_sister and has_direct:
        shape_class = "sister_direct"
    elif has_sister:
        shape_class = "sister_only"
    elif has_direct:
        shape_class = "direct_only"
    else:
        shape_class = "single_only"

    mutations: list[dict[str, Any]] = []
    hamming_ok = True
    for parent, child in edges:
        parent_state = state_of(parent, k)
        child_state = state_of(child, k)
        if len(parent_state) != k or len(child_state) != k:
            hamming_ok = False
            continue
        diffs = [idx for idx, (a, b) in enumerate(zip(parent_state, child_state)) if a != b]
        if len(diffs) != 1 or parent_state[diffs[0]] != "R" or child_state[diffs[0]] != "A":
            hamming_ok = False
        mutations.append(
            {
                "parent": parent,
                "child": child,
                "parent_state": parent_state,
                "child_state": child_state,
                "changed_indices": diffs,
            }
        )

    return {
        "shape_class": shape_class,
        "max_depth": max_depth,
        "branching_nodes": branching_nodes,
        "root_outdegree": len(adjacency.get("ROOT", [])),
        "connected": connected,
        "unit_hamming_edges": hamming_ok,
        "mutations": mutations,
    }


def lineage_profile(lineage: dict[str, Any], k: int) -> dict[str, Any]:
    tree_profiles = [classify_tree(tree, k) for tree in lineage.get("trees", [])]
    shape_classes = sorted({profile["shape_class"] for profile in tree_profiles})
    return {
        "family": str(lineage.get("family")),
        "n_trees": int(lineage.get("n_trees", 0)),
        "n_trees_stored": int(lineage.get("n_trees_stored", 0)),
        "n_distinct_shapes_exact": lineage.get("n_distinct_shapes_exact"),
        "candidate_set_complete": bool(lineage.get("analysis_candidate_set_complete")),
        "verification_status": lineage.get("verification_status"),
        "capped": bool(lineage.get("capped")),
        "shape_classes": shape_classes,
        "shape_invariant": len(shape_classes) == 1,
        "tree_profiles": tree_profiles,
    }


def retained_full_counts(lineage: dict[str, Any]) -> dict[str, int]:
    raw = lineage.get("obs_populations") or {}
    return {str(state): int(count) for state, count in raw.items()}


def observed_full_state_geometry(lineage: dict[str, Any], k: int) -> dict[str, Any]:
    counts = retained_full_counts(lineage)
    mutation_sets = sorted(
        {
            tuple(idx for idx, allele in enumerate(state) if allele == "A")
            for state, count in counts.items()
            if count > 0 and len(state) == k and "X" not in state and "A" in state
        }
    )
    sets = [set(indices) for indices in mutation_sets]
    has_direct = any(a < b or b < a for idx, a in enumerate(sets) for b in sets[idx + 1 :])
    has_sister = any(
        bool(a - b) and bool(b - a) for idx, a in enumerate(sets) for b in sets[idx + 1 :]
    )
    if has_direct and has_sister:
        geometry_class = "sister_direct"
    elif has_direct:
        geometry_class = "direct_only"
    elif has_sister:
        geometry_class = "sister_only"
    else:
        geometry_class = "single_only"
    return {
        "mutation_state_count": len(mutation_sets),
        "mutation_sets": mutation_sets,
        "has_direct": has_direct,
        "has_sister": has_sister,
        "geometry_class": geometry_class,
        "source": "retained full-coverage states only; each state count>=MINREAD upstream",
    }


def alt_count(lineage: dict[str, Any], position: int) -> int:
    counts = (lineage.get("obs_col_coverage") or {}).get(str(position), [0, 0])
    return int(counts[1]) if len(counts) >= 2 else 0


def direct_intermediate(tree: dict[str, Any], k: int) -> str | None:
    profile = classify_tree(tree, k)
    if profile["shape_class"] != "direct_only":
        return None
    root_children = [child for parent, child in tree.get("edges", []) if parent == "ROOT"]
    if len(root_children) != 1:
        return None
    state = state_of(str(root_children[0]), k)
    return state if state.count("A") == 1 else None


def normalized(counts: dict[str, int], states: Iterable[str]) -> dict[str, float] | None:
    wanted = list(states)
    total = sum(counts.get(state, 0) for state in wanted)
    if total <= 0:
        return None
    return {state: counts.get(state, 0) / total for state in wanted}


def infer_two_site_catalog(
    direct: dict[str, Any], sister: dict[str, Any], minread: int
) -> dict[str, Any]:
    """Return a threshold-censored plug-in estimate for the fixed five-state catalog."""
    if int(direct.get("n_trees", 0)) != 1 or int(sister.get("n_trees", 0)) != 1:
        return {"evaluable": False, "reason": "tree_not_unique"}
    direct_trees = direct.get("trees", [])
    sister_trees = sister.get("trees", [])
    if len(direct_trees) != 1 or len(sister_trees) != 1:
        return {"evaluable": False, "reason": "stored_tree_count_mismatch"}

    intermediate = direct_intermediate(direct_trees[0], 2)
    sister_profile = classify_tree(sister_trees[0], 2)
    if intermediate not in {"AR", "RA"} or sister_profile["shape_class"] != "sister_only":
        return {"evaluable": False, "reason": "not_direct_plus_sister"}
    other_singleton = "RA" if intermediate == "AR" else "AR"
    direct_counts = retained_full_counts(direct)
    sister_counts = retained_full_counts(sister)
    required_direct = {"RR", intermediate, "AA"}
    required_sister = {"RR", intermediate, other_singleton}
    observed_direct = {state for state, count in direct_counts.items() if count >= minread}
    observed_sister = {state for state, count in sister_counts.items() if count >= minread}
    missing_direct = sorted(required_direct - observed_direct)
    missing_sister = sorted(required_sister - observed_sister)
    unexpected_direct = sorted(observed_direct - required_direct)
    unexpected_sister = sorted(observed_sister - required_sister)
    catalog_exact = not (missing_direct or missing_sister or unexpected_direct or unexpected_sister)

    result: dict[str, Any] = {
        "evaluable": catalog_exact,
        "reason": "ok" if catalog_exact else "retained_full_state_catalog_mismatch",
        "intermediate_state": intermediate,
        "other_singleton_state": other_singleton,
        "direct_counts": direct_counts,
        "sister_counts": sister_counts,
        "missing_direct": missing_direct,
        "missing_sister": missing_sister,
        "unexpected_direct": unexpected_direct,
        "unexpected_sister": unexpected_sister,
        "catalog_exact_under_minread": catalog_exact,
        "count_semantics": (
            "full-coverage state counts retained only when count>=MINREAD; "
            "obs_subreads are not added because they also contain full reads"
        ),
    }
    if not catalog_exact:
        return result

    p_direct = normalized(direct_counts, ["RR", intermediate, "AA"])
    p_sister = normalized(sister_counts, ["RR", intermediate, other_singleton])
    assert p_direct is not None and p_sister is not None
    estimates = {
        "background_pi0": p_direct["RR"] - p_sister[intermediate],
        "clone1_pi1": p_direct[intermediate] - p_sister[other_singleton],
        "clone2_pi2": p_direct["AA"],
        "clone3_pi3": p_sister[intermediate],
        "clone4_pi4": p_sister[other_singleton],
    }
    tolerance = 1e-12
    compatible = all(value >= -tolerance for value in estimates.values()) and math.isclose(
        sum(estimates.values()), 1.0, abs_tol=1e-9
    )
    result.update(
        {
            "p_direct": p_direct,
            "p_sister": p_sister,
            "plugin_estimates": estimates,
            "plugin_sum": sum(estimates.values()),
            "plugin_nonnegative": compatible,
            "interpretation": (
                "conditional proxy under fixed five-state catalog, equal within-HP exposure, "
                "no HP error, no allelic CN imbalance, and threshold censoring ignored"
            ),
        }
    )
    return result


def load_group_index(sample_dir: Path) -> tuple[dict[tuple[Any, ...], dict[str, Any]], dict[str, Any]]:
    index: dict[tuple[Any, ...], dict[str, Any]] = {}
    diagnostics = {"part_files": 0, "groups": 0, "duplicate_keys": 0, "minread_values": []}
    for path in sorted(sample_dir.glob("mlhp_part_*.json")):
        payload = load_json(path)
        diagnostics["part_files"] += 1
        diagnostics["minread_values"].append(int(payload.get("params", {}).get("MINREAD", MINREAD_DEFAULT)))
        for group in payload.get("groups", []):
            key = (
                str(group["chrom"]),
                int(group["start"]),
                int(group["end"]),
                tuple(int(pos) for pos in group.get("positions", [])),
            )
            if key in index:
                diagnostics["duplicate_keys"] += 1
            index[key] = group
            diagnostics["groups"] += 1
    diagnostics["minread_values"] = sorted(set(diagnostics["minread_values"]))
    return index, diagnostics


def region_key(region: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(region["chrom"]),
        int(region["start"]),
        int(region["end"]),
        tuple(int(pos) for pos in region.get("positions", [])),
    )


def analyze_sample(sample_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    sample = sample_dir.name
    region_view_path = sample_dir / f"layered_region_view_{sample}.json"
    manifest_path = sample_dir / "output_manifest.json"
    region_view = load_json(region_view_path)
    manifest = load_json(manifest_path)
    group_index, group_diagnostics = load_group_index(sample_dir)
    minread_values = group_diagnostics["minread_values"]
    if len(minread_values) != 1:
        raise ValueError(f"{sample}: inconsistent MINREAD values: {minread_values}")
    minread = minread_values[0]
    copy_contract = region_view.get("copy_number_contract", {})
    cn_measured = copy_contract.get("availability") == "measured"

    counts: Counter[str] = Counter()
    records: list[dict[str, Any]] = []
    missing_group_join = 0
    for region in region_view.get("regions", []):
        counts["regions_total"] += 1
        key = region_key(region)
        group = group_index.get(key)
        if group is None:
            missing_group_join += 1
        primary = {
            str(lineage.get("family")): lineage
            for lineage in region.get("lineages", [])
            if lineage.get("is_primary_lineage") and lineage.get("mutation_bearing")
        }
        if "1" not in primary or "2" not in primary:
            continue
        counts["regions_two_primary_hp"] += 1
        k = int(region.get("n_sSNV", len(region.get("positions", []))))
        p1 = lineage_profile(primary["1"], k)
        p2 = lineage_profile(primary["2"], k)
        observed_geometry_hp1 = observed_full_state_geometry(primary["1"], k)
        observed_geometry_hp2 = observed_full_state_geometry(primary["2"], k)
        pair_complete = all(
            profile["candidate_set_complete"]
            and not profile["capped"]
            and profile["verification_status"] == "full_pass"
            for profile in (p1, p2)
        )
        if pair_complete:
            counts["regions_pair_complete"] += 1
        ps_unique = int(group.get("n_unique_phase_sets", -1)) if group else -1
        if ps_unique == 1:
            counts["regions_two_primary_hp_single_ps"] += 1
        if k == 2:
            counts["regions_two_primary_hp_two_site"] += 1
        positions = [int(pos) for pos in region.get("positions", [])]
        collision_positions = [
            pos
            for pos in positions
            if alt_count(primary["1"], pos) >= minread and alt_count(primary["2"], pos) >= minread
        ]
        if collision_positions:
            counts["regions_two_primary_hp_any_cross_hp_same_site_alt"] += 1
        if len(collision_positions) == k and k > 0:
            counts["regions_two_primary_hp_all_sites_cross_hp_same_site_alt"] += 1
            if k == 2:
                counts["regions_two_primary_hp_two_site_all_sites_cross_hp_same_site_alt"] += 1

        shape_pair = {tuple(p1["shape_classes"]), tuple(p2["shape_classes"])}
        direct_sister_invariant = shape_pair == {("direct_only",), ("sister_only",)}
        tree_unique_pair = direct_sister_invariant and p1["n_trees"] == 1 and p2["n_trees"] == 1
        direct_sister_invariant_complete = direct_sister_invariant and pair_complete
        tree_unique_pair_complete = tree_unique_pair and pair_complete
        if direct_sister_invariant:
            counts["regions_direct_sister_shape_invariant"] += 1
        if direct_sister_invariant_complete:
            counts["regions_direct_sister_shape_invariant_pair_complete"] += 1
        if tree_unique_pair:
            counts["regions_direct_sister_tree_unique"] += 1
        if tree_unique_pair_complete:
            counts["regions_direct_sister_tree_unique_pair_complete"] += 1
        if direct_sister_invariant and collision_positions:
            counts["regions_direct_sister_shape_invariant_any_collision"] += 1
        if direct_sister_invariant_complete and collision_positions:
            counts["regions_direct_sister_shape_invariant_pair_complete_any_collision"] += 1
        if tree_unique_pair and collision_positions:
            counts["regions_direct_sister_tree_unique_any_collision"] += 1
        if tree_unique_pair_complete and collision_positions:
            counts["regions_direct_sister_tree_unique_pair_complete_any_collision"] += 1
        if direct_sister_invariant and len(collision_positions) == k and k > 0:
            counts["regions_direct_sister_shape_invariant_all_sites_collision"] += 1
        if tree_unique_pair and len(collision_positions) == k and k > 0:
            counts["regions_direct_sister_tree_unique_all_sites_collision"] += 1

        observed_hp1_direct_hp2_sister_broad = bool(
            observed_geometry_hp1["has_direct"] and observed_geometry_hp2["has_sister"]
        )
        observed_hp1_direct_only_hp2_sister_only = bool(
            observed_geometry_hp1["geometry_class"] == "direct_only"
            and observed_geometry_hp2["geometry_class"] == "sister_only"
        )
        observed_either_direct_only_sister_only = {
            observed_geometry_hp1["geometry_class"], observed_geometry_hp2["geometry_class"]
        } == {"direct_only", "sister_only"}
        if observed_hp1_direct_hp2_sister_broad:
            counts["regions_observed_hp1_direct_hp2_sister_broad"] += 1
        if observed_hp1_direct_only_hp2_sister_only:
            counts["regions_observed_hp1_direct_only_hp2_sister_only"] += 1
        if observed_either_direct_only_sister_only:
            counts["regions_observed_either_direct_only_sister_only"] += 1

        direct_family: str | None = None
        sister_family: str | None = None
        if tree_unique_pair:
            direct_family = "1" if p1["shape_classes"] == ["direct_only"] else "2"
            sister_family = "2" if direct_family == "1" else "1"
        catalog = {"evaluable": False, "reason": "not_two_site_tree_unique_direct_sister"}
        if k == 2 and direct_family and sister_family and pair_complete:
            counts["regions_two_site_direct_sister_tree_unique"] += 1
            catalog = infer_two_site_catalog(primary[direct_family], primary[sister_family], minread)
            if catalog.get("catalog_exact_under_minread"):
                counts["regions_two_site_retained_catalog_exact"] += 1
            if catalog.get("plugin_nonnegative"):
                counts["regions_two_site_plugin_nonnegative"] += 1

        cn_screened_neutral = cn_measured and region.get("cn") == "neutral"
        if cn_screened_neutral:
            counts["regions_two_primary_hp_cn_screened_neutral_proxy"] += 1
        if observed_hp1_direct_only_hp2_sister_only and pair_complete and ps_unique == 1:
            counts["regions_observed_hp1_direct_only_hp2_sister_only_pair_complete_single_ps"] += 1
            if cn_screened_neutral:
                counts[
                    "regions_observed_hp1_direct_only_hp2_sister_only_pair_complete_single_ps_cn_neutral_proxy"
                ] += 1
        pattern_ready = bool(
            k == 2
            and tree_unique_pair
            and pair_complete
            and len(collision_positions) == 2
            and catalog.get("catalog_exact_under_minread")
            and catalog.get("plugin_nonnegative")
        )
        if pattern_ready:
            counts["regions_pattern_level_inverse_ready"] += 1
        cn_ps_screened_probe = pattern_ready and ps_unique == 1 and cn_screened_neutral
        if cn_ps_screened_probe:
            counts["regions_cn_ps_screened_probe"] += 1

        # Keep every direct+sister candidate and every two-site collision candidate for audit.
        if (
            direct_sister_invariant
            or observed_hp1_direct_only_hp2_sister_only
            or (k == 2 and len(collision_positions) == 2)
        ):
            records.append(
                {
                    "dataset": sample,
                    "biological_id": manifest.get("biological_id", sample),
                    "region": region.get("region"),
                    "chrom": region.get("chrom"),
                    "start": int(region.get("start")),
                    "end": int(region.get("end")),
                    "positions": positions,
                    "n_sSNV": k,
                    "cn_label": region.get("cn"),
                    "cn_contract_availability": copy_contract.get("availability"),
                    "cn_contract_source": copy_contract.get("source"),
                    "cn_screened_neutral_proxy": cn_screened_neutral,
                    "allele_specific_cn_1plus1_verified": False,
                    "n_unique_phase_sets": ps_unique,
                    "single_ps": ps_unique == 1,
                    "collision_positions": collision_positions,
                    "all_sites_collision": len(collision_positions) == k and k > 0,
                    "pair_complete": pair_complete,
                    "hp1_shape_classes": p1["shape_classes"],
                    "hp2_shape_classes": p2["shape_classes"],
                    "hp1_n_trees": p1["n_trees"],
                    "hp2_n_trees": p2["n_trees"],
                    "direct_sister_shape_invariant": direct_sister_invariant,
                    "direct_sister_shape_invariant_analysis_complete": (
                        direct_sister_invariant_complete
                    ),
                    "direct_sister_tree_unique": tree_unique_pair,
                    "direct_sister_tree_unique_analysis_complete": tree_unique_pair_complete,
                    "observed_full_geometry_hp1": observed_geometry_hp1,
                    "observed_full_geometry_hp2": observed_geometry_hp2,
                    "observed_hp1_direct_hp2_sister_broad": observed_hp1_direct_hp2_sister_broad,
                    "observed_hp1_direct_only_hp2_sister_only": (
                        observed_hp1_direct_only_hp2_sister_only
                    ),
                    "observed_either_direct_only_sister_only": observed_either_direct_only_sister_only,
                    "direct_family": direct_family,
                    "sister_family": sister_family,
                    "hp1_full_counts_retained": retained_full_counts(primary["1"]),
                    "hp2_full_counts_retained": retained_full_counts(primary["2"]),
                    "hp1_subread_counts_retained": primary["1"].get("obs_subreads") or {},
                    "hp2_subread_counts_retained": primary["2"].get("obs_subreads") or {},
                    "catalog_inference": catalog,
                    "pattern_level_inverse_ready": pattern_ready,
                    "cn_ps_screened_probe": cn_ps_screened_probe,
                    "biological_clone_confirmed": False,
                    "claim_ceiling": (
                        "candidate cross-HP mutation-state pattern; not a confirmed clone pairing or lineage"
                    ),
                }
            )

    if missing_group_join:
        raise ValueError(f"{sample}: {missing_group_join} region-view rows lacked an mlhp group join")
    if group_diagnostics["duplicate_keys"]:
        raise ValueError(f"{sample}: duplicate group keys={group_diagnostics['duplicate_keys']}")
    if counts["regions_total"] != group_diagnostics["groups"]:
        raise ValueError(
            f"{sample}: region conservation failed: view={counts['regions_total']} "
            f"mlhp={group_diagnostics['groups']}"
        )

    summary = {
        "dataset": sample,
        "biological_id": manifest.get("biological_id", sample),
        "run_id": manifest.get("run_id"),
        "minread": minread,
        "copy_number_contract": {
            "availability": copy_contract.get("availability"),
            "source": copy_contract.get("source"),
            "semantics": copy_contract.get("semantics"),
            "allele_specific_cn_1plus1_available": False,
        },
        "input_files": {
            "region_view": str(region_view_path.resolve()),
            "region_view_sha256": sha256_file(region_view_path),
            "output_manifest": str(manifest_path.resolve()),
            "output_manifest_sha256": sha256_file(manifest_path),
            "mlhp_part_count": group_diagnostics["part_files"],
        },
        "join_checks": {
            "region_view_rows": counts["regions_total"],
            "mlhp_groups": group_diagnostics["groups"],
            "missing_group_join": missing_group_join,
            "duplicate_group_keys": group_diagnostics["duplicate_keys"],
            "conservation_pass": counts["regions_total"] == group_diagnostics["groups"],
        },
        "counts": {field: int(counts.get(field, 0)) for field in COUNT_FIELDS},
    }
    return summary, records


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: json.dumps(row.get(field), ensure_ascii=False, sort_keys=True)
                    if isinstance(row.get(field), (dict, list))
                    else row.get(field)
                    for field in fields
                }
            )


def run_audit(run_dir: Path, out_dir: Path) -> dict[str, Any]:
    sample_root = run_dir / "samples"
    sample_dirs = sorted(path for path in sample_root.iterdir() if path.is_dir())
    if not sample_dirs:
        raise ValueError(f"No sample directories found under {sample_root}")
    run_state = load_json(run_dir / "run_state.json")
    if run_state.get("state") != "SUCCEEDED":
        raise ValueError(f"Input run is not SUCCEEDED: {run_state.get('state')}")

    summaries: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []
    for sample_dir in sample_dirs:
        summary, sample_records = analyze_sample(sample_dir)
        summaries.append(summary)
        records.extend(sample_records)

    aggregate: Counter[str] = Counter()
    for summary in summaries:
        aggregate.update(summary["counts"])
    biological_ids = sorted({summary["biological_id"] for summary in summaries})
    strict_ready = sum(
        1
        for record in records
        if record["cn_ps_screened_probe"] and record["allele_specific_cn_1plus1_verified"]
    )
    payload = {
        "schema_name": "intersubmod.cross_hp_clone_state_inverse_audit",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "comprehensive_validation",
        "research_goals": ["G3", "G4", "G5"],
        "input_run": {
            "path": str(run_dir.resolve()),
            "run_id": run_state.get("run_id"),
            "state": run_state.get("state"),
            "run_state_sha256": sha256_file(run_dir / "run_state.json"),
            "scope": "chr1-22; all datasets in succeeded run",
        },
        "dataset_count": len(summaries),
        "biological_sample_count": len(biological_ids),
        "biological_ids": biological_ids,
        "definitions": {
            "cross_hp_same_site_alt": "same sSNV position has ALT count>=MINREAD in both HP1 and HP2",
            "direct_sister_shape_invariant": (
                "all currently enumerated/stored trees for one HP are direct_only and all for the "
                "other HP are sister_only; this is not exhaustive when pair_complete=false"
            ),
            "direct_sister_shape_invariant_analysis_complete": (
                "direct_sister_shape_invariant plus pair_complete=true for both HP candidate sets"
            ),
            "tree_unique": (
                "both HP outputs currently report exactly one tree; only analysis-complete when "
                "pair_complete=true"
            ),
            "tree_unique_analysis_complete": (
                "tree_unique plus pair_complete=true for both HP candidate sets"
            ),
            "pattern_level_inverse_ready": (
                "two sSNVs, unique direct+sister trees, both sites collide across HP, retained full-state "
                "catalog is exact, and the threshold-censored plug-in solution is nonnegative"
            ),
            "cn_ps_screened_probe": (
                "pattern-level ready plus one recorded PS and measured-contract neutral proxy; "
                "still not allele-specific CN=1+1"
            ),
            "strict_biological_inverse_ready": (
                "cn_ps_screened_probe plus allele-specific CN=1+1, HP-error calibration, uncensored "
                "per-molecule state table, and orthogonal validation"
            ),
        },
        "aggregate_counts": {field: int(aggregate.get(field, 0)) for field in COUNT_FIELDS},
        "strict_biological_inverse_ready_count": strict_ready,
        "strict_biological_inverse_ready_verdict": "0 by current canonical contract",
        "claim_ceiling": (
            "Current output can census candidate cross-HP mutation-state patterns. It cannot by itself "
            "confirm cell-level HP pairings, true clone identities, or a biological clone tree."
        ),
        "sample_summaries": summaries,
        "candidate_records": records,
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "cross_hp_candidate_audit.json"
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    sample_fields = [
        "dataset",
        "biological_id",
        "run_id",
        "minread",
        "regions_total",
        "regions_two_primary_hp",
        "regions_pair_complete",
        "regions_two_primary_hp_single_ps",
        "regions_two_primary_hp_two_site",
        "regions_two_primary_hp_any_cross_hp_same_site_alt",
        "regions_two_primary_hp_all_sites_cross_hp_same_site_alt",
        "regions_two_primary_hp_two_site_all_sites_cross_hp_same_site_alt",
        "regions_direct_sister_shape_invariant",
        "regions_direct_sister_shape_invariant_pair_complete",
        "regions_direct_sister_tree_unique",
        "regions_direct_sister_tree_unique_pair_complete",
        "regions_direct_sister_shape_invariant_any_collision",
        "regions_direct_sister_shape_invariant_pair_complete_any_collision",
        "regions_direct_sister_tree_unique_any_collision",
        "regions_direct_sister_tree_unique_pair_complete_any_collision",
        "regions_direct_sister_shape_invariant_all_sites_collision",
        "regions_direct_sister_tree_unique_all_sites_collision",
        "regions_observed_hp1_direct_hp2_sister_broad",
        "regions_observed_hp1_direct_only_hp2_sister_only",
        "regions_observed_either_direct_only_sister_only",
        "regions_observed_hp1_direct_only_hp2_sister_only_pair_complete_single_ps",
        "regions_observed_hp1_direct_only_hp2_sister_only_pair_complete_single_ps_cn_neutral_proxy",
        "regions_two_site_direct_sister_tree_unique",
        "regions_two_site_retained_catalog_exact",
        "regions_two_site_plugin_nonnegative",
        "regions_pattern_level_inverse_ready",
        "regions_two_primary_hp_cn_screened_neutral_proxy",
        "regions_cn_ps_screened_probe",
        "cn_contract_availability",
        "cn_contract_source",
        "allele_specific_cn_1plus1_available",
        "join_conservation_pass",
    ]
    sample_rows: list[dict[str, Any]] = []
    for summary in summaries:
        counts = summary["counts"]
        sample_rows.append(
            {
                "dataset": summary["dataset"],
                "biological_id": summary["biological_id"],
                "run_id": summary["run_id"],
                "minread": summary["minread"],
                **{field: counts.get(field, 0) for field in sample_fields if field.startswith("regions_")},
                "cn_contract_availability": summary["copy_number_contract"]["availability"],
                "cn_contract_source": summary["copy_number_contract"]["source"],
                "allele_specific_cn_1plus1_available": summary["copy_number_contract"][
                    "allele_specific_cn_1plus1_available"
                ],
                "join_conservation_pass": summary["join_checks"]["conservation_pass"],
            }
        )
    write_tsv(out_dir / "sample_summary.tsv", sample_rows, sample_fields)

    candidate_fields = [
        "dataset",
        "biological_id",
        "region",
        "chrom",
        "start",
        "end",
        "positions",
        "n_sSNV",
        "cn_label",
        "cn_contract_availability",
        "cn_screened_neutral_proxy",
        "allele_specific_cn_1plus1_verified",
        "n_unique_phase_sets",
        "single_ps",
        "collision_positions",
        "all_sites_collision",
        "pair_complete",
        "hp1_shape_classes",
        "hp2_shape_classes",
        "hp1_n_trees",
        "hp2_n_trees",
        "direct_sister_shape_invariant",
        "direct_sister_shape_invariant_analysis_complete",
        "direct_sister_tree_unique",
        "direct_sister_tree_unique_analysis_complete",
        "observed_full_geometry_hp1",
        "observed_full_geometry_hp2",
        "observed_hp1_direct_hp2_sister_broad",
        "observed_hp1_direct_only_hp2_sister_only",
        "observed_either_direct_only_sister_only",
        "direct_family",
        "sister_family",
        "hp1_full_counts_retained",
        "hp2_full_counts_retained",
        "catalog_inference",
        "pattern_level_inverse_ready",
        "cn_ps_screened_probe",
        "biological_clone_confirmed",
        "claim_ceiling",
    ]
    write_tsv(out_dir / "candidate_regions.tsv", records, candidate_fields)
    print(
        json.dumps(
            {
                "run_id": run_state.get("run_id"),
                "datasets": len(summaries),
                "biological_samples": len(biological_ids),
                "aggregate_counts": {field: int(aggregate.get(field, 0)) for field in COUNT_FIELDS},
                "candidate_records_written": len(records),
                "strict_biological_inverse_ready": strict_ready,
                "outputs": [
                    str(json_path.resolve()),
                    str((out_dir / "sample_summary.tsv").resolve()),
                    str((out_dir / "candidate_regions.tsv").resolve()),
                ],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_audit(args.run_dir, args.out_dir)
