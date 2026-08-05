#!/usr/bin/env python3
"""Quantify retained regional mutation-state bounds without claiming biological clone K."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "1.0.0"
PRIMARY_FAMILIES = ("1", "2")


def popcount(value: int) -> int:
    """Python 3.9-compatible population count."""
    return bin(value).count("1")


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pattern_conflict(left: str, right: str) -> bool:
    """Return whether two R/A/X partial assignments cannot share one full state."""
    if len(left) != len(right):
        raise ValueError(f"Pattern length mismatch: {left!r}, {right!r}")
    return any(a != "X" and b != "X" and a != b for a, b in zip(left, right))


def compatible(pattern: str, state: str) -> bool:
    if len(pattern) != len(state):
        return False
    return all(observed == "X" or observed == allele for observed, allele in zip(pattern, state))


def _greedy_dsatur(adjacency: list[int]) -> tuple[int, list[int]]:
    n = len(adjacency)
    colors = [-1] * n
    saturation: list[set[int]] = [set() for _ in range(n)]
    degrees = [popcount(mask) for mask in adjacency]
    used = 0
    for _ in range(n):
        vertex = max(
            (idx for idx, color in enumerate(colors) if color < 0),
            key=lambda idx: (len(saturation[idx]), degrees[idx], -idx),
        )
        forbidden = saturation[vertex]
        color = next(candidate for candidate in range(used + 1) if candidate not in forbidden)
        colors[vertex] = color
        used = max(used, color + 1)
        neighbors = adjacency[vertex]
        while neighbors:
            bit = neighbors & -neighbors
            neighbor = bit.bit_length() - 1
            if colors[neighbor] < 0:
                saturation[neighbor].add(color)
            neighbors ^= bit
    return used, colors


def _greedy_clique_lower_bound(adjacency: list[int]) -> int:
    """Return a valid, inexpensive clique lower bound for chromatic search."""
    if not adjacency:
        return 0
    candidates = (1 << len(adjacency)) - 1
    size = 0
    while candidates:
        vertices = [idx for idx in range(len(adjacency)) if candidates & (1 << idx)]
        vertex = max(vertices, key=lambda idx: popcount(adjacency[idx] & candidates))
        size += 1
        candidates &= adjacency[vertex]
    return size


def minimum_compatible_state_count(patterns: Iterable[str]) -> int:
    """Exact minimum number of full R/A states covering unique partial patterns.

    A color class in the conflict graph is pairwise compatible. For coordinate-wise
    R/A/X constraints, pairwise compatibility is also jointly compatible, so the
    exact chromatic number is the exact minimum state-cover size.
    """
    unique = sorted(set(patterns))
    n = len(unique)
    if n == 0:
        return 0
    adjacency = [0] * n
    for left in range(n):
        for right in range(left + 1, n):
            if pattern_conflict(unique[left], unique[right]):
                adjacency[left] |= 1 << right
                adjacency[right] |= 1 << left
    if not any(adjacency):
        return 1

    best, _ = _greedy_dsatur(adjacency)
    lower = _greedy_clique_lower_bound(adjacency)
    if best == lower:
        return best

    colors = [-1] * n
    saturation: list[set[int]] = [set() for _ in range(n)]
    degrees = [popcount(mask) for mask in adjacency]

    def search(colored: int, used: int) -> None:
        nonlocal best
        if used >= best:
            return
        if colored == n:
            best = used
            return
        vertex = max(
            (idx for idx, color in enumerate(colors) if color < 0),
            key=lambda idx: (len(saturation[idx]), degrees[idx], -idx),
        )
        forbidden = saturation[vertex]
        candidate_colors = [color for color in range(used) if color not in forbidden]
        if used + 1 < best:
            candidate_colors.append(used)
        for color in candidate_colors:
            colors[vertex] = color
            changed: list[int] = []
            neighbors = adjacency[vertex]
            while neighbors:
                bit = neighbors & -neighbors
                neighbor = bit.bit_length() - 1
                if colors[neighbor] < 0 and color not in saturation[neighbor]:
                    saturation[neighbor].add(color)
                    changed.append(neighbor)
                neighbors ^= bit
            search(colored + 1, max(used, color + 1))
            for neighbor in changed:
                saturation[neighbor].remove(color)
            colors[vertex] = -1
            if best == lower:
                return

    search(0, 0)
    return best


@dataclass
class _Edge:
    to: int
    rev: int
    cap: int


class _Dinic:
    def __init__(self, n_nodes: int) -> None:
        self.graph: list[list[_Edge]] = [[] for _ in range(n_nodes)]

    def add_edge(self, source: int, target: int, capacity: int) -> None:
        forward = _Edge(target, len(self.graph[target]), capacity)
        reverse = _Edge(source, len(self.graph[source]), 0)
        self.graph[source].append(forward)
        self.graph[target].append(reverse)

    def max_flow(self, source: int, sink: int) -> int:
        total = 0
        while True:
            level = [-1] * len(self.graph)
            level[source] = 0
            queue: deque[int] = deque([source])
            while queue:
                node = queue.popleft()
                for edge in self.graph[node]:
                    if edge.cap > 0 and level[edge.to] < 0:
                        level[edge.to] = level[node] + 1
                        queue.append(edge.to)
            if level[sink] < 0:
                return total
            cursor = [0] * len(self.graph)

            def send(node: int, flow: int) -> int:
                if node == sink:
                    return flow
                while cursor[node] < len(self.graph[node]):
                    edge = self.graph[node][cursor[node]]
                    if edge.cap > 0 and level[edge.to] == level[node] + 1:
                        pushed = send(edge.to, min(flow, edge.cap))
                        if pushed:
                            edge.cap -= pushed
                            self.graph[edge.to][edge.rev].cap += pushed
                            return pushed
                    cursor[node] += 1
                return 0

            while True:
                pushed = send(source, 10**9)
                if not pushed:
                    break
                total += pushed


def full_states(k: int, mutation_only: bool = False) -> list[str]:
    states = []
    for value in range(1 << k):
        state = "".join("A" if value & (1 << idx) else "R" for idx in range(k))
        if mutation_only and "A" not in state:
            continue
        states.append(state)
    return states


def maximum_supported_state_count(
    pattern_counts: dict[str, int], k: int, mutation_only: bool = False
) -> int:
    """Maximum distinct full states assignable to retained aggregate read counts."""
    patterns = [(pattern, count) for pattern, count in sorted(pattern_counts.items()) if count > 0]
    states = full_states(k, mutation_only=mutation_only)
    if not patterns or not states:
        return 0
    source = 0
    pattern_offset = 1
    state_offset = pattern_offset + len(patterns)
    sink = state_offset + len(states)
    flow = _Dinic(sink + 1)
    for idx, (_, count) in enumerate(patterns):
        flow.add_edge(source, pattern_offset + idx, int(count))
    for p_idx, (pattern, _) in enumerate(patterns):
        for s_idx, state in enumerate(states):
            if compatible(pattern, state):
                flow.add_edge(pattern_offset + p_idx, state_offset + s_idx, 1)
    for s_idx in range(len(states)):
        flow.add_edge(state_offset + s_idx, sink, 1)
    return flow.max_flow(source, sink)


def classify_tree_shape(tree: dict[str, Any]) -> str:
    adjacency: dict[str, list[str]] = defaultdict(list)
    depths = {"ROOT": 0}
    for parent, child in tree.get("edges", []):
        adjacency[str(parent)].append(str(child))
    queue: deque[str] = deque(["ROOT"])
    while queue:
        parent = queue.popleft()
        for child in adjacency.get(parent, []):
            if child not in depths or depths[child] > depths[parent] + 1:
                depths[child] = depths[parent] + 1
                queue.append(child)
    has_sister = any(len(children) >= 2 for children in adjacency.values())
    has_direct = max(depths.values(), default=0) >= 2
    if has_sister and has_direct:
        return "sister_direct"
    if has_sister:
        return "sister_only"
    if has_direct:
        return "direct_only"
    return "single_only"


def validate_pattern_counts(lineage: dict[str, Any], k: int) -> tuple[dict[str, int], int]:
    raw = lineage.get("obs_subreads") or {}
    valid: dict[str, int] = {}
    invalid = 0
    for pattern, count in raw.items():
        pattern = str(pattern)
        count = int(count)
        if len(pattern) != k or not set(pattern) <= {"R", "A", "X"} or count <= 0:
            invalid += 1
            continue
        valid[pattern] = valid.get(pattern, 0) + count
    return valid, invalid


def analyze_lineage(lineage: dict[str, Any], k: int) -> dict[str, Any]:
    counts, invalid = validate_pattern_counts(lineage, k)
    patterns = sorted(counts)
    mutation_patterns = [pattern for pattern in patterns if "A" in pattern]
    retained_reads = sum(counts.values())
    reported_reads = int(lineage.get("n_reads", retained_reads))
    censored_reads = max(0, reported_reads - retained_reads)
    all_full = bool(patterns) and all("X" not in pattern for pattern in patterns)
    full_counts = {pattern: count for pattern, count in counts.items() if "X" not in pattern}
    shapes = sorted({classify_tree_shape(tree) for tree in lineage.get("trees", [])})
    return {
        "family": str(lineage.get("family")),
        "retained_pattern_count": len(patterns),
        "retained_pattern_counts": counts,
        "retained_read_count": retained_reads,
        "reported_read_count": reported_reads,
        "minread_censored_read_count": censored_reads,
        "minread_censored_fraction": censored_reads / reported_reads if reported_reads else None,
        "invalid_pattern_count": invalid,
        "full_pattern_count": len(full_counts),
        "has_full_pattern": bool(full_counts),
        "all_retained_patterns_full": all_full,
        "k_state_min_retained": minimum_compatible_state_count(patterns),
        "k_mutation_state_min_retained": minimum_compatible_state_count(mutation_patterns),
        "k_state_max_retained": maximum_supported_state_count(counts, k),
        "k_mutation_state_max_retained": maximum_supported_state_count(
            counts, k, mutation_only=True
        ),
        "full_state_counts_retained": full_counts,
        "catalog_identified_under_retained_table": all_full,
        "n_trees": int(lineage.get("n_trees", 0)),
        "n_trees_stored": int(lineage.get("n_trees_stored", 0)),
        "shape_classes": shapes,
        "shape_invariant": len(shapes) == 1,
        "candidate_set_complete": bool(lineage.get("analysis_candidate_set_complete")),
        "verification_status": lineage.get("verification_status"),
        "capped": bool(lineage.get("capped")),
    }


def load_group_phase_sets(sample_dir: Path) -> dict[tuple[Any, ...], int]:
    index: dict[tuple[Any, ...], int] = {}
    for path in sorted(sample_dir.glob("mlhp_part_*.json")):
        payload = load_json(path)
        for group in payload.get("groups", []):
            key = (
                str(group["chrom"]),
                int(group["start"]),
                int(group["end"]),
                tuple(int(pos) for pos in group.get("positions", [])),
            )
            if key in index:
                raise ValueError(f"Duplicate MLHP group key in {sample_dir.name}: {key}")
            index[key] = int(group.get("n_unique_phase_sets", -1))
    return index


def load_site_ledger(path: Path) -> tuple[dict[str, tuple[tuple[int, str, str], ...]], dict[Any, Any]]:
    group_alleles: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    retained_sites: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("selected_for_read_census") != "true" or row.get("ssnv_branch") != "retained":
                continue
            key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
            group_id = row.get("selected_group_id", "")
            group_alleles[group_id].append((int(row["pos"]), row["ref"], row["alt"]))
            pooled_ref = int(row["pooled_ref_reads"]) if row.get("pooled_ref_reads") else 0
            pooled_alt = int(row["pooled_alt_reads"]) if row.get("pooled_alt_reads") else 0
            pooled_total = pooled_ref + pooled_alt
            caller_af = None
            try:
                caller = json.loads(row.get("sample_values_json") or "{}")
                caller_af = caller.get("SAMPLE", {}).get("AF")
            except json.JSONDecodeError:
                caller_af = None
            record = {
                "caller_af": float(caller_af) if caller_af is not None else None,
                "pooled_read_af": pooled_alt / pooled_total if pooled_total else None,
            }
            if key in retained_sites and retained_sites[key] != record:
                raise ValueError(f"Conflicting retained site rows in {path}: {key}")
            retained_sites[key] = record
    return (
        {group: tuple(sorted(alleles)) for group, alleles in group_alleles.items()},
        retained_sites,
    )


def region_tuple(region: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(region["chrom"]),
        int(region["start"]),
        int(region["end"]),
        tuple(int(pos) for pos in region.get("positions", [])),
    )


def pair_complete(profiles: dict[str, dict[str, Any]]) -> bool:
    return all(
        family in profiles
        and profiles[family]["candidate_set_complete"]
        and not profiles[family]["capped"]
        and profiles[family]["verification_status"] == "full_pass"
        for family in PRIMARY_FAMILIES
    )


def analyze_sample(sample_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]], dict[Any, Any]]:
    sample = sample_dir.name
    region_path = sample_dir / f"layered_region_view_{sample}.json"
    ledger_path = sample_dir / f"ssnv_site_ledger_{sample}.tsv.gz"
    manifest_path = sample_dir / "output_manifest.json"
    region_view = load_json(region_path)
    manifest = load_json(manifest_path)
    phase_sets = load_group_phase_sets(sample_dir)
    group_alleles, sites = load_site_ledger(ledger_path)

    records: list[dict[str, Any]] = []
    aggregate: Counter[str] = Counter()
    state_min_distribution: Counter[int] = Counter()
    mutation_min_distribution: Counter[int] = Counter()
    joint_mutation_min_distribution: Counter[int] = Counter()
    cn_multistate: Counter[str] = Counter()
    missing_phase_join = 0
    allele_position_mismatch = 0
    for region in region_view.get("regions", []):
        aggregate["regions_total"] += 1
        key = region_tuple(region)
        n_phase_sets = phase_sets.get(key)
        if n_phase_sets is None:
            missing_phase_join += 1
            n_phase_sets = -1
        k = int(region.get("n_sSNV", len(region.get("positions", []))))
        primary_raw = {
            str(lineage.get("family")): lineage
            for lineage in region.get("lineages", [])
            if lineage.get("is_primary_lineage") and lineage.get("mutation_bearing")
        }
        profiles = {
            family: analyze_lineage(primary_raw[family], k)
            for family in PRIMARY_FAMILIES
            if family in primary_raw
        }
        for profile in profiles.values():
            aggregate["primary_lineage_units"] += 1
            aggregate["primary_reported_reads"] += profile["reported_read_count"]
            aggregate["primary_retained_reads"] += profile["retained_read_count"]
            aggregate["primary_minread_censored_reads"] += profile["minread_censored_read_count"]
            aggregate["invalid_patterns"] += profile["invalid_pattern_count"]
            state_min_distribution[profile["k_state_min_retained"]] += 1
            mutation_min_distribution[profile["k_mutation_state_min_retained"]] += 1
            if profile["k_mutation_state_min_retained"] >= 2:
                aggregate["primary_lineages_multistate_lower_bound"] += 1
            if profile["catalog_identified_under_retained_table"]:
                aggregate["primary_lineages_exact_retained_catalog"] += 1
            if not profile["has_full_pattern"]:
                aggregate["primary_lineages_no_full_pattern"] += 1
            if profile["minread_censored_read_count"] > 0:
                aggregate["primary_lineages_with_minread_censoring"] += 1

        dual = all(family in profiles for family in PRIMARY_FAMILIES)
        if dual:
            aggregate["regions_dual_primary_hp"] += 1
        current_pair_complete = dual and pair_complete(profiles)
        if current_pair_complete:
            aggregate["regions_pair_complete"] += 1
        mutation_mins = [profile["k_mutation_state_min_retained"] for profile in profiles.values()]
        any_multistate = any(value >= 2 for value in mutation_mins)
        if any_multistate:
            aggregate["regions_any_hp_multistate_lower_bound"] += 1
            cn_multistate[str(region.get("cn"))] += 1
        if any_multistate and current_pair_complete and n_phase_sets == 1:
            aggregate["regions_multistate_pair_complete_single_ps"] += 1
        if any_multistate and str(region.get("cn")) == "neutral":
            aggregate["regions_multistate_cn_neutral_proxy"] += 1
        if any_multistate and current_pair_complete and n_phase_sets == 1 and str(region.get("cn")) == "neutral":
            aggregate["regions_multistate_screened_candidate"] += 1
        if any(profile["catalog_identified_under_retained_table"] for profile in profiles.values()):
            aggregate["regions_any_exact_retained_catalog"] += 1
        if profiles and all(profile["catalog_identified_under_retained_table"] for profile in profiles.values()):
            aggregate["regions_all_available_hp_exact_retained_catalog"] += 1

        joint_min = None
        joint_mutation_min = None
        joint_max = None
        joint_mutation_max = None
        joint_sparse_ceiling = None
        joint_mutation_sparse_ceiling = None
        joint_coupling_degrees_of_freedom = None
        if dual:
            state_supports = [
                profiles[family]["k_state_min_retained"] for family in PRIMARY_FAMILIES
            ]
            mutation_supports = [
                profiles[family]["k_mutation_state_min_retained"] for family in PRIMARY_FAMILIES
            ]
            joint_min = max(state_supports)
            joint_mutation_min = max(
                mutation_supports
            )
            joint_sparse_ceiling = sum(state_supports) - 1
            joint_mutation_sparse_ceiling = sum(mutation_supports) - 1
            joint_coupling_degrees_of_freedom = math.prod(value - 1 for value in state_supports)
            joint_max = math.prod(
                profiles[family]["k_state_max_retained"] for family in PRIMARY_FAMILIES
            )
            joint_mutation_max = math.prod(
                profiles[family]["k_mutation_state_max_retained"] for family in PRIMARY_FAMILIES
            )
            joint_mutation_min_distribution[joint_mutation_min] += 1

        allele_tuple = group_alleles.get(str(region.get("region")), tuple())
        if tuple(pos for pos, _, _ in allele_tuple) != tuple(int(pos) for pos in region.get("positions", [])):
            allele_position_mismatch += 1
        allele_key = (
            str(region["chrom"]),
            tuple((int(pos), ref, alt) for pos, ref, alt in allele_tuple),
        )
        hp_shape_signature = tuple(
            sorted(tuple(profiles[family]["shape_classes"]) for family in PRIMARY_FAMILIES if family in profiles)
        )
        hp_k_signature = tuple(
            sorted(
                (
                    profiles[family]["k_state_min_retained"],
                    profiles[family]["k_mutation_state_min_retained"],
                )
                for family in PRIMARY_FAMILIES
                if family in profiles
            )
        )
        hp_state_k_signature = tuple(
            sorted(
                profiles[family]["k_state_min_retained"]
                for family in PRIMARY_FAMILIES
                if family in profiles
            )
        )
        hp_mutation_k_signature = tuple(
            sorted(
                profiles[family]["k_mutation_state_min_retained"]
                for family in PRIMARY_FAMILIES
                if family in profiles
            )
        )
        records.append(
            {
                "dataset": sample,
                "biological_id": manifest.get("biological_id", sample),
                "region": region.get("region"),
                "chrom": region.get("chrom"),
                "start": int(region.get("start")),
                "end": int(region.get("end")),
                "positions": [int(pos) for pos in region.get("positions", [])],
                "alleles": [list(item) for item in allele_tuple],
                "position_key": [str(region["chrom"]), [int(pos) for pos in region.get("positions", [])]],
                "allele_key_internal": allele_key,
                "n_sSNV": k,
                "cn_label": region.get("cn"),
                "region_determinacy": region.get("region_determinacy"),
                "n_unique_phase_sets": n_phase_sets,
                "single_phase_set": n_phase_sets == 1,
                "has_hp3_auxiliary": any(
                    str(lineage.get("family")) == "3" and lineage.get("mutation_bearing")
                    for lineage in region.get("lineages", [])
                ),
                "primary_hp_count": len(profiles),
                "dual_primary_hp": dual,
                "pair_complete": current_pair_complete,
                "hp_profiles": profiles,
                "hp_shape_signature_flip_invariant": hp_shape_signature,
                "hp_k_signature_flip_invariant": hp_k_signature,
                "hp_state_k_signature_flip_invariant": hp_state_k_signature,
                "hp_mutation_k_signature_flip_invariant": hp_mutation_k_signature,
                "joint_local_state_min_conditional": joint_min,
                "joint_mutation_state_min_conditional": joint_mutation_min,
                "joint_local_state_sparse_ceiling_conditional": joint_sparse_ceiling,
                "joint_mutation_state_sparse_ceiling_conditional": joint_mutation_sparse_ceiling,
                "joint_coupling_degrees_of_freedom_at_min_catalog": (
                    joint_coupling_degrees_of_freedom
                ),
                "joint_local_state_max_retained_combinatorial": joint_max,
                "joint_mutation_state_max_retained_combinatorial": joint_mutation_max,
                "any_hp_multistate_lower_bound": any_multistate,
                "state_evidence_level": (
                    "E2_exact_retained_catalog"
                    if any_multistate
                    and profiles
                    and all(p["catalog_identified_under_retained_table"] for p in profiles.values())
                    else "E1_partial_retained_state_lower_bound"
                    if any_multistate
                    else "E0_no_multistate_lower_bound"
                ),
                "passes_phase_tree_cn_screen": bool(
                    any_multistate
                    and current_pair_complete
                    and n_phase_sets == 1
                    and str(region.get("cn")) == "neutral"
                ),
                "cross_technical_reproduced": False,
                "cross_technical_mutation_k_reproduced": False,
                "biological_clone_k_identified": False,
                "claim_ceiling": (
                    "retained regional mutation-state bounds; not biological clone count or cell-level HP pairing"
                ),
            }
        )

    if missing_phase_join:
        raise ValueError(f"{sample}: {missing_phase_join} regions lack MLHP phase-set join")
    if allele_position_mismatch:
        raise ValueError(f"{sample}: {allele_position_mismatch} regions mismatch site-ledger alleles")
    if aggregate["regions_total"] != len(phase_sets):
        raise ValueError(
            f"{sample}: region conservation failed: view={aggregate['regions_total']} MLHP={len(phase_sets)}"
        )
    summary = {
        "dataset": sample,
        "biological_id": manifest.get("biological_id", sample),
        "input_files": {
            "region_view": str(region_path.resolve()),
            "region_view_sha256": sha256_file(region_path),
            "site_ledger": str(ledger_path.resolve()),
            "site_ledger_sha256": sha256_file(ledger_path),
            "output_manifest": str(manifest_path.resolve()),
            "output_manifest_sha256": sha256_file(manifest_path),
        },
        "counts": dict(sorted(aggregate.items())),
        "distributions": {
            "primary_lineage_k_state_min_retained": {
                str(key): value for key, value in sorted(state_min_distribution.items())
            },
            "primary_lineage_k_mutation_state_min_retained": {
                str(key): value for key, value in sorted(mutation_min_distribution.items())
            },
            "dual_hp_joint_mutation_state_min_conditional": {
                str(key): value for key, value in sorted(joint_mutation_min_distribution.items())
            },
            "cn_label_among_multistate_lower_bound_regions": dict(sorted(cn_multistate.items())),
        },
        "quality": {
            "region_to_mlhp_join_pass": True,
            "region_to_ledger_allele_join_pass": True,
            "invalid_pattern_count": int(aggregate["invalid_patterns"]),
            "minread_censored_read_fraction": (
                aggregate["primary_minread_censored_reads"] / aggregate["primary_reported_reads"]
                if aggregate["primary_reported_reads"]
                else None
            ),
        },
    }
    return summary, records, sites


def concordance_correlation(left: list[float], right: list[float]) -> float | None:
    if len(left) < 2:
        return None
    mean_left = sum(left) / len(left)
    mean_right = sum(right) / len(right)
    var_left = sum((value - mean_left) ** 2 for value in left) / len(left)
    var_right = sum((value - mean_right) ** 2 for value in right) / len(right)
    covariance = sum(
        (a - mean_left) * (b - mean_right) for a, b in zip(left, right)
    ) / len(left)
    denominator = var_left + var_right + (mean_left - mean_right) ** 2
    return 2 * covariance / denominator if denominator else None


def numeric_concordance(
    left_sites: dict[Any, Any], right_sites: dict[Any, Any], field: str
) -> dict[str, Any]:
    from scipy.stats import pearsonr, spearmanr

    common = sorted(set(left_sites) & set(right_sites))
    pairs = [
        (left_sites[key].get(field), right_sites[key].get(field))
        for key in common
        if left_sites[key].get(field) is not None and right_sites[key].get(field) is not None
    ]
    left = [float(pair[0]) for pair in pairs]
    right = [float(pair[1]) for pair in pairs]
    absolute = sorted(abs(a - b) for a, b in zip(left, right))
    median = statistics.median(absolute) if absolute else None
    return {
        "n_common_evaluable": len(pairs),
        "pearson_r": float(pearsonr(left, right).statistic) if len(pairs) >= 2 else None,
        "spearman_rho": float(spearmanr(left, right).statistic) if len(pairs) >= 2 else None,
        "concordance_correlation_coefficient": concordance_correlation(left, right),
        "median_absolute_difference": median,
    }


def compare_samples(
    left_records: list[dict[str, Any]],
    right_records: list[dict[str, Any]],
    left_sites: dict[Any, Any],
    right_sites: dict[Any, Any],
) -> dict[str, Any]:
    left_by_allele = {record["allele_key_internal"]: record for record in left_records}
    right_by_allele = {record["allele_key_internal"]: record for record in right_records}
    left_by_position = {
        (record["chrom"], tuple(record["positions"])): record for record in left_records
    }
    right_by_position = {
        (record["chrom"], tuple(record["positions"])): record for record in right_records
    }
    common_allele = sorted(set(left_by_allele) & set(right_by_allele))
    common_position = sorted(set(left_by_position) & set(right_by_position))
    shared_records = [(left_by_allele[key], right_by_allele[key]) for key in common_allele]
    for left, right in shared_records:
        state_reproduced = (
            left["hp_k_signature_flip_invariant"] == right["hp_k_signature_flip_invariant"]
        )
        left["cross_technical_reproduced"] = state_reproduced
        right["cross_technical_reproduced"] = state_reproduced
        mutation_reproduced = (
            left["hp_mutation_k_signature_flip_invariant"]
            == right["hp_mutation_k_signature_flip_invariant"]
        )
        left["cross_technical_mutation_k_reproduced"] = mutation_reproduced
        right["cross_technical_mutation_k_reproduced"] = mutation_reproduced

    dual_pairs = [(left, right) for left, right in shared_records if left["dual_primary_hp"] and right["dual_primary_hp"]]
    pair_complete_pairs = [(left, right) for left, right in dual_pairs if left["pair_complete"] and right["pair_complete"]]
    left_candidates = {
        key for key, record in left_by_allele.items() if record["any_hp_multistate_lower_bound"]
    }
    right_candidates = {
        key for key, record in right_by_allele.items() if record["any_hp_multistate_lower_bound"]
    }
    candidate_intersection = left_candidates & right_candidates
    left_exact_catalog_candidates = {
        key
        for key, record in left_by_allele.items()
        if record["state_evidence_level"] == "E2_exact_retained_catalog"
    }
    right_exact_catalog_candidates = {
        key
        for key, record in right_by_allele.items()
        if record["state_evidence_level"] == "E2_exact_retained_catalog"
    }
    left_screened_candidates = {
        key for key, record in left_by_allele.items() if record["passes_phase_tree_cn_screen"]
    }
    right_screened_candidates = {
        key for key, record in right_by_allele.items() if record["passes_phase_tree_cn_screen"]
    }
    site_left = set(left_sites)
    site_right = set(right_sites)
    return {
        "exact_retained_site_overlap": {
            "left": len(site_left),
            "right": len(site_right),
            "intersection": len(site_left & site_right),
            "jaccard": len(site_left & site_right) / len(site_left | site_right),
        },
        "exact_position_region_overlap": {
            "left": len(left_by_position),
            "right": len(right_by_position),
            "intersection": len(common_position),
            "jaccard": len(common_position) / len(set(left_by_position) | set(right_by_position)),
        },
        "exact_allele_region_overlap": {
            "left": len(left_by_allele),
            "right": len(right_by_allele),
            "intersection": len(common_allele),
            "jaccard": len(common_allele) / len(set(left_by_allele) | set(right_by_allele)),
        },
        "shared_exact_allele_regions": {
            "n": len(shared_records),
            "region_determinacy_agreement": sum(
                left["region_determinacy"] == right["region_determinacy"] for left, right in shared_records
            ),
            "primary_hp_count_agreement": sum(
                left["primary_hp_count"] == right["primary_hp_count"] for left, right in shared_records
            ),
            "flip_invariant_k_signature_agreement": sum(
                left["hp_k_signature_flip_invariant"] == right["hp_k_signature_flip_invariant"]
                for left, right in shared_records
            ),
            "flip_invariant_state_k_signature_agreement": sum(
                left["hp_state_k_signature_flip_invariant"]
                == right["hp_state_k_signature_flip_invariant"]
                for left, right in shared_records
            ),
            "flip_invariant_mutation_k_signature_agreement": sum(
                left["hp_mutation_k_signature_flip_invariant"]
                == right["hp_mutation_k_signature_flip_invariant"]
                for left, right in shared_records
            ),
            "both_dual_primary_hp": len(dual_pairs),
            "both_pair_complete": len(pair_complete_pairs),
            "flip_invariant_shape_signature_agreement_among_pair_complete": sum(
                left["hp_shape_signature_flip_invariant"]
                == right["hp_shape_signature_flip_invariant"]
                for left, right in pair_complete_pairs
            ),
            "joint_mutation_state_min_agreement_among_dual": sum(
                left["joint_mutation_state_min_conditional"]
                == right["joint_mutation_state_min_conditional"]
                for left, right in dual_pairs
            ),
        },
        "multistate_candidate_overlap": {
            "left": len(left_candidates),
            "right": len(right_candidates),
            "intersection": len(candidate_intersection),
            "jaccard": (
                len(candidate_intersection) / len(left_candidates | right_candidates)
                if left_candidates or right_candidates
                else None
            ),
            "left_reproduced_fraction": (
                len(candidate_intersection) / len(left_candidates) if left_candidates else None
            ),
            "right_reproduced_fraction": (
                len(candidate_intersection) / len(right_candidates) if right_candidates else None
            ),
        },
        "exact_retained_catalog_multistate_overlap": {
            "left": len(left_exact_catalog_candidates),
            "right": len(right_exact_catalog_candidates),
            "intersection": len(left_exact_catalog_candidates & right_exact_catalog_candidates),
            "jaccard": (
                len(left_exact_catalog_candidates & right_exact_catalog_candidates)
                / len(left_exact_catalog_candidates | right_exact_catalog_candidates)
                if left_exact_catalog_candidates or right_exact_catalog_candidates
                else None
            ),
        },
        "phase_tree_cn_screened_multistate_overlap": {
            "left": len(left_screened_candidates),
            "right": len(right_screened_candidates),
            "intersection": len(left_screened_candidates & right_screened_candidates),
            "jaccard": (
                len(left_screened_candidates & right_screened_candidates)
                / len(left_screened_candidates | right_screened_candidates)
                if left_screened_candidates or right_screened_candidates
                else None
            ),
        },
        "caller_af_concordance": numeric_concordance(left_sites, right_sites, "caller_af"),
        "pooled_read_af_concordance": numeric_concordance(
            left_sites, right_sites, "pooled_read_af"
        ),
    }


def serializable_record(record: dict[str, Any]) -> dict[str, Any]:
    result = dict(record)
    result.pop("allele_key_internal", None)
    result["hp_shape_signature_flip_invariant"] = [
        list(item) for item in result["hp_shape_signature_flip_invariant"]
    ]
    result["hp_k_signature_flip_invariant"] = [
        list(item) for item in result["hp_k_signature_flip_invariant"]
    ]
    result["hp_state_k_signature_flip_invariant"] = list(
        result["hp_state_k_signature_flip_invariant"]
    )
    result["hp_mutation_k_signature_flip_invariant"] = list(
        result["hp_mutation_k_signature_flip_invariant"]
    )
    return result


def write_region_tsv(path: Path, records: list[dict[str, Any]]) -> None:
    fields = [
        "dataset",
        "biological_id",
        "region",
        "chrom",
        "start",
        "end",
        "positions",
        "alleles",
        "n_sSNV",
        "cn_label",
        "region_determinacy",
        "n_unique_phase_sets",
        "single_phase_set",
        "has_hp3_auxiliary",
        "primary_hp_count",
        "dual_primary_hp",
        "pair_complete",
        "hp1_k_state_min_retained",
        "hp1_k_mutation_state_min_retained",
        "hp1_k_state_max_retained",
        "hp1_catalog_identified_under_retained_table",
        "hp2_k_state_min_retained",
        "hp2_k_mutation_state_min_retained",
        "hp2_k_state_max_retained",
        "hp2_catalog_identified_under_retained_table",
        "joint_local_state_min_conditional",
        "joint_mutation_state_min_conditional",
        "joint_local_state_sparse_ceiling_conditional",
        "joint_mutation_state_sparse_ceiling_conditional",
        "joint_coupling_degrees_of_freedom_at_min_catalog",
        "joint_local_state_max_retained_combinatorial",
        "joint_mutation_state_max_retained_combinatorial",
        "any_hp_multistate_lower_bound",
        "state_evidence_level",
        "passes_phase_tree_cn_screen",
        "cross_technical_reproduced",
        "cross_technical_mutation_k_reproduced",
        "biological_clone_k_identified",
        "claim_ceiling",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for original in records:
            row = serializable_record(original)
            for family in PRIMARY_FAMILIES:
                profile = row["hp_profiles"].get(family, {})
                prefix = f"hp{family}"
                row[f"{prefix}_k_state_min_retained"] = profile.get("k_state_min_retained")
                row[f"{prefix}_k_mutation_state_min_retained"] = profile.get(
                    "k_mutation_state_min_retained"
                )
                row[f"{prefix}_k_state_max_retained"] = profile.get("k_state_max_retained")
                row[f"{prefix}_catalog_identified_under_retained_table"] = profile.get(
                    "catalog_identified_under_retained_table"
                )
            row["positions"] = json.dumps(row["positions"], ensure_ascii=False)
            row["alleles"] = json.dumps(row["alleles"], ensure_ascii=False)
            writer.writerow({field: row.get(field) for field in fields})


def run(run_dir: Path, out_dir: Path, sample_names: list[str]) -> dict[str, Any]:
    state_path = run_dir / "run_state.json"
    run_state = load_json(state_path)
    if run_state.get("state") != "SUCCEEDED" or not (run_dir / "_SUCCESS").exists():
        raise ValueError(f"Canonical input is not successful: {run_dir}")
    summaries = []
    records_by_sample: dict[str, list[dict[str, Any]]] = {}
    sites_by_sample: dict[str, dict[Any, Any]] = {}
    for sample in sample_names:
        summary, records, sites = analyze_sample(run_dir / "samples" / sample)
        summaries.append(summary)
        records_by_sample[sample] = records
        sites_by_sample[sample] = sites
    comparison = compare_samples(
        records_by_sample[sample_names[0]],
        records_by_sample[sample_names[1]],
        sites_by_sample[sample_names[0]],
        sites_by_sample[sample_names[1]],
    )
    all_records = [record for sample in sample_names for record in records_by_sample[sample]]
    priority_records = [
        serializable_record(record)
        for record in all_records
        if record["passes_phase_tree_cn_screen"]
        or record["state_evidence_level"] == "E2_exact_retained_catalog"
        or (record["joint_mutation_state_min_conditional"] or 0) >= 3
    ]
    priority_records.sort(
        key=lambda record: (
            record["dataset"],
            -(record["joint_mutation_state_min_conditional"] or 0),
            record["chrom"],
            record["start"],
        )
    )
    payload = {
        "schema_name": "intersubmod.hcc1395_unknown_k_state_bounds",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "comprehensive_validation",
        "research_goals": ["G3", "G4", "G5"],
        "input_run": {
            "path": str(run_dir.resolve()),
            "run_id": run_state.get("run_id"),
            "state": run_state.get("state"),
            "run_state_sha256": sha256_file(state_path),
            "scope": "chr1-22; HCC1395 and HCC1395_DORADO",
        },
        "model_contract": {
            "estimand": (
                "minimum and maximum HP-specific full R/A mutation states compatible with the "
                "MINREAD-retained R/A/X pattern table"
            ),
            "k_state_min_retained": (
                "exact chromatic number of the retained-pattern conflict graph"
            ),
            "k_state_max_retained": (
                "maximum distinct full states supportable by retained aggregate pattern counts, "
                "computed as capacitated bipartite max flow"
            ),
            "joint_local_state_min_conditional": (
                "support-only lower bound max(HP1 K_min, HP2 K_min); conditional on error-free "
                "HP states and valid phasing; unknown marginal weights can require more joint states"
            ),
            "joint_sparse_support_note": (
                "for fixed HP support sizes r1,r2 and known marginal weights, minimum coupling "
                "support lies between max(r1,r2) and r1+r2-1"
            ),
            "joint_state_upper_note": (
                "HP1_max*HP2_max is only a retained-table combinatorial ceiling, not a biological "
                "clone upper bound"
            ),
            "biological_clone_k_identified": False,
            "why_not_identified": [
                "HP1 and HP2 are marginal bulk read groups without cell-level pairing",
                "patterns with count below MINREAD=3 are absent from the retained table",
                "many reads contain X and do not span all selected sSNVs",
                "allele-specific CN=1+1 and subclonal CNA are not resolved",
                "HP/base/mapping errors are not calibrated in this inversion",
                "different genome-wide clones can share the same local regional genotype",
            ],
        },
        "evidence_levels": {
            "E0": "no retained-pattern evidence requiring >=2 mutation states",
            "E1": "partial retained patterns require >=2 compatible mutation states",
            "E2": "all retained patterns are full length, so retained state catalog is exact",
            "E3_future": "uncensored per-molecule model + phase/error/allele-specific-CN gates",
            "E4_future": "cross-technical plus orthogonal methyl/SNV validation",
            "E5_future": "cell-resolved or single-cell biological clone confirmation",
        },
        "sample_summaries": summaries,
        "cross_technical_comparison": comparison,
        "priority_candidate_definition": (
            "phase/tree/CN-screened OR exact retained multistate catalog OR dual-HP joint "
            "mutation-state lower bound>=3; prioritization only, not clone confirmation"
        ),
        "priority_candidate_records": priority_records,
        "claim_ceiling": (
            "Current data can prioritize regional mutation-state candidates and quantify conditional "
            "state bounds. They cannot confirm unknown biological clone K, clone fractions, cross-HP "
            "cell pairing, or a unique true clone tree."
        ),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "unknown_k_consistency.json"
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    write_region_tsv(out_dir / "region_clone_state_bounds.tsv", all_records)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--samples", nargs=2, default=["HCC1395", "HCC1395_DORADO"], metavar=("LEFT", "RIGHT")
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = run(args.run_dir.resolve(), args.out_dir.resolve(), list(args.samples))
    print(
        json.dumps(
            {
                "status": "PASS",
                "output": str((args.out_dir / "unknown_k_consistency.json").resolve()),
                "samples": [summary["dataset"] for summary in payload["sample_summaries"]],
                "strict_biological_clone_k_identified": False,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
