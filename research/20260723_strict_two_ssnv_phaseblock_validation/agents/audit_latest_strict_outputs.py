#!/usr/bin/env python3
"""Read-only independent audit of completed strict exact-PS/HP outputs."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
from typing import Iterable


AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
PRIMARY_THRESHOLD = 3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample-root",
        action="append",
        required=True,
        help="DATASET=/absolute/path/to/strict_regions_root",
    )
    return parser.parse_args()


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected JSON object")
    return value


def read_tsv_gz(path: Path) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def sha256(path: Path, cache: dict[Path, str]) -> str:
    path = path.resolve()
    if path not in cache:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        cache[path] = digest.hexdigest()
    return cache[path]


def identity_ok(identity: dict, cache: dict[Path, str]) -> tuple[bool, str]:
    path = Path(str(identity["path"]))
    if not path.is_file():
        return False, f"missing:{path}"
    if path.stat().st_size != int(identity["size_bytes"]):
        return False, f"size:{path}"
    if sha256(path, cache) != identity["sha256"]:
        return False, f"sha256:{path}"
    return True, ""


def is_connected(vertices: set[int], edges: list[tuple[int, int]]) -> bool:
    if not vertices:
        return False
    adjacency: dict[int, set[int]] = {vertex: set() for vertex in vertices}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    seen: set[int] = set()
    stack = [min(vertices)]
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        seen.add(current)
        stack.extend(adjacency[current] - seen)
    return seen == vertices


def audit_sample(
    dataset: str, root: Path, hash_cache: dict[Path, str]
) -> dict[str, object]:
    violations: Counter[str] = Counter()
    examples: dict[str, list[str]] = defaultdict(list)
    totals: Counter[str] = Counter()
    active_unique_loci: set[tuple[str, int]] = set()
    containers: set[tuple[str, str, str]] = set()
    latest_mtime = 0.0
    receipt_created_at: list[str] = []
    observed_chromosomes: list[str] = []

    def fail(kind: str, message: str) -> None:
        violations[kind] += 1
        if len(examples[kind]) < 5:
            examples[kind].append(message)

    for chrom in AUTOSOMES:
        chrom_dir = root / "chromosomes" / chrom
        receipt_path = chrom_dir / "receipt.json"
        if not receipt_path.is_file():
            fail("missing_chromosome_receipt", chrom)
            continue
        observed_chromosomes.append(chrom)
        receipt = read_json(receipt_path)
        receipt_created_at.append(str(receipt.get("created_at_utc", "")))
        latest_mtime = max(
            latest_mtime,
            max(path.stat().st_mtime for path in chrom_dir.iterdir() if path.is_file()),
        )
        if receipt.get("all_pass") is not True:
            fail("receipt_all_pass_false", chrom)
        if any(value is not True for value in (receipt.get("checks") or {}).values()):
            fail("receipt_check_false", chrom)
        if receipt.get("scope", {}).get("dataset") != dataset:
            fail("receipt_dataset_mismatch", chrom)
        if receipt.get("scope", {}).get("chrom") != chrom:
            fail("receipt_chrom_mismatch", chrom)

        sidecar = receipt_path.with_suffix(".json.sha256")
        if not sidecar.is_file():
            fail("receipt_sidecar_missing", chrom)
        else:
            expected = sidecar.read_text(encoding="utf-8").split()[0]
            if expected != sha256(receipt_path, hash_cache):
                fail("receipt_sidecar_hash_mismatch", chrom)
            totals["receipt_sidecars_verified"] += 1

        for group in ("inputs", "outputs"):
            for name, identity in (receipt.get(group) or {}).items():
                ok, reason = identity_ok(identity, hash_cache)
                totals["identities_checked"] += 1
                if not ok:
                    fail("identity_mismatch", f"{chrom}:{group}:{name}:{reason}")

        site_path = Path(receipt["inputs"]["site_catalog"]["path"])
        edge_path = Path(receipt["outputs"]["edges"]["path"])
        component_path = Path(receipt["outputs"]["components"]["path"])
        membership_path = Path(receipt["outputs"]["membership"]["path"])

        catalog: dict[int, tuple[str, int]] = {}
        for row in read_tsv_gz(site_path):
            index = int(row["site_index"])
            value = (row["chrom"], int(row["pos1"]))
            if index in catalog:
                fail("duplicate_catalog_index", f"{chrom}:{index}")
            catalog[index] = value
        totals["candidate_loci_S"] += len(catalog)

        members_by_component: dict[str, list[dict[str, str]]] = defaultdict(list)
        membership_by_node: dict[tuple[str, str, int], str] = {}
        for row in read_tsv_gz(membership_path):
            if int(row["threshold"]) != PRIMARY_THRESHOLD:
                continue
            basis = row["linkage_basis"]
            phase_set = row["phase_set"]
            pos = int(row["pos1"])
            component_id = row["component_id"]
            node = (basis, phase_set, pos)
            if node in membership_by_node:
                fail("duplicate_primary_membership", f"{chrom}:{node}")
            membership_by_node[node] = component_id
            members_by_component[component_id].append(row)
            active_unique_loci.add((chrom, pos))
            containers.add((chrom, basis, phase_set))
            totals["active_node_memberships"] += 1
            if (
                not phase_set
                or basis not in {"PS_HP1", "PS_HP2"}
                or row["dataset"] != dataset
                or row["chrom"] != chrom
            ):
                fail("missing_or_invalid_membership_ps_hp", f"{chrom}:{node}")

        primary_edges_by_component: dict[
            str, list[tuple[int, int]]
        ] = defaultdict(list)
        observed_edge_rows = 0
        observed_edge_state_mass = 0
        retained_edges = 0
        for row in read_tsv_gz(edge_path):
            observed_edge_rows += 1
            support = int(row["support_total"])
            state_mass = sum(
                int(row[field])
                for field in ("support_RR", "support_RA", "support_AR", "support_AA")
            )
            observed_edge_state_mass += support
            if state_mass != support:
                fail("edge_state_mass_mismatch", f"{chrom}:{row['pos_i1']}-{row['pos_j1']}")
            passes = row["passes_primary_threshold"] == "true"
            if passes != (support >= PRIMARY_THRESHOLD):
                fail("primary_threshold_flag_mismatch", f"{chrom}:{support}:{passes}")
            if not passes:
                continue
            retained_edges += 1
            left_index = int(row["site_i_index"])
            right_index = int(row["site_j_index"])
            left = int(row["pos_i1"])
            right = int(row["pos_j1"])
            basis = row["linkage_basis"]
            hp = row["hp_family"]
            phase_set = row["phase_set"]
            if support < PRIMARY_THRESHOLD:
                fail("primary_edge_support_lt3", f"{chrom}:{left}-{right}:{support}")
            if left == right or left_index == right_index:
                fail("primary_edge_same_endpoint", f"{chrom}:{left_index}:{left}")
            if catalog.get(left_index) != (chrom, left):
                fail("left_endpoint_not_catalog", f"{chrom}:{left_index}:{left}")
            if catalog.get(right_index) != (chrom, right):
                fail("right_endpoint_not_catalog", f"{chrom}:{right_index}:{right}")
            if (
                hp not in {"1", "2"}
                or basis != f"PS_HP{hp}"
                or not phase_set
                or row["dataset"] != dataset
                or row["chrom"] != chrom
            ):
                fail("missing_or_invalid_edge_ps_hp", f"{chrom}:{basis}:{hp}:{phase_set}")
            left_component = membership_by_node.get((basis, phase_set, left))
            right_component = membership_by_node.get((basis, phase_set, right))
            if left_component is None or right_component is None:
                fail("edge_endpoint_missing_primary_membership", f"{chrom}:{left}-{right}")
            elif left_component != right_component:
                fail("edge_cross_component", f"{chrom}:{left_component}:{right_component}")
            else:
                primary_edges_by_component[left_component].append((left, right))

        totals["endpoint_pair_rows_observed"] += observed_edge_rows
        totals["endpoint_pair_state_mass"] += observed_edge_state_mass
        totals["retained_endpoint_edges"] += retained_edges
        if observed_edge_rows != int(receipt["counts"]["endpoint_pair_rows_observed"]):
            fail("receipt_edge_rows_mismatch", chrom)
        if observed_edge_state_mass != int(receipt["counts"]["endpoint_pair_state_mass"]):
            fail("receipt_edge_mass_mismatch", chrom)

        observed_components: dict[str, dict[str, str]] = {}
        tree_regions = 0
        singleton_regions = 0
        membership_mass = 0
        max_k = 0
        k_gt12 = 0
        for row in read_tsv_gz(component_path):
            if int(row["threshold"]) != PRIMARY_THRESHOLD:
                continue
            component_id = row["component_id"]
            if component_id in observed_components:
                fail("duplicate_component_id", f"{chrom}:{component_id}")
            observed_components[component_id] = row
            k = int(row["k"])
            membership_mass += k
            max_k = max(max_k, k)
            members = members_by_component.get(component_id, [])
            if len(members) != k:
                fail("component_k_membership_mismatch", f"{chrom}:{component_id}:{k}:{len(members)}")
            basis = row["linkage_basis"]
            phase_set = row["phase_set"]
            if not phase_set or basis not in {"PS_HP1", "PS_HP2"}:
                fail("missing_or_invalid_component_ps_hp", f"{chrom}:{component_id}")
            if any(
                member["linkage_basis"] != basis
                or member["phase_set"] != phase_set
                for member in members
            ):
                fail("component_cross_ps_hp", f"{chrom}:{component_id}")
            vertices = {int(member["pos1"]) for member in members}
            tree_eligible = row["tree_eligible"] == "true"
            if tree_eligible:
                tree_regions += 1
                if k <= 1 or row["inference_role"] != "PRIMARY_PS_AWARE":
                    fail("invalid_tree_eligible_role_or_k", f"{chrom}:{component_id}:{k}")
                if not is_connected(vertices, primary_edges_by_component.get(component_id, [])):
                    fail("tree_component_disconnected", f"{chrom}:{component_id}")
                if k > 12:
                    k_gt12 += 1
            else:
                singleton_regions += 1
                if k != 1 or row["inference_role"] != "ABSTAIN_SINGLETON_UNLINKED":
                    fail("invalid_singleton_role_or_k", f"{chrom}:{component_id}:{k}")

        if set(members_by_component) != set(observed_components):
            fail("component_membership_id_set_mismatch", chrom)
        if membership_mass != len(membership_by_node):
            fail("component_membership_mass_not_conserved", chrom)
        threshold = receipt["summary_by_threshold"][str(PRIMARY_THRESHOLD)]
        expected_pairs = {
            "active_node_memberships": len(membership_by_node),
            "components_all": len(observed_components),
            "singleton_unlinked_components": singleton_regions,
            "tree_eligible_regions": tree_regions,
            "retained_endpoint_edges": retained_edges,
        }
        for key, value in expected_pairs.items():
            if int(threshold[key]) != value:
                fail("receipt_threshold_summary_mismatch", f"{chrom}:{key}:{value}:{threshold[key]}")

        totals["components_all_including_singletons"] += len(observed_components)
        totals["singleton_unlinked_components_excluded_from_tree"] += singleton_regions
        totals["tree_eligible_read_linked_regions"] += tree_regions
        totals["k_gt12_regions"] += k_gt12
        totals["max_k"] = max(totals["max_k"], max_k)
        totals["canonical_molecule_rows"] += int(receipt["counts"]["molecule_rows_total"])
        totals["primary_known_ps_molecule_rows"] += int(receipt["counts"]["primary_known_ps_rows"])

    summary_path = root / "summary" / "summary.json"
    if not summary_path.is_file():
        fail("summary_missing", str(summary_path))
        summary = {}
    else:
        summary = read_json(summary_path)
        latest_mtime = max(latest_mtime, summary_path.stat().st_mtime)
        if summary.get("all_pass") is not True:
            fail("summary_all_pass_false", dataset)
        if any(value is not True for value in (summary.get("checks") or {}).values()):
            fail("summary_check_false", dataset)
        for chrom, identity in (summary.get("inputs") or {}).items():
            ok, reason = identity_ok(identity, hash_cache)
            totals["summary_receipt_identities_checked"] += 1
            if not ok:
                fail("summary_receipt_identity_mismatch", f"{chrom}:{reason}")
        aggregate = summary.get("aggregate") or {}
        totals["active_unique_loci"] = len(active_unique_loci)
        totals["exact_ps_hp_containers_with_active_nodes"] = len(containers)
        aggregate_keys = (
            "candidate_loci_S",
            "active_node_memberships",
            "active_unique_loci",
            "components_all_including_singletons",
            "exact_ps_hp_containers_with_active_nodes",
            "retained_endpoint_edges",
            "singleton_unlinked_components_excluded_from_tree",
            "tree_eligible_read_linked_regions",
            "k_gt12_regions",
            "max_k",
            "canonical_molecule_rows",
            "primary_known_ps_molecule_rows",
        )
        for key in aggregate_keys:
            if int(aggregate.get(key, -1)) != totals[key]:
                fail("genome_summary_aggregate_mismatch", f"{key}:{totals[key]}:{aggregate.get(key)}")

    return {
        "dataset": dataset,
        "root": str(root.resolve()),
        "complete_chr1_22": observed_chromosomes == list(AUTOSOMES),
        "observed_chromosomes": observed_chromosomes,
        "receipt_created_at_utc_min": min(receipt_created_at) if receipt_created_at else None,
        "receipt_created_at_utc_max": max(receipt_created_at) if receipt_created_at else None,
        "latest_artifact_mtime_utc": (
            datetime.fromtimestamp(latest_mtime, tz=timezone.utc).isoformat()
            if latest_mtime
            else None
        ),
        "counts": dict(sorted(totals.items())),
        "violation_counts": dict(sorted(violations.items())),
        "violation_examples": dict(sorted(examples.items())),
        "all_pass": not violations and observed_chromosomes == list(AUTOSOMES),
    }


def main() -> int:
    args = parse_args()
    roots: list[tuple[str, Path]] = []
    for spec in args.sample_root:
        if "=" not in spec:
            raise ValueError(f"invalid --sample-root: {spec}")
        dataset, path = spec.split("=", 1)
        roots.append((dataset, Path(path)))
    cache: dict[Path, str] = {}
    samples = [audit_sample(dataset, root, cache) for dataset, root in roots]
    output = {
        "schema_name": "intersubmod.latest_strict_output_independent_audit",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "parameters": {
            "primary_threshold": PRIMARY_THRESHOLD,
            "expected_chromosomes": list(AUTOSOMES),
        },
        "samples": samples,
        "aggregate": {
            "sample_roots_audited": len(samples),
            "sample_roots_all_pass": sum(sample["all_pass"] for sample in samples),
            "chromosome_receipts_audited": sum(
                len(sample["observed_chromosomes"]) for sample in samples
            ),
            "identity_files_hashed": len(cache),
            "violations_total": sum(
                sum(sample["violation_counts"].values()) for sample in samples
            ),
        },
    }
    print(json.dumps(output, ensure_ascii=False, indent=2, sort_keys=True))
    return 0 if output["aggregate"]["violations_total"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
