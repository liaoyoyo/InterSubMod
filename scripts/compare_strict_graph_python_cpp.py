#!/usr/bin/env python3
"""Compare Python and C++ strict endpoint graphs by scientific semantics."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable, Sequence


SCHEMA_NAME = "intersubmod.strict_endpoint_python_cpp_parity"
SCHEMA_VERSION = "1.0.0"


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_path(path),
    }


def read_tsv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:  # type: ignore[arg-type]
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(set(required) - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{path}: missing columns {','.join(missing)}")
        return [dict(row) for row in reader]


def normalized_hp(value: str) -> str:
    if value in {"1", "HP1", "PS_HP1"}:
        return "HP1"
    if value in {"2", "HP2", "PS_HP2"}:
        return "HP2"
    raise ValueError(f"unexpected primary HP value: {value!r}")


def edge_table(path: Path, *, python: bool, threshold: int) -> dict[tuple[Any, ...], tuple[Any, ...]]:
    if python:
        rows = read_tsv(
            path,
            (
                "dataset", "chrom", "linkage_basis", "phase_set", "site_i_index",
                "site_j_index", "pos_i1", "pos_j1", "support_total", "support_RR",
                "support_RA", "support_AR", "support_AA",
            ),
        )
    else:
        rows = read_tsv(
            path,
            (
                "dataset", "chrom", "phase_set", "hp_family", "left_site_index",
                "right_site_index", "left_pos1", "right_pos1", "support_total",
                "support_RR", "support_RA", "support_AR", "support_AA", "threshold", "qualifies",
            ),
        )
    result = {}
    for row in rows:
        hp = normalized_hp(row["linkage_basis"] if python else row["hp_family"])
        left_index = int(row["site_i_index"] if python else row["left_site_index"])
        right_index = int(row["site_j_index"] if python else row["right_site_index"])
        key = (row["dataset"], row["chrom"], row["phase_set"], hp, left_index, right_index)
        total = int(row["support_total"])
        value = (
            int(row["pos_i1"] if python else row["left_pos1"]),
            int(row["pos_j1"] if python else row["right_pos1"]),
            total,
            int(row["support_RR"]),
            int(row["support_RA"]),
            int(row["support_AR"]),
            int(row["support_AA"]),
            total >= threshold,
        )
        if not python and (
            int(row["threshold"]) != threshold
            or bool(int(row["qualifies"])) != (total >= threshold)
        ):
            raise ValueError(f"{path}: C++ threshold/qualifies mismatch")
        if key in result:
            raise ValueError(f"{path}: duplicate edge key {key}")
        result[key] = value
    return result


def python_components(
    component_path: Path, membership_path: Path, threshold: int
) -> dict[tuple[Any, ...], tuple[Any, ...]]:
    components = read_tsv(
        component_path,
        (
            "dataset", "chrom", "linkage_basis", "phase_set", "threshold", "component_id",
            "k", "retained_edge_count", "component_class", "tree_eligible", "inference_role",
        ),
    )
    memberships = read_tsv(
        membership_path,
        ("threshold", "component_id", "site_index", "pos1"),
    )
    members: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in memberships:
        if int(row["threshold"]) == threshold:
            members[row["component_id"]].append((int(row["site_index"]), int(row["pos1"])))
    result = {}
    for row in components:
        if int(row["threshold"]) != threshold:
            continue
        calls = tuple(sorted(members.get(row["component_id"], ())))
        if len(calls) != int(row["k"]):
            raise ValueError(f"{component_path}: k/membership mismatch")
        key = (
            row["dataset"],
            row["chrom"],
            row["phase_set"],
            normalized_hp(row["linkage_basis"]),
            tuple(index for index, _ in calls),
        )
        value = (
            tuple(position for _, position in calls),
            int(row["retained_edge_count"]),
            row["component_class"],
            row["tree_eligible"],
            row["inference_role"],
        )
        if key in result:
            raise ValueError(f"{component_path}: duplicate semantic component {key}")
        result[key] = value
    return result


def cpp_components(path: Path, threshold: int) -> dict[tuple[Any, ...], tuple[Any, ...]]:
    rows = read_tsv(
        path,
        (
            "dataset", "chrom", "phase_set", "hp_family", "threshold", "k", "site_indices",
            "positions1", "qualifying_edge_count", "component_class", "tree_eligible", "inference_role",
        ),
    )
    result = {}
    for row in rows:
        if int(row["threshold"]) != threshold:
            raise ValueError(f"{path}: C++ component threshold mismatch")
        indices = tuple(int(value) for value in row["site_indices"].split(","))
        positions = tuple(int(value) for value in row["positions1"].split(","))
        if len(indices) != int(row["k"]) or len(indices) != len(positions):
            raise ValueError(f"{path}: C++ k/vector mismatch")
        key = (row["dataset"], row["chrom"], row["phase_set"], normalized_hp(row["hp_family"]), indices)
        value = (
            positions,
            int(row["qualifying_edge_count"]),
            row["component_class"],
            row["tree_eligible"],
            row["inference_role"],
        )
        if key in result:
            raise ValueError(f"{path}: duplicate semantic component {key}")
        result[key] = value
    return result


def mismatch_preview(left: dict[Any, Any], right: dict[Any, Any], limit: int = 20) -> list[dict[str, Any]]:
    result = []
    for key in sorted(set(left) | set(right), key=repr):
        if left.get(key) != right.get(key):
            result.append({"key": repr(key), "python": repr(left.get(key)), "cpp": repr(right.get(key))})
            if len(result) >= limit:
                break
    return result


def compare(
    *,
    python_edges: Path,
    python_components_path: Path,
    python_membership: Path,
    cpp_edges: Path,
    cpp_components_path: Path,
    cpp_receipt: Path,
    threshold: int,
    output: Path,
) -> dict[str, Any]:
    paths = tuple(
        path.resolve(strict=True)
        for path in (
            python_edges, python_components_path, python_membership,
            cpp_edges, cpp_components_path, cpp_receipt,
        )
    )
    py_edges = edge_table(paths[0], python=True, threshold=threshold)
    cxx_edges = edge_table(paths[3], python=False, threshold=threshold)
    py_components = python_components(paths[1], paths[2], threshold)
    cxx_components = cpp_components(paths[4], threshold)
    receipt = json.loads(paths[5].read_text(encoding="utf-8"))
    if receipt.get("schema_name") != "intersubmod.strict_endpoint_graph_receipt":
        raise ValueError("C++ receipt schema mismatch")
    edge_mismatch = mismatch_preview(py_edges, cxx_edges)
    component_mismatch = mismatch_preview(py_components, cxx_components)
    checks = {
        "edge_keys_and_counts_identical": py_edges == cxx_edges,
        "component_partitions_and_roles_identical": py_components == cxx_components,
        "cpp_threshold_matches": receipt.get("threshold") == threshold,
        "cpp_singleton_abstain_contract": receipt.get("invariants", {}).get(
            "singleton_components_are_abstain_not_tree_regions"
        )
        is True,
        "edge_cardinality_matches_receipt": len(cxx_edges)
        == int(receipt.get("counts", {}).get("observed_edges", -1)),
        "component_cardinality_matches_receipt": len(cxx_components)
        == int(receipt.get("counts", {}).get("components", -1)),
    }
    document = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "all_pass": all(checks.values()),
        "threshold": threshold,
        "counts": {
            "python_edges": len(py_edges),
            "cpp_edges": len(cxx_edges),
            "python_components": len(py_components),
            "cpp_components": len(cxx_components),
            "edge_mismatch_preview_count": len(edge_mismatch),
            "component_mismatch_preview_count": len(component_mismatch),
        },
        "checks": checks,
        "edge_mismatch_preview": edge_mismatch,
        "component_mismatch_preview": component_mismatch,
        "inputs": {path.name: identity(path) for path in paths},
    }
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(document, handle, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    return document


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python-edges", required=True, type=Path)
    parser.add_argument("--python-components", required=True, type=Path)
    parser.add_argument("--python-membership", required=True, type=Path)
    parser.add_argument("--cpp-edges", required=True, type=Path)
    parser.add_argument("--cpp-components", required=True, type=Path)
    parser.add_argument("--cpp-receipt", required=True, type=Path)
    parser.add_argument("--threshold", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = compare(
        python_edges=args.python_edges,
        python_components_path=args.python_components,
        python_membership=args.python_membership,
        cpp_edges=args.cpp_edges,
        cpp_components_path=args.cpp_components,
        cpp_receipt=args.cpp_receipt,
        threshold=args.threshold,
        output=args.output,
    )
    print(json.dumps({"all_pass": result["all_pass"], "counts": result["counts"]}, sort_keys=True))
    return 0 if result["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
