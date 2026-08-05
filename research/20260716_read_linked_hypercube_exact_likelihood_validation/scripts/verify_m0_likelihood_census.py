#!/usr/bin/env python3
"""Independent, read-only verifier for the M0 likelihood census.

This verifier intentionally does not import ``run_m0_likelihood_census.py``.
It recomputes the receipt aggregates from the emitted TSV, reconciles every
TSV unit against the canonical-v5 lineage population, and can independently
reconstruct the candidate vertex-set classes used by the census.

The verifier does not refit the likelihood optimizer.  Instead, it verifies
the persisted likelihood ordering fields, status contract, candidate IDs, and
all population/aggregate conservation identities.  This keeps it suitable for
checking the full seven-dataset output without rerunning the expensive census.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import platform
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import scipy


DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)

EXPECTED_ROW_FIELDS = (
    "dataset",
    "region",
    "chrom",
    "start",
    "end",
    "k",
    "family",
    "n_reads_reported",
    "n_scoring_pattern_groups",
    "n_scoring_alignment_exposures",
    "raw_tree_candidates_T",
    "distinct_vertex_sets_V",
    "distinct_edge_sets_E",
    "candidate_generation_status",
    "best_vertex_sets",
    "top_edge_variants",
    "best_log_likelihood",
    "second_log_likelihood",
    "delta_best_second",
    "top_relative_likelihood_weight",
    "selection_status",
    "all_fits_converged",
    "all_fits_monotone",
    "max_emission_rank",
    "vertex_set_ids",
    "top_vertex_set_ids",
)

SELECTION_STATUSES = {
    "T1_CANDIDATE_UNIQUE",
    "T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED",
    "LIKELIHOOD_TIED_VERTEX_SETS",
    "LIKELIHOOD_UNIQUE_VERTEX_EDGE_UNRESOLVED",
    "LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE",
    "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE",
}

GENERATION_STATUSES = {
    "STORED_COMPLETE",
    "RECONSTRUCTED_FROM_FROZEN_SOLVER",
}

HEX64 = re.compile(r"^[0-9a-f]{64}$")


def sha256_path(path: Path, block_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def file_binding(path: Path, *, hash_scope: str = "exact file bytes") -> dict[str, Any]:
    """Return an explicit byte-identity binding for a local file."""
    resolved = path.resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_path(resolved),
        "hash_scope": hash_scope,
    }


def verification_provenance(receipt_path: Path, rows_path: Path) -> dict[str, Any]:
    """Bind the verifier result to inputs, source files, and verifier runtime."""
    verifier_path = Path(__file__).resolve()
    scripts_dir = verifier_path.parent
    return {
        "byte_identity": {
            "m0_receipt": file_binding(receipt_path),
            "m0_rows_gzip": file_binding(
                rows_path,
                hash_scope="exact compressed file bytes, including the gzip container",
            ),
            "m0_census_script_on_disk_at_verification": file_binding(
                scripts_dir / "run_m0_likelihood_census.py"
            ),
            "likelihood_scoring_utility_on_disk_at_verification": file_binding(
                scripts_dir / "hypercube_exact.py"
            ),
            "independent_verifier": file_binding(verifier_path),
        },
        "verification_runtime": {
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "python_executable": sys.executable,
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
            "platform": platform.platform(),
        },
        "interpretation": (
            "source-file hashes and software versions describe the verification-time environment; "
            "because M0 receipt 1.0.0 did not persist them, they do not by themselves prove the "
            "historical census process loaded identical source bytes/runtime"
        ),
    }


def vertex_set_digest(k: int, vertices: Iterable[int]) -> str:
    payload = {"k": int(k), "vertices": sorted(int(vertex) for vertex in vertices)}
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def parse_full(pattern: str) -> int:
    if not pattern or set(pattern) - {"R", "A"}:
        raise ValueError(f"invalid full state: {pattern!r}")
    state = 0
    for bit, symbol in enumerate(pattern):
        if symbol == "A":
            state |= 1 << bit
    return state


def node_to_vertex(node: str, k: int) -> int:
    if node == "ROOT":
        return 0
    pattern = node[2:] if node.startswith("H_") else node
    if len(pattern) != k:
        raise ValueError(f"tree node length mismatch: {node!r}, k={k}")
    return parse_full(pattern)


def tree_vertex_set(tree: dict[str, Any], k: int) -> tuple[int, ...]:
    vertices = {0}
    for edge in tree.get("edges") or []:
        if not isinstance(edge, list) or len(edge) != 2:
            raise ValueError(f"malformed edge: {edge!r}")
        vertices.add(node_to_vertex(str(edge[0]), k))
        vertices.add(node_to_vertex(str(edge[1]), k))
    return tuple(sorted(vertices))


def edge_set_digest(tree: dict[str, Any]) -> str:
    edges = sorted((str(parent), str(child)) for parent, child in (tree.get("edges") or []))
    payload = {
        "edges": edges,
        "recurrence": sorted(str(value) for value in (tree.get("recurrence") or [])),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def parse_bool(value: str, *, label: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"{label}: expected true/false, found {value!r}")


def parse_int(value: str, *, label: str, minimum: int | None = None) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: invalid integer {value!r}") from exc
    if minimum is not None and parsed < minimum:
        raise ValueError(f"{label}: {parsed} < {minimum}")
    return parsed


def parse_optional_float(value: str, *, label: str) -> float | None:
    if value == "":
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: invalid float {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{label}: non-finite float {value!r}")
    return parsed


def parse_id_list(value: str, *, label: str) -> list[str]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label}: invalid JSON") from exc
    if not isinstance(parsed, list) or any(not isinstance(item, str) for item in parsed):
        raise ValueError(f"{label}: expected a JSON string list")
    if any(not HEX64.fullmatch(item) for item in parsed):
        raise ValueError(f"{label}: contains a non-SHA256 identifier")
    if len(set(parsed)) != len(parsed):
        raise ValueError(f"{label}: duplicate identifiers")
    return parsed


def expected_selection_status(row: dict[str, Any]) -> str:
    if not row["all_fits_converged"]:
        return "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE"
    if row["raw_tree_candidates_T"] == 1:
        return "T1_CANDIDATE_UNIQUE"
    if row["distinct_vertex_sets_V"] == 1:
        return "T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED"
    if row["best_vertex_sets"] > 1:
        return "LIKELIHOOD_TIED_VERTEX_SETS"
    if row["top_edge_variants"] > 1:
        return "LIKELIHOOD_UNIQUE_VERTEX_EDGE_UNRESOLVED"
    return "LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE"


def normalize_row(raw: dict[str, str], row_number: int, errors: list[str]) -> dict[str, Any] | None:
    prefix = f"TSV row {row_number}"
    try:
        row: dict[str, Any] = {
            "dataset": raw["dataset"],
            "region": raw["region"],
            "chrom": raw["chrom"],
            "family": raw["family"],
            "candidate_generation_status": raw["candidate_generation_status"],
            "selection_status": raw["selection_status"],
        }
        for key, minimum in (
            ("start", 0),
            ("end", 0),
            ("k", 1),
            ("n_reads_reported", 1),
            ("n_scoring_pattern_groups", 1),
            ("n_scoring_alignment_exposures", 1),
            ("raw_tree_candidates_T", 1),
            ("distinct_vertex_sets_V", 1),
            ("distinct_edge_sets_E", 1),
            ("best_vertex_sets", 1),
            ("top_edge_variants", 1),
            ("max_emission_rank", 0),
        ):
            row[key] = parse_int(raw[key], label=f"{prefix} {key}", minimum=minimum)
        for key in (
            "best_log_likelihood",
            "second_log_likelihood",
            "delta_best_second",
            "top_relative_likelihood_weight",
        ):
            row[key] = parse_optional_float(raw[key], label=f"{prefix} {key}")
        row["all_fits_converged"] = parse_bool(
            raw["all_fits_converged"], label=f"{prefix} all_fits_converged"
        )
        row["all_fits_monotone"] = parse_bool(
            raw["all_fits_monotone"], label=f"{prefix} all_fits_monotone"
        )
        row["vertex_set_ids"] = parse_id_list(raw["vertex_set_ids"], label=f"{prefix} vertex_set_ids")
        row["top_vertex_set_ids"] = parse_id_list(
            raw["top_vertex_set_ids"], label=f"{prefix} top_vertex_set_ids"
        )
    except (KeyError, ValueError) as exc:
        errors.append(str(exc))
        return None

    if row["dataset"] not in DATASETS:
        errors.append(f"{prefix}: unknown dataset {row['dataset']!r}")
    if row["family"] not in {"1", "2"}:
        errors.append(f"{prefix}: family must be 1 or 2")
    if row["start"] > row["end"]:
        errors.append(f"{prefix}: start exceeds end")
    if row["selection_status"] not in SELECTION_STATUSES:
        errors.append(f"{prefix}: unknown selection_status {row['selection_status']!r}")
    if row["candidate_generation_status"] not in GENERATION_STATUSES:
        errors.append(f"{prefix}: unknown candidate_generation_status {row['candidate_generation_status']!r}")
    if row["n_scoring_alignment_exposures"] > row["n_reads_reported"]:
        errors.append(f"{prefix}: scoring alignment exposures exceed reported reads")

    raw_t = row["raw_tree_candidates_T"]
    n_vertex = row["distinct_vertex_sets_V"]
    n_edge = row["distinct_edge_sets_E"]
    best_vertex = row["best_vertex_sets"]
    top_edges = row["top_edge_variants"]
    vertex_ids = row["vertex_set_ids"]
    top_ids = row["top_vertex_set_ids"]
    if not (1 <= n_vertex <= raw_t):
        errors.append(f"{prefix}: expected 1 <= V <= T")
    if n_edge != raw_t:
        errors.append(f"{prefix}: E ({n_edge}) != T ({raw_t})")
    if not (1 <= best_vertex <= n_vertex):
        errors.append(f"{prefix}: expected 1 <= best_vertex_sets <= V")
    if not (best_vertex <= top_edges <= raw_t):
        errors.append(f"{prefix}: expected best_vertex_sets <= top_edge_variants <= T")
    if len(vertex_ids) != n_vertex:
        errors.append(f"{prefix}: vertex_set_ids length != V")
    if len(top_ids) != best_vertex:
        errors.append(f"{prefix}: top_vertex_set_ids length != best_vertex_sets")
    if not set(top_ids).issubset(vertex_ids):
        errors.append(f"{prefix}: top vertex IDs are not a subset of all vertex IDs")

    weight = row["top_relative_likelihood_weight"]
    if weight is None or not 0.0 < weight <= 1.0 + 1e-12:
        errors.append(f"{prefix}: top relative likelihood weight must be in (0,1]")
    elif weight * best_vertex > 1.0 + 2e-6:
        errors.append(f"{prefix}: tied top candidates imply an impossible normalized weight")

    best_ll = row["best_log_likelihood"]
    second_ll = row["second_log_likelihood"]
    delta = row["delta_best_second"]
    if n_vertex == 1:
        if any(value is not None for value in (best_ll, second_ll, delta)):
            errors.append(f"{prefix}: V=1 must leave likelihood comparison cells blank")
        if row["max_emission_rank"] != 0:
            errors.append(f"{prefix}: V=1 must have max_emission_rank=0")
        if weight is not None and not math.isclose(weight, 1.0, abs_tol=1e-12):
            errors.append(f"{prefix}: V=1 must have relative weight 1")
    else:
        if any(value is None for value in (best_ll, second_ll, delta)):
            errors.append(f"{prefix}: V>1 requires best/second/delta likelihood values")
        elif best_ll is not None and second_ll is not None and delta is not None:
            if best_ll + 1e-8 < second_ll:
                errors.append(f"{prefix}: best log likelihood is below second")
            if not math.isclose(delta, best_ll - second_ll, rel_tol=2e-9, abs_tol=2e-9):
                errors.append(f"{prefix}: delta_best_second != best - second")

    expected_status = expected_selection_status(row)
    if row["selection_status"] != expected_status:
        errors.append(
            f"{prefix}: status {row['selection_status']} != independently derived {expected_status}"
        )
    return row


def read_rows(path: Path, errors: list[str]) -> tuple[list[dict[str, Any]], tuple[str, ...]]:
    rows: list[dict[str, Any]] = []
    keys: set[tuple[str, str, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = tuple(reader.fieldnames or ())
        if fields != EXPECTED_ROW_FIELDS:
            errors.append(f"TSV schema mismatch: {fields!r}")
        for row_number, raw in enumerate(reader, start=2):
            normalized = normalize_row(raw, row_number, errors)
            if normalized is None:
                continue
            key = (normalized["dataset"], normalized["region"], normalized["family"])
            if key in keys:
                errors.append(f"TSV duplicate unit key: {key}")
            keys.add(key)
            rows.append(normalized)
    return rows, fields


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    status = Counter(row["selection_status"] for row in rows)
    by_dataset: dict[str, dict[str, Any]] = {}
    for dataset in DATASETS:
        subset = [row for row in rows if row["dataset"] == dataset]
        by_dataset[dataset] = {
            "n_hp_lineage_units": len(subset),
            "n_regions": len({row["region"] for row in subset}),
            "selection_status_counts": dict(Counter(row["selection_status"] for row in subset)),
            "raw_tree_candidates_T": sum(row["raw_tree_candidates_T"] for row in subset),
            "distinct_vertex_sets_V": sum(row["distinct_vertex_sets_V"] for row in subset),
        }
    n = len(rows)
    return {
        "n_hp_lineage_units": n,
        "n_regions_with_at_least_one_eligible_hp_lineage": len(
            {(row["dataset"], row["region"]) for row in rows}
        ),
        "raw_tree_candidates_T": sum(row["raw_tree_candidates_T"] for row in rows),
        "distinct_vertex_sets_V": sum(row["distinct_vertex_sets_V"] for row in rows),
        "selection_status_counts": dict(status),
        "selection_status_fraction_of_eligible_hp_units": {
            key: value / n if n else None for key, value in sorted(status.items())
        },
        "all_fits_converged": all(row["all_fits_converged"] for row in rows),
        "all_fits_monotone": all(row["all_fits_monotone"] for row in rows),
        "n_optimizer_abstain": sum(
            row["selection_status"] == "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE" for row in rows
        ),
        "by_dataset": by_dataset,
    }


def categorical_conservation(rows: list[dict[str, Any]]) -> dict[str, Any]:
    t1_v1 = sum(
        row["raw_tree_candidates_T"] == 1 and row["distinct_vertex_sets_V"] == 1
        for row in rows
    )
    t_gt1_v1 = sum(
        row["raw_tree_candidates_T"] > 1 and row["distinct_vertex_sets_V"] == 1
        for row in rows
    )
    t_gt1_v_gt1 = sum(
        row["raw_tree_candidates_T"] > 1 and row["distinct_vertex_sets_V"] > 1
        for row in rows
    )
    t1_v_gt1 = sum(
        row["raw_tree_candidates_T"] == 1 and row["distinct_vertex_sets_V"] > 1
        for row in rows
    )
    unique_statuses = {
        "LIKELIHOOD_UNIQUE_VERTEX_EDGE_UNRESOLVED",
        "LIKELIHOOD_UNIQUE_VERTEX_SINGLE_EDGE_CANDIDATE",
    }
    likelihood_unique = sum(row["selection_status"] in unique_statuses for row in rows)
    likelihood_tied = sum(
        row["selection_status"] == "LIKELIHOOD_TIED_VERTEX_SETS" for row in rows
    )
    nonconverged = sum(
        row["selection_status"] == "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE" for row in rows
    )
    vertex_equivalent = sum(
        row["selection_status"] == "T_GT1_VERTEX_EQUIVALENT_EDGE_UNRESOLVED" for row in rows
    )
    t1_status = sum(row["selection_status"] == "T1_CANDIDATE_UNIQUE" for row in rows)
    return {
        "T_eq_1_V_eq_1": t1_v1,
        "T_eq_1_V_gt_1_impossible": t1_v_gt1,
        "T_gt_1_V_eq_1": t_gt1_v1,
        "T_gt_1_V_gt_1": t_gt1_v_gt1,
        "likelihood_unique_vertex_set": likelihood_unique,
        "likelihood_tied_vertex_sets": likelihood_tied,
        "optimizer_nonconverged_abstain": nonconverged,
        "T_gt_1_vertex_equivalent_edge_unresolved": vertex_equivalent,
        "T1_status": t1_status,
        "partition_sum": t1_status + vertex_equivalent + likelihood_unique + likelihood_tied + nonconverged,
        "n_hp_lineage_units": len(rows),
        "partition_conserves": (
            t1_status + vertex_equivalent + likelihood_unique + likelihood_tied + nonconverged
            == len(rows)
        ),
        "T_V_partition_conserves": t1_v1 + t1_v_gt1 + t_gt1_v1 + t_gt1_v_gt1 == len(rows),
    }


def eligible_reason(lineage: dict[str, Any]) -> tuple[bool, str]:
    if not lineage.get("is_primary_lineage"):
        return False, "NOT_PRIMARY_HP1_HP2"
    if not lineage.get("mutation_bearing"):
        return False, "REFERENCE_ONLY"
    if str(lineage.get("family")) not in {"1", "2"}:
        return False, "NOT_GERMLINE_HP1_HP2"
    if lineage.get("capped"):
        return False, "CAPPED"
    if lineage.get("analysis_candidate_set_complete") is not True:
        return False, "CANDIDATE_SET_INCOMPLETE"
    if lineage.get("verification_complete") is not True or lineage.get("verify_pass") is not True:
        return False, "SOLVER_VERIFICATION_NOT_FULL_PASS"
    if not lineage.get("obs_subreads"):
        return False, "NO_THRESHOLD_RETAINED_PATTERNS"
    return True, "ELIGIBLE_M0"


def load_frozen_solver(canonical_root: Path):
    path = canonical_root / "source_bundle" / "files" / "imported" / "003_tree_enumeration_solver.py"
    if not path.is_file():
        raise FileNotFoundError(path)
    specification = importlib.util.spec_from_file_location("m0_independent_frozen_solver", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot import frozen solver: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module, path


def candidate_classes(
    dataset: str,
    region: dict[str, Any],
    lineage: dict[str, Any],
    frozen_solver,
) -> tuple[dict[tuple[int, ...], int], str]:
    k = int(region["n_sSNV"])
    raw_t = int(lineage.get("n_trees") or 0)
    if lineage.get("display_trees_complete") is True:
        trees = list(lineage.get("trees") or [])
        if raw_t != len(trees) or int(lineage.get("n_trees_stored") or 0) != raw_t:
            raise ValueError(f"stored T mismatch: {dataset} {region['region']} HP{lineage['family']}")
        edge_ids: dict[tuple[int, ...], set[str]] = defaultdict(set)
        for tree in trees:
            edge_ids[tree_vertex_set(tree, k)].add(edge_set_digest(tree))
        if sum(len(values) for values in edge_ids.values()) != raw_t:
            raise ValueError(f"stored edge candidates do not conserve T: {dataset} {region['region']}")
        return {vertices: len(values) for vertices, values in edge_ids.items()}, "STORED_COMPLETE"

    full = lineage.get("obs_populations") or {}
    partial = list((lineage.get("obs_subreads") or {}).keys())
    rebuilt = frozen_solver.enumerate_min_trees(full, partial, k, tree_cap=1)
    if rebuilt.get("capped"):
        raise ValueError(f"independent frozen reconstruction capped: {dataset} {region['region']}")
    if int(rebuilt.get("n_trees") or 0) != raw_t:
        raise ValueError(f"reconstructed T mismatch: {dataset} {region['region']}")
    if int(rebuilt.get("n_hidden") or 0) != int(lineage.get("n_hidden") or 0):
        raise ValueError(f"reconstructed hidden mismatch: {dataset} {region['region']}")
    classes: dict[tuple[int, ...], int] = {}
    for node_set in rebuilt.get("_feasible_N") or []:
        vertices = tuple(sorted(sum(1 << bit for bit in node) for node in node_set))
        edge_count = 1
        for node in node_set:
            if node:
                edge_count *= sum(1 for bit in node if (node - {bit}) in node_set)
        if vertices in classes:
            raise ValueError(f"duplicate reconstructed vertex set: {dataset} {region['region']}")
        classes[vertices] = edge_count
    if sum(classes.values()) != raw_t:
        raise ValueError(f"reconstructed edge classes do not conserve T: {dataset} {region['region']}")
    return classes, "RECONSTRUCTED_FROM_FROZEN_SOLVER"


def reconcile_canonical(
    rows: list[dict[str, Any]],
    receipt: dict[str, Any],
    canonical_root: Path,
    candidate_mode: str,
    errors: list[str],
) -> dict[str, Any]:
    selected = tuple(receipt.get("selected_datasets") or ())
    row_map = {(row["dataset"], row["region"], row["family"]): row for row in rows}
    canonical_eligible: dict[tuple[str, str, str], tuple[dict[str, Any], dict[str, Any]]] = {}
    primary_regions: set[tuple[str, str]] = set()
    fully_eligible_regions: set[tuple[str, str]] = set()
    exclusion = Counter()
    input_evidence = []
    n_primary_hp_lineage_units = 0
    frozen_solver = None
    frozen_solver_path = canonical_root / "source_bundle" / "files" / "imported" / "003_tree_enumeration_solver.py"
    if candidate_mode == "deep":
        frozen_solver, frozen_solver_path = load_frozen_solver(canonical_root)

    for dataset in selected:
        path = canonical_root / "samples" / dataset / f"layered_region_view_{dataset}.json"
        if not path.is_file():
            errors.append(f"missing canonical input: {path}")
            continue
        evidence = {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_path(path)}
        input_evidence.append(evidence)
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("sample") != dataset:
            errors.append(f"canonical sample identity mismatch: {path}")
        for region in payload.get("regions") or []:
            primary = [
                lineage
                for lineage in (region.get("lineages") or [])
                if lineage.get("is_primary_lineage") and lineage.get("mutation_bearing")
            ]
            if not primary:
                continue
            region_key = (dataset, str(region["region"]))
            primary_regions.add(region_key)
            n_primary_hp_lineage_units += len(primary)
            eligibility = [eligible_reason(lineage) for lineage in primary]
            if all(ok for ok, _ in eligibility):
                fully_eligible_regions.add(region_key)
            for lineage, (eligible, reason) in zip(primary, eligibility):
                if not eligible:
                    exclusion[reason] += 1
                    continue
                key = (dataset, str(region["region"]), str(lineage["family"]))
                if key in canonical_eligible:
                    errors.append(f"duplicate canonical eligible unit: {key}")
                canonical_eligible[key] = (region, lineage)

    missing_rows = sorted(set(canonical_eligible) - set(row_map))
    extra_rows = sorted(set(row_map) - set(canonical_eligible))
    if missing_rows:
        errors.append(f"TSV is missing {len(missing_rows)} canonical eligible units; first={missing_rows[:3]}")
    if extra_rows:
        errors.append(f"TSV has {len(extra_rows)} noneligible/unknown units; first={extra_rows[:3]}")

    deep_checked = 0
    stored_checked = 0
    reconstructed_checked = 0
    for key in sorted(set(row_map) & set(canonical_eligible)):
        row = row_map[key]
        region, lineage = canonical_eligible[key]
        expected_scalar = {
            "chrom": str(region["chrom"]),
            "start": int(region["start"]),
            "end": int(region["end"]),
            "k": int(region["n_sSNV"]),
            "n_reads_reported": int(lineage.get("n_reads") or 0),
            "n_scoring_pattern_groups": len(lineage.get("obs_subreads") or {}),
            "n_scoring_alignment_exposures": sum(
                int(value) for value in (lineage.get("obs_subreads") or {}).values()
            ),
            "raw_tree_candidates_T": int(lineage.get("n_trees") or 0),
            "distinct_edge_sets_E": int(lineage.get("n_trees") or 0),
        }
        for field, expected in expected_scalar.items():
            if row[field] != expected:
                errors.append(f"{key}: {field} TSV={row[field]!r}, canonical={expected!r}")

        expected_generation = (
            "STORED_COMPLETE"
            if lineage.get("display_trees_complete") is True
            else "RECONSTRUCTED_FROM_FROZEN_SOLVER"
        )
        if row["candidate_generation_status"] != expected_generation:
            errors.append(
                f"{key}: generation status {row['candidate_generation_status']} != {expected_generation}"
            )

        if candidate_mode != "deep":
            continue
        try:
            classes, generation = candidate_classes(dataset=key[0], region=region, lineage=lineage, frozen_solver=frozen_solver)
        except Exception as exc:  # keep a full audit instead of stopping at the first bad unit
            errors.append(f"{key}: candidate reconstruction failed: {exc}")
            continue
        deep_checked += 1
        if generation == "STORED_COMPLETE":
            stored_checked += 1
        else:
            reconstructed_checked += 1
        expected_ids = [vertex_set_digest(row["k"], vertices) for vertices in sorted(classes)]
        if row["distinct_vertex_sets_V"] != len(classes):
            errors.append(f"{key}: V={row['distinct_vertex_sets_V']} != reconstructed {len(classes)}")
        if row["vertex_set_ids"] != expected_ids:
            errors.append(f"{key}: ordered vertex_set_ids differ from independent reconstruction")
        top_edge_count = sum(
            edge_count
            for vertices, edge_count in classes.items()
            if vertex_set_digest(row["k"], vertices) in set(row["top_vertex_set_ids"])
        )
        if row["top_edge_variants"] != top_edge_count:
            errors.append(
                f"{key}: top_edge_variants={row['top_edge_variants']} != reconstructed {top_edge_count}"
            )

    declared_inputs = receipt.get("input_files") or []
    if declared_inputs != input_evidence:
        errors.append("receipt input_files metadata/hash differs from canonical files")
    declared_solver = receipt.get("frozen_solver") or {}
    if not frozen_solver_path.is_file():
        errors.append(f"missing frozen solver: {frozen_solver_path}")
        solver_sha = None
    else:
        solver_sha = sha256_path(frozen_solver_path)
        if declared_solver != {"path": str(frozen_solver_path), "sha256": solver_sha}:
            errors.append("receipt frozen_solver path/hash mismatch")

    return {
        "n_primary_mutation_regions": len(primary_regions),
        "n_primary_mutation_hp_lineage_units": n_primary_hp_lineage_units,
        "n_fully_m0_eligible_regions": len(fully_eligible_regions),
        "n_regions_with_any_eligible_hp_lineage": len(
            {(dataset, region) for dataset, region, _family in canonical_eligible}
        ),
        "n_eligible_hp_lineage_units": len(canonical_eligible),
        "excluded_hp_lineage_unit_counts": dict(exclusion),
        "eligible_plus_excluded_equals_primary_hp_units": (
            len(canonical_eligible) + sum(exclusion.values()) == n_primary_hp_lineage_units
        ),
        "missing_tsv_units": len(missing_rows),
        "extra_tsv_units": len(extra_rows),
        "candidate_mode": candidate_mode,
        "n_candidate_units_deep_checked": deep_checked,
        "n_stored_candidate_units_deep_checked": stored_checked,
        "n_reconstructed_candidate_units_deep_checked": reconstructed_checked,
        "input_files": input_evidence,
        "frozen_solver_sha256": solver_sha,
    }


def compare_exact(label: str, actual: Any, expected: Any, errors: list[str]) -> None:
    if actual != expected:
        errors.append(f"{label} mismatch: actual={actual!r}, expected={expected!r}")


def compare_aggregate(
    receipt_aggregate: dict[str, Any],
    computed: dict[str, Any],
    errors: list[str],
) -> None:
    for key in (
        "n_hp_lineage_units",
        "n_regions_with_at_least_one_eligible_hp_lineage",
        "raw_tree_candidates_T",
        "distinct_vertex_sets_V",
        "selection_status_counts",
        "all_fits_converged",
        "all_fits_monotone",
        "n_optimizer_abstain",
        "by_dataset",
    ):
        compare_exact(f"receipt aggregate.{key}", receipt_aggregate.get(key), computed[key], errors)
    actual_fractions = receipt_aggregate.get("selection_status_fraction_of_eligible_hp_units") or {}
    expected_fractions = computed["selection_status_fraction_of_eligible_hp_units"]
    if set(actual_fractions) != set(expected_fractions):
        errors.append("receipt aggregate selection fraction keys mismatch")
    else:
        for key in expected_fractions:
            actual = actual_fractions[key]
            expected = expected_fractions[key]
            if actual is None or expected is None:
                if actual != expected:
                    errors.append(f"receipt aggregate fraction {key} mismatch")
            elif not math.isclose(float(actual), float(expected), rel_tol=0.0, abs_tol=1e-15):
                errors.append(f"receipt aggregate fraction {key} mismatch: {actual} != {expected}")


def resolve_rows_path(output_dir: Path, receipt: dict[str, Any]) -> Path:
    expected = output_dir / "m0_hp_lineage_likelihood_census.tsv.gz"
    if expected.is_file():
        return expected
    declared = Path(str((receipt.get("outputs") or {}).get("rows") or ""))
    if declared.is_absolute() and declared.is_file():
        return declared
    if declared.is_file():
        return declared.resolve()
    raise FileNotFoundError(expected)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--canonical-root", type=Path)
    parser.add_argument(
        "--candidate-mode",
        choices=("metadata", "deep"),
        default="deep",
        help="deep independently reconstructs every candidate vertex-set class",
    )
    parser.add_argument("--json-output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    receipt_path = output_dir / "m0_receipt.json"
    if not receipt_path.is_file():
        raise FileNotFoundError(receipt_path)
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    rows_path = resolve_rows_path(output_dir, receipt)
    canonical_root = (args.canonical_root or Path(str(receipt.get("canonical_root") or ""))).resolve()
    errors: list[str] = []
    warnings: list[str] = []

    compare_exact(
        "schema_name",
        receipt.get("schema_name"),
        "intersubmod.hypercube_m0_likelihood_census_receipt",
        errors,
    )
    compare_exact("schema_version", receipt.get("schema_version"), "1.0.0", errors)
    selected = tuple(receipt.get("selected_datasets") or ())
    if not selected or len(set(selected)) != len(selected) or any(dataset not in DATASETS for dataset in selected):
        errors.append(f"invalid selected_datasets: {selected}")
    if Path(str(receipt.get("canonical_root") or "")).resolve() != canonical_root:
        errors.append("receipt canonical_root differs from verification canonical_root")

    declared_output = receipt.get("outputs") or {}
    actual_rows_size = rows_path.stat().st_size
    actual_rows_sha = sha256_path(rows_path)
    compare_exact("receipt outputs.rows_size_bytes", declared_output.get("rows_size_bytes"), actual_rows_size, errors)
    compare_exact("receipt outputs.rows_sha256", declared_output.get("rows_sha256"), actual_rows_sha, errors)

    rows, fields = read_rows(rows_path, errors)
    computed_aggregate = summarize_rows(rows)
    compare_aggregate(receipt.get("aggregate") or {}, computed_aggregate, errors)
    categories = categorical_conservation(rows)
    if not categories["partition_conserves"] or not categories["T_V_partition_conserves"]:
        errors.append("mutually exclusive status or T/V partition does not conserve eligible units")
    if categories["T_eq_1_V_gt_1_impossible"] != 0:
        errors.append("found impossible T=1 and V>1 rows")

    canonical = reconcile_canonical(
        rows=rows,
        receipt=receipt,
        canonical_root=canonical_root,
        candidate_mode=args.candidate_mode,
        errors=errors,
    )
    population = receipt.get("population") or {}
    for key in (
        "n_primary_mutation_regions",
        "n_fully_m0_eligible_regions",
        "n_regions_with_any_eligible_hp_lineage",
        "excluded_hp_lineage_unit_counts",
    ):
        compare_exact(f"receipt population.{key}", population.get(key), canonical[key], errors)
    compare_exact(
        "canonical eligible units vs TSV rows",
        canonical["n_eligible_hp_lineage_units"],
        len(rows),
        errors,
    )
    if not canonical["eligible_plus_excluded_equals_primary_hp_units"]:
        errors.append("canonical eligible + excluded HP units does not equal primary mutation HP units")

    # Reproduce receipt 1.0.0 fields for tamper/schema consistency only.  In
    # particular, the first legacy field is not accepted as a real invariant;
    # the strong canonical equality is enforced separately above and reported
    # as ``strong_eligible_excluded_conservation`` below.
    legacy_receipt_conservation = {
        "eligible_plus_excluded_hp_units_positive": bool(rows)
        and sum(canonical["excluded_hp_lineage_unit_counts"].values()) >= 0,
        "output_rows_equal_eligible_hp_units": len(rows) == computed_aggregate["n_hp_lineage_units"],
        "selected_datasets_present": all(
            computed_aggregate["by_dataset"][dataset]["n_hp_lineage_units"] > 0 for dataset in selected
        ),
        "all_nonconverged_units_fail_closed": computed_aggregate["n_optimizer_abstain"]
        == computed_aggregate["selection_status_counts"].get(
            "RANK_ABSTAIN_OPTIMIZER_NONCONVERGENCE", 0
        ),
        "all_fits_monotone": computed_aggregate["all_fits_monotone"],
    }
    compare_exact(
        "legacy receipt self-consistency fields",
        receipt.get("conservation"),
        legacy_receipt_conservation,
        errors,
    )
    compare_exact(
        "legacy receipt all_pass self-consistency",
        receipt.get("all_pass"),
        all(legacy_receipt_conservation.values()),
        errors,
    )
    expected_full_scope = set(selected) == set(DATASETS)
    compare_exact("receipt full_task_b_scope", receipt.get("full_task_b_scope"), expected_full_scope, errors)

    if "executing_script" not in receipt and "script_sha256" not in receipt:
        warnings.append(
            "receipt 1.0.0 does not record the census script SHA; input/frozen-solver/rows are hashed, "
            "but the exact scoring implementation is not self-contained in the receipt"
        )
    if "software_versions" not in receipt:
        warnings.append(
            "receipt 1.0.0 does not record Python/NumPy/SciPy versions; optimizer reproducibility is incomplete"
        )
    warnings.append(
        "receipt 1.0.0 field eligible_plus_excluded_hp_units_positive has a naming/semantic weakness: "
        "it is reproduced only as a legacy receipt field, not treated as a real checked invariant; "
        "this verifier separately enforces eligible + excluded == canonical primary mutation HP-lineage units"
    )
    warnings.append(
        "this verifier checks persisted ordering/status consistency but deliberately does not rerun all SLSQP fits"
    )

    provenance = verification_provenance(receipt_path, rows_path)
    result = {
        "schema_name": "intersubmod.hypercube_m0_independent_verification",
        "schema_version": "1.1.0",
        "verdict": "PASS" if not errors else "FAIL",
        "candidate_mode": args.candidate_mode,
        "inputs": {
            "receipt": str(receipt_path),
            "receipt_size_bytes": receipt_path.stat().st_size,
            "receipt_sha256": sha256_path(receipt_path),
            "rows": str(rows_path),
            "rows_size_bytes": actual_rows_size,
            "rows_sha256": actual_rows_sha,
            "canonical_root": str(canonical_root),
        },
        "verification_provenance": provenance,
        "legacy_receipt_self_consistency": {
            "fields_match_original_implementation": receipt.get("conservation")
            == legacy_receipt_conservation,
            "all_pass_matches_original_implementation": receipt.get("all_pass")
            == all(legacy_receipt_conservation.values()),
            "eligible_plus_excluded_positive_is_not_a_strong_invariant": True,
        },
        "scope": {
            "selected_datasets": list(selected),
            "full_task_b_scope": expected_full_scope,
            "row_schema_matches": fields == EXPECTED_ROW_FIELDS,
        },
        "independently_recomputed_aggregate": computed_aggregate,
        "categorical_conservation": categories,
        "canonical_reconciliation": canonical,
        "checks": {
            "n_errors": len(errors),
            "receipt_tsv_hash_matches": declared_output.get("rows_sha256") == actual_rows_sha,
            "receipt_aggregate_matches_tsv": not any(error.startswith("receipt aggregate") for error in errors),
            "tsv_units_match_canonical_eligible_population": canonical["missing_tsv_units"] == 0
            and canonical["extra_tsv_units"] == 0,
            "strong_eligible_excluded_conservation": canonical[
                "eligible_plus_excluded_equals_primary_hp_units"
            ],
            "selection_partition_conserves": categories["partition_conserves"],
            "T_V_partition_conserves": categories["T_V_partition_conserves"],
            "required_files_are_exact_byte_hashed": all(
                item.get("sha256") and item.get("size_bytes", -1) >= 0
                for item in provenance["byte_identity"].values()
            ),
            "verification_runtime_versions_recorded": all(
                provenance["verification_runtime"].get(key)
                for key in ("python_version", "numpy_version", "scipy_version")
            ),
        },
        "warnings": warnings,
        "errors": errors,
    }
    encoded = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        if args.json_output.exists():
            raise FileExistsError(f"refusing to overwrite {args.json_output}")
        args.json_output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    raise SystemExit(0 if not errors else 1)


if __name__ == "__main__":
    main()
