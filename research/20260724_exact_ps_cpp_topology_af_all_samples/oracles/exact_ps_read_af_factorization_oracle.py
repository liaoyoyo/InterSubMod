#!/usr/bin/env python3
"""Replay the frozen Python read-AF ranking by explicit and factorized routes.

This oracle is intentionally read-only with respect to MLHP and layered inputs.
It writes one machine-readable receipt and binds that receipt to input/source
SHA-256 identities.  The exact-PS AF here is the supported-pattern read-AF
stored in ``col_coverage_by_hp``; it is not caller VAF or a calibrated
likelihood.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from types import ModuleType
from typing import Any


SCHEMA_NAME = "intersubmod.exact_ps_read_af_factorization_oracle"
SCHEMA_VERSION = "1.0.0"
REPO_ROOT = Path(__file__).resolve().parents[3]
LEGACY_SCORER = (
    REPO_ROOT
    / "research/20260715_layered_workstation_genome_topology_multiselect"
    / "scripts/build_current_v5_read_af_topology.py"
)
TREE_SOLVER = (
    REPO_ROOT
    / "docs/methodology/_assets/20260627_subclone_4axis_teaching"
    / "scripts/tree_enumeration_solver.py"
)

GOLDEN_REGIONS: dict[str, dict[str, Any]] = {
    "chr1|PS=103318|HP=2|U3b6ee6088350bfcdd65cf9e2:B0001": {
        "read_af_fractions": ["1/1", "12/13"],
        "explicit_candidate_count": 2,
        "factorized_vertex_set_count": 2,
        "top_score_fraction": "1/13",
        "n_top_exact": 1,
        "n_top_shapes_exact": 1,
        "selection_class": "unique_first",
        "top_tree_ids": ["TR-94e0af01c4ad"],
        "top_shape_ids": ["TS-a99c6ae6b3"],
        "candidate_comparison_count_varies": False,
    },
    "chr1|PS=1528513|HP=2|Uc66a96aa84e5a99f332fd3d3:B0001": {
        "read_af_fractions": ["1/1", "1/1"],
        "explicit_candidate_count": 2,
        "factorized_vertex_set_count": 2,
        "top_score_fraction": "0/1",
        "n_top_exact": 2,
        "n_top_shapes_exact": 1,
        "selection_class": "tied_first_same_topology",
        "top_tree_ids": ["TR-2828254f77f0", "TR-94e0af01c4ad"],
        "top_shape_ids": ["TS-a99c6ae6b3"],
        "candidate_comparison_count_varies": False,
    },
    "chr1|PS=653452|HP=1|U206207b58e8c85579c67b04e:B0001": {
        "read_af_fractions": ["1/1", "1/1", "1/1"],
        "explicit_candidate_count": 7,
        "factorized_vertex_set_count": 7,
        "top_score_fraction": "0/1",
        "n_top_exact": 7,
        "n_top_shapes_exact": 2,
        "selection_class": "tied_first_different_topology",
        "top_tree_ids": [
            "TR-025c4567eace",
            "TR-6675f21b589d",
            "TR-688a77edae98",
            "TR-7cd0f5c6f39b",
            "TR-8d31746cf1f0",
            "TR-a107fd1651f0",
            "TR-a2dd238c7a27",
        ],
        "top_shape_ids": ["TS-23ba77f4e8", "TS-30a85152fe"],
        "candidate_comparison_count_varies": True,
    },
}

HCC1395_PILOT_EXPECTED = {
    "complete_ambiguous_t_gt1_units": 4601,
    "af_usable_units": 4505,
    "zero_denominator_units": 96,
    "explicit_candidate_count_all_units": 27687,
    "explicit_candidate_count_scored_units": 26942,
    "factorized_vertex_set_count_all_units": 27687,
    "factorization_mismatch_units": 0,
    "unique_first_units": 2473,
    "tied_first_same_topology_units": 1804,
    "tied_first_different_topology_units": 228,
    "exact_top_candidate_count": 9773,
    "exact_top_shape_count": 4830,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mlhp", type=Path, required=True)
    parser.add_argument("--layered", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pretty", action="store_true")
    return parser.parse_args()


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "sha256": sha256_path(resolved),
    }


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON must be an object")
    return value


def load_legacy_scorer() -> ModuleType:
    spec = importlib.util.spec_from_file_location("frozen_current_v5_read_af", LEGACY_SCORER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import frozen scorer: {LEGACY_SCORER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fraction_string(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def json_digest(value: Any) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def eligible_unit(unit: dict[str, Any]) -> bool:
    return (
        unit.get("analysis_candidate_set_complete") is True
        and unit.get("verification_status") == "full_pass"
        and unit.get("L1_base_class") == "ambiguous"
        and int(unit.get("n_trees") or 0) > 1
    )


def normalized_edges(tree: dict[str, Any]) -> list[tuple[str, str]]:
    output: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for raw_edge in tree.get("edges", []):
        if not isinstance(raw_edge, (list, tuple)) or len(raw_edge) != 2:
            raise ValueError("candidate edge must contain parent and child")
        edge = (str(raw_edge[0]), str(raw_edge[1]))
        if edge in seen:
            raise ValueError(f"candidate contains duplicate edge: {edge}")
        seen.add(edge)
        output.append(edge)
    return output


def factorize_candidates(
    trees: list[dict[str, Any]],
    read_af: dict[int, Fraction] | None,
    scorer: ModuleType,
) -> dict[str, Any]:
    """Count all mappings and, when AF exists, optimize without tree expansion."""

    by_vertex_set: dict[tuple[str, ...], list[list[tuple[str, str]]]] = defaultdict(list)
    for tree in trees:
        edges = normalized_edges(tree)
        nodes = {"ROOT"}
        for parent, child in edges:
            nodes.update((parent, child))
        by_vertex_set[tuple(sorted(nodes))].append(edges)

    candidate_count = 0
    scored_sets: list[tuple[Fraction, int]] = []
    per_set_count_mismatches = 0
    for nodes_key, explicit_mappings in by_vertex_set.items():
        state_by_label = {label: scorer.unlabel(label) for label in nodes_key}
        label_by_state: dict[frozenset[int], str] = {}
        for label, state in state_by_label.items():
            if state in label_by_state and label_by_state[state] != label:
                raise ValueError(
                    f"duplicate labels for one mutation state: "
                    f"{label_by_state[state]} and {label}"
                )
            label_by_state[state] = label

        set_candidate_count = 1
        set_best_count = 1
        set_best_score = Fraction(0, 1)
        for child in nodes_key:
            if child == "ROOT":
                continue
            child_state = state_by_label[child]
            choices: list[tuple[str, Fraction | None]] = []
            for mutation in child_state:
                parent_state = child_state - {mutation}
                parent = label_by_state.get(parent_state)
                if parent is None:
                    continue
                local_score = None
                if read_af is not None:
                    local_score, _ = scorer.exact_ordering_score(
                        [(parent, child)], read_af
                    )
                choices.append((parent, local_score))
            if not choices:
                raise ValueError(f"closed-set violation: {child} has no Hamming-1 parent")
            set_candidate_count *= len(choices)
            if read_af is not None:
                local_best = max(score for _, score in choices if score is not None)
                set_best_score += local_best
                set_best_count *= sum(score == local_best for _, score in choices)

        candidate_count += set_candidate_count
        if set_candidate_count != len(explicit_mappings):
            per_set_count_mismatches += 1
        if read_af is not None:
            scored_sets.append((set_best_score, set_best_count))

    result: dict[str, Any] = {
        "vertex_set_count": len(by_vertex_set),
        "candidate_count": candidate_count,
        "per_vertex_set_candidate_count_mismatches": per_set_count_mismatches,
    }
    if read_af is not None:
        best_score = max(score for score, _ in scored_sets)
        result.update(
            {
                "top_score_fraction": fraction_string(best_score),
                "top_candidate_count": sum(
                    count for score, count in scored_sets if score == best_score
                ),
            }
        )
    return result


def analyze_unit(
    unit: dict[str, Any],
    group: dict[str, Any] | None,
    scorer: ModuleType,
) -> dict[str, Any]:
    region = str(unit.get("region"))
    family = str(unit.get("family"))
    expected_candidates = int(unit.get("n_trees") or 0)
    trees = list(unit.get("trees") or [])
    issues: list[str] = []
    if group is None:
        return {
            "region": region,
            "family": family,
            "af_status": "missing_group",
            "explicit_candidate_count": len(trees),
            "expected_candidate_count": expected_candidates,
            "issues": ["missing_group"],
            "factorization_match": False,
        }
    if len(trees) != expected_candidates:
        issues.append("stored_candidate_count_mismatch")

    positions = [int(position) for position in group.get("positions", [])]
    coverage = (
        (group.get("col_coverage_by_hp") or {}).get(family, {}) or {}
    )
    read_af, read_af_rows, missing_reason = scorer.exact_read_af(coverage, positions)
    factorized = factorize_candidates(trees, read_af, scorer)
    base: dict[str, Any] = {
        "region": region,
        "family": family,
        "positions": positions,
        "expected_candidate_count": expected_candidates,
        "explicit_candidate_count": len(trees),
        "factorized_vertex_set_count": factorized["vertex_set_count"],
        "factorized_candidate_count": factorized["candidate_count"],
        "per_vertex_set_candidate_count_mismatches": factorized[
            "per_vertex_set_candidate_count_mismatches"
        ],
    }
    if factorized["candidate_count"] != len(trees):
        issues.append("factorized_candidate_count_mismatch")
    if factorized["per_vertex_set_candidate_count_mismatches"] != 0:
        issues.append("per_vertex_set_candidate_count_mismatch")

    if read_af is None:
        base.update(
            {
                "af_status": (
                    "zero_denominator"
                    if str(missing_reason).startswith("zero denominator")
                    else "unavailable"
                ),
                "af_unavailable_reason": missing_reason,
                "read_af_site_count_before_failure": len(read_af_rows),
                "issues": issues,
                "factorization_match": not issues,
            }
        )
        return base

    explicit_candidates: list[dict[str, Any]] = []
    for tree in trees:
        edges = normalized_edges(tree)
        score, comparisons = scorer.exact_ordering_score(edges, read_af)
        signature = scorer.canonical_shape(edges)
        explicit_candidates.append(
            {
                "tree_id": scorer.tree_id(edges),
                "shape_id": scorer.shape_id(signature),
                "score": score,
                "score_fraction": fraction_string(score),
                "n_comparisons": len(comparisons),
            }
        )
    explicit_best = max(candidate["score"] for candidate in explicit_candidates)
    explicit_top = [
        candidate
        for candidate in explicit_candidates
        if candidate["score"] == explicit_best
    ]
    top_tree_ids = sorted(candidate["tree_id"] for candidate in explicit_top)
    top_shape_ids = sorted({candidate["shape_id"] for candidate in explicit_top})
    selection_class = (
        "unique_first"
        if len(explicit_top) == 1
        else "tied_first_same_topology"
        if len(top_shape_ids) == 1
        else "tied_first_different_topology"
    )

    explicit_score = fraction_string(explicit_best)
    if factorized.get("top_score_fraction") != explicit_score:
        issues.append("factorized_top_score_mismatch")
    if int(factorized.get("top_candidate_count") or 0) != len(explicit_top):
        issues.append("factorized_top_candidate_count_mismatch")

    read_af_fractions = [
        row["read_af"]["fraction"]
        for row in read_af_rows
    ]
    candidate_digest_rows = sorted(
        [
            candidate["tree_id"],
            candidate["shape_id"],
            candidate["score_fraction"],
            candidate["n_comparisons"],
        ]
        for candidate in explicit_candidates
    )
    base.update(
        {
            "af_status": "usable",
            "read_af_fractions": read_af_fractions,
            "read_af_rows_sha256": json_digest(read_af_rows),
            "top_score_fraction": explicit_score,
            "n_top_exact": len(explicit_top),
            "n_top_shapes_exact": len(top_shape_ids),
            "selection_class": selection_class,
            "top_tree_ids": top_tree_ids,
            "top_shape_ids": top_shape_ids,
            "top_tree_ids_sha256": json_digest(top_tree_ids),
            "candidate_comparison_count_varies": len(
                {candidate["n_comparisons"] for candidate in explicit_candidates}
            )
            > 1,
            "explicit_candidate_score_digest_sha256": json_digest(
                candidate_digest_rows
            ),
            "factorized_top_score_fraction": factorized["top_score_fraction"],
            "factorized_top_candidate_count": factorized["top_candidate_count"],
            "issues": issues,
            "factorization_match": not issues,
        }
    )
    return base


def validate_golden(
    sample: str,
    unit_rows: list[dict[str, Any]],
    summary: dict[str, Any],
) -> dict[str, Any]:
    if sample != "HCC1395":
        return {
            "status": "not_applicable",
            "reason": "golden regions are frozen HCC1395 pilot fixtures",
            "definitions": GOLDEN_REGIONS,
        }

    by_region = {row["region"]: row for row in unit_rows}
    fixture_rows = []
    for region, expected in GOLDEN_REGIONS.items():
        observed = by_region.get(region)
        mismatches = {}
        if observed is None:
            mismatches["region"] = {"expected": "present", "observed": "missing"}
        else:
            for field, expected_value in expected.items():
                observed_value = observed.get(field)
                if observed_value != expected_value:
                    mismatches[field] = {
                        "expected": expected_value,
                        "observed": observed_value,
                    }
        fixture_rows.append(
            {
                "region": region,
                "pass": not mismatches,
                "expected": expected,
                "mismatches": mismatches,
            }
        )

    aggregate_mismatches = {}
    for field, expected_value in HCC1395_PILOT_EXPECTED.items():
        observed_value = summary.get(field)
        if observed_value != expected_value:
            aggregate_mismatches[field] = {
                "expected": expected_value,
                "observed": observed_value,
            }
    all_pass = all(row["pass"] for row in fixture_rows) and not aggregate_mismatches
    return {
        "status": "pass" if all_pass else "fail",
        "all_pass": all_pass,
        "fixtures": fixture_rows,
        "aggregate_expected": HCC1395_PILOT_EXPECTED,
        "aggregate_mismatches": aggregate_mismatches,
    }


def main() -> int:
    args = parse_args()
    scorer = load_legacy_scorer()
    mlhp = load_json(args.mlhp)
    layered = load_json(args.layered)
    sample = str(mlhp.get("sample") or layered.get("sample") or "")
    if not sample or str(layered.get("sample") or "") != sample:
        raise ValueError("MLHP/layered sample mismatch or missing sample")

    groups: dict[str, dict[str, Any]] = {}
    for group in mlhp.get("groups", []):
        region = str(group.get("region_id") or "")
        if not region or region in groups:
            raise ValueError(f"missing or duplicate MLHP region_id: {region}")
        groups[region] = group

    selected_units = [unit for unit in layered.get("detail", []) if eligible_unit(unit)]
    unit_rows = [
        analyze_unit(unit, groups.get(str(unit.get("region"))), scorer)
        for unit in selected_units
    ]
    status_counts = Counter(row["af_status"] for row in unit_rows)
    selection_counts = Counter(
        row["selection_class"]
        for row in unit_rows
        if row.get("af_status") == "usable"
    )
    summary = {
        "complete_ambiguous_t_gt1_units": len(unit_rows),
        "af_usable_units": status_counts["usable"],
        "zero_denominator_units": status_counts["zero_denominator"],
        "other_af_unavailable_units": len(unit_rows) - status_counts["usable"]
        - status_counts["zero_denominator"],
        "explicit_candidate_count_all_units": sum(
            int(row["explicit_candidate_count"]) for row in unit_rows
        ),
        "explicit_candidate_count_scored_units": sum(
            int(row["explicit_candidate_count"])
            for row in unit_rows
            if row.get("af_status") == "usable"
        ),
        "factorized_vertex_set_count_all_units": sum(
            int(row["factorized_vertex_set_count"]) for row in unit_rows
        ),
        "factorized_candidate_count_all_units": sum(
            int(row["factorized_candidate_count"]) for row in unit_rows
        ),
        "factorization_mismatch_units": sum(
            row.get("factorization_match") is not True for row in unit_rows
        ),
        "unique_first_units": selection_counts["unique_first"],
        "tied_first_same_topology_units": selection_counts[
            "tied_first_same_topology"
        ],
        "tied_first_different_topology_units": selection_counts[
            "tied_first_different_topology"
        ],
        "exact_top_candidate_count": sum(
            int(row.get("n_top_exact") or 0) for row in unit_rows
        ),
        "exact_top_shape_count": sum(
            int(row.get("n_top_shapes_exact") or 0) for row in unit_rows
        ),
        "candidate_comparison_count_varies_units": sum(
            row.get("candidate_comparison_count_varies") is True
            for row in unit_rows
        ),
    }
    golden = validate_golden(sample, unit_rows, summary)
    all_pass = (
        summary["factorization_mismatch_units"] == 0
        and (
            golden["status"] == "not_applicable"
            or golden.get("all_pass") is True
        )
    )
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "sample": sample,
        "scope": "complete, full-pass, L1 ambiguous, T>1 exact-PS units",
        "method": {
            "read_af": (
                "exact Fraction ALT/(REF+ALT) from exact-PS "
                "supported-pattern col_coverage_by_hp; not caller VAF"
            ),
            "explicit": "frozen current-v5 exact_ordering_score over every stored tree",
            "factorized": (
                "per minimum vertex set: sum child-local best edge scores and "
                "multiply local argmax counts; sum counts across globally best sets"
            ),
            "tie": "exact Fraction equality",
            "claim_ceiling": (
                "engineering parity oracle; not likelihood, posterior, CCF, "
                "cellular clone confirmation, or true ancestry"
            ),
        },
        "runtime": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
        },
        "inputs": {
            "mlhp": identity(args.mlhp),
            "layered": identity(args.layered),
        },
        "sources": {
            "oracle": identity(Path(__file__)),
            "legacy_read_af_scorer": identity(LEGACY_SCORER),
            "tree_enumeration_solver": identity(TREE_SOLVER),
        },
        "summary": summary,
        "golden_validation": golden,
        "unit_records_sha256": json_digest(unit_rows),
        "units": unit_rows,
        "all_pass": all_pass,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            receipt,
            ensure_ascii=False,
            indent=2 if args.pretty else None,
            separators=None if args.pretty else (",", ":"),
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "sample": sample,
                "summary": summary,
                "golden_status": golden["status"],
                "all_pass": all_pass,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if all_pass else 2


if __name__ == "__main__":
    sys.exit(main())
