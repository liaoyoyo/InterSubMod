#!/usr/bin/env python3
"""Fail-closed HCC1395 parity gate for the exact-PS C++ topology/AF engine.

The gate has three disjoint populations:

1. old-complete units: exact structural parity against the frozen Python
   display-all result;
2. old read-AF oracle units: exact Fraction score and best-tree tie parity;
3. old-capped units: explicitly excluded from old structural parity and
   accepted only when C++ either completes the family or honestly abstains.

The output is a compact receipt.  No candidate tree payload is copied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from collections import Counter
from pathlib import Path
from typing import Any


SCHEMA_NAME = "intersubmod.hcc1395_cpp_topology_af_gate"
SCHEMA_VERSION = "1.0.0"
MLHP_SCHEMA = ("intersubmod.exact_ps_layered_tree_input", "1.0.0")
CPP_UNIT_SCHEMA = ("intersubmod.exact_ps_cpp_topology_af.unit", "1.0.0")
CPP_RECEIPT_SCHEMA = ("intersubmod.exact_ps_cpp_topology_af.receipt", "1.0.0")
AF_ORACLE_SCHEMA = ("intersubmod.exact_ps_read_af_factorization_oracle", "1.0.0")
EXPECTED_SAMPLE = "HCC1395"
EXPECTED_OLD_CAPPED = 420
DEFAULT_CPP_SOURCE = (
    Path(__file__).resolve().parents[1] / "cpp/exact_ps_topology_af.cpp"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mlhp", type=Path, required=True)
    parser.add_argument("--layered", type=Path, required=True)
    parser.add_argument("--cpp-jsonl", type=Path, required=True)
    parser.add_argument("--cpp-receipt", type=Path, required=True)
    parser.add_argument("--af-oracle", type=Path, required=True)
    parser.add_argument("--cpp-source", type=Path, default=DEFAULT_CPP_SOURCE)
    parser.add_argument("--cpp-binary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-mismatch-examples", type=int, default=25)
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


def identity_matches(observed: Any, expected: dict[str, Any]) -> bool:
    return (
        isinstance(observed, dict)
        and observed.get("size_bytes") == expected["size_bytes"]
        and observed.get("sha256") == expected["sha256"]
    )


def json_digest(value: Any) -> str:
    payload = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON must be an object")
    return value


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line.endswith("\n"):
                raise ValueError(f"{path}:{line_number}: JSONL line lacks newline terminator")
            line = raw_line.rstrip("\n")
            if not line:
                raise ValueError(f"{path}:{line_number}: blank JSONL line")
            value = json.loads(line)
            if not isinstance(value, dict):
                raise ValueError(f"{path}:{line_number}: JSONL row must be an object")
            rows.append(value)
    return rows


def schema_is(value: dict[str, Any], expected: tuple[str, str]) -> bool:
    return (
        value.get("schema_name") == expected[0]
        and value.get("schema_version") == expected[1]
    )


def active_indices(group: dict[str, Any]) -> list[int]:
    family = str(group["hp_family"])
    k = len(group["positions"])
    full = ((group.get("populations_by_hp") or {}).get(family, {}) or {})
    partial = ((group.get("subread_groups_by_hp") or {}).get(family, {}) or {})
    patterns = list(full) + list(partial)
    if any(len(pattern) != k for pattern in patterns):
        raise ValueError(f"{group['region_id']}: pattern length differs from positions")
    return [
        index
        for index in range(k)
        if any(pattern[index] == "A" for pattern in patterns)
    ]


def label_mask(label: str, active: list[int], original_k: int) -> int:
    if label == "ROOT":
        return 0
    genotype = label[2:] if label.startswith("H_") else label
    if len(genotype) != original_k or set(genotype) - {"R", "A"}:
        raise ValueError(f"invalid old tree label: {label}")
    mask = 0
    for compact_index, original_index in enumerate(active):
        if genotype[original_index] == "A":
            mask |= 1 << compact_index
    return mask


def old_complete_contract(
    unit: dict[str, Any], group: dict[str, Any]
) -> dict[str, Any]:
    active = active_indices(group)
    original_k = len(group["positions"])
    expected_tree_count = int(unit["n_trees"])
    trees = list(unit.get("trees") or [])
    if len(trees) != expected_tree_count:
        raise ValueError(
            f"{unit['region']}: display-all stores {len(trees)} trees, "
            f"expected {expected_tree_count}"
        )

    family: set[tuple[int, ...]] = set()
    for tree in trees:
        nodes = {"ROOT"}
        for edge in tree.get("edges", []):
            if not isinstance(edge, (list, tuple)) or len(edge) != 2:
                raise ValueError(f"{unit['region']}: malformed old tree edge")
            nodes.update((str(edge[0]), str(edge[1])))
        vertices = tuple(
            sorted(label_mask(label, active, original_k) for label in nodes)
        )
        if len(vertices) != len(set(vertices)):
            raise ValueError(f"{unit['region']}: projected old tree repeats a vertex")
        family.add(vertices)
    if not family:
        raise ValueError(f"{unit['region']}: complete unit has an empty family")

    family_sorted = sorted(family)
    canonical = f"active_bit_count={len(active)}\n" + "".join(
        ",".join(str(vertex) for vertex in vertices) + "\n"
        for vertices in family_sorted
    )
    return {
        "active_bit_count": len(active),
        "objective_h": int(unit["n_hidden"]),
        "minimum_vertex_set_count": len(family_sorted),
        "minimum_family_sha256": hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
        "total_tree_count": expected_tree_count,
    }


def old_complete(unit: dict[str, Any]) -> bool:
    return (
        unit.get("analysis_candidate_set_complete") is True
        and unit.get("capped") is False
        and unit.get("verification_status") == "full_pass"
        and unit.get("verification_complete") is True
        and unit.get("verify_pass") is True
        and int(unit.get("analysis_trees_generated") or -1)
        == int(unit.get("n_trees") or 0)
    )


def append_mismatch(
    counts: Counter[str],
    examples: list[dict[str, Any]],
    limit: int,
    *,
    stage: str,
    region: str,
    field: str,
    expected: Any,
    observed: Any,
) -> None:
    key = f"{stage}:{field}"
    counts[key] += 1
    if len(examples) < limit:
        examples.append(
            {
                "stage": stage,
                "region": region,
                "field": field,
                "expected": expected,
                "observed": observed,
            }
        )


def compare_field(
    counts: Counter[str],
    examples: list[dict[str, Any]],
    limit: int,
    stage: str,
    region: str,
    field: str,
    expected: Any,
    observed: Any,
) -> None:
    if expected != observed:
        append_mismatch(
            counts,
            examples,
            limit,
            stage=stage,
            region=region,
            field=field,
            expected=expected,
            observed=observed,
        )


def parse_nonnegative_decimal(value: Any) -> int | None:
    if not isinstance(value, str) or not value or not value.isdigit():
        return None
    return int(value)


def compare_old_complete(
    old_units: list[dict[str, Any]],
    groups: dict[str, dict[str, Any]],
    cpp_by_region: dict[str, dict[str, Any]],
    mismatch_counts: Counter[str],
    mismatch_examples: list[dict[str, Any]],
    limit: int,
) -> tuple[dict[str, Any], list[list[Any]]]:
    parity_rows = []
    compared = 0
    for unit in old_units:
        region = str(unit["region"])
        group = groups.get(region)
        cpp = cpp_by_region.get(region)
        if group is None or cpp is None:
            append_mismatch(
                mismatch_counts,
                mismatch_examples,
                limit,
                stage="old_complete",
                region=region,
                field="record_presence",
                expected="group and C++ row",
                observed={
                    "group": group is not None,
                    "cpp": cpp is not None,
                },
            )
            continue
        expected = old_complete_contract(unit, group)
        compared += 1
        compare_field(
            mismatch_counts,
            mismatch_examples,
            limit,
            "old_complete",
            region,
            "objective_status",
            "OBJECTIVE_CERTIFIED",
            cpp.get("objective_status"),
        )
        compare_field(
            mismatch_counts,
            mismatch_examples,
            limit,
            "old_complete",
            region,
            "family_status",
            "FAMILY_COMPLETE",
            cpp.get("family_status"),
        )
        for field in (
            "active_bit_count",
            "objective_h",
            "minimum_vertex_set_count",
            "minimum_family_sha256",
        ):
            compare_field(
                mismatch_counts,
                mismatch_examples,
                limit,
                "old_complete",
                region,
                field,
                expected[field],
                cpp.get(field),
            )
        cpp_tree_count = parse_nonnegative_decimal(cpp.get("total_tree_count"))
        compare_field(
            mismatch_counts,
            mismatch_examples,
            limit,
            "old_complete",
            region,
            "total_tree_count",
            expected["total_tree_count"],
            cpp_tree_count,
        )
        parity_rows.append(
            [
                region,
                expected["active_bit_count"],
                expected["objective_h"],
                expected["minimum_vertex_set_count"],
                expected["minimum_family_sha256"],
                expected["total_tree_count"],
            ]
        )
    stage_mismatches = sum(
        count
        for key, count in mismatch_counts.items()
        if key.startswith("old_complete:")
    )
    return {
        "eligible_units": len(old_units),
        "compared_units": compared,
        "mismatch_fields": stage_mismatches,
        "old_family_count_sum": sum(row[3] for row in parity_rows),
    }, parity_rows


def validate_capped_row(row: dict[str, Any]) -> tuple[str, list[str]]:
    issues = []
    objective_status = row.get("objective_status")
    family_status = row.get("family_status")
    if family_status == "FAMILY_COMPLETE":
        if objective_status != "OBJECTIVE_CERTIFIED":
            issues.append("complete_without_certified_objective")
        for field in (
            "objective_h",
            "minimum_vertex_set_count",
            "minimum_family_sha256",
            "total_tree_count",
        ):
            if row.get(field) is None:
                issues.append(f"complete_missing_{field}")
        if parse_nonnegative_decimal(row.get("total_tree_count")) is None:
            issues.append("complete_invalid_total_tree_count")
        return "cpp_complete", issues

    abstain = (
        family_status in {"FAMILY_INCOMPLETE_CAP", "ABSTAIN_RESOURCE_LIMIT"}
        or objective_status == "ABSTAIN_RESOURCE_LIMIT"
    )
    if abstain:
        if row.get("unit_status") != "family_incomplete":
            issues.append("abstain_unit_status")
        if row.get("read_af_status") != "not_evaluable_family_incomplete":
            issues.append("abstain_read_af_status")
        for field in (
            "minimum_vertex_set_count",
            "minimum_family_sha256",
            "total_tree_count",
            "best_score_fraction",
            "best_tree_tie_count",
            "best_tree_unique",
        ):
            if row.get(field) is not None:
                issues.append(f"abstain_leaks_{field}")
        return "cpp_abstain", issues
    return "invalid", ["neither_complete_nor_resource_abstain"]


def check_old_capped(
    capped_units: list[dict[str, Any]],
    groups: dict[str, dict[str, Any]],
    cpp_by_region: dict[str, dict[str, Any]],
    mismatch_counts: Counter[str],
    mismatch_examples: list[dict[str, Any]],
    limit: int,
) -> tuple[dict[str, Any], list[list[Any]]]:
    outcomes = Counter()
    rows_for_digest = []
    for unit in capped_units:
        region = str(unit["region"])
        group = groups.get(region)
        cpp = cpp_by_region.get(region)
        if group is None or cpp is None:
            append_mismatch(
                mismatch_counts,
                mismatch_examples,
                limit,
                stage="old_capped",
                region=region,
                field="record_presence",
                expected="group and C++ row",
                observed={
                    "group": group is not None,
                    "cpp": cpp is not None,
                },
            )
            continue
        expected_active = len(active_indices(group))
        if cpp.get("active_bit_count") != expected_active:
            append_mismatch(
                mismatch_counts,
                mismatch_examples,
                limit,
                stage="old_capped",
                region=region,
                field="active_bit_count",
                expected=expected_active,
                observed=cpp.get("active_bit_count"),
            )
        outcome, issues = validate_capped_row(cpp)
        outcomes[outcome] += 1
        for issue in issues:
            append_mismatch(
                mismatch_counts,
                mismatch_examples,
                limit,
                stage="old_capped",
                region=region,
                field=issue,
                expected="valid C++ complete or fail-closed abstain",
                observed={
                    "objective_status": cpp.get("objective_status"),
                    "family_status": cpp.get("family_status"),
                    "unit_status": cpp.get("unit_status"),
                },
            )
        rows_for_digest.append(
            [
                region,
                expected_active,
                outcome,
                cpp.get("objective_status"),
                cpp.get("family_status"),
                cpp.get("objective_h"),
            ]
        )
    if len(capped_units) != EXPECTED_OLD_CAPPED:
        append_mismatch(
            mismatch_counts,
            mismatch_examples,
            limit,
            stage="old_capped",
            region="__cohort__",
            field="frozen_capped_count",
            expected=EXPECTED_OLD_CAPPED,
            observed=len(capped_units),
        )
    stage_mismatches = sum(
        count
        for key, count in mismatch_counts.items()
        if key.startswith("old_capped:")
    )
    return {
        "old_capped_units": len(capped_units),
        "explicitly_excluded_from_old_parity": len(capped_units),
        "cpp_complete": outcomes["cpp_complete"],
        "cpp_abstain": outcomes["cpp_abstain"],
        "invalid": outcomes["invalid"],
        "mismatch_fields": stage_mismatches,
    }, rows_for_digest


def compare_read_af(
    oracle_units: list[dict[str, Any]],
    cpp_by_region: dict[str, dict[str, Any]],
    mismatch_counts: Counter[str],
    mismatch_examples: list[dict[str, Any]],
    limit: int,
) -> tuple[dict[str, Any], list[list[Any]]]:
    status = Counter()
    parity_rows = []
    for expected in oracle_units:
        region = str(expected["region"])
        cpp = cpp_by_region.get(region)
        if cpp is None:
            append_mismatch(
                mismatch_counts,
                mismatch_examples,
                limit,
                stage="read_af",
                region=region,
                field="record_presence",
                expected="C++ row",
                observed=None,
            )
            continue
        af_status = expected.get("af_status")
        status[str(af_status)] += 1
        if af_status == "usable":
            compare_field(
                mismatch_counts,
                mismatch_examples,
                limit,
                "read_af",
                region,
                "read_af_status",
                "ranked_complete",
                cpp.get("read_af_status"),
            )
            compare_field(
                mismatch_counts,
                mismatch_examples,
                limit,
                "read_af",
                region,
                "best_score_fraction",
                expected.get("top_score_fraction"),
                cpp.get("best_score_fraction"),
            )
            cpp_ties = parse_nonnegative_decimal(cpp.get("best_tree_tie_count"))
            compare_field(
                mismatch_counts,
                mismatch_examples,
                limit,
                "read_af",
                region,
                "best_tree_tie_count",
                int(expected["n_top_exact"]),
                cpp_ties,
            )
            compare_field(
                mismatch_counts,
                mismatch_examples,
                limit,
                "read_af",
                region,
                "best_tree_unique",
                int(expected["n_top_exact"]) == 1,
                cpp.get("best_tree_unique"),
            )
            cpp_fractions = [
                row.get("fraction")
                for row in cpp.get("af_coverage", [])
                if isinstance(row, dict)
            ]
            compare_field(
                mismatch_counts,
                mismatch_examples,
                limit,
                "read_af",
                region,
                "af_coverage_fractions",
                expected.get("read_af_fractions"),
                cpp_fractions,
            )
            parity_rows.append(
                [
                    region,
                    expected.get("top_score_fraction"),
                    int(expected["n_top_exact"]),
                    expected.get("read_af_fractions"),
                ]
            )
        elif af_status == "zero_denominator":
            compare_field(
                mismatch_counts,
                mismatch_examples,
                limit,
                "read_af",
                region,
                "read_af_status",
                "zero_denominator",
                cpp.get("read_af_status"),
            )
            for field in (
                "best_score_fraction",
                "best_tree_tie_count",
                "best_tree_unique",
            ):
                compare_field(
                    mismatch_counts,
                    mismatch_examples,
                    limit,
                    "read_af",
                    region,
                    field,
                    None,
                    cpp.get(field),
                )
            parity_rows.append([region, "zero_denominator", None, None])
        else:
            append_mismatch(
                mismatch_counts,
                mismatch_examples,
                limit,
                stage="read_af",
                region=region,
                field="oracle_af_status",
                expected="usable or zero_denominator",
                observed=af_status,
            )
    stage_mismatches = sum(
        count for key, count in mismatch_counts.items() if key.startswith("read_af:")
    )
    return {
        "oracle_units": len(oracle_units),
        "usable_units": status["usable"],
        "zero_denominator_units": status["zero_denominator"],
        "mismatch_fields": stage_mismatches,
    }, parity_rows


def recompute_census(
    rows: list[dict[str, Any]], field: str
) -> dict[str, int]:
    return dict(sorted(Counter(str(row.get(field)) for row in rows).items()))


def main() -> int:
    args = parse_args()
    if args.max_mismatch_examples < 0:
        raise ValueError("--max-mismatch-examples must be non-negative")

    mlhp_identity = identity(args.mlhp)
    layered_identity = identity(args.layered)
    cpp_jsonl_identity = identity(args.cpp_jsonl)
    cpp_receipt_identity = identity(args.cpp_receipt)
    af_oracle_identity = identity(args.af_oracle)
    cpp_source_identity = identity(args.cpp_source)
    cpp_binary_identity = identity(args.cpp_binary) if args.cpp_binary else None

    mlhp = load_object(args.mlhp)
    layered = load_object(args.layered)
    cpp_rows = load_jsonl(args.cpp_jsonl)
    cpp_receipt = load_object(args.cpp_receipt)
    af_oracle = load_object(args.af_oracle)

    contract_checks: dict[str, bool] = {
        "mlhp_schema": schema_is(mlhp, MLHP_SCHEMA),
        "cpp_receipt_schema": schema_is(cpp_receipt, CPP_RECEIPT_SCHEMA),
        "af_oracle_schema": schema_is(af_oracle, AF_ORACLE_SCHEMA),
        "sample_is_HCC1395": mlhp.get("sample") == EXPECTED_SAMPLE
        and layered.get("sample") == EXPECTED_SAMPLE
        and af_oracle.get("sample") == EXPECTED_SAMPLE,
        "cpp_receipt_all_pass": cpp_receipt.get("all_pass") is True,
        "af_oracle_all_pass": af_oracle.get("all_pass") is True,
        "cpp_receipt_input_identity": identity_matches(
            cpp_receipt.get("input"), mlhp_identity
        ),
        "cpp_receipt_output_identity": identity_matches(
            cpp_receipt.get("output"), cpp_jsonl_identity
        ),
        "af_oracle_mlhp_identity": identity_matches(
            (af_oracle.get("inputs") or {}).get("mlhp"), mlhp_identity
        ),
        "af_oracle_layered_identity": identity_matches(
            (af_oracle.get("inputs") or {}).get("layered"), layered_identity
        ),
    }

    groups_list = list(mlhp.get("groups") or [])
    groups: dict[str, dict[str, Any]] = {}
    duplicate_groups = 0
    for group in groups_list:
        region = str(group.get("region_id") or "")
        if not region or region in groups:
            duplicate_groups += 1
        else:
            groups[region] = group
    cpp_by_region: dict[str, dict[str, Any]] = {}
    duplicate_cpp = 0
    cpp_row_contract_ok = True
    cpp_order_ok = len(cpp_rows) == len(groups_list)
    for index, row in enumerate(cpp_rows):
        cpp_row_contract_ok = cpp_row_contract_ok and schema_is(row, CPP_UNIT_SCHEMA)
        region = str(row.get("region_id") or "")
        if not region or region in cpp_by_region:
            duplicate_cpp += 1
        else:
            cpp_by_region[region] = row
        if index >= len(groups_list):
            cpp_order_ok = False
        else:
            cpp_order_ok = cpp_order_ok and (
                row.get("group_index") == index
                and region == str(groups_list[index].get("region_id"))
            )
    contract_checks.update(
        {
            "unique_mlhp_regions": duplicate_groups == 0,
            "unique_cpp_regions": duplicate_cpp == 0,
            "cpp_unit_schema_all_rows": cpp_row_contract_ok,
            "cpp_row_count_matches_groups": len(cpp_rows) == len(groups_list),
            "cpp_regions_exactly_match_groups": set(cpp_by_region) == set(groups),
            "cpp_group_order_exact": cpp_order_ok,
            "cpp_receipt_total_units": (cpp_receipt.get("counts") or {}).get(
                "total_units"
            )
            == len(cpp_rows),
        }
    )
    for field in ("unit_status", "objective_status", "family_status", "read_af_status"):
        contract_checks[f"cpp_receipt_{field}_census"] = (
            ((cpp_receipt.get("status_census") or {}).get(field) or {})
            == recompute_census(cpp_rows, field)
        )

    mismatch_counts: Counter[str] = Counter()
    mismatch_examples: list[dict[str, Any]] = []
    old_details = list(layered.get("detail") or [])
    old_complete_units = [unit for unit in old_details if old_complete(unit)]
    old_capped_units = [unit for unit in old_details if unit.get("capped") is True]

    structural_summary, structural_rows = compare_old_complete(
        old_complete_units,
        groups,
        cpp_by_region,
        mismatch_counts,
        mismatch_examples,
        args.max_mismatch_examples,
    )
    capped_summary, capped_rows = check_old_capped(
        old_capped_units,
        groups,
        cpp_by_region,
        mismatch_counts,
        mismatch_examples,
        args.max_mismatch_examples,
    )

    oracle_units = list(af_oracle.get("units") or [])
    read_af_summary, read_af_rows = compare_read_af(
        oracle_units,
        cpp_by_region,
        mismatch_counts,
        mismatch_examples,
        args.max_mismatch_examples,
    )
    contract_checks["oracle_unit_count_matches_summary"] = len(oracle_units) == int(
        (af_oracle.get("summary") or {}).get("complete_ambiguous_t_gt1_units")
        or -1
    )
    contract_checks["oracle_units_factorization_pass"] = all(
        unit.get("factorization_match") is True for unit in oracle_units
    )

    all_pass = (
        all(contract_checks.values())
        and not mismatch_counts
        and capped_summary["cpp_complete"] + capped_summary["cpp_abstain"]
        == capped_summary["old_capped_units"]
    )
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "sample": EXPECTED_SAMPLE,
        "all_pass": all_pass,
        "scope": {
            "old_complete_parity": (
                "active bits, exact objective_h, complete minimum family count/digest, "
                "and total parent-tree count"
            ),
            "read_af_parity": (
                "frozen 4,601-unit oracle; exact Fraction top score and exact top-tree count"
            ),
            "old_capped_policy": (
                "excluded from old parity; C++ must return complete family or explicit "
                "resource abstention with ranking fields withheld"
            ),
        },
        "runtime": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
        },
        "inputs": {
            "mlhp": mlhp_identity,
            "layered": layered_identity,
            "cpp_jsonl": cpp_jsonl_identity,
            "cpp_receipt": cpp_receipt_identity,
            "af_oracle": af_oracle_identity,
        },
        "sources": {
            "gate": identity(Path(__file__)),
            "cpp_source": cpp_source_identity,
            "cpp_binary": cpp_binary_identity,
        },
        "contract_checks": contract_checks,
        "summary": {
            "mlhp_groups": len(groups_list),
            "cpp_rows": len(cpp_rows),
            "old_complete": structural_summary,
            "old_capped": capped_summary,
            "read_af": read_af_summary,
            "cpp_status_census": {
                field: recompute_census(cpp_rows, field)
                for field in (
                    "unit_status",
                    "objective_status",
                    "family_status",
                    "read_af_status",
                )
            },
        },
        "parity_digests": {
            "old_complete_expected_sha256": json_digest(structural_rows),
            "old_capped_outcomes_sha256": json_digest(capped_rows),
            "read_af_expected_sha256": json_digest(read_af_rows),
        },
        "mismatch_counts": dict(sorted(mismatch_counts.items())),
        "mismatch_examples": mismatch_examples,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            receipt,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "all_pass": all_pass,
                "contract_checks_pass": sum(contract_checks.values()),
                "contract_checks_total": len(contract_checks),
                "mismatch_fields": sum(mismatch_counts.values()),
                "old_complete": structural_summary,
                "old_capped": capped_summary,
                "read_af": read_af_summary,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if all_pass else 2


if __name__ == "__main__":
    sys.exit(main())
