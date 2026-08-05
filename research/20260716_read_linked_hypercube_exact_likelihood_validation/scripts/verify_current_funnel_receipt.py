#!/usr/bin/env python3
"""Independently verify the current ClairS-to-tree funnel receipt.

This verifier deliberately does not import the receipt builder.  It recomputes
the dataset and cohort totals directly from the seven site-ledger summaries,
checks their conservation identities, reconciles them with the canonical
topology summary, and then compares every reported count, ratio, path, and hash
against the receipt under review.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo


EXPECTED_DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)

EXPECTED_BRANCHES = (
    "excluded_by_longphase_filter",
    "out_of_scope_non_autosomal",
    "max_snv_excluded",
    "positional_singleton",
    "retained",
)

LEDGER_CHECKS = (
    "all_mlhp_retained_groups_joined",
    "all_raw_records_written",
    "autosomal_snv_conservation",
    "longphase_input_equals_recalibrated_all",
    "raw_all_equals_longphase_input",
    "recalibrated_PASS_equals_tree_input",
    "record_keys_are_unique_at_all_four_layers",
)

RECEIPT_CHECKS = (
    "all_7_dataset_ledgers_present",
    "all_ledger_checks_pass",
    "all_record_keys_unique",
    "all_branch_conservation_checks_pass",
    "all_dataset_counts_match_canonical",
    "aggregate_counts_match_canonical",
)


def file_sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest for a local file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    """Read UTF-8 JSON without relying on producer-side helpers."""

    return json.loads(path.read_text(encoding="utf-8"))


def nonnegative_integer(value: Any) -> bool:
    """Reject booleans even though bool is an int subclass in Python."""

    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def exact_ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def ratio_matches(observed: Any, expected: float | None) -> bool:
    if expected is None:
        return observed is None
    return (
        isinstance(observed, (int, float))
        and not isinstance(observed, bool)
        and math.isfinite(observed)
        and math.isclose(float(observed), expected, rel_tol=1e-12, abs_tol=1e-15)
    )


def verify_funnel(canonical_path: Path, ledger_root: Path, receipt_path: Path) -> dict[str, Any]:
    """Recompute the funnel and return a self-contained verification receipt."""

    canonical = read_json(canonical_path)
    receipt = read_json(receipt_path)
    failures: list[str] = []
    check_results: dict[str, bool] = {}

    def check(name: str, condition: bool, failure: str) -> None:
        condition = bool(condition)
        check_results[name] = check_results.get(name, True) and condition
        if not condition:
            failures.append(failure)

    canonical_payload = canonical.get("canonical") if isinstance(canonical, dict) else None
    canonical_payload = canonical_payload if isinstance(canonical_payload, dict) else {}
    canonical_aggregate = canonical_payload.get("aggregate")
    canonical_aggregate = canonical_aggregate if isinstance(canonical_aggregate, dict) else {}
    canonical_samples = canonical_payload.get("samples")
    canonical_samples = canonical_samples if isinstance(canonical_samples, list) else []

    canonical_names = [row.get("sample") for row in canonical_samples if isinstance(row, dict)]
    check("canonical_all_pass", canonical.get("all_pass") is True, "canonical all_pass 不是 true")
    check(
        "canonical_dataset_order",
        canonical_names == list(EXPECTED_DATASETS),
        f"canonical dataset 順序/全集不符：{canonical_names}",
    )
    check(
        "canonical_dataset_unique",
        len(canonical_names) == len(set(canonical_names)) == len(EXPECTED_DATASETS),
        "canonical sample 欄有缺失或重複",
    )
    canonical_by_dataset = {
        row.get("sample"): row
        for row in canonical_samples
        if isinstance(row, dict) and row.get("sample") in EXPECTED_DATASETS
    }

    receipt_inputs = receipt.get("inputs") if isinstance(receipt, dict) else None
    receipt_inputs = receipt_inputs if isinstance(receipt_inputs, dict) else {}
    receipt_canonical = receipt_inputs.get("canonical_json")
    receipt_canonical = receipt_canonical if isinstance(receipt_canonical, dict) else {}
    check(
        "receipt_canonical_path_binding",
        receipt_canonical.get("path") == str(canonical_path.resolve()),
        "receipt canonical path 未綁定本次實際輸入",
    )
    check(
        "receipt_canonical_hash_binding",
        receipt_canonical.get("sha256") == file_sha256(canonical_path),
        "receipt canonical SHA-256 與實際檔案不符",
    )
    check(
        "receipt_ledger_root_binding",
        receipt_inputs.get("site_ledger_root") == str(ledger_root.resolve()),
        "receipt site-ledger root 未綁定本次實際輸入",
    )

    receipt_scope = receipt.get("scope") if isinstance(receipt, dict) else None
    receipt_scope = receipt_scope if isinstance(receipt_scope, dict) else {}
    check(
        "receipt_schema_name",
        receipt.get("schema_name") == "intersubmod_current_ssnv_funnel_receipt",
        "receipt schema_name 不符",
    )
    check("receipt_schema_version", receipt.get("schema_version") == "1.0.0", "receipt schema_version 不符")
    check("receipt_task_type", receipt.get("task_type") == "B_comprehensive_validation", "receipt task_type 不是 Task B")
    check("receipt_all_pass", receipt.get("all_pass") is True, "receipt all_pass 不是 true")
    check(
        "receipt_scope_datasets",
        receipt_scope.get("datasets") == list(EXPECTED_DATASETS),
        "receipt scope datasets 不符",
    )
    check("receipt_scope_dataset_count", receipt_scope.get("dataset_count") == 7, "receipt dataset_count 不是 7")
    check(
        "receipt_scope_biological_count",
        receipt_scope.get("biological_sample_count") == 6,
        "receipt biological_sample_count 不是 6",
    )
    check(
        "receipt_scope_chromosomes",
        receipt_scope.get("chromosomes") == "chr1-chr22 for autosomal_biallelic_sSNV and downstream",
        "receipt chr1-chr22 scope 描述不符",
    )

    reported_sources = receipt_inputs.get("site_ledger_summaries")
    reported_sources = reported_sources if isinstance(reported_sources, list) else []
    source_names = [row.get("dataset") for row in reported_sources if isinstance(row, dict)]
    check(
        "receipt_source_dataset_order",
        source_names == list(EXPECTED_DATASETS),
        f"receipt source 清單順序/全集不符：{source_names}",
    )
    reported_source_by_dataset = {
        row.get("dataset"): row
        for row in reported_sources
        if isinstance(row, dict) and row.get("dataset") in EXPECTED_DATASETS
    }

    receipt_datasets = receipt.get("datasets")
    receipt_datasets = receipt_datasets if isinstance(receipt_datasets, list) else []
    receipt_names = [row.get("dataset") for row in receipt_datasets if isinstance(row, dict)]
    check(
        "receipt_dataset_order",
        receipt_names == list(EXPECTED_DATASETS),
        f"receipt dataset rows 順序/全集不符：{receipt_names}",
    )
    check(
        "receipt_dataset_unique",
        len(receipt_names) == len(set(receipt_names)) == len(EXPECTED_DATASETS),
        "receipt dataset rows 有缺失或重複",
    )
    receipt_by_dataset = {
        row.get("dataset"): row
        for row in receipt_datasets
        if isinstance(row, dict) and row.get("dataset") in EXPECTED_DATASETS
    }

    recomputed_rows: list[dict[str, Any]] = []
    for dataset in EXPECTED_DATASETS:
        ledger_path = ledger_root / "samples" / dataset / f"ssnv_site_ledger_{dataset}.summary.json"
        exists = ledger_path.is_file()
        check(f"{dataset}.ledger_present", exists, f"{dataset}: 缺少 site ledger {ledger_path}")
        if not exists:
            continue

        ledger = read_json(ledger_path)
        check(f"{dataset}.schema", ledger.get("schema_version") == "2.0", f"{dataset}: ledger schema_version 不是 2.0")
        check(f"{dataset}.sample", ledger.get("sample") == dataset, f"{dataset}: ledger sample 欄不符")
        check(f"{dataset}.pass", ledger.get("pass") is True, f"{dataset}: ledger pass 不是 true")
        check(
            f"{dataset}.input_contract",
            ledger.get("longphase_input_contract") == "clairs_raw_all",
            f"{dataset}: longphase input contract 不是 clairs_raw_all",
        )
        check(
            f"{dataset}.tree_contract",
            ledger.get("tree_contract") == "longphase_recalibrated_PASS",
            f"{dataset}: tree contract 不是 longphase_recalibrated_PASS",
        )

        ledger_checks = ledger.get("checks")
        ledger_checks = ledger_checks if isinstance(ledger_checks, dict) else {}
        for check_name in LEDGER_CHECKS:
            check(
                f"{dataset}.source_check.{check_name}",
                ledger_checks.get(check_name) is True,
                f"{dataset}: source check {check_name} 不是 true",
            )

        duplicate_counts = ledger.get("duplicate_record_key_excess")
        duplicate_counts = duplicate_counts if isinstance(duplicate_counts, dict) else {}
        expected_duplicate_layers = {"raw_clairs", "longphase_input", "longphase_recalibrated", "tree_input"}
        check(
            f"{dataset}.duplicate_layers",
            set(duplicate_counts) == expected_duplicate_layers,
            f"{dataset}: duplicate_record_key_excess layer 集合不符",
        )
        check(
            f"{dataset}.record_keys_unique",
            set(duplicate_counts) == expected_duplicate_layers
            and all(value == 0 for value in duplicate_counts.values()),
            f"{dataset}: record key duplicate excess 不全為 0",
        )

        branches = ledger.get("branch_counts")
        branches = branches if isinstance(branches, dict) else {}
        check(
            f"{dataset}.branch_schema",
            set(branches) == set(EXPECTED_BRANCHES),
            f"{dataset}: branch_counts 欄位集合不符：{sorted(branches)}",
        )

        named_counts = {
            "raw_clairs_records": ledger.get("raw_clairs_records"),
            "longphase_input_records": ledger.get("longphase_input_records"),
            "longphase_recalibrated_records": ledger.get("longphase_recalibrated_records"),
            "longphase_s_recalibrated_pass": ledger.get("tree_input_records"),
            "autosomal_biallelic_sSNV": ledger.get("autosomal_biallelic_snvs"),
            "retained_sSNV": branches.get("retained"),
        }
        all_numeric = all(nonnegative_integer(value) for value in [*named_counts.values(), *branches.values()])
        check(f"{dataset}.nonnegative_integer_counts", all_numeric, f"{dataset}: 主要/分支計數含非非負整數")
        if not all_numeric:
            continue

        raw = named_counts["raw_clairs_records"]
        lps_input = named_counts["longphase_input_records"]
        lps_all = named_counts["longphase_recalibrated_records"]
        tree = named_counts["longphase_s_recalibrated_pass"]
        autosomal = named_counts["autosomal_biallelic_sSNV"]
        retained = named_counts["retained_sSNV"]
        assert isinstance(raw, int) and isinstance(lps_input, int) and isinstance(lps_all, int)
        assert isinstance(tree, int) and isinstance(autosomal, int) and isinstance(retained, int)

        check(
            f"{dataset}.raw_input_all_conservation",
            raw == lps_input == lps_all,
            f"{dataset}: raw/input/recalibrated-all 不守恆",
        )
        check(
            f"{dataset}.raw_branch_partition",
            raw == sum(branches[name] for name in EXPECTED_BRANCHES),
            f"{dataset}: 五個互斥 branch 不等於 raw",
        )
        check(
            f"{dataset}.filter_transition_conservation",
            tree == raw - branches["excluded_by_longphase_filter"],
            f"{dataset}: LongPhase-S PASS != raw - excluded_by_longphase_filter",
        )
        check(
            f"{dataset}.autosomal_scope_conservation",
            autosomal == tree - branches["out_of_scope_non_autosomal"],
            f"{dataset}: autosomal != LongPhase-S PASS - out_of_scope_non_autosomal",
        )
        check(
            f"{dataset}.autosomal_branch_partition",
            autosomal == branches["max_snv_excluded"] + branches["positional_singleton"] + retained,
            f"{dataset}: autosomal != MAX_SNV-excluded + singleton + retained",
        )

        canonical_row = canonical_by_dataset.get(dataset)
        canonical_row = canonical_row if isinstance(canonical_row, dict) else {}
        check(
            f"{dataset}.canonical_tree",
            canonical_row.get("tree_input_records") == tree,
            f"{dataset}: tree count 與 canonical 不符",
        )
        check(
            f"{dataset}.canonical_autosomal",
            canonical_row.get("autosomal_biallelic_sSNV") == autosomal,
            f"{dataset}: autosomal count 與 canonical 不符",
        )
        check(
            f"{dataset}.canonical_retained",
            canonical_row.get("retained_sSNV") == retained,
            f"{dataset}: retained count 與 canonical 不符",
        )

        ledger_hash = file_sha256(ledger_path)
        source_row = reported_source_by_dataset.get(dataset)
        source_row = source_row if isinstance(source_row, dict) else {}
        check(
            f"{dataset}.input_source_path",
            source_row.get("path") == str(ledger_path.resolve()),
            f"{dataset}: receipt input source path 不符",
        )
        check(
            f"{dataset}.input_source_hash",
            source_row.get("sha256") == ledger_hash,
            f"{dataset}: receipt input source SHA-256 不符",
        )

        expected_row = {
            "dataset": dataset,
            "raw_clairs_records": raw,
            "longphase_s_recalibrated_pass": tree,
            "autosomal_biallelic_sSNV": autosomal,
            "retained_sSNV": retained,
            "branch_counts": {name: branches[name] for name in EXPECTED_BRANCHES},
            "relative_ratios": {
                "longphase_pass_over_raw": exact_ratio(tree, raw),
                "autosomal_over_longphase_pass": exact_ratio(autosomal, tree),
                "retained_over_autosomal": exact_ratio(retained, autosomal),
            },
            "total_ratios_over_raw": {
                "autosomal": exact_ratio(autosomal, raw),
                "retained": exact_ratio(retained, raw),
            },
            "source_path": str(ledger_path.resolve()),
            "source_sha256": ledger_hash,
        }
        observed_row = receipt_by_dataset.get(dataset)
        observed_row = observed_row if isinstance(observed_row, dict) else {}
        for count_name in (
            "raw_clairs_records",
            "longphase_s_recalibrated_pass",
            "autosomal_biallelic_sSNV",
            "retained_sSNV",
        ):
            check(
                f"{dataset}.receipt_count.{count_name}",
                observed_row.get(count_name) == expected_row[count_name],
                f"{dataset}: receipt {count_name} 不符"
                f"（observed={observed_row.get(count_name)}, expected={expected_row[count_name]}）",
            )
        check(
            f"{dataset}.receipt_branches",
            observed_row.get("branch_counts") == expected_row["branch_counts"],
            f"{dataset}: receipt branch_counts 不符",
        )
        check(
            f"{dataset}.receipt_source_path",
            observed_row.get("source_path") == expected_row["source_path"],
            f"{dataset}: dataset row source_path 不符",
        )
        check(
            f"{dataset}.receipt_source_hash",
            observed_row.get("source_sha256") == ledger_hash,
            f"{dataset}: dataset row source SHA-256 不符",
        )
        for ratio_group in ("relative_ratios", "total_ratios_over_raw"):
            observed_ratios = observed_row.get(ratio_group)
            observed_ratios = observed_ratios if isinstance(observed_ratios, dict) else {}
            for ratio_name, expected_ratio_value in expected_row[ratio_group].items():
                check(
                    f"{dataset}.receipt_ratio.{ratio_group}.{ratio_name}",
                    ratio_matches(observed_ratios.get(ratio_name), expected_ratio_value),
                    f"{dataset}: receipt ratio {ratio_group}.{ratio_name} 不符",
                )

        recomputed_rows.append(expected_row)

    total_fields = (
        "raw_clairs_records",
        "longphase_s_recalibrated_pass",
        "autosomal_biallelic_sSNV",
        "retained_sSNV",
    )
    totals = {field: sum(row[field] for row in recomputed_rows) for field in total_fields}
    aggregate_branches = {
        branch: sum(row["branch_counts"][branch] for row in recomputed_rows) for branch in EXPECTED_BRANCHES
    }
    complete_source_set = len(recomputed_rows) == len(EXPECTED_DATASETS)
    check("all_7_ledgers_recomputed", complete_source_set, f"只重算 {len(recomputed_rows)}/7 datasets")

    if complete_source_set:
        raw = totals["raw_clairs_records"]
        tree = totals["longphase_s_recalibrated_pass"]
        autosomal = totals["autosomal_biallelic_sSNV"]
        retained = totals["retained_sSNV"]
        check("aggregate_raw_partition", raw == sum(aggregate_branches.values()), "aggregate 五個互斥 branches 不等於 raw")
        check(
            "aggregate_filter_conservation",
            tree == raw - aggregate_branches["excluded_by_longphase_filter"],
            "aggregate LongPhase-S PASS != raw - filter-excluded",
        )
        check(
            "aggregate_autosomal_scope_conservation",
            autosomal == tree - aggregate_branches["out_of_scope_non_autosomal"],
            "aggregate autosomal != LongPhase-S PASS - non-autosomal",
        )
        check(
            "aggregate_autosomal_partition",
            autosomal
            == aggregate_branches["max_snv_excluded"]
            + aggregate_branches["positional_singleton"]
            + aggregate_branches["retained"],
            "aggregate autosomal != MAX_SNV-excluded + singleton + retained",
        )
        check("headline_raw", raw == 638259, f"headline raw 應為 638259，實得 {raw}")
        check("headline_lps_pass", tree == 582820, f"headline LongPhase-S PASS 應為 582820，實得 {tree}")
        check("headline_autosomal", autosomal == 469849, f"headline autosomal 應為 469849，實得 {autosomal}")
        check("headline_singleton", aggregate_branches["positional_singleton"] == 50432, "headline singleton 應為 50432")
        check(
            "headline_max_snv",
            aggregate_branches["max_snv_excluded"] == 225268,
            "headline MAX_SNV-excluded 應為 225268",
        )
        check("headline_retained", retained == 194149, f"headline retained 應為 194149，實得 {retained}")

        check(
            "canonical_aggregate_dataset_count",
            canonical_aggregate.get("dataset_count") == 7,
            "canonical aggregate dataset_count 不是 7",
        )
        check(
            "canonical_aggregate_tree",
            canonical_aggregate.get("tree_input_records") == tree,
            "canonical aggregate tree count 不符",
        )
        check(
            "canonical_aggregate_autosomal",
            canonical_aggregate.get("autosomal_biallelic_sSNV") == autosomal,
            "canonical aggregate autosomal count 不符",
        )
        check(
            "canonical_aggregate_retained",
            canonical_aggregate.get("retained_sSNV") == retained,
            "canonical aggregate retained count 不符",
        )

        expected_aggregate_ratios = {
            "relative_ratios": {
                "longphase_pass_over_raw": exact_ratio(tree, raw),
                "autosomal_over_longphase_pass": exact_ratio(autosomal, tree),
                "retained_over_autosomal": exact_ratio(retained, autosomal),
            },
            "total_ratios_over_raw": {
                "autosomal": exact_ratio(autosomal, raw),
                "retained": exact_ratio(retained, raw),
            },
        }
        receipt_aggregate = receipt.get("aggregate")
        receipt_aggregate = receipt_aggregate if isinstance(receipt_aggregate, dict) else {}
        for field, expected_value in totals.items():
            check(
                f"receipt_aggregate_count.{field}",
                receipt_aggregate.get(field) == expected_value,
                f"receipt aggregate {field} 不符（observed={receipt_aggregate.get(field)}, expected={expected_value}）",
            )
        check(
            "receipt_aggregate_branches",
            receipt_aggregate.get("branch_counts") == aggregate_branches,
            "receipt aggregate branch_counts 不符",
        )
        for ratio_group, expected_ratios in expected_aggregate_ratios.items():
            observed_ratios = receipt_aggregate.get(ratio_group)
            observed_ratios = observed_ratios if isinstance(observed_ratios, dict) else {}
            for ratio_name, expected_ratio_value in expected_ratios.items():
                check(
                    f"receipt_aggregate_ratio.{ratio_group}.{ratio_name}",
                    ratio_matches(observed_ratios.get(ratio_name), expected_ratio_value),
                    f"receipt aggregate ratio {ratio_group}.{ratio_name} 不符",
                )

    reported_checks = receipt.get("checks")
    reported_checks = reported_checks if isinstance(reported_checks, dict) else {}
    check(
        "receipt_check_schema",
        set(reported_checks) == set(RECEIPT_CHECKS),
        f"receipt checks 欄位集合不符：{sorted(reported_checks)}",
    )
    for check_name in RECEIPT_CHECKS:
        check(
            f"receipt_reported_check.{check_name}",
            reported_checks.get(check_name) is True,
            f"receipt 自報 check {check_name} 不是 true",
        )

    all_pass = not failures
    return {
        "schema_name": "intersubmod_current_ssnv_funnel_independent_verification",
        "schema_version": "1.0.0",
        "generated_at": dt.datetime.now(ZoneInfo("Asia/Taipei")).isoformat(),
        "task_type": "B_comprehensive_validation",
        "independence_statement": "Direct JSON recomputation; does not import or call build_current_funnel_receipt.py.",
        "inputs": {
            "canonical_json": {"path": str(canonical_path.resolve()), "sha256": file_sha256(canonical_path)},
            "site_ledger_root": str(ledger_root.resolve()),
            "receipt_under_test": {"path": str(receipt_path.resolve()), "sha256": file_sha256(receipt_path)},
            "site_ledger_summaries": [
                {
                    "dataset": row["dataset"],
                    "path": row["source_path"],
                    "sha256": row["source_sha256"],
                }
                for row in recomputed_rows
            ],
        },
        "scope": {
            "datasets": list(EXPECTED_DATASETS),
            "dataset_count": len(recomputed_rows),
            "biological_sample_count": 6,
            "chromosomes": "chr1-chr22 for autosomal_biallelic_sSNV and downstream",
            "unit_note": "Seven technical datasets; HCC1395 and HCC1395_DORADO represent one biological sample.",
        },
        "recomputed": {
            "aggregate": {
                **totals,
                "branch_counts": aggregate_branches,
                "autosomal_partition_identity": (
                    f"{aggregate_branches.get('positional_singleton', 0)} + "
                    f"{aggregate_branches.get('max_snv_excluded', 0)} + "
                    f"{aggregate_branches.get('retained', 0)} = {totals.get('autosomal_biallelic_sSNV', 0)}"
                ),
                "relative_ratios": {
                    "longphase_pass_over_raw": exact_ratio(
                        totals.get("longphase_s_recalibrated_pass", 0), totals.get("raw_clairs_records", 0)
                    ),
                    "autosomal_over_longphase_pass": exact_ratio(
                        totals.get("autosomal_biallelic_sSNV", 0),
                        totals.get("longphase_s_recalibrated_pass", 0),
                    ),
                    "retained_over_autosomal": exact_ratio(
                        totals.get("retained_sSNV", 0), totals.get("autosomal_biallelic_sSNV", 0)
                    ),
                },
                "total_ratios_over_raw": {
                    "autosomal": exact_ratio(
                        totals.get("autosomal_biallelic_sSNV", 0), totals.get("raw_clairs_records", 0)
                    ),
                    "retained": exact_ratio(totals.get("retained_sSNV", 0), totals.get("raw_clairs_records", 0)),
                },
            },
            "datasets": recomputed_rows,
        },
        "checks": check_results,
        "n_checks": len(check_results),
        "n_pass": sum(check_results.values()),
        "n_fail": sum(not value for value in check_results.values()),
        "failures": failures,
        "all_pass": all_pass,
    }


def atomic_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-json", type=Path, required=True)
    parser.add_argument("--ledger-root", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = verify_funnel(args.canonical_json, args.ledger_root, args.receipt)
    except (FileNotFoundError, json.JSONDecodeError) as exc:
        result = {
            "schema_name": "intersubmod_current_ssnv_funnel_independent_verification",
            "schema_version": "1.0.0",
            "all_pass": False,
            "failures": [str(exc)],
        }
    atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "all_pass": result.get("all_pass") is True,
                "output": str(args.output.resolve()),
                "output_sha256": file_sha256(args.output),
                "n_checks": result.get("n_checks"),
                "n_fail": result.get("n_fail"),
                "aggregate": (result.get("recomputed") or {}).get("aggregate"),
                "failures": result.get("failures"),
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0 if result.get("all_pass") is True else 2


if __name__ == "__main__":
    raise SystemExit(main())
