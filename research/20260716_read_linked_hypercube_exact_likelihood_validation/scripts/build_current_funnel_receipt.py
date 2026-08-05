#!/usr/bin/env python3
"""Build a source-backed ClairS-to-tree funnel receipt for the canonical cohort.

The receipt is intentionally strict: all seven dataset ledgers must pass their
own conservation checks and must agree with the canonical layered summary.
No funnel count is accepted from a hand-copied constant.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo


DATASETS = (
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
)

REQUIRED_CHECKS = (
    "all_mlhp_retained_groups_joined",
    "all_raw_records_written",
    "autosomal_snv_conservation",
    "longphase_input_equals_recalibrated_all",
    "raw_all_equals_longphase_input",
    "recalibrated_PASS_equals_tree_input",
    "record_keys_are_unique_at_all_four_layers",
)


class FunnelReceiptError(RuntimeError):
    """Raised when source ledgers cannot support the canonical funnel."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def build_receipt(canonical_json: Path, ledger_root: Path) -> dict[str, Any]:
    errors: list[str] = []
    canonical = load_json(canonical_json)
    require(canonical.get("all_pass") is True, "canonical all_pass 不是 true", errors)
    canonical_payload = canonical.get("canonical") or {}
    aggregate = canonical_payload.get("aggregate") or {}
    sample_rows = canonical_payload.get("samples") or []
    canonical_by_dataset = {row.get("sample"): row for row in sample_rows}
    require(tuple(canonical_by_dataset) == DATASETS, "canonical datasets/order 不是預期 7 datasets", errors)

    dataset_rows: list[dict[str, Any]] = []
    source_files: list[dict[str, Any]] = []
    for dataset in DATASETS:
        path = ledger_root / "samples" / dataset / f"ssnv_site_ledger_{dataset}.summary.json"
        if not path.is_file():
            errors.append(f"缺少 site ledger：{path}")
            continue
        ledger = load_json(path)
        prefix = f"{dataset}: "
        require(ledger.get("schema_version") == "2.0", prefix + "schema_version 不是 2.0", errors)
        require(ledger.get("sample") == dataset, prefix + "sample 欄不符", errors)
        require(ledger.get("pass") is True, prefix + "pass 不是 true", errors)
        require(ledger.get("longphase_input_contract") == "clairs_raw_all", prefix + "非 raw-all input contract", errors)
        require(ledger.get("tree_contract") == "longphase_recalibrated_PASS", prefix + "tree contract 不符", errors)
        checks = ledger.get("checks") or {}
        for check in REQUIRED_CHECKS:
            require(checks.get(check) is True, prefix + f"check {check} 不是 true", errors)
        duplicates = ledger.get("duplicate_record_key_excess") or {}
        require(all(value == 0 for value in duplicates.values()), prefix + "record key 有重複", errors)

        raw = ledger.get("raw_clairs_records")
        lps_input = ledger.get("longphase_input_records")
        lps_all = ledger.get("longphase_recalibrated_records")
        tree = ledger.get("tree_input_records")
        autosomal = ledger.get("autosomal_biallelic_snvs")
        branches = ledger.get("branch_counts") or {}
        retained = branches.get("retained")
        numeric = (raw, lps_input, lps_all, tree, autosomal, retained)
        require(all(isinstance(value, int) and not isinstance(value, bool) and value >= 0 for value in numeric), prefix + "主要計數不是非負整數", errors)
        if all(isinstance(value, int) for value in numeric):
            require(raw == lps_input == lps_all, prefix + "raw/input/recalibrated-all 不守恆", errors)
            require(raw == sum(branches.values()), prefix + "branch_counts 不等於 raw", errors)
            require(
                autosomal == retained + branches.get("positional_singleton", -1) + branches.get("max_snv_excluded", -1),
                prefix + "autosomal retained/singleton/MAX_SNV 分支不守恆",
                errors,
            )
            require(
                autosomal == raw - branches.get("excluded_by_longphase_filter", -1) - branches.get("out_of_scope_non_autosomal", -1),
                prefix + "raw→autosomal 分支不守恆",
                errors,
            )

        canonical_row = canonical_by_dataset.get(dataset) or {}
        require(canonical_row.get("tree_input_records") == tree, prefix + "tree_input 與 canonical 不符", errors)
        require(canonical_row.get("autosomal_biallelic_sSNV") == autosomal, prefix + "autosomal 與 canonical 不符", errors)
        require(canonical_row.get("retained_sSNV") == retained, prefix + "retained 與 canonical 不符", errors)

        dataset_rows.append(
            {
                "dataset": dataset,
                "raw_clairs_records": raw,
                "longphase_s_recalibrated_pass": tree,
                "autosomal_biallelic_sSNV": autosomal,
                "retained_sSNV": retained,
                "relative_ratios": {
                    "longphase_pass_over_raw": ratio(tree, raw),
                    "autosomal_over_longphase_pass": ratio(autosomal, tree),
                    "retained_over_autosomal": ratio(retained, autosomal),
                },
                "total_ratios_over_raw": {
                    "autosomal": ratio(autosomal, raw),
                    "retained": ratio(retained, raw),
                },
                "branch_counts": branches,
                "source_path": str(path.resolve()),
                "source_sha256": sha256(path),
            }
        )
        source_files.append({"dataset": dataset, "path": str(path.resolve()), "sha256": sha256(path)})

    totals: dict[str, int] = {}
    for key in ("raw_clairs_records", "longphase_s_recalibrated_pass", "autosomal_biallelic_sSNV", "retained_sSNV"):
        totals[key] = sum(row[key] for row in dataset_rows)
    branch_keys = (
        "excluded_by_longphase_filter",
        "out_of_scope_non_autosomal",
        "max_snv_excluded",
        "positional_singleton",
        "retained",
    )
    aggregate_branches = {
        key: sum(int(row["branch_counts"].get(key, 0)) for row in dataset_rows)
        for key in branch_keys
    }
    require(len(dataset_rows) == len(DATASETS), "沒有完整讀到 7 datasets", errors)
    require(aggregate.get("dataset_count") == len(DATASETS), "canonical dataset_count 不是 7", errors)
    require(aggregate.get("tree_input_records") == totals.get("longphase_s_recalibrated_pass"), "aggregate tree_input 不符", errors)
    require(aggregate.get("autosomal_biallelic_sSNV") == totals.get("autosomal_biallelic_sSNV"), "aggregate autosomal 不符", errors)
    require(aggregate.get("retained_sSNV") == totals.get("retained_sSNV"), "aggregate retained 不符", errors)
    if len(dataset_rows) == len(DATASETS):
        require(sum(aggregate_branches.values()) == totals.get("raw_clairs_records"), "aggregate raw branch counts 不守恆", errors)
        require(
            aggregate_branches["max_snv_excluded"]
            + aggregate_branches["positional_singleton"]
            + aggregate_branches["retained"]
            == totals.get("autosomal_biallelic_sSNV"),
            "aggregate autosomal branch counts 不守恆",
            errors,
        )

    if errors:
        raise FunnelReceiptError("\n".join(errors))

    raw = totals["raw_clairs_records"]
    tree = totals["longphase_s_recalibrated_pass"]
    autosomal = totals["autosomal_biallelic_sSNV"]
    retained = totals["retained_sSNV"]
    return {
        "schema_name": "intersubmod_current_ssnv_funnel_receipt",
        "schema_version": "1.0.0",
        "generated_at": dt.datetime.now(ZoneInfo("Asia/Taipei")).isoformat(),
        "task_type": "B_comprehensive_validation",
        "scope": {
            "datasets": list(DATASETS),
            "dataset_count": 7,
            "biological_sample_count": 6,
            "chromosomes": "chr1-chr22 for autosomal_biallelic_sSNV and downstream",
            "technical_replication_note": "HCC1395 and HCC1395_DORADO are two technical datasets from one biological sample",
        },
        "inputs": {
            "canonical_json": {"path": str(canonical_json.resolve()), "sha256": sha256(canonical_json)},
            "site_ledger_root": str(ledger_root.resolve()),
            "site_ledger_summaries": source_files,
        },
        "aggregate": {
            **totals,
            "branch_counts": aggregate_branches,
            "relative_ratios": {
                "longphase_pass_over_raw": ratio(tree, raw),
                "autosomal_over_longphase_pass": ratio(autosomal, tree),
                "retained_over_autosomal": ratio(retained, autosomal),
            },
            "total_ratios_over_raw": {
                "autosomal": ratio(autosomal, raw),
                "retained": ratio(retained, raw),
            },
        },
        "datasets": dataset_rows,
        "checks": {
            "all_7_dataset_ledgers_present": True,
            "all_ledger_checks_pass": True,
            "all_record_keys_unique": True,
            "all_branch_conservation_checks_pass": True,
            "all_dataset_counts_match_canonical": True,
            "aggregate_counts_match_canonical": True,
        },
        "all_pass": True,
    }


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-json", type=Path, required=True)
    parser.add_argument("--ledger-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        receipt = build_receipt(args.canonical_json, args.ledger_root)
        atomic_write_json(args.output, receipt)
    except (FunnelReceiptError, FileNotFoundError, json.JSONDecodeError) as exc:
        print(json.dumps({"all_pass": False, "error": str(exc)}, ensure_ascii=False, indent=2))
        return 2
    print(
        json.dumps(
            {
                "all_pass": True,
                "output": str(args.output.resolve()),
                "output_sha256": sha256(args.output),
                "aggregate": receipt["aggregate"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
