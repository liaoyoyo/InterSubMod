#!/usr/bin/env python3
"""Validate conservation, claim ceiling, and independent funnel agreement."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(result_path: Path, region_path: Path, cross_hp_path: Path) -> dict[str, Any]:
    payload = load_json(result_path)
    cross_hp = load_json(cross_hp_path)
    checks: list[dict[str, Any]] = []

    def check(name: str, passed: bool, observed: Any) -> None:
        checks.append({"name": name, "pass": bool(passed), "observed": observed})

    summaries = {summary["dataset"]: summary for summary in payload["sample_summaries"]}
    expected_rows = sum(summary["counts"]["regions_total"] for summary in summaries.values())
    with region_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    check("region_row_conservation", len(rows) == expected_rows, {"rows": len(rows), "expected": expected_rows})
    check(
        "no_biological_clone_k_claim",
        all(row["biological_clone_k_identified"] == "False" for row in rows)
        and payload["model_contract"]["biological_clone_k_identified"] is False,
        {
            "true_rows": sum(row["biological_clone_k_identified"] == "True" for row in rows),
            "contract": payload["model_contract"]["biological_clone_k_identified"],
        },
    )
    for sample, summary in summaries.items():
        counts = summary["counts"]
        state_dist = summary["distributions"]["primary_lineage_k_state_min_retained"]
        mutation_dist = summary["distributions"]["primary_lineage_k_mutation_state_min_retained"]
        joint_dist = summary["distributions"]["dual_hp_joint_mutation_state_min_conditional"]
        check(
            f"{sample}_lineage_distribution_conservation",
            sum(state_dist.values()) == counts["primary_lineage_units"]
            and sum(mutation_dist.values()) == counts["primary_lineage_units"],
            {
                "state": sum(state_dist.values()),
                "mutation": sum(mutation_dist.values()),
                "expected": counts["primary_lineage_units"],
            },
        )
        check(
            f"{sample}_dual_distribution_conservation",
            sum(joint_dist.values()) == counts["regions_dual_primary_hp"],
            {"distribution": sum(joint_dist.values()), "expected": counts["regions_dual_primary_hp"]},
        )
        check(
            f"{sample}_read_count_conservation",
            counts["primary_retained_reads"] + counts["primary_minread_censored_reads"]
            == counts["primary_reported_reads"],
            {
                "retained": counts["primary_retained_reads"],
                "censored": counts["primary_minread_censored_reads"],
                "reported": counts["primary_reported_reads"],
            },
        )
        check(
            f"{sample}_input_joins",
            summary["quality"]["region_to_mlhp_join_pass"]
            and summary["quality"]["region_to_ledger_allele_join_pass"]
            and summary["quality"]["invalid_pattern_count"] == 0,
            summary["quality"],
        )

    independent = {
        summary["dataset"]: summary["counts"] for summary in cross_hp["sample_summaries"]
    }
    for sample, summary in summaries.items():
        current = summary["counts"]
        prior = independent[sample]
        observed = {
            "regions_total": [current["regions_total"], prior["regions_total"]],
            "regions_dual_primary_hp": [
                current["regions_dual_primary_hp"],
                prior["regions_two_primary_hp"],
            ],
            "regions_pair_complete": [current["regions_pair_complete"], prior["regions_pair_complete"]],
        }
        check(
            f"{sample}_independent_funnel_agreement",
            all(left == right for left, right in observed.values()),
            observed,
        )

    comparison = payload["cross_technical_comparison"]
    exact = comparison["exact_allele_region_overlap"]
    check(
        "exact_region_denominators",
        exact["left"] == summaries["HCC1395"]["counts"]["regions_total"]
        and exact["right"] == summaries["HCC1395_DORADO"]["counts"]["regions_total"],
        exact,
    )
    priority = payload["priority_candidate_records"]
    check(
        "priority_definition_enforced",
        all(
            record["passes_phase_tree_cn_screen"]
            or record["state_evidence_level"] == "E2_exact_retained_catalog"
            or (record["joint_mutation_state_min_conditional"] or 0) >= 3
            for record in priority
        ),
        {"priority_records": len(priority)},
    )
    all_pass = all(item["pass"] for item in checks)
    return {
        "schema_name": "intersubmod.unknown_k_output_validation_receipt",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if all_pass else "FAIL",
        "all_pass": all_pass,
        "input_files": {
            "result": str(result_path.resolve()),
            "result_sha256": sha256_file(result_path),
            "region_table": str(region_path.resolve()),
            "region_table_sha256": sha256_file(region_path),
            "independent_cross_hp_audit": str(cross_hp_path.resolve()),
            "independent_cross_hp_audit_sha256": sha256_file(cross_hp_path),
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", required=True, type=Path)
    parser.add_argument("--region-table", required=True, type=Path)
    parser.add_argument("--cross-hp-audit", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    receipt = validate(args.result, args.region_table, args.cross_hp_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(receipt, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"status": receipt["status"], "checks": len(receipt["checks"]), "output": str(args.output.resolve())}, ensure_ascii=False))
    if not receipt["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
