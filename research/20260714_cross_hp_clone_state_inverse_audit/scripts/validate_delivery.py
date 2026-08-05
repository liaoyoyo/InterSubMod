#!/usr/bin/env python3
"""Validate numerical, narrative, and provenance integrity of the final delivery."""

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


def tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--ascn", required=True, type=Path)
    parser.add_argument("--candidate-tsv", required=True, type=Path)
    parser.add_argument("--sample-tsv", required=True, type=Path)
    parser.add_argument("--markdown", required=True, type=Path)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--browser-qa", required=True, type=Path)
    parser.add_argument("--bounded-receipt", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    audit = load_json(args.audit)
    comparison = load_json(args.comparison)
    ascn = load_json(args.ascn)
    browser = load_json(args.browser_qa)
    bounded = load_json(args.bounded_receipt)
    markdown = args.markdown.read_text(encoding="utf-8")
    standalone = args.html.read_text(encoding="utf-8")
    counts = audit["aggregate_counts"]
    summaries = audit["sample_summaries"]
    records = audit["candidate_records"]
    candidate_rows = tsv_rows(args.candidate_tsv)
    sample_rows = tsv_rows(args.sample_tsv)

    checks: dict[str, bool] = {
        "input_run_succeeded": audit["input_run"]["state"] == "SUCCEEDED",
        "dataset_count_7": audit["dataset_count"] == 7 == len(summaries) == len(sample_rows),
        "biological_sample_count_6": audit["biological_sample_count"] == 6,
        "regions_total_reconciles": sum(
            row["counts"]["regions_total"] for row in summaries
        ) == counts["regions_total"] == 51_815,
        "dual_hp_reconciles": sum(
            row["counts"]["regions_two_primary_hp"] for row in summaries
        ) == counts["regions_two_primary_hp"] == 22_779,
        "topology_invariant_recomputed": sum(
            1 for row in records if row["direct_sister_shape_invariant"]
        ) == counts["regions_direct_sister_shape_invariant"] == 54,
        "tree_unique_recomputed": sum(
            1 for row in records if row["direct_sister_tree_unique"]
        ) == counts["regions_direct_sister_tree_unique"] == 17,
        "complete_topology_recomputed": sum(
            1 for row in records if row["direct_sister_shape_invariant_analysis_complete"]
        ) == counts["regions_direct_sister_shape_invariant_pair_complete"] == 45,
        "complete_tree_unique_recomputed": sum(
            1 for row in records if row["direct_sister_tree_unique_analysis_complete"]
        ) == counts["regions_direct_sister_tree_unique_pair_complete"] == 13,
        "observed_strict_recomputed": sum(
            1 for row in records if row.get("observed_hp1_direct_only_hp2_sister_only")
        ) == counts["regions_observed_hp1_direct_only_hp2_sister_only"] == 7,
        "candidate_tsv_row_count": len(candidate_rows) == len(records),
        "pattern_ready_zero": counts["regions_pattern_level_inverse_ready"] == 0,
        "sensitivity_pattern_ready_zero": comparison["fixed_two_site_catalog_conclusion"][
            "sensitivity_pattern_ready"
        ] == 0,
        "ascn_only_one_topology_candidate": ascn["category_stable_pass_counts"][
            "topology_invariant"
        ] == 1,
        "ascn_zero_tree_unique": ascn["category_stable_pass_counts"]["tree_unique"] == 0,
        "h2009_ascn_pass": ascn["H2009_chr1_120007237_120040749"]["stable_pass"],
        "ascn_mother_set_reconciles": ascn["dual_hp_mother_set_screen"]["total"]
        == counts["regions_two_primary_hp"]
        == 22_779,
        "ascn_declared_gate_counts": (
            ascn["dual_hp_mother_set_screen"]["magnitude_pass"] == 1_951
            and ascn["dual_hp_mother_set_screen"]["stable_pass"] == 1_889
            and ascn["dual_hp_mother_set_screen"]["by_dataset"]["H2009"]["magnitude_pass"]
            == 1_779
        ),
        "bounded_receipt_all_pass": bounded["all_pass"],
        "bounded_receipt_binds_current_audit": bounded["input_paths"]["candidate_audit"][
            "sha256_at_validation"
        ]
        == sha256_file(args.audit),
        "bounded_raw_counts_reconcile": (
            bounded["bounded_raw_read_results"]["HCC1937_chr20_58489593_58518912"][
                "alignment_exposures"
            ]
            == 306
            and bounded["bounded_raw_read_results"]["H2009_chr1_120007237_120040749"][
                "alignment_exposures"
            ]
            == 183
        ),
        "bounded_methyl_counts_reconcile": (
            bounded["methyl_results"]["aggregate_54_topology_candidates"][
                "bounded_exact_sidecar_join_regions"
            ]
            == 54
            and bounded["methyl_results"]["aggregate_54_topology_candidates"][
                "both_HP1_HP2_MM_ML_rate_ge_0_95_regions"
            ]
            == 49
            and bounded["methyl_results"]["aggregate_54_topology_candidates"][
                "both_HP1_HP2_have_full_span_MM_ML_reads_regions"
            ]
            == 29
            and bounded["methyl_results"]["H2009_chr1_120007237_120040749"][
                "bounded_5mC_entries_total_HP1_HP2_full_span_reads"
            ]
            == 3_155
        ),
        "browser_qa_all_pass": browser["all_pass"],
        "markdown_metadata_present": markdown.startswith("<!--"),
        "markdown_framework_present": "用 SCQA + Claim–Evidence–Limit" in markdown,
        "markdown_key_counts_present": all(
            token in markdown
            for token in ("51,815", "22,779", "54", "45", "17", "13", "1,951", "1,889")
        ),
        "markdown_completeness_ceiling_present": all(
            token in markdown
            for token in ("analysis-complete", "candidate set 完整", "目前各側只儲存一棵樹")
        ),
        "markdown_old_ascn_counts_absent": all(
            token not in markdown for token in ("1,955", "1,780")
        ),
        "markdown_vaf_ccf_ceiling_present": "不是 cancer cell fraction（CCF）" in markdown,
        "markdown_bounded_receipt_linked": "bounded_read_methyl_validation_receipt.json" in markdown,
        "markdown_claim_ceiling_present": "不支持已確認四個 cellular clones" in markdown,
        "html_title_present": "<title>跨 HP clone-state 反推｜可行性與觀察紀錄</title>" in standalone,
        "html_key_counts_present": all(
            token in standalone
            for token in ("51,815", "22,779", ">54<", ">45<", ">17<", ">13<", "1,951", "1,889")
        ),
        "html_completeness_ceiling_present": all(
            token in standalone for token in ("analysis-complete", "目前各側只儲存一棵樹")
        ),
        "html_old_ascn_counts_absent": all(
            token not in standalone for token in ("1,955", "1,780")
        ),
        "html_vaf_ccf_ceiling_present": "不是 cancer cell fraction" in standalone,
        "html_bounded_receipt_linked": "bounded_read_methyl_validation_receipt.json" in standalone,
        "html_no_external_cdn": all(
            token not in standalone for token in ("https://cdn", "unpkg.com", "fonts.googleapis.com")
        ),
        "html_claim_ceiling_present": "not confirmed cell-level HP pairing" in standalone,
    }
    failed = [name for name, passed in checks.items() if not passed]
    artifacts = [
        args.audit,
        args.comparison,
        args.ascn,
        args.candidate_tsv,
        args.sample_tsv,
        args.markdown,
        args.html,
        args.browser_qa,
        args.bounded_receipt,
    ]
    receipt = {
        "schema_name": "intersubmod.cross_hp_clone_state_delivery_validation",
        "schema_version": "1.1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "overall_assessment": "ready_to_share_as_observation_report" if not failed else "needs_revision",
        "all_pass": not failed,
        "checks": checks,
        "failed_checks": failed,
        "calculation_spot_checks": {
            "regions_total": counts["regions_total"],
            "two_primary_hp": counts["regions_two_primary_hp"],
            "topology_invariant": counts["regions_direct_sister_shape_invariant"],
            "tree_unique": counts["regions_direct_sister_tree_unique"],
            "complete_topology_invariant": counts[
                "regions_direct_sister_shape_invariant_pair_complete"
            ],
            "complete_tree_unique": counts["regions_direct_sister_tree_unique_pair_complete"],
            "observed_strict": counts["regions_observed_hp1_direct_only_hp2_sister_only"],
            "pattern_level_inverse_ready": counts["regions_pattern_level_inverse_ready"],
            "candidate_record_rows": len(records),
            "ascn_stable_topology_overlap": ascn["category_stable_pass_counts"][
                "topology_invariant"
            ],
            "ascn_magnitude_mother_set": ascn["dual_hp_mother_set_screen"]["magnitude_pass"],
            "ascn_stable_mother_set": ascn["dual_hp_mother_set_screen"]["stable_pass"],
            "bounded_receipt_all_pass": bounded["all_pass"],
        },
        "required_caveats": [
            "7 datasets correspond to 6 biological samples",
            "topology and same-site ALT are parallel evidence sets, not a linear funnel",
            "SAVANA 1+1 is a magnitude screen without HP1/HP2 orientation",
            "ranked output row count is not proof of global purity/ploidy solution uniqueness",
            "HP1/HP2 bulk marginals do not provide cell-level homolog pairing",
            "methyl exact join does not confirm a clone or tree edge",
            "VAF/read fraction is not CCF without purity, allele-specific CN, multiplicity, and normal correction",
            "zero inference-ready regions means zero under the fixed thresholded catalog, not biological absence",
        ],
        "artifact_sha256": {
            str(path.resolve()): sha256_file(path) for path in artifacts
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
