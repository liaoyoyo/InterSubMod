import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys

import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from build_report_artifact import (  # noqa: E402
    COMPONENT_REQUIRED_FIELDS,
    EXACT_LOG_REQUIRED_FIELDS,
    HP_PS_PAIR_REQUIRED_FIELDS,
    HP_PS_SOURCE_CHECKS,
    HP_PS_UNIT_REQUIRED_FIELDS,
    ReportInputError,
    SUMMARY_REQUIRED_FIELDS,
    main,
    validate_component_retention,
    validate_summary_row,
)


GENERATED_AT = "2026-07-18T00:00:00Z"


def test_component_retention_blank_is_valid_only_for_zero_denominator():
    zero_denominator = {
        "raw_total_molecule_weight": "0",
        "raw_retained_molecule_weight": "0",
        "old_densest8_retained_molecule_weight": "0",
        "raw_retention_ratio": "",
        "old_densest8_retention_ratio": "",
    }
    assert validate_component_retention(zero_denominator, "zero") == (
        None,
        None,
        None,
    )

    fabricated_zero = dict(zero_denominator, raw_retention_ratio="0")
    with pytest.raises(ReportInputError, match="must be blank exactly"):
        validate_component_retention(fabricated_zero, "fabricated")

    missing_nonzero = dict(
        zero_denominator,
        raw_total_molecule_weight="10",
        raw_retained_molecule_weight="5",
        old_densest8_retained_molecule_weight="4",
    )
    with pytest.raises(ReportInputError, match="must be blank exactly"):
        validate_component_retention(missing_nonzero, "missing")


def _summary_rows():
    chr1 = {
        "chrom": "chr1",
        "partition_stage_status": "FIXTURE_COMPLETED",
        "ssnv_sites": 1000,
        "components": 1,
        "sites": 10,
        "old_selected_sites": 8,
        "old_excluded_sites": 2,
        "new_blocks": 2,
        "new_retained_sites": 10,
        "primary_active_sites": 9,
        "raw_total_molecule_weight": 100,
        "new_retained_molecule_weight": 90,
        "new_lost_molecule_weight": 10,
        "unavoidable_patterns": 2,
        "unavoidable_molecule_weight": 3,
        "unavoidable_n_fixed_ra_gt8_patterns": 1,
        "unavoidable_n_fixed_ra_gt8_molecule_weight": 1,
        "unavoidable_n_fixed_ra_lte8_span_gt8_patterns": 1,
        "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight": 2,
        "new_retention_ratio": 0.9,
        "old_densest8_retained_molecule_weight": 70,
        "old_densest8_retention_ratio": 0.7,
        "weight_stable_components": 1,
        "weight_sensitive_components": 0,
        "zero_evidence_blocks": 0,
        "zero_evidence_block_sites": 0,
        "tree_ready_blocks": 2,
        "tree_ready_block_sites": 10,
        "abstain_blocks": 0,
        "abstain_block_sites": 0,
        "extraction_wall_seconds": 60,
        "partition_wall_seconds": 2,
        "partition_pattern_load_aggregate_seconds": 0.5,
        "partition_ordered_hypergraph_dp_seconds": 0.2,
    }
    chr22 = {
        "chrom": "chr22",
        "partition_stage_status": "FIXTURE_COMPLETED",
        "ssnv_sites": 543,
        "components": 1,
        "sites": 12,
        "old_selected_sites": 8,
        "old_excluded_sites": 4,
        "new_blocks": 2,
        "new_retained_sites": 12,
        "primary_active_sites": 10,
        "raw_total_molecule_weight": 200,
        "new_retained_molecule_weight": 170,
        "new_lost_molecule_weight": 30,
        "unavoidable_patterns": 3,
        "unavoidable_molecule_weight": 5,
        "unavoidable_n_fixed_ra_gt8_patterns": 2,
        "unavoidable_n_fixed_ra_gt8_molecule_weight": 4,
        "unavoidable_n_fixed_ra_lte8_span_gt8_patterns": 1,
        "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight": 1,
        "new_retention_ratio": 0.85,
        "old_densest8_retained_molecule_weight": 120,
        "old_densest8_retention_ratio": 0.6,
        "weight_stable_components": 0,
        "weight_sensitive_components": 1,
        "zero_evidence_blocks": 1,
        "zero_evidence_block_sites": 4,
        "tree_ready_blocks": 1,
        "tree_ready_block_sites": 8,
        "abstain_blocks": 1,
        "abstain_block_sites": 4,
        "extraction_wall_seconds": 120,
        "partition_wall_seconds": 4,
        "partition_pattern_load_aggregate_seconds": 1.0,
        "partition_ordered_hypergraph_dp_seconds": 0.5,
    }
    totals = {
        "chrom": "ALL",
        "partition_stage_status": "FIXTURE_2_CHROMS",
        "ssnv_sites": 1543,
        "components": 2,
        "sites": 22,
        "old_selected_sites": 16,
        "old_excluded_sites": 6,
        "new_blocks": 4,
        "new_retained_sites": 22,
        "primary_active_sites": 19,
        "raw_total_molecule_weight": 300,
        "new_retained_molecule_weight": 260,
        "new_lost_molecule_weight": 40,
        "unavoidable_patterns": 5,
        "unavoidable_molecule_weight": 8,
        "unavoidable_n_fixed_ra_gt8_patterns": 3,
        "unavoidable_n_fixed_ra_gt8_molecule_weight": 5,
        "unavoidable_n_fixed_ra_lte8_span_gt8_patterns": 2,
        "unavoidable_n_fixed_ra_lte8_span_gt8_molecule_weight": 3,
        "new_retention_ratio": 260 / 300,
        "old_densest8_retained_molecule_weight": 190,
        "old_densest8_retention_ratio": 190 / 300,
        "weight_stable_components": 1,
        "weight_sensitive_components": 1,
        "zero_evidence_blocks": 1,
        "zero_evidence_block_sites": 4,
        "tree_ready_blocks": 3,
        "tree_ready_block_sites": 18,
        "abstain_blocks": 1,
        "abstain_block_sites": 4,
        "extraction_wall_seconds": 180,
        "partition_wall_seconds": 6,
        "partition_pattern_load_aggregate_seconds": 1.5,
        "partition_ordered_hypergraph_dp_seconds": 0.7,
    }
    return [chr1, chr22], totals


def _component_rows():
    base = {
        field: "" for field in COMPONENT_REQUIRED_FIELDS
    }
    chr1 = {
        **base,
        "dataset": "HCC1395",
        "chrom": "chr1",
        "legacy_component_id": "fixture_chr1_component_1",
        "start1": "100",
        "end1": "1000",
        "span_bp": "900",
        "pre_cap_k": "10",
        "old_densest8_selected": "8",
        "old_cap_excluded": "2",
        "new_block_count": "2",
        "new_site_retained": "10",
        "primary_active_site_count": "9",
        "primary_active_site_fraction": "0.9",
        "raw_total_molecule_weight": "100",
        "raw_retained_molecule_weight": "90",
        "raw_lost_molecule_weight": "10",
        "raw_retention_ratio": "0.9",
        "old_densest8_retained_molecule_weight": "70",
        "old_densest8_retention_ratio": "0.7",
        "weight_stable": "true",
        "status": "fixture",
        "positions_sha256": "0" * 64,
    }
    chr22 = {
        **base,
        "dataset": "HCC1395",
        "chrom": "chr22",
        "legacy_component_id": "fixture_chr22_component_1",
        "start1": "200",
        "end1": "2200",
        "span_bp": "2000",
        "pre_cap_k": "12",
        "old_densest8_selected": "8",
        "old_cap_excluded": "4",
        "new_block_count": "2",
        "new_site_retained": "12",
        "primary_active_site_count": "10",
        "primary_active_site_fraction": str(10 / 12),
        "raw_total_molecule_weight": "200",
        "raw_retained_molecule_weight": "170",
        "raw_lost_molecule_weight": "30",
        "raw_retention_ratio": "0.85",
        "old_densest8_retained_molecule_weight": "120",
        "old_densest8_retention_ratio": "0.6",
        "weight_stable": "false",
        "status": "fixture",
        "positions_sha256": "1" * 64,
    }
    return [chr1, chr22]


def _write_fixture(tmp_path: Path):
    per_chrom, totals = _summary_rows()
    summary = {
        "schema_name": "intersubmod.hcc1395_full_k_gt8_segmentation.summary",
        "schema_version": "1.1.0",
        "sample": "HCC1395",
        "all_pass": True,
        "comprehensive_all_pass": False,
        "scope": {"chromosomes": ["chr1", "chr22"]},
        "per_chromosome": per_chrom,
        "totals": totals,
        "resources": {"outer": {"elapsed_seconds": 200}},
        "checks": {"fixture_only": True},
        "claim_ceiling": "fixture schema only",
    }
    summary_json = tmp_path / "summary.json"
    summary_json.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    summary_tsv = tmp_path / "summary.tsv"
    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(SUMMARY_REQUIRED_FIELDS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows([*per_chrom, totals])

    component_all = tmp_path / "component_all.tsv.gz"
    with gzip.open(component_all, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(COMPONENT_REQUIRED_FIELDS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(_component_rows())

    audit = tmp_path / "baseline-runtime-audit.md"
    audit.write_text(
        "\n".join(
            [
                "# Fixture copy of audited runtime statement",
                "",
                "高精度時間差為 5,086.484135464 秒。",
                "此值是 filesystem birth-timestamp proxy，不是單一 /usr/bin/time。",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return summary_json, summary_tsv, component_all, audit


def _sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _identity(path: Path):
    return {
        "path": str(path.resolve()),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _write_sha_sidecar(path: Path):
    sidecar = path.with_name(path.name + ".sha256")
    sidecar.write_text(f"{_sha256(path)}  {path.name}\n", encoding="utf-8")
    return sidecar


def _write_exact_log_fixture(tmp_path: Path):
    roots = []
    specifications = (
        (
            "chr1",
            "fixture_chr1_component_1",
            10,
            5,
            True,
            False,
            False,
            True,
        ),
        (
            "chr22",
            "fixture_chr22_component_1",
            12,
            7,
            False,
            False,
            True,
            False,
        ),
    )
    for (
        chrom,
        component_id,
        pre_cap_k,
        pattern_count,
        legacy_stable,
        corrected_stable,
        legacy_matches,
        changed,
    ) in specifications:
        root = tmp_path / f"exact_{chrom}"
        root.mkdir()
        source_identities = {}
        for source_name in (
            "partition_receipt",
            "legacy_components",
            "cut_constraints",
            "site_membership_coordinate_witness",
        ):
            source_path = root / f"source_{source_name}.dat"
            source_path.write_text(
                f"{chrom}:{source_name}:immutable\n", encoding="utf-8"
            )
            source_identities[source_name] = _identity(source_path)

        detail = root / "exact_log_sensitivity.tsv.gz"
        detail_row = {
            field: "" for field in EXACT_LOG_REQUIRED_FIELDS
        }
        detail_row.update(
            {
                "dataset": "HCC1395",
                "chrom": chrom,
                "legacy_component_id": component_id,
                "pre_cap_k": str(pre_cap_k),
                "exact_pattern_count": str(pattern_count),
                "legacy_log_matches_exact": str(legacy_matches).lower(),
                "legacy_weight_stable_reported": str(legacy_stable).lower(),
                "legacy_weight_stable_reconstructed": str(legacy_stable).lower(),
                "corrected_weight_stable": str(corrected_stable).lower(),
                "correction_changed_stability": str(changed).lower(),
                "remediation_class": (
                    "LEGACY_LOG_MATCHES_EXACT"
                    if legacy_matches
                    else "LEGACY_LOG_SUBOPTIMAL_EXACT_PRODUCT"
                ),
            }
        )
        with gzip.open(detail, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(EXACT_LOG_REQUIRED_FIELDS),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(detail_row)

        counts = {
            "components": 1,
            "sites": pre_cap_k,
            "patterns": pattern_count,
            "legacy_log_matches_exact": int(legacy_matches),
            "legacy_log_differs_exact": int(not legacy_matches),
            "legacy_weight_stable": int(legacy_stable),
            "corrected_weight_stable": int(corrected_stable),
            "correction_changed_stability": int(changed),
        }
        semantic_sha = hashlib.sha256(
            f"{chrom}:{component_id}:semantic".encode()
        ).hexdigest()
        summary_payload = {
            "schema_name": (
                "intersubmod.k_gt8_exact_log_sensitivity_remediation.summary"
            ),
            "schema_version": "0.1.0",
            "all_pass": True,
            "mode": "single_chromosome",
            "counts": counts,
            "chromosome_counts": {chrom: counts},
            "semantic_result_sha256": semantic_sha,
            "interpretation": {
                "primary_raw_partition_changed": False,
                "exact_product_scope": "fixture exact product only",
                "claim_ceiling": "fixture only",
            },
        }
        summary = root / "summary.json"
        summary.write_text(
            json.dumps(summary_payload, sort_keys=True) + "\n", encoding="utf-8"
        )
        summary_sha = _write_sha_sidecar(summary)
        receipt_payload = {
            "schema_name": "intersubmod.k_gt8_exact_log_sensitivity_remediation",
            "schema_version": "0.1.0",
            "all_pass": True,
            "mode": "single_chromosome",
            "scope": {
                "datasets": ["HCC1395"],
                "chromosomes": [chrom],
                "components": 1,
                "sites": pre_cap_k,
            },
            "parameters": {
                "max_block_size": 8,
                "primary_raw_partition_changed": False,
            },
            "checks": {
                "source_receipts_and_sha_verified": True,
                "constraint_count_and_mass_conserved": True,
            },
            "counts": counts,
            "chromosome_counts": {chrom: counts},
            "sources": {
                "full_root": None,
                "partitions": [
                    {
                        "dataset": "HCC1395",
                        "chrom": chrom,
                        **source_identities,
                    }
                ],
            },
            "outputs": {
                "exact_log_sensitivity": _identity(detail),
                "summary": _identity(summary),
                "summary_sha256": _identity(summary_sha),
            },
            "script": {
                "path": "/fixture/remediate_exact_log_sensitivity.py",
                "size_bytes": 123,
                "sha256": "a" * 64,
            },
            "semantic_result_sha256": semantic_sha,
        }
        receipt = root / "receipt.json"
        receipt.write_text(
            json.dumps(receipt_payload, sort_keys=True) + "\n", encoding="utf-8"
        )
        _write_sha_sidecar(receipt)
        success = {
            "schema_name": (
                "intersubmod.k_gt8_exact_log_sensitivity_remediation.success_marker"
            ),
            "schema_version": "0.1.0",
            "all_pass": True,
            "receipt": {
                "path": str(receipt.resolve()),
                "sha256": _sha256(receipt),
            },
        }
        (root / "_SUCCESS").write_text(
            json.dumps(success, sort_keys=True) + "\n", encoding="utf-8"
        )
        roots.append(root)
    return roots


def _write_hp_ps_unit_fixture(tmp_path: Path):
    root = tmp_path / "hp_ps_v5"
    root.mkdir()
    source_inputs = []
    for chrom in ("chr1", "chr22"):
        source = {"dataset": "HCC1395", "chrom": chrom}
        for source_name in (
            "partition_receipt",
            "legacy_components",
            "site_membership",
            "cut_constraints",
        ):
            path = root / f"{chrom}.{source_name}.dat"
            path.write_text(f"{chrom}:{source_name}:immutable\n", encoding="utf-8")
            source[source_name] = _identity(path)
        source_inputs.append(source)

    units = [
        {
            "dataset": "HCC1395",
            "chrom": "chr1",
            "legacy_component_id": "fixture_chr1_component_1",
            "hp_family": "1",
            "phase_set": "100",
            "unit_id": "HP1|PS100",
            "component_k": "10",
            "component_primary_active_site_count": "9",
            "unit_active_site_count": "5",
            "unit_active_site_fraction": "0.5",
            "total_pattern_rows": "5",
            "retained_pattern_rows": "1",
            "cut_lost_pattern_rows": "3",
            "unavoidable_pattern_rows": "1",
            "nonretained_pattern_rows": "4",
            "total_molecule_component_incidence_weight": "40",
            "retained_molecule_component_incidence_weight": "10",
            "cut_lost_molecule_component_incidence_weight": "20",
            "unavoidable_molecule_component_incidence_weight": "10",
            "nonretained_molecule_component_incidence_weight": "30",
            "retention_ratio": "0.25",
            "ratio_status": "OBSERVED_CONSTRAINT_DENOMINATOR",
            "support_stratum": "20-49",
            "eligible_headline": "true",
        },
        {
            "dataset": "HCC1395",
            "chrom": "chr1",
            "legacy_component_id": "fixture_chr1_component_1",
            "hp_family": "2",
            "phase_set": "100",
            "unit_id": "HP2|PS100",
            "component_k": "10",
            "component_primary_active_site_count": "9",
            "unit_active_site_count": "4",
            "unit_active_site_fraction": "0.4",
            "total_pattern_rows": "5",
            "retained_pattern_rows": "4",
            "cut_lost_pattern_rows": "1",
            "unavoidable_pattern_rows": "0",
            "nonretained_pattern_rows": "1",
            "total_molecule_component_incidence_weight": "40",
            "retained_molecule_component_incidence_weight": "30",
            "cut_lost_molecule_component_incidence_weight": "5",
            "unavoidable_molecule_component_incidence_weight": "5",
            "nonretained_molecule_component_incidence_weight": "10",
            "retention_ratio": "0.75",
            "ratio_status": "OBSERVED_CONSTRAINT_DENOMINATOR",
            "support_stratum": "20-49",
            "eligible_headline": "true",
        },
    ]
    unit_path = root / "HCC1395.hp_ps_observed_units.tsv.gz"
    with gzip.open(unit_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(HP_PS_UNIT_REQUIRED_FIELDS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(units)

    pair = {
        field: "" for field in HP_PS_PAIR_REQUIRED_FIELDS
    }
    pair.update(
        {
            "dataset": "HCC1395",
            "chrom": "chr1",
            "legacy_component_id": "fixture_chr1_component_1",
            "hp1_total_pattern_rows": "5",
            "hp1_total_molecule_component_incidence_weight": "40",
            "hp1_retained_molecule_component_incidence_weight": "10",
            "hp1_retention_ratio": "0.25",
            "hp2_total_pattern_rows": "5",
            "hp2_total_molecule_component_incidence_weight": "40",
            "hp2_retained_molecule_component_incidence_weight": "30",
            "hp2_retention_ratio": "0.75",
            "hp1_minus_hp2_retention_delta": "-0.5",
            "absolute_retention_delta": "0.5",
            "both_hp_headline_eligible": "true",
        }
    )
    pair_path = root / "HCC1395.hp1_hp2_paired_components.tsv.gz"
    with gzip.open(pair_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(HP_PS_PAIR_REQUIRED_FIELDS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(pair)

    coverage = {
        "components_in_partition_scope": 2,
        "components_with_any_observed_unit": 1,
        "components_without_observed_unit": 1,
        "components_hp1_only": 0,
        "components_hp2_only": 0,
        "components_hp1_and_hp2": 1,
        "by_chromosome": {
            "chr1": {
                "components_in_partition_scope": 1,
                "components_with_any_observed_unit": 1,
                "components_without_observed_unit": 0,
                "components_hp1_only": 0,
                "components_hp2_only": 0,
                "components_hp1_and_hp2": 1,
            },
            "chr22": {
                "components_in_partition_scope": 1,
                "components_with_any_observed_unit": 0,
                "components_without_observed_unit": 1,
                "components_hp1_only": 0,
                "components_hp2_only": 0,
                "components_hp1_and_hp2": 0,
            },
        },
    }
    def _slice(rows):
        total = sum(
            int(row["total_molecule_component_incidence_weight"]) for row in rows
        )
        retained = sum(
            int(row["retained_molecule_component_incidence_weight"])
            for row in rows
        )
        return {
            "observed_constraint_units": len(rows),
            "eligible_headline_units": sum(
                row["eligible_headline"] == "true" for row in rows
            ),
            "molecule_component_incidences": {
                "total": total,
                "retained": retained,
                "cut_lost": sum(
                    int(row["cut_lost_molecule_component_incidence_weight"])
                    for row in rows
                ),
                "unavoidable": sum(
                    int(row["unavoidable_molecule_component_incidence_weight"])
                    for row in rows
                ),
                "weighted_retention_ratio": str(retained / total),
            },
            "pattern_rows": {
                "total": sum(int(row["total_pattern_rows"]) for row in rows),
                "retained": sum(
                    int(row["retained_pattern_rows"]) for row in rows
                ),
                "cut_lost": sum(
                    int(row["cut_lost_pattern_rows"]) for row in rows
                ),
                "unavoidable": sum(
                    int(row["unavoidable_pattern_rows"]) for row in rows
                ),
            },
        }

    summary_payload = {
        "scope_contract": {
            "aggregation_weight": "molecule_x_component_incidence",
            "scope_ceiling": "observed_constraint_units_only",
            "unobserved_opportunity_policy": "no synthetic ratios",
        },
        "counts": {
            "components_in_partition_scope": 2,
            "components_with_observed_constraint_units": 1,
            "components_without_observed_constraint_units": 1,
            "observed_constraint_units": 2,
            "eligible_headline_units": 2,
        },
        "component_hp_unit_coverage": coverage,
        "by_chromosome": {"chr1": _slice(units)},
        "by_hp_family": {
            "HP1": _slice([units[0]]),
            "HP2": _slice([units[1]]),
        },
        "by_support_stratum": {"20-49": _slice(units)},
        "molecule_component_incidence_totals": {
            "total": 80,
            "retained": 40,
            "cut_lost": 25,
            "unavoidable": 15,
            "weighted_retention_ratio": "0.5",
        },
        "pattern_row_totals": {
            "total": 10,
            "retained": 5,
            "cut_lost": 4,
            "unavoidable": 1,
        },
        "retention_distribution_eligible_headline_units": {
            "n_units": 2,
            "quantiles": {"p25": "0.375", "median": "0.5"},
            "cumulative_threshold_counts": {"lt_0_5": 1, "lt_0_8": 2},
        },
        "hp1_hp2_component_paired": {
            "n_pairs": 1,
            "n_both_headline_eligible": 1,
            "absolute_delta_ge_0_25": 1,
            "absolute_delta_quantiles": {"median": "0.5"},
        },
    }
    summary = root / "summary.json"
    summary.write_text(
        json.dumps(summary_payload, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary_tsv = root / "summary.tsv"
    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "scope_type",
                "scope_id",
                "metric",
                "value",
                "unit",
                "denominator_note",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for metric in (
            "components_in_partition_scope",
            "components_with_any_observed_unit",
            "components_without_observed_unit",
            "components_hp1_only",
            "components_hp2_only",
            "components_hp1_and_hp2",
        ):
            writer.writerow(
                {
                    "scope_type": "component_coverage",
                    "scope_id": "overall",
                    "metric": metric,
                    "value": coverage[metric],
                    "unit": "count",
                    "denominator_note": "observed constraints only",
                }
            )

    receipt_payload = {
        "schema_name": "intersubmod.hp_ps_observed_constraint_retention_audit",
        "schema_version": "1.0.0",
        "all_pass": True,
        "scope": {
            "mode": "probe",
            "dataset": "HCC1395",
            "chromosomes_with_partition_targets": ["chr1", "chr22"],
            "skipped_zero_target_chromosomes": [],
            "unit_grain": (
                "dataset x chromosome x legacy_component_id x hp_family x "
                "exact known phase_set"
            ),
            "scope_ceiling": "observed_constraint_units_only",
        },
        "definitions": {
            "aggregation_weight": "molecule_x_component_incidence",
            "headline_eligible": "total_weight >= 20 and total_pattern_rows >= 5",
        },
        "parameters": {
            "eligible_min_total_weight": 20,
            "eligible_min_pattern_rows": 5,
        },
        "summary": summary_payload,
        "source_checks": {
            chrom: {name: True for name in HP_PS_SOURCE_CHECKS}
            for chrom in ("chr1", "chr22")
        },
        "checks": {
            "full_scope_is_exact_autosomes": True,
            "no_unobserved_unit_rows_synthesized": True,
            "source_checks_all_pass": True,
            "unit_molecule_incidence_mass_conserved": True,
            "unit_pattern_mass_conserved": True,
            "unit_rows_have_observed_denominator": True,
        },
        "tool": {
            "path": "/fixture/audit_hp_ps_unit_retention.py",
            "size_bytes": 456,
            "sha256": "b" * 64,
        },
        "inputs": source_inputs,
        "outputs": {
            "unit_table": _identity(unit_path),
            "paired_component_table": _identity(pair_path),
            "summary_json": _identity(summary),
            "summary_tsv": _identity(summary_tsv),
        },
        "semantic_result_sha256": "c" * 64,
    }
    receipt = root / "receipt.json"
    receipt.write_text(
        json.dumps(receipt_payload, sort_keys=True) + "\n", encoding="utf-8"
    )
    _write_sha_sidecar(receipt)
    return root


def _write_span_grid_fixture(tmp_path: Path):
    root = tmp_path / "span_grid"
    root.mkdir()
    settings = {
        50_000: {
            "retained": (70, 140),
            "blocks": (3, 5),
            "wall": (10.0, 20.0),
            "unavoidable": (5, 7),
        },
        100_000: {
            "retained": (80, 150),
            "blocks": (3, 3),
            "wall": (8.0, 12.0),
            "unavoidable": (3, 5),
        },
        200_000: {
            "retained": (85, 165),
            "blocks": (2, 3),
            "wall": (5.0, 10.0),
            "unavoidable": (2, 3),
        },
    }
    chrom_specs = (
        ("chr1", 1000, 1, 10, 10, 100),
        ("chr22", 543, 1, 12, 12, 200),
    )
    rows = []
    cap_totals = []
    for cap, values in settings.items():
        for index, (
            chrom,
            ssnv_sites,
            components,
            sites,
            max_k,
            raw_total,
        ) in enumerate(chrom_specs):
            retained = values["retained"][index]
            unavoidable = values["unavoidable"][index]
            rows.append(
                {
                    "span_cap_bp": cap,
                    "chrom": chrom,
                    "status": "COMPLETED",
                    "ssnv_sites": ssnv_sites,
                    "k_gt8_components": components,
                    "k_gt8_sites": sites,
                    "k_gt8_max_k": max_k,
                    "wall_seconds": values["wall"][index],
                    "new_blocks": values["blocks"][index],
                    "exact_patterns": 10 + index * 10,
                    "raw_total_molecule_weight": raw_total,
                    "raw_retained_molecule_weight": retained,
                    "raw_lost_molecule_weight": raw_total - retained,
                    "unavoidable_patterns": unavoidable,
                    "unavoidable_size_patterns": 1 + index,
                    "unavoidable_span_cap_patterns": unavoidable - 1,
                    "unavoidable_both_limits_patterns": 1,
                }
            )
        cap_rows = rows[-2:]
        cap_totals.append(
            {
                "span_cap_bp": cap,
                "receipt": {
                    "path": f"/fixture/span_{cap}/receipt.json",
                    "size_bytes": 1,
                    "sha256": "a" * 64,
                },
                "totals": {
                    "chromosomes": 2,
                    "completed_partition_chromosomes": 2,
                    "zero_target_skipped_chromosomes": 0,
                    "k_gt8_components": sum(
                        row["k_gt8_components"] for row in cap_rows
                    ),
                    "k_gt8_sites": sum(row["k_gt8_sites"] for row in cap_rows),
                    "cached_partition_wall_seconds": sum(
                        row["wall_seconds"] for row in cap_rows
                    ),
                    "new_blocks": sum(row["new_blocks"] for row in cap_rows),
                    "exact_patterns": sum(row["exact_patterns"] for row in cap_rows),
                    "raw_total_molecule_weight": sum(
                        row["raw_total_molecule_weight"] for row in cap_rows
                    ),
                    "raw_retained_molecule_weight": sum(
                        row["raw_retained_molecule_weight"] for row in cap_rows
                    ),
                    "raw_lost_molecule_weight": sum(
                        row["raw_lost_molecule_weight"] for row in cap_rows
                    ),
                    "unavoidable_patterns": sum(
                        row["unavoidable_patterns"] for row in cap_rows
                    ),
                    "unavoidable_size_patterns": sum(
                        row["unavoidable_size_patterns"] for row in cap_rows
                    ),
                    "unavoidable_span_cap_patterns": sum(
                        row["unavoidable_span_cap_patterns"] for row in cap_rows
                    ),
                    "unavoidable_both_limits_patterns": sum(
                        row["unavoidable_both_limits_patterns"] for row in cap_rows
                    ),
                },
            }
        )

    summary = root / "summary.tsv"
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    wall_total = sum(
        item["totals"]["cached_partition_wall_seconds"] for item in cap_totals
    )
    receipt = {
        "schema_name": "intersubmod.hcc1395_cached_span_sensitivity.run_receipt",
        "schema_version": "1.0.0",
        "all_pass": True,
        "comprehensive_all_pass": False,
        "sample": "HCC1395",
        "scope": {
            "chromosomes": ["chr1", "chr22"],
            "span_caps_bp": [50_000, 100_000, 200_000],
            "test_mode": True,
        },
        "timing_scope": (
            "cached_partition_wall_seconds excludes source identity "
            "validation and excludes upstream BAM extraction"
        ),
        "totals": {
            "span_caps": 3,
            "chromosome_cap_tasks": 6,
            "completed_partition_tasks": 6,
            "zero_target_skipped_tasks": 0,
            "cached_partition_wall_seconds": wall_total,
        },
        "caps": cap_totals,
        "outputs": {
            "summary": {
                "path": str(summary.resolve()),
                "size_bytes": summary.stat().st_size,
                "sha256": _sha256(summary),
            }
        },
        "claim_ceiling": "fixture hard-cap sensitivity only",
    }
    (root / "receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return root


def _argv(tmp_path: Path, *, fixture: bool = True):
    summary_json, summary_tsv, component_all, audit = _write_fixture(tmp_path)
    output = tmp_path / "artifact.json"
    argv = [
        "--summary-json",
        str(summary_json),
        "--summary-tsv",
        str(summary_tsv),
        "--component-all",
        str(component_all),
        "--baseline-runtime-audit",
        str(audit),
        "--output",
        str(output),
        "--generated-at",
        GENERATED_AT,
        "--top-components",
        "2",
    ]
    if fixture:
        argv.append("--fixture")
    return argv, output


def test_chr1_chr22_fixture_builds_canonical_report_shape(tmp_path, capsys):
    argv, output = _argv(tmp_path)

    assert main(argv) == 0
    receipt = json.loads(capsys.readouterr().out)
    artifact = json.loads(output.read_text(encoding="utf-8"))

    assert receipt["snapshot_status"] == "fixture"
    assert artifact["surface"] == "report"
    assert artifact["snapshot"]["status"] == "fixture"
    assert artifact["snapshot"]["generatedAt"] == GENERATED_AT
    assert artifact["manifest"]["title"].startswith("[FIXTURE]")
    assert artifact["manifest"]["blocks"][0] == {
        "id": "title",
        "type": "markdown",
        "body": f"# {artifact['manifest']['title']}",
        "layout": "full",
    }
    assert artifact["manifest"]["blocks"][1]["body"].startswith(
        "## Executive Summary"
    )
    assert "不可當作 HCC1395 chr1–chr22 最終結果" in artifact["manifest"][
        "blocks"
    ][1]["body"]

    datasets = artifact["snapshot"]["datasets"]
    assert all(isinstance(rows, list) for rows in datasets.values())
    assert all(len(rows) <= 2000 for rows in datasets.values())
    assert [row["chrom"] for row in datasets["per_chromosome_metrics"]] == [
        "chr1",
        "chr22",
    ]
    assert len(datasets["retention_by_chrom_method"]) == 4
    assert {
        row["method"] for row in datasets["retention_by_chrom_method"]
    } == {"densest-8 counterfactual", "新 read-supported"}
    runtime_rows = datasets["runtime_by_chrom_stage"]
    assert len(runtime_rows) == 8
    assert {row["stage"] for row in runtime_rows} == {
        "Extract",
        "Partition",
        "Load",
        "Loop",
    }
    assert {row["stage_detail"] for row in runtime_rows} == {
        "read-linkage extraction",
        "partition stage total",
        "pattern load + aggregate",
        "three-weight partition component loop",
    }
    assert all(row["fixture_value"] is True for row in runtime_rows)
    assert all("fixture synthetic" in row["timing_scope"] for row in runtime_rows)
    assert sum(
        row["wall_minutes"]
        for row in runtime_rows
        if row["stage"] == "Load"
    ) == pytest.approx(1.5 / 60)
    assert sum(
        row["wall_minutes"]
        for row in runtime_rows
        if row["stage"] == "Loop"
    ) == pytest.approx(0.7 / 60)
    assert len(datasets["genomic_components"]) == 2
    assert {
        "start_mb",
        "end_mb",
        "span_bp",
        "pre_cap_k",
        "primary_active_sites",
        "weight_stability_status",
        "new_read_retention_ratio",
        "old_read_retention_ratio",
    }.issubset(datasets["genomic_components"][0])
    assert "span_sensitivity" not in datasets

    charts = {chart["id"]: chart for chart in artifact["manifest"]["charts"]}
    assert set(charts) == {
        "chart_read_retention",
        "chart_target_sites",
        "chart_genomic_components",
        "chart_runtime",
    }
    assert charts["chart_read_retention"]["encodings"]["color"]["field"] == "method"
    assert charts["chart_read_retention"]["settings"]["groupMode"] == "grouped"
    assert charts["chart_runtime"]["encodings"]["color"]["field"] == "stage"
    genomic = charts["chart_genomic_components"]
    assert genomic["type"] == "scatter"
    assert genomic["encodings"]["x"]["field"] == "start_mb"
    assert genomic["encodings"]["y"]["field"] == "chromosome_number"
    assert genomic["encodings"]["size"]["field"] == "pre_cap_k"
    assert genomic["encodings"]["color"]["field"] == "weight_stability_status"
    assert any(
        block["id"] == "span_not_provided"
        and "--span-grid-summary" in block["body"]
        for block in artifact["manifest"]["blocks"]
    )
    headline = datasets["headline_metrics"][0]
    assert headline["primary_active_sites"] == 19
    assert headline["tree_ready_local_blocks"] == 3
    assert headline["zero_evidence_blocks"] == 1
    assert headline["zero_evidence_block_sites"] == 4
    assert headline["abstain_blocks"] == 1
    assert headline["partition_pattern_load_aggregate_minutes"] == 1.5 / 60
    assert (
        headline["three_weight_partition_component_loop_minutes"]
        == 0.7 / 60
    )
    expected_stage_wall_seconds = 60.0 + 2.0 + 120.0 + 4.0
    assert headline["runner_overhead_minutes"] == pytest.approx(
        (200.0 - expected_stage_wall_seconds) / 60.0
    )
    assert headline["unavoidable_patterns"] == 5
    assert headline["unavoidable_n_fixed_ra_gt8_patterns"] == 3
    assert headline["unavoidable_n_fixed_ra_lte8_span_gt8_patterns"] == 2
    assert any(
        "Unavoidable 不只有一種原因" in block.get("body", "")
        and "n_fixed_ra>8" in block["body"]
        and "n_fixed_ra≤8" in block["body"]
        for block in artifact["manifest"]["blocks"]
    )
    rendered_text = json.dumps(
        artifact["manifest"], ensure_ascii=False, sort_keys=True
    )
    assert "three-weight partition component loop" in rendered_text
    assert "pure dp" not in rendered_text.lower()
    assert "component-level HP1/HP2 × known-PS" in rendered_text
    assert "不是全基因組 unique reads" in rendered_text
    assert "不是舊 v5 的實測 retention" in rendered_text

    source_ids = {source["id"] for source in artifact["manifest"]["sources"]}
    for collection in ("cards", "charts", "tables"):
        assert all(
            item["sourceId"] in source_ids
            for item in artifact["manifest"][collection]
        )
    assert all(
        not Path(source["path"]).is_absolute()
        for source in artifact["manifest"]["sources"]
        if "path" in source
    )
    for table in artifact["manifest"]["tables"]:
        column_fields = {column["field"] for column in table["columns"]}
        assert table["defaultSort"]["field"] in column_fields


def test_subset_input_fails_closed_without_fixture_flag(tmp_path, capsys):
    argv, output = _argv(tmp_path, fixture=False)

    assert main(argv) == 2
    assert not output.exists()
    assert "formal report requires" in capsys.readouterr().err


def test_builder_refuses_to_overwrite_artifact(tmp_path, capsys):
    argv, output = _argv(tmp_path)

    assert main(argv) == 0
    capsys.readouterr()
    original = output.read_bytes()

    assert main(argv) == 2
    assert output.read_bytes() == original
    assert "refusing to overwrite" in capsys.readouterr().err


def test_builder_fails_closed_when_internal_partition_timing_is_missing(
    tmp_path, capsys
):
    argv, output = _argv(tmp_path)
    summary_path = Path(argv[argv.index("--summary-json") + 1])
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["per_chromosome"][0].pop(
        "partition_ordered_hypergraph_dp_seconds"
    )
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    assert main(argv) == 2
    assert not output.exists()
    assert "partition_ordered_hypergraph_dp_seconds" in capsys.readouterr().err


def test_chr21_skip_requires_explicit_zero_internal_partition_timings():
    row = {field: 0 for field in SUMMARY_REQUIRED_FIELDS}
    row.update(
        {
            "chrom": "chr21",
            "partition_stage_status": "SKIP_NO_K_GT8_TARGET",
            "ssnv_sites": 436,
            "new_retention_ratio": "",
            "old_densest8_retention_ratio": "",
            "extraction_wall_seconds": 12.0,
            "partition_wall_seconds": "",
            "partition_pattern_load_aggregate_seconds": 0.0,
            "partition_ordered_hypergraph_dp_seconds": 0.0,
        }
    )

    validate_summary_row(row, "chr21")

    row["partition_ordered_hypergraph_dp_seconds"] = 0.01
    with pytest.raises(ReportInputError, match="zero internal timings"):
        validate_summary_row(row, "chr21")


def test_builder_fails_closed_on_unavoidable_mechanism_drift(tmp_path, capsys):
    argv, output = _argv(tmp_path)
    summary_path = Path(argv[argv.index("--summary-json") + 1])
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["per_chromosome"][0][
        "unavoidable_n_fixed_ra_gt8_patterns"
    ] += 1
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    assert main(argv) == 2
    assert not output.exists()
    assert "unavoidable pattern-mechanism conservation" in capsys.readouterr().err


def test_optional_span_grid_adds_three_comparison_charts(tmp_path, capsys):
    argv, output = _argv(tmp_path)
    span_root = _write_span_grid_fixture(tmp_path)
    argv.extend(["--span-grid-summary", str(span_root)])

    assert main(argv) == 0
    capsys.readouterr()
    artifact = json.loads(output.read_text(encoding="utf-8"))

    span_rows = artifact["snapshot"]["datasets"]["span_sensitivity"]
    assert [row["cap_label"] for row in span_rows] == [
        "No cap（k≤8 only）",
        "50 kb",
        "100 kb",
        "200 kb",
    ]
    assert [row["raw_total_molecule_weight"] for row in span_rows] == [
        300,
        300,
        300,
        300,
    ]
    assert span_rows[0]["cached_partition_wall_seconds"] == 6
    assert [row["new_blocks"] for row in span_rows] == [4, 8, 6, 5]
    assert [row["read_retention_ratio"] for row in span_rows] == [
        260 / 300,
        210 / 300,
        230 / 300,
        250 / 300,
    ]

    charts = {chart["id"]: chart for chart in artifact["manifest"]["charts"]}
    for chart_id, y_field in (
        ("chart_span_retention", "read_retention_ratio"),
        ("chart_span_blocks", "new_blocks"),
        ("chart_span_runtime", "cached_partition_wall_minutes"),
    ):
        assert charts[chart_id]["dataset"] == "span_sensitivity"
        assert charts[chart_id]["encodings"]["x"]["field"] == "cap_label"
        assert charts[chart_id]["encodings"]["y"]["field"] == y_field
        assert charts[chart_id]["sourceId"] == "src_span_sensitivity"
    assert any(
        source["id"] == "src_span_sensitivity"
        for source in artifact["manifest"]["sources"]
    )
    assert not any(
        block["id"] == "span_not_provided"
        for block in artifact["manifest"]["blocks"]
    )


def test_span_grid_summary_identity_drift_fails_closed(tmp_path, capsys):
    argv, output = _argv(tmp_path)
    span_root = _write_span_grid_fixture(tmp_path)
    with (span_root / "summary.tsv").open("a", encoding="utf-8") as handle:
        handle.write("\n")
    argv.extend(["--span-grid-summary", str(span_root / "receipt.json")])

    assert main(argv) == 2
    assert not output.exists()
    assert "identity mismatch" in capsys.readouterr().err


def test_audited_fixture_uses_corrected_stability_and_hp_ps_tail(
    tmp_path, capsys
):
    argv, output = _argv(tmp_path)
    exact_roots = _write_exact_log_fixture(tmp_path)
    hp_root = _write_hp_ps_unit_fixture(tmp_path)
    for root in exact_roots:
        argv.extend(["--exact-log-audit", str(root)])
    argv.extend(["--hp-ps-unit-audit", str(hp_root)])

    assert main(argv) == 0
    receipt = json.loads(capsys.readouterr().out)
    artifact = json.loads(output.read_text(encoding="utf-8"))
    datasets = artifact["snapshot"]["datasets"]

    headline = datasets["headline_metrics"][0]
    assert headline["legacy_weight_stable_components"] == 1
    assert headline["corrected_weight_stable_components"] == 0
    assert headline["weight_stable_components"] == 0
    assert headline["legacy_log_differs_exact"] == 1
    assert headline["correction_changed_stability"] == 1
    assert headline["primary_raw_partition_changed"] is False
    chr1 = next(
        row for row in datasets["genomic_components"] if row["chrom"] == "chr1"
    )
    assert chr1["legacy_weight_stable"] is True
    assert chr1["corrected_weight_stable"] is False
    assert chr1["weight_stable"] is False
    assert chr1["weight_stability_status"] == "Sensitive"
    assert chr1["evidence_status_detail"] == "Non-ABSTAIN / weight-sensitive"

    hp_metrics = datasets["hp_ps_summary_metrics"][0]
    assert hp_metrics["components_in_partition_scope"] == 2
    assert hp_metrics["components_with_observed_constraint_units"] == 1
    assert hp_metrics["components_without_observed_constraint_units"] == 1
    assert hp_metrics["observed_constraint_units"] == 2
    assert hp_metrics["eligible_headline_units"] == 2
    assert hp_metrics["weighted_retention_ratio"] == pytest.approx(0.5)
    assert hp_metrics["eligible_retention_p10"] == pytest.approx(0.30)
    assert hp_metrics["eligible_retention_p25"] == pytest.approx(0.375)
    assert hp_metrics["eligible_retention_median"] == pytest.approx(0.5)
    assert hp_metrics["eligible_retention_lt_0_5"] == 1
    assert hp_metrics["eligible_retention_lt_0_8"] == 2
    assert hp_metrics["hp1_minus_hp2_weighted_retention_delta"] == pytest.approx(
        -0.5
    )
    assert hp_metrics["paired_component_absolute_delta_median"] == pytest.approx(
        0.5
    )

    distribution = datasets["hp_ps_unit_distribution"]
    assert len(distribution) == 22
    assert sum(row["unit_count"] for row in distribution) == 2
    worst = datasets["hp_ps_worst_units"]
    assert len(worst) == 2
    assert len(worst) <= 25
    assert [row["retention_ratio"] for row in worst] == [0.25, 0.75]
    assert all(row["component_id"] != "fixture_chr22_component_1" for row in worst)

    charts = {chart["id"]: chart for chart in artifact["manifest"]["charts"]}
    hp_chart = charts["chart_hp_ps_unit_distribution"]
    assert hp_chart["dataset"] == "hp_ps_unit_distribution"
    assert hp_chart["encodings"]["x"]["field"] == "retention_bucket"
    assert hp_chart["encodings"]["color"]["field"] == "eligibility_status"
    tables = {table["id"]: table for table in artifact["manifest"]["tables"]}
    hp_table = tables["table_hp_ps_worst_units"]
    assert hp_table["dataset"] == "hp_ps_worst_units"
    assert hp_table["defaultSort"] == {
        "field": "retention_ratio",
        "direction": "asc",
    }
    rendered = json.dumps(artifact["manifest"], ensure_ascii=False)
    assert "Primary raw-molecule cuts 完全未改動" in rendered
    assert "0% 或 100% retention" in rendered
    assert "每個位點 VAF" in rendered

    assert receipt["dataset_rows"]["hp_ps_unit_distribution"] == 22
    assert receipt["dataset_rows"]["hp_ps_worst_units"] == 2
    assert receipt["artifact_profile"] == {
        "datasets": 9,
        "rows": sum(receipt["dataset_rows"].values()),
        "blocks": len(artifact["manifest"]["blocks"]),
        "cards": len(artifact["manifest"]["cards"]),
        "charts": 5,
        "tables": 3,
        "sources": len(artifact["manifest"]["sources"]),
    }
    assert receipt["portable_contract"]["exact_log_audit"] == "authenticated"
    assert receipt["portable_contract"]["hp_ps_unit_audit"] == "authenticated"


def test_exact_log_source_identity_tamper_fails_closed(tmp_path, capsys):
    argv, output = _argv(tmp_path)
    exact_roots = _write_exact_log_fixture(tmp_path)
    source = exact_roots[0] / "source_cut_constraints.dat"
    source.write_text("tampered\n", encoding="utf-8")
    for root in exact_roots:
        argv.extend(["--exact-log-audit", str(root / "summary.json")])

    assert main(argv) == 2
    assert not output.exists()
    assert "recorded" in capsys.readouterr().err


def test_hp_ps_receipt_identity_tamper_fails_closed(tmp_path, capsys):
    argv, output = _argv(tmp_path)
    hp_root = _write_hp_ps_unit_fixture(tmp_path)
    with (hp_root / "summary.tsv").open("a", encoding="utf-8") as handle:
        handle.write("\n")
    argv.extend(["--hp-ps-unit-audit", str(hp_root / "receipt.json")])

    assert main(argv) == 2
    assert not output.exists()
    assert "identity mismatch" in capsys.readouterr().err
