#!/usr/bin/env python3
"""Build the canonical multi-sample drilldown dashboard artifact.

The artifact is a bounded, source-backed snapshot for the shared Data Analytics
dashboard renderer.  This builder deliberately stops at ``artifact.json``; HTML
packaging is a separate delivery step.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parents[1]
DEFAULT_RESULTS = PROJECT_ROOT / "results"
DEFAULT_OUTPUT = PROJECT_ROOT / "multi_bam_dashboard_artifact.json"
ALL_DATASETS_SCOPE = "All"
EXPECTED_DATASET_N = 7
EXPECTED_BIOLOGICAL_N = 6

INPUT_FILES = {
    "cohort": "cohort_topology_metrics.csv",
    "bundle": "bundle_overview.csv",
    "coverage": "methylation_coverage_by_k.csv",
    "axis": "methylation_axis_metrics.csv",
    "audit": "metrics_audit.json",
    "visual": "legacy_browser_visual_audit.generated.json",
}


class ContractError(ValueError):
    """Raised when an input or artifact violates the canonical contract."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS,
        help=f"Audited input directory (default: {DEFAULT_RESULTS})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Canonical artifact JSON path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def read_csv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise ContractError(f"required input is missing: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = sorted(set(required) - fields)
        if missing:
            raise ContractError(f"{path.name} is missing columns: {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ContractError(f"{path.name} has no data rows")
    return rows


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ContractError(f"required input is missing: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ContractError(f"{path.name} must contain a JSON object")
    return value


def integer(value: Any, label: str) -> int:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{label} must be an integer, got {value!r}") from exc
    if not math.isfinite(parsed) or not parsed.is_integer():
        raise ContractError(f"{label} must be an integer, got {value!r}")
    return int(parsed)


def number(value: Any, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{label} must be numeric, got {value!r}") from exc
    if not math.isfinite(parsed):
        raise ContractError(f"{label} must be finite, got {value!r}")
    return parsed


def boolean(value: Any, label: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise ContractError(f"{label} must be true/false, got {value!r}")


def rate_from_percent(value: Any, label: str) -> float:
    rate = number(value, label) / 100.0
    if not 0.0 <= rate <= 1.0:
        raise ContractError(f"{label} must be between 0 and 100, got {value!r}")
    return rate


def ratio(numerator: int, denominator: int, label: str) -> float:
    if denominator <= 0:
        raise ContractError(f"{label} denominator must be positive")
    value = numerator / denominator
    if not 0.0 <= value <= 1.0:
        raise ContractError(f"{label} must be between 0 and 1, got {value}")
    return value


def assert_close(actual: float, expected: float, label: str, tolerance: float = 1e-6) -> None:
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise ContractError(f"{label} mismatch: stored={actual}, recomputed={expected}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return resolved.as_posix()


def scoped_rows(rows: Iterable[dict[str, Any]], sample_field: str = "sample") -> list[dict[str, Any]]:
    """Expose rows in the cohort view and their exact dataset-specific view."""
    output: list[dict[str, Any]] = []
    for row in rows:
        sample = str(row[sample_field])
        output.append({"sample_filter": ALL_DATASETS_SCOPE, **row})
        output.append({"sample_filter": sample, **row})
    return output


def hcc_scoped_rows(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Expose HCC1395-only downstream rows without fabricating other samples."""
    output: list[dict[str, Any]] = []
    for row in rows:
        output.append({"sample_filter": ALL_DATASETS_SCOPE, **row})
        output.append({"sample_filter": "HCC1395", **row})
    return output


def make_snapshot_source(
    source_id: str,
    label: str,
    dataset: str,
    input_paths: list[Path],
    generated_at: str,
    description: str,
    metric_definitions: list[str],
    order_by: str,
    filters: list[str] | None = None,
) -> dict[str, Any]:
    relative_inputs = [source_path(path) for path in input_paths]
    return {
        "id": source_id,
        "label": label,
        "path": relative_inputs[0],
        "inputPaths": relative_inputs,
        "inputSha256": {source_path(path): sha256(path) for path in input_paths},
        "query": {
            "engine": "portable artifact snapshot",
            "language": "sql",
            "sql": f"SELECT * FROM {dataset} ORDER BY {order_by};",
            "description": description,
            "tables_used": [dataset, *relative_inputs],
            "filters": filters or ["No source rows excluded unless stated in the description."],
            "metric_definitions": metric_definitions,
            "id": f"build_multi_bam_dashboard_artifact.py::{dataset}",
            "executed_at": generated_at,
        },
    }


def build_artifact(results_dir: Path) -> dict[str, Any]:
    paths = {key: results_dir / name for key, name in INPUT_FILES.items()}
    cohort_raw = read_csv(
        paths["cohort"],
        {
            "sample", "biological_id", "technical_replicate", "region_n",
            "distinct_ssnv_n", "tree_n", "tree_pct_all_regions", "unique_tree_n",
            "unique_pct_among_tree", "hash_match", "sample_identity_all_match",
            "receipt_all_pass", "family_complete_n", "objective_certified_n",
            "receipt_all_mutation_bearing_families_complete",
            "receipt_all_units_objective_certified",
        },
    )
    bundle_raw = read_csv(
        paths["bundle"],
        {
            "bundle", "sample", "built_at", "validation_status_recomputed", "regions",
            "regions_with_tree", "tree_coverage_pct_all_regions", "unique_best_tree",
            "unique_pct_among_tree", "distinct_ssnv", "methyl_summary_rows",
            "methyl_topology_linked", "methyl_topology_linkage_pct", "selfcheck_pass",
            "selfcheck_fail", "selfcheck_skip", "panel_all_actual_bytes",
            "panel_unreported_bytes", "cluster_dispersion_warning_pct",
        },
    )
    coverage_raw = read_csv(
        paths["coverage"],
        {"bundle", "k_bin", "numerator", "denominator", "percent", "formula", "unit"},
    )
    axis_raw = read_csv(
        paths["axis"],
        {
            "bundle", "scope", "axis", "tested_n", "raw_p_le_0_05_n",
            "raw_p_le_0_05_pct", "bh_fdr_q_le_0_05_n", "bh_fdr_q_le_0_05_pct",
            "effect_median", "effect_q1", "effect_q3", "effect_unit",
            "interpretation_gate", "multiplicity_family", "source_effect_field",
        },
    )
    audit = read_json(paths["audit"])
    visual = read_json(paths["visual"])

    if audit.get("schema_name") != "intersubmod.drilldown.audit.metrics":
        raise ContractError("metrics_audit.json has an unexpected schema_name")
    if visual.get("schema_name") != "intersubmod.legacy_browser_visual_audit":
        raise ContractError("legacy visual audit has an unexpected schema_name")

    cohort: list[dict[str, Any]] = []
    seen_samples: set[str] = set()
    for row in cohort_raw:
        sample = row["sample"]
        if sample in seen_samples:
            raise ContractError(f"duplicate cohort sample: {sample}")
        seen_samples.add(sample)
        regions = integer(row["region_n"], f"{sample}.region_n")
        tree_n = integer(row["tree_n"], f"{sample}.tree_n")
        unique_n = integer(row["unique_tree_n"], f"{sample}.unique_tree_n")
        tree_coverage = rate_from_percent(
            row["tree_pct_all_regions"], f"{sample}.tree_pct_all_regions"
        )
        unique_rate = rate_from_percent(
            row["unique_pct_among_tree"], f"{sample}.unique_pct_among_tree"
        )
        assert_close(tree_coverage, ratio(tree_n, regions, f"{sample}.tree coverage"),
                     f"{sample}.tree coverage")
        assert_close(unique_rate, ratio(unique_n, tree_n, f"{sample}.unique rate"),
                     f"{sample}.unique rate")
        hash_match = boolean(row["hash_match"], f"{sample}.hash_match")
        identity_match = boolean(
            row["sample_identity_all_match"], f"{sample}.sample_identity_all_match"
        )
        receipt_pass = boolean(row["receipt_all_pass"], f"{sample}.receipt_all_pass")
        all_families_complete = boolean(
            row["receipt_all_mutation_bearing_families_complete"],
            f"{sample}.receipt_all_mutation_bearing_families_complete",
        )
        all_objective_certified = boolean(
            row["receipt_all_units_objective_certified"],
            f"{sample}.receipt_all_units_objective_certified",
        )
        cohort.append({
            "sample": sample,
            "biological_id": row["biological_id"],
            "technical_replicate": "yes" if boolean(
                row["technical_replicate"], f"{sample}.technical_replicate"
            ) else "no",
            "regions": regions,
            "distinct_ssnv": integer(row["distinct_ssnv_n"], f"{sample}.distinct_ssnv_n"),
            "tree_n": tree_n,
            "tree_coverage": tree_coverage,
            "unique_tree_n": unique_n,
            "unique_among_tree": unique_rate,
            "family_complete_n": integer(
                row["family_complete_n"], f"{sample}.family_complete_n"
            ),
            "objective_certified_n": integer(
                row["objective_certified_n"], f"{sample}.objective_certified_n"
            ),
            "hash_match": "PASS" if hash_match else "FAIL",
            "sample_identity_match": "PASS" if identity_match else "FAIL",
            "receipt_technical_pass": "PASS" if receipt_pass else "FAIL",
            "technical_gate_pass": receipt_pass,
            "all_families_complete": all_families_complete,
            "all_families_complete_gate": "PASS" if all_families_complete else "FAIL",
            "all_objective_certified": all_objective_certified,
            "all_objective_certified_gate": "PASS" if all_objective_certified else "FAIL",
            "integrity_pass": hash_match and identity_match and receipt_pass,
        })

    biological_ids = {row["biological_id"] for row in cohort}
    if len(cohort) != EXPECTED_DATASET_N or len(biological_ids) != EXPECTED_BIOLOGICAL_N:
        raise ContractError(
            "cohort scope mismatch: "
            f"expected {EXPECTED_DATASET_N} datasets/{EXPECTED_BIOLOGICAL_N} biological, "
            f"got {len(cohort)}/{len(biological_ids)}"
        )

    bundle_samples = {row["sample"] for row in bundle_raw}
    if bundle_samples != {"HCC1395"}:
        raise ContractError(
            "downstream bundle scope changed; expected HCC1395 only, got "
            + ", ".join(sorted(bundle_samples))
        )

    bundle_overview: list[dict[str, Any]] = []
    for row in bundle_raw:
        bundle = row["bundle"]
        bundle_overview.append({
            "sample": row["sample"],
            "bundle": bundle,
            "built_at": row["built_at"],
            "validation_status": row["validation_status_recomputed"],
            "regions": integer(row["regions"], f"{bundle}.regions"),
            "regions_with_tree": integer(
                row["regions_with_tree"], f"{bundle}.regions_with_tree"
            ),
            "tree_coverage": rate_from_percent(
                row["tree_coverage_pct_all_regions"], f"{bundle}.tree_coverage"
            ),
            "unique_tree_n": integer(row["unique_best_tree"], f"{bundle}.unique_best_tree"),
            "unique_among_tree": rate_from_percent(
                row["unique_pct_among_tree"], f"{bundle}.unique_among_tree"
            ),
            "distinct_ssnv": integer(row["distinct_ssnv"], f"{bundle}.distinct_ssnv"),
            "methyl_summary_rows": integer(
                row["methyl_summary_rows"], f"{bundle}.methyl_summary_rows"
            ),
            "ism_linked": integer(
                row["methyl_topology_linked"], f"{bundle}.methyl_topology_linked"
            ),
            "ism_linkage": rate_from_percent(
                row["methyl_topology_linkage_pct"], f"{bundle}.methyl_topology_linkage"
            ),
            "selfcheck_pass": integer(row["selfcheck_pass"], f"{bundle}.selfcheck_pass"),
            "selfcheck_fail": integer(row["selfcheck_fail"], f"{bundle}.selfcheck_fail"),
            "selfcheck_skip": integer(row["selfcheck_skip"], f"{bundle}.selfcheck_skip"),
            "panel_actual_bytes": integer(
                row["panel_all_actual_bytes"], f"{bundle}.panel_all_actual_bytes"
            ),
            "panel_unreported_bytes": integer(
                row["panel_unreported_bytes"], f"{bundle}.panel_unreported_bytes"
            ),
            "cluster_warning_rate": rate_from_percent(
                row["cluster_dispersion_warning_pct"], f"{bundle}.cluster_warning_rate"
            ),
        })

    k_order = {str(i): i for i in range(1, 8)} | {"8+": 8}
    coverage: list[dict[str, Any]] = []
    for row in coverage_raw:
        if row["k_bin"] not in k_order:
            raise ContractError(f"unexpected active-k bin: {row['k_bin']}")
        numerator = integer(row["numerator"], f"coverage.{row['bundle']}.{row['k_bin']}.n")
        denominator = integer(
            row["denominator"], f"coverage.{row['bundle']}.{row['k_bin']}.denominator"
        )
        stored_rate = rate_from_percent(
            row["percent"], f"coverage.{row['bundle']}.{row['k_bin']}.percent"
        )
        assert_close(
            stored_rate,
            ratio(numerator, denominator, "ISM coverage"),
            f"coverage.{row['bundle']}.{row['k_bin']}",
        )
        coverage.append({
            "sample": "HCC1395",
            "bundle": row["bundle"],
            "k_bin": row["k_bin"],
            "k_order": k_order[row["k_bin"]],
            "linked": numerator,
            "denominator": denominator,
            "coverage": stored_rate,
            "formula": row["formula"],
            "unit": row["unit"],
        })
    coverage.sort(key=lambda item: (item["k_order"], item["bundle"]))

    axis_details: list[dict[str, Any]] = []
    axis_rates_base: list[dict[str, Any]] = []
    for row in axis_raw:
        label = f"axis.{row['bundle']}.{row['scope']}.{row['axis']}"
        tested_n = integer(row["tested_n"], f"{label}.tested_n")
        raw_n = integer(row["raw_p_le_0_05_n"], f"{label}.raw_n")
        bh_n = integer(row["bh_fdr_q_le_0_05_n"], f"{label}.bh_n")
        raw_rate = rate_from_percent(row["raw_p_le_0_05_pct"], f"{label}.raw_rate")
        bh_rate = rate_from_percent(row["bh_fdr_q_le_0_05_pct"], f"{label}.bh_rate")
        assert_close(raw_rate, ratio(raw_n, tested_n, f"{label}.raw"), f"{label}.raw")
        assert_close(bh_rate, ratio(bh_n, tested_n, f"{label}.bh"), f"{label}.bh")
        item = {
            "sample": "HCC1395",
            "bundle": row["bundle"],
            "analysis_scope": row["scope"],
            "axis": row["axis"],
            "tested_n": tested_n,
            "raw_sig_n": raw_n,
            "raw_sig_rate": raw_rate,
            "bh_sig_n": bh_n,
            "bh_sig_rate": bh_rate,
            "effect_median": number(row["effect_median"], f"{label}.effect_median"),
            "effect_q1": number(row["effect_q1"], f"{label}.effect_q1"),
            "effect_q3": number(row["effect_q3"], f"{label}.effect_q3"),
            "effect_unit": row["effect_unit"],
            "effect_field": row["source_effect_field"],
            "multiplicity_family": row["multiplicity_family"],
            "interpretation_gate": row["interpretation_gate"],
        }
        axis_details.append(item)
        if row["scope"] == "topology_linked_distinct_ssnv":
            for rate_type, significant_n, significant_rate in (
                ("raw p≤0.05", raw_n, raw_rate),
                ("BH q≤0.05", bh_n, bh_rate),
            ):
                axis_rates_base.append({
                    "sample": "HCC1395",
                    "bundle": row["bundle"],
                    "axis": row["axis"],
                    "rate_type": rate_type,
                    "series": f"{row['bundle']} · {rate_type}",
                    "significant_n": significant_n,
                    "tested_n": tested_n,
                    "rate": significant_rate,
                    "interpretation_gate": row["interpretation_gate"],
                })

    lca_raw = audit.get("lca_ab_summary")
    if not isinstance(lca_raw, dict):
        raise ContractError("metrics_audit.json is missing lca_ab_summary")
    shared_chrom = integer(lca_raw.get("shared_chrom_n"), "lca.shared_chrom_n")
    same_in_bam = integer(lca_raw.get("same_in_bam_n"), "lca.same_in_bam_n")
    same_threads = integer(lca_raw.get("same_threads_n"), "lca.same_threads_n")
    lca_summary = [{
        "sample": "HCC1395",
        "shared_chromosomes": shared_chrom,
        "same_in_bam_n": same_in_bam,
        "same_in_bam_rate": ratio(same_in_bam, shared_chrom, "LCA in_bam identity"),
        "same_threads_n": same_threads,
        "same_threads_rate": ratio(same_threads, shared_chrom, "LCA thread identity"),
        "pre_lv_written": integer(lca_raw.get("pre_lv_written"), "lca.pre_lv_written"),
        "post_lv_written": integer(lca_raw.get("post_lv_written"), "lca.post_lv_written"),
        "net_new_lv_written": integer(
            lca_raw.get("net_new_lv_written"), "lca.net_new_lv_written"
        ),
        "lca_resolved": integer(lca_raw.get("lca_resolved"), "lca.lca_resolved"),
        "causal_ab_gate": str(lca_raw.get("causal_ab_gate", "UNKNOWN")),
    }]

    visual_rows: list[dict[str, Any]] = []
    for surface in ("legacy", "current"):
        for viewport in ("desktop", "mobile"):
            key = f"{surface}_{viewport}"
            raw = visual.get(key)
            if not isinstance(raw, dict):
                raise ContractError(f"visual audit is missing {key}")
            widths = raw.get("widths")
            if not isinstance(widths, dict):
                raise ContractError(f"visual audit {key} is missing widths")
            browser_errors = raw.get("browserErrors", [])
            if not isinstance(browser_errors, list):
                raise ContractError(f"visual audit {key}.browserErrors must be a list")
            visual_rows.append({
                "sample": "HCC1395",
                "surface": surface,
                "viewport": viewport,
                "evidence_status": visual.get("evidence_status", "UNKNOWN"),
                "viewport_width": integer(
                    raw.get("viewport", {}).get("width"), f"visual.{key}.viewport.width"
                ),
                "horizontal_overflow": "yes" if boolean(
                    widths.get("horizontalOverflow"), f"visual.{key}.horizontalOverflow"
                ) else "no",
                "browser_error_n": len(browser_errors),
                "screenshot": raw.get("screenshot", ""),
            })

    topology_details = scoped_rows(cohort)
    topology_rates: list[dict[str, Any]] = []
    for item in cohort:
        for metric, numerator, denominator, value in (
            ("Tree coverage", item["tree_n"], item["regions"], item["tree_coverage"]),
            (
                "Unique among tree",
                item["unique_tree_n"],
                item["tree_n"],
                item["unique_among_tree"],
            ),
        ):
            topology_rates.extend(scoped_rows([{
                "sample": item["sample"],
                "biological_id": item["biological_id"],
                "technical_replicate": item["technical_replicate"],
                "metric": metric,
                "numerator": numerator,
                "denominator": denominator,
                "rate": value,
            }]))

    downstream_samples = {row["sample"] for row in bundle_overview}
    scope_summary: list[dict[str, Any]] = []
    scope_groups = [(ALL_DATASETS_SCOPE, cohort)] + [
        (item["sample"], [item]) for item in cohort
    ]
    availability: list[dict[str, Any]] = []
    for sample_filter, group in scope_groups:
        dataset_count = len(group)
        biological_count = len({item["biological_id"] for item in group})
        integrity_pass_n = sum(bool(item["integrity_pass"]) for item in group)
        technical_gate_pass_n = sum(bool(item["technical_gate_pass"]) for item in group)
        all_families_complete_n = sum(bool(item["all_families_complete"]) for item in group)
        all_objective_certified_n = sum(bool(item["all_objective_certified"]) for item in group)
        biological_rows = [
            item for item in group if item["technical_replicate"] == "no"
        ]
        biological_macro_n = len(biological_rows)
        biological_macro_tree_coverage_median = (
            median(item["tree_coverage"] for item in biological_rows)
            if biological_rows else None
        )
        biological_macro_unique_median = (
            median(item["unique_among_tree"] for item in biological_rows)
            if biological_rows else None
        )
        downstream_available_n = sum(
            item["sample"] in downstream_samples for item in group
        )
        downstream_available_value = (
            downstream_available_n if downstream_available_n else None
        )
        scope_summary.append({
            "sample_filter": sample_filter,
            "dataset_count": dataset_count,
            "biological_count": biological_count,
            "integrity_pass_n": integrity_pass_n,
            "integrity_pass_rate": ratio(
                integrity_pass_n, dataset_count, f"{sample_filter}.integrity pass rate"
            ),
            "technical_gate_pass_n": technical_gate_pass_n,
            "technical_gate_total_n": dataset_count,
            "all_families_complete_n": all_families_complete_n,
            "all_families_complete_total_n": dataset_count,
            "all_objective_certified_n": all_objective_certified_n,
            "all_objective_certified_total_n": dataset_count,
            "biological_macro_n": biological_macro_n,
            "biological_macro_tree_coverage_median": biological_macro_tree_coverage_median,
            "biological_macro_unique_among_tree_median": biological_macro_unique_median,
            "downstream_available_count": downstream_available_value,
            "downstream_total_count": dataset_count,
        })
        if sample_filter == ALL_DATASETS_SCOPE:
            ism_available = downstream_available_value
            ism_status = "PARTIAL_HCC1395_ONLY"
            ism_note = (
                "Only the HCC1395 primary dataset has an audited bundle; v1 source is missing "
                "and v3 source is not content-addressed."
            )
            lineage_available = downstream_available_value
            lineage_status = "PARTIAL_HCC1395_ONLY_AB_GATE_FAIL"
            lineage_note = (
                "Only HCC1395 primary has lineage/LCA receipts, and its controlled A/B identity gate fails."
            )
        elif sample_filter == "HCC1395":
            ism_available = 1
            ism_status = "AVAILABLE_BLOCKED_OBSERVATION"
            ism_note = (
                "HCC1395 v1/v3 observations exist, but v1 source is missing and v3 source is not content-addressed."
            )
            lineage_available = 1
            lineage_status = "AVAILABLE_AB_GATE_FAIL"
            lineage_note = "HCC1395 lineage/LCA receipts exist, but the controlled A/B identity gate fails."
        else:
            ism_available = None
            ism_status = "ABSENT_NO_EQUIVALENT_BUNDLE"
            ism_note = (
                f"No sample-matched ISM bundle exists for {sample_filter}; HCC1395 primary data are not borrowed."
            )
            lineage_available = None
            lineage_status = "ABSENT_NO_EQUIVALENT_BUNDLE"
            lineage_note = (
                f"No sample-matched lineage/LCA bundle exists for {sample_filter}; HCC1395 primary data are not borrowed."
            )
        component_rows = [
            (
                "Topology + MLHP",
                dataset_count,
                dataset_count,
                "AVAILABLE",
                "Technical identity/hash/receipt gates only; not truth-set validation.",
            ),
            (
                "ISM drilldown",
                ism_available,
                dataset_count,
                ism_status,
                ism_note,
            ),
            (
                "Lineage / LCA",
                lineage_available,
                dataset_count,
                lineage_status,
                lineage_note,
            ),
            (
                "Truth benchmark",
                None,
                dataset_count,
                "UNAVAILABLE_NOT_MEASURED",
                "SEQC2 truth calls, HC BED, and som.py benchmark are not integrated.",
            ),
            (
                "Per-BAM mapping rate",
                None,
                dataset_count,
                "UNAVAILABLE_NOT_MEASURED",
                "No reviewed per-BAM mapping-rate extract is present in this snapshot.",
            ),
            (
                "Per-BAM depth",
                None,
                dataset_count,
                "UNAVAILABLE_NOT_MEASURED",
                "No reviewed depth or coverage distribution is present in this snapshot.",
            ),
            (
                "Read N50",
                None,
                dataset_count,
                "UNAVAILABLE_NOT_MEASURED",
                "No reviewed per-BAM read-length N50 is present in this snapshot.",
            ),
            (
                "HP/PS tagging rate",
                None,
                dataset_count,
                "UNAVAILABLE_NOT_MEASURED",
                "No reviewed HP/PS tagged-read numerator and denominator are present.",
            ),
            (
                "MM/ML completeness",
                None,
                dataset_count,
                "UNAVAILABLE_NOT_MEASURED",
                "No reviewed MM/ML tag completeness numerator and denominator are present.",
            ),
            (
                "KDE provenance",
                None,
                dataset_count,
                "UNAVAILABLE_NOT_MEASURED",
                "The inputs do not declare whether methylation values are KDE-corrected or identify a KDE receipt.",
            ),
        ]
        for component, available_n, expected_n, status, note in component_rows:
            coverage_rate = (
                None
                if available_n is None
                else ratio(available_n, expected_n, f"{sample_filter}.{component}")
            )
            availability.append({
                "sample_filter": sample_filter,
                "component": component,
                "available_datasets": available_n,
                "expected_datasets": expected_n,
                "coverage": coverage_rate,
                "status": status,
                "note": note,
            })

    generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    sources = [
        make_snapshot_source(
            "q_scope_summary",
            "Cohort scope summary",
            "scope_summary",
            [paths["cohort"], paths["bundle"]],
            generated_at,
            "Counts dataset and biological-sample scope, technical integrity, and audited downstream availability.",
            [
                "dataset_count = topology dataset rows in the selected view.",
                "biological_count = distinct biological_id values; HCC1395_DORADO remains a technical replicate.",
                "integrity_pass_rate = datasets passing hash, sample identity, and technical receipt gates / datasets in view.",
                "technical_gate_pass_n = datasets with receipt_all_pass=true.",
                "all_families_complete_n = datasets with receipt_all_mutation_bearing_families_complete=true.",
                "all_objective_certified_n = datasets with receipt_all_units_objective_certified=true.",
                "biological macro medians exclude rows marked technical_replicate=true.",
                "downstream_available_count = datasets with an audited drilldown bundle; availability does not imply scientific validation.",
            ],
            "sample_filter",
        ),
        make_snapshot_source(
            "q_availability",
            "Cross-sample data availability matrix",
            "availability",
            [paths["cohort"], paths["bundle"], paths["audit"]],
            generated_at,
            "Shows the measured availability boundary for topology, ISM, lineage/LCA, and truth benchmarking.",
            ["coverage = available datasets / expected topology datasets in the selected view."],
            "sample_filter, component",
        ),
        make_snapshot_source(
            "q_topology_rates",
            "Per-dataset topology comparison",
            "topology_rates",
            [paths["cohort"]],
            generated_at,
            "Long-form per-dataset topology coverage and unique-among-tree rates with exact numerators and denominators.",
            [
                "tree coverage = tree_n / region_n.",
                "unique among tree = unique_tree_n / tree_n.",
            ],
            "sample_filter, sample, metric",
        ),
        make_snapshot_source(
            "q_coverage_by_k",
            "HCC1395 ISM linkage by active-k bin",
            "coverage_by_k",
            [paths["coverage"]],
            generated_at,
            "HCC1395 v1/v3 ISM linkage by topology active-k bin; rows are not a multi-sample comparison.",
            ["coverage = linked distinct genomic sSNV / topology distinct genomic sSNV."],
            "sample_filter, k_order, bundle",
        ),
        make_snapshot_source(
            "q_axis_rates",
            "HCC1395 topology-linked methylation-axis rates",
            "axis_rates",
            [paths["axis"]],
            generated_at,
            "Exploratory raw-p and within-family BH rates for topology-linked HCC1395 rows.",
            [
                "raw rate = rows with raw p <= 0.05 / tested rows.",
                "BH rate = rows with recomputed within-axis-and-scope BH q <= 0.05 / tested rows.",
                "Cluster is circular double-dipping and is invalid as evidence.",
            ],
            "sample_filter, axis, series",
            ["analysis_scope = topology_linked_distinct_ssnv"],
        ),
        make_snapshot_source(
            "q_bundle_overview",
            "HCC1395 v1/v3 bundle audit",
            "bundle_overview",
            [paths["bundle"]],
            generated_at,
            "Exact HCC1395 v1/v3 bundle metrics, selfcheck status, and asset-accounting gaps.",
            [
                "ISM linkage = methyl_topology_linked / distinct topology sSNV.",
                "All rate fields are stored as decimal fractions in this artifact.",
            ],
            "sample_filter, bundle",
        ),
        make_snapshot_source(
            "q_topology_details",
            "Seven-dataset topology detail",
            "topology_details",
            [paths["cohort"]],
            generated_at,
            "Exact topology denominators and technical provenance gates for seven datasets representing six biological samples.",
            [
                "Technical PASS does not prove a unique biological clone tree.",
                "HCC1395_DORADO is a technical replicate of biological_id HCC1395.",
            ],
            "sample_filter, sample",
        ),
        make_snapshot_source(
            "q_axis_details",
            "HCC1395 methylation-axis detail",
            "axis_details",
            [paths["axis"]],
            generated_at,
            "Full v1/v3 axis metrics for both all-summary and topology-linked scopes.",
            [
                "BH correction is within one bundle x axis x declared scope, not cohort-wide.",
                "Effect units differ by axis and must not be compared across incompatible units.",
            ],
            "sample_filter, bundle, analysis_scope, axis",
        ),
        make_snapshot_source(
            "q_lca_summary",
            "HCC1395 LCA A/B identity audit",
            "lca_summary",
            [paths["audit"]],
            generated_at,
            "Summarizes the 22-chromosome HCC1395 LCA comparison and its failed causal A/B gate.",
            [
                "same_in_bam_rate = chromosomes with identical in_bam / shared chromosomes.",
                "same_threads_rate = chromosomes with identical thread count / shared chromosomes.",
            ],
            "sample_filter, sample",
        ),
        make_snapshot_source(
            "q_visual_qa",
            "Legacy/current browser visual QA",
            "visual_qa",
            [paths["visual"]],
            generated_at,
            "Desktop/mobile overflow and browser-error receipts for the legacy and direct-generated HCC1395 interfaces.",
            ["browser_error_n = captured browserErrors array length for one surface and viewport."],
            "sample_filter, surface, viewport",
        ),
        {
            "id": "builder_source",
            "label": "Canonical dashboard artifact builder",
            "path": source_path(Path(__file__)),
        },
    ]

    cards = [
        {
            "id": "c_dataset_count",
            "dataset": "scope_summary",
            "sourceId": "q_scope_summary",
            "description": "目前篩選範圍內的 topology 技術資料集數；資料集不等同生物樣本。",
            "metrics": [
                {"label": "目前資料集", "field": "dataset_count", "format": "number"}
            ],
        },
        {
            "id": "c_biological_count",
            "dataset": "scope_summary",
            "sourceId": "q_scope_summary",
            "description": "不同 biological_id 的數量；HCC1395_DORADO 不重複計為另一個生物樣本。",
            "metrics": [
                {"label": "生物樣本", "field": "biological_count", "format": "number"}
            ],
        },
        {
            "id": "c_technical_gate",
            "dataset": "scope_summary",
            "sourceId": "q_scope_summary",
            "description": "receipt_all_pass=true 的資料集；技術 PASS 不等於生物真值驗證。",
            "metrics": [
                {
                    "label": "技術檢核 PASS",
                    "field": "technical_gate_pass_n",
                    "format": "number",
                },
                {"label": "目前資料集", "field": "technical_gate_total_n", "format": "number"},
            ],
        },
        {
            "id": "c_all_families",
            "dataset": "scope_summary",
            "sourceId": "q_scope_summary",
            "description": "所有帶突變 family 均完整的資料集；0 是實測 gate 結果，不是缺值。",
            "metrics": [
                {
                    "label": "全 family 完整",
                    "field": "all_families_complete_n",
                    "format": "number",
                },
                {"label": "目前資料集", "field": "all_families_complete_total_n", "format": "number"},
            ],
        },
        {
            "id": "c_all_objective",
            "dataset": "scope_summary",
            "sourceId": "q_scope_summary",
            "description": "所有分析單位均具 objective certification 的資料集；0 是實測 gate 結果。",
            "metrics": [
                {
                    "label": "全單位 objective-certified",
                    "field": "all_objective_certified_n",
                    "format": "number",
                },
                {"label": "目前資料集", "field": "all_objective_certified_total_n", "format": "number"},
            ],
        },
        {
            "id": "c_macro_tree",
            "dataset": "scope_summary",
            "sourceId": "q_scope_summary",
            "description": "排除技術重複後，各生物樣本 tree coverage 的不加權中位數；不 pooled loci。",
            "metrics": [
                {
                    "label": "生物樣本 macro tree coverage",
                    "field": "biological_macro_tree_coverage_median",
                    "format": "percent",
                },
                {"label": "生物樣本 n", "field": "biological_macro_n", "format": "number"},
            ],
        },
        {
            "id": "c_macro_unique",
            "dataset": "scope_summary",
            "sourceId": "q_scope_summary",
            "description": "排除技術重複後，各生物樣本 unique/tree 比例的不加權中位數。",
            "metrics": [
                {
                    "label": "生物樣本 macro unique / tree",
                    "field": "biological_macro_unique_among_tree_median",
                    "format": "percent",
                },
                {"label": "生物樣本 n", "field": "biological_macro_n", "format": "number"},
            ],
        },
        {
            "id": "c_downstream_count",
            "dataset": "scope_summary",
            "sourceId": "q_scope_summary",
            "description": "具有已稽核 drilldown bundle 的資料集；目前沒有可比較的多樣本 ISM/lineage 矩陣。",
            "metrics": [
                {
                    "label": "已稽核下游資料集",
                    "field": "downstream_available_count",
                    "format": "number",
                },
                {
                    "label": "目前資料集",
                    "field": "downstream_total_count",
                    "format": "number",
                },
            ],
        },
    ]

    charts = [
        {
            "id": "ch_topology_rates",
            "title": "Topology coverage and uniqueness by dataset",
            "subtitle": "Seven technical datasets represent six biological samples; exact numerators and denominators remain in the chart data.",
            "type": "bar",
            "intent": "comparison",
            "dataset": "topology_rates",
            "sourceId": "q_topology_rates",
            "valueFormat": "percent",
            "encodings": {
                "x": {"field": "sample", "type": "nominal", "label": "Dataset"},
                "y": {
                    "field": "rate",
                    "type": "quantitative",
                    "label": "Rate",
                    "format": "percent",
                },
                "color": {"field": "metric", "type": "nominal", "label": "Metric"},
                "tooltip": [
                    {"field": "biological_id", "type": "nominal", "label": "Biological sample"},
                    {"field": "numerator", "type": "quantitative", "label": "Numerator"},
                    {"field": "denominator", "type": "quantitative", "label": "Denominator"},
                    {"field": "technical_replicate", "type": "nominal", "label": "Technical replicate"},
                ],
            },
            "maxRows": 40,
        },
        {
            "id": "ch_coverage_by_k",
            "title": "HCC1395 ISM linkage by active-k bin",
            "subtitle": "v1 and v3 use different windows, so this is a descriptive availability diagnostic rather than a controlled A/B comparison.",
            "type": "line",
            "intent": "trend",
            "dataset": "coverage_by_k",
            "sourceId": "q_coverage_by_k",
            "valueFormat": "percent",
            "encodings": {
                "x": {"field": "k_bin", "type": "ordinal", "label": "Active k"},
                "y": {
                    "field": "coverage",
                    "type": "quantitative",
                    "label": "ISM linkage",
                    "format": "percent",
                },
                "color": {"field": "bundle", "type": "nominal", "label": "Bundle"},
                "tooltip": [
                    {"field": "linked", "type": "quantitative", "label": "Linked"},
                    {"field": "denominator", "type": "quantitative", "label": "Denominator"},
                ],
            },
            "maxRows": 40,
        },
        {
            "id": "ch_axis_rates",
            "title": "HCC1395 topology-linked methylation-axis rates",
            "subtitle": "Raw p and within-family BH rates are exploratory; cluster is circular double-dipping and invalid as evidence.",
            "type": "bar",
            "intent": "comparison",
            "dataset": "axis_rates",
            "sourceId": "q_axis_rates",
            "valueFormat": "percent",
            "encodings": {
                "x": {"field": "axis", "type": "nominal", "label": "Methylation axis"},
                "y": {
                    "field": "rate",
                    "type": "quantitative",
                    "label": "Significant-row rate",
                    "format": "percent",
                },
                "color": {"field": "series", "type": "nominal", "label": "Bundle · threshold"},
                "tooltip": [
                    {"field": "significant_n", "type": "quantitative", "label": "Significant rows"},
                    {"field": "tested_n", "type": "quantitative", "label": "Tested rows"},
                    {"field": "interpretation_gate", "type": "nominal", "label": "Interpretation gate"},
                ],
            },
            "maxRows": 50,
        },
    ]

    tables = [
        {
            "id": "t_availability",
            "title": "Cross-sample data availability",
            "subtitle": "Unavailable values are null (blank), never zero; this table separates measured coverage from data gaps.",
            "dataset": "availability",
            "sourceId": "q_availability",
            "defaultSort": {"field": "component", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "component", "label": "Component", "type": "text"},
                {"field": "available_datasets", "label": "Available", "format": "number"},
                {"field": "expected_datasets", "label": "Expected", "format": "number"},
                {"field": "coverage", "label": "Coverage", "format": "percent"},
                {"field": "status", "label": "Status", "type": "text"},
                {"field": "note", "label": "Boundary", "type": "text"},
            ],
        },
        {
            "id": "t_bundle_overview",
            "title": "HCC1395 v1/v3 bundle diagnostics",
            "subtitle": "Both bundles remain blocked observations; rates are decimal fractions in the snapshot.",
            "dataset": "bundle_overview",
            "sourceId": "q_bundle_overview",
            "defaultSort": {"field": "bundle", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "bundle", "label": "Bundle", "type": "text"},
                {"field": "validation_status", "label": "Validation", "type": "text"},
                {"field": "regions", "label": "Regions", "format": "number"},
                {"field": "tree_coverage", "label": "Tree coverage", "format": "percent"},
                {"field": "unique_among_tree", "label": "Unique / tree", "format": "percent"},
                {"field": "ism_linked", "label": "ISM linked", "format": "number"},
                {"field": "ism_linkage", "label": "ISM linkage", "format": "percent"},
                {"field": "selfcheck_fail", "label": "FAIL", "format": "number"},
                {"field": "selfcheck_skip", "label": "SKIP", "format": "number"},
                {"field": "panel_actual_bytes", "label": "Panel bytes", "format": "number"},
                {"field": "panel_unreported_bytes", "label": "Unreported bytes", "format": "number"},
            ],
        },
        {
            "id": "t_topology_details",
            "title": "Seven-dataset topology detail",
            "subtitle": "Per-dataset denominators are retained; HCC1395_DORADO is a technical replicate.",
            "dataset": "topology_details",
            "sourceId": "q_topology_details",
            "defaultSort": {"field": "tree_coverage", "direction": "desc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "sample", "label": "Dataset", "type": "text"},
                {"field": "biological_id", "label": "Biological sample", "type": "text"},
                {"field": "technical_replicate", "label": "Technical replicate", "type": "text"},
                {"field": "regions", "label": "Regions", "format": "number"},
                {"field": "distinct_ssnv", "label": "Distinct sSNV", "format": "number"},
                {"field": "tree_n", "label": "Trees", "format": "number"},
                {"field": "tree_coverage", "label": "Tree coverage", "format": "percent"},
                {"field": "unique_tree_n", "label": "Unique trees", "format": "number"},
                {"field": "unique_among_tree", "label": "Unique / tree", "format": "percent"},
                {"field": "hash_match", "label": "Hash", "type": "text"},
                {"field": "sample_identity_match", "label": "Identity", "type": "text"},
                {"field": "receipt_technical_pass", "label": "Receipt", "type": "text"},
                {"field": "all_families_complete_gate", "label": "All families", "type": "text"},
                {"field": "all_objective_certified_gate", "label": "All objective", "type": "text"},
            ],
        },
        {
            "id": "t_axis_details",
            "title": "HCC1395 methylation-axis detail",
            "subtitle": "Effect units differ by axis; the cluster axis is not valid evidence.",
            "dataset": "axis_details",
            "sourceId": "q_axis_details",
            "defaultSort": {"field": "tested_n", "direction": "desc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "bundle", "label": "Bundle", "type": "text"},
                {"field": "analysis_scope", "label": "Analysis scope", "type": "text"},
                {"field": "axis", "label": "Axis", "type": "text"},
                {"field": "tested_n", "label": "Tested", "format": "number"},
                {"field": "raw_sig_rate", "label": "Raw p≤.05", "format": "percent"},
                {"field": "bh_sig_rate", "label": "BH q≤.05", "format": "percent"},
                {"field": "effect_median", "label": "Effect median", "format": "number"},
                {"field": "effect_unit", "label": "Effect unit", "type": "text"},
                {"field": "interpretation_gate", "label": "Interpretation gate", "type": "text"},
            ],
        },
        {
            "id": "t_lca_summary",
            "title": "HCC1395 lineage/LCA controlled-comparison gate",
            "subtitle": "Input identity is not held constant, so the observed delta is not a causal LCA effect.",
            "dataset": "lca_summary",
            "sourceId": "q_lca_summary",
            "defaultSort": {"field": "sample", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "sample", "label": "Sample", "type": "text"},
                {"field": "shared_chromosomes", "label": "Shared chr", "format": "number"},
                {"field": "same_in_bam_n", "label": "Same BAM", "format": "number"},
                {"field": "same_in_bam_rate", "label": "Same BAM rate", "format": "percent"},
                {"field": "same_threads_n", "label": "Same threads", "format": "number"},
                {"field": "same_threads_rate", "label": "Same threads rate", "format": "percent"},
                {"field": "pre_lv_written", "label": "Pre lv", "format": "number"},
                {"field": "post_lv_written", "label": "Post lv", "format": "number"},
                {"field": "net_new_lv_written", "label": "Net new lv", "format": "number"},
                {"field": "causal_ab_gate", "label": "Causal A/B gate", "type": "text"},
            ],
        },
        {
            "id": "t_visual_qa",
            "title": "Legacy/current browser QA receipt",
            "subtitle": "Visual QA verifies interface behavior; it does not raise the scientific claim ceiling.",
            "dataset": "visual_qa",
            "sourceId": "q_visual_qa",
            "defaultSort": {"field": "surface", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "surface", "label": "Surface", "type": "text"},
                {"field": "viewport", "label": "Viewport", "type": "text"},
                {"field": "viewport_width", "label": "Width", "format": "number"},
                {"field": "horizontal_overflow", "label": "Overflow", "type": "text"},
                {"field": "browser_error_n", "label": "Browser errors", "format": "number"},
                {"field": "evidence_status", "label": "Evidence status", "type": "text"},
                {"field": "screenshot", "label": "Screenshot receipt", "type": "text"},
            ],
        },
    ]

    used_datasets = sorted({
        *(card["dataset"] for card in cards),
        *(chart["dataset"] for chart in charts),
        *(table["dataset"] for table in tables),
    })
    filter_targets = [
        {"dataset": dataset, "field": "sample_filter"}
        for dataset in used_datasets
        if dataset != "scope_summary"
    ]

    cohort_scope = scope_summary[0]
    blocks = [
        {
            "id": "b_summary_header",
            "type": "markdown",
            "body": (
                "## Summary — technical cohort coverage is available; downstream evidence is partial\n\n"
                f"At the default All scope, the technical gate is {cohort_scope['technical_gate_pass_n']}/"
                f"{cohort_scope['technical_gate_total_n']}, while all-families-complete is "
                f"{cohort_scope['all_families_complete_n']}/{cohort_scope['all_families_complete_total_n']} "
                f"and all-objective-certified is {cohort_scope['all_objective_certified_n']}/"
                f"{cohort_scope['all_objective_certified_total_n']}. Biological macro medians exclude "
                f"the technical replicate (n={cohort_scope['biological_macro_n']}). Use **Sample filter** "
                "to isolate one technical dataset. Missing ISM, lineage, BAM QC, KDE provenance, or "
                "truth-benchmark values are null/unavailable, never zero signal. The claim ceiling remains "
                "descriptive observation and internal data-product QA, not truth-set validation."
            ),
        },
        {
            "id": "b_summary_metrics",
            "type": "metric-strip",
            "cardIds": [
                "c_dataset_count", "c_biological_count", "c_technical_gate",
                "c_all_families", "c_all_objective", "c_macro_tree",
                "c_macro_unique", "c_downstream_count",
            ],
        },
        {"id": "b_availability", "type": "table", "tableId": "t_availability", "layout": "full"},
        {"id": "b_topology_chart", "type": "chart", "chartId": "ch_topology_rates", "layout": "full"},
        {
            "id": "b_diagnostics_header",
            "type": "markdown",
            "body": "## Diagnostics — HCC1395 downstream selection and multiplicity\n\nThese diagnostics retain linked/denominator counts and distinguish raw-p discovery from within-family BH correction. v1/v3 window differences and the circular cluster axis remain explicit gates.",
        },
        {"id": "b_coverage_chart", "type": "chart", "chartId": "ch_coverage_by_k", "layout": "half"},
        {"id": "b_axis_chart", "type": "chart", "chartId": "ch_axis_rates", "layout": "half"},
        {"id": "b_bundle_table", "type": "table", "tableId": "t_bundle_overview", "layout": "full"},
        {
            "id": "b_details_header",
            "type": "markdown",
            "body": "## Details — exact denominators, gates, and QA receipts\n\nThe tables below are audit surfaces: they preserve technical replicates, analysis scope, effect units, causal gates, and browser receipts rather than collapsing them into pooled cohort claims.",
        },
        {"id": "b_topology_table", "type": "table", "tableId": "t_topology_details", "layout": "full"},
        {"id": "b_axis_table", "type": "table", "tableId": "t_axis_details", "layout": "full"},
        {"id": "b_lca_table", "type": "table", "tableId": "t_lca_summary", "layout": "full"},
        {"id": "b_visual_table", "type": "table", "tableId": "t_visual_qa", "layout": "full"},
        {
            "id": "b_method",
            "type": "markdown",
            "body": "## Interpretation guardrails\n\n- Cohort topology remains one row per technical dataset; no pooled-locus inference is performed.\n- HCC1395_DORADO is retained for reproducibility and is not a seventh biological sample.\n- All fields rendered as percent are stored on a 0–1 scale.\n- Missing per-BAM mapping rate, depth, read N50, HP/PS tagging rate, MM/ML completeness, KDE provenance, and truth benchmarking remain null with explicit access status.\n- ISM availability is selection-sensitive across active-k; untested is not nonsignificant.\n- Lineage/LCA output is descriptive until input identity is held constant.\n- Browser screenshots validate presentation only, not biological truth.",
        },
    ]

    artifact = {
        "surface": "dashboard",
        "manifest": {
            "version": 1,
            "surface": "dashboard",
            "title": "Canonical multi-sample drilldown readiness",
            "description": "Seven topology datasets / six biological samples with HCC1395-only downstream diagnostics, explicit missing-access gates, and source-backed details.",
            "generatedAt": generated_at,
            "filters": [
                {
                    "id": "sample_filter",
                    "label": "Sample filter",
                    "dataset": "scope_summary",
                    "field": "sample_filter",
                    "defaultValue": ALL_DATASETS_SCOPE,
                    "includeAll": False,
                    "targets": filter_targets,
                }
            ],
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": sources,
            "blocks": blocks,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "partial",
            "datasets": {
                "scope_summary": scope_summary,
                "availability": availability,
                "topology_rates": topology_rates,
                "coverage_by_k": hcc_scoped_rows(coverage),
                "axis_rates": hcc_scoped_rows(axis_rates_base),
                "bundle_overview": hcc_scoped_rows(bundle_overview),
                "topology_details": topology_details,
                "axis_details": hcc_scoped_rows(axis_details),
                "lca_summary": hcc_scoped_rows(lca_summary),
                "visual_qa": hcc_scoped_rows(visual_rows),
            },
            "accessIssues": [
                {
                    "id": "multisample_ism_unavailable",
                    "scope": "G4 multi-sample ISM",
                    "sourceId": "q_availability",
                    "message": "Comparable ISM is available only as HCC1395 v1/v3 observations. The v1 source root is missing and the v3 source directory is not content-addressed; no 7-dataset/6-biological-sample ISM matrix exists.",
                },
                {
                    "id": "multisample_lineage_unavailable",
                    "scope": "G4 multi-sample lineage",
                    "sourceId": "q_lca_summary",
                    "message": "Lineage/LCA receipts cover HCC1395 only, and the 22-chromosome comparison fails the controlled A/B gate because input BAM and thread identities are not held constant.",
                },
                {
                    "id": "biological_replicate_policy_unset",
                    "scope": "cohort aggregation",
                    "sourceId": "q_topology_details",
                    "message": "HCC1395_DORADO is a technical replicate and no frozen biological-replicate aggregation policy is declared; the dashboard keeps per-dataset rows and does not pool loci or report a six-sample macro estimate.",
                },
                {
                    "id": "truth_benchmark_unavailable",
                    "scope": "scientific validation",
                    "sourceId": "q_availability",
                    "message": "SEQC2 truth calls, high-confidence BED, and som.py metrics are absent, so technical PASS cannot be promoted to truth-set validation.",
                },
                {
                    "id": "per_bam_qc_metrics_unavailable",
                    "scope": "per-BAM data quality",
                    "sourceId": "q_availability",
                    "message": "This snapshot has no reviewed per-BAM mapping rate, depth, read N50, HP/PS tagging rate, or MM/ML completeness. These values are null/unavailable and are not represented as zero.",
                },
                {
                    "id": "kde_provenance_unavailable",
                    "scope": "methylation provenance",
                    "sourceId": "q_availability",
                    "message": "The reviewed inputs do not declare KDE-corrected status or identify a KDE provenance receipt. KDE status remains unavailable rather than false or zero.",
                },
            ],
        },
        "sources": sources,
    }
    validate_artifact_contract(artifact)
    return artifact


def validate_artifact_contract(artifact: dict[str, Any]) -> None:
    """Run bounded checks that complement the shared validate_artifact tool."""
    manifest = artifact["manifest"]
    snapshot = artifact["snapshot"]
    if artifact.get("surface") != "dashboard" or manifest.get("surface") != "dashboard":
        raise ContractError("artifact surface must be dashboard")
    if snapshot.get("status") != "partial":
        raise ContractError("snapshot status must remain partial")

    source_ids = {source.get("id") for source in manifest.get("sources", [])}
    if None in source_ids or len(source_ids) != len(manifest.get("sources", [])):
        raise ContractError("manifest source IDs must be unique and non-empty")
    for kind in ("cards", "charts", "tables"):
        for item in manifest.get(kind, []):
            if item.get("sourceId") not in source_ids:
                raise ContractError(
                    f"{kind[:-1]} {item.get('id')} has unresolved sourceId {item.get('sourceId')!r}"
                )
            if item.get("dataset") not in snapshot["datasets"]:
                raise ContractError(
                    f"{kind[:-1]} {item.get('id')} references missing dataset {item.get('dataset')!r}"
                )

    card_ids = {item["id"] for item in manifest["cards"]}
    chart_ids = {item["id"] for item in manifest["charts"]}
    table_ids = {item["id"] for item in manifest["tables"]}
    for block in manifest["blocks"]:
        if block["type"] == "metric-strip" and not set(block["cardIds"]).issubset(card_ids):
            raise ContractError(f"block {block['id']} has unresolved cardIds")
        if block["type"] == "chart" and block.get("chartId") not in chart_ids:
            raise ContractError(f"block {block['id']} has unresolved chartId")
        if block["type"] == "table" and block.get("tableId") not in table_ids:
            raise ContractError(f"block {block['id']} has unresolved tableId")

    percent_fields: dict[str, set[str]] = {}
    for card in manifest["cards"]:
        for metric in card["metrics"]:
            if metric.get("format") == "percent":
                percent_fields.setdefault(card["dataset"], set()).add(metric["field"])
    for chart in manifest["charts"]:
        for encoding in chart.get("encodings", {}).values():
            encodings = encoding if isinstance(encoding, list) else [encoding]
            for item in encodings:
                if isinstance(item, dict) and item.get("format") == "percent":
                    percent_fields.setdefault(chart["dataset"], set()).add(item["field"])
    for table in manifest["tables"]:
        for column in table["columns"]:
            if column.get("format") == "percent":
                percent_fields.setdefault(table["dataset"], set()).add(column["field"])
        sort_field = table["defaultSort"]["field"]
        if sort_field not in {column["field"] for column in table["columns"]}:
            raise ContractError(f"table {table['id']} defaultSort is not a declared column")

    for dataset, fields in percent_fields.items():
        for row_index, row in enumerate(snapshot["datasets"][dataset]):
            for field in fields:
                value = row.get(field)
                if value is None:
                    continue
                if not isinstance(value, (int, float)) or isinstance(value, bool) or not 0 <= value <= 1:
                    raise ContractError(
                        f"{dataset}[{row_index}].{field} must be a 0-1 percent value, got {value!r}"
                    )

    asset_datasets = {
        *(item["dataset"] for item in manifest["cards"]),
        *(item["dataset"] for item in manifest["charts"]),
        *(item["dataset"] for item in manifest["tables"]),
    }
    dashboard_filter = manifest["filters"][0]
    filter_datasets = {dashboard_filter["dataset"], *(
        target["dataset"] for target in dashboard_filter.get("targets", [])
    )}
    if filter_datasets != asset_datasets:
        raise ContractError(
            "dataset filter must target every card/chart/table dataset: "
            f"missing={sorted(asset_datasets - filter_datasets)}"
        )
    for dataset in asset_datasets:
        if not all("sample_filter" in row for row in snapshot["datasets"][dataset]):
            raise ContractError(f"dataset {dataset} has rows without sample_filter")

    scope_all = next(
        row for row in snapshot["datasets"]["scope_summary"]
        if row["sample_filter"] == ALL_DATASETS_SCOPE
    )
    if (
        scope_all["dataset_count"] != EXPECTED_DATASET_N
        or scope_all["biological_count"] != EXPECTED_BIOLOGICAL_N
    ):
        raise ContractError("artifact headline must preserve 7 datasets / 6 biological samples")
    expected_gate_counts = {
        "technical_gate_pass_n": 7,
        "all_families_complete_n": 0,
        "all_objective_certified_n": 0,
        "biological_macro_n": 6,
    }
    for field, expected in expected_gate_counts.items():
        if scope_all[field] != expected:
            raise ContractError(
                f"All-scope {field} must be {expected}, got {scope_all[field]}"
            )
    biological_rows = [
        row for row in snapshot["datasets"]["topology_details"]
        if row["sample_filter"] == ALL_DATASETS_SCOPE
        and row["technical_replicate"] == "no"
    ]
    assert_close(
        scope_all["biological_macro_tree_coverage_median"],
        median(row["tree_coverage"] for row in biological_rows),
        "biological macro tree coverage median",
    )
    assert_close(
        scope_all["biological_macro_unique_among_tree_median"],
        median(row["unique_among_tree"] for row in biological_rows),
        "biological macro unique-among-tree median",
    )

    for sample in sorted(
        row["sample"] for row in biological_rows
        if row["sample"] != "HCC1395"
    ) + ["HCC1395_DORADO"]:
        downstream_rows = [
            row for row in snapshot["datasets"]["availability"]
            if row["sample_filter"] == sample
            and row["component"] in {"ISM drilldown", "Lineage / LCA"}
        ]
        if len(downstream_rows) != 2 or any(
            row["status"] != "ABSENT_NO_EQUIVALENT_BUNDLE"
            or row["available_datasets"] is not None
            or row["coverage"] is not None
            for row in downstream_rows
        ):
            raise ContractError(
                f"{sample} must not borrow HCC1395 ISM/lineage data"
            )

    issue_ids = {issue["id"] for issue in snapshot.get("accessIssues", [])}
    required_issues = {"multisample_ism_unavailable", "multisample_lineage_unavailable"}
    if not required_issues.issubset(issue_ids):
        raise ContractError("partial snapshot must explicitly declare missing ISM and lineage access")


def main() -> int:
    args = parse_args()
    artifact = build_artifact(args.results_dir.resolve())
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    summary = {
        "output": source_path(output),
        "surface": artifact["surface"],
        "status": artifact["snapshot"]["status"],
        "topology_datasets": EXPECTED_DATASET_N,
        "biological_samples": EXPECTED_BIOLOGICAL_N,
        "snapshot_datasets": {
            key: len(rows) for key, rows in artifact["snapshot"]["datasets"].items()
        },
        "cards": len(artifact["manifest"]["cards"]),
        "charts": len(artifact["manifest"]["charts"]),
        "tables": len(artifact["manifest"]["tables"]),
        "access_issues": [
            issue["id"] for issue in artifact["snapshot"]["accessIssues"]
        ],
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
