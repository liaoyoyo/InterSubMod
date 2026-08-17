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
import html
import importlib.util
import json
import math
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parents[1]
DEFAULT_RESULTS = PROJECT_ROOT / "results"
DEFAULT_OUTPUT = PROJECT_ROOT / "multi_bam_dashboard_artifact.json"
DEFAULT_BAM_MANIFEST = PROJECT_ROOT / "multi_bam_input_manifest.json"
DEFAULT_BAM_MANIFEST_RECEIPT = DEFAULT_RESULTS / "multi_bam_input_manifest_validation.json"
DEFAULT_BAM_MANIFEST_SCHEMA = PROJECT_ROOT / "multi_bam_input_manifest.schema.json"
BAM_MANIFEST_BUILDER = PROJECT_ROOT / "scripts" / "build_multi_bam_input_manifest.py"
DEFAULT_SOURCE_SCHEMA = (
    REPO_ROOT
    / "docs"
    / "methodology"
    / "_assets"
    / "20260627_subclone_4axis_teaching"
    / "schemas"
    / "layered_input_manifest_v3.schema.json"
)
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


def validate_bam_manifest_contract(manifest: dict[str, Any]) -> None:
    spec = importlib.util.spec_from_file_location(
        "intersubmod_multi_bam_manifest_builder", BAM_MANIFEST_BUILDER
    )
    if spec is None or spec.loader is None:
        raise ContractError(f"cannot import BAM manifest validator: {BAM_MANIFEST_BUILDER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    try:
        module.validate_output(manifest)
    except module.ContractError as exc:
        raise ContractError(f"multi-BAM manifest contract failed: {exc}") from exc


def atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=path.parent
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except BaseException:
        if temporary_path.exists():
            temporary_path.rename(temporary_path.with_suffix(".failed"))
        raise


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
    parser.add_argument(
        "--bam-manifest",
        type=Path,
        default=DEFAULT_BAM_MANIFEST,
        help=f"Verified multi-BAM input manifest (default: {DEFAULT_BAM_MANIFEST})",
    )
    parser.add_argument(
        "--bam-manifest-receipt",
        type=Path,
        default=DEFAULT_BAM_MANIFEST_RECEIPT,
        help=(
            "Validation receipt for --bam-manifest "
            f"(default: {DEFAULT_BAM_MANIFEST_RECEIPT})"
        ),
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


def percentile_cont(values: Iterable[float], quantile: float, label: str) -> float:
    """Return SQL PERCENTILE_CONT-style linear interpolation."""
    ordered = sorted(values)
    if not ordered:
        raise ContractError(f"{label} requires at least one value")
    if not 0.0 <= quantile <= 1.0:
        raise ContractError(f"{label} quantile must be between 0 and 1")
    position = (len(ordered) - 1) * quantile
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] + (ordered[upper] - ordered[lower]) * fraction


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


def build_artifact(
    results_dir: Path,
    bam_manifest_path: Path,
    bam_manifest_receipt_path: Path,
) -> dict[str, Any]:
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
    bam_manifest = read_json(bam_manifest_path)
    bam_manifest_receipt = read_json(bam_manifest_receipt_path)

    if audit.get("schema_name") != "intersubmod.drilldown.audit.metrics":
        raise ContractError("metrics_audit.json has an unexpected schema_name")
    if visual.get("schema_name") != "intersubmod.legacy_browser_visual_audit":
        raise ContractError("legacy visual audit has an unexpected schema_name")
    if (
        bam_manifest.get("schema_name")
        != "intersubmod.multi_bam_dashboard_input_manifest"
        or bam_manifest.get("schema_version") != "1.1.0"
        or bam_manifest.get("dataset_count") != EXPECTED_DATASET_N
        or bam_manifest.get("biological_sample_count") != EXPECTED_BIOLOGICAL_N
    ):
        raise ContractError("multi-BAM input manifest schema or scope mismatch")
    validate_bam_manifest_contract(bam_manifest)
    manifest_summary = bam_manifest.get("verification_summary")
    receipt_summary = bam_manifest_receipt.get("verification_summary")
    if (
        bam_manifest_receipt.get("schema_name")
        != "intersubmod.multi_bam_dashboard_input_manifest_validation"
        or bam_manifest_receipt.get("schema_version") != "1.0.0"
        or bam_manifest_receipt.get("status")
        not in {"PASS_BOUNDED", "PASS_BOUNDED_WITH_METADATA_DRIFT"}
        or bam_manifest_receipt.get("status") != manifest_summary.get("status")
        or receipt_summary != manifest_summary
        or bam_manifest_receipt.get("claim_ceiling") != bam_manifest.get("claim_ceiling")
        or bam_manifest_receipt.get("output", {}).get("sha256") != sha256(bam_manifest_path)
        or bam_manifest_receipt.get("schema", {}).get("path")
        != source_path(DEFAULT_BAM_MANIFEST_SCHEMA)
        or bam_manifest_receipt.get("schema", {}).get("sha256")
        != sha256(DEFAULT_BAM_MANIFEST_SCHEMA)
        or bam_manifest_receipt.get("source_schema", {}).get("path")
        != source_path(DEFAULT_SOURCE_SCHEMA)
        or bam_manifest_receipt.get("source_schema", {}).get("sha256")
        != sha256(DEFAULT_SOURCE_SCHEMA)
        or bam_manifest_receipt.get("builder", {}).get("path")
        != source_path(BAM_MANIFEST_BUILDER)
        or bam_manifest_receipt.get("builder", {}).get("sha256")
        != sha256(BAM_MANIFEST_BUILDER)
        or bam_manifest_receipt.get("inputs", {}).get("source_manifest", {}).get("path")
        != bam_manifest.get("source_manifest", {}).get("path")
        or bam_manifest_receipt.get("inputs", {}).get("source_manifest", {}).get("sha256")
        != bam_manifest.get("source_manifest", {}).get("sha256")
        or bam_manifest_receipt.get("inputs", {}).get("topology_csv", {}).get("path")
        != source_path(paths["cohort"])
        or bam_manifest_receipt.get("inputs", {}).get("topology_csv", {}).get("sha256")
        != sha256(paths["cohort"])
    ):
        raise ContractError("multi-BAM input manifest validation receipt mismatch")

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

    cohort_by_sample = {row["sample"]: row for row in cohort}
    bam_rows_raw = bam_manifest.get("datasets")
    if not isinstance(bam_rows_raw, list) or len(bam_rows_raw) != EXPECTED_DATASET_N:
        raise ContractError("multi-BAM input manifest must contain seven dataset rows")
    bam_by_sample = {
        str(row.get("dataset_id")): row
        for row in bam_rows_raw
        if isinstance(row, dict)
    }
    if set(bam_by_sample) != set(cohort_by_sample):
        raise ContractError(
            "BAM/topology dataset set mismatch: "
            f"bam={sorted(bam_by_sample)} topology={sorted(cohort_by_sample)}"
        )

    topology_counts: list[dict[str, Any]] = []
    bam_input_identity: list[dict[str, Any]] = []
    producer_tag_details: list[dict[str, Any]] = []
    read_tag_rates: list[dict[str, Any]] = []
    for sample in (row["sample"] for row in cohort):
        cohort_row = cohort_by_sample[sample]
        bam_row = bam_by_sample[sample]
        if (
            bam_row.get("biological_id") != cohort_row["biological_id"]
            or bool(bam_row.get("technical_replicate"))
            != (cohort_row["technical_replicate"] == "yes")
        ):
            raise ContractError(f"{sample}: BAM/topology identity contract mismatch")
        for metric, count in (
            ("Analysis regions", cohort_row["regions"]),
            ("Distinct sSNV loci", cohort_row["distinct_ssnv"]),
        ):
            topology_counts.append({
                "sample": sample,
                "sample_label": (
                    f"{sample} · technical"
                    if cohort_row["technical_replicate"] == "yes"
                    else sample
                ),
                "biological_id": cohort_row["biological_id"],
                "technical_replicate": cohort_row["technical_replicate"],
                "metric": metric,
                "count": count,
            })

        tumor = bam_row.get("bam")
        normal = bam_row.get("paired_normal_bam")
        reference = bam_row.get("reference")
        compatibility = bam_row.get("input_compatibility")
        tags = bam_row.get("producer_read_tags")
        if not all(isinstance(value, dict) for value in (tumor, normal, reference, compatibility, tags)):
            raise ContractError(f"{sample}: incomplete BAM intake row")
        bam_input_identity.append({
            "sample": sample,
            "biological_id": cohort_row["biological_id"],
            "technical_replicate": cohort_row["technical_replicate"],
            "platform": bam_row.get("platform"),
            "tumor_bam_path": tumor.get("path"),
            "tumor_bam_gib": number(tumor.get("size_bytes"), f"{sample}.tumor_size") / (1024 ** 3),
            "tumor_quickcheck": tumor.get("quickcheck_status"),
            "tumor_bounded_payload": tumor.get("bounded_payload_status"),
            "tumor_strict_identity": tumor.get("strict_storage_identity_status"),
            "tumor_identity_differences": ", ".join(tumor.get("strict_differing_fields", [])),
            "tumor_bai_sha": tumor.get("index_sha256_status"),
            "normal_bam_path": normal.get("path"),
            "normal_bam_gib": number(normal.get("size_bytes"), f"{sample}.normal_size") / (1024 ** 3),
            "normal_quickcheck": normal.get("quickcheck_status"),
            "normal_bounded_payload": normal.get("bounded_payload_status"),
            "normal_strict_identity": normal.get("strict_storage_identity_status"),
            "normal_bai_sha": normal.get("index_sha256_status"),
            "reference_path": reference.get("path"),
            "reference_bounded_payload": reference.get("bounded_payload_status"),
            "reference_strict_identity": reference.get("strict_storage_identity_status"),
            "reference_dictionary": compatibility.get("tumor_reference_dictionary_status"),
            "read_group_identity": compatibility.get("read_group_identity_status"),
            "caller_input_path": compatibility.get("caller_input_path"),
            "source_directory_pairing": compatibility.get("source_directory_pairing_status"),
            "full_bam_sha256": "NOT_COLLECTED" if tumor.get("full_content_sha256") is None else "AVAILABLE",
        })

        total_alignments = integer(tags.get("total_alignment_records"), f"{sample}.tag_total")
        hp_alignments = integer(tags.get("hp_assigned_alignment_records"), f"{sample}.tag_hp")
        hp_ps_alignments = integer(tags.get("hp_and_ps_alignment_records"), f"{sample}.tag_hp_ps")
        stored_hp_rate = number(
            tags.get("hp_assigned_rate_all_alignment_records"), f"{sample}.hp_rate"
        )
        stored_hp_ps_all_rate = number(
            tags.get("hp_and_ps_rate_all_alignment_records"), f"{sample}.hp_ps_all_rate"
        )
        assert_close(stored_hp_rate, ratio(hp_alignments, total_alignments, f"{sample}.HP/all"), f"{sample}.HP/all")
        assert_close(
            stored_hp_ps_all_rate,
            ratio(hp_ps_alignments, total_alignments, f"{sample}.HP+PS/all"),
            f"{sample}.HP+PS/all",
        )
        hp_ps_among_hp = ratio(hp_ps_alignments, hp_alignments, f"{sample}.HP+PS/HP")
        stored_hp_ps_among_hp = number(
            tags.get("hp_and_ps_rate_among_hp_assigned_alignment_records"),
            f"{sample}.hp_ps_among_hp",
        )
        assert_close(
            stored_hp_ps_among_hp, hp_ps_among_hp, f"{sample}.HP+PS/HP"
        )
        duplicate_identity_rows = integer(
            tags.get("duplicate_identity_rows"), f"{sample}.duplicates"
        )
        duplicate_rate = ratio(
            duplicate_identity_rows, total_alignments, f"{sample}.duplicate/all"
        )
        stored_duplicate_rate = number(
            tags.get("duplicate_identity_rate_all_alignment_records"),
            f"{sample}.duplicate_rate",
        )
        assert_close(
            stored_duplicate_rate, duplicate_rate, f"{sample}.duplicate/all"
        )
        producer_tag_details.append({
            "sample": sample,
            "biological_id": cohort_row["biological_id"],
            "technical_replicate": cohort_row["technical_replicate"],
            "validation_status": tags.get("validation_status"),
            "alignment_population": tags.get("denominator_population"),
            "total_alignments": total_alignments,
            "primary_alignments": integer(tags.get("primary_alignment_records"), f"{sample}.primary"),
            "secondary_alignments": integer(tags.get("secondary_alignment_records"), f"{sample}.secondary"),
            "supplementary_alignments": integer(tags.get("supplementary_alignment_records"), f"{sample}.supplementary"),
            "hp_alignments": hp_alignments,
            "hp_rate_all_alignments": stored_hp_rate,
            "hp_ps_alignments": hp_ps_alignments,
            "hp_ps_rate_all_alignments": stored_hp_ps_all_rate,
            "hp_ps_rate_among_hp": stored_hp_ps_among_hp,
            "duplicate_identity_rows": duplicate_identity_rows,
            "duplicate_identity_rate_all_alignments": stored_duplicate_rate,
            "record_key_missing": integer(tags.get("record_key_missing"), f"{sample}.missing"),
            "record_key_extra": integer(tags.get("record_key_extra"), f"{sample}.extra"),
            "identity_conflicts": integer(tags.get("duplicate_identity_conflicts"), f"{sample}.conflicts"),
            "unknown_hp_category_n": integer(tags.get("unknown_hp_category_n"), f"{sample}.unknown_hp"),
        })
        for metric, numerator, denominator, value in (
            ("HP / all alignments", hp_alignments, total_alignments, stored_hp_rate),
            (
                "HP+PS / HP-tagged",
                hp_ps_alignments,
                hp_alignments,
                stored_hp_ps_among_hp,
            ),
            (
                "Duplicate identities / all alignments",
                duplicate_identity_rows,
                total_alignments,
                stored_duplicate_rate,
            ),
        ):
            read_tag_rates.append({
                "sample": sample,
                "sample_label": (
                    f"{sample} · technical"
                    if cohort_row["technical_replicate"] == "yes"
                    else sample
                ),
                "biological_id": cohort_row["biological_id"],
                "technical_replicate": cohort_row["technical_replicate"],
                "metric": metric,
                "numerator": numerator,
                "denominator": denominator,
                "rate": value,
                "population": (
                    "mapped alignment records; includes primary, secondary, supplementary"
                    if metric == "HP / all alignments"
                    else (
                        "HP-tagged mapped alignment records"
                        if metric == "HP+PS / HP-tagged"
                        else "all mapped alignment records; denominator-composition guard only"
                    )
                ),
            })

    region_counts = [
        row for row in topology_counts if row["metric"] == "Analysis regions"
    ]
    ssnv_counts = [
        row for row in topology_counts if row["metric"] == "Distinct sSNV loci"
    ]
    if len(region_counts) != EXPECTED_DATASET_N or len(ssnv_counts) != EXPECTED_DATASET_N:
        raise ContractError("topology opportunity datasets must each contain seven rows")

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
    coverage_by_key = {
        (item["k_bin"], item["bundle"]): item for item in coverage
    }
    coverage_denominator_rail: list[dict[str, Any]] = []
    for k_bin, order in sorted(k_order.items(), key=lambda item: item[1]):
        v1 = coverage_by_key[(k_bin, "v1")]
        v3 = coverage_by_key[(k_bin, "v3")]
        if v1["denominator"] != v3["denominator"]:
            raise ContractError(f"coverage denominator differs across bundles at active-k {k_bin}")
        coverage_denominator_rail.append({
            "sample": "HCC1395",
            "k_bin": k_bin,
            "k_order": order,
            "denominator": v1["denominator"],
            "v1_linked": v1["linked"],
            "v1_coverage": v1["coverage"],
            "v3_linked": v3["linked"],
            "v3_coverage": v3["coverage"],
        })

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
        if (
            row["scope"] == "topology_linked_distinct_ssnv"
            and row["interpretation_gate"] != "DOUBLE_DIPPING_INVALID_AS_EVIDENCE"
        ):
            for rate_type, significant_n, significant_rate in (
                ("raw p≤0.05", raw_n, raw_rate),
                ("BH q≤0.05", bh_n, bh_rate),
            ):
                axis_rates_base.append({
                    "sample": "HCC1395",
                    "bundle": row["bundle"],
                    "axis": row["axis"],
                    "axis_display": (
                        "cluster · INVALID (double-dipping)"
                        if row["interpretation_gate"]
                        == "DOUBLE_DIPPING_INVALID_AS_EVIDENCE"
                        else row["axis"]
                    ),
                    "rate_type": rate_type,
                    "series": f"{row['bundle']} · {rate_type}",
                    "significant_n": significant_n,
                    "tested_n": tested_n,
                    "rate": significant_rate,
                    "interpretation_gate": row["interpretation_gate"],
                })

    topology_axis = {
        (item["axis"], item["bundle"]): item
        for item in axis_details
        if item["analysis_scope"] == "topology_linked_distinct_ssnv"
    }
    axis_denominator_rail: list[dict[str, Any]] = []
    axis_order = []
    for item in axis_details:
        if (
            item["analysis_scope"] == "topology_linked_distinct_ssnv"
            and item["axis"] not in axis_order
        ):
            axis_order.append(item["axis"])
    for axis in axis_order:
        v1 = topology_axis[(axis, "v1")]
        v3 = topology_axis[(axis, "v3")]
        gate = v1["interpretation_gate"]
        if v3["interpretation_gate"] != gate:
            raise ContractError(f"axis interpretation gate differs across bundles: {axis}")
        axis_denominator_rail.append({
            "sample": "HCC1395",
            "axis": axis,
            "display_status": (
                "INVALID · DOUBLE-DIPPING"
                if gate == "DOUBLE_DIPPING_INVALID_AS_EVIDENCE"
                else "EXPLORATORY"
            ),
            "v1_tested": v1["tested_n"],
            "v1_raw_n": v1["raw_sig_n"],
            "v1_raw_rate": v1["raw_sig_rate"],
            "v1_bh_n": v1["bh_sig_n"],
            "v1_bh_rate": v1["bh_sig_rate"],
            "v3_tested": v3["tested_n"],
            "v3_raw_n": v3["raw_sig_n"],
            "v3_raw_rate": v3["raw_sig_rate"],
            "v3_bh_n": v3["bh_sig_n"],
            "v3_bh_rate": v3["bh_sig_rate"],
            "v1_raw_nd": f"{v1['raw_sig_n']:,} / {v1['tested_n']:,} · {v1['raw_sig_rate']:.1%}",
            "v1_bh_nd": f"{v1['bh_sig_n']:,} / {v1['tested_n']:,} · {v1['bh_sig_rate']:.1%}",
            "v3_raw_nd": f"{v3['raw_sig_n']:,} / {v3['tested_n']:,} · {v3['raw_sig_rate']:.1%}",
            "v3_bh_nd": f"{v3['bh_sig_n']:,} / {v3['tested_n']:,} · {v3['bh_sig_rate']:.1%}",
            "interpretation_gate": gate,
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
    bam_scope_summary: list[dict[str, Any]] = []
    scope_groups = [(ALL_DATASETS_SCOPE, cohort)] + [
        (item["sample"], [item]) for item in cohort
    ]
    availability: list[dict[str, Any]] = []
    for sample_filter, group in scope_groups:
        bam_group = [bam_by_sample[item["sample"]] for item in group]
        dataset_count = len(group)
        biological_count = len({item["biological_id"] for item in group})
        technical_replicate_n = sum(
            item["technical_replicate"] == "yes" for item in group
        )
        hash_match_n = sum(item["hash_match"] == "PASS" for item in group)
        sample_identity_match_n = sum(
            item["sample_identity_match"] == "PASS" for item in group
        )
        receipt_technical_pass_n = sum(
            item["receipt_technical_pass"] == "PASS" for item in group
        )
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
        biological_macro_tree_q1 = (
            percentile_cont(
                (item["tree_coverage"] for item in biological_rows),
                0.25,
                f"{sample_filter}.tree coverage P25",
            )
            if biological_rows else None
        )
        biological_macro_tree_q3 = (
            percentile_cont(
                (item["tree_coverage"] for item in biological_rows),
                0.75,
                f"{sample_filter}.tree coverage P75",
            )
            if biological_rows else None
        )
        biological_macro_unique_q1 = (
            percentile_cont(
                (item["unique_among_tree"] for item in biological_rows),
                0.25,
                f"{sample_filter}.unique/tree P25",
            )
            if biological_rows else None
        )
        biological_macro_unique_q3 = (
            percentile_cont(
                (item["unique_among_tree"] for item in biological_rows),
                0.75,
                f"{sample_filter}.unique/tree P75",
            )
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
            "technical_replicate_n": technical_replicate_n,
            "hash_match_n": hash_match_n,
            "sample_identity_match_n": sample_identity_match_n,
            "receipt_technical_pass_n": receipt_technical_pass_n,
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
            "biological_macro_tree_coverage_q1": biological_macro_tree_q1,
            "biological_macro_tree_coverage_q3": biological_macro_tree_q3,
            "biological_macro_unique_among_tree_median": biological_macro_unique_median,
            "biological_macro_unique_among_tree_q1": biological_macro_unique_q1,
            "biological_macro_unique_among_tree_q3": biological_macro_unique_q3,
            # This is a reviewed inventory count: zero means no audited bundle
            # in the selected scope, not zero biological signal.
            "downstream_available_count": downstream_available_n,
            "downstream_total_count": dataset_count,
        })
        tumor_quickcheck_pass_n = sum(
            row["bam"]["quickcheck_status"] == "PASS" for row in bam_group
        )
        normal_quickcheck_pass_n = sum(
            row["paired_normal_bam"]["quickcheck_status"] == "PASS"
            for row in bam_group
        )
        tumor_bounded_payload_match_n = sum(
            row["bam"]["bounded_payload_status"] == "MATCH" for row in bam_group
        )
        normal_bounded_payload_match_n = sum(
            row["paired_normal_bam"]["bounded_payload_status"] == "MATCH"
            for row in bam_group
        )
        reference_compatible_n = sum(
            row["input_compatibility"]["tumor_reference_dictionary_status"] == "PASS"
            and row["input_compatibility"]["normal_reference_dictionary_status"] == "PASS"
            and row["input_compatibility"]["tumor_normal_dictionary_status"] == "PASS"
            for row in bam_group
        )
        strict_storage_identity_match_n = sum(
            row["bam"]["strict_storage_identity_status"] == "MATCH"
            and row["paired_normal_bam"]["strict_storage_identity_status"] == "MATCH"
            and row["reference"]["strict_storage_identity_status"] == "MATCH"
            for row in bam_group
        )
        mount_device_drift_n = sum(
            all(
                record["strict_storage_identity_status"]
                in {"MATCH", "MOUNT_DEVICE_DRIFT_ONLY"}
                for record in (row["bam"], row["paired_normal_bam"], row["reference"])
            )
            and any(
                record["strict_storage_identity_status"]
                == "MOUNT_DEVICE_DRIFT_ONLY"
                for record in (row["bam"], row["paired_normal_bam"], row["reference"])
            )
            for row in bam_group
        )
        tag_receipt_pass_n = sum(
            row["producer_read_tags"]["validation_status"] == "PASS"
            for row in bam_group
        )
        full_bam_sha256_n = sum(
            row["bam"]["full_content_sha256"] is not None for row in bam_group
        )
        no_read_group_n = sum(
            row["input_compatibility"]["read_group_identity_status"]
            == "UNAVAILABLE_NO_RG"
            for row in bam_group
        )
        cross_directory_pairing_n = sum(
            row["input_compatibility"]["source_directory_pairing_status"]
            == "CROSS_DIRECTORY_REVIEW_REQUIRED"
            for row in bam_group
        )
        bam_scope_summary.append({
            "sample_filter": sample_filter,
            "dataset_count": dataset_count,
            "tumor_quickcheck_pass_n": tumor_quickcheck_pass_n,
            "normal_quickcheck_pass_n": normal_quickcheck_pass_n,
            "tumor_bounded_payload_match_n": tumor_bounded_payload_match_n,
            "normal_bounded_payload_match_n": normal_bounded_payload_match_n,
            "reference_compatible_n": reference_compatible_n,
            "strict_storage_identity_match_n": strict_storage_identity_match_n,
            "mount_device_drift_n": mount_device_drift_n,
            "tag_receipt_pass_n": tag_receipt_pass_n,
            "full_bam_sha256_n": full_bam_sha256_n,
            "no_read_group_n": no_read_group_n,
            "cross_directory_pairing_n": cross_directory_pairing_n,
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
        if strict_storage_identity_match_n == dataset_count:
            bounded_identity_status = "AVAILABLE_BOUNDED"
            bounded_identity_note = (
                "Tumor/normal/reference strict storage identity is MATCH for every selected "
                "dataset; quickcheck and three frozen 1 MiB chunks also pass. Full BAM "
                "SHA-256 remains uncollected, so this is not full-content attestation."
            )
        elif strict_storage_identity_match_n + mount_device_drift_n == dataset_count:
            bounded_identity_status = "AVAILABLE_BOUNDED_WITH_METADATA_DRIFT"
            bounded_identity_note = (
                f"{strict_storage_identity_match_n}/{dataset_count} selected datasets have "
                f"all-role strict MATCH; {mount_device_drift_n}/{dataset_count} are limited "
                "to mount-device metadata drift while frozen payload chunks match. Full BAM "
                "SHA-256 remains uncollected."
            )
        else:
            bounded_identity_status = "REVIEW_REQUIRED"
            bounded_identity_note = (
                "At least one selected dataset is neither all-role strict MATCH nor bounded "
                "mount-device-only drift; inspect the per-dataset identity table."
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
                "Tumor + normal bounded payload identity",
                min(
                    tumor_quickcheck_pass_n,
                    normal_quickcheck_pass_n,
                    tumor_bounded_payload_match_n,
                    normal_bounded_payload_match_n,
                ),
                dataset_count,
                bounded_identity_status,
                bounded_identity_note,
            ),
            (
                "Reference dictionary compatibility",
                reference_compatible_n,
                dataset_count,
                "AVAILABLE_BOUNDED",
                "Tumor, paired normal, and reference FAI dictionaries reconcile for every selected dataset.",
            ),
            (
                "LongPhase-S output alignment tags",
                tag_receipt_pass_n,
                dataset_count,
                "AVAILABLE_EXISTING_PRODUCER_RECEIPT",
                "Rates use all captured mapped alignment records, including primary, secondary, and supplementary records; this is not a primary-read metric.",
            ),
            (
                "Full BAM content SHA-256",
                None if full_bam_sha256_n == 0 else full_bam_sha256_n,
                dataset_count,
                "UNAVAILABLE_NOT_COLLECTED",
                "No full-content BAM SHA-256 scan was run; bounded chunk identity must not be promoted to full identity.",
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
                "Primary-read HP/PS tagging rate",
                None,
                dataset_count,
                "UNAVAILABLE_NOT_MEASURED",
                "Existing producer receipts use all mapped alignment records; a primary-read-only numerator and denominator are not present.",
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
                "technical_replicate_n = datasets marked technical_replicate=true in the selected view.",
                "hash_match_n, sample_identity_match_n, and receipt_technical_pass_n are separate technical gate counts.",
                "integrity_pass_rate = datasets passing hash, sample identity, and technical receipt gates / datasets in view.",
                "technical_gate_pass_n = datasets with receipt_all_pass=true.",
                "all_families_complete_n = datasets with receipt_all_mutation_bearing_families_complete=true.",
                "all_objective_certified_n = datasets with receipt_all_units_objective_certified=true.",
                "biological macro median/P25/P75 use SQL PERCENTILE_CONT-style linear interpolation and exclude rows marked technical_replicate=true.",
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
            "q_bam_scope_summary",
            "Bounded BAM intake summary",
            "bam_scope_summary",
            [bam_manifest_path, bam_manifest_receipt_path],
            generated_at,
            "Selected-scope tumor/normal/reference readiness, strict-identity categories, missing read groups, and caller-directory review counts.",
            [
                "payload sample match = first/middle/last frozen 1 MiB chunks plus size/index checks; it is not a full-content BAM hash.",
                "strict storage identity match = all three roles MATCH; mount-device drift = all three roles are MATCH or MOUNT_DEVICE_DRIFT_ONLY and at least one role has that bounded drift status.",
                "tag receipt pass uses existing producer sidecar validation, not a new primary-read scan.",
            ],
            "sample_filter",
        ),
        make_snapshot_source(
            "q_bam_input_identity",
            "Per-dataset BAM/normal/reference intake",
            "bam_input_identity",
            [bam_manifest_path, bam_manifest_receipt_path],
            generated_at,
            "Per-dataset bounded payload, BAM index, header dictionary, source-directory pairing, and full-hash availability details.",
            [
                "tumor_bam_gib and normal_bam_gib = file size bytes / 1024^3.",
                "UNAVAILABLE_NO_RG means BAM headers contain no @RG SM field, so header identity cannot be asserted.",
                "CROSS_DIRECTORY_REVIEW_REQUIRED is a transfer/confound flag, not evidence of corruption.",
            ],
            "sample_filter, sample",
        ),
        make_snapshot_source(
            "q_read_tag_rates",
            "LongPhase-S output alignment tag rates",
            "read_tag_rates",
            [bam_manifest_path, bam_manifest_receipt_path],
            generated_at,
            "Two exact tag rates from existing producer receipts, with denominator population retained for every dataset.",
            [
                "HP / all alignments = HP-assigned mapped alignment records / all captured mapped alignment records, including primary, secondary, and supplementary.",
                "HP+PS / HP-tagged = alignment records carrying both HP and PS / HP-assigned alignment records.",
                "These are alignment-record metrics and must not be relabeled as unique-read or primary-read rates.",
            ],
            "sample_filter, sample, metric",
        ),
        make_snapshot_source(
            "q_producer_tag_details",
            "LongPhase-S tag conservation detail",
            "producer_tag_details",
            [bam_manifest_path, bam_manifest_receipt_path],
            generated_at,
            "Exact alignment-class and HP/PS counts plus missing/extra/conflict conservation fields from producer validation receipts.",
            [
                "total alignments = primary + secondary + supplementary captured mapped alignment records.",
                "HP+PS among HP = HP-and-PS records / HP-assigned records.",
                "duplicate identity rows are reported separately from duplicate identity conflicts.",
            ],
            "sample_filter, sample",
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
            "q_region_counts",
            "Per-dataset analysis-region opportunity",
            "region_counts",
            [paths["cohort"]],
            generated_at,
            "Analysis-region counts are shown separately from distinct sSNV counts to avoid a dual-axis comparison.",
            ["count = region_n for one technical dataset; it is not variant burden or recall."],
            "sample_filter, sample",
        ),
        make_snapshot_source(
            "q_ssnv_counts",
            "Per-dataset distinct sSNV opportunity",
            "ssnv_counts",
            [paths["cohort"]],
            generated_at,
            "Distinct topology sSNV counts are shown separately from region opportunity counts.",
            ["count = distinct_ssnv_n for one technical dataset; it is not truth-set sensitivity."],
            "sample_filter, sample",
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
            "q_coverage_denominator_rail",
            "HCC1395 active-k exact denominator rail",
            "coverage_denominator_rail",
            [paths["coverage"]],
            generated_at,
            "Compact exact N/D companion to the active-k chart, grouped into one row per k bin.",
            ["v1/v3 coverage = linked / the displayed denominator."],
            "sample_filter, k_order",
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
            "q_axis_denominator_rail",
            "HCC1395 axis exact denominator and validity rail",
            "axis_denominator_rail",
            [paths["axis"]],
            generated_at,
            "One row per topology-linked axis with v1/v3 exact tested, raw, and BH counts; invalid cluster is retained here but excluded from the chart.",
            [
                "raw rate = raw-significant / tested.",
                "BH rate = BH-significant / tested.",
                "cluster is invalid as evidence because of circular double-dipping.",
            ],
            "sample_filter, axis",
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
                {"label": "目前資料集", "field": "dataset_count", "format": "number"},
                {"label": "技術 replicate", "field": "technical_replicate_n", "format": "number"},
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
                    "label": "Hash PASS",
                    "field": "hash_match_n",
                    "format": "number",
                },
                {
                    "label": "Identity PASS",
                    "field": "sample_identity_match_n",
                    "format": "number",
                },
                {
                    "label": "Receipt PASS",
                    "field": "receipt_technical_pass_n",
                    "format": "number",
                },
                {
                    "label": "目前資料集",
                    "field": "technical_gate_total_n",
                    "format": "number",
                },
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
                {
                    "label": "P25",
                    "field": "biological_macro_tree_coverage_q1",
                    "format": "percent",
                },
                {
                    "label": "P75",
                    "field": "biological_macro_tree_coverage_q3",
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
                {
                    "label": "P25",
                    "field": "biological_macro_unique_among_tree_q1",
                    "format": "percent",
                },
                {
                    "label": "P75",
                    "field": "biological_macro_unique_among_tree_q3",
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
        {
            "id": "c_bam_payload_identity",
            "dataset": "bam_scope_summary",
            "sourceId": "q_bam_scope_summary",
            "description": "Tumor／normal／reference bounded payload只涵蓋固定三段 1 MiB／index／metadata；strict MATCH與MOUNT_DEVICE_DRIFT_ONLY依各資料集動態重算，兩者都不是full-content attestation。",
            "metrics": [
                {
                    "label": "Tumor payload MATCH",
                    "field": "tumor_bounded_payload_match_n",
                    "format": "number",
                },
                {
                    "label": "Normal payload MATCH",
                    "field": "normal_bounded_payload_match_n",
                    "format": "number",
                },
                {
                    "label": "Strict identity MATCH",
                    "field": "strict_storage_identity_match_n",
                    "format": "number",
                },
                {
                    "label": "Mount metadata drift",
                    "field": "mount_device_drift_n",
                    "format": "number",
                },
            ],
        },
        {
            "id": "c_bam_input_readiness",
            "dataset": "bam_scope_summary",
            "sourceId": "q_bam_scope_summary",
            "description": "Tumor/normal quickcheck、reference dictionary 與既有 producer tag receipt；不包含 depth、N50、MM/ML 或 truth benchmark。",
            "metrics": [
                {
                    "label": "Tumor quickcheck",
                    "field": "tumor_quickcheck_pass_n",
                    "format": "number",
                },
                {
                    "label": "Normal quickcheck",
                    "field": "normal_quickcheck_pass_n",
                    "format": "number",
                },
                {
                    "label": "Reference compatible",
                    "field": "reference_compatible_n",
                    "format": "number",
                },
                {
                    "label": "Tag receipt PASS",
                    "field": "tag_receipt_pass_n",
                    "format": "number",
                },
            ],
        },
    ]

    charts = [
        {
            "id": "ch_topology_rates",
            "title": "各資料集 topology coverage 與唯一性",
            "subtitle": "7 組技術資料代表 6 個生物樣本；圖表資料保留精確分子與分母。",
            "type": "bar",
            "intent": "comparison",
            "dataset": "topology_rates",
            "sourceId": "q_topology_rates",
            "valueFormat": "percent",
            "encodings": {
                "x": {"field": "sample", "type": "nominal", "label": "資料集"},
                "y": {
                    "field": "rate",
                    "type": "quantitative",
                    "label": "比例",
                    "format": "percent",
                },
                "color": {"field": "metric", "type": "nominal", "label": "指標"},
                "tooltip": [
                    {"field": "biological_id", "type": "nominal", "label": "生物樣本"},
                    {"field": "numerator", "type": "quantitative", "label": "分子"},
                    {"field": "denominator", "type": "quantitative", "label": "分母"},
                    {"field": "technical_replicate", "type": "nominal", "label": "技術重複"},
                ],
            },
            "maxRows": 40,
        },
        {
            "id": "ch_region_counts",
            "title": "分析區域機會量",
            "subtitle": "水平列出全部 7 個資料集；technical replicate 直接寫入標籤，與 sSNV 分圖。",
            "type": "bar",
            "intent": "comparison",
            "dataset": "region_counts",
            "sourceId": "q_region_counts",
            "valueFormat": "number",
            "encodings": {
                "x": {"field": "sample_label", "type": "nominal", "label": "資料集"},
                "y": {
                    "field": "count",
                    "type": "quantitative",
                    "label": "分析區域數",
                    "format": "number",
                },
                "color": {
                    "field": "technical_replicate",
                    "type": "nominal",
                    "label": "技術重複",
                },
                "tooltip": [
                    {"field": "biological_id", "type": "nominal", "label": "生物樣本"},
                    {"field": "count", "type": "quantitative", "label": "區域數"},
                ],
            },
            "maxRows": 20,
        },
        {
            "id": "ch_ssnv_counts",
            "title": "Distinct sSNV 機會量",
            "subtitle": "水平列出全部 7 個資料集；不是 mutation burden、recall 或樣本品質排名。",
            "type": "bar",
            "intent": "comparison",
            "dataset": "ssnv_counts",
            "sourceId": "q_ssnv_counts",
            "valueFormat": "number",
            "encodings": {
                "x": {"field": "sample_label", "type": "nominal", "label": "資料集"},
                "y": {
                    "field": "count",
                    "type": "quantitative",
                    "label": "Distinct sSNV",
                    "format": "number",
                },
                "color": {
                    "field": "technical_replicate",
                    "type": "nominal",
                    "label": "技術重複",
                },
                "tooltip": [
                    {"field": "biological_id", "type": "nominal", "label": "生物樣本"},
                    {"field": "count", "type": "quantitative", "label": "Distinct sSNV"},
                ],
            },
            "maxRows": 20,
        },
        {
            "id": "ch_read_tag_rates",
            "title": "LongPhase-S 輸出 alignment tag 比例",
            "subtitle": "HP/all、HP+PS/HP-tagged 與 duplicate/all 分開呈現；duplicate 是分母組成 guard，不是 phasing accuracy。",
            "type": "bar",
            "intent": "comparison",
            "dataset": "read_tag_rates",
            "sourceId": "q_read_tag_rates",
            "valueFormat": "percent",
            "encodings": {
                "x": {"field": "sample_label", "type": "nominal", "label": "資料集"},
                "y": {
                    "field": "rate",
                    "type": "quantitative",
                    "label": "Alignment record 比例",
                    "format": "percent",
                },
                "color": {"field": "metric", "type": "nominal", "label": "分母口徑"},
                "tooltip": [
                    {"field": "numerator", "type": "quantitative", "label": "分子"},
                    {"field": "denominator", "type": "quantitative", "label": "分母"},
                    {"field": "population", "type": "nominal", "label": "母體"},
                ],
            },
            "maxRows": 40,
        },
        {
            "id": "ch_coverage_by_k",
            "title": "HCC1395：各 active-k 區間的 ISM 連結率",
            "subtitle": "v1 與 v3 使用不同視窗；此圖是描述性可用性診斷，不是受控 A/B 比較。",
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
                    "label": "ISM 連結率",
                    "format": "percent",
                },
                "color": {"field": "bundle", "type": "nominal", "label": "Bundle"},
                "tooltip": [
                    {"field": "linked", "type": "quantitative", "label": "已連結"},
                    {"field": "denominator", "type": "quantitative", "label": "分母"},
                ],
            },
            "maxRows": 40,
        },
        {
            "id": "ch_axis_rates",
            "title": "HCC1395：有效 topology-linked 甲基化軸比例",
            "subtitle": "raw p 與 family 內 BH 僅供探索；cluster 因 double-dipping 已自圖排除，精確值保留在 validity rail。",
            "type": "bar",
            "intent": "comparison",
            "dataset": "axis_rates",
            "sourceId": "q_axis_rates",
            "valueFormat": "percent",
            "encodings": {
                "x": {"field": "axis_display", "type": "nominal", "label": "甲基化軸"},
                "y": {
                    "field": "rate",
                    "type": "quantitative",
                    "label": "顯著列比例",
                    "format": "percent",
                },
                "color": {"field": "series", "type": "nominal", "label": "Bundle · 門檻"},
                "tooltip": [
                    {"field": "significant_n", "type": "quantitative", "label": "顯著列"},
                    {"field": "tested_n", "type": "quantitative", "label": "受測列"},
                    {"field": "interpretation_gate", "type": "nominal", "label": "解讀 gate"},
                ],
            },
            "maxRows": 50,
        },
    ]

    tables = [
        {
            "id": "t_availability",
            "title": "跨樣本資料可用性",
            "subtitle": "不可用值一律為 null（空白），絕不以 0 代替；本表分開呈現實測覆蓋與資料缺口。",
            "dataset": "availability",
            "sourceId": "q_availability",
            "defaultSort": {"field": "component", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "component", "label": "資料構面", "type": "text"},
                {"field": "available_datasets", "label": "可用", "format": "number"},
                {"field": "expected_datasets", "label": "預期", "format": "number"},
                {"field": "coverage", "label": "覆蓋率", "format": "percent"},
                {"field": "status", "label": "狀態", "type": "text"},
                {"field": "note", "label": "證據界線", "type": "text"},
            ],
        },
        {
            "id": "t_bam_input_identity",
            "title": "BAM／paired-normal／reference 輸入綁定與 bounded payload 身分",
            "subtitle": "Bounded payload MATCH 為固定三段 1 MiB + index；不是 biological identity 或 full-content hash，且 14 BAM 都沒有 @RG。",
            "dataset": "bam_input_identity",
            "sourceId": "q_bam_input_identity",
            "defaultSort": {"field": "sample", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "sample", "label": "資料集", "type": "text"},
                {"field": "biological_id", "label": "生物樣本", "type": "text"},
                {"field": "technical_replicate", "label": "技術重複", "type": "text"},
                {"field": "platform", "label": "平台", "type": "text"},
                {"field": "tumor_bam_gib", "label": "Tumor GiB", "format": "number"},
                {"field": "normal_bam_gib", "label": "Normal GiB", "format": "number"},
                {"field": "tumor_quickcheck", "label": "Tumor QC", "type": "text"},
                {"field": "normal_quickcheck", "label": "Normal QC", "type": "text"},
                {"field": "tumor_bounded_payload", "label": "Tumor bounded", "type": "text"},
                {"field": "normal_bounded_payload", "label": "Normal bounded", "type": "text"},
                {"field": "tumor_strict_identity", "label": "Strict identity", "type": "text"},
                {"field": "normal_strict_identity", "label": "Normal strict", "type": "text"},
                {"field": "reference_strict_identity", "label": "Reference strict", "type": "text"},
                {"field": "tumor_identity_differences", "label": "Strict 差異", "type": "text"},
                {"field": "tumor_bai_sha", "label": "Tumor BAI", "type": "text"},
                {"field": "normal_bai_sha", "label": "Normal BAI", "type": "text"},
                {"field": "reference_dictionary", "label": "Reference dict", "type": "text"},
                {"field": "read_group_identity", "label": "@RG identity", "type": "text"},
                {"field": "source_directory_pairing", "label": "BAM ↔ caller", "type": "text"},
                {"field": "full_bam_sha256", "label": "Full BAM SHA", "type": "text"},
                {"field": "tumor_bam_path", "label": "Tumor BAM", "type": "text"},
                {"field": "normal_bam_path", "label": "Normal BAM", "type": "text"},
                {"field": "reference_path", "label": "Reference", "type": "text"},
                {"field": "caller_input_path", "label": "Caller VCF", "type": "text"},
            ],
        },
        {
            "id": "t_producer_tag_details",
            "title": "LongPhase-S alignment tag 分母與守恆明細",
            "subtitle": "每列均保留 alignment class、HP、HP+PS 與 identity conservation；duplicate rows 不等同 conflicts。",
            "dataset": "producer_tag_details",
            "sourceId": "q_producer_tag_details",
            "defaultSort": {"field": "sample", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "sample", "label": "資料集", "type": "text"},
                {"field": "validation_status", "label": "Receipt", "type": "text"},
                {"field": "total_alignments", "label": "All alignments", "format": "number"},
                {"field": "primary_alignments", "label": "Primary", "format": "number"},
                {"field": "secondary_alignments", "label": "Secondary", "format": "number"},
                {"field": "supplementary_alignments", "label": "Supplementary", "format": "number"},
                {"field": "hp_alignments", "label": "HP", "format": "number"},
                {"field": "hp_rate_all_alignments", "label": "HP / all", "format": "percent"},
                {"field": "hp_ps_alignments", "label": "HP+PS", "format": "number"},
                {"field": "hp_ps_rate_all_alignments", "label": "HP+PS / all", "format": "percent"},
                {"field": "hp_ps_rate_among_hp", "label": "HP+PS / HP", "format": "percent"},
                {"field": "duplicate_identity_rows", "label": "Duplicate rows", "format": "number"},
                {"field": "duplicate_identity_rate_all_alignments", "label": "Duplicate / all", "format": "percent"},
                {"field": "record_key_missing", "label": "Missing", "format": "number"},
                {"field": "record_key_extra", "label": "Extra", "format": "number"},
                {"field": "identity_conflicts", "label": "Conflicts", "format": "number"},
                {"field": "unknown_hp_category_n", "label": "Unknown HP", "format": "number"},
                {"field": "alignment_population", "label": "分母母體", "type": "text"},
            ],
        },
        {
            "id": "t_bundle_overview",
            "title": "HCC1395 v1/v3 bundle 診斷",
            "subtitle": "兩個 bundle 仍是 blocked observation；snapshot 內的比例以 0–1 小數儲存。",
            "dataset": "bundle_overview",
            "sourceId": "q_bundle_overview",
            "defaultSort": {"field": "bundle", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "bundle", "label": "Bundle", "type": "text"},
                {"field": "validation_status", "label": "驗證狀態", "type": "text"},
                {"field": "regions", "label": "區域數", "format": "number"},
                {"field": "tree_coverage", "label": "Tree coverage", "format": "percent"},
                {"field": "unique_among_tree", "label": "Unique / tree", "format": "percent"},
                {"field": "ism_linked", "label": "ISM 已連結", "format": "number"},
                {"field": "ism_linkage", "label": "ISM 連結率", "format": "percent"},
                {"field": "selfcheck_fail", "label": "FAIL", "format": "number"},
                {"field": "selfcheck_skip", "label": "SKIP", "format": "number"},
                {"field": "panel_actual_bytes", "label": "Panel bytes", "format": "number"},
                {"field": "panel_unreported_bytes", "label": "未回報 bytes", "format": "number"},
            ],
        },
        {
            "id": "t_coverage_denominator_rail",
            "title": "Active-k exact N/D",
            "subtitle": "圖表旁的 8-bin 精確分母；v1/v3 共用 denominator，但 linked numerator 分開保存。",
            "dataset": "coverage_denominator_rail",
            "sourceId": "q_coverage_denominator_rail",
            "defaultSort": {"field": "k_bin", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "k_bin", "label": "Active k", "type": "text"},
                {"field": "denominator", "label": "Denominator", "format": "number"},
                {"field": "v1_linked", "label": "v1 linked", "format": "number"},
                {"field": "v1_coverage", "label": "v1 rate", "format": "percent"},
                {"field": "v3_linked", "label": "v3 linked", "format": "number"},
                {"field": "v3_coverage", "label": "v3 rate", "format": "percent"},
            ],
        },
        {
            "id": "t_axis_denominator_rail",
            "title": "甲基化軸 exact N/D 與 validity rail",
            "subtitle": "cluster 保留精確數值但標為 INVALID 並自有效 axis 圖排除；其餘仍僅 exploratory。",
            "dataset": "axis_denominator_rail",
            "sourceId": "q_axis_denominator_rail",
            "defaultSort": {"field": "axis", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "axis", "label": "Axis", "type": "text"},
                {"field": "display_status", "label": "Validity", "type": "text"},
                {"field": "v1_raw_nd", "label": "v1 raw N / D · %", "type": "text"},
                {"field": "v1_bh_nd", "label": "v1 BH N / D · %", "type": "text"},
                {"field": "v3_raw_nd", "label": "v3 raw N / D · %", "type": "text"},
                {"field": "v3_bh_nd", "label": "v3 BH N / D · %", "type": "text"},
            ],
        },
        {
            "id": "t_topology_details",
            "title": "7 組資料的 topology 明細",
            "subtitle": "保留每組資料的分母；HCC1395_DORADO 明確標為技術重複。",
            "dataset": "topology_details",
            "sourceId": "q_topology_details",
            "defaultSort": {"field": "tree_coverage", "direction": "desc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "sample", "label": "資料集", "type": "text"},
                {"field": "biological_id", "label": "生物樣本", "type": "text"},
                {"field": "technical_replicate", "label": "技術重複", "type": "text"},
                {"field": "regions", "label": "區域數", "format": "number"},
                {"field": "distinct_ssnv", "label": "Distinct sSNV", "format": "number"},
                {"field": "tree_n", "label": "Tree 數", "format": "number"},
                {"field": "tree_coverage", "label": "Tree coverage", "format": "percent"},
                {"field": "unique_tree_n", "label": "Unique tree 數", "format": "number"},
                {"field": "unique_among_tree", "label": "Unique / tree", "format": "percent"},
                {"field": "hash_match", "label": "Hash", "type": "text"},
                {"field": "sample_identity_match", "label": "樣本身分", "type": "text"},
                {"field": "receipt_technical_pass", "label": "技術收據", "type": "text"},
                {"field": "all_families_complete_gate", "label": "全 family", "type": "text"},
                {"field": "all_objective_certified_gate", "label": "全 objective", "type": "text"},
            ],
        },
        {
            "id": "t_axis_details",
            "title": "HCC1395 甲基化軸明細",
            "subtitle": "不同軸的 effect unit 不同；cluster 軸不可作為有效證據。",
            "dataset": "axis_details",
            "sourceId": "q_axis_details",
            "defaultSort": {"field": "tested_n", "direction": "desc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "bundle", "label": "Bundle", "type": "text"},
                {"field": "analysis_scope", "label": "分析範圍", "type": "text"},
                {"field": "axis", "label": "甲基化軸", "type": "text"},
                {"field": "tested_n", "label": "受測列", "format": "number"},
                {"field": "raw_sig_rate", "label": "Raw p≤.05", "format": "percent"},
                {"field": "bh_sig_rate", "label": "BH q≤.05", "format": "percent"},
                {"field": "effect_median", "label": "Effect 中位數", "format": "number"},
                {"field": "effect_unit", "label": "Effect unit", "type": "text"},
                {"field": "interpretation_gate", "label": "解讀 gate", "type": "text"},
            ],
        },
        {
            "id": "t_lca_summary",
            "title": "HCC1395 lineage/LCA 受控比較 gate",
            "subtitle": "輸入身分未固定，因此觀察到的差值不能解讀為 LCA 因果效果。",
            "dataset": "lca_summary",
            "sourceId": "q_lca_summary",
            "defaultSort": {"field": "sample", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "sample", "label": "樣本", "type": "text"},
                {"field": "shared_chromosomes", "label": "共同 chr", "format": "number"},
                {"field": "same_in_bam_n", "label": "相同 BAM", "format": "number"},
                {"field": "same_in_bam_rate", "label": "相同 BAM 比例", "format": "percent"},
                {"field": "same_threads_n", "label": "相同 threads", "format": "number"},
                {"field": "same_threads_rate", "label": "相同 threads 比例", "format": "percent"},
                {"field": "pre_lv_written", "label": "Pre lv", "format": "number"},
                {"field": "post_lv_written", "label": "Post lv", "format": "number"},
                {"field": "net_new_lv_written", "label": "Net new lv", "format": "number"},
                {"field": "causal_ab_gate", "label": "因果 A/B gate", "type": "text"},
            ],
        },
        {
            "id": "t_visual_qa",
            "title": "Legacy／目前 browser QA 收據",
            "subtitle": "視覺 QA 僅驗證介面行為，不會提高科學 claim ceiling。",
            "dataset": "visual_qa",
            "sourceId": "q_visual_qa",
            "defaultSort": {"field": "surface", "direction": "asc"},
            "density": "dense",
            "layout": "full",
            "columns": [
                {"field": "surface", "label": "介面", "type": "text"},
                {"field": "viewport", "label": "視窗", "type": "text"},
                {"field": "viewport_width", "label": "寬度", "format": "number"},
                {"field": "horizontal_overflow", "label": "橫向溢位", "type": "text"},
                {"field": "browser_error_n", "label": "瀏覽器錯誤", "format": "number"},
                {"field": "evidence_status", "label": "證據狀態", "type": "text"},
                {"field": "screenshot", "label": "截圖收據", "type": "text"},
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
    bam_cohort_scope = bam_scope_summary[0]
    bam_scope_n = bam_cohort_scope["dataset_count"]
    bam_strict_n = bam_cohort_scope["strict_storage_identity_match_n"]
    bam_drift_n = bam_cohort_scope["mount_device_drift_n"]
    if bam_strict_n == bam_scope_n:
        bam_identity_detail_title = (
            f"BAM strict identity 已達 {bam_strict_n}/{bam_scope_n}；還缺什麼？"
        )
        bam_identity_detail_body = (
            f"三段 frozen 1 MiB payload 與 BAI 是 "
            f"{bam_cohort_scope['tumor_bounded_payload_match_n']}/{bam_scope_n} MATCH，"
            f"tumor／normal／reference all-role strict storage identity 也是 "
            f"{bam_strict_n}/{bam_scope_n} MATCH。但 full BAM SHA-256 尚未執行；另有 "
            f"{bam_cohort_scope['no_read_group_n']}/{bam_scope_n} tumor BAM 沒有 "
            f"<code>@RG</code>，且有 {bam_cohort_scope['cross_directory_pairing_n']} 筆 "
            "BAM↔caller cross-directory review，因此仍不可宣稱 full-content、完整生物樣本身分或 platform effect。"
        )
    else:
        bam_identity_detail_title = (
            f"BAM identity 為何只有 {bam_strict_n}/{bam_scope_n} strict MATCH？"
        )
        bam_identity_detail_body = (
            f"三段 frozen 1 MiB payload 與 BAI 是 "
            f"{bam_cohort_scope['tumor_bounded_payload_match_n']}/{bam_scope_n} MATCH；"
            f"all-role strict identity 為 {bam_strict_n}/{bam_scope_n}，另有 "
            f"{bam_drift_n}/{bam_scope_n} 被明確分類為 "
            "<code>MOUNT_DEVICE_DRIFT_ONLY</code>，逐 role 差異保留在 identity table。"
            f"此外 {bam_cohort_scope['no_read_group_n']}/{bam_scope_n} tumor BAM 沒有 "
            f"<code>@RG</code>，且有 {bam_cohort_scope['cross_directory_pairing_n']} 筆 "
            "BAM↔caller cross-directory review；因此不可宣稱 full-content、完整生物樣本身分或 platform effect。"
        )
    canonical_hcc = next(item for item in cohort if item["sample"] == "HCC1395")
    technical_hcc = next(
        item for item in cohort if item["sample"] == "HCC1395_DORADO"
    )
    canonical_biological = [
        item for item in cohort if item["technical_replicate"] == "no"
    ]
    extreme_observed = min(
        canonical_biological,
        key=lambda item: (item["tree_coverage"], item["unique_among_tree"]),
    )
    if extreme_observed["sample"] != "H2009" or extreme_observed != min(
        canonical_biological,
        key=lambda item: (item["unique_among_tree"], item["tree_coverage"]),
    ):
        raise ContractError("four-layer extreme selection no longer resolves to H2009")
    four_layer_strip = f"""
<style>
  :root {{ color-scheme: light dark; }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; font-family: "Noto Sans TC", "PingFang TC", "Microsoft JhengHei", sans-serif; }}
  .anchor {{ padding: 16px; border: 1px solid color-mix(in srgb, CanvasText 14%, transparent);
    border-radius: 14px; background: color-mix(in srgb, Canvas 97%, #0f766e 3%); }}
  .head {{ display: flex; flex-wrap: wrap; gap: 8px 12px; align-items: baseline; margin-bottom: 10px; }}
  .eyebrow {{ color: #0f766e; font-size: 11px; font-weight: 750; letter-spacing: .1em; }}
  h2 {{ margin: 0; font-family: "Iowan Old Style", "Palatino Linotype", serif;
    font-size: 18px; line-height: 1.25; letter-spacing: .01em; }}
  .status {{ margin-left: auto; padding: 2px 8px; border: 1px solid #b45309; border-radius: 999px;
    color: #92400e; font-size: 11px; font-weight: 750; }}
  .grid {{ display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 8px; }}
  article {{ min-width: 0; padding: 11px; border-radius: 10px;
    background: color-mix(in srgb, CanvasText 5%, Canvas); }}
  h3 {{ margin: 0 0 5px; font-size: 12px; color: #0f766e; }}
  .value {{ margin: 0; font-size: 14px; font-weight: 700; font-variant-numeric: tabular-nums; }}
  .rule {{ margin: 5px 0 0; color: color-mix(in srgb, CanvasText 65%, transparent);
    font-size: 10.5px; line-height: 1.35; }}
  @media (max-width: 720px) {{ .grid {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }} }}
  @media (max-width: 420px) {{ .anchor {{ padding: 12px; }} .grid {{ grid-template-columns: 1fr; }}
    .status {{ margin-left: 0; }} }}
</style>
<section class="anchor" lang="zh-Hant" aria-labelledby="four-layer-title">
  <header class="head">
    <span class="eyebrow">ALL-SCOPE EVIDENCE ANCHOR</span>
    <h2 id="four-layer-title">見林、見樹、見極端，也看技術對照</h2>
    <span class="status">PARTIAL · 描述性</span>
  </header>
  <div class="grid">
    <article>
      <h3>Aggregate · n={cohort_scope['biological_macro_n']}</h3>
      <p class="value">Tree {cohort_scope['biological_macro_tree_coverage_median']:.1%} · Unique/tree {cohort_scope['biological_macro_unique_among_tree_median']:.1%}</p>
      <p class="rule">6 biological rows 的 macro median；IQR 分別 {cohort_scope['biological_macro_tree_coverage_q1']:.1%}–{cohort_scope['biological_macro_tree_coverage_q3']:.1%}、{cohort_scope['biological_macro_unique_among_tree_q1']:.1%}–{cohort_scope['biological_macro_unique_among_tree_q3']:.1%}。</p>
    </article>
    <article>
      <h3>Canonical · HCC1395</h3>
      <p class="value">{canonical_hcc['regions']:,} regions · Tree {canonical_hcc['tree_coverage']:.1%} · Unique/tree {canonical_hcc['unique_among_tree']:.1%}</p>
      <p class="rule">唯一具有 v1/v3 drilldown 的主資料集；可觀察但 bundle 仍 BLOCKED。</p>
    </article>
    <article>
      <h3>Extreme observed · {extreme_observed['sample']}</h3>
      <p class="value">{extreme_observed['regions']:,} regions · Tree {extreme_observed['tree_coverage']:.1%} · Unique/tree {extreme_observed['unique_among_tree']:.1%}</p>
      <p class="rule">在 6 個 biological rows 同時具有最低兩比例；稱 observed extreme，不稱統計 outlier。</p>
    </article>
    <article>
      <h3>Technical control · HCC1395_DORADO</h3>
      <p class="value">Δ Tree {technical_hcc['tree_coverage'] - canonical_hcc['tree_coverage']:+.1%} · Δ Unique/tree {technical_hcc['unique_among_tree'] - canonical_hcc['unique_among_tree']:+.1%}</p>
      <p class="rule">paired biological_id=HCC1395；排除於 biological macro，差值不宣稱 platform effect。</p>
    </article>
  </div>
</section>
""".strip()
    def opportunity_bar_rows(field: str) -> str:
        maximum = max(item[field] for item in cohort)
        rows: list[str] = []
        for item in cohort:
            value = item[field]
            width = 100 * value / maximum
            sample_label = item["sample"] + (
                " · technical" if item["technical_replicate"] == "yes" else ""
            )
            tone = " technical" if item["technical_replicate"] == "yes" else ""
            rows.append(
                f'<div class="bar-row{tone}" role="img" '
                f'aria-label="{html.escape(sample_label)}: {value:,}">'
                f'<span class="sample">{html.escape(sample_label)}</span>'
                f'<span class="track"><span class="fill" style="width:{width:.5f}%"></span></span>'
                f'<strong>{value:,}</strong></div>'
            )
        return "".join(rows)

    opportunity_comparison = f"""
<style>
  :root {{ color-scheme: light dark; }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; font-family: "Noto Sans TC", "PingFang TC", "Microsoft JhengHei", sans-serif; }}
  .panel {{ padding: 16px; border: 1px solid color-mix(in srgb, CanvasText 14%, transparent);
    border-radius: 14px; background: color-mix(in srgb, Canvas 97%, #0f766e 3%); }}
  header {{ display: flex; flex-wrap: wrap; align-items: baseline; gap: 7px 12px; margin-bottom: 12px; }}
  .eyebrow {{ color: #0f766e; font-size: 11px; font-weight: 750; letter-spacing: .1em; }}
  h2 {{ margin: 0; font: 650 18px/1.25 "Iowan Old Style", "Palatino Linotype", serif; }}
  .note {{ margin: 0; color: color-mix(in srgb, CanvasText 62%, transparent); font-size: 11px; }}
  .grid {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 20px; }}
  h3 {{ margin: 0 0 9px; font-size: 13px; }}
  .bar-row {{ display: grid; grid-template-columns: minmax(112px, 1.25fr) minmax(80px, 2fr) auto;
    gap: 8px; align-items: center; min-height: 27px; }}
  .sample {{ overflow-wrap: anywhere; font-size: 11px; font-weight: 650; }}
  .track {{ height: 10px; overflow: hidden; border-radius: 999px;
    background: color-mix(in srgb, CanvasText 9%, Canvas); }}
  .fill {{ display: block; height: 100%; border-radius: inherit; background: #0f766e; }}
  .technical .fill {{ background: #b7791f; }}
  strong {{ min-width: 58px; font-size: 11px; font-variant-numeric: tabular-nums; text-align: right; }}
  @media (max-width: 760px) {{ .grid {{ grid-template-columns: 1fr; }} }}
  @media (max-width: 420px) {{ .panel {{ padding: 12px; }}
    .bar-row {{ grid-template-columns: minmax(105px, 1.25fr) minmax(54px, 1fr) auto; }} }}
</style>
<section class="panel" lang="zh-Hant" aria-labelledby="opportunity-fixed-title">
  <header>
    <span class="eyebrow">ALL 7 DATASETS · FIXED LABELS</span>
    <h2 id="opportunity-fixed-title">機會量比較：名稱、bar 與精確值同列</h2>
    <p class="note">各 panel 獨立縮放；長 bar 不代表較高品質。</p>
  </header>
  <div class="grid">
    <section aria-labelledby="region-bars-title">
      <h3 id="region-bars-title">Analysis regions</h3>
      {opportunity_bar_rows('regions')}
    </section>
    <section aria-labelledby="ssnv-bars-title">
      <h3 id="ssnv-bars-title">Distinct sSNV</h3>
      {opportunity_bar_rows('distinct_ssnv')}
    </section>
  </div>
</section>
""".strip()
    folded_detail_guide = f"""
<style>
  :root {{ color-scheme: light dark; }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; font-family: "Noto Sans TC", "PingFang TC", "Microsoft JhengHei", sans-serif; }}
  .guide {{ padding: 18px; border: 1px solid color-mix(in srgb, CanvasText 14%, transparent);
    border-radius: 14px; background: color-mix(in srgb, Canvas 96%, #0f766e 4%); }}
  .eyebrow {{ color: #0f766e; font-size: 11px; font-weight: 750; letter-spacing: .12em; }}
  h2 {{ margin: 4px 0 6px; font-family: "Iowan Old Style", "Palatino Linotype", serif;
    font-size: 19px; line-height: 1.25; letter-spacing: .01em; }}
  .lede {{ margin: 0 0 12px; color: color-mix(in srgb, CanvasText 70%, transparent); }}
  details {{ margin-top: 8px; border-top: 1px solid color-mix(in srgb, CanvasText 12%, transparent); }}
  summary {{ display: flex; gap: 10px; align-items: baseline; padding: 10px 2px; cursor: pointer;
    font-weight: 650; list-style: none; }}
  summary::-webkit-details-marker {{ display: none; }}
  summary::after {{ content: '+'; margin-left: auto; color: #0f766e; font-size: 18px; }}
  details[open] summary::after {{ content: '−'; }}
  summary:focus-visible {{ outline: 2px solid #0f766e; outline-offset: 3px; border-radius: 4px; }}
  .index {{ color: #0f766e; font-variant-numeric: tabular-nums; }}
  .detail {{ padding: 0 2px 12px 32px; color: color-mix(in srgb, CanvasText 78%, transparent); }}
  .detail p {{ margin: 0 0 7px; }}
  .detail code {{ font: 12px ui-monospace, SFMono-Regular, Menlo, monospace; }}
  .source {{ color: color-mix(in srgb, CanvasText 55%, transparent); font-size: 11px; }}
  @media (max-width: 520px) {{ .guide {{ padding: 14px; }} .detail {{ padding-left: 0; }} }}
</style>
<section class="guide" lang="zh-Hant" aria-labelledby="fold-guide-title">
  <div class="eyebrow">DETAILS ON DEMAND</div>
  <h2 id="fold-guide-title">摺疊式細節與證據界線</h2>
  <p class="lede">主畫面先顯示決策級訊號；需要分母、口徑或缺失欄位時再展開。</p>
  <details>
    <summary><span class="index">01</span><span>預設 All 範圍目前代表什麼？</span></summary>
    <div class="detail">
      <p>7 datasets = 6 biological + {cohort_scope['technical_replicate_n']} technical replicate；hash、identity、receipt 三個技術 gate 都是 7/7 PASS。全 family 完整 {cohort_scope['all_families_complete_n']}/{cohort_scope['all_families_complete_total_n']}；全 objective-certified {cohort_scope['all_objective_certified_n']}/{cohort_scope['all_objective_certified_total_n']}。排除技術重複後的 macro tree coverage 為 {cohort_scope['biological_macro_tree_coverage_median']:.1%}（IQR {cohort_scope['biological_macro_tree_coverage_q1']:.1%}–{cohort_scope['biological_macro_tree_coverage_q3']:.1%}），macro unique/tree 為 {cohort_scope['biological_macro_unique_among_tree_median']:.1%}（IQR {cohort_scope['biological_macro_unique_among_tree_q1']:.1%}–{cohort_scope['biological_macro_unique_among_tree_q3']:.1%}）。</p>
      <p class="source">來源：scope_summary；切換資料集後，以頁面上方動態卡片為準。</p>
    </div>
  </details>
  <details>
    <summary><span class="index">02</span><span>比例如何計算，哪些不能互相比？</span></summary>
    <div class="detail">
      <p><code>tree coverage = tree_n / region_n</code>；<code>unique/tree = unique_tree_n / tree_n</code>；<code>ISM linkage = linked sSNV / topology distinct sSNV</code>。raw p 與 family 內 BH 需分開閱讀；不同甲基化軸的 effect unit 不可直接比較。</p>
      <p class="source">來源：topology_rates、bundle_overview、axis_rates。</p>
    </div>
  </details>
  <details>
    <summary><span class="index">03</span><span>目前仍缺哪些 multi-BAM QC 與真值資料？</span></summary>
    <div class="detail">
      <p>已收集 tumor／paired-normal quickcheck、bounded payload chunks、index SHA、reference dictionary 與既有 LongPhase-S alignment-tag receipt。仍缺 mapping rate、depth／IQR、read N50、primary-read HP／PS、MM／ML 與 CpG 完整度、KDE provenance，以及 truth VCF／high-confidence BED／som.py 指標；缺值維持 null／UNAVAILABLE，絕不以 0 代替。</p>
      <p class="source">來源：availability；詳細欄位契約見同目錄的 multi_bam_dashboard_metric_contract.md。</p>
    </div>
  </details>
  <details>
    <summary><span class="index">04</span><span>{bam_identity_detail_title}</span></summary>
    <div class="detail">
      <p>{bam_identity_detail_body}</p>
      <p class="source">來源：bam_scope_summary、bam_input_identity；full BAM SHA-256 尚未執行。</p>
    </div>
  </details>
</section>
""".strip()
    blocks = [
        {
            "id": "b_summary_header",
            "type": "markdown",
            "body": (
                "## 總覽 — 已接入的 topology／bounded-input gates 完整，下游科學證據仍為 partial\n\n"
                "預設 **All** 是 7 datasets＝6 biological＋1 technical replicate；切換 **選擇資料集** 後，"
                "指標、圖表與表格會同步隔離。已接入 topology、bounded tumor／normal／reference 與既有 tag receipt；"
                "ISM／lineage／truth 等未完成構面維持 unavailable，claim ceiling 仍是描述性觀察與資料產品 QA。"
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
        {"id": "b_topology_chart", "type": "chart", "chartId": "ch_topology_rates", "layout": "full"},
        {
            "id": "b_opportunity_header",
            "type": "markdown",
            "body": "## 機會量 — 先看每組資料能被分析多少\n\nRegion 與 distinct sSNV 分開成同層小圖；避免雙軸掩蓋尺度差異，也不把較大的 opportunity 誤稱為較好的樣本。",
        },
        {"id": "b_region_counts", "type": "chart", "chartId": "ch_region_counts", "layout": "half"},
        {"id": "b_ssnv_counts", "type": "chart", "chartId": "ch_ssnv_counts", "layout": "half"},
        {"id": "b_four_layer", "type": "html", "body": four_layer_strip, "layout": "full"},
        {"id": "b_opportunity_fixed_labels", "type": "html", "body": opportunity_comparison, "layout": "full"},
        {
            "id": "b_bam_header",
            "type": "markdown",
            "body": "## 輸入與 phasing readiness — bounded identity，非 full-content attestation\n\n下方兩張卡與 tag 圖會隨資料集選擇器同步更新。Payload MATCH只代表frozen chunks／index／metadata；strict MATCH與MOUNT_DEVICE_DRIFT_ONLY依scope動態重算，缺少@RG、caller配對檢查與full BAM SHA缺口都保留在明細。",
        },
        {
            "id": "b_bam_metrics",
            "type": "metric-strip",
            "cardIds": ["c_bam_payload_identity", "c_bam_input_readiness"],
        },
        {"id": "b_read_tag_chart", "type": "chart", "chartId": "ch_read_tag_rates", "layout": "full"},
        {"id": "b_availability", "type": "table", "tableId": "t_availability", "layout": "full"},
        {
            "id": "b_diagnostics_header",
            "type": "markdown",
            "body": "## 診斷 — HCC1395 下游選擇與多重比較\n\n保留 linked／denominator 計數，並分開呈現 raw-p discovery 與 family 內 BH correction。v1/v3 視窗不同與 cluster 軸循環性仍是明示 gate。",
        },
        {"id": "b_coverage_chart", "type": "chart", "chartId": "ch_coverage_by_k", "layout": "half"},
        {"id": "b_axis_chart", "type": "chart", "chartId": "ch_axis_rates", "layout": "half"},
        {"id": "b_coverage_denominator_rail", "type": "table", "tableId": "t_coverage_denominator_rail", "layout": "half"},
        {"id": "b_axis_denominator_rail", "type": "table", "tableId": "t_axis_denominator_rail", "layout": "full"},
        {"id": "b_bundle_table", "type": "table", "tableId": "t_bundle_overview", "layout": "full"},
        {
            "id": "b_details_header",
            "type": "markdown",
            "body": "## 明細 — 精確分母、gate 與 QA 收據\n\n下方是稽核介面：保留技術重複、分析範圍、effect unit、因果 gate 與瀏覽器收據，不壓縮成 pooled cohort claim。",
        },
        {"id": "b_folded_detail_guide", "type": "html", "body": folded_detail_guide, "layout": "full"},
        {"id": "b_bam_input_table", "type": "table", "tableId": "t_bam_input_identity", "layout": "full"},
        {"id": "b_tag_table", "type": "table", "tableId": "t_producer_tag_details", "layout": "full"},
        {"id": "b_topology_table", "type": "table", "tableId": "t_topology_details", "layout": "full"},
        {"id": "b_axis_table", "type": "table", "tableId": "t_axis_details", "layout": "full"},
        {"id": "b_lca_table", "type": "table", "tableId": "t_lca_summary", "layout": "full"},
        {"id": "b_visual_table", "type": "table", "tableId": "t_visual_qa", "layout": "full"},
        {
            "id": "b_method",
            "type": "markdown",
            "body": "## 解讀 guardrails\n\n- Cohort topology 維持每個技術資料集一列，不執行 pooled-locus inference。\n- HCC1395_DORADO 為 reproducibility 保留，但不是第 7 個生物樣本；BAM↔caller 跨目錄差異不宣稱 platform effect。\n- 頁面顯示為百分比的欄位，在 snapshot 中皆以 0–1 儲存，並保留 exact numerator／denominator。\n- Bounded payload MATCH 是 frozen chunks／index／metadata 證據，不是 full BAM content SHA；strict identity category依每個input role重算，不把bounded status升格為full-content PASS。\n- LongPhase-S tag 圖使用 all captured mapped alignment records（含 primary／secondary／supplementary），不是 primary-read 或 unique-read rate。\n- 缺少的 mapping rate、depth、read N50、primary-read HP/PS、MM/ML completeness、KDE provenance 與 truth benchmark 維持 null，並附明確 access status。\n- ISM 可用性會隨 active-k 選擇改變；untested 不等於 nonsignificant。\n- Lineage/LCA 在受控輸入身分固定前只能作描述性觀察。\n- Browser 截圖只驗證呈現，不驗證生物真值。",
        },
    ]

    artifact = {
        "surface": "dashboard",
        "manifest": {
            "version": 1,
            "surface": "dashboard",
            "title": "InterSubMod 多樣本分析總覽",
            "description": "PARTIAL snapshot：7 組 topology/BAM 資料＝6 個生物樣本＋1 技術 replicate；整合 bounded tumor-normal-reference intake、LongPhase-S alignment tags、HCC1395 下游診斷、缺失 gate、可選資料集與可追溯明細。",
            "generatedAt": generated_at,
            "filters": [
                {
                    "id": "sample_filter",
                    "label": "選擇資料集",
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
                "bam_scope_summary": bam_scope_summary,
                "availability": availability,
                "topology_rates": topology_rates,
                "region_counts": scoped_rows(region_counts),
                "ssnv_counts": scoped_rows(ssnv_counts),
                "read_tag_rates": scoped_rows(read_tag_rates),
                "bam_input_identity": scoped_rows(bam_input_identity),
                "producer_tag_details": scoped_rows(producer_tag_details),
                "coverage_by_k": hcc_scoped_rows(coverage),
                "coverage_denominator_rail": hcc_scoped_rows(
                    coverage_denominator_rail
                ),
                "axis_rates": hcc_scoped_rows(axis_rates_base),
                "axis_denominator_rail": hcc_scoped_rows(axis_denominator_rail),
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
                    "message": "僅 HCC1395 有 ISM 觀察；v1 source 缺失，v3 尚未 content-addressed。",
                },
                {
                    "id": "multisample_lineage_unavailable",
                    "scope": "G4 multi-sample lineage",
                    "sourceId": "q_lca_summary",
                    "message": "僅 HCC1395 有 lineage/LCA；輸入 BAM／threads 未固定，受控 A/B gate FAIL。",
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
        "technical_replicate_n": 1,
        "hash_match_n": 7,
        "sample_identity_match_n": 7,
        "receipt_technical_pass_n": 7,
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

    bam_scope_rows = snapshot["datasets"]["bam_scope_summary"]
    if len(bam_scope_rows) != EXPECTED_DATASET_N + 1:
        raise ContractError("BAM scope summary must contain All plus seven dataset rows")
    bam_scope_all = next(
        row for row in bam_scope_rows if row["sample_filter"] == ALL_DATASETS_SCOPE
    )
    all_bam_rows = [
        row for row in snapshot["datasets"]["bam_input_identity"]
        if row["sample_filter"] == ALL_DATASETS_SCOPE
    ]
    if len(all_bam_rows) != EXPECTED_DATASET_N:
        raise ContractError("All-scope BAM identity table must contain seven rows")
    strict_match_n = 0
    mount_drift_n = 0
    for row in all_bam_rows:
        role_statuses = (
            row["tumor_strict_identity"],
            row["normal_strict_identity"],
            row["reference_strict_identity"],
        )
        if all(status == "MATCH" for status in role_statuses):
            strict_match_n += 1
        elif (
            all(status in {"MATCH", "MOUNT_DEVICE_DRIFT_ONLY"} for status in role_statuses)
            and any(status == "MOUNT_DEVICE_DRIFT_ONLY" for status in role_statuses)
        ):
            mount_drift_n += 1
        else:
            raise ContractError(
                f"{row['sample']}: unsupported strict storage identity combination {role_statuses}"
            )
    expected_bam_counts = {
        "dataset_count": 7,
        "tumor_quickcheck_pass_n": 7,
        "normal_quickcheck_pass_n": 7,
        "tumor_bounded_payload_match_n": 7,
        "normal_bounded_payload_match_n": 7,
        "reference_compatible_n": 7,
        "strict_storage_identity_match_n": strict_match_n,
        "mount_device_drift_n": mount_drift_n,
        "tag_receipt_pass_n": 7,
        "full_bam_sha256_n": 0,
        "no_read_group_n": 7,
        "cross_directory_pairing_n": 1,
    }
    for field, expected in expected_bam_counts.items():
        if bam_scope_all[field] != expected:
            raise ContractError(
                f"All-scope BAM {field} must be {expected}, got {bam_scope_all[field]}"
            )

    expected_scoped_row_counts = {
        "region_counts": EXPECTED_DATASET_N * 2,
        "ssnv_counts": EXPECTED_DATASET_N * 2,
        "bam_input_identity": EXPECTED_DATASET_N * 2,
        "producer_tag_details": EXPECTED_DATASET_N * 2,
        "read_tag_rates": EXPECTED_DATASET_N * 6,
    }
    for dataset, expected in expected_scoped_row_counts.items():
        observed = len(snapshot["datasets"][dataset])
        if observed != expected:
            raise ContractError(
                f"{dataset} must contain {expected} All+sample scoped rows, got {observed}"
            )

    for row in all_bam_rows:
        if any(
            row[field] != expected
            for field, expected in (
                ("tumor_quickcheck", "PASS"),
                ("normal_quickcheck", "PASS"),
                ("tumor_bounded_payload", "MATCH"),
                ("normal_bounded_payload", "MATCH"),
                ("reference_dictionary", "PASS"),
                ("read_group_identity", "UNAVAILABLE_NO_RG"),
                ("full_bam_sha256", "NOT_COLLECTED"),
            )
        ):
            raise ContractError(f"{row['sample']}: BAM identity detail contract changed")
        expected_pairing = (
            "CROSS_DIRECTORY_REVIEW_REQUIRED"
            if row["sample"] == "HCC1395"
            else "SAME_DIRECTORY_FAMILY"
        )
        if row["source_directory_pairing"] != expected_pairing:
            raise ContractError(f"{row['sample']}: BAM/caller pairing gate changed")

    all_tag_rows = [
        row for row in snapshot["datasets"]["producer_tag_details"]
        if row["sample_filter"] == ALL_DATASETS_SCOPE
    ]
    if len(all_tag_rows) != EXPECTED_DATASET_N:
        raise ContractError("All-scope producer tag detail must contain seven rows")
    for row in all_tag_rows:
        if (
            row["primary_alignments"]
            + row["secondary_alignments"]
            + row["supplementary_alignments"]
            != row["total_alignments"]
            or row["record_key_missing"] != 0
            or row["record_key_extra"] != 0
            or row["identity_conflicts"] != 0
            or row["unknown_hp_category_n"] != 0
        ):
            raise ContractError(f"{row['sample']}: producer tag conservation changed")
        assert_close(
            row["hp_rate_all_alignments"],
            ratio(row["hp_alignments"], row["total_alignments"], f"{row['sample']}.HP/all"),
            f"{row['sample']}.HP/all",
        )
        assert_close(
            row["hp_ps_rate_among_hp"],
            ratio(row["hp_ps_alignments"], row["hp_alignments"], f"{row['sample']}.HP+PS/HP"),
            f"{row['sample']}.HP+PS/HP",
        )
        assert_close(
            row["duplicate_identity_rate_all_alignments"],
            ratio(
                row["duplicate_identity_rows"],
                row["total_alignments"],
                f"{row['sample']}.duplicate/all",
            ),
            f"{row['sample']}.duplicate/all",
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
        scope_all["biological_macro_tree_coverage_q1"],
        percentile_cont(
            (row["tree_coverage"] for row in biological_rows),
            0.25,
            "validation tree coverage P25",
        ),
        "biological macro tree coverage P25",
    )
    assert_close(
        scope_all["biological_macro_tree_coverage_q3"],
        percentile_cont(
            (row["tree_coverage"] for row in biological_rows),
            0.75,
            "validation tree coverage P75",
        ),
        "biological macro tree coverage P75",
    )
    assert_close(
        scope_all["biological_macro_unique_among_tree_median"],
        median(row["unique_among_tree"] for row in biological_rows),
        "biological macro unique-among-tree median",
    )
    assert_close(
        scope_all["biological_macro_unique_among_tree_q1"],
        percentile_cont(
            (row["unique_among_tree"] for row in biological_rows),
            0.25,
            "validation unique/tree P25",
        ),
        "biological macro unique-among-tree P25",
    )
    assert_close(
        scope_all["biological_macro_unique_among_tree_q3"],
        percentile_cont(
            (row["unique_among_tree"] for row in biological_rows),
            0.75,
            "validation unique/tree P75",
        ),
        "biological macro unique-among-tree P75",
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


def validate_io_path_separation(args: argparse.Namespace) -> None:
    """Prevent the atomic dashboard output from replacing any explicit or implicit input."""

    inputs = {
        "bam_manifest": args.bam_manifest.resolve(),
        "bam_manifest_receipt": args.bam_manifest_receipt.resolve(),
        "bam_manifest_schema": DEFAULT_BAM_MANIFEST_SCHEMA.resolve(),
        "bam_manifest_builder": BAM_MANIFEST_BUILDER.resolve(),
        "source_manifest_schema": DEFAULT_SOURCE_SCHEMA.resolve(),
        "dashboard_builder": Path(__file__).resolve(),
        **{
            f"results_{key}": (args.results_dir / filename).resolve()
            for key, filename in INPUT_FILES.items()
        },
    }
    output = args.output.resolve()
    collisions = [name for name, input_path in inputs.items() if output == input_path]
    if collisions:
        raise ContractError(f"dashboard output collides with input path(s): {collisions}")


def main() -> int:
    args = parse_args()
    validate_io_path_separation(args)
    output = args.output.resolve()
    artifact = build_artifact(
        args.results_dir.resolve(),
        args.bam_manifest.resolve(),
        args.bam_manifest_receipt.resolve(),
    )
    atomic_write_json(output, artifact)
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
