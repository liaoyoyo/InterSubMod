#!/usr/bin/env python3
"""Build the canonical seven-dataset layered reconstruction workstation.

The repository machine summary is the only authority for sample paths and all
landing-page measurements.  The command fails closed before writing index.html
when that authority, any bound artifact, or any sample-page provenance marker
is stale or inconsistent.

Default mode delegates the seven large sample pages to
``build_layered_workstation_v5.py`` and then writes the cohort index.  ``--index-only``
never rebuilds sample pages; it succeeds only when all seven existing pages are
already bound to the current summary and region-view hashes.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import subprocess
import sys
from datetime import datetime
from html.parser import HTMLParser
from pathlib import Path
from string import Template
from typing import Any


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
BUILD = HERE / "build_layered_workstation_v5.py"
OUTDIR = (HERE / ".." / ".." / "layered_workstation").resolve()
CURRENT_SUMMARY = (
    REPO_ROOT
    / "research"
    / "20260710_layered_reconstruction_v2"
    / "current_layered_topology_v3_raw_all_v1.json"
)
READ_AF_TOPOLOGY_DIR = (
    REPO_ROOT
    / "research"
    / "20260715_layered_workstation_genome_topology_multiselect"
    / "data"
    / "current_v5_read_af_topology"
)
READ_AF_TOPOLOGY_INDEX = READ_AF_TOPOLOGY_DIR / "current_v5_read_af_topology.index.json"
READ_AF_TOPOLOGY_BUILD = (
    REPO_ROOT
    / "research"
    / "20260715_layered_workstation_genome_topology_multiselect"
    / "scripts"
    / "build_current_v5_read_af_topology.py"
)
SAMPLE_COMPARISON_DIR = (
    REPO_ROOT
    / "research"
    / "20260715_sample_topology_comparison"
)
SAMPLE_COMPARISON_JSON = SAMPLE_COMPARISON_DIR / "artifacts" / "sample_topology_comparison.json"
SAMPLE_COMPARISON_RECEIPT = SAMPLE_COMPARISON_DIR / "artifacts" / "validation_receipt.json"
SAMPLE_COMPARISON_BUILD = (
    SAMPLE_COMPARISON_DIR / "scripts" / "analyze_sample_topology_comparison.py"
)
HCC_EXACT_SIGNATURE_JSON = SAMPLE_COMPARISON_DIR / "artifacts" / "hcc1395_exact_signature_validation.json"
HCC_EXACT_SIGNATURE_BUILD = SAMPLE_COMPARISON_DIR / "scripts" / "analyze_hcc_exact_signature.py"
PYTHON = sys.executable or "python3"

EXPECTED_SAMPLES = {
    "COLO829",
    "H1437",
    "H2009",
    "HCC1395",
    "HCC1395_DORADO",
    "HCC1937",
    "HCC1954",
}
EXPECTED_SCOPE = "7 datasets / 6 biological samples / chr1-22"
EXPECTED_BACKBONE = "longphase_s_recalibrated_FILTER_PASS"
EXPECTED_TREE_SOURCE = "longphase_s_recalibrated_filter_pass"
EXPECTED_REGION_SCOPE = "chr1-22 primary; chrX/chrY out-of-scope census only"
EXPECTED_UI_CONTRACT = "layered-workstation-v5-grch38-topology-multiselect-2"
PAGE_META_NAMES = {
    "summary_sha256": "intersubmod-current-summary-sha256",
    "region_sha256": "intersubmod-region-view-sha256",
    "read_af_sha256": "intersubmod-read-af-topology-sha256",
    "sample": "intersubmod-canonical-sample",
    "renderer_sha256": "intersubmod-renderer-sha256",
    "ui_contract": "intersubmod-ui-contract",
}
TOPOLOGY_CLASSES = (
    ("exact_and_topology_unique", "C=1 / Topo=1", "exact 與拓撲皆唯一", "unique"),
    ("topology_unique_exact_multiple", "C>1 / Topo=1", "exact 多解、拓撲唯一", "shape"),
    ("topology_multiple_exact_multiple", "C>1 / Topo>1", "exact 與拓撲皆多解", "multiple"),
    ("incomplete", "Incomplete", "候選集合未完整", "incomplete"),
)
READ_AF_SELECTION_CLASSES = (
    (
        "structural_exact_unique",
        "結構已 exact 唯一",
        "C=1 / Topo=1；不需 read-AF 排序",
        "read-structural",
    ),
    (
        "read_af_unique_first",
        "read-AF 唯一第一順位",
        "Exact Fraction score 只有一個第一順位 candidate",
        "read-unique",
    ),
    (
        "read_af_tied_same_topology",
        "並列第一 · 同一 Topo",
        "Exact candidates 並列，但 canonical shape 一致",
        "read-same",
    ),
    (
        "read_af_tied_different_topology",
        "並列第一 · 多 Topo",
        "Exact 第一順位集合仍跨多個 canonical shapes",
        "read-multiple",
    ),
    (
        "read_af_unavailable",
        "read-AF 不可排序",
        "Coverage／recurrence 等條件不足，維持不可評估",
        "read-unavailable",
    ),
    (
        "incomplete",
        "候選未完整",
        "不可在 capped／incomplete candidate prefix 上排序",
        "incomplete",
    ),
)
MORPHOLOGY_CLASSES = (
    (
        "single_no_within_hp_relation",
        "單支／無 HP 內關係",
        "各 primary HP 均未見 depth≥2 或 outdegree≥2",
        "morph-single",
    ),
    (
        "direct_chain",
        "直系鏈",
        "至少一個 primary HP depth≥2；subclone-compatible",
        "morph-direct",
    ),
    (
        "sister_branch",
        "旁系／分支",
        "至少一個 primary HP outdegree≥2；branch-compatible",
        "morph-sister",
    ),
    (
        "direct_and_sister",
        "直系＋旁系",
        "Nesting 與 branching 並存；可能分屬不同 HP",
        "morph-mixed",
    ),
    (
        "unresolved",
        "形態未解",
        "候選未完整、ranking 不可用或第一順位仍跨多 shape",
        "incomplete",
    ),
)
COMPARISON_DIMENSIONS = {
    "structural": {
        "label": "結構可辨識度",
        "short": "Structural",
        "classes": TOPOLOGY_CLASSES,
    },
    "read_af_selection": {
        "label": "read-AF 第一順位",
        "short": "Read-AF",
        "classes": READ_AF_SELECTION_CLASSES,
    },
    "morphology": {
        "label": "拓撲形態",
        "short": "Morphology",
        "classes": MORPHOLOGY_CLASSES,
    },
}


class AuthorityError(RuntimeError):
    """Raised when canonical provenance cannot be proven exactly."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuthorityError(message)


def escaped(value: Any) -> str:
    return html.escape(str(value), quote=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AuthorityError(f"cannot read canonical JSON {path}: {exc}") from exc
    require(isinstance(payload, dict), f"canonical JSON is not an object: {path}")
    return payload


def parse_timestamp(value: str, label: str) -> datetime:
    try:
        timestamp = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (TypeError, ValueError) as exc:
        raise AuthorityError(f"invalid {label} timestamp: {value!r}") from exc
    require(timestamp.tzinfo is not None, f"{label} timestamp has no timezone")
    return timestamp


def resolve_authority_path(value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (REPO_ROOT / path).resolve()


def require_under(path: Path, parent: Path, label: str) -> None:
    try:
        path.relative_to(parent)
    except ValueError as exc:
        raise AuthorityError(f"{label} escapes canonical run root: {path}") from exc


def verify_bound_file(path: Path, expected_sha256: str, label: str) -> None:
    require(path.is_file(), f"missing {label}: {path}")
    actual = sha256_file(path)
    require(actual == expected_sha256, f"{label} hash drift: expected={expected_sha256} actual={actual}")


def validate_aggregate(aggregate: dict[str, Any], samples: list[dict[str, Any]]) -> None:
    require(aggregate.get("dataset_count") == 7, "canonical aggregate dataset_count must be 7")
    require(aggregate.get("all_invariants_pass") is True, "canonical aggregate invariants are not all PASS")

    additive = (
        "tree_input_records",
        "autosomal_biallelic_sSNV",
        "retained_sSNV",
        "W_tree",
        "W_primary",
        "no_primary_lineage",
        "primary_units",
        "primary_units_full_and_partial",
        "primary_units_partial_only",
        "primary_units_full_only",
        "complete_regions",
        "incomplete_regions",
        "read_tag_exposures",
        "read_tag_exact_matches",
        "mixed_PS_regions",
    )
    for field in additive:
        observed = sum(int(sample[field]) for sample in samples)
        require(observed == int(aggregate[field]), f"canonical aggregate mismatch for {field}")

    topology = aggregate["topology_classes"]
    for key, _, _, _ in TOPOLOGY_CLASSES:
        observed = sum(int(sample["topology_classes"][key]) for sample in samples)
        require(observed == int(topology[key]), f"canonical topology aggregate mismatch for {key}")
    require(int(topology.get("impossible_exact_unique_topology_multiple", -1)) == 0,
            "impossible C=1 / Topo>1 state is non-zero")
    require(int(aggregate["W_tree"]) == int(aggregate["W_primary"]) + int(aggregate["no_primary_lineage"]),
            "W_tree does not conserve W_primary + no-primary")
    require(int(aggregate["W_primary"]) == int(aggregate["complete_regions"]) + int(aggregate["incomplete_regions"]),
            "W_primary does not conserve complete + incomplete")
    require(int(aggregate["W_primary"]) == sum(int(topology[key]) for key, _, _, _ in TOPOLOGY_CLASSES),
            "topology classes do not conserve W_primary")
    require(int(aggregate["read_tag_exposures"]) == int(aggregate["read_tag_exact_matches"]),
            "read-tag exposures do not equal exact matches")


def load_authority() -> dict[str, Any]:
    """Load and exhaustively bind the current canonical machine summary."""
    require(CURRENT_SUMMARY.is_file(), f"missing current summary: {CURRENT_SUMMARY}")
    summary_sha256 = sha256_file(CURRENT_SUMMARY)
    summary = load_json(CURRENT_SUMMARY)
    require(summary.get("schema_name") == "intersubmod.current_layered_topology_summary",
            "unexpected current-summary schema_name")
    require(summary.get("schema_version") == "1.0.0", "unexpected current-summary schema_version")
    require(summary.get("all_pass") is True, "current summary all_pass is not true")
    require(summary.get("task_type") == "B_comprehensive_validation", "current summary is not Task B")
    require(summary.get("scope") == EXPECTED_SCOPE, f"current summary scope is not {EXPECTED_SCOPE}")
    require(summary.get("claim_scope") == "regional mutation-state candidate trees; not confirmed cell clones",
            "current summary claim boundary drifted")
    definitions = summary.get("definitions", {})
    require(definitions.get("PS") == "phase-block QC context only; not a topology edge or lineage label",
            "current summary PS contract drifted")

    canonical = summary.get("canonical")
    require(isinstance(canonical, dict), "current summary has no canonical object")
    require(canonical.get("label") == EXPECTED_BACKBONE, "canonical backbone is not LongPhase-S recalibrated PASS")
    run_root = resolve_authority_path(canonical["run_root"])
    require(run_root.is_dir(), f"missing canonical run root: {run_root}")

    success_path = run_root / "_SUCCESS"
    verify_bound_file(success_path, canonical["success_sha256"], "canonical _SUCCESS")
    success = load_json(success_path)
    require(success.get("run_id") == run_root.name, "_SUCCESS run_id/root mismatch")
    require(success.get("extra", {}).get("dataset_count") == 7, "_SUCCESS dataset_count is not 7")
    require(success.get("extra", {}).get("biological_sample_count") == 6,
            "_SUCCESS biological_sample_count is not 6")
    require(success.get("extra", {}).get("scope") == "chr1-22", "_SUCCESS scope is not chr1-22")
    require(success.get("extra", {}).get("mode") == "full", "_SUCCESS mode is not full")

    verification_path = resolve_authority_path(canonical["verification_path"])
    require_under(verification_path, run_root, "verification summary")
    verify_bound_file(verification_path, canonical["verification_sha256"], "verification summary")
    verification = load_json(verification_path)
    require(verification.get("all_pass") is True, "verification all_pass is not true")
    require(verification.get("dataset_count") == 7, "verification dataset_count is not 7")
    require(verification.get("n_pass") == 7 and verification.get("n_fail") == 0,
            "verification is not 7/7 PASS")

    records = canonical.get("samples", [])
    require(isinstance(records, list) and len(records) == 7, "canonical sample list is not exactly seven records")
    names = [record.get("sample") for record in records]
    require(len(names) == len(set(names)), "canonical sample list contains duplicates")
    require(set(names) == EXPECTED_SAMPLES, f"canonical sample set drifted: {names}")
    verification_names = {record.get("sample") for record in verification.get("samples", [])}
    require(verification_names == EXPECTED_SAMPLES, "verification sample set drifted")

    require(READ_AF_TOPOLOGY_INDEX.is_file(), f"missing read-AF topology index: {READ_AF_TOPOLOGY_INDEX}")
    read_af_index = load_json(READ_AF_TOPOLOGY_INDEX)
    require(
        read_af_index.get("schema_name") == "intersubmod.current_v5_read_af_topology_index"
        and read_af_index.get("schema_version") == "1.0.0",
        "unexpected read-AF topology index schema",
    )
    require(read_af_index.get("all_checks_pass") is True, "read-AF topology index checks failed")
    require(read_af_index.get("dataset_count") == 7, "read-AF topology index dataset_count is not 7")
    read_af_provenance = read_af_index.get("provenance") or {}
    require(
        Path(read_af_provenance.get("current_summary", "")).resolve() == CURRENT_SUMMARY.resolve(),
        "read-AF topology current-summary path drifted",
    )
    require(
        read_af_provenance.get("current_summary_sha256") == summary_sha256,
        "read-AF topology current-summary hash drifted",
    )
    require(
        Path(read_af_provenance.get("run_root", "")).resolve() == run_root,
        "read-AF topology run root drifted",
    )
    read_af_manifest = Path(read_af_provenance.get("input_manifest", "")).resolve()
    require(
        read_af_manifest == (run_root / "input_manifest.snapshot.json").resolve(),
        "read-AF topology input manifest differs from canonical run manifest",
    )
    verify_bound_file(
        read_af_manifest,
        read_af_provenance.get("input_manifest_sha256", ""),
        "read-AF topology input manifest",
    )
    read_af_builder = Path(read_af_provenance.get("builder", "")).resolve()
    require(
        read_af_builder == READ_AF_TOPOLOGY_BUILD.resolve(),
        "read-AF topology builder path differs from repository builder",
    )
    verify_bound_file(
        read_af_builder,
        read_af_provenance.get("builder_sha256", ""),
        "read-AF topology builder",
    )
    read_af_solver = Path(read_af_provenance.get("solver", "")).resolve()
    require(
        read_af_solver == (HERE / "tree_enumeration_solver.py").resolve(),
        "read-AF topology solver path differs from canonical solver",
    )
    verify_bound_file(
        read_af_solver,
        read_af_provenance.get("solver_sha256", ""),
        "read-AF topology solver",
    )
    read_af_records = {row.get("sample"): row for row in read_af_index.get("samples", [])}
    require(set(read_af_records) == EXPECTED_SAMPLES, "read-AF topology sample set drifted")

    summary_generated = parse_timestamp(summary["generated_at"], "current summary generated_at")
    max_source_mtime = 0.0
    enriched_records: list[dict[str, Any]] = []
    for record in records:
        sample = record["sample"]
        require(record.get("pass") is True, f"canonical sample is not PASS: {sample}")
        require(record.get("all_invariants_pass") is True, f"sample invariants are not PASS: {sample}")
        require(record.get("tree_backbone_source") == EXPECTED_TREE_SOURCE,
                f"tree source drifted for {sample}")
        require(int(record["W_tree"]) == int(record["W_primary"]) + int(record["no_primary_lineage"]),
                f"W conservation failed for {sample}")
        require(int(record["W_primary"]) == int(record["complete_regions"]) + int(record["incomplete_regions"]),
                f"complete/incomplete conservation failed for {sample}")
        sample_topology = record["topology_classes"]
        require(int(record["W_primary"]) == sum(int(sample_topology[key]) for key, _, _, _ in TOPOLOGY_CLASSES),
                f"topology conservation failed for {sample}")
        require(int(sample_topology.get("impossible_exact_unique_topology_multiple", -1)) == 0,
                f"impossible topology state present for {sample}")

        region_path = resolve_authority_path(record["paths"]["layered_region_view"])
        layered_path = resolve_authority_path(record["paths"]["layered_reconstruction"])
        require_under(region_path, run_root, f"{sample} region view")
        require_under(layered_path, run_root, f"{sample} layered reconstruction")
        verify_bound_file(region_path, record["sha256"]["layered_region_view"], f"{sample} region view")
        verify_bound_file(layered_path, record["sha256"]["layered_reconstruction"],
                          f"{sample} layered reconstruction")
        max_source_mtime = max(max_source_mtime, region_path.stat().st_mtime, layered_path.stat().st_mtime)

        region_payload = load_json(region_path)
        require(region_payload.get("sample") == sample, f"region-view sample mismatch: {sample}")
        require(region_payload.get("schema_version") == "2.0", f"region-view schema drifted: {sample}")
        census = region_payload.get("census", {})
        require(census.get("U1_backbone_source") == EXPECTED_TREE_SOURCE,
                f"region-view backbone drifted: {sample}")
        require(census.get("analysis_scope") == EXPECTED_REGION_SCOPE,
                f"region-view scope drifted: {sample}")
        require(int(census.get("n_regions", -1)) == int(record["W_tree"]),
                f"region-view W_tree drifted: {sample}")
        require(int(census.get("U1_sSNV_somatic_total", -1)) == int(record["tree_input_records"]),
                f"region-view tree-input count drifted: {sample}")
        l3 = region_payload.get("L3_methyl", {})
        require(l3.get("status") == "not_evaluated" and l3.get("role") == "bounded_auxiliary",
                f"L3 contract drifted: {sample}")
        ps_contract = region_payload.get("analysis_contract", {}).get("PS", "")
        require(ps_contract == "preserved for phase-block QC; not used as a topology edge or lineage label",
                f"PS contract drifted: {sample}")

        enriched = dict(record)
        enriched["region_path"] = region_path
        enriched["layered_path"] = layered_path
        enriched["region_sha256"] = record["sha256"]["layered_region_view"]
        enriched["l3_status"] = l3["status"]
        enriched["canonical_record"] = record
        read_af_record = read_af_records[sample]
        require(read_af_record.get("all_checks_pass") is True, f"read-AF topology checks failed: {sample}")
        require(int(read_af_record.get("W_tree", -1)) == int(record["W_tree"]),
                f"read-AF topology W_tree drifted: {sample}")
        require(int(read_af_record.get("W_primary", -1)) == int(record["W_primary"]),
                f"read-AF topology W_primary drifted: {sample}")
        read_af_path = Path(read_af_record["output"]).resolve()
        require_under(read_af_path, READ_AF_TOPOLOGY_DIR.resolve(), f"{sample} read-AF topology")
        verify_bound_file(read_af_path, read_af_record["output_sha256"], f"{sample} read-AF topology")
        enriched["read_af_path"] = read_af_path
        enriched["read_af_sha256"] = read_af_record["output_sha256"]
        enriched["read_af_summary"] = read_af_record
        enriched_records.append(enriched)

    newest_source = datetime.fromtimestamp(max_source_mtime, tz=summary_generated.tzinfo)
    require(summary_generated >= newest_source,
            f"current summary is stale: generated={summary_generated.isoformat()} newest_source={newest_source.isoformat()}")
    validate_aggregate(canonical["aggregate"], enriched_records)

    require(SAMPLE_COMPARISON_JSON.is_file(), f"missing sample comparison: {SAMPLE_COMPARISON_JSON}")
    sample_comparison = load_json(SAMPLE_COMPARISON_JSON)
    require(
        sample_comparison.get("schema_name") == "intersubmod.sample_topology_comparison"
        and sample_comparison.get("schema_version") == "1.0.0",
        "unexpected sample-comparison schema",
    )
    require(
        sample_comparison.get("task_type") == "B_comprehensive_validation"
        and sample_comparison.get("all_checks_pass") is True,
        "sample comparison is not a passing Task B artifact",
    )
    require(
        sample_comparison.get("scope")
        == "7 datasets / 6 biological samples / GRCh38 chr1-22 / current canonical v5",
        "sample-comparison scope drifted",
    )
    sample_comparison_provenance = sample_comparison.get("provenance") or {}
    require(
        Path(sample_comparison_provenance.get("sidecar_index", "")).resolve()
        == READ_AF_TOPOLOGY_INDEX.resolve(),
        "sample comparison is bound to a different read-AF index",
    )
    require(
        sample_comparison_provenance.get("sidecar_index_sha256")
        == sha256_file(READ_AF_TOPOLOGY_INDEX),
        "sample-comparison read-AF index hash drifted",
    )
    require(
        Path(sample_comparison_provenance.get("current_summary", "")).resolve()
        == CURRENT_SUMMARY.resolve()
        and sample_comparison_provenance.get("current_summary_sha256") == summary_sha256,
        "sample comparison is bound to a different current summary",
    )
    require(
        Path(sample_comparison_provenance.get("analysis_script", "")).resolve()
        == SAMPLE_COMPARISON_BUILD.resolve(),
        "sample-comparison analysis script path drifted",
    )
    verify_bound_file(
        SAMPLE_COMPARISON_BUILD,
        sample_comparison_provenance.get("analysis_script_sha256", ""),
        "sample-comparison analysis script",
    )
    comparison_generated = parse_timestamp(sample_comparison["generated_at"], "sample comparison generated_at")
    sidecar_generated = parse_timestamp(read_af_index["generated_at"], "read-AF index generated_at")
    require(comparison_generated >= sidecar_generated, "sample comparison predates its sidecar index")

    comparison_datasets = {
        record.get("dataset"): record for record in sample_comparison.get("datasets", [])
    }
    require(set(comparison_datasets) == EXPECTED_SAMPLES, "sample-comparison dataset set drifted")
    comparison_dimension_keys = {
        "structural": "structural_classes",
        "read_af_selection": "selection_classes",
        "morphology": "morphology_classes",
    }
    for row in enriched_records:
        sample = row["sample"]
        comparison_record = comparison_datasets[sample]
        require(
            int(comparison_record.get("W_primary", -1)) == int(row["W_primary"])
            and int(comparison_record.get("W_tree", -1)) == int(row["W_tree"]),
            f"sample-comparison W drifted for {sample}",
        )
        require(
            comparison_record.get("source_sha256") == row["read_af_sha256"],
            f"sample-comparison source hash drifted for {sample}",
        )
        for dimension, summary_key in comparison_dimension_keys.items():
            expected_counts = row["read_af_summary"][summary_key]
            observed_counts = comparison_record["dimensions"][dimension]["counts"]
            require(
                observed_counts == expected_counts,
                f"sample-comparison {dimension} counts drifted for {sample}",
            )
    pairwise = sample_comparison.get("pairwise_composition") or {}
    require(len(pairwise.get("records", [])) == 21, "sample comparison does not contain 21 pairs")
    hcc_pair = sample_comparison.get("hcc1395_technical_pair") or {}
    require(
        hcc_pair.get("relationship")
        == "same biological sample; technical/cross-platform pair; not biological replication",
        "HCC1395 technical-pair relationship drifted",
    )
    primary_match_n = int(hcc_pair["exact_coordinate_overlap"]["primary_both"]["intersection"])
    require(
        all(
            int(values.get("n", -1)) == primary_match_n
            for values in hcc_pair.get("matched_primary_dimensions", {}).values()
        ),
        "HCC1395 matched confusion denominators drifted",
    )
    exact_signature = sample_comparison.get("hcc_exact_signature_evidence") or {}
    require(
        Path(exact_signature.get("source", "")).resolve() == HCC_EXACT_SIGNATURE_JSON.resolve(),
        "HCC exact-signature artifact path drifted",
    )
    verify_bound_file(
        HCC_EXACT_SIGNATURE_JSON,
        exact_signature.get("source_sha256", ""),
        "HCC exact-signature artifact",
    )
    require(
        Path(exact_signature.get("analysis_script", "")).resolve()
        == HCC_EXACT_SIGNATURE_BUILD.resolve(),
        "HCC exact-signature analysis-script path drifted",
    )
    verify_bound_file(
        HCC_EXACT_SIGNATURE_BUILD,
        exact_signature.get("analysis_script_sha256", ""),
        "HCC exact-signature analysis script",
    )
    require(exact_signature.get("all_checks_pass") is True, "HCC exact-signature checks failed")
    require(SAMPLE_COMPARISON_RECEIPT.is_file(), f"missing sample-comparison receipt: {SAMPLE_COMPARISON_RECEIPT}")
    sample_comparison_receipt = load_json(SAMPLE_COMPARISON_RECEIPT)
    require(sample_comparison_receipt.get("status") == "PASS", "sample-comparison receipt is not PASS")
    receipt_outputs = sample_comparison_receipt.get("outputs") or {}
    receipt_record = receipt_outputs.get(str(SAMPLE_COMPARISON_JSON.resolve()))
    require(
        isinstance(receipt_record, dict)
        and receipt_record.get("sha256") == sha256_file(SAMPLE_COMPARISON_JSON),
        "sample-comparison receipt hash drifted",
    )

    for role, provenance in summary.get("code_provenance", {}).items():
        source_path = resolve_authority_path(provenance["path"])
        verify_bound_file(source_path, provenance["sha256"], f"code provenance {role}")

    comparison = summary.get("backbone_comparison", {})
    comparison_path = resolve_authority_path(comparison["path"])
    verify_bound_file(comparison_path, comparison["sha256"], "backbone comparison")
    require(comparison.get("aggregate", {}).get("verdict") == "backbone_sensitive",
            "backbone comparison verdict drifted")

    return {
        "summary": summary,
        "summary_path": CURRENT_SUMMARY,
        "summary_sha256": summary_sha256,
        "summary_generated": summary_generated,
        "canonical": canonical,
        "aggregate": canonical["aggregate"],
        "samples": enriched_records,
        "run_root": run_root,
        "success": success,
        "verification": verification,
        "verification_path": verification_path,
        "read_af_index": read_af_index,
        "read_af_index_path": READ_AF_TOPOLOGY_INDEX,
        "read_af_index_sha256": sha256_file(READ_AF_TOPOLOGY_INDEX),
        "sample_comparison": sample_comparison,
        "sample_comparison_path": SAMPLE_COMPARISON_JSON,
        "sample_comparison_sha256": sha256_file(SAMPLE_COMPARISON_JSON),
        "sample_comparison_receipt_path": SAMPLE_COMPARISON_RECEIPT,
        "sample_comparison_receipt_sha256": sha256_file(SAMPLE_COMPARISON_RECEIPT),
        "hcc_exact_signature_path": HCC_EXACT_SIGNATURE_JSON,
        "hcc_exact_signature_sha256": sha256_file(HCC_EXACT_SIGNATURE_JSON),
    }


class _MetaParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.values: dict[str, str] = {}

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "meta":
            return
        values = {key.lower(): value for key, value in attrs if value is not None}
        name = values.get("name")
        if name in PAGE_META_NAMES.values():
            self.values[name] = values.get("content", "")


def read_page_meta(path: Path) -> dict[str, str]:
    parser = _MetaParser()
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            parser.feed(handle.read(256 * 1024))
    except OSError as exc:
        raise AuthorityError(f"cannot inspect sample page {path}: {exc}") from exc
    return parser.values


def page_freshness(output: Path, sample: dict[str, Any], authority: dict[str, Any]) -> tuple[bool, str]:
    if not output.is_file():
        return False, "missing page"
    if output.stat().st_size < 1024:
        return False, "page is unexpectedly small"
    markers = read_page_meta(output)
    expected = {
        PAGE_META_NAMES["summary_sha256"]: authority["summary_sha256"],
        PAGE_META_NAMES["region_sha256"]: sample["region_sha256"],
        PAGE_META_NAMES["read_af_sha256"]: sample["read_af_sha256"],
        PAGE_META_NAMES["sample"]: sample["sample"],
        PAGE_META_NAMES["renderer_sha256"]: sha256_file(BUILD),
        PAGE_META_NAMES["ui_contract"]: EXPECTED_UI_CONTRACT,
    }
    for name, value in expected.items():
        if markers.get(name) != value:
            return False, f"marker {name} mismatch"
    newest_input = max(
        authority["summary_path"].stat().st_mtime,
        sample["region_path"].stat().st_mtime,
        sample["read_af_path"].stat().st_mtime,
        BUILD.stat().st_mtime,
    )
    if output.stat().st_mtime < newest_input:
        return False, "page mtime predates its summary or region view"
    return True, "hash-bound"


def collect_rows(authority: dict[str, Any], build_samples: bool) -> list[dict[str, Any]]:
    """Build or verify seven hash-bound sample pages without fallback data."""
    rows: list[dict[str, Any]] = []
    for sample in authority["samples"]:
        name = sample["sample"]
        output = OUTDIR / f"{name}.html"
        if build_samples:
            env = dict(
                os.environ,
                SM_RV=str(sample["region_path"]),
                SM_OUT=str(output),
                SM_CURRENT_SUMMARY=str(authority["summary_path"]),
                SM_CURRENT_SUMMARY_SHA256=authority["summary_sha256"],
                SM_CANONICAL_SAMPLE=json.dumps(
                    sample["canonical_record"], ensure_ascii=False, sort_keys=True, separators=(",", ":")
                ),
                SM_REGION_VIEW_SHA256=sample["region_sha256"],
                SM_READ_AF_TOPOLOGY=str(sample["read_af_path"]),
                SM_READ_AF_TOPOLOGY_SHA256=sample["read_af_sha256"],
            )
            result = subprocess.run(
                [PYTHON, str(BUILD)],
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            if result.returncode != 0:
                output_excerpt = result.stdout.decode(errors="replace")[-1200:]
                raise AuthorityError(f"sample build failed for {name} (rc={result.returncode}):\n{output_excerpt}")

        ready, reason = page_freshness(output, sample, authority)
        require(ready, f"{'built' if build_samples else 'existing'} sample page is stale for {name}: {reason}")
        row = dict(sample)
        row.update({"output": output, "page_ready": True, "page_size": output.stat().st_size})
        rows.append(row)
        print(f"[{name}] {'BUILT' if build_samples else 'VERIFIED'} hash-bound page")
    require(len(rows) == 7, "did not collect exactly seven canonical sample pages")
    return rows


def percentage(value: int, total: int) -> float:
    return 100.0 * value / total if total else 0.0


def topology_segments(topology: dict[str, Any], total: int) -> str:
    segments = []
    for key, short_label, long_label, class_name in TOPOLOGY_CLASSES:
        count = int(topology[key])
        share = percentage(count, total)
        segments.append(
            f'<span class="bar-segment segment-{class_name}" style="width:{share:.4f}%" '
            f'title="{escaped(short_label)} · {count:,} · {share:.1f}%">'
            f'<span class="sr-only">{escaped(long_label)} {count:,}，{share:.1f}%</span></span>'
        )
    return "".join(segments)


def topology_legend(aggregate: dict[str, Any]) -> str:
    total = int(aggregate["W_primary"])
    topology = aggregate["topology_classes"]
    items = []
    for key, short_label, long_label, class_name in TOPOLOGY_CLASSES:
        count = int(topology[key])
        items.append(
            f'<li><span class="legend-swatch segment-{class_name}" aria-hidden="true"></span>'
            f'<span><strong>{escaped(short_label)}</strong><small>{escaped(long_label)}</small></span>'
            f'<b>{count:,}<small>{percentage(count, total):.1f}%</small></b></li>'
        )
    return "".join(items)


def aggregate_sidecar_dimension(
    rows: list[dict[str, Any]],
    summary_key: str,
    classes: tuple[tuple[str, str, str, str], ...],
) -> dict[str, int]:
    """Aggregate one mutually exclusive current-v5 read-AF sidecar partition."""
    counts = {
        key: sum(
            int((row["read_af_summary"].get(summary_key) or {}).get(key, 0))
            for row in rows
        )
        for key, _, _, _ in classes
    }
    total = sum(counts.values())
    expected = sum(int(row["W_primary"]) for row in rows)
    require(total == expected, f"{summary_key} cohort partition mismatch ({total} != {expected})")
    return counts


def distribution_segments(
    counts: dict[str, int],
    total: int,
    classes: tuple[tuple[str, str, str, str], ...],
) -> str:
    segments = []
    for key, short_label, _, class_name in classes:
        count = int(counts[key])
        share = percentage(count, total)
        segments.append(
            f'<span class="bar-segment segment-{class_name}" style="width:{share:.4f}%" '
            f'title="{escaped(short_label)} · {count:,} · {share:.1f}%">'
            f'<span class="sr-only">{escaped(short_label)} {count:,}，{share:.1f}%</span></span>'
        )
    return "".join(segments)


def distribution_legend(
    counts: dict[str, int],
    total: int,
    classes: tuple[tuple[str, str, str, str], ...],
) -> str:
    items = []
    for key, short_label, long_label, class_name in classes:
        count = int(counts[key])
        items.append(
            f'<li><span class="legend-swatch segment-{class_name}" aria-hidden="true"></span>'
            f'<span><strong>{escaped(short_label)}</strong><small>{escaped(long_label)}</small></span>'
            f'<b>{count:,}<small>{percentage(count, total):.1f}%</small></b></li>'
        )
    return "".join(items)


def cohort_read_af_morphology(rows: list[dict[str, Any]]) -> str:
    total = sum(int(row["W_primary"]) for row in rows)
    selection = aggregate_sidecar_dimension(
        rows, "selection_classes", READ_AF_SELECTION_CLASSES
    )
    morphology = aggregate_sidecar_dimension(rows, "morphology_classes", MORPHOLOGY_CLASSES)
    selection_label = ", ".join(
        f"{short_label} {selection[key]:,}"
        for key, short_label, _, _ in READ_AF_SELECTION_CLASSES
    )
    morphology_label = ", ".join(
        f"{short_label} {morphology[key]:,}"
        for key, short_label, _, _ in MORPHOLOGY_CLASSES
    )
    return (
        '<div class="interpretation-board">'
        '<article class="distribution-card"><div class="distribution-head">'
        '<div><span>01 · Read-AF selection</span><h3>第一順位縮到哪一層？</h3></div>'
        f'<b>7-dataset W_primary {total:,}</b></div><p>先保留 structural determinacy，再以同一 primary HP family 的 '
        'ALT/(REF+ALT) 與 score=Σ(AF_ancestor−AF_newly-acquired) 做描述性排序；只在同一 unit 內比較，不是 VCF caller AF。</p>'
        f'<div class="topology-bar" role="img" aria-label="全 cohort read-AF selection：{escaped(selection_label)}">'
        f'{distribution_segments(selection, total, READ_AF_SELECTION_CLASSES)}</div>'
        '<div class="bar-scale"><span>0</span><span>read-AF selection composition</span><span>100%</span></div>'
        f'<ul class="distribution-legend">{distribution_legend(selection, total, READ_AF_SELECTION_CLASSES)}</ul>'
        '<aside>「結構 exact 唯一」與「read-AF 唯一第一順位」都不是正交生物確認；score 也不是 posterior、CCF 或 calibrated probability。</aside>'
        '</article>'
        '<article class="distribution-card"><div class="distribution-head">'
        '<div><span>02 · Mutation-state morphology</span><h3>單支、直系、旁系如何分布？</h3></div>'
        f'<b>7-dataset W_primary {total:,}</b></div><p>由 structural Topo=1 或 read-AF co-top 唯一 shape 的 primary-HP graph pattern 重算。</p>'
        f'<div class="topology-bar" role="img" aria-label="全 cohort morphology：{escaped(morphology_label)}">'
        f'{distribution_segments(morphology, total, MORPHOLOGY_CLASSES)}</div>'
        '<div class="bar-scale"><span>0</span><span>morphology composition</span><span>100%</span></div>'
        f'<ul class="distribution-legend">{distribution_legend(morphology, total, MORPHOLOGY_CLASSES)}</ul>'
        '<aside>直系／旁系只描述 mutation-state nesting／branching；HP 間不建立親緣，不能當 clone census 或 confirmed subclone。</aside>'
        '</article></div>'
    )


def comparison_class_contract(dimension: str) -> list[dict[str, str]]:
    contract = COMPARISON_DIMENSIONS[dimension]
    return [
        {
            "key": key,
            "label": short_label,
            "description": long_label,
            "class_name": class_name,
        }
        for key, short_label, long_label, class_name in contract["classes"]
    ]


def comparison_legend(dimension: str) -> str:
    return "".join(
        f'<span class="comparison-key-item"><i class="legend-swatch segment-{escaped(item["class_name"])}" '
        f'aria-hidden="true"></i><b>{escaped(item["label"])}</b></span>'
        for item in comparison_class_contract(dimension)
    )


def comparison_profile_rows(comparison: dict[str, Any], dimension: str) -> str:
    classes = comparison_class_contract(dimension)
    rows = []
    for record in comparison["datasets"]:
        sample = record["dataset"]
        values = record["dimensions"][dimension]
        segments = "".join(
            f'<span class="bar-segment segment-{escaped(item["class_name"])}" '
            f'style="width:{float(values["proportions"][item["key"]]) * 100:.6f}%" '
            f'title="{escaped(item["label"])} · {int(values["counts"][item["key"]]):,} · '
            f'{float(values["proportions"][item["key"]]) * 100:.1f}%"></span>'
            for item in classes
        )
        technical = sample in {"HCC1395", "HCC1395_DORADO"}
        role = "5kHz canonical" if sample == "HCC1395" else "Dorado technical pair"
        rows.append(
            f'<div class="comparison-profile-row{" technical-pair" if technical else ""}" '
            f'data-dataset="{escaped(sample)}"><div class="comparison-profile-name">'
            f'<a href="{escaped(sample)}.html"><strong>{escaped(sample)}</strong></a>'
            f'{f"<small>{escaped(role)}</small>" if technical else ""}</div>'
            f'<div class="topology-bar" role="img" aria-label="{escaped(sample)} {escaped(COMPARISON_DIMENSIONS[dimension]["label"])} composition">'
            f'{segments}</div><b class="comparison-profile-n">{int(record["W_primary"]):,}<small>W_primary</small></b></div>'
        )
    return "".join(rows)


def comparison_ledger_rows(comparison: dict[str, Any], dimension: str) -> str:
    classes = comparison_class_contract(dimension)
    rendered = []
    for record in comparison["datasets"]:
        values = record["dimensions"][dimension]
        cells = "".join(
            f'<td class="numeric" data-label="{escaped(item["label"])}">'
            f'<strong>{int(values["counts"][item["key"]]):,}</strong>'
            f'<small>{float(values["proportions"][item["key"]]) * 100:.1f}%</small></td>'
            for item in classes
        )
        rendered.append(
            f'<tr><th scope="row" class="sample-cell"><a href="{escaped(record["dataset"])}.html">'
            f'<strong>{escaped(record["dataset"])}</strong></a><small>{escaped(record["biological_id"])} biological ID</small></th>'
            f'<td class="numeric" data-label="W_primary">{int(record["W_primary"]):,}</td>{cells}</tr>'
        )
    return "".join(rendered)


def comparison_ledger_head(dimension: str) -> str:
    classes = comparison_class_contract(dimension)
    return (
        '<tr><th scope="col">Dataset</th><th scope="col" class="numeric">W_primary</th>'
        + "".join(
            f'<th scope="col" class="numeric">{escaped(item["label"])}</th>'
            for item in classes
        )
        + "</tr>"
    )


def comparison_matrix_table(comparison: dict[str, Any], dimension: str) -> str:
    order = comparison["pairwise_composition"]["matrices"]["order"]
    matrix = comparison["pairwise_composition"]["matrices"]["tvd"][dimension]
    short = {
        "COLO829": "COLO829",
        "H1437": "H1437",
        "H2009": "H2009",
        "HCC1395": "HCC1395",
        "HCC1395_DORADO": "DORADO",
        "HCC1937": "HCC1937",
        "HCC1954": "HCC1954",
    }
    header = "".join(f'<th scope="col">{escaped(short[sample])}</th>' for sample in order)
    body = []
    for left in order:
        cells = []
        for right in order:
            value = float(matrix[left][right])
            if left == right:
                cells.append('<td class="matrix-diagonal"><span aria-label="same dataset">0</span></td>')
                continue
            hcc_pair = {left, right} == {"HCC1395", "HCC1395_DORADO"}
            alpha = 0.08 + min(value / 0.5, 1.0) * 0.64
            cells.append(
                f'<td><button type="button" class="matrix-cell{" is-hcc" if hcc_pair else ""}" '
                f'data-left="{escaped(left)}" data-right="{escaped(right)}" '
                f'style="background:rgba(22,139,141,{alpha:.4f})" '
                f'aria-label="{escaped(left)} 與 {escaped(right)}，TVD {value:.3f}">{value:.3f}</button></td>'
            )
        body.append(
            f'<tr><th scope="row">{escaped(short[left])}</th>{"".join(cells)}</tr>'
        )
    return (
        '<table class="distance-matrix"><caption class="sr-only">7 datasets pairwise composition TVD</caption>'
        f'<thead><tr><th scope="col">TVD</th>{header}</tr></thead><tbody>{"".join(body)}</tbody></table>'
    )


def comparison_pair_inspector(comparison: dict[str, Any], dimension: str) -> str:
    record = next(
        row
        for row in comparison["pairwise_composition"]["records"]
        if {row["left"], row["right"]} == {"HCC1395", "HCC1395_DORADO"}
    )
    metric = record["dimensions"][dimension]
    return (
        '<span class="pair-chip">Same biological sample · technical pair</span>'
        '<h3>HCC1395 × HCC1395_DORADO</h3>'
        f'<strong class="pair-distance">{float(metric["tvd"]):.3f}<small>TVD · rank '
        f'{int(metric["tvd_rank_among_21_pairs"])}/21</small></strong>'
        '<p>此格只比較類別比例；下方 exact-region panel 才檢查相同區域是否得到相同標籤。</p>'
    )


def sample_signature_cards(comparison: dict[str, Any]) -> str:
    records = comparison["datasets"]

    def maximum(dimension: str, keys: tuple[str, ...]) -> tuple[dict[str, Any], float]:
        candidates = [
            (
                record,
                sum(float(record["dimensions"][dimension]["proportions"][key]) for key in keys),
            )
            for record in records
        ]
        return max(candidates, key=lambda item: item[1])

    signatures = (
        (
            "可評估性負擔",
            *maximum("structural", ("incomplete",)),
            "incomplete",
        ),
        (
            "Exact 可辨識度",
            *maximum("structural", ("exact_and_topology_unique",)),
            "exact + Topo unique",
        ),
        (
            "Read-AF 並列負擔",
            *maximum("read_af_selection", ("read_af_tied_same_topology", "read_af_tied_different_topology")),
            "co-top ties",
        ),
        (
            "分枝相容形態",
            *maximum("morphology", ("sister_branch", "direct_and_sister")),
            "sister + mixed",
        ),
    )
    return "".join(
        f'<article class="signature-card"><span>{escaped(title)}</span><strong>{escaped(record["dataset"])}</strong>'
        f'<b>{value * 100:.1f}%</b><small>{escaped(note)} · W_primary {int(record["W_primary"]):,}</small></article>'
        for title, record, value, note in signatures
    )


def hcc_matched_rows(comparison: dict[str, Any]) -> str:
    hcc = comparison["hcc1395_technical_pair"]
    rows = []
    for dimension, contract in COMPARISON_DIMENSIONS.items():
        values = hcc["matched_primary_dimensions"][dimension]
        agreement_ci = values["bootstrap"]["raw_agreement_ci95"]
        kappa_ci = values["bootstrap"]["cohen_kappa_ci95"]
        rows.append(
            '<tr>'
            f'<th scope="row">{escaped(contract["label"])}</th>'
            f'<td class="numeric"><strong>{float(values["raw_agreement"]) * 100:.2f}%</strong>'
            f'<small>{int(values["agreement_count"]):,} / {int(values["n"]):,}</small></td>'
            f'<td class="numeric">{float(agreement_ci[0]) * 100:.2f}–{float(agreement_ci[1]) * 100:.2f}%</td>'
            f'<td class="numeric"><strong>{float(values["cohen_kappa"]):.3f}</strong>'
            f'<small>{float(kappa_ci[0]):.3f}–{float(kappa_ci[1]):.3f}</small></td>'
            f'<td><span class="verdict-badge verdict-{escaped(values["raw_agreement_band"])}">'
            f'{escaped(values["raw_agreement_band"])}</span></td></tr>'
        )
    return "".join(rows)


def hcc_confusion_table(comparison: dict[str, Any], dimension: str) -> str:
    classes = comparison_class_contract(dimension)
    matrix = comparison["hcc1395_technical_pair"]["matched_primary_dimensions"][dimension]["matrix"]
    header = "".join(f'<th scope="col">{escaped(item["label"])}</th>' for item in classes)
    body = []
    for left in classes:
        cells = "".join(
            f'<td class="numeric{" confusion-diagonal" if left["key"] == right["key"] else ""}">'
            f'{int(matrix[left["key"]][right["key"]]):,}</td>'
            for right in classes
        )
        body.append(f'<tr><th scope="row">{escaped(left["label"])}</th>{cells}</tr>')
    return (
        '<table class="confusion-table"><caption>Rows = HCC1395 5kHz；columns = HCC1395_DORADO</caption>'
        f'<thead><tr><th scope="col">5kHz ↓ / DORADO →</th>{header}</tr></thead>'
        f'<tbody>{"".join(body)}</tbody></table>'
    )


def compact_comparison_payload(comparison: dict[str, Any]) -> dict[str, Any]:
    hcc = comparison["hcc1395_technical_pair"]
    compact_hcc_dimensions = {}
    for dimension, values in hcc["matched_primary_dimensions"].items():
        compact_hcc_dimensions[dimension] = {
            key: value for key, value in values.items() if key != "by_chromosome"
        }
    return {
        "dimensions": {
            dimension: {
                "label": contract["label"],
                "short": contract["short"],
                "classes": comparison_class_contract(dimension),
            }
            for dimension, contract in COMPARISON_DIMENSIONS.items()
        },
        "datasets": comparison["datasets"],
        "aggregates": comparison["aggregates"],
        "pairwise": comparison["pairwise_composition"],
        "hcc": {
            "marginal_profiles": hcc["marginal_profiles"],
            "exact_coordinate_overlap": hcc["exact_coordinate_overlap"],
            "matched_primary_dimensions": compact_hcc_dimensions,
        },
    }


def sample_comparison_dashboard(comparison: dict[str, Any]) -> str:
    hcc = comparison["hcc1395_technical_pair"]
    marginal = hcc["marginal_profiles"]
    overlap = hcc["exact_coordinate_overlap"]["primary_both"]
    prior = comparison["independent_current_v5_hcc_evidence"]
    exact = comparison["hcc_exact_signature_evidence"]
    retained = prior["exact_retained_site_overlap"]
    same_internal = exact["internal_ssnv_pairing"]["same_set_within_exact_common_region"]
    shape = exact["shape_agreement"]["same_internal_ssnv_both_shape_resolved"]
    exact_edges = exact["exact_labeled_edge_agreement"]["same_internal_ssnv_both_candidate_unique"]
    bootstrap = marginal["chromosome_block_bootstrap"]
    embedded = json.dumps(
        compact_comparison_payload(comparison),
        ensure_ascii=False,
        separators=(",", ":"),
    ).replace("</", "<\\/")
    return f'''
  <section class="section" id="cohort-comparison" aria-labelledby="comparison-title">
    <div class="section-head"><div><p class="section-kicker">Cross-dataset comparison</p><h2 id="comparison-title">樣本狀況、比例與差異放在同一張桌上</h2></div><p class="section-note">先看 7-dataset operational profile，再以 6-biological-sample macro 防止 HCC1395 technical pair 被當成兩個生物樣本。TVD 是類別比例距離，不是 clone／演化距離。</p></div>
    <div class="comparison-scope-grid">
      <article><span>Operational micro</span><strong>7 datasets</strong><b>{int(comparison["aggregates"]["dataset_micro"]["W_primary"]):,}</b><small>W_primary region rows；HCC pair 各占一列</small></article>
      <article><span>Biological macro</span><strong>6 biological IDs</strong><b>Equal weight</b><small>HCC pair 先平均，再與其他五個 ID 平均</small></article>
      <article class="technical"><span>Validation pair</span><strong>HCC1395 × DORADO</strong><b>1 biological ID</b><small>technical／cross-platform pair；不是 biological replication</small></article>
    </div>
    <div class="signature-grid" aria-label="樣本 profile 的自動極值摘要">{sample_signature_cards(comparison)}</div>
    <div class="comparison-controls" role="group" aria-label="切換樣本比較維度">
      {''.join(f'<button type="button" class="comparison-tab{" is-active" if dimension == "structural" else ""}" data-comparison-dimension="{escaped(dimension)}" aria-pressed="{"true" if dimension == "structural" else "false"}"><span>{escaped(contract["short"])}</span>{escaped(contract["label"])}</button>' for dimension, contract in COMPARISON_DIMENSIONS.items())}
    </div>
    <div class="comparison-key" id="comparison-key" aria-live="polite">{comparison_legend("structural")}</div>
    <article class="comparison-profile-panel">
      <div class="comparison-panel-head"><div><span>01 · Composition profiles</span><h3 id="comparison-profile-title">結構可辨識度 · 逐 dataset</h3></div><p>每列獨立以該 dataset 的 W_primary=100%；長條適合看 pattern，精確 count／% 在下方表格。</p></div>
      <div class="comparison-profile-list" id="comparison-profile-chart" aria-live="polite">{comparison_profile_rows(comparison, "structural")}</div>
    </article>
    <details class="comparison-ledger">
      <summary>查看目前維度的 7 datasets 精確 count 與比例</summary>
      <div class="table-scroll" role="region" aria-label="樣本拓撲比例精確表" tabindex="0"><table class="comparison-table"><thead id="comparison-ledger-head">{comparison_ledger_head("structural")}</thead><tbody id="comparison-ledger-body">{comparison_ledger_rows(comparison, "structural")}</tbody></table></div>
    </details>
    <div class="comparison-matrix-layout">
      <article class="matrix-panel"><div class="comparison-panel-head"><div><span>02 · Pairwise distance</span><h3 id="comparison-matrix-title">結構可辨識度 · TVD</h3></div><p>0=比例完全相同；數值越高差異越大。每格直接印數值，色深只作輔助。</p></div><div class="heatmap-scroll" role="region" aria-label="7乘7樣本組成距離矩陣" tabindex="0" id="comparison-matrix">{comparison_matrix_table(comparison, "structural")}</div></article>
      <aside class="pair-inspector" id="pair-inspector" aria-live="polite">{comparison_pair_inspector(comparison, "structural")}</aside>
    </div>
  </section>

  <section class="section" id="hcc1395-technical-validation" aria-labelledby="hcc-validation-title">
    <div class="section-head"><div><p class="section-kicker">HCC1395 technical reproducibility</p><h2 id="hcc-validation-title">相對第 2 近，但不是可互換的同一份結果</h2></div><p class="section-note">兩列資料屬同一 biological sample；驗證的是 current pipeline 跨技術／輸入重現性。結論必同時滿足 marginal composition 與 matched-region 兩個視角。</p></div>
    <article class="hcc-verdict">
      <div><span class="verdict-badge verdict-moderate">PARTIAL TECHNICAL REPRODUCIBILITY</span><h3>粗粒度拓撲中度重現；exact ancestry 尚未確認</h3><p>三維 profile 在 21 組中排名第 {int(marginal["rank_among_21_pairs"])} 近，但 full-profile TVD={float(marginal["profile_mean_tvd"]):.3f} 已跨過事先鎖定的 0.10 明顯差異門檻。共同 primary regions 的三類 κ 都只落在 moderate。</p></div>
      <dl><div><dt>Full profile</dt><dd>{float(marginal["profile_mean_tvd"]):.3f}<small>95% chr-block {float(bootstrap["full_profile_mean_tvd_ci95"][0]):.3f}–{float(bootstrap["full_profile_mean_tvd_ci95"][1]):.3f}</small></dd></div><div><dt>Conditional evaluable</dt><dd>{float(marginal["conditional_evaluable_mean_tvd"]):.3f}<small>rank {int(marginal["conditional_evaluable_rank_among_21_pairs"])}/21</small></dd></div></dl>
    </article>
    <div class="hcc-metric-grid">
      <article><span>Exact primary region</span><strong>{int(overlap["intersection"]):,}</strong><b>Jaccard {float(overlap["jaccard"]):.3f}</b><small>{float(overlap["left_coverage"])*100:.1f}%／{float(overlap["right_coverage"])*100:.1f}% coverage</small></article>
      <article><span>Exact retained sSNV</span><strong>{int(retained["intersection"]):,}</strong><b>Jaccard {float(retained["jaccard"]):.3f}</b><small>site backbone 高；不等於 tree 相同</small></article>
      <article><span>Unlabeled rooted shape</span><strong>{float(shape["rate"])*100:.1f}%</strong><b>{int(shape["numerator"]):,} / {int(shape["denominator"]):,}</b><small>same internal sSNV；忽略 HP label</small></article>
      <article><span>Exact labeled edges</span><strong>{float(exact_edges["rate"])*100:.1f}%</strong><b>{int(exact_edges["numerator"]):,} / {int(exact_edges["denominator"]):,}</b><small>unique candidate；exact ancestry 重現偏弱</small></article>
    </div>
    <div class="hcc-analysis-grid">
      <article class="hcc-table-panel"><div class="comparison-panel-head"><div><span>03 · Exact-coordinate agreement</span><h3>共同 primary region 的分類是否一致？</h3></div><p>95% CI 以 22 條 autosome block bootstrap，避免把 6,402 regions 當獨立 replicate。</p></div><div class="table-scroll" role="region" aria-label="HCC1395 matched region agreement" tabindex="0"><table class="hcc-agreement-table"><thead><tr><th scope="col">維度</th><th scope="col" class="numeric">Agreement</th><th scope="col" class="numeric">95% block CI</th><th scope="col" class="numeric">κ（95% CI）</th><th scope="col">判讀</th></tr></thead><tbody>{hcc_matched_rows(comparison)}</tbody></table></div></article>
      <aside class="evidence-ladder"><span>Independent current-v5 layers</span><ol><li><b>{float(retained["jaccard"])*100:.1f}%</b><small>retained-site Jaccard · high backbone overlap</small></li><li><b>{float(same_internal["rate"])*100:.1f}%</b><small>exact common regions with the same internal sSNV set</small></li><li><b>{float(shape["rate"])*100:.1f}%</b><small>same-site unlabeled rooted-shape agreement</small></li><li><b>{float(exact_edges["rate"])*100:.1f}%</b><small>same-site exact labeled-edge agreement</small></li></ol></aside>
    </div>
    <article class="confusion-panel"><div class="comparison-panel-head"><div><span>04 · Confusion matrix</span><h3 id="hcc-confusion-title">結構可辨識度 · 哪些類別互相轉換？</h3></div><p>對角線是相同 label；非對角線是 technical divergence，不稱「錯誤」，因兩側都沒有 ground truth 身分。</p></div><div class="heatmap-scroll" role="region" aria-label="HCC1395 技術配對 confusion matrix" tabindex="0" id="hcc-confusion">{hcc_confusion_table(comparison, "structural")}</div></article>
    <aside class="claim-boundary-box"><strong>判定界線</strong><p>可寫：「同一 HCC1395 的粗粒度 regional mutation-state profile 具有中度跨技術重現性。」不可寫：「兩個 dataset 得到相同真樹／clone tree」「已驗證 biological clone」或「可互換」。若要驗樣本身分，需另做 germline SNP fingerprint；topology comparison 只驗分析輸出。</p></aside>
  </section>
  <script type="application/json" id="sample-comparison-data">{embedded}</script>
'''
def sample_topology_rows(rows: list[dict[str, Any]]) -> str:
    rendered = []
    for row in rows:
        total = int(row["W_primary"])
        topology = row["topology_classes"]
        label = ", ".join(
            f"{short} {int(topology[key]):,}"
            for key, short, _, _ in TOPOLOGY_CLASSES
        )
        rendered.append(
            f'<article class="sample-topology"><div class="sample-topology-head">'
            f'<a href="{escaped(row["sample"])}.html"><strong>{escaped(row["sample"])}</strong></a>'
            f'<span>W_primary {total:,}</span></div>'
            f'<div class="topology-bar" role="img" aria-label="{escaped(row["sample"])}：{escaped(label)}">'
            f'{topology_segments(topology, total)}</div>'
            f'<div class="bar-scale"><span>0</span><span>候選拓撲狀態</span><span>100%</span></div></article>'
        )
    return "".join(rendered)


def genome_launchers(rows: list[dict[str, Any]]) -> str:
    links = []
    for index, row in enumerate(rows, start=1):
        links.append(
            f'<a class="genome-link" href="{escaped(row["sample"])}.html">'
            f'<span class="launcher-index">{index:02d}</span>'
            f'<span><strong>{escaped(row["sample"])}</strong>'
            f'<small>GRCh38 座標分布 · structural / read-AF / morphology 七組觀察 · {int(row["W_tree"]):,} regions</small></span>'
            f'<span class="launcher-action">進入巡覽</span></a>'
        )
    return "".join(links)


def sample_position_leaders(rows: list[dict[str, Any]]) -> str:
    """Render per-dataset W_primary chromosome leaders from bound region views."""
    rendered = []
    for row in rows:
        region_view = load_json(row["region_path"])
        counts = {f"chr{number}": 0 for number in range(1, 23)}
        regions = region_view.get("regions") or []
        for region in regions:
            chrom = region.get("chrom")
            if chrom not in counts:
                continue
            if any(line.get("is_primary_lineage") is True for line in (region.get("lineages") or [])):
                counts[chrom] += 1
        observed_primary = sum(counts.values())
        expected_primary = int(row["W_primary"])
        require(
            observed_primary == expected_primary,
            f'{row["sample"]}: position leader W_primary mismatch ({observed_primary} != {expected_primary})',
        )
        leader, count = min(counts.items(), key=lambda item: (-item[1], int(item[0][3:])))
        rendered.append(
            f'<a href="{escaped(row["sample"])}.html#chr={escaped(leader)}">'
            f'<span>{escaped(row["sample"])}</span><b>{escaped(leader)} · {count:,}</b></a>'
        )
    return "".join(rendered)


def cohort_rows(rows: list[dict[str, Any]]) -> str:
    rendered = []
    for row in rows:
        complete = int(row["complete_regions"])
        primary = int(row["W_primary"])
        rendered.append(
            '<tr>'
            f'<th scope="row" class="sample-cell"><a href="{escaped(row["sample"])}.html">'
            f'<strong>{escaped(row["sample"])}</strong></a><small>canonical PASS</small></th>'
            f'<td class="numeric" data-label="LPS PASS">{int(row["tree_input_records"]):,}</td>'
            f'<td class="numeric" data-label="Retained sSNV">{int(row["retained_sSNV"]):,}</td>'
            f'<td class="numeric" data-label="W tree">{int(row["W_tree"]):,}</td>'
            f'<td class="numeric" data-label="W primary">{primary:,}</td>'
            f'<td class="numeric" data-label="Complete"><strong>{complete:,}</strong>'
            f'<small>{percentage(complete, primary):.1f}% of W_primary</small></td>'
            f'<td class="numeric" data-label="Incomplete">{int(row["incomplete_regions"]):,}</td>'
            f'<td class="numeric" data-label="Primary units">{int(row["primary_units"]):,}</td>'
            f'<td><a class="table-action" href="{escaped(row["sample"])}.html">開啟全基因頁</a></td>'
            '</tr>'
        )
    return "".join(rendered)


def evidence_links(authority: dict[str, Any], rows: list[dict[str, Any]]) -> str:
    summary_href = os.path.relpath(authority["summary_path"], OUTDIR)
    comparison_meta = authority["summary"]["backbone_comparison"]
    comparison_path = resolve_authority_path(comparison_meta["path"])
    comparison_href = os.path.relpath(comparison_path, OUTDIR)
    links = [
        f'<li><a href="{escaped(summary_href)}">Current canonical machine summary JSON</a>'
        f'<code>{escaped(authority["summary_sha256"])}</code></li>',
        f'<li><a href="{escaped(comparison_href)}">Backbone sensitivity comparison JSON</a>'
        f'<code>{escaped(comparison_meta["sha256"])}</code></li>',
        f'<li><a href="{escaped(os.path.relpath(authority["read_af_index_path"], OUTDIR))}">'
        f'Current-v5 read-AF／morphology index JSON</a>'
        f'<code>{escaped(authority["read_af_index_sha256"])}</code></li>',
        f'<li><a href="{escaped(os.path.relpath(authority["sample_comparison_path"], OUTDIR))}">'
        f'Seven-dataset topology comparison JSON</a>'
        f'<code>{escaped(authority["sample_comparison_sha256"])}</code></li>',
        f'<li><a href="{escaped(os.path.relpath(authority["sample_comparison_receipt_path"], OUTDIR))}">'
        f'Sample-comparison validation receipt JSON</a>'
        f'<code>{escaped(authority["sample_comparison_receipt_sha256"])}</code></li>',
        f'<li><a href="{escaped(os.path.relpath(authority["hcc_exact_signature_path"], OUTDIR))}">'
        f'HCC1395 exact-signature validation JSON</a>'
        f'<code>{escaped(authority["hcc_exact_signature_sha256"])}</code></li>',
    ]
    for row in rows:
        links.append(
            f'<li><a href="{escaped(row["region_path"].as_uri())}">{escaped(row["sample"])} region-view JSON</a>'
            f'<code>{escaped(row["region_sha256"])}</code></li>'
        )
        links.append(
            f'<li><a href="{escaped(row["read_af_path"].as_uri())}">{escaped(row["sample"])} '
            f'read-AF／morphology sidecar JSON</a><code>{escaped(row["read_af_sha256"])}</code></li>'
        )
    return "".join(links)


def build_index(authority: dict[str, Any], rows: list[dict[str, Any]]) -> str:
    aggregate = authority["aggregate"]
    comparison = authority["summary"]["backbone_comparison"]["aggregate"]
    generated_at = datetime.now().astimezone().isoformat(timespec="minutes")
    verified_at = parse_timestamp(authority["verification"]["verified_at_utc"], "verification").astimezone()
    run_id = authority["run_root"].name
    topology = aggregate["topology_classes"]
    exact_unique = int(topology["exact_and_topology_unique"])
    shape_unique_exact_multiple = int(topology["topology_unique_exact_multiple"])
    multi_shape = int(topology["topology_multiple_exact_multiple"])
    incomplete = int(topology["incomplete"])
    shape_unique = exact_unique + shape_unique_exact_multiple
    complete_rate = percentage(int(aggregate["complete_regions"]), int(aggregate["W_primary"]))
    incomplete_rate = percentage(int(aggregate["incomplete_regions"]), int(aggregate["W_primary"]))
    no_primary_rate = percentage(int(aggregate["no_primary_lineage"]), int(aggregate["W_tree"]))

    template = Template("""<!DOCTYPE html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light">
<link rel="icon" href="data:,">
<meta name="intersubmod-current-summary-sha256" content="$summary_sha256">
<meta name="intersubmod-backbone-comparison-sha256" content="$comparison_sha256">
<meta name="intersubmod-read-af-topology-index-sha256" content="$read_af_index_sha256">
<meta name="intersubmod-sample-topology-comparison-sha256" content="$sample_comparison_sha256">
<meta name="intersubmod-canonical-run" content="$run_id">
<title>Layered reconstruction · 全基因 cohort command center</title>
<style>
:root{--ink:#122b32;--ink-soft:#345058;--muted:#596a6e;--paper:#f2efe7;--panel:#fffdf7;--line:#c8c6bc;--rail:#0f3138;--teal:#168b8d;--cyan:#53c4bd;--amber:#d47a2b;--brick:#b84f3d;--sand:#dfd5bc;--blue:#426f9e;--focus:#f2a900;--shadow:0 15px 45px rgba(19,45,52,.09)}
*{box-sizing:border-box}
html{background:var(--paper);scroll-behavior:smooth}
body{margin:0;color:var(--ink);background:linear-gradient(90deg,rgba(15,49,56,.035) 1px,transparent 1px),linear-gradient(180deg,#e8e4d9 0,var(--paper) 360px);background-size:48px 100%,100% 100%;font-family:"IBM Plex Sans","Noto Sans TC","PingFang TC","Microsoft JhengHei",sans-serif;line-height:1.58}
a{color:var(--rail);text-decoration-color:rgba(15,49,56,.42);text-underline-offset:3px}
a:hover{text-decoration-thickness:2px}
a:focus-visible,summary:focus-visible,.table-scroll:focus-visible{outline:3px solid var(--focus);outline-offset:3px}
.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}
.skip-link{position:fixed;z-index:100;left:14px;top:10px;transform:translateY(-180%);padding:9px 13px;border:2px solid var(--ink);background:#fff;color:#111;font-weight:800}
.skip-link:focus{transform:none}
.shell{width:min(100%,1360px);margin:0 auto;padding:22px clamp(14px,3vw,42px) 54px}
.hero{position:relative;display:grid;grid-template-columns:minmax(0,1.55fr) minmax(290px,.75fr);overflow:hidden;border:1px solid #789095;background:var(--rail);color:#f7f3e9;box-shadow:var(--shadow)}
.hero:after{content:"CHR 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22";position:absolute;left:0;right:0;bottom:0;padding:7px 24px;border-top:1px solid rgba(255,255,255,.18);color:rgba(255,255,255,.55);font:600 9px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.18em;word-spacing:.55em;white-space:nowrap;overflow:hidden}
.hero-main{padding:clamp(28px,5vw,66px) clamp(22px,5vw,64px) 60px}
.eyebrow{display:flex;flex-wrap:wrap;gap:8px;margin:0 0 20px;font:700 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.11em;text-transform:uppercase}
.eyebrow span{padding:5px 8px;border:1px solid rgba(255,255,255,.28)}
.eyebrow .canonical{border-color:var(--cyan);background:rgba(83,196,189,.12);color:#aaf0e8}
h1{max-width:830px;margin:0;font-family:"Iowan Old Style","Noto Serif TC","Songti TC",serif;font-size:clamp(34px,5.3vw,72px);font-weight:600;letter-spacing:-.045em;line-height:.98}
.lede{max-width:780px;margin:22px 0 0;color:#cbd9d9;font-size:clamp(15px,1.6vw,19px)}
.hero-aside{display:flex;flex-direction:column;justify-content:space-between;padding:30px 28px 54px;border-left:1px solid rgba(255,255,255,.18);background:rgba(0,0,0,.12)}
.authority-label{margin:0;color:#9bb1b4;font:700 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.1em;text-transform:uppercase}
.authority-value{display:block;margin-top:7px;overflow-wrap:anywhere;font:700 13px/1.45 "IBM Plex Mono","Cascadia Code",monospace}
.authority-list{display:grid;gap:17px;margin:0}.authority-list div{padding-bottom:14px;border-bottom:1px solid rgba(255,255,255,.13)}.authority-list div:last-child{border-bottom:0}
.sensitivity-banner{display:grid;grid-template-columns:minmax(235px,.8fr) minmax(0,1.7fr);gap:20px;margin-top:12px;padding:17px 19px;border:1px solid #c88a4f;border-left:7px solid var(--amber);background:#fff7e9}.sensitivity-banner h2{font-family:"IBM Plex Sans","Noto Sans TC",sans-serif;font-size:18px;font-weight:800;letter-spacing:0}.sensitivity-banner p{margin:6px 0 0;color:#664522;font-size:12px}.sensitivity-metrics{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:7px}.sensitivity-metric{padding:9px 10px;border:1px solid #dec4a5;background:#fffdf8}.sensitivity-metric span{display:block;color:#765737;font-size:9.5px}.sensitivity-metric strong{display:block;margin-top:5px;font:800 15px/1.1 "IBM Plex Mono","Cascadia Code",monospace}
.section{margin-top:34px}
.section-head{display:flex;align-items:end;justify-content:space-between;gap:24px;margin-bottom:13px;padding-bottom:10px;border-bottom:1px solid var(--line)}
.section-kicker{margin:0 0 4px;color:var(--teal);font:800 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.12em;text-transform:uppercase}
h2{margin:0;font-family:"Iowan Old Style","Noto Serif TC","Songti TC",serif;font-size:clamp(23px,3vw,37px);font-weight:600;letter-spacing:-.025em;line-height:1.1;text-wrap:balance}
.section-note{max-width:610px;margin:0;color:var(--muted);font-size:12.5px}
.metric-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px}
.metric-card{position:relative;min-height:170px;padding:20px;border:1px solid var(--line);background:var(--panel);box-shadow:0 6px 20px rgba(19,45,52,.035)}
.metric-card:before{content:"";position:absolute;left:-1px;right:-1px;top:-1px;height:4px;background:var(--teal)}
.metric-card:nth-child(2):before{background:var(--blue)}.metric-card:nth-child(3):before{background:var(--cyan)}.metric-card:nth-child(4):before{background:var(--brick)}
.metric-code{color:var(--muted);font:800 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.1em;text-transform:uppercase}
.metric-number{display:block;margin:18px 0 7px;font:700 clamp(29px,3vw,43px)/1 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:-.05em;font-variant-numeric:tabular-nums}
.metric-card p{margin:0;color:var(--ink-soft);font-size:12.5px}
.metric-sub{display:block;margin-top:8px;color:var(--muted);font:700 10.5px/1.4 "IBM Plex Mono","Cascadia Code",monospace}
.topology-board{display:grid;grid-template-columns:minmax(0,1.55fr) minmax(320px,.75fr);gap:12px}
.topology-main,.topology-legend,.sample-topology,.dimension-card,.claim-card,.layer-card,.table-shell,.evidence-drawer{border:1px solid var(--line);background:var(--panel)}
.topology-main{padding:24px}.topology-main h3{margin:0 0 5px;font-size:16px}.topology-main>p{margin:0 0 24px;color:var(--muted);font-size:12.5px}
.topology-bar{display:flex;width:100%;height:30px;overflow:hidden;border:1px solid #829195;background:#e8e4d9}
.bar-segment{display:block;min-width:0;height:100%}.segment-unique{background:var(--teal)}.segment-shape{background:var(--blue)}.segment-multiple{background:var(--amber)}.segment-incomplete{background:var(--brick)}
.segment-read-structural{background:#138b78}.segment-read-unique{background:#426f9e}.segment-read-same{background:#53bcb8}.segment-read-multiple{background:#8a4f9e}.segment-read-unavailable{background:#d08a31}
.segment-morph-single{background:#809196}.segment-morph-direct{background:#426f9e}.segment-morph-sister{background:#35a9a5}.segment-morph-mixed{background:#7e50a0}
.bar-scale{display:flex;justify-content:space-between;margin-top:6px;color:var(--muted);font:700 9.5px/1.2 "IBM Plex Mono","Cascadia Code",monospace;text-transform:uppercase}
.topology-legend{padding:18px}.topology-legend ul{display:grid;gap:12px;margin:0;padding:0;list-style:none}.topology-legend li{display:grid;grid-template-columns:10px 1fr auto;gap:10px;align-items:start}.legend-swatch{width:10px;height:30px}.topology-legend strong,.topology-legend small,.topology-legend b{display:block}.topology-legend strong{font:750 11px/1.3 "IBM Plex Mono","Cascadia Code",monospace}.topology-legend small{margin-top:2px;color:var(--muted);font-size:10.5px}.topology-legend b{text-align:right;font:750 13px/1.2 "IBM Plex Mono","Cascadia Code",monospace}.topology-legend b small{font-size:9.5px}
.sample-bars{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px;margin-top:12px}
.sample-topology{padding:15px 16px}.sample-topology-head{display:flex;justify-content:space-between;gap:12px;margin-bottom:9px}.sample-topology-head a{font-size:12.5px}.sample-topology-head span{color:var(--muted);font:700 10px/1.3 "IBM Plex Mono","Cascadia Code",monospace}.sample-topology .topology-bar{height:16px}
.interpretation-board{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}.distribution-card{min-width:0;padding:22px;border:1px solid var(--line);background:var(--panel)}.distribution-card:nth-child(1){border-top:4px solid var(--blue)}.distribution-card:nth-child(2){border-top:4px solid #7e50a0}.distribution-head{display:flex;justify-content:space-between;gap:16px;align-items:start}.distribution-head span{color:var(--muted);font:800 9.5px/1.25 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.06em;text-transform:uppercase}.distribution-head h3{margin:6px 0 0;font-size:17px}.distribution-head>b{color:var(--muted);font:750 10px/1.4 "IBM Plex Mono","Cascadia Code",monospace;white-space:nowrap}.distribution-card>p{min-height:40px;margin:10px 0 18px;color:var(--ink-soft);font-size:12px}.distribution-card .topology-bar{height:24px}.distribution-legend{display:grid;gap:7px;margin:17px 0 0;padding:0;list-style:none}.distribution-legend li{display:grid;grid-template-columns:10px minmax(0,1fr) auto;gap:9px;align-items:start;padding-top:7px;border-top:1px solid #e3e0d7}.distribution-legend .legend-swatch{height:29px}.distribution-legend strong,.distribution-legend small,.distribution-legend b{display:block}.distribution-legend strong{font-size:11.5px}.distribution-legend small{margin-top:2px;color:var(--muted);font-size:10px}.distribution-legend b{text-align:right;font:750 12px/1.2 "IBM Plex Mono","Cascadia Code",monospace}.distribution-card aside{margin-top:16px;padding:10px 12px;border-left:4px solid var(--amber);background:#f8efe2;color:#624a2e;font-size:11px}
.comparison-scope-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}.comparison-scope-grid article{min-height:146px;padding:18px;border:1px solid var(--line);border-top:4px solid var(--teal);background:var(--panel)}.comparison-scope-grid article:nth-child(2){border-top-color:var(--blue)}.comparison-scope-grid article.technical{border-top-color:var(--amber);background:#fff9ed}.comparison-scope-grid span,.comparison-panel-head span,.evidence-ladder>span{display:block;color:var(--muted);font:800 9.5px/1.25 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.07em;text-transform:uppercase}.comparison-scope-grid strong,.comparison-scope-grid b,.comparison-scope-grid small{display:block}.comparison-scope-grid strong{margin-top:14px;font-size:17px}.comparison-scope-grid b{margin-top:4px;font:800 14px/1.3 "IBM Plex Mono","Cascadia Code",monospace}.comparison-scope-grid small{margin-top:9px;color:var(--muted);font-size:10.5px}
.signature-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px;margin-top:10px}.signature-card{min-width:0;padding:13px 14px;border:1px solid #d2cfc5;background:#f7f4ec}.signature-card span,.signature-card strong,.signature-card b,.signature-card small{display:block}.signature-card span{color:var(--muted);font-size:9.5px}.signature-card strong{margin-top:6px;overflow-wrap:anywhere;font:800 11.5px/1.25 "IBM Plex Mono","Cascadia Code",monospace}.signature-card b{margin-top:5px;color:var(--teal);font:850 20px/1 "IBM Plex Mono","Cascadia Code",monospace}.signature-card small{margin-top:5px;color:var(--muted);font-size:9px}
.comparison-controls{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:7px;margin:16px 0 8px;padding:7px;border:1px solid #9da9a7;background:#e4e8e4}.comparison-tab{min-height:48px;padding:8px 12px;border:1px solid transparent;background:transparent;color:var(--ink);cursor:pointer;font:750 11px/1.3 "IBM Plex Sans","Noto Sans TC",sans-serif}.comparison-tab span{display:block;color:var(--muted);font:750 8.5px/1.2 "IBM Plex Mono","Cascadia Code",monospace;text-transform:uppercase}.comparison-tab.is-active{border-color:var(--rail);background:var(--rail);color:#fff}.comparison-tab.is-active span{color:#b5d7d4}.comparison-tab:focus-visible,.matrix-cell:focus-visible{outline:3px solid var(--focus);outline-offset:2px}
.comparison-key{display:flex;flex-wrap:wrap;gap:6px 14px;min-height:40px;padding:9px 11px;border:1px solid #d8d5cc;background:#faf8f2}.comparison-key-item{display:inline-grid;grid-template-columns:9px auto;gap:6px;align-items:center}.comparison-key-item .legend-swatch{width:9px;height:18px}.comparison-key-item b{font-size:9.5px}
.comparison-profile-panel,.matrix-panel,.hcc-table-panel,.confusion-panel{min-width:0;margin-top:10px;padding:19px;border:1px solid var(--line);background:var(--panel)}.comparison-panel-head{display:flex;justify-content:space-between;gap:24px;align-items:end;margin-bottom:15px}.comparison-panel-head h3{margin:5px 0 0;font-size:17px}.comparison-panel-head p{max-width:520px;margin:0;color:var(--muted);font-size:11px}.comparison-profile-list{display:grid;gap:7px}.comparison-profile-row{display:grid;grid-template-columns:minmax(150px,190px) minmax(0,1fr) 86px;gap:12px;align-items:center;min-height:43px;padding:6px 8px;border-left:3px solid transparent}.comparison-profile-row.technical-pair{border-left-color:var(--amber);background:#fff9ed}.comparison-profile-name strong,.comparison-profile-name small,.comparison-profile-n small{display:block}.comparison-profile-name strong{overflow-wrap:anywhere;font:800 11px/1.25 "IBM Plex Mono","Cascadia Code",monospace}.comparison-profile-name small{margin-top:2px;color:#8a5a24;font-size:8.5px}.comparison-profile-row .topology-bar{height:21px}.comparison-profile-n{text-align:right;font:800 11px/1.25 "IBM Plex Mono","Cascadia Code",monospace}.comparison-profile-n small{color:var(--muted);font-size:8px}.comparison-ledger{min-width:0;margin-top:8px;border:1px solid var(--line);background:var(--panel)}.comparison-ledger summary{min-height:48px;padding:13px 16px;cursor:pointer;font-weight:800}.comparison-ledger[open] summary{border-bottom:1px solid var(--line)}.comparison-table{min-width:1050px}.comparison-table strong,.comparison-table small{display:block}
.comparison-matrix-layout{display:grid;grid-template-columns:minmax(0,1.6fr) minmax(250px,.6fr);gap:10px;align-items:stretch}.heatmap-scroll{width:100%;min-width:0;max-width:100%;overflow-x:auto;overscroll-behavior-inline:contain}.heatmap-scroll:focus-visible{outline:3px solid var(--focus);outline-offset:2px}.distance-matrix{min-width:760px;table-layout:fixed;border-collapse:separate;border-spacing:3px}.distance-matrix th,.distance-matrix td{position:static;padding:0;border:0;background:transparent;text-align:center}.distance-matrix thead th,.distance-matrix tbody th{padding:5px;background:transparent;color:var(--muted);font:750 8.5px/1.15 "IBM Plex Mono","Cascadia Code",monospace}.distance-matrix tbody th{text-align:left}.matrix-cell,.matrix-diagonal span{display:grid;width:100%;min-height:44px;place-items:center;border:1px solid rgba(15,49,56,.22);color:#0a2429;font:800 10px/1 "IBM Plex Mono","Cascadia Code",monospace;font-variant-numeric:tabular-nums}.matrix-cell{cursor:pointer}.matrix-cell:hover{border-color:var(--rail);box-shadow:inset 0 0 0 2px #fff}.matrix-cell.is-hcc{outline:2px solid var(--amber);outline-offset:-2px}.matrix-cell.is-selected{box-shadow:inset 0 0 0 3px #fff;outline:3px solid var(--rail);outline-offset:-1px}.matrix-diagonal span{border-style:dashed;background:#ece9e0;color:#7a8584}
.pair-inspector{min-width:0;margin-top:10px;padding:21px;border:1px solid #9db0b2;border-top:5px solid var(--amber);background:#f8f5ed}.pair-chip,.verdict-badge{display:inline-block;padding:4px 7px;border:1px solid currentColor;color:#83501f;font:800 8.5px/1.2 "IBM Plex Mono","Cascadia Code",monospace;text-transform:uppercase}.pair-inspector h3{margin:16px 0 10px;overflow-wrap:anywhere;font-size:16px}.pair-distance{display:block;font:850 34px/1 "IBM Plex Mono","Cascadia Code",monospace}.pair-distance small{display:block;margin-top:6px;color:var(--muted);font-size:9px}.pair-inspector p{margin:18px 0 0;color:var(--muted);font-size:11px}.pair-delta-list{display:grid;gap:5px;margin-top:14px}.pair-delta-list div{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:8px;padding-top:5px;border-top:1px solid #ddd8cc}.pair-delta-list span{font-size:9.5px}.pair-delta-list b{font:750 9.5px/1.2 "IBM Plex Mono","Cascadia Code",monospace}
.hcc-verdict{display:grid;grid-template-columns:minmax(0,1.4fr) minmax(260px,.6fr);gap:18px;padding:22px;border:1px solid #c38a4c;border-left:7px solid var(--amber);background:#fff8e8}.hcc-verdict h3{margin:15px 0 7px;font-size:20px}.hcc-verdict p{margin:0;color:#654722;font-size:12px}.hcc-verdict dl{display:grid;grid-template-columns:1fr 1fr;gap:8px;margin:0}.hcc-verdict dl div{padding:12px;border:1px solid #dec5a6;background:#fffdf8}.hcc-verdict dt{color:var(--muted);font-size:9px}.hcc-verdict dd{margin:7px 0 0;font:850 24px/1 "IBM Plex Mono","Cascadia Code",monospace}.hcc-verdict dd small{display:block;margin-top:7px;color:var(--muted);font-size:8.5px;line-height:1.4}.verdict-moderate{color:#8a5720;background:#fff4d9}.verdict-high{color:#176d5b;background:#e5f3ee}.verdict-limited{color:#8a3f34;background:#f8e8e4}
.hcc-metric-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px;margin-top:9px}.hcc-metric-grid article{min-height:150px;padding:16px;border:1px solid var(--line);background:var(--panel)}.hcc-metric-grid span,.hcc-metric-grid strong,.hcc-metric-grid b,.hcc-metric-grid small{display:block}.hcc-metric-grid span{color:var(--muted);font-size:9px}.hcc-metric-grid strong{margin-top:13px;font:850 25px/1 "IBM Plex Mono","Cascadia Code",monospace}.hcc-metric-grid b{margin-top:7px;color:var(--teal);font:800 11px/1.3 "IBM Plex Mono","Cascadia Code",monospace}.hcc-metric-grid small{margin-top:10px;color:var(--muted);font-size:9.5px}
.hcc-analysis-grid{display:grid;grid-template-columns:minmax(0,1.55fr) minmax(250px,.45fr);gap:10px}.hcc-agreement-table{min-width:760px}.hcc-agreement-table strong,.hcc-agreement-table small{display:block}.evidence-ladder{margin-top:10px;padding:19px;border:1px solid #9dacac;background:#e9efed}.evidence-ladder ol{display:grid;gap:8px;margin:14px 0 0;padding:0;list-style:none}.evidence-ladder li{padding:10px;border-left:4px solid var(--teal);background:#f9fbf8}.evidence-ladder b,.evidence-ladder small{display:block}.evidence-ladder b{font:850 20px/1 "IBM Plex Mono","Cascadia Code",monospace}.evidence-ladder small{margin-top:5px;color:var(--muted);font-size:9.5px}.confusion-table{min-width:820px}.confusion-table caption{padding:0 0 10px;text-align:left;color:var(--muted);font-size:10px}.confusion-table th,.confusion-table td{position:static;padding:8px;border:2px solid var(--panel);background:#f3e7d3;text-align:center}.confusion-table thead th,.confusion-table tbody th{background:#e5e9e5;font-size:9px}.confusion-table .confusion-diagonal{background:#d7ebe4;color:#154d42;font-weight:850}.claim-boundary-box{display:grid;grid-template-columns:auto 1fr;gap:15px;align-items:start;margin-top:9px;padding:15px 17px;border-left:6px solid var(--brick);background:#f8eae6}.claim-boundary-box strong{font:850 9px/1.3 "IBM Plex Mono","Cascadia Code",monospace}.claim-boundary-box p{margin:0;color:#613c35;font-size:11px}
.dimension-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}.dimension-card{display:flex;min-width:0;min-height:306px;flex-direction:column;padding:20px;border-top:4px solid var(--teal)}.dimension-card:nth-child(2){border-top-color:var(--amber)}.dimension-card:nth-child(3){border-top-color:var(--blue)}
.dimension-code{display:flex;justify-content:space-between;gap:10px;color:var(--teal);font:800 10px/1.2 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.07em;text-transform:uppercase}.dimension-card:nth-child(2) .dimension-code{color:var(--amber)}.dimension-card:nth-child(3) .dimension-code{color:var(--blue)}
.dimension-card h3{margin:18px 0 8px;font-size:17px}.dimension-answer{margin:0;color:var(--ink-soft);font-size:12.5px}.dimension-answer strong{color:var(--ink);font:800 22px/1.1 "IBM Plex Mono","Cascadia Code",monospace}.dimension-list{display:grid;gap:5px;margin:15px 0}.dimension-row{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:10px;padding:7px 8px;border:1px solid #ddd9cf;background:#f6f3eb}.dimension-row span{color:var(--muted);font-size:10.5px}.dimension-row b{font:750 10.5px/1.3 "IBM Plex Mono","Cascadia Code",monospace;text-align:right}.dimension-boundary{margin:auto 0 0;color:var(--muted);font-size:11.5px}.dimension-link{display:inline-block;margin-top:13px;padding:8px 10px;border:1px solid var(--rail);background:var(--rail);color:#fff;text-decoration:none;font-size:11px;font-weight:800}.dimension-link:hover{background:var(--teal)}
.position-leaders{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:4px;margin:0 0 15px}.position-leaders a{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:6px;padding:6px 7px;border:1px solid #ddd9cf;background:#f6f3eb;text-decoration:none}.position-leaders span{min-width:0;overflow-wrap:anywhere;color:var(--ink-soft);font-size:9.5px}.position-leaders b{font:750 9.5px/1.35 "IBM Plex Mono","Cascadia Code",monospace;white-space:nowrap}
.genome-launchers{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:9px}
.genome-link{display:grid;grid-template-columns:auto 1fr;gap:11px;align-items:start;min-height:122px;padding:16px;border:1px solid #98a5a4;background:var(--panel);color:var(--ink);text-decoration:none;transition:transform .18s ease,box-shadow .18s ease,border-color .18s ease}
.genome-link:hover{transform:translateY(-2px);border-color:var(--teal);box-shadow:0 10px 24px rgba(19,45,52,.09)}
.launcher-index{display:grid;place-items:center;width:29px;height:29px;border-radius:50%;background:var(--rail);color:#fff;font:700 10px/1 "IBM Plex Mono","Cascadia Code",monospace}.genome-link strong,.genome-link small{display:block}.genome-link strong{overflow-wrap:anywhere;font:800 13px/1.25 "IBM Plex Mono","Cascadia Code",monospace}.genome-link small{margin-top:6px;color:var(--muted);font-size:10.5px}.launcher-action{grid-column:2;margin-top:auto;color:var(--teal);font-size:11px;font-weight:800}
.table-shell{max-width:100%;box-shadow:var(--shadow)}.scroll-cue{display:flex;justify-content:space-between;gap:12px;padding:9px 12px;border-bottom:1px solid var(--line);background:#e9e7df;color:var(--muted);font-size:11px}.table-scroll{width:100%;overflow-x:auto;overscroll-behavior-inline:contain}.table-scroll:focus-visible{outline-offset:-3px}
table{width:100%;min-width:1020px;border-collapse:collapse;font-size:12px}th,td{padding:11px 10px;text-align:left;vertical-align:middle;border-bottom:1px solid #e2dfd6}thead th{position:sticky;top:0;z-index:2;background:#dfe5e3;color:#2f494f;font:800 9.5px/1.3 "IBM Plex Mono","Cascadia Code",monospace;letter-spacing:.03em}tbody tr:last-child>*{border-bottom:0}tbody tr:hover>*{background:#f4f5ef}th:first-child,td:first-child{position:sticky;left:0;z-index:1}thead th:first-child{z-index:3}tbody th:first-child{background:var(--panel);box-shadow:7px 0 12px -13px #000}.sample-cell{min-width:166px}.sample-cell strong,.sample-cell small,td small{display:block}.sample-cell small,td small{margin-top:3px;color:var(--muted);font-size:10px}.numeric{text-align:right;font-family:"IBM Plex Mono","Cascadia Code",monospace;font-variant-numeric:tabular-nums}.table-action{display:inline-block;padding:7px 9px;border:1px solid var(--rail);background:var(--rail);color:#fff;text-decoration:none;white-space:nowrap;font-weight:750}.table-action:hover{background:var(--teal)}
.claim-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}.claim-card{min-height:185px;padding:20px}.claim-card:nth-child(1){border-top:4px solid var(--teal)}.claim-card:nth-child(2){border-top:4px solid var(--amber)}.claim-card:nth-child(3){border-top:4px solid var(--brick)}.claim-index{color:var(--muted);font:800 10px/1 "IBM Plex Mono","Cascadia Code",monospace}.claim-card h3{margin:25px 0 9px;font-size:16px}.claim-card p{margin:0;color:var(--ink-soft);font-size:12.5px}
.layer-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:0;border:1px solid var(--line)}.layer-card{min-height:168px;padding:19px;border:0;border-right:1px solid var(--line)}.layer-card:last-child{border-right:0}.layer-code{color:var(--teal);font:850 11px/1 "IBM Plex Mono","Cascadia Code",monospace}.layer-card h3{margin:17px 0 7px;font-size:14px}.layer-card p{margin:0;color:var(--muted);font-size:11.5px}.layer-card.pending{background:#f8eee9}.layer-card.pending .layer-code{color:var(--brick)}
.qc-note{display:grid;grid-template-columns:auto 1fr;gap:14px;align-items:start;margin-top:10px;padding:15px 18px;border-left:5px solid var(--blue);background:#e8edf0}.qc-note strong{font:800 11px/1.3 "IBM Plex Mono","Cascadia Code",monospace}.qc-note p{margin:0;color:var(--ink-soft);font-size:12px}
.evidence-drawer{margin-top:34px}.evidence-drawer summary{cursor:pointer;padding:14px 17px;font-weight:800}.evidence-drawer[open] summary{border-bottom:1px solid var(--line)}.evidence-inner{padding:16px 18px}.evidence-inner p{margin:0 0 12px;color:var(--muted);font-size:11.5px}.evidence-list{margin:0;padding-left:20px}.evidence-list li+li{margin-top:9px}.evidence-list a{font-size:11.5px}.evidence-list code{display:block;margin-top:2px;overflow-wrap:anywhere;color:var(--muted);font-size:9.5px}
footer{display:flex;justify-content:space-between;gap:20px;margin-top:18px;padding-top:13px;border-top:1px solid var(--line);color:var(--muted);font-size:10.5px}.mono{font-family:"IBM Plex Mono","Cascadia Code",monospace;overflow-wrap:anywhere}
@media(max-width:1050px){.metric-grid{grid-template-columns:repeat(2,1fr)}.genome-launchers{grid-template-columns:repeat(3,1fr)}.topology-board{grid-template-columns:1fr}.dimension-grid{grid-template-columns:1fr}.dimension-card{min-height:0}.layer-grid{grid-template-columns:repeat(2,1fr)}.layer-card:nth-child(2){border-right:0}.layer-card:nth-child(-n+2){border-bottom:1px solid var(--line)}.signature-grid,.hcc-metric-grid{grid-template-columns:repeat(2,1fr)}.comparison-matrix-layout,.hcc-analysis-grid{grid-template-columns:1fr}}
@media(max-width:760px){.shell{padding:12px 11px 36px}.hero{grid-template-columns:1fr}.hero-aside{border-left:0;border-top:1px solid rgba(255,255,255,.18)}.sensitivity-banner{grid-template-columns:1fr}.sensitivity-metrics{grid-template-columns:1fr 1fr}.section-head{display:block}.section-note{margin-top:7px}.sample-bars,.claim-grid,.comparison-scope-grid{grid-template-columns:1fr}.genome-launchers{grid-template-columns:repeat(2,1fr)}footer{display:block}footer span{display:block;margin-top:5px}.comparison-panel-head{display:block}.comparison-panel-head p{margin-top:7px}.comparison-profile-row{grid-template-columns:1fr;gap:6px;padding:9px}.comparison-profile-n{text-align:left}.hcc-verdict{grid-template-columns:1fr}.claim-boundary-box{grid-template-columns:1fr}}
@media(max-width:760px){.interpretation-board{grid-template-columns:1fr}.distribution-card>p{min-height:0}}
@media(max-width:500px){.metric-grid,.genome-launchers,.layer-grid,.signature-grid,.hcc-metric-grid,.comparison-controls{grid-template-columns:1fr}.hero-main{padding:28px 20px 58px}.hero-aside{padding:24px 20px 48px}.metric-card{min-height:0}.position-leaders{grid-template-columns:1fr}.layer-card{min-height:0;border-right:0;border-bottom:1px solid var(--line)}.layer-card:last-child{border-bottom:0}.scroll-cue span:last-child{display:none}.topology-main{padding:18px}.topology-legend{padding:15px}.distribution-card{padding:18px}.distribution-head{display:block}.distribution-head>b{display:block;margin-top:7px}.distribution-legend li{grid-template-columns:9px minmax(0,1fr) auto}.hcc-verdict dl{grid-template-columns:1fr}.comparison-profile-panel,.matrix-panel,.hcc-table-panel,.confusion-panel{padding:14px}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}.genome-link{transition:none}}
@media print{body{background:#fff}.shell{width:100%;padding:0}.hero,.table-shell{box-shadow:none}.genome-link{break-inside:avoid}.table-scroll{overflow:visible}table{min-width:0;font-size:8px}.skip-link,.scroll-cue,.table-action{display:none}}
</style>
</head>
<body>
<a class="skip-link" href="#main-content">跳至主要內容</a>
<main class="shell" id="main-content">
  <header class="hero">
    <div class="hero-main">
      <p class="eyebrow"><span class="canonical">Canonical / 7 of 7 PASS</span><span>chr1–22 全基因範圍</span></p>
      <h1>分層重建<br>全基因指揮中心</h1>
      <p class="lede">以 LongPhase-S recalibrated FILTER=PASS 為唯一主骨幹，從 7-dataset 總覽進入 GRCh38 座標比例分布、current-v5 family-specific read-AF 描述性第一順位，以及 mutation-state graph morphology（僅 clone/subclone-compatible，不是 clone census）；每個數字保留明確 grain、分母與候選完整性。</p>
    </div>
    <aside class="hero-aside" aria-label="Canonical authority">
      <p class="authority-label">Machine-bound authority</p>
      <dl class="authority-list">
        <div><dt class="authority-label">Run</dt><dd class="authority-value">$run_id</dd></div>
        <div><dt class="authority-label">Verified</dt><dd class="authority-value">$verified_at</dd></div>
        <div><dt class="authority-label">Summary SHA-256</dt><dd class="authority-value">$summary_short</dd></div>
      </dl>
    </aside>
  </header>

  <aside class="sensitivity-banner" aria-labelledby="sensitivity-title">
    <div><p class="section-kicker">Backbone sensitivity · 不可混合分母</p><h2 id="sensitivity-title">$backbone_verdict</h2><p>LongPhase-S recalibrated FILTER=PASS 仍是 canonical 主結果；ClairS FILTER=PASS 僅作 sensitivity。低 overlap 指標表示所有解讀都必須標示 backbone-sensitive。</p></div>
    <div class="sensitivity-metrics" aria-label="Backbone sensitivity minimum concordance metrics">
      <div class="sensitivity-metric"><span>Retained-position Jaccard · min</span><strong>$retained_jaccard</strong></div>
      <div class="sensitivity-metric"><span>Primary-unit-key Jaccard · min</span><strong>$primary_jaccard</strong></div>
      <div class="sensitivity-metric"><span>Shared topology digest · min</span><strong>$topology_concordance</strong></div>
      <div class="sensitivity-metric"><span>Max absolute proportion delta</span><strong>$max_delta_pp pp</strong></div>
    </div>
  </aside>

  <section class="section" aria-labelledby="pulse-title">
    <div class="section-head"><div><p class="section-kicker">Cohort pulse</p><h2 id="pulse-title">先鎖定 W 與候選完整性</h2></div><p class="section-note">W_tree 是 retained multilocus regions；C/Topo 只以具有 mutation-bearing HP1/HP2 primary units 的 W_primary 為分母。</p></div>
    <div class="metric-grid">
      <article class="metric-card"><span class="metric-code">W_tree · region</span><strong class="metric-number">$w_tree</strong><p>chr1–22 retained regional groups。</p><small class="metric-sub">no-primary $no_primary · $no_primary_rate%</small></article>
      <article class="metric-card"><span class="metric-code">W_primary · region</span><strong class="metric-number">$w_primary</strong><p>至少一個 mutation-bearing HP1/HP2 primary unit。</p><small class="metric-sub">C / Topo denominator</small></article>
      <article class="metric-card"><span class="metric-code">Complete · region</span><strong class="metric-number">$complete</strong><p>所有 primary units 均 non-capped、full-pass、candidate-complete。</p><small class="metric-sub">$complete_rate% of W_primary</small></article>
      <article class="metric-card"><span class="metric-code">Incomplete · region</span><strong class="metric-number">$incomplete</strong><p>任一 primary unit 未完成時，C 與 Topo 必須保持 unavailable。</p><small class="metric-sub">$incomplete_rate% of W_primary</small></article>
    </div>
  </section>

  <section class="section" aria-labelledby="dimensions-title">
    <div class="section-head"><div><p class="section-kicker">Three reading dimensions</p><h2 id="dimensions-title">先用三個問題定位每個數字</h2></div><p class="section-note">舊頁的三軸教學結構保留；所有名詞與分母改綁 canonical v5，不沿用 clone-first 類別或 archive 數字。</p></div>
    <div class="dimension-grid">
      <article class="dimension-card"><div class="dimension-code"><span>01 · 拓樸型態</span><span>Topology shape</span></div><h3>分子累積形狀有幾種？</h3><p class="dimension-answer"><strong>$shape_unique</strong> 個 complete regions 的 Topo=1。</p><div class="dimension-list"><div class="dimension-row"><span>形狀唯一 / complete</span><b>$shape_unique / $complete · $shape_unique_rate%</b></div><div class="dimension-row"><span>多種形狀相容</span><b>$multi_shape · $multi_shape_rate%</b></div></div><p class="dimension-boundary">Topo 比較無節點標籤的 regional mutation-state shape；不等於 biological ancestry 或真實時間順序。</p></article>
      <article class="dimension-card"><div class="dimension-code"><span>02 · 可辨識度</span><span>Determinacy</span></div><h3>目前能唯一辨識到哪一層？</h3><p class="dimension-answer"><strong>$exact_unique</strong> 個 regions 可辨識到 exact candidate。</p><div class="dimension-list"><div class="dimension-row"><span>Exact 唯一</span><b>$exact_unique</b></div><div class="dimension-row"><span>只到 shape 唯一</span><b>$shape_only</b></div><div class="dimension-row"><span>Multi-shape</span><b>$multi_shape</b></div><div class="dimension-row"><span>尚未評估</span><b>$incomplete</b></div></div><p class="dimension-boundary">Incomplete 是 candidate set 未完整，所以 C／Topo 不可評估；不是 0，也不是已證明多解。</p></article>
      <article class="dimension-card"><div class="dimension-code"><span>03 · 基因體位置</span><span>Genome position</span></div><h3>這些 regions 位於哪裡？</h3><p class="dimension-answer"><strong>7 × chr1–22</strong> 的 canonical 全基因巡覽。</p><div class="dimension-list"><div class="dimension-row"><span>W_tree regions</span><b>$w_tree</b></div><div class="dimension-row"><span>Dataset pages</span><b>7</b></div><div class="dimension-row"><span>Autosomes / page</span><b>22</b></div></div><div class="position-leaders" aria-label="各 dataset 的 W primary count leader chromosome">$position_leaders</div><p class="dimension-boundary">內頁 ideogram 依 GRCh38 bp 長度與 region midpoint 定位；上列仍是各 dataset 的 raw W_primary count leader，未除以 chromosome 長度或輸入 sSNV density，因此只供巡覽，不作 hotspot／enrichment 判定。</p><a class="dimension-link" href="#launch-title">選擇 dataset 看位置</a></article>
    </div>
  </section>

  <section class="section" aria-labelledby="topology-title">
    <div class="section-head"><div><p class="section-kicker">Candidate topology</p><h2 id="topology-title">Exact 組合與拓撲形狀分開讀</h2></div><p class="section-note">C_region 是 exact candidate-tree 組合乘積；Topo_region 是 exact unlabeled-shape 組合乘積。Incomplete 不以 0 代替。</p></div>
    <div class="topology-board">
      <article class="topology-main"><h3>7-dataset aggregate · W_primary $w_primary</h3><p>四類互斥且守恆；不可能狀態 C=1 / Topo&gt;1 為 0。HCC1395 與 DORADO 是同一 biological sample 的兩個 dataset。</p><div class="topology-bar" role="img" aria-label="7-dataset 候選拓撲四類">$aggregate_segments</div><div class="bar-scale"><span>0</span><span>W_primary composition</span><span>100%</span></div></article>
      <aside class="topology-legend" aria-label="Topology legend"><ul>$topology_legend</ul></aside>
    </div>
    <div class="sample-bars" aria-label="7 dataset topology distributions">$sample_topology_rows</div>
  </section>

  <section class="section" aria-labelledby="read-af-morphology-title">
    <div class="section-head"><div><p class="section-kicker">Current-v5 interpretation layers</p><h2 id="read-af-morphology-title">順位與形態是兩個不同問題</h2></div><p class="section-note">這是 7 datasets／6 biological samples 的 dataset aggregate（HCC1395 含兩個 dataset），不是 biological-sample 加權。兩張圖都以 W_primary 為分母且各自互斥守恆。</p></div>
    $cohort_read_af_morphology
  </section>

$sample_comparison_dashboard

  <section class="section" aria-labelledby="launch-title">
    <div class="section-head"><div><p class="section-kicker">Genome launchpad</p><h2 id="launch-title">進入 GRCh38 chr1–22 全基因巡覽</h2></div><p class="section-note">每頁都含座標比例 ideogram、read-AF 第一順位與 clone-compatible morphology 等七組 sample-wide 觀察，並由 current summary、sample region-view、read-AF sidecar、renderer SHA 與 UI contract 綁定；stale 頁面不會進入 index。</p></div>
    <nav class="genome-launchers" aria-label="開啟各 dataset 全基因頁">$genome_launchers</nav>
  </section>

  <section class="section" id="cohort-table" aria-labelledby="table-title">
    <div class="section-head"><div><p class="section-kicker">Cohort ledger</p><h2 id="table-title">7 dataset 簡潔盤點</h2></div><p class="section-note">LPS PASS 是全 contig tree-input records；retained sSNV 與 W metrics 是 chr1–22。不同 grain 不直接互除。</p></div>
    <div class="table-shell"><div class="scroll-cue"><span>表格可水平捲動；樣本欄固定。</span><span>Canonical main only</span></div><div class="table-scroll" role="region" aria-label="Canonical cohort table" tabindex="0"><table>
      <caption class="sr-only">Canonical v5 七個 dataset 的輸入、retained sSNV、W_tree、W_primary、候選完整性與 primary units</caption>
      <thead><tr><th scope="col">Dataset</th><th scope="col" class="numeric">LPS PASS<br>records</th><th scope="col" class="numeric">Retained<br>sSNV</th><th scope="col" class="numeric">W_tree<br>regions</th><th scope="col" class="numeric">W_primary<br>regions</th><th scope="col" class="numeric">Complete<br>regions</th><th scope="col" class="numeric">Incomplete<br>regions</th><th scope="col" class="numeric">Primary<br>units</th><th scope="col">巡覽</th></tr></thead>
      <tbody>$cohort_rows</tbody>
    </table></div></div>
  </section>

  <section class="section" aria-labelledby="boundary-title">
    <div class="section-head"><div><p class="section-kicker">Claim boundary</p><h2 id="boundary-title">觀察、推論、限制分開</h2></div><p class="section-note">這是 regional mutation-state candidate sets；不是細胞層級真相或已確認的全腫瘤演化關係。</p></div>
    <div class="claim-grid">
      <article class="claim-card"><span class="claim-index">01 / 回答</span><h3>哪些區域候選與資料相容？</h3><p>Analysis producer 對完整 minimal candidate universe 枚舉與計數；頁面可能只展開 structural stored subset，並另顯示 exhaustive read-AF ranking 的 preview 與 co-top representatives。</p></article>
      <article class="claim-card"><span class="claim-index">02 / 證據</span><h3>read-level sSNV 共現</h3><p>LongPhase-S recalibrated PASS 提供 tree backbone；L0 family partition 與 L1 read-state constraints 承重，CN 只作 post-tree context。</p></article>
      <article class="claim-card"><span class="claim-index">03 / 限制</span><h3>不可越過資料識別度</h3><p>單一 bulk 的 regional candidates 不能單獨證明細胞比例、真實 ancestry 或正交生物學確認；多解與 incomplete 都是應保留的結果。</p></article>
    </div>
  </section>

  <section class="section" aria-labelledby="layers-title">
    <div class="section-head"><div><p class="section-kicker">Evidence stack</p><h2 id="layers-title">L0 → L3 的承重順序</h2></div><p class="section-note">後續 auxiliary layer 不得回頭重寫 read-level observation，也不得偷偷選樹。</p></div>
    <div class="layer-grid">
      <article class="layer-card"><span class="layer-code">L0</span><h3>HP family partition</h3><p>HP1/HP2 mutation-bearing units 是 primary；H3、unphased 與 reference-only 都排除於 primary denominator。</p></article>
      <article class="layer-card"><span class="layer-code">L1</span><h3>Candidate enumeration</h3><p>以 read-state constraints 枚舉 minimal candidates；capped 代表候選集合未完成。</p></article>
      <article class="layer-card"><span class="layer-code">L2</span><h3>Copy-number context</h3><p>只對 recurrence-like outcome 提供 post-tree annotation；不重新定義 observed state。</p></article>
      <article class="layer-card pending"><span class="layer-code">L3 · NOT EVALUATED</span><h3>Bounded methyl auxiliary</h3><p>本 canonical run 未評估；未來僅限 negative screen / residual flag，禁止排序或確認候選。</p></article>
    </div>
    <aside class="qc-note"><strong>PS / QC ONLY</strong><p>Phase-set 只保留作 phase-block 品質脈絡；不是 topology edge，也不是 lineage label。全 cohort 有 $mixed_ps 個 multi-PS regions。</p></aside>
  </section>

  <details class="evidence-drawer">
    <summary>機器證據與原始 JSON（預設收合）</summary>
    <div class="evidence-inner"><p>下列連結只供 provenance readback；W／structural C／Topo 來自 hash-verified canonical machine summary，read-AF selection／morphology 來自另行 hash 綁定且守恆驗證的 current-v5 sidecar index。</p><ul class="evidence-list">$evidence_links</ul></div>
  </details>

  <footer><span>Canonical main · LongPhase-S recalibrated FILTER=PASS · chr1–22 · 7/7 PASS</span><span class="mono">Page generated $generated_at</span></footer>
</main>
<script>
(() => {
  const source = document.getElementById('sample-comparison-data');
  if (!source) return;
  const data = JSON.parse(source.textContent);
  let dimension = 'structural';
  let selectedPair = ['HCC1395', 'HCC1395_DORADO'];
  const safe = value => String(value).replace(/[&<>"']/g, character => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[character]));
  const pct = value => (Number(value) * 100).toFixed(1) + '%';
  const pairKey = (left, right) => [left, right].sort().join('|');
  const findPair = (left, right) => data.pairwise.records.find(record => pairKey(record.left, record.right) === pairKey(left, right));
  const shortName = sample => sample === 'HCC1395_DORADO' ? 'DORADO' : sample;

  function renderKey() {
    const contract = data.dimensions[dimension];
    document.getElementById('comparison-key').innerHTML = contract.classes.map(item =>
      '<span class="comparison-key-item"><i class="legend-swatch segment-' + safe(item.class_name) + '" aria-hidden="true"></i><b>' + safe(item.label) + '</b></span>'
    ).join('');
  }

  function renderProfiles() {
    const contract = data.dimensions[dimension];
    document.getElementById('comparison-profile-title').textContent = contract.label + ' · 逐 dataset';
    document.getElementById('comparison-profile-chart').innerHTML = data.datasets.map(record => {
      const values = record.dimensions[dimension];
      const technical = record.dataset === 'HCC1395' || record.dataset === 'HCC1395_DORADO';
      const role = record.dataset === 'HCC1395' ? '5kHz canonical' : 'Dorado technical pair';
      const segments = contract.classes.map(item =>
        '<span class="bar-segment segment-' + safe(item.class_name) + '" style="width:' + (Number(values.proportions[item.key]) * 100).toFixed(6) + '%" title="' + safe(item.label) + ' · ' + Number(values.counts[item.key]).toLocaleString() + ' · ' + pct(values.proportions[item.key]) + '"></span>'
      ).join('');
      return '<div class="comparison-profile-row' + (technical ? ' technical-pair' : '') + '" data-dataset="' + safe(record.dataset) + '"><div class="comparison-profile-name"><a href="' + safe(record.dataset) + '.html"><strong>' + safe(record.dataset) + '</strong></a>' + (technical ? '<small>' + safe(role) + '</small>' : '') + '</div><div class="topology-bar" role="img" aria-label="' + safe(record.dataset) + ' ' + safe(contract.label) + ' composition">' + segments + '</div><b class="comparison-profile-n">' + Number(record.W_primary).toLocaleString() + '<small>W_primary</small></b></div>';
    }).join('');
  }

  function renderLedger() {
    const contract = data.dimensions[dimension];
    document.getElementById('comparison-ledger-head').innerHTML = '<tr><th scope="col">Dataset</th><th scope="col" class="numeric">W_primary</th>' + contract.classes.map(item => '<th scope="col" class="numeric">' + safe(item.label) + '</th>').join('') + '</tr>';
    document.getElementById('comparison-ledger-body').innerHTML = data.datasets.map(record => {
      const values = record.dimensions[dimension];
      const cells = contract.classes.map(item => '<td class="numeric" data-label="' + safe(item.label) + '"><strong>' + Number(values.counts[item.key]).toLocaleString() + '</strong><small>' + pct(values.proportions[item.key]) + '</small></td>').join('');
      return '<tr><th scope="row" class="sample-cell"><a href="' + safe(record.dataset) + '.html"><strong>' + safe(record.dataset) + '</strong></a><small>' + safe(record.biological_id) + ' biological ID</small></th><td class="numeric" data-label="W_primary">' + Number(record.W_primary).toLocaleString() + '</td>' + cells + '</tr>';
    }).join('');
  }

  function renderPairInspector() {
    const record = findPair(selectedPair[0], selectedPair[1]);
    if (!record) return;
    const left = data.datasets.find(item => item.dataset === selectedPair[0]);
    const right = data.datasets.find(item => item.dataset === selectedPair[1]);
    const metric = record.dimensions[dimension];
    const technical = record.same_biological_sample === true;
    const deltas = data.dimensions[dimension].classes.map(item => ({
      label: item.label,
      delta: Number(right.dimensions[dimension].proportions[item.key]) - Number(left.dimensions[dimension].proportions[item.key])
    })).sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta)).slice(0, 3);
    document.getElementById('pair-inspector').innerHTML = '<span class="pair-chip">' + (technical ? 'Same biological sample · technical pair' : 'Different biological samples · composition only') + '</span><h3>' + safe(selectedPair[0]) + ' × ' + safe(selectedPair[1]) + '</h3><strong class="pair-distance">' + Number(metric.tvd).toFixed(3) + '<small>TVD · rank ' + Number(metric.tvd_rank_among_21_pairs) + '/21</small></strong><div class="pair-delta-list">' + deltas.map(item => '<div><span>' + safe(item.label) + '</span><b>' + (item.delta >= 0 ? '+' : '') + (item.delta * 100).toFixed(1) + ' pp</b></div>').join('') + '</div><p>' + (technical ? '此格只比較類別比例；下方 exact-region panel 才檢查相同區域是否得到相同標籤。' : '不同 biological samples 只比較 profile；不以 exact-coordinate overlap 當一般跨細胞株品質指標。') + '</p>';
  }

  function renderMatrix() {
    const order = data.pairwise.matrices.order;
    const matrix = data.pairwise.matrices.tvd[dimension];
    document.getElementById('comparison-matrix-title').textContent = data.dimensions[dimension].label + ' · TVD';
    let html = '<table class="distance-matrix"><caption class="sr-only">7 datasets pairwise composition TVD</caption><thead><tr><th scope="col">TVD</th>' + order.map(sample => '<th scope="col">' + safe(shortName(sample)) + '</th>').join('') + '</tr></thead><tbody>';
    order.forEach(left => {
      html += '<tr><th scope="row">' + safe(shortName(left)) + '</th>';
      order.forEach(right => {
        const value = Number(matrix[left][right]);
        if (left === right) {
          html += '<td class="matrix-diagonal"><span aria-label="same dataset">0</span></td>';
          return;
        }
        const isHcc = pairKey(left, right) === pairKey('HCC1395', 'HCC1395_DORADO');
        const isSelected = pairKey(left, right) === pairKey(selectedPair[0], selectedPair[1]);
        const alpha = 0.08 + Math.min(value / 0.5, 1) * 0.64;
        html += '<td><button type="button" class="matrix-cell' + (isHcc ? ' is-hcc' : '') + (isSelected ? ' is-selected' : '') + '" data-left="' + safe(left) + '" data-right="' + safe(right) + '" style="background:rgba(22,139,141,' + alpha.toFixed(4) + ')" aria-label="' + safe(left) + ' 與 ' + safe(right) + '，TVD ' + value.toFixed(3) + '">' + value.toFixed(3) + '</button></td>';
      });
      html += '</tr>';
    });
    html += '</tbody></table>';
    const container = document.getElementById('comparison-matrix');
    container.innerHTML = html;
    container.querySelectorAll('.matrix-cell').forEach(button => button.addEventListener('click', () => {
      selectedPair = [button.dataset.left, button.dataset.right];
      renderMatrix();
      renderPairInspector();
    }));
  }

  function renderHccConfusion() {
    const contract = data.dimensions[dimension];
    const matrix = data.hcc.matched_primary_dimensions[dimension].matrix;
    document.getElementById('hcc-confusion-title').textContent = contract.label + ' · 哪些類別互相轉換？';
    let html = '<table class="confusion-table"><caption>Rows = HCC1395 5kHz；columns = HCC1395_DORADO</caption><thead><tr><th scope="col">5kHz ↓ / DORADO →</th>' + contract.classes.map(item => '<th scope="col">' + safe(item.label) + '</th>').join('') + '</tr></thead><tbody>';
    contract.classes.forEach(left => {
      html += '<tr><th scope="row">' + safe(left.label) + '</th>';
      contract.classes.forEach(right => {
        html += '<td class="numeric' + (left.key === right.key ? ' confusion-diagonal' : '') + '">' + Number(matrix[left.key][right.key]).toLocaleString() + '</td>';
      });
      html += '</tr>';
    });
    document.getElementById('hcc-confusion').innerHTML = html + '</tbody></table>';
  }

  function renderAll() {
    document.querySelectorAll('.comparison-tab').forEach(button => {
      const active = button.dataset.comparisonDimension === dimension;
      button.classList.toggle('is-active', active);
      button.setAttribute('aria-pressed', active ? 'true' : 'false');
    });
    renderKey();
    renderProfiles();
    renderLedger();
    renderMatrix();
    renderPairInspector();
    renderHccConfusion();
  }

  document.querySelectorAll('.comparison-tab').forEach(button => button.addEventListener('click', () => {
    dimension = button.dataset.comparisonDimension;
    renderAll();
  }));
  renderAll();
})();
</script>
</body>
</html>
""")
    return template.substitute(
        summary_sha256=escaped(authority["summary_sha256"]),
        comparison_sha256=escaped(authority["summary"]["backbone_comparison"]["sha256"]),
        read_af_index_sha256=escaped(authority["read_af_index_sha256"]),
        sample_comparison_sha256=escaped(authority["sample_comparison_sha256"]),
        summary_short=escaped(authority["summary_sha256"][:16] + "…"),
        run_id=escaped(run_id),
        backbone_verdict=escaped(str(comparison["verdict"]).replace("_", " ").upper()),
        retained_jaccard=f'{float(comparison["min_retained_position_jaccard"]):.3f}',
        primary_jaccard=f'{float(comparison["min_primary_unit_key_jaccard"]):.3f}',
        topology_concordance=f'{float(comparison["min_shared_topology_digest_concordance"]):.3f}',
        max_delta_pp=f'{float(comparison["max_abs_delta_pp"]):.2f}',
        verified_at=escaped(verified_at.isoformat(timespec="minutes")),
        generated_at=escaped(generated_at),
        w_tree=f'{int(aggregate["W_tree"]):,}',
        w_primary=f'{int(aggregate["W_primary"]):,}',
        no_primary=f'{int(aggregate["no_primary_lineage"]):,}',
        no_primary_rate=f'{no_primary_rate:.1f}',
        complete=f'{int(aggregate["complete_regions"]):,}',
        incomplete=f'{int(aggregate["incomplete_regions"]):,}',
        complete_rate=f'{complete_rate:.1f}',
        incomplete_rate=f'{incomplete_rate:.1f}',
        exact_unique=f'{exact_unique:,}',
        shape_only=f'{shape_unique_exact_multiple:,}',
        multi_shape=f'{multi_shape:,}',
        shape_unique=f'{shape_unique:,}',
        shape_unique_rate=f'{percentage(shape_unique, int(aggregate["complete_regions"])):.1f}',
        multi_shape_rate=f'{percentage(multi_shape, int(aggregate["complete_regions"])):.1f}',
        mixed_ps=f'{int(aggregate["mixed_PS_regions"]):,}',
        aggregate_segments=topology_segments(topology, int(aggregate["W_primary"])),
        topology_legend=topology_legend(aggregate),
        sample_topology_rows=sample_topology_rows(rows),
        cohort_read_af_morphology=cohort_read_af_morphology(rows),
        sample_comparison_dashboard=sample_comparison_dashboard(authority["sample_comparison"]),
        position_leaders=sample_position_leaders(rows),
        genome_launchers=genome_launchers(rows),
        cohort_rows=cohort_rows(rows),
        evidence_links=evidence_links(authority, rows),
    )


def write_index_atomic(document: str) -> Path:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    output = OUTDIR / "index.html"
    staging = OUTDIR / f".index.html.stage.{os.getpid()}"
    with staging.open("w", encoding="utf-8") as handle:
        handle.write(document)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(staging, output)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--index-only",
        action="store_true",
        help="Verify seven existing hash-bound sample pages and rebuild only index.html.",
    )
    args = parser.parse_args()
    try:
        authority = load_authority()
        rows = collect_rows(authority, build_samples=not args.index_only)
        index_path = write_index_atomic(build_index(authority, rows))
    except AuthorityError as exc:
        print(f"FAIL CLOSED: {exc}", file=sys.stderr)
        return 2
    print(f"OK canonical index + 7/7 hash-bound sample pages -> {index_path}")
    print(f"SUMMARY_SHA256={authority['summary_sha256']}")
    print(f"CANONICAL_RUN={authority['run_root']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
