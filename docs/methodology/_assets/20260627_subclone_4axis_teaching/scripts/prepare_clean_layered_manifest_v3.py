#!/usr/bin/env python3
"""Build a strict layered v3 manifest from completed production coordinate sidecars."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import validate_layered_v3_inputs as contract


MANIFEST_SCHEMA_REF = (
    "InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/"
    "schemas/layered_input_manifest_v3.schema.json"
)


def _require_exact_source_set(document: dict[str, Any], label: str) -> dict[str, dict[str, Any]]:
    samples = document.get("samples")
    if not isinstance(samples, list):
        raise contract.ContractError("E_SCHEMA_TYPE", f"{label}.samples must be an array")
    result: dict[str, dict[str, Any]] = {}
    for item in samples:
        if not isinstance(item, dict) or not isinstance(item.get("sample"), str):
            raise contract.ContractError("E_REQUIRED_FIELD", f"{label} sample entry lacks sample ID")
        sample = item["sample"]
        if sample in result:
            raise contract.ContractError("E_SAMPLE_DUPLICATE", f"duplicate {label} sample: {sample}")
        result[sample] = item
    if set(result) != set(contract.EXPECTED_DATASETS):
        raise contract.ContractError(
            "E_DATASET_SET_MISMATCH",
            f"{label} must contain the canonical seven datasets",
            expected=sorted(contract.EXPECTED_DATASETS),
            observed=sorted(result),
        )
    return result


def _path(value: Any, label: str, sample: str) -> Path:
    if not isinstance(value, str) or not value or not Path(value).is_absolute():
        raise contract.ContractError(
            "E_REQUIRED_ARTIFACT", f"{sample} {label} must be an absolute path",
            exit_code=3, stage="artifact_identity", sample=sample,
        )
    return Path(value)


def _optional_artifact(value: Any, label: str, sample: str) -> dict[str, Any] | None:
    if value is None:
        return None
    return contract.artifact(_path(value, label, sample))


def _copy_number(base: dict[str, Any], sample: str) -> dict[str, Any]:
    cn_bed = base.get("cn_bed")
    availability = "measured" if cn_bed is not None else "unavailable"
    if availability == "unavailable" and any(
        base.get(key) is not None for key in ("cn_int_gain", "cn_int_loss", "integration_json")
    ):
        raise contract.ContractError("E_CN_CONTRACT", f"{sample} unavailable CN has measured artifacts")
    source = base.get("cn_source")
    semantics = base.get("cn_semantics")
    if not isinstance(source, str) or not source or not isinstance(semantics, str) or not semantics:
        raise contract.ContractError("E_CN_CONTRACT", f"{sample} CN source/semantics are not explicit")
    if availability == "unavailable" and (source != "unavailable" or semantics != "missing; never interpreted neutral"):
        raise contract.ContractError("E_CN_CONTRACT", f"{sample} unavailable CN semantics are ambiguous")
    return {
        "availability": availability,
        "source": source,
        "semantics": semantics,
        "coordinate_system": "0_based_half_open" if availability == "measured" else None,
        "unlisted_position_semantics": "neutral" if availability == "measured" else "unavailable",
        "allowed_states": ["gain", "loss", "loh", "neutral"] if availability == "measured" else [],
        "overlap_policy": "forbid" if availability == "measured" else "not_applicable",
        "reason": None if availability == "measured" else "No reviewed measured CN source is available",
        "cn_bed": _optional_artifact(cn_bed, "cn_bed", sample),
        "cn_int_gain": _optional_artifact(base.get("cn_int_gain"), "cn_int_gain", sample),
        "cn_int_loss": _optional_artifact(base.get("cn_int_loss"), "cn_int_loss", sample),
        "integration_json": _optional_artifact(base.get("integration_json"), "integration_json", sample),
    }


def _build_sample(
    base: dict[str, Any], source: dict[str, Any], production_root: Path, tree_input_contract: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    sample = base["sample"]
    sample_dir = production_root / "samples" / sample
    bam = _path(source.get("tumor_bam"), "tumor_bam", sample)
    bam_index = Path(f"{bam}.bai")
    storage = contract.storage_identity_v1(bam, bam_index)
    raw = contract.indexed_artifact(_path(source.get("caller_raw_vcf"), "caller_raw_vcf", sample))
    normalized = contract.indexed_artifact(sample_dir / f"{sample}.clairs.raw_all.normalized.vcf.gz")
    caller_pass_baseline = contract.indexed_artifact(
        _path(source.get("caller_pass_vcf"), "caller_pass_vcf", sample)
    )
    recal_all = contract.indexed_artifact(sample_dir / f"{sample}.longphase_s.recalibrated.all.vcf.gz")
    recal_pass = contract.indexed_artifact(sample_dir / f"{sample}.longphase_s.recalibrated.pass.vcf.gz")
    sidecar = contract.artifact(sample_dir / f"{sample}.read_tags.tsv.gz")
    sidecar_index = contract.artifact(sample_dir / f"{sample}.read_tags.tsv.gz.tbi")
    validation_artifact = contract.artifact(sample_dir / "sidecar_validation.json")
    stream_capture_artifact = contract.artifact(sample_dir / "stream_capture_summary.json")
    normalization_artifact = contract.artifact(sample_dir / "normalization_audit.json")
    transition_artifact = contract.artifact(sample_dir / "filter_transition_audit.json")
    sample_verification_artifact = contract.artifact(sample_dir / "sample_verification.json")
    producer_receipt_artifact = contract.artifact(sample_dir / "producer_capture_receipt_v2.json")
    validation = contract.strict_json_load(Path(validation_artifact["path"]))
    contract._validate_native_sidecar_document(validation, sample)
    producer_receipt = contract.strict_json_load(Path(producer_receipt_artifact["path"]))
    expected_receipt = {
        "tumor_bam": {"path": str(bam), "index_path": str(bam_index), "storage_identity_v1": storage},
        "caller_raw_vcf": raw,
        "longphase_input_vcf": normalized,
        "caller_pass_baseline_vcf": caller_pass_baseline,
        "capture_outputs": {
            "sidecar": sidecar,
            "sidecar_index": sidecar_index,
            "native_validation": validation_artifact,
            "stream_capture_summary": stream_capture_artifact,
            "normalization_audit": normalization_artifact,
            "filter_transition_audit": transition_artifact,
            "sample_verification": sample_verification_artifact,
            "longphase_recalibrated_all_vcf": recal_all,
            "longphase_recalibrated_pass_vcf": recal_pass,
        },
        "producer_receipt_sha256": producer_receipt_artifact["identity"]["sha256"],
    }
    evidence = contract._validate_raw_all_producer_receipt(
        producer_receipt,
        sample,
        Path(__file__).resolve().parent.parent / "schemas" / "longphase_raw_all_capture_receipt_v2.schema.json",
        expected_receipt,
    )
    binding = {
        "schema_name": contract.SIDECAR_SUBJECT_SCHEMA,
        "schema_version": contract.SIDECAR_SUBJECT_VERSION,
        "sample": sample,
        "duplicate_identity_policy": contract.REDUNDANT_IDENTITY_POLICY,
        "coordinate_identity_columns": contract.COORDINATE_IDENTITY_COLUMNS,
        "sidecar_sha256": sidecar["identity"]["sha256"],
        "sidecar_index_sha256": sidecar_index["identity"]["sha256"],
        "validation_sha256": validation_artifact["identity"]["sha256"],
        "producer_capture_receipt_sha256": producer_receipt_artifact["identity"]["sha256"],
        "alignment_payload_storage_identity_sha256": storage["identity_sha256"],
        "producer_command_argv_sha256": evidence["command_argv_sha256"],
        "producer_input_binding_sha256": evidence["input_binding_sha256"],
        "producer_effective_options_sha256": evidence["effective_options_sha256"],
        "caller_raw_vcf_sha256": raw["identity"]["sha256"],
        "longphase_input_vcf_sha256": normalized["identity"]["sha256"],
        "caller_pass_baseline_vcf_sha256": caller_pass_baseline["identity"]["sha256"],
        "longphase_recalibrated_all_vcf_sha256": recal_all["identity"]["sha256"],
        "longphase_recalibrated_pass_vcf_sha256": recal_pass["identity"]["sha256"],
        "mapped_alignment_count": evidence["mapped_alignment_count"],
        "identity_unique_count": evidence["identity_unique_count"],
        "duplicate_count": evidence["duplicate_count"],
        "conflict_count": evidence["conflict_count"],
    }
    if tree_input_contract == contract.TREE_INPUT_CONTRACT:
        backbone_role = "longphase_s_recalibrated_filter_pass"
        selected_tree = recal_pass
    elif tree_input_contract == contract.SENSITIVITY_TREE_INPUT_CONTRACT:
        backbone_role = "clairs_filter_pass_sensitivity"
        selected_tree = caller_pass_baseline
    else:
        raise contract.ContractError(
            "E_TREE_VCF_ROLE", f"unsupported selected tree contract: {tree_input_contract}",
            exit_code=4, stage="producer_scope", sample=sample,
        )
    item = {
        "sample": sample,
        "biological_id": base["biological_id"],
        "platform": source.get("platform") or ("ONT_DORADO" if sample.endswith("_DORADO") else "ONT"),
        "replicate_role": "platform_replica" if sample == "HCC1395_DORADO" else "canonical",
        "alignment_payload": {
            "path": str(bam),
            "index_path": str(bam_index),
            "embedded_hp_ps_policy": "ignore",
            "identity_policy": "storage_identity_v1",
            "storage_identity_v1": storage,
        },
        "somatic": {
            "backbone_role": backbone_role,
            "caller_raw_vcf": raw,
            "longphase_input_vcf": normalized,
            "caller_pass_baseline_vcf": caller_pass_baseline,
            "longphase_recalibrated_all_vcf": recal_all,
            "longphase_recalibrated_pass_vcf": recal_pass,
            "tree_vcf": selected_tree,
        },
        "read_tags": {
            "authority": "external_coordinate_sidecar",
            "identity_schema": "coordinate_join_v1",
            "assurance": "bounded_coordinate_equivalence_not_global_payload_identity",
            "tagging_semantics": contract.RAW_ALL_TAGGING_SEMANTICS,
            "duplicate_identity_policy": contract.REDUNDANT_IDENTITY_POLICY,
            "require_unique_identity": False,
            "fallback_policy": "forbidden",
            "sidecar": sidecar,
            "index": sidecar_index,
            "validation": validation_artifact,
            "producer_capture_receipt": producer_receipt_artifact,
            "producer_policy": contract.PRODUCTION_POLICY,
            "subject_binding": binding,
        },
        "copy_number": _copy_number(base, sample),
    }
    summary = {
        "sample": sample,
        "pass": True,
        "filter_policy": "production_default_filter",
        "truth_flags_present": False,
        "identity_schema": "coordinate_join_v1",
        "longphase_input_contract": contract.RAW_ALL_INPUT_CONTRACT,
        "tree_backbone_role": "longphase_s_recalibrated_filter_pass",
        "validation_sha256": validation_artifact["identity"]["sha256"],
    }
    return item, summary


def build_manifest(
    base_manifest: Path,
    longphase_manifest: Path,
    production_root: Path,
    real_data_receipt: Path,
    synthetic_receipt: Path,
    manifest_id: str,
    *,
    tree_input_contract: str = contract.TREE_INPUT_CONTRACT,
    created_at_utc: str | None = None,
) -> dict[str, Any]:
    base = contract.strict_json_load(base_manifest)
    longphase = contract.strict_json_load(longphase_manifest)
    base_by_sample = _require_exact_source_set(base, "base_manifest")
    source_by_sample = _require_exact_source_set(longphase, "longphase_manifest")
    if base.get("dataset_count") != 7 or base.get("biological_sample_count") != 6:
        raise contract.ContractError("E_COUNT_MISMATCH", "base manifest does not declare 7 datasets/6 biological IDs")
    join_method = {
        "assurance": "bounded_coordinate_equivalence_not_global_payload_identity",
        "claim_limit": "join_method_only_not_per_sample_global_payload_identity",
        "real_data_bounded_receipt": contract.artifact(real_data_receipt.resolve()),
        "synthetic_three_case_receipt": contract.artifact(synthetic_receipt.resolve()),
    }
    samples: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for sample in contract.EXPECTED_DATASETS:
        item, summary = _build_sample(
            base_by_sample[sample], source_by_sample[sample], production_root, tree_input_contract
        )
        samples.append(item)
        summaries.append(summary)
    manifest: dict[str, Any] = {
        "$schema": MANIFEST_SCHEMA_REF,
        "schema_name": contract.SCHEMA_NAME,
        "schema_version": contract.SCHEMA_VERSION,
        "manifest_id": manifest_id,
        "created_at_utc": created_at_utc or contract.utc_now(),
        "dataset_count": 7,
        "biological_sample_count": 6,
        "analysis_contract": {
            "task_type": (
                "comprehensive_validation"
                if tree_input_contract == contract.TREE_INPUT_CONTRACT
                else "backbone_sensitivity"
            ),
            "scope": {"name": "whole_autosomes_chr1_22", "contigs": contract.AUTOSOMES},
            "read_tag_mode": "external_sidecar",
            "embedded_tag_policy": "ignore",
            "require_exact_join": True,
            "sidecar_identity_schema": "coordinate_join_v1",
            "sidecar_assurance": "bounded_coordinate_equivalence_not_global_payload_identity",
            "bam_identity_policy": "storage_identity_v1",
            "bam_identity_assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
            "production_filter_policy": contract.PRODUCTION_POLICY,
            "longphase_input_contract": contract.RAW_ALL_INPUT_CONTRACT,
            "tree_input_contract": tree_input_contract,
            "tagging_semantics": contract.RAW_ALL_TAGGING_SEMANTICS,
            "duplicate_identity_policy": contract.REDUNDANT_IDENTITY_POLICY,
            "join_method_validation": join_method,
        },
        "production_summary": {
            "expected_dataset_count": 7,
            "completed_dataset_count": 7,
            "passed_dataset_count": 7,
            "all_pass": True,
            "datasets": summaries,
        },
        "samples": samples,
    }
    schema_path = Path(__file__).resolve().parent.parent / "schemas" / "layered_input_manifest_v3.schema.json"
    contract._pre_schema_checks(manifest)
    contract.apply_json_schema(manifest, schema_path)
    contract.validate_semantics(manifest)
    contract.verify_join_method(manifest["analysis_contract"])
    for sample in manifest["samples"]:
        contract.verify_sample(sample)
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-manifest", required=True, type=Path)
    parser.add_argument("--longphase-manifest", required=True, type=Path)
    parser.add_argument("--production-root", required=True, type=Path)
    parser.add_argument("--real-data-receipt", required=True, type=Path)
    parser.add_argument("--synthetic-receipt", required=True, type=Path)
    parser.add_argument("--manifest-id", required=True)
    parser.add_argument(
        "--tree-input-contract",
        choices=(contract.TREE_INPUT_CONTRACT, contract.SENSITIVITY_TREE_INPUT_CONTRACT),
        default=contract.TREE_INPUT_CONTRACT,
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--failure-report", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.output.exists() or args.failure_report.exists():
        print("E_OUTPUT_EXISTS: manifest/failure targets must both be new", file=sys.stderr)
        return 7
    try:
        manifest = build_manifest(
            args.base_manifest,
            args.longphase_manifest,
            args.production_root,
            args.real_data_receipt,
            args.synthetic_receipt,
            args.manifest_id,
            tree_input_contract=args.tree_input_contract,
        )
        contract.atomic_write_json(args.output, manifest, mode=0o444)
        print(f"LAYERED V3 MANIFEST: PASS 7/7 -> {args.output}")
        return 0
    except contract.ContractError as error:
        failure = contract.failure_document(error, args.base_manifest, args.output)
        failure["schema_name"] = "intersubmod.layered_manifest_build_failure"
        contract.atomic_write_json(args.failure_report, failure)
        print(f"{error.code}: {error} -> {args.failure_report}", file=sys.stderr)
        return error.exit_code
    except Exception as error:  # pragma: no cover
        internal = contract.ContractError("E_INTERNAL", repr(error), exit_code=70, stage="internal")
        failure = contract.failure_document(internal, args.base_manifest, args.output)
        failure["schema_name"] = "intersubmod.layered_manifest_build_failure"
        contract.atomic_write_json(args.failure_report, failure)
        print(f"{internal.code}: {internal} -> {args.failure_report}", file=sys.stderr)
        return internal.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
