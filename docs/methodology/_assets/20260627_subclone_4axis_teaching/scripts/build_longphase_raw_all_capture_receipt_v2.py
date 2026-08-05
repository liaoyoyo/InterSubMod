#!/usr/bin/env python3
"""Normalize a completed raw-all LongPhase-S sample into a strict receipt v2."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import platform
import shlex
import stat
import sys
from pathlib import Path
from typing import Any

import audit_longphase_filter_transitions as filter_audit
import build_longphase_production_capture_receipt_v3 as legacy
import validate_layered_v3_inputs as contract


_SCRIPT_DIR = Path(__file__).resolve().parent
_SCHEMA_CANDIDATES = (
    _SCRIPT_DIR.parent / "schemas" / "longphase_raw_all_capture_receipt_v2.schema.json",
    _SCRIPT_DIR / "longphase_raw_all_capture_receipt_v2.schema.json",
)
RECEIPT_SCHEMA = next((path for path in _SCHEMA_CANDIDATES if path.is_file()), _SCHEMA_CANDIDATES[0])
RECEIPT_SCHEMA_NAME = "intersubmod.longphase_raw_all_capture_receipt"
RECEIPT_SCHEMA_VERSION = "2.0.0"
NATIVE_CHECKS = legacy.NATIVE_CHECKS


def _verify_inventory(path: Path, expected: dict[str, Path]) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = {row["role"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(rows) != set(expected):
        raise contract.ContractError("E_INPUT_INVENTORY", "raw-all input inventory roles differ", exit_code=3)
    for role, expected_path in expected.items():
        row = rows[role]
        observed = expected_path.stat()
        if Path(row["path"]) != expected_path or int(row["size_bytes"]) != observed.st_size:
            raise contract.ContractError(
                "E_INPUT_BINDING_MISMATCH",
                f"raw-all input inventory mismatch: {role}",
                exit_code=4,
                stage="producer_scope",
            )
        try:
            date_part, zone = row["mtime"].rsplit(" ", 1)
            whole, fraction = date_part.rsplit(".", 1)
            seconds = int(dt.datetime.strptime(f"{whole} {zone}", "%Y-%m-%d %H:%M:%S %z").timestamp())
            recorded_ns = seconds * 1_000_000_000 + int(fraction[:9].ljust(9, "0"))
        except (TypeError, ValueError) as error:
            raise contract.ContractError(
                "E_INPUT_INVENTORY", f"unparseable raw-all inventory mtime: {role}", exit_code=3
            ) from error
        if recorded_ns != observed.st_mtime_ns:
            raise contract.ContractError(
                "E_STORAGE_IDENTITY", f"raw-all input inventory mtime changed: {role}", exit_code=3
            )


def _parse_command(
    path: Path,
    binary: Path,
    sample_dir: Path,
    sample_meta: dict[str, Any],
    threads: int,
) -> list[str]:
    lines = [line for line in path.read_text(encoding="utf-8", errors="strict").splitlines() if line.strip()]
    if len(lines) != 1:
        raise contract.ContractError(
            "E_PRODUCER_AMBIGUOUS",
            "command.sh.txt must contain exactly one non-empty command",
            exit_code=4,
            stage="producer_scope",
        )
    try:
        argv = shlex.split(lines[0], posix=True)
    except ValueError as error:
        raise contract.ContractError(
            "E_PRODUCER_AMBIGUOUS", f"cannot parse command argv: {error}", exit_code=4, stage="producer_scope"
        ) from error
    if not argv or any(token in {"|", "||", "&&", ";", ">", ">>", "<"} for token in argv):
        raise contract.ContractError(
            "E_PRODUCER_AMBIGUOUS", "command contains shell control operators", exit_code=4, stage="producer_scope"
        )
    if Path(argv[0]).resolve(strict=True) != binary.resolve(strict=True) or argv[1:2] != ["somatic_haplotag"]:
        raise contract.ContractError(
            "E_PRODUCER_AMBIGUOUS",
            "command binary/subcommand differs from frozen LongPhase-S",
            exit_code=4,
            stage="producer_scope",
        )
    normalized = sample_dir / f"{sample_meta['sample']}.clairs.raw_all.normalized.vcf.gz"
    expected_options = {
        ("-s",): sample_meta["germline_phased_vcf"],
        ("-b",): sample_meta["normal_bam"],
        ("--tumor-snv-file",): str(normalized),
        ("--tumor-bam-file",): sample_meta["tumor_bam"],
        ("-r",): sample_meta["reference"],
        ("-t", "--threads"): str(threads),
        ("-q", "--mapq"): "20",
        ("-o",): str(sample_dir / f"{sample_meta['sample']}_production"),
    }
    for options, expected in expected_options.items():
        observed = legacy._option_value(argv, options)
        if observed != str(expected):
            raise contract.ContractError(
                "E_INPUT_BINDING_MISMATCH",
                f"argv {options} differs from frozen raw-all contract",
                exit_code=4,
                stage="producer_scope",
                expected=str(expected),
                observed=observed,
            )
    if not contract._argv_contains(argv, "--tagSupplementary") or not contract._argv_contains(
        argv, "--output-somatic-vcf"
    ):
        raise contract.ContractError(
            "E_PRODUCER_OPTION_MISMATCH", "required output/tag options are absent", exit_code=4, stage="producer_scope"
        )
    if any(contract._argv_contains(argv, option) for option in ("--truth-vcf", "--truth-bed")):
        raise contract.ContractError(
            "E_TRUTH_FLAG_PRESENT", "truth option found in raw-all command", exit_code=4, stage="producer_scope"
        )
    if contract._argv_contains(argv, "--disableFilter"):
        raise contract.ContractError(
            "E_PRODUCTION_FILTER_POLICY", "--disableFilter found in raw-all command", exit_code=4, stage="producer_scope"
        )
    return argv


def _require_json(path: Path, label: str) -> dict[str, Any]:
    legacy._required_file(path, label)
    return contract.strict_json_load(path)


def _audit_without_paths(value: dict[str, Any]) -> dict[str, Any]:
    projected = dict(value)
    for role in ("input", "output"):
        artifact = dict(projected.get(role, {}))
        artifact.pop("path", None)
        projected[role] = artifact
    return projected


def _verify_audits(
    raw_vcf: Path,
    normalized_vcf: Path,
    recalibrated_vcf: Path,
    normalization_path: Path,
    transition_path: Path,
    sample_verification_path: Path,
    expected_records: int,
    sample: str,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    normalization = _require_json(normalization_path, "normalization audit")
    transitions = _require_json(transition_path, "FILTER transition audit")
    sample_verification = _require_json(sample_verification_path, "sample verification")
    fresh_normalization = filter_audit.audit(raw_vcf, normalized_vcf)
    fresh_transitions = filter_audit.audit(normalized_vcf, recalibrated_vcf)
    if normalization != fresh_normalization:
        raise contract.ContractError(
            "E_NORMALIZATION_AUDIT", f"{sample} stored normalization audit differs from fresh audit", exit_code=4
        )
    native_recalibrated_vcf = transition_path.parent / f"{sample}_production_sc.vcf"
    legacy._required_file(native_recalibrated_vcf, "native recalibrated VCF")
    stored_paths = {
        "input": transitions.get("input", {}).get("path"),
        "output": transitions.get("output", {}).get("path"),
    }
    fresh_paths = {
        "input": fresh_transitions.get("input", {}).get("path"),
        "output": fresh_transitions.get("output", {}).get("path"),
    }
    if stored_paths != {"input": str(normalized_vcf), "output": str(native_recalibrated_vcf)}:
        raise contract.ContractError(
            "E_FILTER_TRANSITION_AUDIT",
            f"{sample} stored FILTER audit paths do not bind the native producer artifacts",
            exit_code=4,
        )
    if fresh_paths != {"input": str(normalized_vcf), "output": str(recalibrated_vcf)}:
        raise contract.ContractError(
            "E_FILTER_TRANSITION_AUDIT",
            f"{sample} fresh FILTER audit paths do not bind the packaged artifacts",
            exit_code=4,
        )
    if _audit_without_paths(transitions) != _audit_without_paths(fresh_transitions):
        raise contract.ContractError(
            "E_FILTER_TRANSITION_AUDIT",
            f"{sample} native and packaged FILTER audits differ semantically",
            exit_code=4,
        )
    if (
        normalization.get("pass") is not True
        or normalization.get("input", {}).get("record_count") != expected_records
        or normalization.get("output", {}).get("record_count") != expected_records
        or normalization.get("rescued_nonpass_to_pass") != 0
        or normalization.get("removed_pass_to_nonpass") != 0
    ):
        raise contract.ContractError(
            "E_NORMALIZATION_AUDIT", f"{sample} normalization changed the raw-all record universe", exit_code=4
        )
    if (
        transitions.get("pass") is not True
        or transitions.get("input", {}).get("record_count") != expected_records
        or transitions.get("output", {}).get("record_count") != expected_records
    ):
        raise contract.ContractError(
            "E_FILTER_TRANSITION_AUDIT", f"{sample} LongPhase-S did not conserve every raw-all record", exit_code=4
        )
    if (
        sample_verification.get("sample") != sample
        or sample_verification.get("expected_raw_records") != expected_records
        or sample_verification.get("pass") is not True
        or sample_verification.get("normalization") != normalization
        or sample_verification.get("filter_transitions") != transitions
    ):
        raise contract.ContractError(
            "E_SAMPLE_VERIFICATION", f"{sample} sample verification does not bind both audits", exit_code=4
        )
    return normalization, transitions, sample_verification


def _validate_native_redundancy(
    document: dict[str, Any], scan: dict[str, int], sample: str
) -> None:
    if document.get("region") != "all" or document.get("pass") is not True:
        raise contract.ContractError(
            "E_SIDECAR_VALIDATION", "native validation is not all-region PASS", exit_code=5
        )
    failed = [name for name in NATIVE_CHECKS if document.get("checks", {}).get(name) is not True]
    if failed:
        raise contract.ContractError(
            "E_SIDECAR_VALIDATION", f"native validation checks missing/false: {failed}", exit_code=5
        )
    if (
        document.get("duplicate_exact_alignment_rows") != scan["duplicate_count"]
        or document.get("duplicate_exact_alignment_conflicts") != scan["conflict_count"]
        or document.get("record_key_missing") != 0
        or document.get("record_key_extra") != 0
        or document.get("capture", {}).get("rows_mapped") != scan["mapped_alignment_count"]
    ):
        raise contract.ContractError(
            "E_SIDECAR_VALIDATION",
            "native validation does not match redundant/conflicting coordinate counts",
            exit_code=5,
        )
    contract.validate_redundant_coordinate_counts(scan, sample)


def build_receipt(
    sample_dir: Path,
    production_manifest: Path,
    sample: str,
    run_root: Path,
    longphase_binary: Path,
    normalizer_argv: list[str],
) -> dict[str, Any]:
    if sample_dir.resolve() != (run_root / "samples" / sample).resolve():
        raise contract.ContractError("E_PATH_CONTAINMENT", "sample directory is outside run root", exit_code=7)
    if production_manifest.resolve() != (run_root / "input_manifest.json").resolve():
        raise contract.ContractError(
            "E_INPUT_BINDING_MISMATCH", "production manifest is not the frozen run-root copy", exit_code=4
        )
    root_success = _require_json(run_root / "_SUCCESS", "production root success")
    sample_success = _require_json(sample_dir / "_SUCCESS", "sample success")
    if root_success.get("status") != "SUCCESS" or sample_success.get("status") != "SUCCESS":
        raise contract.ContractError("E_PRODUCTION_INCOMPLETE", "production root/sample is not SUCCESS", exit_code=4)
    manifest = _require_json(production_manifest, "production manifest")
    if manifest.get("schema_version") != "2.0":
        raise contract.ContractError("E_SCHEMA_VERSION", "raw-all production manifest must be 2.0")
    matches = [item for item in manifest.get("samples", []) if item.get("sample") == sample]
    if len(matches) != 1:
        raise contract.ContractError("E_DATASET_SET_MISMATCH", f"sample must occur exactly once: {sample}")
    meta = matches[0]
    if (
        meta.get("biological_id") != contract.EXPECTED_DATASETS.get(sample)
        or meta.get("longphase_input_contract") != "normalized_ClairS_raw_all"
    ):
        raise contract.ContractError("E_DATASET_SET_MISMATCH", f"unexpected raw-all binding for {sample}")

    params_path = legacy._required_file(run_root / "params.json", "run params")
    code_manifest = legacy._required_file(run_root / "code.sha256", "code hash manifest")
    inventory = legacy._required_file(sample_dir / "input_files.tsv", "sample input inventory")
    input_hash_manifest = legacy._required_file(sample_dir / "input.sha256", "sample input hash manifest")
    output_hash_manifest = legacy._required_file(sample_dir / "output.sha256", "sample output hash manifest")
    command_path = legacy._required_file(sample_dir / "command.sh.txt", "producer command")
    native_path = legacy._required_file(sample_dir / "sidecar_validation.json", "native validation")
    stream_summary = legacy._required_file(sample_dir / "stream_capture_summary.json", "stream capture summary")
    normalization_path = legacy._required_file(sample_dir / "normalization_audit.json", "normalization audit")
    transition_path = legacy._required_file(sample_dir / "filter_transition_audit.json", "FILTER audit")
    sample_verification_path = legacy._required_file(sample_dir / "sample_verification.json", "sample verification")
    sidecar = legacy._required_file(sample_dir / f"{sample}.read_tags.tsv.gz", "read-tag sidecar")
    sidecar_index = legacy._required_file(Path(f"{sidecar}.tbi"), "read-tag sidecar index")
    plain_sidecar = legacy._required_file(sample_dir / f"{sample}.read_tags.tsv", "plain read-tag sidecar")
    normalized = legacy._required_file(
        sample_dir / f"{sample}.clairs.raw_all.normalized.vcf.gz", "normalized raw-all VCF"
    )
    recal_all = legacy._required_file(
        sample_dir / f"{sample}.longphase_s.recalibrated.all.vcf.gz", "recalibrated all VCF"
    )
    recal_pass = legacy._required_file(
        sample_dir / f"{sample}.longphase_s.recalibrated.pass.vcf.gz", "recalibrated PASS VCF"
    )
    binary = legacy._required_file(longphase_binary, "frozen LongPhase-S binary")
    params = _require_json(params_path, "run params")
    expected_params = {
        "threads": 12,
        "truth_flags": False,
        "mapq": 20,
        "tag_supplementary": True,
        "output_mode": "read_tag_sidecar",
        "longphase_input_contract": "normalized_ClairS_raw_all",
        "filter_contract": "bidirectional_recalibration",
    }
    if any(params.get(key) != value for key, value in expected_params.items()):
        raise contract.ContractError("E_PRODUCTION_FILTER_POLICY", "run params differ from raw-all profile", exit_code=4)
    if Path(params.get("run_root", "")).resolve() != run_root.resolve():
        raise contract.ContractError("E_INPUT_BINDING_MISMATCH", "params.run_root differs", exit_code=4)
    argv = _parse_command(command_path, binary, sample_dir, meta, int(params["threads"]))

    paths = {
        key: Path(meta[key])
        for key in (
            "germline_phased_vcf",
            "normal_bam",
            "tumor_bam",
            "caller_raw_vcf",
            "caller_pass_vcf",
            "reference",
        )
    }
    paths["longphase_input_vcf"] = normalized
    for value in paths.values():
        legacy._required_file(value, "manifest input")
    inventory_expected = {
        "germline": paths["germline_phased_vcf"],
        "normal_bam": paths["normal_bam"],
        "tumor_bam": paths["tumor_bam"],
        "reference": paths["reference"],
        "caller_raw": paths["caller_raw_vcf"],
        "longphase_input": normalized,
        "caller_pass_baseline": paths["caller_pass_vcf"],
    }
    _verify_inventory(inventory, inventory_expected)
    germline_index = contract.find_index(paths["germline_phased_vcf"])
    raw_index = contract.find_index(paths["caller_raw_vcf"])
    normalized_index = contract.find_index(normalized)
    caller_pass_index = contract.find_index(paths["caller_pass_vcf"])
    normal_index = legacy._required_file(Path(f"{paths['normal_bam']}.bai"), "normal BAM index")
    tumor_index = legacy._required_file(Path(f"{paths['tumor_bam']}.bai"), "tumor BAM index")
    expected_input_hashes = {
        paths["caller_raw_vcf"], raw_index, normalized, normalized_index,
        paths["caller_pass_vcf"], caller_pass_index,
        paths["germline_phased_vcf"], germline_index, normal_index, tumor_index, binary,
    }
    legacy._verify_hash_manifest(input_hash_manifest, expected_input_hashes)
    expected_output_hashes = {
        plain_sidecar,
        sidecar,
        sidecar_index,
        recal_all,
        contract.find_index(recal_all),
        recal_pass,
        contract.find_index(recal_pass),
        normalization_path,
        transition_path,
        native_path,
        stream_summary,
        sample_verification_path,
    }
    legacy._verify_hash_manifest(output_hash_manifest, expected_output_hashes)
    code_rows = legacy._verify_hash_manifest(code_manifest)
    code_by_basename = {Path(path).name: Path(path) for path in code_rows}
    required_code = {
        "runner": "run_longphase_raw_all_production_sidecars.sh",
        "capture": "capture_longphase_tagged_bam_sidecar.py",
        "validator": "validate_streamed_longphase_sidecar.py",
        "filter_auditor": "audit_longphase_filter_transitions.py",
    }
    if any(name not in code_by_basename for name in required_code.values()):
        raise contract.ContractError("E_SOURCE_BUNDLE_MISMATCH", "code hash lacks raw-all sources", exit_code=6)

    expected_records = int(meta["caller_raw_scan"]["records"])
    normalization, transitions, sample_verification = _verify_audits(
        paths["caller_raw_vcf"], normalized, recal_all, normalization_path, transition_path,
        sample_verification_path, expected_records, sample,
    )
    native = _require_json(native_path, "native validation")

    fifo = sample_dir / "consumed_tagged_bam.fifo"
    if not fifo.exists() or not stat.S_ISFIFO(fifo.stat().st_mode):
        raise contract.ContractError("E_BAM_OUTPUT_POLICY", f"consumed FIFO is missing: {fifo}", exit_code=4)
    regular_bams = sorted(path for path in sample_dir.glob("*.bam") if path.is_file())
    if regular_bams:
        raise contract.ContractError(
            "E_BAM_OUTPUT_POLICY", f"persisted BAM payloads are forbidden: {regular_bams}", exit_code=4
        )

    frozen_patch_receipt = run_root / "source_bundle" / Path(manifest["longphase_binary"]["patch_receipt"]).name
    patch_document = _require_json(frozen_patch_receipt, "frozen patch receipt")
    frozen_patch = run_root / "source_bundle" / Path(patch_document["patch"]).name
    legacy._required_file(frozen_patch, "frozen source patch")
    if (
        patch_document.get("status") != "APPROVED_FOR_FAIL_CLOSED_7_DATASET_VALIDATION"
        or patch_document.get("approval", {}).get("scope") != "FAIL_CLOSED_7_DATASET_VALIDATION_ONLY"
        or patch_document.get("patched_binary_sha256") != contract.file_sha256(binary)
    ):
        raise contract.ContractError("E_PATCH_EVIDENCE", "frozen patch approval/binary binding mismatch", exit_code=4)

    producer_inputs = {
        "germline_phased_vcf": legacy._indexed(paths["germline_phased_vcf"]),
        "normal_bam": {
            "path": str(paths["normal_bam"]),
            "index_path": str(normal_index),
            "storage_identity_v1": contract.storage_identity_v1(paths["normal_bam"], normal_index),
        },
        "tumor_bam": {
            "path": str(paths["tumor_bam"]),
            "index_path": str(tumor_index),
            "storage_identity_v1": contract.storage_identity_v1(paths["tumor_bam"], tumor_index),
        },
        "caller_raw_vcf": legacy._indexed(paths["caller_raw_vcf"]),
        "longphase_input_vcf": legacy._indexed(normalized),
        "caller_pass_baseline_vcf": legacy._indexed(paths["caller_pass_vcf"]),
        "reference": {
            "path": str(paths["reference"]),
            "fai_path": str(Path(f"{paths['reference']}.fai")),
            "storage_identity_v1": contract.storage_identity_v1(
                paths["reference"], Path(f"{paths['reference']}.fai")
            ),
        },
    }
    scan = contract.inspect_coordinate_sidecar(sidecar, sample)
    _validate_native_redundancy(native, scan, sample)
    receipt = {
        "schema_name": RECEIPT_SCHEMA_NAME,
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "sample": sample,
        "evidence_origin": "post_run_normalization_from_frozen_raw_all_execution_artifacts",
        "identity_schema": "coordinate_join_v1",
        "assurance": "bounded_coordinate_equivalence_not_global_payload_identity",
        "longphase_input_contract": "normalized_ClairS_raw_all",
        "tree_input_contract": "longphase_s_recalibrated_FILTER_PASS",
        "duplicate_identity_policy": contract.REDUNDANT_IDENTITY_POLICY,
        "production_policy": contract.PRODUCTION_POLICY,
        "effective_options": contract.PRODUCTION_EFFECTIVE_OPTIONS,
        "command_argv": argv,
        "command_argv_sha256": contract.canonical_sha256(argv),
        "producer_inputs": producer_inputs,
        "producer_input_binding_sha256": contract.canonical_sha256(producer_inputs),
        "longphase_executable": {**contract.artifact(binary), "version": legacy._longphase_version(binary)},
        "patch_evidence": {
            "approval_scope": "FAIL_CLOSED_7_DATASET_VALIDATION_ONLY",
            "patch_receipt": contract.artifact(frozen_patch_receipt),
            "source_patch": contract.artifact(frozen_patch),
        },
        "producer_code": {
            role: contract.artifact(code_by_basename[name]) for role, name in required_code.items()
        },
        "environment_lock": {
            "production_manifest": contract.artifact(production_manifest),
            "run_params": contract.artifact(params_path),
            "code_hash_manifest": contract.artifact(code_manifest),
            "sample_input_inventory": contract.artifact(inventory),
            "sample_input_hash_manifest": contract.artifact(input_hash_manifest),
            "sample_output_hash_manifest": contract.artifact(output_hash_manifest),
            "python_executable": contract.artifact(Path(sys.executable).resolve()),
            "python_version": sys.version,
            "platform": platform.platform(),
            "normalizer_source": contract.artifact(Path(__file__).resolve()),
            "normalizer_argv": normalizer_argv,
        },
        "capture_outputs": {
            "sidecar": contract.artifact(sidecar),
            "sidecar_index": contract.artifact(sidecar_index),
            "native_validation": contract.artifact(native_path),
            "stream_capture_summary": contract.artifact(stream_summary),
            "normalization_audit": contract.artifact(normalization_path),
            "filter_transition_audit": contract.artifact(transition_path),
            "sample_verification": contract.artifact(sample_verification_path),
            "longphase_recalibrated_all_vcf": legacy._indexed(recal_all),
            "longphase_recalibrated_pass_vcf": legacy._indexed(recal_pass),
        },
        "global_coordinate_counts": scan,
        "filter_transition_summary": {
            "input_record_count": transitions["input"]["record_count"],
            "output_record_count": transitions["output"]["record_count"],
            "rescued_nonpass_to_pass": transitions["rescued_nonpass_to_pass"],
            "removed_pass_to_nonpass": transitions["removed_pass_to_nonpass"],
            "transition_counts": transitions["filter_transition_counts"],
            "pass": transitions["pass"],
        },
        "bam_output_policy": {
            "transport": "named_fifo",
            "persisted_bam": False,
            "regular_bam_count": 0,
            "consumed_fifo_path": str(fifo),
            "is_fifo_at_closeout": True,
        },
    }
    contract.apply_json_schema(receipt, RECEIPT_SCHEMA)
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--production-manifest", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--longphase-binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--failure-report", required=True, type=Path)
    args = parser.parse_args(argv)
    effective_argv = [str(Path(__file__).resolve()), *(argv if argv is not None else sys.argv[1:])]
    if args.output.exists() or args.failure_report.exists():
        print("E_OUTPUT_EXISTS: receipt/failure targets must both be new", file=sys.stderr)
        return 7
    try:
        receipt = build_receipt(
            args.sample_dir.resolve(),
            args.production_manifest.resolve(),
            args.sample,
            args.run_root.resolve(),
            args.longphase_binary.resolve(),
            effective_argv,
        )
        contract.atomic_write_json(args.output, receipt, mode=0o444)
        readback = contract.strict_json_load(args.output)
        contract.apply_json_schema(readback, RECEIPT_SCHEMA)
        if contract.canonical_sha256(readback) != contract.canonical_sha256(receipt):
            raise contract.ContractError("E_RECEIPT_READBACK", "receipt readback changed", exit_code=7)
        print(f"RAW-ALL CAPTURE RECEIPT: PASS {args.sample} -> {args.output}")
        return 0
    except contract.ContractError as error:
        failure = contract.failure_document(error, args.production_manifest, args.output)
        failure["schema_name"] = "intersubmod.longphase_raw_all_capture_receipt_failure"
        contract.atomic_write_json(args.failure_report, failure)
        print(f"{error.code}: {error} -> {args.failure_report}", file=sys.stderr)
        return error.exit_code
    except Exception as error:  # pragma: no cover
        internal = contract.ContractError("E_INTERNAL", repr(error), exit_code=70, stage="internal")
        failure = contract.failure_document(internal, args.production_manifest, args.output)
        failure["schema_name"] = "intersubmod.longphase_raw_all_capture_receipt_failure"
        contract.atomic_write_json(args.failure_report, failure)
        print(f"{internal.code}: {internal} -> {args.failure_report}", file=sys.stderr)
        return internal.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
