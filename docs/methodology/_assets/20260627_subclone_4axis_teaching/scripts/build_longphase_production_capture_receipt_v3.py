#!/usr/bin/env python3
"""Normalize immutable LongPhase production evidence into a strict v3 capture receipt."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import platform
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

import validate_layered_v3_inputs as contract


RECEIPT_SCHEMA = Path(__file__).resolve().parent.parent / "schemas" / "longphase_production_capture_receipt_v1.schema.json"
NATIVE_CHECKS = (
    "truth_flags_absent",
    "parser_count_matches_input",
    "capture_pass",
    "execution_alignment_count_matches_capture",
    "sidecar_row_count_matches_capture",
    "tagged_count_matches_execution",
    "sidecar_no_malformed_rows",
    "sidecar_no_unknown_HP",
    "sidecar_no_exact_identity_conflicts",
    "recalibrated_preserves_all_input_keys",
)


def _required_file(path: Path, label: str) -> Path:
    contract.ensure_regular_file(path, label)
    return path


def _indexed(path: Path) -> dict[str, Any]:
    return contract.indexed_artifact(_required_file(path, "indexed artifact"))


def _option_value(argv: list[str], options: tuple[str, ...]) -> str:
    values: list[str] = []
    for index, token in enumerate(argv):
        for option in options:
            if token == option:
                if index + 1 >= len(argv):
                    raise contract.ContractError(
                        "E_PRODUCER_OPTION_MISMATCH", f"missing value after {option}",
                        exit_code=4, stage="producer_scope",
                    )
                values.append(argv[index + 1])
            elif token.startswith(f"{option}="):
                values.append(token.split("=", 1)[1])
    if len(values) != 1:
        raise contract.ContractError(
            "E_PRODUCER_OPTION_MISMATCH", f"requires exactly one of {options}; observed {values}",
            exit_code=4, stage="producer_scope",
        )
    return values[0]


def parse_command(path: Path, binary: Path, sample_dir: Path, sample_meta: dict[str, Any]) -> list[str]:
    lines = [line for line in path.read_text(encoding="utf-8", errors="strict").splitlines() if line.strip()]
    if len(lines) != 1:
        raise contract.ContractError(
            "E_PRODUCER_AMBIGUOUS", "command.sh.txt must contain exactly one non-empty command",
            exit_code=4, stage="producer_scope", artifact=str(path),
        )
    try:
        argv = shlex.split(lines[0], posix=True)
    except ValueError as error:
        raise contract.ContractError(
            "E_PRODUCER_AMBIGUOUS", f"cannot parse command argv: {error}",
            exit_code=4, stage="producer_scope", artifact=str(path),
        ) from error
    if not argv or any(token in {"|", "||", "&&", ";", ">", ">>", "<"} for token in argv):
        raise contract.ContractError(
            "E_PRODUCER_AMBIGUOUS", "command contains shell control operators",
            exit_code=4, stage="producer_scope", artifact=str(path),
        )
    if Path(argv[0]).resolve(strict=True) != binary.resolve(strict=True) or argv[1:2] != ["somatic_haplotag"]:
        raise contract.ContractError(
            "E_PRODUCER_AMBIGUOUS", "command binary/subcommand differs from reviewed LongPhase somatic_haplotag",
            exit_code=4, stage="producer_scope",
        )
    expected_options = {
        ("-s",): sample_meta["germline_phased_vcf"],
        ("-b",): sample_meta["normal_bam"],
        ("--tumor-snv-file",): sample_meta["caller_pass_vcf"],
        ("--tumor-bam-file",): sample_meta["tumor_bam"],
        ("-r",): sample_meta["reference"],
        ("-t", "--threads"): "12",
        ("-q", "--mapq"): "20",
        ("-o",): str(sample_dir / f"{sample_meta['sample']}_production"),
    }
    for options, expected in expected_options.items():
        observed = _option_value(argv, options)
        if observed != str(expected):
            raise contract.ContractError(
                "E_INPUT_BINDING_MISMATCH", f"argv {options} differs from frozen manifest",
                exit_code=4, stage="producer_scope", expected=str(expected), observed=observed,
            )
    if not contract._argv_contains(argv, "--tagSupplementary") or not contract._argv_contains(argv, "--output-somatic-vcf"):
        raise contract.ContractError(
            "E_PRODUCER_OPTION_MISMATCH", "required explicit output/tag options are absent",
            exit_code=4, stage="producer_scope",
        )
    if any(contract._argv_contains(argv, option) for option in ("--truth-vcf", "--truth-bed")):
        raise contract.ContractError(
            "E_TRUTH_FLAG_PRESENT", "truth option found in command argv",
            exit_code=4, stage="producer_scope",
        )
    if contract._argv_contains(argv, "--disableFilter"):
        raise contract.ContractError(
            "E_PRODUCTION_FILTER_POLICY", "--disableFilter found in command argv",
            exit_code=4, stage="producer_scope",
        )
    return argv


def _parse_hash_file(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        match = re.fullmatch(r"([0-9a-f]{64})\s+(.+)", line)
        if not match:
            raise contract.ContractError("E_HASH_MANIFEST", f"malformed hash row {line_number}: {path}", exit_code=3)
        digest, value = match.groups()
        if value in rows:
            raise contract.ContractError("E_HASH_MANIFEST", f"duplicate hash path: {value}", exit_code=3)
        rows[value] = digest
    return rows


def _verify_hash_manifest(path: Path, expected_paths: set[Path] | None = None) -> dict[str, str]:
    rows = _parse_hash_file(path)
    if expected_paths is not None and {Path(value) for value in rows} != expected_paths:
        raise contract.ContractError(
            "E_HASH_MANIFEST", "hash manifest path set differs from expected",
            exit_code=3, stage="artifact_identity",
            expected=sorted(str(value) for value in expected_paths), observed=sorted(rows),
        )
    for value, recorded in rows.items():
        artifact_path = _required_file(Path(value), "hash-manifest artifact")
        observed = contract.file_sha256(artifact_path)
        if observed != recorded:
            raise contract.ContractError(
                "E_HASH_MISMATCH", f"recorded hash changed: {value}",
                exit_code=3, stage="artifact_identity", expected=recorded, observed=observed,
            )
    return rows


def _verify_inventory(path: Path, expected: dict[str, Path]) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = {row["role"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(rows) != set(expected):
        raise contract.ContractError("E_INPUT_INVENTORY", "input inventory roles differ from expected", exit_code=3)
    for role, expected_path in expected.items():
        row = rows[role]
        # The active producer records GNU `stat -c` metadata for the requested
        # path.  GNU stat does not dereference a symlink unless `-L` is used, so
        # compare the logical-path lstat here.  The opened target is bound
        # separately by storage_identity_v1 (realpath, target metadata, sampled
        # chunks and the full index hash).
        logical = expected_path.lstat()
        if Path(row["path"]) != expected_path or int(row["size_bytes"]) != logical.st_size:
            raise contract.ContractError(
                "E_INPUT_BINDING_MISMATCH", f"input inventory mismatch: {role}",
                exit_code=4, stage="producer_scope",
            )
        try:
            date_part, zone = row["mtime"].rsplit(" ", 1)
            whole, fraction = date_part.rsplit(".", 1)
            seconds = int(dt.datetime.strptime(f"{whole} {zone}", "%Y-%m-%d %H:%M:%S %z").timestamp())
            recorded_ns = seconds * 1_000_000_000 + int(fraction[:9].ljust(9, "0"))
        except (ValueError, TypeError) as error:
            raise contract.ContractError("E_INPUT_INVENTORY", f"unparseable mtime: {role}", exit_code=3) from error
        if recorded_ns != logical.st_mtime_ns:
            raise contract.ContractError("E_STORAGE_IDENTITY", f"input inventory mtime changed: {role}", exit_code=3)


def _reference_dict(reference: Path) -> Path:
    candidates = [reference.with_suffix(".dict"), Path(f"{reference}.dict")]
    existing = [path for path in dict.fromkeys(candidates) if path.is_file()]
    if len(existing) != 1:
        raise contract.ContractError(
            "E_REQUIRED_ARTIFACT", f"expected exactly one reference dict; found {[str(x) for x in existing]}",
            exit_code=3, stage="artifact_identity",
        )
    return existing[0]


def _longphase_version(binary: Path) -> str:
    result = subprocess.run([str(binary), "--version"], text=True, capture_output=True, timeout=10, check=False)
    if result.returncode != 0:
        raise contract.ContractError(
            "E_PRODUCER_AMBIGUOUS", f"LongPhase --version exit {result.returncode}",
            exit_code=4, stage="producer_scope",
        )
    match = re.search(r"Version:\s*([^\s]+)", f"{result.stdout}\n{result.stderr}")
    if not match:
        raise contract.ContractError("E_PRODUCER_AMBIGUOUS", "LongPhase version is not parseable", exit_code=4)
    return match.group(1)


def _validate_native(document: dict[str, Any], scan: dict[str, int], sample: str) -> None:
    if document.get("region") != "all" or document.get("pass") is not True:
        raise contract.ContractError("E_SIDECAR_VALIDATION", "native validation is not all-region PASS", exit_code=5)
    failed = [name for name in NATIVE_CHECKS if document.get("checks", {}).get(name) is not True]
    if failed:
        raise contract.ContractError("E_SIDECAR_VALIDATION", f"native validation checks missing/false: {failed}", exit_code=5)
    if (
        document.get("duplicate_exact_alignment_rows") != 0
        or document.get("duplicate_exact_alignment_conflicts") != 0
        or document.get("record_key_missing") != 0
        or document.get("record_key_extra") != 0
    ):
        raise contract.ContractError("E_SIDECAR_VALIDATION", "native validation reports duplicates/conflicts/key loss", exit_code=5)
    capture_rows = document.get("capture", {}).get("rows_mapped")
    if capture_rows != scan["mapped_alignment_count"]:
        raise contract.ContractError(
            "E_SIDECAR_SUBJECT_MISMATCH", f"{sample} native/sidecar row counts differ",
            exit_code=5, expected=scan["mapped_alignment_count"], observed=capture_rows,
        )


def build_receipt(
    sample_dir: Path,
    production_manifest: Path,
    sample: str,
    run_root: Path,
    longphase_binary: Path,
    normalizer_argv: list[str],
) -> dict[str, Any]:
    if sample_dir.resolve() != (run_root / "samples" / sample).resolve():
        raise contract.ContractError(
            "E_PATH_CONTAINMENT", "sample directory is not run_root/samples/<sample>",
            exit_code=7, stage="lifecycle",
        )
    if production_manifest.resolve() != (run_root / "input_manifest.json").resolve():
        raise contract.ContractError(
            "E_INPUT_BINDING_MISMATCH", "production manifest is not the frozen run-root copy",
            exit_code=4, stage="producer_scope",
        )
    manifest = contract.strict_json_load(_required_file(production_manifest, "production manifest"))
    matches = [item for item in manifest.get("samples", []) if item.get("sample") == sample]
    if len(matches) != 1:
        raise contract.ContractError("E_DATASET_SET_MISMATCH", f"sample must occur exactly once: {sample}")
    meta = matches[0]
    if meta.get("biological_id") != contract.EXPECTED_DATASETS.get(sample):
        raise contract.ContractError("E_DATASET_SET_MISMATCH", f"unexpected biological ID for {sample}")
    params_path = _required_file(run_root / "params.json", "run params")
    code_manifest = _required_file(run_root / "code.sha256", "code hash manifest")
    inventory = _required_file(sample_dir / "input_files.tsv", "sample input inventory")
    input_hash_manifest = _required_file(sample_dir / "input.sha256", "sample input hash manifest")
    output_hash_manifest = _required_file(sample_dir / "output.sha256", "sample output hash manifest")
    command_path = _required_file(sample_dir / "command.sh.txt", "producer command")
    native_path = _required_file(sample_dir / "sidecar_validation.json", "native validation")
    sidecar = _required_file(sample_dir / f"{sample}.read_tags.tsv.gz", "read-tag sidecar")
    sidecar_index = _required_file(Path(f"{sidecar}.tbi"), "read-tag sidecar index")
    recal_all = _required_file(sample_dir / f"{sample}.longphase_s.recalibrated.all.vcf.gz", "recalibrated all VCF")
    recal_pass = _required_file(sample_dir / f"{sample}.longphase_s.recalibrated.pass.vcf.gz", "recalibrated PASS VCF")
    binary = _required_file(longphase_binary, "LongPhase binary")
    params = contract.strict_json_load(params_path)
    expected_params = {"threads": 12, "truth_flags": False, "mapq": 20, "tag_supplementary": True,
                       "output_mode": "read_tag_sidecar"}
    if any(params.get(key) != value for key, value in expected_params.items()):
        raise contract.ContractError("E_PRODUCTION_FILTER_POLICY", "run params differ from production profile", exit_code=4)
    if Path(params.get("run_root", "")).resolve() != run_root.resolve():
        raise contract.ContractError(
            "E_INPUT_BINDING_MISMATCH", "params.run_root differs from receipt run root",
            exit_code=4, stage="producer_scope",
        )
    argv = parse_command(command_path, binary, sample_dir, meta)
    paths = {key: Path(meta[key]) for key in (
        "germline_phased_vcf", "normal_bam", "tumor_bam", "caller_pass_vcf", "caller_raw_vcf", "reference"
    )}
    inventory_expected = {
        "germline": paths["germline_phased_vcf"], "normal_bam": paths["normal_bam"],
        "tumor_bam": paths["tumor_bam"], "reference": paths["reference"], "caller_pass": paths["caller_pass_vcf"],
    }
    for value in paths.values():
        _required_file(value, "manifest input")
    _verify_inventory(inventory, inventory_expected)
    germline_index = contract.find_index(paths["germline_phased_vcf"])
    caller_pass_index = contract.find_index(paths["caller_pass_vcf"])
    normal_index = _required_file(Path(f"{paths['normal_bam']}.bai"), "normal BAM index")
    tumor_index = _required_file(Path(f"{paths['tumor_bam']}.bai"), "tumor BAM index")
    expected_input_hashes = {
        paths["caller_pass_vcf"], caller_pass_index, paths["germline_phased_vcf"], germline_index,
        normal_index, tumor_index,
    }
    _verify_hash_manifest(input_hash_manifest, expected_input_hashes)
    plain_sidecar = _required_file(sample_dir / f"{sample}.read_tags.tsv", "plain read-tag sidecar")
    expected_output_hashes = {plain_sidecar, sidecar, recal_all, recal_pass, native_path}
    _verify_hash_manifest(output_hash_manifest, expected_output_hashes)
    code_rows = _verify_hash_manifest(code_manifest)
    code_by_basename = {Path(path).name: Path(path) for path in code_rows}
    required_code = {
        "runner": "run_longphase_production_sidecars.sh",
        "capture": "capture_longphase_tagged_bam_sidecar.py",
        "validator": "validate_streamed_longphase_sidecar.py",
    }
    if any(name not in code_by_basename for name in required_code.values()):
        raise contract.ContractError("E_SOURCE_BUNDLE_MISMATCH", "code.sha256 lacks required producer sources", exit_code=6)
    native = contract.strict_json_load(native_path)
    scan = contract.inspect_coordinate_sidecar(sidecar, sample)
    _validate_native(native, scan, sample)
    if scan["duplicate_count"] != 0 or scan["conflict_count"] != 0 or scan["identity_unique_count"] != scan["mapped_alignment_count"]:
        raise contract.ContractError("E_JOIN_DUPLICATE", "global sidecar scan is not one-to-one", exit_code=5)
    producer_inputs = {
        "germline_phased_vcf": _indexed(paths["germline_phased_vcf"]),
        "normal_bam": {
            "path": str(paths["normal_bam"]), "index_path": str(normal_index),
            "storage_identity_v1": contract.storage_identity_v1(paths["normal_bam"], normal_index),
        },
        "tumor_bam": {
            "path": str(paths["tumor_bam"]), "index_path": str(tumor_index),
            "storage_identity_v1": contract.storage_identity_v1(paths["tumor_bam"], tumor_index),
        },
        "caller_pass_vcf": _indexed(paths["caller_pass_vcf"]),
        "reference": {
            "path": str(paths["reference"]), "fai_path": str(Path(f"{paths['reference']}.fai")),
            "storage_identity_v1": contract.storage_identity_v1(paths["reference"], Path(f"{paths['reference']}.fai")),
            "dict": contract.artifact(_reference_dict(paths["reference"])),
        },
    }
    receipt = {
        "schema_name": contract.PRODUCER_RECEIPT_SCHEMA,
        "schema_version": contract.PRODUCER_RECEIPT_VERSION,
        "sample": sample,
        "evidence_origin": "post_run_normalization_from_frozen_execution_artifacts",
        "identity_schema": "coordinate_join_v1",
        "assurance": "bounded_coordinate_equivalence_not_global_payload_identity",
        "production_policy": contract.PRODUCTION_POLICY,
        "effective_options": contract.PRODUCTION_EFFECTIVE_OPTIONS,
        "command_argv": argv,
        "command_argv_sha256": contract.canonical_sha256(argv),
        "producer_inputs": producer_inputs,
        "producer_input_binding_sha256": contract.canonical_sha256(producer_inputs),
        "layered_ledger_artifacts": {"caller_raw_vcf": _indexed(paths["caller_raw_vcf"])},
        "longphase_executable": {**contract.artifact(binary), "version": _longphase_version(binary)},
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
            "longphase_recalibrated_all_vcf": _indexed(recal_all),
            "longphase_recalibrated_pass_vcf": _indexed(recal_pass),
        },
        "global_coordinate_counts": scan,
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
            args.sample_dir.resolve(), args.production_manifest.resolve(), args.sample,
            args.run_root.resolve(), args.longphase_binary.resolve(), effective_argv,
        )
        contract.atomic_write_json(args.output, receipt, mode=0o444)
        readback = contract.strict_json_load(args.output)
        contract.apply_json_schema(readback, RECEIPT_SCHEMA)
        if contract.canonical_sha256(readback) != contract.canonical_sha256(receipt):
            raise contract.ContractError("E_RECEIPT_READBACK", "receipt readback changed", exit_code=7, stage="lifecycle")
        print(f"PRODUCER CAPTURE RECEIPT: PASS {args.sample} -> {args.output}")
        return 0
    except contract.ContractError as error:
        failure = contract.failure_document(error, args.production_manifest, args.output)
        failure["schema_name"] = "intersubmod.longphase_production_capture_receipt_failure"
        contract.atomic_write_json(args.failure_report, failure)
        print(f"{error.code}: {error} -> {args.failure_report}", file=sys.stderr)
        return error.exit_code
    except Exception as error:  # pragma: no cover
        internal = contract.ContractError("E_INTERNAL", repr(error), exit_code=70, stage="internal")
        failure = contract.failure_document(internal, args.production_manifest, args.output)
        failure["schema_name"] = "intersubmod.longphase_production_capture_receipt_failure"
        contract.atomic_write_json(args.failure_report, failure)
        print(f"{internal.code}: {internal} -> {args.failure_report}", file=sys.stderr)
        return internal.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
