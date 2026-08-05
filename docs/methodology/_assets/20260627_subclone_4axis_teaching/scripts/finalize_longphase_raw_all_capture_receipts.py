#!/usr/bin/env python3
"""Create restart-safe raw-all capture receipts after a completed production run."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import build_longphase_raw_all_capture_receipt_v2 as builder
import validate_layered_v3_inputs as contract


EXPECTED_SAMPLES = tuple(contract.EXPECTED_DATASETS)


def _builder_argv(run_root: Path, sample: str, binary: Path, output: Path, failure: Path) -> list[str]:
    return [
        str(Path(builder.__file__).resolve()),
        "--sample-dir", str(run_root / "samples" / sample),
        "--production-manifest", str(run_root / "input_manifest.json"),
        "--sample", sample,
        "--run-root", str(run_root),
        "--longphase-binary", str(binary),
        "--output", str(output),
        "--failure-report", str(failure),
    ]


def _write_or_verify(path: Path, value: dict[str, Any]) -> None:
    if path.exists():
        observed = contract.strict_json_load(path)
        if contract.canonical_sha256(observed) != contract.canonical_sha256(value):
            raise contract.ContractError(
                "E_RECEIPT_DRIFT", f"existing receipt differs from fresh normalization: {path}",
                exit_code=7, stage="lifecycle",
            )
        return
    contract.atomic_write_json(path, value, mode=0o444)


def verify_existing_closeout(run_root: Path, closeout: Path, marker: Path) -> dict[str, Any]:
    closeout_document = contract.strict_json_load(closeout)
    marker_document = contract.strict_json_load(marker)
    if (
        marker_document.get("status") != "SUCCESS"
        or marker_document.get("dataset_count") != 7
        or Path(marker_document.get("closeout", "")).resolve() != closeout.resolve()
        or marker_document.get("closeout_sha256") != contract.file_sha256(closeout)
    ):
        raise contract.ContractError(
            "E_RECEIPT_CLOSEOUT", "existing raw-all receipt marker does not bind its closeout",
            exit_code=7, stage="lifecycle",
        )
    rows = closeout_document.get("receipts")
    if (
        closeout_document.get("schema_name") != "intersubmod.longphase_raw_all_receipt_closeout"
        or closeout_document.get("schema_version") != "1.0.0"
        or closeout_document.get("run_root") != str(run_root)
        or closeout_document.get("dataset_count") != 7
        or closeout_document.get("all_pass") is not True
        or not isinstance(rows, list)
        or tuple(row.get("sample") for row in rows) != EXPECTED_SAMPLES
    ):
        raise contract.ContractError(
            "E_RECEIPT_CLOSEOUT", "existing raw-all receipt closeout is not exact 7/7 PASS",
            exit_code=7, stage="lifecycle",
        )
    contract.verify_artifact(closeout_document["builder"], "ALL", "receipt closeout builder")
    contract.verify_artifact(closeout_document["receipt_schema"], "ALL", "receipt closeout schema")
    for row in rows:
        sample = row["sample"]
        verified = contract.verify_artifact(row["receipt"], sample, "raw-all producer receipt")
        receipt = contract.strict_json_load(Path(verified["path"]))
        contract.apply_json_schema(receipt, builder.RECEIPT_SCHEMA)
        if (
            receipt.get("sample") != sample
            or receipt.get("bam_output_policy", {}).get("persisted_bam") is not False
            or receipt.get("global_coordinate_counts", {}).get("mapped_alignment_count")
            != row.get("mapped_alignment_count")
            or receipt.get("filter_transition_summary", {}).get("input_record_count")
            != row.get("input_record_count")
            or row.get("persisted_bam") is not False
        ):
            raise contract.ContractError(
                "E_RECEIPT_CLOSEOUT", f"existing receipt row does not bind {sample}",
                exit_code=7, stage="lifecycle", sample=sample,
            )
    return closeout_document


def finalize(run_root: Path) -> dict[str, Any]:
    root_success = contract.strict_json_load(run_root / "_SUCCESS")
    manifest = contract.strict_json_load(run_root / "input_manifest.json")
    if root_success.get("status") != "SUCCESS":
        raise contract.ContractError("E_PRODUCTION_INCOMPLETE", "production root is not SUCCESS", exit_code=4)
    if (
        manifest.get("schema_version") != "2.0"
        or manifest.get("dataset_count") != 7
        or tuple(item.get("sample") for item in manifest.get("samples", [])) != EXPECTED_SAMPLES
    ):
        raise contract.ContractError("E_DATASET_SET_MISMATCH", "raw-all manifest is not the canonical 7 datasets")
    binary = run_root / "source_bundle" / "longphase-s"
    contract.ensure_regular_file(binary, "frozen LongPhase-S")
    rows = []
    for sample in EXPECTED_SAMPLES:
        sample_dir = run_root / "samples" / sample
        output = sample_dir / "producer_capture_receipt_v2.json"
        failure = sample_dir / "producer_capture_receipt_v2.failure.json"
        if failure.exists():
            raise contract.ContractError(
                "E_PRIOR_RECEIPT_FAILURE", f"prior receipt failure requires audit: {failure}",
                exit_code=7, stage="lifecycle", sample=sample,
            )
        argv = _builder_argv(run_root, sample, binary, output, failure)
        receipt = builder.build_receipt(
            sample_dir,
            run_root / "input_manifest.json",
            sample,
            run_root,
            binary,
            argv,
        )
        _write_or_verify(output, receipt)
        readback = contract.strict_json_load(output)
        contract.apply_json_schema(readback, builder.RECEIPT_SCHEMA)
        rows.append({
            "sample": sample,
            "receipt": contract.artifact(output),
            "mapped_alignment_count": readback["global_coordinate_counts"]["mapped_alignment_count"],
            "input_record_count": readback["filter_transition_summary"]["input_record_count"],
            "rescued_nonpass_to_pass": readback["filter_transition_summary"]["rescued_nonpass_to_pass"],
            "removed_pass_to_nonpass": readback["filter_transition_summary"]["removed_pass_to_nonpass"],
            "persisted_bam": readback["bam_output_policy"]["persisted_bam"],
        })
        print(f"RAW-ALL RECEIPT {sample}: PASS -> {output}", flush=True)
    return {
        "schema_name": "intersubmod.longphase_raw_all_receipt_closeout",
        "schema_version": "1.0.0",
        "run_root": str(run_root),
        "dataset_count": len(rows),
        "all_pass": len(rows) == 7 and all(item["persisted_bam"] is False for item in rows),
        "total_mapped_alignments": sum(item["mapped_alignment_count"] for item in rows),
        "total_input_records": sum(item["input_record_count"] for item in rows),
        "total_rescued_nonpass_to_pass": sum(item["rescued_nonpass_to_pass"] for item in rows),
        "total_removed_pass_to_nonpass": sum(item["removed_pass_to_nonpass"] for item in rows),
        "builder": contract.artifact(Path(builder.__file__).resolve()),
        "receipt_schema": contract.artifact(builder.RECEIPT_SCHEMA.resolve()),
        "receipts": rows,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args(argv)
    run_root = args.run_root.resolve()
    closeout = run_root / "raw_all_receipt_closeout.json"
    marker = run_root / "_RAW_ALL_RECEIPTS_SUCCESS"
    failure = run_root / "raw_all_receipt_closeout.failure.json"
    if failure.exists():
        print(f"E_PRIOR_RECEIPT_FAILURE: {failure}", file=sys.stderr)
        return 7
    try:
        if closeout.exists() or marker.exists():
            if not closeout.is_file() or not marker.is_file():
                raise contract.ContractError(
                    "E_RECEIPT_CLOSEOUT", "existing closeout and success marker must appear together",
                    exit_code=7, stage="lifecycle",
                )
            verify_existing_closeout(run_root, closeout, marker)
            print(f"RAW-ALL RECEIPT CLOSEOUT: EXISTING 7/7 VERIFIED -> {closeout}")
            return 0
        value = finalize(run_root)
        if value["all_pass"] is not True:
            raise contract.ContractError("E_RECEIPT_CLOSEOUT", "receipt closeout is not 7/7 PASS", exit_code=7)
        _write_or_verify(closeout, value)
        marker_value = {
            "status": "SUCCESS",
            "closeout": str(closeout),
            "closeout_sha256": contract.file_sha256(closeout),
            "dataset_count": 7,
        }
        _write_or_verify(marker, marker_value)
        print(f"RAW-ALL RECEIPT CLOSEOUT: 7/7 PASS -> {closeout}")
        return 0
    except contract.ContractError as error:
        document = {
            "schema_name": "intersubmod.longphase_raw_all_receipt_closeout_failure",
            "schema_version": "1.0.0",
            "run_root": str(run_root),
            "error": error.record(),
        }
        if not failure.exists():
            contract.atomic_write_json(failure, document)
        print(f"{error.code}: {error} -> {failure}", file=sys.stderr)
        return error.exit_code
    except Exception as error:  # pragma: no cover
        internal = contract.ContractError("E_INTERNAL", repr(error), exit_code=70, stage="internal")
        document = {
            "schema_name": "intersubmod.longphase_raw_all_receipt_closeout_failure",
            "schema_version": "1.0.0",
            "run_root": str(run_root),
            "error": internal.record(),
        }
        if not failure.exists():
            contract.atomic_write_json(failure, document)
        print(f"{internal.code}: {internal} -> {failure}", file=sys.stderr)
        return internal.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
