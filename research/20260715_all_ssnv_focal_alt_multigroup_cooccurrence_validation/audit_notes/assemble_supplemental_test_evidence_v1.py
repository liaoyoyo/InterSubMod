#!/usr/bin/env python3
"""Assemble the one-time signed supplemental report test-evidence manifest."""

from __future__ import annotations

from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import stat
from types import ModuleType
from typing import Any


HERE = Path(__file__).resolve().parent
TOPIC_ROOT = HERE.parent
BUILDER_PATH = TOPIC_ROOT / "scripts" / "build_positional_singleton_report.py"


def load_builder() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "positional_singleton_report_builder", BUILDER_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load positional-singleton report builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_new_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    path.chmod(0o444)
    descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def main() -> None:
    builder = load_builder()
    authority = builder.FINAL_RELEASE.SOURCE_AUTHORITY
    if (
        builder.TEST_EVIDENCE_MANIFEST.exists()
        or builder.TEST_EVIDENCE_SIGNATURE.exists()
    ):
        raise RuntimeError("Test-evidence output already exists")
    if builder.CANONICAL_PYTEST_XML.resolve() != Path(
        os.environ.get(
            "INTERSUBMOD_SUPPLEMENTAL_TEST_XML",
            str(builder.CANONICAL_PYTEST_XML),
        )
    ).resolve():
        raise RuntimeError("Canonical pytest XML environment override differs")
    verification = builder.validate_pytest_xml(builder.CANONICAL_PYTEST_XML)
    artifact_paths = {
        "canonical_pytest_xml": builder.CANONICAL_PYTEST_XML,
        "report_builder": builder.TEST_EVIDENCE_REPORT_BUILDER_PATH,
        "report_builder_tests": builder.TEST_EVIDENCE_TEST_SOURCE_PATH,
    }
    artifacts = {
        name: authority.artifact(path, include_mode=True)
        for name, path in artifact_paths.items()
    }
    if any(record["mode"] != "0o444" for record in artifacts.values()):
        raise RuntimeError("Canonical test evidence inputs are not mode 0444")
    public_key = authority.artifact(
        builder.TEST_EVIDENCE_PUBLIC_KEY, include_mode=True
    )
    if (
        public_key["sha256"] != builder.TEST_EVIDENCE_PUBLIC_KEY_SHA256
        or public_key["mode"] != "0o444"
    ):
        raise RuntimeError("Supplemental test-evidence public key drift")
    if (
        stat.S_IMODE(builder.TEST_EVIDENCE_PRIVATE_KEY.stat().st_mode)
        != 0o400
    ):
        raise RuntimeError("One-time private key is not prepared at mode 0400")
    payload = {
        "schema_name": "intersubmod.supplemental_report_test_evidence",
        "schema_version": "1.1.0",
        "authority_id": builder.TEST_EVIDENCE_AUTHORITY_ID,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_comprehensive_validation",
        "scope": "positional_singleton_supplemental_report_release",
        "artifacts": artifacts,
        "canonical_test_summary": verification,
        "required_test_names": sorted(builder.REQUIRED_SUPPLEMENTAL_TESTS),
        "signature_contract": {
            "algorithm": "Ed25519",
            "public_key": public_key,
            "signed_artifact": str(
                builder.TEST_EVIDENCE_MANIFEST.resolve()
            ),
            "signature": str(builder.TEST_EVIDENCE_SIGNATURE.resolve()),
            "private_key_lifecycle": (
                "encrypted_one_time_key_chmod_000_after_signing"
            ),
            "private_key_path": str(
                builder.TEST_EVIDENCE_PRIVATE_KEY.resolve()
            ),
        },
        "signature_state_at_manifest_creation": {
            "state": "UNSIGNED_PENDING_OUT_OF_BAND",
            "private_key_mode": "0o400",
            "required_private_key_mode_after_signing": "0o0",
        },
        "checks": {
            "canonical_pytest_pass": True,
            "report_builder_identity_bound": True,
            "report_builder_tests_identity_bound": True,
        },
        "pass": True,
    }
    write_new_json(builder.TEST_EVIDENCE_MANIFEST, payload)
    print(
        json.dumps(
            {
                "manifest": str(builder.TEST_EVIDENCE_MANIFEST),
                "canonical_pytest_xml": artifacts["canonical_pytest_xml"],
                "report_builder": artifacts["report_builder"],
                "report_builder_tests": artifacts["report_builder_tests"],
                "public_key": public_key,
                "canonical_test_summary": verification,
                "pass": True,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
