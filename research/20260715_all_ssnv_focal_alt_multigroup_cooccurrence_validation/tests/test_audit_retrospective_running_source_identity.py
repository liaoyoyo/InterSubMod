from __future__ import annotations

import argparse
import importlib.util
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "audit_retrospective_running_source_identity.py"
)
SPEC = importlib.util.spec_from_file_location("audit_retrospective_source_identity", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def make_receipts(tmp_path: Path) -> tuple[dict[str, Path], Path, Path]:
    sources = {
        "analyzer": tmp_path / "producer.py",
        "focal_alt_cluster_lib": tmp_path / "focal_alt_cluster_lib.py",
    }
    sources["analyzer"].write_text("print('fixed')\n", encoding="utf-8")
    sources["focal_alt_cluster_lib"].write_text("MIN_SIZE = 3\n", encoding="utf-8")
    source_identities = {
        role: MODULE.identity(path) for role, path in sources.items()
    }
    base = datetime.now(timezone.utc) + timedelta(seconds=2)
    command = [str(sources["analyzer"].resolve()), "--fixed"]
    snapshot = {
        "schema_name": MODULE.SNAPSHOT_SCHEMA,
        "schema_version": MODULE.SCHEMA_VERSION,
        "created_at_utc": (base + timedelta(seconds=1)).isoformat(),
        "task_type": "B_comprehensive_validation",
        "audit_class": MODULE.AUDIT_CLASS,
        "process": {
            "started_at_utc": (base - timedelta(seconds=1)).isoformat(),
            "cmdline": ["/usr/bin/python3", *command],
        },
        "expected_command_fragment": source_identities["analyzer"]["path"],
        "source_identity_during_execution": source_identities,
        "snapshot_creator_source_identity": MODULE.identity(SCRIPT),
        "checks": {
            "process_alive_at_snapshot": True,
            "source_role_set_exact": True,
            "analyzer_role_equals_live_process_script": True,
            "interpreter_token_matches_process_executable": True,
            "all_sources_predate_process": True,
        },
        "pass": True,
    }
    manifest = {
        "schema_name": MODULE.EXPECTED_MANIFEST_SCHEMA,
        "schema_version": MODULE.EXPECTED_MANIFEST_VERSION,
        "started_at_utc": base.isoformat(),
        "finished_at_utc": (base + timedelta(seconds=2)).isoformat(),
        "command": command,
        "source_code": {
            role: {
                key: source_identity[key]
                for key in ("path", "size_bytes", "sha256")
            }
            for role, source_identity in source_identities.items()
        },
        "counts": {
            "primary_stable_sites": MODULE.EXPECTED_PROCESSED_SITES,
            "processed_sites": MODULE.EXPECTED_PROCESSED_SITES,
        },
        "pass": True,
    }
    snapshot_path = tmp_path / "snapshot.json"
    manifest_path = tmp_path / "run_manifest.json"
    snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return sources, snapshot_path, manifest_path


def test_verify_hash_binds_unchanged_sources_and_manifest(tmp_path: Path) -> None:
    _, snapshot_path, manifest_path = make_receipts(tmp_path)
    output = tmp_path / "receipt.json"
    payload = MODULE.verify_snapshot(
        argparse.Namespace(
            snapshot=snapshot_path,
            run_manifest=manifest_path,
            output=output,
        )
    )
    assert payload["pass"] is True
    assert payload["checks"]["all_source_identities_unchanged_after_snapshot"] is True
    assert payload["checks"]["manifest_source_artifacts_match_snapshot"] is True
    assert payload["tumor_ref_run_manifest"]["sha256"] == MODULE.sha256(manifest_path)
    assert output.is_file()


def test_verify_rejects_source_drift_after_snapshot(tmp_path: Path) -> None:
    sources, snapshot_path, manifest_path = make_receipts(tmp_path)
    sources["analyzer"].write_text("print('changed')\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="changed after"):
        MODULE.verify_snapshot(
            argparse.Namespace(
                snapshot=snapshot_path,
                run_manifest=manifest_path,
                output=tmp_path / "receipt.json",
            )
        )


def test_verify_rejects_missing_source_role(tmp_path: Path) -> None:
    _, snapshot_path, manifest_path = make_receipts(tmp_path)
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    snapshot["source_identity_during_execution"].pop("focal_alt_cluster_lib")
    snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
    with pytest.raises(RuntimeError, match="role set"):
        MODULE.verify_snapshot(
            argparse.Namespace(
                snapshot=snapshot_path,
                run_manifest=manifest_path,
                output=tmp_path / "receipt.json",
            )
        )


def test_verify_rejects_snapshot_creator_identity_drift(tmp_path: Path) -> None:
    _, snapshot_path, manifest_path = make_receipts(tmp_path)
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    snapshot["snapshot_creator_source_identity"]["sha256"] = "0" * 64
    snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
    with pytest.raises(RuntimeError, match="creator/verifier"):
        MODULE.verify_snapshot(
            argparse.Namespace(
                snapshot=snapshot_path,
                run_manifest=manifest_path,
                output=tmp_path / "receipt.json",
            )
        )


@pytest.mark.parametrize(
    ("artifact_name", "mutate", "match"),
    [
        (
            "snapshot",
            lambda payload: payload.update(schema_version="0.0.0"),
            "passing supported receipt",
        ),
        (
            "snapshot",
            lambda payload: payload["checks"].pop("all_sources_predate_process"),
            "passing supported receipt",
        ),
        (
            "manifest",
            lambda payload: payload["counts"].update(processed_sites=1),
            "does not cover all",
        ),
        (
            "manifest",
            lambda payload: payload["command"].append("--drift"),
            "not exactly bound",
        ),
        (
            "manifest",
            lambda payload: payload.update(schema_version="0.0.0"),
            "schema mismatch",
        ),
    ],
)
def test_verify_rejects_contract_drift(
    tmp_path: Path,
    artifact_name: str,
    mutate: object,
    match: str,
) -> None:
    _, snapshot_path, manifest_path = make_receipts(tmp_path)
    path = snapshot_path if artifact_name == "snapshot" else manifest_path
    payload = json.loads(path.read_text(encoding="utf-8"))
    mutate(payload)  # type: ignore[operator]
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(RuntimeError, match=match):
        MODULE.verify_snapshot(
            argparse.Namespace(
                snapshot=snapshot_path,
                run_manifest=manifest_path,
                output=tmp_path / "receipt.json",
            )
        )
