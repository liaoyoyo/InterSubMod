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
    / "verify_retrospective_running_source_identity_v2.py"
)
SPEC = importlib.util.spec_from_file_location("verify_source_identity_v2", SCRIPT)
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
    identities = {role: MODULE.BASE.identity(path) for role, path in sources.items()}
    base = datetime.now(timezone.utc) + timedelta(seconds=2)
    command = [str(sources["analyzer"].resolve()), "--fixed"]
    snapshot = {
        "schema_name": MODULE.BASE.SNAPSHOT_SCHEMA,
        "schema_version": MODULE.BASE.SCHEMA_VERSION,
        "created_at_utc": (base + timedelta(seconds=1)).isoformat(),
        "task_type": "B_comprehensive_validation",
        "audit_class": MODULE.BASE.AUDIT_CLASS,
        "process": {
            "started_at_utc": (base - timedelta(seconds=1)).isoformat(),
            "cmdline": ["/usr/bin/python3", *command],
        },
        "expected_command_fragment": identities["analyzer"]["path"],
        "source_identity_during_execution": identities,
        "snapshot_creator_source_identity": MODULE.BASE.identity(MODULE.BASE_PATH),
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
        "schema_name": MODULE.BASE.EXPECTED_MANIFEST_SCHEMA,
        "schema_version": MODULE.BASE.EXPECTED_MANIFEST_VERSION,
        "started_at_utc": base.isoformat(),
        "finished_at_utc": (base + timedelta(seconds=2)).isoformat(),
        "command": command,
        "source_code": {
            role: {
                key: identity[key] for key in ("path", "size_bytes", "sha256")
            }
            for role, identity in identities.items()
        },
        "counts": {
            "primary_stable_sites": MODULE.BASE.EXPECTED_PROCESSED_SITES,
            "processed_sites": MODULE.BASE.EXPECTED_PROCESSED_SITES,
        },
        "pass": True,
    }
    snapshot_path = tmp_path / "snapshot.json"
    manifest_path = tmp_path / "run_manifest.json"
    snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return sources, snapshot_path, manifest_path


def verify(tmp_path: Path, snapshot: Path, manifest: Path) -> dict[str, object]:
    return MODULE.verify_snapshot(
        argparse.Namespace(
            snapshot=snapshot,
            run_manifest=manifest,
            output=tmp_path / "receipt.json",
        )
    )


def test_absolute_script_token_is_bound_to_attested_source(tmp_path: Path) -> None:
    _, snapshot, manifest = make_receipts(tmp_path)
    payload = verify(tmp_path, snapshot, manifest)
    assert payload["schema_version"] == "1.2.0"
    assert payload["pass"] is True
    assert payload["command_binding"]["manifest_script_token_mode"] == "absolute"
    assert payload["checks"]["manifest_script_token_matches_attested_source"] is True
    assert payload["post_run_verifier_source_identity"]["sha256"] == MODULE.BASE.sha256(
        SCRIPT
    )


def test_only_exact_repo_relative_token_is_accepted() -> None:
    relative = SCRIPT.resolve().relative_to(MODULE.REPO_ROOT).as_posix()
    assert MODULE.script_token_matches_attested_source(relative, str(SCRIPT)) is True
    assert MODULE.script_token_matches_attested_source(SCRIPT.name, str(SCRIPT)) is False
    assert (
        MODULE.script_token_matches_attested_source(
            f"scripts/{SCRIPT.name}", str(SCRIPT)
        )
        is False
    )


@pytest.mark.parametrize("token", ["other.py", "../producer.py", "./producer.py"])
def test_exact_live_and_manifest_tokens_still_reject_wrong_source_binding(
    tmp_path: Path, token: str
) -> None:
    _, snapshot_path, manifest_path = make_receipts(tmp_path)
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    snapshot["process"]["cmdline"][1] = token
    manifest["command"][0] = token
    snapshot_path.write_text(json.dumps(snapshot), encoding="utf-8")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(RuntimeError, match="not bound to the attested"):
        verify(tmp_path, snapshot_path, manifest_path)


def test_manifest_token_drift_rejected_before_normalization(tmp_path: Path) -> None:
    _, snapshot_path, manifest_path = make_receipts(tmp_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["command"].append("--drift")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(RuntimeError, match="not exactly bound"):
        verify(tmp_path, snapshot_path, manifest_path)


def test_producer_source_drift_rejected(tmp_path: Path) -> None:
    sources, snapshot_path, manifest_path = make_receipts(tmp_path)
    sources["analyzer"].write_text("print('changed')\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="changed after"):
        verify(tmp_path, snapshot_path, manifest_path)
