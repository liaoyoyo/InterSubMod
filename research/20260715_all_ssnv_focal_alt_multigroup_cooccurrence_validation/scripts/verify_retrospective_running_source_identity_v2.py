#!/usr/bin/env python3
"""Verify a v1.1 running-source snapshot with explicit command-token normalization."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


BASE_PATH = Path(__file__).with_name("audit_retrospective_running_source_identity.py")
BASE_SPEC = importlib.util.spec_from_file_location(
    "retrospective_running_source_identity_v1", BASE_PATH
)
if BASE_SPEC is None or BASE_SPEC.loader is None:
    raise RuntimeError(f"Unable to load snapshot creator module: {BASE_PATH}")
BASE = importlib.util.module_from_spec(BASE_SPEC)
BASE_SPEC.loader.exec_module(BASE)

RECEIPT_SCHEMA_VERSION = "1.2.0"
REPO_ROOT = Path(__file__).resolve().parents[3]


def script_token_matches_attested_source(token: str, source_path: str) -> bool:
    """Bind an absolute or exact repository-relative script token."""
    if not token or "\x00" in token:
        return False
    token_path = Path(token)
    source = Path(source_path)
    if token_path.is_absolute():
        try:
            return token_path.resolve(strict=True) == source.resolve(strict=True)
        except FileNotFoundError:
            return False
    raw_parts = token.split("/")
    if any(part in {"", ".", ".."} for part in raw_parts):
        return False
    if not token_path.parts:
        return False
    try:
        expected_relative = source.resolve(strict=True).relative_to(REPO_ROOT)
    except (FileNotFoundError, ValueError):
        return False
    return token_path.as_posix() == expected_relative.as_posix()


def _source_artifact_from_identity(source: Mapping[str, Any]) -> dict[str, Any]:
    return {key: source[key] for key in ("path", "size_bytes", "sha256")}


def verify_snapshot(args: argparse.Namespace) -> dict[str, Any]:
    snapshot, snapshot_artifact = BASE.read_json_artifact(args.snapshot)
    manifest, manifest_artifact = BASE.read_json_artifact(args.run_manifest)
    if (
        snapshot.get("schema_name") != BASE.SNAPSHOT_SCHEMA
        or snapshot.get("schema_version") != BASE.SCHEMA_VERSION
        or snapshot.get("audit_class") != BASE.AUDIT_CLASS
        or snapshot.get("task_type") != "B_comprehensive_validation"
        or snapshot.get("pass") is not True
        or not isinstance(snapshot.get("checks"), dict)
        or set(snapshot["checks"]) != BASE.EXPECTED_SNAPSHOT_CHECKS
        or not all(value is True for value in snapshot["checks"].values())
    ):
        raise RuntimeError("Source identity snapshot is not a passing supported receipt")
    if (
        manifest.get("schema_name") != BASE.EXPECTED_MANIFEST_SCHEMA
        or manifest.get("schema_version") != BASE.EXPECTED_MANIFEST_VERSION
    ):
        raise RuntimeError("Tumor-REF run-manifest schema mismatch")
    if manifest.get("pass") is not True:
        raise RuntimeError("Tumor-REF run manifest is not passing")
    counts = manifest.get("counts") or {}
    if counts.get("primary_stable_sites") != BASE.EXPECTED_PROCESSED_SITES or counts.get(
        "processed_sites"
    ) != BASE.EXPECTED_PROCESSED_SITES:
        raise RuntimeError("Tumor-REF run manifest does not cover all 102842 M1 sites")

    during = snapshot.get("source_identity_during_execution")
    if not isinstance(during, dict) or set(during) != BASE.EXPECTED_SOURCE_ROLES:
        raise RuntimeError("Snapshot source role set is not exact")
    creator_during = snapshot.get("snapshot_creator_source_identity")
    if not isinstance(creator_during, dict):
        raise RuntimeError("Snapshot creator source identity is missing")
    creator_path = Path(str(creator_during.get("path", "")))
    creator_after = BASE.identity(creator_path)
    if creator_during != creator_after:
        raise RuntimeError("Snapshot creator source identity changed after snapshot")
    post_run_verifier = BASE.identity(Path(__file__))

    after = {role: BASE.identity(Path(item["path"])) for role, item in during.items()}
    if after != during:
        raise RuntimeError("A source file changed after the during-execution snapshot")
    manifest_sources = manifest.get("source_code")
    if not isinstance(manifest_sources, dict) or set(
        manifest_sources
    ) != BASE.EXPECTED_SOURCE_ROLES:
        raise RuntimeError("Tumor-REF manifest source role set is not exact")
    for role, source in during.items():
        if manifest_sources.get(role) != _source_artifact_from_identity(source):
            raise RuntimeError(f"Manifest source identity mismatch for role: {role}")

    started = BASE.parse_utc(manifest["started_at_utc"])
    finished = BASE.parse_utc(manifest["finished_at_utc"])
    observed = BASE.parse_utc(snapshot["created_at_utc"])
    process_started = BASE.parse_utc(snapshot["process"]["started_at_utc"])
    process_start_clock_tolerance_seconds = 2.0
    if not (
        abs(process_started.timestamp() - started.timestamp())
        <= process_start_clock_tolerance_seconds
        and started <= observed <= finished
    ):
        raise RuntimeError("Process, manifest, snapshot, and finish timestamps are inconsistent")
    started_ns = int(started.timestamp() * 1e9)
    if not all(
        source["mtime_ns"] <= started_ns and source["ctime_ns"] <= started_ns
        for source in during.values()
    ):
        raise RuntimeError("A source identity does not predate manifest start")

    manifest_command = manifest.get("command")
    live_cmdline = snapshot.get("process", {}).get("cmdline")
    if (
        not isinstance(manifest_command, list)
        or not manifest_command
        or not all(isinstance(token, str) for token in manifest_command)
        or not isinstance(live_cmdline, list)
        or len(live_cmdline) < 2
        or live_cmdline[1:] != manifest_command
    ):
        raise RuntimeError("Live process command and manifest command are not exactly bound")
    if not script_token_matches_attested_source(
        manifest_command[0], during["analyzer"]["path"]
    ):
        raise RuntimeError("Manifest script token is not bound to the attested analyzer source")
    if snapshot.get("expected_command_fragment") != during["analyzer"]["path"]:
        raise RuntimeError("Snapshot expected analyzer path is not bound to source identity")

    payload = {
        "schema_name": BASE.RECEIPT_SCHEMA,
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "created_at_utc": BASE.now_utc(),
        "task_type": "B_comprehensive_validation",
        "audit_class": BASE.AUDIT_CLASS,
        "snapshot": snapshot_artifact,
        "tumor_ref_run_manifest": manifest_artifact,
        "source_identity_during_execution": during,
        "source_identity_after_execution": after,
        "snapshot_creator_source_identity": creator_during,
        "snapshot_creator_source_identity_after_execution": creator_after,
        "post_run_verifier_source_identity": post_run_verifier,
        "auditor_source_identity_after_execution": post_run_verifier,
        "command_binding": {
            "live_python_launcher_token": live_cmdline[0],
            "manifest_script_token": manifest_command[0],
            "manifest_script_token_mode": (
                "absolute"
                if Path(manifest_command[0]).is_absolute()
                else "repo_relative_exact"
            ),
            "attested_analyzer_path": during["analyzer"]["path"],
            "live_after_launcher_exactly_equals_manifest": True,
            "relative_token_rejects_dot_and_parent_segments": True,
            "repo_relative_token_must_equal_attested_source_relative_to_repo_root": True,
        },
        "checks": {
            "snapshot_pass": True,
            "producer_manifest_pass": True,
            "producer_processed_all_102842_sites": True,
            "source_role_sets_exact": True,
            "snapshot_creator_source_unchanged": True,
            "post_run_verifier_source_identity_recorded": True,
            "process_start_within_manifest_start_clock_tolerance": True,
            "process_start_clock_tolerance_seconds": process_start_clock_tolerance_seconds,
            "snapshot_observed_between_manifest_start_and_finish": True,
            "all_sources_predate_manifest_start": True,
            "all_source_identities_unchanged_after_snapshot": True,
            "manifest_source_artifacts_match_snapshot": True,
            "manifest_command_exactly_matches_live_process_command": True,
            "manifest_script_token_matches_attested_source": True,
            "snapshot_and_manifest_parsed_from_hash_bound_bytes": True,
            "manifest_hash_bound_to_receipt": True,
        },
        "pass_semantics": (
            "Named producer source files were observed during execution, predated the run, "
            "remained identity-equal afterward, and match the passing producer manifest."
        ),
        "limitation": (
            "Retrospective bounded source-file attestation only; not a prelaunch lock or a "
            "complete interpreter, package, kernel, hardware, or environment attestation."
        ),
        "pass": True,
    }
    BASE.write_exclusive(args.output, payload)
    return payload


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)
    verify = subparsers.add_parser("verify")
    verify.add_argument("--snapshot", type=Path, required=True)
    verify.add_argument("--run-manifest", type=Path, required=True)
    verify.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    payload = verify_snapshot(args)
    print(
        json.dumps(
            {
                "schema_name": payload["schema_name"],
                "schema_version": payload["schema_version"],
                "output": str(args.output.resolve()),
                "pass": payload["pass"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
