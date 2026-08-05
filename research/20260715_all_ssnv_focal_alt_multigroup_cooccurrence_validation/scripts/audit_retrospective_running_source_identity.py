#!/usr/bin/env python3
"""Hash-bind source files observed during a running producer to its final manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat as stat_module
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence


SNAPSHOT_SCHEMA = "intersubmod.retrospective_running_source_identity.snapshot"
RECEIPT_SCHEMA = "intersubmod.retrospective_running_source_identity.receipt"
SCHEMA_VERSION = "1.1.0"
EXPECTED_MANIFEST_SCHEMA = "intersubmod.all_ssnv_tumor_ref_controls.run_manifest"
EXPECTED_MANIFEST_VERSION = "1.0.0"
EXPECTED_SOURCE_ROLES = frozenset({"analyzer", "focal_alt_cluster_lib"})
EXPECTED_SNAPSHOT_CHECKS = frozenset(
    {
        "process_alive_at_snapshot",
        "source_role_set_exact",
        "analyzer_role_equals_live_process_script",
        "interpreter_token_matches_process_executable",
        "all_sources_predate_process",
    }
)
EXPECTED_PROCESSED_SITES = 102_842
AUDIT_CLASS = "bounded_retrospective_source_file_identity"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def parse_utc(value: str) -> datetime:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise ValueError(f"Timestamp lacks timezone: {value}")
    return parsed.astimezone(timezone.utc)


def _stat_identity(stat: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        stat.st_dev,
        stat.st_ino,
        stat.st_size,
        stat.st_mode,
        stat.st_mtime_ns,
        stat.st_ctime_ns,
    )


def read_stable_bytes(path: Path) -> tuple[Path, bytes, os.stat_result]:
    resolved = path.resolve(strict=True)
    with resolved.open("rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat_module.S_ISREG(before.st_mode):
            raise FileNotFoundError(f"Not a regular file: {resolved}")
        content = handle.read()
        after = os.fstat(handle.fileno())
    if _stat_identity(before) != _stat_identity(after) or len(content) != before.st_size:
        raise RuntimeError(f"File identity changed while reading: {resolved}")
    return resolved, content, before


def sha256(path: Path) -> str:
    _, content, _ = read_stable_bytes(path)
    return hashlib.sha256(content).hexdigest()


def identity(path: Path) -> dict[str, Any]:
    resolved, content, stat = read_stable_bytes(path)
    return {
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "sha256": hashlib.sha256(content).hexdigest(),
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "mode": oct(stat.st_mode & 0o777),
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
    }


def artifact(path: Path) -> dict[str, Any]:
    resolved, content, stat = read_stable_bytes(path)
    return {
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def read_json_artifact(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    resolved, content, stat = read_stable_bytes(path)
    payload = json.loads(content.decode("utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"JSON artifact is not an object: {resolved}")
    return payload, {
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def process_info(pid: int) -> dict[str, Any]:
    proc = Path("/proc") / str(pid)
    cmdline = [
        token.decode("utf-8", errors="surrogateescape")
        for token in (proc / "cmdline").read_bytes().split(b"\0")
        if token
    ]
    stat_text = (proc / "stat").read_text(encoding="utf-8")
    right_paren = stat_text.rfind(")")
    if right_paren < 0:
        raise RuntimeError(f"Malformed /proc/{pid}/stat")
    fields_from_three = stat_text[right_paren + 2 :].split()
    parent_pid = int(fields_from_three[1])
    start_ticks = int(fields_from_three[19])
    clock_ticks = os.sysconf("SC_CLK_TCK")
    boot_epoch = time.time() - float(Path("/proc/uptime").read_text().split()[0])
    started_epoch = boot_epoch + (start_ticks / clock_ticks)
    return {
        "pid": pid,
        "parent_pid": parent_pid,
        "executable": str((proc / "exe").resolve(strict=True)),
        "cmdline": cmdline,
        "parent_cmdline": [
            token.decode("utf-8", errors="surrogateescape")
            for token in (Path("/proc") / str(parent_pid) / "cmdline")
            .read_bytes()
            .split(b"\0")
            if token
        ],
        "linux_start_ticks": start_ticks,
        "clock_ticks_per_second": clock_ticks,
        "started_at_utc": datetime.fromtimestamp(
            started_epoch, tz=timezone.utc
        ).isoformat(timespec="microseconds"),
    }


def parse_sources(values: Sequence[str]) -> dict[str, Path]:
    sources: dict[str, Path] = {}
    for value in values:
        role, separator, raw_path = value.partition("=")
        if not separator or not role or not raw_path:
            raise ValueError(f"--source must use ROLE=PATH: {value}")
        if role in sources:
            raise ValueError(f"Duplicate source role: {role}")
        sources[role] = Path(raw_path)
    if set(sources) != EXPECTED_SOURCE_ROLES:
        raise ValueError(
            "Source roles must be exact: "
            f"expected={sorted(EXPECTED_SOURCE_ROLES)} observed={sorted(sources)}"
        )
    return sources


def write_exclusive(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    path.chmod(0o444)


def is_contiguous_subsequence(container: Sequence[str], expected: Sequence[str]) -> bool:
    if not expected or len(expected) > len(container):
        return False
    return any(
        list(container[index : index + len(expected)]) == list(expected)
        for index in range(len(container) - len(expected) + 1)
    )


def create_snapshot(args: argparse.Namespace) -> dict[str, Any]:
    sources = parse_sources(args.source)
    analyzer = sources["analyzer"].resolve(strict=True)
    expected_fragment = Path(args.expected_command_fragment).resolve(strict=True)
    if expected_fragment != analyzer:
        raise RuntimeError("Expected command fragment must equal the analyzer role path")
    process = process_info(args.pid)
    if len(process["cmdline"]) < 2 or Path(process["cmdline"][1]).resolve(strict=True) != analyzer:
        raise RuntimeError(
            "The live process script token does not equal the analyzer role path"
        )
    if Path(process["cmdline"][0]).resolve(strict=True) != Path(
        process["executable"]
    ).resolve(strict=True):
        raise RuntimeError("The live process interpreter token/executable identity differs")
    source_identities = {role: identity(path) for role, path in sources.items()}
    creator_identity = identity(Path(__file__))
    process_start_ns = int(parse_utc(process["started_at_utc"]).timestamp() * 1e9)
    sources_predate_process = all(
        item["mtime_ns"] <= process_start_ns and item["ctime_ns"] <= process_start_ns
        for item in source_identities.values()
    )
    if not sources_predate_process:
        raise RuntimeError("At least one source identity does not predate the producer process")
    payload = {
        "schema_name": SNAPSHOT_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "task_type": "B_comprehensive_validation",
        "audit_class": AUDIT_CLASS,
        "process": process,
        "expected_command_fragment": str(expected_fragment),
        "source_identity_during_execution": source_identities,
        "snapshot_creator_source_identity": creator_identity,
        "checks": {
            "process_alive_at_snapshot": True,
            "source_role_set_exact": True,
            "analyzer_role_equals_live_process_script": True,
            "interpreter_token_matches_process_executable": True,
            "all_sources_predate_process": True,
        },
        "limitation": (
            "This snapshot attests named source-file identity during execution; it is not "
            "a prelaunch lock or a complete environment attestation."
        ),
        "pass": True,
    }
    write_exclusive(args.output, payload)
    return payload


def verify_snapshot(args: argparse.Namespace) -> dict[str, Any]:
    snapshot, snapshot_artifact = read_json_artifact(args.snapshot)
    manifest, manifest_artifact = read_json_artifact(args.run_manifest)
    if (
        snapshot.get("schema_name") != SNAPSHOT_SCHEMA
        or snapshot.get("schema_version") != SCHEMA_VERSION
        or snapshot.get("audit_class") != AUDIT_CLASS
        or snapshot.get("task_type") != "B_comprehensive_validation"
        or snapshot.get("pass") is not True
        or not isinstance(snapshot.get("checks"), dict)
        or set(snapshot["checks"]) != EXPECTED_SNAPSHOT_CHECKS
        or not all(value is True for value in snapshot["checks"].values())
    ):
        raise RuntimeError("Source identity snapshot is not a passing supported receipt")
    if (
        manifest.get("schema_name") != EXPECTED_MANIFEST_SCHEMA
        or manifest.get("schema_version") != EXPECTED_MANIFEST_VERSION
    ):
        raise RuntimeError("Tumor-REF run-manifest schema mismatch")
    if manifest.get("pass") is not True:
        raise RuntimeError("Tumor-REF run manifest is not passing")
    counts = manifest.get("counts") or {}
    if counts.get("primary_stable_sites") != EXPECTED_PROCESSED_SITES or counts.get(
        "processed_sites"
    ) != EXPECTED_PROCESSED_SITES:
        raise RuntimeError("Tumor-REF run manifest does not cover all 102842 M1 sites")

    during = snapshot.get("source_identity_during_execution")
    if not isinstance(during, dict) or set(during) != EXPECTED_SOURCE_ROLES:
        raise RuntimeError("Snapshot source role set is not exact")
    creator_during = snapshot.get("snapshot_creator_source_identity")
    creator_after = identity(Path(__file__))
    if creator_during != creator_after:
        raise RuntimeError("Snapshot creator/verifier source identity changed")
    after = {
        role: identity(Path(item["path"]))
        for role, item in during.items()
    }
    if after != during:
        raise RuntimeError("A source file changed after the during-execution snapshot")

    manifest_sources = manifest.get("source_code")
    if not isinstance(manifest_sources, dict) or set(manifest_sources) != EXPECTED_SOURCE_ROLES:
        raise RuntimeError("Tumor-REF manifest source role set is not exact")
    for role, source in during.items():
        expected = {
            key: source[key] for key in ("path", "size_bytes", "sha256")
        }
        if manifest_sources.get(role) != expected:
            raise RuntimeError(f"Manifest source identity mismatch for role: {role}")

    started = parse_utc(manifest["started_at_utc"])
    finished = parse_utc(manifest["finished_at_utc"])
    observed = parse_utc(snapshot["created_at_utc"])
    process_started = parse_utc(snapshot["process"]["started_at_utc"])
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
        or not isinstance(live_cmdline, list)
        or live_cmdline[1:] != manifest_command
        or manifest_command[0] != during["analyzer"]["path"]
        or snapshot.get("expected_command_fragment") != during["analyzer"]["path"]
    ):
        raise RuntimeError("Live process command and manifest command are not exactly bound")

    payload = {
        "schema_name": RECEIPT_SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": now_utc(),
        "task_type": "B_comprehensive_validation",
        "audit_class": AUDIT_CLASS,
        "snapshot": snapshot_artifact,
        "tumor_ref_run_manifest": manifest_artifact,
        "source_identity_during_execution": during,
        "source_identity_after_execution": after,
        "snapshot_creator_source_identity": creator_during,
        "auditor_source_identity_after_execution": creator_after,
        "checks": {
            "snapshot_pass": True,
            "producer_manifest_pass": True,
            "producer_processed_all_102842_sites": True,
            "source_role_sets_exact": True,
            "snapshot_creator_source_unchanged": True,
            "process_start_within_manifest_start_clock_tolerance": True,
            "process_start_clock_tolerance_seconds": process_start_clock_tolerance_seconds,
            "snapshot_observed_between_manifest_start_and_finish": True,
            "all_sources_predate_manifest_start": True,
            "all_source_identities_unchanged_after_snapshot": True,
            "manifest_source_artifacts_match_snapshot": True,
            "manifest_command_exactly_matches_live_process_command": True,
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
    write_exclusive(args.output, payload)
    return payload


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)

    snapshot = subparsers.add_parser("snapshot")
    snapshot.add_argument("--pid", type=int, required=True)
    snapshot.add_argument("--source", action="append", default=[])
    snapshot.add_argument("--expected-command-fragment", required=True)
    snapshot.add_argument("--output", type=Path, required=True)

    verify = subparsers.add_parser("verify")
    verify.add_argument("--snapshot", type=Path, required=True)
    verify.add_argument("--run-manifest", type=Path, required=True)
    verify.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    payload = create_snapshot(args) if args.mode == "snapshot" else verify_snapshot(args)
    print(
        json.dumps(
            {
                "schema_name": payload["schema_name"],
                "output": str(args.output.resolve()),
                "pass": payload["pass"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
