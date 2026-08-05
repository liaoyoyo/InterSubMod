#!/usr/bin/env python3
"""Run a screen producer while proving its source files stayed byte-identical."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


SCHEMA_NAME = "intersubmod.source_locked_screen_run"
SCHEMA_VERSION = "1.1.0"
THREAD_ENV_KEYS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "sha256": sha256(resolved),
        "inode": stat.st_ino,
        "device": stat.st_dev,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
        "mode": oct(stat.st_mode & 0o777),
    }


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--producer", type=Path, required=True)
    parser.add_argument("--locked-file", type=Path, action="append", default=[])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if len(args.command) < 2:
        parser.error("command must contain the Python interpreter and producer path")
    return args


def validate_child_receipt(
    output_dir: Path,
    command: Sequence[str],
    producer_before: dict[str, Any],
) -> tuple[Path, dict[str, Any]]:
    receipt_path = output_dir / "run_manifest.json"
    if not receipt_path.is_file() or receipt_path.stat().st_size == 0:
        raise RuntimeError("Screen producer did not create a non-empty run_manifest.json")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("schema_name") != "intersubmod.all_ssnv_focal_alt_multigroup_run_manifest":
        raise RuntimeError("Screen run-manifest schema mismatch")
    if receipt.get("pass") is not True or receipt.get("status") != "EXECUTION_PASS":
        raise RuntimeError("Screen run-manifest is not passing")
    if receipt.get("command") != list(command[1:]):
        raise RuntimeError("Screen run-manifest command differs from the executed command")
    analyzer = (receipt.get("source_code") or {}).get("analyzer")
    if not isinstance(analyzer, dict):
        raise RuntimeError("Screen run-manifest lacks analyzer identity")
    expected_artifact = {
        key: producer_before[key] for key in ("path", "size_bytes", "sha256")
    }
    if analyzer != expected_artifact:
        raise RuntimeError("Screen run-manifest analyzer identity differs from start lock")
    return receipt_path, receipt


def run(args: argparse.Namespace) -> Path:
    output_dir = args.output_dir.resolve()
    if os.path.lexists(output_dir):
        raise FileExistsError(f"Refusing to reuse output directory: {output_dir}")
    workdir = args.workdir.resolve(strict=True)
    producer = args.producer.resolve(strict=True)
    command = [str(value) for value in args.command]
    if Path(command[1]).resolve(strict=True) != producer:
        raise RuntimeError("Executed producer path differs from --producer")
    if Path(command[0]).resolve(strict=True) != Path(sys.executable).resolve(strict=True):
        raise RuntimeError("Source-lock wrapper must use the current Python interpreter")
    locked_paths = [producer, *(path.resolve(strict=True) for path in args.locked_file)]
    if len(set(locked_paths)) != len(locked_paths):
        raise RuntimeError("Locked source file list contains duplicates")
    before = {str(path): identity(path) for path in locked_paths}
    wrapper_before = identity(Path(__file__))
    started_at = now_utc()
    completed = subprocess.run(command, cwd=workdir, check=False)
    finished_at = now_utc()
    after = {str(path): identity(path) for path in locked_paths}
    wrapper_after = identity(Path(__file__))
    if completed.returncode != 0:
        raise RuntimeError(f"Screen producer exited with code {completed.returncode}")
    if before != after:
        raise RuntimeError("A locked source file changed while the screen was running")
    if wrapper_before != wrapper_after:
        raise RuntimeError("Source-lock wrapper changed while the screen was running")
    receipt_path, child_receipt = validate_child_receipt(
        output_dir, command, before[str(producer)]
    )
    source_lock_path = output_dir / "source_lock_receipt.json"
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": finished_at,
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "task_type": "B_comprehensive_validation",
        "scope": child_receipt["execution"]["selected_samples"],
        "command": command,
        "workdir": str(workdir),
        "computational_thread_environment": {
            key: os.environ.get(key) for key in THREAD_ENV_KEYS
        },
        "source_identity_before": before,
        "source_identity_after": after,
        "wrapper_identity_before": wrapper_before,
        "wrapper_identity_after": wrapper_after,
        "child_run_manifest": artifact(receipt_path),
        "checks": {
            "child_exit_zero": True,
            "all_locked_sources_unchanged": True,
            "wrapper_unchanged": True,
            "child_receipt_pass": True,
            "child_receipt_command_exact": True,
            "child_receipt_analyzer_matches_start_lock": True,
        },
        "pass_semantics": "execution_and_source_identity_only_not_scientific_confirmation",
        "pass": True,
    }
    with source_lock_path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(
        json.dumps(
            {
                "source_lock_receipt": str(source_lock_path),
                "selected_samples": payload["scope"],
                "analyzer_sha256": before[str(producer)]["sha256"],
                "pass": True,
            },
            indent=2,
        )
    )
    return source_lock_path


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
