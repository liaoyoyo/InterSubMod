#!/usr/bin/env python3
"""Run the packaged Data Analytics portable builder with no-clobber receipts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


TOPIC_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PLUGIN_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/"
    "data-analytics/0.2.8-13ceeea1f599"
)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_json(path: Path, payload: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def ensure_topic_path(path: Path) -> None:
    try:
        path.resolve().relative_to(TOPIC_ROOT.resolve())
    except ValueError as exc:
        raise ValueError(f"output path must be within {TOPIC_ROOT}") from exc


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--plugin-root", type=Path, default=DEFAULT_PLUGIN_ROOT)
    parser.add_argument("--failure-screenshot", type=Path)
    parser.add_argument("--ready-timeout-ms", type=int, default=5_000)
    parser.add_argument("--action-timeout-ms", type=int, default=2_500)
    parser.add_argument("--timeout-ms", type=int, default=10_000)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    for path in (args.html, args.receipt):
        ensure_topic_path(path)
        if path.exists():
            raise FileExistsError(path)
    if args.failure_screenshot is not None:
        ensure_topic_path(args.failure_screenshot)
        if args.failure_screenshot.exists():
            raise FileExistsError(args.failure_screenshot)
    if not args.artifact.is_file():
        raise FileNotFoundError(args.artifact)
    for label, value in (
        ("ready-timeout-ms", args.ready_timeout_ms),
        ("action-timeout-ms", args.action_timeout_ms),
        ("timeout-ms", args.timeout_ms),
    ):
        if not 1_000 <= value <= 60_000:
            raise ValueError(f"{label} must be between 1000 and 60000")

    builder = args.plugin_root / "skills" / "build-report" / "scripts" / "deliver_portable_artifact.mjs"
    if not builder.is_file():
        raise FileNotFoundError(builder)
    command = [
        "node",
        str(builder.resolve()),
        "--input",
        str(args.artifact.resolve()),
        "--output",
        str(args.html.resolve()),
        "--ready-timeout-ms",
        str(args.ready_timeout_ms),
        "--action-timeout-ms",
        str(args.action_timeout_ms),
        "--timeout-ms",
        str(args.timeout_ms),
    ]
    if args.failure_screenshot is not None:
        command.extend(("--screenshot", str(args.failure_screenshot.resolve())))
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    created_at = datetime.now(timezone.utc).isoformat()
    try:
        builder_payload = json.loads(completed.stdout) if completed.stdout.strip() else None
    except json.JSONDecodeError:
        builder_payload = None
    receipt = {
        "schema_name": "intersubmod.singleton_sidecar_portable_builder_receipt",
        "schema_version": "1.0.0",
        "created_at_utc": created_at,
        "ok": completed.returncode == 0,
        "command": command,
        "returncode": completed.returncode,
        "builder": {
            "path": str(builder.resolve()),
            "sha256": sha256_path(builder),
            "plugin_root": str(args.plugin_root.resolve()),
        },
        "browser_environment": {
            "CHROMIUM_EXECUTABLE_PATH": os.environ.get("CHROMIUM_EXECUTABLE_PATH"),
        },
        "artifact": {
            "path": str(args.artifact.resolve()),
            "sha256": sha256_path(args.artifact),
        },
        "stdout_payload": builder_payload,
        "stdout": completed.stdout if builder_payload is None else "",
        "stderr": completed.stderr,
    }
    if completed.returncode == 0:
        if not args.html.is_file():
            receipt["ok"] = False
            receipt["error"] = "portable builder returned zero without creating HTML"
        else:
            receipt["html"] = {
                "path": str(args.html.resolve()),
                "sha256": sha256_path(args.html),
                "size_bytes": args.html.stat().st_size,
            }
    safe_json(args.receipt, receipt)
    args.receipt.chmod(0o444)
    if args.html.exists():
        args.html.chmod(0o444)
    if args.failure_screenshot is not None and args.failure_screenshot.exists():
        args.failure_screenshot.chmod(0o444)
    print(json.dumps(receipt, ensure_ascii=False, sort_keys=True))
    return 0 if receipt["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
