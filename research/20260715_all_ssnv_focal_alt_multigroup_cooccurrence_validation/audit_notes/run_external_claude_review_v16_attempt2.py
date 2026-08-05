#!/usr/bin/env python3
"""Run one read-only external Claude attempt2 review with bounded shell usage."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import types
import uuid
from typing import Any


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
HERE = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "audit_notes"
)
BASE_RUNNER = HERE / "run_external_claude_review_v15.py"
BASE_RUNNER_SHA256 = (
    "97f497695f869f5f1308220ffde67de75b81711d52c0b380b003d5d9716b8556"
)
SCHEMA = HERE / "external_source_review_v15_output.schema.json"
SCHEMA_SHA256 = "9a51390058a157afc1355c766d007cef51696b6f83e0dfb37c8921b50f233d85"
CLAUDE = Path("/bip7_disk/liaoyoyo2001/.local/bin/claude").resolve(strict=True)
REVIEW_SPECS = {
    "A": {
        "request": (
            HERE
            / "20260718_external_claude_reviewer_A_v15_v4_signer_request.md"
        ),
        "request_sha256": (
            "21079007398f13d43648221439641a75ff2bf3adf2d0b1a90ffe8d9e451802fd"
        ),
        "output_prefix": (
            "20260718_external_claude_reviewer_A_v15_v4_signer_output.attempt2"
        ),
        "prior_attempt": {
            "status": "TRANSPORT_FAILURE_SESSION_LIMIT_NOT_A_REVIEW_VERDICT",
            "returncode": 1,
            "result": (
                "You've hit your session limit · resets 12:40am (Asia/Taipei)"
            ),
            "command_sha256": (
                "f304dfa9693dcc1ea39bb701921a92ac098ba20c7d9d8ff9ac3dcdc5a2fe4aaa"
            ),
            "envelope_sha256": (
                "608c3689fa40ece5976e3eab3484568234f679c19b4a3249e8e1dcd5e11d3811"
            ),
            "stderr_sha256": (
                "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a"
            ),
        },
    },
    "B": {
        "request": (
            HERE
            / "20260718_external_claude_reviewer_B_v15_adversarial_request.md"
        ),
        "request_sha256": (
            "cf00effd120069ace81e4af770609aa3d2c3022569f5cdd858ce7aec2030e345"
        ),
        "output_prefix": (
            "20260718_external_claude_reviewer_B_v15_v4_signer_output.attempt2"
        ),
        "prior_attempt": {
            "status": "TRANSPORT_FAILURE_SESSION_LIMIT_NOT_A_REVIEW_VERDICT",
            "returncode": 1,
            "result": (
                "You've hit your session limit · resets 12:40am (Asia/Taipei)"
            ),
            "command_sha256": (
                "af3f4077bcd1aba51a905699916351f2678bb3c18543eab22b5d1b48eb315268"
            ),
            "envelope_sha256": (
                "a44e204ae9b75c413de16e821bedd8b370810c56c99b8b6b584c51b630ab530f"
            ),
            "stderr_sha256": (
                "222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a"
            ),
        },
    },
}
ALLOWED_TOOLS = ",".join(
    (
        "Read",
        "Glob",
        "Grep",
        "Bash(git *)",
        "Bash(sha256sum *)",
        "Bash(stat *)",
        "Bash(realpath *)",
        "Bash(ls *)",
        "Bash(find *)",
        "Bash(rg *)",
        "Bash(sed *)",
        "Bash(wc *)",
        "Bash(printf *)",
        "Bash(test *)",
        "Bash(head *)",
        "Bash(tail *)",
        "Bash(sort *)",
        "Bash(cut *)",
    )
)
DISALLOWED_TOOLS = ",".join(
    (
        "Edit",
        "Write",
        "NotebookEdit",
        "WebFetch",
        "WebSearch",
    )
)


class ReviewRunnerError(RuntimeError):
    """Raised when an attempt2 external review cannot be preserved safely."""


def load_bound_base_runner() -> tuple[types.ModuleType, dict[str, Any]]:
    resolved = BASE_RUNNER.resolve(strict=True)
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(resolved, flags)
    try:
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            raise ReviewRunnerError(f"Base runner is not regular: {resolved}")
        source = b"".join(
            os.pread(descriptor, min(1024 * 1024, opened.st_size - offset), offset)
            for offset in range(0, opened.st_size, 1024 * 1024)
        )
        closed_read = os.fstat(descriptor)
        live = resolved.stat()
        identity = (
            opened.st_dev,
            opened.st_ino,
            opened.st_mode,
            opened.st_size,
            opened.st_mtime_ns,
            opened.st_ctime_ns,
        )
        if (
            len(source) != opened.st_size
            or identity
            != (
                closed_read.st_dev,
                closed_read.st_ino,
                closed_read.st_mode,
                closed_read.st_size,
                closed_read.st_mtime_ns,
                closed_read.st_ctime_ns,
            )
            or identity
            != (
                live.st_dev,
                live.st_ino,
                live.st_mode,
                live.st_size,
                live.st_mtime_ns,
                live.st_ctime_ns,
            )
        ):
            raise ReviewRunnerError(f"Base runner changed while reading: {resolved}")
        digest = hashlib.sha256(source).hexdigest()
        if digest != BASE_RUNNER_SHA256 or stat.S_IMODE(opened.st_mode) != 0o444:
            raise ReviewRunnerError(f"Base runner identity drifted: {resolved}")
    finally:
        os.close(descriptor)

    module = types.ModuleType("external_review_v15_bound")
    module.__file__ = str(resolved)
    exec(compile(source, str(resolved), "exec"), module.__dict__)
    return module, {
        "path": str(resolved),
        "size_bytes": len(source),
        "sha256": digest,
        "mode": oct(stat.S_IMODE(opened.st_mode)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reviewer", choices=("A", "B"), required=True)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base, base_record = load_bound_base_runner()
    spec = REVIEW_SPECS[args.reviewer]
    schema_bytes, schema_record = base.read_bound(SCHEMA, SCHEMA_SHA256)
    request_bytes, request_record = base.read_bound(
        spec["request"],
        spec["request_sha256"],
    )
    try:
        schema_text = schema_bytes.decode("utf-8")
        request_text = request_bytes.decode("utf-8")
        json.loads(schema_text)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ReviewRunnerError("Schema or request encoding is invalid") from exc

    session_id = str(uuid.uuid4())
    prefix = str(spec["output_prefix"])
    envelope_path = HERE / f"{prefix}.claude_cli.envelope.json"
    stderr_path = HERE / f"{prefix}.claude_cli.stderr.txt"
    command_path = HERE / f"{prefix}.claude_cli.command.json"
    for output_path in (envelope_path, stderr_path, command_path):
        if os.path.lexists(output_path):
            raise ReviewRunnerError(f"Review output slot is occupied: {output_path}")

    system_prompt = (
        "Act only as an independent read-only reviewer. Never mutate files, "
        "run tests or producers, install software, or use web/network tools. "
        "Use only bounded reads and explicitly allowed non-mutating checks. "
        "Do not use Python, cd, shell loops, pipelines, redirection, command "
        "substitution, semicolons, &&, or ||. Every Bash tool call must contain "
        "exactly one primitive permitted command and use absolute paths. Prefer "
        "Read, Glob, and Grep for inspection. Return the exact structured JSON "
        "contract requested."
    )
    command = [
        str(CLAUDE),
        "--print",
        "--safe-mode",
        "--disable-slash-commands",
        "--no-chrome",
        "--strict-mcp-config",
        "--mcp-config",
        "{\"mcpServers\":{}}",
        "--permission-mode",
        "dontAsk",
        "--tools",
        "Read,Glob,Grep,Bash",
        "--allowedTools",
        ALLOWED_TOOLS,
        "--disallowedTools",
        DISALLOWED_TOOLS,
        "--add-dir",
        "/bip7_disk/liaoyoyo2001/.config",
        "--model",
        "opus",
        "--effort",
        "max",
        "--max-budget-usd",
        "10",
        "--output-format",
        "json",
        "--json-schema",
        schema_text,
        "--session-id",
        session_id,
        "--name",
        f"InterSubMod Task-B external Reviewer {args.reviewer} v15 attempt2",
        "--append-system-prompt",
        system_prompt,
        request_text,
    ]
    command_record = {
        "reviewer": args.reviewer,
        "attempt": 2,
        "session_id": session_id,
        "cwd": str(REPO_ROOT),
        "claude_path": str(CLAUDE),
        "claude_sha256": hashlib.sha256(CLAUDE.read_bytes()).hexdigest(),
        "base_runner": base_record,
        "request": request_record,
        "schema": schema_record,
        "output_envelope": str(envelope_path),
        "output_stderr": str(stderr_path),
        "model": "opus",
        "effort": "max",
        "permission_mode": "dontAsk",
        "allowed_tools": ALLOWED_TOOLS,
        "disallowed_tools": DISALLOWED_TOOLS,
        "system_prompt": system_prompt,
        "formal_artifact": False,
        "requires_transport_extraction_and_v15_normalization": True,
        "prior_attempt": spec["prior_attempt"],
    }
    if args.dry_run:
        print(json.dumps(command_record, ensure_ascii=False, indent=2, sort_keys=True))
        return 0

    command_bytes = (
        json.dumps(command_record, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    command_output_record = base.write_exclusive_readonly(command_path, command_bytes)
    environment = dict(os.environ)
    environment.pop("PYTHONPATH", None)
    environment["NO_COLOR"] = "1"
    environment["TERM"] = "dumb"
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        check=False,
        timeout=3600,
    )
    envelope_record = base.write_exclusive_readonly(envelope_path, completed.stdout)
    stderr_record = base.write_exclusive_readonly(stderr_path, completed.stderr)

    envelope_summary: dict[str, Any] = {}
    try:
        envelope = json.loads(completed.stdout.decode("utf-8"))
        envelope_summary = {
            "type": envelope.get("type"),
            "subtype": envelope.get("subtype"),
            "is_error": envelope.get("is_error"),
            "session_id": envelope.get("session_id"),
            "structured_output_type": type(
                envelope.get("structured_output")
            ).__name__,
            "permission_denials": len(envelope.get("permission_denials", [])),
        }
    except (UnicodeDecodeError, json.JSONDecodeError):
        envelope_summary = {"parseable_json_envelope": False}
    print(
        json.dumps(
            {
                "reviewer": args.reviewer,
                "returncode": completed.returncode,
                "command_record": command_output_record,
                "envelope": envelope_record,
                "stderr": stderr_record,
                "envelope_summary": envelope_summary,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return completed.returncode


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.TimeoutExpired as exc:
        print(f"EXTERNAL_REVIEW_TIMEOUT error={exc}", file=sys.stderr)
        raise SystemExit(124)
