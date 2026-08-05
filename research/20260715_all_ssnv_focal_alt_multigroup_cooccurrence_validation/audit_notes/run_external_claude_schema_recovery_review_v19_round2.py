#!/usr/bin/env python3
"""Run and preserve the read-only external Claude v19 round-2 review."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
AUDIT_ROOT = (
    REPO_ROOT
    / "research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
    / "audit_notes"
)
CLAUDE = Path("/bip7_disk/liaoyoyo2001/.local/bin/claude")
PROMPT = AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v19_round2_prompt.md"
SCHEMA = AUDIT_ROOT / "20260723_external_claude_schema_recovery_review_v19_round2.schema.json"
ENVELOPE = (
    AUDIT_ROOT
    / "20260723_external_claude_schema_recovery_review_v19_round2.claude_cli.envelope.json"
)
STDERR = (
    AUDIT_ROOT
    / "20260723_external_claude_schema_recovery_review_v19_round2.claude_cli.stderr.txt"
)
SESSION_ID = "db7885d1-6d8c-4bdb-8fcc-7bb3f6f03124"


def write_exclusive_readonly(path: Path, data: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
    try:
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise RuntimeError(f"Short write: {path}")
            view = view[written:]
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
    finally:
        os.close(descriptor)


def main() -> int:
    for output in (ENVELOPE, STDERR):
        if os.path.lexists(output):
            raise RuntimeError(f"External-review output already exists: {output}")
    prompt = PROMPT.read_text(encoding="utf-8")
    schema = SCHEMA.read_text(encoding="utf-8")
    json.loads(schema)
    allowed_tools = ",".join(
        (
            "Read",
            "Glob",
            "Grep",
            "Bash(/usr/bin/env *)",
            "Bash(sha256sum *)",
            "Bash(stat *)",
            "Bash(realpath *)",
            "Bash(ls *)",
            "Bash(find *)",
            "Bash(rg *)",
            "Bash(sed *)",
            "Bash(wc *)",
            "Bash(head *)",
            "Bash(tail *)",
            "Bash(sort *)",
            "Bash(cut *)",
        )
    )
    system_prompt = (
        "Act only as an independent read-only reviewer. Never mutate files, run "
        "producers, install software, or use network tools. Use only bounded reads "
        "and the explicitly requested read-only probe. Return only schema-valid JSON."
    )
    command = [
        str(CLAUDE),
        "--print",
        "--safe-mode",
        "--disable-slash-commands",
        "--no-chrome",
        "--strict-mcp-config",
        "--mcp-config",
        '{"mcpServers":{}}',
        "--permission-mode",
        "dontAsk",
        "--tools",
        "Read,Glob,Grep,Bash",
        "--allowedTools",
        allowed_tools,
        "--disallowedTools",
        "Edit,Write,NotebookEdit,WebFetch,WebSearch",
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
        schema,
        "--session-id",
        SESSION_ID,
        "--name",
        "InterSubMod Task-B schema recovery external reviewer v19 round2",
        "--append-system-prompt",
        system_prompt,
        prompt,
    ]
    result = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env=os.environ.copy(),
        capture_output=True,
        check=False,
        timeout=1800,
    )
    write_exclusive_readonly(ENVELOPE, result.stdout)
    write_exclusive_readonly(STDERR, result.stderr)
    if result.returncode != 0:
        return result.returncode
    envelope = json.loads(result.stdout)
    structured = envelope.get("structured_output")
    if (
        envelope.get("session_id") != SESSION_ID
        or not isinstance(structured, dict)
        or structured.get("reviewer_agent_id") != SESSION_ID
    ):
        raise RuntimeError("External Claude envelope/session binding drift")
    print(
        json.dumps(
            {
                "envelope": str(ENVELOPE),
                "stderr": str(STDERR),
                "session_id": envelope.get("session_id"),
                "verdict": structured.get("verdict"),
                "pass": structured.get("pass"),
                "high_findings": len(structured.get("high_findings", [])),
                "medium_findings": len(structured.get("medium_findings", [])),
                "unresolved_conditions": len(
                    structured.get("unresolved_conditions", [])
                ),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
