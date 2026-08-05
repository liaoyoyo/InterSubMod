#!/usr/bin/env python3
"""Validate signed test-run and external-review process evidence."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping
import uuid
import xml.etree.ElementTree as ET


TOPIC_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = TOPIC_ROOT.parents[1]
HERE = TOPIC_ROOT / "audit_notes"
TEST_ROOT = TOPIC_ROOT / "tests"
TEST_XML_PATH = (
    TOPIC_ROOT
    / "logs"
    / "pytest_full_pre_authority_v27_git_env_alignment_attested_canonical.xml"
)
RAW_TEST_XML_PATH = (
    TOPIC_ROOT
    / "logs"
    / "pytest_full_pre_authority_v27_git_env_alignment_raw_evidence"
    / "pytest.xml"
)
BOUND_TEST_DIR = TOPIC_ROOT / ".pytest_full_pre_authority_v27_git_env_alignment_bound_test_fds"
TEST_STDOUT_PATH = TEST_XML_PATH.with_suffix(".stdout.txt")
TEST_STDERR_PATH = TEST_XML_PATH.with_suffix(".stderr.txt")
TEST_MANIFEST_PATH = (
    TOPIC_ROOT
    / "logs"
    / "pytest_full_pre_authority_v27_git_env_alignment_test_run_manifest.v1.json"
)
TEST_MANIFEST_SIGNATURE_PATH = Path(str(TEST_MANIFEST_PATH) + ".ed25519.sig")
TEST_MANIFEST_PUBLIC_KEY_PATH = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_test_run_authority/"
    "20260722_all_ssnv_v11_git_env_alignment_bootstrap/ed25519_public.pem"
)
TEST_MANIFEST_PRIVATE_KEY_PATH = TEST_MANIFEST_PUBLIC_KEY_PATH.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
TEST_MANIFEST_PUBLIC_KEY_SHA256 = (
    "92255c453329b8b21dda721048259e8cd567addb5b3f0ee2aa03fb414a9d5015"
)
TEST_RUNNER_PATH = HERE / "run_attested_canonical_pytest_v5.py"
TEST_PYTHON_PATH = Path(
    "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"
)
TEST_PYTHON_SHA256 = (
    "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e"
)
EXPECTED_TESTCASE_COUNT = 743
EXPECTED_TESTCASE_INVENTORY_SHA256 = (
    "831f0f9a0ea4b0b1490fadcedc95bcfeb1e458ea159aa71f79bb064843447c22"
)
AUTHORITY_ASSEMBLER_PATH = HERE / "assemble_release_source_authority_v15.py"
REVIEW_NORMALIZER_PATH = HERE / "normalize_external_source_review_v20.py"
REVIEW_RUNNER_PATH = HERE / "run_external_claude_review_v22_attested.py"
SIGNER_PATH = HERE / "run_one_time_ed25519_signer_v6.py"
EVIDENCE_VALIDATOR_PATH = Path(__file__).resolve()
CLAUDE_CLI_PATH = Path(
    "/bip7_disk/liaoyoyo2001/.local/share/claude/versions/2.1.202"
)
CLAUDE_CLI_SHA256 = (
    "71590202249892db3805ecd5b867f831f04b8129eaabd3f9a5bd4ba16b52c839"
)
EXPECTED_CLAUDE_MODEL_ID = "claude-opus-4-8"
REQUIRED_REVIEW_READ_PATHS = (
    TOPIC_ROOT / "scripts" / "release_source_authority.py",
    HERE / "attested_release_evidence_v5.py",
    HERE / "run_attested_canonical_pytest_v5.py",
    HERE / "run_external_claude_review_v22_attested.py",
    HERE / "normalize_external_source_review_v20.py",
    HERE / "assemble_release_source_authority_v15.py",
    HERE / "run_one_time_ed25519_signer_v6.py",
)
REVIEW_ATTESTATION_SPECS = {
    "A": {
        "path": (
            HERE
            / "20260722_external_claude_reviewer_A_v21_strict_command_parity_process_attestation.v1.json"
        ),
        "public_key": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_external_review_attestation/"
            "20260722_reviewer_A_v22_strict_command_parity_bootstrap/ed25519_public.pem"
        ),
        "public_key_sha256": (
            "cbf3ad997f9b4f0a174d2ba0482e35334091ea717464e89ac22175481581bd80"
        ),
    },
    "B": {
        "path": (
            HERE
            / "20260722_external_claude_reviewer_B_v21_strict_command_parity_process_attestation.v1.json"
        ),
        "public_key": Path(
            "/bip7_disk/liaoyoyo2001/.config/"
            "intersubmod_external_review_attestation/"
            "20260722_reviewer_B_v22_strict_command_parity_bootstrap/ed25519_public.pem"
        ),
        "public_key_sha256": (
            "e4396b37bf5c04b9f61b91028bdd06fa6639503abe58def649cbd60656d8b74a"
        ),
    },
}
for _spec in REVIEW_ATTESTATION_SPECS.values():
    _spec["signature"] = Path(str(_spec["path"]) + ".ed25519.sig")
    _spec["private_key"] = _spec["public_key"].with_name(
        "ed25519_private_encrypted_unrecoverable_after_signing.pem"
    )

ARTIFACT_KEYS = frozenset({"path", "size_bytes", "sha256", "mode"})
JUNIT_COUNT_KEYS = frozenset({"tests", "failures", "errors", "skipped"})
REVIEW_PAYLOAD_FIELD_ORDER = (
    "schema_name",
    "schema_version",
    "reviewer_label",
    "reviewer_id",
    "model",
    "verdict",
    "findings_closed",
    "f1_status",
    "f2_status",
    "reviewed_source_set_sha256",
    "reviewed_git_head",
    "reviewed_assembler",
    "canonical_junit",
    "review_scope",
    "blocking_findings",
    "nonblocking_findings",
    "evidence",
)
REVIEW_PAYLOAD_KEYS = frozenset(REVIEW_PAYLOAD_FIELD_ORDER)
TEST_MANIFEST_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "task_type",
        "scope",
        "run_id",
        "started_at_utc",
        "finished_at_utc",
        "cwd",
        "command",
        "environment",
        "git_head",
        "protected_sources",
        "source_set_sha256",
        "test_sources",
        "supporting_sources",
        "runtime",
        "junit",
        "captured_output",
        "signature_contract",
        "checks",
        "pass",
    }
)
TEST_ENVIRONMENT = {
    "PATH": "/usr/bin:/bin",
    "HOME": "/tmp",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PYTHONHASHSEED": "0",
    "PYTHONNOUSERSITE": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1",
    "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1",
}
ATTESTATION_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "reviewer_token",
        "attestation_id",
        "created_at_utc",
        "process",
        "review_payload",
        "review_payload_sha256",
        "reviewed_test_run",
        "signature_contract",
        "checks",
        "pass",
    }
)
PROCESS_KEYS = frozenset(
    {
        "session_id",
        "entrypoint",
        "claude_cli",
        "command_record",
        "request",
        "transcript",
        "stderr",
        "tool_trace",
        "returncode",
        "model",
        "effort",
        "permission_mode",
        "allowed_tools",
        "disallowed_tools",
        "environment",
    }
)
COMMAND_RECORD_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "reviewer_token",
        "session_id",
        "cwd",
        "entrypoint",
        "claude_cli",
        "runner",
        "evidence_validator",
        "request",
        "model",
        "effort",
        "permission_mode",
        "allowed_tools",
        "disallowed_tools",
        "environment",
        "schema_sha256",
        "command",
    }
)
BOOTSTRAP_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "argv0",
        "script",
        "python",
        "isolated_mode",
        "no_user_site",
        "environment",
    }
)
ALLOWED_REVIEW_TOOLS = "Read,Glob,Grep,StructuredOutput"
DISALLOWED_REVIEW_TOOLS = (
    "Bash,Edit,Write,NotebookEdit,WebFetch,WebSearch"
)
CLEAN_REVIEW_ENVIRONMENT = {
    "PATH": "/usr/bin:/bin",
    "HOME": "/bip7_disk/liaoyoyo2001",
    "USER": "liaoyoyo2001",
    "LOGNAME": "liaoyoyo2001",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "TERM": "dumb",
    "NO_COLOR": "1",
    "PYTHONHASHSEED": "0",
    "PYTHONNOUSERSITE": "1",
}
REVIEW_SYSTEM_PROMPT = (
    "Act only as an independent read-only source and release-chain reviewer. "
    "Never mutate files, run tests/producers, install software, or use network "
    "tools. Use only Read, Glob, Grep, and the required StructuredOutput "
    "finalizer. "
    "Return exactly the requested structured JSON."
)


class EvidenceValidationError(RuntimeError):
    """Raised when signed process evidence is incomplete or inconsistent."""


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise EvidenceValidationError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_nonfinite_constant(value: str) -> None:
    raise EvidenceValidationError(f"Non-finite JSON constant: {value}")


def parse_json_bytes(data: bytes, *, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(
            data.decode("utf-8"),
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=reject_nonfinite_constant,
        )
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        EvidenceValidationError,
    ) as error:
        raise EvidenceValidationError(f"{label} is malformed") from error
    if not isinstance(payload, dict):
        raise EvidenceValidationError(f"{label} must be a JSON object")
    return payload


def encode_json(payload: Any) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        + "\n"
    ).encode("utf-8")


def artifact_schema(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "type": "object",
        "additionalProperties": False,
        "required": ["path", "size_bytes", "sha256", "mode"],
        "properties": {
            field: {"const": record[field]}
            for field in ("path", "size_bytes", "sha256", "mode")
        },
    }


def expected_review_schema(
    *,
    source_set_sha256: str,
    git_head: str,
    assembler: Mapping[str, Any],
    junit: Mapping[str, Any],
    counts: Mapping[str, int],
) -> dict[str, Any]:
    # Claude Code can demote StructuredOutput when a transport schema is deeply
    # nested. Exact identities and approval semantics remain enforced below by
    # validate_clean_review_payload(), independent of this shallow transport.
    _ = (source_set_sha256, git_head, assembler, junit, counts)
    field_types = {
        "schema_name": "string",
        "schema_version": "string",
        "reviewer_label": "string",
        "reviewer_id": "string",
        "model": "string",
        "verdict": "string",
        "findings_closed": "boolean",
        "f1_status": "string",
        "f2_status": "string",
        "reviewed_source_set_sha256": "string",
        "reviewed_git_head": "string",
        "reviewed_assembler": "object",
        "canonical_junit": "object",
        "review_scope": "string",
        "blocking_findings": "array",
        "nonblocking_findings": "array",
        "evidence": "object",
    }
    return {
        "type": "object",
        "additionalProperties": False,
        "required": list(REVIEW_PAYLOAD_FIELD_ORDER),
        "properties": {
            field: {"type": field_types[field]}
            for field in REVIEW_PAYLOAD_FIELD_ORDER
        },
    }


def expected_review_request(
    reviewer_token: str,
    *,
    contract: Mapping[str, Any],
) -> bytes:
    emphasis = (
        "Trace every authorization claim to executable checks and look for "
        "coherent substitution, TOCTOU, environment injection, or false-pass paths."
        if reviewer_token == "A"
        else
        "Act adversarially: attempt to falsify source, test, signer, review-process, "
        "private-key retirement, and old-v4-revocation claims."
    )
    contract_text = json.dumps(
        contract,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    required_reads = "\n".join(
        f"- {path.resolve()}" for path in REQUIRED_REVIEW_READ_PATHS
    )
    return (
        f"# External Reviewer {reviewer_token} v21 final git-env-alignment review\n\n"
        "Read-only independent review. Do not modify files or run tests/producers.\n\n"
        f"{emphasis}\n\n"
        "Review all protected sources, the attested canonical pytest manifest, "
        "the v15 assembler, v20 normalizer/v22 runner, signed-evidence validation, "
        "FD lifetime, clean subprocess environment, supplemental ABA/post-sign "
        "verification, the corrected v7 authority consumers, and block-specific "
        "Git/OpenSSL subprocess environment isolation.\n\n"
        "Security boundary: this local evidence chain covers persistent or "
        "before/after-observable drift, no-clobber publication, exact identities, "
        "and reproducible process contracts. Concurrent same-UID/root/ptrace "
        "compromise or live mutation of the Python interpreter/site-packages is "
        "out of scope and must be reported as a residual limitation, not treated "
        "as a release blocker by itself.\n\n"
        "Return APPROVE only if no HIGH/MEDIUM false-authority path remains and "
        "F1/F2 are RESOLVED_VERIFIED. Otherwise return REQUEST_CHANGES with "
        "specific file:line evidence. StructuredOutput only transports the "
        "17 top-level fields; the independent local consumer revalidates every "
        "exact identity and approval rule.\n\n"
        "StructuredOutput must contain exactly these top-level fields:\n"
        f"{','.join(REVIEW_PAYLOAD_FIELD_ORDER)}\n"
        "For APPROVE use schema_name=intersubmod.external_claude_source_review, "
        "schema_version=1.0.0, model=claude-opus-4-8, findings_closed=true, "
        "f1_status=RESOLVED_VERIFIED, f2_status=RESOLVED_VERIFIED, and an empty "
        "blocking_findings array. reviewer_label must identify this Reviewer "
        f"{reviewer_token}; reviewer_id must be a fresh canonical UUIDv4. Copy "
        "reviewed_source_set_sha256, reviewed_git_head, reviewed_assembler, and "
        "canonical_junit exactly from REVIEW_CONTRACT_JSON. review_scope must be "
        "at least 80 characters and evidence must contain at least 3 properties. "
        "Any mismatch is rejected locally even if transport validation passes.\n\n"
        "You MUST invoke Read on every path below before deciding. Read each whole "
        "file with only file_path: do not set offset or limit. Use Glob/Grep for "
        "broader coverage as needed. The signed transcript is rejected if any "
        "required Read is absent, bounded, empty, or returns an error:\n"
        f"{required_reads}\n\n"
        f"REVIEW_CONTRACT_JSON={contract_text}\n"
    ).encode("utf-8")


def canonical_uuid4(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    try:
        parsed = uuid.UUID(value)
    except ValueError:
        return False
    return parsed.version == 4 and str(parsed) == value


def parse_utc_timestamp(value: Any, *, label: str) -> datetime:
    if not isinstance(value, str):
        raise EvidenceValidationError(f"{label} is not an ISO timestamp")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as error:
        raise EvidenceValidationError(f"{label} is not an ISO timestamp") from error
    if parsed.tzinfo is None or parsed.utcoffset() != timezone.utc.utcoffset(parsed):
        raise EvidenceValidationError(f"{label} is not UTC")
    return parsed


def require_artifact_shape(
    value: Any,
    *,
    expected_path: Path | None = None,
    label: str,
) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != set(ARTIFACT_KEYS):
        raise EvidenceValidationError(f"{label} is not an exact artifact identity")
    record = dict(value)
    sha256 = record.get("sha256")
    if (
        (expected_path is not None and record.get("path") != str(expected_path.resolve()))
        or type(record.get("size_bytes")) is not int
        or record["size_bytes"] < 0
        or not isinstance(sha256, str)
        or len(sha256) != 64
        or any(character not in "0123456789abcdef" for character in sha256)
        or record.get("mode") != "0o444"
    ):
        raise EvidenceValidationError(f"{label} identity is malformed")
    return record


def source_set_digest(sources: Mapping[str, Mapping[str, Any]]) -> str:
    lines = [
        "|".join(
            (
                role,
                str(record["path"]),
                str(record["size_bytes"]),
                str(record["sha256"]),
                str(record["mode"]),
            )
        )
        for role, record in sorted(sources.items())
    ]
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def parse_junit_counts(data: bytes) -> dict[str, int]:
    try:
        root = ET.fromstring(data)
    except ET.ParseError as error:
        raise EvidenceValidationError("Canonical JUnit XML is malformed") from error
    suites = [node for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "testsuite"]
    if not suites:
        raise EvidenceValidationError("Canonical JUnit XML has no testsuite")
    attributes = {
        key: sum(int(suite.attrib.get(key, "0")) for suite in suites)
        for key in JUNIT_COUNT_KEYS
    }
    elements = {
        "tests": sum(
            1 for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "testcase"
        ),
        "failures": sum(
            1 for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "failure"
        ),
        "errors": sum(
            1 for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "error"
        ),
        "skipped": sum(
            1 for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "skipped"
        ),
    }
    if attributes != elements:
        raise EvidenceValidationError(
            f"JUnit attributes differ from elements: {attributes} != {elements}"
        )
    return attributes


def parse_junit_inventory(data: bytes) -> dict[str, Any]:
    try:
        root = ET.fromstring(data)
    except ET.ParseError as error:
        raise EvidenceValidationError("Canonical JUnit XML is malformed") from error
    cases = []
    for node in root.iter():
        if node.tag.rsplit("}", 1)[-1] != "testcase":
            continue
        cases.append(
            {
                "classname": node.attrib.get("classname", ""),
                "name": node.attrib.get("name", ""),
                "file": node.attrib.get("file", ""),
                "line": node.attrib.get("line", ""),
            }
        )
    canonical = json.dumps(
        sorted(cases, key=lambda item: tuple(item.values())),
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return {
        "count": len(cases),
        "sha256": hashlib.sha256(canonical).hexdigest(),
    }


def validate_entrypoint_record(
    value: Any,
    *,
    expected_script: Mapping[str, Any],
    expected_python: Mapping[str, Any],
    expected_environment: Mapping[str, str],
    label: str,
) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != set(BOOTSTRAP_KEYS):
        raise EvidenceValidationError(f"{label} entrypoint schema is not exact")
    record = dict(value)
    argv0 = record.get("argv0")
    if (
        record.get("schema_name") != "intersubmod.bound_python_entrypoint"
        or record.get("schema_version") != "1.0.0"
        or not isinstance(argv0, str)
        or not argv0.startswith("/proc/self/fd/")
        or not argv0.removeprefix("/proc/self/fd/").isdigit()
        or record.get("script") != dict(expected_script)
        or record.get("python") != dict(expected_python)
        or record.get("isolated_mode") is not True
        or record.get("no_user_site") is not True
        or record.get("environment") != dict(expected_environment)
    ):
        raise EvidenceValidationError(f"{label} bound entrypoint drift")
    return record


def expected_test_runtime_record(
    *,
    live_python: Mapping[str, Any],
    live_git: Mapping[str, Any],
    entrypoint: Mapping[str, Any],
    pytest_process: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "python": dict(live_python),
        "git": dict(live_git),
        "entrypoint": dict(entrypoint),
        "pytest_process": dict(pytest_process),
    }


def parse_claude_stream_transcript(
    data: bytes,
    *,
    reviewer_token: str,
    session_id: str,
) -> dict[str, Any]:
    events: list[dict[str, Any]] = []
    for line_number, raw_line in enumerate(data.splitlines(), 1):
        if not raw_line.strip():
            continue
        event = parse_json_bytes(
            raw_line,
            label=f"Claude transcript line {line_number}",
        )
        if not isinstance(event, Mapping):
            raise EvidenceValidationError(
                f"Claude transcript line {line_number} is not an object"
            )
        event = dict(event)
        events.append(event)
    if not events:
        raise EvidenceValidationError("Claude transcript is empty")

    init_events = [
        event
        for event in events
        if event.get("type") == "system" and event.get("subtype") == "init"
    ]
    result_events = [event for event in events if event.get("type") == "result"]
    if len(init_events) != 1 or len(result_events) != 1:
        raise EvidenceValidationError("Claude transcript init/result cardinality drift")
    init = init_events[0]
    result = result_events[0]
    expected_tools = {"Read", "Glob", "Grep", "StructuredOutput"}
    if (
        init.get("cwd") != str(REPO_ROOT)
        or init.get("session_id") != session_id
        or set(init.get("tools", [])) != expected_tools
        or init.get("mcp_servers") != []
        or init.get("model") != EXPECTED_CLAUDE_MODEL_ID
        or init.get("permissionMode") != "dontAsk"
        or init.get("claude_code_version") != "2.1.202"
        or result.get("session_id") != session_id
        or result.get("subtype") != "success"
        or result.get("is_error") is not False
        or result.get("terminal_reason") != "completed"
        or result.get("permission_denials") != []
        or not isinstance(result.get("structured_output"), Mapping)
    ):
        raise EvidenceValidationError("Claude transcript process identity drift")

    tool_uses: list[dict[str, Any]] = []
    tool_results: dict[str, dict[str, Any]] = {}
    for event in events:
        if event.get("type") == "assistant":
            message = event.get("message")
            if not isinstance(message, Mapping):
                continue
            if message.get("model") != EXPECTED_CLAUDE_MODEL_ID:
                raise EvidenceValidationError("Claude transcript model changed mid-review")
            content = message.get("content")
            if not isinstance(content, list):
                continue
            for block in content:
                if isinstance(block, Mapping) and block.get("type") == "tool_use":
                    tool_uses.append(dict(block))
        elif event.get("type") == "user":
            message = event.get("message")
            if not isinstance(message, Mapping):
                continue
            content = message.get("content")
            if not isinstance(content, list):
                continue
            for block in content:
                if (
                    isinstance(block, Mapping)
                    and block.get("type") == "tool_result"
                    and isinstance(block.get("tool_use_id"), str)
                ):
                    tool_use_id = str(block["tool_use_id"])
                    if tool_use_id in tool_results:
                        raise EvidenceValidationError(
                            "Claude transcript has duplicate tool results"
                        )
                    tool_results[tool_use_id] = dict(block)

    allowed_tool_events = {"Read", "Glob", "Grep", "StructuredOutput"}
    if not tool_uses or any(
        block.get("name") not in allowed_tool_events for block in tool_uses
    ):
        raise EvidenceValidationError("Claude transcript contains invalid tool events")
    tool_use_ids = [block.get("id") for block in tool_uses]
    if (
        any(not isinstance(tool_use_id, str) for tool_use_id in tool_use_ids)
        or len(tool_use_ids) != len(set(tool_use_ids))
    ):
        raise EvidenceValidationError("Claude transcript has duplicate tool-use ids")
    tool_use_id_set = {str(tool_use_id) for tool_use_id in tool_use_ids}
    missing_results = tool_use_id_set.difference(tool_results)
    if missing_results:
        raise EvidenceValidationError("Claude transcript has unmatched tool events")
    unexpected_results = set(tool_results).difference(tool_use_id_set)
    if unexpected_results:
        raise EvidenceValidationError("Claude transcript has unexpected tool results")

    structured_tool_uses = [
        block for block in tool_uses if block.get("name") == "StructuredOutput"
    ]
    if len(structured_tool_uses) != 1:
        raise EvidenceValidationError(
            "Claude transcript must contain exactly one StructuredOutput tool use"
        )
    structured_tool_use = structured_tool_uses[0]
    structured_tool_input = structured_tool_use.get("input")
    structured_tool_use_id = structured_tool_use.get("id")
    if (
        not isinstance(structured_tool_input, Mapping)
        or not isinstance(structured_tool_use_id, str)
        or dict(structured_tool_input) != dict(result["structured_output"])
    ):
        raise EvidenceValidationError(
            "Claude StructuredOutput tool input/result payload binding drift"
        )
    structured_result_block = tool_results.get(structured_tool_use_id)
    if (
        structured_result_block is None
        or structured_result_block.get("is_error") not in (None, False)
        or structured_result_block.get("content") in (None, "", [])
    ):
        raise EvidenceValidationError(
            "Claude StructuredOutput tool did not return successful nonempty content"
        )

    required_reads = {str(path.resolve()) for path in REQUIRED_REVIEW_READ_PATHS}
    read_paths: set[str] = set()
    tool_counts = {name: 0 for name in sorted(allowed_tool_events)}
    for block in tool_uses:
        name = str(block.get("name"))
        tool_counts[name] += 1
        if name != "Read":
            continue
        inputs = block.get("input")
        tool_use_id = block.get("id")
        if (
            not isinstance(inputs, Mapping)
            or not isinstance(inputs.get("file_path"), str)
            or not isinstance(tool_use_id, str)
        ):
            raise EvidenceValidationError(
                "Claude transcript contains a malformed Read tool use"
            )
        resolved_read_path = str(
            Path(str(inputs["file_path"])).resolve(strict=False)
        )
        if resolved_read_path not in required_reads:
            continue
        result_block = tool_results.get(tool_use_id)
        if (
            set(inputs) != {"file_path"}
            or result_block is None
            or result_block.get("is_error") not in (None, False)
            or result_block.get("content") in (None, "", [])
        ):
            continue
        read_paths.add(resolved_read_path)
    if not required_reads.issubset(read_paths):
        raise EvidenceValidationError(
            "Claude review did not complete an unbounded successful nonempty "
            "Read for every required trust-chain source"
        )
    structured = dict(result["structured_output"])
    if f"Reviewer {reviewer_token}" not in str(structured.get("reviewer_label", "")):
        raise EvidenceValidationError("Claude structured-output reviewer drift")
    return {
        "schema_name": "intersubmod.claude_tool_trace",
        "schema_version": "1.0.0",
        "session_id": session_id,
        "model": EXPECTED_CLAUDE_MODEL_ID,
        "claude_code_version": "2.1.202",
        "tool_counts": tool_counts,
        "read_paths": sorted(read_paths),
        "required_read_paths": sorted(required_reads),
        "event_count": len(events),
        "structured_output": structured,
    }


def verify_bound_signature(
    *,
    reader: Any,
    release_module: Any,
    payload_path: Path,
    signature_path: Path,
    public_key_path: Path,
    private_key_path: Path,
    public_key_sha256: str,
    label: str,
) -> tuple[bytes, dict[str, dict[str, Any]]]:
    openssl_fd, _, openssl_record = reader.open(
        release_module.OPENSSL_PATH, include_mode=True
    )
    payload_fd, payload_bytes, payload_record = reader.open(
        payload_path, include_mode=True
    )
    signature_fd, _, signature_record = reader.open(
        signature_path, include_mode=True
    )
    public_key_fd, _, public_key_record = reader.open(
        public_key_path, include_mode=True
    )
    _, private_key_record = reader.open_metadata(
        private_key_path, include_mode=True
    )
    if openssl_record["sha256"] != release_module.OPENSSL_SHA256:
        raise EvidenceValidationError(f"{label} OpenSSL verifier drift")
    if (
        payload_record["mode"] != "0o444"
        or signature_record["mode"] != "0o444"
        or public_key_record["mode"] != "0o444"
        or public_key_record["sha256"] != public_key_sha256
        or private_key_record["mode"] != "0o0"
    ):
        raise EvidenceValidationError(f"{label} signature/key lifecycle drift")
    if not release_module.verify_ed25519_signature_fds(
        openssl_fd=openssl_fd,
        data_fd=payload_fd,
        public_key_fd=public_key_fd,
        signature_fd=signature_fd,
    ):
        raise EvidenceValidationError(f"{label} detached signature failed")
    return payload_bytes, {
        "manifest": payload_record,
        "signature": signature_record,
        "public_key": public_key_record,
        "private_key": private_key_record,
    }


def validate_signed_test_run(
    *,
    reader: Any,
    release_module: Any,
    expected_head: str,
    expected_support_paths: Mapping[str, Path],
) -> dict[str, Any]:
    manifest_bytes, signature_records = verify_bound_signature(
        reader=reader,
        release_module=release_module,
        payload_path=TEST_MANIFEST_PATH,
        signature_path=TEST_MANIFEST_SIGNATURE_PATH,
        public_key_path=TEST_MANIFEST_PUBLIC_KEY_PATH,
        private_key_path=TEST_MANIFEST_PRIVATE_KEY_PATH,
        public_key_sha256=TEST_MANIFEST_PUBLIC_KEY_SHA256,
        label="Attested canonical test run",
    )
    manifest = parse_json_bytes(manifest_bytes, label="test-run manifest")
    started_at = parse_utc_timestamp(
        manifest.get("started_at_utc"),
        label="test-run started_at_utc",
    )
    finished_at = parse_utc_timestamp(
        manifest.get("finished_at_utc"),
        label="test-run finished_at_utc",
    )
    manifest_command = manifest.get("command")
    if (
        not isinstance(manifest_command, list)
        or not manifest_command
        or not isinstance(manifest_command[0], str)
        or manifest_command[0] != str(TEST_PYTHON_PATH.resolve())
    ):
        raise EvidenceValidationError(
            "Attested pytest argv0 is not the canonical Python runtime"
        )
    expected_command = [
        str(TEST_PYTHON_PATH.resolve()),
        "-I",
        "-X",
        f"pycache_prefix={RAW_TEST_XML_PATH.parent / 'python_cache'}",
        "-m",
        "pytest",
        "-q",
        "-c",
        "/dev/null",
        "--noconftest",
        "--strict-config",
        "--strict-markers",
        "-p",
        "no:cacheprovider",
    ]
    expected_test_paths = sorted(TEST_ROOT.glob("test_*.py"))
    expected_command.extend(str(BOUND_TEST_DIR / path.name) for path in expected_test_paths)
    expected_command.append(f"--junitxml={RAW_TEST_XML_PATH}")
    if (
        set(manifest) != set(TEST_MANIFEST_KEYS)
        or manifest.get("schema_name") != "intersubmod.attested_pytest_run"
        or manifest.get("schema_version") != "1.0.0"
        or manifest.get("task_type") != "B_comprehensive_validation"
        or manifest.get("scope")
        != "all_protected_release_sources_and_canonical_topic_tests"
        or manifest.get("git_head") != expected_head
        or manifest.get("cwd") != str(REPO_ROOT)
        or manifest.get("command") != expected_command
        or manifest.get("environment") != TEST_ENVIRONMENT
        or manifest.get("pass") is not True
        or not canonical_uuid4(manifest.get("run_id"))
        or finished_at < started_at
    ):
        raise EvidenceValidationError("Attested canonical test-run contract drift")

    live_sources = {
        role: reader.open(path, include_mode=True)[2]
        for role, path in release_module.EXPECTED_SOURCE_PATHS.items()
    }
    if (
        manifest.get("protected_sources") != live_sources
        or manifest.get("source_set_sha256") != source_set_digest(live_sources)
    ):
        raise EvidenceValidationError("Test run protected source-set drift")

    test_records = manifest.get("test_sources")
    if not isinstance(test_records, list) or not test_records:
        raise EvidenceValidationError("Test-run source list is empty")
    if len(test_records) != len(expected_test_paths):
        raise EvidenceValidationError("Canonical test source cardinality drift")
    observed_test_records: list[dict[str, Any]] = []
    seen_paths: set[str] = set()
    for index, (record, expected_path) in enumerate(
        zip(test_records, expected_test_paths),
        1,
    ):
        checked = require_artifact_shape(
            record,
            expected_path=expected_path,
            label=f"test source {index}",
        )
        path = Path(checked["path"]).resolve(strict=True)
        if TEST_ROOT.resolve() not in path.parents or path.suffix != ".py":
            raise EvidenceValidationError(f"Test source is outside canonical root: {path}")
        if str(path) in seen_paths:
            raise EvidenceValidationError(f"Duplicate test source: {path}")
        seen_paths.add(str(path))
        live = reader.open(path, include_mode=True)[2]
        if live != checked:
            raise EvidenceValidationError(f"Test source identity drift: {path}")
        observed_test_records.append(live)
    if observed_test_records != test_records:
        raise EvidenceValidationError("Test source ordering or identity drift")

    supporting = manifest.get("supporting_sources")
    if not isinstance(supporting, Mapping) or set(supporting) != set(
        expected_support_paths
    ):
        raise EvidenceValidationError("Test-run supporting source role set drift")
    live_support = {
        role: reader.open(path, include_mode=True)[2]
        for role, path in expected_support_paths.items()
    }
    if dict(supporting) != live_support:
        raise EvidenceValidationError("Test-run supporting source identity drift")

    _, junit_bytes, junit_record = reader.open(
        TEST_XML_PATH, include_mode=True
    )
    _, raw_junit_bytes, raw_junit_record = reader.open(
        RAW_TEST_XML_PATH,
        include_mode=True,
    )
    counts = parse_junit_counts(junit_bytes)
    inventory = parse_junit_inventory(junit_bytes)
    if (
        manifest.get("junit")
        != {
            "artifact": junit_record,
            "counts": counts,
            "testcase_inventory": inventory,
        }
        or raw_junit_bytes != junit_bytes
        or raw_junit_record["sha256"] != junit_record["sha256"]
        or raw_junit_record["size_bytes"] != junit_record["size_bytes"]
        or counts["tests"] != EXPECTED_TESTCASE_COUNT
        or inventory["count"] != EXPECTED_TESTCASE_COUNT
        or inventory["sha256"] != EXPECTED_TESTCASE_INVENTORY_SHA256
        or counts["failures"] != 0
        or counts["errors"] != 0
        or counts["skipped"] != 0
    ):
        raise EvidenceValidationError("Attested canonical JUnit evidence drift")
    captured = manifest.get("captured_output")
    if not isinstance(captured, Mapping) or set(captured) != {
        "stdout",
        "stderr",
        "raw_junit",
    }:
        raise EvidenceValidationError("Test-run captured output schema drift")
    live_captured = {
        "stdout": reader.open(TEST_STDOUT_PATH, include_mode=True)[2],
        "stderr": reader.open(TEST_STDERR_PATH, include_mode=True)[2],
        "raw_junit": raw_junit_record,
    }
    if dict(captured) != live_captured:
        raise EvidenceValidationError("Test-run captured output identity drift")
    runtime = manifest.get("runtime")
    if not isinstance(runtime, Mapping) or set(runtime) != {
        "python",
        "git",
        "entrypoint",
        "pytest_process",
    }:
        raise EvidenceValidationError("Test-run runtime schema drift")
    _, _, live_python = reader.open(TEST_PYTHON_PATH, include_mode=True)
    _, _, live_git = reader.open(release_module.GIT_PATH, include_mode=True)
    if (
        live_python["sha256"] != TEST_PYTHON_SHA256
        or live_python["mode"] != "0o775"
        or live_git["sha256"] != release_module.GIT_SHA256
    ):
        raise EvidenceValidationError("Test-run runtime identity drift")
    validate_entrypoint_record(
        runtime.get("entrypoint"),
        expected_script=live_support["canonical_test_runner"],
        expected_python=live_python,
        expected_environment=TEST_ENVIRONMENT,
        label="canonical pytest runner",
    )
    pytest_process = runtime.get("pytest_process")
    pytest_executable = (
        pytest_process.get("executable")
        if isinstance(pytest_process, Mapping)
        else None
    )
    if (
        not isinstance(pytest_process, Mapping)
        or set(pytest_process) != {"argv0", "executable", "python"}
        or pytest_process.get("argv0") != live_python["path"]
        or pytest_process.get("python") != live_python
        or not isinstance(pytest_executable, str)
        or not pytest_executable.startswith("/proc/self/fd/")
        or not pytest_executable.removeprefix("/proc/self/fd/").isdigit()
    ):
        raise EvidenceValidationError("Bound pytest process identity drift")
    if dict(runtime) != expected_test_runtime_record(
        live_python=live_python,
        live_git=live_git,
        entrypoint=runtime["entrypoint"],
        pytest_process=pytest_process,
    ):
        raise EvidenceValidationError("Test-run runtime record drift")
    if manifest.get("signature_contract") != {
        "algorithm": "Ed25519",
        "public_key": signature_records["public_key"],
        "expected_signature_path": str(TEST_MANIFEST_SIGNATURE_PATH),
        "private_key_path": str(TEST_MANIFEST_PRIVATE_KEY_PATH.resolve()),
        "required_private_key_mode_after_signing": "0o0",
    }:
        raise EvidenceValidationError("Test-run signature contract drift")
    checks = manifest.get("checks")
    expected_checks = {
        "inputs_bound_before_pytest": True,
        "inputs_unchanged_after_pytest": True,
        "pytest_returncode_zero": True,
        "junit_attributes_match_testcase_elements": True,
        "pytest_config_conftest_and_plugin_autoload_disabled": True,
        "pytest_test_sources_executed_through_bound_fds": True,
        "testcase_inventory_matches_pinned_release_contract": True,
        "junit_zero_failures_errors_skips": True,
        "all_outputs_no_clobber_mode_0444": True,
        "manifest_requires_detached_signature": True,
    }
    if checks != expected_checks:
        raise EvidenceValidationError("Attested test-run checks are incomplete")
    return {
        "payload": manifest,
        "records": signature_records,
        "junit": junit_record,
        "counts": counts,
        "testcase_inventory": inventory,
        "source_set_sha256": manifest["source_set_sha256"],
        "supporting_sources": live_support,
    }


def extract_review_contract(request_bytes: bytes) -> dict[str, Any]:
    try:
        text = request_bytes.decode("utf-8")
    except UnicodeDecodeError as error:
        raise EvidenceValidationError("External review request is not UTF-8") from error
    matches = [
        line.removeprefix("REVIEW_CONTRACT_JSON=")
        for line in text.splitlines()
        if line.startswith("REVIEW_CONTRACT_JSON=")
    ]
    if len(matches) != 1:
        raise EvidenceValidationError("External review request contract marker drift")
    return parse_json_bytes(matches[0].encode("utf-8"), label="review request contract")


def validate_clean_review_payload(
    payload: Any,
    *,
    reviewer_token: str,
    source_set_sha256: str,
    expected_head: str,
    assembler_record: Mapping[str, Any],
    junit_record: Mapping[str, Any],
    junit_counts: Mapping[str, int],
) -> dict[str, Any]:
    if not isinstance(payload, Mapping) or set(payload) != set(REVIEW_PAYLOAD_KEYS):
        raise EvidenceValidationError("External review payload schema is not exact")
    result = dict(payload)
    review_scope = result.get("review_scope")
    evidence = result.get("evidence")
    if (
        result["schema_name"] != "intersubmod.external_claude_source_review"
        or result["schema_version"] != "1.0.0"
        or not isinstance(result["reviewer_label"], str)
        or f"Reviewer {reviewer_token}" not in result["reviewer_label"]
        or not canonical_uuid4(result["reviewer_id"])
        or result["model"] != EXPECTED_CLAUDE_MODEL_ID
        or result["verdict"] != "APPROVE"
        or result["findings_closed"] is not True
        or result["f1_status"] != "RESOLVED_VERIFIED"
        or result["f2_status"] != "RESOLVED_VERIFIED"
        or result["blocking_findings"] != []
        or not isinstance(result["nonblocking_findings"], list)
        or not isinstance(evidence, Mapping)
        or len(evidence) < 3
        or not isinstance(review_scope, str)
        or len(review_scope) < 80
        or result["reviewed_source_set_sha256"] != source_set_sha256
        or result["reviewed_git_head"] != expected_head
        or result["reviewed_assembler"] != dict(assembler_record)
        or result["canonical_junit"]
        != {"artifact": dict(junit_record), "counts": dict(junit_counts)}
    ):
        raise EvidenceValidationError("External review is not a clean bound approval")
    return result


def expected_review_process_paths(
    reviewer_token: str,
) -> dict[str, Path]:
    attestation = REVIEW_ATTESTATION_SPECS[reviewer_token]["path"]
    prefix = attestation.with_suffix("")
    return {
        "command_record": Path(str(prefix) + ".command.json"),
        "request": Path(str(prefix) + ".request.md"),
        "transcript": Path(str(prefix) + ".stream.jsonl"),
        "stderr": Path(str(prefix) + ".stderr.txt"),
    }


def validate_review_command_record(
    *,
    command: Mapping[str, Any],
    reviewer_token: str,
    session_id: str,
    raw_request: bytes,
    request_record: Mapping[str, Any],
    claude_record: Mapping[str, Any],
    assembler_record: Mapping[str, Any],
    source_set_sha256: str,
    expected_head: str,
    test_run: Mapping[str, Any],
) -> None:
    if set(command) != set(COMMAND_RECORD_KEYS):
        raise EvidenceValidationError("Claude command-record schema is not exact")
    expected_contract = {
        "reviewer_token": reviewer_token,
        "source_set_sha256": source_set_sha256,
        "git_head": expected_head,
        "assembler": dict(assembler_record),
        "test_run_manifest": test_run["records"]["manifest"],
        "test_run_signature": test_run["records"]["signature"],
        "test_run_public_key": test_run["records"]["public_key"],
        "canonical_junit": {
            "artifact": test_run["junit"],
            "counts": test_run["counts"],
        },
    }
    expected_request = expected_review_request(
        reviewer_token,
        contract=expected_contract,
    )
    if raw_request != expected_request:
        raise EvidenceValidationError("External review request bytes drifted")
    command_line = command.get("command")
    if (
        not isinstance(command_line, list)
        or not command_line
        or not isinstance(command_line[0], str)
        or not command_line[0].startswith("/proc/self/fd/")
        or not command_line[0].removeprefix("/proc/self/fd/").isdigit()
    ):
        raise EvidenceValidationError("Claude command line is malformed")
    try:
        schema_index = command_line.index("--json-schema") + 1
    except (ValueError, IndexError) as error:
        raise EvidenceValidationError("Claude command has no JSON schema") from error
    schema_text = command_line[schema_index]
    if not isinstance(schema_text, str):
        raise EvidenceValidationError("Claude command JSON schema is malformed")
    schema = parse_json_bytes(
        schema_text.encode("utf-8"),
        label="Claude command JSON schema",
    )
    expected_schema = expected_review_schema(
        source_set_sha256=source_set_sha256,
        git_head=expected_head,
        assembler=assembler_record,
        junit=test_run["junit"],
        counts=test_run["counts"],
    )
    canonical_schema_text = json.dumps(
        expected_schema,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    expected_command = [
        command_line[0],
        "--print",
        "--safe-mode",
        "--disable-slash-commands",
        "--no-chrome",
        "--no-session-persistence",
        "--strict-mcp-config",
        "--mcp-config",
        "{\"mcpServers\":{}}",
        "--permission-mode",
        "dontAsk",
        "--tools",
        ALLOWED_REVIEW_TOOLS,
        "--allowedTools",
        ALLOWED_REVIEW_TOOLS,
        "--disallowedTools",
        DISALLOWED_REVIEW_TOOLS,
        "--model",
        EXPECTED_CLAUDE_MODEL_ID,
        "--effort",
        "max",
        "--max-budget-usd",
        "12",
        "--output-format",
        "stream-json",
        "--verbose",
        "--json-schema",
        canonical_schema_text,
        "--session-id",
        session_id,
        "--name",
        f"InterSubMod Task-B external Reviewer {reviewer_token} v21 final git env alignment",
        "--append-system-prompt",
        REVIEW_SYSTEM_PROMPT,
        expected_request.decode("utf-8"),
    ]
    expected_support = test_run["supporting_sources"]
    if (
        schema != expected_schema
        or command_line != expected_command
        or command.get("schema_name")
        != "intersubmod.external_claude_command_record"
        or command.get("schema_version") != "1.0.0"
        or command.get("reviewer_token") != reviewer_token
        or command.get("session_id") != session_id
        or command.get("cwd") != str(REPO_ROOT)
        or command.get("claude_cli") != dict(claude_record)
        or command.get("runner") != expected_support["external_review_runner"]
        or command.get("evidence_validator")
        != expected_support["attested_evidence_validator"]
        or command.get("request") != dict(request_record)
        or command.get("model") != EXPECTED_CLAUDE_MODEL_ID
        or command.get("effort") != "max"
        or command.get("permission_mode") != "dontAsk"
        or command.get("allowed_tools") != ALLOWED_REVIEW_TOOLS
        or command.get("disallowed_tools") != DISALLOWED_REVIEW_TOOLS
        or command.get("environment") != CLEAN_REVIEW_ENVIRONMENT
        or command.get("schema_sha256")
        != hashlib.sha256(encode_json(expected_schema)).hexdigest()
    ):
        raise EvidenceValidationError("Claude command/process contract drift")
    validate_entrypoint_record(
        command.get("entrypoint"),
        expected_script=expected_support["external_review_runner"],
        expected_python=test_run["payload"]["runtime"]["python"],
        expected_environment=CLEAN_REVIEW_ENVIRONMENT,
        label="external Claude review runner",
    )


def validate_signed_review_attestation(
    *,
    reader: Any,
    release_module: Any,
    reviewer_token: str,
    expected_head: str,
    source_set_sha256: str,
    assembler_record: Mapping[str, Any],
    test_run: Mapping[str, Any],
) -> dict[str, Any]:
    if reviewer_token not in REVIEW_ATTESTATION_SPECS:
        raise EvidenceValidationError(f"Unknown reviewer token: {reviewer_token}")
    spec = REVIEW_ATTESTATION_SPECS[reviewer_token]
    attestation_bytes, signature_records = verify_bound_signature(
        reader=reader,
        release_module=release_module,
        payload_path=spec["path"],
        signature_path=spec["signature"],
        public_key_path=spec["public_key"],
        private_key_path=spec["private_key"],
        public_key_sha256=spec["public_key_sha256"],
        label=f"Reviewer {reviewer_token} process attestation",
    )
    attestation = parse_json_bytes(
        attestation_bytes,
        label=f"Reviewer {reviewer_token} process attestation",
    )
    if (
        set(attestation) != set(ATTESTATION_KEYS)
        or
        attestation.get("schema_name")
        != "intersubmod.external_claude_process_attestation"
        or attestation.get("schema_version") != "1.0.0"
        or attestation.get("reviewer_token") != reviewer_token
        or not canonical_uuid4(attestation.get("attestation_id"))
        or attestation.get("pass") is not True
    ):
        raise EvidenceValidationError(
            f"Reviewer {reviewer_token} process-attestation contract drift"
        )
    process = attestation.get("process")
    if not isinstance(process, Mapping) or set(process) != set(PROCESS_KEYS):
        raise EvidenceValidationError("External review process record is malformed")
    raw_records: dict[str, dict[str, Any]] = {}
    raw_bytes: dict[str, bytes] = {}
    process_paths = expected_review_process_paths(reviewer_token)
    for field in ("command_record", "request", "transcript", "stderr"):
        record = require_artifact_shape(
            process.get(field),
            expected_path=process_paths[field],
            label=f"review process {field}",
        )
        path = Path(record["path"]).resolve(strict=True)
        _, data, live = reader.open(path, include_mode=True)
        if live != record:
            raise EvidenceValidationError(f"External review {field} identity drift")
        raw_records[field] = live
        raw_bytes[field] = data
    _, _, claude_record = reader.open(CLAUDE_CLI_PATH, include_mode=True)
    if (
        process.get("claude_cli") != claude_record
        or claude_record["sha256"] != CLAUDE_CLI_SHA256
        or claude_record["mode"] != "0o755"
    ):
        raise EvidenceValidationError("Pinned Claude CLI identity drift")
    command = parse_json_bytes(raw_bytes["command_record"], label="Claude command record")
    request_contract = extract_review_contract(raw_bytes["request"])
    session_id = process.get("session_id")
    if not canonical_uuid4(session_id):
        raise EvidenceValidationError("External Claude session id drift")
    tool_trace = parse_claude_stream_transcript(
        raw_bytes["transcript"],
        reviewer_token=reviewer_token,
        session_id=session_id,
    )
    if (
        command.get("session_id") != session_id
        or process.get("returncode") != 0
        or process.get("model") != EXPECTED_CLAUDE_MODEL_ID
        or process.get("effort") != "max"
        or process.get("permission_mode") != "dontAsk"
        or process.get("allowed_tools") != ALLOWED_REVIEW_TOOLS
        or process.get("disallowed_tools") != DISALLOWED_REVIEW_TOOLS
        or process.get("environment") != CLEAN_REVIEW_ENVIRONMENT
        or process.get("tool_trace") != tool_trace
    ):
        raise EvidenceValidationError("External Claude process/session linkage drift")
    validate_entrypoint_record(
        process.get("entrypoint"),
        expected_script=test_run["supporting_sources"]["external_review_runner"],
        expected_python=test_run["payload"]["runtime"]["python"],
        expected_environment=CLEAN_REVIEW_ENVIRONMENT,
        label="external Claude review runner",
    )
    expected_contract = {
        "reviewer_token": reviewer_token,
        "source_set_sha256": source_set_sha256,
        "git_head": expected_head,
        "assembler": dict(assembler_record),
        "test_run_manifest": test_run["records"]["manifest"],
        "test_run_signature": test_run["records"]["signature"],
        "test_run_public_key": test_run["records"]["public_key"],
        "canonical_junit": {
            "artifact": test_run["junit"],
            "counts": test_run["counts"],
        },
    }
    if request_contract != expected_contract:
        raise EvidenceValidationError("External review request binding drift")
    validate_review_command_record(
        command=command,
        reviewer_token=reviewer_token,
        session_id=session_id,
        raw_request=raw_bytes["request"],
        request_record=raw_records["request"],
        claude_record=claude_record,
        assembler_record=assembler_record,
        source_set_sha256=source_set_sha256,
        expected_head=expected_head,
        test_run=test_run,
    )
    review_payload = validate_clean_review_payload(
        attestation.get("review_payload"),
        reviewer_token=reviewer_token,
        source_set_sha256=source_set_sha256,
        expected_head=expected_head,
        assembler_record=assembler_record,
        junit_record=test_run["junit"],
        junit_counts=test_run["counts"],
    )
    if (
        tool_trace["structured_output"] != review_payload
        or attestation.get("review_payload_sha256")
        != hashlib.sha256(encode_json(review_payload)).hexdigest()
        or attestation.get("reviewed_test_run")
        != {
            "manifest": test_run["records"]["manifest"],
            "signature": test_run["records"]["signature"],
            "public_key": test_run["records"]["public_key"],
            "junit": test_run["junit"],
            "counts": test_run["counts"],
            "source_set_sha256": source_set_sha256,
            "git_head": expected_head,
        }
    ):
        raise EvidenceValidationError("External review attestation evidence drift")
    if attestation.get("signature_contract") != {
        "algorithm": "Ed25519",
        "public_key": signature_records["public_key"],
        "expected_signature_path": str(spec["signature"]),
        "private_key_path": str(spec["private_key"].resolve()),
        "required_private_key_mode_after_signing": "0o0",
    }:
        raise EvidenceValidationError("External review signature contract drift")
    checks = attestation.get("checks")
    expected_checks = {
        "pinned_claude_cli_identity_verified": True,
        "clean_environment_and_read_only_tools_enforced": True,
        "command_request_session_transcript_linked": True,
        "required_read_tool_events_verified": True,
        "bound_runner_entrypoint_verified": True,
        "signed_test_run_reverified": True,
        "review_payload_schema_and_clean_approval_verified": True,
        "attestation_requires_detached_signature": True,
    }
    if checks != expected_checks:
        raise EvidenceValidationError("External review attestation checks are incomplete")
    return {
        "review_payload": review_payload,
        "attestation": signature_records["manifest"],
        "signature": signature_records["signature"],
        "public_key": signature_records["public_key"],
        "private_key": signature_records["private_key"],
        "process_records": raw_records,
        "session_id": session_id,
    }


__all__ = [
    "AUTHORITY_ASSEMBLER_PATH",
    "CLAUDE_CLI_PATH",
    "CLAUDE_CLI_SHA256",
    "EVIDENCE_VALIDATOR_PATH",
    "EvidenceValidationError",
    "REVIEW_ATTESTATION_SPECS",
    "REVIEW_NORMALIZER_PATH",
    "REVIEW_RUNNER_PATH",
    "SIGNER_PATH",
    "TEST_MANIFEST_PATH",
    "TEST_MANIFEST_PRIVATE_KEY_PATH",
    "TEST_MANIFEST_PUBLIC_KEY_PATH",
    "TEST_MANIFEST_PUBLIC_KEY_SHA256",
    "TEST_MANIFEST_SIGNATURE_PATH",
    "TEST_ROOT",
    "TEST_RUNNER_PATH",
    "TEST_XML_PATH",
    "encode_json",
    "parse_junit_counts",
    "source_set_digest",
    "validate_clean_review_payload",
    "validate_signed_review_attestation",
    "validate_signed_test_run",
]
