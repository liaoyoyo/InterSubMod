from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = TOPIC_ROOT / "audit_notes" / "attested_release_evidence_v5.py"
RELEASE_PATH = TOPIC_ROOT / "scripts" / "release_source_authority.py"
EXTERNAL_RUNNER_PATH = (
    TOPIC_ROOT / "audit_notes" / "run_external_claude_review_v22_attested.py"
)


def load(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def module() -> ModuleType:
    return load(MODULE_PATH, "attested_release_evidence_v5_test")


def record(path: str, token: str) -> dict[str, object]:
    return {
        "path": path,
        "size_bytes": 10,
        "sha256": token * 64,
        "mode": "0o444",
    }


def test_strict_json_rejects_duplicate_and_nonfinite(module: ModuleType) -> None:
    with pytest.raises(module.EvidenceValidationError, match="malformed"):
        module.parse_json_bytes(b'{"a":1,"a":2}', label="duplicate")
    with pytest.raises(module.EvidenceValidationError, match="malformed"):
        module.parse_json_bytes(b'{"a":NaN}', label="nonfinite")


def test_junit_requires_attribute_element_agreement(module: ModuleType) -> None:
    valid = (
        b'<testsuites><testsuite tests="1" failures="0" errors="0" skipped="0">'
        b'<testcase name="ok"/></testsuite></testsuites>'
    )
    assert module.parse_junit_counts(valid) == {
        "tests": 1,
        "failures": 0,
        "errors": 0,
        "skipped": 0,
    }
    with pytest.raises(module.EvidenceValidationError, match="differ"):
        module.parse_junit_counts(
            valid.replace(b'tests="1"', b'tests="2"')
        )
    inventory = module.parse_junit_inventory(valid)
    assert inventory["count"] == 1
    assert len(inventory["sha256"]) == 64
    python = record("/tmp/python", "1")
    git = record("/tmp/git", "2")
    entrypoint = {"argv0": "/proc/self/fd/3"}
    pytest_process = {
        "argv0": "/tmp/python",
        "executable": "/proc/self/fd/4",
        "python": python,
    }
    assert module.expected_test_runtime_record(
        live_python=python,
        live_git=git,
        entrypoint=entrypoint,
        pytest_process=pytest_process,
    ) == {
        "python": python,
        "git": git,
        "entrypoint": entrypoint,
        "pytest_process": pytest_process,
    }


def test_review_request_has_one_exact_contract_marker(module: ModuleType) -> None:
    contract = {
        "reviewer_token": "A",
        "source_set_sha256": "d" * 64,
    }
    request = module.expected_review_request("A", contract=contract)
    assert request.count(b"REVIEW_CONTRACT_JSON=") == 1
    assert module.extract_review_contract(request) == contract
    runner_source = EXTERNAL_RUNNER_PATH.read_text(encoding="utf-8")
    assert "evidence_module.expected_review_request(" in runner_source
    assert "def build_request(" not in runner_source
    assert "evidence_module.expected_review_schema(" in runner_source
    assert "def review_schema(" not in runner_source


def valid_command_inputs(module: ModuleType) -> tuple[
    dict[str, object],
    bytes,
    dict[str, object],
    dict[str, object],
    dict[str, object],
    dict[str, object],
]:
    assembler = record("/tmp/assembler.py", "a")
    junit = record("/tmp/junit.xml", "b")
    request_record = record("/tmp/request.md", "c")
    claude = {
        **record(str(module.CLAUDE_CLI_PATH), "d"),
        "mode": "0o755",
    }
    test_run = {
        "records": {
            "manifest": record("/tmp/test-manifest.json", "e"),
            "signature": record("/tmp/test-manifest.sig", "f"),
            "public_key": record("/tmp/test-public.pem", "1"),
        },
        "junit": junit,
        "counts": {
            "tests": 5,
            "failures": 0,
            "errors": 0,
            "skipped": 0,
        },
        "supporting_sources": {
            "external_review_runner": record("/tmp/review-runner.py", "2"),
            "attested_evidence_validator": record("/tmp/evidence.py", "3"),
        },
        "payload": {
            "runtime": {
                "python": {
                    **record(str(module.TEST_PYTHON_PATH), "6"),
                    "mode": "0o775",
                }
            }
        },
    }
    contract = {
        "reviewer_token": "A",
        "source_set_sha256": "4" * 64,
        "git_head": "5" * 40,
        "assembler": assembler,
        "test_run_manifest": test_run["records"]["manifest"],
        "test_run_signature": test_run["records"]["signature"],
        "test_run_public_key": test_run["records"]["public_key"],
        "canonical_junit": {
            "artifact": junit,
            "counts": test_run["counts"],
        },
    }
    request = module.expected_review_request("A", contract=contract)
    schema = module.expected_review_schema(
        source_set_sha256="4" * 64,
        git_head="5" * 40,
        assembler=assembler,
        junit=junit,
        counts=test_run["counts"],
    )
    schema_text = json.dumps(schema, separators=(",", ":"), sort_keys=True)
    command_line = [
        "/proc/self/fd/31",
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
        module.ALLOWED_REVIEW_TOOLS,
        "--allowedTools",
        module.ALLOWED_REVIEW_TOOLS,
        "--disallowedTools",
        module.DISALLOWED_REVIEW_TOOLS,
        "--model",
        module.EXPECTED_CLAUDE_MODEL_ID,
        "--effort",
        "max",
        "--max-budget-usd",
        "12",
        "--output-format",
        "stream-json",
        "--verbose",
        "--json-schema",
        schema_text,
        "--session-id",
        "45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
        "--name",
        "InterSubMod Task-B external Reviewer A v21 final git env alignment",
        "--append-system-prompt",
        module.REVIEW_SYSTEM_PROMPT,
        request.decode("utf-8"),
    ]
    command = {
        "schema_name": "intersubmod.external_claude_command_record",
        "schema_version": "1.0.0",
        "reviewer_token": "A",
        "session_id": "45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
        "cwd": str(module.REPO_ROOT),
        "entrypoint": {
            "schema_name": "intersubmod.bound_python_entrypoint",
            "schema_version": "1.0.0",
            "argv0": "/proc/self/fd/30",
            "script": test_run["supporting_sources"]["external_review_runner"],
            "python": test_run["payload"]["runtime"]["python"],
            "isolated_mode": True,
            "no_user_site": True,
            "environment": module.CLEAN_REVIEW_ENVIRONMENT,
        },
        "claude_cli": claude,
        "runner": test_run["supporting_sources"]["external_review_runner"],
        "evidence_validator": test_run["supporting_sources"][
            "attested_evidence_validator"
        ],
        "request": request_record,
        "model": module.EXPECTED_CLAUDE_MODEL_ID,
        "effort": "max",
        "permission_mode": "dontAsk",
        "allowed_tools": module.ALLOWED_REVIEW_TOOLS,
        "disallowed_tools": module.DISALLOWED_REVIEW_TOOLS,
        "environment": module.CLEAN_REVIEW_ENVIRONMENT,
        "schema_sha256": hashlib.sha256(module.encode_json(schema)).hexdigest(),
        "command": command_line,
    }
    return command, request, request_record, claude, assembler, test_run


def test_command_consumer_accepts_only_fd_readonly_contract(
    module: ModuleType,
) -> None:
    command, request, request_record, claude, assembler, test_run = (
        valid_command_inputs(module)
    )
    module.validate_review_command_record(
        command=command,
        reviewer_token="A",
        session_id="45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
        raw_request=request,
        request_record=request_record,
        claude_record=claude,
        assembler_record=assembler,
        source_set_sha256="4" * 64,
        expected_head="5" * 40,
        test_run=test_run,
    )
    schema_index = command["command"].index("--json-schema") + 1
    transport_schema = json.loads(command["command"][schema_index])
    assert set(transport_schema) == {
        "type",
        "additionalProperties",
        "required",
        "properties",
    }
    assert transport_schema["required"] == list(
        module.REVIEW_PAYLOAD_FIELD_ORDER
    )
    assert transport_schema["properties"]["reviewed_assembler"] == {
        "type": "object"
    }
    assert transport_schema["properties"]["blocking_findings"] == {
        "type": "array"
    }

    valid_review = {
        "schema_name": "intersubmod.external_claude_source_review",
        "schema_version": "1.0.0",
        "reviewer_label": "External Reviewer A v21",
        "reviewer_id": "45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
        "model": module.EXPECTED_CLAUDE_MODEL_ID,
        "verdict": "APPROVE",
        "findings_closed": True,
        "f1_status": "RESOLVED_VERIFIED",
        "f2_status": "RESOLVED_VERIFIED",
        "reviewed_source_set_sha256": "4" * 64,
        "reviewed_git_head": "5" * 40,
        "reviewed_assembler": assembler,
        "canonical_junit": {
            "artifact": test_run["junit"],
            "counts": test_run["counts"],
        },
        "review_scope": "x" * 80,
        "blocking_findings": [],
        "nonblocking_findings": [],
        "evidence": {"a": 1, "b": 2, "c": 3},
    }
    assert module.validate_clean_review_payload(
        valid_review,
        reviewer_token="A",
        source_set_sha256="4" * 64,
        expected_head="5" * 40,
        assembler_record=assembler,
        junit_record=test_run["junit"],
        junit_counts=test_run["counts"],
    ) == valid_review
    for field, bad_value in (
        ("review_scope", "too short"),
        ("evidence", {}),
    ):
        invalid_review = copy.deepcopy(valid_review)
        invalid_review[field] = bad_value
        with pytest.raises(
            module.EvidenceValidationError,
            match="not a clean bound approval",
        ):
            module.validate_clean_review_payload(
                invalid_review,
                reviewer_token="A",
                source_set_sha256="4" * 64,
                expected_head="5" * 40,
                assembler_record=assembler,
                junit_record=test_run["junit"],
                junit_counts=test_run["counts"],
            )

    tampered = copy.deepcopy(command)
    tampered["command"][0] = str(module.CLAUDE_CLI_PATH)
    with pytest.raises(module.EvidenceValidationError, match="command line"):
        module.validate_review_command_record(
            command=tampered,
            reviewer_token="A",
            session_id="45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
            raw_request=request,
            request_record=request_record,
            claude_record=claude,
            assembler_record=assembler,
            source_set_sha256="4" * 64,
            expected_head="5" * 40,
            test_run=test_run,
        )


def test_external_runner_source_exposes_no_bash_or_config_tree(
    module: ModuleType,
) -> None:
    source = module.REVIEW_RUNNER_PATH.read_text(encoding="utf-8")
    assert 'ALLOWED_TOOLS = "Read,Glob,Grep,StructuredOutput"' in source
    runner = load(EXTERNAL_RUNNER_PATH, "structured_output_tool_parity_runner")
    assert runner.ALLOWED_TOOLS == module.ALLOWED_REVIEW_TOOLS
    assert set(runner.ALLOWED_TOOLS.split(",")) == {
        "Read",
        "Glob",
        "Grep",
        "StructuredOutput",
    }
    assert "system_prompt = evidence_module.REVIEW_SYSTEM_PROMPT" in source
    assert '"--add-dir"' not in source
    assert "paths[\"request\"].read_text" not in source
    assert "pass_fds=(claude_fd,)" in source


def test_external_runner_binds_modules_and_rejects_path_replacement(
    tmp_path: Path,
) -> None:
    runner = load(EXTERNAL_RUNNER_PATH, "external_review_v20_binding_test")
    source_text = EXTERNAL_RUNNER_PATH.read_text(encoding="utf-8")
    assert "execute_module_from_bytes(" in source_text
    assert "release_bytes" in source_text
    assert "evidence_bytes" in source_text
    assert "importlib" not in source_text
    assert "spec_from_file_location" not in source_text

    source = tmp_path / "review_input.json"
    displaced = tmp_path / "review_input.displaced.json"
    source.write_bytes(b'{"trusted":true}\n')
    source.chmod(0o444)
    with runner.BoundArtifactReader() as reader:
        reader.open(source, include_mode=True)
        source.rename(displaced)
        source.write_bytes(b'{"trusted":false}\n')
        source.chmod(0o444)
        with pytest.raises(runner.ReviewRunnerError, match="binding changed"):
            reader.require_paths_still_bound()


def synthetic_claude_transcript(
    module: ModuleType,
    *,
    include_all_required_reads: bool,
    limited_read: bool = False,
    error_read: bool = False,
    duplicate_tool_use_ids: bool = False,
    structured_output_count: int = 1,
    structured_input_mismatch: bool = False,
    structured_result_error: bool = False,
    omit_structured_result: bool = False,
    unexpected_tool_result: bool = False,
    optional_bounded_read: bool = False,
    optional_error_read: bool = False,
) -> bytes:
    session_id = "45df7edf-681f-4a4c-b69f-d0ea2c3ef527"
    required_paths = list(module.REQUIRED_REVIEW_READ_PATHS)
    if not include_all_required_reads:
        required_paths.pop()
    events: list[dict[str, object]] = [
        {
            "type": "system",
            "subtype": "init",
            "cwd": str(module.REPO_ROOT),
            "session_id": session_id,
            "tools": ["Read", "Glob", "Grep", "StructuredOutput"],
            "mcp_servers": [],
            "model": module.EXPECTED_CLAUDE_MODEL_ID,
            "permissionMode": "dontAsk",
            "claude_code_version": "2.1.202",
        }
    ]
    for index, path in enumerate(required_paths, 1):
        tool_id = "tool-1" if duplicate_tool_use_ids else f"tool-{index}"
        read_input: dict[str, object] = {"file_path": str(path)}
        if limited_read and index == 1:
            read_input["limit"] = 1
        events.append(
            {
                "type": "assistant",
                "message": {
                    "model": module.EXPECTED_CLAUDE_MODEL_ID,
                    "content": [
                        {
                            "type": "tool_use",
                            "id": tool_id,
                            "name": "Read",
                            "input": read_input,
                        }
                    ],
                },
            }
        )
        if not duplicate_tool_use_ids or index == 1:
            events.append(
                {
                    "type": "user",
                    "message": {
                        "content": [
                            {
                                "type": "tool_result",
                                "tool_use_id": tool_id,
                                "content": "reviewed",
                                "is_error": error_read and index == 1,
                            }
                        ]
                    },
                }
            )
    optional_index = len(required_paths) + 1
    if optional_bounded_read:
        events.extend(
            [
                {
                    "type": "assistant",
                    "message": {
                        "model": module.EXPECTED_CLAUDE_MODEL_ID,
                        "content": [
                            {
                                "type": "tool_use",
                                "id": f"tool-{optional_index}",
                                "name": "Read",
                                "input": {
                                    "file_path": str(module.REPO_ROOT / "README.md"),
                                    "offset": 1,
                                    "limit": 1,
                                },
                            }
                        ],
                    },
                },
                {
                    "type": "user",
                    "message": {
                        "content": [
                            {
                                "type": "tool_result",
                                "tool_use_id": f"tool-{optional_index}",
                                "content": "optional excerpt",
                                "is_error": False,
                            }
                        ]
                    },
                },
            ]
        )
        optional_index += 1
    if optional_error_read:
        events.extend(
            [
                {
                    "type": "assistant",
                    "message": {
                        "model": module.EXPECTED_CLAUDE_MODEL_ID,
                        "content": [
                            {
                                "type": "tool_use",
                                "id": f"tool-{optional_index}",
                                "name": "Read",
                                "input": {
                                    "file_path": str(
                                        module.REPO_ROOT / "expected_absent.json"
                                    )
                                },
                            }
                        ],
                    },
                },
                {
                    "type": "user",
                    "message": {
                        "content": [
                            {
                                "type": "tool_result",
                                "tool_use_id": f"tool-{optional_index}",
                                "content": "File does not exist",
                                "is_error": True,
                            }
                        ]
                    },
                },
            ]
        )
    structured = {"reviewer_label": "External Claude Reviewer A"}
    for index in range(structured_output_count):
        structured_tool_input = dict(structured)
        if structured_input_mismatch and index == 0:
            structured_tool_input["unexpected"] = True
        structured_tool_id = f"structured-output-{index + 1}"
        events.append(
            {
                "type": "assistant",
                "message": {
                    "model": module.EXPECTED_CLAUDE_MODEL_ID,
                    "content": [
                        {
                            "type": "tool_use",
                            "id": structured_tool_id,
                            "name": "StructuredOutput",
                            "input": structured_tool_input,
                        }
                    ],
                },
            }
        )
        if not omit_structured_result:
            events.append(
                {
                    "type": "user",
                    "message": {
                        "content": [
                            {
                                "type": "tool_result",
                                "tool_use_id": structured_tool_id,
                                "content": "Structured output provided successfully",
                                "is_error": structured_result_error and index == 0,
                            }
                        ]
                    },
                }
            )
    if unexpected_tool_result:
        events.append(
            {
                "type": "user",
                "message": {
                    "content": [
                        {
                            "type": "tool_result",
                            "tool_use_id": "unexpected-tool-use",
                            "content": "unexpected",
                            "is_error": False,
                        }
                    ]
                },
            }
        )
    events.append(
        {
            "type": "result",
            "subtype": "success",
            "is_error": False,
            "session_id": session_id,
            "terminal_reason": "completed",
            "permission_denials": [],
            "structured_output": structured,
        }
    )
    return b"".join(
        json.dumps(event, separators=(",", ":"), sort_keys=True).encode("utf-8")
        + b"\n"
        for event in events
    )


def test_claude_stream_requires_actual_read_events(module: ModuleType) -> None:
    session_id = "45df7edf-681f-4a4c-b69f-d0ea2c3ef527"
    trace = module.parse_claude_stream_transcript(
        synthetic_claude_transcript(
            module,
            include_all_required_reads=True,
        ),
        reviewer_token="A",
        session_id=session_id,
    )
    assert trace["tool_counts"]["Read"] == len(module.REQUIRED_REVIEW_READ_PATHS)
    assert trace["tool_counts"]["StructuredOutput"] == 1
    assert trace["model"] == module.EXPECTED_CLAUDE_MODEL_ID

    exploratory_trace = module.parse_claude_stream_transcript(
        synthetic_claude_transcript(
            module,
            include_all_required_reads=True,
            optional_bounded_read=True,
            optional_error_read=True,
        ),
        reviewer_token="A",
        session_id=session_id,
    )
    assert exploratory_trace["tool_counts"]["Read"] == (
        len(module.REQUIRED_REVIEW_READ_PATHS) + 2
    )
    assert exploratory_trace["read_paths"] == sorted(
        str(path.resolve()) for path in module.REQUIRED_REVIEW_READ_PATHS
    )

    with pytest.raises(module.EvidenceValidationError, match="did not complete"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=False,
            ),
            reviewer_token="A",
            session_id=session_id,
        )
    with pytest.raises(module.EvidenceValidationError, match="did not complete"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=True,
                limited_read=True,
            ),
            reviewer_token="A",
            session_id=session_id,
        )
    with pytest.raises(module.EvidenceValidationError, match="did not complete"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=True,
                error_read=True,
            ),
            reviewer_token="A",
            session_id=session_id,
        )
    with pytest.raises(module.EvidenceValidationError, match="unmatched tool events"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=True,
                omit_structured_result=True,
            ),
            reviewer_token="A",
            session_id=session_id,
        )
    with pytest.raises(module.EvidenceValidationError, match="unexpected tool results"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=True,
                unexpected_tool_result=True,
            ),
            reviewer_token="A",
            session_id=session_id,
        )
    with pytest.raises(module.EvidenceValidationError, match="duplicate tool-use"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=True,
                duplicate_tool_use_ids=True,
            ),
            reviewer_token="A",
            session_id=session_id,
        )
    with pytest.raises(module.EvidenceValidationError, match="exactly one"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=True,
                structured_output_count=0,
            ),
            reviewer_token="A",
            session_id=session_id,
        )
    with pytest.raises(module.EvidenceValidationError, match="exactly one"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=True,
                structured_output_count=2,
            ),
            reviewer_token="A",
            session_id=session_id,
        )
    with pytest.raises(module.EvidenceValidationError, match="payload binding drift"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=True,
                structured_input_mismatch=True,
            ),
            reviewer_token="A",
            session_id=session_id,
        )
    with pytest.raises(module.EvidenceValidationError, match="successful nonempty"):
        module.parse_claude_stream_transcript(
            synthetic_claude_transcript(
                module,
                include_all_required_reads=True,
                structured_result_error=True,
            ),
            reviewer_token="A",
            session_id=session_id,
        )


def test_bound_signature_requires_retired_private_key(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    release = load(RELEASE_PATH, "release_for_attested_signature_test")
    private_key = tmp_path / "private.pem"
    public_key = tmp_path / "public.pem"
    payload = tmp_path / "payload.json"
    signature = tmp_path / "payload.json.ed25519.sig"
    payload.write_bytes(b'{"signed":true}\n')
    subprocess.run(
        [
            str(release.OPENSSL_PATH),
            "genpkey",
            "-algorithm",
            "ED25519",
            "-out",
            str(private_key),
        ],
        check=True,
    )
    subprocess.run(
        [
            str(release.OPENSSL_PATH),
            "pkey",
            "-in",
            str(private_key),
            "-pubout",
            "-out",
            str(public_key),
        ],
        check=True,
    )
    subprocess.run(
        [
            str(release.OPENSSL_PATH),
            "pkeyutl",
            "-sign",
            "-rawin",
            "-inkey",
            str(private_key),
            "-in",
            str(payload),
            "-out",
            str(signature),
        ],
        check=True,
    )
    for path in (payload, public_key, signature):
        path.chmod(0o444)
    private_key.chmod(0o000)
    public_sha = hashlib.sha256(public_key.read_bytes()).hexdigest()
    try:
        with release.BoundArtifactReader() as reader:
            data, records = module.verify_bound_signature(
                reader=reader,
                release_module=release,
                payload_path=payload,
                signature_path=signature,
                public_key_path=public_key,
                private_key_path=private_key,
                public_key_sha256=public_sha,
                label="unit-test attestation",
            )
            assert data == b'{"signed":true}\n'
            assert records["private_key"]["mode"] == "0o0"
    finally:
        private_key.chmod(0o600)
