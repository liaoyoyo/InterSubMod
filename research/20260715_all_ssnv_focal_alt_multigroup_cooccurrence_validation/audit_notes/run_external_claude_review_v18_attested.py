#!/usr/bin/env python3
"""Run one read-only Claude review and create a signable process attestation."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
from types import ModuleType
from typing import Any, Mapping
import uuid


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
HERE = TOPIC_ROOT / "audit_notes"
RELEASE_MODULE_PATH = TOPIC_ROOT / "scripts" / "release_source_authority.py"
EVIDENCE_MODULE_PATH = HERE / "attested_release_evidence_v1.py"
ASSEMBLER_PATH = HERE / "assemble_release_source_authority_v11.py"
ALLOWED_TOOLS = "Read,Glob,Grep,StructuredOutput"
DISALLOWED_TOOLS = ",".join(
    (
        "Bash",
        "Edit",
        "Write",
        "NotebookEdit",
        "WebFetch",
        "WebSearch",
    )
)


class ReviewRunnerError(RuntimeError):
    """Raised when an external-review process cannot be attested."""


_PROCESS_LIFETIME_READERS: list["BoundArtifactReader"] = []


class BoundArtifactReader:
    """Bind review inputs before executing any trust-bearing module."""

    def __init__(self) -> None:
        self._opened: list[tuple[Path, int, os.stat_result]] = []
        self._retained = False

    def __enter__(self) -> "BoundArtifactReader":
        return self

    def __exit__(self, *_: Any) -> None:
        if not self._retained:
            self.close()

    @staticmethod
    def _same(left: os.stat_result, right: os.stat_result) -> bool:
        return (
            left.st_dev,
            left.st_ino,
            left.st_mode,
            left.st_nlink,
            left.st_size,
            left.st_mtime_ns,
            left.st_ctime_ns,
        ) == (
            right.st_dev,
            right.st_ino,
            right.st_mode,
            right.st_nlink,
            right.st_size,
            right.st_mtime_ns,
            right.st_ctime_ns,
        )

    def open(
        self,
        path: Path,
        *,
        include_mode: bool = False,
    ) -> tuple[int, bytes, dict[str, Any]]:
        resolved = path.expanduser().resolve(strict=True)
        flags = os.O_RDONLY | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(resolved, flags)
        first = os.fstat(descriptor)
        if not stat.S_ISREG(first.st_mode):
            os.close(descriptor)
            raise ReviewRunnerError(f"Review input is not regular: {resolved}")
        data = b"".join(
            os.pread(
                descriptor,
                min(8 * 1024 * 1024, first.st_size - offset),
                offset,
            )
            for offset in range(0, first.st_size, 8 * 1024 * 1024)
        )
        second = os.fstat(descriptor)
        live = resolved.stat()
        if (
            len(data) != first.st_size
            or not self._same(first, second)
            or not self._same(second, live)
        ):
            os.close(descriptor)
            raise ReviewRunnerError(f"Review input changed while binding: {resolved}")
        record: dict[str, Any] = {
            "path": str(resolved),
            "size_bytes": second.st_size,
            "sha256": hashlib.sha256(data).hexdigest(),
        }
        if include_mode:
            record["mode"] = oct(stat.S_IMODE(second.st_mode))
        self._opened.append((resolved, descriptor, second))
        return descriptor, data, record

    def open_metadata(
        self,
        path: Path,
        *,
        include_mode: bool = False,
    ) -> tuple[int, dict[str, Any]]:
        if not hasattr(os, "O_PATH"):
            raise ReviewRunnerError("O_PATH is required for signer-key binding")
        resolved = path.expanduser().resolve(strict=True)
        flags = os.O_PATH | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(resolved, flags)
        opened = os.fstat(descriptor)
        live = resolved.stat()
        if (
            not stat.S_ISREG(opened.st_mode)
            or not self._same(opened, live)
        ):
            os.close(descriptor)
            raise ReviewRunnerError(f"Signer-key metadata binding drift: {resolved}")
        record: dict[str, Any] = {"path": str(resolved)}
        if include_mode:
            record["mode"] = oct(stat.S_IMODE(opened.st_mode))
        self._opened.append((resolved, descriptor, opened))
        return descriptor, record

    def require_paths_still_bound(self) -> None:
        for path, descriptor, opened in self._opened:
            if (
                not self._same(os.fstat(descriptor), opened)
                or not self._same(path.stat(), opened)
            ):
                raise ReviewRunnerError(f"Review input binding changed: {path}")

    def close(self) -> None:
        while self._opened:
            _, descriptor, _ = self._opened.pop()
            os.close(descriptor)

    def retain_until_process_exit(self) -> None:
        if self._retained:
            raise ReviewRunnerError("Review input reader was already retained")
        self.require_paths_still_bound()
        _PROCESS_LIFETIME_READERS.append(self)
        self._retained = True


def execute_module_from_bytes(
    name: str,
    path: Path,
    data: bytes,
) -> ModuleType:
    module = ModuleType(name)
    module.__file__ = str(path)
    module.__package__ = ""
    exec(compile(data, str(path), "exec"), module.__dict__)
    return module


def artifact_from_fd(path: Path, descriptor: int) -> dict[str, Any]:
    opened = os.fstat(descriptor)
    live = path.resolve(strict=True).stat()
    if (
        not stat.S_ISREG(opened.st_mode)
        or not BoundArtifactReader._same(opened, live)
    ):
        raise ReviewRunnerError(f"Bound entrypoint identity drift: {path}")
    data = b"".join(
        os.pread(
            descriptor,
            min(8 * 1024 * 1024, opened.st_size - offset),
            offset,
        )
        for offset in range(0, opened.st_size, 8 * 1024 * 1024)
    )
    if len(data) != opened.st_size:
        raise ReviewRunnerError(f"Bound entrypoint read was incomplete: {path}")
    return {
        "path": str(path.resolve()),
        "size_bytes": opened.st_size,
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(stat.S_IMODE(opened.st_mode)),
    }


def bound_entrypoint_record(
    *,
    expected_script: Mapping[str, Any],
    expected_python: Mapping[str, Any],
    environment: Mapping[str, str],
) -> dict[str, Any]:
    argv0 = sys.argv[0]
    if (
        not argv0.startswith("/proc/self/fd/")
        or not argv0.removeprefix("/proc/self/fd/").isdigit()
        or sys.flags.isolated != 1
        or sys.flags.no_user_site != 1
        or dict(os.environ) != dict(environment)
    ):
        raise ReviewRunnerError("External review runner lacks bound isolated bootstrap")
    script_fd = int(argv0.removeprefix("/proc/self/fd/"))
    observed_script = artifact_from_fd(Path(__file__).resolve(), script_fd)
    python_fd = os.open("/proc/self/exe", os.O_RDONLY | os.O_CLOEXEC)
    try:
        observed_python = artifact_from_fd(
            Path(str(expected_python["path"])),
            python_fd,
        )
    finally:
        os.close(python_fd)
    if (
        observed_script != dict(expected_script)
        or observed_python != dict(expected_python)
    ):
        raise ReviewRunnerError("External review bound bootstrap identity mismatch")
    return {
        "schema_name": "intersubmod.bound_python_entrypoint",
        "schema_version": "1.0.0",
        "argv0": argv0,
        "script": observed_script,
        "python": observed_python,
        "isolated_mode": True,
        "no_user_site": True,
        "environment": dict(environment),
    }


def encode_json(payload: Any, *, pretty: bool = True) -> bytes:
    if pretty:
        text = json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True)
    else:
        text = json.dumps(
            payload,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    return (text + "\n").encode("utf-8")


def write_exclusive_readonly(path: Path, data: bytes) -> dict[str, Any]:
    output = path.parent.resolve(strict=True) / path.name
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(output, flags, 0o444)
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise ReviewRunnerError(f"Unable to write complete output: {output}")
            offset += written
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
        opened = os.fstat(descriptor)
        live = output.stat()
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_dev != live.st_dev
            or opened.st_ino != live.st_ino
            or opened.st_size != len(data)
            or opened.st_mode != live.st_mode
        ):
            raise ReviewRunnerError(f"Published output binding drift: {output}")
        return {
            "path": str(output),
            "size_bytes": opened.st_size,
            "sha256": hashlib.sha256(data).hexdigest(),
            "mode": oct(stat.S_IMODE(opened.st_mode)),
        }
    finally:
        os.close(descriptor)


def run_bound_git_head(reader: Any, release_module: Any) -> str:
    git_fd, _, git_record = reader.open(
        release_module.GIT_PATH, include_mode=True
    )
    if git_record["sha256"] != release_module.GIT_SHA256:
        raise ReviewRunnerError("Pinned Git verifier drift")
    completed = subprocess.run(
        [
            f"/proc/self/fd/{git_fd}",
            "-C",
            str(REPO_ROOT),
            "rev-parse",
            "HEAD",
        ],
        check=False,
        capture_output=True,
        text=True,
        pass_fds=(git_fd,),
        env={
            "PATH": "/usr/bin:/bin",
            "LANG": "C",
            "LC_ALL": "C",
            "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_CONFIG_GLOBAL": "/dev/null",
        },
    )
    if completed.returncode != 0:
        raise ReviewRunnerError(f"Bound Git HEAD query failed: {completed.stderr}")
    return completed.stdout.strip()


def output_paths(evidence_module: Any, reviewer_token: str) -> dict[str, Path]:
    attestation = evidence_module.REVIEW_ATTESTATION_SPECS[reviewer_token]["path"]
    prefix = attestation.with_suffix("")
    return {
        "request": Path(str(prefix) + ".request.md"),
        "command_record": Path(str(prefix) + ".command.json"),
        "transcript": Path(str(prefix) + ".stream.jsonl"),
        "stderr": Path(str(prefix) + ".stderr.txt"),
        "attestation": attestation,
        "signature": Path(str(attestation) + ".ed25519.sig"),
    }


def require_empty_output_slots(paths: dict[str, Path]) -> None:
    occupied = [str(path) for path in paths.values() if os.path.lexists(path)]
    if occupied:
        raise ReviewRunnerError(f"External review output slots are occupied: {occupied}")


def clean_environment() -> dict[str, str]:
    return {
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


def main() -> int:
    if len(sys.argv) != 3 or sys.argv[1] != "--reviewer" or sys.argv[2] not in (
        "A",
        "B",
    ):
        raise ReviewRunnerError("Usage: run_external...py --reviewer A|B")
    reviewer_token = sys.argv[2]
    with BoundArtifactReader() as reader:
        _, release_bytes, _ = reader.open(
            RELEASE_MODULE_PATH,
            include_mode=True,
        )
        _, evidence_bytes, _ = reader.open(
            EVIDENCE_MODULE_PATH,
            include_mode=True,
        )
        release_module = execute_module_from_bytes(
            "external_review_release_source_authority_bound",
            RELEASE_MODULE_PATH,
            release_bytes,
        )
        evidence_module = execute_module_from_bytes(
            "external_review_attested_evidence_bound",
            EVIDENCE_MODULE_PATH,
            evidence_bytes,
        )
        paths = output_paths(evidence_module, reviewer_token)
        require_empty_output_slots(paths)
        _, _, assembler_record = reader.open(
            ASSEMBLER_PATH, include_mode=True
        )
        _, _, runner_record = reader.open(Path(__file__), include_mode=True)
        _, _, evidence_record = reader.open(
            EVIDENCE_MODULE_PATH, include_mode=True
        )
        claude_fd, _, claude_record = reader.open(
            evidence_module.CLAUDE_CLI_PATH, include_mode=True
        )
        if (
            assembler_record["mode"] != "0o444"
            or runner_record["mode"] != "0o444"
            or evidence_record["mode"] != "0o444"
            or claude_record["sha256"] != evidence_module.CLAUDE_CLI_SHA256
            or claude_record["mode"] != "0o755"
        ):
            raise ReviewRunnerError("External-review executable source identity drift")
        git_head = run_bound_git_head(reader, release_module)
        support_paths = {
            "canonical_test_runner": evidence_module.TEST_RUNNER_PATH,
            "authority_assembler": evidence_module.AUTHORITY_ASSEMBLER_PATH,
            "review_normalizer": evidence_module.REVIEW_NORMALIZER_PATH,
            "external_review_runner": evidence_module.REVIEW_RUNNER_PATH,
            "attested_evidence_validator": evidence_module.EVIDENCE_VALIDATOR_PATH,
            "one_time_signer": evidence_module.SIGNER_PATH,
        }
        test_run = evidence_module.validate_signed_test_run(
            reader=reader,
            release_module=release_module,
            expected_head=git_head,
            expected_support_paths=support_paths,
        )
        environment = clean_environment()
        entrypoint = bound_entrypoint_record(
            expected_script=runner_record,
            expected_python=test_run["payload"]["runtime"]["python"],
            environment=environment,
        )
        contract = {
            "reviewer_token": reviewer_token,
            "source_set_sha256": test_run["source_set_sha256"],
            "git_head": git_head,
            "assembler": assembler_record,
            "test_run_manifest": test_run["records"]["manifest"],
            "test_run_signature": test_run["records"]["signature"],
            "test_run_public_key": test_run["records"]["public_key"],
            "canonical_junit": {
                "artifact": test_run["junit"],
                "counts": test_run["counts"],
            },
        }
        request_text = evidence_module.expected_review_request(
            reviewer_token,
            contract=contract,
        ).decode("utf-8")
        request_record = write_exclusive_readonly(
            paths["request"],
            request_text.encode("utf-8"),
        )
        reader.open(paths["request"], include_mode=True)
        schema = evidence_module.expected_review_schema(
            source_set_sha256=test_run["source_set_sha256"],
            git_head=git_head,
            assembler=assembler_record,
            junit=test_run["junit"],
            counts=test_run["counts"],
        )
        session_id = str(uuid.uuid4())
        system_prompt = evidence_module.REVIEW_SYSTEM_PROMPT
        command = [
            f"/proc/self/fd/{claude_fd}",
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
            ALLOWED_TOOLS,
            "--allowedTools",
            ALLOWED_TOOLS,
            "--disallowedTools",
            DISALLOWED_TOOLS,
            "--model",
            evidence_module.EXPECTED_CLAUDE_MODEL_ID,
            "--effort",
            "max",
            "--max-budget-usd",
            "12",
            "--output-format",
            "stream-json",
            "--verbose",
            "--json-schema",
            json.dumps(schema, separators=(",", ":"), sort_keys=True),
            "--session-id",
            session_id,
            "--name",
            f"InterSubMod Task-B external Reviewer {reviewer_token} v18",
            "--append-system-prompt",
            system_prompt,
            request_text,
        ]
        command_record_payload = {
            "schema_name": "intersubmod.external_claude_command_record",
            "schema_version": "1.0.0",
            "reviewer_token": reviewer_token,
            "session_id": session_id,
            "cwd": str(REPO_ROOT),
            "entrypoint": entrypoint,
            "claude_cli": claude_record,
            "runner": runner_record,
            "evidence_validator": evidence_record,
            "request": request_record,
            "model": evidence_module.EXPECTED_CLAUDE_MODEL_ID,
            "effort": "max",
            "permission_mode": "dontAsk",
            "allowed_tools": ALLOWED_TOOLS,
            "disallowed_tools": DISALLOWED_TOOLS,
            "environment": environment,
            "schema_sha256": hashlib.sha256(
                encode_json(schema, pretty=False)
            ).hexdigest(),
            "command": command,
        }
        command_record = write_exclusive_readonly(
            paths["command_record"],
            encode_json(command_record_payload),
        )
        reader.open(paths["command_record"], include_mode=True)
        reader.require_paths_still_bound()

        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=environment,
            capture_output=True,
            check=False,
            pass_fds=(claude_fd,),
            timeout=3600,
        )
        transcript_record = write_exclusive_readonly(
            paths["transcript"], completed.stdout
        )
        stderr_record = write_exclusive_readonly(
            paths["stderr"], completed.stderr
        )
        reader.open(paths["transcript"], include_mode=True)
        reader.open(paths["stderr"], include_mode=True)
        if completed.returncode != 0:
            raise ReviewRunnerError(
                f"Claude review process failed: returncode={completed.returncode}"
            )
        tool_trace = evidence_module.parse_claude_stream_transcript(
            completed.stdout,
            reviewer_token=reviewer_token,
            session_id=session_id,
        )
        review_payload = evidence_module.validate_clean_review_payload(
            tool_trace["structured_output"],
            reviewer_token=reviewer_token,
            source_set_sha256=test_run["source_set_sha256"],
            expected_head=git_head,
            assembler_record=assembler_record,
            junit_record=test_run["junit"],
            junit_counts=test_run["counts"],
        )
        key_spec = evidence_module.REVIEW_ATTESTATION_SPECS[reviewer_token]
        _, _, public_key = reader.open(
            key_spec["public_key"], include_mode=True
        )
        _, private_key = reader.open_metadata(
            key_spec["private_key"], include_mode=True
        )
        if (
            public_key["sha256"] != key_spec["public_key_sha256"]
            or public_key["mode"] != "0o444"
            or private_key["mode"] != "0o400"
        ):
            raise ReviewRunnerError("External-review signer key state drift")
        created_at = datetime.now(timezone.utc).isoformat()
        attestation = {
            "schema_name": "intersubmod.external_claude_process_attestation",
            "schema_version": "1.0.0",
            "reviewer_token": reviewer_token,
            "attestation_id": str(uuid.uuid4()),
            "created_at_utc": created_at,
            "process": {
                "session_id": session_id,
                "entrypoint": entrypoint,
                "claude_cli": claude_record,
                "command_record": command_record,
                "request": request_record,
                "transcript": transcript_record,
                "stderr": stderr_record,
                "tool_trace": tool_trace,
                "returncode": completed.returncode,
                "model": evidence_module.EXPECTED_CLAUDE_MODEL_ID,
                "effort": "max",
                "permission_mode": "dontAsk",
                "allowed_tools": ALLOWED_TOOLS,
                "disallowed_tools": DISALLOWED_TOOLS,
                "environment": environment,
            },
            "review_payload": review_payload,
            "review_payload_sha256": hashlib.sha256(
                evidence_module.encode_json(review_payload)
            ).hexdigest(),
            "reviewed_test_run": {
                "manifest": test_run["records"]["manifest"],
                "signature": test_run["records"]["signature"],
                "public_key": test_run["records"]["public_key"],
                "junit": test_run["junit"],
                "counts": test_run["counts"],
                "source_set_sha256": test_run["source_set_sha256"],
                "git_head": git_head,
            },
            "signature_contract": {
                "algorithm": "Ed25519",
                "public_key": public_key,
                "expected_signature_path": str(key_spec["signature"]),
                "private_key_path": str(key_spec["private_key"].resolve()),
                "required_private_key_mode_after_signing": "0o0",
            },
            "checks": {
                "pinned_claude_cli_identity_verified": True,
                "clean_environment_and_read_only_tools_enforced": True,
                "command_request_session_transcript_linked": True,
                "required_read_tool_events_verified": True,
                "bound_runner_entrypoint_verified": True,
                "signed_test_run_reverified": True,
                "review_payload_schema_and_clean_approval_verified": True,
                "attestation_requires_detached_signature": True,
            },
            "pass": True,
        }
        reader.require_paths_still_bound()
        attestation_record = write_exclusive_readonly(
            paths["attestation"],
            encode_json(attestation),
        )
        reader.open(paths["attestation"], include_mode=True)
        reader.retain_until_process_exit()

    signing_instruction = "SIGN " + json.dumps(
        attestation_record,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )
    print(
        json.dumps(
            {
                "reviewer": reviewer_token,
                "session_id": session_id,
                "attestation": attestation_record,
                "review_payload": review_payload,
                "signing_instruction": signing_instruction,
                "pass": True,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.TimeoutExpired as error:
        print(f"EXTERNAL_REVIEW_TIMEOUT error={error}", file=sys.stderr)
        raise SystemExit(124)
