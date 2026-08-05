#!/usr/bin/env python3
"""Run the final Git-environment-alignment canonical pytest suite."""

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
TEST_ROOT = TOPIC_ROOT / "tests"
RELEASE_MODULE_PATH = TOPIC_ROOT / "scripts" / "release_source_authority.py"
EVIDENCE_MODULE_PATH = HERE / "attested_release_evidence_v5.py"
PYTHON = Path(
    "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python"
)
XML_PATH = (
    TOPIC_ROOT
    / "logs"
    / "pytest_full_pre_authority_v27_git_env_alignment_attested_canonical.xml"
)
RAW_OUTPUT_DIR = (
    TOPIC_ROOT / "logs" / "pytest_full_pre_authority_v27_git_env_alignment_raw_evidence"
)
RAW_XML_PATH = RAW_OUTPUT_DIR / "pytest.xml"
BOUND_TEST_DIR = TOPIC_ROOT / ".pytest_full_pre_authority_v27_git_env_alignment_bound_test_fds"
PYTHON_CACHE_DIR = RAW_OUTPUT_DIR / "python_cache"
STDOUT_PATH = XML_PATH.with_suffix(".stdout.txt")
STDERR_PATH = XML_PATH.with_suffix(".stderr.txt")
MANIFEST_PATH = (
    TOPIC_ROOT
    / "logs"
    / "pytest_full_pre_authority_v27_git_env_alignment_test_run_manifest.v1.json"
)
MANIFEST_SIGNATURE_PATH = Path(str(MANIFEST_PATH) + ".ed25519.sig")
TEST_PUBLIC_KEY = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_test_run_authority/"
    "20260722_all_ssnv_v11_git_env_alignment_bootstrap/ed25519_public.pem"
)
TEST_PRIVATE_KEY = TEST_PUBLIC_KEY.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
TEST_PUBLIC_KEY_SHA256 = (
    "92255c453329b8b21dda721048259e8cd567addb5b3f0ee2aa03fb414a9d5015"
)
SUPPORT_PATHS = {
    "canonical_test_runner": Path(__file__).resolve(),
    "authority_assembler": HERE / "assemble_release_source_authority_v15.py",
    "review_normalizer": HERE / "normalize_external_source_review_v20.py",
    "external_review_runner": HERE / "run_external_claude_review_v22_attested.py",
    "attested_evidence_validator": EVIDENCE_MODULE_PATH,
    "one_time_signer": HERE / "run_one_time_ed25519_signer_v6.py",
}
ENVIRONMENT = {
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


class TestRunError(RuntimeError):
    """Raised when the canonical test-run evidence cannot be trusted."""


_PROCESS_LIFETIME_READERS: list["BoundArtifactReader"] = []


class BoundArtifactReader:
    """Bind test inputs before executing any trust-bearing module."""

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
            raise TestRunError(f"Test input is not regular: {resolved}")
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
            raise TestRunError(f"Test input changed while binding: {resolved}")
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
            raise TestRunError("O_PATH is required for signer-key binding")
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
            raise TestRunError(f"Signer-key metadata binding drift: {resolved}")
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
                raise TestRunError(f"Test input binding changed: {path}")

    def close(self) -> None:
        while self._opened:
            _, descriptor, _ = self._opened.pop()
            os.close(descriptor)

    def retain_until_process_exit(self) -> None:
        if self._retained:
            raise TestRunError("Test input reader was already retained")
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
        raise TestRunError(f"Bound entrypoint identity drift: {path}")
    data = b"".join(
        os.pread(
            descriptor,
            min(8 * 1024 * 1024, opened.st_size - offset),
            offset,
        )
        for offset in range(0, opened.st_size, 8 * 1024 * 1024)
    )
    if len(data) != opened.st_size:
        raise TestRunError(f"Bound entrypoint read was incomplete: {path}")
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
) -> dict[str, Any]:
    argv0 = sys.argv[0]
    if (
        not argv0.startswith("/proc/self/fd/")
        or not argv0.removeprefix("/proc/self/fd/").isdigit()
        or sys.flags.isolated != 1
        or sys.flags.no_user_site != 1
        or dict(os.environ) != ENVIRONMENT
    ):
        raise TestRunError("Canonical pytest runner lacks bound isolated bootstrap")
    script_fd = int(argv0.removeprefix("/proc/self/fd/"))
    observed_script = artifact_from_fd(
        SUPPORT_PATHS["canonical_test_runner"],
        script_fd,
    )
    python_fd = os.open("/proc/self/exe", os.O_RDONLY | os.O_CLOEXEC)
    try:
        observed_python = artifact_from_fd(PYTHON, python_fd)
    finally:
        os.close(python_fd)
    if (
        observed_script != dict(expected_script)
        or observed_python != dict(expected_python)
    ):
        raise TestRunError("Canonical pytest bound bootstrap identity mismatch")
    return {
        "schema_name": "intersubmod.bound_python_entrypoint",
        "schema_version": "1.0.0",
        "argv0": argv0,
        "script": observed_script,
        "python": observed_python,
        "isolated_mode": True,
        "no_user_site": True,
        "environment": ENVIRONMENT,
    }


def encode_json(payload: Any) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def write_exclusive_readonly(path: Path, data: bytes) -> dict[str, Any]:
    resolved_parent = path.parent.resolve(strict=True)
    output = resolved_parent / path.name
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(output, flags, 0o444)
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise TestRunError(f"Unable to write complete output: {output}")
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
            raise TestRunError(f"Published output binding drift: {output}")
        return {
            "path": str(output),
            "size_bytes": opened.st_size,
            "sha256": hashlib.sha256(data).hexdigest(),
            "mode": oct(stat.S_IMODE(opened.st_mode)),
        }
    finally:
        os.close(descriptor)


def require_empty_output_slots() -> None:
    slots = (
        XML_PATH,
        STDOUT_PATH,
        STDERR_PATH,
        MANIFEST_PATH,
        MANIFEST_SIGNATURE_PATH,
        RAW_OUTPUT_DIR,
        BOUND_TEST_DIR,
    )
    occupied = [str(path) for path in slots if os.path.lexists(path)]
    if occupied:
        raise TestRunError(f"Canonical test output slots are occupied: {occupied}")


def canonical_command(
    python_executable: str,
    test_entrypoints: list[str] | None = None,
) -> list[str]:
    command = [
        python_executable,
        "-I",
        "-X",
        f"pycache_prefix={PYTHON_CACHE_DIR}",
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
    command.extend(test_entrypoints or [str(TEST_ROOT)])
    command.append(f"--junitxml={RAW_XML_PATH}")
    return command


def run_bound_git_head(reader: Any, release_module: Any) -> tuple[str, dict[str, Any]]:
    git_fd, _, git_record = reader.open(
        release_module.GIT_PATH, include_mode=True
    )
    if git_record["sha256"] != release_module.GIT_SHA256:
        raise TestRunError("Pinned Git verifier drift")
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
        raise TestRunError(f"Bound Git HEAD query failed: {completed.stderr}")
    return completed.stdout.strip(), git_record


def build_manifest_payload(
    *,
    run_id: str,
    started_at: str,
    finished_at: str,
    git_head: str,
    sources: Mapping[str, Mapping[str, Any]],
    test_sources: list[dict[str, Any]],
    support: Mapping[str, Mapping[str, Any]],
    runtime: Mapping[str, Any],
    git_record: Mapping[str, Any],
    public_key: Mapping[str, Any],
    junit: Mapping[str, Any],
    counts: Mapping[str, int],
    stdout_record: Mapping[str, Any],
    stderr_record: Mapping[str, Any],
    raw_junit_record: Mapping[str, Any],
    evidence_module: Any,
    command: list[str],
    entrypoint: Mapping[str, Any],
    pytest_process: Mapping[str, Any],
    inventory: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema_name": "intersubmod.attested_pytest_run",
        "schema_version": "1.0.0",
        "task_type": "B_comprehensive_validation",
        "scope": "all_protected_release_sources_and_canonical_topic_tests",
        "run_id": run_id,
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "cwd": str(REPO_ROOT),
        "command": command,
        "environment": ENVIRONMENT,
        "git_head": git_head,
        "protected_sources": dict(sources),
        "source_set_sha256": evidence_module.source_set_digest(sources),
        "test_sources": test_sources,
        "supporting_sources": dict(support),
        "runtime": {
            "python": dict(runtime),
            "git": dict(git_record),
            "entrypoint": dict(entrypoint),
            "pytest_process": dict(pytest_process),
        },
        "junit": {
            "artifact": dict(junit),
            "counts": dict(counts),
            "testcase_inventory": dict(inventory),
        },
        "captured_output": {
            "stdout": dict(stdout_record),
            "stderr": dict(stderr_record),
            "raw_junit": dict(raw_junit_record),
        },
        "signature_contract": {
            "algorithm": "Ed25519",
            "public_key": dict(public_key),
            "expected_signature_path": str(MANIFEST_SIGNATURE_PATH),
            "private_key_path": str(TEST_PRIVATE_KEY.resolve()),
            "required_private_key_mode_after_signing": "0o0",
        },
        "checks": {
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
        },
        "pass": True,
    }


def main() -> int:
    require_empty_output_slots()
    os.mkdir(RAW_OUTPUT_DIR, mode=0o700)
    os.mkdir(BOUND_TEST_DIR, mode=0o700)
    os.mkdir(PYTHON_CACHE_DIR, mode=0o700)
    raw_directory_stat = RAW_OUTPUT_DIR.stat()
    if (
        not stat.S_ISDIR(raw_directory_stat.st_mode)
        or stat.S_IMODE(raw_directory_stat.st_mode) != 0o700
    ):
        raise TestRunError("Raw pytest evidence directory is not mode 0700")
    for label, directory in (
        ("Bound-test FD", BOUND_TEST_DIR),
        ("Python cache", PYTHON_CACHE_DIR),
    ):
        observed = directory.stat()
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_IMODE(observed.st_mode) != 0o700
        ):
            raise TestRunError(f"{label} directory is not mode 0700")
    run_id = str(uuid.uuid4())
    started_at = datetime.now(timezone.utc).isoformat()

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
            "attested_pytest_release_source_authority_bound",
            RELEASE_MODULE_PATH,
            release_bytes,
        )
        evidence_module = execute_module_from_bytes(
            "attested_pytest_evidence_validator_bound",
            EVIDENCE_MODULE_PATH,
            evidence_bytes,
        )
        sources = {
            role: reader.open(path, include_mode=True)[2]
            for role, path in release_module.EXPECTED_SOURCE_PATHS.items()
        }
        if any(record["mode"] != "0o444" for record in sources.values()):
            raise TestRunError("Protected source is not mode 0444")
        test_paths = sorted(TEST_ROOT.glob("test_*.py"))
        if not test_paths:
            raise TestRunError("Canonical test root contains no test files")
        test_bindings = [
            reader.open(path, include_mode=True) for path in test_paths
        ]
        test_sources = [binding[2] for binding in test_bindings]
        if any(record["mode"] != "0o444" for record in test_sources):
            raise TestRunError("Canonical test source is not mode 0444")
        test_entrypoints: list[str] = []
        for path, binding in zip(test_paths, test_bindings):
            entrypoint = BOUND_TEST_DIR / path.name
            os.symlink(f"/proc/self/fd/{binding[0]}", entrypoint)
            test_entrypoints.append(str(entrypoint))
        support = {
            role: reader.open(path, include_mode=True)[2]
            for role, path in SUPPORT_PATHS.items()
        }
        if any(record["mode"] != "0o444" for record in support.values()):
            raise TestRunError("Test-run supporting source is not mode 0444")
        runtime_fd, _, runtime = reader.open(PYTHON, include_mode=True)
        entrypoint = bound_entrypoint_record(
            expected_script=support["canonical_test_runner"],
            expected_python=runtime,
        )
        pytest_executable = f"/proc/self/fd/{runtime_fd}"
        command = canonical_command(str(runtime["path"]), test_entrypoints)
        pytest_process = {
            "argv0": command[0],
            "executable": pytest_executable,
            "python": runtime,
        }
        _, _, public_key = reader.open(TEST_PUBLIC_KEY, include_mode=True)
        _, private_key = reader.open_metadata(
            TEST_PRIVATE_KEY, include_mode=True
        )
        if (
            public_key["sha256"] != TEST_PUBLIC_KEY_SHA256
            or public_key["mode"] != "0o444"
            or private_key["mode"] != "0o400"
        ):
            raise TestRunError("Test-run signer key state drift")
        git_head, git_record = run_bound_git_head(reader, release_module)
        reader.require_paths_still_bound()

        completed = subprocess.run(
            command,
            executable=pytest_executable,
            cwd=REPO_ROOT,
            env=ENVIRONMENT,
            check=False,
            capture_output=True,
            pass_fds=(runtime_fd, *(binding[0] for binding in test_bindings)),
        )
        stdout_record = write_exclusive_readonly(STDOUT_PATH, completed.stdout)
        stderr_record = write_exclusive_readonly(STDERR_PATH, completed.stderr)
        if not RAW_XML_PATH.exists():
            raise TestRunError("pytest did not produce raw JUnit XML")
        RAW_XML_PATH.chmod(0o444)
        _, xml_bytes, raw_junit = reader.open(
            RAW_XML_PATH,
            include_mode=True,
        )
        counts = evidence_module.parse_junit_counts(xml_bytes)
        inventory = evidence_module.parse_junit_inventory(xml_bytes)
        if (
            completed.returncode != 0
            or counts["tests"] != evidence_module.EXPECTED_TESTCASE_COUNT
            or inventory["count"] != evidence_module.EXPECTED_TESTCASE_COUNT
            or inventory["sha256"]
            != evidence_module.EXPECTED_TESTCASE_INVENTORY_SHA256
            or counts["failures"] != 0
            or counts["errors"] != 0
            or counts["skipped"] != 0
        ):
            raise TestRunError(
                f"Canonical pytest failed: returncode={completed.returncode}, "
                f"counts={counts}"
            )
        junit = write_exclusive_readonly(XML_PATH, xml_bytes)
        reader.open(XML_PATH, include_mode=True)
        if (
            junit["sha256"] != raw_junit["sha256"]
            or junit["size_bytes"] != raw_junit["size_bytes"]
        ):
            raise TestRunError("Canonical JUnit differs from raw pytest evidence")
        finished_at = datetime.now(timezone.utc).isoformat()
        manifest = build_manifest_payload(
            run_id=run_id,
            started_at=started_at,
            finished_at=finished_at,
            git_head=git_head,
            sources=sources,
            test_sources=test_sources,
            support=support,
            runtime=runtime,
            git_record=git_record,
            public_key=public_key,
            junit=junit,
            counts=counts,
            stdout_record=stdout_record,
            stderr_record=stderr_record,
            raw_junit_record=raw_junit,
            evidence_module=evidence_module,
            command=command,
            entrypoint=entrypoint,
            pytest_process=pytest_process,
            inventory=inventory,
        )
        live_raw_directory = RAW_OUTPUT_DIR.stat()
        if (
            live_raw_directory.st_dev != raw_directory_stat.st_dev
            or live_raw_directory.st_ino != raw_directory_stat.st_ino
            or live_raw_directory.st_mode != raw_directory_stat.st_mode
            or live_raw_directory.st_nlink != raw_directory_stat.st_nlink
        ):
            raise TestRunError("Raw pytest evidence directory binding drifted")
        reader.require_paths_still_bound()
        predicted = {
            "path": str(MANIFEST_PATH),
            "size_bytes": len(encode_json(manifest)),
            "sha256": hashlib.sha256(encode_json(manifest)).hexdigest(),
            "mode": "0o444",
        }
        observed = write_exclusive_readonly(MANIFEST_PATH, encode_json(manifest))
        if observed != predicted:
            raise TestRunError("Published test-run manifest differs from prediction")
        reader.open(MANIFEST_PATH, include_mode=True)
        reader.retain_until_process_exit()

    signing_instruction = "SIGN " + json.dumps(
        observed,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )
    print(
        json.dumps(
            {
                "manifest": observed,
                "junit": junit,
                "counts": counts,
                "source_set_sha256": manifest["source_set_sha256"],
                "signing_instruction": signing_instruction,
                "pass": True,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
