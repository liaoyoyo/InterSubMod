#!/usr/bin/env python3
"""Assemble v7 authority after final Git subprocess environment alignment."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
from types import ModuleType
from typing import Any, Mapping


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
HERE = TOPIC_ROOT / "audit_notes"
RELEASE_MODULE_PATH = TOPIC_ROOT / "scripts" / "release_source_authority.py"
EVIDENCE_MODULE_PATH = HERE / "attested_release_evidence_v5.py"
REVIEW_NORMALIZER_PATH = HERE / "normalize_external_source_review_v20.py"
TEST_RUNNER_PATH = HERE / "run_attested_canonical_pytest_v5.py"
REVIEW_RUNNER_PATH = HERE / "run_external_claude_review_v22_attested.py"
SIGNER_SOURCE = HERE / "run_one_time_ed25519_signer_v6.py"
AUTHORITY_PATH = (
    REPO_ROOT
    / "docs"
    / "provenance"
    / "source_authorities"
    / "20260722_all_ssnv_focal_alt_release_source_authority.v7.json"
)
APPROVAL_PATH = AUTHORITY_PATH.with_name(
    "20260722_all_ssnv_focal_alt_release_source_authority.v7.approval.json"
)
APPROVAL_SIGNATURE_PATH = Path(str(APPROVAL_PATH) + ".ed25519.sig")
PUBLIC_KEY_PATH = Path(
    "/bip7_disk/liaoyoyo2001/.config/intersubmod_release_authority/"
    "20260722_all_ssnv_v10_strict_command_parity_bootstrap/ed25519_public.pem"
)
PRIVATE_KEY_PATH = PUBLIC_KEY_PATH.with_name(
    "ed25519_private_encrypted_unrecoverable_after_signing.pem"
)
EXPECTED_PUBLIC_KEY_SHA256 = (
    "cecb2287cf87f8bd948c7390bd5f7059578966e0fc917de420572ddd82d8f21b"
)
EXPECTED_SIGNER_SOURCE_SHA256 = (
    "c9c32222d4e8ca6daad8f6b3f8a91746fc99dd200d24cf3d77e500234c0f8bee"
)
EXPECTED_AUTHORITY_ID = "20260722_all_ssnv_focal_alt_task_b_release_v7_strict_command_parity"
EXPECTED_APPROVAL_STATUS = "APPROVED_FOR_FULL_TASK_B_RUN"
EXPECTED_SCHEMA_NAME = "intersubmod.release_source_authority"
EXPECTED_SCHEMA_VERSION = "1.0.0"
EXPECTED_SCOPE = "all_7_datasets_469849_sites_and_release_chain"
EXPECTED_CLAIM_CONTRACT_SIZE = 9_144
EXPECTED_CLAIM_CONTRACT_SHA256 = (
    "da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af"
)
ENTRYPOINT_ENVIRONMENT = {
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
EXPECTED_ROLE_NAMES = frozenset(
    {
        "release_source_authority_validator",
        "cooccurrence_analyzer",
        "cooccurrence_preflight",
        "cooccurrence_release_finalizer",
        "cooccurrence_release_runner",
        "ssnv_cooccurrence_lib",
        "latest_tag_join",
        "m2_screen_gate",
        "independent_m2_auditor",
        "primary_artifact_auditor",
        "strict_producer",
        "focal_alt_cluster_lib",
        "completion_runner",
        "final_dataset_builder",
        "final_result_release_finalizer",
        "matched_normal_runner",
        "matched_normal_analyzer",
        "cn_ccf_annotator",
        "frozen_input_auditor",
        "tumor_ref_source_verifier",
        "report_builder",
        "portable_report_deliverer",
        "portable_layout_qa",
    }
)


class AssemblyError(RuntimeError):
    """Raised when a release authority cannot be assembled fail-closed."""


_PROCESS_LIFETIME_READERS: list["BoundArtifactReader"] = []


class BoundArtifactReader:
    """Bind every authority input by descriptor for the full assembly."""

    def __init__(self) -> None:
        self._opened: list[tuple[Path, int, os.stat_result]] = []
        self._directories: list[tuple[Path, int, os.stat_result]] = []
        self._retained = False

    def __enter__(self) -> "BoundArtifactReader":
        return self

    def __exit__(self, *_: Any) -> None:
        if not self._retained:
            self.close()

    def close(self) -> None:
        while self._opened:
            _, descriptor, _ = self._opened.pop()
            os.close(descriptor)
        while self._directories:
            _, descriptor, _ = self._directories.pop()
            os.close(descriptor)

    @staticmethod
    def _same_identity(left: os.stat_result, right: os.stat_result) -> bool:
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
            raise AssemblyError(f"Authority input is not regular: {resolved}")
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
            or not self._same_identity(first, second)
            or not self._same_identity(second, live)
        ):
            os.close(descriptor)
            raise AssemblyError(f"Authority input changed while binding: {resolved}")
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
            raise AssemblyError("O_PATH is required for retired-key binding")
        resolved = path.expanduser().resolve(strict=True)
        flags = os.O_PATH | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(resolved, flags)
        opened = os.fstat(descriptor)
        live = resolved.stat()
        if (
            not stat.S_ISREG(opened.st_mode)
            or not self._same_identity(opened, live)
        ):
            os.close(descriptor)
            raise AssemblyError(f"Metadata-only input binding drift: {resolved}")
        record: dict[str, Any] = {"path": str(resolved)}
        if include_mode:
            record["mode"] = oct(stat.S_IMODE(opened.st_mode))
        self._opened.append((resolved, descriptor, opened))
        return descriptor, record

    def open_directory(self, path: Path) -> tuple[int, os.stat_result]:
        resolved = path.expanduser().resolve(strict=True)
        flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(resolved, flags)
        opened = os.fstat(descriptor)
        live = resolved.stat()
        if (
            not stat.S_ISDIR(opened.st_mode)
            or not self._same_identity(opened, live)
        ):
            os.close(descriptor)
            raise AssemblyError(f"Authority directory binding drift: {resolved}")
        self._directories.append((resolved, descriptor, opened))
        return descriptor, opened

    def require_paths_still_bound(self) -> None:
        for path, descriptor, opened in self._opened:
            if (
                not self._same_identity(os.fstat(descriptor), opened)
                or not self._same_identity(path.stat(), opened)
            ):
                raise AssemblyError(f"Authority input binding changed: {path}")
        for path, descriptor, opened in self._directories:
            if (
                not self._same_identity(os.fstat(descriptor), opened)
                or not self._same_identity(path.stat(), opened)
            ):
                raise AssemblyError(f"Authority directory binding changed: {path}")

    def retain_until_process_exit(self) -> None:
        if self._retained:
            raise AssemblyError("Authority reader was already retained")
        self.require_paths_still_bound()
        _PROCESS_LIFETIME_READERS.append(self)
        self._retained = True


def execute_module_from_bytes(name: str, path: Path, data: bytes) -> ModuleType:
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
        or not BoundArtifactReader._same_identity(opened, live)
    ):
        raise AssemblyError(f"Bound entrypoint identity drift: {path}")
    data = b"".join(
        os.pread(
            descriptor,
            min(8 * 1024 * 1024, opened.st_size - offset),
            offset,
        )
        for offset in range(0, opened.st_size, 8 * 1024 * 1024)
    )
    if len(data) != opened.st_size:
        raise AssemblyError(f"Bound entrypoint read was incomplete: {path}")
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
        or dict(os.environ) != ENTRYPOINT_ENVIRONMENT
    ):
        raise AssemblyError("Authority assembler lacks bound isolated bootstrap")
    observed_script = artifact_from_fd(
        Path(__file__).resolve(),
        int(argv0.removeprefix("/proc/self/fd/")),
    )
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
        raise AssemblyError("Authority assembler bound bootstrap identity mismatch")
    return {
        "schema_name": "intersubmod.bound_python_entrypoint",
        "schema_version": "1.0.0",
        "argv0": argv0,
        "script": observed_script,
        "python": observed_python,
        "isolated_mode": True,
        "no_user_site": True,
        "environment": ENTRYPOINT_ENVIRONMENT,
    }


def encode_json(payload: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def predicted_readonly_record(path: Path, data: bytes) -> dict[str, Any]:
    output = path.parent.resolve(strict=True) / path.name
    return {
        "path": str(output),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": "0o444",
    }


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
                raise AssemblyError(f"Unable to publish complete artifact: {output}")
            offset += written
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
        opened = os.fstat(descriptor)
        live = output.stat()
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_dev != live.st_dev
            or opened.st_ino != live.st_ino
            or opened.st_mode != live.st_mode
            or opened.st_size != len(data)
        ):
            raise AssemblyError(f"Published authority binding drift: {output}")
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
        AUTHORITY_PATH,
        APPROVAL_PATH,
        APPROVAL_SIGNATURE_PATH,
        Path(str(APPROVAL_SIGNATURE_PATH) + ".pending"),
    )
    occupied = [str(path) for path in slots if os.path.lexists(path)]
    if occupied:
        raise AssemblyError(f"Authority output slots are occupied: {occupied}")


def run_bound_git_head(
    reader: Any,
    release_module: ModuleType,
) -> tuple[str, dict[str, Any]]:
    git_fd, _, git_record = reader.open(
        release_module.GIT_PATH,
        include_mode=True,
    )
    if git_record["sha256"] != release_module.GIT_SHA256:
        raise AssemblyError("Pinned Git verifier identity drifted")
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
        raise AssemblyError(f"Bound Git HEAD query failed: {completed.stderr}")
    return completed.stdout.strip(), git_record


def support_paths() -> dict[str, Path]:
    return {
        "canonical_test_runner": TEST_RUNNER_PATH,
        "authority_assembler": Path(__file__).resolve(),
        "review_normalizer": REVIEW_NORMALIZER_PATH,
        "external_review_runner": REVIEW_RUNNER_PATH,
        "attested_evidence_validator": EVIDENCE_MODULE_PATH,
        "one_time_signer": SIGNER_SOURCE,
    }


def require_constant_bindings(
    release_module: ModuleType,
    evidence_module: ModuleType,
) -> None:
    expected = {
        "AUTHORITY_PATH": AUTHORITY_PATH,
        "APPROVAL_PATH": APPROVAL_PATH,
        "APPROVAL_SIGNATURE_PATH": APPROVAL_SIGNATURE_PATH,
        "PUBLIC_KEY_PATH": PUBLIC_KEY_PATH,
        "PRIVATE_KEY_PATH": PRIVATE_KEY_PATH,
        "PUBLIC_KEY_SHA256": EXPECTED_PUBLIC_KEY_SHA256,
        "AUTHORITY_ID": EXPECTED_AUTHORITY_ID,
        "APPROVAL_STATUS": EXPECTED_APPROVAL_STATUS,
        "AUTHORITY_SCHEMA_NAME": EXPECTED_SCHEMA_NAME,
        "AUTHORITY_SCHEMA_VERSION": EXPECTED_SCHEMA_VERSION,
        "CLAIM_CONTRACT_SIZE": EXPECTED_CLAIM_CONTRACT_SIZE,
        "CLAIM_CONTRACT_SHA256": EXPECTED_CLAIM_CONTRACT_SHA256,
        "AUTHORITY_ASSEMBLER_PATH": Path(__file__).resolve(),
        "CANONICAL_PYTEST_XML_PATH": evidence_module.TEST_XML_PATH,
    }
    for field, expected_value in expected.items():
        observed = getattr(release_module, field, None)
        if isinstance(expected_value, Path):
            if not isinstance(observed, Path) or observed.resolve(
                strict=False
            ) != expected_value.resolve(strict=False):
                raise AssemblyError(f"Release-module path binding drifted: {field}")
        elif observed != expected_value:
            raise AssemblyError(f"Release-module constant drifted: {field}")
    if frozenset(release_module.EXPECTED_SOURCE_PATHS) != EXPECTED_ROLE_NAMES:
        raise AssemblyError("Release-module protected role set drifted")
    expected_reviews = tuple(
        evidence_module.REVIEW_ATTESTATION_SPECS[token]["path"]
        for token in ("A", "B")
    )
    if tuple(release_module.EXPECTED_EXTERNAL_REVIEW_PATHS) != expected_reviews:
        raise AssemblyError("Release-module signed review paths drifted")


def approval_from_review(
    verified: Mapping[str, Any],
) -> dict[str, Any]:
    payload = verified["review_payload"]
    return {
        "reviewer_id": payload["reviewer_id"],
        "reviewer_label": payload["reviewer_label"],
        "model": payload["model"],
        "verdict": payload["verdict"],
        "findings_closed": payload["findings_closed"],
        "review_scope": payload["review_scope"],
        "reviewed_source_set_sha256": payload["reviewed_source_set_sha256"],
        "reviewed_git_head": payload["reviewed_git_head"],
        "reviewed_assembler": payload["reviewed_assembler"],
        "canonical_junit": payload["canonical_junit"],
        "process_attestation": verified["attestation"],
        "process_attestation_signature": verified["signature"],
        "process_attestation_public_key": verified["public_key"],
        "process_attestation_session_id": verified["session_id"],
        "process_records": verified["process_records"],
    }


def main() -> int:
    require_empty_output_slots()
    with BoundArtifactReader() as reader:
        _, assembler_bytes, assembler_record = reader.open(
            Path(__file__),
            include_mode=True,
        )
        _, release_bytes, release_record = reader.open(
            RELEASE_MODULE_PATH,
            include_mode=True,
        )
        _, evidence_bytes, evidence_record = reader.open(
            EVIDENCE_MODULE_PATH,
            include_mode=True,
        )
        _, normalizer_bytes, normalizer_record = reader.open(
            REVIEW_NORMALIZER_PATH,
            include_mode=True,
        )
        _, signer_bytes, signer_record = reader.open(
            SIGNER_SOURCE,
            include_mode=True,
        )
        if any(
            record["mode"] != "0o444"
            for record in (
                assembler_record,
                release_record,
                evidence_record,
                normalizer_record,
                signer_record,
            )
        ):
            raise AssemblyError("Authority support source is not mode 0444")
        if signer_record["sha256"] != EXPECTED_SIGNER_SOURCE_SHA256:
            raise AssemblyError("One-time signer implementation identity drifted")

        release_module = execute_module_from_bytes(
            "all_ssnv_release_source_authority_v13_bound",
            RELEASE_MODULE_PATH,
            release_bytes,
        )
        evidence_module = execute_module_from_bytes(
            "all_ssnv_attested_release_evidence_v5_bound",
            EVIDENCE_MODULE_PATH,
            evidence_bytes,
        )
        normalizer = execute_module_from_bytes(
            "all_ssnv_external_review_normalizer_v18_bound",
            REVIEW_NORMALIZER_PATH,
            normalizer_bytes,
        )
        signer = execute_module_from_bytes(
            "all_ssnv_digest_bound_signer_v6",
            SIGNER_SOURCE,
            signer_bytes,
        )
        require_constant_bindings(release_module, evidence_module)

        current_head, git_record = run_bound_git_head(reader, release_module)
        test_run = evidence_module.validate_signed_test_run(
            reader=reader,
            release_module=release_module,
            expected_head=current_head,
            expected_support_paths=support_paths(),
        )
        assembler_entrypoint = bound_entrypoint_record(
            expected_script=assembler_record,
            expected_python=test_run["payload"]["runtime"]["python"],
        )
        evidence_module.validate_entrypoint_record(
            assembler_entrypoint,
            expected_script=assembler_record,
            expected_python=test_run["payload"]["runtime"]["python"],
            expected_environment=ENTRYPOINT_ENVIRONMENT,
            label="source authority assembler",
        )
        sources = {
            role: reader.open(path, include_mode=True)[2]
            for role, path in release_module.EXPECTED_SOURCE_PATHS.items()
        }
        if (
            set(sources) != set(EXPECTED_ROLE_NAMES)
            or any(record["mode"] != "0o444" for record in sources.values())
            or evidence_module.source_set_digest(sources)
            != test_run["source_set_sha256"]
        ):
            raise AssemblyError("Live protected sources differ from signed test run")
        if (
            test_run["supporting_sources"]["authority_assembler"]
            != assembler_record
            or test_run["supporting_sources"]["attested_evidence_validator"]
            != evidence_record
            or test_run["supporting_sources"]["review_normalizer"]
            != normalizer_record
            or test_run["supporting_sources"]["one_time_signer"]
            != signer_record
        ):
            raise AssemblyError("Live authority support differs from signed test run")

        reviews = [
            normalizer.validate_attested_review(
                reader=reader,
                release_module=release_module,
                reviewer_token=token,
                expected_head=current_head,
                source_set_sha256=test_run["source_set_sha256"],
                assembler_record=assembler_record,
                test_run=test_run,
                evidence_module=evidence_module,
            )
            for token in ("A", "B")
        ]
        reviewer_ids = [
            review["review_payload"]["reviewer_id"] for review in reviews
        ]
        session_ids = [review["session_id"] for review in reviews]
        public_keys = [
            review["public_key"]["sha256"] for review in reviews
        ]
        if (
            len(set(reviewer_ids)) != 2
            or len(set(session_ids)) != 2
            or len(set(public_keys)) != 2
        ):
            raise AssemblyError(
                "External reviews do not have distinct reviewer/session/key identities"
            )

        _, _, claim_contract = reader.open(
            release_module.CLAIM_CONTRACT_PATH,
            include_mode=True,
        )
        if (
            claim_contract["size_bytes"] != EXPECTED_CLAIM_CONTRACT_SIZE
            or claim_contract["sha256"] != EXPECTED_CLAIM_CONTRACT_SHA256
            or claim_contract["mode"] != "0o444"
        ):
            raise AssemblyError("Canonical claim contract identity drifted")
        _, key_directory = reader.open_directory(PUBLIC_KEY_PATH.parent)
        _, _, public_key = reader.open(PUBLIC_KEY_PATH, include_mode=True)
        _, private_key = reader.open_metadata(
            PRIVATE_KEY_PATH,
            include_mode=True,
        )
        if (
            public_key["sha256"] != EXPECTED_PUBLIC_KEY_SHA256
            or public_key["mode"] != "0o444"
            or private_key["mode"] != "0o400"
            or stat.S_IMODE(key_directory.st_mode) != 0o700
        ):
            raise AssemblyError("Source-authority signer key state drifted")

        review_approvals = [approval_from_review(review) for review in reviews]
        review_evidence = {
            "authority_assembler": assembler_record,
            "authority_assembler_entrypoint": assembler_entrypoint,
            "attested_test_run": {
                "manifest": test_run["records"]["manifest"],
                "signature": test_run["records"]["signature"],
                "public_key": test_run["records"]["public_key"],
                "retired_private_key": test_run["records"]["private_key"],
                "junit": test_run["junit"],
                "counts": test_run["counts"],
                "supporting_sources": test_run["supporting_sources"],
            },
            "external_review_normalizer": normalizer_record,
            "attested_evidence_validator": evidence_record,
            "canonical_test_runner": test_run["supporting_sources"][
                "canonical_test_runner"
            ],
            "external_review_runner": test_run["supporting_sources"][
                "external_review_runner"
            ],
            "one_time_signer_implementation": signer_record,
            "repository_head_verifier": git_record,
            "external_review_attestations": [
                review["attestation"] for review in reviews
            ],
            "external_review_attestation_signatures": [
                review["signature"] for review in reviews
            ],
            "external_review_public_keys": [
                review["public_key"] for review in reviews
            ],
            "external_review_process_records": [
                review["process_records"] for review in reviews
            ],
        }
        authority = {
            "schema_name": EXPECTED_SCHEMA_NAME,
            "schema_version": EXPECTED_SCHEMA_VERSION,
            "authority_id": EXPECTED_AUTHORITY_ID,
            "task_type": "B_comprehensive_validation",
            "scope": EXPECTED_SCOPE,
            "approval_status": EXPECTED_APPROVAL_STATUS,
            "repository": {
                "path": str(REPO_ROOT.resolve()),
                "git_head_at_authorization": current_head,
                "worktree_policy": (
                    "DIRTY_WORKTREE_ALLOWED_ONLY_WITH_EXACT_CONTENT_HASH_AUTHORITY"
                ),
            },
            "sources": sources,
            "source_set_sha256": test_run["source_set_sha256"],
            "claim_contract": claim_contract,
            "review_evidence": review_evidence,
        }
        authority_encoded = encode_json(authority)
        predicted_authority = predicted_readonly_record(
            AUTHORITY_PATH,
            authority_encoded,
        )
        approval = {
            "schema_name": "intersubmod.release_source_authority.approval",
            "schema_version": EXPECTED_SCHEMA_VERSION,
            "authority_id": EXPECTED_AUTHORITY_ID,
            "approval_status": EXPECTED_APPROVAL_STATUS,
            "public_key": public_key,
            "private_key_retirement": {
                "path": str(PRIVATE_KEY_PATH.resolve()),
                "required_post_signature_mode": "0o0",
                "passphrase_lifecycle": "discarded_after_single_signature",
            },
            "authority_manifest": predicted_authority,
            "review_approvals": review_approvals,
        }
        approval_encoded = encode_json(approval)

        reader.require_paths_still_bound()
        observed_authority = write_exclusive_readonly(
            AUTHORITY_PATH,
            authority_encoded,
        )
        if observed_authority != predicted_authority:
            raise AssemblyError("Authority output differs from prediction")
        observed_approval = write_exclusive_readonly(
            APPROVAL_PATH,
            approval_encoded,
        )
        reader.open(AUTHORITY_PATH, include_mode=True)
        reader.open(APPROVAL_PATH, include_mode=True)
        signing_instruction = signer.format_signing_instruction(observed_approval)
        try:
            parsed = signer.parse_signing_instruction(
                signing_instruction,
                expected_path=APPROVAL_PATH.resolve(strict=False),
            )
        except signer.SignerError as error:
            raise AssemblyError("Signer instruction is not canonical") from error
        if parsed != observed_approval:
            raise AssemblyError("Signer instruction target drifted")
        reader.require_paths_still_bound()
        reader.retain_until_process_exit()

    print(
        json.dumps(
            {
                "authority": observed_authority,
                "approval": observed_approval,
                "git_head": current_head,
                "source_roles": len(sources),
                "source_set_sha256": test_run["source_set_sha256"],
                "pytest": test_run["counts"],
                "reviewer_ids": reviewer_ids,
                "review_session_ids": session_ids,
                "signing_instruction": signing_instruction,
                "pass": True,
                "pass_semantics": (
                    "assembly_complete_signature_and_consumer_verification_pending"
                ),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
