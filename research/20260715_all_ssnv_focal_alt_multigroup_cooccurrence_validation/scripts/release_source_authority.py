#!/usr/bin/env python3
"""Validate the reviewer-approved source authority for the all-sSNV release."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
from types import ModuleType, SimpleNamespace
from typing import Any, Iterable, Mapping
import uuid
import xml.etree.ElementTree as ET


SCRIPT_DIR = Path(__file__).resolve().parent
TOPIC_ROOT = SCRIPT_DIR.parent
REPO_ROOT = TOPIC_ROOT.parents[1]
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
PUBLIC_KEY_SHA256 = "cecb2287cf87f8bd948c7390bd5f7059578966e0fc917de420572ddd82d8f21b"
OPENSSL_PATH = Path("/usr/bin/openssl")
OPENSSL_SHA256 = "38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08"
GIT_PATH = Path("/usr/bin/git")
GIT_SHA256 = "587ef21868c948b883993e23209b86a72a6ddc06aab1545c697ffc31075acd4a"
CLEAN_SUBPROCESS_ENV = {
    "PATH": "/usr/bin:/bin",
    "LC_ALL": "C",
    "LANG": "C",
}
AUTHORITY_SCHEMA_NAME = "intersubmod.release_source_authority"
AUTHORITY_SCHEMA_VERSION = "1.0.0"
AUTHORITY_ID = "20260722_all_ssnv_focal_alt_task_b_release_v7_strict_command_parity"
APPROVAL_STATUS = "APPROVED_FOR_FULL_TASK_B_RUN"
CLAIM_CONTRACT_PATH = TOPIC_ROOT / "claim-contract-v5.md"
CLAIM_CONTRACT_SIZE = 9_144
CLAIM_CONTRACT_SHA256 = (
    "da8778af6bcc9b1e8b2887eb2bcc87eca83d909be32a819ec5eb8f5f12c9f2af"
)
AUTHORITY_ASSEMBLER_PATH = (
    TOPIC_ROOT / "audit_notes" / "assemble_release_source_authority_v15.py"
)
CANONICAL_PYTEST_XML_PATH = (
    TOPIC_ROOT
    / "logs"
    / "pytest_full_pre_authority_v27_git_env_alignment_attested_canonical.xml"
)
ATTESTED_EVIDENCE_PATH = (
    TOPIC_ROOT / "audit_notes" / "attested_release_evidence_v5.py"
)
CANONICAL_TEST_RUNNER_PATH = (
    TOPIC_ROOT / "audit_notes" / "run_attested_canonical_pytest_v5.py"
)
REVIEW_NORMALIZER_PATH = (
    TOPIC_ROOT / "audit_notes" / "normalize_external_source_review_v20.py"
)
EXTERNAL_REVIEW_RUNNER_PATH = (
    TOPIC_ROOT / "audit_notes" / "run_external_claude_review_v22_attested.py"
)
ONE_TIME_SIGNER_PATH = (
    TOPIC_ROOT / "audit_notes" / "run_one_time_ed25519_signer_v6.py"
)
EXPECTED_SUPPORT_SHA256 = {
    "authority_assembler": (
        "5dcec8c2a8b3dd9143a306ca3fe6302ec71a3676cca005ce36e30b2b20421ba5"
    ),
    "attested_evidence_validator": (
        "ba75e6b454863536ac56f08f1630adedf39c7f155be2758d13e7d434ec606fc8"
    ),
    "canonical_test_runner": (
        "10822d820bf362770bc424ddb558833638fc9ca0db0a19602469d9f52397e2fb"
    ),
    "review_normalizer": (
        "8e14b74c133353e73fb61b1b274b991603f63866c017dcab9e5fc09ba276c273"
    ),
    "external_review_runner": (
        "380dd4cc63729d087ecd152e74ddef222b1b8cf470469e3c7492a1fff5cbdad4"
    ),
    "one_time_signer": (
        "c9c32222d4e8ca6daad8f6b3f8a91746fc99dd200d24cf3d77e500234c0f8bee"
    ),
}
ARTIFACT_IDENTITY_KEYS = frozenset({"path", "size_bytes", "sha256", "mode"})
JUNIT_COUNT_KEYS = frozenset({"tests", "failures", "errors", "skipped"})
AUTHORITY_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "authority_id",
        "task_type",
        "scope",
        "approval_status",
        "repository",
        "sources",
        "source_set_sha256",
        "claim_contract",
        "review_evidence",
    }
)
APPROVAL_KEYS = frozenset(
    {
        "schema_name",
        "schema_version",
        "authority_id",
        "approval_status",
        "public_key",
        "private_key_retirement",
        "authority_manifest",
        "review_approvals",
    }
)

EXPECTED_SOURCE_PATHS = {
    "release_source_authority_validator": Path(__file__).resolve(),
    "cooccurrence_analyzer": SCRIPT_DIR / "analyze_methyl_ssnv_cooccurrence.py",
    "cooccurrence_preflight": SCRIPT_DIR / "audit_cooccurrence_task_contract_preflight.py",
    "cooccurrence_release_finalizer": SCRIPT_DIR / "finalize_cooccurrence_release_receipt.py",
    "cooccurrence_release_runner": SCRIPT_DIR / "run_cooccurrence_v6_source_locked.sh",
    "ssnv_cooccurrence_lib": REPO_ROOT / "scripts" / "ssnv_cooccurrence_lib.py",
    "latest_tag_join": SCRIPT_DIR / "latest_tag_join.py",
    "m2_screen_gate": SCRIPT_DIR / "m2_screen_gate.py",
    "independent_m2_auditor": SCRIPT_DIR / "audit_independent_m2_gate_recount.py",
    "primary_artifact_auditor": SCRIPT_DIR / "audit_stable_primary_artifacts.py",
    "strict_producer": SCRIPT_DIR / "run_strict_methyl_candidate_confirmation.py",
    "focal_alt_cluster_lib": (
        TOPIC_ROOT.parent
        / "20260715_single_fp_alt_multicluster_subclone_validation"
        / "scripts"
        / "focal_alt_cluster_lib.py"
    ),
    "completion_runner": SCRIPT_DIR / "run_m2v5_recovered_completion_chain.sh",
    "final_dataset_builder": SCRIPT_DIR / "build_all_ssnv_final_report_dataset.py",
    "final_result_release_finalizer": SCRIPT_DIR / "finalize_task_b_result_release.py",
    "matched_normal_runner": SCRIPT_DIR / "run_matched_normal_candidate_controls.py",
    "matched_normal_analyzer": SCRIPT_DIR / "analyze_matched_normal_candidate_controls.py",
    "cn_ccf_annotator": SCRIPT_DIR / "annotate_candidate_cn_ccf.py",
    "frozen_input_auditor": SCRIPT_DIR / "audit_frozen_input_immutability.py",
    "tumor_ref_source_verifier": (
        SCRIPT_DIR / "verify_retrospective_running_source_identity_v2.py"
    ),
    "report_builder": SCRIPT_DIR / "build_all_ssnv_report_artifact.py",
    "portable_report_deliverer": (
        SCRIPT_DIR / "deliver_portable_report_scrollbar_safe.mjs"
    ),
    "portable_layout_qa": SCRIPT_DIR / "qa_portable_report_layout.py",
}
EXPECTED_EXTERNAL_REVIEW_PATHS = (
    TOPIC_ROOT
    / "audit_notes"
    / "20260722_external_claude_reviewer_A_v21_strict_command_parity_process_attestation.v1.json",
    TOPIC_ROOT
    / "audit_notes"
    / "20260722_external_claude_reviewer_B_v21_strict_command_parity_process_attestation.v1.json",
)
EXPECTED_REVIEW_KEYS = frozenset(
    {
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
    }
)
EXPECTED_REVIEW_APPROVAL_KEYS = frozenset(
    {
        "reviewer_id",
        "reviewer_label",
        "model",
        "verdict",
        "findings_closed",
        "review_scope",
        "reviewed_source_set_sha256",
        "reviewed_git_head",
        "reviewed_assembler",
        "canonical_junit",
        "process_attestation",
        "process_attestation_signature",
        "process_attestation_public_key",
        "process_attestation_session_id",
        "process_records",
    }
)
_PROCESS_LIFETIME_BOUND_READERS: list["BoundArtifactReader"] = []


class SourceAuthorityError(RuntimeError):
    """Raised when the source authority or a protected source has drifted."""


class BoundArtifactReader:
    """Read each trust artifact once and retain its inode through validation."""

    def __init__(self) -> None:
        self._opened: list[tuple[Path, int, os.stat_result]] = []
        self._retained_until_process_exit = False

    def __enter__(self) -> "BoundArtifactReader":
        return self

    def __exit__(self, *_: Any) -> None:
        if not self._retained_until_process_exit:
            self.close()

    @property
    def opened_descriptor_count(self) -> int:
        return len(self._opened)

    def close(self) -> None:
        while self._opened:
            _, descriptor, _ = self._opened.pop()
            os.close(descriptor)

    def retain_until_process_exit(self) -> None:
        """Keep successful validation descriptors open until process teardown."""
        if self._retained_until_process_exit:
            raise SourceAuthorityError("Bound artifact reader was already retained")
        self.require_paths_still_bound()
        _PROCESS_LIFETIME_BOUND_READERS.append(self)
        self._retained_until_process_exit = True

    def open(
        self, path: Path, *, include_mode: bool = False
    ) -> tuple[int, bytes, dict[str, Any]]:
        resolved = path.expanduser().resolve(strict=True)
        flags = os.O_RDONLY | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(resolved, flags)
        first = os.fstat(descriptor)
        if not stat.S_ISREG(first.st_mode):
            os.close(descriptor)
            raise SourceAuthorityError(f"Trust artifact is not a regular file: {resolved}")
        data = b"".join(
            os.pread(descriptor, min(8 * 1024 * 1024, first.st_size - offset), offset)
            for offset in range(0, first.st_size, 8 * 1024 * 1024)
        )
        second = os.fstat(descriptor)
        if (
            len(data) != first.st_size
            or first.st_dev != second.st_dev
            or first.st_ino != second.st_ino
            or first.st_mode != second.st_mode
            or first.st_size != second.st_size
            or first.st_mtime_ns != second.st_mtime_ns
            or first.st_ctime_ns != second.st_ctime_ns
        ):
            os.close(descriptor)
            raise SourceAuthorityError(f"Trust artifact changed while reading: {resolved}")
        live = resolved.stat()
        if (
            live.st_dev != second.st_dev
            or live.st_ino != second.st_ino
            or live.st_mode != second.st_mode
            or live.st_size != second.st_size
            or live.st_mtime_ns != second.st_mtime_ns
            or live.st_ctime_ns != second.st_ctime_ns
        ):
            os.close(descriptor)
            raise SourceAuthorityError(f"Trust artifact path was replaced: {resolved}")
        record: dict[str, Any] = {
            "path": str(resolved),
            "size_bytes": second.st_size,
            "sha256": hashlib.sha256(data).hexdigest(),
        }
        if include_mode:
            record["mode"] = oct(second.st_mode & 0o777)
        self._opened.append((resolved, descriptor, second))
        return descriptor, data, record

    def open_metadata(
        self, path: Path, *, include_mode: bool = False
    ) -> tuple[int, dict[str, Any]]:
        """Bind an unreadable regular file by inode without opening its contents."""
        if not hasattr(os, "O_PATH"):
            raise SourceAuthorityError("O_PATH is required for metadata-only trust binding")
        resolved = path.expanduser().resolve(strict=True)
        flags = os.O_PATH | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(resolved, flags)
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            os.close(descriptor)
            raise SourceAuthorityError(f"Trust artifact is not a regular file: {resolved}")
        live = resolved.stat()
        if (
            live.st_dev != opened.st_dev
            or live.st_ino != opened.st_ino
            or live.st_mode != opened.st_mode
            or live.st_size != opened.st_size
            or live.st_mtime_ns != opened.st_mtime_ns
            or live.st_ctime_ns != opened.st_ctime_ns
        ):
            os.close(descriptor)
            raise SourceAuthorityError(f"Trust artifact path was replaced: {resolved}")
        record: dict[str, Any] = {"path": str(resolved)}
        if include_mode:
            record["mode"] = oct(opened.st_mode & 0o777)
        self._opened.append((resolved, descriptor, opened))
        return descriptor, record

    def require_paths_still_bound(self) -> None:
        for path, descriptor, opened_stat in self._opened:
            descriptor_stat = os.fstat(descriptor)
            live = path.stat()
            if (
                descriptor_stat.st_dev != opened_stat.st_dev
                or descriptor_stat.st_ino != opened_stat.st_ino
                or descriptor_stat.st_mode != opened_stat.st_mode
                or descriptor_stat.st_size != opened_stat.st_size
                or descriptor_stat.st_mtime_ns != opened_stat.st_mtime_ns
                or descriptor_stat.st_ctime_ns != opened_stat.st_ctime_ns
                or live.st_dev != opened_stat.st_dev
                or live.st_ino != opened_stat.st_ino
                or live.st_mode != opened_stat.st_mode
                or live.st_size != opened_stat.st_size
                or live.st_mtime_ns != opened_stat.st_mtime_ns
                or live.st_ctime_ns != opened_stat.st_ctime_ns
            ):
                raise SourceAuthorityError(
                    f"Trust artifact binding changed during validation: {path}"
                )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path, *, include_mode: bool = False) -> dict[str, Any]:
    resolved = path.expanduser().resolve(strict=True)
    result: dict[str, Any] = {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }
    if include_mode:
        result["mode"] = oct(resolved.stat().st_mode & 0o777)
    return result


def _require_exact_artifact(
    record: Any,
    path: Path,
    *,
    label: str,
    include_mode: bool = False,
) -> dict[str, Any]:
    expected = artifact(path, include_mode=include_mode)
    if not isinstance(record, Mapping) or dict(record) != expected:
        raise SourceAuthorityError(f"{label} identity differs from live artifact")
    return expected


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


def verify_ed25519_signature_fds(
    *,
    openssl_fd: int,
    data_fd: int,
    public_key_fd: int,
    signature_fd: int,
) -> bool:
    """Verify an Ed25519 signature from seekable, already-bound file descriptors."""
    verification = subprocess.run(
        [
            f"/proc/self/fd/{openssl_fd}",
            "pkeyutl",
            "-verify",
            "-rawin",
            "-pubin",
            "-inkey",
            f"/proc/self/fd/{public_key_fd}",
            "-sigfile",
            f"/proc/self/fd/{signature_fd}",
            "-in",
            f"/proc/self/fd/{data_fd}",
        ],
        check=False,
        pass_fds=(openssl_fd, data_fd, public_key_fd, signature_fd),
        capture_output=True,
        env=CLEAN_SUBPROCESS_ENV,
    )
    return verification.returncode == 0


def _reject_duplicate_review_keys(
    pairs: list[tuple[str, Any]],
) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for key, value in pairs:
        if key in payload:
            raise SourceAuthorityError(f"Duplicate review JSON key: {key}")
        payload[key] = value
    return payload


def _reject_nonfinite_review_constant(value: str) -> None:
    raise SourceAuthorityError(f"Non-finite review JSON constant: {value}")


def _is_canonical_uuid4(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    try:
        parsed = uuid.UUID(value)
    except ValueError:
        return False
    return parsed.version == 4 and str(parsed) == value


def _is_exact_artifact_identity(value: Any, expected_path: Path) -> bool:
    if not isinstance(value, Mapping) or set(value) != set(ARTIFACT_IDENTITY_KEYS):
        return False
    sha256 = value.get("sha256")
    return (
        value.get("path") == str(expected_path.resolve(strict=False))
        and type(value.get("size_bytes")) is int
        and value["size_bytes"] >= 0
        and isinstance(sha256, str)
        and len(sha256) == 64
        and all(character in "0123456789abcdef" for character in sha256)
        and value.get("mode") == "0o444"
    )


def _parse_junit_counts(data: bytes) -> dict[str, int]:
    try:
        root = ET.fromstring(data)
    except ET.ParseError as exc:
        raise SourceAuthorityError("Canonical pytest XML is malformed") from exc
    suites = [node for node in root.iter() if node.tag.rsplit("}", 1)[-1] == "testsuite"]
    if not suites:
        raise SourceAuthorityError("Canonical pytest XML contains no testsuite")
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
        raise SourceAuthorityError(
            f"JUnit aggregate attributes differ from testcase elements: "
            f"{attributes} != {elements}"
        )
    return attributes


def execute_module_from_bytes(name: str, path: Path, data: bytes) -> ModuleType:
    """Execute only bytes already bound by the caller's descriptor lease."""
    module = ModuleType(name)
    module.__file__ = str(path)
    module.__package__ = ""
    exec(compile(data, str(path), "exec"), module.__dict__)
    return module


def authority_support_paths() -> dict[str, Path]:
    return {
        "canonical_test_runner": CANONICAL_TEST_RUNNER_PATH,
        "authority_assembler": AUTHORITY_ASSEMBLER_PATH,
        "review_normalizer": REVIEW_NORMALIZER_PATH,
        "external_review_runner": EXTERNAL_REVIEW_RUNNER_PATH,
        "attested_evidence_validator": ATTESTED_EVIDENCE_PATH,
        "one_time_signer": ONE_TIME_SIGNER_PATH,
    }


def approval_from_verified_review(
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


def validate_live_review_approvals(
    reader: BoundArtifactReader,
    approvals: Any,
    review_evidence: Any,
    *,
    observed_digest: str,
    expected_head: str,
) -> list[dict[str, Any]]:
    """Re-verify signed tests and both signed Claude process attestations."""
    if not isinstance(approvals, list) or len(approvals) != 2:
        raise SourceAuthorityError("Source authority requires exactly two reviews")
    if not isinstance(review_evidence, Mapping):
        raise SourceAuthorityError("Source authority review evidence is malformed")

    support_records: dict[str, dict[str, Any]] = {}
    support_bytes: dict[str, bytes] = {}
    for role, path in authority_support_paths().items():
        _, data, record = reader.open(path, include_mode=True)
        if (
            record["mode"] != "0o444"
            or record["sha256"] != EXPECTED_SUPPORT_SHA256[role]
        ):
            raise SourceAuthorityError(
                f"Authority support source {role} identity drifted"
            )
        support_records[role] = record
        support_bytes[role] = data

    evidence_module = execute_module_from_bytes(
        "all_ssnv_attested_release_evidence_v5_consumer",
        ATTESTED_EVIDENCE_PATH,
        support_bytes["attested_evidence_validator"],
    )
    normalizer = execute_module_from_bytes(
        "all_ssnv_external_review_normalizer_v18_consumer",
        REVIEW_NORMALIZER_PATH,
        support_bytes["review_normalizer"],
    )
    expected_review_paths = tuple(
        evidence_module.REVIEW_ATTESTATION_SPECS[token]["path"]
        for token in ("A", "B")
    )
    if expected_review_paths != EXPECTED_EXTERNAL_REVIEW_PATHS:
        raise SourceAuthorityError("Signed external-review path contract drifted")
    release_api = SimpleNamespace(
        OPENSSL_PATH=OPENSSL_PATH,
        OPENSSL_SHA256=OPENSSL_SHA256,
        GIT_PATH=GIT_PATH,
        GIT_SHA256=GIT_SHA256,
        EXPECTED_SOURCE_PATHS=EXPECTED_SOURCE_PATHS,
        verify_ed25519_signature_fds=verify_ed25519_signature_fds,
    )
    test_run = evidence_module.validate_signed_test_run(
        reader=reader,
        release_module=release_api,
        expected_head=expected_head,
        expected_support_paths=authority_support_paths(),
    )
    if (
        test_run["source_set_sha256"] != observed_digest
        or test_run["junit"]["path"] != str(CANONICAL_PYTEST_XML_PATH.resolve())
        or test_run["supporting_sources"] != support_records
    ):
        raise SourceAuthorityError(
            "Signed canonical test run differs from live release inputs"
        )

    assembler_record = support_records["authority_assembler"]
    reviews = [
        normalizer.validate_attested_review(
            reader=reader,
            release_module=release_api,
            reviewer_token=token,
            expected_head=expected_head,
            source_set_sha256=observed_digest,
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
    public_keys = [review["public_key"]["sha256"] for review in reviews]
    if (
        len(set(reviewer_ids)) != 2
        or len(set(session_ids)) != 2
        or len(set(public_keys)) != 2
    ):
        raise SourceAuthorityError(
            "External review identities, sessions, and keys must be distinct"
        )

    expected_approvals = [
        approval_from_verified_review(review) for review in reviews
    ]
    if any(
        not isinstance(approval, Mapping)
        or set(approval) != set(EXPECTED_REVIEW_APPROVAL_KEYS)
        for approval in approvals
    ) or approvals != expected_approvals:
        raise SourceAuthorityError(
            "Detached review approvals differ from signed process attestations"
        )
    assembler_entrypoint = review_evidence.get("authority_assembler_entrypoint")
    evidence_module.validate_entrypoint_record(
        assembler_entrypoint,
        expected_script=assembler_record,
        expected_python=test_run["payload"]["runtime"]["python"],
        expected_environment=evidence_module.TEST_ENVIRONMENT,
        label="source authority assembler",
    )
    expected_evidence = {
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
        "external_review_normalizer": support_records["review_normalizer"],
        "attested_evidence_validator": support_records[
            "attested_evidence_validator"
        ],
        "canonical_test_runner": support_records["canonical_test_runner"],
        "external_review_runner": support_records["external_review_runner"],
        "one_time_signer_implementation": support_records["one_time_signer"],
        "repository_head_verifier": test_run["payload"]["runtime"]["git"],
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
    if dict(review_evidence) != expected_evidence:
        raise SourceAuthorityError(
            "Signed authority review evidence differs from live evidence"
        )
    return expected_approvals


def _validate_release_source_authority(
    reader: BoundArtifactReader,
    required_roles: Iterable[str] | None = None,
) -> dict[str, Any]:
    _, authority_bytes, expected_authority_artifact = reader.open(
        AUTHORITY_PATH, include_mode=True
    )
    authority_mode = str(expected_authority_artifact["mode"])
    if authority_mode != "0o444":
        raise SourceAuthorityError(
            f"Source authority must be mode 0444, observed {authority_mode}"
        )
    try:
        payload = json.loads(
            authority_bytes.decode("utf-8"),
            object_pairs_hook=_reject_duplicate_review_keys,
            parse_constant=_reject_nonfinite_review_constant,
        )
    except (
        OSError,
        UnicodeDecodeError,
        json.JSONDecodeError,
        SourceAuthorityError,
    ) as error:
        raise SourceAuthorityError("Source authority is unreadable or malformed") from error
    if not isinstance(payload, Mapping) or set(payload) != set(AUTHORITY_KEYS):
        raise SourceAuthorityError("Source authority root schema is not exact")
    approval_fd, approval_bytes, observed_approval = reader.open(
        APPROVAL_PATH, include_mode=True
    )
    approval_mode = str(observed_approval["mode"])
    if approval_mode != "0o444":
        raise SourceAuthorityError(
            f"Source authority approval must be mode 0444, observed {approval_mode}"
        )
    try:
        approval_payload = json.loads(
            approval_bytes.decode("utf-8"),
            object_pairs_hook=_reject_duplicate_review_keys,
            parse_constant=_reject_nonfinite_review_constant,
        )
    except (
        OSError,
        UnicodeDecodeError,
        json.JSONDecodeError,
        SourceAuthorityError,
    ) as error:
        raise SourceAuthorityError(
            "Source authority approval is unreadable or malformed"
        ) from error
    if not isinstance(approval_payload, Mapping) or set(approval_payload) != set(
        APPROVAL_KEYS
    ):
        raise SourceAuthorityError("Source authority approval schema is not exact")
    public_key_fd, _, observed_public_key = reader.open(
        PUBLIC_KEY_PATH, include_mode=True
    )
    if (
        observed_public_key["sha256"] != PUBLIC_KEY_SHA256
        or observed_public_key["mode"] != "0o444"
    ):
        raise SourceAuthorityError("External Ed25519 public-key trust anchor drifted")
    _, private_key_record = reader.open_metadata(
        PRIVATE_KEY_PATH, include_mode=True
    )
    private_key_mode = str(private_key_record["mode"])
    if private_key_mode != "0o0":
        raise SourceAuthorityError(
            "One-time source-authority private key was not retired to mode 000"
        )
    openssl_fd, _, observed_openssl = reader.open(
        OPENSSL_PATH, include_mode=True
    )
    if observed_openssl["sha256"] != OPENSSL_SHA256:
        raise SourceAuthorityError("Pinned OpenSSL verifier binary drifted")
    git_fd, _, observed_git = reader.open(GIT_PATH, include_mode=True)
    if observed_git["sha256"] != GIT_SHA256:
        raise SourceAuthorityError("Pinned Git binary drifted")
    signature_fd, _, observed_signature = reader.open(
        APPROVAL_SIGNATURE_PATH, include_mode=True
    )
    signature_mode = str(observed_signature["mode"])
    if signature_mode != "0o444":
        raise SourceAuthorityError(
            f"Detached approval signature must be mode 0444, observed {signature_mode}"
        )
    if not verify_ed25519_signature_fds(
        openssl_fd=openssl_fd,
        data_fd=approval_fd,
        public_key_fd=public_key_fd,
        signature_fd=signature_fd,
    ):
        raise SourceAuthorityError("Detached approval Ed25519 signature verification failed")
    fixed_fields = {
        "schema_name": AUTHORITY_SCHEMA_NAME,
        "schema_version": AUTHORITY_SCHEMA_VERSION,
        "authority_id": AUTHORITY_ID,
        "task_type": "B_comprehensive_validation",
        "scope": "all_7_datasets_469849_sites_and_release_chain",
        "approval_status": APPROVAL_STATUS,
    }
    for field, expected in fixed_fields.items():
        if payload.get(field) != expected:
            raise SourceAuthorityError(f"Source authority {field} drifted")
    if (
        approval_payload.get("schema_name")
        != "intersubmod.release_source_authority.approval"
        or approval_payload.get("schema_version") != AUTHORITY_SCHEMA_VERSION
        or approval_payload.get("authority_id") != AUTHORITY_ID
        or approval_payload.get("approval_status") != APPROVAL_STATUS
    ):
        raise SourceAuthorityError("Source authority detached approval drifted")
    if approval_payload.get("public_key") != observed_public_key:
        raise SourceAuthorityError("Detached approval public-key binding drifted")
    if approval_payload.get("private_key_retirement") != {
        "path": str(private_key_record["path"]),
        "required_post_signature_mode": "0o0",
        "passphrase_lifecycle": "discarded_after_single_signature",
    }:
        raise SourceAuthorityError("Detached approval private-key retirement drifted")
    if approval_payload.get("authority_manifest") != expected_authority_artifact:
        raise SourceAuthorityError(
            "Detached approval does not bind the current authority manifest"
        )
    repository = payload.get("repository")
    if (
        not isinstance(repository, Mapping)
        or repository.get("path") != str(REPO_ROOT.resolve())
        or repository.get("worktree_policy")
        != "DIRTY_WORKTREE_ALLOWED_ONLY_WITH_EXACT_CONTENT_HASH_AUTHORITY"
    ):
        raise SourceAuthorityError("Source authority repository policy drifted")
    expected_head = str(repository.get("git_head_at_authorization", ""))
    if len(expected_head) != 40 or any(
        character not in "0123456789abcdef" for character in expected_head
    ):
        raise SourceAuthorityError("Source authority git commit is malformed")
    current_head = subprocess.run(
        [f"/proc/self/fd/{git_fd}", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
        pass_fds=(git_fd,),
        env={
            **CLEAN_SUBPROCESS_ENV,
            "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_CONFIG_GLOBAL": "/dev/null",
        },
    ).stdout.strip()
    if current_head != expected_head:
        raise SourceAuthorityError("Repository HEAD differs from source authorization")
    sources = payload.get("sources")
    if not isinstance(sources, Mapping) or set(sources) != set(EXPECTED_SOURCE_PATHS):
        raise SourceAuthorityError("Source authority role set is not exact")
    verified_sources: dict[str, dict[str, Any]] = {}
    for role, path in EXPECTED_SOURCE_PATHS.items():
        _, _, live_source = reader.open(path, include_mode=True)
        if not isinstance(sources.get(role), Mapping) or dict(
            sources[role]
        ) != live_source:
            raise SourceAuthorityError(
                f"protected source {role} identity differs from live artifact"
            )
        verified_sources[role] = live_source
        if verified_sources[role]["mode"] != "0o444":
            raise SourceAuthorityError(f"Protected source {role} is not mode 0444")
    requested = set(required_roles or EXPECTED_SOURCE_PATHS)
    unknown = requested.difference(EXPECTED_SOURCE_PATHS)
    if unknown:
        raise SourceAuthorityError(f"Unknown required authority roles: {sorted(unknown)}")
    observed_digest = source_set_digest(verified_sources)
    if payload.get("source_set_sha256") != observed_digest:
        raise SourceAuthorityError("Source authority aggregate digest drifted")
    claim_contract = payload.get("claim_contract")
    _, _, observed_claim = reader.open(CLAIM_CONTRACT_PATH, include_mode=True)
    if not isinstance(claim_contract, Mapping) or dict(claim_contract) != observed_claim:
        raise SourceAuthorityError(
            "claim-contract-v5 identity differs from live artifact"
        )
    if (
        observed_claim["size_bytes"] != CLAIM_CONTRACT_SIZE
        or observed_claim["sha256"] != CLAIM_CONTRACT_SHA256
        or observed_claim["mode"] != "0o444"
    ):
        raise SourceAuthorityError("Canonical claim-contract-v5 identity drifted")
    approvals = approval_payload.get("review_approvals")
    verified_approvals = validate_live_review_approvals(
        reader,
        approvals,
        payload.get("review_evidence"),
        observed_digest=observed_digest,
        expected_head=expected_head,
    )
    reader.require_paths_still_bound()
    return {
        "authority_id": AUTHORITY_ID,
        "authority_manifest": expected_authority_artifact,
        "detached_approval": observed_approval,
        "detached_approval_signature": observed_signature,
        "external_public_key": observed_public_key,
        "retired_private_key": {
            "path": str(private_key_record["path"]),
            "mode": private_key_mode,
            "passphrase_lifecycle": "discarded_after_single_signature",
        },
        "signature_verifier": observed_openssl,
        "repository_head_verifier": observed_git,
        "authorized_git_head": expected_head,
        "source_set_sha256": observed_digest,
        "required_roles_verified": sorted(verified_sources),
        "all_source_roles_verified": sorted(verified_sources),
        "claim_contract": observed_claim,
        "review_approvals": verified_approvals,
        "pass": True,
        "pass_semantics": (
            "local_reproducible_process_attestation_for_signed_tests_signed_"
            "claude_reviews_and_exact_sources_with_process_lifetime_fd_lease"
        ),
    }


def validate_release_source_authority(
    required_roles: Iterable[str] | None = None,
) -> dict[str, Any]:
    with BoundArtifactReader() as reader:
        result = _validate_release_source_authority(reader, required_roles)
        result["validation_binding_lease"] = {
            "policy": "BOUND_FILE_DESCRIPTORS_RETAINED_UNTIL_PROCESS_EXIT",
            "descriptor_count": reader.opened_descriptor_count,
        }
        reader.retain_until_process_exit()
        return result


__all__ = [
    "APPROVAL_STATUS",
    "APPROVAL_PATH",
    "APPROVAL_SIGNATURE_PATH",
    "AUTHORITY_ASSEMBLER_PATH",
    "AUTHORITY_ID",
    "AUTHORITY_PATH",
    "CLAIM_CONTRACT_PATH",
    "CLAIM_CONTRACT_SHA256",
    "CLAIM_CONTRACT_SIZE",
    "CANONICAL_PYTEST_XML_PATH",
    "EXPECTED_SOURCE_PATHS",
    "GIT_PATH",
    "GIT_SHA256",
    "OPENSSL_PATH",
    "OPENSSL_SHA256",
    "PUBLIC_KEY_PATH",
    "PUBLIC_KEY_SHA256",
    "PRIVATE_KEY_PATH",
    "SourceAuthorityError",
    "artifact",
    "sha256",
    "source_set_digest",
    "verify_ed25519_signature_fds",
    "validate_release_source_authority",
]
