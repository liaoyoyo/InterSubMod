#!/usr/bin/env python3
"""Publish the reviewed transition-relation recovery v3 authority."""

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


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
VALIDATOR = AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v3.py"
EXPECTED_VALIDATOR_SHA256 = (
    "f216d4c1d18fc364b3393308c015ce73d06fc335774bf4d5d5e8a75b8c765e9b"
)

REVIEWER_IDS = {
    "mendel": "019f8af7-757a-7c22-905c-b416760f7e34",
    "nash": "019f8af7-a9bd-7fd0-b6e5-c88ec1fe41b7",
    "external_claude_opus": "c1bc858d-7bfb-4289-a286-91b4e7d05f9a",
}


def stable_json(payload: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True) + "\n"
    ).encode("ascii")


def basic_artifact(path: Path, *, expected_mode: str | None = None) -> dict[str, Any]:
    if path.is_symlink() or path.resolve(strict=True) != path:
        raise RuntimeError(f"Artifact path is not canonical: {path}")
    before = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
        raise RuntimeError(f"Artifact is not a single-link regular file: {path}")
    data = path.read_bytes()
    after = os.stat(path, follow_symlinks=False)
    fields = (
        "st_dev",
        "st_ino",
        "st_mode",
        "st_nlink",
        "st_size",
        "st_mtime_ns",
        "st_ctime_ns",
    )
    if any(getattr(before, field) != getattr(after, field) for field in fields):
        raise RuntimeError(f"Artifact changed while reading: {path}")
    record = {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(stat.S_IMODE(after.st_mode)),
    }
    if expected_mode is not None and record["mode"] != expected_mode:
        raise RuntimeError(f"Artifact mode drift: {path}")
    return record


def predicted_artifact(path: Path, data: bytes) -> dict[str, Any]:
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": "0o444",
    }


def load_validator() -> ModuleType:
    record = basic_artifact(VALIDATOR, expected_mode="0o444")
    if record["sha256"] != EXPECTED_VALIDATOR_SHA256:
        raise RuntimeError("Frozen recovery validator SHA-256 drift")
    data = VALIDATOR.read_bytes()
    module = ModuleType("bound_tumor_ref_schema_recovery_authority_builder_v2")
    module.__file__ = str(VALIDATOR)
    module.__package__ = ""
    exec(compile(data, str(VALIDATOR), "exec", dont_inherit=True), module.__dict__)
    return module


def require_private_pre_signature(path: Path) -> tuple[int, os.stat_result]:
    if path.is_symlink() or path.resolve(strict=True) != path:
        raise RuntimeError("Recovery private-key path is not canonical")
    before = os.stat(path, follow_symlinks=False)
    if (
        not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or stat.S_IMODE(before.st_mode) != 0o400
    ):
        raise RuntimeError("Recovery private key is not single-link mode 0400")
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    opened = os.fstat(descriptor)
    if (
        opened.st_dev != before.st_dev
        or opened.st_ino != before.st_ino
        or opened.st_mode != before.st_mode
        or opened.st_nlink != before.st_nlink
        or opened.st_size != before.st_size
    ):
        os.close(descriptor)
        raise RuntimeError("Recovery private key changed while binding")
    return descriptor, opened


def require_slots_absent(module: ModuleType) -> None:
    paths = [module.AUTHORITY, module.AUTHORITY_SIGNATURE]
    occupied = [str(path) for path in paths if os.path.lexists(path)]
    if occupied:
        raise RuntimeError(f"Recovery authority output slots are occupied: {occupied}")


def run_probe(module: ModuleType) -> dict[str, Any]:
    result = subprocess.run(
        module.READONLY_PROBE_COMMAND,
        cwd=REPO_ROOT,
        capture_output=True,
        check=False,
        timeout=900,
    )
    if result.returncode != 0 or result.stderr:
        raise RuntimeError(f"Recovery read-only probe failed: exit={result.returncode}")
    payload = json.loads(result.stdout)
    if (
        payload.get("pass") is not True
        or payload.get("no_output_writes") is not True
        or payload.get("forbidden_output_slots_checked") != 30
        or payload.get("prior_inputs_verified") is not True
        or payload.get("prior_recorded_verifier_recheck", {}).get("pass") is not True
        or payload.get("static_contract", {}).get("pass") is not True
        or payload.get("regression_tests", {}).get("summary") != "25 passed"
        or payload.get("review_evidence_state") != "all_present"
    ):
        raise RuntimeError("Recovery read-only probe contract drift")
    return payload


def build_payloads(module: ModuleType) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    recovery_sources = module._records(module.RECOVERY_SOURCE_PATHS)
    legacy_sources = module._records(module.LEGACY_SOURCE_PATHS)
    original_chain = module._records(module.ORIGINAL_CHAIN_PATHS)
    module._validate_pinned_records(legacy_sources, original_chain)
    prior_recovery_records = module._records(module.PRIOR_RECOVERY_CHAIN_PATHS)
    prior_recovery_chain = module._validate_prior_recovery_chain(
        prior_recovery_records
    )
    public_key = basic_artifact(module.PUBLIC_KEY, expected_mode="0o444")
    if public_key["sha256"] != module.EXPECTED_PUBLIC_KEY_SHA256:
        raise RuntimeError("Recovery public-key SHA-256 drift")

    review_records: dict[str, dict[str, Any]] = {}
    for role, path in module.REVIEW_PATHS.items():
        data = path.read_bytes()
        payload = module._load_json(data, f"pre-existing recovery review {role}")
        if payload.get("reviewer_agent_id") != REVIEWER_IDS[role]:
            raise RuntimeError(f"Recovery reviewer identity drift: {role}")
        module._validate_review(
            role,
            payload,
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            public_key,
        )
        review_records[role] = basic_artifact(path, expected_mode="0o444")
    authority = {
        "schema_name": module.SCHEMA_NAME,
        "schema_version": module.SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "authority_id": module.AUTHORITY_ID,
        "task_type": "B_comprehensive_validation",
        "status": "APPROVED_FOR_TRANSITION_RELATION_RECOVERY_ONLY",
        "scope": module.RECOVERY_SCOPE,
        "public_key": public_key,
        "retired_private_key": {"path": str(module.PRIVATE_KEY), "mode": "0o0"},
        "original_signed_chain": original_chain,
        "prior_recovery_chain": prior_recovery_chain,
        "legacy_sources": legacy_sources,
        "recovery_sources": recovery_sources,
        "legacy_commands": module.LEGACY_COMMANDS,
        "recovery_commands": module.RECOVERY_COMMANDS,
        "schemas": module.SCHEMAS,
        "receipt_paths": module.RECEIPT_PATHS,
        "review_artifacts": review_records,
        "checks": module.AUTHORITY_CHECKS,
        "pass": True,
        "pass_semantics": (
            "authorizes only transition-context-aware treatment of the exact historical "
            "eight-key observation while retaining live binding of the current nine-key record; "
            "does not alter scientific payload, canonical receipt bytes, downstream "
            "gates, or biological claim ceiling"
        ),
    }
    return review_records, authority


def memfd(name: str, data: bytes) -> int:
    descriptor = os.memfd_create(name, flags=os.MFD_CLOEXEC)
    offset = 0
    while offset < len(data):
        written = os.write(descriptor, data[offset:])
        if written <= 0:
            os.close(descriptor)
            raise RuntimeError(f"Short memfd write: {name}")
        offset += written
    os.lseek(descriptor, 0, os.SEEK_SET)
    return descriptor


def sign_and_verify(
    module: ModuleType, authority_data: bytes, private_fd: int
) -> bytes:
    authority_fd = memfd("schema_recovery_authority", authority_data)
    signature_fd = os.memfd_create("schema_recovery_signature", flags=os.MFD_CLOEXEC)
    public_flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        public_flags |= os.O_NOFOLLOW
    public_fd = os.open(module.PUBLIC_KEY, public_flags)
    environment = {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C"}
    try:
        sign = subprocess.run(
            [
                str(module.OPENSSL),
                "pkeyutl",
                "-sign",
                "-rawin",
                "-inkey",
                f"/proc/self/fd/{private_fd}",
                "-in",
                f"/proc/self/fd/{authority_fd}",
                "-out",
                f"/proc/self/fd/{signature_fd}",
            ],
            env=environment,
            pass_fds=(private_fd, authority_fd, signature_fd),
            capture_output=True,
            check=False,
            timeout=60,
        )
        if sign.returncode != 0 or sign.stderr:
            raise RuntimeError(f"Recovery authority signing failed: exit={sign.returncode}")
        signature = os.pread(signature_fd, 128, 0)
        if len(signature) != 64:
            raise RuntimeError("Recovery authority signature is not 64 bytes")
        verify = subprocess.run(
            [
                str(module.OPENSSL),
                "pkeyutl",
                "-verify",
                "-rawin",
                "-pubin",
                "-inkey",
                f"/proc/self/fd/{public_fd}",
                "-in",
                f"/proc/self/fd/{authority_fd}",
                "-sigfile",
                f"/proc/self/fd/{signature_fd}",
            ],
            env=environment,
            pass_fds=(public_fd, authority_fd, signature_fd),
            capture_output=True,
            check=False,
            timeout=60,
        )
        if verify.returncode != 0 or verify.stderr:
            raise RuntimeError("In-memory recovery authority signature verification failed")
        return signature
    finally:
        os.close(public_fd)
        os.close(signature_fd)
        os.close(authority_fd)


def write_exclusive(path: Path, data: bytes) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o400)
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise RuntimeError(f"Short publication write: {path}")
            offset += written
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    parent_fd = os.open(path.parent, os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY)
    try:
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)
    if path.read_bytes() != data or basic_artifact(path)["mode"] != "0o444":
        raise RuntimeError(f"Published artifact failed reopen verification: {path}")


def preflight(module: ModuleType) -> tuple[dict[str, dict[str, Any]], bytes, dict[str, Any]]:
    require_slots_absent(module)
    private_fd, _ = require_private_pre_signature(module.PRIVATE_KEY)
    os.close(private_fd)
    probe = run_probe(module)
    require_slots_absent(module)
    reviews, authority = build_payloads(module)
    authority_data = stable_json(authority)
    return reviews, authority_data, probe


def publish(module: ModuleType) -> dict[str, Any]:
    reviews, authority_data, probe = preflight(module)
    private_fd, private_stat = require_private_pre_signature(module.PRIVATE_KEY)
    try:
        signature = sign_and_verify(module, authority_data, private_fd)
        require_slots_absent(module)
        for role, expected in reviews.items():
            if basic_artifact(
                module.REVIEW_PATHS[role], expected_mode="0o444"
            ) != expected:
                raise RuntimeError(f"Recovery review changed before publication: {role}")
        write_exclusive(module.AUTHORITY, authority_data)
        write_exclusive(module.AUTHORITY_SIGNATURE, signature)
        live_private = os.stat(module.PRIVATE_KEY, follow_symlinks=False)
        if (
            live_private.st_dev != private_stat.st_dev
            or live_private.st_ino != private_stat.st_ino
            or live_private.st_mode != private_stat.st_mode
            or live_private.st_nlink != private_stat.st_nlink
        ):
            raise RuntimeError("Recovery private-key path changed before retirement")
        os.fchmod(private_fd, 0o000)
        os.fsync(private_fd)
        retired = os.fstat(private_fd)
        if stat.S_IMODE(retired.st_mode) != 0:
            raise RuntimeError("Recovery private-key retirement failed")
    finally:
        os.close(private_fd)

    validation = module.validate_recovery_authority(
        expected_runtime_role="continuation_verifier"
    )
    if validation.get("pass") is not True:
        raise RuntimeError("Published recovery authority failed runtime validation")
    return {
        "schema_name": "intersubmod.tumor_ref_schema_recovery_authority_build",
        "schema_version": "1.0.0",
        "authority": basic_artifact(module.AUTHORITY, expected_mode="0o444"),
        "signature": basic_artifact(
            module.AUTHORITY_SIGNATURE, expected_mode="0o444"
        ),
        "reviews": {
            role: basic_artifact(path, expected_mode="0o444")
            for role, path in module.REVIEW_PATHS.items()
        },
        "private_key_mode": oct(
            stat.S_IMODE(os.stat(module.PRIVATE_KEY, follow_symlinks=False).st_mode)
        ),
        "probe_pass": probe["pass"],
        "runtime_validation": validation,
        "pass": True,
    }


def main() -> int:
    module = load_validator()
    if sys.argv[1:] == ["--preflight"]:
        reviews, authority_data, probe = preflight(module)
        result = {
            "mode": "preflight",
            "review_count": len(reviews),
            "predicted_authority_size_bytes": len(authority_data),
            "predicted_authority_sha256": hashlib.sha256(authority_data).hexdigest(),
            "probe_pass": probe["pass"],
            "no_output_writes": probe["no_output_writes"],
            "pass": True,
        }
    elif sys.argv[1:] == ["--publish-after-three-approvals"]:
        result = publish(module)
    else:
        raise RuntimeError(
            "Usage: builder --preflight | --publish-after-three-approvals"
        )
    print(json.dumps(result, ensure_ascii=True, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
