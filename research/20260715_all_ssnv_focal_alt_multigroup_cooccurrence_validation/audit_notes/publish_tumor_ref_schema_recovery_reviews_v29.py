#!/usr/bin/env python3
"""Validate three independent v29 reviews and publish immutable records."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import runpy
import stat
import sys
from typing import Any, Mapping


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
REVIEW_ROOT = TOPIC_ROOT / "reviews"
VALIDATOR = AUDIT_ROOT / "validate_tumor_ref_schema_recovery_authority_v29.py"
EXTERNAL_ENVELOPE = (
    AUDIT_ROOT
    / "20260724_external_claude_schema_recovery_review_v29.claude_cli.envelope.json"
)
INTERNAL_TRANSPORT_PATHS = {
    "mendel": (
        AUDIT_ROOT
        / "20260724_tumor_ref_schema_recovery_mendel.v29.multi_agent.transport.json"
    ),
    "nash": (
        AUDIT_ROOT
        / "20260724_tumor_ref_schema_recovery_nash.v29.multi_agent.transport.json"
    ),
}


class ReviewPublicationError(RuntimeError):
    """Raised when the formal review publication contract is not exact."""


def reject_constant(value: str) -> None:
    raise ReviewPublicationError(f"Non-finite JSON constant: {value}")


def strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ReviewPublicationError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def strict_loads(data: bytes, label: str) -> Any:
    try:
        return json.loads(
            data.decode("utf-8"),
            object_pairs_hook=strict_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ReviewPublicationError(f"Invalid {label}: {error}") from error


def encode_json(payload: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(payload, ensure_ascii=True, indent=2, sort_keys=True) + "\n"
    ).encode("ascii")


def write_exclusive_readonly(path: Path, data: bytes) -> dict[str, Any]:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
        0o444,
    )
    try:
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise ReviewPublicationError(f"Short review write: {path}")
            view = view[written:]
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or stat.S_IMODE(opened.st_mode) != 0o444
            or opened.st_nlink != 1
            or opened.st_size != len(data)
        ):
            raise ReviewPublicationError(f"Published review identity drift: {path}")
    finally:
        os.close(descriptor)
    parent = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(parent)
    finally:
        os.close(parent)
    return {
        "mode": "0o444",
        "path": str(path),
        "sha256": hashlib.sha256(data).hexdigest(),
        "size_bytes": len(data),
    }


def require_absent(paths: list[Path]) -> None:
    occupied = [str(path) for path in paths if os.path.lexists(path)]
    if occupied:
        raise ReviewPublicationError(f"Review output already exists: {occupied}")


def load_internal_reviews() -> dict[str, Mapping[str, Any]]:
    data = sys.stdin.buffer.read(256 * 1024 + 1)
    if len(data) > 256 * 1024:
        raise ReviewPublicationError("Internal-review transport exceeds size limit")
    payload = strict_loads(data, "internal-review transport")
    if not isinstance(payload, Mapping) or set(payload) != {"mendel", "nash"}:
        raise ReviewPublicationError("Internal-review role set is not exact")
    reviews: dict[str, Mapping[str, Any]] = {}
    for role in ("mendel", "nash"):
        review = payload.get(role)
        if not isinstance(review, Mapping):
            raise ReviewPublicationError(f"Internal review is not an object: {role}")
        reviews[role] = review
    return reviews


def load_external_review() -> Mapping[str, Any]:
    envelope = strict_loads(EXTERNAL_ENVELOPE.read_bytes(), "Claude CLI envelope")
    if not isinstance(envelope, Mapping):
        raise ReviewPublicationError("Claude CLI envelope is not an object")
    structured = envelope.get("structured_output")
    if (
        envelope.get("session_id") != "8e6d6be7-d3bd-400d-ae0a-43a7960deec3"
        or not isinstance(structured, Mapping)
        or structured.get("reviewer_agent_id") != envelope.get("session_id")
    ):
        raise ReviewPublicationError("External review session binding drift")
    return structured


def validate_reviews(
    payloads: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, Any], Mapping[str, Any]]:
    module = runpy.run_path(str(VALIDATOR), run_name="v29_review_publisher")
    module["_validate_ceremony_absence_contract"]()
    recovery_sources = module["_records"](module["RECOVERY_SOURCE_PATHS"])
    legacy_sources = module["_records"](module["LEGACY_SOURCE_PATHS"])
    prior_recovery_chain = module["_validate_prior_recovery_chain"](
        module["_records"](module["PRIOR_RECOVERY_CHAIN_PATHS"])
    )
    rejected_generations = module["_validate_rejected_generations"]()
    terminal_key_rotation, _ = module["_validate_terminal_key_rotation"](
        expected_recovery_private_mode="0o400"
    )
    _, public_key, _ = module["_open_bound"](module["PUBLIC_KEY"])
    reviewer_ids = {
        module["_validate_review"](
            role,
            payloads[role],
            recovery_sources,
            legacy_sources,
            prior_recovery_chain,
            rejected_generations,
            public_key,
            terminal_key_rotation,
        )
        for role in module["REVIEW_PATHS"]
    }
    if len(reviewer_ids) != 3:
        raise ReviewPublicationError("Reviewer transport identities are not distinct")
    return module, public_key


def verify_live(path: Path, expected: bytes) -> None:
    opened = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
    try:
        observed = os.fstat(opened)
        data = os.pread(opened, observed.st_size, 0)
        if (
            not stat.S_ISREG(observed.st_mode)
            or stat.S_IMODE(observed.st_mode) != 0o444
            or observed.st_nlink != 1
            or data != expected
        ):
            raise ReviewPublicationError(f"Post-publication review drift: {path}")
    finally:
        os.close(opened)


def main() -> int:
    internal = load_internal_reviews()
    payloads: dict[str, Mapping[str, Any]] = {
        **internal,
        "external_claude_opus": load_external_review(),
    }
    formal_paths = {
        "mendel": REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_mendel.v29.json",
        "nash": REVIEW_ROOT / "20260724_tumor_ref_schema_recovery_nash.v29.json",
        "external_claude_opus": (
            REVIEW_ROOT
            / "20260724_tumor_ref_schema_recovery_external_claude_opus.v29.json"
        ),
    }
    require_absent([*INTERNAL_TRANSPORT_PATHS.values(), *formal_paths.values()])
    module, public_key = validate_reviews(payloads)
    encoded = {role: encode_json(payload) for role, payload in payloads.items()}

    transport_records = {
        role: write_exclusive_readonly(INTERNAL_TRANSPORT_PATHS[role], encoded[role])
        for role in ("mendel", "nash")
    }
    formal_records = {
        role: write_exclusive_readonly(formal_paths[role], encoded[role])
        for role in ("mendel", "nash", "external_claude_opus")
    }
    for role, path in formal_paths.items():
        verify_live(path, encoded[role])

    print(
        json.dumps(
            {
                "schema_name": "intersubmod.tumor_ref_schema_recovery_review_publication",
                "schema_version": "1.0.0",
                "authority_id": module["AUTHORITY_ID"],
                "formal_reviews": formal_records,
                "internal_review_transports": transport_records,
                "trusted_recovery_public_key": public_key,
                "verdicts": {
                    role: payloads[role]["verdict"] for role in payloads
                },
                "pass": True,
            },
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
