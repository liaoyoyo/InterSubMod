#!/usr/bin/env python3
"""Normalize one final signed external Claude v21 process attestation."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping


TOPIC_ROOT = Path(__file__).resolve().parents[1]
HERE = TOPIC_ROOT / "audit_notes"
EVIDENCE_MODULE_PATH = HERE / "attested_release_evidence_v5.py"


class ReviewValidationError(RuntimeError):
    """Raised when signed external-review evidence is not a clean approval."""


def validate_attested_review(
    *,
    reader: Any,
    release_module: Any,
    reviewer_token: str,
    expected_head: str,
    source_set_sha256: str,
    assembler_record: Mapping[str, Any],
    test_run: Mapping[str, Any],
    evidence_module: Any,
) -> dict[str, Any]:
    evidence = evidence_module
    try:
        verified = evidence.validate_signed_review_attestation(
            reader=reader,
            release_module=release_module,
            reviewer_token=reviewer_token,
            expected_head=expected_head,
            source_set_sha256=source_set_sha256,
            assembler_record=assembler_record,
            test_run=test_run,
        )
    except evidence.EvidenceValidationError as error:
        raise ReviewValidationError(
            f"Reviewer {reviewer_token} signed process evidence failed: {error}"
        ) from error
    review = verified.get("review_payload")
    if (
        not isinstance(review, dict)
        or review.get("verdict") != "APPROVE"
        or review.get("findings_closed") is not True
        or review.get("f1_status") != "RESOLVED_VERIFIED"
        or review.get("f2_status") != "RESOLVED_VERIFIED"
        or review.get("blocking_findings") != []
    ):
        raise ReviewValidationError(
            f"Reviewer {reviewer_token} is not a clean signed approval"
        )
    return verified


__all__ = [
    "EVIDENCE_MODULE_PATH",
    "ReviewValidationError",
    "validate_attested_review",
]
