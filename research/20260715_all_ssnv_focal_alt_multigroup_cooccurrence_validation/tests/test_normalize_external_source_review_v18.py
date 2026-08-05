from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    TOPIC_ROOT / "audit_notes" / "normalize_external_source_review_v20.py"
)


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "normalize_external_source_review_v20_test",
        MODULE_PATH,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def module() -> ModuleType:
    return load_module()


class EvidenceFailure:
    class EvidenceValidationError(RuntimeError):
        pass

    @staticmethod
    def validate_signed_review_attestation(**_kwargs: object) -> object:
        raise EvidenceFailure.EvidenceValidationError("unsigned review")


class EvidenceNonclean:
    class EvidenceValidationError(RuntimeError):
        pass

    @staticmethod
    def validate_signed_review_attestation(**_kwargs: object) -> object:
        return {
            "review_payload": {
                "verdict": "REQUEST_CHANGES",
                "findings_closed": False,
                "f1_status": "OPEN",
                "f2_status": "OPEN",
                "blocking_findings": ["still open"],
            }
        }


def call(module: ModuleType, evidence: object) -> object:
    return module.validate_attested_review(
        reader=object(),
        release_module=object(),
        reviewer_token="A",
        expected_head="a" * 40,
        source_set_sha256="b" * 64,
        assembler_record={},
        test_run={},
        evidence_module=evidence,
    )


def test_unsigned_or_invalid_process_evidence_is_rejected(
    module: ModuleType,
) -> None:
    with pytest.raises(module.ReviewValidationError, match="signed process"):
        call(module, EvidenceFailure)


def test_signed_but_nonclean_review_is_rejected(module: ModuleType) -> None:
    with pytest.raises(module.ReviewValidationError, match="not a clean"):
        call(module, EvidenceNonclean)


def test_v18_normalizer_has_no_v17_or_raw_json_fallback(
    module: ModuleType,
) -> None:
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "validate_signed_review_attestation" in source
    assert "parse_single_json_object" not in source
    assert "v17" not in source
