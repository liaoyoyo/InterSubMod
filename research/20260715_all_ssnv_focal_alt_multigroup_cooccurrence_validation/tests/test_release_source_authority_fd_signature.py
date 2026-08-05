from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = TOPIC_ROOT / "scripts" / "release_source_authority.py"
SPEC = importlib.util.spec_from_file_location("release_source_authority_fd_test", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_protected_release_sources_contain_no_superseded_path_tokens() -> None:
    superseded_tokens = (
        "methyl_ssnv_cooccurrence_v6_m2v5_raw_identity_contract_source_locked",
        "cooccurrence_task_contract_preflight.v6_analysis_payload_release_full_runtime.json",
        "stable_primary_artifact_audit.v2_source_authorized_pre_downstream.json",
        "strict_methyl_candidate_confirmation_v2_m2v5",
        "matched_normal_candidate_controls_v2_m2v5",
        "matched_normal_candidate_control_analysis_v2_m2v5",
        "matched_normal_methyl_tag_audit.v2.json",
        "candidate_cn_ccf_annotations_v2_m2v5",
        "all_ssnv_final_report_dataset_v4_m2v5_source_attested",
        "all_ssnv_final_report_v4_m2v5_source_attested",
    )
    offenders = {
        role: sorted(
            token
            for token in superseded_tokens
            if token in path.read_text(encoding="utf-8")
        )
        for role, path in MODULE.EXPECTED_SOURCE_PATHS.items()
    }
    assert {role: tokens for role, tokens in offenders.items() if tokens} == {}


def generate_signature(tmp_path: Path) -> tuple[Path, Path, Path]:
    private_key = tmp_path / "private.pem"
    public_key = tmp_path / "public.pem"
    payload = tmp_path / "payload.json"
    signature = tmp_path / "payload.json.ed25519.sig"
    payload.write_bytes(b'{"scope":"real-fd-integration-test"}\n')
    subprocess.run(
        [
            str(MODULE.OPENSSL_PATH),
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
            str(MODULE.OPENSSL_PATH),
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
            str(MODULE.OPENSSL_PATH),
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
    return payload, public_key, signature


def verify(payload: Path, public_key: Path, signature: Path) -> bool:
    descriptors = [
        os.open(path, os.O_RDONLY | os.O_CLOEXEC)
        for path in (MODULE.OPENSSL_PATH, payload, public_key, signature)
    ]
    try:
        return MODULE.verify_ed25519_signature_fds(
            openssl_fd=descriptors[0],
            data_fd=descriptors[1],
            public_key_fd=descriptors[2],
            signature_fd=descriptors[3],
        )
    finally:
        for descriptor in descriptors:
            os.close(descriptor)


def test_real_regular_fd_ed25519_signature_verifies(tmp_path: Path) -> None:
    payload, public_key, signature = generate_signature(tmp_path)
    assert verify(payload, public_key, signature) is True


def test_real_regular_fd_ed25519_signature_rejects_payload_tamper(
    tmp_path: Path,
) -> None:
    payload, public_key, signature = generate_signature(tmp_path)
    payload.write_bytes(b'{"scope":"tampered-after-signing"}\n')
    assert verify(payload, public_key, signature) is False


def test_signature_verifier_ignores_inherited_loader_and_openssl_environment(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload, public_key, signature = generate_signature(tmp_path)
    monkeypatch.setenv("LD_PRELOAD", "/nonexistent/attacker.so")
    monkeypatch.setenv("LD_LIBRARY_PATH", "/nonexistent/attacker-libraries")
    monkeypatch.setenv("OPENSSL_CONF", "/nonexistent/attacker.cnf")
    monkeypatch.setenv("OPENSSL_MODULES", "/nonexistent/attacker-modules")

    assert verify(payload, public_key, signature) is True
    payload.write_bytes(b'{"scope":"tampered-with-poisoned-parent-env"}\n')
    assert verify(payload, public_key, signature) is False


def fake_verified_review() -> dict[str, object]:
    artifact = {
        "path": "/tmp/evidence.json",
        "size_bytes": 1,
        "sha256": "a" * 64,
        "mode": "0o444",
    }
    payload = {
        "reviewer_id": "45df7edf-681f-4a4c-b69f-d0ea2c3ef527",
        "reviewer_label": "External Reviewer A v19",
        "model": "claude-opus-4-8",
        "verdict": "APPROVE",
        "findings_closed": True,
        "review_scope": "All protected sources and signed process evidence.",
        "reviewed_source_set_sha256": "d" * 64,
        "reviewed_git_head": "a" * 40,
        "reviewed_assembler": artifact,
        "canonical_junit": {
            "artifact": artifact,
            "counts": {
                "tests": 2,
                "failures": 0,
                "errors": 0,
                "skipped": 0,
            },
        },
    }
    return {
        "review_payload": payload,
        "attestation": artifact,
        "signature": artifact,
        "public_key": artifact,
        "session_id": "7b543513-8f45-41af-a1df-02ec98f518d4",
        "process_records": {
            "command_record": artifact,
            "request": artifact,
            "envelope": artifact,
            "stderr": artifact,
        },
    }


def test_attested_release_routes_bind_strict_parity_v25_and_v21_reviews() -> None:
    assert MODULE.AUTHORITY_ASSEMBLER_PATH.name == (
        "assemble_release_source_authority_v15.py"
    )
    assert MODULE.CANONICAL_PYTEST_XML_PATH.name == (
        "pytest_full_pre_authority_v27_git_env_alignment_attested_canonical.xml"
    )
    assert all("v21_strict_command_parity_process_attestation" in path.name for path in (
        MODULE.EXPECTED_EXTERNAL_REVIEW_PATHS
    ))
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "v17_v6_two_stage_output" not in source
    assert "assemble_release_source_authority_v10.py" not in source
    assert "assemble_release_source_authority_v12.py" not in source
    assert "assemble_release_source_authority_v13.py" not in source
    assert "assemble_release_source_authority_v14.py" not in source


def test_approval_is_derived_from_signed_process_attestation() -> None:
    review = fake_verified_review()
    approval = MODULE.approval_from_verified_review(review)
    assert set(approval) == set(MODULE.EXPECTED_REVIEW_APPROVAL_KEYS)
    assert approval["process_attestation"] == review["attestation"]
    assert approval["process_attestation_signature"] == review["signature"]
    assert approval["process_attestation_public_key"] == review["public_key"]
    assert approval["process_attestation_session_id"] == review["session_id"]
    assert approval["process_records"] == review["process_records"]


def test_support_hash_contract_has_exact_roles_and_no_placeholder() -> None:
    assert set(MODULE.EXPECTED_SUPPORT_SHA256) == set(
        MODULE.authority_support_paths()
    )
    assert all(
        len(value) == 64
        and set(value) <= set("0123456789abcdef")
        and value != "0" * 64
        for value in MODULE.EXPECTED_SUPPORT_SHA256.values()
    )


def test_bound_module_execution_uses_supplied_bytes(tmp_path: Path) -> None:
    path = tmp_path / "module.py"
    path.write_text("VALUE = 2\n", encoding="utf-8")
    loaded = MODULE.execute_module_from_bytes(
        "bound_release_test_module",
        path,
        b"VALUE = 1\n",
    )
    assert loaded.VALUE == 1


def test_bound_reader_detects_mode_change_after_open(tmp_path: Path) -> None:
    review = tmp_path / "review.json"
    review.write_text("{}\n", encoding="utf-8")
    review.chmod(0o444)
    with MODULE.BoundArtifactReader() as reader:
        reader.open(review, include_mode=True)
        review.chmod(0o644)
        with pytest.raises(MODULE.SourceAuthorityError, match="binding changed"):
            reader.require_paths_still_bound()
    assert stat.S_IMODE(review.stat().st_mode) == 0o644


def test_bound_reader_tracks_mode_zero_metadata_file(tmp_path: Path) -> None:
    private_key = tmp_path / "retired-private-key.pem"
    private_key.write_text("test-only encrypted key material\n", encoding="utf-8")
    private_key.chmod(0o000)
    try:
        with MODULE.BoundArtifactReader() as reader:
            _, record = reader.open_metadata(private_key, include_mode=True)
            assert record["mode"] == "0o0"
            private_key.chmod(0o600)
            with pytest.raises(MODULE.SourceAuthorityError, match="binding changed"):
                reader.require_paths_still_bound()
    finally:
        private_key.chmod(0o600)


def test_successful_public_validation_does_not_close_bound_descriptors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    protected = tmp_path / "protected.json"
    protected.write_text("{}\n", encoding="utf-8")
    protected.chmod(0o444)
    close_calls: list[int] = []
    real_close = os.close

    def fake_validate(
        reader: object, _required_roles: object
    ) -> dict[str, object]:
        assert isinstance(reader, MODULE.BoundArtifactReader)
        reader.open(protected, include_mode=True)
        return {"pass": True}

    def injected_close(descriptor: int) -> None:
        close_calls.append(descriptor)
        protected.chmod(0o644)
        real_close(descriptor)

    monkeypatch.setattr(MODULE, "_validate_release_source_authority", fake_validate)
    monkeypatch.setattr(os, "close", injected_close)
    result = MODULE.validate_release_source_authority()

    assert result["pass"] is True
    assert result["validation_binding_lease"] == {
        "policy": "BOUND_FILE_DESCRIPTORS_RETAINED_UNTIL_PROCESS_EXIT",
        "descriptor_count": 1,
    }
    assert close_calls == []
    assert stat.S_IMODE(protected.stat().st_mode) == 0o444
