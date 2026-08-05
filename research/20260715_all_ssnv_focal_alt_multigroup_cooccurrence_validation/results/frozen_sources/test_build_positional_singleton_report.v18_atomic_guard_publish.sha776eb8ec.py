from __future__ import annotations

import ast
import csv
import gzip
import importlib.util
import inspect
import json
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "build_positional_singleton_report.py"
)
SPEC = importlib.util.spec_from_file_location("singleton_report", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_select_cases_includes_dataset_maxima_and_required_caveat_cases() -> None:
    rows = []
    for index, dataset in enumerate(MODULE.EXPECTED_DATASETS):
        if dataset == "COLO829":
            continue
        maximum = {
            "dataset": dataset,
            "chrom": "chr1",
            "pos": 200 + index,
            "ref": "G",
            "alt": "A",
            "core_read_n": 30,
            "methyl_group_count": 2,
        }
        if dataset == "HCC1395_DORADO":
            maximum = {
                "dataset": dataset,
                "chrom": "chr22",
                "pos": 47_518_662,
                "ref": "A",
                "alt": "G",
                "core_read_n": 30,
                "methyl_group_count": 2,
            }
        rows.extend(
            [
                {
                    "dataset": dataset,
                    "chrom": "chr1",
                    "pos": 100 + index,
                    "ref": "C",
                    "alt": "T",
                    "core_read_n": 20,
                    "methyl_group_count": 2,
                },
                maximum,
            ]
        )
    rows.extend(
        [
            {
                "dataset": "HCC1395_DORADO",
                "chrom": "chr1",
                "pos": 20_467_811,
                "ref": "G",
                "alt": "C",
                "core_read_n": 18,
                "methyl_group_count": 2,
            },
            {
                "dataset": "HCC1937",
                "chrom": "chr5",
                "pos": 43_849_776,
                "ref": "T",
                "alt": "C",
                "core_read_n": 19,
                "methyl_group_count": 2,
            },
        ]
    )
    selected = MODULE.select_cases(rows)
    assert len(selected) == 8
    selected_keys = {
        (
            row["dataset"],
            row["chrom"],
            row["pos"],
            row["ref"],
            row["alt"],
        )
        for row in selected
    }
    assert set(MODULE.REQUIRED_VISUAL_CASES).issubset(selected_keys)
    assert {row["dataset"] for row in selected} == set(
        MODULE.EXPECTED_DATASETS
    ) - {"COLO829"}


def test_select_cases_rejects_missing_positive_dataset() -> None:
    rows = [
        {
            "dataset": "H2009",
            "chrom": "chr1",
            "pos": 100,
            "ref": "C",
            "alt": "T",
            "core_read_n": 20,
            "methyl_group_count": 2,
        }
    ]
    with pytest.raises(
        MODULE.ReportBuildError, match="Expected illustrated datasets"
    ):
        MODULE.select_cases(rows)


def test_select_cases_rejects_missing_required_caveat_case() -> None:
    rows = []
    for index, dataset in enumerate(MODULE.EXPECTED_POSITIVE_DATASETS):
        rows.append(
            {
                "dataset": dataset,
                "chrom": "chr1",
                "pos": 100 + index,
                "ref": "C",
                "alt": "T",
                "core_read_n": 20,
                "methyl_group_count": 2,
            }
        )
    with pytest.raises(MODULE.ReportBuildError, match="Required illustrated"):
        MODULE.select_cases(rows)


def make_phase_assignment(
    key: tuple[str, str, int, str, str],
    ps_values: list[int | None],
    *,
    include_latest_ps: bool = True,
) -> dict[str, object]:
    reads = []
    for index, value in enumerate(ps_values):
        read = {"read_id": str(index), "label": "1-1" if index < 2 else "1-2"}
        if include_latest_ps:
            read["latest_ps"] = value
        reads.append(read)
    return {
        "dataset": key[0],
        "chrom": key[1],
        "pos": key[2],
        "core_reads": reads,
    }


def test_summarize_phase_sets_counts_missing_and_multiple_ps() -> None:
    first = ("HCC1937", "chr5", 43_849_776, "T", "C")
    second = ("H2009", "chr1", 10, "A", "G")
    candidates = [
        {
            "dataset": key[0],
            "chrom": key[1],
            "pos": key[2],
            "ref": key[3],
            "alt": key[4],
            "core_read_n": count,
        }
        for key, count in ((first, 4), (second, 2))
    ]
    assignments = {
        first: make_phase_assignment(first, [None, 10, 20, None]),
        second: make_phase_assignment(second, [30, 30]),
    }
    audit = MODULE.summarize_phase_sets(candidates, assignments)
    assert audit["site_count"] == 2
    assert audit["total_core_reads"] == 6
    assert audit["total_missing_ps"] == 2
    assert audit["sites_with_missing_ps"] == 1
    assert audit["sites_with_multiple_nonmissing_ps"] == 1
    assert audit["sites_with_missing_and_multiple_nonmissing_ps"] == 1
    assert audit["m2_phase_set_axis_evaluated"] is False


def test_summarize_phase_sets_rejects_core_count_drift() -> None:
    key = ("H2009", "chr1", 10, "A", "G")
    candidate = {
        "dataset": key[0],
        "chrom": key[1],
        "pos": key[2],
        "ref": key[3],
        "alt": key[4],
        "core_read_n": 3,
    }
    with pytest.raises(MODULE.ReportBuildError, match="core-read count drift"):
        MODULE.summarize_phase_sets(
            [candidate], {key: make_phase_assignment(key, [10, 10])}
        )


def test_summarize_phase_sets_rejects_missing_latest_ps() -> None:
    key = ("H2009", "chr1", 10, "A", "G")
    candidate = {
        "dataset": key[0],
        "chrom": key[1],
        "pos": key[2],
        "ref": key[3],
        "alt": key[4],
        "core_read_n": 1,
    }
    with pytest.raises(MODULE.ReportBuildError, match="lacks latest_ps"):
        MODULE.summarize_phase_sets(
            [candidate],
            {
                key: make_phase_assignment(
                    key, [None], include_latest_ps=False
                )
            },
        )


def test_validate_phase_set_audit_rejects_expected_count_drift() -> None:
    with pytest.raises(MODULE.ReportBuildError, match="audit drift"):
        MODULE.validate_phase_set_audit(
            {
                **MODULE.EXPECTED_PHASE_SET_AUDIT,
                "total_missing_ps": (
                    MODULE.EXPECTED_PHASE_SET_AUDIT["total_missing_ps"] - 1
                ),
                "sites": [],
            }
        )


def make_expected_phase_set_audit() -> dict[str, object]:
    key = MODULE.EXPECTED_PHASE_SET_CASE["key"]
    return {
        **MODULE.EXPECTED_PHASE_SET_AUDIT,
        "sites": [
            {
                "dataset": key[0],
                "chrom": key[1],
                "pos": key[2],
                "ref": key[3],
                "alt": key[4],
                **{
                    field: MODULE.EXPECTED_PHASE_SET_CASE[field]
                    for field in (
                        "core_read_n",
                        "ps_missing_n",
                        "ps_nonmissing_n",
                        "ps_distinct_nonmissing_n",
                        "ps_values",
                    )
                },
            }
        ],
    }


def test_validate_phase_set_audit_accepts_expected_hcc1937_case() -> None:
    MODULE.validate_phase_set_audit(make_expected_phase_set_audit())


def test_validate_phase_set_audit_rejects_hcc1937_case_drift() -> None:
    audit = make_expected_phase_set_audit()
    audit["sites"][0]["ps_missing_n"] -= 1
    with pytest.raises(MODULE.ReportBuildError, match="case drift"):
        MODULE.validate_phase_set_audit(audit)


def test_html_table_escapes_source_values() -> None:
    output = MODULE.html_table(["A&B"], [["<script>", '"quoted"']])
    assert '<th scope="col">' in output
    assert '<th scope="row">&lt;script&gt;</th>' in output
    assert "A&amp;B" in output
    assert "&lt;script&gt;" in output
    assert "&quot;quoted&quot;" in output
    assert "<script>" not in output


def test_report_rename_no_replace_preserves_target(tmp_path: Path) -> None:
    source = tmp_path / "source"
    target = tmp_path / "target"
    source.mkdir()
    target.mkdir()
    marker = target / "marker.txt"
    marker.write_text("first publisher\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="No-replace"):
        MODULE.rename_no_replace(source, target)
    assert marker.read_text(encoding="utf-8") == "first publisher\n"


def make_signed_receipt(tmp_path: Path) -> tuple[Path, Path, Path]:
    private_key = tmp_path / "private.pem"
    public_key = tmp_path / "public.pem"
    receipt = tmp_path / "receipt.json"
    signature = tmp_path / "receipt.json.ed25519.sig"
    subprocess.run(
        [
            str(MODULE.OPENSSL),
            "genpkey",
            "-algorithm",
            "ED25519",
            "-out",
            str(private_key),
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            str(MODULE.OPENSSL),
            "pkey",
            "-in",
            str(private_key),
            "-pubout",
            "-out",
            str(public_key),
        ],
        check=True,
        capture_output=True,
    )
    receipt.write_text('{"pass":true}\n', encoding="utf-8")
    subprocess.run(
        [
            str(MODULE.OPENSSL),
            "pkeyutl",
            "-sign",
            "-rawin",
            "-inkey",
            str(private_key),
            "-in",
            str(receipt),
            "-out",
            str(signature),
        ],
        check=True,
        capture_output=True,
    )
    for path in (receipt, signature, public_key):
        path.chmod(0o444)
    return receipt, signature, public_key


def write_passing_pytest_xml(path: Path, tests: int | None = None) -> None:
    total = tests or MODULE.MINIMUM_CANONICAL_TEST_COUNT
    required = sorted(MODULE.REQUIRED_SUPPLEMENTAL_TESTS)
    if total < len(required):
        raise ValueError("Test total is smaller than the required named set")
    names = required + [
        f"test_filler_{index}" for index in range(total - len(required))
    ]
    path.write_text(
        (
            f'<testsuites><testsuite tests="{total}" failures="0" '
            'errors="0" skipped="0">'
            + "".join(
                f'<testcase classname="current" name="{name}" />'
                for name in names
            )
            + "</testsuite></testsuites>"
        ),
        encoding="utf-8",
    )


def make_signed_test_evidence_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> dict[str, Path]:
    paths = {
        "xml": tmp_path / "canonical.xml",
        "builder": tmp_path / "report_builder.py",
        "tests": tmp_path / "test_report_builder.py",
        "manifest": tmp_path / "test_evidence.json",
        "signature": tmp_path / "test_evidence.json.ed25519.sig",
        "private_key": tmp_path / "test_evidence.private.pem",
        "public_key": tmp_path / "test_evidence.public.pem",
    }
    write_passing_pytest_xml(paths["xml"])
    paths["builder"].write_text("BUILDER = 'current'\n", encoding="utf-8")
    paths["tests"].write_text("TESTS = 'current'\n", encoding="utf-8")
    subprocess.run(
        [
            str(MODULE.OPENSSL),
            "genpkey",
            "-algorithm",
            "ED25519",
            "-out",
            str(paths["private_key"]),
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            str(MODULE.OPENSSL),
            "pkey",
            "-in",
            str(paths["private_key"]),
            "-pubout",
            "-out",
            str(paths["public_key"]),
        ],
        check=True,
        capture_output=True,
    )
    for key in ("xml", "builder", "tests", "public_key"):
        paths[key].chmod(0o444)
    paths["private_key"].chmod(0o400)
    for attribute, key in (
        ("CANONICAL_PYTEST_XML", "xml"),
        ("TEST_EVIDENCE_MANIFEST", "manifest"),
        ("TEST_EVIDENCE_SIGNATURE", "signature"),
        ("TEST_EVIDENCE_PUBLIC_KEY", "public_key"),
        ("TEST_EVIDENCE_PRIVATE_KEY", "private_key"),
        ("TEST_EVIDENCE_REPORT_BUILDER_PATH", "builder"),
        ("TEST_EVIDENCE_TEST_SOURCE_PATH", "tests"),
    ):
        monkeypatch.setattr(MODULE, attribute, paths[key])
    monkeypatch.setattr(
        MODULE,
        "TEST_EVIDENCE_PUBLIC_KEY_SHA256",
        MODULE.sha256(paths["public_key"]),
    )
    authority = MODULE.FINAL_RELEASE.SOURCE_AUTHORITY
    artifacts = {
        "canonical_pytest_xml": authority.artifact(
            paths["xml"], include_mode=True
        ),
        "report_builder": authority.artifact(
            paths["builder"], include_mode=True
        ),
        "report_builder_tests": authority.artifact(
            paths["tests"], include_mode=True
        ),
    }
    public_key = authority.artifact(paths["public_key"], include_mode=True)
    payload = {
        "schema_name": "intersubmod.supplemental_report_test_evidence",
        "schema_version": "1.0.0",
        "authority_id": MODULE.TEST_EVIDENCE_AUTHORITY_ID,
        "task_type": "B_comprehensive_validation",
        "scope": "positional_singleton_supplemental_report_release",
        "artifacts": artifacts,
        "canonical_test_summary": MODULE.validate_pytest_xml(paths["xml"]),
        "required_test_names": sorted(MODULE.REQUIRED_SUPPLEMENTAL_TESTS),
        "signature_contract": {
            "algorithm": "Ed25519",
            "public_key": public_key,
            "signed_artifact": str(paths["manifest"].resolve()),
            "signature": str(paths["signature"].resolve()),
            "private_key_lifecycle": (
                "encrypted_one_time_key_chmod_000_after_signing"
            ),
            "private_key_path": str(paths["private_key"].resolve()),
        },
        "checks": {
            "canonical_pytest_pass": True,
            "one_time_private_key_retired_after_signature": True,
            "report_builder_identity_bound": True,
            "report_builder_tests_identity_bound": True,
        },
        "pass": True,
    }
    paths["manifest"].write_text(
        json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8"
    )
    subprocess.run(
        [
            str(MODULE.OPENSSL),
            "pkeyutl",
            "-sign",
            "-rawin",
            "-inkey",
            str(paths["private_key"]),
            "-in",
            str(paths["manifest"]),
            "-out",
            str(paths["signature"]),
        ],
        check=True,
        capture_output=True,
    )
    for key in ("manifest", "signature"):
        paths[key].chmod(0o444)
    paths["private_key"].chmod(0o000)
    return paths


def test_validate_canonical_pytest_release_evidence_accepts_signed_current(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = make_signed_test_evidence_fixture(tmp_path, monkeypatch)
    observed = MODULE.validate_canonical_pytest_release_evidence(paths["xml"])
    assert observed["tests"] == MODULE.MINIMUM_CANONICAL_TEST_COUNT
    assert observed["testcases"] == MODULE.MINIMUM_CANONICAL_TEST_COUNT
    assert observed["signed_evidence_verified"] is True


def test_validate_canonical_pytest_release_evidence_rejects_noncanonical_xml(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    make_signed_test_evidence_fixture(tmp_path, monkeypatch)
    other = tmp_path / "synthetic.xml"
    write_passing_pytest_xml(other)
    with pytest.raises(MODULE.ReportBuildError, match="signed release path"):
        MODULE.validate_canonical_pytest_release_evidence(other)


def test_validate_canonical_pytest_release_evidence_rejects_report_source_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = make_signed_test_evidence_fixture(tmp_path, monkeypatch)
    paths["builder"].chmod(0o644)
    paths["builder"].write_text("BUILDER = 'drifted'\n", encoding="utf-8")
    paths["builder"].chmod(0o444)
    with pytest.raises(MODULE.ReportBuildError, match="contract drift"):
        MODULE.validate_canonical_pytest_release_evidence(paths["xml"])


def test_validate_canonical_pytest_release_evidence_rejects_unretired_private_key(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = make_signed_test_evidence_fixture(tmp_path, monkeypatch)
    paths["private_key"].chmod(0o400)
    with pytest.raises(MODULE.ReportBuildError, match="not retired"):
        MODULE.validate_canonical_pytest_release_evidence(paths["xml"])


def make_source_authority_fixture(
    tmp_path: Path,
    *,
    review_verdicts: tuple[str, ...] = ("APPROVE", "APPROVE"),
) -> tuple[dict[str, object], Path, Path]:
    private_key = tmp_path / "authority.private.pem"
    public_key = tmp_path / "authority.public.pem"
    subprocess.run(
        [
            str(MODULE.OPENSSL),
            "genpkey",
            "-algorithm",
            "ED25519",
            "-out",
            str(private_key),
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            str(MODULE.OPENSSL),
            "pkey",
            "-in",
            str(private_key),
            "-pubout",
            "-out",
            str(public_key),
        ],
        check=True,
        capture_output=True,
    )
    public_key.chmod(0o444)
    supplemental_auditor = tmp_path / "supplemental_auditor.py"
    supplemental_auditor.write_text("SOURCE = 'fixture'\n", encoding="utf-8")
    supplemental_auditor.chmod(0o444)
    source_set = "a" * 64
    authorized_head = "b" * 40
    authority = tmp_path / "authority.json"
    authority.write_text(
        json.dumps(
            {
                "schema_name": "intersubmod.release_source_authority",
                "schema_version": "1.0.0",
                "authority_id": MODULE.SOURCE_AUTHORITY_ID,
                "approval_status": "APPROVED_FOR_FULL_TASK_B_RUN",
                "source_set_sha256": source_set,
                "repository": {
                    "git_head_at_authorization": authorized_head
                },
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    authority.chmod(0o444)
    authority_record = MODULE.artifact(authority)
    public_key_record = MODULE.artifact(public_key)
    approval = tmp_path / "approval.json"
    approval.write_text(
        json.dumps(
            {
                "schema_name": (
                    "intersubmod.release_source_authority.approval"
                ),
                "schema_version": "1.0.0",
                "authority_id": MODULE.SOURCE_AUTHORITY_ID,
                "approval_status": "APPROVED_FOR_FULL_TASK_B_RUN",
                "authority_manifest": {
                    field: authority_record[field]
                    for field in ("path", "size_bytes", "sha256", "mode")
                },
                "public_key": {
                    field: public_key_record[field]
                    for field in ("path", "size_bytes", "sha256", "mode")
                },
                "review_approvals": [
                    {"reviewer_id": str(index), "verdict": verdict}
                    for index, verdict in enumerate(review_verdicts)
                ],
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    approval.chmod(0o444)
    signature = tmp_path / "approval.json.ed25519.sig"
    subprocess.run(
        [
            str(MODULE.OPENSSL),
            "pkeyutl",
            "-sign",
            "-rawin",
            "-inkey",
            str(private_key),
            "-in",
            str(approval),
            "-out",
            str(signature),
        ],
        check=True,
        capture_output=True,
    )
    signature.chmod(0o444)
    summary = {
        "source_chain": {
            "v4_source_authority": {
                "pass": True,
                "authority_id": MODULE.SOURCE_AUTHORITY_ID,
                "authority_manifest": MODULE.artifact(authority),
                "detached_approval": MODULE.artifact(approval),
                "detached_approval_signature": MODULE.artifact(signature),
                "public_key": MODULE.artifact(public_key),
                "source_set_sha256": source_set,
                "authorized_git_head": authorized_head,
            }
        },
        "code": {
            "supplemental_auditor": MODULE.artifact(
                supplemental_auditor
            )
        },
    }
    return summary, approval, signature


def test_verify_signature_accepts_valid_ed25519_fixture(
    tmp_path: Path,
) -> None:
    receipt, signature, public_key = make_signed_receipt(tmp_path)
    MODULE.verify_signature(
        receipt_path=receipt,
        signature_path=signature,
        public_key_path=public_key,
        expected_public_key_sha256=MODULE.sha256(public_key),
        label="test release",
    )


def test_verify_signature_rejects_tampered_receipt(tmp_path: Path) -> None:
    receipt, signature, public_key = make_signed_receipt(tmp_path)
    receipt.chmod(0o644)
    receipt.write_text('{"pass":false}\n', encoding="utf-8")
    receipt.chmod(0o444)
    with pytest.raises(MODULE.ReportBuildError, match="detached signature"):
        MODULE.verify_signature(
            receipt_path=receipt,
            signature_path=signature,
            public_key_path=public_key,
            expected_public_key_sha256=MODULE.sha256(public_key),
            label="test release",
        )


def test_verify_signature_rejects_mode_drift(tmp_path: Path) -> None:
    receipt, signature, public_key = make_signed_receipt(tmp_path)
    signature.chmod(0o644)
    with pytest.raises(MODULE.ReportBuildError, match="protected mode drift"):
        MODULE.verify_signature(
            receipt_path=receipt,
            signature_path=signature,
            public_key_path=public_key,
            expected_public_key_sha256=MODULE.sha256(public_key),
            label="test release",
        )


def test_verify_signature_rejects_public_key_sha_drift(
    tmp_path: Path,
) -> None:
    receipt, signature, public_key = make_signed_receipt(tmp_path)
    with pytest.raises(MODULE.ReportBuildError, match="public-key SHA-256"):
        MODULE.verify_signature(
            receipt_path=receipt,
            signature_path=signature,
            public_key_path=public_key,
            expected_public_key_sha256="0" * 64,
            label="test release",
        )


def make_formal_release_verifications() -> tuple[
    dict[str, object], dict[str, object]
]:
    source_authority = {"authority_id": "test-v5", "pass": True}
    dataset = {
        "schema_name": (
            "intersubmod.task_b_final_dataset_release_receipt.verification"
        ),
        "signature_verified": True,
        "all_final_outputs_reverified": True,
        "source_authority": source_authority,
        "pass": True,
    }
    report = {
        "schema_name": (
            "intersubmod.task_b_final_report_release_receipt.verification"
        ),
        "signed_dataset_release_reverified": True,
        "signature_verified": True,
        "all_report_outputs_reverified": True,
        "source_authority": source_authority,
        "pass": True,
    }
    return dataset, report


def configure_formal_release_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> dict[str, Path]:
    paths = {
        "dataset_receipt": tmp_path / "dataset.json",
        "dataset_signature": tmp_path / "dataset.sig",
        "dataset_public_key": tmp_path / "dataset.pem",
        "report_receipt": tmp_path / "report.json",
        "report_signature": tmp_path / "report.sig",
        "report_public_key": tmp_path / "report.pem",
    }
    for attribute, key in (
        ("RELEASE_RECEIPT", "dataset_receipt"),
        ("RELEASE_SIGNATURE", "dataset_signature"),
        ("RESULT_PUBLIC_KEY", "dataset_public_key"),
        ("REPORT_RELEASE_RECEIPT", "report_receipt"),
        ("REPORT_RELEASE_SIGNATURE", "report_signature"),
        ("REPORT_PUBLIC_KEY", "report_public_key"),
    ):
        monkeypatch.setattr(MODULE.FINAL_RELEASE, attribute, paths[key])
    return paths


def call_formal_release_validation(
    paths: dict[str, Path],
) -> tuple[dict[str, object], dict[str, object]]:
    return MODULE.validate_formal_release_chain(
        final_dataset_receipt_path=paths["dataset_receipt"],
        final_dataset_signature_path=paths["dataset_signature"],
        final_dataset_public_key_path=paths["dataset_public_key"],
        final_report_receipt_path=paths["report_receipt"],
        final_report_signature_path=paths["report_signature"],
        final_report_public_key_path=paths["report_public_key"],
    )


def test_validate_formal_release_chain_uses_canonical_fd_validators(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = configure_formal_release_paths(tmp_path, monkeypatch)
    dataset, report = make_formal_release_verifications()
    monkeypatch.setattr(
        MODULE.FINAL_RELEASE,
        "validate_release_signature_artifacts",
        lambda: dataset,
    )
    monkeypatch.setattr(
        MODULE.FINAL_RELEASE,
        "validate_report_release_signature_artifacts",
        lambda: report,
    )
    assert call_formal_release_validation(paths) == (dataset, report)


def test_validate_formal_release_chain_rejects_noncanonical_path(
    tmp_path: Path,
) -> None:
    with pytest.raises(MODULE.ReportBuildError, match="not canonical"):
        MODULE.validate_formal_release_chain(
            final_dataset_receipt_path=tmp_path / "other-dataset.json",
            final_dataset_signature_path=MODULE.FINAL_RELEASE.RELEASE_SIGNATURE,
            final_dataset_public_key_path=MODULE.FINAL_RELEASE.RESULT_PUBLIC_KEY,
            final_report_receipt_path=MODULE.FINAL_RELEASE.REPORT_RELEASE_RECEIPT,
            final_report_signature_path=(
                MODULE.FINAL_RELEASE.REPORT_RELEASE_SIGNATURE
            ),
            final_report_public_key_path=MODULE.FINAL_RELEASE.REPORT_PUBLIC_KEY,
        )


def test_validate_formal_release_chain_rejects_report_contract_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = configure_formal_release_paths(tmp_path, monkeypatch)
    dataset, report = make_formal_release_verifications()
    report["all_report_outputs_reverified"] = False
    monkeypatch.setattr(
        MODULE.FINAL_RELEASE,
        "validate_release_signature_artifacts",
        lambda: dataset,
    )
    monkeypatch.setattr(
        MODULE.FINAL_RELEASE,
        "validate_report_release_signature_artifacts",
        lambda: report,
    )
    with pytest.raises(MODULE.ReportBuildError, match="contract drift"):
        call_formal_release_validation(paths)


def test_validate_pre_publish_state_detects_late_input_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tracked = tmp_path / "tracked.json"
    tracked.write_text('{"pass":true}\n', encoding="utf-8")
    input_paths = {
        "canonical_pytest_xml": tracked,
        "final_dataset_receipt": tracked,
        "final_dataset_signature": tracked,
        "final_dataset_public_key": tracked,
        "final_report_receipt": tracked,
        "final_report_signature": tracked,
        "final_report_public_key": tracked,
    }
    identities = {
        name: MODULE.artifact(path) for name, path in input_paths.items()
    }
    dataset, report = make_formal_release_verifications()
    monkeypatch.setattr(
        MODULE, "repository_metadata", lambda _summary: {"head": "test"}
    )
    monkeypatch.setattr(
        MODULE,
        "validate_canonical_pytest_release_evidence",
        lambda _path: {"tests": MODULE.MINIMUM_CANONICAL_TEST_COUNT},
    )
    monkeypatch.setattr(
        MODULE, "validate_source_authority_chain", lambda _summary: None
    )
    monkeypatch.setattr(
        MODULE,
        "validate_formal_release_chain",
        lambda **_kwargs: (dataset, report),
    )
    source = MODULE.artifact(MODULE.Path(MODULE.__file__).resolve())
    tracked.write_text('{"pass":false}\n', encoding="utf-8")
    with pytest.raises(MODULE.ReportBuildError, match="before atomic publish"):
        MODULE.validate_pre_publish_state(
            input_paths=input_paths,
            summary={},
            identities_expected=identities,
            source_expected=source,
            repository_expected={"head": "test"},
            verification_expected={
                "tests": MODULE.MINIMUM_CANONICAL_TEST_COUNT
            },
            dataset_verification_expected=dataset,
            report_verification_expected=report,
        )


def test_validate_pre_publish_state_detects_transitive_authority_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    summary, approval, _ = make_source_authority_fixture(tmp_path)
    public_key = Path(
        summary["source_chain"]["v4_source_authority"]["public_key"]["path"]
    )
    monkeypatch.setattr(
        MODULE,
        "SOURCE_AUTHORITY_PUBLIC_KEY_SHA256",
        MODULE.sha256(public_key),
    )
    tracked = tmp_path / "tracked.json"
    tracked.write_text('{"pass":true}\n', encoding="utf-8")
    input_paths = {
        "canonical_pytest_xml": tracked,
        "final_dataset_receipt": tracked,
        "final_dataset_signature": tracked,
        "final_dataset_public_key": tracked,
        "final_report_receipt": tracked,
        "final_report_signature": tracked,
        "final_report_public_key": tracked,
    }
    input_paths.update(MODULE.resolve_source_authority_input_paths(summary))
    identities = {
        name: MODULE.artifact(path) for name, path in input_paths.items()
    }
    dataset, report = make_formal_release_verifications()
    monkeypatch.setattr(
        MODULE, "repository_metadata", lambda _summary: {"head": "test"}
    )
    monkeypatch.setattr(
        MODULE,
        "validate_canonical_pytest_release_evidence",
        lambda _path: {"tests": MODULE.MINIMUM_CANONICAL_TEST_COUNT},
    )
    monkeypatch.setattr(
        MODULE,
        "validate_formal_release_chain",
        lambda **_kwargs: (dataset, report),
    )
    source = MODULE.artifact(MODULE.Path(MODULE.__file__).resolve())
    approval.chmod(0o644)
    approval.write_text('{"pass":false}\n', encoding="utf-8")
    approval.chmod(0o444)
    with pytest.raises(MODULE.ReportBuildError, match="identity drift"):
        MODULE.validate_pre_publish_state(
            input_paths=input_paths,
            summary=summary,
            identities_expected=identities,
            source_expected=source,
            repository_expected={"head": "test"},
            verification_expected={
                "tests": MODULE.MINIMUM_CANONICAL_TEST_COUNT
            },
            dataset_verification_expected=dataset,
            report_verification_expected=report,
        )


def test_main_orders_publish_guard_after_outputs_and_before_rename() -> None:
    source = inspect.getsource(MODULE.main)
    assert source.index("success_path.write_text") < source.index(
        "fsync_directory(staging)"
    )
    assert source.index("fsync_directory(staging)") < source.index(
        "staging.chmod(0o555)"
    )
    assert source.index("staging.chmod(0o555)") < source.index(
        "validate_and_publish("
    )
    assert "rename_no_replace(" not in source


def test_validate_and_publish_contains_only_guard_then_atomic_rename(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = inspect.getsource(MODULE.validate_and_publish)
    function = ast.parse(source).body[0]
    calls = [
        node.value.func.id
        for node in function.body
        if isinstance(node, ast.Expr)
        and isinstance(node.value, ast.Call)
        and isinstance(node.value.func, ast.Name)
    ]
    assert calls == ["validate_pre_publish_state", "rename_no_replace"]

    events: list[str] = []
    monkeypatch.setattr(
        MODULE,
        "validate_pre_publish_state",
        lambda **_kwargs: events.append("guard"),
    )
    monkeypatch.setattr(
        MODULE,
        "rename_no_replace",
        lambda _source, _target: events.append("rename"),
    )
    dataset, report = make_formal_release_verifications()
    MODULE.validate_and_publish(
        staging=Path("/staging"),
        output_dir=Path("/output"),
        input_paths={},
        summary={},
        identities_expected={},
        source_expected={},
        repository_expected={},
        verification_expected={},
        dataset_verification_expected=dataset,
        report_verification_expected=report,
    )
    assert events == ["guard", "rename"]


def test_validate_source_authority_chain_accepts_signed_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    summary, _, _ = make_source_authority_fixture(tmp_path)
    public_key = Path(
        summary["source_chain"]["v4_source_authority"]["public_key"]["path"]
    )
    monkeypatch.setattr(
        MODULE,
        "SOURCE_AUTHORITY_PUBLIC_KEY_SHA256",
        MODULE.sha256(public_key),
    )
    MODULE.validate_source_authority_chain(summary)


@pytest.mark.parametrize(
    "review_verdicts",
    [("APPROVE",), ("APPROVE", "FIX")],
)
def test_validate_source_authority_chain_rejects_review_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    review_verdicts: tuple[str, ...],
) -> None:
    summary, _, _ = make_source_authority_fixture(
        tmp_path, review_verdicts=review_verdicts
    )
    public_key = Path(
        summary["source_chain"]["v4_source_authority"]["public_key"]["path"]
    )
    monkeypatch.setattr(
        MODULE,
        "SOURCE_AUTHORITY_PUBLIC_KEY_SHA256",
        MODULE.sha256(public_key),
    )
    with pytest.raises(MODULE.ReportBuildError, match="signed contract drift"):
        MODULE.validate_source_authority_chain(summary)


def test_validate_source_authority_chain_rejects_tampered_approval(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    summary, approval, _ = make_source_authority_fixture(tmp_path)
    public_key = Path(
        summary["source_chain"]["v4_source_authority"]["public_key"]["path"]
    )
    monkeypatch.setattr(
        MODULE,
        "SOURCE_AUTHORITY_PUBLIC_KEY_SHA256",
        MODULE.sha256(public_key),
    )
    approval.chmod(0o644)
    value = json.loads(approval.read_text(encoding="utf-8"))
    value["approval_status"] = "REVOKED"
    approval.write_text(json.dumps(value) + "\n", encoding="utf-8")
    approval.chmod(0o444)
    summary["source_chain"]["v4_source_authority"][
        "detached_approval"
    ] = MODULE.artifact(approval)
    with pytest.raises(MODULE.ReportBuildError, match="detached signature"):
        MODULE.validate_source_authority_chain(summary)


def test_validate_source_authority_chain_rejects_source_set_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    summary, _, _ = make_source_authority_fixture(tmp_path)
    public_key = Path(
        summary["source_chain"]["v4_source_authority"]["public_key"]["path"]
    )
    monkeypatch.setattr(
        MODULE,
        "SOURCE_AUTHORITY_PUBLIC_KEY_SHA256",
        MODULE.sha256(public_key),
    )
    summary["source_chain"]["v4_source_authority"][
        "source_set_sha256"
    ] = "0" * 64
    with pytest.raises(MODULE.ReportBuildError, match="signed contract drift"):
        MODULE.validate_source_authority_chain(summary)


def test_validate_source_authority_chain_rejects_authority_id_drift(
    tmp_path: Path,
) -> None:
    summary, _, _ = make_source_authority_fixture(tmp_path)
    summary["source_chain"]["v4_source_authority"][
        "authority_id"
    ] = "wrong"
    with pytest.raises(MODULE.ReportBuildError, match="incomplete"):
        MODULE.validate_source_authority_chain(summary)


def test_require_identity_rejects_content_drift(tmp_path: Path) -> None:
    path = tmp_path / "artifact.txt"
    path.write_text("before\n", encoding="utf-8")
    record = MODULE.artifact(path)
    path.write_text("after\n", encoding="utf-8")
    with pytest.raises(MODULE.ReportBuildError, match="identity drift"):
        MODULE.require_identity(path, record, label="test artifact")


def test_require_identity_rejects_protected_mode_drift(
    tmp_path: Path,
) -> None:
    path = tmp_path / "artifact.txt"
    path.write_text("content\n", encoding="utf-8")
    path.chmod(0o444)
    record = MODULE.artifact(path)
    path.chmod(0o644)
    with pytest.raises(MODULE.ReportBuildError, match="mode"):
        MODULE.require_identity(path, record, label="test artifact")


def make_summary_contract_fixture() -> dict[str, object]:
    counts = {
        "singleton_sites": MODULE.EXPECTED_COUNTS["singleton_sites"],
        "m1_evaluable": MODULE.EXPECTED_COUNTS["m1_evaluable"],
        "m1_not_evaluable": (
            MODULE.EXPECTED_COUNTS["singleton_sites"]
            - MODULE.EXPECTED_COUNTS["m1_evaluable"]
        ),
        "m1_flagged": MODULE.EXPECTED_COUNTS["m1_flagged"],
        "m2_pass": MODULE.EXPECTED_COUNTS["m2_pass"],
        "m2_fail": MODULE.EXPECTED_COUNTS["m2_fail"],
        "m2_not_evaluable": MODULE.EXPECTED_COUNTS["m2_not_evaluable"],
        "m2_not_run": MODULE.EXPECTED_COUNTS["m2_not_run"],
    }
    total_row = {
        "sites": counts["singleton_sites"],
        "m1_evaluable": counts["m1_evaluable"],
        "m1_flagged": counts["m1_flagged"],
        "m2_pass": counts["m2_pass"],
        "m2_fail": counts["m2_fail"],
        "m2_not_evaluable": counts["m2_not_evaluable"],
        "m2_not_run": counts["m2_not_run"],
    }
    zero_row = {field: 0 for field in total_row}
    first_dataset = MODULE.EXPECTED_DATASETS[0]
    datasets = {
        dataset: dict(total_row if dataset == first_dataset else zero_row)
        for dataset in MODULE.EXPECTED_DATASETS
    }
    truths = {
        truth: dict(total_row if truth == "TP" else zero_row)
        for truth in ("TP", "FP", "UNASSESSED")
    }
    dataset_truth = {
        f"{dataset}|{truth}": dict(
            total_row
            if dataset == first_dataset and truth == "TP"
            else zero_row
        )
        for dataset in MODULE.EXPECTED_DATASETS
        for truth in ("TP", "FP", "UNASSESSED")
    }
    all_rows = MODULE.EXPECTED_COUNTS["all_dataset_sites"]
    dataset_counts = {}
    for dataset in MODULE.EXPECTED_DATASETS:
        is_first = dataset == first_dataset
        dataset_counts[dataset] = {
            "all_ssnv": all_rows if is_first else 0,
            "branch_max_snv_excluded": (
                all_rows - counts["singleton_sites"] if is_first else 0
            ),
            "branch_positional_singleton": (
                counts["singleton_sites"] if is_first else 0
            ),
            "branch_retained": 0,
            "truth_fp": 0,
            "truth_tp": all_rows if is_first else 0,
            "truth_unassessed": 0,
        }
    biological_ids = {
        dataset: ("HCC1395" if dataset == "HCC1395_DORADO" else dataset)
        for dataset in MODULE.EXPECTED_DATASETS
    }
    return {
        "counts": counts,
        "status_census": {
            "m1": {
                "FLAGGED": counts["m1_flagged"],
                "NOT_FLAGGED": counts["m2_not_run"],
            },
            "m2": {
                "PASS": counts["m2_pass"],
                "FAIL": counts["m2_fail"],
                "NOT_EVALUABLE": counts["m2_not_evaluable"],
                "NOT_RUN": counts["m2_not_run"],
            },
            "g1": {"NOT_RUN": counts["singleton_sites"]},
            "g2": {"NOT_RUN": counts["singleton_sites"]},
            "r1": {"NOT_RUN": counts["singleton_sites"]},
        },
        "methyl_group_count_distribution": {
            "2": counts["m1_flagged"]
        },
        "screen_metadata": {
            "all_rows": all_rows,
            "datasets": list(MODULE.EXPECTED_DATASETS),
            "unique_keys": all_rows,
            "duplicate_keys": 0,
            "branch_counts": {
                "max_snv_excluded": all_rows - counts["singleton_sites"],
                "positional_singleton": counts["singleton_sites"],
                "retained": 0,
            },
            "dataset_biological_ids": biological_ids,
            "dataset_counts": dataset_counts,
            "truth_counts": {
                "TP": all_rows,
                "FP": 0,
                "UNASSESSED": 0,
            },
        },
        "positional_recomputation": {
            "component_rows": all_rows,
            "recomputed_singletons": counts["singleton_sites"],
            "component_identity_mismatch": 0,
            "singleton_branch_mismatch": 0,
            "singleton_positional_contract_failure": 0,
            "minimum_finite_singleton_nearest_gap_bp": 50_003,
        },
        "breakdowns": {
            "dataset": datasets,
            "truth": truths,
            "dataset_truth": dataset_truth,
        },
    }


def test_validate_summary_contract_accepts_conserved_fixture() -> None:
    MODULE.validate_summary_contract(make_summary_contract_fixture())


def test_validate_summary_contract_rejects_m2_conservation_drift() -> None:
    summary = make_summary_contract_fixture()
    summary["counts"]["m2_fail"] += 1
    with pytest.raises(MODULE.ReportBuildError, match="within-M1"):
        MODULE.validate_summary_contract(summary)


def test_validate_summary_contract_rejects_status_census_drift() -> None:
    summary = make_summary_contract_fixture()
    summary["status_census"]["m2"]["PASS"] -= 1
    summary["status_census"]["m2"]["FAIL"] += 1
    with pytest.raises(MODULE.ReportBuildError, match="status census"):
        MODULE.validate_summary_contract(summary)


def test_validate_summary_contract_rejects_visual_distribution_drift() -> None:
    summary = make_summary_contract_fixture()
    summary["methyl_group_count_distribution"]["2"] -= 1
    with pytest.raises(MODULE.ReportBuildError, match="distribution"):
        MODULE.validate_summary_contract(summary)


def test_validate_summary_contract_rejects_breakdown_drift() -> None:
    summary = make_summary_contract_fixture()
    summary["breakdowns"]["dataset"]["COLO829"]["m2_pass"] -= 1
    with pytest.raises(MODULE.ReportBuildError, match="Dataset breakdown"):
        MODULE.validate_summary_contract(summary)


def test_validate_summary_contract_rejects_positional_boundary_drift() -> None:
    summary = make_summary_contract_fixture()
    summary["positional_recomputation"][
        "minimum_finite_singleton_nearest_gap_bp"
    ] = 50_000
    with pytest.raises(MODULE.ReportBuildError, match="Positional singleton"):
        MODULE.validate_summary_contract(summary)


def test_validate_inputs_rejects_count_drift_before_release(
    tmp_path: Path,
) -> None:
    summary = {
        "schema_name": MODULE.AUDIT_SCHEMA,
        "schema_version": "1.0.0",
        "pass": True,
        "checks": {"all_rows": True},
        "counts": {
            **MODULE.EXPECTED_COUNTS,
            "singleton_sites": MODULE.EXPECTED_COUNTS["singleton_sites"] - 1,
        },
        "screen_metadata": {
            "all_rows": MODULE.EXPECTED_COUNTS["all_dataset_sites"]
        },
    }
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    missing = tmp_path / "missing"
    with pytest.raises(MODULE.ReportBuildError, match="count drift"):
        MODULE.validate_inputs(
            summary_path=summary_path,
            site_audit_path=missing,
            candidate_path=missing,
            assignments_path=missing,
            final_dataset_receipt_path=missing,
            final_dataset_signature_path=missing,
            final_dataset_public_key_path=missing,
            final_report_receipt_path=missing,
            final_report_signature_path=missing,
            final_report_public_key_path=missing,
            tumor_ref_controls_path=missing,
        )


@pytest.mark.parametrize(
    (
        "tests",
        "failures",
        "errors",
        "skipped",
        "include_required",
        "passes",
    ),
    [
        (488, 0, 0, 0, True, True),
        (487, 0, 0, 0, True, False),
        (488, 1, 0, 0, True, False),
        (488, 0, 1, 0, True, False),
        (488, 0, 0, 1, True, False),
        (488, 0, 0, 0, False, False),
    ],
)
def test_validate_pytest_xml_contract(
    tmp_path: Path,
    tests: int,
    failures: int,
    errors: int,
    skipped: int,
    include_required: bool,
    passes: bool,
) -> None:
    path = tmp_path / "pytest.xml"
    required_names = (
        sorted(MODULE.REQUIRED_SUPPLEMENTAL_TESTS)
        if include_required
        else []
    )
    filler_names = [
        f"test_filler_{index}"
        for index in range(tests - len(required_names))
    ]
    cases = "".join(
        f'<testcase classname="current" name="{name}" />'
        for name in required_names + filler_names
    )
    path.write_text(
        (
            '<testsuites><testsuite tests="'
            f'{tests}" failures="{failures}" errors="{errors}" '
            f'skipped="{skipped}">'
            f"{cases}</testsuite></testsuites>"
        ),
        encoding="utf-8",
    )
    if passes:
        observed = MODULE.validate_pytest_xml(path)
        assert observed["tests"] == tests
        assert observed["required_supplemental_tests"] == len(
            MODULE.REQUIRED_SUPPLEMENTAL_TESTS
        )
    else:
        with pytest.raises(MODULE.ReportBuildError, match="did not pass"):
            MODULE.validate_pytest_xml(path)


@pytest.mark.parametrize(
    ("payload", "message"),
    [
        ("<testsuites>", "Invalid canonical pytest XML"),
        ("<testsuites></testsuites>", "contains no test suites"),
    ],
)
def test_validate_pytest_xml_rejects_malformed_contract(
    tmp_path: Path, payload: str, message: str
) -> None:
    path = tmp_path / "pytest.xml"
    path.write_text(payload, encoding="utf-8")
    with pytest.raises(MODULE.ReportBuildError, match=message):
        MODULE.validate_pytest_xml(path)


def write_candidate_fixture(
    path: Path, rows: list[dict[str, str]]
) -> None:
    fields = [
        "dataset",
        "chrom",
        "pos",
        "ref",
        "alt",
        "methyl_group_count",
        "core_read_n",
        "min_group_n",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_load_candidates_rejects_count_drift(tmp_path: Path) -> None:
    path = tmp_path / "candidates.tsv"
    write_candidate_fixture(path, [])
    with pytest.raises(MODULE.ReportBuildError, match="Expected 30"):
        MODULE.load_candidates(path)


def test_load_candidates_rejects_duplicate_key(tmp_path: Path) -> None:
    path = tmp_path / "candidates.tsv"
    rows = []
    for index in range(MODULE.EXPECTED_COUNTS["m2_pass"]):
        position = 1 if index < 2 else index
        rows.append(
            {
                "dataset": "S",
                "chrom": "chr1",
                "pos": str(position),
                "ref": "A",
                "alt": "C",
                "methyl_group_count": "2",
                "core_read_n": "20",
                "min_group_n": "10",
            }
        )
    write_candidate_fixture(path, rows)
    with pytest.raises(MODULE.ReportBuildError, match="Duplicate M2 PASS"):
        MODULE.load_candidates(path)


def test_repository_metadata_rejects_authorized_head_drift() -> None:
    summary = {
        "source_chain": {
            "v4_source_authority": {"authorized_git_head": "0" * 40}
        }
    }
    with pytest.raises(MODULE.ReportBuildError, match="Repository HEAD drift"):
        MODULE.repository_metadata(summary)


def write_singleton_site_fixture(
    path: Path,
    rows: list[dict[str, str]],
    *,
    fields: list[str] | None = None,
) -> None:
    fieldnames = fields or [
        "dataset",
        "chrom",
        "pos",
        "ref",
        "alt",
        "m1_status",
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)


def test_load_m1_singleton_keys_rejects_missing_fields(
    tmp_path: Path,
) -> None:
    path = tmp_path / "sites.tsv.gz"
    write_singleton_site_fixture(path, [], fields=["dataset", "chrom"])
    with pytest.raises(MODULE.ReportBuildError, match="lacks M1 key fields"):
        MODULE.load_m1_singleton_keys(path)


def test_load_m1_singleton_keys_rejects_duplicate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "sites.tsv.gz"
    row = {
        "dataset": "S",
        "chrom": "chr1",
        "pos": "1",
        "ref": "A",
        "alt": "C",
        "m1_status": "FLAGGED",
    }
    write_singleton_site_fixture(path, [row, row])
    monkeypatch.setitem(MODULE.EXPECTED_COUNTS, "m1_flagged", 2)
    with pytest.raises(MODULE.ReportBuildError, match="Duplicate singleton"):
        MODULE.load_m1_singleton_keys(path)


def test_load_m1_singleton_keys_rejects_count_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "sites.tsv.gz"
    write_singleton_site_fixture(
        path,
        [
            {
                "dataset": "S",
                "chrom": "chr1",
                "pos": "1",
                "ref": "A",
                "alt": "C",
                "m1_status": "FLAGGED",
            }
        ],
    )
    monkeypatch.setitem(MODULE.EXPECTED_COUNTS, "m1_flagged", 2)
    with pytest.raises(MODULE.ReportBuildError, match="key count drift"):
        MODULE.load_m1_singleton_keys(path)


def write_tumor_ref_fixture(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "dataset",
        "chrom",
        "pos",
        "ref",
        "alt",
        "n_tumor_ref",
        "ref_status",
        "ref_evaluable",
        "ref_stable_null_multigroup",
        "joint_stable_null_multigroup",
        "joint_allele_axis_aligned",
        "background_control_interpretation",
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_tumor_ref_control_join_and_conservation(tmp_path: Path) -> None:
    keys = {
        ("S", "chr1", 1, "A", "C"),
        ("S", "chr1", 2, "G", "T"),
        ("S", "chr1", 3, "C", "A"),
    }
    rows = [
        {
            "dataset": "S",
            "chrom": "chr1",
            "pos": str(index),
            "ref": ref,
            "alt": alt,
            "n_tumor_ref": "20" if index < 3 else "0",
            "ref_status": "evaluable" if index < 3 else "insufficient_reads",
            "ref_evaluable": "true" if index < 3 else "false",
            "ref_stable_null_multigroup": "true" if index == 1 else "false",
            "joint_stable_null_multigroup": "false",
            "joint_allele_axis_aligned": "true" if index == 2 else "false",
            "background_control_interpretation": (
                "REF_REPLICATION_WEAKENS_ALT_SPECIFICITY"
                if index == 1
                else (
                    "REF_NONREPLICATION_DOES_NOT_PROVE_SUBCLONE"
                    if index == 2
                    else "REF_CONTROL_NOT_TESTABLE"
                )
            ),
        }
        for index, ref, alt in (
            (1, "A", "C"),
            (2, "G", "T"),
            (3, "C", "A"),
        )
    ]
    path = tmp_path / "tumor_ref.tsv.gz"
    write_tumor_ref_fixture(path, rows)
    controls, summary = MODULE.load_tumor_ref_controls(
        path,
        m1_keys=keys,
        m2_keys={("S", "chr1", 2, "G", "T")},
    )
    assert len(controls) == 3
    assert summary["m1"] == {
        "sites": 3,
        "ref_evaluable": 2,
        "ref_stable_multigroup": 1,
        "ref_nonreplication": 1,
        "ref_not_evaluable": 1,
        "joint_stable_multigroup": 0,
        "joint_allele_axis_aligned": 1,
    }
    assert summary["m2_pass"]["ref_nonreplication"] == 1


def test_tumor_ref_control_rejects_interpretation_drift(
    tmp_path: Path,
) -> None:
    key = ("S", "chr1", 1, "A", "C")
    path = tmp_path / "tumor_ref.tsv.gz"
    write_tumor_ref_fixture(
        path,
        [
            {
                "dataset": "S",
                "chrom": "chr1",
                "pos": "1",
                "ref": "A",
                "alt": "C",
                "n_tumor_ref": "20",
                "ref_status": "evaluable",
                "ref_evaluable": "true",
                "ref_stable_null_multigroup": "false",
                "joint_stable_null_multigroup": "false",
                "joint_allele_axis_aligned": "false",
                "background_control_interpretation": "REF_CONTROL_NOT_TESTABLE",
            }
        ],
    )
    with pytest.raises(MODULE.ReportBuildError, match="interpretation drift"):
        MODULE.load_tumor_ref_controls(
            path, m1_keys={key}, m2_keys={key}
        )


def test_tumor_ref_control_rejects_m2_outside_m1(
    tmp_path: Path,
) -> None:
    path = tmp_path / "unused.tsv.gz"
    with pytest.raises(MODULE.ReportBuildError, match="not a subset"):
        MODULE.load_tumor_ref_controls(
            path,
            m1_keys={("S", "chr1", 1, "A", "C")},
            m2_keys={("S", "chr1", 2, "G", "T")},
        )


def test_tumor_ref_control_rejects_missing_key(tmp_path: Path) -> None:
    path = tmp_path / "tumor_ref.tsv.gz"
    write_tumor_ref_fixture(path, [])
    key = ("S", "chr1", 1, "A", "C")
    with pytest.raises(MODULE.ReportBuildError, match="lack tumor-REF"):
        MODULE.load_tumor_ref_controls(
            path, m1_keys={key}, m2_keys=set()
        )


@pytest.mark.parametrize(
    ("evaluable", "stable", "message"),
    [
        ("TRUE", "false", "Invalid tumor-REF boolean"),
        ("false", "true", "Non-evaluable tumor-REF control is stable"),
    ],
)
def test_tumor_ref_control_rejects_boolean_contract_drift(
    tmp_path: Path,
    evaluable: str,
    stable: str,
    message: str,
) -> None:
    path = tmp_path / "tumor_ref.tsv.gz"
    key = ("S", "chr1", 1, "A", "C")
    write_tumor_ref_fixture(
        path,
        [
            {
                "dataset": "S",
                "chrom": "chr1",
                "pos": "1",
                "ref": "A",
                "alt": "C",
                "n_tumor_ref": "0",
                "ref_status": "insufficient_reads",
                "ref_evaluable": evaluable,
                "ref_stable_null_multigroup": stable,
                "joint_stable_null_multigroup": "false",
                "joint_allele_axis_aligned": "false",
                "background_control_interpretation": (
                    "REF_CONTROL_NOT_TESTABLE"
                ),
            }
        ],
    )
    with pytest.raises(MODULE.ReportBuildError, match=message):
        MODULE.load_tumor_ref_controls(
            path, m1_keys={key}, m2_keys=set()
        )


def test_technical_pair_overlap_exact_contract(tmp_path: Path) -> None:
    path = tmp_path / "sites.tsv.gz"
    fields = [
        "dataset",
        "chrom",
        "pos",
        "ref",
        "alt",
        "m1_status",
        "m2_status",
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for dataset, maximum in (
            ("HCC1395", 7_484),
            ("HCC1395_DORADO", 9_116),
        ):
            for position in range(1, maximum + 1):
                if dataset == "HCC1395":
                    m1 = position <= 407
                    m2 = position in {1, 2}
                else:
                    m1 = position <= 1_289
                    m2 = position in {3, 4}
                writer.writerow(
                    {
                        "dataset": dataset,
                        "chrom": "chr1",
                        "pos": position,
                        "ref": "A",
                        "alt": "C",
                        "m1_status": "FLAGGED" if m1 else "NOT_FLAGGED",
                        "m2_status": "PASS" if m2 else "NOT_RUN",
                    }
                )
    result = MODULE.technical_pair_overlap(path)
    assert result["all"]["intersection"] == 7_484
    assert result["all"]["union"] == 9_116
    assert result["m1"]["intersection"] == 407
    assert result["m1"]["union"] == 1_289
    assert result["m2"]["intersection"] == 0
    assert result["m2"]["union"] == 4


def test_technical_pair_overlap_rejects_contract_drift(
    tmp_path: Path,
) -> None:
    path = tmp_path / "sites.tsv.gz"
    fields = [
        "dataset",
        "chrom",
        "pos",
        "ref",
        "alt",
        "m1_status",
        "m2_status",
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "dataset": "HCC1395",
                "chrom": "chr1",
                "pos": "1",
                "ref": "A",
                "alt": "C",
                "m1_status": "NOT_FLAGGED",
                "m2_status": "NOT_RUN",
            }
        )
    with pytest.raises(MODULE.ReportBuildError, match="overlap drift"):
        MODULE.technical_pair_overlap(path)


def test_main_rejects_existing_output_before_reading_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "existing"
    output.mkdir()
    names = [
        "--audit-summary",
        "--site-audit",
        "--candidate-cases",
        "--stable-assignments",
        "--tumor-ref-controls",
        "--canonical-pytest-xml",
        "--final-dataset-receipt",
        "--final-dataset-signature",
        "--final-dataset-public-key",
        "--final-report-receipt",
        "--final-report-signature",
        "--final-report-public-key",
    ]
    arguments = [str(SCRIPT)]
    for name in names:
        arguments.extend([name, str(tmp_path / name.removeprefix("--"))])
    arguments.extend(["--output-dir", str(output)])
    monkeypatch.setattr(sys, "argv", arguments)
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.main()
