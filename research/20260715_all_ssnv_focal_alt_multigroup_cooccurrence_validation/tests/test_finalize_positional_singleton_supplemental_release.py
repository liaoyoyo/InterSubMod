from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import stat
import subprocess
import sys

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "finalize_positional_singleton_supplemental_release.py"
)
SPEC = importlib.util.spec_from_file_location("supplemental_release", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_json(path: Path, value: dict[str, object]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    path.chmod(0o444)


@pytest.mark.parametrize(
    "payload",
    [
        b'{"pass":true,"pass":false}',
        b'{"value":NaN}',
    ],
)
def test_supplemental_json_is_strict(
    payload: bytes,
    tmp_path: Path,
) -> None:
    with pytest.raises(MODULE.ReleaseError, match="Unable to read"):
        MODULE.parse_json_bytes(
            payload,
            path=tmp_path / "signed.json",
            label="signed supplemental receipt",
        )


def create_key(tmp_path: Path, name: str) -> tuple[Path, Path]:
    private_key = tmp_path / f"{name}.private.pem"
    public_key = tmp_path / f"{name}.public.pem"
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
    return private_key, public_key


def sign_receipt(receipt: Path, private_key: Path) -> Path:
    signature = receipt.with_suffix(receipt.suffix + ".ed25519.sig")
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
    signature.chmod(0o444)
    return signature


def test_atomic_create_refuses_existing_output(tmp_path: Path) -> None:
    output = tmp_path / "receipt.json"
    MODULE.atomic_create(output, {"pass": True})
    assert stat.S_IMODE(output.stat().st_mode) == 0o444
    with pytest.raises(FileExistsError, match="Refusing to overwrite"):
        MODULE.atomic_create(output, {"pass": False})
    assert json.loads(output.read_text())["pass"] is True


def create_release_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> dict[str, Path | list[str]]:
    audit_summary = tmp_path / "audit_summary.json"
    write_json(
        audit_summary,
        {
            "schema_name": (
                "intersubmod.positional_singleton_methyl_multigroup_audit"
            ),
            "schema_version": "2.0.0",
            "pass": True,
            "checks": {"complete": True},
            "counts": {
                "singleton_sites": 50_432,
                "m1_evaluable": 48_347,
                "m1_flagged": 5_961,
                "m2_pass": 30,
            },
        },
    )
    audit_success = tmp_path / "audit_success.json"
    write_json(
        audit_success,
        {
            "schema_name": "intersubmod.atomic_release_marker",
            "schema_version": "1.0.0",
            "summary": MODULE.artifact(audit_summary),
            "pass": True,
        },
    )
    report_receipt = tmp_path / "report_receipt.json"
    write_json(
        report_receipt,
        {
            "schema_name": (
                "intersubmod.positional_singleton_methyl_multigroup_report"
            ),
            "schema_version": "1.0.0",
            "inputs": {"audit_summary": MODULE.artifact(audit_summary)},
            "formal_release_links": {
                "final_dataset_receipt_pass": True,
                "final_report_receipt_pass": True,
                "final_dataset_builder_receipt_pass": True,
            },
            "claim_ceiling": (
                "M2_read_level_residual_epigenetic_partition"
            ),
            "pass": True,
        },
    )
    report_success = tmp_path / "report_success.json"
    write_json(
        report_success,
        {
            "schema_name": "intersubmod.atomic_release_marker",
            "schema_version": "1.0.0",
            "receipt": MODULE.artifact(report_receipt),
            "pass": True,
        },
    )

    dataset_receipt = tmp_path / "dataset_receipt.json"
    write_json(
        dataset_receipt,
        {
            "schema_name": (
                "intersubmod.task_b_final_dataset_release_receipt"
            ),
            "pass": True,
        },
    )
    dataset_private, dataset_public = create_key(tmp_path, "dataset")
    dataset_signature = sign_receipt(dataset_receipt, dataset_private)

    formal_report_receipt = tmp_path / "formal_report_receipt.json"
    write_json(
        formal_report_receipt,
        {
            "schema_name": (
                "intersubmod.task_b_final_report_release_receipt"
            ),
            "pass": True,
        },
    )
    report_private, report_public = create_key(tmp_path, "report")
    formal_report_signature = sign_receipt(
        formal_report_receipt, report_private
    )

    supplemental_private, supplemental_public = create_key(
        tmp_path, "supplemental"
    )
    monkeypatch.setattr(
        MODULE,
        "FINAL_DATASET_PUBLIC_KEY_SHA256",
        MODULE.sha256(dataset_public),
    )
    monkeypatch.setattr(
        MODULE,
        "FINAL_REPORT_PUBLIC_KEY_SHA256",
        MODULE.sha256(report_public),
    )
    monkeypatch.setattr(
        MODULE,
        "SUPPLEMENTAL_PUBLIC_KEY_SHA256",
        MODULE.sha256(supplemental_public),
    )
    output = tmp_path / "supplemental_release.json"
    argv = [
        str(SCRIPT),
        "--audit-summary",
        str(audit_summary),
        "--audit-success",
        str(audit_success),
        "--report-receipt",
        str(report_receipt),
        "--report-success",
        str(report_success),
        "--final-dataset-receipt",
        str(dataset_receipt),
        "--final-dataset-signature",
        str(dataset_signature),
        "--final-dataset-public-key",
        str(dataset_public),
        "--final-report-receipt",
        str(formal_report_receipt),
        "--final-report-signature",
        str(formal_report_signature),
        "--final-report-public-key",
        str(report_public),
        "--supplemental-public-key",
        str(supplemental_public),
        "--output",
        str(output),
    ]
    return {
        "argv": argv,
        "audit_summary": audit_summary,
        "audit_success": audit_success,
        "report_receipt": report_receipt,
        "report_success": report_success,
        "dataset_receipt": dataset_receipt,
        "dataset_signature": dataset_signature,
        "formal_report_receipt": formal_report_receipt,
        "formal_report_signature": formal_report_signature,
        "supplemental_public": supplemental_public,
        "supplemental_private": supplemental_private,
        "output": output,
    }


def test_supplemental_release_main_with_real_signatures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = create_release_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(sys, "argv", fixture["argv"])
    MODULE.main()
    output = fixture["output"]
    assert isinstance(output, Path)
    receipt = json.loads(output.read_text(encoding="utf-8"))
    assert receipt["pass"] is True
    assert (
        receipt["schema_name"]
        == "intersubmod.positional_singleton_supplemental_release_receipt"
    )
    assert receipt["counts"]["singleton_sites"] == 50_432
    assert (
        receipt["signature_contract"]["state_at_receipt_creation"]
        == "UNSIGNED_PENDING_OUT_OF_BAND"
    )
    assert stat.S_IMODE(output.stat().st_mode) == 0o444


def test_supplemental_release_rejects_formal_link_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = create_release_fixture(tmp_path, monkeypatch)
    report_receipt = fixture["report_receipt"]
    assert isinstance(report_receipt, Path)
    value = json.loads(report_receipt.read_text(encoding="utf-8"))
    value["formal_release_links"]["final_report_receipt_pass"] = False
    report_receipt.chmod(0o644)
    write_json(report_receipt, value)
    monkeypatch.setattr(sys, "argv", fixture["argv"])
    with pytest.raises(MODULE.ReleaseError, match="formal-release link"):
        MODULE.main()


def test_supplemental_release_rejects_audit_identity_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = create_release_fixture(tmp_path, monkeypatch)
    audit_summary = fixture["audit_summary"]
    assert isinstance(audit_summary, Path)
    value = json.loads(audit_summary.read_text(encoding="utf-8"))
    value["extra"] = "tampered"
    audit_summary.chmod(0o644)
    write_json(audit_summary, value)
    monkeypatch.setattr(sys, "argv", fixture["argv"])
    with pytest.raises(MODULE.ReleaseError, match="identity drift"):
        MODULE.main()


def test_supplemental_release_rejects_signature_tamper(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = create_release_fixture(tmp_path, monkeypatch)
    signature = fixture["dataset_signature"]
    assert isinstance(signature, Path)
    signature.chmod(0o644)
    signature.write_bytes(b"\x00" * 64)
    signature.chmod(0o444)
    monkeypatch.setattr(sys, "argv", fixture["argv"])
    with pytest.raises(MODULE.ReleaseError, match="detached signature"):
        MODULE.main()


def test_supplemental_release_rejects_public_key_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = create_release_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(
        MODULE, "SUPPLEMENTAL_PUBLIC_KEY_SHA256", "0" * 64
    )
    monkeypatch.setattr(sys, "argv", fixture["argv"])
    with pytest.raises(
        MODULE.ReleaseError, match="Supplemental integrity public-key drift"
    ):
        MODULE.main()


def test_supplemental_release_rejects_aba_even_after_path_is_restored(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = create_release_fixture(tmp_path, monkeypatch)
    audit_summary = fixture["audit_summary"]
    supplemental_public = fixture["supplemental_public"]
    output = fixture["output"]
    assert isinstance(audit_summary, Path)
    assert isinstance(supplemental_public, Path)
    assert isinstance(output, Path)
    good_backup = tmp_path / "audit_summary.good.json"
    bad_pending = tmp_path / "audit_summary.bad.pending.json"
    bad_after = tmp_path / "audit_summary.bad.after.json"
    write_json(
        bad_pending,
        {
            "schema_name": (
                "intersubmod.positional_singleton_methyl_multigroup_audit"
            ),
            "schema_version": "1.0.0",
            "pass": True,
            "checks": {"complete": True},
            "counts": {
                "singleton_sites": 999,
                "m1_evaluable": 999,
                "m1_flagged": 999,
                "m2_pass": 999,
            },
        },
    )
    original_open = MODULE.SOURCE_AUTHORITY.BoundArtifactReader.open
    state = {"switched": False, "restored": False}

    def injected_open(
        reader: object,
        path: Path,
        *,
        include_mode: bool = False,
    ) -> tuple[int, bytes, dict[str, object]]:
        result = original_open(reader, path, include_mode=include_mode)
        resolved = path.resolve()
        if resolved == audit_summary.resolve() and not state["switched"]:
            audit_summary.replace(good_backup)
            bad_pending.replace(audit_summary)
            state["switched"] = True
        elif (
            resolved == supplemental_public.resolve()
            and state["switched"]
            and not state["restored"]
        ):
            audit_summary.replace(bad_after)
            good_backup.replace(audit_summary)
            state["restored"] = True
        return result

    monkeypatch.setattr(
        MODULE.SOURCE_AUTHORITY.BoundArtifactReader,
        "open",
        injected_open,
    )
    monkeypatch.setattr(sys, "argv", fixture["argv"])
    with pytest.raises(
        MODULE.SOURCE_AUTHORITY.SourceAuthorityError,
        match="binding changed",
    ):
        MODULE.main()

    assert state == {"switched": True, "restored": True}
    assert not output.exists()


def run_post_signature_verification(
    fixture: dict[str, Path | list[str]],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> Path:
    output = fixture["output"]
    private_key = fixture["supplemental_private"]
    assert isinstance(output, Path)
    assert isinstance(private_key, Path)
    signature = sign_receipt(output, private_key)
    private_key.chmod(0o000)
    verification_output = tmp_path / "post_signature_verification.json"
    argv = list(fixture["argv"])
    argv.extend(
        [
            "--mode",
            "verify",
            "--supplemental-signature",
            str(signature),
            "--supplemental-private-key",
            str(private_key),
            "--verification-output",
            str(verification_output),
        ]
    )
    monkeypatch.setattr(sys, "argv", argv)
    MODULE.main()
    return verification_output


def test_post_signature_verification_recomputes_full_live_contract(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = create_release_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(sys, "argv", fixture["argv"])
    MODULE.main()

    verification_output = run_post_signature_verification(
        fixture, monkeypatch, tmp_path
    )
    verification = json.loads(
        verification_output.read_text(encoding="utf-8")
    )
    assert verification["pass"] is True
    assert (
        verification["verification_status"]
        == "SIGNED_SUPPLEMENTAL_RELEASE_REVERIFIED"
    )
    assert all(verification["checks"].values())
    assert verification["checks"][
        "all_bound_direct_input_identities_reverified"
    ] is True
    assert "all_transitive_input_identities_reverified" not in verification["checks"]
    assert verification["reverification_scope"] == {
        "identity_count": len(MODULE.INPUT_NAMES),
        "scope": "direct_receipt_inputs_only",
        "transitive_artifacts_cryptographically_bound_but_not_reopened": True,
    }
    assert stat.S_IMODE(verification_output.stat().st_mode) == 0o444


def test_post_signature_verification_rejects_tampered_receipt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = create_release_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(sys, "argv", fixture["argv"])
    MODULE.main()
    output = fixture["output"]
    private_key = fixture["supplemental_private"]
    assert isinstance(output, Path)
    assert isinstance(private_key, Path)
    signature = sign_receipt(output, private_key)
    private_key.chmod(0o000)
    output.chmod(0o644)
    value = json.loads(output.read_text(encoding="utf-8"))
    value["counts"]["m2_pass"] = 999
    write_json(output, value)
    verification_output = tmp_path / "post_signature_verification.json"
    argv = list(fixture["argv"])
    argv.extend(
        [
            "--mode",
            "verify",
            "--supplemental-signature",
            str(signature),
            "--supplemental-private-key",
            str(private_key),
            "--verification-output",
            str(verification_output),
        ]
    )
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(MODULE.ReleaseError, match="signature failed"):
        MODULE.main()


def test_post_signature_verification_rejects_unretired_private_key(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = create_release_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(sys, "argv", fixture["argv"])
    MODULE.main()
    output = fixture["output"]
    private_key = fixture["supplemental_private"]
    assert isinstance(output, Path)
    assert isinstance(private_key, Path)
    signature = sign_receipt(output, private_key)
    verification_output = tmp_path / "post_signature_verification.json"
    argv = list(fixture["argv"])
    argv.extend(
        [
            "--mode",
            "verify",
            "--supplemental-signature",
            str(signature),
            "--supplemental-private-key",
            str(private_key),
            "--verification-output",
            str(verification_output),
        ]
    )
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(MODULE.ReleaseError, match="not retired"):
        MODULE.main()
