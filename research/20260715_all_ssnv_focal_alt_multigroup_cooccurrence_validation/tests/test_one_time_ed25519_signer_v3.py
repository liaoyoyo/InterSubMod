from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = TOPIC_ROOT / "audit_notes" / "run_one_time_ed25519_signer_v3.py"


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "run_one_time_ed25519_signer_v3", MODULE_PATH
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def module() -> ModuleType:
    return load_module()


def file_mode(path: Path) -> int:
    return stat.S_IMODE(path.stat().st_mode)


def run_cli(
    *,
    key_dir: Path,
    target: Path,
    command_input: str,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "-I",
            str(MODULE_PATH),
            str(key_dir),
            "ed25519_private.pem",
            str(target),
        ],
        input=command_input,
        text=True,
        capture_output=True,
        check=False,
    )


def verify_signature(target: Path, public_key: Path, signature: Path) -> None:
    completed = subprocess.run(
        [
            "/usr/bin/openssl",
            "pkeyutl",
            "-verify",
            "-rawin",
            "-pubin",
            "-inkey",
            str(public_key),
            "-sigfile",
            str(signature),
            "-in",
            str(target),
        ],
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr.decode(
        "utf-8", errors="replace"
    )


def test_success_uses_bound_key_and_retires_only_after_publication(
    tmp_path: Path,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    key_dir = tmp_path / "keys"

    completed = run_cli(
        key_dir=key_dir,
        target=target,
        command_input=f"SIGN {target.resolve()}\n",
    )

    assert completed.returncode == 0, completed.stderr
    records = [json.loads(line) for line in completed.stdout.splitlines()]
    assert [record["event"] for record in records] == ["SIGNER_READY", "SIGNED"]
    assert records[0]["openssl"]["sha256"] == (
        "38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08"
    )
    private_key = key_dir / "ed25519_private.pem"
    public_key = key_dir / "ed25519_public.pem"
    signature = Path(str(target) + ".ed25519.sig")
    assert file_mode(key_dir) == 0o700
    assert file_mode(private_key) == 0o000
    assert file_mode(public_key) == 0o444
    assert file_mode(signature) == 0o444
    assert signature.stat().st_size == 64
    assert not Path(str(signature) + ".pending").exists()
    verify_signature(target, public_key, signature)


def test_existing_signature_is_never_overwritten_and_private_key_stays_live(
    tmp_path: Path,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    signature = Path(str(target) + ".ed25519.sig")
    signature.write_bytes(b"preexisting")
    original_signature = signature.read_bytes()
    key_dir = tmp_path / "keys"

    completed = run_cli(
        key_dir=key_dir,
        target=target,
        command_input=f"SIGN {target.resolve()}\n",
    )

    assert completed.returncode == 8
    assert "SIGN_FAILED_PRIVATE_KEY_RETAINED" in completed.stderr
    assert signature.read_bytes() == original_signature
    assert file_mode(key_dir / "ed25519_private.pem") == 0o400


def test_publication_race_collision_is_no_replace(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    signature = Path(str(target) + ".ed25519.sig")
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_verify = signer._verify_signature
    calls = 0

    def verify_then_create_collision(*, target_fd: int, signature_fd: int) -> None:
        nonlocal calls
        real_verify(target_fd=target_fd, signature_fd=signature_fd)
        calls += 1
        if calls == 1:
            signature.write_bytes(b"racing-output")

    monkeypatch.setattr(signer, "_verify_signature", verify_then_create_collision)
    try:
        with pytest.raises(module.SignerError, match="File exists"):
            signer.sign()
        assert signature.read_bytes() == b"racing-output"
        assert file_mode(key_dir / "ed25519_private.pem") == 0o400
    finally:
        signer.close()


def test_target_replacement_before_publication_is_rejected(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    displaced = tmp_path / "authority.original.json"
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_verify = signer._verify_signature
    calls = 0

    def verify_then_replace_target(*, target_fd: int, signature_fd: int) -> None:
        nonlocal calls
        real_verify(target_fd=target_fd, signature_fd=signature_fd)
        calls += 1
        if calls == 1:
            os.rename(target, displaced)
            target.write_text('{"authorized":false}\n', encoding="utf-8")

    monkeypatch.setattr(signer, "_verify_signature", verify_then_replace_target)
    try:
        with pytest.raises(module.SignerError, match="Path binding changed"):
            signer.sign()
        assert not Path(str(target) + ".ed25519.sig").exists()
        assert file_mode(key_dir / "ed25519_private.pem") == 0o400
    finally:
        signer.close()


def test_retirement_failure_keeps_private_key_live_after_valid_publication(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )

    def fail_retirement() -> None:
        raise module.SignerError("injected retirement failure")

    monkeypatch.setattr(signer, "_retire_private_key", fail_retirement)
    try:
        with pytest.raises(module.SignerError, match="injected retirement failure"):
            signer.sign()
        signature = Path(str(target) + ".ed25519.sig")
        assert signature.is_file()
        assert file_mode(signature) == 0o444
        assert file_mode(key_dir / "ed25519_private.pem") == 0o400
        verify_signature(target, key_dir / "ed25519_public.pem", signature)
    finally:
        signer.close()


def test_signature_chmod_failure_prevents_publication_and_retirement(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_fchmod = module.os.fchmod

    def fail_signature_chmod(fd: int, mode: int) -> None:
        if mode == 0o444 and fd not in {signer.public_fd}:
            raise OSError("injected signature chmod failure")
        real_fchmod(fd, mode)

    monkeypatch.setattr(module.os, "fchmod", fail_signature_chmod)
    try:
        with pytest.raises(OSError, match="injected signature chmod failure"):
            signer.sign()
        assert not Path(str(target) + ".ed25519.sig").exists()
        assert file_mode(key_dir / "ed25519_private.pem") == 0o400
    finally:
        signer.close()


def test_directory_fsync_failure_after_publication_is_not_reported_as_success(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_fsync = module.os.fsync

    def fail_directory_fsync(fd: int) -> None:
        if stat.S_ISDIR(os.fstat(fd).st_mode):
            raise OSError("injected directory fsync failure")
        real_fsync(fd)

    monkeypatch.setattr(module.os, "fsync", fail_directory_fsync)
    try:
        with pytest.raises(OSError, match="injected directory fsync failure"):
            signer.sign()
        signature = Path(str(target) + ".ed25519.sig")
        assert signature.is_file()
        assert file_mode(key_dir / "ed25519_private.pem") == 0o400
        verify_signature(target, key_dir / "ed25519_public.pem", signature)
    finally:
        signer.close()


def test_wrong_command_then_eof_does_not_sign_or_retire(tmp_path: Path) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    key_dir = tmp_path / "keys"

    completed = run_cli(
        key_dir=key_dir,
        target=target,
        command_input="SIGN wrong-target\n",
    )

    assert completed.returncode == 7
    assert "WAITING_FOR_EXACT_SIGN_TARGET" in completed.stdout
    assert "SIGNER_INPUT_CLOSED_WITHOUT_SIGNATURE" in completed.stderr
    assert not Path(str(target) + ".ed25519.sig").exists()
    assert file_mode(key_dir / "ed25519_private.pem") == 0o400


def test_existing_key_directory_is_a_setup_failure_without_mutation(
    tmp_path: Path,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    key_dir = tmp_path / "keys"
    key_dir.mkdir(mode=0o700)
    sentinel = key_dir / "sentinel"
    sentinel.write_text("unchanged\n", encoding="utf-8")

    completed = run_cli(
        key_dir=key_dir,
        target=target,
        command_input=f"SIGN {target.resolve()}\n",
    )

    assert completed.returncode == 3
    assert "SIGNER_SETUP_FAILED" in completed.stderr
    assert sentinel.read_text(encoding="utf-8") == "unchanged\n"
    assert sorted(path.name for path in key_dir.iterdir()) == ["sentinel"]
