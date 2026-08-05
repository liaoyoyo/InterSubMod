from __future__ import annotations

import ast
import importlib.util
import hashlib
import inspect
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import textwrap
from types import ModuleType

import pytest


TOPIC_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = TOPIC_ROOT / "audit_notes" / "run_one_time_ed25519_signer_v6.py"


def load_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "run_one_time_ed25519_signer_v6", MODULE_PATH
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


def freeze_target(path: Path) -> dict[str, object]:
    path.chmod(0o444)
    data = path.read_bytes()
    return {
        "path": str(path.resolve()),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": "0o444",
    }


def signing_instruction(record: dict[str, object]) -> str:
    return "SIGN " + json.dumps(
        record,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def run_cli(
    *,
    key_dir: Path,
    target: Path,
    command_input: str,
) -> subprocess.CompletedProcess[str]:
    python_fd = os.open(sys.executable, os.O_RDONLY | os.O_CLOEXEC)
    script_fd = os.open(MODULE_PATH, os.O_RDONLY | os.O_CLOEXEC)
    try:
        return subprocess.run(
            [
                f"/proc/self/fd/{python_fd}",
                "-I",
                f"/proc/self/fd/{script_fd}",
                str(key_dir),
                "ed25519_private.pem",
                str(target),
            ],
            input=command_input,
            text=True,
            capture_output=True,
            check=False,
            pass_fds=(python_fd, script_fd),
            env={
                "PATH": "/usr/bin:/bin",
                "HOME": "/tmp",
                "LANG": "C.UTF-8",
                "LC_ALL": "C.UTF-8",
                "PYTHONHASHSEED": "0",
                "PYTHONNOUSERSITE": "1",
            },
        )
    finally:
        os.close(script_fd)
        os.close(python_fd)


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
    target_record = freeze_target(target)
    key_dir = tmp_path / "keys"

    completed = run_cli(
        key_dir=key_dir,
        target=target,
        command_input=signing_instruction(target_record) + "\n",
    )

    assert completed.returncode == 0, completed.stderr
    records = [json.loads(line) for line in completed.stdout.splitlines()]
    assert [record["event"] for record in records] == [
        "SIGNER_READY",
        "CEREMONY_OUTPUT_AVAILABLE_REQUIRES_INDEPENDENT_VERIFICATION",
    ]
    assert "grants no release authority" in records[0]["publication_semantics"]
    assert records[1]["point_in_time_checks"] == (
        "PASS_BEFORE_DESCRIPTOR_TEARDOWN"
    )
    assert records[1]["consumer_verification_required"] is True
    assert records[1]["release_authority_granted"] is False
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
    target_record = freeze_target(target)
    signature = Path(str(target) + ".ed25519.sig")
    signature.write_bytes(b"preexisting")
    original_signature = signature.read_bytes()
    key_dir = tmp_path / "keys"

    completed = run_cli(
        key_dir=key_dir,
        target=target,
        command_input=signing_instruction(target_record) + "\n",
    )

    assert completed.returncode == 8
    failure = json.loads(completed.stderr)
    assert failure["event"] == "SIGN_FAILED"
    assert failure["private_key_state"] == "RETAINED_MODE_400"
    assert failure["signature_published_by_signer"] is False
    assert signature.read_bytes() == original_signature
    assert file_mode(key_dir / "ed25519_private.pem") == 0o400


def test_publication_race_collision_is_no_replace(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
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
            signer.sign(target_record)
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
    target_record = freeze_target(target)
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
            target.chmod(0o444)

    monkeypatch.setattr(signer, "_verify_signature", verify_then_replace_target)
    try:
        with pytest.raises(module.SignerError, match="Path binding changed"):
            signer.sign(target_record)
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
    target_record = freeze_target(target)
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
            signer.sign(target_record)
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
    target_record = freeze_target(target)
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
            signer.sign(target_record)
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
    target_record = freeze_target(target)
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
            signer.sign(target_record)
        signature = Path(str(target) + ".ed25519.sig")
        assert signature.is_file()
        assert file_mode(key_dir / "ed25519_private.pem") == 0o400
        verify_signature(target, key_dir / "ed25519_public.pem", signature)
    finally:
        signer.close()


def test_wrong_command_then_eof_does_not_sign_or_retire(tmp_path: Path) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    freeze_target(target)
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
    target_record = freeze_target(target)
    key_dir = tmp_path / "keys"
    key_dir.mkdir(mode=0o700)
    sentinel = key_dir / "sentinel"
    sentinel.write_text("unchanged\n", encoding="utf-8")

    completed = run_cli(
        key_dir=key_dir,
        target=target,
        command_input=signing_instruction(target_record) + "\n",
    )

    assert completed.returncode == 3
    assert "SIGNER_SETUP_FAILED" in completed.stderr
    assert sentinel.read_text(encoding="utf-8") == "unchanged\n"
    assert sorted(path.name for path in key_dir.iterdir()) == ["sentinel"]


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("sha256", "0" * 64, "does not match"),
        ("size_bytes", 999, "does not match"),
        ("mode", "0o644", "mode must be 0o444"),
    ],
)
def test_wrong_approved_identity_never_signs_or_retires(
    module: ModuleType,
    tmp_path: Path,
    field: str,
    value: object,
    message: str,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    target_record[field] = value
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    try:
        with pytest.raises(module.SignerError, match=message):
            signer.sign(target_record)
        assert not Path(str(target) + ".ed25519.sig").exists()
        assert file_mode(key_dir / "ed25519_private.pem") == 0o400
    finally:
        signer.close()


def test_target_replacement_after_approval_identity_is_rejected(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    displaced = tmp_path / "authority.approved.json"
    target.rename(displaced)
    target.write_text('{"authorized":false}\n', encoding="utf-8")
    target.chmod(0o444)
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    try:
        with pytest.raises(module.SignerError, match="does not match"):
            signer.sign(target_record)
        assert not Path(str(target) + ".ed25519.sig").exists()
        assert file_mode(key_dir / "ed25519_private.pem") == 0o400
    finally:
        signer.close()


def test_key_directory_chmod_after_ready_blocks_signing(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    signer.ready_record()
    key_dir.chmod(0o775)
    try:
        with pytest.raises(module.SignerError, match="Directory binding changed"):
            signer.sign(target_record)
        assert not Path(str(target) + ".ed25519.sig").exists()
        assert file_mode(key_dir / "ed25519_private.pem") == 0o400
    finally:
        key_dir.chmod(0o700)
        signer.close()


def test_private_fsync_failure_after_chmod_reports_retired_state(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_fsync = module.os.fsync

    def fail_retired_private_fsync(fd: int) -> None:
        if fd == signer.private_fd and file_mode(key_dir / "ed25519_private.pem") == 0:
            raise OSError("injected retired-private fsync failure")
        real_fsync(fd)

    monkeypatch.setattr(module.os, "fsync", fail_retired_private_fsync)
    try:
        with pytest.raises(OSError, match="injected retired-private fsync failure") as caught:
            signer.sign(target_record)
        failure = signer.failure_record(caught.value)
        signature = Path(str(target) + ".ed25519.sig")
        assert failure["signature_published_by_signer"] is True
        assert failure["signature_path_exists"] is True
        assert failure["private_retirement_started"] is True
        assert failure["private_key_state"] == "RETIRED_MODE_000"
        assert failure["private_key_mode"] == "0o0"
        assert signature.is_file()
        verify_signature(target, key_dir / "ed25519_public.pem", signature)
    finally:
        signer.close()


def test_signature_path_replacement_after_retirement_never_reports_signed(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    signature = Path(str(target) + ".ed25519.sig")
    displaced = Path(str(signature) + ".displaced")
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_verify = signer._verify_signature
    calls = 0

    def verify_then_replace_published_path(
        *,
        target_fd: int,
        signature_fd: int,
    ) -> None:
        nonlocal calls
        real_verify(target_fd=target_fd, signature_fd=signature_fd)
        calls += 1
        if calls == 3:
            signature.rename(displaced)
            signature.write_bytes(b"x" * 64)
            signature.chmod(0o444)

    monkeypatch.setattr(
        signer,
        "_verify_signature",
        verify_then_replace_published_path,
    )
    try:
        with pytest.raises(module.SignerError, match="Path binding changed") as caught:
            signer.sign(target_record)
        failure = signer.failure_record(caught.value)
        assert calls == 3
        assert failure["event"] == "SIGN_FAILED"
        assert failure["signature_published_by_signer"] is True
        assert failure["private_retirement_started"] is True
        assert failure["private_key_state"] == "RETIRED_MODE_000"
        assert failure["private_key_mode"] == "0o0"
        assert signature.read_bytes() == b"x" * 64
        verify_signature(target, key_dir / "ed25519_public.pem", displaced)
    finally:
        signer.close()


def test_target_path_replacement_after_terminal_verify_never_reports_signed(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    displaced = tmp_path / "authority.approved.json"
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_verify = signer._verify_signature
    calls = 0

    def verify_then_replace_target(
        *,
        target_fd: int,
        signature_fd: int,
    ) -> None:
        nonlocal calls
        real_verify(target_fd=target_fd, signature_fd=signature_fd)
        calls += 1
        if calls == 3:
            target.rename(displaced)
            target.write_text('{"authorized":false}\n', encoding="utf-8")
            target.chmod(0o444)

    monkeypatch.setattr(signer, "_verify_signature", verify_then_replace_target)
    try:
        with pytest.raises(module.SignerError, match="Path binding changed") as caught:
            signer.sign(target_record)
        failure = signer.failure_record(caught.value)
        assert calls == 3
        assert failure["event"] == "SIGN_FAILED"
        assert failure["signature_published_by_signer"] is True
        assert failure["private_key_state"] == "RETIRED_MODE_000"
    finally:
        signer.close()


def test_same_inode_signature_pwrite_after_terminal_verify_never_reports_signed(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    signature = Path(str(target) + ".ed25519.sig")
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_verify = signer._verify_signature
    real_link = module.link_fd_no_replace
    calls = 0
    writable_signature_fd = -1

    def capture_writable_fd_then_link(
        source_fd: int,
        target_directory_fd: int,
        target_name: str,
    ) -> None:
        nonlocal writable_signature_fd
        writable_signature_fd = os.dup(source_fd)
        real_link(source_fd, target_directory_fd, target_name)

    def verify_then_pwrite_signature(
        *,
        target_fd: int,
        signature_fd: int,
    ) -> None:
        nonlocal calls
        real_verify(target_fd=target_fd, signature_fd=signature_fd)
        calls += 1
        if calls == 3:
            assert writable_signature_fd >= 0
            os.pwrite(writable_signature_fd, b"x" * 64, 0)
            os.fsync(writable_signature_fd)

    monkeypatch.setattr(module, "link_fd_no_replace", capture_writable_fd_then_link)
    monkeypatch.setattr(signer, "_verify_signature", verify_then_pwrite_signature)
    try:
        with pytest.raises(
            module.SignerError,
            match=(
                "Opened input changed during signing|"
                "Signature identity changed during terminal"
            ),
        ) as caught:
            signer.sign(target_record)
        failure = signer.failure_record(caught.value)
        assert calls == 3
        assert failure["event"] == "SIGN_FAILED"
        assert failure["signature_published_by_signer"] is True
        assert failure["private_key_state"] == "RETIRED_MODE_000"
        assert signature.read_bytes() == b"x" * 64
    finally:
        if writable_signature_fd >= 0:
            os.close(writable_signature_fd)
        signer.close()


def test_post_retirement_key_directory_fsync_replacement_is_revalidated(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    signature = Path(str(target) + ".ed25519.sig")
    displaced = Path(str(signature) + ".displaced")
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_fsync = module.os.fsync
    injected = False

    def fsync_then_replace_signature(fd: int) -> None:
        nonlocal injected
        real_fsync(fd)
        if (
            not injected
            and fd == signer.key_directory_fd
            and file_mode(key_dir / "ed25519_private.pem") == 0o000
        ):
            injected = True
            signature.rename(displaced)
            signature.write_bytes(b"x" * 64)
            signature.chmod(0o444)

    monkeypatch.setattr(module.os, "fsync", fsync_then_replace_signature)
    try:
        with pytest.raises(module.SignerError, match="Path binding changed") as caught:
            signer.sign(target_record)
        failure = signer.failure_record(caught.value)
        assert injected is True
        assert failure["event"] == "SIGN_FAILED"
        assert failure["signature_published_by_signer"] is True
        assert failure["private_key_state"] == "RETIRED_MODE_000"
        verify_signature(target, key_dir / "ed25519_public.pem", displaced)
    finally:
        signer.close()


def test_point_in_time_validation_is_last_pre_teardown_operation(
    module: ModuleType,
) -> None:
    tree = ast.parse(textwrap.dedent(inspect.getsource(module.OneTimeSigner.sign)))
    function = tree.body[0]
    assert isinstance(function, ast.FunctionDef)
    try_statement = next(
        statement for statement in function.body if isinstance(statement, ast.Try)
    )
    validation_index = next(
        index
        for index, statement in enumerate(try_statement.body)
        if isinstance(statement, (ast.Assign, ast.AnnAssign))
        and any(
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "_point_in_time_validate_published_state"
            for node in ast.walk(statement)
        )
    )
    assert len(try_statement.body) == validation_index + 2
    provisional_return = try_statement.body[validation_index + 1]
    assert isinstance(provisional_return, ast.Return)
    assert not any(isinstance(node, ast.Call) for node in ast.walk(provisional_return))
    assert not any(
        isinstance(node, ast.Constant) and node.value == "SIGNED"
        for node in ast.walk(ast.parse(MODULE_PATH.read_text(encoding="utf-8")))
    )


def test_teardown_race_requires_independent_consumer_rejection(
    module: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    signature = Path(str(target) + ".ed25519.sig")
    displaced = Path(str(signature) + ".displaced")
    key_dir = tmp_path / "keys"
    signer = module.OneTimeSigner(
        key_dir=key_dir,
        private_filename="ed25519_private.pem",
        expected_target=target,
    )
    real_close = module.os.close
    injected = False

    def close_then_replace_live_signature(fd: int) -> None:
        nonlocal injected
        closed_signature_inode = False
        if not injected and signature.is_file():
            try:
                closed_signature_inode = (
                    os.fstat(fd).st_dev,
                    os.fstat(fd).st_ino,
                ) == (
                    signature.stat().st_dev,
                    signature.stat().st_ino,
                )
            except OSError:
                pass
        real_close(fd)
        if (
            not injected
            and closed_signature_inode
            and signature.is_file()
            and file_mode(key_dir / "ed25519_private.pem") == 0o000
        ):
            injected = True
            signature.rename(displaced)
            signature.write_bytes(b"x" * 64)
            signature.chmod(0o444)

    monkeypatch.setattr(module.os, "close", close_then_replace_live_signature)
    try:
        result = signer.sign(target_record)
        assert injected is True
        assert result["event"] == (
            "CEREMONY_OUTPUT_AVAILABLE_REQUIRES_INDEPENDENT_VERIFICATION"
        )
        assert result["consumer_verification_required"] is True
        assert result["release_authority_granted"] is False
        completed = subprocess.run(
            [
                "/usr/bin/openssl",
                "pkeyutl",
                "-verify",
                "-rawin",
                "-pubin",
                "-inkey",
                str(key_dir / "ed25519_public.pem"),
                "-sigfile",
                str(signature),
                "-in",
                str(target),
            ],
            capture_output=True,
            check=False,
        )
        assert completed.returncode != 0
        verify_signature(target, key_dir / "ed25519_public.pem", displaced)
    finally:
        monkeypatch.setattr(module.os, "close", real_close)
        signer.close()


def test_noncanonical_signing_instruction_is_rejected(
    module: ModuleType,
    tmp_path: Path,
) -> None:
    target = tmp_path / "authority.json"
    target.write_text('{"authorized":true}\n', encoding="utf-8")
    target_record = freeze_target(target)
    noncanonical = "SIGN " + json.dumps(target_record, sort_keys=False)

    with pytest.raises(module.SignerError, match="not canonical"):
        module.parse_signing_instruction(
            noncanonical,
            expected_path=target.resolve(),
        )
