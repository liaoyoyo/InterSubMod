#!/usr/bin/env python3
"""Publish one digest-bound Ed25519 signature for independent verification."""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import secrets
import stat
import subprocess
import sys
from typing import Any, Mapping, Sequence


OPENSSL = Path("/usr/bin/openssl")
OPENSSL_SHA256 = "38c064f53a6619364c7947fc11cad09e1fc4218fc2c3e9016a7786b1341ecd08"
PYTHON = Path("/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python")
PYTHON_SHA256 = "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e"
BOOTSTRAP_ENVIRONMENT = {
    "PATH": "/usr/bin:/bin",
    "HOME": "/tmp",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PYTHONHASHSEED": "0",
    "PYTHONNOUSERSITE": "1",
}
PUBLIC_FILENAME = "ed25519_public.pem"
AT_FDCWD = -100
AT_SYMLINK_FOLLOW = 0x400


class SignerError(RuntimeError):
    """Raised when a one-time signing ceremony cannot complete safely."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_fd(fd: int) -> bytes:
    opened = os.fstat(fd)
    if not stat.S_ISREG(opened.st_mode):
        raise SignerError("Descriptor is not a regular file")
    chunks: list[bytes] = []
    offset = 0
    while offset < opened.st_size:
        chunk = os.pread(fd, min(1024 * 1024, opened.st_size - offset), offset)
        if not chunk:
            break
        chunks.append(chunk)
        offset += len(chunk)
    data = b"".join(chunks)
    closed_read = os.fstat(fd)
    if (
        len(data) != opened.st_size
        or opened.st_dev != closed_read.st_dev
        or opened.st_ino != closed_read.st_ino
        or opened.st_size != closed_read.st_size
        or opened.st_mtime_ns != closed_read.st_mtime_ns
    ):
        raise SignerError("Descriptor changed while reading")
    return data


def open_bound_regular(path: Path, *, writable: bool = False) -> tuple[int, os.stat_result]:
    flags = (os.O_RDWR if writable else os.O_RDONLY) | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags)
    opened = os.fstat(fd)
    if not stat.S_ISREG(opened.st_mode):
        os.close(fd)
        raise SignerError(f"Not a regular file: {path}")
    live = path.stat()
    if live.st_dev != opened.st_dev or live.st_ino != opened.st_ino:
        os.close(fd)
        raise SignerError(f"Path is not bound to the opened file: {path}")
    return fd, opened


def open_bound_directory(path: Path) -> tuple[int, os.stat_result]:
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0)
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags)
    opened = os.fstat(fd)
    if not stat.S_ISDIR(opened.st_mode):
        os.close(fd)
        raise SignerError(f"Not a directory: {path}")
    live = path.stat()
    if live.st_dev != opened.st_dev or live.st_ino != opened.st_ino:
        os.close(fd)
        raise SignerError(f"Directory path is not bound to the opened directory: {path}")
    return fd, opened


def require_path_bound(
    path: Path,
    fd: int,
    *,
    opened: os.stat_result | None = None,
    expected_mode: int | None = None,
) -> os.stat_result:
    descriptor_stat = os.fstat(fd)
    live = path.stat()
    if (
        descriptor_stat.st_dev != live.st_dev
        or descriptor_stat.st_ino != live.st_ino
        or descriptor_stat.st_mode != live.st_mode
        or descriptor_stat.st_size != live.st_size
        or descriptor_stat.st_mtime_ns != live.st_mtime_ns
        or descriptor_stat.st_ctime_ns != live.st_ctime_ns
    ):
        raise SignerError(f"Path binding changed: {path}")
    if opened is not None and (
        descriptor_stat.st_dev != opened.st_dev
        or descriptor_stat.st_ino != opened.st_ino
        or descriptor_stat.st_mode != opened.st_mode
        or descriptor_stat.st_size != opened.st_size
        or descriptor_stat.st_mtime_ns != opened.st_mtime_ns
        or descriptor_stat.st_ctime_ns != opened.st_ctime_ns
    ):
        raise SignerError(f"Opened input changed during signing: {path}")
    if expected_mode is not None and stat.S_IMODE(descriptor_stat.st_mode) != expected_mode:
        raise SignerError(
            f"Unexpected mode for {path}: {oct(stat.S_IMODE(descriptor_stat.st_mode))}"
        )
    return descriptor_stat


def require_directory_bound(
    path: Path,
    fd: int,
    *,
    opened: os.stat_result,
    expected_mode: int | None = None,
) -> os.stat_result:
    descriptor_stat = os.fstat(fd)
    live = path.stat()
    if (
        descriptor_stat.st_dev != live.st_dev
        or descriptor_stat.st_ino != live.st_ino
        or descriptor_stat.st_mode != live.st_mode
        or descriptor_stat.st_nlink != live.st_nlink
        or descriptor_stat.st_dev != opened.st_dev
        or descriptor_stat.st_ino != opened.st_ino
        or descriptor_stat.st_mode != opened.st_mode
        or descriptor_stat.st_nlink != opened.st_nlink
    ):
        raise SignerError(f"Directory binding changed: {path}")
    if expected_mode is not None and stat.S_IMODE(descriptor_stat.st_mode) != expected_mode:
        raise SignerError(
            f"Unexpected directory mode for {path}: "
            f"{oct(stat.S_IMODE(descriptor_stat.st_mode))}"
        )
    return descriptor_stat


def create_exclusive_file(path: Path, mode: int) -> int:
    flags = os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    return os.open(path, flags, mode)


def create_anonymous_file(directory_fd: int, mode: int) -> int:
    if not hasattr(os, "O_TMPFILE"):
        raise SignerError("O_TMPFILE is required for anonymous signature creation")
    flags = os.O_RDWR | os.O_TMPFILE | os.O_CLOEXEC
    return os.open(".", flags, mode, dir_fd=directory_fd)


def write_all(fd: int, data: bytes) -> None:
    offset = 0
    while offset < len(data):
        written = os.write(fd, data[offset:])
        if written <= 0:
            raise SignerError("Unable to complete exclusive file write")
        offset += written
    os.fsync(fd)


def fsync_directory(path: Path) -> None:
    fd = os.open(path, os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def link_fd_no_replace(source_fd: int, target_directory_fd: int, target_name: str) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    linkat = getattr(libc, "linkat", None)
    if linkat is None:
        raise SignerError("linkat is required for FD-bound no-replace publication")
    linkat.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
    )
    linkat.restype = ctypes.c_int
    result = linkat(
        AT_FDCWD,
        os.fsencode(f"/proc/self/fd/{source_fd}"),
        target_directory_fd,
        os.fsencode(target_name),
        AT_SYMLINK_FOLLOW,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise SignerError(
            f"FD-bound no-replace publication failed: {os.strerror(error_number)}"
        )


def minimal_openssl_environment(passphrase: str) -> dict[str, str]:
    environment = {
        key: value
        for key, value in os.environ.items()
        if key
        not in {
            "CONDA_PREFIX",
            "LD_LIBRARY_PATH",
            "OPENSSL_CONF",
            "OPENSSL_MODULES",
            "PYTHONPATH",
        }
    }
    environment["PATH"] = "/usr/bin:/bin"
    environment["PASS"] = passphrase
    return environment


def run_openssl(
    arguments: Sequence[str],
    *,
    executable_fd: int,
    environment: dict[str, str],
    pass_fds: Sequence[int] = (),
) -> bytes:
    for fd in pass_fds:
        os.lseek(fd, 0, os.SEEK_SET)
    inherited_fds = tuple(dict.fromkeys((executable_fd, *pass_fds)))
    completed = subprocess.run(
        [f"/proc/self/fd/{executable_fd}", *arguments],
        check=False,
        capture_output=True,
        env=environment,
        pass_fds=inherited_fds,
    )
    if completed.returncode != 0:
        raise SignerError(
            f"OpenSSL command failed ({completed.returncode}): "
            f"{completed.stderr.decode('utf-8', errors='replace').strip()}"
        )
    return completed.stdout


def artifact_from_fd(
    path: Path,
    fd: int,
    *,
    opened: os.stat_result | None = None,
) -> dict[str, Any]:
    data = read_fd(fd)
    observed = require_path_bound(path, fd, opened=opened)
    return {
        "path": str(path),
        "size_bytes": observed.st_size,
        "sha256": sha256_bytes(data),
        "mode": oct(stat.S_IMODE(observed.st_mode)),
    }


def validate_expected_target_record(
    record: Mapping[str, Any],
    *,
    expected_path: Path,
) -> dict[str, Any]:
    expected_keys = {"path", "size_bytes", "sha256", "mode"}
    if set(record) != expected_keys:
        raise SignerError("Signing instruction target identity keys are not exact")
    path = record["path"]
    size_bytes = record["size_bytes"]
    sha256 = record["sha256"]
    mode = record["mode"]
    if (
        not isinstance(path, str)
        or Path(path).expanduser().resolve(strict=False) != expected_path
    ):
        raise SignerError("Signing instruction targets the wrong path")
    if not isinstance(size_bytes, int) or isinstance(size_bytes, bool) or size_bytes < 0:
        raise SignerError("Signing instruction size_bytes is invalid")
    if (
        not isinstance(sha256, str)
        or len(sha256) != 64
        or any(character not in "0123456789abcdef" for character in sha256)
    ):
        raise SignerError("Signing instruction sha256 is invalid")
    if mode != "0o444":
        raise SignerError("Signing instruction mode must be 0o444")
    return {
        "path": str(expected_path),
        "size_bytes": size_bytes,
        "sha256": sha256,
        "mode": mode,
    }


def format_signing_instruction(target: Mapping[str, Any]) -> str:
    return "SIGN " + json.dumps(
        dict(target),
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def parse_signing_instruction(line: str, *, expected_path: Path) -> dict[str, Any]:
    if not line.startswith("SIGN "):
        raise SignerError("Signing instruction does not start with SIGN")
    raw_payload = line[len("SIGN ") :]
    try:
        payload = json.loads(raw_payload)
    except json.JSONDecodeError as exc:
        raise SignerError("Signing instruction JSON is invalid") from exc
    if not isinstance(payload, dict):
        raise SignerError("Signing instruction target identity must be an object")
    target = validate_expected_target_record(payload, expected_path=expected_path)
    if format_signing_instruction(target) != line:
        raise SignerError("Signing instruction is not canonical")
    return target


class OneTimeSigner:
    """Publish with one encrypted key pair, then irreversibly retire that key."""

    def __init__(
        self,
        *,
        key_dir: Path,
        private_filename: str,
        expected_target: Path,
    ) -> None:
        if Path(private_filename).name != private_filename:
            raise SignerError("Private key filename must not contain a path")
        self.key_dir = key_dir.expanduser().resolve(strict=False)
        self.private_path = self.key_dir / private_filename
        self.public_path = self.key_dir / PUBLIC_FILENAME
        self.expected_target = expected_target.expanduser().resolve(strict=False)
        self.signature_path = Path(str(self.expected_target) + ".ed25519.sig")
        self.pending_path = Path(str(self.signature_path) + ".pending")
        self._passphrase = secrets.token_hex(48)
        self._environment = minimal_openssl_environment(self._passphrase)
        self.openssl_fd = -1
        self.openssl_opened: os.stat_result | None = None
        self.key_directory_fd = -1
        self.key_directory_opened: os.stat_result | None = None
        self.private_fd = -1
        self.private_opened: os.stat_result | None = None
        self.public_fd = -1
        self.public_opened: os.stat_result | None = None
        self.signature_published = False
        self.private_retirement_started = False
        try:
            self.openssl_fd, self.openssl_opened = open_bound_regular(OPENSSL)
            if sha256_bytes(read_fd(self.openssl_fd)) != OPENSSL_SHA256:
                raise SignerError("Pinned OpenSSL binary identity drifted")
            require_path_bound(
                OPENSSL,
                self.openssl_fd,
                opened=self.openssl_opened,
            )
            self._create_key_material()
        except Exception:
            self.close()
            raise

    def _run_openssl(
        self,
        arguments: Sequence[str],
        *,
        pass_fds: Sequence[int] = (),
    ) -> bytes:
        if self.openssl_fd < 0 or self.openssl_opened is None:
            raise SignerError("Pinned OpenSSL descriptor is unavailable")
        require_path_bound(
            OPENSSL,
            self.openssl_fd,
            opened=self.openssl_opened,
        )
        self._require_key_directory()
        output = run_openssl(
            arguments,
            executable_fd=self.openssl_fd,
            environment=self._environment,
            pass_fds=pass_fds,
        )
        require_path_bound(
            OPENSSL,
            self.openssl_fd,
            opened=self.openssl_opened,
        )
        self._require_key_directory()
        return output

    def _require_key_directory(self) -> os.stat_result:
        if self.key_directory_fd < 0 or self.key_directory_opened is None:
            raise SignerError("Bound key directory descriptor is unavailable")
        return require_directory_bound(
            self.key_dir,
            self.key_directory_fd,
            opened=self.key_directory_opened,
            expected_mode=0o700,
        )

    def _require_live_key_material(self) -> None:
        self._require_key_directory()
        if self.private_opened is None or self.public_opened is None:
            raise SignerError("Bound key material identity is unavailable")
        require_path_bound(
            self.private_path,
            self.private_fd,
            opened=self.private_opened,
            expected_mode=0o400,
        )
        require_path_bound(
            self.public_path,
            self.public_fd,
            opened=self.public_opened,
            expected_mode=0o444,
        )

    def _create_key_material(self) -> None:
        if not self.key_dir.parent.is_dir():
            raise SignerError(f"Key parent directory is missing: {self.key_dir.parent}")
        if os.path.lexists(self.key_dir):
            raise SignerError(f"Refusing existing key directory: {self.key_dir}")
        old_umask = os.umask(0o077)
        try:
            os.mkdir(self.key_dir, 0o700)
        finally:
            os.umask(old_umask)
        self.key_directory_fd, self.key_directory_opened = open_bound_directory(
            self.key_dir
        )
        self._require_key_directory()

        private_pem = self._run_openssl(
            [
                "genpkey",
                "-algorithm",
                "ED25519",
                "-aes-256-cbc",
                "-pass",
                "env:PASS",
            ],
        )
        if not private_pem:
            raise SignerError("OpenSSL produced an empty encrypted private key")
        self.private_fd = create_exclusive_file(self.private_path, 0o400)
        write_all(self.private_fd, private_pem)
        os.fchmod(self.private_fd, 0o400)
        os.fsync(self.private_fd)
        self.private_opened = os.fstat(self.private_fd)

        public_pem = self._run_openssl(
            [
                "pkey",
                "-in",
                f"/proc/self/fd/{self.private_fd}",
                "-passin",
                "env:PASS",
                "-pubout",
            ],
            pass_fds=(self.private_fd,),
        )
        if not public_pem:
            raise SignerError("OpenSSL produced an empty public key")
        self.public_fd = create_exclusive_file(self.public_path, 0o444)
        write_all(self.public_fd, public_pem)
        os.fchmod(self.public_fd, 0o444)
        os.fsync(self.public_fd)
        self.public_opened = os.fstat(self.public_fd)
        os.fsync(self.key_directory_fd)
        self._require_live_key_material()

    def ready_record(self) -> dict[str, Any]:
        self._require_live_key_material()
        return {
            "event": "SIGNER_READY",
            "openssl": artifact_from_fd(
                OPENSSL,
                self.openssl_fd,
                opened=self.openssl_opened,
            ),
            "public_key": artifact_from_fd(
                self.public_path,
                self.public_fd,
                opened=self.public_opened,
            ),
            "private_key_mode": oct(
                stat.S_IMODE(os.fstat(self.private_fd).st_mode)
            ),
            "key_directory_mode": oct(
                stat.S_IMODE(self._require_key_directory().st_mode)
            ),
            "expected_target": str(self.expected_target),
            "signing_protocol": (
                "canonical SIGN JSON with exact path,size_bytes,sha256,mode"
            ),
            "publication_semantics": (
                "exit 0 publishes ceremony output but grants no release authority; "
                "an independent consumer must reopen and verify every live path"
            ),
        }

    def _verify_signature(
        self,
        *,
        target_fd: int,
        signature_fd: int,
    ) -> None:
        self._run_openssl(
            [
                "pkeyutl",
                "-verify",
                "-rawin",
                "-pubin",
                "-inkey",
                f"/proc/self/fd/{self.public_fd}",
                "-sigfile",
                f"/proc/self/fd/{signature_fd}",
                "-in",
                f"/proc/self/fd/{target_fd}",
            ],
            pass_fds=(self.public_fd, signature_fd, target_fd),
        )

    def _retire_private_key(self) -> None:
        self._require_live_key_material()
        self.private_retirement_started = True
        os.fchmod(self.private_fd, 0o000)
        os.fsync(self.private_fd)

    def _point_in_time_validate_published_state(
        self,
        *,
        expected_target: Mapping[str, Any],
        expected_signature: Mapping[str, Any],
        target_fd: int,
        target_opened: os.stat_result,
        signature_fd: int,
        signature_directory: Path,
        signature_directory_fd: int,
        signature_directory_opened: os.stat_result,
    ) -> tuple[dict[str, Any], dict[str, Any], str]:
        """Capture a verified point-in-time snapshot before descriptor teardown."""
        require_path_bound(self.private_path, self.private_fd, expected_mode=0o000)
        private_mode = oct(stat.S_IMODE(os.fstat(self.private_fd).st_mode))
        self._require_key_directory()
        require_path_bound(
            self.public_path,
            self.public_fd,
            opened=self.public_opened,
            expected_mode=0o444,
        )
        require_directory_bound(
            signature_directory,
            signature_directory_fd,
            opened=signature_directory_opened,
        )
        require_path_bound(
            self.expected_target,
            target_fd,
            opened=target_opened,
            expected_mode=0o444,
        )
        require_path_bound(
            self.signature_path,
            signature_fd,
            expected_mode=0o444,
        )

        live_target_fd, live_target_opened = open_bound_regular(self.expected_target)
        live_signature_fd, live_signature_opened = open_bound_regular(
            self.signature_path
        )
        try:
            live_target = artifact_from_fd(
                self.expected_target,
                live_target_fd,
                opened=live_target_opened,
            )
            live_signature = artifact_from_fd(
                self.signature_path,
                live_signature_fd,
                opened=live_signature_opened,
            )
            if live_target != dict(expected_target):
                raise SignerError(
                    "Terminal target identity differs from the approved artifact"
                )
            if live_signature != dict(expected_signature):
                raise SignerError(
                    "Terminal signature identity differs from the published signature"
                )

            self._verify_signature(
                target_fd=live_target_fd,
                signature_fd=live_signature_fd,
            )

            # The verification callback may return after a concurrent mutation.
            # Reread both byte streams and path bindings before forming success.
            terminal_target = artifact_from_fd(
                self.expected_target,
                live_target_fd,
                opened=live_target_opened,
            )
            terminal_signature = artifact_from_fd(
                self.signature_path,
                live_signature_fd,
                opened=live_signature_opened,
            )
            if terminal_target != dict(expected_target):
                raise SignerError(
                    "Target identity changed during terminal signature verification"
                )
            if terminal_signature != dict(expected_signature):
                raise SignerError(
                    "Signature identity changed during terminal signature verification"
                )

            self._require_key_directory()
            require_path_bound(
                self.public_path,
                self.public_fd,
                opened=self.public_opened,
                expected_mode=0o444,
            )
            require_path_bound(
                self.private_path,
                self.private_fd,
                expected_mode=0o000,
            )
            require_directory_bound(
                signature_directory,
                signature_directory_fd,
                opened=signature_directory_opened,
            )

            held_target = artifact_from_fd(
                self.expected_target,
                target_fd,
                opened=target_opened,
            )
            held_signature = artifact_from_fd(
                self.signature_path,
                signature_fd,
            )
            if held_target != terminal_target:
                raise SignerError(
                    "Held and reopened target identities differ at terminal validation"
                )
            if held_signature != terminal_signature:
                raise SignerError(
                    "Held and reopened signature identities differ at terminal validation"
                )
            return terminal_target, terminal_signature, private_mode
        finally:
            os.close(live_signature_fd)
            os.close(live_target_fd)

    def failure_record(self, exc: BaseException) -> dict[str, Any]:
        try:
            private_mode = stat.S_IMODE(os.fstat(self.private_fd).st_mode)
            if private_mode == 0o000:
                private_state = "RETIRED_MODE_000"
            elif private_mode == 0o400:
                private_state = "RETAINED_MODE_400"
            else:
                private_state = f"UNEXPECTED_MODE_{oct(private_mode)}"
        except OSError:
            private_mode = None
            private_state = "UNAVAILABLE"
        return {
            "event": "SIGN_FAILED",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "signature_published_by_signer": self.signature_published,
            "signature_path_exists": os.path.lexists(self.signature_path),
            "private_retirement_started": self.private_retirement_started,
            "private_key_mode": None if private_mode is None else oct(private_mode),
            "private_key_state": private_state,
        }

    def sign(self, expected_target_record: Mapping[str, Any]) -> dict[str, Any]:
        expected_target = validate_expected_target_record(
            expected_target_record,
            expected_path=self.expected_target,
        )
        if (
            not self.expected_target.is_file()
            or os.path.lexists(self.signature_path)
            or os.path.lexists(self.pending_path)
        ):
            raise SignerError("Signing precondition failed")

        self._require_live_key_material()
        target_fd, target_opened = open_bound_regular(self.expected_target)
        signature_directory = self.signature_path.parent
        signature_directory_fd, signature_directory_opened = open_bound_directory(
            signature_directory
        )
        signature_fd = -1
        try:
            target_record = artifact_from_fd(
                self.expected_target,
                target_fd,
                opened=target_opened,
            )
            if target_record != expected_target:
                raise SignerError(
                    "Opened target identity does not match the approved signing instruction"
                )
            self._require_live_key_material()
            signature = self._run_openssl(
                [
                    "pkeyutl",
                    "-sign",
                    "-rawin",
                    "-inkey",
                    f"/proc/self/fd/{self.private_fd}",
                    "-passin",
                    "env:PASS",
                    "-in",
                    f"/proc/self/fd/{target_fd}",
                ],
                pass_fds=(self.private_fd, target_fd),
            )
            if len(signature) != 64:
                raise SignerError(
                    f"Unexpected Ed25519 signature length: {len(signature)}"
                )

            signature_fd = create_anonymous_file(signature_directory_fd, 0o400)
            write_all(signature_fd, signature)
            self._verify_signature(target_fd=target_fd, signature_fd=signature_fd)
            require_path_bound(
                self.expected_target,
                target_fd,
                opened=target_opened,
                expected_mode=0o444,
            )
            self._require_live_key_material()
            require_directory_bound(
                signature_directory,
                signature_directory_fd,
                opened=signature_directory_opened,
            )

            os.fchmod(signature_fd, 0o444)
            os.fsync(signature_fd)
            link_fd_no_replace(
                signature_fd,
                signature_directory_fd,
                self.signature_path.name,
            )
            self.signature_published = True
            os.fsync(signature_directory_fd)
            require_directory_bound(
                signature_directory,
                signature_directory_fd,
                opened=signature_directory_opened,
            )
            require_path_bound(
                self.signature_path,
                signature_fd,
                expected_mode=0o444,
            )
            self._verify_signature(target_fd=target_fd, signature_fd=signature_fd)
            require_path_bound(
                self.signature_path,
                signature_fd,
                expected_mode=0o444,
            )
            require_directory_bound(
                signature_directory,
                signature_directory_fd,
                opened=signature_directory_opened,
            )
            require_path_bound(
                self.expected_target,
                target_fd,
                opened=target_opened,
                expected_mode=0o444,
            )
            self._require_live_key_material()

            expected_signature_record = artifact_from_fd(
                self.signature_path,
                signature_fd,
            )
            target_record = artifact_from_fd(
                self.expected_target,
                target_fd,
                opened=target_opened,
            )
            if target_record != expected_target:
                raise SignerError(
                    "Published signature target identity drifted from approval"
                )
            self._retire_private_key()
            require_path_bound(self.private_path, self.private_fd, expected_mode=0o000)
            self._require_key_directory()
            require_path_bound(
                self.public_path,
                self.public_fd,
                opened=self.public_opened,
                expected_mode=0o444,
            )
            target_record = artifact_from_fd(
                self.expected_target,
                target_fd,
                opened=target_opened,
            )
            if target_record != expected_target:
                raise SignerError(
                    "Target identity changed after private-key retirement"
                )
            os.fsync(self.key_directory_fd)
            self._passphrase = ""
            self._environment.pop("PASS", None)
            target_record, signature_record, private_mode = (
                self._point_in_time_validate_published_state(
                    expected_target=expected_target,
                    expected_signature=expected_signature_record,
                    target_fd=target_fd,
                    target_opened=target_opened,
                    signature_fd=signature_fd,
                    signature_directory=signature_directory,
                    signature_directory_fd=signature_directory_fd,
                    signature_directory_opened=signature_directory_opened,
                )
            )
            return {
                "event": (
                    "CEREMONY_OUTPUT_AVAILABLE_REQUIRES_INDEPENDENT_VERIFICATION"
                ),
                "target": target_record,
                "signature": signature_record,
                "private_key_mode": private_mode,
                "point_in_time_checks": "PASS_BEFORE_DESCRIPTOR_TEARDOWN",
                "consumer_verification_required": True,
                "release_authority_granted": False,
            }
        finally:
            if signature_fd >= 0:
                os.close(signature_fd)
            os.close(signature_directory_fd)
            os.close(target_fd)

    def close(self) -> None:
        if self.public_fd >= 0:
            os.close(self.public_fd)
            self.public_fd = -1
        if self.private_fd >= 0:
            os.close(self.private_fd)
            self.private_fd = -1
        if self.openssl_fd >= 0:
            os.close(self.openssl_fd)
            self.openssl_fd = -1
        if self.key_directory_fd >= 0:
            os.close(self.key_directory_fd)
            self.key_directory_fd = -1
        self._passphrase = ""
        self._environment.pop("PASS", None)


def require_bound_bootstrap() -> dict[str, Any]:
    argv0 = sys.argv[0]
    if (
        not argv0.startswith("/proc/self/fd/")
        or not argv0.removeprefix("/proc/self/fd/").isdigit()
        or sys.flags.isolated != 1
        or sys.flags.no_user_site != 1
        or dict(os.environ) != BOOTSTRAP_ENVIRONMENT
    ):
        raise SignerError("Signer requires clean python -I /proc/self/fd bootstrap")
    script_path = Path(__file__).resolve(strict=True)
    script_fd = int(argv0.removeprefix("/proc/self/fd/"))
    script_record = artifact_from_fd(script_path, script_fd)
    if script_record["mode"] != "0o444":
        raise SignerError("Signer source must be mode 0444 before key generation")
    python_path = PYTHON.resolve(strict=True)
    python_fd, python_opened = open_bound_regular(python_path)
    try:
        python_record = artifact_from_fd(
            python_path,
            python_fd,
            opened=python_opened,
        )
    finally:
        os.close(python_fd)
    if (
        python_record["sha256"] != PYTHON_SHA256
        or python_record["mode"] != "0o775"
    ):
        raise SignerError("Signer Python runtime identity drifted")
    return {
        "schema_name": "intersubmod.bound_python_entrypoint",
        "schema_version": "1.0.0",
        "argv0": argv0,
        "script": script_record,
        "python": python_record,
        "isolated_mode": True,
        "no_user_site": True,
        "environment": BOOTSTRAP_ENVIRONMENT,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("key_dir", type=Path)
    parser.add_argument("private_filename")
    parser.add_argument("expected_target", type=Path)
    return parser.parse_args()


def main() -> int:
    signer: OneTimeSigner | None = None
    try:
        bootstrap = require_bound_bootstrap()
        args = parse_args()
        signer = OneTimeSigner(
            key_dir=args.key_dir,
            private_filename=args.private_filename,
            expected_target=args.expected_target,
        )
        ready = signer.ready_record()
        ready["bootstrap"] = bootstrap
        print(json.dumps(ready, sort_keys=True), flush=True)
        for raw_line in sys.stdin:
            line = raw_line.rstrip("\r\n")
            try:
                expected_target_record = parse_signing_instruction(
                    line,
                    expected_path=signer.expected_target,
                )
            except SignerError:
                print("WAITING_FOR_EXACT_SIGN_TARGET", flush=True)
                continue
            try:
                result = signer.sign(expected_target_record)
            except Exception as exc:
                print(
                    json.dumps(signer.failure_record(exc), sort_keys=True),
                    file=sys.stderr,
                    flush=True,
                )
                return 8
            print(json.dumps(result, sort_keys=True), flush=True)
            return 0
        print("SIGNER_INPUT_CLOSED_WITHOUT_SIGNATURE", file=sys.stderr, flush=True)
        return 7
    except Exception as exc:
        print(f"SIGNER_SETUP_FAILED error={exc}", file=sys.stderr, flush=True)
        return 3
    finally:
        if signer is not None:
            signer.close()


if __name__ == "__main__":
    raise SystemExit(main())
