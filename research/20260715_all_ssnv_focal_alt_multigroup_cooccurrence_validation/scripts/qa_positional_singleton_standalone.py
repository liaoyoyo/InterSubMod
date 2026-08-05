#!/usr/bin/env python3
"""Run fail-closed browser QA for the positional-singleton standalone report."""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import importlib.metadata
import json
import os
import re
import stat
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

import playwright
import PIL
import PIL._imaging
from PIL import Image, ImageChops
from playwright.sync_api import BrowserContext, Page, Route, WebSocketRoute, sync_playwright


SCHEMA_NAME = "intersubmod.positional_singleton_standalone_qa"
SCHEMA_VERSION = "2.4.0"
SCREEN_PROFILES = (
    ("desktop_1440x1000", 1440, 1000),
    ("desktop_1024x768", 1024, 768),
    ("mobile_390x844", 390, 844),
    ("mobile_320x568", 320, 568),
)
PRINT_PROFILE = ("print_a4_794x1123", 794, 1123)
EXPECTED_IMAGE_COUNT = 5
EXPECTED_TABLE_COUNT = 8
EXPECTED_ANCHOR_COUNT = 9
IMAGE_CONTENT_RATIO_RELATIVE_ERROR_MAX = 0.005
IMAGE_DENSE_NON_DOMINANT_PIXEL_FRACTION_MIN = 0.05
IMAGE_SPARSE_NON_DOMINANT_PIXEL_FRACTION_MIN = 0.005
IMAGE_SPARSE_NON_DOMINANT_BBOX_FRACTION_MIN = 0.20
IMAGE_SPARSE_NON_DOMINANT_GRID_FRACTION_MIN = 0.10
IMAGE_SPARSE_MINIMUM_UNIQUE_OPAQUE_COLORS = 4
IMAGE_SPATIAL_GRID_SIZE = 8
TOPIC_ROOT = Path(__file__).resolve().parents[1]
FORMAL_RELEASE_RECEIPT = (
    TOPIC_ROOT / "results" / "positional_singleton_supplemental_release_receipt.v1.json"
)
FORMAL_RELEASE_VERIFICATION = (
    TOPIC_ROOT
    / "results"
    / "positional_singleton_supplemental_release_post_signature_verification.v1.json"
)
FORMAL_RELEASE_PUBLIC_KEY_SHA256 = (
    "c8490467e7ad9c0477429d37c799dd4143773b781b6fc33ab08199a67dd99b91"
)
FORMAL_RELEASE_FINALIZER = TOPIC_ROOT / "scripts" / "finalize_positional_singleton_supplemental_release.py"
PDFINFO = Path("/usr/bin/pdfinfo")
PDFTOPPM = Path("/usr/bin/pdftoppm")
PDFTOTEXT = Path("/usr/bin/pdftotext")
OPENSSL = Path("/usr/bin/openssl")
FC_MATCH = Path("/usr/bin/fc-match")
LDD = Path("/usr/bin/ldd")


class QaError(RuntimeError):
    """Raised when the QA contract or publication contract fails."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--executable-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--evidence-kind",
        choices=("formal_report", "self_test"),
        required=True,
    )
    parser.add_argument("--release-receipt", type=Path)
    parser.add_argument("--release-signature", type=Path)
    parser.add_argument("--release-public-key", type=Path)
    parser.add_argument("--release-private-key", type=Path)
    parser.add_argument("--release-verification", type=Path)
    parser.add_argument("--wait-ms", type=int, default=1_500)
    return parser.parse_args()


def path_exists(path: Path) -> bool:
    return os.path.lexists(path)


def artifact(path: Path, *, public_path: Path | None = None) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(resolved, flags)
    try:
        status_before = os.fstat(descriptor)
        if not stat.S_ISREG(status_before.st_mode):
            raise QaError(f"Expected regular file: {resolved}")
        digest = hashlib.sha256()
        size_bytes = 0
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            size_bytes += len(chunk)
        status_after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    stable_fields = (
        "st_dev",
        "st_ino",
        "st_mode",
        "st_size",
        "st_mtime_ns",
        "st_ctime_ns",
    )
    if any(
        getattr(status_before, field) != getattr(status_after, field)
        for field in stable_fields
    ) or size_bytes != status_before.st_size:
        raise QaError(f"File changed while its identity was read: {resolved}")
    return {
        "path": str((public_path or resolved).resolve(strict=False)),
        "size_bytes": size_bytes,
        "sha256": digest.hexdigest(),
        "mode": oct(stat.S_IMODE(status_before.st_mode)),
    }


def read_file_snapshot(path: Path) -> tuple[bytes, dict[str, Any]]:
    resolved = path.resolve(strict=True)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(resolved, flags)
    try:
        status_before = os.fstat(descriptor)
        if not stat.S_ISREG(status_before.st_mode):
            raise QaError(f"Expected regular file: {resolved}")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        status_after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    data = b"".join(chunks)
    stable_fields = (
        "st_dev",
        "st_ino",
        "st_mode",
        "st_size",
        "st_mtime_ns",
        "st_ctime_ns",
    )
    if any(
        getattr(status_before, field) != getattr(status_after, field)
        for field in stable_fields
    ) or len(data) != status_before.st_size:
        raise QaError(f"File changed while its bytes were snapshotted: {resolved}")
    return data, {
        "path": str(resolved),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(stat.S_IMODE(status_before.st_mode)),
    }


def directory_manifest(root: Path) -> dict[str, Any]:
    resolved = root.resolve(strict=True)
    if not resolved.is_dir():
        raise QaError(f"Expected runtime bundle directory: {resolved}")
    entries: list[dict[str, Any]] = []
    total_size = 0
    for candidate in sorted(resolved.rglob("*"), key=lambda path: str(path.relative_to(resolved))):
        relative = str(candidate.relative_to(resolved))
        status = candidate.lstat()
        mode = oct(stat.S_IMODE(status.st_mode))
        if stat.S_ISREG(status.st_mode):
            identity = artifact(candidate)
            total_size += identity["size_bytes"]
            entries.append(
                {
                    "relative_path": relative,
                    "type": "regular",
                    "size_bytes": identity["size_bytes"],
                    "sha256": identity["sha256"],
                    "mode": identity["mode"],
                }
            )
        elif stat.S_ISLNK(status.st_mode):
            entries.append(
                {
                    "relative_path": relative,
                    "type": "symlink",
                    "target": os.readlink(candidate),
                    "mode": mode,
                }
            )
    canonical = json.dumps(
        entries,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return {
        "path": str(resolved),
        "root_mode": oct(stat.S_IMODE(resolved.stat().st_mode)),
        "entry_count": len(entries),
        "regular_file_count": sum(entry["type"] == "regular" for entry in entries),
        "symlink_count": sum(entry["type"] == "symlink" for entry in entries),
        "total_regular_file_size_bytes": total_size,
        "manifest_sha256": hashlib.sha256(canonical).hexdigest(),
        "entries": entries,
    }


class StaticSecurityParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=False)
        self.violations: list[dict[str, Any]] = []

    def inspect(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        normalized_tag = tag.lower()
        if normalized_tag == "script":
            self.violations.append({"tag": normalized_tag, "reason": "script_tag"})
        for name, value in attrs:
            normalized_name = name.lower()
            normalized_value = (value or "").strip()
            if normalized_name.startswith("on"):
                self.violations.append(
                    {
                        "tag": normalized_tag,
                        "attribute": normalized_name,
                        "reason": "inline_event_handler",
                    }
                )
            if normalized_value.lower().startswith(("javascript:", "vbscript:")):
                self.violations.append(
                    {
                        "tag": normalized_tag,
                        "attribute": normalized_name,
                        "reason": "active_url_scheme",
                    }
                )

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.inspect(tag, attrs)

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.inspect(tag, attrs)


def static_html_security_scan(data: bytes) -> dict[str, Any]:
    try:
        source = data.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise QaError("Standalone HTML is not strict UTF-8") from error
    parser = StaticSecurityParser()
    parser.feed(source)
    parser.close()
    return {
        "parser": "python.html.parser.HTMLParser",
        "source_size_bytes": len(data),
        "source_sha256": hashlib.sha256(data).hexdigest(),
        "violations": parser.violations,
        "pass": not parser.violations,
    }


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise QaError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_nonfinite(value: str) -> None:
    raise QaError(f"Non-finite JSON value: {value}")


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=reject_duplicate_keys,
            parse_constant=reject_nonfinite,
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, QaError) as error:
        raise QaError(f"Unable to read {label}: {path}") from error
    if not isinstance(value, dict):
        raise QaError(f"{label} must be a JSON object")
    return value


def require_identity(observed: dict[str, Any], expected: Any, *, label: str) -> None:
    if not isinstance(expected, dict) or set(expected) != {
        "path",
        "size_bytes",
        "sha256",
        "mode",
    }:
        raise QaError(f"{label} identity schema is not exact")
    if observed != expected:
        raise QaError(f"{label} identity drift: {observed!r} != {expected!r}")


def verify_ed25519(
    *, receipt: Path, signature: Path, public_key: Path, label: str
) -> None:
    completed = subprocess.run(
        [
            str(OPENSSL),
            "pkeyutl",
            "-verify",
            "-rawin",
            "-pubin",
            "-inkey",
            str(public_key),
            "-sigfile",
            str(signature),
            "-in",
            str(receipt),
        ],
        env={"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C"},
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise QaError(f"{label} detached Ed25519 signature failed")


def resolve_formal_release_paths(args: argparse.Namespace) -> dict[str, Path] | None:
    values = {
        "receipt": args.release_receipt,
        "signature": args.release_signature,
        "public_key": args.release_public_key,
        "private_key": args.release_private_key,
        "verification": args.release_verification,
    }
    if args.evidence_kind == "self_test":
        if any(value is not None for value in values.values()):
            raise QaError("Self-test mode must not accept formal release arguments")
        return None
    missing = sorted(name for name, value in values.items() if value is None)
    if missing:
        raise QaError(
            "Formal report mode requires signed release arguments: " + ", ".join(missing)
        )
    resolved = {
        name: value.expanduser().resolve(strict=True)
        for name, value in values.items()
        if value is not None
    }
    if resolved["receipt"] != FORMAL_RELEASE_RECEIPT.resolve(strict=True):
        raise QaError("Formal QA release receipt is not the canonical Task-B path")
    if resolved["signature"] != Path(str(resolved["receipt"]) + ".ed25519.sig"):
        raise QaError("Formal QA signature path differs from the receipt contract")
    if resolved["verification"] != FORMAL_RELEASE_VERIFICATION.resolve(strict=True):
        raise QaError("Formal QA verification is not the canonical Task-B path")
    return resolved


def validate_formal_release_binding(
    paths: dict[str, Path] | None,
    *,
    html_identity: dict[str, Any],
) -> dict[str, Any] | None:
    if paths is None:
        return None
    identities = {
        name: artifact(path)
        for name, path in paths.items()
        if name != "private_key"
    }
    private_status = paths["private_key"].stat()
    if not stat.S_ISREG(private_status.st_mode):
        raise QaError("Formal release private-key path is not a regular file")
    identities["private_key"] = {
        "path": str(paths["private_key"]),
        "mode": oct(stat.S_IMODE(private_status.st_mode)),
    }
    if any(
        identities[name]["mode"] != "0o444"
        for name in ("receipt", "signature", "public_key", "verification")
    ):
        raise QaError("Formal release evidence is not protected mode 0444")
    if identities["private_key"]["mode"] != "0o0":
        raise QaError("Formal release one-time private key is not retired")
    if identities["public_key"]["sha256"] != FORMAL_RELEASE_PUBLIC_KEY_SHA256:
        raise QaError("Formal release public-key SHA-256 drift")
    verify_ed25519(
        receipt=paths["receipt"],
        signature=paths["signature"],
        public_key=paths["public_key"],
        label="formal supplemental release",
    )

    receipt = load_json(paths["receipt"], label="formal supplemental release receipt")
    if (
        receipt.get("schema_name")
        != "intersubmod.positional_singleton_supplemental_release_receipt"
        or receipt.get("schema_version") != "1.0.0"
        or receipt.get("task_type") != "B_comprehensive_validation"
        or receipt.get("scope") != "all_50432_positional_singleton_dataset_sites"
        or receipt.get("claim_ceiling")
        != "M2_read_level_residual_epigenetic_partition"
        or receipt.get("pass") is not True
    ):
        raise QaError("Formal supplemental release receipt contract drift")
    signature_contract = receipt.get("signature_contract", {})
    require_identity(
        identities["public_key"],
        signature_contract.get("public_key"),
        label="formal release public key",
    )
    if signature_contract.get("expected_signature_path") != str(paths["signature"]):
        raise QaError("Formal release signature contract path drift")

    finalizer_identity = artifact(FORMAL_RELEASE_FINALIZER)
    require_identity(
        finalizer_identity,
        receipt.get("code", {}).get("supplemental_release_finalizer"),
        label="formal release finalizer",
    )
    report_receipt_record = receipt.get("inputs", {}).get("report_receipt")
    if not isinstance(report_receipt_record, dict):
        raise QaError("Formal release lacks supplemental report receipt identity")
    report_receipt_path = Path(str(report_receipt_record.get("path", ""))).resolve(
        strict=True
    )
    report_receipt_identity = artifact(report_receipt_path)
    require_identity(
        report_receipt_identity,
        report_receipt_record,
        label="supplemental report receipt",
    )
    report_receipt = load_json(report_receipt_path, label="supplemental report receipt")
    if (
        report_receipt.get("schema_name")
        != "intersubmod.positional_singleton_methyl_multigroup_report"
        or report_receipt.get("schema_version") != "1.0.0"
        or report_receipt.get("task_type") != "B_comprehensive_validation"
        or report_receipt.get("scope")
        != "all_50432_positional_singleton_dataset_sites"
        or report_receipt.get("claim_ceiling")
        != "M2_read_level_residual_epigenetic_partition"
        or report_receipt.get("pass") is not True
    ):
        raise QaError("Supplemental report receipt contract drift")
    require_identity(
        html_identity,
        report_receipt.get("outputs", {}).get("portable_html"),
        label="formal portable HTML",
    )
    report_builder_path = Path(
        str(report_receipt.get("code", {}).get("report_builder", {}).get("path", ""))
    ).resolve(strict=True)
    require_identity(
        artifact(report_builder_path),
        report_receipt.get("code", {}).get("report_builder"),
        label="supplemental report builder",
    )

    report_success_record = receipt.get("inputs", {}).get("report_success")
    if not isinstance(report_success_record, dict):
        raise QaError("Formal release lacks supplemental report success identity")
    report_success_path = Path(str(report_success_record.get("path", ""))).resolve(
        strict=True
    )
    require_identity(
        artifact(report_success_path),
        report_success_record,
        label="supplemental report success marker",
    )
    report_success = load_json(report_success_path, label="supplemental report success marker")
    if (
        report_success.get("schema_name") != "intersubmod.atomic_release_marker"
        or report_success.get("pass") is not True
    ):
        raise QaError("Supplemental report success-marker contract drift")
    require_identity(
        report_receipt_identity,
        report_success.get("receipt"),
        label="supplemental report success binding",
    )

    verification = load_json(paths["verification"], label="post-sign verification")
    if (
        verification.get("schema_name")
        != "intersubmod.positional_singleton_supplemental_post_signature_verification"
        or verification.get("schema_version") != "1.0.0"
        or verification.get("task_type") != "B_comprehensive_validation"
        or verification.get("scope")
        != "all_50432_positional_singleton_dataset_sites"
        or verification.get("verification_status")
        != "SIGNED_SUPPLEMENTAL_RELEASE_REVERIFIED"
        or verification.get("pass") is not True
        or not all(
            value is True for value in verification.get("checks", {}).values()
        )
    ):
        raise QaError("Formal post-sign verification contract drift")
    for name, field in (
        ("receipt", "verified_receipt"),
        ("signature", "detached_signature"),
        ("public_key", "public_key"),
    ):
        require_identity(
            identities[name], verification.get(field), label=f"post-sign {field}"
        )
    if identities["private_key"] != verification.get("retired_private_key"):
        raise QaError("Post-sign retired private-key metadata drift")
    require_identity(
        finalizer_identity,
        verification.get("code", {}).get("supplemental_release_finalizer"),
        label="post-sign finalizer",
    )
    return {
        "release_receipt": identities["receipt"],
        "release_signature": identities["signature"],
        "release_public_key": identities["public_key"],
        "retired_private_key": identities["private_key"],
        "post_signature_verification": identities["verification"],
        "report_receipt": report_receipt_identity,
        "report_success": artifact(report_success_path),
        "report_builder": artifact(report_builder_path),
        "release_finalizer": finalizer_identity,
        "portable_html_identity_bound": True,
        "detached_signature_verified": True,
        "pass": True,
    }


def repository_head(repo_root: Path) -> str:
    result = subprocess.run(
        ["/usr/bin/git", "rev-parse", "HEAD"],
        cwd=repo_root,
        env={"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C"},
        check=True,
        capture_output=True,
        text=True,
    )
    head = result.stdout.strip()
    if len(head) != 40 or any(character not in "0123456789abcdef" for character in head):
        raise QaError("Unable to resolve a canonical Git HEAD")
    return head


def command_snapshot(command: list[str], *, label: str) -> dict[str, Any]:
    completed = subprocess.run(
        command,
        env={"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C"},
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise QaError(f"{label} command failed with exit {completed.returncode}")
    return {
        "command": command,
        "returncode": completed.returncode,
        "stdout_size_bytes": len(completed.stdout),
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stdout": completed.stdout.decode("utf-8", errors="replace"),
        "stderr_size_bytes": len(completed.stderr),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "stderr": completed.stderr.decode("utf-8", errors="replace"),
    }


def ldd_snapshot(path: Path) -> dict[str, Any]:
    snapshot = command_snapshot([str(LDD), str(path)], label=f"ldd {path}")
    libraries: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for match in re.finditer(r"(?:=>\s+)?(/[^\s(]+)", snapshot["stdout"]):
        library_path = Path(match.group(1))
        if library_path in seen or not library_path.is_file():
            continue
        seen.add(library_path)
        libraries.append(artifact(library_path))
    normalized_stdout = re.sub(r"\(0x[0-9a-fA-F]+\)", "(ASLR_ADDRESS)", snapshot["stdout"])
    normalized_bytes = normalized_stdout.encode("utf-8")
    snapshot["stdout"] = normalized_stdout
    snapshot["stdout_size_bytes"] = len(normalized_bytes)
    snapshot["stdout_sha256"] = hashlib.sha256(normalized_bytes).hexdigest()
    snapshot["normalization"] = "replace_parenthesized_ASLR_hex_addresses"
    snapshot["resolved_libraries"] = libraries
    return snapshot


def runtime_identity(executable_path: Path) -> dict[str, Any]:
    required_tools = {
        "pdfinfo": PDFINFO,
        "pdftoppm": PDFTOPPM,
        "pdftotext": PDFTOTEXT,
        "openssl": OPENSSL,
        "fc_match": FC_MATCH,
        "ldd": LDD,
    }
    for label, tool_path in required_tools.items():
        if not tool_path.is_file():
            raise QaError(f"Missing required QA runtime tool {label}: {tool_path}")
    font_snapshot = command_snapshot(
        [str(FC_MATCH), "-f", "%{file}\n", "sans-serif"],
        label="fontconfig sans-serif resolution",
    )
    font_lines = [line for line in font_snapshot["stdout"].splitlines() if line]
    if len(font_lines) != 1:
        raise QaError("Fontconfig did not resolve exactly one sans-serif font")
    font_path = Path(font_lines[0]).resolve(strict=True)
    playwright_root = Path(playwright.__file__).resolve().parent
    pillow_root = Path(PIL.__file__).resolve().parent
    playwright_driver = playwright_root / "driver" / "node"
    chromium_root = executable_path.resolve(strict=True).parent
    return {
        "python": {
            "version": sys.version,
            "executable": artifact(Path(sys.executable)),
            "ldd": ldd_snapshot(Path(sys.executable)),
        },
        "playwright": {
            "version": importlib.metadata.version("playwright"),
            "module": artifact(Path(playwright.__file__)),
            "package_bundle": directory_manifest(playwright_root),
            "driver_node": artifact(playwright_driver),
            "driver_node_ldd": ldd_snapshot(playwright_driver),
        },
        "pillow": {
            "version": PIL.__version__,
            "module": artifact(Path(PIL.__file__)),
            "imaging_extension": artifact(Path(PIL._imaging.__file__)),
            "package_bundle": directory_manifest(pillow_root),
        },
        "chromium": {
            "executable": artifact(executable_path),
            "ldd": ldd_snapshot(executable_path),
            "runtime_bundle": directory_manifest(chromium_root),
        },
        "tools": {label: artifact(path) for label, path in required_tools.items()},
        "font_resolution": {
            "snapshot": font_snapshot,
            "resolved_font": artifact(font_path),
        },
    }


def fsync_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def rename_no_replace(source: Path, target: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    if renameat2 is None:
        raise QaError("renameat2 is required for no-replace publication")
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        -100,
        os.fsencode(source),
        -100,
        os.fsencode(target),
        1,
    )
    if result != 0:
        error_number = ctypes.get_errno()
        raise FileExistsError(
            error_number,
            f"No-replace publication failed: {os.strerror(error_number)}",
            str(target),
        )


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    path.chmod(0o444)
    fsync_file(path)


def write_bytes_exclusive(path: Path, payload: bytes) -> None:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0),
        0o600,
    )
    try:
        offset = 0
        while offset < len(payload):
            offset += os.write(descriptor, payload[offset:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    path.chmod(0o444)
    fsync_file(path)


def collect_layout(
    page: Page,
    *,
    mobile: bool,
    print_media: bool,
    allow_svg_fixture: bool,
) -> dict[str, Any]:
    result = page.evaluate(
        """
        ({ expectedImages, expectedTables, expectedAnchors, mobile, printMedia,
           allowSvgFixture, imageThresholds }) => {
          const round = (value) => Math.round(value * 1000) / 1000;
          const visible = (element) => {
            const rect = element.getBoundingClientRect();
            const style = getComputedStyle(element);
            return rect.width > 0 && rect.height > 0 &&
              style.display !== 'none' && style.visibility !== 'hidden';
          };
          const label = (element) => {
            if (element.id) return `${element.tagName.toLowerCase()}#${element.id}`;
            const classes = Array.from(element.classList || []).slice(0, 3);
            return `${element.tagName.toLowerCase()}${classes.length ? '.' + classes.join('.') : ''}`;
          };
          const allowedOverflow = (element) => Boolean(
            element.closest('.table-wrap') || element.closest('nav .inner')
          );
          const documentElement = document.documentElement;
          const viewportWidth = documentElement.clientWidth;
          const outOfBounds = Array.from(document.querySelectorAll('body *'))
            .filter(visible)
            .filter((element) => {
              const rect = element.getBoundingClientRect();
              return !allowedOverflow(element) &&
                (rect.left < -1 || rect.right > viewportWidth + 1);
            })
            .map((element) => {
              const rect = element.getBoundingClientRect();
              return {
                element: label(element),
                left: round(rect.left),
                right: round(rect.right),
                width: round(rect.width),
              };
            });
          const unexpectedPositionedElements = Array.from(
            document.querySelectorAll('body *')
          ).filter(visible).filter((element) => {
            const position = getComputedStyle(element).position;
            return ['absolute', 'fixed'].includes(position) &&
              !element.matches('.skip-link');
          }).map((element) => ({
            element: label(element),
            position: getComputedStyle(element).position,
            zIndex: getComputedStyle(element).zIndex,
          }));
          const unexpectedPositionedPseudoElements = [];
          for (const element of [document.body, ...document.querySelectorAll('body *')]) {
            for (const pseudo of ['::before', '::after']) {
              const style = getComputedStyle(element, pseudo);
              const declared = !['none', 'normal'].includes(style.content) &&
                style.display !== 'none';
              if (declared) {
                unexpectedPositionedPseudoElements.push({
                  element: label(element),
                  pseudo,
                  content: style.content,
                  position: style.position,
                  zIndex: style.zIndex,
                  pointerEvents: style.pointerEvents,
                  backgroundColor: style.backgroundColor,
                  backgroundImage: style.backgroundImage,
                });
              }
            }
          }

          const overlapPairs = [];
          const containers = [
            document.body,
            document.querySelector('main'),
            ...document.querySelectorAll('section > .inner, .kpis, figure, details'),
          ].filter(Boolean);
          for (const container of containers) {
            const children = Array.from(container.children).filter(visible);
            for (let leftIndex = 0; leftIndex < children.length; leftIndex += 1) {
              const left = children[leftIndex];
              const leftRect = left.getBoundingClientRect();
              for (let rightIndex = leftIndex + 1; rightIndex < children.length; rightIndex += 1) {
                const right = children[rightIndex];
                const rightRect = right.getBoundingClientRect();
                const overlapWidth = Math.min(leftRect.right, rightRect.right) -
                  Math.max(leftRect.left, rightRect.left);
                const overlapHeight = Math.min(leftRect.bottom, rightRect.bottom) -
                  Math.max(leftRect.top, rightRect.top);
                const intentionalStickyOverlap =
                  [left, right].some((element) =>
                    getComputedStyle(element).position === 'sticky' &&
                    element.matches('nav, .ribbon, .scope-ribbon')
                  ) && [left, right].some((element) => element.matches('main'));
                if (overlapWidth > 1 && overlapHeight > 1 &&
                    !intentionalStickyOverlap) {
                  overlapPairs.push({
                    container: label(container),
                    left: label(left),
                    right: label(right),
                    overlapWidth: round(overlapWidth),
                    overlapHeight: round(overlapHeight),
                  });
                }
              }
            }
          }

          const images = Array.from(document.images).map((image, index) => {
            const canvas = document.createElement('canvas');
            const rasterScale = Math.min(
              1,
              512 / image.naturalWidth,
              512 / image.naturalHeight
            );
            canvas.width = Math.max(1, Math.round(image.naturalWidth * rasterScale));
            canvas.height = Math.max(1, Math.round(image.naturalHeight * rasterScale));
            const context = canvas.getContext('2d', { willReadFrequently: true });
            context.imageSmoothingEnabled = false;
            let pixelCheckError = null;
            let nonTransparentPixels = 0;
            let uniqueOpaqueColors = 0;
            let alphaCoverage = 0;
            let dominantOpaqueFraction = 1;
            let nonDominantPixelFraction = 0;
            let nonDominantBboxFraction = 0;
            let nonDominantOccupiedGridCells = 0;
            let nonDominantGridFraction = 0;
            try {
              context.drawImage(image, 0, 0, canvas.width, canvas.height);
              const pixels = context.getImageData(0, 0, canvas.width, canvas.height).data;
              const colors = new Map();
              for (let offset = 0; offset < pixels.length; offset += 4) {
                if (pixels[offset + 3] > 16) {
                  nonTransparentPixels += 1;
                  const color = `${pixels[offset]},${pixels[offset + 1]},${pixels[offset + 2]}`;
                  colors.set(color, (colors.get(color) || 0) + 1);
                }
              }
              uniqueOpaqueColors = colors.size;
              alphaCoverage = nonTransparentPixels / (canvas.width * canvas.height);
              let dominantColor = null;
              let dominantCount = 0;
              for (const [color, count] of colors.entries()) {
                if (count > dominantCount) {
                  dominantColor = color;
                  dominantCount = count;
                }
              }
              dominantOpaqueFraction = nonTransparentPixels > 0
                ? dominantCount / nonTransparentPixels
                : 1;
              nonDominantPixelFraction = nonTransparentPixels > 0
                ? (nonTransparentPixels - dominantCount) / nonTransparentPixels
                : 0;
              let minX = canvas.width;
              let minY = canvas.height;
              let maxX = -1;
              let maxY = -1;
              const occupiedGridCells = new Set();
              for (let offset = 0; offset < pixels.length; offset += 4) {
                if (pixels[offset + 3] <= 16) continue;
                const color = `${pixels[offset]},${pixels[offset + 1]},${pixels[offset + 2]}`;
                if (color === dominantColor) continue;
                const pixelIndex = offset / 4;
                const x = pixelIndex % canvas.width;
                const y = Math.floor(pixelIndex / canvas.width);
                minX = Math.min(minX, x);
                minY = Math.min(minY, y);
                maxX = Math.max(maxX, x);
                maxY = Math.max(maxY, y);
                const gridX = Math.min(
                  imageThresholds.spatialGridSize - 1,
                  Math.floor(x * imageThresholds.spatialGridSize / canvas.width)
                );
                const gridY = Math.min(
                  imageThresholds.spatialGridSize - 1,
                  Math.floor(y * imageThresholds.spatialGridSize / canvas.height)
                );
                occupiedGridCells.add(`${gridX},${gridY}`);
              }
              nonDominantOccupiedGridCells = occupiedGridCells.size;
              nonDominantGridFraction = occupiedGridCells.size /
                (imageThresholds.spatialGridSize * imageThresholds.spatialGridSize);
              if (maxX >= minX && maxY >= minY) {
                nonDominantBboxFraction =
                  ((maxX - minX + 1) * (maxY - minY + 1)) /
                  (canvas.width * canvas.height);
              }
            } catch (error) {
              pixelCheckError = String(error);
            }
            const style = getComputedStyle(image);
            const rect = image.getBoundingClientRect();
            const intrinsicRatio = image.naturalWidth / image.naturalHeight;
            const contentRatio = image.clientWidth / image.clientHeight;
            const boxRatio = rect.width / rect.height;
            const contentRatioRelativeError = Math.abs(contentRatio - intrinsicRatio) /
              intrinsicRatio;
            const boxRatioRelativeError = Math.abs(boxRatio - intrinsicRatio) /
              intrinsicRatio;
            const sourceAllowed = allowSvgFixture
              ? image.currentSrc.startsWith('data:image/')
              : image.currentSrc.startsWith('data:image/png;base64,');
            const deferredByClosedDetails = Boolean(
              image.closest('details:not([open])')
            );
            const ancestorViolations = [];
            const imageMaskImage = style.maskImage || style.webkitMaskImage || 'none';
            let ancestor = image.parentElement;
            while (ancestor && ancestor !== document.documentElement) {
              const ancestorStyle = getComputedStyle(ancestor);
              const ancestorMaskImage = ancestorStyle.maskImage ||
                ancestorStyle.webkitMaskImage || 'none';
              if (
                ancestorStyle.display === 'none' ||
                ['hidden', 'collapse'].includes(ancestorStyle.visibility) ||
                Number(ancestorStyle.opacity) < 0.99 ||
                ancestorStyle.filter !== 'none' ||
                ancestorStyle.clipPath !== 'none' ||
                ancestorMaskImage !== 'none' ||
                ancestorStyle.mixBlendMode !== 'normal'
              ) {
                ancestorViolations.push({
                  element: label(ancestor),
                  display: ancestorStyle.display,
                  visibility: ancestorStyle.visibility,
                  opacity: ancestorStyle.opacity,
                  filter: ancestorStyle.filter,
                  clipPath: ancestorStyle.clipPath,
                  maskImage: ancestorMaskImage,
                  mixBlendMode: ancestorStyle.mixBlendMode,
                });
              }
              ancestor = ancestor.parentElement;
            }
            const minimumUniqueColors = 2;
            const densePixelContentPass = nonDominantPixelFraction >=
              imageThresholds.denseNonDominantPixelFractionMin;
            const sparseScientificContentPass = nonDominantPixelFraction >=
                imageThresholds.sparseNonDominantPixelFractionMin &&
              nonDominantBboxFraction >=
                imageThresholds.sparseNonDominantBboxFractionMin &&
              nonDominantGridFraction >=
                imageThresholds.sparseNonDominantGridFractionMin &&
              uniqueOpaqueColors >= imageThresholds.sparseMinimumUniqueOpaqueColors;
            const pixelContentPass = densePixelContentPass || sparseScientificContentPass;
            const visualPass = image.clientWidth > 0 && image.clientHeight > 0 &&
              contentRatioRelativeError <= imageThresholds.contentRatioRelativeErrorMax &&
              Number(style.opacity) >= 0.99 && style.filter === 'none' &&
              style.clipPath === 'none' && style.transform === 'none' &&
              imageMaskImage === 'none' && style.mixBlendMode === 'normal' &&
              style.display !== 'none' && style.visibility === 'visible' &&
              alphaCoverage >= 0.25 && uniqueOpaqueColors >= minimumUniqueColors &&
              pixelContentPass &&
              ancestorViolations.length === 0 && pixelCheckError === null;
            return {
              index,
              sourceKind: image.currentSrc.split(',', 1)[0],
              sourceLength: image.currentSrc.length,
              sourceAllowed,
              complete: image.complete,
              naturalWidth: image.naturalWidth,
              naturalHeight: image.naturalHeight,
              rasterWidth: canvas.width,
              rasterHeight: canvas.height,
              rasterScale: round(rasterScale),
              imageSmoothingEnabled: context.imageSmoothingEnabled,
              clientWidth: image.clientWidth,
              clientHeight: image.clientHeight,
              boxWidth: round(rect.width),
              boxHeight: round(rect.height),
              intrinsicRatio: round(intrinsicRatio),
              contentRatio: round(contentRatio),
              boxRatio: round(boxRatio),
              contentRatioRelativeError: round(contentRatioRelativeError),
              boxRatioRelativeError: round(boxRatioRelativeError),
              opacity: style.opacity,
              filter: style.filter,
              clipPath: style.clipPath,
              maskImage: imageMaskImage,
              mixBlendMode: style.mixBlendMode,
              transform: style.transform,
              display: style.display,
              visibility: style.visibility,
              nonTransparentPixels,
              uniqueOpaqueColors,
              minimumUniqueColors,
              alphaCoverage: round(alphaCoverage),
              dominantOpaqueFraction: round(dominantOpaqueFraction),
              nonDominantPixelFraction: round(nonDominantPixelFraction),
              nonDominantBboxFraction: round(nonDominantBboxFraction),
              nonDominantOccupiedGridCells,
              nonDominantGridFraction: round(nonDominantGridFraction),
              densePixelContentPass,
              sparseScientificContentPass,
              pixelContentPass,
              pixelContentMode: densePixelContentPass
                ? 'dense'
                : sparseScientificContentPass
                ? 'sparse_spatially_distributed'
                : 'insufficient',
              deferredByClosedDetails,
              ancestorViolations,
              pixelCheckError,
              pass: sourceAllowed && image.complete && image.naturalWidth > 0 &&
                image.naturalHeight > 0 && (deferredByClosedDetails || visualPass),
            };
          });

          const allTables = Array.from(document.querySelectorAll('table'));
          const tables = Array.from(document.querySelectorAll('.table-wrap > table'))
            .map((table, index) => {
              const wrapper = table.parentElement;
              const overflowX = getComputedStyle(wrapper).overflowX;
              const before = wrapper.scrollLeft;
              wrapper.scrollLeft = Math.max(0, wrapper.scrollWidth - wrapper.clientWidth);
              const after = wrapper.scrollLeft;
              wrapper.scrollLeft = before;
              const hasOverflow = wrapper.scrollWidth > wrapper.clientWidth + 1;
              const deferredByClosedDetails = Boolean(
                table.closest('details:not([open])')
              );
              const visualPass = printMedia
                ? overflowX === 'visible' && !hasOverflow
                : ['auto', 'scroll'].includes(overflowX) &&
                  (!mobile || (hasOverflow && after > before + 1));
              return {
                index,
                overflowX,
                clientWidth: wrapper.clientWidth,
                scrollWidth: wrapper.scrollWidth,
                hasOverflow,
                scrollMoved: after > before + 1,
                deferredByClosedDetails,
                pass: deferredByClosedDetails || visualPass,
              };
            });
          const unwrappedTables = allTables
            .filter((table) => !table.matches('.table-wrap > table'))
            .map(label);

          const anchors = Array.from(document.querySelectorAll('a[href^="#"]'))
            .map((anchor, index) => {
              const href = anchor.getAttribute('href');
              let target = null;
              let targetId = null;
              try {
                targetId = href && href.length > 1
                  ? decodeURIComponent(href.slice(1))
                  : null;
                target = targetId ? document.getElementById(targetId) : null;
              } catch (error) {
                target = null;
              }
              const idCount = targetId
                ? Array.from(document.querySelectorAll('[id]'))
                    .filter((element) => element.id === targetId).length
                : 0;
              return { index, href, targetId, targetExists: Boolean(target), idCount };
            });
          const scriptCount = document.querySelectorAll('script').length;
          const resourceViolations = [];
          for (const element of document.querySelectorAll('*')) {
            for (const attribute of Array.from(element.attributes || [])) {
              const name = attribute.name.toLowerCase();
              const value = attribute.value.trim().toLowerCase();
              if (name.startsWith('on')) {
                resourceViolations.push({
                  element: label(element),
                  attribute: name,
                  reason: 'inline_event_handler',
                });
              }
              if (value.startsWith('javascript:') || value.startsWith('vbscript:')) {
                resourceViolations.push({
                  element: label(element),
                  attribute: name,
                  reason: 'active_url_scheme',
                });
              }
            }
          }
          for (const image of images) {
            if (!image.sourceAllowed) {
              resourceViolations.push({ element: `img[${image.index}]`, reason: 'image_source' });
            }
          }
          const resourceSelectors = [
            ['link[href]', 'href'],
            ['script[src]', 'src'],
            ['source[src]', 'src'],
            ['video[src]', 'src'],
            ['audio[src]', 'src'],
            ['iframe[src]', 'src'],
            ['iframe[srcdoc]', 'srcdoc'],
            ['object[data]', 'data'],
            ['embed[src]', 'src'],
            ['track[src]', 'src'],
            ['base[href]', 'href'],
            ['[srcset]', 'srcset'],
            ['[imagesrcset]', 'imagesrcset'],
            ['form[action]', 'action'],
            ['meta[http-equiv="refresh"]', 'content'],
            ['a[target]:not([target="_self"])', 'target'],
          ];
          for (const [selector, attribute] of resourceSelectors) {
            for (const element of document.querySelectorAll(selector)) {
              resourceViolations.push({
                element: label(element),
                attribute,
                value: element.getAttribute(attribute),
              });
            }
          }
          for (const element of document.querySelectorAll('style, [style]')) {
            const cssText = element.matches('style')
              ? element.textContent || ''
              : element.getAttribute('style') || '';
            if (/url\s*\(|@import/i.test(cssText)) {
              resourceViolations.push({
                element: label(element),
                reason: 'css_external_reference_syntax',
              });
            }
          }
          for (const sheet of Array.from(document.styleSheets)) {
            try {
              for (const rule of Array.from(sheet.cssRules || [])) {
                if (/url\s*\(|@import/i.test(rule.cssText || '')) {
                  resourceViolations.push({
                    element: 'stylesheet.cssRules',
                    reason: 'css_external_reference_rule',
                  });
                }
              }
            } catch (error) {
              resourceViolations.push({
                element: 'stylesheet.cssRules',
                reason: 'css_rules_unreadable',
              });
            }
          }
          const distinctAnchorTargets = new Set(
            anchors.map((entry) => entry.targetId).filter(Boolean)
          ).size;

          return {
            viewport: {
              width: window.innerWidth,
              height: window.innerHeight,
              documentClientWidth: documentElement.clientWidth,
              documentScrollWidth: documentElement.scrollWidth,
              bodyClientWidth: document.body.clientWidth,
              bodyScrollWidth: document.body.scrollWidth,
            },
            expected: {
              images: expectedImages,
              tables: expectedTables,
              wrappedTables: expectedTables,
              anchors: expectedAnchors,
              distinctAnchorTargets: expectedAnchors,
              scripts: 0,
            },
            counts: {
              images: images.length,
              tables: allTables.length,
              wrappedTables: tables.length,
              anchors: anchors.length,
              distinctAnchorTargets,
              scripts: scriptCount,
            },
            outOfBounds,
            unexpectedPositionedElements,
            unexpectedPositionedPseudoElements,
            overlapPairs,
            images,
            tables,
            unwrappedTables,
            anchors,
            resourceViolations,
            media: {
              screen: matchMedia('screen').matches,
              print: matchMedia('print').matches,
            },
          };
        }
        """,
        {
            "expectedImages": EXPECTED_IMAGE_COUNT,
            "expectedTables": EXPECTED_TABLE_COUNT,
            "expectedAnchors": EXPECTED_ANCHOR_COUNT,
            "mobile": mobile,
            "printMedia": print_media,
            "allowSvgFixture": allow_svg_fixture,
            "imageThresholds": {
                "contentRatioRelativeErrorMax": (
                    IMAGE_CONTENT_RATIO_RELATIVE_ERROR_MAX
                ),
                "denseNonDominantPixelFractionMin": (
                    IMAGE_DENSE_NON_DOMINANT_PIXEL_FRACTION_MIN
                ),
                "sparseNonDominantPixelFractionMin": (
                    IMAGE_SPARSE_NON_DOMINANT_PIXEL_FRACTION_MIN
                ),
                "sparseNonDominantBboxFractionMin": (
                    IMAGE_SPARSE_NON_DOMINANT_BBOX_FRACTION_MIN
                ),
                "sparseNonDominantGridFractionMin": (
                    IMAGE_SPARSE_NON_DOMINANT_GRID_FRACTION_MIN
                ),
                "sparseMinimumUniqueOpaqueColors": (
                    IMAGE_SPARSE_MINIMUM_UNIQUE_OPAQUE_COLORS
                ),
                "spatialGridSize": IMAGE_SPATIAL_GRID_SIZE,
            },
        },
    )
    return result


def validate_anchors(
    page: Page,
    *,
    mobile: bool,
    print_media: bool,
    allow_svg_fixture: bool,
) -> list[dict[str, Any]]:
    outcomes: list[dict[str, Any]] = []
    anchor_count = page.locator('a[href^="#"]').count()
    for index in range(anchor_count):
        page.evaluate("window.scrollTo(0, 0)")
        locator = page.locator('a[href^="#"]').nth(index)
        href = locator.get_attribute("href")
        if not href:
            outcomes.append({"index": index, "href": href, "pass": False})
            continue
        locator.focus()
        page.keyboard.press("Enter")
        page.wait_for_timeout(700)
        outcome = page.evaluate(
            """
            ({ href }) => {
              let target = null;
              let targetId = null;
              try {
                targetId = href.length > 1 ? decodeURIComponent(href.slice(1)) : null;
                target = targetId ? document.getElementById(targetId) : null;
              } catch (error) {
                target = null;
              }
              if (!target) {
                return { href, targetExists: false, pass: false };
              }
              const heading = target.matches('h1,h2,h3')
                ? target
                : target.querySelector('h1,h2,h3');
              const hitTarget = heading || target;
              const rect = hitTarget.getBoundingClientRect();
              const obstructionBottom = Math.max(0, ...Array.from(document.querySelectorAll('body *'))
                .filter((element) => element !== hitTarget && !element.contains(hitTarget))
                .filter((element) => {
                  const style = getComputedStyle(element);
                  const elementRect = element.getBoundingClientRect();
                  return ['fixed', 'sticky'].includes(style.position) &&
                    style.display !== 'none' && style.visibility !== 'hidden' &&
                    Number(style.opacity) > 0 && elementRect.width > 0 &&
                    elementRect.height > 0 && elementRect.top <= 1 &&
                    elementRect.bottom > 0;
                })
                .map((element) => element.getBoundingClientRect().bottom));
              const expectedHash = `#${encodeURIComponent(targetId)}`;
              const hashMatches = decodeURIComponent(location.hash) ===
                decodeURIComponent(expectedHash);
              const targetVisible = rect.bottom > obstructionBottom - 1 &&
                rect.top < window.innerHeight - 1 && rect.width > 0 && rect.height > 0;
              const notOccluded = rect.top >= obstructionBottom - 1;
              const x = Math.min(window.innerWidth - 1, Math.max(0, rect.left + rect.width / 2));
              const y = Math.min(window.innerHeight - 1,
                Math.max(obstructionBottom + 1, rect.top + Math.min(rect.height / 2, 12)));
              const hitStack = document.elementsFromPoint(x, y);
              const topHit = hitStack[0] || null;
              const hitTestPass = Boolean(topHit) &&
                (topHit === hitTarget || hitTarget.contains(topHit));
              return {
                href,
                targetId,
                targetExists: true,
                hash: location.hash,
                hashMatches,
                targetTop: Math.round(rect.top * 1000) / 1000,
                targetBottom: Math.round(rect.bottom * 1000) / 1000,
                obstructionBottom: Math.round(obstructionBottom * 1000) / 1000,
                targetVisible,
                notOccluded,
                hitTestPoint: { x: Math.round(x * 1000) / 1000, y: Math.round(y * 1000) / 1000 },
                hitStack: hitStack.slice(0, 5).map((element) =>
                  element.id ? `${element.tagName.toLowerCase()}#${element.id}` :
                    element.tagName.toLowerCase()
                ),
                topHit: topHit
                  ? (topHit.id ? `${topHit.tagName.toLowerCase()}#${topHit.id}` :
                    topHit.tagName.toLowerCase())
                  : null,
                hitTestPass,
                pass: hashMatches && targetVisible && notOccluded && hitTestPass,
              };
            }
            """,
            {"href": href},
        )
        outcome["index"] = index
        outcome["layout"] = collect_layout(
            page,
            mobile=mobile,
            print_media=print_media,
            allow_svg_fixture=allow_svg_fixture,
        )
        outcome["layout_pass"] = layout_contract_pass(outcome["layout"])
        outcome["pass"] = outcome.get("pass") is True and outcome["layout_pass"]
        outcomes.append(outcome)
    page.evaluate("window.scrollTo(0, 0)")
    page.wait_for_timeout(700)
    return outcomes


def validate_image_hit_tests(page: Page) -> list[dict[str, Any]]:
    outcomes: list[dict[str, Any]] = []
    images = page.locator("img")
    for index in range(images.count()):
        locator = images.nth(index)
        locator.scroll_into_view_if_needed()
        page.wait_for_timeout(100)
        outcome = locator.evaluate(
            """
            (image) => {
              const rect = image.getBoundingClientRect();
              const points = [];
              for (const xFraction of [0.25, 0.5, 0.75]) {
                for (const yFraction of [0.25, 0.5, 0.75]) {
                  const x = Math.min(window.innerWidth - 1,
                    Math.max(0, rect.left + rect.width * xFraction));
                  const y = Math.min(window.innerHeight - 1,
                    Math.max(0, rect.top + rect.height * yFraction));
                  const stack = document.elementsFromPoint(x, y);
                  const topHit = stack[0] || null;
                  const hit = Boolean(topHit) &&
                    (topHit === image || image.contains(topHit));
                  points.push({
                    xFraction,
                    yFraction,
                    point: {
                      x: Math.round(x * 1000) / 1000,
                      y: Math.round(y * 1000) / 1000,
                    },
                    topHit: topHit
                      ? (topHit.id ? `${topHit.tagName.toLowerCase()}#${topHit.id}` :
                        topHit.tagName.toLowerCase())
                      : null,
                    hitStack: stack.slice(0, 5).map((element) =>
                      element.id ? `${element.tagName.toLowerCase()}#${element.id}` :
                        element.tagName.toLowerCase()
                    ),
                    pass: hit,
                  });
                }
              }
              return {
                rect: {
                  left: Math.round(rect.left * 1000) / 1000,
                  top: Math.round(rect.top * 1000) / 1000,
                  width: Math.round(rect.width * 1000) / 1000,
                  height: Math.round(rect.height * 1000) / 1000,
                },
                points,
                pass: rect.width > 0 && rect.height > 0 &&
                  points.every((point) => point.pass),
              };
            }
            """
        )
        outcome["index"] = index
        outcomes.append(outcome)
    page.evaluate("window.scrollTo(0, 0)")
    page.wait_for_timeout(100)
    return outcomes


def validate_details_states(
    page: Page,
    *,
    mobile: bool,
    print_media: bool,
    allow_svg_fixture: bool,
) -> list[dict[str, Any]]:
    outcomes: list[dict[str, Any]] = []
    details = page.locator("details")
    for index in range(details.count()):
        locator = details.nth(index)
        initial_open = locator.evaluate("element => element.open")
        initial_layout = collect_layout(
            page,
            mobile=mobile,
            print_media=print_media,
            allow_svg_fixture=allow_svg_fixture,
        )
        locator.scroll_into_view_if_needed()
        locator.locator("summary").evaluate("summary => summary.click()")
        page.wait_for_timeout(250)
        toggled_open = locator.evaluate("element => element.open")
        toggled_layout = collect_layout(
            page,
            mobile=mobile,
            print_media=print_media,
            allow_svg_fixture=allow_svg_fixture,
        )
        locator.evaluate("(element, value) => { element.open = value; }", initial_open)
        page.wait_for_timeout(250)
        restored_open = locator.evaluate("element => element.open")
        restored_layout = collect_layout(
            page,
            mobile=mobile,
            print_media=print_media,
            allow_svg_fixture=allow_svg_fixture,
        )
        outcomes.append(
            {
                "index": index,
                "initial_open": initial_open,
                "toggled_open": toggled_open,
                "restored_open": restored_open,
                "initial_layout": initial_layout,
                "toggled_layout": toggled_layout,
                "restored_layout": restored_layout,
                "pass": all(
                    (
                        toggled_open is not initial_open,
                        restored_open is initial_open,
                        layout_contract_pass(initial_layout),
                        layout_contract_pass(toggled_layout),
                        layout_contract_pass(restored_layout),
                    )
                ),
            }
        )
    page.evaluate("window.scrollTo(0, 0)")
    page.wait_for_timeout(250)
    return outcomes


def layout_contract_pass(result: dict[str, Any]) -> bool:
    viewport = result["viewport"]
    return all(
        (
            result["counts"] == result["expected"],
            viewport["documentScrollWidth"] <= viewport["documentClientWidth"] + 1,
            viewport["bodyScrollWidth"] <= viewport["bodyClientWidth"] + 1,
            not result["outOfBounds"],
            not result["unexpectedPositionedElements"],
            not result["unexpectedPositionedPseudoElements"],
            not result["overlapPairs"],
            not result["unwrappedTables"],
            not result["resourceViolations"],
            all(entry["pass"] for entry in result["images"]),
            all(entry["pass"] for entry in result["tables"]),
            all(
                entry["targetExists"] and entry["idCount"] == 1
                for entry in result["anchors"]
            ),
        )
    )


def profile_pass(
    result: dict[str, Any],
    *,
    console_errors: list[str],
    page_errors: list[str],
    network_attempts: list[str],
    websocket_attempts: list[str],
    popup_attempts: list[str],
    anchor_outcomes: list[dict[str, Any]] | None,
    details_outcomes: list[dict[str, Any]],
    image_hit_tests: list[dict[str, Any]],
    open_result: dict[str, Any],
    post_interaction_result: dict[str, Any],
    pdf_verification: dict[str, Any] | None,
    expect_print: bool,
) -> bool:
    anchor_interaction_pass = anchor_outcomes is None or all(
        entry.get("pass") is True for entry in anchor_outcomes
    )
    return all(
        (
            layout_contract_pass(result),
            layout_contract_pass(open_result),
            layout_contract_pass(post_interaction_result),
            anchor_interaction_pass,
            all(entry.get("pass") is True for entry in details_outcomes),
            all(entry.get("pass") is True for entry in image_hit_tests),
            pdf_verification is None or pdf_verification.get("pass") is True,
            result["media"]["print"] is expect_print,
            open_result["media"]["print"] is expect_print,
            post_interaction_result["media"]["print"] is expect_print,
            not console_errors,
            not page_errors,
            not network_attempts,
            not websocket_attempts,
            not popup_attempts,
        )
    )


def abort_route(route: Route) -> None:
    route.abort("blockedbyclient")


def reject_websocket(route: WebSocketRoute, attempts: list[str]) -> None:
    attempts.append(route.url)
    route.close(code=1008, reason="offline QA policy")


def verify_pdf(path: Path, staging: Path) -> dict[str, Any]:
    info_result = subprocess.run(
        [str(PDFINFO), str(path)],
        env={"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C"},
        check=True,
        capture_output=True,
        text=True,
    )
    fields: dict[str, str] = {}
    for line in info_result.stdout.splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            fields[key.strip()] = value.strip()
    pages = int(fields.get("Pages", "0"))
    page_size = fields.get("Page size", "")
    size_match = re.search(
        r"([0-9]+(?:\.[0-9]+)?)\s+x\s+([0-9]+(?:\.[0-9]+)?)\s+pts",
        page_size,
    )
    width_points = float(size_match.group(1)) if size_match else 0.0
    height_points = float(size_match.group(2)) if size_match else 0.0
    page_size_a4 = (
        abs(width_points - 595.28) <= 2.0
        and abs(height_points - 841.89) <= 2.0
    )

    raster_prefix = staging / "print_a4_page"
    raster_command = [
        str(PDFTOPPM),
        "-png",
        "-r",
        "96",
        str(path),
        str(raster_prefix),
    ]
    raster_result = subprocess.run(
        raster_command,
        env={"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C"},
        check=True,
        capture_output=True,
        text=True,
    )

    def page_number(page_path: Path) -> int:
        match = re.search(r"-(\d+)\.png$", page_path.name)
        return int(match.group(1)) if match else 0

    raster_pages: list[dict[str, Any]] = []
    for raster_path in sorted(staging.glob("print_a4_page-*.png"), key=page_number):
        current_page = page_number(raster_path)
        page_text_command = [
            str(PDFTOTEXT),
            "-f",
            str(current_page),
            "-l",
            str(current_page),
            str(path),
            "-",
        ]
        page_text_result = subprocess.run(
            page_text_command,
            env={"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C"},
            check=True,
            capture_output=True,
        )
        page_text = page_text_result.stdout.decode("utf-8", errors="replace")
        page_nonwhitespace_characters = len(re.sub(r"\s+", "", page_text))
        page_text_pass = page_nonwhitespace_characters >= 20
        raster_path.chmod(0o444)
        fsync_file(raster_path)
        with Image.open(raster_path) as source:
            source.load()
            rgb = source.convert("RGB")
            width, height = rgb.size
            white = Image.new("RGB", rgb.size, (255, 255, 255))
            difference = ImageChops.difference(rgb, white).convert("L")
            threshold = difference.point(lambda value: 255 if value > 4 else 0)
            histogram = threshold.histogram()
            nonwhite_pixels = histogram[255]
            nonwhite_fraction = nonwhite_pixels / (width * height)
            bbox = threshold.getbbox()
            grid_columns = 8
            grid_rows = 12
            occupied_grid_cells: list[list[int]] = []
            for grid_y in range(grid_rows):
                for grid_x in range(grid_columns):
                    cell_box = (
                        grid_x * width // grid_columns,
                        grid_y * height // grid_rows,
                        (grid_x + 1) * width // grid_columns,
                        (grid_y + 1) * height // grid_rows,
                    )
                    cell_nonwhite = threshold.crop(cell_box).histogram()[255]
                    if cell_nonwhite >= 20:
                        occupied_grid_cells.append([grid_x, grid_y])
            bbox_width_fraction = (
                (bbox[2] - bbox[0]) / width if bbox is not None else 0.0
            )
            bbox_height_fraction = (
                (bbox[3] - bbox[1]) / height if bbox is not None else 0.0
            )
            spatial_coverage_pass = (
                len(occupied_grid_cells) >= 12
                and bbox_width_fraction >= 0.4
                and bbox_height_fraction >= 0.05
            )
            margins = (
                {
                    "left": bbox[0],
                    "top": bbox[1],
                    "right": width - bbox[2],
                    "bottom": height - bbox[3],
                }
                if bbox
                else {"left": 0, "top": 0, "right": 0, "bottom": 0}
            )
            dimensions_pass = 790 <= width <= 798 and 1118 <= height <= 1128
            content_pass = (
                bbox is not None
                and nonwhite_fraction >= 0.005
                and spatial_coverage_pass
            )
            edge_pass = bbox is not None and min(margins.values()) >= 8
            raster_pages.append(
                {
                    "page": current_page,
                    "artifact": artifact(
                        raster_path,
                        public_path=Path("PLACEHOLDER") / raster_path.name,
                    ),
                    "width_px": width,
                    "height_px": height,
                    "nonwhite_pixels": nonwhite_pixels,
                    "nonwhite_fraction": round(nonwhite_fraction, 8),
                    "nonwhite_bbox": list(bbox) if bbox else None,
                    "bbox_width_fraction": round(bbox_width_fraction, 8),
                    "bbox_height_fraction": round(bbox_height_fraction, 8),
                    "grid_columns": grid_columns,
                    "grid_rows": grid_rows,
                    "occupied_grid_cells": occupied_grid_cells,
                    "occupied_grid_cell_count": len(occupied_grid_cells),
                    "spatial_coverage_pass": spatial_coverage_pass,
                    "content_margins_px": margins,
                    "dimensions_pass": dimensions_pass,
                    "content_pass": content_pass,
                    "edge_pass": edge_pass,
                    "text_command": page_text_command,
                    "text_sha256": hashlib.sha256(page_text_result.stdout).hexdigest(),
                    "text_size_bytes": len(page_text_result.stdout),
                    "nonwhitespace_text_characters": page_nonwhitespace_characters,
                    "text_pass": page_text_pass,
                    "pass": dimensions_pass and content_pass and edge_pass and page_text_pass,
                }
            )

    text_command = [str(PDFTOTEXT), str(path), "-"]
    text_result = subprocess.run(
        text_command,
        env={"PATH": "/usr/bin:/bin", "LC_ALL": "C", "LANG": "C"},
        check=True,
        capture_output=True,
    )
    extracted_text = text_result.stdout.decode("utf-8", errors="replace")
    nonwhitespace_characters = len(re.sub(r"\s+", "", extracted_text))
    text_pass = nonwhitespace_characters >= 100
    raster_count_pass = pages > 0 and len(raster_pages) == pages
    return {
        "command": [str(PDFINFO), str(path)],
        "pages": pages,
        "page_size": page_size,
        "page_width_points": width_points,
        "page_height_points": height_points,
        "page_size_a4": page_size_a4,
        "stdout": info_result.stdout,
        "raster_command": raster_command,
        "raster_stdout": raster_result.stdout,
        "raster_stderr": raster_result.stderr,
        "raster_pages": raster_pages,
        "raster_page_count_matches_pdf": raster_count_pass,
        "text_command": text_command,
        "text_sha256": hashlib.sha256(text_result.stdout).hexdigest(),
        "text_size_bytes": len(text_result.stdout),
        "nonwhitespace_text_characters": nonwhitespace_characters,
        "text_pass": text_pass,
        "pass": all(
            (
                pages > 0,
                page_size_a4,
                raster_count_pass,
                all(page["pass"] for page in raster_pages),
                text_pass,
            )
        ),
    }


def run_profile(
    browser: Any,
    *,
    html_path: Path,
    staging: Path,
    name: str,
    width: int,
    height: int,
    wait_ms: int,
    print_media: bool,
    evidence_kind: str,
) -> dict[str, Any]:
    context: BrowserContext = browser.new_context(
        viewport={"width": width, "height": height},
        service_workers="block",
    )
    context.set_offline(True)
    page = context.new_page()
    console_errors: list[str] = []
    page_errors: list[str] = []
    network_attempts: list[str] = []
    websocket_attempts: list[str] = []
    popup_attempts: list[str] = []
    page.on(
        "console",
        lambda message: console_errors.append(message.text)
        if message.type == "error"
        else None,
    )
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    main_document_url = html_path.as_uri()
    context.on(
        "request",
        lambda request: network_attempts.append(request.url)
        if request.url != main_document_url
        else None,
    )
    context.on("page", lambda popup: popup_attempts.append(popup.url) if popup is not page else None)

    def route_resource(route: Route) -> None:
        if route.request.url == main_document_url:
            route.continue_()
        else:
            route.abort("blockedbyclient")

    context.route("**/*", route_resource)
    context.route_web_socket(
        "**/*", lambda route: reject_websocket(route, websocket_attempts)
    )
    if print_media:
        page.emulate_media(media="print")
    page.goto(main_document_url, wait_until="load")
    page.wait_for_function(
        "Array.from(document.images).every((image) => image.complete)",
        timeout=120_000,
    )
    page.wait_for_timeout(wait_ms)
    mobile = width <= 390 and not print_media
    if print_media:
        page.eval_on_selector_all(
            "details", "elements => elements.forEach((element) => element.open = true)"
        )
    page.evaluate("window.scrollTo(0, 0)")
    allow_svg_fixture = evidence_kind == "self_test"
    diagnostics = collect_layout(
        page,
        mobile=mobile,
        print_media=print_media,
        allow_svg_fixture=allow_svg_fixture,
    )
    details_states = (
        []
        if print_media
        else validate_details_states(
            page,
            mobile=mobile,
            print_media=print_media,
            allow_svg_fixture=allow_svg_fixture,
        )
    )
    page.eval_on_selector_all(
        "details", "elements => elements.forEach((element) => element.open = true)"
    )
    page.evaluate("window.scrollTo(0, 0)")
    page.wait_for_timeout(250)
    open_diagnostics = collect_layout(
        page,
        mobile=mobile,
        print_media=print_media,
        allow_svg_fixture=allow_svg_fixture,
    )
    anchors = (
        None
        if print_media
        else validate_anchors(
            page,
            mobile=mobile,
            print_media=print_media,
            allow_svg_fixture=allow_svg_fixture,
        )
    )
    image_hit_tests = validate_image_hit_tests(page)
    page.evaluate("window.scrollTo(0, 0)")
    page.wait_for_timeout(100)
    post_interaction_diagnostics = collect_layout(
        page,
        mobile=mobile,
        print_media=print_media,
        allow_svg_fixture=allow_svg_fixture,
    )

    screenshot_path = staging / f"{name}.png"
    page.screenshot(path=str(screenshot_path), full_page=True)
    screenshot_path.chmod(0o444)
    fsync_file(screenshot_path)
    if screenshot_path.stat().st_size <= 10_000:
        raise QaError(f"Screenshot is unexpectedly small: {screenshot_path}")

    pdf_path: Path | None = None
    pdf_verification: dict[str, Any] | None = None
    if print_media:
        pdf_path = staging / "print_a4.pdf"
        page.pdf(
            path=str(pdf_path),
            format="A4",
            print_background=True,
            prefer_css_page_size=False,
        )
        pdf_path.chmod(0o444)
        fsync_file(pdf_path)
        if pdf_path.stat().st_size <= 10_000 or pdf_path.read_bytes()[:5] != b"%PDF-":
            raise QaError(f"Print PDF is missing or invalid: {pdf_path}")
        pdf_verification = verify_pdf(pdf_path, staging)

    passed = profile_pass(
        diagnostics,
        console_errors=console_errors,
        page_errors=page_errors,
        network_attempts=network_attempts,
        websocket_attempts=websocket_attempts,
        popup_attempts=popup_attempts,
        anchor_outcomes=anchors,
        details_outcomes=details_states,
        image_hit_tests=image_hit_tests,
        open_result=open_diagnostics,
        post_interaction_result=post_interaction_diagnostics,
        pdf_verification=pdf_verification,
        expect_print=print_media,
    )
    result = {
        "name": name,
        "width": width,
        "height": height,
        "media": "print" if print_media else "screen",
        "diagnostics": diagnostics,
        "details_state_interactions": details_states,
        "open_diagnostics": open_diagnostics,
        "post_interaction_diagnostics": post_interaction_diagnostics,
        "anchor_interactions": anchors,
        "image_hit_tests": image_hit_tests,
        "console_errors": console_errors,
        "page_errors": page_errors,
        "network_attempts": network_attempts,
        "websocket_attempts": websocket_attempts,
        "popup_attempts": popup_attempts,
        "screenshot": artifact(
            screenshot_path,
            public_path=Path("PLACEHOLDER") / screenshot_path.name,
        ),
        "pdf": artifact(
            pdf_path,
            public_path=Path("PLACEHOLDER") / pdf_path.name,
        )
        if pdf_path
        else None,
        "pdf_verification": pdf_verification,
        "pass": passed,
    }
    context.close()
    return result


def replace_public_paths(payload: dict[str, Any], output_dir: Path) -> None:
    for key in ("executed_html", "executed_html_after"):
        executed_html = payload.get(key)
        if isinstance(executed_html, dict):
            executed_html["path"] = str(
                (output_dir / Path(executed_html["path"]).name).resolve(strict=False)
            )
    for profile in payload["profiles"]:
        profile["screenshot"]["path"] = str(
            (output_dir / Path(profile["screenshot"]["path"]).name).resolve(strict=False)
        )
        if profile["pdf"]:
            profile["pdf"]["path"] = str(
                (output_dir / Path(profile["pdf"]["path"]).name).resolve(strict=False)
            )
        pdf_verification = profile.get("pdf_verification")
        if pdf_verification:
            pdf_path = str((output_dir / "print_a4.pdf").resolve(strict=False))
            pdf_verification["command"][1] = pdf_path
            pdf_verification["raster_command"][4] = pdf_path
            pdf_verification["raster_command"][5] = str(
                (output_dir / "print_a4_page").resolve(strict=False)
            )
            pdf_verification["text_command"][1] = pdf_path
            for raster_page in pdf_verification["raster_pages"]:
                record = raster_page["artifact"]
                record["path"] = str(
                    (output_dir / Path(record["path"]).name).resolve(strict=False)
                )
                raster_page["text_command"][5] = pdf_path


def publish_failure(source_dir: Path, output_dir: Path, payload: dict[str, Any]) -> Path:
    source_dir.chmod(0o700)
    failure_name = (
        f"{output_dir.name}.failed."
        f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}."
        f"{uuid.uuid4().hex[:8]}"
    )
    failure_dir = output_dir.with_name(failure_name)
    replace_public_paths(payload, failure_dir)
    payload["pass"] = False
    quarantine_dir: Path | None = None
    quarantined: list[dict[str, Any]] = []
    success_artifact_names = ("qa_receipt.json", "_SUCCESS.json")
    present_success_artifacts = [
        source_dir / name
        for name in success_artifact_names
        if path_exists(source_dir / name)
    ]
    if present_success_artifacts:
        quarantine_dir = output_dir.with_name(
            f"{failure_name}.quarantined_success_artifacts"
        )
        quarantine_dir.mkdir(mode=0o700)
        for source_path in present_success_artifacts:
            observed = artifact(source_path)
            target_path = quarantine_dir / source_path.name
            rename_no_replace(source_path, target_path)
            observed["path"] = str(target_path.resolve(strict=False))
            quarantined.append(observed)
        fsync_directory(source_dir)
        quarantine_dir.chmod(0o555)
        fsync_directory(quarantine_dir)
        fsync_directory(quarantine_dir.parent)
    receipt_name = "failure_receipt.json"
    if path_exists(source_dir / receipt_name):
        raise QaError(f"Failure receipt already exists: {source_dir / receipt_name}")
    payload["failure_publication"] = {
        "directory": str(failure_dir.resolve(strict=False)),
        "receipt": str((failure_dir / receipt_name).resolve(strict=False)),
        "authoritative_receipt": receipt_name,
        "quarantined_success_artifacts_directory": (
            str(quarantine_dir.resolve(strict=False)) if quarantine_dir else None
        ),
        "quarantined_success_artifacts": quarantined,
        "success_markers_present_in_failure_directory": False,
    }
    write_json(source_dir / receipt_name, payload)
    source_dir.chmod(0o555)
    fsync_directory(source_dir)
    rename_no_replace(source_dir, failure_dir)
    fsync_directory(failure_dir.parent)
    if path_exists(output_dir):
        raise QaError(f"Failed QA left an authoritative output directory: {output_dir}")
    if any(path_exists(failure_dir / name) for name in success_artifact_names):
        raise QaError("Failed QA directory retained a success artifact")
    return failure_dir


def evidence_descriptor(evidence_kind: str) -> dict[str, Any]:
    if evidence_kind == "formal_report":
        return {
            "evidence_kind": evidence_kind,
            "task_type": "B_comprehensive_validation",
            "scope": "positional_singleton_supplemental_standalone_html",
            "validation_eligibility": "formal_report_QA_evidence",
        }
    return {
        "evidence_kind": evidence_kind,
        "task_type": "F_demo_illustration",
        "scope": "QA_harness_self_test_fixture_only",
        "validation_eligibility": "self_test_only_not_formal_report_evidence",
    }


def base_payload(
    *,
    args: argparse.Namespace,
    html_before: dict[str, Any],
    executable_before: dict[str, Any],
    script_before: dict[str, Any],
    git_head: str,
    profiles: list[dict[str, Any]],
    browser_version: str,
    runtime_before: dict[str, Any],
    runtime_after: dict[str, Any] | None,
    executed_html: dict[str, Any] | None,
    static_html_security: dict[str, Any] | None,
    formal_release_binding: dict[str, Any] | None,
) -> dict[str, Any]:
    return {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        **evidence_descriptor(args.evidence_kind),
        "command": [str(Path(sys.executable).resolve()), *sys.argv],
        "repository": {
            "root": str(Path(__file__).resolve().parents[3]),
            "git_head": git_head,
        },
        "qa_script": script_before,
        "input_html": html_before,
        "executed_html": executed_html,
        "static_html_security": static_html_security,
        "formal_release_binding": formal_release_binding,
        "runtime_identity": runtime_before,
        "runtime_identity_after": runtime_after,
        "runtime_identity_stable": runtime_after == runtime_before,
        "verification_trust_model": {
            "environment": "trusted_local_scientific_QA_without_malicious_same_uid_concurrent_mutation",
            "executed_html_scope": (
                "browser_reads_a_private_read_only_copy_whose_size_and_sha256_equal_the_"
                "release_bound_input_HTML"
            ),
            "runtime_scope": (
                "complete_direct_runtime_bundle_identity_before_after_stability_check;_"
                "not_a_claim_of_kernel_enforced_executed_byte_attestation_or_resistance_"
                "to_malicious_same_uid_ABA_replacement"
            ),
            "independent_review_required_for_formal_release": True,
        },
        "browser": {
            "executable": executable_before,
            "version": browser_version,
            "launch": {
                "headless": True,
                "ignore_default_args": ["--hide-scrollbars"],
                "additional_args": [],
            },
            "context_policy": {
                "offline": True,
                "service_workers": "block",
                "http_https_routes": "abort_blockedbyclient",
                "websocket_routes": "close_1008",
                "popup_attempts_allowed": 0,
            },
        },
        "profiles": profiles,
    }


def main() -> None:
    args = parse_args()
    if args.wait_ms < 0:
        raise SystemExit("--wait-ms must be non-negative")
    html_path = args.html.resolve(strict=True)
    executable_path = args.executable_path.resolve(strict=True)
    output_dir = args.output_dir.resolve(strict=False)
    if not html_path.is_file():
        raise SystemExit(f"Missing HTML: {html_path}")
    if not executable_path.is_file():
        raise SystemExit(f"Missing Chromium executable: {executable_path}")
    if path_exists(output_dir):
        raise SystemExit(f"Refusing to overwrite output: {output_dir}")
    repo_root = Path(__file__).resolve().parents[3]
    html_bytes, html_before = read_file_snapshot(html_path)
    executable_before = artifact(executable_path)
    script_before = artifact(Path(__file__).resolve())
    git_head = repository_head(repo_root)
    formal_paths = resolve_formal_release_paths(args)
    formal_binding_before = validate_formal_release_binding(
        formal_paths,
        html_identity=html_before,
    )
    runtime_before = runtime_identity(executable_path)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = output_dir.with_name(
        f".{output_dir.name}.staging.{os.getpid()}.{uuid.uuid4().hex}"
    )
    started = time.monotonic()
    profiles: list[dict[str, Any]] = []
    browser_version = ""
    runtime_after: dict[str, Any] | None = None
    executed_html_identity: dict[str, Any] | None = None
    static_html_security: dict[str, Any] | None = None
    html_after: dict[str, Any] | None = None
    executed_html_after: dict[str, Any] | None = None
    executable_after: dict[str, Any] | None = None
    script_after: dict[str, Any] | None = None
    formal_binding_after: dict[str, Any] | None = None
    payload: dict[str, Any] | None = None
    executed_html_path = staging / "executed_input.html"
    try:
        staging.mkdir(mode=0o700)
        write_bytes_exclusive(executed_html_path, html_bytes)
        executed_html_identity = artifact(
            executed_html_path,
            public_path=output_dir / executed_html_path.name,
        )
        if (
            executed_html_identity["size_bytes"] != html_before["size_bytes"]
            or executed_html_identity["sha256"] != html_before["sha256"]
        ):
            raise QaError("Private browser input copy differs from the signed source bytes")
        static_html_security = static_html_security_scan(html_bytes)
        if static_html_security["pass"] is not True:
            raise QaError("Static HTML security scan rejected active content")
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(
                headless=True,
                executable_path=str(executable_path),
                ignore_default_args=["--hide-scrollbars"],
                args=[],
            )
            browser_version = browser.version
            for name, width, height in SCREEN_PROFILES:
                profiles.append(
                    run_profile(
                        browser,
                        html_path=executed_html_path,
                        staging=staging,
                        name=name,
                        width=width,
                        height=height,
                        wait_ms=args.wait_ms,
                        print_media=False,
                        evidence_kind=args.evidence_kind,
                    )
                )
            name, width, height = PRINT_PROFILE
            profiles.append(
                run_profile(
                    browser,
                    html_path=executed_html_path,
                    staging=staging,
                    name=name,
                    width=width,
                    height=height,
                    wait_ms=args.wait_ms,
                    print_media=True,
                    evidence_kind=args.evidence_kind,
                )
            )
            browser.close()
        _, html_after = read_file_snapshot(html_path)
        if html_after != html_before:
            raise QaError("HTML input identity changed during QA")
        executed_html_after = artifact(
            executed_html_path,
            public_path=output_dir / executed_html_path.name,
        )
        if executed_html_after != executed_html_identity:
            raise QaError("Private browser input copy changed during QA")
        executable_after = artifact(executable_path)
        if executable_after != executable_before:
            raise QaError("Chromium executable identity changed during QA")
        script_after = artifact(Path(__file__).resolve())
        if script_after != script_before:
            raise QaError("QA script identity changed during QA")
        formal_binding_after = validate_formal_release_binding(
            formal_paths,
            html_identity=html_after,
        )
        if formal_binding_after != formal_binding_before:
            raise QaError("Formal release binding changed during QA")
        runtime_after = runtime_identity(executable_path)
        if runtime_after != runtime_before:
            raise QaError("QA runtime bundle identity changed during execution")

        passed = all(profile["pass"] for profile in profiles)
        payload = base_payload(
            args=args,
            html_before=html_before,
            executable_before=executable_before,
            script_before=script_before,
            git_head=git_head,
            profiles=profiles,
            browser_version=browser_version,
            runtime_before=runtime_before,
            runtime_after=runtime_after,
            executed_html=executed_html_identity,
            static_html_security=static_html_security,
            formal_release_binding=formal_binding_before,
        )
        payload.update(
            {
                "contracts": {
                    "screen_profiles": [name for name, _, _ in SCREEN_PROFILES],
                    "print_profile": PRINT_PROFILE[0],
                    "expected_images": EXPECTED_IMAGE_COUNT,
                    "expected_tables": EXPECTED_TABLE_COUNT,
                    "expected_internal_anchors": EXPECTED_ANCHOR_COUNT,
                    "image_content_ratio_relative_error_max": (
                        IMAGE_CONTENT_RATIO_RELATIVE_ERROR_MAX
                    ),
                    "image_box_ratio": "informational_only_includes_CSS_border",
                    "image_dense_non_dominant_pixel_fraction_min": (
                        IMAGE_DENSE_NON_DOMINANT_PIXEL_FRACTION_MIN
                    ),
                    "image_sparse_non_dominant_pixel_fraction_min": (
                        IMAGE_SPARSE_NON_DOMINANT_PIXEL_FRACTION_MIN
                    ),
                    "image_sparse_non_dominant_bbox_fraction_min": (
                        IMAGE_SPARSE_NON_DOMINANT_BBOX_FRACTION_MIN
                    ),
                    "image_sparse_non_dominant_grid_fraction_min": (
                        IMAGE_SPARSE_NON_DOMINANT_GRID_FRACTION_MIN
                    ),
                    "image_sparse_minimum_unique_opaque_colors": (
                        IMAGE_SPARSE_MINIMUM_UNIQUE_OPAQUE_COLORS
                    ),
                    "image_spatial_grid_size": IMAGE_SPATIAL_GRID_SIZE,
                    "horizontal_overflow_tolerance_px": 1,
                    "overlap_tolerance_px": 1,
                    "network_policy": (
                        "offline_context_plus_zero_http_https_ws_wss_popup_attempts"
                    ),
                    "post_interaction_layout_revalidation": True,
                    "per_anchor_layout_revalidation": True,
                    "details_default_toggled_restored_revalidation": True,
                    "print_pdf_pdfinfo_a4_validation": True,
                    "print_pdf_raster_all_pages_validation": True,
                    "print_pdf_per_page_text_validation": True,
                    "private_executed_html_snapshot": True,
                    "static_active_content_rejection": True,
                    "full_runtime_bundle_before_after_stability_check": True,
                    "formal_release_signature_binding": (
                        args.evidence_kind == "formal_report"
                    ),
                },
                "input_html_after": html_after,
                "executed_html_after": executed_html_after,
                "browser_executable_after": executable_after,
                "qa_script_after": script_after,
                "formal_release_binding_after": formal_binding_after,
                "elapsed_seconds": round(time.monotonic() - started, 3),
                "pass_semantics": (
                    "all_screen_and_print_profiles_pass_pre_and_post_interaction_"
                    "layout_image_pixel_hit_test_table_anchor_error_offline_resource_"
                    "popup_and_pdf_structure_gates"
                ),
                "pass": passed,
            }
        )
        if not passed:
            raise QaError("One or more QA profiles failed")

        replace_public_paths(payload, output_dir)
        receipt_path = staging / "qa_receipt.json"
        write_json(receipt_path, payload)
        success_path = staging / "_SUCCESS.json"
        write_json(
            success_path,
            {
                "schema_name": "intersubmod.atomic_release_marker",
                "schema_version": "1.0.0",
                "receipt": artifact(
                    receipt_path,
                    public_path=output_dir / receipt_path.name,
                ),
                "pass": True,
            },
        )
        staging.chmod(0o555)
        fsync_directory(staging)
        rename_no_replace(staging, output_dir)
        fsync_directory(output_dir.parent)
        output_message = json.dumps(
            {
                "output_dir": str(output_dir),
                "receipt": str(output_dir / receipt_path.name),
                "profiles": len(profiles),
                "pass": True,
            },
            ensure_ascii=False,
            indent=2,
        )
        sys.stdout.write(output_message + "\n")
        sys.stdout.flush()
    except Exception as error:
        after_identity_capture_errors: dict[str, str] = {}

        def capture_after(label: str, callback: Any) -> Any:
            try:
                return callback()
            except Exception as capture_error:
                after_identity_capture_errors[label] = (
                    f"{type(capture_error).__name__}: {capture_error}"
                )
                return None

        if html_after is None:
            html_snapshot = capture_after(
                "input_html_after", lambda: read_file_snapshot(html_path)
            )
            if html_snapshot is not None:
                _, html_after = html_snapshot
        if executed_html_after is None and path_exists(executed_html_path):
            executed_html_after = capture_after(
                "executed_html_after",
                lambda: artifact(
                    executed_html_path,
                    public_path=output_dir / executed_html_path.name,
                ),
            )
        if executable_after is None:
            executable_after = capture_after(
                "browser_executable_after", lambda: artifact(executable_path)
            )
        if script_after is None:
            script_after = capture_after(
                "qa_script_after", lambda: artifact(Path(__file__).resolve())
            )
        if runtime_after is None:
            runtime_after = capture_after(
                "runtime_identity_after", lambda: runtime_identity(executable_path)
            )
        if formal_binding_after is None and html_after is not None:
            formal_binding_after = capture_after(
                "formal_release_binding_after",
                lambda: validate_formal_release_binding(
                    formal_paths,
                    html_identity=html_after,
                ),
            )
        if payload is None:
            payload = base_payload(
                args=args,
                html_before=html_before,
                executable_before=executable_before,
                script_before=script_before,
                git_head=git_head,
                profiles=profiles,
                browser_version=browser_version,
                runtime_before=runtime_before,
                runtime_after=runtime_after,
                executed_html=executed_html_identity,
                static_html_security=static_html_security,
                formal_release_binding=formal_binding_before,
            )
        payload.update(
            {
                "input_html_after": html_after,
                "executed_html_after": executed_html_after,
                "browser_executable_after": executable_after,
                "qa_script_after": script_after,
                "formal_release_binding_after": formal_binding_after,
                "runtime_identity_after": runtime_after,
                "runtime_identity_stable": runtime_after == runtime_before,
                "after_identity_capture_errors": after_identity_capture_errors,
                "elapsed_seconds": round(time.monotonic() - started, 3),
                "exception": {"type": type(error).__name__, "message": str(error)},
                "pass": False,
            }
        )
        failure_source = (
            staging
            if path_exists(staging)
            else output_dir
            if path_exists(output_dir)
            else None
        )
        if failure_source is not None:
            try:
                failure_dir = publish_failure(failure_source, output_dir, payload)
            except Exception as publication_error:
                raise SystemExit(
                    "QA failed and failure publication also failed; "
                    f"source preserved at {failure_source}: {error}; "
                    f"publication: {publication_error}"
                ) from error
            raise SystemExit(
                f"QA failed; evidence preserved at {failure_dir}: {error}"
            ) from error
        raise SystemExit(f"QA failed before staging creation: {error}") from error

if __name__ == "__main__":
    main()
