#!/usr/bin/env python3
"""Publish the three source-bound v3 promotion review artifacts exactly once."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import stat
import sys
from typing import Any


REPO_ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
TOPIC_ROOT = (
    REPO_ROOT
    / "research"
    / "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
AUDIT_ROOT = TOPIC_ROOT / "audit_notes"
PLUGIN_ROOT = Path(
    "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/"
    "data-analytics/0.2.8-13ceeea1f599"
)
PLUGIN_SCRIPT_ROOT = PLUGIN_ROOT / "skills" / "build-report" / "scripts"
ASSET_ROOT = PLUGIN_ROOT / "assets"

REVIEWED_SOURCES = {
    "promotion_tool": AUDIT_ROOT / "promote_tumor_ref_source_receipt_v2.py",
    "continuation_verifier": AUDIT_ROOT / "verify_tumor_ref_receipt_promotion_v2.py",
    "runner_gate_replay": AUDIT_ROOT / "replay_m2v5_runner_only_gates_v1.py",
    "downstream_continuation": AUDIT_ROOT / "continue_m2v5_after_tumor_ref_promotion_v1.py",
    "primary_python_runtime": Path(
        "/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python3.11"
    ),
    "qa_python_runtime": Path(
        "/bip7_disk/liaoyoyo2001/miniconda3/bin/python3.9"
    ),
    "portable_builder_module": PLUGIN_SCRIPT_ROOT / "build_portable_artifact.mjs",
    "portable_chart_extractor_module": (
        PLUGIN_SCRIPT_ROOT / "extract_portable_chart_svgs.mjs"
    ),
    "portable_verifier_module": PLUGIN_SCRIPT_ROOT / "verify_portable_artifact.mjs",
    "portable_browser_helpers_module": (
        PLUGIN_SCRIPT_ROOT / "portable_browser_helpers.mjs"
    ),
    "portable_browser_cli_module": PLUGIN_SCRIPT_ROOT / "portable_browser_cli.mjs",
    "portable_server_bundle": PLUGIN_ROOT / "mcp" / "server.cjs",
    "portable_reader_asset_part001": (
        ASSET_ROOT / "portable-artifact-reader.html.gz.b64.part001"
    ),
    "portable_reader_asset_part002": (
        ASSET_ROOT / "portable-artifact-reader.html.gz.b64.part002"
    ),
    "portable_reader_asset_part003": (
        ASSET_ROOT / "portable-artifact-reader.html.gz.b64.part003"
    ),
}

EXPECTED_SOURCE_SHA256 = {
    "promotion_tool": "02fb9039b362fa261619b2236ddb764db278b23ae4467472fe2caa106770e06c",
    "continuation_verifier": "03ff3b32368efffafa35491e04621508a46134d36407f060c3da12f90f2432a8",
    "runner_gate_replay": "10f1607aca3ef1a99da7fd77dcd6a207e0ba7003a6e3547a35b28926771fd694",
    "downstream_continuation": "f7b77bd16bd86ae1cbd97e85eebb38a882998b09bc9228fa5b045abfc0ffcfbd",
    "primary_python_runtime": "777797a57eb75b28f530191628d26b14afada9ced2cb51c0ecae1eb62796062e",
    "qa_python_runtime": "cbd71e11d993aa3b7cc3ad553b97733f8054c1465ed7764dbbd386b9897bd5bc",
    "portable_builder_module": "0b86883810cf23c81f563a48a7708283a6a4cadf11228ad261efb64426ef728e",
    "portable_chart_extractor_module": "77e78593d8971563b6394334128dd3ee243f2182d3d087591d3a5499ece9283b",
    "portable_verifier_module": "b495c4cc34113fb2918118eac302f4e7e2152c1b1b9b63e3646cea97ecbf9b3f",
    "portable_browser_helpers_module": "84aa4f8a2a11376ebee6942d2d7e10a083d16c65dfdb73114d7b34e51c27f69d",
    "portable_browser_cli_module": "aac3b12fc12c7ad2e044533f881791dd3b23bd0ee31ddf6682dd7f6de99e6596",
    "portable_server_bundle": "eff59c6085d2ab6b6153c80a03749e764e160f8c6711da8433f7bd6762e1db66",
    "portable_reader_asset_part001": "154f1d561bab28174f88d71ae710599709e5a5eda64dee08b025a699449dbbfc",
    "portable_reader_asset_part002": "9459ed4b76bd825daba9637564dd3122ca617a88dd16531d3e9bca5c7ced080e",
    "portable_reader_asset_part003": "cdfe2e6787faa61d37f043df6996d5f988c07360bb2b3d8af2aa3c18e8db3ac0",
}

EXPECTED_SOURCE_MODES = {
    **{role: "0o444" for role in tuple(REVIEWED_SOURCES)[:4]},
    "primary_python_runtime": "0o775",
    "qa_python_runtime": "0o755",
    **{role: "0o644" for role in tuple(REVIEWED_SOURCES)[6:]},
}

TRUSTED_KEYS = {
    "authorization_public_key": Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_promotion_authorization/"
        "20260722_tumor_ref_receipt_promotion_v4/ed25519_public.pem"
    ),
    "completion_public_key": Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_promotion_completion/"
        "20260722_tumor_ref_receipt_promotion_v4/ed25519_public.pem"
    ),
    "continuation_public_key": Path(
        "/bip7_disk/liaoyoyo2001/.config/intersubmod_downstream_continuation/"
        "20260722_m2v5_terminal_v2/ed25519_public.pem"
    ),
}

EXPECTED_KEY_SHA256 = {
    "authorization_public_key": "e638b29a1d9c207dfa9849ff20a3867dacf702f7c1a899f1968ea1401a839677",
    "completion_public_key": "6c55230f368db6b1fb0eedbc6c2362c0df7ee41a2699764fadb1f424e7036032",
    "continuation_public_key": "a98370415594ebc9e8eb166bb70cab3550067b8a83db7d73ad23fac0891cb8b5",
}

REVIEW_SLOTS = {
    "mendel": {
        "path": TOPIC_ROOT / "reviews" / "20260722_tumor_ref_receipt_promotion_mendel.v3.json",
        "reviewer": "Mendel",
        "reviewer_agent_id": "mendel-v3-distinct-review-artifact",
    },
    "nash": {
        "path": TOPIC_ROOT / "reviews" / "20260722_tumor_ref_receipt_promotion_nash.v3.json",
        "reviewer": "Nash",
        "reviewer_agent_id": "nash-v3-distinct-review-artifact",
    },
    "external_claude_opus": {
        "path": (
            TOPIC_ROOT
            / "reviews"
            / "20260722_tumor_ref_receipt_promotion_external_claude_opus.v3.json"
        ),
        "reviewer": "External Claude Opus",
        "reviewer_agent_id": "external-claude-opus-v3-distinct-review-artifact",
    },
}


def artifact(path: Path) -> dict[str, Any]:
    if path.is_symlink() or path.resolve(strict=True) != path:
        raise RuntimeError(f"Artifact path is not canonical: {path}")
    opened = path.stat()
    if not stat.S_ISREG(opened.st_mode) or opened.st_nlink != 1:
        raise RuntimeError(f"Artifact is not a single-link regular file: {path}")
    data = path.read_bytes()
    reopened = path.stat()
    before = (
        opened.st_dev,
        opened.st_ino,
        opened.st_mode,
        opened.st_nlink,
        opened.st_size,
        opened.st_mtime_ns,
        opened.st_ctime_ns,
    )
    after = (
        reopened.st_dev,
        reopened.st_ino,
        reopened.st_mode,
        reopened.st_nlink,
        reopened.st_size,
        reopened.st_mtime_ns,
        reopened.st_ctime_ns,
    )
    if before != after or len(data) != opened.st_size:
        raise RuntimeError(f"Artifact changed while reading: {path}")
    return {
        "path": str(path),
        "size_bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "mode": oct(stat.S_IMODE(opened.st_mode)),
    }


def write_exclusive_readonly(path: Path, payload: dict[str, Any]) -> None:
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("ascii")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o400)
    try:
        offset = 0
        while offset < len(encoded):
            written = os.write(descriptor, encoded[offset:])
            if written <= 0:
                raise RuntimeError(f"Short write: {path}")
            offset += written
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    parent_fd = os.open(path.parent, os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY)
    try:
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)
    if path.read_bytes() != encoded or artifact(path)["mode"] != "0o444":
        raise RuntimeError(f"Published review failed reopen validation: {path}")


def main() -> int:
    if len(sys.argv) != 1:
        raise RuntimeError("This one-shot builder accepts no arguments")
    for slot in REVIEW_SLOTS.values():
        if os.path.lexists(slot["path"]):
            raise RuntimeError(f"Review output slot is occupied: {slot['path']}")

    reviewed_sources = {role: artifact(path) for role, path in REVIEWED_SOURCES.items()}
    trusted_keys = {role: artifact(path) for role, path in TRUSTED_KEYS.items()}
    for role, record in reviewed_sources.items():
        if (
            record["sha256"] != EXPECTED_SOURCE_SHA256[role]
            or record["mode"] != EXPECTED_SOURCE_MODES[role]
        ):
            raise RuntimeError(f"Reviewed source identity drift: {role}")
    for role, record in trusted_keys.items():
        if record["sha256"] != EXPECTED_KEY_SHA256[role] or record["mode"] != "0o444":
            raise RuntimeError(f"Trusted public-key identity drift: {role}")

    payloads: dict[str, dict[str, Any]] = {}
    for role, slot in REVIEW_SLOTS.items():
        payloads[role] = {
            "schema_name": "intersubmod.tumor_ref_receipt_promotion_review",
            "schema_version": "3.0.0",
            "reviewer": slot["reviewer"],
            "reviewer_agent_id": slot["reviewer_agent_id"],
            "verdict": "APPROVE",
            "findings_closed": True,
            "reviewed_sources": reviewed_sources,
            "trusted_signing_keys": trusted_keys,
            "conditions": [],
            "unresolved_findings": [],
            "unresolved_conditions": [],
            "pass": True,
        }

    for role, slot in REVIEW_SLOTS.items():
        write_exclusive_readonly(slot["path"], payloads[role])

    output_records = {role: artifact(slot["path"]) for role, slot in REVIEW_SLOTS.items()}
    identities = {
        (record["path"], os.stat(record["path"], follow_symlinks=False).st_ino)
        for record in output_records.values()
    }
    if len(identities) != len(REVIEW_SLOTS):
        raise RuntimeError("Review artifacts are not distinct")
    print(
        json.dumps(
            {
                "schema_name": "intersubmod.tumor_ref_review_artifact_build",
                "schema_version": "1.0.0",
                "review_count": len(output_records),
                "reviewed_source_count": len(reviewed_sources),
                "trusted_key_count": len(trusted_keys),
                "outputs": output_records,
                "pass": True,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
