#!/usr/bin/env python3
"""Extract one schema-validated review object from a bound Claude CLI envelope."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Any


class EnvelopeExtractionError(RuntimeError):
    """Raised when a Claude CLI envelope cannot be bound to its review command."""


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise EnvelopeExtractionError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_nonfinite_constant(value: str) -> None:
    raise EnvelopeExtractionError(f"Non-finite JSON constant: {value}")


def parse_json_object(data: bytes, *, label: str) -> dict[str, Any]:
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise EnvelopeExtractionError(f"{label} is not UTF-8") from exc
    decoder = json.JSONDecoder(
        object_pairs_hook=reject_duplicate_keys,
        parse_constant=reject_nonfinite_constant,
    )
    stripped = text.lstrip()
    try:
        payload, end = decoder.raw_decode(stripped)
    except (json.JSONDecodeError, EnvelopeExtractionError) as exc:
        raise EnvelopeExtractionError(f"{label} is invalid JSON: {exc}") from exc
    if stripped[end:].strip():
        raise EnvelopeExtractionError(f"{label} contains trailing data")
    if not isinstance(payload, dict):
        raise EnvelopeExtractionError(f"{label} top level must be an object")
    return payload


def read_bound_regular(path: Path) -> tuple[bytes, dict[str, Any]]:
    resolved = path.expanduser().resolve(strict=True)
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(resolved, flags)
    try:
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            raise EnvelopeExtractionError(f"Input is not regular: {resolved}")
        chunks: list[bytes] = []
        offset = 0
        while offset < opened.st_size:
            chunk = os.pread(descriptor, min(1024 * 1024, opened.st_size - offset), offset)
            if not chunk:
                break
            chunks.append(chunk)
            offset += len(chunk)
        data = b"".join(chunks)
        closed_read = os.fstat(descriptor)
        live = resolved.stat()
        identity = (
            opened.st_dev,
            opened.st_ino,
            opened.st_mode,
            opened.st_size,
            opened.st_mtime_ns,
            opened.st_ctime_ns,
        )
        if (
            len(data) != opened.st_size
            or identity
            != (
                closed_read.st_dev,
                closed_read.st_ino,
                closed_read.st_mode,
                closed_read.st_size,
                closed_read.st_mtime_ns,
                closed_read.st_ctime_ns,
            )
            or identity
            != (
                live.st_dev,
                live.st_ino,
                live.st_mode,
                live.st_size,
                live.st_mtime_ns,
                live.st_ctime_ns,
            )
        ):
            raise EnvelopeExtractionError(f"Input binding changed: {resolved}")
        return data, {
            "path": str(resolved),
            "size_bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "mode": oct(stat.S_IMODE(opened.st_mode)),
        }
    finally:
        os.close(descriptor)


def extract_structured_review(
    envelope: dict[str, Any],
    command_record: dict[str, Any],
    *,
    reviewer_token: str,
) -> dict[str, Any]:
    if (
        envelope.get("type") != "result"
        or envelope.get("subtype") != "success"
        or envelope.get("is_error") is not False
    ):
        raise EnvelopeExtractionError("Claude CLI envelope is not a clean success")
    expected_session = command_record.get("session_id")
    if (
        not isinstance(expected_session, str)
        or not expected_session
        or envelope.get("session_id") != expected_session
    ):
        raise EnvelopeExtractionError("Envelope session does not match command record")
    if command_record.get("reviewer") != reviewer_token:
        raise EnvelopeExtractionError("Command record reviewer token is wrong")
    structured = envelope.get("structured_output")
    if not isinstance(structured, dict):
        raise EnvelopeExtractionError("Envelope structured_output must be an object")
    reviewer_label = structured.get("reviewer_label")
    if (
        not isinstance(reviewer_label, str)
        or f"Reviewer {reviewer_token}" not in reviewer_label
    ):
        raise EnvelopeExtractionError("Structured output reviewer identity is wrong")
    return structured


def write_new_readonly_json(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    parent = path.expanduser().parent.resolve(strict=True)
    output = parent / path.name
    encoded = (
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(output, flags, 0o444)
    try:
        offset = 0
        while offset < len(encoded):
            written = os.write(descriptor, encoded[offset:])
            if written <= 0:
                raise EnvelopeExtractionError(f"Unable to write output: {output}")
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        written_stat = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    parent_fd = os.open(
        parent,
        os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_DIRECTORY", 0),
    )
    try:
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)
    return {
        "path": str(output),
        "size_bytes": written_stat.st_size,
        "sha256": hashlib.sha256(encoded).hexdigest(),
        "mode": oct(stat.S_IMODE(written_stat.st_mode)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--envelope", type=Path, required=True)
    parser.add_argument("--command-record", type=Path, required=True)
    parser.add_argument("--reviewer-token", choices=("A", "B"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    envelope_data, envelope_record = read_bound_regular(args.envelope)
    command_data, command_record_identity = read_bound_regular(args.command_record)
    if envelope_record["mode"] != "0o444" or command_record_identity["mode"] != "0o444":
        raise EnvelopeExtractionError("Envelope and command record must be mode 0444")
    envelope = parse_json_object(envelope_data, label="Claude CLI envelope")
    command_record = parse_json_object(command_data, label="Claude CLI command record")
    structured = extract_structured_review(
        envelope,
        command_record,
        reviewer_token=args.reviewer_token,
    )
    output_record = write_new_readonly_json(args.output, structured)
    print(
        json.dumps(
            {
                "pass": True,
                "reviewer_token": args.reviewer_token,
                "session_id": envelope["session_id"],
                "verdict": structured.get("verdict"),
                "reviewer_id": structured.get("reviewer_id"),
                "permission_denials": envelope.get("permission_denials", []),
                "command_record": command_record_identity,
                "envelope": envelope_record,
                "extracted_payload": output_record,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
