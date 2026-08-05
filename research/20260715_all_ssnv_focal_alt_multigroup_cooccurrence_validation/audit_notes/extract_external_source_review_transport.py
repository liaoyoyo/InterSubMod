#!/usr/bin/env python3
"""Extract one EOF-terminated JSON object from an external-review text transport."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
from typing import Any


class TransportExtractionError(RuntimeError):
    """Raised when a review transport is ambiguous or unsafe to normalize."""


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Build a JSON object while rejecting duplicate keys at every depth."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise TransportExtractionError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_nonfinite_constant(value: str) -> None:
    """Reject non-standard JSON constants such as NaN and Infinity."""
    raise TransportExtractionError(f"Non-finite JSON constant: {value}")


def extract_eof_json_object(raw_text: str) -> tuple[str, dict[str, Any], str]:
    """Return the sole JSON object when only non-braced prose precedes it."""
    object_start = raw_text.find("{")
    if object_start < 0:
        raise TransportExtractionError("Transport contains no JSON object")

    prefix = raw_text[:object_start]
    if "}" in prefix:
        raise TransportExtractionError("Transport prefix contains a closing brace")

    decoder = json.JSONDecoder(
        object_pairs_hook=reject_duplicate_keys,
        parse_constant=reject_nonfinite_constant,
    )
    candidate = raw_text[object_start:]
    try:
        payload, object_end = decoder.raw_decode(candidate)
    except (json.JSONDecodeError, TransportExtractionError) as exc:
        raise TransportExtractionError(f"Invalid transport JSON object: {exc}") from exc
    if candidate[object_end:].strip():
        raise TransportExtractionError(
            "Transport contains trailing non-whitespace data after the JSON object"
        )
    if not isinstance(payload, dict):
        raise TransportExtractionError("Transport JSON top level must be an object")
    return candidate[:object_end], payload, prefix


def read_bound_regular_file(path: Path) -> tuple[bytes, dict[str, Any]]:
    """Read one regular file and fail if its identity changes during the read."""
    resolved = path.expanduser().resolve(strict=True)
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(resolved, flags)
    try:
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            raise TransportExtractionError(f"Input is not a regular file: {resolved}")
        data = bytearray()
        offset = 0
        while offset < opened.st_size:
            chunk = os.pread(
                descriptor, min(1024 * 1024, opened.st_size - offset), offset
            )
            if not chunk:
                break
            data.extend(chunk)
            offset += len(chunk)
        closed_read = os.fstat(descriptor)
        live = resolved.stat()
        identity = (
            opened.st_dev,
            opened.st_ino,
            opened.st_size,
            opened.st_mtime_ns,
        )
        if len(data) != opened.st_size or identity != (
            closed_read.st_dev,
            closed_read.st_ino,
            closed_read.st_size,
            closed_read.st_mtime_ns,
        ):
            raise TransportExtractionError(f"Input changed while reading: {resolved}")
        if (live.st_dev, live.st_ino) != (opened.st_dev, opened.st_ino):
            raise TransportExtractionError(f"Input path was replaced: {resolved}")
        encoded = bytes(data)
        return encoded, {
            "path": str(resolved),
            "size_bytes": opened.st_size,
            "sha256": hashlib.sha256(encoded).hexdigest(),
            "mode": oct(opened.st_mode & 0o777),
        }
    finally:
        os.close(descriptor)


def write_new_readonly_payload(path: Path, payload_text: str) -> dict[str, Any]:
    """Create a new mode-0444 JSON payload without permitting overwrite."""
    parent = path.expanduser().parent.resolve(strict=True)
    output_path = parent / path.name
    encoded = (payload_text.rstrip() + "\n").encode("utf-8")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(output_path, flags, 0o444)
    try:
        offset = 0
        while offset < len(encoded):
            written = os.write(descriptor, encoded[offset:])
            if written <= 0:
                raise TransportExtractionError(
                    f"Unable to complete payload write: {output_path}"
                )
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        output_stat = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    parent_descriptor = os.open(parent, os.O_RDONLY | os.O_CLOEXEC)
    try:
        os.fsync(parent_descriptor)
    finally:
        os.close(parent_descriptor)
    return {
        "path": str(output_path),
        "size_bytes": output_stat.st_size,
        "sha256": hashlib.sha256(encoded).hexdigest(),
        "mode": oct(output_stat.st_mode & 0o777),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw_data, raw_record = read_bound_regular_file(args.raw_input)
    try:
        raw_text = raw_data.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise TransportExtractionError("Transport is not valid UTF-8") from exc
    payload_text, payload, prefix = extract_eof_json_object(raw_text)
    output_record = write_new_readonly_payload(args.output, payload_text)
    prefix_bytes = prefix.encode("utf-8")
    print(
        json.dumps(
            {
                "pass": True,
                "schema_name": payload.get("schema_name"),
                "reviewer_id": payload.get("reviewer_id"),
                "verdict": payload.get("verdict"),
                "prefix": {
                    "size_bytes": len(prefix_bytes),
                    "sha256": hashlib.sha256(prefix_bytes).hexdigest(),
                },
                "raw_input": raw_record,
                "extracted_payload": output_record,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
