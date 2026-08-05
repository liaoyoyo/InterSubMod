#!/usr/bin/env python3
"""Capture or verify sampled identities for canonical LongPhase-S tagged BAMs."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import stat
import tempfile
from pathlib import Path
from typing import Any


EXPECTED_SAMPLES = (
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
)
CHUNK_SIZE = 1024 * 1024


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def strict_json_load(path: Path) -> Any:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key {key!r}: {path}")
            result[key] = value
        return result

    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=reject_duplicates)


def chunk_identity(path: Path, label: str, offset: int, length: int) -> dict[str, Any]:
    with path.open("rb") as handle:
        handle.seek(offset)
        payload = handle.read(length)
    if len(payload) != length:
        raise RuntimeError(f"short read for {path} at offset {offset}: {len(payload)} != {length}")
    return {
        "label": label,
        "offset": offset,
        "length": length,
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def bam_identity(requested: Path) -> dict[str, Any]:
    logical = requested.lstat()
    resolved = requested.resolve(strict=True)
    target = resolved.stat()
    if not stat.S_ISREG(target.st_mode):
        raise RuntimeError(f"BAM target is not a regular file: {requested} -> {resolved}")
    size = target.st_size
    length = min(CHUNK_SIZE, size)
    middle = max(0, (size - length) // 2)
    last = max(0, size - length)
    return {
        "requested_path": str(requested),
        "realpath": str(resolved),
        "logical_is_symlink": stat.S_ISLNK(logical.st_mode),
        "logical_size_bytes": logical.st_size,
        "logical_mtime_ns": logical.st_mtime_ns,
        "logical_ctime_ns": logical.st_ctime_ns,
        "st_dev": target.st_dev,
        "st_ino": target.st_ino,
        "size_bytes": size,
        "mtime_ns": target.st_mtime_ns,
        "ctime_ns": target.st_ctime_ns,
        "chunk_size_bytes": CHUNK_SIZE,
        "chunks": [
            chunk_identity(resolved, "first", 0, length),
            chunk_identity(resolved, "middle", middle, length),
            chunk_identity(resolved, "last", last, length),
        ],
    }


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def atomic_write_json(path: Path, value: Any) -> None:
    if path.exists() or path.is_symlink():
        raise RuntimeError(f"refusing to overwrite audit receipt: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    temporary_path = Path(temporary)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2).encode("utf-8"))
            handle.write(b"\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except BaseException:
        if temporary_path.exists():
            temporary_path.rename(temporary_path.with_suffix(temporary_path.suffix + ".failed"))
        raise


def manifest_bams(manifest_path: Path) -> list[tuple[str, Path]]:
    manifest = strict_json_load(manifest_path)
    samples = manifest.get("samples")
    if not isinstance(samples, list):
        raise RuntimeError("manifest samples must be an array")
    observed = [item.get("sample") for item in samples if isinstance(item, dict)]
    if tuple(observed) != EXPECTED_SAMPLES:
        raise RuntimeError(f"expected ordered canonical seven datasets; observed {observed}")
    result: list[tuple[str, Path]] = []
    for item in samples:
        raw = item.get("tumor_bam")
        if not isinstance(raw, str) or not Path(raw).is_absolute():
            raise RuntimeError(f"{item['sample']} tumor_bam must be an absolute path")
        result.append((item["sample"], Path(raw)))
    if len({str(path) for _, path in result}) != len(EXPECTED_SAMPLES):
        raise RuntimeError("canonical tagged BAM paths are not unique")
    return result


def capture(manifest_path: Path) -> dict[str, Any]:
    samples = [
        {"sample": sample, "identity": bam_identity(path)}
        for sample, path in manifest_bams(manifest_path)
    ]
    payload = {
        "schema_name": "intersubmod.canonical_longphase_bam_immutability",
        "schema_version": "1.0.0",
        "mode": "baseline",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "manifest_path": str(manifest_path.resolve(strict=True)),
        "manifest_sha256": sha256_file(manifest_path),
        "dataset_count": len(samples),
        "samples": samples,
    }
    payload["identity_set_sha256"] = hashlib.sha256(canonical_bytes(samples)).hexdigest()
    return payload


def verify(manifest_path: Path, baseline_path: Path) -> tuple[dict[str, Any], bool]:
    baseline = strict_json_load(baseline_path)
    if (
        baseline.get("schema_name") != "intersubmod.canonical_longphase_bam_immutability"
        or baseline.get("schema_version") != "1.0.0"
        or baseline.get("mode") != "baseline"
        or baseline.get("dataset_count") != 7
    ):
        raise RuntimeError("baseline receipt contract mismatch")
    observed = capture(manifest_path)
    expected_by_sample = {item["sample"]: item["identity"] for item in baseline["samples"]}
    comparisons = []
    for item in observed["samples"]:
        sample = item["sample"]
        expected = expected_by_sample.get(sample)
        differences = []
        if expected is None:
            differences.append("sample_missing_from_baseline")
        else:
            differences = sorted(key for key in set(expected) | set(item["identity"]) if expected.get(key) != item["identity"].get(key))
        comparisons.append({
            "sample": sample,
            "match": not differences,
            "differing_fields": differences,
            "expected_identity": expected,
            "observed_identity": item["identity"],
        })
    all_match = len(expected_by_sample) == 7 and all(item["match"] for item in comparisons)
    result = {
        "schema_name": "intersubmod.canonical_longphase_bam_immutability_verification",
        "schema_version": "1.0.0",
        "mode": "verification",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "baseline_path": str(baseline_path.resolve(strict=True)),
        "baseline_sha256": sha256_file(baseline_path),
        "manifest_path": observed["manifest_path"],
        "manifest_sha256": observed["manifest_sha256"],
        "dataset_count": len(comparisons),
        "all_match": all_match,
        "comparisons": comparisons,
    }
    return result, all_match


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--baseline", type=Path)
    args = parser.parse_args()
    try:
        if args.baseline is None:
            result = capture(args.manifest)
            all_match = True
        else:
            result, all_match = verify(args.manifest, args.baseline)
        atomic_write_json(args.output, result)
        label = "BASELINE" if args.baseline is None else "VERIFY"
        print(f"CANONICAL LONGPHASE BAM {label}: {'PASS' if all_match else 'FAIL'} -> {args.output}")
        return 0 if all_match else 7
    except BaseException as exc:
        print(f"CANONICAL LONGPHASE BAM AUDIT ERROR: {exc}", file=os.sys.stderr)
        return 7


if __name__ == "__main__":
    raise SystemExit(main())
