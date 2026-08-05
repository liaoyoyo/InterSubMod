#!/usr/bin/env python3
"""Verify every stable-site primary artifact against the screen assignment contract."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
from collections import deque
from concurrent.futures import Future, ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence, TextIO

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import release_source_authority as SOURCE_AUTHORITY


SCHEMA_NAME = "intersubmod.stable_primary_artifact_audit"
SCHEMA_VERSION = "2.0.0"
ASSIGNMENT_SCHEMA = "intersubmod.all_ssnv_stable_multigroup_read_assignment"
SCREEN_CONTRACT = "phylo-v4.1_column_null95_modal_K10_RNULL40_min_group3"
WORKSPACE_ROOT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/"
    "20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation"
)
CANONICAL_PYTHON_CACHE_DIRNAME = (
    ".python_cache_m2v5_completion_v2_bound_bootstrap"
)
CANONICAL_PYTHON_CACHE_ROOT = WORKSPACE_ROOT / CANONICAL_PYTHON_CACHE_DIRNAME
ARTIFACT_RELATIVE_PATHS = {
    "reads": Path("reads") / "reads.tsv",
    "distance_matrix": Path("distance") / "BERNOULLI" / "matrix.csv",
    "methylation_matrix": Path("methylation") / "methylation.csv",
}
CANONICAL_PARAMETERS = {
    "workers": 40,
    "chunk_size": 64,
    "max_pending_chunks": 80,
    "progress_every": 1_000,
}


class AuditContractError(RuntimeError):
    """Raised when an assignment or primary artifact violates the frozen contract."""


SiteKey = tuple[str, str, int, str, str]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256(resolved),
    }


def source_paths() -> dict[str, Path]:
    return {
        "primary_artifact_auditor": Path(__file__).resolve(),
        "source_authority_validator": Path(SOURCE_AUTHORITY.__file__).resolve(),
    }


def capture_source_identity() -> dict[str, dict[str, Any]]:
    return {role: artifact(path) for role, path in source_paths().items()}


def capture_source_modes() -> dict[str, str]:
    return {
        role: oct(path.stat().st_mode & 0o777)
        for role, path in source_paths().items()
    }


def canonical_python_prefix() -> list[str]:
    return [
        sys.executable,
        "-I",
        "-X",
        f"pycache_prefix={CANONICAL_PYTHON_CACHE_ROOT}",
    ]


def canonical_command(
    *,
    site_results: Path,
    assignments: Path,
    consumer_receipts: Sequence[Path],
    output: Path,
    workers: int,
    chunk_size: int,
    max_pending_chunks: int,
    progress_every: int,
) -> list[str]:
    command = [
        *canonical_python_prefix(),
        str(Path(__file__).resolve()),
        "--site-results",
        str(site_results.resolve()),
        "--stable-assignments",
        str(assignments.resolve()),
    ]
    for receipt in consumer_receipts:
        command.extend(("--consumer-receipt", str(receipt.resolve())))
    command.extend(
        (
            "--output",
            str(output.resolve()),
            "--workers",
            str(workers),
            "--chunk-size",
            str(chunk_size),
            "--max-pending-chunks",
            str(max_pending_chunks),
            "--progress-every",
            str(progress_every),
        )
    )
    return command


def observed_process_command() -> list[str]:
    raw = Path("/proc/self/cmdline").read_bytes()
    if not raw.endswith(b"\0"):
        raise AuditContractError("Process command line is unavailable or malformed")
    return [os.fsdecode(token) for token in raw[:-1].split(b"\0")]


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def parse_bool(value: Any, *, field_name: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise AuditContractError(f"Invalid boolean for {field_name}: {value!r}")


def site_key(row: dict[str, Any], *, source: str) -> SiteKey:
    try:
        sample = str(row.get("sample") or row.get("dataset") or "").strip()
        chrom = str(row["chrom"]).strip()
        pos = int(row["pos"])
        ref = str(row["ref"]).strip().upper()
        alt = str(row["alt"]).strip().upper()
    except (KeyError, TypeError, ValueError) as error:
        raise AuditContractError(f"Malformed site identity in {source}") from error
    if not sample or not chrom or pos <= 0 or not ref or not alt:
        raise AuditContractError(f"Invalid site identity in {source}: {(sample, chrom, pos, ref, alt)}")
    return sample, chrom, pos, ref, alt


def load_stable_sites(path: Path) -> set[SiteKey]:
    stable: set[SiteKey] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "sample",
            "chrom",
            "pos",
            "ref",
            "alt",
            "stable_null_multigroup",
            "screen_contract",
        }
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise AuditContractError(f"Site table missing fields: {missing}")
        for row in reader:
            if not parse_bool(row["stable_null_multigroup"], field_name="stable_null_multigroup"):
                continue
            if row["screen_contract"] != SCREEN_CONTRACT:
                raise AuditContractError("Stable site screen contract drift")
            key = site_key(row, source="site_results")
            if key in stable:
                raise AuditContractError(f"Duplicate stable site: {key}")
            stable.add(key)
    return stable


def load_assignment_tasks(path: Path) -> dict[SiteKey, dict[str, Any]]:
    tasks: dict[SiteKey, dict[str, Any]] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                assignment = json.loads(line)
            except json.JSONDecodeError as error:
                raise AuditContractError(f"Malformed assignment JSON at line {line_number}") from error
            if not isinstance(assignment, dict):
                raise AuditContractError(f"Assignment line {line_number} is not an object")
            if (
                assignment.get("schema_name") != ASSIGNMENT_SCHEMA
                or assignment.get("schema_version") != "1.0.0"
                or assignment.get("screen_contract") != SCREEN_CONTRACT
                or assignment.get("artifact_identity_contract") != "sha256_size_path_v1"
                or assignment.get("strict_confirm_candidate") is not True
            ):
                raise AuditContractError(f"Assignment contract drift at line {line_number}")
            posthoc = assignment.get("posthoc")
            if not isinstance(posthoc, dict):
                raise AuditContractError(f"Assignment posthoc missing at line {line_number}")
            identity = {
                "sample": assignment.get("sample") or assignment.get("dataset"),
                "chrom": assignment.get("chrom"),
                "pos": assignment.get("pos"),
                "ref": posthoc.get("ref"),
                "alt": posthoc.get("alt"),
            }
            key = site_key(identity, source=f"assignment line {line_number}")
            if key in tasks:
                raise AuditContractError(f"Duplicate assignment site: {key}")
            region_raw = assignment.get("region_dir")
            if not isinstance(region_raw, str) or not region_raw:
                raise AuditContractError(f"Assignment region_dir missing at {key}")
            contracts = assignment.get("primary_artifacts")
            if not isinstance(contracts, dict) or set(contracts) != set(ARTIFACT_RELATIVE_PATHS):
                raise AuditContractError(f"Assignment primary_artifacts contract drift at {key}")
            tasks[key] = {
                "key": key,
                "region_dir": str(Path(region_raw).expanduser().resolve()),
                "contracts": contracts,
            }
    return tasks


def verify_task(task: dict[str, Any]) -> dict[str, Any]:
    region_dir = Path(task["region_dir"])
    verified: list[dict[str, Any]] = []
    for role, relative_path in ARTIFACT_RELATIVE_PATHS.items():
        record = task["contracts"].get(role)
        if not isinstance(record, dict):
            raise AuditContractError(f"Missing {role} identity at {task['key']}")
        missing = sorted({"path", "size_bytes", "sha256"}.difference(record))
        if missing:
            raise AuditContractError(f"{role} identity missing {missing} at {task['key']}")
        expected_path = (region_dir / relative_path).resolve()
        declared_path = Path(str(record["path"])).expanduser().resolve()
        if declared_path != expected_path:
            raise AuditContractError(f"{role} path mismatch at {task['key']}")
        if not expected_path.is_file():
            raise AuditContractError(f"Missing {role} artifact at {task['key']}: {expected_path}")
        before = expected_path.stat()
        observed_sha = sha256(expected_path)
        after = expected_path.stat()
        if before.st_size != after.st_size or before.st_mtime_ns != after.st_mtime_ns:
            raise AuditContractError(f"{role} changed during hashing at {task['key']}")
        if int(record["size_bytes"]) != after.st_size:
            raise AuditContractError(f"{role} size mismatch at {task['key']}")
        if str(record["sha256"]) != observed_sha:
            raise AuditContractError(f"{role} sha256 mismatch at {task['key']}")
        verified.append(
            {
                "role": role,
                "path": str(expected_path),
                "size_bytes": after.st_size,
                "mtime_ns": after.st_mtime_ns,
                "sha256": observed_sha,
            }
        )
    return {"key": list(task["key"]), "artifacts": verified}


def chunked(values: Iterable[dict[str, Any]], size: int) -> Iterator[list[dict[str, Any]]]:
    chunk: list[dict[str, Any]] = []
    for value in values:
        chunk.append(value)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def verify_chunk(chunk: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [verify_task(task) for task in chunk]


def bounded_results(
    executor: ProcessPoolExecutor,
    chunks: Iterable[list[dict[str, Any]]],
    max_pending: int,
) -> Iterator[dict[str, Any]]:
    iterator = iter(chunks)
    pending: deque[Future[list[dict[str, Any]]]] = deque()

    def fill() -> None:
        while len(pending) < max_pending:
            try:
                pending.append(executor.submit(verify_chunk, next(iterator)))
            except StopIteration:
                return

    fill()
    while pending:
        yield from pending.popleft().result()
        fill()


def artifact_set_digest(results: Iterable[dict[str, Any]]) -> str:
    lines: list[str] = []
    for result in results:
        key_text = "|".join(str(value) for value in result["key"])
        for record in result["artifacts"]:
            lines.append(
                "|".join(
                    (
                        key_text,
                        record["role"],
                        record["path"],
                        str(record["size_bytes"]),
                        record["sha256"],
                    )
                )
            )
    return hashlib.sha256("\n".join(sorted(lines)).encode()).hexdigest()


def run_audit(
    site_results: Path,
    assignments: Path,
    output: Path,
    *,
    consumer_receipts: Sequence[Path] = (),
    workers: int,
    chunk_size: int,
    max_pending_chunks: int,
    progress_every: int,
    command: list[str] | None = None,
) -> dict[str, Any]:
    if os.path.lexists(output):
        raise FileExistsError(f"Refusing to overwrite existing audit receipt: {output}")
    execution = {
        "workers": workers,
        "chunk_size": chunk_size,
        "max_pending_chunks": max_pending_chunks,
        "progress_every": progress_every,
    }
    expected_command = canonical_command(
        site_results=site_results,
        assignments=assignments,
        consumer_receipts=consumer_receipts,
        output=output,
        **execution,
    )
    formal_release = execution == CANONICAL_PARAMETERS
    if formal_release and (
        command != expected_command or observed_process_command() != expected_command
    ):
        raise AuditContractError(
            "Formal primary artifact audit is direct-CLI canonical-process only"
        )
    source_authority = (
        SOURCE_AUTHORITY.validate_release_source_authority(
            {"release_source_authority_validator", "primary_artifact_auditor"}
        )
        if formal_release
        else {"status": "NOT_APPLICABLE_NON_RELEASE_TEST_OR_DEVELOPMENT_RUN"}
    )
    source_identity_before = capture_source_identity()
    source_modes_before = capture_source_modes()
    if formal_release and set(source_modes_before.values()) != {"0o444"}:
        raise AuditContractError(
            f"Primary audit sources are not mode 0444: {source_modes_before}"
        )
    started_at = now_utc()
    input_artifacts_before = {
        "site_results": artifact(site_results),
        "stable_assignments": artifact(assignments),
        "consumer_receipts": [artifact(path) for path in consumer_receipts],
    }
    stable_keys = load_stable_sites(site_results)
    tasks_by_key = load_assignment_tasks(assignments)
    assignment_keys = set(tasks_by_key)
    if stable_keys != assignment_keys:
        raise AuditContractError(
            "Stable site/assignment key-set mismatch: "
            f"stable={len(stable_keys)} assignments={len(assignment_keys)} "
            f"missing={sorted(stable_keys - assignment_keys)[:3]} "
            f"extra={sorted(assignment_keys - stable_keys)[:3]}"
        )
    ordered_tasks = [tasks_by_key[key] for key in sorted(tasks_by_key)]
    results: list[dict[str, Any]] = []
    if ordered_tasks and workers == 1:
        for task in ordered_tasks:
            results.append(verify_task(task))
            if len(results) % progress_every == 0 or len(results) == len(ordered_tasks):
                print(f"verified={len(results)}/{len(ordered_tasks)}", flush=True)
    elif ordered_tasks:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            for result in bounded_results(
                executor,
                chunked(ordered_tasks, chunk_size),
                max_pending_chunks,
            ):
                results.append(result)
                if len(results) % progress_every == 0 or len(results) == len(ordered_tasks):
                    print(f"verified={len(results)}/{len(ordered_tasks)}", flush=True)
    verification_digest = artifact_set_digest(results)
    input_artifacts = {
        "site_results": artifact(site_results),
        "stable_assignments": artifact(assignments),
        "consumer_receipts": [artifact(path) for path in consumer_receipts],
    }
    if input_artifacts != input_artifacts_before:
        raise AuditContractError("Audit input artifacts changed during verification")
    source_identity_after_compute = capture_source_identity()
    source_modes_after_compute = capture_source_modes()
    if (
        source_identity_after_compute != source_identity_before
        or source_modes_after_compute != source_modes_before
    ):
        raise AuditContractError("Primary audit source identity changed during execution")
    finished_at = now_utc()
    recorded_command = command or canonical_command(
        site_results=site_results,
        assignments=assignments,
        consumer_receipts=consumer_receipts,
        output=output,
        workers=workers,
        chunk_size=chunk_size,
        max_pending_chunks=max_pending_chunks,
        progress_every=progress_every,
    )
    receipt = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": finished_at,
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "task_type": (
            "B_comprehensive_validation"
            if formal_release
            else "NON_RELEASE_TEST_OR_DEVELOPMENT_AUDIT"
        ),
        "scope": (
            "complete_primary_stable_null_multigroup_set"
            if formal_release
            else "explicit_non_release_input_set"
        ),
        "formal_task_b_release_eligible": formal_release,
        "command": recorded_command,
        "source_authority": source_authority,
        "code": source_identity_after_compute,
        "source_lock": {
            "source_identity_before": source_identity_before,
            "source_identity_after_compute": source_identity_after_compute,
            "source_modes_before": source_modes_before,
            "source_modes_after_compute": source_modes_after_compute,
            "all_sources_read_only_and_unchanged": True,
        },
        "inputs": input_artifacts,
        "verification": {
            "stable_site_assignment_key_sets_exact": True,
            "artifact_roles_exact": sorted(ARTIFACT_RELATIVE_PATHS),
            "path_size_sha256_verified": True,
            "artifact_set_sha256": verification_digest,
            "errors": 0,
        },
        "counts": {
            "stable_sites": len(stable_keys),
            "assignment_records": len(tasks_by_key),
            "primary_artifacts_expected": len(stable_keys) * len(ARTIFACT_RELATIVE_PATHS),
            "primary_artifacts_verified": sum(len(result["artifacts"]) for result in results),
        },
        "execution": execution,
        "pass_semantics": "execution_integrity_and_primary_artifact_identity_only",
        "scientific_interpretation": {
            "is_negative_result": False,
            "statement": "Artifact identity verification does not establish a biological claim.",
        },
        "pass": True,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(receipt, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    output.chmod(0o444)
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-results", type=Path, required=True)
    parser.add_argument("--stable-assignments", type=Path, required=True)
    parser.add_argument("--consumer-receipt", type=Path, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=max(1, min(16, (os.cpu_count() or 4) - 2)))
    parser.add_argument("--chunk-size", type=int, default=64)
    parser.add_argument("--max-pending-chunks", type=int, default=32)
    parser.add_argument("--progress-every", type=int, default=1000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if min(
        args.workers,
        args.chunk_size,
        args.max_pending_chunks,
        args.progress_every,
    ) < 1:
        raise ValueError("workers, chunk-size, max-pending-chunks, and progress-every must be positive")
    actual_command = [
        *canonical_python_prefix(),
        str(Path(__file__).resolve()),
        *sys.argv[1:],
    ]
    expected_command = canonical_command(
        site_results=args.site_results,
        assignments=args.stable_assignments,
        consumer_receipts=args.consumer_receipt,
        output=args.output,
        workers=args.workers,
        chunk_size=args.chunk_size,
        max_pending_chunks=args.max_pending_chunks,
        progress_every=args.progress_every,
    )
    if actual_command != expected_command:
        raise AuditContractError("Primary audit command is not canonical")
    receipt = run_audit(
        args.site_results,
        args.stable_assignments,
        args.output,
        consumer_receipts=args.consumer_receipt,
        workers=args.workers,
        chunk_size=args.chunk_size,
        max_pending_chunks=args.max_pending_chunks,
        progress_every=args.progress_every,
        command=actual_command,
    )
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "stable_sites": receipt["counts"]["stable_sites"],
                "primary_artifacts_verified": receipt["counts"]["primary_artifacts_verified"],
                "artifact_set_sha256": receipt["verification"]["artifact_set_sha256"],
                "pass": receipt["pass"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
