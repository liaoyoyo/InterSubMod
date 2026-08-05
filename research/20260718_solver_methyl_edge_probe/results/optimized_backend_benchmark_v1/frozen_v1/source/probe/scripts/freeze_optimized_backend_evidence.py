#!/usr/bin/env python3
"""Create and verify a no-clobber frozen evidence bundle for the bounded probe."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import platform
import stat
import subprocess
import sys
import time
from typing import Any, Iterable


REPO = pathlib.Path("/big7_disk/liaoyoyo2001/InterSubMod")
PROBE_ROOT = REPO / "research/20260718_solver_methyl_edge_probe"
BENCHMARK_ROOT = PROBE_ROOT / "results/optimized_backend_benchmark_v1"
CURRENT_ROOT = (
    REPO / "research/20260716_read_linked_hypercube_exact_likelihood_validation"
)


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_read_only(command: list[str]) -> str:
    result = subprocess.run(
        command,
        cwd=REPO,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout.strip()


def iter_source_files(root: pathlib.Path, suffixes: set[str]) -> Iterable[pathlib.Path]:
    for path in sorted(root.rglob("*")):
        if (
            path.is_file()
            and "__pycache__" not in path.parts
            and path.suffix in suffixes
        ):
            yield path


class Freezer:
    """Copy bytes once, record their identity, and reject destination collisions."""

    def __init__(self, output_dir: pathlib.Path) -> None:
        self.output_dir = output_dir
        self.entries: list[dict[str, Any]] = []
        self.destinations: set[pathlib.PurePosixPath] = set()

    def copy(self, source: pathlib.Path, destination: pathlib.PurePosixPath) -> None:
        source = source.resolve(strict=True)
        if not source.is_file():
            raise ValueError(f"source is not a regular file: {source}")
        if destination.is_absolute() or ".." in destination.parts:
            raise ValueError(f"unsafe destination: {destination}")
        if destination in self.destinations:
            raise FileExistsError(f"duplicate bundle destination: {destination}")
        self.destinations.add(destination)

        target = self.output_dir.joinpath(*destination.parts)
        target.parent.mkdir(mode=0o755, parents=True, exist_ok=True)
        descriptor = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
        try:
            with source.open("rb") as input_handle, os.fdopen(
                descriptor, "wb", closefd=True
            ) as output_handle:
                while True:
                    chunk = input_handle.read(1024 * 1024)
                    if not chunk:
                        break
                    output_handle.write(chunk)
                output_handle.flush()
                os.fsync(output_handle.fileno())
        except Exception:
            # The top-level directory is new and a failed freeze remains auditable.
            raise

        source_hash = sha256_file(source)
        target_hash = sha256_file(target)
        if target_hash != source_hash:
            raise RuntimeError(f"copy hash mismatch: {source} -> {target}")
        self.entries.append(
            {
                "role": destination.parts[0],
                "source_path": str(source),
                "bundle_path": destination.as_posix(),
                "size_bytes": source.stat().st_size,
                "sha256": source_hash,
                "source_mode": stat.S_IMODE(source.stat().st_mode),
                "frozen_mode": 0o444,
            }
        )

    def copy_tree(
        self,
        source_root: pathlib.Path,
        destination_root: pathlib.PurePosixPath,
        suffixes: set[str] | None = None,
    ) -> None:
        for source in sorted(source_root.rglob("*")):
            if not source.is_file() or "__pycache__" in source.parts:
                continue
            if suffixes is not None and source.suffix not in suffixes:
                continue
            self.copy(source, destination_root / source.relative_to(source_root))


def write_exclusive(path: pathlib.Path, payload: bytes, mode: int = 0o444) -> None:
    path.parent.mkdir(mode=0o755, parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, mode)
    with os.fdopen(descriptor, "wb", closefd=True) as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def verify_bundle(output_dir: pathlib.Path, manifest: dict[str, Any]) -> dict[str, Any]:
    mismatches: list[str] = []
    for entry in manifest["entries"]:
        path = output_dir / entry["bundle_path"]
        if not path.is_file() or path.is_symlink():
            mismatches.append(f"not_regular:{entry['bundle_path']}")
            continue
        if path.stat().st_size != entry["size_bytes"]:
            mismatches.append(f"size:{entry['bundle_path']}")
        if sha256_file(path) != entry["sha256"]:
            mismatches.append(f"sha256:{entry['bundle_path']}")
        if stat.S_IMODE(path.stat().st_mode) != 0o444:
            mismatches.append(f"mode:{entry['bundle_path']}")

    manifest_path = output_dir / "manifest.json"
    manifest_sidecar = output_dir / "manifest.json.sha256"
    expected_sidecar = f"{sha256_file(manifest_path)}  manifest.json\n"
    if manifest_sidecar.read_text(encoding="utf-8") != expected_sidecar:
        mismatches.append("manifest_sidecar")

    sums_path = output_dir / "SHA256SUMS"
    expected_sums = []
    for path in sorted(output_dir.rglob("*")):
        if (
            path.is_file()
            and not path.is_symlink()
            and path.name not in {"SHA256SUMS", "SHA256SUMS.sha256"}
        ):
            expected_sums.append(
                f"{sha256_file(path)}  {path.relative_to(output_dir).as_posix()}"
            )
    if sums_path.read_text(encoding="utf-8") != "\n".join(expected_sums) + "\n":
        mismatches.append("SHA256SUMS")
    sums_sidecar = output_dir / "SHA256SUMS.sha256"
    expected_sums_sidecar = f"{sha256_file(sums_path)}  SHA256SUMS\n"
    if sums_sidecar.read_text(encoding="utf-8") != expected_sums_sidecar:
        mismatches.append("SHA256SUMS_sidecar")

    return {
        "all_pass": not mismatches,
        "mismatches": mismatches,
        "verified_entries": len(manifest["entries"]),
        "manifest_sha256": sha256_file(manifest_path),
        "sha256sums_sha256": sha256_file(sums_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report",
        type=pathlib.Path,
        default=PROBE_ROOT
        / "20260718_HypercubeOptimizedBackend完整流程與效能驗證_01.html",
    )
    parser.add_argument("--browser-qa", type=pathlib.Path, required=True)
    parser.add_argument("--screenshot", type=pathlib.Path, required=True)
    parser.add_argument(
        "--current-baseline-dir",
        type=pathlib.Path,
        default=pathlib.Path("/tmp/intersubmod_solver_baseline_20260718"),
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=BENCHMARK_ROOT / "frozen_v1",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(
            f"no-clobber freeze refused because output exists: {output_dir}"
        )
    output_dir.mkdir(mode=0o755, parents=False, exist_ok=False)
    freezer = Freezer(output_dir)

    # Freeze the bounded implementation and every directly exercised probe source.
    freezer.copy_tree(
        PROBE_ROOT / "scripts",
        pathlib.PurePosixPath("source/probe/scripts"),
        {".py", ".json"},
    )
    freezer.copy_tree(
        PROBE_ROOT / "tests",
        pathlib.PurePosixPath("source/probe/tests"),
        {".py", ".json"},
    )
    freezer.copy(
        CURRENT_ROOT / "scripts/hypercube_exact.py",
        pathlib.PurePosixPath("source/current/scripts/hypercube_exact.py"),
    )
    freezer.copy(
        CURRENT_ROOT / "tests/test_hypercube_exact.py",
        pathlib.PurePosixPath("source/current/tests/test_hypercube_exact.py"),
    )

    # Freeze decisions and the mutable plan as observed at handoff time.
    for name in (
        "20260718_Hypercube邊與subcube改良研究計畫_01.md",
        "pre-decision-audit.md",
        "implementation-notes.md",
    ):
        freezer.copy(
            PROBE_ROOT / name,
            pathlib.PurePosixPath(f"context/{name}"),
        )
    state_audit = (
        REPO
        / "state/cycles/cycle_20260718-solver-methyl-edge-probe/audit.json"
    )
    if state_audit.exists():
        freezer.copy(state_audit, pathlib.PurePosixPath("context/cycle_audit.json"))

    # Freeze receipts and both sides of the measured raw-run evidence.
    freezer.copy(
        BENCHMARK_ROOT / "benchmark_receipt.json",
        pathlib.PurePosixPath("evidence/benchmark_receipt.json"),
    )
    freezer.copy(
        BENCHMARK_ROOT / "oracle_receipt.json",
        pathlib.PurePosixPath("evidence/oracle_receipt.json"),
    )
    freezer.copy_tree(
        BENCHMARK_ROOT / "runs",
        pathlib.PurePosixPath("evidence/optimized_raw_runs"),
    )
    freezer.copy_tree(
        args.current_baseline_dir.resolve(strict=True),
        pathlib.PurePosixPath("evidence/current_baseline_raw_runs"),
    )

    # Freeze reader-facing and visual QA artifacts.
    freezer.copy(
        args.report.resolve(strict=True),
        pathlib.PurePosixPath(f"report/{args.report.name}"),
    )
    freezer.copy(
        args.browser_qa.resolve(strict=True),
        pathlib.PurePosixPath(f"qa/{args.browser_qa.name}"),
    )
    freezer.copy(
        args.screenshot.resolve(strict=True),
        pathlib.PurePosixPath(f"qa/{args.screenshot.name}"),
    )

    scoped_status = run_read_only(
        [
            "git",
            "status",
            "--porcelain=v1",
            "--",
            str(PROBE_ROOT.relative_to(REPO)),
            str((CURRENT_ROOT / "scripts/hypercube_exact.py").relative_to(REPO)),
        ]
    )
    manifest = {
        "schema": "intersubmod.optimized_backend_frozen_evidence.v1",
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "task_type": "A_EXPLORATORY_IMPLEMENTATION_WITH_EXACTNESS_GATES",
        "scope": {
            "tier": "BOUNDED_DUAL_PILOT",
            "partial": True,
            "production_claim_allowed": False,
            "canonical_promotion_allowed": False,
            "claim_ceiling": (
                "Exact structural equivalence and local performance for three "
                "declared complete cases; finite-oracle and fail-closed behavior."
            ),
        },
        "provenance": {
            "repo_root": str(REPO),
            "git_head": run_read_only(["git", "rev-parse", "HEAD"]),
            "scoped_git_status_porcelain": scoped_status,
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "thread_environment": {
                key: os.environ.get(key)
                for key in (
                    "PYTHONHASHSEED",
                    "OPENBLAS_NUM_THREADS",
                    "OMP_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
            },
            "freeze_command": " ".join(sys.argv),
        },
        "promotion_gate": {
            "overall": "PASS_FOR_BOUNDED_DUAL_PILOT",
            "authorized_next_action": (
                "Run a separately frozen 33-tail plus H2009/H1437 stress panel."
            ),
            "canonical_or_production_promotion_allowed": False,
        },
        "entries": sorted(freezer.entries, key=lambda row: row["bundle_path"]),
    }
    manifest_payload = (
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    manifest_path = output_dir / "manifest.json"
    write_exclusive(manifest_path, manifest_payload)
    write_exclusive(
        output_dir / "manifest.json.sha256",
        f"{sha256_file(manifest_path)}  manifest.json\n".encode("utf-8"),
    )

    checksum_lines = []
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and not path.is_symlink():
            checksum_lines.append(
                f"{sha256_file(path)}  {path.relative_to(output_dir).as_posix()}"
            )
    sums_path = output_dir / "SHA256SUMS"
    write_exclusive(sums_path, ("\n".join(checksum_lines) + "\n").encode("utf-8"))
    write_exclusive(
        output_dir / "SHA256SUMS.sha256",
        f"{sha256_file(sums_path)}  SHA256SUMS\n".encode("utf-8"),
    )

    for path in sorted(output_dir.rglob("*"), reverse=True):
        os.chmod(path, 0o444 if path.is_file() else 0o555)
    os.chmod(output_dir, 0o555)

    verification = verify_bundle(output_dir, manifest)
    print(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "entries": len(manifest["entries"]),
                "verification": verification,
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if verification["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
