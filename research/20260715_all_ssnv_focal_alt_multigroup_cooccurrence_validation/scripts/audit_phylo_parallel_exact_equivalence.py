#!/usr/bin/env python3
"""Create a machine receipt for serial versus seed-parallel screen equivalence."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


SCHEMA_NAME = "intersubmod.phylo_parallel_exact_equivalence"
SCHEMA_VERSION = "1.0.0"


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def load_analyzer(path: Path):
    resolved = path.resolve(strict=True)
    module_name = f"_intersubmod_equivalence_analyzer_{sha256(resolved)[:16]}"
    if str(resolved.parent) not in sys.path:
        sys.path.insert(0, str(resolved.parent))
    spec = importlib.util.spec_from_file_location(module_name, resolved)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load analyzer: {resolved}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def payload_digest(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def exact_payload_equal(left: Any, right: Any) -> bool:
    return canonical_json(left) == canonical_json(right)


def find_task(
    analyzer: Any,
    manifest: Mapping[str, Any],
    intersubmod_root: Path,
    sample: str,
    chrom: str,
    pos: int,
) -> dict[str, Any]:
    entries = [entry for entry in manifest["samples"] if entry["sample"] == sample]
    if len(entries) != 1:
        raise RuntimeError(f"Expected exactly one manifest entry for {sample}")
    found = None
    for task in analyzer.iter_sample_tasks(entries[0], intersubmod_root):
        screen = task["screen"]
        if screen["chrom"] == chrom and int(screen["pos"]) == pos:
            if found is not None:
                raise RuntimeError("Duplicate requested real-site task")
            found = task
    if found is None:
        raise RuntimeError(f"Real-site fixture not found: {sample} {chrom}:{pos}")
    return found


def synthetic_equivalence(analyzer: Any, workers: int) -> dict[str, Any]:
    methylation = np.asarray(
        [
            [0.05, 0.10, 0.15, 0.10, 0.05, 0.10],
            [0.10, 0.05, 0.10, 0.15, 0.10, 0.05],
            [0.15, 0.10, 0.05, 0.10, 0.15, 0.10],
            [0.08, 0.12, 0.10, 0.08, 0.12, 0.10],
            [0.90, 0.85, 0.95, 0.90, 0.85, 0.95],
            [0.85, 0.95, 0.90, 0.85, 0.95, 0.90],
            [0.95, 0.90, 0.85, 0.95, 0.90, 0.85],
            [0.88, 0.92, 0.90, 0.88, 0.92, 0.90],
        ],
        dtype=float,
    )
    distance = analyzer.F.bernoulli_distance(methylation)
    serial = analyzer.F.analyze_phylo(distance, methylation, seeds=10, rnull=analyzer.F.RNULL)
    parallel = analyzer.analyze_phylo_parallel(
        distance,
        methylation,
        workers=workers,
        seeds=10,
        rnull=analyzer.F.RNULL,
    )
    exact = exact_payload_equal(serial, parallel)
    return {
        "matrix_shape": list(methylation.shape),
        "seeds": 10,
        "rnull": analyzer.F.RNULL,
        "serial_payload_sha256": payload_digest(serial),
        "parallel_payload_sha256": payload_digest(parallel),
        "exact_payload_equal": exact,
    }


def real_nested_equivalence(
    analyzer: Any,
    task: dict[str, Any],
    workers: int,
    min_reads: int,
) -> dict[str, Any]:
    analyzer.configure_phylo_parallel(1, 0)
    serial_started = time.monotonic()
    serial = analyzer.process_site_task(task)
    serial_seconds = time.monotonic() - serial_started
    parallel_started = time.monotonic()
    with ProcessPoolExecutor(
        max_workers=1,
        initializer=analyzer.configure_phylo_parallel,
        initargs=(workers, min_reads),
    ) as executor:
        parallel = executor.submit(analyzer.process_site_task, task).result()
    parallel_seconds = time.monotonic() - parallel_started
    exact = exact_payload_equal(serial, parallel)
    serial_row, serial_detail = serial
    peeled = int(serial_row["n_alt_after_peel"])
    if peeled < min_reads:
        raise RuntimeError(
            f"Real fixture did not trigger nested parallelism: {peeled} < {min_reads}"
        )
    return {
        "sample": serial_row["sample"],
        "chrom": serial_row["chrom"],
        "pos": int(serial_row["pos"]),
        "ref": serial_row["ref"],
        "alt": serial_row["alt"],
        "n_reads_total": int(serial_row["n_reads_total"]),
        "n_alt_after_peel": peeled,
        "parallel_min_reads": min_reads,
        "outer_worker_processes": 1,
        "inner_seed_workers": workers,
        "nested_outer_inner_execution_triggered": True,
        "serial_seconds": serial_seconds,
        "parallel_seconds": parallel_seconds,
        "serial_row_sha256": payload_digest(serial_row),
        "parallel_row_sha256": payload_digest(parallel[0]),
        "serial_detail_sha256": payload_digest(serial_detail),
        "parallel_detail_sha256": payload_digest(parallel[1]),
        "full_row_and_assignment_payload_exact": exact,
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analyzer", type=Path, required=True)
    parser.add_argument("--expected-analyzer-sha256", required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--intersubmod-root", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--pos", type=int, required=True)
    parser.add_argument("--parallel-workers", type=int, default=11)
    parser.add_argument("--parallel-min-reads", type=int, default=200)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> dict[str, Any]:
    if os.path.lexists(args.output):
        raise FileExistsError(f"Refusing to overwrite equivalence receipt: {args.output}")
    if args.parallel_workers < 2 or args.parallel_min_reads < 2:
        raise ValueError("Parallel workers and minimum reads are too small")
    analyzer_path = args.analyzer.resolve(strict=True)
    if sha256(analyzer_path) != args.expected_analyzer_sha256:
        raise RuntimeError("Analyzer SHA-256 differs from the expected pinned identity")
    analyzer = load_analyzer(analyzer_path)
    source_before = {
        "analyzer": artifact(analyzer_path),
        "focal_alt_cluster_lib": artifact(Path(analyzer.F.__file__)),
        "latest_tag_join": artifact(Path(analyzer.TAGS.__file__)),
        "claim_contract_v2": artifact(analyzer.TOPIC_ROOT / "claim-contract-v2.md"),
        "equivalence_auditor": artifact(Path(__file__)),
    }
    started_at = now_utc()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    analyzer.validate_manifest(manifest)
    task = find_task(
        analyzer,
        manifest,
        args.intersubmod_root.resolve(strict=True),
        args.sample,
        args.chrom,
        args.pos,
    )
    synthetic = synthetic_equivalence(analyzer, args.parallel_workers)
    real = real_nested_equivalence(
        analyzer,
        task,
        args.parallel_workers,
        args.parallel_min_reads,
    )
    finished_at = now_utc()
    source_after = {
        "analyzer": artifact(analyzer_path),
        "focal_alt_cluster_lib": artifact(Path(analyzer.F.__file__)),
        "latest_tag_join": artifact(Path(analyzer.TAGS.__file__)),
        "claim_contract_v2": artifact(analyzer.TOPIC_ROOT / "claim-contract-v2.md"),
        "equivalence_auditor": artifact(Path(__file__)),
    }
    checks = {
        "pinned_analyzer_sha256_exact": source_before["analyzer"]["sha256"]
        == args.expected_analyzer_sha256,
        "source_identity_unchanged": source_before == source_after,
        "synthetic_full_payload_exact": synthetic["exact_payload_equal"] is True,
        "real_nested_full_payload_exact": real[
            "full_row_and_assignment_payload_exact"
        ]
        is True,
        "real_nested_parallel_triggered": real["nested_outer_inner_execution_triggered"]
        is True,
    }
    if not all(checks.values()):
        raise RuntimeError(f"Serial/parallel equivalence gate failed: {checks}")
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": finished_at,
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "task_type": "B_comprehensive_validation",
        "scope": "algorithm_and_real_nested_high_read_fixture",
        "command": sys.argv,
        "inputs": {
            "manifest": artifact(args.manifest),
            "intersubmod_root": str(args.intersubmod_root.resolve(strict=True)),
            "source_identity_before": source_before,
            "source_identity_after": source_after,
        },
        "synthetic_fixture": synthetic,
        "real_fixture": real,
        "checks": checks,
        "limitations": [
            "Exact equality is demonstrated for the complete synthetic phylo payload and one real nested outer/inner screen payload.",
            "Identity of the pinned implementation plus deterministic per-run RNG streams is the whole-run equivalence basis; every genome site is not duplicated serially.",
        ],
        "pass_semantics": "algorithmic_execution_equivalence_only_not_scientific_confirmation",
        "pass": True,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(
        json.dumps(
            {
                "output": str(args.output.resolve()),
                "real_fixture": real,
                "checks": checks,
                "pass": True,
            },
            indent=2,
        )
    )
    return payload


def main() -> None:
    run(parse_args())


if __name__ == "__main__":
    main()
