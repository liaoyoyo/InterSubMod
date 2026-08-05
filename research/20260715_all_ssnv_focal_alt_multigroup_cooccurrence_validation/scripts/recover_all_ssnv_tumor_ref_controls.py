#!/usr/bin/env python3
"""Recover a canonical tumor-REF prefix and recompute only its missing suffix."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import multiprocessing
import os
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


TOPIC_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = Path(__file__).resolve().parent
BASE_ANALYZER_PATH = SCRIPT_ROOT / "analyze_all_ssnv_tumor_ref_controls.py"
PINNED_PHYLO_PATH = SCRIPT_ROOT / "analyze_all_ssnv_focal_alt_multigroup.pinned_390228ce.py"
EXPECTED_BASE_ANALYZER_SHA256 = "95bf7cdca5b636eb0905693ee3c35f1bab699cf4872370fe6f22157cac4c8b87"
EXPECTED_PINNED_PHYLO_SHA256 = "390228ce1fb98b59409f2517341334ab94329ef6c73a4838555fca85878027b2"
EXPECTED_FOCAL_LIB_SHA256 = "ed5f3c99461248248b20a9f49597ab5de7340a1e2055d77a1d83dbcc2799b72a"
SCHEMA_NAME = "intersubmod.all_ssnv_tumor_ref_controls"
SCHEMA_VERSION = "1.0.0"


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))
B = load_module("tumor_ref_serial_base_for_recovery", BASE_ANALYZER_PATH)
P = load_module("pinned_phylo_for_tumor_ref_recovery", PINNED_PHYLO_PATH)
F = B.F


INNER_PHYLO_WORKERS = 11
INNER_PHYLO_MIN_READS = 200


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    stat = resolved.stat()
    return {
        **artifact(resolved),
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "mtime_ns": stat.st_mtime_ns,
        "ctime_ns": stat.st_ctime_ns,
    }


def verify_dependencies() -> dict[str, dict[str, Any]]:
    expected = {
        "serial_base_analyzer": (BASE_ANALYZER_PATH, EXPECTED_BASE_ANALYZER_SHA256),
        "pinned_seed_parallel_phylo": (PINNED_PHYLO_PATH, EXPECTED_PINNED_PHYLO_SHA256),
        "focal_alt_cluster_lib": (Path(F.__file__).resolve(), EXPECTED_FOCAL_LIB_SHA256),
    }
    observed = {role: artifact(path) for role, (path, _sha) in expected.items()}
    drift = {
        role: {"expected": expected_sha, "observed": observed[role]["sha256"]}
        for role, (_path, expected_sha) in expected.items()
        if observed[role]["sha256"] != expected_sha
    }
    if drift:
        raise RuntimeError(f"Recovery dependency SHA-256 drift: {drift}")
    return observed


def parse_bool(value: Any, *, field: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"Invalid boolean for {field}: {value!r}")


def typed_prefix_row(row: Mapping[str, str]) -> dict[str, Any]:
    output: dict[str, Any] = dict(row)
    boolean_fields = {
        "ref_evaluable", "ref_stable_null_multigroup",
        "joint_evaluable", "joint_stable_null_multigroup",
        "joint_allele_testable", "joint_allele_axis_aligned",
    }
    optional_boolean_fields = {"ref_unstable", "joint_unstable"}
    integer_fields = {
        "primary_assignment_n_core_reads", "n_tumor_alt", "n_tumor_ref",
        "n_tumor_alt_matrix", "n_tumor_ref_matrix", "ref_n_after_peel", "ref_n_peeled",
        "joint_n_after_peel", "joint_n_peeled", "joint_allele_n",
    }
    optional_integer_fields = {"ref_coarse_ng", "joint_coarse_ng"}
    optional_float_fields = {
        "ref_modal_fraction", "joint_modal_fraction", "joint_allele_v", "joint_allele_p_perm"
    }
    for field in boolean_fields:
        output[field] = parse_bool(row[field], field=field)
    for field in optional_boolean_fields:
        output[field] = None if row[field] == "" else parse_bool(row[field], field=field)
    for field in integer_fields:
        output[field] = int(row[field])
    for field in optional_integer_fields:
        output[field] = None if row[field] == "" else int(row[field])
    for field in optional_float_fields:
        output[field] = None if row[field] == "" else float(row[field])
    return output


def salvage_canonical_prefix(
    path: Path,
    fields: Sequence[str],
    tasks: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    accepted: list[dict[str, Any]] = []
    reached_truncated_eof = False
    source_before = file_identity(path)
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if list(reader.fieldnames or ()) != list(fields):
                raise RuntimeError("Live prefix header differs from canonical tumor-REF fields")
            for row_number, row in enumerate(reader, 1):
                if None in row or set(row) != set(fields) or any(row[field] is None for field in fields):
                    raise RuntimeError(f"Malformed complete prefix row at data row {row_number}")
                if row_number > len(tasks):
                    raise RuntimeError("Live prefix contains more rows than canonical task set")
                expected_site = tasks[row_number - 1]["site"]
                observed_key = B.site_row_key(row)
                expected_key = B.site_row_key(expected_site)
                if observed_key != expected_key:
                    raise RuntimeError(
                        f"Live prefix is not canonical at row {row_number}: "
                        f"observed={observed_key} expected={expected_key}"
                    )
                source_drift = [field for field in expected_site if row[field] != expected_site[field]]
                if source_drift:
                    raise RuntimeError(
                        f"Live prefix/source-screen field drift at row {row_number}: {source_drift[:5]}"
                    )
                if int(row["primary_assignment_n_core_reads"]) != int(
                    tasks[row_number - 1]["primary_assignment_n_core_reads"]
                ):
                    raise RuntimeError(f"Live prefix assignment count drift at row {row_number}")
                if row["primary_assignment_labels_sha256"] != tasks[row_number - 1][
                    "primary_assignment_labels_sha256"
                ]:
                    raise RuntimeError(f"Live prefix assignment digest drift at row {row_number}")
                accepted.append(typed_prefix_row(row))
    except (EOFError, gzip.BadGzipFile, zlib.error):
        reached_truncated_eof = True
    source_after = file_identity(path)
    keys_payload = "\n".join(":\t".join(map(str, B.site_row_key(row))) for row in accepted)
    return accepted, {
        "source_before": source_before,
        "source_after": source_after,
        "source_changed_while_reading_allowed_for_live_producer": source_before != source_after,
        "accepted_complete_rows": len(accepted),
        "accepted_prefix_keys_sha256": hashlib.sha256(keys_payload.encode()).hexdigest(),
        "truncated_gzip_eof_observed": reached_truncated_eof,
        "complete_rows_only": True,
        "canonical_task_order_exact": True,
        "source_screen_fields_exact": True,
        "assignment_count_and_digest_exact": True,
    }


def analyze_subset_parallel(
    read_ids: Sequence[str],
    distance_ids: Sequence[str],
    distance: np.ndarray,
    methylation_ids: Sequence[str],
    methylation: np.ndarray,
    seed: int,
) -> dict[str, Any]:
    distance_index = {read_id: index for index, read_id in enumerate(distance_ids)}
    methylation_index = {read_id: index for index, read_id in enumerate(methylation_ids)}
    usable = [
        read_id for read_id in read_ids if read_id in distance_index and read_id in methylation_index
    ]
    if len(read_ids) < 2 * F.MIN_SIZE:
        return B._not_evaluable(
            "insufficient_reads",
            f"fewer_than_2x_MIN_SIZE_raw_reads:{len(read_ids)}<{2 * F.MIN_SIZE}",
            len(read_ids),
            len(usable),
        )
    if len(usable) < 2 * F.MIN_SIZE:
        return B._not_evaluable(
            "insufficient_matrix_reads",
            f"fewer_than_2x_MIN_SIZE_matrix_reads:{len(usable)}<{2 * F.MIN_SIZE}",
            len(read_ids),
            len(usable),
        )
    distance_rows = [distance_index[read_id] for read_id in usable]
    sub_distance = distance[np.ix_(distance_rows, distance_rows)]
    kept = F.peel_complete(sub_distance)
    kept_ids = [usable[index] for index in kept]
    if len(kept_ids) < 2 * F.MIN_SIZE:
        return B._not_evaluable(
            "incomplete_distance_below_minimum",
            f"fewer_than_2x_MIN_SIZE_after_peel:{len(kept_ids)}<{2 * F.MIN_SIZE}",
            len(read_ids),
            len(usable),
            len(kept_ids),
        )
    sub_distance = sub_distance[np.ix_(kept, kept)]
    sub_methylation = methylation[[methylation_index[read_id] for read_id in kept_ids]]
    if len(kept_ids) >= INNER_PHYLO_MIN_READS and INNER_PHYLO_WORKERS > 1:
        phylo = P.analyze_phylo_parallel(
            sub_distance,
            sub_methylation,
            workers=INNER_PHYLO_WORKERS,
            base_seed=seed,
            seeds=10,
            rnull=F.RNULL,
            null_mode="column",
            empirical_alpha=None,
        )
    else:
        phylo = F.analyze_phylo(
            sub_distance,
            sub_methylation,
            base_seed=seed,
            seeds=10,
            rnull=F.RNULL,
            null_mode="column",
            empirical_alpha=None,
        )
    labels = phylo["coarse_labels"]
    return {
        "status": "evaluable",
        "evaluable": True,
        "not_testable_reason": None,
        "n_raw": len(read_ids),
        "n_matrix": len(usable),
        "n_after_peel": len(kept_ids),
        "n_peeled": len(usable) - len(kept_ids),
        "kept_ids": kept_ids,
        "labels": labels,
        "coarse_ng": phylo["coarse_ng"],
        "modal_fraction": phylo["modal_fraction"],
        "unstable": phylo["unstable"],
        "stable_null_multigroup": phylo["coarse_ng"] >= 2 and not phylo["unstable"],
        "cluster_sizes": dict(
            B.Counter(label for label in labels if label not in {"other", "outlier"})
        ),
    }


B.analyze_subset = analyze_subset_parallel


def configure_inner_workers(workers: int, min_reads: int) -> None:
    if workers < 1 or min_reads < 2 * F.MIN_SIZE:
        raise ValueError("Invalid inner phylo worker configuration")
    global INNER_PHYLO_WORKERS, INNER_PHYLO_MIN_READS
    INNER_PHYLO_WORKERS = workers
    INNER_PHYLO_MIN_READS = min_reads


def process_suffix(
    tasks: Sequence[dict[str, Any]], workers: int, chunk_size: int, max_pending: int
) -> Iterable[list[dict[str, Any]]]:
    if "fork" not in multiprocessing.get_all_start_methods():
        raise RuntimeError("Tumor-REF suffix recovery requires Linux fork")
    context = multiprocessing.get_context("fork")
    chunks = B.chunked(tasks, chunk_size)
    with ProcessPoolExecutor(max_workers=workers, mp_context=context) as executor:
        yield from B.bounded_chunk_results(executor, chunks, max_pending)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-results", type=Path, required=True)
    parser.add_argument("--stable-assignments", type=Path, required=True)
    parser.add_argument("--primary-artifact-audit-pre", type=Path, required=True)
    parser.add_argument("--live-prefix-sites", type=Path, required=True)
    parser.add_argument("--live-prefix-source-snapshot", type=Path, required=True)
    parser.add_argument("--serial-parallel-equivalence-receipt", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-prefix-rows", type=int, default=102_000)
    parser.add_argument("--expected-total-rows", type=int, default=102_842)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--chunk-size", type=int, default=1)
    parser.add_argument("--max-pending-chunks", type=int, default=6)
    parser.add_argument("--inner-phylo-workers", type=int, default=11)
    parser.add_argument("--inner-phylo-min-reads", type=int, default=200)
    parser.add_argument("--progress-every", type=int, default=100)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if any(
        value < 1
        for value in (
            args.workers, args.chunk_size, args.max_pending_chunks,
            args.inner_phylo_workers, args.progress_every,
        )
    ):
        raise ValueError("Worker, chunk, pending, and progress values must be positive")
    configure_inner_workers(args.inner_phylo_workers, args.inner_phylo_min_reads)
    started_at = B.now_utc()
    dependencies_before = verify_dependencies()
    self_before = artifact(Path(__file__).resolve())
    immutable_inputs_before = {
        "site_results": artifact(args.site_results),
        "stable_assignments": artifact(args.stable_assignments),
        "primary_artifact_audit_pre": artifact(args.primary_artifact_audit_pre),
        "live_prefix_source_snapshot": artifact(args.live_prefix_source_snapshot),
        "serial_parallel_equivalence_receipt": artifact(args.serial_parallel_equivalence_receipt),
    }
    equivalence = json.loads(args.serial_parallel_equivalence_receipt.read_text(encoding="utf-8"))
    if equivalence.get("pass") is not True:
        raise RuntimeError("Serial/parallel exact-equivalence receipt is not passing")
    source_fields, tasks = B.load_primary_inputs(args.site_results, args.stable_assignments)
    if len(tasks) != args.expected_total_rows:
        raise RuntimeError(f"Canonical task count drift: {len(tasks)} != {args.expected_total_rows}")
    fields = [*source_fields, *B.CONTROL_FIELDS]
    prefix_rows, prefix_audit = salvage_canonical_prefix(
        args.live_prefix_sites, fields, tasks
    )
    if len(prefix_rows) < args.minimum_prefix_rows or len(prefix_rows) >= len(tasks):
        raise RuntimeError(
            f"Recovered prefix length outside expected range: {len(prefix_rows)}"
        )
    suffix_tasks = tasks[len(prefix_rows):]
    B.create_output_dir(args.output_dir)
    site_path = args.output_dir / "all_ssnv_tumor_ref_control_site_results.tsv.gz"
    summary_path = args.output_dir / "all_ssnv_tumor_ref_control_summary.json"
    manifest_path = args.output_dir / "run_manifest.json"
    prefix_receipt_path = args.output_dir / "recovered_prefix_receipt.json"
    summaries = B.SummaryBook()
    dedup = B.DedupAnnotator()
    processed = 0
    suffix_processed = 0
    with gzip.open(site_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in prefix_rows:
            dedup.annotate(row)
            writer.writerow({field: B._tsv_value(row[field]) for field in fields})
            summaries.add(row)
            processed += 1
        for result_chunk in process_suffix(
            suffix_tasks, args.workers, args.chunk_size, args.max_pending_chunks
        ):
            for row in result_chunk:
                expected_key = B.site_row_key(tasks[processed]["site"])
                if B.site_row_key(row) != expected_key:
                    raise RuntimeError(
                        f"Suffix result order drift: {B.site_row_key(row)} != {expected_key}"
                    )
                dedup.annotate(row)
                writer.writerow({field: B._tsv_value(row[field]) for field in fields})
                summaries.add(row)
                processed += 1
                suffix_processed += 1
                if suffix_processed % args.progress_every == 0 or processed == len(tasks):
                    print(
                        f"suffix_processed={suffix_processed}/{len(suffix_tasks)} "
                        f"total_processed={processed}/{len(tasks)}",
                        flush=True,
                    )
    if processed != len(tasks) or suffix_processed != len(suffix_tasks):
        raise RuntimeError("Recovered tumor-REF row-count conservation failed")
    prefix_receipt = {
        "schema_name": "intersubmod.tumor_ref_live_prefix_recovery",
        "schema_version": "1.0.0",
        "created_at_utc": B.now_utc(),
        "source": prefix_audit,
        "accepted_rows": len(prefix_rows),
        "recomputed_suffix_rows": len(suffix_tasks),
        "final_rows": processed,
        "pass": True,
    }
    prefix_receipt_path.write_text(
        json.dumps(prefix_receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    summary = {
        "schema_name": f"{SCHEMA_NAME}.summary",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": B.now_utc(),
        "task_type": "B_comprehensive_validation",
        "scope": "complete_primary_stable_null_multigroup_set",
        "selection_contract": {
            "site_results_filter": "stable_null_multigroup=true",
            "site_key": ["sample", "chrom", "pos", "ref", "alt"],
            "site_assignment_key_sets_exact": True,
            "site_assignment_posthoc_metadata_exact": True,
            "selected_site_duplicates": 0,
            "assignment_duplicates": 0,
            "primary_alt_labels_rewritten": False,
        },
        "clustering_contract": {
            "screen": B.SCREEN_CONTRACT,
            "minimum_group_size": F.MIN_SIZE,
            "rnull": F.RNULL,
            "seeds": 10,
            "null_mode": "column",
            "ref_selection": "reads.tsv is_tumor=true and alt_support=REF",
            "joint_selection": "reads.tsv is_tumor=true and alt_support in {ALT,REF}",
            "joint_allele_association": (
                "posthoc Cramer's V with 499 label permutations after stable joint clustering"
            ),
            "seed_parallel_exact_equivalent": True,
        },
        "denominator_contract": {
            "site_weighted": "one row per dataset site",
            "component_deduplicated": ["biological_id", "component_id"],
            "readset_deduplicated": ["biological_id", "alt_readset_sha256"],
            "missing_component_or_readset": "site-specific fail-open identity; never merged",
        },
        "recovery": {
            "canonical_prefix_rows": len(prefix_rows),
            "recomputed_suffix_rows": len(suffix_tasks),
            "final_rows": processed,
            "prefix_receipt": artifact(prefix_receipt_path),
            "serial_parallel_equivalence_receipt": artifact(
                args.serial_parallel_equivalence_receipt
            ),
            "prefix_and_suffix_disjoint_by_canonical_order": True,
        },
        **summaries.payload(),
        "guardrail": {
            "ref_nonreplication": "Tumor-REF nonreplication does not prove a subclone.",
            "ref_replication": "Tumor-REF replication weakens ALT-specificity.",
        },
        "pass": True,
    }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    immutable_inputs_after = {
        "site_results": artifact(args.site_results),
        "stable_assignments": artifact(args.stable_assignments),
        "primary_artifact_audit_pre": artifact(args.primary_artifact_audit_pre),
        "live_prefix_source_snapshot": artifact(args.live_prefix_source_snapshot),
        "serial_parallel_equivalence_receipt": artifact(args.serial_parallel_equivalence_receipt),
    }
    dependencies_after = verify_dependencies()
    self_after = artifact(Path(__file__).resolve())
    if immutable_inputs_after != immutable_inputs_before:
        raise RuntimeError("Immutable recovery inputs changed during execution")
    if dependencies_after != dependencies_before or self_after != self_before:
        raise RuntimeError("Recovery source code changed during execution")
    finished_at = B.now_utc()
    manifest = {
        "schema_name": f"{SCHEMA_NAME}.run_manifest",
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": finished_at,
        "started_at_utc": started_at,
        "finished_at_utc": finished_at,
        "command": sys.argv,
        "inputs": {
            **immutable_inputs_after,
            "live_prefix_sites_observed": prefix_audit["source_after"],
        },
        "source_code": {
            "analyzer": self_after,
            "focal_alt_cluster_lib": dependencies_after["focal_alt_cluster_lib"],
        },
        "source_dependencies": {
            "serial_base_analyzer": dependencies_after["serial_base_analyzer"],
            "pinned_seed_parallel_phylo": dependencies_after["pinned_seed_parallel_phylo"],
        },
        "execution": {
            "workers": args.workers,
            "chunk_size": args.chunk_size,
            "max_pending_chunks": args.max_pending_chunks,
            "inner_phylo_workers": args.inner_phylo_workers,
            "inner_phylo_min_reads": args.inner_phylo_min_reads,
            "per_site_future_submission": args.chunk_size == 1,
            "site_exceptions_hard_fail": True,
            "linux_fork_required": True,
        },
        "recovery": prefix_receipt,
        "contracts": {
            "complete_primary_stable_set_only": True,
            "site_assignment_keys_exact": True,
            "site_assignment_posthoc_metadata_exact": True,
            "primary_alt_labels_rewritten": False,
            "existing_results_overwritten": False,
            "canonical_prefix_only": True,
            "prefix_suffix_key_sets_disjoint": True,
            "serial_parallel_exact_equivalence_required": True,
            "guardrail": summary["guardrail"],
        },
        "outputs": {
            "site_results": artifact(site_path),
            "summary": artifact(summary_path),
            "prefix_receipt": artifact(prefix_receipt_path),
        },
        "counts": {
            "primary_stable_sites": len(tasks),
            "processed_sites": processed,
            "recovered_prefix_sites": len(prefix_rows),
            "recomputed_suffix_sites": suffix_processed,
        },
        "pass": True,
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output_dir": str(args.output_dir.resolve()),
                "site_results": str(site_path.resolve()),
                "summary": str(summary_path.resolve()),
                "run_manifest": str(manifest_path.resolve()),
                "recovered_prefix_sites": len(prefix_rows),
                "recomputed_suffix_sites": suffix_processed,
                "processed_sites": processed,
                "pass": True,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
