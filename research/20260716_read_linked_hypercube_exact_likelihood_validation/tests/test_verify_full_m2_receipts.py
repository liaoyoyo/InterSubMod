#!/usr/bin/env python3
import argparse
import csv
import copy
import gzip
import hashlib
import inspect
import json
import math
import os
import pathlib
import shutil
import sys
import tempfile
import time
import unittest
from unittest import mock


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))

import verify_full_m2_receipts as verifier  # noqa: E402
import run_full_m2_extraction as production_extraction  # noqa: E402
import run_full_m2_ranking as production_ranking  # noqa: E402
from test_full_m2_extraction import make_release_contract  # noqa: E402


def make_verifier_release_contract(root: pathlib.Path):
    manifest, paths = make_release_contract(root)
    document = json.loads(manifest.read_text(encoding="utf-8"))
    repo_root = pathlib.Path(document["source_snapshot"]["repo_root"])
    for entry in document["source_snapshot"]["entries"]:
        source_path = repo_root / entry["repo_relative_path"]
        source_path.parent.mkdir(parents=True, exist_ok=True)
        snapshot_path = root / entry["snapshot"]["path"]
        source_path.write_bytes(snapshot_path.read_bytes())
        source_path.chmod(0o444)
        observed = source_path.stat()
        entry["source"] = {
            "path": str(source_path.resolve()),
            "st_dev": int(observed.st_dev), "st_ino": int(observed.st_ino),
            "st_nlink": int(observed.st_nlink), "size_bytes": int(observed.st_size),
            "mtime_ns": int(observed.st_mtime_ns), "ctime_ns": int(observed.st_ctime_ns),
            "mode_octal": "0444", "sha256": verifier.sha256_path(source_path),
        }
    by_role = {
        row["role"]: row for row in document["source_snapshot"]["entries"]
    }
    input_verifier = by_role["input_identity_verifier"]["source"]
    document["fresh_input_identity_verification"]["verifier_path"] = input_verifier["path"]
    document["fresh_input_identity_verification"]["verifier_sha256"] = input_verifier["sha256"]
    freezer = by_role["release_contract_freezer"]
    document["producer"].update({
        "source_sha256": freezer["source"]["sha256"],
        "immutable_copy_path": freezer["snapshot"]["path"],
        "immutable_copy_sha256": freezer["snapshot"]["sha256"],
    })
    _rewrite_immutable_json(manifest, document)
    return manifest, paths


def _source_identity(path: pathlib.Path) -> dict:
    return {"path": str(path.resolve()), "sha256": verifier.sha256_path(path)}


def _rewrite_immutable_json(path: pathlib.Path, payload: dict) -> None:
    sidecar = path.with_name(f"{path.name}.sha256")
    path.chmod(0o644)
    sidecar.chmod(0o644)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    sidecar.write_text(f"{verifier.sha256_path(path)}  {path.name}\n", encoding="ascii")
    path.chmod(0o444)
    sidecar.chmod(0o444)


def _minimal_release_binding(root: pathlib.Path) -> tuple[dict, pathlib.Path]:
    release_root = root / "release"
    release_root.mkdir(parents=True)
    release_path = release_root / "m2_run_manifest.json"
    release = {
        "runtime": {
            "python": {
                "executable": str(pathlib.Path(sys.executable).resolve()),
                "version": "test", "implementation": "CPython",
            },
            "packages": {"numpy": "test", "scipy": "test", "pysam": "test"},
            "samtools": {
                "path": "/test/samtools", "version_line": "test",
                "htslib_version_line": "test",
            },
            "platform": "test",
        },
        "receipt_integrity": {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{release_path.name}.sha256", "covers": release_path.name,
        },
    }
    production_extraction.write_immutable_json_exclusive(release_path, release)
    canonical = release_root / "canonical_manifest.json"
    canonical.write_text('{"canonical":true}\n', encoding="utf-8")
    canonical.chmod(0o444)
    release_sidecar = release_path.with_name(f"{release_path.name}.sha256")
    sources = {
        "full_extraction_runner": _source_identity(production_extraction.RUNNER),
        "extractor": _source_identity(production_extraction.EXTRACTOR),
        "lossless_read_contract": _source_identity(production_extraction.LOSSLESS_READ_CONTRACT),
        "full_ranking_runner": _source_identity(production_ranking.RUNNER),
        "ranker": _source_identity(production_ranking.RANKER),
        "hypercube_solver": _source_identity(production_ranking.HYPERCUBE_SOLVER),
    }
    binding = {
        "schema_name": "intersubmod.m2_release_binding",
        "schema_version": "1.0.0",
        "release_manifest": {
            "path": str(release_path.resolve()),
            "sha256": verifier.sha256_path(release_path),
            "semantic_sha256": verifier.semantic_json_sha256(release),
            "sidecar": {
                "path": str(release_sidecar.resolve()),
                "sha256": verifier.sha256_path(release_sidecar),
            },
        },
        "authority_mode": "CANONICAL_V5_FROZEN",
        "validation_evidence_eligible": True,
        "canonical_input_manifest": {
            "path": str(canonical.resolve()), "sha256": verifier.sha256_path(canonical),
        },
        "snapshot_sources": sources,
        "frozen_parameters": {},
        "frozen_parameters_semantic_sha256": verifier.semantic_json_sha256({}),
        "deep_release_verification": {"all_pass": True},
    }
    return binding, canonical


def build_production_extraction_orchestration_fixture(root: pathlib.Path) -> dict:
    """Generate the complete 154-task chain through production runner helpers."""
    binding, manifest = _minimal_release_binding(root)
    output_root = root / "extraction"
    thresholds = [1, 2, 3, 5]
    parameters = {
        "mapq_min": 20, "baseq_min": 20, "samtools_threads": 1,
        "bridge_thresholds": thresholds,
        "component_linkage_bases": list(production_extraction.LINKAGE_BASES),
    }
    run_contract = {
        "manifest_sha256": verifier.sha256_path(manifest),
        "extractor_sha256": verifier.sha256_path(production_extraction.EXTRACTOR),
        "extractor": _source_identity(production_extraction.EXTRACTOR),
        "lossless_read_contract": _source_identity(production_extraction.LOSSLESS_READ_CONTRACT),
        "runner": _source_identity(production_extraction.RUNNER),
        "workers": 2, "samtools_threads": 1, "mapq_min": 20, "baseq_min": 20,
        "bridge_thresholds": thresholds,
        "component_linkage_bases": list(production_extraction.LINKAGE_BASES),
        "task_timeout_seconds": 28800, "timeout_grace_seconds": 300,
        "release_binding": binding,
        "orchestration_policy": dict(production_extraction.ORCHESTRATION_POLICY),
    }
    output_root.mkdir(parents=True)
    output_stat = output_root.stat()
    session_target = {
        "output_root": {
            "path": str(output_root.resolve()), "st_dev": int(output_stat.st_dev),
            "st_ino": int(output_stat.st_ino),
        },
        "release_manifest": dict(binding["release_manifest"]),
        "run_contract_semantic_sha256": production_extraction.semantic_json_sha256(run_contract),
    }
    _, gate = production_extraction.create_resource_gate_receipt(
        output_root, stage="extraction", gate_scope="session", target=session_target,
        producer_source=production_extraction.release_producer_sources(binding)["runner"],
        conflicts={"process_count": 0, "root_count": 0, "representatives": []},
    )
    session = production_extraction.ensure_release_orchestration_session(
        output_root, binding, run_contract, gate
    )
    args = argparse.Namespace(
        mapq_min=20, baseq_min=20, bridge_thresholds="1,2,3,5", samtools_threads=1,
    )
    pairs = list(verifier.expected_keys(verifier.DATASETS, verifier.AUTOSOMES))
    results: list[dict] = []
    previous = None
    previous_count = 0
    child_receipts: dict[tuple[str, str], dict] = {}
    for batch_index, count in enumerate(verifier.ORCHESTRATION_COMPLETED_COUNTS, start=1):
        specs = []
        for dataset, chrom in pairs[previous_count:count]:
            task_dir = output_root / "samples" / dataset / chrom
            command = production_extraction.task_command(
                manifest, dataset, chrom, task_dir, args
            )
            specs.append((dataset, chrom, task_dir, command, parameters))
        batch_target = {
            "output_root": dict(session["output_root"]),
            "session_id": session["session_id"],
            "session_sha256": verifier.sha256_path(output_root / "_orchestration/session.json"),
            "batch_index": batch_index,
            "before_count": previous_count,
            "max_new_tasks": 8 if batch_index == 1 else 16,
            "effective_count": len(specs),
            "selected_task_ids": [f"{spec[0]}:{spec[1]}" for spec in specs],
            "previous_chain_head": previous,
        }
        _, batch_gate = production_extraction.create_resource_gate_receipt(
            output_root, stage="extraction", gate_scope="batch", batch_index=batch_index,
            target=batch_target, producer_source=session["producer_sources"]["runner"],
            conflicts={"process_count": 0, "root_count": 0, "representatives": []},
        )
        batch, grants = production_extraction.create_batch_start_and_grants(
            output_root, session, run_contract, specs,
            before_count=previous_count, previous_chain_head=previous,
            batch_index=batch_index, max_new_tasks=(8 if batch_index == 1 else 16),
            gate=batch_gate,
        )
        current_results = []
        for dataset, chrom, task_dir, _command, _parameters in specs:
            task_dir.mkdir(parents=True, exist_ok=True)
            output = task_dir / "tiny.tsv.gz"
            output.write_bytes(f"{dataset}\t{chrom}\n".encode("ascii"))
            child_path = task_dir / "receipt.json"
            child = {
                "outputs": {
                    output.name: {
                        "path": str(output.resolve()), "size_bytes": output.stat().st_size,
                        "sha256": verifier.sha256_path(output),
                    },
                },
                "receipt_integrity": {
                    "scheme": "external_sha256_sidecar_v1",
                    "sidecar_name": f"{child_path.name}.sha256", "covers": child_path.name,
                },
            }
            production_extraction.write_immutable_json_exclusive(child_path, child)
            child_receipts[(dataset, chrom)] = child
            result = {
                "dataset": dataset, "chrom": chrom, "status": "PASS", "returncode": 0,
                "timed_out": False, "process_group_id": os.getpid(),
                "started_monotonic_ns": time.monotonic_ns(), "receipt": child,
            }
            task_id = f"{dataset}:{chrom}"
            result["orchestration_completion"] = production_extraction.write_child_completion(
                output_root, session, batch, grants[task_id], result
            )
            current_results.append(result)
        # Production workers may complete out of order; publication must be canonical.
        current_results.sort(
            key=lambda row: pairs.index((row["dataset"], row["chrom"]))
        )
        results.extend(current_results)
        path = (
            output_root / "full_extraction_receipt.json" if count == 154
            else output_root / "checkpoints" / f"checkpoint_{count:03d}_of_154.json"
        )
        chain = {
            "schema_name": "intersubmod.m2_full_extraction_receipt",
            "schema_version": "1.2.0", "run_contract": run_contract,
            "results": copy.deepcopy(results), "all_pass": count == 154,
            "orchestration": {
                "session_identity": production_extraction._session_identity(output_root, session),
                "batch_start_identity": {
                    "path": batch["path"], "sha256": batch["sha256"],
                    "batch_id": batch["batch_id"], "batch_index": batch["batch_index"],
                },
                "previous_chain_head": previous,
                "batch_completion_attestations": [
                    {"task_id": f"{row['dataset']}:{row['chrom']}",
                     **row["orchestration_completion"]}
                    for row in current_results
                ],
                "cumulative_attested_tasks": count,
            },
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": f"{path.name}.sha256", "covers": path.name,
            },
        }
        production_extraction.write_immutable_json_exclusive(path, chain)
        previous = {"path": str(path.resolve()), "sha256": verifier.sha256_path(path)}
        previous_count = count
    return {
        "root": output_root, "binding": binding, "full": chain,
        "children": child_receipts, "manifest": manifest,
    }


def build_production_ranking_orchestration_fixture(
    root: pathlib.Path, extraction: dict
) -> dict:
    output_root = root / "ranking"
    args = argparse.Namespace(
        structural_exact_pattern_minreads=(1, 2, 3, 5),
        primary_structural_exact_pattern_minread=3,
        exact_k_max=12, max_vertex_sets=256,
        solver_time_limit_seconds=30.0,
        fixed_error_grid_values=(0.005, 0.01, 0.02, 0.05),
        minimum_bq_error_rate=0.000001, maximum_bq_error_rate=0.25,
        conditional_candidate_ranking_bootstrap_replicates=20,
        conditional_candidate_ranking_bootstrap_seed=20260716,
        tie_tolerance=0.000001,
        thresholds=(1, 2, 3, 5), component_bases=("PS_HP1", "PS_HP2"),
        hp_families=("1", "2"),
    )
    parameters = production_ranking._rank_parameters(args)
    extraction_terminal = extraction["root"] / "full_extraction_receipt.json"
    run_contract = {
        "full_extraction_receipt": {
            "path": str(extraction_terminal.resolve()),
            "sha256": verifier.sha256_path(extraction_terminal),
        },
        "ranker": _source_identity(production_ranking.RANKER),
        "hypercube_solver": _source_identity(production_ranking.HYPERCUBE_SOLVER),
        "runner": _source_identity(production_ranking.RUNNER),
        "method_contract": dict(production_ranking.EXPECTED_METHOD_CONTRACT),
        "thresholds": list(args.thresholds),
        "component_bases": list(args.component_bases),
        "hp_families": list(args.hp_families),
        "parameters": parameters, "workers": 2,
        "task_timeout_seconds": 28800, "timeout_grace_seconds": 300,
        "release_binding": extraction["binding"],
        "orchestration_policy": dict(production_ranking.ORCHESTRATION_POLICY),
    }
    parent = production_ranking.extraction_parent_identity(
        extraction["root"], extraction_terminal
    )
    output_root.mkdir(parents=True)
    output_stat = output_root.stat()
    session_target = {
        "output_root": {
            "path": str(output_root.resolve()), "st_dev": int(output_stat.st_dev),
            "st_ino": int(output_stat.st_ino),
        },
        "release_manifest": dict(extraction["binding"]["release_manifest"]),
        "run_contract_semantic_sha256": production_ranking.semantic_json_sha256(run_contract),
        "parent_extraction": dict(parent),
    }
    _, session_gate = production_ranking.create_resource_gate_receipt(
        output_root, stage="ranking", gate_scope="session", target=session_target,
        producer_source=production_ranking.release_producer_sources(
            extraction["binding"]
        )["runner"],
        conflicts={"process_count": 0, "root_count": 0, "representatives": []},
    )
    session = production_ranking.ensure_release_orchestration_session(
        output_root, extraction["binding"], run_contract, parent, session_gate
    )
    pairs = list(verifier.expected_keys(verifier.DATASETS, verifier.AUTOSOMES))
    results: list[dict] = []
    children: dict[tuple[str, str], dict] = {}
    previous = None
    previous_count = 0
    for batch_index, count in enumerate(verifier.ORCHESTRATION_COMPLETED_COUNTS, start=1):
        specs = []
        for dataset, chrom in pairs[previous_count:count]:
            extraction_dir = extraction["root"] / "samples" / dataset / chrom
            output_dir = output_root / "samples" / dataset / chrom
            command = production_ranking.task_command(extraction_dir, output_dir, args)
            specs.append((
                dataset, chrom, extraction_dir, output_dir, command,
                production_ranking.RANKER, verifier.sha256_path(production_ranking.RANKER),
                parameters,
            ))
        batch_target = {
            "output_root": dict(session["output_root"]),
            "session_id": session["session_id"],
            "session_sha256": verifier.sha256_path(output_root / "_orchestration/session.json"),
            "batch_index": batch_index,
            "before_count": previous_count,
            "max_new_tasks": 8 if batch_index == 1 else 16,
            "effective_count": len(specs),
            "selected_task_ids": [f"{spec[0]}:{spec[1]}" for spec in specs],
            "previous_chain_head": previous,
        }
        _, batch_gate = production_ranking.create_resource_gate_receipt(
            output_root, stage="ranking", gate_scope="batch", batch_index=batch_index,
            target=batch_target, producer_source=session["producer_sources"]["runner"],
            conflicts={"process_count": 0, "root_count": 0, "representatives": []},
        )
        batch, grants = production_ranking.create_batch_start_and_grants(
            output_root, session, run_contract, specs,
            before_count=previous_count, previous_chain_head=previous,
            batch_index=batch_index, max_new_tasks=(8 if batch_index == 1 else 16),
            gate=batch_gate,
        )
        current_results = []
        for dataset, chrom, extraction_dir, output_dir, _command, *_rest in specs:
            output_dir.mkdir(parents=True, exist_ok=True)
            output = output_dir / "rank.tiny.tsv.gz"
            output.write_bytes(f"rank\t{dataset}\t{chrom}\n".encode("ascii"))
            child_path = output_dir / "receipt.json"
            extraction_path = extraction_dir / "receipt.json"
            child = {
                "parameters": dict(parameters),
                "provenance": {
                    "upstream_extraction_receipt": {
                        "path": str(extraction_path.resolve()),
                        "sha256": verifier.sha256_path(extraction_path),
                    },
                },
                "outputs": {
                    output.name: {
                        "path": str(output.resolve()), "size_bytes": output.stat().st_size,
                        "sha256": verifier.sha256_path(output),
                    },
                },
                "receipt_integrity": {
                    "scheme": "external_sha256_sidecar_v1",
                    "sidecar_name": f"{child_path.name}.sha256", "covers": child_path.name,
                },
            }
            production_ranking.write_immutable_json_exclusive(child_path, child)
            children[(dataset, chrom)] = child
            result = {
                "dataset": dataset, "chrom": chrom, "status": "PASS", "returncode": 0,
                "timed_out": False, "process_group_id": os.getpid(),
                "started_monotonic_ns": time.monotonic_ns(),
                "receipt_path": str(child_path.resolve()),
                "rank_receipt": {
                    "path": str(child_path.resolve()), "sha256": verifier.sha256_path(child_path),
                },
            }
            task_id = f"{dataset}:{chrom}"
            result["orchestration_completion"] = production_ranking.write_child_completion(
                output_root, session, batch, grants[task_id], result
            )
            current_results.append(result)
        current_results = production_ranking.canonical_sort_results(current_results)
        results.extend(current_results)
        path = (
            output_root / "full_ranking_receipt.json" if count == 154
            else output_root / "checkpoints" / f"checkpoint_{count:03d}_of_154.json"
        )
        chain = {
            "schema_name": "intersubmod.m2_full_ranking_receipt",
            "schema_version": "2.0.0", "run_contract": run_contract,
            "results": copy.deepcopy(results), "all_pass": count == 154,
            "orchestration": {
                "session_identity": production_ranking._session_identity(output_root, session),
                "batch_start_identity": {
                    "path": batch["path"], "sha256": batch["sha256"],
                    "batch_id": batch["batch_id"], "batch_index": batch["batch_index"],
                },
                "previous_chain_head": previous,
                "batch_completion_attestations": [
                    {"task_id": f"{row['dataset']}:{row['chrom']}",
                     **row["orchestration_completion"]}
                    for row in current_results
                ],
                "cumulative_attested_tasks": count,
            },
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": f"{path.name}.sha256", "covers": path.name,
            },
        }
        production_ranking.write_immutable_json_exclusive(path, chain)
        previous = {"path": str(path.resolve()), "sha256": verifier.sha256_path(path)}
        previous_count = count
    return {"root": output_root, "full": chain, "children": children, "parent": parent}


def build_post_input_identity_fixture(root: pathlib.Path) -> dict:
    release_root = root / "post-release"
    input_root = release_root / "input_contract"
    code_root = release_root / "code_snapshot"
    input_root.mkdir(parents=True)
    code_root.mkdir(parents=True)
    canonical = input_root / "canonical_manifest.json"
    canonical.write_text('{"canonical":true}\n', encoding="utf-8")
    canonical.chmod(0o444)
    schema = code_root / "schema.json"
    schema.write_text('{"type":"object"}\n', encoding="utf-8")
    schema.chmod(0o444)
    input_verifier = code_root / "verify_m2_input_identity.py"
    input_verifier.write_text("# frozen verifier\n", encoding="utf-8")
    input_verifier.chmod(0o444)
    artifacts = [
        {
            "dataset": verifier.DATASETS[index % len(verifier.DATASETS)],
            "role": f"role_{index:02d}", "policy": "full_sha256",
            "path": f"/input/{index:02d}", "observed": {"sha256": f"{index:064x}"},
        }
        for index in range(42)
    ]
    snapshot = {
        "manifest_sha256": verifier.sha256_path(canonical),
        "schema_sha256": verifier.sha256_path(schema),
        "datasets": list(verifier.DATASETS), "artifacts": artifacts,
    }
    snapshot_sha = verifier.semantic_json_sha256(snapshot)
    pre_path = input_root / "pre.json"
    pre = {
        "schema_name": verifier.POST_INPUT_IDENTITY_SCHEMA_NAME,
        "schema_version": verifier.POST_INPUT_IDENTITY_SCHEMA_VERSION,
        "task_type": "B_COMPREHENSIVE_VALIDATION", "mode": "PRE",
        "authority_mode": "CANONICAL_V5_FROZEN",
        "validation_evidence_eligible": True,
        "input_identity_snapshot": snapshot,
        "input_identity_snapshot_sha256": snapshot_sha,
        "checks": {"fixture": True}, "all_pass": True,
        "receipt_integrity": {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{pre_path.name}.sha256", "covers": pre_path.name,
        },
    }
    production_extraction.write_immutable_json_exclusive(pre_path, pre)
    release_path = release_root / "m2_run_manifest.json"
    release = {
        "pre_input_identity_receipt": {
            "origin": {},
            "immutable_copy": {
                "path": pre_path.relative_to(release_root).as_posix(),
                "sha256": verifier.sha256_path(pre_path),
            },
            "input_identity_snapshot_sha256": snapshot_sha,
            "receipt_semantic_sha256": verifier.semantic_json_sha256(pre),
            "authority_mode": "CANONICAL_V5_FROZEN",
            "validation_evidence_eligible": True, "artifact_roles": 42,
        },
        "receipt_integrity": {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{release_path.name}.sha256", "covers": release_path.name,
        },
    }
    production_extraction.write_immutable_json_exclusive(release_path, release)
    release_sidecar = release_path.with_name(f"{release_path.name}.sha256")
    binding = {
        "release_manifest": {
            "path": str(release_path.resolve()), "sha256": verifier.sha256_path(release_path),
            "semantic_sha256": verifier.semantic_json_sha256(release),
            "sidecar": {
                "path": str(release_sidecar.resolve()),
                "sha256": verifier.sha256_path(release_sidecar),
            },
        },
        "canonical_input_manifest": _source_identity(canonical),
        "snapshot_sources": {
            "canonical_manifest_schema": _source_identity(schema),
            "input_identity_verifier": _source_identity(input_verifier),
        },
    }
    post_path = root / "post.json"
    post = {
        "schema_name": verifier.POST_INPUT_IDENTITY_SCHEMA_NAME,
        "schema_version": verifier.POST_INPUT_IDENTITY_SCHEMA_VERSION,
        "task_type": "B_COMPREHENSIVE_VALIDATION", "mode": "POST_COMPARE",
        "authority_mode": "CANONICAL_V5_FROZEN",
        "validation_evidence_eligible": True,
        "authority": {
            "canonical_manifest_sha256": verifier.sha256_path(canonical),
            "canonical_schema_sha256": verifier.sha256_path(schema),
            "canonical_schema_path": str(schema.resolve()),
            "selected_authority_is_canonical": True, "test_only_override": False,
        },
        "scope": {
            "technical_datasets": 7, "biological_samples": 6,
            "chromosomes": "chr1-chr22", "tasks": 154,
            "datasets": list(verifier.DATASETS), "direct_input_artifacts": 42,
        },
        "manifest": {
            "path": str(canonical.resolve()), "sha256": verifier.sha256_path(canonical),
            "expected_sha256": verifier.sha256_path(canonical),
        },
        "canonical_schema": {
            "path": str(schema.resolve()), "sha256": verifier.sha256_path(schema),
            "expected_sha256": verifier.sha256_path(schema),
        },
        "verifier": _source_identity(input_verifier),
        "assurance": {
            "bam": "fixture", "other_direct_inputs": "fixture",
            "pre_post": "persistent equality only", "temporal_immutability_proven": False,
        },
        "input_identity_snapshot": snapshot,
        "input_identity_snapshot_sha256": snapshot_sha,
        "artifact_audit": {
            "n_artifacts": 42, "n_unique_resolved_files": 42,
            "n_storage_identity_v1": 7, "n_full_sha256": 35,
            "n_sampled_bam_chunks": 21, "n_mismatches": 0,
        },
        "compare_to": {
            "path": str(pre_path.resolve()), "sha256": verifier.sha256_path(pre_path),
            "pre_snapshot_sha256": snapshot_sha, "exact_snapshot_equal": True,
        },
        "checks": {name: True for name in verifier.POST_INPUT_IDENTITY_CHECK_NAMES},
        "all_pass": True,
        "receipt_integrity": {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{post_path.name}.sha256", "covers": post_path.name,
        },
    }
    production_extraction.write_immutable_json_exclusive(post_path, post)
    return {"path": post_path, "post": post, "binding": binding}


class ReleaseFinalVerificationTests(unittest.TestCase):
    @staticmethod
    def ranking_parameters(frozen):
        return {
            "structural_exact_pattern_minread_grid": frozen["structural_exact_pattern_minread_grid"],
            "primary_structural_exact_pattern_minread": 3, "scoring_minread": 1,
            "exact_k_max": 12, "max_vertex_sets": 256,
            "solver_time_limit_seconds_per_milp": 30.0,
            "minimum_bq_error_rate": 1e-6, "maximum_bq_error_rate": 0.25,
            "fixed_error_grid_conditional_binary_flip_probability": [0.005, 0.01, 0.02, 0.05],
            "conditional_candidate_ranking_bootstrap_replicates": 20,
            "conditional_candidate_ranking_bootstrap_base_seed": 20260716,
            "tie_tolerance_log_likelihood": 1e-6,
        }

    def test_final_cli_requires_release_contract_manifest(self):
        with mock.patch.object(sys, "argv", ["verify_full_m2_receipts.py"]):
            with self.assertRaises(SystemExit) as caught:
                verifier.parse_args()
        self.assertEqual(caught.exception.code, 2)

    def test_verifier_reauthenticates_release_and_full_receipt_bindings(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, paths = make_verifier_release_contract(pathlib.Path(tmp))
            binding = verifier.load_release_contract_binding(
                manifest, required_sources={"full_verifier": paths["full_verifier"]},
                _skip_deep_verification_for_test=True,
            )
            extraction = {
                "run_contract": {
                    "release_binding": binding,
                    "runner": {
                        "path": str(paths["full_extraction_runner"]),
                        "sha256": verifier.sha256_path(paths["full_extraction_runner"]),
                    },
                    "extractor": {
                        "path": str(paths["extractor"]),
                        "sha256": verifier.sha256_path(paths["extractor"]),
                    },
                    "lossless_read_contract": {
                        "path": str(paths["lossless_read_contract"]),
                        "sha256": verifier.sha256_path(paths["lossless_read_contract"]),
                    },
                    "extractor_sha256": verifier.sha256_path(paths["extractor"]),
                    "manifest_sha256": binding["canonical_input_manifest"]["sha256"],
                    "mapq_min": 20, "baseq_min": 20,
                    "bridge_thresholds": [1, 2, 3, 5], "workers": 2,
                    "samtools_threads": 1, "task_timeout_seconds": 28800,
                    "timeout_grace_seconds": 300,
                },
                "invocation": {"max_new_tasks": 16},
            }
            verifier.verify_extraction_release_binding(extraction, binding)
            correct_lossless = dict(extraction["run_contract"]["lossless_read_contract"])
            extraction["run_contract"]["lossless_read_contract"]["sha256"] = "0" * 64
            with self.assertRaisesRegex(verifier.VerificationError, "declared SHA differs"):
                verifier.verify_extraction_release_binding(extraction, binding)
            extraction["run_contract"]["lossless_read_contract"] = correct_lossless
            correct_runner = dict(extraction["run_contract"]["runner"])
            extraction["run_contract"]["runner"] = {
                "path": str(paths["full_ranking_runner"]),
                "sha256": verifier.sha256_path(paths["full_ranking_runner"]),
            }
            with self.assertRaisesRegex(verifier.VerificationError, "source path differs"):
                verifier.verify_extraction_release_binding(extraction, binding)
            extraction["run_contract"]["runner"] = correct_runner
            ranking_frozen = binding["frozen_parameters"]["ranking"]
            ranking = {
                "run_contract": {
                    "release_binding": binding,
                    "runner": {
                        "path": str(paths["full_ranking_runner"]),
                        "sha256": verifier.sha256_path(paths["full_ranking_runner"]),
                    },
                    "ranker": {
                        "path": str(paths["ranker"]),
                        "sha256": verifier.sha256_path(paths["ranker"]),
                    },
                    "hypercube_solver": {
                        "path": str(paths["hypercube_solver"]),
                        "sha256": verifier.sha256_path(paths["hypercube_solver"]),
                    },
                    "parameters": self.ranking_parameters(ranking_frozen),
                    "thresholds": [1, 2, 3, 5], "component_bases": ["PS_HP1", "PS_HP2"],
                    "hp_families": ["1", "2"], "workers": 2,
                    "task_timeout_seconds": 28800, "timeout_grace_seconds": 300,
                },
                "invocation": {"max_new_tasks": 16, "aggregate_only": False},
            }
            verifier.verify_ranking_release_binding(ranking, extraction, binding)
            correct_solver = dict(ranking["run_contract"]["hypercube_solver"])
            ranking["run_contract"]["hypercube_solver"]["sha256"] = "0" * 64
            with self.assertRaisesRegex(verifier.VerificationError, "declared SHA differs"):
                verifier.verify_ranking_release_binding(ranking, extraction, binding)
            ranking["run_contract"]["hypercube_solver"] = correct_solver

            ranking["run_contract"]["parameters"][
                "conditional_candidate_ranking_bootstrap_replicates"
            ] = 19
            with self.assertRaisesRegex(verifier.VerificationError, "scientific parameters"):
                verifier.verify_ranking_release_binding(ranking, extraction, binding)

    def test_verifier_rejects_sidecar_source_and_cross_stage_contract_tamper(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, paths = make_verifier_release_contract(pathlib.Path(tmp))
            sidecar = manifest.with_name(f"{manifest.name}.sha256")
            sidecar.chmod(0o644)
            sidecar.write_text(f"{'0' * 64}  {manifest.name}\n", encoding="ascii")
            sidecar.chmod(0o444)
            with self.assertRaisesRegex(verifier.VerificationError, "sidecar mismatch"):
                verifier.load_release_contract_binding(
                    manifest, required_sources={"full_verifier": paths["full_verifier"]},
                    _skip_deep_verification_for_test=True,
                )

        with tempfile.TemporaryDirectory() as tmp:
            manifest, paths = make_verifier_release_contract(pathlib.Path(tmp))
            binding = verifier.load_release_contract_binding(
                manifest, required_sources={"full_verifier": paths["full_verifier"]},
                _skip_deep_verification_for_test=True,
            )
            extraction = {"run_contract": {"release_binding": dict(binding)}}
            ranking = {"run_contract": {"release_binding": binding}}
            extraction["run_contract"]["release_binding"] = json.loads(json.dumps(binding))
            extraction["run_contract"]["release_binding"]["release_manifest"]["sha256"] = "0" * 64
            with self.assertRaisesRegex(verifier.VerificationError, "contracts differ"):
                verifier.verify_ranking_release_binding(ranking, extraction, binding)

    def test_strict_loader_rejects_before_malicious_freezer_and_schema_resigning(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest, paths = make_verifier_release_contract(root)
            document = json.loads(manifest.read_text(encoding="utf-8"))
            sentinel = root / "MALICIOUS_FREEZER_EXECUTED"
            freezer = paths["release_contract_freezer"]
            freezer.chmod(0o644)
            freezer.write_text(
                f"from pathlib import Path\nPath({str(sentinel)!r}).write_text('executed')\n",
                encoding="utf-8",
            )
            freezer.chmod(0o444)
            document["authority_mode"] = "TEST_ONLY_UNFROZEN"
            _rewrite_immutable_json(manifest, document)
            with self.assertRaisesRegex(verifier.VerificationError, "authority"):
                verifier.load_release_contract_binding(
                    manifest, required_sources={"full_verifier": paths["full_verifier"]}
                )
            self.assertFalse(sentinel.exists())

        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest, paths = make_verifier_release_contract(root)
            document = json.loads(manifest.read_text(encoding="utf-8"))
            sentinel = root / "CORRECT_AUTHORITY_FREEZER_EXECUTED"
            freezer = paths["release_contract_freezer"]
            freezer.chmod(0o644)
            freezer.write_text(
                f"from pathlib import Path\nPath({str(sentinel)!r}).write_text('executed')\n",
                encoding="utf-8",
            )
            freezer.chmod(0o444)
            entry = next(
                row for row in document["source_snapshot"]["entries"]
                if row["role"] == "release_contract_freezer"
            )
            observed = freezer.stat()
            entry["snapshot"] = {
                "path": freezer.relative_to(root).as_posix(),
                "st_dev": int(observed.st_dev), "st_ino": int(observed.st_ino),
                "st_nlink": int(observed.st_nlink), "size_bytes": int(observed.st_size),
                "mtime_ns": int(observed.st_mtime_ns), "ctime_ns": int(observed.st_ctime_ns),
                "mode_octal": "0444", "sha256": verifier.sha256_path(freezer),
            }
            document["producer"]["immutable_copy_sha256"] = entry["snapshot"]["sha256"]
            _rewrite_immutable_json(manifest, document)
            with self.assertRaisesRegex(verifier.VerificationError, "source/copy"):
                verifier.load_release_contract_binding(
                    manifest, required_sources={"full_verifier": paths["full_verifier"]}
                )
            self.assertFalse(sentinel.exists())

        mutations = (
            ("missing top key", lambda doc: doc.pop("producer"), "exact-key"),
            (
                "extra role",
                lambda doc: (
                    doc["source_snapshot"]["entries"].append(
                        copy.deepcopy(doc["source_snapshot"]["entries"][0])
                    ),
                    doc["source_snapshot"].__setitem__("n_files", 12),
                ),
                "eleven",
            ),
            (
                "parameter drift",
                lambda doc: doc["parameters"]["ranking"].__setitem__("exact_k_max", 99),
                "parameters",
            ),
            (
                "PRE escape",
                lambda doc: doc["pre_input_identity_receipt"]["immutable_copy"].__setitem__(
                    "path", "../escaped-pre.json"
                ),
                "escapes",
            ),
        )
        for label, mutate, message in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                manifest, paths = make_verifier_release_contract(pathlib.Path(tmp))
                document = json.loads(manifest.read_text(encoding="utf-8"))
                mutate(document)
                _rewrite_immutable_json(manifest, document)
                with self.assertRaisesRegex(verifier.VerificationError, message):
                    verifier.load_release_contract_binding(
                        manifest,
                        required_sources={"full_verifier": paths["full_verifier"]},
                        _skip_deep_verification_for_test=True,
                    )

    def test_exclusive_final_writer_never_overwrites_and_freezes_both_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "final.json"
            receipt = {
                "all_pass": True,
                "receipt_integrity": {
                    "scheme": "external_sha256_sidecar_v1",
                    "sidecar_name": "final.json.sha256", "covers": "final.json",
                },
            }
            verifier.write_immutable_receipt_exclusive(path, receipt)
            sidecar = path.with_name("final.json.sha256")
            self.assertEqual(path.stat().st_mode & 0o777, 0o444)
            self.assertEqual(sidecar.stat().st_mode & 0o777, 0o444)
            original = path.read_bytes()
            with self.assertRaisesRegex(verifier.VerificationError, "already exists"):
                verifier.write_immutable_receipt_exclusive(path, {"all_pass": False})
            self.assertEqual(path.read_bytes(), original)

        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "final.json"
            sidecar = path.with_name("final.json.sha256")
            sidecar.write_text("preexisting\n", encoding="ascii")
            with self.assertRaisesRegex(verifier.VerificationError, "already exists"):
                verifier.write_immutable_receipt_exclusive(path, {"all_pass": True})
            self.assertFalse(path.exists())
            self.assertEqual(sidecar.read_text(encoding="ascii"), "preexisting\n")

    def test_mid_verification_release_or_post_drift_blocks_publication_boundary(self):
        orchestration = {
            "n_attested_children": 154, "session_id": "1" * 64,
            "session_receipt": {"path": "/session", "sha256": "2" * 64},
            "terminal_receipt": {"path": "/terminal", "sha256": "3" * 64},
        }
        extraction = {
            "n_tasks": 154, "children": {}, "full_receipt": {},
            "orchestration": orchestration,
        }
        ranking = {
            "n_tasks": 154, "all_aggregate_cells_conserved": True,
            "candidate_table": {
                "n_units": 1, "n_rows": 1,
                "all_rows_match_independent_child_reconstruction": True,
                "winner_partitions_conserved": True,
            },
            "runtime_diagnostics": {
                "all_child_and_full_runtime_summaries_independently_recomputed": True,
            },
            "method_contract_verification": {
                "all_children_exactly_matched_and_source_bound": True,
            },
            "profile_likelihood_independent_recomputation": {
                "n_units": 1, "n_candidates": 1,
                "all_profile_likelihoods_and_certificates_match": True,
                "all_relative_weights_match": True,
                "all_winner_tie_partitions_match": True,
            },
            "orchestration": orchestration,
        }
        start_binding = {"binding": "start"}
        drifted_binding = {"binding": "drifted"}
        post = {"exact_snapshot_equal": True}
        with (
            mock.patch.object(
                verifier, "load_release_contract_binding",
                side_effect=[start_binding, drifted_binding],
            ) as load_binding,
            mock.patch.object(verifier, "verify_post_input_identity_receipt", return_value=post),
            mock.patch.object(verifier, "verify_extraction_bundle", return_value=extraction),
            mock.patch.object(verifier, "verify_ranking_bundle", return_value=ranking),
        ):
            with self.assertRaisesRegex(verifier.VerificationError, "drifted"):
                verifier.verify_bundle(
                    pathlib.Path("/extraction"), pathlib.Path("/ranking"),
                    release_contract_manifest=pathlib.Path("/release"),
                    post_input_identity_receipt=pathlib.Path("/post"),
                )
            self.assertTrue(load_binding.call_args_list[-1].kwargs["_force_deep_reverification"])

        with (
            mock.patch.object(
                verifier, "load_release_contract_binding",
                side_effect=[start_binding, start_binding],
            ),
            mock.patch.object(
                verifier, "verify_post_input_identity_receipt",
                side_effect=[post, {"exact_snapshot_equal": True, "drift": True}],
            ),
            mock.patch.object(verifier, "verify_extraction_bundle", return_value=extraction),
            mock.patch.object(verifier, "verify_ranking_bundle", return_value=ranking),
        ):
            with self.assertRaisesRegex(verifier.VerificationError, "POST"):
                verifier.verify_bundle(
                    pathlib.Path("/extraction"), pathlib.Path("/ranking"),
                    release_contract_manifest=pathlib.Path("/release"),
                    post_input_identity_receipt=pathlib.Path("/post"),
                )


class OrchestrationCrossModuleGoldenTests(unittest.TestCase):
    def test_production_extraction_helpers_generate_verifier_accepted_154_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = build_production_extraction_orchestration_fixture(pathlib.Path(tmp))
            audit = verifier.verify_orchestration_stage(
                fixture["root"], fixture["full"], stage="extraction",
                datasets=verifier.DATASETS, chromosomes=verifier.AUTOSOMES,
                release_binding=fixture["binding"],
            )
            self.assertEqual(audit["n_attested_children"], 154)
            self.assertEqual(audit["n_authenticated_resource_gates"], 12)
            self.assertEqual(
                audit["legal_cumulative_counts"],
                [8, 24, 40, 56, 72, 88, 104, 120, 136, 152, 154],
            )
            final_batch = json.loads((
                fixture["root"] / "_orchestration/batches/batch_011_start.json"
            ).read_text(encoding="utf-8"))
            self.assertEqual(
                (final_batch["before_count"], final_batch["max_new_tasks"],
                 final_batch["effective_count"]),
                (152, 16, 2),
            )
            session_path = fixture["root"] / "_orchestration/session.json"
            session = json.loads(session_path.read_text(encoding="utf-8"))
            drifted_session = copy.deepcopy(session)
            drifted_session["created_at_utc"] = "2099-01-01T00:00:00Z"
            _rewrite_immutable_json(session_path, drifted_session)
            with self.assertRaisesRegex(verifier.VerificationError, "drifted"):
                verifier._reauthenticate_orchestration_publication_boundary(
                    audit, "extraction"
                )
            _rewrite_immutable_json(session_path, session)

            terminal_path = fixture["root"] / "full_extraction_receipt.json"
            terminal = json.loads(terminal_path.read_text(encoding="utf-8"))
            drifted_terminal = copy.deepcopy(terminal)
            drifted_terminal["all_pass"] = False
            _rewrite_immutable_json(terminal_path, drifted_terminal)
            with self.assertRaisesRegex(verifier.VerificationError, "drifted"):
                verifier._reauthenticate_orchestration_publication_boundary(
                    audit, "extraction"
                )
            _rewrite_immutable_json(terminal_path, terminal)

    def test_resigned_tamper_or_orphan_symlink_hardlink_and_root_swap_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = build_production_extraction_orchestration_fixture(pathlib.Path(tmp))
            root = fixture["root"]

            def verify_now():
                return verifier.verify_orchestration_stage(
                    root, fixture["full"], stage="extraction",
                    datasets=verifier.DATASETS, chromosomes=verifier.AUTOSOMES,
                    release_binding=fixture["binding"],
                )

            grant_path = root / "_orchestration/tasks/COLO829/chr1/grant.json"
            original_grant = json.loads(grant_path.read_text(encoding="utf-8"))
            tampered_grant = copy.deepcopy(original_grant)
            tampered_grant["parameters_semantic_sha256"] = "0" * 64
            _rewrite_immutable_json(grant_path, tampered_grant)
            with self.assertRaises(verifier.VerificationError):
                verify_now()
            _rewrite_immutable_json(grant_path, original_grant)

            completion_path = root / "_orchestration/tasks/COLO829/chr1/completion.json"
            original_completion = json.loads(completion_path.read_text(encoding="utf-8"))
            tampered_completion = copy.deepcopy(original_completion)
            tampered_completion["child_receipt_identity"]["sha256"] = "0" * 64
            _rewrite_immutable_json(completion_path, tampered_completion)
            with self.assertRaises(verifier.VerificationError):
                verify_now()
            _rewrite_immutable_json(completion_path, original_completion)

            gate_path = root / "_orchestration/resource_gates/batch_001.json"
            original_gate = json.loads(gate_path.read_text(encoding="utf-8"))
            tampered_gate = copy.deepcopy(original_gate)
            tampered_gate["process_snapshot"]["semantic_sha256"] = "f" * 64
            _rewrite_immutable_json(gate_path, tampered_gate)
            with self.assertRaises(verifier.VerificationError):
                verify_now()
            _rewrite_immutable_json(gate_path, original_gate)

            checkpoint_path = root / "checkpoints/checkpoint_024_of_154.json"
            original_checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            tampered_checkpoint = copy.deepcopy(original_checkpoint)
            tampered_checkpoint["orchestration"]["previous_chain_head"]["sha256"] = "0" * 64
            _rewrite_immutable_json(checkpoint_path, tampered_checkpoint)
            with self.assertRaisesRegex(verifier.VerificationError, "chain"):
                verify_now()
            _rewrite_immutable_json(checkpoint_path, original_checkpoint)

            terminal_path = root / "full_extraction_receipt.json"
            original_terminal = json.loads(terminal_path.read_text(encoding="utf-8"))
            tampered_terminal = copy.deepcopy(original_terminal)
            tampered_terminal["orchestration"]["session_identity"]["sha256"] = "0" * 64
            _rewrite_immutable_json(terminal_path, tampered_terminal)
            with self.assertRaisesRegex(verifier.VerificationError, "orchestration"):
                verify_now()
            _rewrite_immutable_json(terminal_path, original_terminal)

            marker = root / "samples/COLO829/chr1/runner_task_status.json"
            marker.write_text('{"timed_out":true}\n', encoding="utf-8")
            with self.assertRaisesRegex(verifier.VerificationError, "marker"):
                verify_now()
            marker.unlink()

            orphan = root / "_orchestration/batches/batch_999_start.json"
            orphan.write_text("{}\n", encoding="utf-8")
            with self.assertRaises(verifier.VerificationError):
                verify_now()
            orphan.unlink()

            child_path = root / "samples/COLO829/chr1/receipt.json"
            child_sidecar = child_path.with_name(f"{child_path.name}.sha256")
            saved_child = child_path.read_bytes()
            saved_sidecar = child_sidecar.read_bytes()
            child_path.unlink()
            child_path.symlink_to(root / "samples/COLO829/chr2/receipt.json")
            with self.assertRaisesRegex(verifier.VerificationError, "non-symlink"):
                verify_now()
            child_path.unlink()
            child_path.write_bytes(saved_child)
            child_path.chmod(0o444)
            child_sidecar.chmod(0o644)
            child_sidecar.write_bytes(saved_sidecar)
            child_sidecar.chmod(0o444)

            alias = root / "samples/COLO829/chr1/receipt.alias.json"
            os.link(child_path, alias)
            with self.assertRaisesRegex(verifier.VerificationError, "single-link"):
                verify_now()
            alias.unlink()

            output_path = root / "samples/COLO829/chr1/tiny.tsv.gz"
            output_bytes = output_path.read_bytes()
            output_path.unlink()
            output_path.symlink_to(root / "samples/COLO829/chr2/tiny.tsv.gz")
            with self.assertRaisesRegex(verifier.VerificationError, "non-symlink"):
                verify_now()
            output_path.unlink()
            output_path.write_bytes(output_bytes)
            output_alias = output_path.with_name("tiny.alias.tsv.gz")
            os.link(output_path, output_alias)
            with self.assertRaisesRegex(verifier.VerificationError, "single-link"):
                verify_now()
            output_alias.unlink()

            attestation_ancestor = root / "_orchestration/tasks/COLO829"
            moved_ancestor = pathlib.Path(tmp) / "COLO829-attestation-real"
            attestation_ancestor.rename(moved_ancestor)
            attestation_ancestor.symlink_to(moved_ancestor, target_is_directory=True)
            with self.assertRaisesRegex(verifier.VerificationError, "ancestor"):
                verify_now()
            attestation_ancestor.unlink()
            moved_ancestor.rename(attestation_ancestor)

            backup = pathlib.Path(tmp) / "extraction-original"
            root.rename(backup)
            shutil.copytree(backup, root, copy_function=shutil.copy2)
            with self.assertRaisesRegex(verifier.VerificationError, "output-root"):
                verify_now()
            shutil.rmtree(root)
            backup.rename(root)

    def test_production_ranking_helpers_generate_parent_bound_154_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            extraction = build_production_extraction_orchestration_fixture(root)
            extraction_audit = verifier.verify_orchestration_stage(
                extraction["root"], extraction["full"], stage="extraction",
                datasets=verifier.DATASETS, chromosomes=verifier.AUTOSOMES,
                release_binding=extraction["binding"],
            )
            ranking = build_production_ranking_orchestration_fixture(root, extraction)
            ranking_audit = verifier.verify_orchestration_stage(
                ranking["root"], ranking["full"], stage="ranking",
                datasets=verifier.DATASETS, chromosomes=verifier.AUTOSOMES,
                release_binding=extraction["binding"], extraction_root=extraction["root"],
                extraction_children=extraction["children"],
                parent_extraction=extraction_audit,
            )
            self.assertEqual(ranking_audit["n_attested_children"], 154)
            self.assertEqual(ranking_audit["n_authenticated_resource_gates"], 12)
            self.assertTrue(ranking_audit["all_grants_and_completions_session_bound"])

            session = json.loads((
                ranking["root"] / "_orchestration/session.json"
            ).read_text(encoding="utf-8"))
            batch1 = json.loads((
                ranking["root"] / "_orchestration/batches/batch_001_start.json"
            ).read_text(encoding="utf-8"))
            batch2 = json.loads((
                ranking["root"] / "_orchestration/batches/batch_002_start.json"
            ).read_text(encoding="utf-8"))
            gate1_path = ranking["root"] / "_orchestration/resource_gates/batch_001.json"
            gate1 = json.loads(gate1_path.read_text(encoding="utf-8"))
            with self.assertRaisesRegex(verifier.VerificationError, "path swap"):
                verifier._verify_resource_gate(
                    batch2["resource_gate"], root=ranking["root"], stage="ranking",
                    gate_scope="batch", batch_index=1,
                    expected_target=gate1["target"],
                    expected_producer_source=session["producer_sources"]["runner"],
                    label="ranking reused/path-swapped gate",
                )
            with self.assertRaises(verifier.VerificationError):
                verifier._verify_resource_gate(
                    None, root=ranking["root"], stage="ranking",
                    gate_scope="batch", batch_index=1,
                    expected_target=gate1["target"],
                    expected_producer_source=session["producer_sources"]["runner"],
                    label="ranking null gate",
                )
            with self.assertRaisesRegex(verifier.VerificationError, "producer source"):
                verifier._verify_resource_gate(
                    batch1["resource_gate"], root=ranking["root"], stage="ranking",
                    gate_scope="batch", batch_index=1,
                    expected_target=gate1["target"],
                    expected_producer_source={"path": "/swapped", "sha256": "0" * 64},
                    label="ranking producer-swapped gate",
                )
            tampered_gate = copy.deepcopy(gate1)
            filesystem = tampered_gate["filesystem_snapshot"]
            filesystem.update({
                "f_bavail": verifier.RESOURCE_GATE_REQUIRED_RESERVE_BYTES - 1,
                "f_frsize": 1,
                "available_bytes": verifier.RESOURCE_GATE_REQUIRED_RESERVE_BYTES - 1,
                "disk_pass": True,
            })
            filesystem["semantic_sha256"] = verifier.semantic_json_sha256({
                key: filesystem[key]
                for key in (
                    "probe_path", "target_output_root", "st_dev", "f_bavail", "f_frsize",
                    "available_bytes", "required_reserve_bytes", "disk_pass",
                )
            })
            tampered_gate["gate_id"] = verifier.semantic_json_sha256({
                key: value for key, value in tampered_gate.items()
                if key not in {"gate_id", "receipt_integrity"}
            })
            _rewrite_immutable_json(gate1_path, tampered_gate)
            gate1_sidecar = gate1_path.with_name(f"{gate1_path.name}.sha256")
            tampered_identity = {
                "path": str(gate1_path.resolve()),
                "sha256": verifier.sha256_path(gate1_path),
                "semantic_sha256": verifier.semantic_json_sha256(tampered_gate),
                "gate_id": tampered_gate["gate_id"],
                "sidecar_path": str(gate1_sidecar.resolve()),
                "sidecar_sha256": verifier.sha256_path(gate1_sidecar),
            }
            with self.assertRaisesRegex(verifier.VerificationError, "300 GiB"):
                verifier._verify_resource_gate(
                    tampered_identity, root=ranking["root"], stage="ranking",
                    gate_scope="batch", batch_index=1,
                    expected_target=gate1["target"],
                    expected_producer_source=session["producer_sources"]["runner"],
                    label="ranking low-disk gate",
                )
            _rewrite_immutable_json(gate1_path, gate1)

            grant_path = (
                ranking["root"] /
                "_orchestration/tasks/COLO829/chr1/grant.json"
            )
            grant = json.loads(grant_path.read_text(encoding="utf-8"))
            grant["input_identity"]["extraction_child_receipt_sha256"] = "0" * 64
            _rewrite_immutable_json(grant_path, grant)
            with self.assertRaises(verifier.VerificationError):
                verifier.verify_orchestration_stage(
                    ranking["root"], ranking["full"], stage="ranking",
                    datasets=verifier.DATASETS, chromosomes=verifier.AUTOSOMES,
                    release_binding=extraction["binding"],
                    extraction_root=extraction["root"],
                    extraction_children=extraction["children"],
                    parent_extraction=extraction_audit,
                )

        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            extraction = build_production_extraction_orchestration_fixture(root)
            extraction_audit = verifier.verify_orchestration_stage(
                extraction["root"], extraction["full"], stage="extraction",
                datasets=verifier.DATASETS, chromosomes=verifier.AUTOSOMES,
                release_binding=extraction["binding"],
            )
            ranking = build_production_ranking_orchestration_fixture(root, extraction)
            wrong_parent = copy.deepcopy(extraction_audit)
            wrong_parent["session_id"] = "0" * 64
            with self.assertRaisesRegex(verifier.VerificationError, "not bound"):
                verifier.verify_orchestration_stage(
                    ranking["root"], ranking["full"], stage="ranking",
                    datasets=verifier.DATASETS, chromosomes=verifier.AUTOSOMES,
                    release_binding=extraction["binding"],
                    extraction_root=extraction["root"],
                    extraction_children=extraction["children"],
                    parent_extraction=wrong_parent,
                )


class PostInputIdentityPublicationGateTests(unittest.TestCase):
    def test_post_exactly_matches_immutable_pre_and_frozen_verifier(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture = build_post_input_identity_fixture(pathlib.Path(tmp))
            audit = verifier.verify_post_input_identity_receipt(
                fixture["path"], fixture["binding"]
            )
            self.assertTrue(audit["exact_snapshot_equal"])
            self.assertEqual(audit["n_artifact_roles"], 42)
            self.assertIn("transient mutation", audit["claim_boundary"])

    def test_post_resigned_snapshot_verifier_check_and_sidecar_tamper_fail(self):
        mutations = (
            (
                "snapshot",
                lambda post: post["input_identity_snapshot"]["artifacts"][0].__setitem__(
                    "path", "/changed"
                ),
                "not exactly equal",
            ),
            (
                "verifier",
                lambda post: post["verifier"].__setitem__("sha256", "0" * 64),
                "verifier",
            ),
            (
                "check",
                lambda post: post["checks"].__setitem__(
                    "post_snapshot_exactly_matches_authenticated_pre", False
                ),
                "named check",
            ),
            (
                "compare",
                lambda post: post["compare_to"].__setitem__("exact_snapshot_equal", False),
                "compare_to",
            ),
        )
        for label, mutate, message in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                fixture = build_post_input_identity_fixture(pathlib.Path(tmp))
                post = copy.deepcopy(fixture["post"])
                mutate(post)
                _rewrite_immutable_json(fixture["path"], post)
                with self.assertRaisesRegex(verifier.VerificationError, message):
                    verifier.verify_post_input_identity_receipt(
                        fixture["path"], fixture["binding"]
                    )

        with tempfile.TemporaryDirectory() as tmp:
            fixture = build_post_input_identity_fixture(pathlib.Path(tmp))
            sidecar = fixture["path"].with_name(f"{fixture['path'].name}.sha256")
            sidecar.chmod(0o644)
            sidecar.write_text(f"{'0' * 64}  {fixture['path'].name}\n", encoding="ascii")
            sidecar.chmod(0o444)
            with self.assertRaisesRegex(verifier.VerificationError, "sidecar mismatch"):
                verifier.verify_post_input_identity_receipt(
                    fixture["path"], fixture["binding"]
                )


def component_cell(maximum: int) -> dict:
    distribution = {"1": 1, str(maximum): 1} if maximum != 1 else {"1": 2}
    return {
        "n_components": 2,
        "n_singletons_k1": distribution.get("1", 0),
        "n_multisite_k_gt1": sum(value for key, value in distribution.items() if int(key) > 1),
        "max_k_component_sites": maximum,
        "max_k": maximum,
        "k_component_sites_distribution": distribution,
        "k_distribution": distribution,
    }


def extraction_child(maximum: int) -> dict:
    return {
        "counts": {
            "raw_overlapping_alignments": 10,
            "alignment_class_primary": 10,
            "alignment_class_secondary": 0,
            "alignment_class_supplementary": 0,
            "alignment_class_unmapped": 0,
            "excluded_by_flag": 1,
            "mapq_rejected_after_flag": 1,
            "canonical_eligible_alignments": 8,
            "molecule_sparse_rows_written": 8,
            "unique_molecule_ids": 8,
            "sidecar_exact_matches": 8,
            "sidecar_missing": 0,
            "site_call_rows_sparse": 20,
            "fixed_ra_calls": 15,
            "alt_calls": 5,
        },
        "phase_set_contract_counts": {
            "known_phase_sets_by_hp_family": {"1": 2, "2": 1},
            "known_ps_active_site_memberships_by_hp_family": {"1": 3, "2": 2},
            "missing_ps_active_sites_by_hp_family": {"1": 0, "2": 1},
        },
        "legacy_cross_phase_set_aggregation_audit": {
            "1": {"1": {"legacy_components": 2}}
        },
        "component_summary_by_linkage_basis": {
            basis: {"1": component_cell(maximum)} for basis in verifier.LINKAGE_BASES
        },
    }


def valid_rank_summary(n_units: int = 1) -> dict:
    return {
        "n_component_hp_units": n_units,
        "n_components": n_units,
        "molecule_component_projections": 5 * n_units,
        "informative_scoring_molecules": 4 * n_units,
        "all_x_excluded_molecules": n_units,
        "structural_retained_molecules": 3 * n_units,
        "below_minread_scoring_molecules": n_units,
        "bq_scoring_molecules": 4 * n_units,
        "raw_tree_candidates_T_complete_units": 2 * n_units,
        "distinct_vertex_sets_V_complete_units": n_units,
        "solver_complete_units": n_units,
        "solver_incomplete_or_not_run_units": 0,
        "quality_primary_unique_vertex_units": n_units,
        "quality_primary_tied_vertex_units": 0,
        "rank_abstain_units": 0,
        "coarse_topology_class_unique_units": n_units,
        "coarse_topology_multiple_class_units": 0,
        "topology_evaluated_units": n_units,
        "k_component_sites_total": 2 * n_units,
        "k_observed_alt_active_total": n_units,
        "not_structural_alt_active_sites_total": n_units,
        "partial_group_coverage_denominator": n_units,
        "partial_groups_covered": n_units,
        "partial_groups_unsatisfied": 0,
        "selection_status_counts": {"T1_CANDIDATE_UNIQUE": n_units},
        "candidate_generation_status_counts": {"EXACT_COMPLETE": n_units},
        "k_route_counts": {"EXACT_COMPLETE": n_units},
        "projected_call_class_counts": {"R_or_A": 4 * n_units},
        "conditional_candidate_ranking_bootstrap_status_counts": {"NOT_RUN": n_units},
        "topology_class_inclusion_counts": {"single-only": n_units},
        "coarse_topology_unique_class_counts": {"single-only": n_units},
        "coarse_topology_ambiguous_class_set_counts": {},
        "topology_derivation_status_counts": {"COMPLETE": n_units},
        "exact_topology_uniqueness_status_counts": {"UNIQUE": n_units},
    }


class IndependenceAndScopeTests(unittest.TestCase):
    def test_verifier_does_not_import_production_aggregators_or_ranker(self):
        source = inspect.getsource(verifier)
        self.assertNotIn("from run_full_m2_extraction import", source)
        self.assertNotIn("from run_full_m2_ranking import", source)
        self.assertNotIn("from build_m2_patterns_and_rank import", source)
        self.assertNotIn("from hypercube_exact import", source)

    def test_duplicate_task_is_rejected_even_if_count_matches(self):
        rows = [
            {"dataset": "D", "chrom": "chr1"},
            {"dataset": "D", "chrom": "chr1"},
        ]
        with self.assertRaises(verifier.VerificationError):
            verifier.validate_result_index(rows, ("D",), ("chr1", "chr2"), "fixture")


class IndependentMethodContractTests(unittest.TestCase):
    def fixture(self, root: pathlib.Path) -> tuple[pathlib.Path, dict]:
        ranker = root / "ranker.py"
        ranker.write_text("# independent verifier fixture\n", encoding="utf-8")
        child = {
            "parameters": {
                "method_contract": copy.deepcopy(verifier.EXPECTED_METHOD_CONTRACT)
            },
            "provenance": {
                "ranker": {
                    "path": str(ranker.resolve()),
                    "sha256": verifier.sha256_path(ranker),
                }
            },
        }
        return ranker, child

    def verify(self, ranker: pathlib.Path, child: dict) -> None:
        verifier.verify_child_method_contract_and_ranker_source(
            child, ranker, verifier.sha256_path(ranker), "fixture"
        )

    def test_exact_contract_and_actual_source_bytes_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            ranker, child = self.fixture(pathlib.Path(tmp))
            self.verify(ranker, child)

    def test_old_key_only_and_contract_drift_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            ranker, child = self.fixture(pathlib.Path(tmp))
            child["parameters"].pop("method_contract")
            child.setdefault("checks", {})[
                "no_" + "partial_completions_materialized"
            ] = True
            with self.assertRaisesRegex(verifier.VerificationError, "method contract"):
                self.verify(ranker, child)
        with tempfile.TemporaryDirectory() as tmp:
            ranker, child = self.fixture(pathlib.Path(tmp))
            child["parameters"]["method_contract"]["structural_group_scope"] = "DRIFT"
            with self.assertRaisesRegex(verifier.VerificationError, "method contract"):
                self.verify(ranker, child)

    def test_vaf_or_edge_true_is_rejected(self):
        for field in (
            "same_molecule_vaf_added_as_separate_term",
            "parent_edges_or_transitions_scored",
        ):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmp:
                ranker, child = self.fixture(pathlib.Path(tmp))
                child["parameters"]["method_contract"][field] = True
                with self.assertRaisesRegex(verifier.VerificationError, "method contract"):
                    self.verify(ranker, child)

    def test_source_mutation_is_rejected_even_if_receipt_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            ranker, child = self.fixture(pathlib.Path(tmp))
            expected_sha = verifier.sha256_path(ranker)
            ranker.write_text("# mutated after receipt\n", encoding="utf-8")
            with self.assertRaisesRegex(verifier.VerificationError, "actual ranker source bytes"):
                verifier.verify_child_method_contract_and_ranker_source(
                    child, ranker, expected_sha, "fixture"
                )


class IndependentExtractionAggregateTests(unittest.TestCase):
    def test_counts_sum_while_max_k_and_distributions_are_recomputed(self):
        children = {
            ("D", "chr1"): extraction_child(5),
            ("D", "chr2"): extraction_child(7),
        }
        result = verifier.aggregate_extraction_children(
            children,
            ("D",),
            (1,),
            {("D", "chr1"): "PASS", ("D", "chr2"): "REUSED_PASS"},
        )
        self.assertEqual(result["counts"]["raw_overlapping_alignments"], 20)
        cell = result["component_summary_by_linkage_basis"]["PS_HP1"]["1"]
        self.assertEqual(cell["n_components"], 4)
        self.assertEqual(cell["max_k_component_sites"], 7)
        self.assertEqual(cell["k_component_sites_distribution"], {"1": 2, "5": 1, "7": 1})
        self.assertEqual(result["by_dataset"]["D"]["task_status_counts"], {"PASS": 1, "REUSED_PASS": 1})

    def test_child_funnel_violation_is_detected_independently(self):
        child = extraction_child(3)
        child["counts"]["canonical_eligible_alignments"] = 9
        with self.assertRaises(verifier.VerificationError):
            verifier.extraction_count_checks(child["counts"], "tampered")


class IndependentRankingConservationTests(unittest.TestCase):
    def test_summary_is_recomputed_and_conserves_units(self):
        cell = verifier.new_rank_cell()
        verifier.add_rank_summary(cell, valid_rank_summary(2), {"molecule_projections": {"denominator": 10}})
        frozen = verifier.freeze_rank_cell(cell)
        checks = verifier.rank_conservation(frozen)
        self.assertTrue(all(checks.values()), checks)
        self.assertEqual(frozen["n_component_hp_units"], 2)

    def test_forged_full_aggregate_and_bad_unit_partition_are_rejected(self):
        cell = verifier.new_rank_cell()
        verifier.add_rank_summary(cell, valid_rank_summary(), None)
        frozen = verifier.freeze_rank_cell(cell)
        forged = dict(frozen)
        forged["solver_complete_units"] = 2
        with self.assertRaises(verifier.VerificationError):
            verifier.compare_rank_cell(forged, frozen, "forged")
        bad = dict(frozen)
        bad["selection_status_counts"] = {}
        self.assertFalse(verifier.rank_conservation(bad)["selection_status_sum_equals_units"])


class IndependentRuntimeDiagnosticsTests(unittest.TestCase):
    def test_nearest_rank_and_child_rows_are_independently_recomputed(self):
        self.assertEqual(
            verifier.independent_runtime_summary([4.0, 1.0, 5.0, 2.0, 3.0]),
            {"n": 5, "sum": 15.0, "max": 5.0, "p50": 3.0, "p95": 5.0, "p99": 5.0},
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / verifier.RUNTIME_SOURCE_NAME
            rows = [
                {
                    "dataset": "D", "chrom": "chr1", "threshold": 3,
                    "component_basis": "PS_HP1", "phase_set": "100", "component_id": "C",
                    "family": "1", "structural_exact_pattern_minread": minread,
                    "structural_minread_role": role,
                    "candidate_generation_invoked": "true",
                    "likelihood_fit_invoked": "true",
                    "candidate_generation_elapsed_seconds": candidate,
                    "likelihood_fit_elapsed_seconds": likelihood,
                    "unit_total_elapsed_seconds": total,
                }
                for minread, role, candidate, likelihood, total in (
                    (3, "PRIMARY", 0.1, 0.2, 0.5),
                    (5, "SENSITIVITY", 0.3, 0.4, 0.9),
                )
            ]
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, verifier.RUNTIME_DIAGNOSTIC_FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            def scope(selected):
                return {
                    "n_unit_evaluations": len(selected),
                    **{
                        metric: verifier.independent_runtime_summary(
                            float(row[metric]) for row in selected
                        )
                        for metric in verifier.RUNTIME_METRICS
                    },
                }

            child = {
                "outputs": {verifier.RUNTIME_SOURCE_NAME: {"path": str(path)}},
                "runtime_diagnostics": {
                    "schema_name": "intersubmod.m2_unit_runtime_diagnostics",
                    "schema_version": "1.0.0",
                    "clock": "time.monotonic_ns",
                    "unit": "seconds",
                    "per_unit_output": verifier.RUNTIME_SOURCE_NAME,
                    "scopes": {
                        "primary_unit_evaluations": scope(rows[:1]),
                        "all_structural_minread_unit_evaluations": scope(rows),
                    },
                    "primary_invoked_segment_scopes": {
                        metric: scope(rows[:1])[metric]
                        for metric in (
                            "candidate_generation_elapsed_seconds",
                            "likelihood_fit_elapsed_seconds",
                        )
                    },
                },
            }
            values = verifier.verify_child_runtime_diagnostics(child, "D", "chr1")
            self.assertEqual(values["all_primary"]["candidate_generation_elapsed_seconds"], [0.1])
            self.assertEqual(values["when_invoked"]["likelihood_fit_elapsed_seconds"], [0.2])
            child["runtime_diagnostics"]["scopes"]["primary_unit_evaluations"][
                "unit_total_elapsed_seconds"
            ]["p95"] = 99.0
            with self.assertRaises(verifier.VerificationError):
                verifier.verify_child_runtime_diagnostics(child, "D", "chr1")


class IndependentProfileLikelihoodTests(unittest.TestCase):
    PARAMETERS = {
        "minimum_bq_error_rate": 1e-6,
        "maximum_bq_error_rate": 0.25,
        "tie_tolerance_log_likelihood": 1e-6,
        "primary_structural_exact_pattern_minread": 3,
        "max_vertex_sets": 256,
    }

    @staticmethod
    def vertex_id(k: int, states: list[int]) -> str:
        return hashlib.sha256(json.dumps(
            {"k": k, "vertices": sorted(states)},
            sort_keys=True, separators=(",", ":"),
        ).encode()).hexdigest()

    @staticmethod
    def state_records(k: int, states: list[int]) -> list[dict]:
        return [
            {
                "bitmask": state,
                "state_rax": "".join(
                    "A" if state & (1 << bit) else "R" for bit in range(k)
                ),
                "roles": ["root"] if state == 0 else ["full-observed"],
            }
            for state in states
        ]

    @staticmethod
    def manual_fit(patterns, states, weights):
        emissions_by_pattern = []
        counts = []
        for pattern, qualities, count in patterns:
            emissions = []
            for state in states:
                probability = 1.0
                for bit, (symbol, quality) in enumerate(zip(pattern, qualities)):
                    if symbol == "X":
                        continue
                    raw_error = 10 ** (-quality / 10)
                    error = min(max(raw_error, 1e-6), 0.25)
                    denominator = 1 - 2 * error / 3
                    same = (1 - error) / denominator
                    changed = (error / 3) / denominator
                    probability *= same if bool(state & (1 << bit)) == (symbol == "A") else changed
                emissions.append(probability)
            emissions_by_pattern.append(emissions)
            counts.append(count)
        denominators = [
            math.fsum(weight * value for weight, value in zip(weights, emissions))
            for emissions in emissions_by_pattern
        ]
        ll = math.fsum(count * math.log(value) for count, value in zip(counts, denominators))
        gradients = [
            math.fsum(
                count * emissions[index] / denominator
                for emissions, count, denominator in zip(emissions_by_pattern, counts, denominators)
            )
            for index in range(len(states))
        ]
        gap = max(0.0, max(gradients) - math.fsum(
            weight * gradient for weight, gradient in zip(weights, gradients)
        ))
        simplex = max(abs(math.fsum(weights) - 1), max(0.0, -min(weights)))
        return ll, gap, simplex

    def make_fixture(self, root: pathlib.Path):
        pattern_path = root / verifier.PATTERN_SOURCE_NAME
        candidate_path = root / verifier.CANDIDATE_SOURCE_NAME
        common = {
            "dataset": "D", "chrom": "chr1", "threshold": "3",
            "component_basis": "PS_HP1", "phase_set": "100",
            "component_id": "C", "family": "1",
            "structural_exact_pattern_minread": "3",
        }
        pattern_rows = [
            {
                **common, "structural_minread_role": "PRIMARY", "k": "2",
                "pattern_rax": pattern, "fixed_base_qualities": "30,",
                "n_molecules": str(count), "scoring_eligible": "true",
            }
            for pattern, count in (("RX", 40), ("AX", 60))
        ]
        with gzip.open(pattern_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, verifier.PROFILE_PATTERN_REQUIRED_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(pattern_rows)

        error = 10 ** (-30 / 10)
        denominator = 1 - 2 * error / 3
        match = (1 - error) / denominator
        flip = (error / 3) / denominator
        alt_weight = (0.60 - flip) / (match - flip)
        specifications = [
            ([0, 1], [1 - alt_weight, alt_weight]),
            ([0], [1.0]),
        ]
        manual_patterns = [("RX", (30, -1), 40), ("AX", (30, -1), 60)]
        fitted = [self.manual_fit(manual_patterns, states, weights) for states, weights in specifications]
        best_ll = max(result[0] for result in fitted)
        relative_terms = [math.exp(max(-745, result[0] - best_ll)) for result in fitted]
        normalizer = math.fsum(relative_terms)
        candidate_rows = []
        for index, ((states, weights), (ll, gap, simplex), relative) in enumerate(
            zip(specifications, fitted, relative_terms)
        ):
            candidate_rows.append({
                **common,
                "vertex_set_id": self.vertex_id(2, states),
                "states_json": json.dumps(self.state_records(2, states), separators=(",", ":")),
                "primary_log_likelihood": format(ll, ".12g"),
                "relative_likelihood_weight": format(relative / normalizer, ".12g"),
                "mixture_weights_json": json.dumps(
                    {str(state): weight for state, weight in zip(states, weights)},
                    separators=(",", ":"),
                ),
                "fit_converged": "true", "fit_monotone": "true",
                "optimizer_status": "CERTIFIED_SLSQP_WARM_START",
                "global_log_likelihood_gap_bound": format(gap, ".12g"),
                "simplex_residual": format(simplex, ".12g"),
                "is_winner": "true" if index == 0 else "false",
                "is_tied_winner": "false",
            })
        with gzip.open(candidate_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, verifier.PROFILE_CANDIDATE_REQUIRED_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(candidate_rows)
        child = {
            "parameters": {
                **self.PARAMETERS,
                "optimizer_contract": {
                    "global_log_likelihood_gap_tolerance": 1e-8,
                    "simplex_residual_tolerance": 1e-12,
                },
            },
            "outputs": {
                verifier.PATTERN_SOURCE_NAME: {"path": str(pattern_path)},
                verifier.CANDIDATE_SOURCE_NAME: {"path": str(candidate_path)},
            },
        }
        return child, pattern_path, candidate_path, pattern_rows, candidate_rows

    @staticmethod
    def rewrite(path: pathlib.Path, fields, rows):
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_positive_child_audit_is_directly_callable_and_bounded_by_one_unit(self):
        with tempfile.TemporaryDirectory() as tmp:
            child, _, _, _, _ = self.make_fixture(pathlib.Path(tmp))
            result = verifier.verify_child_profile_likelihood(
                child, "D", "chr1", self.PARAMETERS
            )
            self.assertEqual(result["n_units"], 1)
            self.assertEqual(result["n_candidates"], 2)
            self.assertEqual(result["n_pattern_rows"], 2)
            self.assertEqual(result["n_scoring_molecules"], 100)
            self.assertEqual(result["peak_pattern_rows_per_unit"], 2)
            self.assertEqual(result["peak_candidates_per_unit"], 2)
            self.assertEqual(result["peak_states_per_candidate"], 2)

    def test_true_ll_tie_with_distinct_state_sets_is_reconstructed(self):
        with tempfile.TemporaryDirectory() as tmp:
            child, _, candidate_path, _, rows = self.make_fixture(pathlib.Path(tmp))
            first_weights = json.loads(rows[0]["mixture_weights_json"])
            tied_states = [0, 1, 2]
            tied_weights = [first_weights["0"], first_weights["1"], 0.0]
            tied_ll, tied_gap, tied_simplex = self.manual_fit(
                [("RX", (30, -1), 40), ("AX", (30, -1), 60)],
                tied_states, tied_weights,
            )
            rows[1].update({
                "vertex_set_id": self.vertex_id(2, tied_states),
                "states_json": json.dumps(
                    self.state_records(2, tied_states), separators=(",", ":")
                ),
                "primary_log_likelihood": format(tied_ll, ".12g"),
                "mixture_weights_json": json.dumps(
                    {str(state): weight for state, weight in zip(tied_states, tied_weights)},
                    separators=(",", ":"),
                ),
                "global_log_likelihood_gap_bound": format(tied_gap, ".12g"),
                "simplex_residual": format(tied_simplex, ".12g"),
            })
            for row in rows:
                row["relative_likelihood_weight"] = "0.5"
                row["is_winner"] = "true"
                row["is_tied_winner"] = "true"
            self.rewrite(candidate_path, verifier.PROFILE_CANDIDATE_REQUIRED_COLUMNS, rows)
            result = verifier.verify_child_profile_likelihood(
                child, "D", "chr1", self.PARAMETERS
            )
            self.assertEqual(result["n_candidates"], 2)

    def test_ll_or_vaf_like_additive_score_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            child, _, candidate_path, _, rows = self.make_fixture(pathlib.Path(tmp))
            rows[0]["primary_log_likelihood"] = str(float(rows[0]["primary_log_likelihood"]) + 0.25)
            self.rewrite(candidate_path, verifier.PROFILE_CANDIDATE_REQUIRED_COLUMNS, rows)
            with self.assertRaisesRegex(verifier.VerificationError, "primary_log_likelihood"):
                verifier.verify_child_profile_likelihood(child, "D", "chr1", self.PARAMETERS)

    def test_bq_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            child, pattern_path, _, rows, _ = self.make_fixture(pathlib.Path(tmp))
            rows[0]["fixed_base_qualities"] = "10,"
            self.rewrite(pattern_path, verifier.PROFILE_PATTERN_REQUIRED_COLUMNS, rows)
            with self.assertRaisesRegex(verifier.VerificationError, "primary_log_likelihood"):
                verifier.verify_child_profile_likelihood(child, "D", "chr1", self.PARAMETERS)

    def test_pi_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            child, _, candidate_path, _, rows = self.make_fixture(pathlib.Path(tmp))
            rows[0]["mixture_weights_json"] = '{"0":0.5,"1":0.5}'
            self.rewrite(candidate_path, verifier.PROFILE_CANDIDATE_REQUIRED_COLUMNS, rows)
            with self.assertRaisesRegex(verifier.VerificationError, "primary_log_likelihood"):
                verifier.verify_child_profile_likelihood(child, "D", "chr1", self.PARAMETERS)

    def test_state_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            child, _, candidate_path, _, rows = self.make_fixture(pathlib.Path(tmp))
            states = json.loads(rows[0]["states_json"])
            states[1]["state_rax"] = "R"
            rows[0]["states_json"] = json.dumps(states, separators=(",", ":"))
            self.rewrite(candidate_path, verifier.PROFILE_CANDIDATE_REQUIRED_COLUMNS, rows)
            with self.assertRaisesRegex(verifier.VerificationError, "state_rax"):
                verifier.verify_child_profile_likelihood(child, "D", "chr1", self.PARAMETERS)

    def test_kkt_gap_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            child, _, candidate_path, _, rows = self.make_fixture(pathlib.Path(tmp))
            rows[0]["global_log_likelihood_gap_bound"] = "0.001"
            self.rewrite(candidate_path, verifier.PROFILE_CANDIDATE_REQUIRED_COLUMNS, rows)
            with self.assertRaisesRegex(
                verifier.VerificationError, "global_log_likelihood_gap_bound"
            ):
                verifier.verify_child_profile_likelihood(child, "D", "chr1", self.PARAMETERS)

    def test_winner_or_tie_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            child, _, candidate_path, _, rows = self.make_fixture(pathlib.Path(tmp))
            rows[1]["is_winner"] = "true"
            rows[1]["is_tied_winner"] = "true"
            self.rewrite(candidate_path, verifier.PROFILE_CANDIDATE_REQUIRED_COLUMNS, rows)
            with self.assertRaisesRegex(verifier.VerificationError, "winner flag"):
                verifier.verify_child_profile_likelihood(child, "D", "chr1", self.PARAMETERS)

    def test_aggregate_exposes_full_receipt_contract(self):
        child = {
            "dataset": "D", "chrom": "chr1", "n_units": 1, "n_candidates": 2,
            "n_pattern_rows": 2, "n_scoring_molecules": 100,
            "max_abs_ll_delta": 1e-12, "max_abs_relative_weight_delta": 2e-12,
            "max_abs_gap_delta": 3e-12, "max_abs_simplex_residual_delta": 4e-12,
            "peak_pattern_rows_per_unit": 2, "peak_candidates_per_unit": 2,
            "peak_states_per_candidate": 2,
            "all_profile_likelihoods_and_certificates_match": True,
            "all_relative_weights_match": True, "all_winner_tie_partitions_match": True,
        }
        result = verifier.aggregate_profile_likelihood_audits([child])
        self.assertEqual(result["n_children"], 1)
        self.assertTrue(result["all_profile_likelihoods_and_certificates_match"])
        self.assertIn("one child opened at a time", result["streaming_memory_contract"])

    def test_new_profile_checks_are_part_of_all_pass(self):
        extraction = {"children": {}, "n_tasks": 1}
        ranking = {
            "n_tasks": 1,
            "all_aggregate_cells_conserved": True,
            "candidate_table": {
                "n_units": 1, "n_rows": 1,
                "all_rows_match_independent_child_reconstruction": True,
                "winner_partitions_conserved": True,
            },
            "runtime_diagnostics": {
                "all_child_and_full_runtime_summaries_independently_recomputed": True,
            },
            "method_contract_verification": {
                "all_children_exactly_matched_and_source_bound": True,
            },
            "profile_likelihood_independent_recomputation": {
                "n_units": 1, "n_candidates": 1,
                "all_profile_likelihoods_and_certificates_match": False,
                "all_relative_weights_match": True,
                "all_winner_tie_partitions_match": True,
            },
        }
        with (
            mock.patch.object(verifier, "verify_extraction_bundle", return_value=extraction),
            mock.patch.object(verifier, "verify_ranking_bundle", return_value=ranking),
        ):
            receipt = verifier.verify_bundle(
                pathlib.Path("/unused/extraction"), pathlib.Path("/unused/ranking"),
                ("D",), ("chr1",),
            )
        self.assertFalse(receipt["checks"]["profile_likelihood_certificates_match"])
        self.assertFalse(receipt["all_pass"])
        self.assertTrue(receipt["verification_independence"][
            "profile_likelihood_recomputed_from_patterns_states_pi"
        ])


class CandidateTableReconstructionTests(unittest.TestCase):
    SOURCE_FIELDS = (
        "dataset", "chrom", "component_basis", "phase_set", "threshold", "component_id", "family",
        "structural_exact_pattern_minread", "vertex_set_id", "states_json", "parent_choice_count",
        "primary_log_likelihood", "mixture_weights_json", "fit_converged", "fit_monotone",
        "coarse_topology_classes_json", "is_winner", "is_tied_winner",
    )

    def make_fixture(self, root: pathlib.Path):
        source = root / verifier.CANDIDATE_SOURCE_NAME
        states = [
            {"bitmask": 0, "state_rax": "R", "roles": ["root"]},
            {"bitmask": 1, "state_rax": "A", "roles": ["full-observed"]},
        ]
        source_row = {
            "dataset": "D", "chrom": "chr1", "component_basis": "PS_HP1", "phase_set": "100",
            "threshold": "3", "component_id": "C", "family": "1",
            "structural_exact_pattern_minread": "3", "vertex_set_id": "v1",
            "states_json": json.dumps(states), "parent_choice_count": "1",
            "primary_log_likelihood": "-2.0", "mixture_weights_json": '{"0":0.4,"1":0.6}',
            "fit_converged": "true", "fit_monotone": "true",
            "coarse_topology_classes_json": '["single-only"]', "is_winner": "true",
            "is_tied_winner": "false",
        }
        with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, self.SOURCE_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerow(source_row)
        unit_key = "D|chr01|PS_HP1|PS=100|B003|C|HP1|M3"
        expected = {
            "unit_key": unit_key, "dataset": "D", "chrom": "chr1", "component_id": "C",
            "threshold": "3", "hp_family": "1", "ps": "100", "candidate_id": "C000001",
            "vertex_set_id": "v1", "vertex_states": '{"0":"R","1":"A"}',
            "vertex_roles": '{"0":["root"],"1":["full-observed"]}', "parent_choice_count": "1",
            "profile_log_likelihood": "-2.0", "relative_log_likelihood": "0",
            "mixture_weights_pi": '{"0":0.4,"1":0.6}', "winner_status": "UNIQUE_WINNER",
            "tie_group": "", "coarse_topology_class": '["single-only"]',
            "candidate_set_complete": "true",
        }
        full = root / "m2_ps_aware_candidate_table.tsv.gz"
        with gzip.open(full, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, verifier.CANDIDATE_TABLE_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerow(expected)
        semantic = hashlib.sha256(
            json.dumps(expected, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode() + b"\n"
        ).hexdigest()
        metadata = {
            "path": str(full.resolve()), "size_bytes": full.stat().st_size,
            "sha256": verifier.sha256_path(full), "semantic_sha256": semantic,
            "columns": list(verifier.CANDIDATE_TABLE_COLUMNS), "n_rows": 1, "n_units": 1,
        }
        child = {"outputs": {verifier.CANDIDATE_SOURCE_NAME: {"path": str(source.resolve())}}}
        return full, metadata, child, expected

    def test_full_candidate_table_matches_independent_child_reconstruction(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            _, metadata, child, _ = self.make_fixture(root)
            audit = verifier.verify_candidate_table(
                root, metadata, {("D", "chr1"): child}, ("D",), ("chr1",), 3
            )
            self.assertEqual(audit["n_rows"], 1)
            self.assertEqual(audit["n_units"], 1)

    def test_self_consistent_but_wrong_full_row_is_rejected_against_child_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            full, metadata, child, row = self.make_fixture(root)
            row["ps"] = "999"
            with gzip.open(full, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, verifier.CANDIDATE_TABLE_COLUMNS, delimiter="\t")
                writer.writeheader()
                writer.writerow(row)
            metadata["size_bytes"] = full.stat().st_size
            metadata["sha256"] = verifier.sha256_path(full)
            metadata["semantic_sha256"] = hashlib.sha256(
                json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode() + b"\n"
            ).hexdigest()
            with self.assertRaises(verifier.VerificationError):
                verifier.verify_candidate_table(
                    root, metadata, {("D", "chr1"): child}, ("D",), ("chr1",), 3
                )


if __name__ == "__main__":
    unittest.main()
