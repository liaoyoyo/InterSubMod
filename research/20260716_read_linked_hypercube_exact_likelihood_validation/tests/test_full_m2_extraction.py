#!/usr/bin/env python3
from __future__ import annotations

import json
import argparse
import hashlib
import pathlib
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_full_m2_extraction as extraction_runner  # noqa: E402
from run_full_m2_extraction import (  # noqa: E402
    AUTOSOMES,
    DATASETS,
    LINKAGE_BASES,
    active_conflicts,
    aggregate_results,
    build_extraction_receipt,
    canonical_sort_results,
    create_batch_start_and_grants,
    ensure_release_orchestration_session,
    ensure_preflight_contract,
    embedded_extractor_producers_match_release,
    load_passing_receipt,
    load_release_contract_binding,
    load_release_orchestration_state,
    run_command_with_process_group_timeout,
    run_task,
    run_specs_rolling,
    scan_existing_specs,
    sha256_path,
    validate_release_extraction_parameters,
    write_child_completion,
    verify_sha256_sidecar,
    write_json_and_sha256_exclusive,
    write_sha256_sidecar,
)


EXPECTED_PARAMETERS = {
    "mapq_min": 20,
    "baseq_min": 20,
    "samtools_threads": 1,
    "bridge_thresholds": [1, 2, 3, 5],
    "component_linkage_bases": list(LINKAGE_BASES),
}


def make_release_contract(root: pathlib.Path):
    def identity(path: pathlib.Path, stored_path: str | None = None):
        observed = path.stat()
        return {
            "path": str(path if stored_path is None else stored_path),
            "st_dev": int(observed.st_dev), "st_ino": int(observed.st_ino),
            "st_nlink": int(observed.st_nlink), "size_bytes": int(observed.st_size),
            "mtime_ns": int(observed.st_mtime_ns), "ctime_ns": int(observed.st_ctime_ns),
            "mode_octal": format(observed.st_mode & 0o7777, "04o"),
            "sha256": sha256_path(path),
        }

    snapshot = root / "code_snapshot"
    origin_root = root / "origin_repo"
    entries = []
    role_paths = {}
    for role, repo_relative in extraction_runner.RELEASE_SOURCE_PATHS.items():
        path = snapshot / repo_relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if role == "canonical_manifest_schema":
            payload = (
                ROOT.parents[1] / extraction_runner.CANONICAL_SCHEMA_RELATIVE
            ).read_bytes()
        else:
            payload = f"# {role}\n".encode()
        path.write_bytes(payload)
        path.chmod(0o444)
        relative = path.relative_to(root).as_posix()
        role_paths[role] = path
        copied = identity(path, relative)
        source = dict(copied)
        source["path"] = str(origin_root / repo_relative)
        entries.append({
            "role": role,
            "repo_relative_path": repo_relative,
            "source": source,
            "snapshot": copied,
        })
    canonical_copy = root / "input_contract" / "canonical_manifest.json"
    canonical_copy.parent.mkdir(parents=True)
    canonical_source = pathlib.Path(
        "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
        "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/"
        "input_manifest.snapshot.json"
    )
    canonical_copy.write_bytes(canonical_source.read_bytes())
    canonical_copy.chmod(0o444)
    role_paths["canonical_manifest_copy"] = canonical_copy
    canonical_sha = sha256_path(canonical_copy)
    assert canonical_sha == extraction_runner.CANONICAL_MANIFEST_SHA256
    pre_snapshot = {"schema_name": "fixture_snapshot", "roles": 42}
    pre_snapshot_sha = extraction_runner.semantic_json_sha256(pre_snapshot)
    pre_document = {
        "schema_name": "intersubmod.m2_input_identity_verification",
        "schema_version": "1.0.0", "task_type": "B_COMPREHENSIVE_VALIDATION",
        "mode": "PRE", "authority_mode": "CANONICAL_V5_FROZEN",
        "validation_evidence_eligible": True,
        "input_identity_snapshot": pre_snapshot,
        "input_identity_snapshot_sha256": pre_snapshot_sha,
        "artifact_audit": {"n_artifacts": 42}, "checks": {"fixture": True},
        "all_pass": True,
        "receipt_integrity": {"scheme": "external_sha256_sidecar_v1", "sidecar_name": "pre.json.sha256", "covers": "pre.json"},
    }
    pre_copy_path = root / "input_contract" / "pre.json"
    pre_copy_path.write_text(json.dumps(pre_document), encoding="utf-8")
    pre_copy_path.chmod(0o444)
    pre_sidecar_path = root / "input_contract" / "pre.json.sha256"
    pre_sidecar_path.write_text(f"{sha256_path(pre_copy_path)}  pre.json\n", encoding="ascii")
    pre_sidecar_path.chmod(0o444)
    pre_copy = identity(pre_copy_path, "input_contract/pre.json")
    pre_copy["sidecar"] = identity(pre_sidecar_path, "input_contract/pre.json.sha256")
    pre_origin = dict(identity(pre_copy_path, str(origin_root / "pre.json")))
    pre_origin["sidecar"] = identity(pre_sidecar_path, str(origin_root / "pre.json.sha256"))
    runtime = {
        "python": {"executable": "/python", "version": "3.test", "implementation": "CPython"},
        "packages": {"numpy": "1", "scipy": "1", "pysam": "1"},
        "samtools": {"path": "/samtools", "version_line": "samtools 1", "htslib_version_line": "Using htslib 1"},
        "platform": "test-platform",
    }
    code_entrypoints = {
        role: (pathlib.Path("code_snapshot") / relative).as_posix()
        for role, relative in extraction_runner.RELEASE_SOURCE_PATHS.items()
        if role != "canonical_manifest_schema"
    }
    freezer_entry = next(row for row in entries if row["role"] == "release_contract_freezer")
    verifier_entry = next(row for row in entries if row["role"] == "input_identity_verifier")
    document = {
        "schema_name": "intersubmod.m2_release_run_manifest",
        "schema_version": "1.0.0",
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "authority_mode": "CANONICAL_V5_FROZEN",
        "validation_evidence_eligible": True,
        "created_at_utc": "2026-07-16T00:00:00Z",
        "scope": {
            "technical_datasets": 7, "biological_samples": 6, "chromosomes": 22,
            "chromosome_names": list(AUTOSOMES), "tasks": 154, "datasets": list(DATASETS),
        },
        "canonical_manifest": {
            "expected_sha256": canonical_sha,
            "origin": identity(canonical_copy, str(origin_root / "canonical_manifest.json")),
            "immutable_copy": identity(canonical_copy, "input_contract/canonical_manifest.json"),
        },
        "pre_input_identity_receipt": {
            "origin": pre_origin, "immutable_copy": pre_copy,
            "input_identity_snapshot_sha256": pre_snapshot_sha,
            "receipt_semantic_sha256": extraction_runner.semantic_json_sha256(pre_document),
            "authority_mode": "CANONICAL_V5_FROZEN", "validation_evidence_eligible": True,
            "artifact_roles": 42,
        },
        "fresh_input_identity_verification": {
            "execution_mode": "production_subprocess_required",
            "verifier_path": verifier_entry["source"]["path"],
            "verifier_sha256": verifier_entry["source"]["sha256"],
            "receipt_sha256": "1" * 64, "receipt_semantic_sha256": "2" * 64,
            "checks_semantic_sha256": "3" * 64, "artifact_audit_semantic_sha256": "4" * 64,
            "input_identity_snapshot_sha256": pre_snapshot_sha, "all_pass": True,
            "validation_evidence_eligible": True, "exactly_equals_supplied_pre_snapshot": True,
        },
        "entrypoints": {**code_entrypoints,
            "canonical_manifest_copy": "input_contract/canonical_manifest.json",
            "pre_input_identity_receipt_copy": "input_contract/pre.json",
            "canonical_schema_copy": (
                pathlib.Path("code_snapshot") / extraction_runner.CANONICAL_SCHEMA_RELATIVE
            ).as_posix(),
        },
        "source_snapshot": {
            "repo_root": str(origin_root),
            "snapshot_root": "code_snapshot", "n_files": len(entries), "entries": entries,
            "entrypoints": code_entrypoints,
            "exact_allowlist_semantic_sha256": extraction_runner.semantic_json_sha256([
                {"role": role, "repo_relative_path": relative}
                for role, relative in extraction_runner.RELEASE_SOURCE_PATHS.items()
            ]),
        },
        "canonical_schema": {"role": "canonical_manifest_schema", "repo_relative_path": extraction_runner.CANONICAL_SCHEMA_RELATIVE, "sha256": extraction_runner.CANONICAL_SCHEMA_SHA256},
        "runtime": runtime,
        "runtime_semantic_sha256": extraction_runner.semantic_json_sha256(runtime),
        "parameters": extraction_runner.FROZEN_RELEASE_PARAMETERS,
        "producer": {
            "role": "release_contract_freezer",
            "repo_relative_path": extraction_runner.RELEASE_SOURCE_PATHS["release_contract_freezer"],
            "source_sha256": freezer_entry["source"]["sha256"],
            "immutable_copy_path": freezer_entry["snapshot"]["path"],
            "immutable_copy_sha256": freezer_entry["snapshot"]["sha256"],
        },
        "assurance": extraction_runner.RELEASE_ASSURANCE,
        "checks": {name: True for name in extraction_runner.RELEASE_CHECK_NAMES},
        "all_pass": True,
        "receipt_integrity": {"scheme": "external_sha256_sidecar_v1", "sidecar_name": "m2_run_manifest.json.sha256", "covers": "m2_run_manifest.json"},
    }
    manifest = root / "m2_run_manifest.json"
    manifest.write_text(json.dumps(document), encoding="utf-8")
    manifest.chmod(0o444)
    sidecar = write_sha256_sidecar(manifest)
    sidecar.chmod(0o444)
    return manifest, role_paths


def rewrite_release_manifest(manifest: pathlib.Path, document: dict) -> None:
    manifest.chmod(0o644)
    manifest.write_text(json.dumps(document), encoding="utf-8")
    manifest.chmod(0o444)
    sidecar = manifest.with_name(f"{manifest.name}.sha256")
    sidecar.chmod(0o644)
    sidecar.write_text(
        f"{sha256_path(manifest)}  {manifest.name}\n", encoding="ascii"
    )
    sidecar.chmod(0o444)


class ReleaseBindingTests(unittest.TestCase):
    def test_authenticates_sidecar_and_exact_runner_dependency_paths_and_shas(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, roles = make_release_contract(pathlib.Path(tmp))
            binding = load_release_contract_binding(
                manifest,
                required_sources={
                    "full_extraction_runner": roles["full_extraction_runner"],
                    "extractor": roles["extractor"],
                    "lossless_read_contract": roles["lossless_read_contract"],
                },
                input_manifest=roles["canonical_manifest_copy"],
                _skip_deep_verification_for_test=True,
            )
            self.assertEqual(
                binding["snapshot_sources"]["extractor"]["sha256"],
                sha256_path(roles["extractor"]),
            )
            self.assertFalse(binding["validation_evidence_eligible"])
            child = [{
                "receipt": {"provenance": {
                    "extractor": dict(binding["snapshot_sources"]["extractor"]),
                    "manifest": dict(binding["canonical_input_manifest"]),
                }}
            }]
            self.assertTrue(embedded_extractor_producers_match_release(child, binding))
            child[0]["receipt"]["provenance"]["extractor"]["sha256"] = "0" * 64
            self.assertFalse(embedded_extractor_producers_match_release(child, binding))

            sidecar = manifest.with_name(f"{manifest.name}.sha256")
            sidecar.chmod(0o644)
            sidecar.write_text(
                f"{'0' * 64}  {manifest.name}\n", encoding="ascii"
            )
            sidecar.chmod(0o444)
            with self.assertRaisesRegex(RuntimeError, "sidecar mismatch"):
                load_release_contract_binding(
                    manifest,
                    required_sources={"full_extraction_runner": roles["full_extraction_runner"]},
                    input_manifest=roles["canonical_manifest_copy"],
                    _skip_deep_verification_for_test=True,
                )

    def test_wrong_source_path_or_sha_and_parameter_or_gate_drift_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, roles = make_release_contract(pathlib.Path(tmp))
            wrong = pathlib.Path(tmp) / "wrong.py"
            wrong.write_text("# wrong\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "current source"):
                load_release_contract_binding(
                    manifest,
                    required_sources={"extractor": wrong},
                    input_manifest=roles["canonical_manifest_copy"],
                    _skip_deep_verification_for_test=True,
                )
            roles["lossless_read_contract"].chmod(0o644)
            roles["lossless_read_contract"].write_text("# tampered\n", encoding="utf-8")
            roles["lossless_read_contract"].chmod(0o444)
            with self.assertRaisesRegex(
                RuntimeError,
                "snapshot/lossless_read_contract identity/stat/SHA mismatch",
            ):
                load_release_contract_binding(
                    manifest,
                    required_sources={"extractor": roles["extractor"]},
                    input_manifest=roles["canonical_manifest_copy"],
                    _skip_deep_verification_for_test=True,
                )

        with tempfile.TemporaryDirectory() as tmp:
            manifest, roles = make_release_contract(pathlib.Path(tmp))
            binding = load_release_contract_binding(
                manifest,
                required_sources={"extractor": roles["extractor"]},
                input_manifest=roles["canonical_manifest_copy"],
                _skip_deep_verification_for_test=True,
            )
            args = argparse.Namespace(
                mapq_min=20, baseq_min=20, workers=2, samtools_threads=1,
                task_timeout_seconds=28800.0, timeout_grace_seconds=300.0,
                max_new_tasks=8, ignore_resource_gate=False,
            )
            validate_release_extraction_parameters(binding, args, (1, 2, 3, 5))
            args.baseq_min = 21
            with self.assertRaisesRegex(RuntimeError, "CLI parameters drift"):
                validate_release_extraction_parameters(binding, args, (1, 2, 3, 5))
            args.baseq_min = 20
            args.ignore_resource_gate = True
            with self.assertRaisesRegex(RuntimeError, "forbidden"):
                validate_release_extraction_parameters(binding, args, (1, 2, 3, 5))

    def test_minimal_resigned_self_asserted_contract_is_rejected_by_exact_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, roles = make_release_contract(pathlib.Path(tmp))
            document = json.loads(manifest.read_text(encoding="utf-8"))
            document["parameters"]["extraction"]["mapq_min"] = 0
            document["checks"] = {"trust_me": True}
            rewrite_release_manifest(manifest, document)
            with self.assertRaisesRegex(RuntimeError, "release checks exact-key schema mismatch"):
                load_release_contract_binding(
                    manifest,
                    required_sources={"full_extraction_runner": roles["full_extraction_runner"]},
                    input_manifest=roles["canonical_manifest_copy"],
                    _skip_deep_verification_for_test=True,
                )

    def test_missing_exact_entrypoint_and_pre_escape_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, roles = make_release_contract(pathlib.Path(tmp))
            document = json.loads(manifest.read_text(encoding="utf-8"))
            del document["source_snapshot"]["entrypoints"]["extractor"]
            rewrite_release_manifest(manifest, document)
            with self.assertRaisesRegex(RuntimeError, "exact entrypoint map mismatch"):
                load_release_contract_binding(
                    manifest,
                    required_sources={"full_extraction_runner": roles["full_extraction_runner"]},
                    input_manifest=roles["canonical_manifest_copy"],
                    _skip_deep_verification_for_test=True,
                )

        with tempfile.TemporaryDirectory() as tmp:
            manifest, roles = make_release_contract(pathlib.Path(tmp))
            document = json.loads(manifest.read_text(encoding="utf-8"))
            document["pre_input_identity_receipt"]["immutable_copy"]["path"] = "../pre.json"
            rewrite_release_manifest(manifest, document)
            with self.assertRaisesRegex(RuntimeError, "PRE copy path escapes"):
                load_release_contract_binding(
                    manifest,
                    required_sources={"full_extraction_runner": roles["full_extraction_runner"]},
                    input_manifest=roles["canonical_manifest_copy"],
                    _skip_deep_verification_for_test=True,
                )

    def test_manifest_leaf_symlink_is_rejected_before_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest, roles = make_release_contract(root)
            symlink_manifest = root / "linked_manifest.json"
            symlink_manifest.symlink_to(manifest.name)
            with self.assertRaisesRegex(RuntimeError, "immutable single-link 0444"):
                load_release_contract_binding(
                    symlink_manifest,
                    required_sources={"full_extraction_runner": roles["full_extraction_runner"]},
                    input_manifest=roles["canonical_manifest_copy"],
                    _skip_deep_verification_for_test=True,
                )

    def test_malicious_freezer_never_executes_before_static_authentication(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            manifest, roles = make_release_contract(root)
            sentinel = root / "MALICIOUS_FREEZER_EXECUTED"
            freezer = roles["release_contract_freezer"]
            freezer.chmod(0o644)
            freezer.write_text(
                "from pathlib import Path\n"
                f"Path({str(sentinel)!r}).write_text('executed')\n",
                encoding="utf-8",
            )
            freezer.chmod(0o444)
            with self.assertRaisesRegex(
                RuntimeError, "snapshot/release_contract_freezer identity/stat/SHA mismatch"
            ):
                load_release_contract_binding(
                    manifest,
                    required_sources={"full_extraction_runner": roles["full_extraction_runner"]},
                    input_manifest=roles["canonical_manifest_copy"],
                )
            self.assertFalse(sentinel.exists())

    def test_publication_boundary_force_flag_bypasses_deep_verify_cache_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, roles = make_release_contract(pathlib.Path(tmp))
            summary = {
                "mode": "FROZEN_FREEZER_VERIFY_RELEASE_CONTRACT", "all_pass": True,
                "freezer_path": str(roles["release_contract_freezer"]),
                "freezer_sha256": sha256_path(roles["release_contract_freezer"]),
                "release_manifest_sha256": sha256_path(manifest),
                "verified_snapshot_files": 11,
            }
            with mock.patch.object(
                extraction_runner, "_deep_verify_frozen_release", return_value=summary
            ) as deep:
                common = {
                    "required_sources": {
                        "full_extraction_runner": roles["full_extraction_runner"],
                        "extractor": roles["extractor"],
                        "lossless_read_contract": roles["lossless_read_contract"],
                    },
                    "input_manifest": roles["canonical_manifest_copy"],
                }
                load_release_contract_binding(manifest, **common)
                load_release_contract_binding(
                    manifest, _force_deep_reverification=True, **common
                )
            self.assertFalse(deep.call_args_list[0].kwargs["force_fresh"])
            self.assertTrue(deep.call_args_list[1].kwargs["force_fresh"])


class ReleaseOrchestrationTests(unittest.TestCase):
    @staticmethod
    def fixture(root: pathlib.Path):
        manifest = root / "canonical_manifest.json"
        manifest.write_text("{}\n", encoding="utf-8")
        source = lambda name: {
            "path": str(root / name),
            "sha256": hashlib.sha256(name.encode()).hexdigest(),
        }
        binding = {
            "release_manifest": {
                "path": str(root / "m2_run_manifest.json"), "sha256": "a" * 64,
                "semantic_sha256": "b" * 64,
                "sidecar": {"path": str(root / "m2_run_manifest.json.sha256"), "sha256": "c" * 64},
            },
            "canonical_input_manifest": {"path": str(manifest.resolve()), "sha256": sha256_path(manifest)},
            "snapshot_sources": {
                "full_extraction_runner": source("runner.py"),
                "extractor": source("extractor.py"),
                "lossless_read_contract": source("lossless.py"),
            },
        }
        run_contract = {
            "manifest_sha256": sha256_path(manifest), "release_binding": binding,
            "orchestration_policy": dict(extraction_runner.ORCHESTRATION_POLICY),
            "bridge_thresholds": [1, 2, 3, 5],
        }
        return manifest, binding, run_contract

    @staticmethod
    def create_session(root, binding, run_contract):
        root.mkdir(parents=True)
        observed = root.stat()
        target = {
            "output_root": {
                "path": str(root.resolve()), "st_dev": int(observed.st_dev),
                "st_ino": int(observed.st_ino),
            },
            "release_manifest": dict(binding["release_manifest"]),
            "run_contract_semantic_sha256": extraction_runner.semantic_json_sha256(run_contract),
        }
        _, gate = extraction_runner.create_resource_gate_receipt(
            root, stage="extraction", gate_scope="session", target=target,
            producer_source=extraction_runner.release_producer_sources(binding)["runner"],
            conflicts={"process_count": 0, "root_count": 0, "representatives": []},
        )
        session = ensure_release_orchestration_session(root, binding, run_contract, gate)
        return session

    @staticmethod
    def create_batch_gate(root, session, state, specs):
        batch_index = state["next_batch_index"]
        before_count = state["count"]
        max_new_tasks = 8 if before_count == 0 else 16
        selected_ids = [f"{spec[0]}:{spec[1]}" for spec in specs]
        target = {
            "output_root": dict(session["output_root"]),
            "session_id": session["session_id"],
            "session_sha256": sha256_path(root / "_orchestration/session.json"),
            "batch_index": batch_index,
            "before_count": before_count,
            "max_new_tasks": max_new_tasks,
            "effective_count": len(specs),
            "selected_task_ids": selected_ids,
            "previous_chain_head": state["previous_chain_head"],
        }
        _, gate = extraction_runner.create_resource_gate_receipt(
            root, stage="extraction", gate_scope="batch", batch_index=batch_index,
            target=target, producer_source=session["producer_sources"]["runner"],
            conflicts={"process_count": 0, "root_count": 0, "representatives": []},
        )
        return gate

    @staticmethod
    def specs(root: pathlib.Path, manifest: pathlib.Path, start: int, count: int):
        parameters = dict(EXPECTED_PARAMETERS)
        values = []
        pairs = [(dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES]
        for dataset, chrom in pairs[start:start + count]:
            task_dir = root / "samples" / dataset / chrom
            command = [
                "/python", "/extractor.py", "--manifest", str(manifest.resolve()),
                "--sample", dataset, "--chrom", chrom, "--output-dir", str(task_dir.resolve()),
                "--mapq-min", "20", "--baseq-min", "20", "--bridge-thresholds", "1,2,3,5",
                "--samtools-threads", "1",
            ]
            values.append((
                dataset, chrom, task_dir, command, parameters,
                sha256_path(manifest), "e" * 64, 28800.0, 300.0,
            ))
        return values

    @staticmethod
    def write_child(task_dir: pathlib.Path, dataset: str, chrom: str, manifest: pathlib.Path):
        task_dir.mkdir(parents=True)
        output = task_dir / f"{dataset}.{chrom}.calls.tsv.gz"
        output.write_bytes(b"calls\n")
        receipt_path = task_dir / "receipt.json"
        receipt = {
            "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
            "schema_version": "1.2.0", "scope": {"dataset": dataset, "chrom": chrom},
            "provenance": {
                "manifest": {"path": str(manifest.resolve()), "sha256": sha256_path(manifest)},
                "extractor": {"path": "/extractor.py", "sha256": "e" * 64},
            },
            "parameters": {
                **dict(EXPECTED_PARAMETERS),
                "excluded_flags": 3844,
                "bridge_definition": "production-like extra receipt metadata",
            },
            "outputs": {output.name: {
                "path": str(output.resolve()), "size_bytes": output.stat().st_size,
                "sha256": sha256_path(output),
            }},
            "counts": {"synthetic_units": 1},
            "phase_set_contract_counts": {},
            "legacy_cross_phase_set_aggregation_audit": {},
            "component_summary_by_linkage_basis": {
                basis: {
                    str(threshold): {
                        "n_components": 0,
                        "max_k_component_sites": 0,
                        "max_k": 0,
                        "k_component_sites_distribution": {},
                        "k_distribution": {},
                    }
                    for threshold in (1, 2, 3, 5)
                }
                for basis in LINKAGE_BASES
            },
            "all_pass": True,
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": "receipt.json.sha256", "covers": "receipt.json",
            },
        }
        receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
        write_sha256_sidecar(receipt_path)
        return receipt

    def publish_batch(
        self, root, manifest, session, run_contract, state, count,
        *, mutate_child=None,
    ):
        specs = self.specs(root, manifest, state["count"], count)
        gate = self.create_batch_gate(root, session, state, specs)
        batch, grants = create_batch_start_and_grants(
            root, session, run_contract, specs,
            before_count=state["count"], previous_chain_head=state["previous_chain_head"],
            batch_index=state["next_batch_index"],
            max_new_tasks=8 if state["count"] == 0 else 16, gate=gate,
        )
        raw_results = []
        for spec in reversed(specs):
            dataset, chrom, task_dir = spec[:3]
            receipt = self.write_child(task_dir, dataset, chrom, manifest)
            if mutate_child is not None:
                mutate_child(task_dir, dataset, chrom, receipt)
            raw_results.append({
                "dataset": dataset, "chrom": chrom, "status": "PASS", "returncode": 0,
                "process_group_id": 1000 + len(raw_results),
                "started_monotonic_ns": time.monotonic_ns(),
                "receipt": receipt,
            })
        results = canonical_sort_results(raw_results)
        for result in results:
            task_id = f"{result['dataset']}:{result['chrom']}"
            result["orchestration_completion"] = write_child_completion(
                root, session, batch, grants[task_id], result
            )
        cumulative_results = [
            {
                "dataset": task_id.split(":")[0], "chrom": task_id.split(":")[1],
                "status": "REUSED_PASS",
                "receipt": state["child_receipts"][task_id],
                "orchestration_completion": identity,
            }
            for task_id, identity in state["completions"].items()
        ] + results
        cumulative_results = canonical_sort_results(cumulative_results)
        cumulative = len(cumulative_results)
        path = (
            root / "full_extraction_receipt.json"
            if cumulative == extraction_runner.EXPECTED_TASKS
            else root / "checkpoints" / f"checkpoint_{cumulative:03d}_of_154.json"
        )
        receipt = build_extraction_receipt(
            cumulative_results,
            (1, 2, 3, 5),
            run_contract,
            {
                "max_new_tasks": 8 if state["count"] == 0 else 16,
                "reused_tasks": state["count"],
                "selected_task_ids": [
                    f"{row['dataset']}:{row['chrom']}" for row in results
                ],
            },
            1.0,
        )
        receipt["orchestration"] = {
            "session_identity": extraction_runner._session_identity(root, session),
            "batch_start_identity": {
                "path": batch["path"], "sha256": batch["sha256"],
                "batch_id": batch["batch_id"], "batch_index": batch["batch_index"],
            },
            "previous_chain_head": state["previous_chain_head"],
            "batch_completion_attestations": [
                {"task_id": f"{row['dataset']}:{row['chrom']}", **row["orchestration_completion"]}
                for row in results
            ],
            "cumulative_attested_tasks": cumulative,
        }
        receipt["receipt_integrity"] = {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{path.name}.sha256", "covers": path.name,
        }
        extraction_runner.write_immutable_json_exclusive(path, receipt)
        return load_release_orchestration_state(root, session), raw_results, results

    def test_existing_root_without_session_and_symlink_parent_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            _, binding, run_contract = self.fixture(base)
            existing = base / "existing"
            existing.mkdir()
            with self.assertRaisesRegex(RuntimeError, "requires.*resource gate"):
                ensure_release_orchestration_session(existing, binding, run_contract, None)

            real = base / "real"
            real.mkdir()
            linked = base / "linked"
            linked.symlink_to(real, target_is_directory=True)
            with self.assertRaisesRegex(RuntimeError, "symlink component"):
                ensure_release_orchestration_session(
                    linked / "release", binding, run_contract, None
                )

    def test_first_batch_is_exactly_8_open_batch_rejected_and_last_batch_is_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "release"
            session = self.create_session(root, binding, run_contract)
            state = load_release_orchestration_state(root, session)
            specs = self.specs(root, manifest, 0, 8)
            gate = self.create_batch_gate(root, session, state, specs)
            with self.assertRaisesRegex(RuntimeError, "illegal release batch"):
                create_batch_start_and_grants(
                    root, session, run_contract, specs,
                    before_count=0, previous_chain_head=None, batch_index=1,
                    max_new_tasks=16, gate=gate,
                )
            create_batch_start_and_grants(
                root, session, run_contract, specs,
                before_count=0, previous_chain_head=None, batch_index=1,
                max_new_tasks=8, gate=gate,
            )
            with self.assertRaisesRegex(RuntimeError, "open/orphan"):
                load_release_orchestration_state(root, session)

        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "last"
            session = self.create_session(root, binding, run_contract)
            specs = self.specs(root, manifest, 152, 2)
            state = {
                "count": 152, "next_batch_index": 11,
                "previous_chain_head": {"path": "/previous", "sha256": "d" * 64},
            }
            gate = self.create_batch_gate(root, session, state, specs)
            batch, _ = create_batch_start_and_grants(
                root, session, run_contract, specs, before_count=152,
                previous_chain_head={"path": "/previous", "sha256": "d" * 64},
                batch_index=11, max_new_tasks=16, gate=gate,
            )
            self.assertEqual(batch["effective_count"], 2)

    def test_legal_8_to_24_revalidates_first_batch_and_canonicalizes_reverse_completion(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "release"
            session = self.create_session(root, binding, run_contract)
            state = load_release_orchestration_state(root, session)
            state, raw, ordered = self.publish_batch(
                root, manifest, session, run_contract, state, 8
            )
            self.assertNotEqual(
                [f"{r['dataset']}:{r['chrom']}" for r in raw],
                [f"{r['dataset']}:{r['chrom']}" for r in ordered],
            )
            self.assertEqual(state["count"], 8)
            state, _, _ = self.publish_batch(
                root, manifest, session, run_contract, state, 16
            )
            self.assertEqual(state["count"], 24)
            self.assertEqual(state["next_batch_index"], 3)

    def test_wrong_session_completion_and_timeout_marker_with_late_receipt_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "release"
            session = self.create_session(root, binding, run_contract)
            state, _, _ = self.publish_batch(
                root, manifest, session, run_contract,
                load_release_orchestration_state(root, session), 8,
            )
            completion = root / "_orchestration/tasks/COLO829/chr1/completion.json"
            payload = json.loads(completion.read_text(encoding="utf-8"))
            payload["session_id"] = "0" * 64
            rewrite_release_manifest(completion, payload)
            with self.assertRaisesRegex(RuntimeError, "session/task/status"):
                load_release_orchestration_state(root, session)

        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "release"
            session = self.create_session(root, binding, run_contract)
            self.publish_batch(
                root, manifest, session, run_contract,
                load_release_orchestration_state(root, session), 8,
            )
            marker = root / "samples/COLO829/chr1" / extraction_runner.TASK_STATUS_NAME
            marker.write_text("late timeout\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "receipt/output identity"):
                load_release_orchestration_state(root, session)

    def test_attested_parameter_subset_tamper_and_wrong_receipt_identity_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "parameter_tamper"
            session = self.create_session(root, binding, run_contract)

            def tamper_parameter(task_dir, dataset, chrom, receipt):
                if (dataset, chrom) != ("COLO829", "chr1"):
                    return
                receipt["parameters"]["mapq_min"] = 21
                rewrite_release_manifest(task_dir / "receipt.json", receipt)

            with self.assertRaisesRegex(RuntimeError, "parameter identity"):
                self.publish_batch(
                    root, manifest, session, run_contract,
                    load_release_orchestration_state(root, session), 8,
                    mutate_child=tamper_parameter,
                )

        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "receipt_identity_tamper"
            session = self.create_session(root, binding, run_contract)
            self.publish_batch(
                root, manifest, session, run_contract,
                load_release_orchestration_state(root, session), 8,
            )
            completion = root / "_orchestration/tasks/COLO829/chr1/completion.json"
            payload = json.loads(completion.read_text(encoding="utf-8"))
            payload["child_receipt_identity"]["sha256"] = "0" * 64
            rewrite_release_manifest(completion, payload)
            with self.assertRaisesRegex(RuntimeError, "receipt/output identity"):
                load_release_orchestration_state(root, session)

    def test_attested_declared_output_with_second_hardlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "hardlink_tamper"
            session = self.create_session(root, binding, run_contract)
            self.publish_batch(
                root, manifest, session, run_contract,
                load_release_orchestration_state(root, session), 8,
            )
            output = root / "samples/COLO829/chr1/COLO829.chr1.calls.tsv.gz"
            output.link_to(output.with_name("second_link.tsv.gz"))
            with self.assertRaisesRegex(RuntimeError, "declared output identity"):
                load_release_orchestration_state(root, session)

    def test_renamed_and_resigned_checkpoint_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "renamed_checkpoint"
            session = self.create_session(root, binding, run_contract)
            self.publish_batch(
                root, manifest, session, run_contract,
                load_release_orchestration_state(root, session), 8,
            )
            original = root / "checkpoints/checkpoint_008_of_154.json"
            renamed = root / "checkpoints/checkpoint_999_of_154.json"
            original.rename(renamed)
            original.with_name(f"{original.name}.sha256").rename(
                renamed.with_name(f"{renamed.name}.sha256")
            )
            payload = json.loads(renamed.read_text(encoding="utf-8"))
            payload["receipt_integrity"] = extraction_runner._receipt_integrity(renamed)
            rewrite_release_manifest(renamed, payload)
            with self.assertRaisesRegex(RuntimeError, "exact path/integrity"):
                load_release_orchestration_state(root, session)

    def test_resigned_checkpoint_aggregate_and_checks_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            manifest, binding, run_contract = self.fixture(base)
            root = base / "checkpoint_aggregate_tamper"
            session = self.create_session(root, binding, run_contract)
            self.publish_batch(
                root, manifest, session, run_contract,
                load_release_orchestration_state(root, session), 8,
            )
            checkpoint = root / "checkpoints/checkpoint_008_of_154.json"
            payload = json.loads(checkpoint.read_text(encoding="utf-8"))
            payload["aggregate"]["counts"]["synthetic_units"] = 10**18
            payload["checks"]["passing_count_matches_results"] = False
            rewrite_release_manifest(checkpoint, payload)
            with self.assertRaisesRegex(RuntimeError, "aggregate, checks"):
                load_release_orchestration_state(root, session)


class ResourceGateTests(unittest.TestCase):
    def test_declared_comprehensive_scope_is_154_dataset_chromosome_tasks(self):
        self.assertEqual(len(DATASETS), 7)
        self.assertEqual(len(AUTOSOMES), 22)
        self.assertEqual(len(DATASETS) * len(AUTOSOMES), 154)

    def test_ancestors_and_descendants_are_excluded_but_pinned_jobs_are_grouped(self):
        process_table = "\n".join(
            [
                "1 0 9999 0.0 /sbin/init",
                "100 1 40 0.0 /usr/bin/time -v python run_full_m2_extraction.py --workers 1",
                "101 100 39 1.0 python run_full_m2_extraction.py --workers 1",
                "102 101 1 0.0 samtools view -u input.bam",
                "200 1 500 3.0 python run_source_locked_screen.py --producer /x/analyze_all_ssnv_focal_alt_multigroup.pinned_abc.py",
                "201 200 490 80.0 python /x/analyze_all_ssnv_focal_alt_multigroup.pinned_abc.py --workers 44",
                "202 1 300 2.0 python /x/audit_cooccurrence_task_contract_preflight.py --workers 40",
                "203 202 290 95.0 python /x/audit_cooccurrence_task_contract_preflight.py --worker",
                "300 1 20 4.0 python /x/extract_lossless_read_linkage.py --chrom chr22",
                "500 1 10 2.0 python /x/run_full_m2_ranking.py --max-new-tasks 8",
                "501 500 9 90.0 python /x/build_m2_patterns_and_rank.py --chrom chr1",
                "400 1 1 0.0 rg analyze_all_ssnv_focal_alt_multigroup|run_full_m2_extraction",
            ]
        )
        result = active_conflicts(process_table=process_table, current_pid=101)
        self.assertEqual(result["process_count"], 7)
        self.assertEqual(result["root_count"], 4)
        self.assertEqual(
            [row["pid"] for row in result["representatives"]],
            [200, 202, 300, 500],
        )
        self.assertEqual(result["representatives"][0]["group_process_count"], 2)
        self.assertEqual(result["representatives"][0]["member_pids"], [200, 201])
        self.assertEqual(
            result["representatives"][1]["group_kinds"],
            ["all_ssnv_cooccurrence_audit"],
        )

    def test_pilot_gate_is_exclusive_immutable_and_requires_exact_300_gib(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp) / "pilot_output"
            root.mkdir()
            manifest = root.parent / "manifest.json"
            manifest.write_text("{}\n", encoding="utf-8")
            observed = root.stat()
            target = {
                "task_type": "B_COMPREHENSIVE_VALIDATION_RELEASE_PILOT",
                "dataset": "COLO829", "chrom": "chr1", "gate_label": "extraction",
                "output_root": {
                    "path": str(root.resolve()), "st_dev": int(observed.st_dev),
                    "st_ino": int(observed.st_ino),
                },
                "manifest": {
                    "path": str(manifest.resolve()), "sha256": sha256_path(manifest),
                },
                "release_manifest": {"path": "/release.json", "sha256": "a" * 64},
            }
            receipt_path = root.parent / "resource_gates/extraction.json"
            payload, identity = extraction_runner.create_resource_gate_receipt(
                root, stage="extraction", gate_scope="pilot", target=target,
                producer_source={
                    "path": str(extraction_runner.RUNNER.resolve()),
                    "sha256": sha256_path(extraction_runner.RUNNER),
                },
                conflicts={"process_count": 0, "root_count": 0, "representatives": []},
                receipt_path=receipt_path,
            )
            self.assertEqual(
                payload["filesystem_snapshot"]["required_reserve_bytes"],
                300 * 1024**3,
            )
            self.assertEqual(identity["path"], str(receipt_path.resolve()))
            self.assertEqual(receipt_path.stat().st_mode & 0o777, 0o444)
            self.assertEqual(receipt_path.stat().st_nlink, 1)
            with self.assertRaises(FileExistsError):
                extraction_runner.create_resource_gate_receipt(
                    root, stage="extraction", gate_scope="pilot", target=target,
                    producer_source=payload["producer_source"],
                    conflicts={"process_count": 0, "root_count": 0, "representatives": []},
                    receipt_path=receipt_path,
                )

        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            with mock.patch.object(
                extraction_runner.os,
                "statvfs",
                return_value=mock.Mock(f_bavail=1, f_frsize=4096),
            ):
                with self.assertRaisesRegex(RuntimeError, "required_reserve_bytes"):
                    extraction_runner.create_resource_gate_receipt(
                        root, stage="ranking", gate_scope="pilot", target={},
                        producer_source={"path": "/runner", "sha256": "a" * 64},
                        conflicts={
                            "process_count": 0, "root_count": 0, "representatives": [],
                        },
                        receipt_path=root / "low_disk.json",
                    )


class ReceiptReuseTests(unittest.TestCase):
    def make_receipt(self, directory: pathlib.Path) -> tuple[pathlib.Path, pathlib.Path]:
        output = directory / "calls.tsv.gz"
        output.write_bytes(b"verified output\n")
        receipt_path = directory / "receipt.json"
        receipt = {
            "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
            "schema_version": "1.2.0",
            "scope": {"dataset": "COLO829", "chrom": "chr1", "n_sSNV": 2},
            "provenance": {
                "manifest": {"path": "/manifest.json", "sha256": "m" * 64},
                "extractor": {"path": "/extractor.py", "sha256": "e" * 64},
            },
            "parameters": dict(EXPECTED_PARAMETERS),
            "outputs": {
                output.name: {
                    "path": str(output.resolve()),
                    "size_bytes": output.stat().st_size,
                    "sha256": sha256_path(output),
                }
            },
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": "receipt.json.sha256",
                "covers": "receipt.json",
            },
            "checks": {"synthetic": True},
            "all_pass": True,
        }
        receipt_path.write_text(json.dumps(receipt, sort_keys=True) + "\n", encoding="utf-8")
        write_sha256_sidecar(receipt_path)
        return receipt_path, output

    def load(self, receipt_path: pathlib.Path, **overrides):
        arguments = {
            "dataset": "COLO829",
            "chrom": "chr1",
            "expected_parameters": EXPECTED_PARAMETERS,
            "manifest_sha256": "m" * 64,
            "extractor_sha256": "e" * 64,
        }
        arguments.update(overrides)
        return load_passing_receipt(receipt_path, **arguments)

    def test_verified_receipt_is_reusable(self):
        with tempfile.TemporaryDirectory() as tmp:
            receipt_path, _ = self.make_receipt(pathlib.Path(tmp))
            self.assertIsNotNone(self.load(receipt_path))

    def test_output_tamper_invalidates_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            receipt_path, output = self.make_receipt(pathlib.Path(tmp))
            output.write_bytes(b"tampered\n")
            self.assertIsNone(self.load(receipt_path))

    def test_scope_parameter_and_producer_mismatch_invalidate_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            receipt_path, _ = self.make_receipt(pathlib.Path(tmp))
            self.assertIsNone(self.load(receipt_path, chrom="chr2"))
            changed = dict(EXPECTED_PARAMETERS) | {"baseq_min": 30}
            self.assertIsNone(self.load(receipt_path, expected_parameters=changed))
            self.assertIsNone(self.load(receipt_path, extractor_sha256="x" * 64))

    def test_receipt_tamper_invalidates_checksum(self):
        with tempfile.TemporaryDirectory() as tmp:
            receipt_path, _ = self.make_receipt(pathlib.Path(tmp))
            with receipt_path.open("a", encoding="utf-8") as handle:
                handle.write(" \n")
            self.assertIsNone(self.load(receipt_path))


class AggregateSemanticsTests(unittest.TestCase):
    def test_component_counts_sum_but_max_k_uses_max(self):
        thresholds = (1,)

        def result(dataset: str, maximum: int):
            summary = {
                basis: {
                    "1": {
                        "n_components": 2,
                        "n_singletons_k1": 1,
                        "n_multisite_k_gt1": 1,
                        "max_k_component_sites": maximum,
                        "max_k": maximum,
                        "k_component_sites_distribution": {str(maximum): 1, "1": 1},
                        "k_distribution": {str(maximum): 1, "1": 1},
                    }
                }
                for basis in LINKAGE_BASES
            }
            return {
                "dataset": dataset,
                "chrom": "chr1",
                "status": "PASS",
                "receipt": {
                    "counts": {},
                    "phase_set_contract_counts": {},
                    "legacy_cross_phase_set_aggregation_audit": {},
                    "component_summary_by_linkage_basis": summary,
                },
            }

        aggregate = aggregate_results([result("D1", 5), result("D2", 7)], thresholds)
        for basis in LINKAGE_BASES:
            cell = aggregate["component_summary_by_linkage_basis"][basis]["1"]
            self.assertEqual(cell["n_components"], 4)
            self.assertEqual(cell["max_k_component_sites"], 7)
            self.assertEqual(cell["max_k"], 7)
            self.assertEqual(cell["k_component_sites_distribution"], {"1": 2, "5": 1, "7": 1})


class DeterministicBatchAndTimeoutTests(unittest.TestCase):
    @staticmethod
    def write_pass(directory: pathlib.Path, dataset: str, chrom: str) -> None:
        directory.mkdir(parents=True)
        output = directory / "calls.tsv.gz"
        output.write_bytes(f"{dataset}\t{chrom}\n".encode())
        receipt_path = directory / "receipt.json"
        receipt = {
            "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
            "schema_version": "1.2.0",
            "scope": {"dataset": dataset, "chrom": chrom},
            "provenance": {
                "manifest": {"sha256": "m" * 64},
                "extractor": {"sha256": "e" * 64},
            },
            "parameters": dict(EXPECTED_PARAMETERS),
            "outputs": {
                output.name: {
                    "path": str(output.resolve()),
                    "size_bytes": output.stat().st_size,
                    "sha256": sha256_path(output),
                }
            },
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": "receipt.json.sha256",
            },
            "all_pass": True,
        }
        receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
        write_sha256_sidecar(receipt_path)

    @staticmethod
    def specs(root: pathlib.Path, n: int = 20):
        identities = [
            (dataset, chrom) for dataset in DATASETS for chrom in AUTOSOMES
        ][:n]
        return [
            (
                dataset,
                chrom,
                root / "samples" / dataset / chrom,
                ["unused"],
                EXPECTED_PARAMETERS,
                "m" * 64,
                "e" * 64,
                28800.0,
                300.0,
            )
            for dataset, chrom in identities
        ]

    def test_resume_8_to_16_validates_reused_and_selects_canonical_next_tasks(self):
        with tempfile.TemporaryDirectory() as tmp:
            specs = self.specs(pathlib.Path(tmp))
            reused, pending, invalid = scan_existing_specs(specs)
            self.assertEqual((len(reused), len(pending), invalid), (0, 20, []))
            self.assertEqual([(row[0], row[1]) for row in pending[:8]], [
                ("COLO829", f"chr{index}") for index in range(1, 9)
            ])
            for dataset, chrom, task_dir, *_ in pending[:8]:
                self.write_pass(task_dir, dataset, chrom)
            reused, pending, invalid = scan_existing_specs(specs)
            self.assertEqual((len(reused), len(pending), invalid), (8, 12, []))
            self.assertEqual([(row[0], row[1]) for row in pending[:8]], [
                ("COLO829", f"chr{index}") for index in range(9, 17)
            ])
            for dataset, chrom, task_dir, *_ in pending[:8]:
                self.write_pass(task_dir, dataset, chrom)
            reused, pending, invalid = scan_existing_specs(specs)
            self.assertEqual((len(reused), len(pending), invalid), (16, 4, []))

    def test_existing_nonpass_directory_fails_scan_before_any_new_task(self):
        with tempfile.TemporaryDirectory() as tmp:
            specs = self.specs(pathlib.Path(tmp), n=2)
            specs[0][2].mkdir(parents=True)
            reused, pending, invalid = scan_existing_specs(specs)
            self.assertEqual(reused, [])
            self.assertEqual([(row[0], row[1]) for row in pending], [("COLO829", "chr2")])
            self.assertEqual(invalid[0]["status"], "FAIL_EXISTING_NONPASS_DIRECTORY")

    def test_checkpoint_and_sidecar_are_exclusive(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "checkpoints" / "checkpoint_008_of_154.json"
            payload = {"all_pass": False, "checkpoint_complete": True, "remaining_tasks": 146}
            write_json_and_sha256_exclusive(path, payload)
            original = path.read_bytes()
            self.assertTrue(verify_sha256_sidecar(path))
            with self.assertRaises(FileExistsError):
                write_json_and_sha256_exclusive(path, payload | {"remaining_tasks": 145})
            self.assertEqual(path.read_bytes(), original)

    def test_runner_sha_drift_rejected_but_operational_batch_size_can_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "preflight.json"
            contract = {
                "runner": {"path": "/runner.py", "sha256": "a" * 64},
                "task_timeout_seconds": 28800.0,
            }
            ensure_preflight_contract(
                path,
                {"run_contract": contract, "invocation": {"max_new_tasks": 8}},
                contract,
            )
            ensure_preflight_contract(
                path,
                {"run_contract": contract, "invocation": {"max_new_tasks": 16}},
                contract,
            )
            drift = contract | {"runner": {"path": "/runner.py", "sha256": "b" * 64}}
            with self.assertRaises(RuntimeError):
                ensure_preflight_contract(path, {"run_contract": drift}, drift)

    def test_timeout_terminates_parent_and_grandchild_process_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            pid_path = pathlib.Path(tmp) / "grandchild.pid"
            child_code = "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(60)"
            parent_code = (
                "import pathlib,subprocess,sys,time; "
                f"p=subprocess.Popen([sys.executable,'-c',{child_code!r}], "
                "stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL); "
                f"pathlib.Path({str(pid_path)!r}).write_text(str(p.pid)); "
                "time.sleep(60)"
            )
            result = run_command_with_process_group_timeout(
                [sys.executable, "-c", parent_code],
                task_timeout_seconds=0.5,
                timeout_grace_seconds=0.2,
            )
            self.assertTrue(result["timed_out"])
            self.assertEqual(result["termination_stage"], "KILL")
            grandchild_pid = int(pid_path.read_text())
            deadline = time.monotonic() + 2.0
            while time.monotonic() < deadline:
                stat = pathlib.Path(f"/proc/{grandchild_pid}/stat")
                if not stat.exists() or stat.read_text().split()[2] == "Z":
                    break
                time.sleep(0.02)
            else:
                self.fail(f"grandchild process {grandchild_pid} survived process-group timeout")

    def test_timeout_retains_fail_closed_task_directory_and_status_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            task_dir = pathlib.Path(tmp) / "samples" / "COLO829" / "chr1"
            code = "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(60)"
            spec = (
                "COLO829",
                "chr1",
                task_dir,
                [sys.executable, "-c", code],
                EXPECTED_PARAMETERS,
                "m" * 64,
                "e" * 64,
                0.2,
                0.1,
            )
            result = run_task(spec)
            self.assertEqual(result["status"], "TIMEOUT")
            marker = task_dir / "runner_task_status.json"
            self.assertTrue(marker.is_file())
            marker_payload = json.loads(marker.read_text())
            self.assertEqual(marker_payload["status"], "TIMEOUT")
            self.assertEqual(marker_payload["returncode"], result["returncode"])
            self.assertGreaterEqual(marker_payload["elapsed_seconds"], 0.2)
            self.assertIn("stdout_tail", marker_payload)
            self.assertIn("stderr_tail", marker_payload)
            reused, pending, invalid = scan_existing_specs([spec])
            self.assertEqual((reused, pending), ([], []))
            self.assertEqual(invalid[0]["status"], "FAIL_EXISTING_NONPASS_DIRECTORY")

    def test_non_timeout_failure_retains_marker_and_is_not_implicitly_retried(self):
        with tempfile.TemporaryDirectory() as tmp:
            task_dir = pathlib.Path(tmp) / "samples" / "COLO829" / "chr1"
            spec = (
                "COLO829", "chr1", task_dir,
                [sys.executable, "-c", "import sys; print('bad'); sys.exit(7)"],
                EXPECTED_PARAMETERS, "m" * 64, "e" * 64, 10.0, 1.0,
            )
            result = run_task(spec)
            self.assertEqual((result["status"], result["returncode"]), ("FAIL", 7))
            marker = json.loads((task_dir / "runner_task_status.json").read_text())
            self.assertEqual(marker["status"], "CHILD_FAILED_OR_INVALID_RECEIPT")
            self.assertIn("bad", marker["stdout_tail"])
            reused, pending, invalid = scan_existing_specs([spec])
            self.assertEqual((reused, pending), ([], []))
            self.assertEqual(invalid[0]["status"], "FAIL_EXISTING_NONPASS_DIRECTORY")

    def test_same_done_batch_failure_does_not_submit_replacement(self):
        started = []
        barrier = threading.Barrier(2)
        specs = [(index,) for index in range(3)]

        def task_runner(spec):
            index = spec[0]
            started.append(index)
            if index < 2:
                barrier.wait(timeout=2)
            return {"status": "FAIL" if index == 1 else "PASS", "dataset": "D", "chrom": str(index)}

        original_wait = extraction_runner.concurrent.futures.wait

        def wait_all(futures, **_kwargs):
            return original_wait(futures, return_when=extraction_runner.concurrent.futures.ALL_COMPLETED)

        with mock.patch.object(extraction_runner.concurrent.futures, "wait", side_effect=wait_all):
            results, max_inflight = run_specs_rolling(specs, workers=2, task_runner=task_runner)
        self.assertEqual(set(started), {0, 1})
        self.assertNotIn(2, started)
        self.assertEqual(max_inflight, 2)
        self.assertTrue(any(row["status"] == "FAIL" for row in results))

    def test_terminal_154_payload_preserves_full_aggregate_semantics(self):
        summary = {
            basis: {
                "1": {
                    "n_components": 1,
                    "max_k_component_sites": 1,
                    "max_k": 1,
                    "k_component_sites_distribution": {"1": 1},
                    "k_distribution": {"1": 1},
                }
            }
            for basis in LINKAGE_BASES
        }
        results = [
            {
                "dataset": dataset,
                "chrom": chrom,
                "status": "REUSED_PASS",
                "receipt": {
                    "schema_version": "1.2.0",
                    "counts": {"synthetic_units": 1},
                    "phase_set_contract_counts": {},
                    "legacy_cross_phase_set_aggregation_audit": {},
                    "component_summary_by_linkage_basis": summary,
                },
            }
            for dataset in DATASETS for chrom in AUTOSOMES
        ]
        checkpoint = build_extraction_receipt(
            results[:16],
            (1,),
            {"runner": {"sha256": "a" * 64}},
            {"max_new_tasks": 8, "reused_tasks": 8},
            0.5,
        )
        self.assertEqual(
            checkpoint["schema_name"], "intersubmod.m2_full_extraction_checkpoint"
        )
        self.assertFalse(checkpoint["all_pass"])
        self.assertTrue(checkpoint["checkpoint_complete"])
        self.assertEqual((checkpoint["passing_tasks"], checkpoint["remaining_tasks"]), (16, 138))
        receipt = build_extraction_receipt(
            results,
            (1,),
            {"runner": {"sha256": "a" * 64}},
            {"max_new_tasks": 8},
            1.0,
        )
        self.assertEqual(receipt["schema_name"], "intersubmod.m2_full_extraction_receipt")
        self.assertTrue(receipt["all_pass"])
        self.assertEqual(receipt["remaining_tasks"], 0)
        self.assertEqual(receipt["aggregate"]["counts"]["synthetic_units"], 154)
        self.assertEqual(
            receipt["aggregate"]["component_summary_by_linkage_basis"]["pooled"]["1"]["n_components"],
            154,
        )

        duplicate_scope = list(results)
        duplicate_scope[-1] = dict(duplicate_scope[0])
        rejected = build_extraction_receipt(
            duplicate_scope,
            (1,),
            {"runner": {"sha256": "a" * 64}},
            {"max_new_tasks": 154},
            1.0,
        )
        self.assertEqual(rejected["schema_name"], "intersubmod.m2_full_extraction_checkpoint")
        self.assertFalse(rejected["all_pass"])
        self.assertFalse(rejected["checks"]["recorded_task_pairs_are_unique_and_canonical"])


class TerminalResumeDeepValidationTests(unittest.TestCase):
    @staticmethod
    def make_terminal(root: pathlib.Path):
        summary = {
            basis: {
                "1": {
                    "n_components": 1,
                    "max_k_component_sites": 1,
                    "max_k": 1,
                    "k_component_sites_distribution": {"1": 1},
                    "k_distribution": {"1": 1},
                }
            }
            for basis in LINKAGE_BASES
        }
        reusable = [
            {
                "dataset": dataset,
                "chrom": chrom,
                "status": "REUSED_PASS",
                "receipt": {
                    "schema_version": "1.2.0",
                    "counts": {"synthetic_units": 1},
                    "phase_set_contract_counts": {},
                    "legacy_cross_phase_set_aggregation_audit": {},
                    "component_summary_by_linkage_basis": summary,
                },
            }
            for dataset in DATASETS for chrom in AUTOSOMES
        ]
        run_contract = {"runner": {"sha256": "a" * 64}}
        receipt = build_extraction_receipt(
            reusable,
            (1,),
            run_contract,
            {"max_new_tasks": 16, "reused_tasks": 152},
            1.0,
        )
        path = root / "full_extraction_receipt.json"
        receipt["receipt_integrity"] = {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{path.name}.sha256",
            "covers": path.name,
        }
        path.write_text(json.dumps(receipt), encoding="utf-8")
        write_sha256_sidecar(path)
        baseline = extraction_runner.validated_existing_final(
            path,
            run_contract,
            reusable_results=reusable,
            thresholds=(1,),
        )
        if baseline is None:
            raise AssertionError("synthetic terminal receipt did not pass the baseline deep validator")
        return path, receipt, reusable, run_contract

    @staticmethod
    def resign(path: pathlib.Path, receipt: dict) -> None:
        path.write_text(json.dumps(receipt), encoding="utf-8")
        write_sha256_sidecar(path)

    def test_terminal_resume_rejects_duplicate_for_missing_result_after_resign(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, receipt, reusable, run_contract = self.make_terminal(pathlib.Path(tmp))
            tampered = json.loads(json.dumps(receipt))
            tampered["results"][-1] = dict(tampered["results"][0])
            self.resign(path, tampered)
            self.assertIsNone(extraction_runner.validated_existing_final(
                path,
                run_contract,
                reusable_results=reusable,
                thresholds=(1,),
            ))

    def test_terminal_resume_rejects_aggregate_tamper_after_resign(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, receipt, reusable, run_contract = self.make_terminal(pathlib.Path(tmp))
            tampered = json.loads(json.dumps(receipt))
            tampered["aggregate"]["counts"]["synthetic_units"] = 10**18
            self.resign(path, tampered)
            self.assertIsNone(extraction_runner.validated_existing_final(
                path,
                run_contract,
                reusable_results=reusable,
                thresholds=(1,),
            ))


if __name__ == "__main__":
    unittest.main()
