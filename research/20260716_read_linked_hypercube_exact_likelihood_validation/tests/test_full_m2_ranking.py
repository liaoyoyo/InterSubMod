#!/usr/bin/env python3
import json
import argparse
import hashlib
import csv
import gzip
import pathlib
import sys
import tempfile
import threading
import time
import unittest
from unittest import mock


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))

from build_m2_patterns_and_rank import (  # noqa: E402
    CANDIDATE_FIELDS,
    RUNTIME_DIAGNOSTIC_FIELDS,
    sha256_path,
    summarize_runtime_values,
    write_sha256_sidecar,
)
import run_full_m2_ranking as ranking_runner  # noqa: E402
from test_full_m2_extraction import (  # noqa: E402
    make_release_contract,
    rewrite_release_manifest,
)
from run_full_m2_ranking import (  # noqa: E402
    CANDIDATE_TABLE_COLUMNS,
    EXPECTED_METHOD_CONTRACT,
    aggregate_conservation_audit,
    aggregate_primary_runtime_diagnostics,
    aggregate_rank_receipts,
    build_ranking_checkpoint,
    build_full_candidate_table,
    canonical_sort_results,
    create_batch_start_and_grants,
    ensure_release_orchestration_session,
    ensure_preflight_contract,
    load_full_extraction_receipt,
    load_release_contract_binding,
    load_release_orchestration_state,
    load_passing_rank_receipt,
    run_command_with_process_group_timeout,
    run_task,
    run_specs_rolling,
    require_matching_upstream_release_binding,
    scan_existing_rank_specs,
    semantic_json_sha256,
    validate_release_ranking_parameters,
    write_child_completion,
)


def make_ranking_release_contract(root: pathlib.Path):
    return make_release_contract(root)


class RankingResourceGateConflictTests(unittest.TestCase):
    def test_cooccurrence_audit_parent_and_worker_are_a_single_conflict_root(self):
        process_table = "\n".join(
            [
                "1 0 9999 0.0 /sbin/init",
                "100 1 40 0.0 python run_full_m2_ranking.py --workers 1",
                "101 100 39 1.0 python run_full_m2_ranking.py --workers 1",
                "200 1 300 2.0 python /x/audit_cooccurrence_task_contract_preflight.py --workers 40",
                "201 200 290 95.0 python /x/audit_cooccurrence_task_contract_preflight.pinned_abc.py --worker",
            ]
        )
        result = ranking_runner.active_conflicts(
            process_table=process_table,
            current_pid=101,
        )
        self.assertEqual(result["process_count"], 2)
        self.assertEqual(result["root_count"], 1)
        self.assertEqual(result["representatives"][0]["pid"], 200)
        self.assertEqual(result["representatives"][0]["member_pids"], [200, 201])
        self.assertEqual(
            result["representatives"][0]["group_kinds"],
            ["all_ssnv_cooccurrence_audit"],
        )


class ReleaseRankingBindingTests(unittest.TestCase):
    def args(self):
        return argparse.Namespace(
            thresholds=(1, 2, 3, 5), component_bases=("PS_HP1", "PS_HP2"),
            hp_families=("1", "2"), structural_exact_pattern_minreads=(1, 2, 3, 5),
            primary_structural_exact_pattern_minread=3, exact_k_max=12,
            max_vertex_sets=256, solver_time_limit_seconds=30.0,
            fixed_error_grid_values=(0.005, 0.01, 0.02, 0.05),
            minimum_bq_error_rate=1e-6, maximum_bq_error_rate=0.25,
            conditional_candidate_ranking_bootstrap_replicates=20,
            conditional_candidate_ranking_bootstrap_seed=20260716,
            tie_tolerance=1e-6, workers=2, task_timeout_seconds=28800.0,
            timeout_grace_seconds=300.0, max_new_tasks=16, aggregate_only=False,
        )

    def test_frozen_ranking_params_and_upstream_contract_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, paths = make_ranking_release_contract(pathlib.Path(tmp))
            binding = load_release_contract_binding(
                manifest,
                required_sources={
                    "full_ranking_runner": paths["full_ranking_runner"],
                    "ranker": paths["ranker"], "hypercube_solver": paths["hypercube_solver"],
                },
                _skip_deep_verification_for_test=True,
            )
            args = self.args()
            validate_release_ranking_parameters(binding, args)
            require_matching_upstream_release_binding(
                {"run_contract": {"release_binding": binding}}, binding
            )
            args.conditional_candidate_ranking_bootstrap_replicates = 19
            with self.assertRaisesRegex(RuntimeError, "parameters drift"):
                validate_release_ranking_parameters(binding, args)

    def test_extraction_and_ranking_different_contract_or_missing_binding_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest, paths = make_ranking_release_contract(pathlib.Path(tmp))
            binding = load_release_contract_binding(
                manifest, required_sources={"ranker": paths["ranker"]},
                _skip_deep_verification_for_test=True,
            )
            with self.assertRaisesRegex(RuntimeError, "requires a release-bound"):
                require_matching_upstream_release_binding({}, binding)
            other = json.loads(json.dumps(binding))
            other["release_manifest"]["sha256"] = "0" * 64
            with self.assertRaisesRegex(RuntimeError, "different release contracts"):
                require_matching_upstream_release_binding(
                    {"run_contract": {"release_binding": other}}, binding
                )


class ReleaseRankingOrchestrationTests(unittest.TestCase):
    @staticmethod
    def fixture(root: pathlib.Path):
        digest = lambda value: hashlib.sha256(value.encode()).hexdigest()
        binding = {
            "release_manifest": {
                "path": str(root / "m2_run_manifest.json"), "sha256": digest("release"),
                "semantic_sha256": digest("semantic"),
                "sidecar": {"path": str(root / "m2_run_manifest.json.sha256"), "sha256": digest("sidecar")},
            },
            "snapshot_sources": {
                "full_ranking_runner": {"path": str(root / "runner.py"), "sha256": digest("runner")},
                "ranker": {"path": str(root / "ranker.py"), "sha256": digest("ranker")},
                "hypercube_solver": {"path": str(root / "solver.py"), "sha256": digest("solver")},
            },
        }
        parent = {
            "session_id": digest("parent"),
            "terminal_receipt_path": str(root / "parent/full_extraction_receipt.json"),
            "terminal_receipt_sha256": digest("parent-terminal"),
        }
        run_contract = {
            "release_binding": binding,
            "orchestration_policy": dict(ranking_runner.ORCHESTRATION_POLICY),
        }
        parameters = {
            "structural_exact_pattern_minread_grid": [1, 2, 3, 5],
            "primary_structural_exact_pattern_minread": 3, "scoring_minread": 1,
            "exact_k_max": 12, "max_vertex_sets": 256,
            "solver_time_limit_seconds_per_milp": 30.0,
            "minimum_bq_error_rate": 1e-6, "maximum_bq_error_rate": 0.25,
            "fixed_error_grid_conditional_binary_flip_probability": [0.005, 0.01, 0.02, 0.05],
            "conditional_candidate_ranking_bootstrap_replicates": 20,
            "conditional_candidate_ranking_bootstrap_base_seed": 20260716,
            "tie_tolerance_log_likelihood": 1e-6,
        }
        return binding, parent, run_contract, parameters

    @staticmethod
    def create_session(root, binding, parent, run_contract):
        root.mkdir(parents=True)
        observed = root.stat()
        target = {
            "output_root": {
                "path": str(root.resolve()), "st_dev": int(observed.st_dev),
                "st_ino": int(observed.st_ino),
            },
            "release_manifest": dict(binding["release_manifest"]),
            "run_contract_semantic_sha256": semantic_json_sha256(run_contract),
            "parent_extraction": dict(parent),
        }
        _, gate = ranking_runner.create_resource_gate_receipt(
            root, stage="ranking", gate_scope="session", target=target,
            producer_source=ranking_runner.release_producer_sources(binding)["runner"],
            conflicts={"process_count": 0, "root_count": 0, "representatives": []},
        )
        return ensure_release_orchestration_session(
            root, binding, run_contract, parent, gate
        )

    @staticmethod
    def create_batch_gate(root, session, state, specs):
        batch_index = state["next_batch_index"]
        before_count = state["count"]
        max_new_tasks = 8 if before_count == 0 else 16
        target = {
            "output_root": dict(session["output_root"]),
            "session_id": session["session_id"],
            "session_sha256": sha256_path(root / "_orchestration/session.json"),
            "batch_index": batch_index,
            "before_count": before_count,
            "max_new_tasks": max_new_tasks,
            "effective_count": len(specs),
            "selected_task_ids": [f"{spec[0]}:{spec[1]}" for spec in specs],
            "previous_chain_head": state["previous_chain_head"],
        }
        _, gate = ranking_runner.create_resource_gate_receipt(
            root, stage="ranking", gate_scope="batch", batch_index=batch_index,
            target=target, producer_source=session["producer_sources"]["runner"],
            conflicts={"process_count": 0, "root_count": 0, "representatives": []},
        )
        return gate

    @staticmethod
    def specs(root: pathlib.Path, parameters: dict, start: int, count: int):
        pairs = [(dataset, chrom) for dataset in ranking_runner.DATASETS for chrom in ranking_runner.AUTOSOMES]
        rows = []
        for dataset, chrom in pairs[start:start + count]:
            extraction_dir = root / "extraction" / dataset / chrom
            extraction_dir.mkdir(parents=True, exist_ok=True)
            extraction_output = extraction_dir / "source.tsv.gz"
            extraction_output.write_bytes(b"source\n")
            extraction_receipt = {
                "outputs": {extraction_output.name: {
                    "path": str(extraction_output.resolve()),
                    "size_bytes": extraction_output.stat().st_size,
                    "sha256": sha256_path(extraction_output),
                }}
            }
            extraction_receipt_path = extraction_dir / "receipt.json"
            extraction_receipt_path.write_text(json.dumps(extraction_receipt), encoding="utf-8")
            write_sha256_sidecar(extraction_receipt_path)
            output_dir = root / "ranking" / "samples" / dataset / chrom
            command = ["/python", "/ranker.py", "--input-dir", str(extraction_dir.resolve()), "--output-dir", str(output_dir.resolve())]
            rows.append((
                dataset, chrom, extraction_dir, output_dir, command,
                root / "ranker.py", "r" * 64, parameters,
                ("PS_HP1", "PS_HP2"), (1, 2, 3, 5), ("1", "2"),
                "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY",
                semantic_json_sha256(extraction_receipt), False, 28800.0, 300.0,
            ))
        return rows

    @staticmethod
    def write_rank_child(spec):
        dataset, chrom, extraction_dir, output_dir = spec[:4]
        output_dir.mkdir(parents=True)
        output = output_dir / "rank.tsv.gz"
        output.write_bytes(b"rank\n")
        extraction_path = extraction_dir / "receipt.json"
        receipt = {
            "parameters": dict(spec[7]),
            "provenance": {"upstream_extraction_receipt": {
                "path": str(extraction_path.resolve()), "sha256": sha256_path(extraction_path),
            }},
            "aggregate_by_linkage_basis_threshold": {"PS_HP1": {"3": {}}},
            "outputs": {output.name: {
                "path": str(output.resolve()), "size_bytes": output.stat().st_size,
                "sha256": sha256_path(output),
            }},
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1", "sidecar_name": "receipt.json.sha256",
                "covers": "receipt.json",
            },
        }
        path = output_dir / "receipt.json"
        path.write_text(json.dumps(receipt), encoding="utf-8")
        write_sha256_sidecar(path)
        return receipt, path

    def publish_batch(self, root, session, run_contract, parameters, state, count):
        specs = self.specs(root, parameters, state["count"], count)
        gate = self.create_batch_gate(root / "ranking", session, state, specs)
        batch, grants = create_batch_start_and_grants(
            root / "ranking", session, run_contract, specs,
            before_count=state["count"], previous_chain_head=state["previous_chain_head"],
            batch_index=state["next_batch_index"], max_new_tasks=8 if state["count"] == 0 else 16,
            gate=gate,
        )
        raw = []
        for spec in reversed(specs):
            receipt, receipt_path = self.write_rank_child(spec)
            raw.append({
                "dataset": spec[0], "chrom": spec[1], "status": "PASS", "returncode": 0,
                "receipt": receipt, "receipt_path": str(receipt_path),
                "process_group_id": 10, "started_monotonic_ns": time.monotonic_ns(),
            })
        ordered = canonical_sort_results(raw)
        for result in ordered:
            task_id = f"{result['dataset']}:{result['chrom']}"
            result["orchestration_completion"] = write_child_completion(
                root / "ranking", session, batch, grants[task_id], result
            )
        cumulative = [
            {"dataset": task.split(":")[0], "chrom": task.split(":")[1],
             "status": "REUSED_PASS", "receipt": state["child_receipts"][task],
             "receipt_path": state["child_receipt_paths"][task],
             "orchestration_completion": identity}
            for task, identity in state["completions"].items()
        ] + ordered
        cumulative = canonical_sort_results(cumulative)
        n = len(cumulative)
        path = (
            root / "ranking/full_ranking_receipt.json"
            if n == ranking_runner.EXPECTED_TASKS
            else root / "ranking/checkpoints" / f"checkpoint_{n:03d}_of_154.json"
        )
        receipt = build_ranking_checkpoint(
            cumulative,
            run_contract=run_contract,
            invocation={
                "max_new_tasks": 8 if state["count"] == 0 else 16,
                "reused_tasks": state["count"],
                "selected_task_ids": [
                    f"{row['dataset']}:{row['chrom']}" for row in ordered
                ],
            },
            elapsed_seconds=1.0,
        )
        receipt["orchestration"] = {
            "session_identity": ranking_runner._session_identity(root / "ranking", session),
            "batch_start_identity": {"path": batch["path"], "sha256": batch["sha256"], "batch_id": batch["batch_id"], "batch_index": batch["batch_index"]},
            "previous_chain_head": state["previous_chain_head"],
            "batch_completion_attestations": [
                {"task_id": f"{row['dataset']}:{row['chrom']}", **row["orchestration_completion"]}
                for row in ordered
            ],
            "cumulative_attested_tasks": n,
        }
        receipt["receipt_integrity"] = {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{path.name}.sha256", "covers": path.name,
        }
        ranking_runner.write_immutable_json_exclusive(path, receipt)
        return load_release_orchestration_state(root / "ranking", session), raw, ordered

    def test_ranking_root_session_batch_policy_and_reverse_completion_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            binding, parent, run_contract, parameters = self.fixture(root)
            ranking_root = root / "ranking"
            session = self.create_session(ranking_root, binding, parent, run_contract)
            state = load_release_orchestration_state(ranking_root, session)
            bad_specs = self.specs(root, parameters, 0, 8)
            with self.assertRaisesRegex(RuntimeError, "illegal ranking release batch"):
                create_batch_start_and_grants(
                    ranking_root, session, run_contract, bad_specs, before_count=0,
                    previous_chain_head=None, batch_index=1, max_new_tasks=16,
                    gate={},
                )
            state, raw, ordered = self.publish_batch(root, session, run_contract, parameters, state, 8)
            self.assertEqual(state["count"], 8)
            self.assertNotEqual(
                [f"{r['dataset']}:{r['chrom']}" for r in raw],
                [f"{r['dataset']}:{r['chrom']}" for r in ordered],
            )
            state, _, _ = self.publish_batch(root, session, run_contract, parameters, state, 16)
            self.assertEqual(state["count"], 24)

    def test_production_ranking_helpers_generate_complete_154_task_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            binding, parent, run_contract, parameters = self.fixture(root)
            ranking_root = root / "ranking"
            session = self.create_session(ranking_root, binding, parent, run_contract)
            state = load_release_orchestration_state(ranking_root, session)
            for batch_size in (8, 16, 16, 16, 16, 16, 16, 16, 16, 16, 2):
                state, _, _ = self.publish_batch(
                    root, session, run_contract, parameters, state, batch_size
                )
            self.assertEqual(state["count"], ranking_runner.EXPECTED_TASKS)
            self.assertEqual(state["next_batch_index"], 12)
            self.assertTrue((ranking_root / "full_ranking_receipt.json").is_file())
            self.assertFalse(
                (ranking_root / "checkpoints/checkpoint_154_of_154.json").exists()
            )

    def test_ranking_resume_rejects_parent_change_and_upstream_child_drift(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            binding, parent, run_contract, parameters = self.fixture(root)
            ranking_root = root / "ranking"
            session = self.create_session(ranking_root, binding, parent, run_contract)
            state, _, _ = self.publish_batch(
                root, session, run_contract, parameters,
                load_release_orchestration_state(ranking_root, session), 8,
            )
            changed_parent = dict(parent)
            changed_parent["session_id"] = "0" * 64
            with self.assertRaisesRegex(RuntimeError, "parent mismatch"):
                ensure_release_orchestration_session(
                    ranking_root, binding, run_contract, changed_parent, None
                )
            extraction_receipt = root / "extraction/COLO829/chr1/receipt.json"
            extraction_receipt.write_text('{"tampered":true}\n', encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "upstream extraction identity"):
                load_release_orchestration_state(ranking_root, session)

    def test_ranking_existing_root_without_session_and_open_grant_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            binding, parent, run_contract, parameters = self.fixture(root)
            existing = root / "existing"
            existing.mkdir()
            with self.assertRaisesRegex(RuntimeError, "requires.*resource gate"):
                ensure_release_orchestration_session(
                    existing, binding, run_contract, parent, None
                )

            ranking_root = root / "ranking"
            session = self.create_session(ranking_root, binding, parent, run_contract)
            specs = self.specs(root, parameters, 0, 8)
            state = load_release_orchestration_state(ranking_root, session)
            gate = self.create_batch_gate(ranking_root, session, state, specs)
            create_batch_start_and_grants(
                ranking_root, session, run_contract, specs, before_count=0,
                previous_chain_head=None, batch_index=1, max_new_tasks=8,
                gate=gate,
            )
            with self.assertRaisesRegex(RuntimeError, "open/orphan"):
                load_release_orchestration_state(ranking_root, session)

    def test_ranking_renamed_and_resigned_checkpoint_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            binding, parent, run_contract, parameters = self.fixture(root)
            ranking_root = root / "ranking"
            session = self.create_session(ranking_root, binding, parent, run_contract)
            self.publish_batch(
                root, session, run_contract, parameters,
                load_release_orchestration_state(ranking_root, session), 8,
            )
            original = ranking_root / "checkpoints/checkpoint_008_of_154.json"
            renamed = ranking_root / "checkpoints/checkpoint_999_of_154.json"
            original.rename(renamed)
            original.with_name(f"{original.name}.sha256").rename(
                renamed.with_name(f"{renamed.name}.sha256")
            )
            payload = json.loads(renamed.read_text(encoding="utf-8"))
            payload["receipt_integrity"] = ranking_runner._receipt_integrity(renamed)
            renamed.chmod(0o644)
            renamed.write_text(json.dumps(payload), encoding="utf-8")
            renamed.chmod(0o444)
            sidecar = renamed.with_name(f"{renamed.name}.sha256")
            sidecar.chmod(0o644)
            sidecar.write_text(
                f"{sha256_path(renamed)}  {renamed.name}\n", encoding="ascii"
            )
            sidecar.chmod(0o444)
            with self.assertRaisesRegex(RuntimeError, "exact path/integrity"):
                load_release_orchestration_state(ranking_root, session)

    def test_ranking_resigned_checkpoint_aggregate_and_checks_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            binding, parent, run_contract, parameters = self.fixture(root)
            ranking_root = root / "ranking"
            session = self.create_session(ranking_root, binding, parent, run_contract)
            self.publish_batch(
                root, session, run_contract, parameters,
                load_release_orchestration_state(ranking_root, session), 8,
            )
            checkpoint = ranking_root / "checkpoints/checkpoint_008_of_154.json"
            payload = json.loads(checkpoint.read_text(encoding="utf-8"))
            payload["aggregate"]["ATTACKER_CONTROLLED"] = 10**18
            payload["checks"]["passing_count_matches_results"] = False
            rewrite_release_manifest(checkpoint, payload)
            with self.assertRaisesRegex(RuntimeError, "aggregate, checks"):
                load_release_orchestration_state(ranking_root, session)


class FullRankingReuseTests(unittest.TestCase):
    @staticmethod
    def write_full_extraction_receipt(root: pathlib.Path, results: list[dict]) -> pathlib.Path:
        root.mkdir(parents=True, exist_ok=True)
        path = root / "full_extraction_receipt.json"
        receipt = {
            "schema_name": "intersubmod.m2_full_extraction_receipt",
            "schema_version": "1.2.0",
            "scope": {
                "datasets": list(ranking_runner.DATASETS),
                "chromosomes": list(ranking_runner.AUTOSOMES),
                "expected_tasks": 154,
            },
            "n_results": 154,
            "results": results,
            "all_pass": True,
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": "full_extraction_receipt.json.sha256",
            },
        }
        path.write_text(json.dumps(receipt), encoding="utf-8")
        write_sha256_sidecar(path)
        return path

    def test_full_extraction_receipt_requires_exact_embedded_154_child_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            results = []
            for dataset in ranking_runner.DATASETS:
                for chrom in ranking_runner.AUTOSOMES:
                    results.append({
                        "dataset": dataset,
                        "chrom": chrom,
                        "status": "PASS",
                        "receipt": {
                            "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
                            "schema_version": "1.2.0",
                            "scope": {"dataset": dataset, "chrom": chrom},
                            "all_pass": True,
                        },
                    })
            self.write_full_extraction_receipt(root, results)
            loaded, _ = load_full_extraction_receipt(root)
            self.assertEqual(len(loaded["results"]), 154)

            self.write_full_extraction_receipt(root, [])
            with self.assertRaisesRegex(RuntimeError, "exactly 154 child"):
                load_full_extraction_receipt(root)

            duplicated = list(results)
            duplicated[-1] = dict(duplicated[0])
            self.write_full_extraction_receipt(root, duplicated)
            with self.assertRaisesRegex(RuntimeError, "exact canonical 154"):
                load_full_extraction_receipt(root)

    def fixture(self, root: pathlib.Path):
        extraction = root / "extraction"
        ranking = root / "ranking"
        ranker_source = root / "ranker.py"
        ranker_source.write_text("# fixture ranker source\n", encoding="utf-8")
        extraction.mkdir()
        ranking.mkdir()
        extraction_outputs = {}
        for name in (
            "D.chr1.molecule_sparse_calls.tsv.gz",
            "D.chr1.site_catalog.tsv.gz",
            "D.chr1.components.tsv.gz",
            "D.chr1.site_component_membership.tsv.gz",
        ):
            path = extraction / name
            path.write_bytes((name + "\n").encode())
            extraction_outputs[name] = {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_path(path),
            }
        extraction_receipt_path = extraction / "receipt.json"
        extraction_receipt = {
            "schema_name": "intersubmod.lossless_read_linkage_chromosome_receipt",
            "schema_version": "1.2.0",
            "scope": {"dataset": "COLO829", "chrom": "chr1"},
            "outputs": extraction_outputs,
            "checks": {"synthetic": True},
            "all_pass": True,
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": "receipt.json.sha256",
            },
        }
        extraction_receipt_path.write_text(json.dumps(extraction_receipt), encoding="utf-8")
        write_sha256_sidecar(extraction_receipt_path)

        rank_output = ranking / "units.tsv.gz"
        rank_output.write_bytes(b"ranked\n")
        runtime_output = ranking / "m2_unit_runtime_diagnostics.tsv.gz"
        runtime_row = {
            "dataset": "COLO829", "chrom": "chr1", "threshold": 3,
            "component_basis": "PS_HP1", "phase_set": "100", "component_id": "C",
            "family": "1", "structural_exact_pattern_minread": 3,
            "structural_minread_role": "PRIMARY",
            "candidate_generation_invoked": "true",
            "likelihood_fit_invoked": "true",
            "candidate_generation_elapsed_seconds": 0.1,
            "likelihood_fit_elapsed_seconds": 0.2,
            "unit_total_elapsed_seconds": 0.5,
        }
        with gzip.open(runtime_output, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, RUNTIME_DIAGNOSTIC_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerow(runtime_row)
        timing_scope = {
            "n_unit_evaluations": 1,
            **{
                metric: summarize_runtime_values([runtime_row[metric]])
                for metric in (
                    "candidate_generation_elapsed_seconds",
                    "likelihood_fit_elapsed_seconds",
                    "unit_total_elapsed_seconds",
                )
            },
        }
        rank_receipt_path = ranking / "receipt.json"
        rank_receipt = {
            "schema_name": "intersubmod.m2_symbolic_patterns_vertex_rank_receipt",
            "schema_version": "2.0.0",
            "scope": {
                "dataset": ["COLO829"],
                "chrom": ["chr1"],
                "component_bases": ["PS_HP1", "PS_HP2"],
                "thresholds": [1, 2, 3, 5],
                "hp_families": ["1", "2"],
                "component_basis_mode": "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY",
            },
            "parameters": {
                "primary_structural_exact_pattern_minread": 3,
                "method_contract": EXPECTED_METHOD_CONTRACT,
            },
            "provenance": {
                "ranker": {
                    "path": str(ranker_source.resolve()),
                    "sha256": sha256_path(ranker_source),
                },
                "upstream_extraction_receipt": {
                    "path": str(extraction_receipt_path.resolve()),
                    "sha256": sha256_path(extraction_receipt_path),
                },
            },
            "outputs": {
                rank_output.name: {
                    "path": str(rank_output.resolve()),
                    "size_bytes": rank_output.stat().st_size,
                    "sha256": sha256_path(rank_output),
                },
                runtime_output.name: {
                    "path": str(runtime_output.resolve()),
                    "size_bytes": runtime_output.stat().st_size,
                    "sha256": sha256_path(runtime_output),
                },
            },
            "runtime_diagnostics": {
                "schema_name": "intersubmod.m2_unit_runtime_diagnostics",
                "schema_version": "1.0.0",
                "clock": "time.monotonic_ns",
                "unit": "seconds",
                "per_unit_output": runtime_output.name,
                "scopes": {
                    "primary_unit_evaluations": timing_scope,
                    "all_structural_minread_unit_evaluations": timing_scope,
                },
                "primary_invoked_segment_scopes": {
                    metric: timing_scope[metric]
                    for metric in (
                        "candidate_generation_elapsed_seconds",
                        "likelihood_fit_elapsed_seconds",
                    )
                },
            },
            "all_pass": True,
            "receipt_integrity": {
                "scheme": "external_sha256_sidecar_v1",
                "sidecar_name": "receipt.json.sha256",
            },
        }
        rank_receipt_path.write_text(json.dumps(rank_receipt), encoding="utf-8")
        write_sha256_sidecar(rank_receipt_path)
        return extraction_receipt_path, extraction_receipt, rank_receipt_path, extraction / "D.chr1.molecule_sparse_calls.tsv.gz"

    def load(self, extraction_path, extraction_receipt, rank_path, **overrides):
        ranker_path = rank_path.parent.parent / "ranker.py"
        kwargs = {
            "dataset": "COLO829",
            "chrom": "chr1",
            "extraction_receipt_path": extraction_path,
            "ranker_path": ranker_path,
            "ranker_sha256": sha256_path(ranker_path),
            "expected_parameters": {"primary_structural_exact_pattern_minread": 3},
            "expected_bases": ("PS_HP1", "PS_HP2"),
            "expected_thresholds": (1, 2, 3, 5),
            "expected_hp_families": ("1", "2"),
            "expected_component_basis_mode": "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY",
            "expected_extraction_receipt_semantic_sha256": semantic_json_sha256(extraction_receipt),
        }
        kwargs.update(overrides)
        return load_passing_rank_receipt(rank_path, **kwargs)

    @staticmethod
    def rewrite_receipt(path: pathlib.Path, receipt: dict) -> None:
        path.write_text(json.dumps(receipt), encoding="utf-8")
        write_sha256_sidecar(path)

    def test_valid_reuse_and_scope_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            extraction_path, extraction_receipt, rank_path, _ = self.fixture(pathlib.Path(tmp))
            self.assertIsNotNone(self.load(extraction_path, extraction_receipt, rank_path))
            self.assertIsNone(self.load(
                extraction_path, extraction_receipt, rank_path, expected_hp_families=("1",)
            ))
            self.assertIsNone(self.load(
                extraction_path,
                extraction_receipt,
                rank_path,
                expected_component_basis_mode="EXPLICIT_USER_SELECTION_NONPRIMARY_OR_SENSITIVITY",
            ))

    def test_tampered_extraction_tsv_invalidates_reuse_even_when_receipt_is_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            extraction_path, extraction_receipt, rank_path, calls = self.fixture(pathlib.Path(tmp))
            calls.write_bytes(b"tampered upstream without changing receipt\n")
            self.assertIsNone(self.load(extraction_path, extraction_receipt, rank_path))

    def test_ranker_source_mutation_invalidates_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            extraction_path, extraction_receipt, rank_path, _ = self.fixture(pathlib.Path(tmp))
            ranker_path = rank_path.parent.parent / "ranker.py"
            declared_sha = sha256_path(ranker_path)
            ranker_path.write_text("# source mutated after child receipt\n", encoding="utf-8")
            self.assertIsNone(
                self.load(
                    extraction_path,
                    extraction_receipt,
                    rank_path,
                    ranker_sha256=declared_sha,
                )
            )

    def test_missing_or_tampered_runtime_diagnostic_invalidates_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            extraction_path, extraction_receipt, rank_path, _ = self.fixture(pathlib.Path(tmp))
            runtime_path = rank_path.parent / "m2_unit_runtime_diagnostics.tsv.gz"
            runtime_path.write_bytes(b"not a gzip runtime table")
            self.assertIsNone(self.load(extraction_path, extraction_receipt, rank_path))

    def test_old_key_only_cannot_substitute_for_method_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            extraction_path, extraction_receipt, rank_path, _ = self.fixture(pathlib.Path(tmp))
            receipt = json.loads(rank_path.read_text(encoding="utf-8"))
            receipt["parameters"].pop("method_contract")
            receipt.setdefault("checks", {})[
                "no_" + "partial_completions_materialized"
            ] = True
            self.rewrite_receipt(rank_path, receipt)
            self.assertIsNone(self.load(extraction_path, extraction_receipt, rank_path))

    def test_method_contract_drift_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            extraction_path, extraction_receipt, rank_path, _ = self.fixture(pathlib.Path(tmp))
            receipt = json.loads(rank_path.read_text(encoding="utf-8"))
            receipt["parameters"]["method_contract"]["structural_group_scope"] = "DRIFT"
            self.rewrite_receipt(rank_path, receipt)
            self.assertIsNone(self.load(extraction_path, extraction_receipt, rank_path))

    def test_vaf_or_parent_edge_scoring_true_is_rejected(self):
        for field in (
            "same_molecule_vaf_added_as_separate_term",
            "parent_edges_or_transitions_scored",
        ):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmp:
                extraction_path, extraction_receipt, rank_path, _ = self.fixture(pathlib.Path(tmp))
                receipt = json.loads(rank_path.read_text(encoding="utf-8"))
                receipt["parameters"]["method_contract"][field] = True
                self.rewrite_receipt(rank_path, receipt)
                self.assertIsNone(self.load(extraction_path, extraction_receipt, rank_path))

    def test_full_runtime_aggregate_streams_rows_and_uses_exact_nearest_rank(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            results = []
            expected_candidate = []
            expected_likelihood = []
            expected_total = []
            for index, values in enumerate(((0.1, 0.2, 0.5), (0.4, 0.8, 1.5)), start=1):
                fixture_root = root / str(index)
                fixture_root.mkdir()
                extraction_path, extraction_receipt, rank_path, _ = self.fixture(fixture_root)
                receipt = json.loads(rank_path.read_text(encoding="utf-8"))
                runtime_path = rank_path.parent / "m2_unit_runtime_diagnostics.tsv.gz"
                with gzip.open(runtime_path, "rt", encoding="utf-8", newline="") as handle:
                    row = next(csv.DictReader(handle, delimiter="\t"))
                for metric, value in zip(
                    (
                        "candidate_generation_elapsed_seconds",
                        "likelihood_fit_elapsed_seconds",
                        "unit_total_elapsed_seconds",
                    ),
                    values,
                ):
                    row[metric] = value
                with gzip.open(runtime_path, "wt", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, RUNTIME_DIAGNOSTIC_FIELDS, delimiter="\t")
                    writer.writeheader()
                    writer.writerow(row)
                identity = receipt["outputs"][runtime_path.name]
                identity.update({
                    "size_bytes": runtime_path.stat().st_size,
                    "sha256": sha256_path(runtime_path),
                })
                scope = {
                    "n_unit_evaluations": 1,
                    **{
                        metric: summarize_runtime_values([float(row[metric])])
                        for metric in (
                            "candidate_generation_elapsed_seconds",
                            "likelihood_fit_elapsed_seconds",
                            "unit_total_elapsed_seconds",
                        )
                    },
                }
                receipt["runtime_diagnostics"]["scopes"] = {
                    "primary_unit_evaluations": scope,
                    "all_structural_minread_unit_evaluations": scope,
                }
                receipt["runtime_diagnostics"]["primary_invoked_segment_scopes"] = {
                    metric: scope[metric]
                    for metric in (
                        "candidate_generation_elapsed_seconds",
                        "likelihood_fit_elapsed_seconds",
                    )
                }
                results.append({
                    "dataset": f"D{index}", "chrom": "chr1", "status": "PASS",
                    "receipt": receipt,
                })
                expected_candidate.append(values[0])
                expected_likelihood.append(values[1])
                expected_total.append(values[2])
            aggregate = aggregate_primary_runtime_diagnostics(results)
            self.assertEqual(aggregate["n_child_runtime_files"], 2)
            self.assertEqual(aggregate["n_unit_evaluations"], 2)
            self.assertEqual(
                aggregate["metrics"]["candidate_generation_elapsed_seconds"],
                summarize_runtime_values(expected_candidate),
            )
            self.assertEqual(
                aggregate["metrics"]["likelihood_fit_elapsed_seconds"],
                summarize_runtime_values(expected_likelihood),
            )
            self.assertEqual(
                aggregate["metrics"]["unit_total_elapsed_seconds"],
                summarize_runtime_values(expected_total),
            )
            self.assertEqual(
                aggregate["metrics_when_invoked"]["candidate_generation_elapsed_seconds"],
                summarize_runtime_values(expected_candidate),
            )

    def test_full_candidate_table_is_primary_ps_aware_and_canonical(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            child = root / "child.tsv.gz"
            source_row = {
                "dataset": "COLO829", "chrom": "chr1", "component_basis": "PS_HP1",
                "phase_set": "100", "threshold": "3", "component_id": "C", "family": "1",
                "structural_exact_pattern_minread": "3", "vertex_set_id": "v1",
                "states_json": json.dumps([
                    {"bitmask": 0, "state_rax": "R", "roles": ["root"]},
                    {"bitmask": 1, "state_rax": "A", "roles": ["full-observed"]},
                ]),
                "parent_choice_count": "1", "unique_parent_edges_json": "[]",
                "primary_log_likelihood": "-2.0", "relative_likelihood_weight": "1.0",
                "mixture_weights_json": json.dumps({"0": 0.4, "1": 0.6}),
                "fit_converged": "true", "fit_monotone": "true",
                "coarse_topology_classes_json": json.dumps(["single-only"]),
                "is_winner": "true", "is_tied_winner": "false",
                "candidate_semantic_sha256": "x" * 64,
            }
            with gzip.open(child, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, CANDIDATE_FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerow(source_row)
            results = [{
                "dataset": "COLO829", "chrom": "chr1", "status": "PASS",
                "receipt": {"outputs": {"m2_compressed_vertex_set_candidates.tsv.gz": {
                    "path": str(child), "size_bytes": child.stat().st_size,
                    "sha256": sha256_path(child),
                }}},
            }]
            metadata = build_full_candidate_table(results, root, 3)
            self.assertEqual(metadata["n_rows"], 1)
            self.assertEqual(metadata["n_units"], 1)
            self.assertTrue(ranking_runner._verify_existing_candidate_table(
                metadata, root, results=results, primary_minread=3,
                require_full_scope=False,
            ))
            with gzip.open(metadata["path"], "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(tuple(reader.fieldnames), CANDIDATE_TABLE_COLUMNS)
                row = next(reader)
            self.assertEqual(row["ps"], "100")
            self.assertEqual(row["winner_status"], "UNIQUE_WINNER")

    def test_candidate_table_marks_entire_unit_abstain_when_any_fit_is_invalid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            child = root / "child.tsv.gz"
            common = {
                "dataset": "COLO829", "chrom": "chr1", "component_basis": "PS_HP1",
                "phase_set": "100", "threshold": "3", "component_id": "C", "family": "1",
                "structural_exact_pattern_minread": "3",
                "states_json": json.dumps([
                    {"bitmask": 0, "state_rax": "R", "roles": ["root"]},
                    {"bitmask": 1, "state_rax": "A", "roles": ["full-observed"]},
                ]),
                "parent_choice_count": "1", "unique_parent_edges_json": "[]",
                "relative_likelihood_weight": "0.5",
                "mixture_weights_json": json.dumps({"0": 0.4, "1": 0.6}),
                "coarse_topology_classes_json": json.dumps(["single-only"]),
                "is_winner": "false", "is_tied_winner": "false",
                "candidate_semantic_sha256": "x" * 64,
            }
            source_rows = [
                common | {
                    "vertex_set_id": "v1", "primary_log_likelihood": "-2.0",
                    "fit_converged": "true", "fit_monotone": "true",
                },
                common | {
                    "vertex_set_id": "v2", "primary_log_likelihood": "-3.0",
                    "fit_converged": "false", "fit_monotone": "true",
                },
            ]
            with gzip.open(child, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, CANDIDATE_FIELDS, delimiter="\t")
                writer.writeheader()
                writer.writerows(source_rows)
            results = [{
                "dataset": "COLO829", "chrom": "chr1", "status": "PASS",
                "receipt": {"outputs": {"m2_compressed_vertex_set_candidates.tsv.gz": {
                    "path": str(child), "size_bytes": child.stat().st_size,
                    "sha256": sha256_path(child),
                }}},
            }]
            metadata = build_full_candidate_table(results, root, 3)
            with gzip.open(metadata["path"], "rt", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["winner_status"] for row in rows}, {"ABSTAIN_UNIT_OPTIMIZER"})

    def test_ranking_rolling_queue_never_submits_more_than_workers(self):
        lock = threading.Lock()
        running = 0
        max_running = 0

        def task_runner(spec):
            nonlocal running, max_running
            with lock:
                running += 1
                max_running = max(max_running, running)
            time.sleep(0.01)
            with lock:
                running -= 1
            return {"dataset": "D", "chrom": f"chr{spec}", "status": "PASS"}

        results, max_inflight = run_specs_rolling(
            list(range(20)), workers=3, task_runner=task_runner
        )
        self.assertEqual(len(results), 20)
        self.assertLessEqual(max_inflight, 3)
        self.assertLessEqual(max_running, 3)

    def test_same_done_batch_failure_does_not_submit_replacement(self):
        started = []
        barrier = threading.Barrier(2)
        specs = [(index,) for index in range(3)]

        def task_runner(spec):
            index = spec[0]
            started.append(index)
            if index < 2:
                barrier.wait(timeout=2)
            return {"dataset": "D", "chrom": f"chr{index}", "status": "FAIL" if index == 1 else "PASS"}

        original_wait = ranking_runner.concurrent.futures.wait

        def wait_all(futures, **_kwargs):
            return original_wait(futures, return_when=ranking_runner.concurrent.futures.ALL_COMPLETED)

        with mock.patch.object(ranking_runner.concurrent.futures, "wait", side_effect=wait_all):
            results, max_inflight = run_specs_rolling(specs, workers=2, task_runner=task_runner)
        self.assertEqual(set(started), {0, 1})
        self.assertNotIn(2, started)
        self.assertEqual(max_inflight, 2)
        self.assertTrue(any(row["status"] == "FAIL" for row in results))

    def test_ranking_timeout_uses_new_process_group_and_kill_fallback(self):
        code = "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(60)"
        result = run_command_with_process_group_timeout(
            [sys.executable, "-c", code],
            task_timeout_seconds=0.2,
            timeout_grace_seconds=0.1,
        )
        self.assertTrue(result["timed_out"])
        self.assertEqual(result["termination_stage"], "KILL")

    def test_ranking_timeout_marker_preserves_diagnostic_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            extraction_path, extraction_receipt, _, _ = self.fixture(root)
            output_dir = root / "timed_out_ranking"
            ranker_path = root / "ranker.py"
            code = "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(60)"
            spec = (
                "COLO829", "chr1", extraction_path.parent, output_dir,
                [sys.executable, "-c", code], ranker_path, sha256_path(ranker_path),
                {"primary_structural_exact_pattern_minread": 3},
                ("PS_HP1", "PS_HP2"), (1, 2, 3, 5), ("1", "2"),
                "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY",
                semantic_json_sha256(extraction_receipt), False, 0.2, 0.1,
            )
            result = run_task(spec)
            self.assertEqual(result["status"], "TIMEOUT")
            marker = json.loads((output_dir / "runner_task_status.json").read_text())
            self.assertEqual(marker["returncode"], result["returncode"])
            self.assertGreaterEqual(marker["elapsed_seconds"], 0.2)
            self.assertIn("stdout_tail", marker)
            self.assertIn("stderr_tail", marker)

    def test_ranking_runner_sha_drift_is_resume_hard_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "preflight.json"
            contract = {
                "runner": {"path": "/rank-runner.py", "sha256": "a" * 64},
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
            drift = contract | {
                "runner": {"path": "/rank-runner.py", "sha256": "b" * 64}
            }
            with self.assertRaises(RuntimeError):
                ensure_preflight_contract(path, {"run_contract": drift}, drift)

    def test_ranking_existing_timeout_marker_is_nonpass_and_checkpoint_is_not_final(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            extraction_path, extraction_receipt, rank_path, _ = self.fixture(root)
            marker = rank_path.parent / "runner_task_status.json"
            marker.write_text(json.dumps({"status": "TIMEOUT"}), encoding="utf-8")
            ranker_path = root / "ranker.py"
            spec = (
                "COLO829",
                "chr1",
                extraction_path.parent,
                rank_path.parent,
                ["unused"],
                ranker_path,
                sha256_path(ranker_path),
                {"primary_structural_exact_pattern_minread": 3},
                ("PS_HP1", "PS_HP2"),
                (1, 2, 3, 5),
                ("1", "2"),
                "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY",
                semantic_json_sha256(extraction_receipt),
                False,
                28800.0,
                300.0,
            )
            reused, pending, invalid = scan_existing_rank_specs([spec])
            self.assertEqual((reused, pending), ([], []))
            self.assertEqual(invalid[0]["status"], "FAIL_EXISTING_NONPASS_DIRECTORY")

    def test_non_timeout_ranking_failure_is_durable_and_not_implicitly_retried(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            extraction_path, extraction_receipt, _, _ = self.fixture(root)
            output_dir = root / "failed_ranking"
            ranker_path = root / "ranker.py"
            spec = (
                "COLO829", "chr1", extraction_path.parent, output_dir,
                [sys.executable, "-c", "import sys; print('rank bad'); sys.exit(9)"],
                ranker_path, sha256_path(ranker_path),
                {"primary_structural_exact_pattern_minread": 3},
                ("PS_HP1", "PS_HP2"), (1, 2, 3, 5), ("1", "2"),
                "PS_AWARE_HP_FAMILY_X_KNOWN_PS_PRIMARY",
                semantic_json_sha256(extraction_receipt), False, 10.0, 1.0,
            )
            result = run_task(spec)
            self.assertEqual((result["status"], result["returncode"]), ("FAIL", 9))
            marker = json.loads((output_dir / "runner_task_status.json").read_text())
            self.assertEqual(marker["status"], "CHILD_FAILED_OR_INVALID_RECEIPT")
            self.assertIn("rank bad", marker["stdout_tail"])
            reused, pending, invalid = scan_existing_rank_specs([spec])
            self.assertEqual((reused, pending), ([], []))
            self.assertEqual(invalid[0]["status"], "FAIL_EXISTING_NONPASS_DIRECTORY")

    def test_partial_ranking_payload_is_checkpoint_only_and_includes_reused(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            results = []
            for index in range(2):
                fixture_root = root / str(index)
                fixture_root.mkdir()
                _, _, rank_path, _ = self.fixture(fixture_root)
                receipt = json.loads(rank_path.read_text(encoding="utf-8"))
                results.append({
                    "dataset": "COLO829",
                    "chrom": f"chr{index + 1}",
                    "status": "REUSED_PASS",
                    "receipt": receipt,
                    "receipt_path": str(rank_path),
                    "orchestration_completion": {
                        "path": f"/completion/{index}.json", "sha256": str(index) * 64,
                    },
                })
            with mock.patch.object(
                ranking_runner,
                "aggregate_conservation_audit",
                return_value={"all_aggregate_cells_conserved": True, "failed_cells": []},
            ):
                checkpoint = build_ranking_checkpoint(
                    results,
                    run_contract={"runner": {"sha256": "a" * 64}},
                    invocation={"max_new_tasks": 2, "reused_tasks": 2},
                    elapsed_seconds=0.1,
                )
            self.assertEqual(
                checkpoint["schema_name"], "intersubmod.m2_full_ranking_checkpoint"
            )
            self.assertFalse(checkpoint["all_pass"])
            self.assertTrue(checkpoint["checkpoint_complete"])
            self.assertEqual(checkpoint["passing_tasks"], 2)
            self.assertEqual(checkpoint["remaining_tasks"], 152)
            self.assertEqual(checkpoint["task_status_counts"], {"REUSED_PASS": 2})
            self.assertEqual(
                checkpoint["results"][0]["orchestration_completion"],
                results[0]["orchestration_completion"],
            )

            with mock.patch.object(
                ranking_runner,
                "aggregate_conservation_audit",
                return_value={"all_aggregate_cells_conserved": False, "failed_cells": ["x"]},
            ):
                rejected = build_ranking_checkpoint(
                    results,
                    run_contract={"runner": {"sha256": "a" * 64}},
                    invocation={"max_new_tasks": 2},
                    elapsed_seconds=0.1,
                )
            self.assertFalse(rejected["checkpoint_complete"])
            self.assertFalse(rejected["checks"]["recorded_aggregate_cells_conserve"])

    def test_154_child_rank_aggregate_conservation_does_not_regress(self):
        summary = {
            "n_component_hp_units": 1,
            "molecule_component_projections": 1,
            "informative_scoring_molecules": 1,
            "all_x_excluded_molecules": 0,
            "structural_retained_molecules": 1,
            "below_minread_scoring_molecules": 0,
            "bq_scoring_molecules": 1,
            "raw_tree_candidates_T_complete_units": 1,
            "distinct_vertex_sets_V_complete_units": 1,
            "solver_complete_units": 1,
            "solver_incomplete_or_not_run_units": 0,
            "quality_primary_unique_vertex_units": 1,
            "quality_primary_tied_vertex_units": 0,
            "rank_abstain_units": 0,
            "fixed_error_grid_stable_units": 1,
            "fixed_error_grid_evaluated_units": 1,
            "coarse_topology_class_unique_units": 1,
            "coarse_topology_multiple_class_units": 0,
            "topology_evaluated_units": 1,
            "partial_group_coverage_denominator": 0,
            "partial_groups_covered": 0,
            "partial_groups_unsatisfied": 0,
            "selection_status_counts": {"UNIQUE_WINNER": 1},
            "candidate_generation_status_counts": {"COMPLETE": 1},
            "k_route_counts": {"EXACT_K_LE_12": 1},
        }
        results = []
        datasets = (
            "COLO829", "H1437", "H2009", "HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
        )
        for dataset in datasets:
            for chrom_index in range(1, 23):
                results.append({
                    "dataset": dataset,
                    "chrom": f"chr{chrom_index}",
                    "status": "REUSED_PASS",
                    "receipt": {
                        "input_counts": {
                            "sparse_molecule_rows_total": 1,
                            "sparse_molecule_rows_included_in_at_least_one_selected_unit": 1,
                        },
                        "aggregate_by_linkage_basis_threshold": {
                            "PS_HP1": {"3": summary}
                        },
                        "partial_pattern_funnel_by_linkage_basis_threshold": {},
                        "sensitivity_by_structural_exact_pattern_minread": {},
                    },
                })
        aggregate = aggregate_rank_receipts(results)
        audit = aggregate_conservation_audit(aggregate)
        self.assertEqual(aggregate["task_status_counts"], {"REUSED_PASS": 154})
        self.assertEqual(set(aggregate["by_dataset"]), set(datasets))
        self.assertEqual(
            aggregate["by_linkage_basis_threshold"]["PS_HP1"]["3"]["n_component_hp_units"],
            154,
        )
        self.assertTrue(audit["all_aggregate_cells_conserved"])


class TerminalRankingResumeDeepValidationTests(unittest.TestCase):
    @staticmethod
    def synthetic_aggregate():
        return {
            "task_status_counts": {"REUSED_PASS": 154},
            "input_call_funnel": {},
            "input_sparse_call_code_counts": {},
            "by_linkage_basis_threshold": {},
            "by_linkage_basis_and_threshold": {
                "PS_HP1": {"3": {
                    "solver_complete_units": 1,
                    "n_component_hp_units": 1,
                }}
            },
            "by_dataset": {dataset: {} for dataset in ranking_runner.DATASETS},
            "by_structural_exact_pattern_minread": {},
        }

    @staticmethod
    def synthetic_runtime():
        return {
            "n_child_runtime_files": 154,
            "n_unit_evaluations": 1,
            "all_child_summaries_recomputed_from_per_unit_rows": True,
        }

    @staticmethod
    def validate(path, run_contract, reusable, output_root, extraction_path):
        aggregate = TerminalRankingResumeDeepValidationTests.synthetic_aggregate()
        conservation = {"all_aggregate_cells_conserved": True, "failed_cells": []}
        runtime = TerminalRankingResumeDeepValidationTests.synthetic_runtime()
        with (
            mock.patch.object(ranking_runner, "aggregate_rank_receipts", return_value=aggregate),
            mock.patch.object(
                ranking_runner, "aggregate_conservation_audit", return_value=conservation
            ),
            mock.patch.object(
                ranking_runner,
                "aggregate_primary_runtime_diagnostics",
                return_value=runtime,
            ),
            mock.patch.object(ranking_runner, "child_runtime_contract_valid", return_value=True),
            mock.patch.object(
                ranking_runner,
                "child_method_contract_and_ranker_source_bound",
                return_value=True,
            ),
        ):
            return ranking_runner.validated_existing_final(
                path,
                run_contract,
                reusable_results=reusable,
                output_root=output_root,
                full_extraction_path=extraction_path,
                ranker_sha256=sha256_path(ranking_runner.RANKER),
            )

    @classmethod
    def make_terminal(cls, root: pathlib.Path):
        output_root = root / "ranking"
        extraction_path = root / "extraction" / "full_extraction_receipt.json"
        extraction_path.parent.mkdir(parents=True)
        extraction_path.write_text("{}\n", encoding="utf-8")
        reusable = []
        child_receipt_base = {
            "schema_version": "2.0.0",
            "parameters": {"method_contract": EXPECTED_METHOD_CONTRACT},
            "provenance": {
                "ranker": {
                    "path": str(ranking_runner.RANKER.resolve()),
                    "sha256": sha256_path(ranking_runner.RANKER),
                }
            },
            "checks": {
                "no_cross_ps_pattern_pooling": True,
                "known_ps_never_mixed": True,
                "missing_ps_separate_diagnostic": True,
            },
        }
        source_row = {field: "" for field in CANDIDATE_FIELDS}
        source_row.update({
            "dataset": "COLO829", "chrom": "chr1",
            "component_basis": "PS_HP1", "phase_set": "100",
            "threshold": "3", "component_id": "C1", "family": "1",
            "structural_exact_pattern_minread": "3", "vertex_set_id": "v1",
            "states_json": json.dumps([
                {"bitmask": 0, "state_rax": "R", "roles": ["root"]},
                {"bitmask": 1, "state_rax": "A", "roles": ["full-observed"]},
            ]),
            "parent_choice_count": "1", "primary_log_likelihood": "-2.0",
            "mixture_weights_json": json.dumps({"0": 0.4, "1": 0.6}),
            "fit_converged": "true", "fit_monotone": "true",
            "coarse_topology_classes_json": json.dumps(["single-only"]),
            "is_winner": "true", "is_tied_winner": "false",
        })
        for dataset in ranking_runner.DATASETS:
            for chrom in ranking_runner.AUTOSOMES:
                receipt_path = output_root / "samples" / dataset / chrom / "receipt.json"
                receipt_path.parent.mkdir(parents=True)
                child_candidate = receipt_path.parent / "m2_compressed_vertex_set_candidates.tsv.gz"
                with gzip.open(child_candidate, "wt", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, CANDIDATE_FIELDS, delimiter="\t")
                    writer.writeheader()
                    if dataset == "COLO829" and chrom == "chr1":
                        writer.writerow(source_row)
                child_receipt = json.loads(json.dumps(child_receipt_base))
                child_receipt["outputs"] = {
                    "m2_compressed_vertex_set_candidates.tsv.gz": {
                        "path": str(child_candidate.resolve()),
                        "size_bytes": child_candidate.stat().st_size,
                        "sha256": sha256_path(child_candidate),
                    }
                }
                receipt_path.write_text(json.dumps(child_receipt), encoding="utf-8")
                reusable.append({
                    "dataset": dataset,
                    "chrom": chrom,
                    "status": "REUSED_PASS",
                    "receipt": child_receipt,
                    "receipt_path": str(receipt_path),
                })

        candidate_table = build_full_candidate_table(reusable, output_root, 3)
        aggregate = cls.synthetic_aggregate()
        conservation = {"all_aggregate_cells_conserved": True, "failed_cells": []}
        runtime = cls.synthetic_runtime()
        run_contract = {
            "full_extraction_receipt": {
                "path": str(extraction_path.resolve()),
                "sha256": sha256_path(extraction_path),
            },
            "parameters": {"primary_structural_exact_pattern_minread": 3},
        }
        with (
            mock.patch.object(ranking_runner, "child_runtime_contract_valid", return_value=True),
            mock.patch.object(
                ranking_runner,
                "child_method_contract_and_ranker_source_bound",
                return_value=True,
            ),
        ):
            receipt = ranking_runner.build_ranking_terminal_receipt(
                reusable,
                output_root=output_root,
                full_extraction_path=extraction_path,
                run_contract=run_contract,
                invocation={"max_new_tasks": 16, "reused_tasks": 152},
                elapsed_seconds=1.0,
                aggregate=aggregate,
                conservation_audit=conservation,
                candidate_table=candidate_table,
                runtime_diagnostics=runtime,
                ranker_sha256=sha256_path(ranking_runner.RANKER),
                release_orchestration=None,
            )
        path = output_root / "full_ranking_receipt.json"
        receipt["receipt_integrity"] = {
            "scheme": "external_sha256_sidecar_v1",
            "sidecar_name": f"{path.name}.sha256",
            "covers": path.name,
        }
        path.write_text(json.dumps(receipt), encoding="utf-8")
        write_sha256_sidecar(path)
        if cls.validate(path, run_contract, reusable, output_root, extraction_path) is None:
            raise AssertionError("synthetic ranking terminal did not pass baseline deep validation")
        return path, receipt, reusable, run_contract, output_root, extraction_path

    @staticmethod
    def resign(path: pathlib.Path, receipt: dict) -> None:
        path.write_text(json.dumps(receipt), encoding="utf-8")
        write_sha256_sidecar(path)

    def test_terminal_resume_rejects_duplicate_for_missing_result_after_resign(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, receipt, reusable, contract, output_root, extraction = self.make_terminal(
                pathlib.Path(tmp)
            )
            tampered = json.loads(json.dumps(receipt))
            tampered["results"][-1] = dict(tampered["results"][0])
            self.resign(path, tampered)
            self.assertIsNone(self.validate(path, contract, reusable, output_root, extraction))

    def test_terminal_resume_rejects_aggregate_tamper_after_resign(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, receipt, reusable, contract, output_root, extraction = self.make_terminal(
                pathlib.Path(tmp)
            )
            tampered = json.loads(json.dumps(receipt))
            tampered["aggregate"]["ATTACKER_CONTROLLED"] = 10**18
            self.resign(path, tampered)
            self.assertIsNone(self.validate(path, contract, reusable, output_root, extraction))

    def test_terminal_resume_recomputes_candidate_semantic_sha(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, receipt, reusable, contract, output_root, extraction = self.make_terminal(
                pathlib.Path(tmp)
            )
            tampered = json.loads(json.dumps(receipt))
            tampered["candidate_table"]["semantic_sha256"] = "0" * 64
            self.resign(path, tampered)
            self.assertIsNone(self.validate(path, contract, reusable, output_root, extraction))

    def test_terminal_resume_rejects_resigned_candidate_content_drift(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, receipt, reusable, contract, output_root, extraction = self.make_terminal(
                pathlib.Path(tmp)
            )
            candidate_path = pathlib.Path(receipt["candidate_table"]["path"])
            with gzip.open(candidate_path, "rt", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            rows[0]["profile_log_likelihood"] = "-999.0"
            with gzip.open(candidate_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, CANDIDATE_TABLE_COLUMNS, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            semantic = hashlib.sha256()
            for row in rows:
                semantic.update(
                    json.dumps(
                        row, sort_keys=True, separators=(",", ":"), ensure_ascii=False
                    ).encode()
                    + b"\n"
                )
            tampered = json.loads(json.dumps(receipt))
            tampered["candidate_table"].update({
                "size_bytes": candidate_path.stat().st_size,
                "sha256": sha256_path(candidate_path),
                "semantic_sha256": semantic.hexdigest(),
            })
            self.resign(path, tampered)
            self.assertIsNone(
                self.validate(path, contract, reusable, output_root, extraction)
            )

    def test_main_does_not_emit_reused_final_for_resigned_candidate_drift(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, receipt, reusable, _, output_root, extraction = self.make_terminal(
                pathlib.Path(tmp)
            )
            candidate_path = pathlib.Path(receipt["candidate_table"]["path"])
            with gzip.open(candidate_path, "rt", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            rows[0]["profile_log_likelihood"] = "-999.0"
            with gzip.open(candidate_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, CANDIDATE_TABLE_COLUMNS, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            semantic = hashlib.sha256()
            for row in rows:
                semantic.update(
                    json.dumps(
                        row, sort_keys=True, separators=(",", ":"), ensure_ascii=False
                    ).encode()
                    + b"\n"
                )

            args = argparse.Namespace(
                extraction_root=extraction.parent,
                release_contract_manifest=None,
                output_root=output_root,
                workers=1,
                thresholds="1,2,3,5",
                component_bases="PS_HP1,PS_HP2",
                hp_families="1,2",
                structural_exact_pattern_minread_grid="1,2,3,5",
                primary_structural_exact_pattern_minread=3,
                exact_k_max=12,
                max_vertex_sets=256,
                solver_time_limit_seconds=30.0,
                fixed_error_grid="0.005,0.01,0.02,0.05",
                minimum_bq_error_rate=1e-6,
                maximum_bq_error_rate=0.25,
                conditional_candidate_ranking_bootstrap_replicates=20,
                conditional_candidate_ranking_bootstrap_seed=20260716,
                tie_tolerance=1e-6,
                max_new_tasks=16,
                task_timeout_seconds=28800.0,
                timeout_grace_seconds=300.0,
                aggregate_only=False,
                preflight_only=False,
            )
            processed = argparse.Namespace(**vars(args))
            processed.thresholds = (1, 2, 3, 5)
            processed.component_bases = ("PS_HP1", "PS_HP2")
            processed.hp_families = ("1", "2")
            processed.structural_exact_pattern_minreads = (1, 2, 3, 5)
            processed.fixed_error_grid_values = (0.005, 0.01, 0.02, 0.05)
            run_contract = {
                "full_extraction_receipt": {
                    "path": str(extraction.resolve()),
                    "sha256": sha256_path(extraction),
                },
                "ranker": {
                    "path": str(ranking_runner.RANKER.resolve()),
                    "sha256": sha256_path(ranking_runner.RANKER),
                },
                "hypercube_solver": {
                    "path": str(ranking_runner.HYPERCUBE_SOLVER.resolve()),
                    "sha256": sha256_path(ranking_runner.HYPERCUBE_SOLVER),
                },
                "runner": {
                    "path": str(ranking_runner.RUNNER),
                    "sha256": sha256_path(ranking_runner.RUNNER),
                },
                "method_contract": EXPECTED_METHOD_CONTRACT,
                "thresholds": [1, 2, 3, 5],
                "component_bases": ["PS_HP1", "PS_HP2"],
                "hp_families": ["1", "2"],
                "parameters": ranking_runner._rank_parameters(processed),
                "workers": 1,
                "task_timeout_seconds": 28800.0,
                "timeout_grace_seconds": 300.0,
            }
            tampered = json.loads(json.dumps(receipt))
            tampered["run_contract"] = run_contract
            tampered["candidate_table"].update({
                "size_bytes": candidate_path.stat().st_size,
                "sha256": sha256_path(candidate_path),
                "semantic_sha256": semantic.hexdigest(),
            })
            self.resign(path, tampered)
            full_extraction = {
                "all_pass": True,
                "results": [
                    {
                        "dataset": dataset, "chrom": chrom, "status": "PASS",
                        "receipt": {"scope": {"dataset": dataset, "chrom": chrom}},
                    }
                    for dataset in ranking_runner.DATASETS
                    for chrom in ranking_runner.AUTOSOMES
                ],
            }
            aggregate = self.synthetic_aggregate()
            conservation = {"all_aggregate_cells_conserved": True, "failed_cells": []}
            runtime = self.synthetic_runtime()
            with (
                mock.patch.object(ranking_runner, "parse_args", return_value=args),
                mock.patch.object(
                    ranking_runner, "load_full_extraction_receipt",
                    return_value=(full_extraction, extraction),
                ),
                mock.patch.object(
                    ranking_runner, "active_conflicts",
                    return_value={"process_count": 0, "root_count": 0, "representatives": []},
                ),
                mock.patch.object(
                    ranking_runner, "resource_gate_preview", return_value={"disk_pass": True}
                ),
                mock.patch.object(ranking_runner, "ensure_preflight_contract"),
                mock.patch.object(
                    ranking_runner, "scan_existing_rank_specs",
                    return_value=(reusable, [], []),
                ),
                mock.patch.object(
                    ranking_runner, "aggregate_rank_receipts", return_value=aggregate
                ),
                mock.patch.object(
                    ranking_runner, "aggregate_conservation_audit",
                    return_value=conservation,
                ),
                mock.patch.object(
                    ranking_runner, "aggregate_primary_runtime_diagnostics",
                    return_value=runtime,
                ),
                mock.patch.object(
                    ranking_runner, "child_runtime_contract_valid", return_value=True
                ),
                mock.patch.object(
                    ranking_runner, "child_method_contract_and_ranker_source_bound",
                    return_value=True,
                ),
                mock.patch("builtins.print") as printed,
            ):
                with self.assertRaisesRegex(RuntimeError, "terminal full ranking receipt is invalid"):
                    ranking_runner.main()
            self.assertFalse(any(
                "REUSED_FINAL" in str(call.args) for call in printed.call_args_list
            ))


if __name__ == "__main__":
    unittest.main()
