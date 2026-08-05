#!/usr/bin/env python3
"""Synthetic, non-launching tests for the layered-v3 handoff supervisor."""

from __future__ import annotations

import fcntl
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPT = Path(__file__).resolve().parent / "continue_layered_v3_after_producer.py"
SPEC = importlib.util.spec_from_file_location("layered_v3_handoff", SCRIPT)
assert SPEC and SPEC.loader
handoff = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = handoff
SPEC.loader.exec_module(handoff)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def pin(path: Path) -> dict[str, str]:
    return {
        "path": str(path),
        "resolved_path": str(path.resolve()),
        "sha256": digest(path.resolve()),
    }


class FakeProcess:
    def __init__(self, pid: int, exit_code: int):
        self.pid = pid
        self.exit_code = exit_code

    def wait(self) -> int:
        return self.exit_code


class MockExecutor:
    def __init__(self, case: "HandoffSupervisorTest", *, fail_stage: str | None = None, runner_exit: int = 0):
        self.case = case
        self.fail_stage = fail_stage
        self.runner_exit = runner_exit
        self.commands: list[list[str]] = []
        self.popen_commands: list[list[str]] = []

    @staticmethod
    def option(argv: list[str], name: str) -> str:
        return argv[argv.index(name) + 1]

    def run(self, argv: list[str], log_path: Path, env: dict[str, str]) -> int:
        argv = list(argv)
        self.commands.append(argv)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("mock command\n", encoding="utf-8")
        stage = Path(argv[1]).name
        sample = self.option(argv, "--sample") if "--sample" in argv else None
        label = f"receipt:{sample}" if sample else stage
        if self.fail_stage in {label, stage}:
            return 5
        if stage == "finalize_longphase_production_sidecars.py":
            self.case.publish_closeout()
        elif stage == "build_longphase_production_capture_receipt_v3.py":
            write_json(Path(self.option(argv, "--output")), {
                "schema_name": "intersubmod.longphase_production_capture_receipt",
                "schema_version": "1.0.0",
                "sample": sample,
            })
        elif stage == "prepare_clean_layered_manifest_v3.py":
            write_json(Path(self.option(argv, "--output")), {
                "schema_name": "intersubmod.layered_input_manifest",
                "schema_version": "3.0.0",
                "manifest_id": self.option(argv, "--manifest-id"),
                "dataset_count": 7,
                "biological_sample_count": 6,
                "samples": [{"sample": sample} for sample in handoff.EXPECTED_SAMPLES],
            })
        return 0

    def popen(self, argv: list[str], log_path: Path, env: dict[str, str]):
        argv = list(argv)
        self.popen_commands.append(argv)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log = log_path.open("xb")
        log.write(b"mock foreground runner\n")
        if self.runner_exit == 0:
            run_root = self.case.run_parent / self.case.run_id
            run_root.mkdir()
            verification = run_root / "verification_summary.json"
            write_json(verification, {"dataset_count": 7, "all_pass": True})
            write_json(run_root / "_SUCCESS", {
                "schema_name": "intersubmod.layered_success_marker",
                "schema_version": "1.0.0",
                "run_id": self.case.run_id,
                "verification_path": str(verification),
                "verification_sha256": digest(verification),
                "extra": {"mode": "full", "dataset_count": 7, "scope": "chr1-22"},
            })
        return FakeProcess(424242, self.runner_exit), log


class HandoffSupervisorTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix="layered-v3-handoff-")
        self.tmp = Path(self.temp.name)
        self.run_parent = self.tmp / "research_rounds"
        self.run_parent.mkdir()
        self.production = self.run_parent / "producer"
        (self.production / "samples").mkdir(parents=True)
        self.reviewed = self.tmp / "longphase_production_sidecar_manifest.json"
        self.snapshot = self.production / "input_manifest.json"
        self.base = self.tmp / "layered_v2_input_manifest.json"
        self.real = self.tmp / "equivalence_probe.json"
        self.synthetic = self.tmp / "synthetic_contract_receipt.json"
        self.longphase = self.tmp / "longphase-s"
        self.longphase.write_bytes(b"synthetic longphase\n")
        self.longphase.chmod(0o755)
        self.longphase_sha = digest(self.longphase)
        self.run_id = "20260711_layered_reconstruction_v3_full_no_truth_v1"
        self.plan_id = "20260711_handoff_test"
        self.workspace = self.run_parent / ".layered_v3_handoff" / self.plan_id
        samples = [
            {"sample": sample, "biological_id": "HCC1395" if sample == "HCC1395_DORADO" else sample}
            for sample in handoff.EXPECTED_SAMPLES
        ]
        producer_manifest = {
            "dataset_count": 7,
            "biological_sample_count": 6,
            "samples": samples,
        }
        write_json(self.reviewed, producer_manifest)
        self.snapshot.write_bytes(self.reviewed.read_bytes())
        write_json(self.base, {"dataset_count": 7, "biological_sample_count": 6, "samples": samples})
        write_json(self.real, {"pass": True, "kind": "real"})
        write_json(self.synthetic, {"pass": True, "kind": "synthetic"})
        write_json(self.production / "params.json", {
            "truth_flags": False,
            "output_mode": "read_tag_sidecar",
        })
        with (self.production / "run_status.tsv").open("w", encoding="utf-8") as handle:
            handle.write("timestamp\tsample\tstage\tstatus\tdetail\n")
            for sample in handoff.EXPECTED_SAMPLES:
                handle.write(f"t0\t{sample}\tproduction_tagging\tSTART\t\n")
                handle.write(f"t1\t{sample}\tproduction_tagging\tPASS\trows=1\n")
                (self.production / "samples" / sample).mkdir()
            handle.write("t2\tALL\tverify\tPASS\t7/7\n")
        write_json(self.production / "verification_summary.json", {
            "dataset_count": 7,
            "n_pass": 7,
            "all_pass": True,
            "samples": [
                {
                    "pass": True,
                    "sidecar": str(self.production / "samples" / sample / f"{sample}.read_tags.tsv"),
                }
                for sample in handoff.EXPECTED_SAMPLES
            ],
        })
        producer_source = self.tmp / "producer_source.py"
        producer_source.write_text("# producer\n", encoding="utf-8")
        (self.production / "code.sha256").write_text(
            f"{digest(self.reviewed)}  {self.reviewed}\n{digest(producer_source)}  {producer_source}\n",
            encoding="utf-8",
        )
        self.tool_dir = self.tmp / "tools"
        self.tool_dir.mkdir()
        self.tools = {}
        for name in ("stat", "samtools", "bcftools", "bgzip", "tabix", "ps"):
            path = self.tool_dir / name
            path.write_bytes(f"{name}-v1\n".encode())
            path.chmod(0o755)
            self.tools[name] = path

        self.patch = mock.patch.multiple(
            handoff,
            EXPECTED_PRODUCTION_ROOT=self.production,
            EXPECTED_RUN_PARENT=self.run_parent,
            EXPECTED_REVIEWED_PRODUCER_MANIFEST=self.reviewed,
            EXPECTED_BASE_MANIFEST=self.base,
            EXPECTED_REAL_RECEIPT=self.real,
            EXPECTED_SYNTHETIC_RECEIPT=self.synthetic,
            EXPECTED_LONGPHASE_BINARY=self.longphase,
            EXPECTED_LONGPHASE_SHA256=self.longphase_sha,
            EXPECTED_PS_BINARY=self.tools["ps"],
        )
        self.patch.start()
        self.plan_path = self.tmp / "launch_plan.json"
        self.write_plan(execute=True)
        self.plan = handoff.load_launch_plan(self.plan_path)

    def tearDown(self) -> None:
        self.patch.stop()
        self.temp.cleanup()

    def executable_paths(self) -> dict[str, Path]:
        actual = {
            "python": Path(sys.executable),
            "handoff_supervisor": SCRIPT,
            "longphase_binary": self.longphase,
            "producer_finalizer": handoff.METHOD_SCRIPTS / "finalize_longphase_production_sidecars.py",
            "receipt_normalizer": handoff.METHOD_SCRIPTS / "build_longphase_production_capture_receipt_v3.py",
            "v3_preparer": handoff.METHOD_SCRIPTS / "prepare_clean_layered_manifest_v3.py",
            "v3_runner": handoff.REPO_ROOT / "scripts/run_layered_v3.py",
            "v3_lifecycle": handoff.REPO_ROOT / "scripts/layered_v3_lifecycle.py",
            "v3_validator": handoff.METHOD_SCRIPTS / "validate_layered_v3_inputs.py",
            "v3_verifier": handoff.REPO_ROOT / "scripts/verify_layered_v3.py",
            "schema_manifest": handoff.METHOD_SCHEMAS / "layered_input_manifest_v3.schema.json",
            "schema_lock": handoff.METHOD_SCHEMAS / "layered_input_lock_v1.schema.json",
            "schema_receipt": handoff.METHOD_SCHEMAS / "longphase_production_capture_receipt_v1.schema.json",
            "science_sm_linkage_genomewide": handoff.METHOD_SCRIPTS / "sm_linkage_genomewide.py",
            "science_sm_multilocus_combinations": handoff.METHOD_SCRIPTS / "sm_multilocus_combinations.py",
            "science_tree_enumeration_solver": handoff.METHOD_SCRIPTS / "tree_enumeration_solver.py",
            "science_layered_tree_reconstruction": handoff.METHOD_SCRIPTS / "layered_tree_reconstruction.py",
            "science_build_region_view": handoff.METHOD_SCRIPTS / "build_region_view.py",
            "science_build_ssnv_site_ledger": handoff.METHOD_SCRIPTS / "build_ssnv_site_ledger.py",
            "tool_stat": self.tools["stat"],
            "tool_samtools": self.tools["samtools"],
            "tool_bcftools": self.tools["bcftools"],
            "tool_bgzip": self.tools["bgzip"],
            "tool_tabix": self.tools["tabix"],
            "tool_ps": self.tools["ps"],
        }
        self.assertEqual(set(actual), set(handoff.REQUIRED_EXECUTABLE_ROLES))
        return actual

    def write_plan(self, *, execute: bool) -> None:
        inputs = {
            "producer_manifest_reviewed": self.reviewed,
            "producer_manifest_snapshot": self.snapshot,
            "base_manifest": self.base,
            "real_method_receipt": self.real,
            "synthetic_method_receipt": self.synthetic,
        }
        value = {
            "schema_name": "intersubmod.layered_v3_handoff_launch_plan",
            "schema_version": "1.0.0",
            "plan_id": self.plan_id,
            "execution_authorization": {
                "execute": execute,
                "reviewed_by": "synthetic-reviewer" if execute else None,
                "reviewed_at_utc": "2026-07-11T00:00:00Z" if execute else None,
            },
            "production_root": str(self.production),
            "handoff": {
                "lock": str(self.production.parent / f".{self.production.name}.handoff.lock"),
                "workspace": str(self.workspace),
                "source_manifest": str(self.workspace / "layered_input_manifest_v3.json"),
                "source_manifest_failure": str(self.workspace / "layered_input_manifest_v3.failure.json"),
                "manifest_id": "20260711_layered_v3_full_test_manifest",
            },
            "layered_run": {"run_parent": str(self.run_parent), "run_id": self.run_id},
            "inputs": {role: pin(path) for role, path in inputs.items()},
            "executables": {role: pin(path) for role, path in self.executable_paths().items()},
        }
        write_json(self.plan_path, value)

    def publish_closeout(self, *, tracked: Path | None = None) -> None:
        closeout = self.production / "closeout"
        closeout.mkdir()
        receipt = closeout / "production_closeout.json"
        write_json(receipt, {
            "status": "PASS",
            "dataset_count": 7,
            "n_pass": 7,
            "truth_flags": False,
            "binary_sha256": self.longphase_sha,
            "samples": [{"sample": sample} for sample in handoff.EXPECTED_SAMPLES],
        })
        tracked = tracked or (self.production / "run_status.tsv")
        hashes = closeout / "artifacts.final.sha256"
        hashes.write_text(f"{digest(tracked)}  {tracked}\n", encoding="utf-8")
        write_json(self.production / "_SUCCESS", {
            "status": "SUCCESS",
            "closeout_receipt": str(receipt),
            "closeout_receipt_sha256": digest(receipt),
            "artifacts_manifest": str(hashes),
            "artifacts_manifest_sha256": digest(hashes),
        })

    def supervisor(self, *, executor=None, processes=(), process_provider=None, sleeper=lambda _: None):
        return handoff.HandoffSupervisor(
            self.plan,
            process_provider=process_provider or (lambda: list(processes)),
            executor=executor,
            sleeper=sleeper,
            poll_seconds=0,
        )

    def assert_code(self, expected: str, callable_) -> None:
        with self.assertRaises(handoff.HandoffError) as caught:
            callable_()
        self.assertEqual(caught.exception.code, expected)

    def test_default_is_dry_run_and_starts_zero_subprocesses(self) -> None:
        args = handoff.build_parser().parse_args(["--launch-plan", str(self.plan_path)])
        self.assertFalse(args.execute_reviewed_plan)
        result = self.supervisor(executor=MockExecutor(self)).dry_run()
        self.assertEqual(result["subprocesses_started"], 0)
        self.assertEqual(len(result["commands"]), 10)

    def test_dry_run_commands_bind_reviewed_vs_snapshot_authorities(self) -> None:
        result = self.supervisor(executor=MockExecutor(self)).dry_run()
        commands = {item["stage"]: item["argv"] for item in result["commands"]}
        self.assertIn(str(self.reviewed), commands["producer_closeout"])
        self.assertIn(str(self.snapshot), commands["receipt:HCC1395"])
        self.assertIn(str(self.reviewed), commands["prepare_v3_manifest"])
        self.assertEqual(commands["run_layered_v3"][-len(handoff.RUNNER_OPTIONS):], list(handoff.RUNNER_OPTIONS))
        self.assertFalse(any("--truth-vcf" in token or "--truth-bed" in token for argv in commands.values() for token in argv))
        self.assertTrue(any("no_truth" in token for token in commands["run_layered_v3"]))

    def test_token_level_truth_gate_does_not_reject_no_truth_run_id(self) -> None:
        commands = handoff.build_commands(self.plan, handoff.EXPECTED_SAMPLES)
        self.assertTrue(any(self.run_id in argv for _, argv in commands))
        with self.assertRaises(handoff.HandoffError) as caught:
            handoff.validate_command_policy([("bad", ["tool", "--truth-bed=/tmp/forbidden.bed"])])
        self.assertEqual(caught.exception.code, "E_COMMAND_POLICY")

    def test_active_producer_is_rejected(self) -> None:
        rows = [{"pid": 10, "ppid": 1, "argv": f"longphase-s -o {self.production}/samples/HCC1395/x"}]
        self.assert_code("E_PRODUCER_PENDING", lambda: self.supervisor(processes=rows).dry_run())

    def test_external_v2_launcher_is_rejected(self) -> None:
        rows = [{"pid": 11, "ppid": 1, "argv": "bash /repo/run_layered_7samples_newbb.sh"}]
        self.assert_code("E_PROCESS_ACTIVE", lambda: self.supervisor(processes=rows).dry_run())

    def test_start_only_status_is_rejected(self) -> None:
        lines = (self.production / "run_status.tsv").read_text().splitlines()
        (self.production / "run_status.tsv").write_text("\n".join(line for line in lines if "\tPASS\t" not in line) + "\n")
        self.assert_code("E_PRODUCER_ABANDONED", lambda: self.supervisor().dry_run())

    def test_pass_before_start_is_fatal(self) -> None:
        (self.production / "run_status.tsv").write_text(
            "timestamp\tsample\tstage\tstatus\tdetail\n"
            "t0\tHCC1395\tproduction_tagging\tPASS\trows=1\n",
            encoding="utf-8",
        )
        self.assert_code("E_STATUS_ORDER", lambda: self.supervisor().dry_run())

    def test_wait_progresses_start_to_7_of_7_then_quiescent(self) -> None:
        status = self.production / "run_status.tsv"
        terminal = status.read_bytes()
        initial = (
            "timestamp\tsample\tstage\tstatus\tdetail\n"
            "t0\tHCC1395\tproduction_tagging\tSTART\t\n"
        ).encode()
        status.write_bytes(initial)
        state = {"sleep": 0}

        def provider():
            if state["sleep"] < 2:
                return [{"pid": 10, "ppid": 1, "argv": f"longphase-s -o {self.production}/samples/HCC1395/x"}]
            return []

        def sleeper(_seconds):
            if state["sleep"] == 0:
                status.write_bytes(terminal)
            state["sleep"] += 1

        result = self.supervisor(process_provider=provider, sleeper=sleeper).dry_run(wait_for_producer=True)
        self.assertTrue(result["pass"])
        self.assertEqual(state["sleep"], 2)

    def test_wait_aborts_immediately_on_fail_status(self) -> None:
        status = self.production / "run_status.tsv"
        status.write_text(
            "timestamp\tsample\tstage\tstatus\tdetail\n"
            "t0\tHCC1395\tproduction_tagging\tSTART\t\n",
            encoding="utf-8",
        )
        state = {"sleep": 0}

        def sleeper(_seconds):
            state["sleep"] += 1
            with status.open("a", encoding="utf-8") as handle:
                handle.write("t1\tHCC1395\tproduction_tagging\tFAIL\trc=1\n")

        rows = [{"pid": 10, "ppid": 1, "argv": f"longphase-s -o {self.production}/samples/HCC1395/x"}]
        self.assert_code(
            "E_PRODUCER_FAILED",
            lambda: self.supervisor(processes=rows, sleeper=sleeper).dry_run(wait_for_producer=True),
        )
        self.assertEqual(state["sleep"], 1)

    def test_wait_aborts_on_pinned_hash_drift(self) -> None:
        status = self.production / "run_status.tsv"
        status.write_text(
            "timestamp\tsample\tstage\tstatus\tdetail\n"
            "t0\tHCC1395\tproduction_tagging\tSTART\t\n",
            encoding="utf-8",
        )

        def sleeper(_seconds):
            self.tools["tabix"].write_bytes(b"drift-during-wait\n")

        rows = [{"pid": 10, "ppid": 1, "argv": f"longphase-s -o {self.production}/samples/HCC1395/x"}]
        self.assert_code(
            "E_PIN_HASH_DRIFT",
            lambda: self.supervisor(processes=rows, sleeper=sleeper).dry_run(wait_for_producer=True),
        )

    def test_wait_rejects_status_regression(self) -> None:
        status = self.production / "run_status.tsv"
        status.write_text(
            "timestamp\tsample\tstage\tstatus\tdetail\n"
            "t0\tHCC1395\tproduction_tagging\tSTART\t\n",
            encoding="utf-8",
        )

        def sleeper(_seconds):
            status.write_text("timestamp\tsample\tstage\tstatus\tdetail\n", encoding="utf-8")

        rows = [{"pid": 10, "ppid": 1, "argv": f"longphase-s -o {self.production}/samples/HCC1395/x"}]
        self.assert_code(
            "E_STATUS_REGRESSION",
            lambda: self.supervisor(processes=rows, sleeper=sleeper).dry_run(wait_for_producer=True),
        )

    def test_non_7_of_7_verification_is_rejected(self) -> None:
        value = json.loads((self.production / "verification_summary.json").read_text())
        value["n_pass"] = 6
        write_json(self.production / "verification_summary.json", value)
        self.assert_code("E_VERIFICATION_INCOMPLETE", lambda: self.supervisor().dry_run())

    def test_producer_code_hash_drift_is_rejected(self) -> None:
        self.reviewed.write_bytes(self.reviewed.read_bytes() + b"drift")
        self.assert_code("E_PIN_HASH_DRIFT", lambda: self.supervisor().dry_run())

    def test_post_closeout_artifact_drift_is_rejected(self) -> None:
        tracked = self.tmp / "tracked.dat"
        tracked.write_text("before\n", encoding="utf-8")
        self.publish_closeout(tracked=tracked)
        tracked.write_text("after\n", encoding="utf-8")
        self.assert_code("E_HASH_DRIFT", lambda: self.supervisor().dry_run())

    def test_partial_closeout_is_rejected(self) -> None:
        (self.production / "closeout").mkdir()
        self.assert_code("E_CLOSEOUT_PARTIAL", lambda: self.supervisor().dry_run())

    def test_partial_receipts_are_rejected(self) -> None:
        write_json(
            self.production / "samples/HCC1395/producer_capture_receipt_v3.json",
            {"sample": "HCC1395"},
        )
        self.assert_code("E_PARTIAL_RECEIPTS", lambda: self.supervisor().dry_run())

    def test_existing_layered_run_id_is_rejected(self) -> None:
        (self.run_parent / self.run_id).mkdir()
        self.assert_code("E_RUN_ROOT_EXISTS", lambda: self.supervisor().dry_run())

    def test_source_hash_drift_is_rejected(self) -> None:
        self.tools["tabix"].write_bytes(b"tabix-v2\n")
        self.assert_code("E_PIN_HASH_DRIFT", lambda: self.supervisor().dry_run())

    def test_controlled_path_resolves_every_pinned_tool(self) -> None:
        env, resolutions = handoff.controlled_environment(self.plan)
        self.assertEqual(Path(env["PATH"]), self.tool_dir)
        self.assertEqual(set(resolutions), {"stat", "samtools", "bcftools", "bgzip", "tabix"})
        self.assertTrue(all(Path(path).parent == self.tool_dir for path in resolutions.values()))

    def test_atomic_json_never_replaces_existing_target(self) -> None:
        target = self.tmp / "immutable.json"
        handoff.atomic_write_json(target, {"version": 1})
        before = target.read_bytes()
        self.assert_code("E_OUTPUT_EXISTS", lambda: handoff.atomic_write_json(target, {"version": 2}))
        self.assertEqual(target.read_bytes(), before)

    def test_second_supervisor_is_rejected_by_flock(self) -> None:
        lock_path = self.production.parent / f".{self.production.name}.handoff.lock"
        held = lock_path.open("a+b")
        fcntl.flock(held.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            self.assert_code("E_HANDOFF_LOCKED", lambda: self.supervisor().dry_run())
        finally:
            held.close()

    def test_execution_requires_out_of_band_plan_hash(self) -> None:
        rc = handoff.main(["--launch-plan", str(self.plan_path), "--execute-reviewed-plan"])
        self.assertEqual(rc, 7)

    def test_execution_rejects_wrong_out_of_band_plan_hash(self) -> None:
        rc = handoff.main([
            "--launch-plan", str(self.plan_path),
            "--execute-reviewed-plan",
            "--expected-launch-plan-sha256", "0" * 64,
        ])
        self.assertEqual(rc, 7)

    def test_execution_rejects_unreviewed_json_authorization(self) -> None:
        self.write_plan(execute=False)
        plan_sha = digest(self.plan_path)
        rc = handoff.main([
            "--launch-plan", str(self.plan_path),
            "--execute-reviewed-plan",
            "--expected-launch-plan-sha256", plan_sha,
        ])
        self.assertEqual(rc, 7)

    def test_happy_execute_uses_mocked_foreground_runner(self) -> None:
        executor = MockExecutor(self)
        result = self.supervisor(executor=executor).execute()
        self.assertTrue(result["pass"])
        self.assertEqual(len(result["receipt_sha256"]), 7)
        self.assertEqual(len(executor.commands), 9)  # closeout + seven receipts + preparer
        self.assertEqual(len(executor.popen_commands), 1)
        self.assertTrue((self.workspace / "runner.started.json").is_file())
        self.assertTrue((self.workspace / "runner.exit.json").is_file())
        self.assertTrue((self.workspace / "handoff.success.json").is_file())

    def test_normalizer_failure_stops_without_runner_or_success(self) -> None:
        executor = MockExecutor(self, fail_stage="receipt:COLO829")
        self.assert_code("E_STAGE_EXIT", lambda: self.supervisor(executor=executor).execute())
        self.assertEqual(executor.popen_commands, [])
        self.assertFalse((self.workspace / "handoff.success.json").exists())

    def test_preparer_failure_stops_without_runner(self) -> None:
        executor = MockExecutor(self, fail_stage="prepare_clean_layered_manifest_v3.py")
        self.assert_code("E_STAGE_EXIT", lambda: self.supervisor(executor=executor).execute())
        self.assertEqual(executor.popen_commands, [])

    def test_runner_nonzero_is_preserved_and_no_success_is_fabricated(self) -> None:
        executor = MockExecutor(self, runner_exit=23)
        with self.assertRaises(handoff.HandoffError) as caught:
            self.supervisor(executor=executor).execute()
        self.assertEqual(caught.exception.code, "E_LAYERED_RUNNER_EXIT")
        self.assertEqual(caught.exception.exit_code, 23)
        self.assertFalse((self.workspace / "handoff.success.json").exists())
        exit_receipt = json.loads((self.workspace / "runner.exit.json").read_text())
        self.assertEqual(exit_receipt["exit_code"], 23)


if __name__ == "__main__":
    unittest.main(verbosity=2)
