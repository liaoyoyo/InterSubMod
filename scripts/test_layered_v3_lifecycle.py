#!/usr/bin/env python3
"""Tiny, production-independent tests for layered_v3_lifecycle.py."""

from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

from layered_v3_lifecycle import (
    ChildProcessFailure,
    EnvironmentContractError,
    HeartbeatError,
    LifecycleSignal,
    RunLifecycle,
    RunLockError,
    SourceBundleError,
    StateTransitionError,
    atomic_write_json,
    build_source_bundle,
    capture_environment_lock,
    sha256_file,
    validate_heartbeat,
    verify_source_bundle,
)


class LayeredV3LifecycleTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="layered-v3-lifecycle-test-")
        self.root = Path(self.temporary.name)
        self.run_parent = self.root / "runs"
        self.run_parent.mkdir()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def make_sources(self) -> dict:
        source_dir = self.root / "source"
        source_dir.mkdir(exist_ok=True)
        sources = {}
        for name in ("runner", "validator", "verifier", "consumer", "solver"):
            path = source_dir / f"{name}.py"
            path.write_text(f"# {name}\nVALUE = {name!r}\n", encoding="utf-8")
            path.chmod(0o750 if name == "runner" else 0o640)
            sources[name] = path
        return sources

    @staticmethod
    def locked_environment(extra: dict = None) -> dict:
        environment = {
            "LC_ALL": "C",
            "TZ": "UTC",
            "PYTHONHASHSEED": "0",
            "SM_WORKERS": "2",
        }
        environment.update(extra or {})
        return environment

    def freeze_preflight(self, lifecycle: RunLifecycle, sources: dict) -> None:
        lifecycle.begin_preflight()
        lifecycle.write_frozen_lock(
            {
                "schema_name": "intersubmod.layered_input_lock",
                "schema_version": "1.0.0",
                "dataset_count": 7,
                "sample_ids": ["S1", "S2", "S3", "S4", "S5", "S6", "S7"],
            }
        )
        lifecycle.build_source_bundle(
            sources["runner"],
            sources["validator"],
            sources["verifier"],
            [sources["consumer"], sources["solver"]],
        )
        lifecycle.capture_environment_lock(
            sm_allowlist={"SM_WORKERS"},
            required_sm={"SM_WORKERS"},
            environment=self.locked_environment(),
            distributions=(),
            tools=(),
            storage_path=self.run_parent,
        )

    def test_preflight_failure_never_creates_published_root(self) -> None:
        staging = None
        with RunLifecycle(self.run_parent, "preflight-fail", install_traps=False) as lifecycle:
            staging = lifecycle.staging_root
            lifecycle.begin_preflight()
            lifecycle.fail("E_FIXTURE_PREFLIGHT", "intentional fixture failure", grace_seconds=0.01)
            self.assertEqual(lifecycle.state, "FAILED")
            self.assertFalse(lifecycle.published)
            self.assertFalse((self.run_parent / "preflight-fail").exists())
            self.assertFalse((lifecycle.root / "_SUCCESS").exists())
        self.assertIsNotNone(staging)
        self.assertTrue(staging.exists())
        state = json.loads((staging / "run_state.json").read_text(encoding="utf-8"))
        self.assertEqual(state["state"], "FAILED")

    def test_parent_flock_rejects_second_launcher_for_same_run_id(self) -> None:
        first = RunLifecycle(self.run_parent, "same-id", install_traps=False)
        try:
            with self.assertRaises(RunLockError) as captured:
                RunLifecycle(self.run_parent, "same-id", install_traps=False)
            self.assertEqual(captured.exception.code, "E_RUN_LOCKED")
            candidates = list(self.run_parent.glob(".same-id.staging.*"))
            self.assertEqual(len(candidates), 1)
        finally:
            first.fail("E_TEST_END", "test cleanup", grace_seconds=0.01)
            first.close()

    def test_illegal_state_transition_is_rejected(self) -> None:
        with self.assertRaises(StateTransitionError) as captured:
            with RunLifecycle(self.run_parent, "bad-edge", install_traps=False) as lifecycle:
                lifecycle.transition("RUNNING", "running", "skip preflight")
        self.assertEqual(captured.exception.code, "E_STATE_TRANSITION")
        staging = next(self.run_parent.glob(".bad-edge.staging.*"))
        state = json.loads((staging / "run_state.json").read_text(encoding="utf-8"))
        self.assertEqual(state["state"], "FAILED")
        self.assertFalse((self.run_parent / "bad-edge").exists())

    def test_incomplete_preflight_cannot_publish(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "E_PREFLIGHT_INCOMPLETE"):
            with RunLifecycle(self.run_parent, "incomplete-preflight", install_traps=False) as lifecycle:
                lifecycle.begin_preflight()
                lifecycle.write_frozen_lock({"schema_name": "fixture", "schema_version": "1"})
                lifecycle.publish_ready()
        self.assertFalse((self.run_parent / "incomplete-preflight").exists())
        staging = next(self.run_parent.glob(".incomplete-preflight.staging.*"))
        state = json.loads((staging / "run_state.json").read_text(encoding="utf-8"))
        self.assertEqual(state["state"], "FAILED")
        self.assertFalse((staging / "_SUCCESS").exists())

    def test_source_bundle_copies_bytes_and_detects_tamper(self) -> None:
        sources = self.make_sources()
        original_consumer = sources["consumer"].read_bytes()
        bundle = build_source_bundle(
            self.root / "bundle",
            sources["runner"],
            sources["validator"],
            sources["verifier"],
            [sources["consumer"], sources["solver"]],
        )
        self.assertEqual(bundle.role_paths["imported:000:consumer.py"].read_bytes(), original_consumer)

        # Working-tree edits do not change the frozen execution authority.
        sources["consumer"].write_text("# changed after freeze\n", encoding="utf-8")
        verified = verify_source_bundle(
            bundle.manifest_path,
            expected_manifest_sha256=bundle.manifest_sha256,
            expected_content_sha256=bundle.content_sha256,
        )
        self.assertEqual(verified.role_paths["imported:000:consumer.py"].read_bytes(), original_consumer)

        # A one-byte edit inside the bundle invalidates the frozen authority.
        bundled_runner = bundle.role_paths["runner"]
        tampered = bytearray(bundled_runner.read_bytes())
        tampered[0] ^= 1
        bundled_runner.write_bytes(bytes(tampered))
        with self.assertRaises(SourceBundleError) as captured:
            verify_source_bundle(
                bundle.manifest_path,
                expected_manifest_sha256=bundle.manifest_sha256,
                expected_content_sha256=bundle.content_sha256,
            )
        self.assertEqual(captured.exception.code, "E_SOURCE_BUNDLE_MISMATCH")

    def test_source_bundle_requires_imported_scripts(self) -> None:
        sources = self.make_sources()
        with self.assertRaises(SourceBundleError) as captured:
            build_source_bundle(
                self.root / "missing-imports",
                sources["runner"],
                sources["validator"],
                sources["verifier"],
                [],
            )
        self.assertEqual(captured.exception.code, "E_SOURCE_BUNDLE_IMPORTS_EMPTY")

    def test_unknown_sm_environment_is_rejected_before_publish(self) -> None:
        with RunLifecycle(self.run_parent, "bad-env", install_traps=False) as lifecycle:
            lifecycle.begin_preflight()
            with self.assertRaises(EnvironmentContractError) as captured:
                lifecycle.capture_environment_lock(
                    sm_allowlist={"SM_WORKERS"},
                    environment=self.locked_environment({"SM_UNREVIEWED": "1"}),
                    distributions=(),
                    tools=(),
                    storage_path=self.run_parent,
                )
            self.assertEqual(captured.exception.code, "E_UNKNOWN_SM_ENV")
            lifecycle.fail(captured.exception.code, captured.exception.message, grace_seconds=0.01)
        self.assertFalse((self.run_parent / "bad-env").exists())
        staging = next(self.run_parent.glob(".bad-env.staging.*"))
        self.assertFalse((staging / "environment_lock.json").exists())

    def test_environment_lock_records_allowlisted_value_and_origin(self) -> None:
        lock_path = self.root / "environment.json"
        digest = capture_environment_lock(
            lock_path,
            sm_allowlist={"SM_WORKERS", "SM_OPTIONAL"},
            sm_defaults={"SM_OPTIONAL": "default-value"},
            required_sm={"SM_WORKERS"},
            environment=self.locked_environment(),
            distributions=(),
            tools=(),
            storage_path=self.root,
        )
        self.assertEqual(digest, sha256_file(lock_path, reject_symlink=True))
        document = json.loads(lock_path.read_text(encoding="utf-8"))
        self.assertEqual(document["sm_variables"]["SM_WORKERS"], {"origin": "environment", "value": "2"})
        self.assertEqual(
            document["sm_variables"]["SM_OPTIONAL"],
            {"origin": "default", "value": "default-value"},
        )
        self.assertEqual(document["determinism"], {"LC_ALL": "C", "PYTHONHASHSEED": "0", "TZ": "UTC"})

    def test_publish_heartbeat_verify_and_success_marker(self) -> None:
        sources = self.make_sources()
        with RunLifecycle(self.run_parent, "happy-path", install_traps=False) as lifecycle:
            self.freeze_preflight(lifecycle, sources)
            self.assertFalse((self.run_parent / "happy-path").exists())
            published = lifecycle.publish_ready({"fixture": "tiny"})
            self.assertEqual(published, self.run_parent / "happy-path")
            self.assertFalse(lifecycle.staging_root.exists())
            self.assertFalse((published / "_SUCCESS").exists())
            lifecycle.begin_running()
            lifecycle.set_active_samples(["S1"])
            heartbeat = lifecycle.start_heartbeat(interval_seconds=0.02)
            time.sleep(0.06)
            payload = validate_heartbeat(
                published / "heartbeat.json",
                max_age_seconds=1.0,
                expected_launcher_pid_start_time=lifecycle.launcher_pid_start_time,
            )
            self.assertGreaterEqual(payload["seq"], 2)
            self.assertEqual(payload["active_samples"], ["S1"])
            self.assertEqual(payload["frozen_lock_sha256"], lifecycle.frozen_lock_sha256)
            self.assertEqual(payload["source_bundle_sha256"], lifecycle.source_bundle_content_sha256)
            lifecycle.stop_heartbeat()
            self.assertIsNone(heartbeat.error)
            verification = published / "verification.json"
            atomic_write_json(verification, {"all_pass": True, "n_pass": 7, "n_fail": 0})
            lifecycle.begin_verifying()
            marker = lifecycle.succeed(verification, {"dataset_count": 7})
            self.assertEqual(marker, published / "_SUCCESS")
            success = json.loads(marker.read_text(encoding="utf-8"))
            self.assertEqual(success["verification_sha256"], sha256_file(verification, reject_symlink=True))
            self.assertEqual(success["frozen_lock_sha256"], lifecycle.frozen_lock_sha256)
            self.assertEqual(lifecycle.state, "SUCCEEDED")
        state = json.loads((published / "run_state.json").read_text(encoding="utf-8"))
        self.assertEqual(state["state"], "SUCCEEDED")
        self.assertTrue((published / "_SUCCESS").exists())

    def test_stale_heartbeat_and_pid_reuse_are_rejected(self) -> None:
        path = self.root / "heartbeat.json"
        old_time = (dt.datetime.now(dt.timezone.utc) - dt.timedelta(seconds=181)).isoformat().replace("+00:00", "Z")
        payload = {
            "schema_name": "intersubmod.layered_heartbeat",
            "schema_version": "1.0.0",
            "seq": 1,
            "wall_time_utc": old_time,
            "monotonic_seconds": 1.0,
            "host": "fixture",
            "launcher_pid": 1,
            "launcher_pid_start_time": "100",
            "state": "RUNNING",
            "stage": "running",
            "active_samples": [],
            "child_pid_start_times": {},
            "last_event_seq": 4,
            "frozen_lock_sha256": "a" * 64,
            "source_bundle_sha256": "b" * 64,
        }
        atomic_write_json(path, payload)
        with self.assertRaises(HeartbeatError) as stale:
            validate_heartbeat(path, max_age_seconds=180)
        self.assertEqual(stale.exception.code, "E_HEARTBEAT_STALE")

        payload["wall_time_utc"] = dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")
        atomic_write_json(path, payload)
        with self.assertRaises(HeartbeatError) as reused:
            validate_heartbeat(path, expected_launcher_pid_start_time="different")
        self.assertEqual(reused.exception.code, "E_HEARTBEAT_STALE")

    def test_fail_fast_child_registry_cancels_and_reaps_sibling(self) -> None:
        sources = self.make_sources()
        with RunLifecycle(self.run_parent, "child-fail", install_traps=False) as lifecycle:
            self.freeze_preflight(lifecycle, sources)
            lifecycle.publish_ready()
            lifecycle.begin_running()
            failure = lifecycle.children.launch(
                [sys.executable, "-c", "raise SystemExit(7)"],
                "failing-child",
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            sibling = lifecycle.children.launch(
                [sys.executable, "-c", "import time; time.sleep(5)"],
                "sleeping-sibling",
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            with self.assertRaises(ChildProcessFailure) as captured:
                lifecycle.children.wait_all_fail_fast(poll_seconds=0.005, grace_seconds=0.05)
            self.assertEqual(captured.exception.code, "E_CHILD_FAILED")
            self.assertEqual(failure.process.returncode, 7)
            self.assertIsNotNone(sibling.process.returncode)
            self.assertTrue(lifecycle.children.all_reaped())
            lifecycle.fail(captured.exception.code, captured.exception.message, grace_seconds=0.01)
            self.assertFalse((lifecycle.root / "_SUCCESS").exists())
        state = json.loads((self.run_parent / "child-fail" / "run_state.json").read_text(encoding="utf-8"))
        self.assertEqual(state["state"], "FAILED")

    def test_signal_trap_aborts_and_reaps_child_without_success(self) -> None:
        sources = self.make_sources()
        published = self.run_parent / "signal-abort"
        with self.assertRaises(LifecycleSignal):
            with RunLifecycle(self.run_parent, "signal-abort", install_traps=False) as lifecycle:
                self.freeze_preflight(lifecycle, sources)
                lifecycle.publish_ready()
                lifecycle.begin_running()
                lifecycle.children.launch(
                    [sys.executable, "-c", "import time; time.sleep(5)"],
                    "dummy-worker",
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                lifecycle.handle_signal(signal.SIGTERM)
        state = json.loads((published / "run_state.json").read_text(encoding="utf-8"))
        children = json.loads((published / "children.json").read_text(encoding="utf-8"))["children"]
        self.assertEqual(state["state"], "ABORTED")
        self.assertTrue(all(child["returncode"] is not None for child in children))
        self.assertFalse((published / "_SUCCESS").exists())

    def test_installed_sigterm_trap_aborts_real_subprocess(self) -> None:
        sources = self.make_sources()
        ready_path = self.root / "signal-ready.json"
        helper = self.root / "signal_helper.py"
        helper.write_text(
            "\n".join(
                [
                    "import subprocess",
                    "import sys",
                    "import time",
                    "from layered_v3_lifecycle import RunLifecycle, atomic_write_json",
                    f"with RunLifecycle({str(self.run_parent)!r}, 'real-signal', install_traps=True) as lifecycle:",
                    "    lifecycle.begin_preflight()",
                    "    lifecycle.write_frozen_lock({'schema_name': 'fixture', 'schema_version': '1'})",
                    "    lifecycle.build_source_bundle(",
                    f"        {str(sources['runner'])!r},",
                    f"        {str(sources['validator'])!r},",
                    f"        {str(sources['verifier'])!r},",
                    f"        [{str(sources['consumer'])!r}, {str(sources['solver'])!r}],",
                    "    )",
                    "    lifecycle.capture_environment_lock(",
                    "        sm_allowlist={'SM_WORKERS'},",
                    "        required_sm={'SM_WORKERS'},",
                    "        environment={",
                    "            'LC_ALL': 'C',",
                    "            'TZ': 'UTC',",
                    "            'PYTHONHASHSEED': '0',",
                    "            'SM_WORKERS': '1',",
                    "        },",
                    "        distributions=(),",
                    "        tools=(),",
                    f"        storage_path={str(self.run_parent)!r},",
                    "    )",
                    "    lifecycle.publish_ready()",
                    "    lifecycle.begin_running()",
                    "    lifecycle.children.launch(",
                    "        [sys.executable, '-c', 'import time; time.sleep(5)'],",
                    "        'real-signal-child',",
                    "        stdout=subprocess.DEVNULL,",
                    "        stderr=subprocess.DEVNULL,",
                    "    )",
                    f"    atomic_write_json({str(ready_path)!r}, {{'ready': True}})",
                    "    while True:",
                    "        time.sleep(1)",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        environment = dict(os.environ)
        module_dir = str(Path(__file__).resolve().parent)
        environment["PYTHONPATH"] = module_dir + os.pathsep + environment.get("PYTHONPATH", "")
        process = subprocess.Popen(
            [sys.executable, str(helper)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=environment,
        )
        try:
            deadline = time.monotonic() + 5
            while not ready_path.exists() and process.poll() is None and time.monotonic() < deadline:
                time.sleep(0.02)
            self.assertTrue(ready_path.exists(), "signal fixture did not reach RUNNING")
            os.kill(process.pid, signal.SIGTERM)
            returncode = process.wait(timeout=5)
            self.assertNotEqual(returncode, 0)
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=5)
        published = self.run_parent / "real-signal"
        state = json.loads((published / "run_state.json").read_text(encoding="utf-8"))
        children = json.loads((published / "children.json").read_text(encoding="utf-8"))["children"]
        self.assertEqual(state["state"], "ABORTED")
        self.assertTrue(all(child["returncode"] is not None for child in children))
        self.assertFalse((published / "_SUCCESS").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
