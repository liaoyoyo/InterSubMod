#!/usr/bin/env python3
"""Synthetic fail-closed tests for the presentation snapshot contract."""

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import stat
import sys
import tempfile
import unittest
from unittest import mock


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import freeze_presentation_snapshot as freezer  # noqa: E402


def identity(path: str, payload: bytes) -> dict[str, object]:
    return {"path": path, "size_bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()}


SYNTHETIC_RUNTIME = {
    "python": {
        "executable": identity("/synthetic/python", b"python"),
        "version": "0.0",
        "implementation": "synthetic",
    },
    "playwright": {
        "version": "0.0",
        "module": identity("/synthetic/playwright.py", b"playwright"),
    },
    "chromium": {
        "executable": identity("/synthetic/chromium", b"chromium"),
    },
    "platform": "synthetic-test-platform",
}


def make_writable(root: pathlib.Path) -> None:
    if not root.exists() or root.is_symlink():
        return
    for path in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        if path.is_symlink():
            continue
        try:
            path.chmod(0o755 if path.is_dir() else 0o644)
        except FileNotFoundError:
            pass
    root.chmod(0o755)


class PresentationSnapshotFixture(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.cleanup)
        self.base = pathlib.Path(self.temporary.name)
        self.repo = self.base / "repo"
        self.repo.mkdir()
        for index, (role, relative) in enumerate(freezer.SOURCE_ALLOWLIST):
            path = self.repo / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"# synthetic {role}\nVALUE = {index}\n", encoding="utf-8")
        self.evidence: dict[str, pathlib.Path] = {}
        evidence_dir = self.base / "evidence"
        evidence_dir.mkdir()
        for index, role in enumerate(freezer.EVIDENCE_ROLES):
            path = evidence_dir / f"{index:02d}_{role}.dat"
            path.write_text(f"synthetic evidence {role}\n", encoding="utf-8")
            self.evidence[role] = path
        self.output = self.base / "presentation_snapshot"

    def cleanup(self) -> None:
        make_writable(self.base)
        self.temporary.cleanup()

    def freeze(self) -> dict:
        return freezer.freeze_presentation_snapshot(
            self.evidence,
            self.output,
            _repo_root=self.repo,
            _runtime=SYNTHETIC_RUNTIME,
        )

    def verify(self, runtime: dict | None = None) -> dict:
        return freezer.verify_presentation_snapshot(
            self.output,
            _runtime=SYNTHETIC_RUNTIME if runtime is None else runtime,
        )


class PositiveTests(PresentationSnapshotFixture):
    def test_freeze_and_verify_exact_contract(self) -> None:
        result = self.freeze()
        self.assertTrue(result["all_pass"])
        self.assertEqual(result["programs_frozen"], 3)
        self.assertEqual(result["evidence_inputs_bound"], 12)
        verification = self.verify()
        self.assertTrue(verification["all_pass"])
        self.assertEqual(verification["programs_verified"], 3)
        self.assertEqual(verification["evidence_inputs_verified"], 12)
        self.assertEqual(stat.S_IMODE(self.output.stat().st_mode), 0o555)
        for path in self.output.rglob("*"):
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o555 if path.is_dir() else 0o444)
        manifest = json.loads((self.output / freezer.MANIFEST_NAME).read_text(encoding="utf-8"))
        self.assertEqual(manifest["artifact_class"], freezer.ARTIFACT_CLASS)
        self.assertFalse(manifest["assurance"]["scientific_evidence_authority"])
        self.assertEqual([row["role"] for row in manifest["evidence_inputs"]["entries"]], list(freezer.EVIDENCE_ROLES))

    def test_verify_only_does_not_modify_snapshot(self) -> None:
        self.freeze()
        before = {
            path.relative_to(self.output).as_posix(): (
                path.stat().st_size, path.stat().st_mtime_ns, path.stat().st_ctime_ns
            )
            for path in self.output.rglob("*") if path.is_file()
        }
        self.verify()
        after = {
            path.relative_to(self.output).as_posix(): (
                path.stat().st_size, path.stat().st_mtime_ns, path.stat().st_ctime_ns
            )
            for path in self.output.rglob("*") if path.is_file()
        }
        self.assertEqual(before, after)


class FailClosedTests(PresentationSnapshotFixture):
    def test_output_root_must_not_exist(self) -> None:
        self.output.mkdir()
        with self.assertRaisesRegex(freezer.PresentationSnapshotError, "must not exist"):
            self.freeze()

    def test_wrong_evidence_role_set_is_rejected(self) -> None:
        evidence = dict(self.evidence)
        evidence.pop("final_numeric_summary")
        with self.assertRaisesRegex(freezer.PresentationSnapshotError, "role allowlist"):
            freezer.freeze_presentation_snapshot(
                evidence, self.output, _repo_root=self.repo, _runtime=SYNTHETIC_RUNTIME
            )

    def test_symlink_evidence_is_rejected(self) -> None:
        original = self.evidence["canonical_json"]
        link = self.base / "canonical_link"
        link.symlink_to(original)
        self.evidence["canonical_json"] = link
        with self.assertRaisesRegex(freezer.PresentationSnapshotError, "symlink"):
            self.freeze()

    def test_changed_evidence_fails_verify(self) -> None:
        self.freeze()
        self.evidence["final_numeric_summary"].write_text("changed\n", encoding="utf-8")
        with self.assertRaisesRegex(freezer.PresentationSnapshotError, "identity changed"):
            self.verify()

    def test_changed_runtime_fails_verify(self) -> None:
        self.freeze()
        runtime = json.loads(json.dumps(SYNTHETIC_RUNTIME))
        runtime["playwright"]["version"] = "different"
        with self.assertRaisesRegex(freezer.PresentationSnapshotError, "runtime differs"):
            self.verify(runtime)

    def test_snapshot_program_tamper_fails_verify(self) -> None:
        self.freeze()
        program = self.output / "code_snapshot" / "build_validated_html_report.py"
        program.chmod(0o644)
        program.write_text("tampered\n", encoding="utf-8")
        with self.assertRaisesRegex(freezer.PresentationSnapshotError, "snapshot"):
            self.verify()

    def test_extra_snapshot_file_fails_verify(self) -> None:
        self.freeze()
        self.output.chmod(0o755)
        extra = self.output / "extra.txt"
        extra.write_text("extra\n", encoding="utf-8")
        extra.chmod(0o444)
        self.output.chmod(0o555)
        with self.assertRaisesRegex(freezer.PresentationSnapshotError, "file set"):
            self.verify()

    def test_atomic_publish_failure_preserves_named_readonly_staging(self) -> None:
        real_rename = freezer._rename_noreplace
        calls = 0

        def fail_publish_once(source: pathlib.Path, destination: pathlib.Path) -> None:
            nonlocal calls
            calls += 1
            if calls == 1:
                raise freezer.PresentationSnapshotError("synthetic publish failure")
            real_rename(source, destination)

        with mock.patch.object(freezer, "_rename_noreplace", side_effect=fail_publish_once):
            with self.assertRaisesRegex(freezer.PresentationSnapshotError, "synthetic publish failure"):
                self.freeze()
        self.assertFalse(os.path.lexists(self.output))
        preserved = list(self.base.glob(f"{self.output.name}.failed-staging.*"))
        self.assertEqual(len(preserved), 1)
        self.assertEqual(list(self.base.glob(f".{self.output.name}.staging.*")), [])
        self.assertEqual(stat.S_IMODE(preserved[0].stat().st_mode), 0o555)
        verification = freezer.verify_presentation_snapshot(
            preserved[0], _runtime=SYNTHETIC_RUNTIME
        )
        self.assertTrue(verification["all_pass"])

    def test_archive_rename_failure_leaves_original_readonly_staging(self) -> None:
        with mock.patch.object(
            freezer,
            "_rename_noreplace",
            side_effect=freezer.PresentationSnapshotError("all renames unavailable"),
        ):
            with self.assertRaisesRegex(
                freezer.PresentationSnapshotError,
                r"failed staging preserved at .*\.staging\..*preservation note",
            ):
                self.freeze()
        self.assertFalse(os.path.lexists(self.output))
        staging = list(self.base.glob(f".{self.output.name}.staging.*"))
        self.assertEqual(len(staging), 1)
        self.assertEqual(stat.S_IMODE(staging[0].stat().st_mode), 0o555)


if __name__ == "__main__":
    unittest.main()
