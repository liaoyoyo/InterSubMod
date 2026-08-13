import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "handoff" / "repo_hygiene.py"
SPEC = importlib.util.spec_from_file_location("repo_hygiene", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class RepoHygieneTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.repo = self.root / "repo"
        self.archive = self.root / "archive"
        self.repo.mkdir()
        subprocess.run(["git", "init", "-q", self.repo], check=True)
        subprocess.run(["git", "-C", self.repo, "config", "user.email", "fixture@example.invalid"], check=True)
        subprocess.run(["git", "-C", self.repo, "config", "user.name", "Fixture"], check=True)

        (self.repo / "source").mkdir()
        (self.repo / "source" / "tool.py").write_text("print('fixture')\n", encoding="utf-8")
        (self.repo / "external").mkdir()
        (self.repo / "external" / "figure.png").write_bytes(b"fixture-image")
        (self.repo / "links").mkdir()
        os.symlink("/proc/self/fd/999", self.repo / "links" / "runtime")
        os.symlink(str(self.repo / "source" / "tool.py"), self.repo / "links" / "tool.py")
        os.symlink(str(self.repo / "external"), self.repo / "links" / "figures")
        os.symlink("../missing", self.repo / "links" / "broken")
        (self.repo / ".claude").mkdir()
        (self.repo / ".claude" / "settings.local.json").write_text('{"local": true}\n', encoding="utf-8")
        (self.repo / ".claude" / "settings.local.example.json").write_text("{}\n", encoding="utf-8")
        subprocess.run(["git", "-C", self.repo, "add", "."], check=True)
        subprocess.run(
            ["git", "-C", self.repo, "add", "-f", ".claude/settings.local.json"], check=True
        )
        subprocess.run(["git", "-C", self.repo, "commit", "-qm", "fixture"], check=True)

        self.rules = {
            "links/tool.py": MODULE.LinkRule(
                action="convert_to_relative",
                expected_target=str(self.repo / "source" / "tool.py"),
                replacement_target="../source/tool.py",
            ),
            "links/figures": MODULE.LinkRule(
                action="replace_with_pointer",
                expected_target=str(self.repo / "external"),
                semantic_description="Fixture figures.",
                availability="LOCAL_DERIVED_UNTRACKED",
                claim_ceiling="Fixture only.",
            ),
            "links/broken": MODULE.LinkRule(
                action="remove_broken",
                expected_target="../missing",
                availability="MISSING",
                claim_ceiling="No payload.",
            ),
        }

    def tearDown(self):
        self.temporary.cleanup()

    def make_plan(self):
        return MODULE.build_plan(
            self.repo,
            rules=self.rules,
            expected_counts={
                "tracked_symlinks": 4,
                "remove_proc_self_fd": 1,
                "remove_broken": 1,
                "convert_to_relative": 1,
                "replace_with_pointer": 1,
            },
            proc_source_prefix="links/",
        )

    def test_plan_classifies_exact_inventory_without_hashing_external_payload(self):
        plan = self.make_plan()
        actions = {operation["path"]: operation["action"] for operation in plan}
        self.assertEqual(actions["links/runtime"], "remove_proc_self_fd")
        self.assertEqual(actions["links/tool.py"], "convert_to_relative")
        pointer = next(operation for operation in plan if operation["path"] == "links/figures")
        self.assertEqual(pointer["target_inventory"]["file_count"], 1)
        self.assertIn("file contents are intentionally not read", pointer["target_inventory"]["hash_strategy"])

    def test_unreviewed_symlink_fails_closed(self):
        (self.repo / "links" / "unexpected").symlink_to("../source/tool.py")
        subprocess.run(["git", "-C", self.repo, "add", "links/unexpected"], check=True)
        with self.assertRaises(MODULE.HygieneError):
            self.make_plan()

    def test_archive_apply_and_restore_round_trip(self):
        baseline = subprocess.check_output(["git", "-C", self.repo, "rev-parse", "HEAD"], text=True).strip()
        plan = self.make_plan()
        before = MODULE.make_before_receipt(self.repo, baseline, plan)
        _, manifest = MODULE.build_archive(self.repo, self.archive, before)
        self.assertTrue(MODULE.verify_archive(self.archive)["pass"])

        MODULE.apply_operations(self.repo, manifest)
        self.assertFalse((self.repo / "links" / "runtime").is_symlink())
        self.assertEqual(os.readlink(self.repo / "links" / "tool.py"), "../source/tool.py")
        self.assertTrue((self.repo / "links" / "figures" / "ARTIFACT_POINTER.json").is_file())
        self.assertFalse((self.repo / ".claude" / "settings.local.json").exists())

        restored = MODULE.restore_from_archive(self.repo, self.archive, confirm=True)
        self.assertTrue(restored["pass"])
        self.assertEqual(os.readlink(self.repo / "links" / "runtime"), "/proc/self/fd/999")
        self.assertEqual(os.readlink(self.repo / "links" / "tool.py"), str(self.repo / "source" / "tool.py"))
        restored_settings = self.repo / ".claude" / "settings.local.json"
        self.assertEqual(MODULE.sha256_file(restored_settings), before["settings"]["sha256"])

    def test_archive_refuses_overwrite(self):
        baseline = subprocess.check_output(["git", "-C", self.repo, "rev-parse", "HEAD"], text=True).strip()
        before = MODULE.make_before_receipt(self.repo, baseline, self.make_plan())
        MODULE.build_archive(self.repo, self.archive, before)
        with self.assertRaises(MODULE.HygieneError):
            MODULE.build_archive(self.repo, self.archive, before)

    def test_restore_recovers_a_partially_applied_cleanup(self):
        baseline = subprocess.check_output(["git", "-C", self.repo, "rev-parse", "HEAD"], text=True).strip()
        before = MODULE.make_before_receipt(self.repo, baseline, self.make_plan())
        MODULE.build_archive(self.repo, self.archive, before)
        (self.repo / "links" / "runtime").unlink()

        restored = MODULE.restore_from_archive(self.repo, self.archive, confirm=True)
        self.assertTrue(restored["pass"])
        self.assertEqual(restored["restored_symlinks"], 1)
        self.assertEqual(restored["already_original_symlinks"], 3)
        self.assertEqual(os.readlink(self.repo / "links" / "runtime"), "/proc/self/fd/999")

    def test_pointer_payload_is_machine_readable(self):
        operation = next(item for item in self.make_plan() if item["action"] == "replace_with_pointer")
        payload = MODULE.pointer_payload(operation, "fixture-head")
        encoded = json.dumps(payload)
        self.assertIn("metadata_sha256_v1", encoded)
        self.assertEqual(payload["claim_ceiling"], "Fixture only.")


if __name__ == "__main__":
    unittest.main()
