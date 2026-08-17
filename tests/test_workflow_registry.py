import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "handoff" / "build_workflow_registry.py"
REGISTRY = (
    ROOT
    / "docs"
    / "handoff"
    / "20260813_完整研究資料與軟體交接_01"
    / "registries"
    / "workflow_registry.json"
)
SPEC = importlib.util.spec_from_file_location("build_workflow_registry", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class WorkflowRegistryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
        cls.tracked = MODULE.git_tracked_regular_files(ROOT)
        cls.by_path = {row["path"]: row for row in cls.registry["records"]}

    def test_all_tracked_regular_scripts_are_registered_exactly_once(self):
        tracked_paths = [entry.path for entry in self.tracked]
        registered_paths = [row["path"] for row in self.registry["records"]]
        self.assertEqual(registered_paths, sorted(tracked_paths))
        self.assertEqual(len(registered_paths), len(set(registered_paths)))
        self.assertFalse(any(MODULE.is_pycache_path(path) for path in registered_paths))

    def test_record_metadata_matches_git_index_blobs(self):
        required = {
            "path",
            "type",
            "classification",
            "reason",
            "shebang",
            "suffix",
            "git_mode",
            "size_bytes",
            "sha256",
            "absolute_path_token_count",
            "known_limits",
        }
        for entry in self.tracked:
            row = self.by_path[entry.path]
            self.assertTrue(required.issubset(row))
            self.assertEqual(row["git_mode"], entry.mode)
            self.assertEqual(row["size_bytes"], len(entry.content))
            self.assertEqual(row["sha256"], hashlib.sha256(entry.content).hexdigest())
            self.assertEqual(
                row["absolute_path_token_count"], len(MODULE.ABSOLUTE_PATH_TOKEN_RE.findall(entry.content))
            )
            self.assertTrue(row["shebang"] or row["suffix"])
            self.assertIsInstance(row["known_limits"], list)
            self.assertTrue(row["known_limits"])

    def test_fixed_enum_and_supported_allowlist_fail_closed(self):
        self.assertEqual(self.registry["classification_enum"], list(MODULE.CLASSIFICATIONS))
        classifications = {row["classification"] for row in self.registry["records"]}
        self.assertLessEqual(classifications, set(MODULE.CLASSIFICATIONS))
        supported = {row["path"] for row in self.registry["records"] if row["classification"] == "SUPPORTED"}
        tracked_paths = {entry.path for entry in self.tracked}
        self.assertEqual(supported, set(MODULE.SUPPORTED_ALLOWLIST) & tracked_paths)
        self.assertNotIn(MODULE.UNSAFE_CLEANUP_PATH, supported)
        self.assertEqual(
            self.by_path[MODULE.UNSAFE_CLEANUP_PATH]["classification"], "REPRODUCIBLE_LEGACY"
        )
        self.assertIn("directly deletes", self.by_path[MODULE.UNSAFE_CLEANUP_PATH]["reason"])
        self.assertIn("--skip-cleanup", " ".join(self.by_path[MODULE.UNSAFE_CLEANUP_PATH]["known_limits"]))

    def test_legacy_and_archive_meanings_are_explicit(self):
        for row in self.registry["records"]:
            if row["path"].startswith(MODULE.ARCHIVED_PREFIX):
                self.assertEqual(row["classification"], "ARCHIVED")
            if row["classification"] == "REPRODUCIBLE_LEGACY":
                self.assertIn("does not mean broken", " ".join(row["known_limits"]))
        self.assertGreater(self.registry["summary"]["legacy_or_archived_records"], 200)

    def test_summary_and_hashes_match_records(self):
        records = self.registry["records"]
        counts = Counter(row["classification"] for row in records)
        expected_counts = {name: counts[name] for name in MODULE.CLASSIFICATIONS}
        self.assertEqual(self.registry["summary"]["total_records"], len(records))
        self.assertEqual(self.registry["summary"]["by_classification"], expected_counts)
        self.assertEqual(
            self.registry["summary"]["legacy_or_archived_records"],
            counts["REPRODUCIBLE_LEGACY"] + counts["ARCHIVED"],
        )
        self.assertEqual(self.registry["tree_content_hash"]["value"], MODULE.tree_content_sha256(records))
        self.assertEqual(self.registry["records_hash"]["value"], MODULE.records_sha256(records))

    def test_rebuild_is_byte_deterministic(self):
        first = MODULE.build_registry(ROOT, MODULE.DEFAULT_GENERATED_AT)
        second = MODULE.build_registry(ROOT, MODULE.DEFAULT_GENERATED_AT)
        self.assertEqual(MODULE.render_registry(first), MODULE.render_registry(second))
        self.assertEqual(MODULE.render_registry(first), REGISTRY.read_bytes())

        changed_timestamp = MODULE.build_registry(ROOT, "2026-08-13T01:00:00+08:00")
        self.assertNotEqual(first["generated_at"], changed_timestamp["generated_at"])
        self.assertEqual(first["tree_content_hash"], changed_timestamp["tree_content_hash"])
        self.assertEqual(first["records_hash"], changed_timestamp["records_hash"])

    def test_dirty_and_untracked_files_do_not_affect_inventory(self):
        with tempfile.TemporaryDirectory(prefix="workflow-registry-test-") as directory:
            repo = Path(directory)
            subprocess.run(["git", "init", "-q", str(repo)], check=True)
            script_dir = repo / "scripts"
            script_dir.mkdir()
            tracked = script_dir / "tracked.py"
            indexed_bytes = b"#!/usr/bin/env python3\nprint('indexed')\n"
            tracked.write_bytes(indexed_bytes)
            subprocess.run(["git", "-C", str(repo), "add", "scripts/tracked.py"], check=True)

            tracked.write_text("#!/usr/bin/env python3\nprint('dirty')\n", encoding="utf-8")
            (script_dir / "untracked.py").write_text("print('untracked')\n", encoding="utf-8")
            cache = script_dir / "__pycache__"
            cache.mkdir()
            (cache / "tracked.cpython-39.pyc").write_bytes(b"untracked bytecode")

            entries = MODULE.git_tracked_regular_files(repo)
            self.assertEqual([entry.path for entry in entries], ["scripts/tracked.py"])
            self.assertEqual(entries[0].content, indexed_bytes)


if __name__ == "__main__":
    unittest.main()
