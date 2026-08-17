import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "handoff" / "replay_authority.py"
RUNBOOK = (
    Path(__file__).parents[1]
    / "docs/handoff/20260813_完整研究資料與軟體交接_01/20260813_bip7_bip8操作與驗證_01.md"
)
SPEC = importlib.util.spec_from_file_location("replay_authority", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class AuthorityReplayTest(unittest.TestCase):
    def make_manifest(self, root: Path, expected_hash: str) -> Path:
        manifest = {
            "as_of_date": "2026-08-01",
            "artifacts": [
                {
                    "artifact_id": "fixture",
                    "path": str(root / "fixture.txt"),
                    "sha256": expected_hash,
                }
            ],
            "implementation": {"source_snapshots": []},
        }
        path = root / "authority_manifest.json"
        path.write_text(json.dumps(manifest), encoding="utf-8")
        return path

    def test_match(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = root / "fixture.txt"
            fixture.write_text("frozen\n", encoding="utf-8")
            manifest = self.make_manifest(root, MODULE.sha256_file(fixture))
            receipt = MODULE.replay(manifest)
            self.assertTrue(receipt["pass"])
            self.assertEqual(receipt["tally"], {"MATCH": 1, "MISSING": 0, "HASH_MISMATCH": 0})

    def test_mismatch_is_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "fixture.txt").write_text("changed\n", encoding="utf-8")
            manifest = self.make_manifest(root, "0" * 64)
            receipt = MODULE.replay(manifest)
            self.assertFalse(receipt["pass"])
            self.assertEqual(receipt["tally"]["HASH_MISMATCH"], 1)

    def test_missing_is_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = self.make_manifest(root, "0" * 64)
            receipt = MODULE.replay(manifest)
            self.assertFalse(receipt["pass"])
            self.assertEqual(receipt["tally"]["MISSING"], 1)

    def test_public_runbook_uses_manifest_cli_flag(self):
        text = RUNBOOK.read_text(encoding="utf-8")
        replay_command = text.split("python3 scripts/handoff/replay_authority.py", 1)[1].split("```", 1)[0]
        self.assertIn("--manifest ", replay_command)
        self.assertNotIn("--authority ", replay_command)


if __name__ == "__main__":
    unittest.main()
