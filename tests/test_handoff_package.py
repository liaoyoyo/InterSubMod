import contextlib
import importlib.util
import io
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).parents[1]
SCRIPT = REPO_ROOT / "scripts" / "handoff" / "validate_handoff_package.py"
PACKAGE_ROOT = REPO_ROOT / "docs" / "handoff" / "20260813_完整研究資料與軟體交接_01"
SPEC = importlib.util.spec_from_file_location("validate_handoff_package", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class HandoffPackageValidationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.valid_receipt = MODULE.validate_package(PACKAGE_ROOT)

    def copy_package(self, directory: str) -> Path:
        fixture_root = Path(directory) / "repo"
        destination = fixture_root / "docs" / "handoff" / PACKAGE_ROOT.name
        destination.parent.mkdir(parents=True)
        shutil.copytree(PACKAGE_ROOT, destination)
        return destination

    @staticmethod
    def failed_check(receipt, check_id):
        return next(row for row in receipt["checks"] if row["check_id"] == check_id)

    def test_complete_package_passes_all_gates(self):
        self.assertTrue(self.valid_receipt["pass"], self.valid_receipt)
        self.assertEqual(self.valid_receipt["tally"], {"PASS": 13, "FAIL": 0})

    def test_receipt_exposes_required_counts(self):
        by_id = {row["check_id"]: row for row in self.valid_receipt["checks"]}
        self.assertEqual(by_id["dataset_registry"]["details"]["technical_datasets"], 7)
        self.assertEqual(by_id["dataset_registry"]["details"]["biological_samples"], 6)
        self.assertEqual(by_id["run_registry"]["details"]["physical_runs"], 51)
        self.assertEqual(by_id["artifact_registry"]["details"]["artifacts"], 36)
        self.assertEqual(by_id["artifact_registry"]["details"]["final_for_scope"], 20)
        self.assertEqual(by_id["artifact_registry"]["details"]["tagged_bam_partial_non_final"], 14)
        self.assertEqual(by_id["claim_registries"]["details"]["claims"], 158)
        self.assertEqual(by_id["authority_replay"]["details"]["match"], 19)

    def test_broken_markdown_link_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            index = package / "00_INDEX.md"
            index.write_text(index.read_text(encoding="utf-8") + "\n[broken](missing.json)\n", encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "markdown_links")["status"], "FAIL")
            self.assertFalse(receipt["pass"])

    def test_invalid_json_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            (package / "registries/run_registry.json").write_text("{", encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "json_parse")["status"], "FAIL")
            self.assertFalse(receipt["pass"])

    def test_evidence_hash_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            evidence = package / "evidence/denominator_registry.tsv"
            evidence.write_text(evidence.read_text(encoding="utf-8") + "\nmutation\n", encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "evidence_manifest")["status"], "FAIL")

    def test_summary_source_hash_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            source = package / "evidence/repo_hygiene_full_receipt.json"
            source.write_text(source.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "evidence_manifest")["status"], "FAIL")

    def test_duplicate_physical_run_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            path = package / "registries/run_registry.json"
            registry = json.loads(path.read_text(encoding="utf-8"))
            registry["records"][1]["physical_run_id"] = registry["records"][0]["physical_run_id"]
            path.write_text(json.dumps(registry), encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "run_registry")["status"], "FAIL")

    def test_final_artifact_without_producer_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            path = package / "registries/artifact_registry.json"
            registry = json.loads(path.read_text(encoding="utf-8"))
            final = next(row for row in registry if row["finality"] == "FINAL_FOR_SCOPE")
            final["producer"] = ""
            path.write_text(json.dumps(registry), encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "artifact_registry")["status"], "FAIL")

    def test_tagged_bam_promoted_to_final_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            path = package / "registries/artifact_registry.json"
            registry = json.loads(path.read_text(encoding="utf-8"))
            tagged_bam = next(row for row in registry if row["artifact_type"] == "longphase_s_tagged_bam")
            tagged_bam["finality"] = "FINAL_FOR_SCOPE"
            path.write_text(json.dumps(registry), encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "artifact_registry")["status"], "FAIL")

    def test_claim_pointer_hash_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            path = package / "registries/claim_registry.json"
            pointer = json.loads(path.read_text(encoding="utf-8"))
            pointer["records_sha256"] = "0" * 64
            path.write_text(json.dumps(pointer), encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "claim_registries")["status"], "FAIL")

    def test_unregistered_site_profile_alias_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            path = package / "registries/dataset_alias_registry.json"
            aliases = json.loads(path.read_text(encoding="utf-8"))
            aliases["aliases"] = []
            path.write_text(json.dumps(aliases), encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "site_profile_join")["status"], "FAIL")

    def test_external_html_runtime_resource_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            path = package / "20260813_完整研究交接總覽_01.standalone.html"
            text = path.read_text(encoding="utf-8")
            text = text.replace("</body>", '<script src="https://example.invalid/runtime.js"></script></body>')
            path.write_text(text, encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "standalone_html")["status"], "FAIL")

    def test_authority_replay_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            path = package / "evidence/authority_replay_receipt.json"
            replay = json.loads(path.read_text(encoding="utf-8"))
            replay["results"][0]["status"] = "HASH_MISMATCH"
            replay["tally"] = {"MATCH": 18, "MISSING": 0, "HASH_MISMATCH": 1}
            replay["pass"] = False
            path.write_text(json.dumps(replay), encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "authority_replay")["status"], "FAIL")

    def test_no_receipt_argument_does_not_write_package(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            before = {path.relative_to(package): path.stat().st_mtime_ns for path in package.rglob("*") if path.is_file()}
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                exit_code = MODULE.main([str(package)])
            after = {path.relative_to(package): path.stat().st_mtime_ns for path in package.rglob("*") if path.is_file()}
            self.assertEqual(exit_code, 0)
            self.assertEqual(before, after)
            self.assertTrue(json.loads(output.getvalue())["pass"])

    def test_reader_contract_rejects_p0_denominator_promotion(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            index = package / "00_INDEX.md"
            index.write_text(index.read_text(encoding="utf-8").replace("33/33", "34/34", 1), encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "reader_contract")["status"], "FAIL")

    def test_reader_contract_rejects_invalid_bootstrap_interface(self):
        with tempfile.TemporaryDirectory() as directory:
            package = self.copy_package(directory)
            html = package / "20260813_完整研究交接總覽_01.standalone.html"
            text = html.read_text(encoding="utf-8").replace("scripts/site/bootstrap --template", "scripts/site/bootstrap --profile", 1)
            html.write_text(text, encoding="utf-8")
            receipt = MODULE.validate_package(package)
            self.assertEqual(self.failed_check(receipt, "reader_contract")["status"], "FAIL")

    def test_explicit_receipt_path_is_written(self):
        with tempfile.TemporaryDirectory() as directory:
            receipt_path = Path(directory) / "receipts/handoff.json"
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                exit_code = MODULE.main([str(PACKAGE_ROOT), "--receipt", str(receipt_path)])
            self.assertEqual(exit_code, 0)
            self.assertTrue(receipt_path.is_file())
            self.assertTrue(json.loads(receipt_path.read_text(encoding="utf-8"))["pass"])


if __name__ == "__main__":
    unittest.main()
