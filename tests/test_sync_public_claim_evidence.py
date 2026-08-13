import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).parents[1]
SCRIPT = REPO_ROOT / "scripts" / "handoff" / "sync_public_claim_evidence.py"
SPEC = importlib.util.spec_from_file_location("sync_public_claim_evidence", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class PublicClaimEvidenceSyncTest(unittest.TestCase):
    def write_json(self, path: Path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    def build_fixture(self, root: Path):
        claims = [{"claim_id": f"C{index:03d}"} for index in range(1, 159)]
        registry = {
            "schema_name": "intersubmod.public_claim_remediation_registry",
            "schema_version": "2.0.0",
            "counts": {
                "claims": 158,
                "by_priority": {"P0": 34, "P1": 20, "P2": 35, "P3": 69},
                "by_current_verdict": {
                    "CONFIRMED": 69,
                    "CONFIRMED_WITH_LIMITS": 65,
                    "UNVERIFIED": 24,
                },
                "by_source_status": {
                    "EXTERNAL_LIVE_CONFIRMED_WITH_LIMITS": 1,
                    "FROZEN_CONFIRMED_WITH_LIMITS": 31,
                    "FROZEN_EVIDENCE_CONFIRMED": 69,
                    "SOURCE_EDITED_REVIEW_REQUIRED": 24,
                    "SOURCE_READY": 33,
                },
            },
            "gates": {
                gate: {"status": status, "reason": f"bounded {gate}"}
                for gate, status in MODULE.EXPECTED_GATE_STATUSES.items()
            },
            "claims": claims,
        }
        browser_receipt = {
            "schema_version": 3,
            "verdict": "PASS_FOR_LOCAL_SOURCE_QA",
        }
        runner_source = b"# validation runner fixture\n"
        browser_bytes = (json.dumps(browser_receipt, ensure_ascii=False, indent=2) + "\n").encode()
        registry_bytes = (json.dumps(registry, ensure_ascii=False, indent=2) + "\n").encode()
        import hashlib
        receipt = {
            "schema_name": "intersubmod.public_claim_validation_receipt",
            "schema_version": "2.0.0",
            "created_at": "2026-08-13T15:49:30+08:00",
            "scope": "LOCAL_SOURCE_PLUS_C108_LIVE_RECEIPT_PLUS_CHROMIUM_QA",
            "verdict": "PASS",
            "publication_status": "BLOCKED_DEFAULT_BRANCH_WIKI_PAGES_NOT_PUBLISHED_OR_REFETCHED__ABOUT_C108_CONFIRMED",
            "release_status": "BLOCKED",
            "runner": {
                "path": MODULE.SOURCE_VALIDATION_RUNNER_RELATIVE.as_posix(),
                "sha256": hashlib.sha256(runner_source).hexdigest(),
                "python": "3.10.12",
            },
            "claim_registry_contract": {
                "registry_sha256": hashlib.sha256(registry_bytes).hexdigest(),
                "public_source_files": 34,
                "verdict": "PASS",
            },
            "browser_qa_contract": {
                "receipt_sha256": hashlib.sha256(browser_bytes).hexdigest(),
                "html_files": 21,
                "standalone_svg_files": 21,
                "browser_cases": 84,
                "verdict": "PASS",
            },
        }
        pointer = {
            "registry_type": "claim_registry_snapshot_pointer",
            "schema_version": "1.0.0",
            "snapshot_id": "intersubmod-public-claims-remediation-20260813",
            "records_uri": MODULE.EXPECTED_RECORDS_URI,
            "records_sha256": "0" * 64,
            "records_count": 0,
            "materialization": {"preserve_me": True},
            "counts": {},
            "gates": {},
            "claim_ceiling": "bounded",
            "verified_at": "old",
        }
        manifest = {
            "manifest_type": "research_handoff_embedded_evidence",
            "schema_version": "1.0.0",
            "records": [
                {
                    "evidence_id": "unrelated",
                    "path": "evidence/unrelated.json",
                    "sha256": "a" * 64,
                    "copy_status": "NEW_PACKAGE_RECEIPT",
                    "sentinel": "preserve",
                },
                {
                    "evidence_id": "public_claim_validation_20260813",
                    "path": "evidence/public_claim_validation_receipt.json",
                    "sha256": "b" * 64,
                    "copy_status": "EXACT_COPY",
                    "source_path": "InterSubMod/research/20260813_public_docs_p0_correction/validation_receipt.json",
                },
                {
                    "evidence_id": "full_claim_registry_20260813",
                    "path": "evidence/full_claim_registry.json",
                    "sha256": "c" * 64,
                    "copy_status": "EXACT_COPY",
                    "source_path": "InterSubMod/research/20260813_public_docs_p0_correction/claim_remediation_registry.json",
                },
            ],
        }

        source_registry = root / MODULE.SOURCE_REGISTRY_RELATIVE
        source_receipt = root / MODULE.SOURCE_RECEIPT_RELATIVE
        package = root / MODULE.PACKAGE_RELATIVE
        self.write_json(source_registry, registry)
        self.write_json(source_receipt, receipt)
        browser_path = root / MODULE.SOURCE_BROWSER_RECEIPT_RELATIVE
        validation_runner = root / MODULE.SOURCE_VALIDATION_RUNNER_RELATIVE
        self.write_json(browser_path, browser_receipt)
        validation_runner.parent.mkdir(parents=True, exist_ok=True)
        validation_runner.write_bytes(runner_source)
        self.write_json(package / MODULE.POINTER_RELATIVE, pointer)
        self.write_json(package / MODULE.MANIFEST_RELATIVE, manifest)
        return source_registry, source_receipt, package

    def test_sync_exact_copies_and_derives_pointer_and_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_registry, source_receipt, package = self.build_fixture(root)

            result = MODULE.sync_evidence(root)

            target_registry = package / MODULE.TARGET_REGISTRY_RELATIVE
            target_receipt = package / MODULE.TARGET_RECEIPT_RELATIVE
            self.assertEqual(target_registry.read_bytes(), source_registry.read_bytes())
            self.assertEqual(target_receipt.read_bytes(), source_receipt.read_bytes())
            self.assertTrue(result["outputs"][MODULE.TARGET_REGISTRY_RELATIVE.as_posix()]["exact_copy"])

            pointer = json.loads((package / MODULE.POINTER_RELATIVE).read_text(encoding="utf-8"))
            registry = json.loads(source_registry.read_text(encoding="utf-8"))
            self.assertEqual(pointer["records_sha256"], MODULE.sha256_bytes(source_registry.read_bytes()))
            self.assertEqual(pointer["records_count"], 158)
            self.assertEqual(pointer["counts"], {key: value for key, value in registry["counts"].items() if key != "claims"})
            self.assertEqual(pointer["gates"], MODULE.EXPECTED_GATE_STATUSES)
            self.assertEqual(pointer["verified_at"], "2026-08-13T15:49:30+08:00")
            self.assertTrue(pointer["materialization"]["preserve_me"])
            self.assertEqual(
                pointer["materialization"]["source_patch_delivery"],
                "PUBLIC_CLAIM_CORRECTIONS_WORKING_STACK_VALIDATED",
            )
            self.assertIn("P0_SOURCE_READY is PASS", pointer["materialization"]["known_limit"])
            self.assertIn("RELEASE_READY BLOCKED", pointer["materialization"]["known_limit"])

            manifest = json.loads((package / MODULE.MANIFEST_RELATIVE).read_text(encoding="utf-8"))
            by_id = {row["evidence_id"]: row for row in manifest["records"]}
            self.assertEqual(by_id["unrelated"]["sentinel"], "preserve")
            for evidence_id, source in (
                ("full_claim_registry_20260813", source_registry),
                ("public_claim_validation_20260813", source_receipt),
            ):
                self.assertEqual(by_id[evidence_id]["sha256"], MODULE.sha256_bytes(source.read_bytes()))
                self.assertEqual(by_id[evidence_id]["size_bytes"], len(source.read_bytes()))
                self.assertEqual(by_id[evidence_id]["copy_status"], "EXACT_COPY")

    def test_invalid_receipt_fails_before_any_package_write(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, source_receipt, package = self.build_fixture(root)
            receipt = json.loads(source_receipt.read_text(encoding="utf-8"))
            receipt["release_status"] = "PASS"
            self.write_json(source_receipt, receipt)
            pointer_before = (package / MODULE.POINTER_RELATIVE).read_bytes()
            manifest_before = (package / MODULE.MANIFEST_RELATIVE).read_bytes()

            with self.assertRaisesRegex(MODULE.SyncError, "release_status must remain BLOCKED"):
                MODULE.sync_evidence(root)

            self.assertEqual((package / MODULE.POINTER_RELATIVE).read_bytes(), pointer_before)
            self.assertEqual((package / MODULE.MANIFEST_RELATIVE).read_bytes(), manifest_before)
            self.assertFalse((package / MODULE.TARGET_REGISTRY_RELATIVE).exists())
            self.assertFalse((package / MODULE.TARGET_RECEIPT_RELATIVE).exists())

    def test_non_exact_manifest_record_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, _, package = self.build_fixture(root)
            manifest_path = package / MODULE.MANIFEST_RELATIVE
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            row = next(item for item in manifest["records"] if item["evidence_id"] == "full_claim_registry_20260813")
            row["copy_status"] = "SUMMARY_HASH_BOUND"
            self.write_json(manifest_path, manifest)

            with self.assertRaisesRegex(MODULE.SyncError, "not EXACT_COPY"):
                MODULE.sync_evidence(root)

    def test_gate_promotion_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_registry, _, _ = self.build_fixture(root)
            registry = json.loads(source_registry.read_text(encoding="utf-8"))
            registry["gates"]["RELEASE_READY"]["status"] = "PASS"
            self.write_json(source_registry, registry)

            with self.assertRaisesRegex(MODULE.SyncError, "not fail-closed"):
                MODULE.sync_evidence(root)


if __name__ == "__main__":
    unittest.main()
