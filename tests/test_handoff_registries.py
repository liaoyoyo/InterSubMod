import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "handoff" / "build_registries.py"
REGISTRY_ROOT = (
    Path(__file__).parents[1]
    / "docs"
    / "handoff"
    / "20260813_完整研究資料與軟體交接_01"
    / "registries"
)
SPEC = importlib.util.spec_from_file_location("build_registries", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class HandoffRegistryTest(unittest.TestCase):
    def test_generated_registry_contract(self):
        run_registry = json.loads((REGISTRY_ROOT / "run_registry.json").read_text(encoding="utf-8"))
        artifacts = json.loads((REGISTRY_ROOT / "artifact_registry.json").read_text(encoding="utf-8"))
        machine_paths = json.loads((REGISTRY_ROOT / "machine_path_registry.json").read_text(encoding="utf-8"))

        expected_reconciliation = {
            "manifest_file_lines": 19,
            "logical_manifest_rows": 18,
            "physical_directories_total": 51,
            "current_physical_directories": 35,
            "pending_archive_physical_directories": 16,
            "logical_rows_merged_current": 9,
            "logical_rows_merged_pending_archive": 9,
            "current_physical_unregistered": 26,
            "pending_archive_extras": 7,
        }
        reconciliation = run_registry["reconciliation"]
        for key, value in expected_reconciliation.items():
            self.assertEqual(reconciliation[key], value)

        run_records = run_registry["records"]
        self.assertEqual(len(run_records), 51)
        self.assertEqual(len({row["physical_run_id"] for row in run_records}), 51)
        self.assertEqual(sum(row["logical_manifest_member"] for row in run_records), 18)

        tagged_bams = [row for row in artifacts if row["artifact_type"] == "longphase_s_tagged_bam"]
        self.assertEqual(len(tagged_bams), 14)
        self.assertEqual(sum(row["size_bytes"] for row in tagged_bams), 3_709_322_840_333)
        self.assertTrue(all(row["sha256"] is None for row in tagged_bams))
        self.assertTrue(all(row["filesystem_mtime"] for row in tagged_bams))
        for row in tagged_bams:
            bam_path = Path(row["machine_locations"][0]["path"])
            self.assertEqual(row["filesystem_mtime"], MODULE.iso_from_stat(bam_path))
        self.assertTrue(all(row["scope"] == "PARTIAL" for row in tagged_bams))
        self.assertTrue(all(row["evidence_status"] == "PARTIAL" for row in tagged_bams))
        self.assertTrue(all(row["finality"] == "NON_FINAL" for row in tagged_bams))

        machine_by_id = {row["path_id"]: row for row in machine_paths["records"]}
        self.assertEqual(machine_by_id["legacy.runbook"]["observation_status"], "MISSING")
        self.assertEqual(machine_by_id["legacy.big8_output"]["observation_status"], "NFS_VISIBLE_NOT_HOST_VERIFIED")
        self.assertEqual(machine_by_id["legacy.bip8_output"]["observation_status"], "NFS_VISIBLE_NOT_HOST_VERIFIED")

        artifacts_by_id = {row["artifact_id"]: row for row in artifacts}
        source4 = artifacts_by_id["authority.source_snapshot_4"]
        frozen = artifacts_by_id["authority.frozen_binary"]
        self.assertEqual(source4["used_by"], [])
        self.assertNotIn("authority.source_snapshot_4", frozen["inputs"])

    def test_complete_matrix_classification(self):
        self.assertEqual(MODULE.classify_run(Path("/x/20260212_H1437_paired_full_complete_matrix")), "COMPLETE_MATRIX")

    def test_pending_path_does_not_duplicate_canonical_component(self):
        canonical_root = Path("/output/canonical")
        pending_root = Path("/output/_ARCHIVE_pending_cleanup_202606/canonical")
        canonical_path = canonical_root / "HCC1395/paired_full/run1"
        self.assertEqual(
            MODULE.pending_path_for_manifest(canonical_path, canonical_root, pending_root),
            pending_root / "HCC1395/paired_full/run1",
        )

    def test_invalidated_claim_is_not_final(self):
        row = MODULE.artifact_row(
            artifact_id="claim.invalid",
            artifact_type="invalidated_claim",
            semantic_description="fixture",
            producer="test",
            claim_ceiling="invalid",
            evidence_status="INVALIDATED",
            finality="SUPERSEDED",
        )
        self.assertEqual(row["finality"], "SUPERSEDED")
        self.assertIsNone(row["sha256"])

    def test_final_artifact_without_hash_fails(self):
        artifact = MODULE.artifact_row(
            artifact_id="bad.final",
            artifact_type="fixture",
            semantic_description="fixture",
            producer="test",
            claim_ceiling="fixture",
            finality="FINAL_FOR_SCOPE",
        )
        datasets = [
            {"technical_dataset_id": name, "biological_id": "HCC1395" if name == "HCC1395_DORADO" else name}
            for name in MODULE.TECHNICAL_TO_BIOLOGICAL
        ]
        run_registry = {
            "reconciliation": {
                "manifest_file_lines": 19,
                "logical_manifest_rows": 18,
                "physical_directories_total": 51,
                "current_physical_directories": 35,
                "pending_archive_physical_directories": 16,
                "logical_rows_merged_current": 9,
                "logical_rows_merged_pending_archive": 9,
                "current_physical_unregistered": 26,
                "pending_archive_extras": 7,
            },
            "records": [
                {"physical_run_id": f"p{i}", "logical_manifest_member": i < 18} for i in range(51)
            ],
        }
        machine = {"records": [{"path_id": "legacy.runbook", "observation_status": "MISSING"}]}
        storage = {"records": [{"root_id": "fixture"}]}
        crosswalk = {"crosswalks": [{"crosswalk_id": "fixture"}]}
        source4 = MODULE.artifact_row(
            artifact_id="authority.source_snapshot_4",
            artifact_type="source_snapshot",
            semantic_description="fixture",
            producer="test",
            claim_ceiling="fixture",
            finality="FINAL_FOR_SCOPE",
            sha256="0" * 64,
        )
        frozen = MODULE.artifact_row(
            artifact_id="authority.frozen_binary",
            artifact_type="binary",
            semantic_description="fixture",
            producer="test",
            claim_ceiling="fixture",
            finality="FINAL_FOR_SCOPE",
            sha256="1" * 64,
        )
        errors = MODULE.validate_registries(datasets, run_registry, [artifact, source4, frozen], machine, storage, crosswalk)
        self.assertTrue(any("lacks hash" in error for error in errors))

    def test_source_snapshot_4_cannot_feed_frozen_binary(self):
        source4 = MODULE.artifact_row(
            artifact_id="authority.source_snapshot_4",
            artifact_type="source_snapshot",
            semantic_description="fixture",
            producer="test",
            claim_ceiling="fixture",
            finality="FINAL_FOR_SCOPE",
            sha256="0" * 64,
            used_by=["authority.frozen_binary"],
        )
        frozen = MODULE.artifact_row(
            artifact_id="authority.frozen_binary",
            artifact_type="binary",
            semantic_description="fixture",
            producer="test",
            claim_ceiling="fixture",
            finality="FINAL_FOR_SCOPE",
            sha256="1" * 64,
            inputs=["authority.source_snapshot_4"],
        )
        datasets = [
            {"technical_dataset_id": name, "biological_id": "HCC1395" if name == "HCC1395_DORADO" else name}
            for name in MODULE.TECHNICAL_TO_BIOLOGICAL
        ]
        run_registry = {
            "reconciliation": {
                "manifest_file_lines": 19, "logical_manifest_rows": 18, "physical_directories_total": 51,
                "current_physical_directories": 35, "pending_archive_physical_directories": 16,
                "logical_rows_merged_current": 9, "logical_rows_merged_pending_archive": 9,
                "current_physical_unregistered": 26, "pending_archive_extras": 7,
            },
            "records": [{"physical_run_id": f"p{i}", "logical_manifest_member": i < 18} for i in range(51)],
        }
        machine = {"records": [{"path_id": "legacy.runbook", "observation_status": "MISSING"}]}
        errors = MODULE.validate_registries(
            datasets, run_registry, [source4, frozen], machine, {"records": []}, {"crosswalks": []}
        )
        self.assertTrue(any("source_snapshot_4" in error for error in errors))

    def test_large_bam_inventory_uses_stat_without_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "HCC1395/paired_full/run/longphase_s/HCC1395_tagged.bam"
            bam.parent.mkdir(parents=True)
            bam.write_bytes(b"not-a-real-bam")
            paths = MODULE.tagged_bam_paths(root)
            self.assertEqual(paths, [bam])
            self.assertEqual(paths[0].stat().st_size, len(b"not-a-real-bam"))

    def test_storage_manifest_is_immediate_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "child/deep").mkdir(parents=True)
            (root / "child/deep/payload").write_bytes(b"payload")
            machine = {
                "records": [
                    {
                        "path_id": "storage.big7_output",
                        "path": str(root),
                        "observation_status": "PRESENT",
                    }
                ]
            }
            manifest = MODULE.build_storage_manifest(machine)
            row = manifest["records"][0]
            self.assertEqual(row["immediate_child_count"], 1)
            self.assertIsNone(row["recursive_size_bytes"])
            self.assertEqual(row["size_measurement_status"], "NOT_ATTEMPTED_TIB_SCALE_SAFETY")

    def test_bounded_tree_size_measures_only_small_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "a").write_bytes(b"123")
            (root / "nested").mkdir()
            (root / "nested/b").write_bytes(b"4567")
            size, status, blocker = MODULE.bounded_tree_size(root, max_entries=10)
            self.assertEqual((size, status, blocker), (7, "MEASURED", None))

    def test_bounded_tree_size_fails_closed_on_entry_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for index in range(3):
                (root / f"{index}.txt").write_text("x", encoding="utf-8")
            size, status, blocker = MODULE.bounded_tree_size(root, max_entries=2)
            self.assertIsNone(size)
            self.assertEqual(status, "BLOCKED_SAFETY_BOUND")
            self.assertIn("safe bound", blocker)

    def test_schema_validation_rejects_invalid_document_when_available(self):
        # Invalid input or an unavailable validator must both fail closed.
        errors = MODULE.load_and_validate_json_schema({}, {"type": "array"})
        self.assertTrue(errors)


if __name__ == "__main__":
    unittest.main()
