#!/usr/bin/env python3
"""Contract tests for the bounded multi-BAM input manifest."""
from __future__ import annotations

import copy
import importlib.util
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
BUILDER_PATH = SCRIPT_DIR / "build_multi_bam_input_manifest.py"
DASHBOARD_BUILDER_PATH = SCRIPT_DIR / "build_multi_bam_dashboard_artifact.py"
MANIFEST_PATH = PROJECT_ROOT / "multi_bam_input_manifest.json"
SCHEMA_PATH = PROJECT_ROOT / "multi_bam_input_manifest.schema.json"
RECEIPT_PATH = PROJECT_ROOT / "results" / "multi_bam_input_manifest_validation.json"
SOURCE_MANIFEST_PATH = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json"
)


def load_builder():
    spec = importlib.util.spec_from_file_location("multi_bam_input_builder", BUILDER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {BUILDER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_dashboard_builder():
    spec = importlib.util.spec_from_file_location(
        "multi_bam_dashboard_builder", DASHBOARD_BUILDER_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {DASHBOARD_BUILDER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BUILDER = load_builder()
DASHBOARD_BUILDER = load_dashboard_builder()


class MultiBamInputManifestContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.valid = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        cls.schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))

    def mutated(self) -> dict:
        return copy.deepcopy(self.valid)

    def assert_rejected(self, manifest: dict, pattern: str) -> None:
        with self.assertRaisesRegex(BUILDER.ContractError, pattern):
            BUILDER.validate_output(manifest)

    def test_current_manifest_is_accepted(self) -> None:
        BUILDER.validate_output(self.valid)

    def test_canonical_source_schema_binding_is_accepted(self) -> None:
        schema, observed_sha = BUILDER.strict_json_load_with_sha256(
            BUILDER.DEFAULT_SOURCE_SCHEMA
        )
        BUILDER.validate_canonical_source_schema(
            schema,
            observed_sha,
            path=BUILDER.DEFAULT_SOURCE_SCHEMA,
        )

    def test_permissive_source_schema_is_rejected(self) -> None:
        permissive = {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            "$id": BUILDER.CANONICAL_SOURCE_SCHEMA_ID,
            "type": "object",
        }
        with self.assertRaisesRegex(BUILDER.ContractError, "canonical source schema SHA256"):
            BUILDER.validate_canonical_source_schema(permissive, "0" * 64)

    def test_dataset_set_is_exact(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["dataset_id"] = "NOT_A_DATASET"
        self.assert_rejected(manifest, "dataset IDs mismatch")

    def test_technical_replicate_identity_is_fixed(self) -> None:
        manifest = self.mutated()
        dorado = next(
            row for row in manifest["datasets"] if row["dataset_id"] == "HCC1395_DORADO"
        )
        dorado["technical_replicate"] = False
        self.assert_rejected(manifest, "expected one technical replicate")

    def test_tumor_payload_mismatch_is_rejected(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["bam"]["bounded_payload_status"] = "MISMATCH"
        self.assert_rejected(manifest, "sampled tumor BAM payload mismatch")

    def test_normal_quickcheck_failure_is_rejected(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["paired_normal_bam"]["quickcheck_status"] = "FAIL"
        self.assert_rejected(manifest, "normal samtools quickcheck failed")

    def test_reference_dictionary_failure_is_rejected(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["input_compatibility"][
            "tumor_reference_dictionary_status"
        ] = "FAIL"
        self.assert_rejected(manifest, "const mismatch.*tumor_reference_dictionary_status")

    def test_strict_identity_can_regenerate_as_full_match(self) -> None:
        manifest = self.mutated()
        for row in manifest["datasets"]:
            for role in ("bam", "paired_normal_bam", "reference"):
                row[role]["strict_storage_identity_status"] = "MATCH"
                row[role]["strict_differing_fields"] = []
        summary = manifest["verification_summary"]
        summary["status"] = "PASS_BOUNDED"
        summary["tumor_strict_storage_identity_match_n"] = 7
        summary["tumor_mount_device_drift_only_n"] = 0
        summary["all_input_roles_strict_storage_identity_match_n"] = 7
        summary["all_input_roles_mount_device_drift_only_n"] = 0
        BUILDER.validate_output(manifest)

    def test_dashboard_accepts_full_match_manifest_dynamically(self) -> None:
        manifest = self.mutated()
        for row in manifest["datasets"]:
            for role in ("bam", "paired_normal_bam", "reference"):
                row[role]["strict_storage_identity_status"] = "MATCH"
                row[role]["strict_differing_fields"] = []
        summary = manifest["verification_summary"]
        summary["status"] = "PASS_BOUNDED"
        summary["tumor_strict_storage_identity_match_n"] = 7
        summary["tumor_mount_device_drift_only_n"] = 0
        summary["all_input_roles_strict_storage_identity_match_n"] = 7
        summary["all_input_roles_mount_device_drift_only_n"] = 0
        BUILDER.validate_output(manifest)
        receipt = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
        receipt["status"] = summary["status"]
        receipt["verification_summary"] = copy.deepcopy(summary)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = root / "manifest.json"
            receipt_path = root / "receipt.json"
            manifest_path.write_text(
                json.dumps(manifest, ensure_ascii=False, sort_keys=True),
                encoding="utf-8",
            )
            receipt["output"] = {
                "path": str(manifest_path),
                "sha256": DASHBOARD_BUILDER.sha256(manifest_path),
            }
            receipt_path.write_text(
                json.dumps(receipt, ensure_ascii=False, sort_keys=True),
                encoding="utf-8",
            )
            artifact = DASHBOARD_BUILDER.build_artifact(
                PROJECT_ROOT / "results",
                manifest_path,
                receipt_path,
            )
        all_scope = next(
            row
            for row in artifact["snapshot"]["datasets"]["bam_scope_summary"]
            if row["sample_filter"] == "All"
        )
        self.assertEqual(all_scope["strict_storage_identity_match_n"], 7)
        self.assertEqual(all_scope["mount_device_drift_n"], 0)
        availability = next(
            row
            for row in artifact["snapshot"]["datasets"]["availability"]
            if row["sample_filter"] == "All"
            and row["component"] == "Tumor + normal bounded payload identity"
        )
        self.assertEqual(availability["status"], "AVAILABLE_BOUNDED")
        serialized = json.dumps(artifact, ensure_ascii=False)
        self.assertIn("BAM strict identity 已達 7/7", serialized)
        for stale_text in (
            "strict storage identity remains metadata-drifted",
            "strict storage identity match remains zero",
            "strict match 因 st_dev 掛載漂移保持 0",
            "目前僅 <code>st_dev</code>",
        ):
            self.assertNotIn(stale_text, serialized)

    def test_dashboard_aggregates_mixed_match_and_mount_drift(self) -> None:
        manifest = self.mutated()
        by_sample = {row["dataset_id"]: row for row in manifest["datasets"]}
        for role in ("bam", "paired_normal_bam", "reference"):
            by_sample["HCC1395"][role]["strict_storage_identity_status"] = "MATCH"
            by_sample["HCC1395"][role]["strict_differing_fields"] = []
        by_sample["H1437"]["bam"]["strict_storage_identity_status"] = "MATCH"
        by_sample["H1437"]["bam"]["strict_differing_fields"] = []
        summary = manifest["verification_summary"]
        summary["tumor_strict_storage_identity_match_n"] = 2
        summary["tumor_mount_device_drift_only_n"] = 5
        summary["all_input_roles_strict_storage_identity_match_n"] = 1
        summary["all_input_roles_mount_device_drift_only_n"] = 6
        BUILDER.validate_output(manifest)
        receipt = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
        receipt["verification_summary"] = copy.deepcopy(summary)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = root / "manifest.json"
            receipt_path = root / "receipt.json"
            manifest_path.write_text(
                json.dumps(manifest, ensure_ascii=False, sort_keys=True),
                encoding="utf-8",
            )
            receipt["output"] = {
                "path": str(manifest_path),
                "sha256": DASHBOARD_BUILDER.sha256(manifest_path),
            }
            receipt_path.write_text(
                json.dumps(receipt, ensure_ascii=False, sort_keys=True),
                encoding="utf-8",
            )
            artifact = DASHBOARD_BUILDER.build_artifact(
                PROJECT_ROOT / "results",
                manifest_path,
                receipt_path,
            )
        scopes = {
            row["sample_filter"]: row
            for row in artifact["snapshot"]["datasets"]["bam_scope_summary"]
        }
        self.assertEqual(
            (scopes["All"]["strict_storage_identity_match_n"], scopes["All"]["mount_device_drift_n"]),
            (1, 6),
        )
        self.assertEqual(
            (scopes["HCC1395"]["strict_storage_identity_match_n"], scopes["HCC1395"]["mount_device_drift_n"]),
            (1, 0),
        )
        self.assertEqual(
            (scopes["H1437"]["strict_storage_identity_match_n"], scopes["H1437"]["mount_device_drift_n"]),
            (0, 1),
        )

    def test_dashboard_rejects_incoherent_receipt_summary(self) -> None:
        receipt = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
        receipt["verification_summary"] = {"status": "FABRICATED", "dataset_count": 999}
        with tempfile.TemporaryDirectory() as directory:
            receipt_path = Path(directory) / "receipt.json"
            receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaisesRegex(
                DASHBOARD_BUILDER.ContractError,
                "validation receipt mismatch",
            ):
                DASHBOARD_BUILDER.build_artifact(
                    PROJECT_ROOT / "results",
                    MANIFEST_PATH,
                    receipt_path,
                )

    def test_full_bam_sha_cannot_be_claimed_without_scan(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["bam"]["full_content_sha256"] = "0" * 64
        self.assert_rejected(manifest, "full BAM SHA must remain null")

    def test_hp_count_cannot_exceed_all_alignments(self) -> None:
        manifest = self.mutated()
        tags = manifest["datasets"][0]["producer_read_tags"]
        tags["hp_assigned_alignment_records"] = tags["total_alignment_records"] + 1
        self.assert_rejected(manifest, "tagged alignments exceed total")

    def test_hp_ps_count_cannot_exceed_hp_count(self) -> None:
        manifest = self.mutated()
        tags = manifest["datasets"][0]["producer_read_tags"]
        tags["hp_and_ps_alignment_records"] = tags["hp_assigned_alignment_records"] + 1
        self.assert_rejected(manifest, "HP\+PS alignments exceed HP-assigned")

    def test_tag_rate_must_reconcile_with_exact_counts(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["producer_read_tags"][
            "hp_assigned_rate_all_alignment_records"
        ] = 0.999
        self.assert_rejected(manifest, "hp_assigned_rate_all_alignment_records does not reconcile")

    def test_duplicate_rate_must_reconcile_with_exact_counts(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["producer_read_tags"][
            "duplicate_identity_rate_all_alignment_records"
        ] = 0.999
        self.assert_rejected(manifest, "duplicate_identity_rate_all_alignment_records")

    def test_hcc1395_cross_directory_pairing_is_not_silenced(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["input_compatibility"][
            "source_directory_pairing_status"
        ] = "SAME_DIRECTORY_FAMILY"
        self.assert_rejected(manifest, "source-directory pairing changed")

    def test_missing_read_group_cannot_be_reported_as_identity_pass(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["input_compatibility"][
            "read_group_identity_status"
        ] = "AVAILABLE_UNVERIFIED"
        self.assert_rejected(manifest, "read-group identity status does not match BAM headers")

    def test_verification_summary_is_fixed_to_observed_scope(self) -> None:
        manifest = self.mutated()
        manifest["verification_summary"]["normal_quickcheck_pass_n"] = 6
        self.assert_rejected(manifest, "verification_summary.normal_quickcheck_pass_n")

    def test_schema_const_claim_ceiling_is_enforced(self) -> None:
        manifest = self.mutated()
        manifest["claim_ceiling"] = "scientifically_validated"
        self.assert_rejected(manifest, "const mismatch.*claim_ceiling")

    def test_schema_const_task_type_is_enforced(self) -> None:
        manifest = self.mutated()
        manifest["task_type"] = "production_deployment"
        self.assert_rejected(manifest, "const mismatch.*task_type")

    def test_schema_authority_and_assurance_are_enforced(self) -> None:
        for field, invalid in (("authority", "fabricated"), ("assurance", "global_full_bam")):
            with self.subTest(field=field):
                manifest = self.mutated()
                manifest["datasets"][0]["producer_read_tags"][field] = invalid
                self.assert_rejected(manifest, f"const mismatch.*{field}")

    def test_schema_denominator_is_enforced(self) -> None:
        manifest = self.mutated()
        manifest["datasets"][0]["producer_read_tags"][
            "denominator_population"
        ] = "primary_reads"
        self.assert_rejected(manifest, "const mismatch.*denominator_population")

    def test_schema_required_fields_are_enforced(self) -> None:
        for field in ("source_manifest", "generated_at_utc"):
            with self.subTest(field=field):
                manifest = self.mutated()
                del manifest[field]
                self.assert_rejected(manifest, "required field")

    def test_schema_additional_properties_are_rejected(self) -> None:
        manifest = self.mutated()
        manifest["unexpected_field"] = True
        self.assert_rejected(manifest, "additionalProperties")

    def synthetic_frozen_identity(self) -> tuple[dict, Path, Path]:
        payload = Path("/bounded/example.bam")
        index = Path("/bounded/example.bam.bai")
        size = BUILDER.FROZEN_CHUNK_SIZE_BYTES * 4
        chunk_size = BUILDER.FROZEN_CHUNK_SIZE_BYTES
        frozen = {
            "policy": "storage_identity_v1",
            "assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
            "is_full_content_hash": False,
            "requested_path": str(payload),
            "realpath": str(payload),
            "logical_is_symlink": False,
            "logical_size_bytes": size,
            "logical_mtime_ns": 1,
            "st_dev": 1,
            "st_ino": 1,
            "size_bytes": size,
            "mtime_ns": 1,
            "ctime_ns": 1,
            "chunk_size_bytes": chunk_size,
            "chunks": [
                {"label": "first", "offset": 0, "length": chunk_size, "sha256": "1" * 64},
                {
                    "label": "middle",
                    "offset": (size - chunk_size) // 2,
                    "length": chunk_size,
                    "sha256": "2" * 64,
                },
                {
                    "label": "last",
                    "offset": size - chunk_size,
                    "length": chunk_size,
                    "sha256": "3" * 64,
                },
            ],
            "index": {
                "path": str(index),
                "identity": {"policy": "full_sha256", "size_bytes": 1, "sha256": "4" * 64},
            },
        }
        frozen["identity_sha256"] = BUILDER.canonical_sha256(frozen)
        return frozen, payload, index

    def test_chunk_policy_accepts_exact_fixed_layout(self) -> None:
        frozen, payload, index = self.synthetic_frozen_identity()
        BUILDER.validate_frozen_storage_identity(
            frozen, label="test", payload_path=payload, index_path=index
        )

    def test_chunk_policy_rejects_two_mib_and_arbitrary_label(self) -> None:
        for mutation, pattern in (("length", "exactly"), ("label", "must be 'first'")):
            with self.subTest(mutation=mutation):
                frozen, payload, index = self.synthetic_frozen_identity()
                if mutation == "length":
                    frozen["chunks"][0]["length"] = 2 * BUILDER.FROZEN_CHUNK_SIZE_BYTES
                else:
                    frozen["chunks"][0]["label"] = "arbitrary"
                frozen["identity_sha256"] = BUILDER.canonical_sha256(
                    {key: value for key, value in frozen.items() if key != "identity_sha256"}
                )
                with self.assertRaisesRegex(BUILDER.ContractError, pattern):
                    BUILDER.validate_frozen_storage_identity(
                        frozen, label="test", payload_path=payload, index_path=index
                    )

    def test_io_output_and_receipt_collision_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = Namespace(
                source_manifest=root / "source.json",
                topology_csv=root / "topology.csv",
                source_schema=root / "source.schema.json",
                schema=root / "output.schema.json",
                output=root / "same.json",
                receipt=root / "same.json",
            )
            with self.assertRaisesRegex(BUILDER.ContractError, "different paths"):
                BUILDER.validate_io_path_separation(args)

    def test_io_output_cannot_overwrite_an_input(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            args = Namespace(
                source_manifest=source,
                topology_csv=root / "topology.csv",
                source_schema=root / "source.schema.json",
                schema=root / "output.schema.json",
                output=source,
                receipt=root / "receipt.json",
            )
            with self.assertRaisesRegex(BUILDER.ContractError, "collides with input"):
                BUILDER.validate_io_path_separation(args)

    def test_dashboard_output_cannot_overwrite_implicit_input(self) -> None:
        args = Namespace(
            results_dir=PROJECT_ROOT / "results",
            bam_manifest=MANIFEST_PATH,
            bam_manifest_receipt=RECEIPT_PATH,
            output=DASHBOARD_BUILDER.DEFAULT_BAM_MANIFEST_SCHEMA,
        )
        with self.assertRaisesRegex(
            DASHBOARD_BUILDER.ContractError,
            "dashboard output collides with input",
        ):
            DASHBOARD_BUILDER.validate_io_path_separation(args)

    @unittest.skipUnless(SOURCE_MANIFEST_PATH.is_file(), "canonical external source unavailable")
    def test_cross_sample_read_tag_receipt_substitution_is_rejected(self) -> None:
        source = json.loads(SOURCE_MANIFEST_PATH.read_text(encoding="utf-8"))
        h1437 = next(row for row in source["samples"] if row["sample"] == "H1437")
        with self.assertRaisesRegex(BUILDER.ContractError, "sample mismatch"):
            BUILDER.producer_inputs("HCC1395", h1437["read_tags"])


if __name__ == "__main__":
    unittest.main()
