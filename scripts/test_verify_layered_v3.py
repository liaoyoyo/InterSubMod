#!/usr/bin/env python3
"""Cross-component tiny tests for the layered-v3 verifier.

Every fixture uses the real contract builder to create the frozen lock and the
real ``RunLifecycle`` to create/publish source/environment/receipt/state
artifacts through VERIFYING.  All biological/scientific files remain tiny and
inside ``TemporaryDirectory``; no production path is opened.
"""

from __future__ import annotations

import gzip
import fcntl
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
VERIFIER_PATH = HERE / "verify_layered_v3.py"
DOC_SCRIPTS = (
    REPO / "docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts"
).resolve()
sys.path.insert(0, str(DOC_SCRIPTS))
sys.path.insert(0, str(HERE))

import layered_v3_lifecycle as lifecycle_contract  # noqa: E402
import sm_multilocus_combinations as combinations_contract  # noqa: E402
import test_layered_v3_contract as contract_fixture  # noqa: E402
import validate_layered_v3_inputs as input_contract  # noqa: E402


SPEC = importlib.util.spec_from_file_location("verify_layered_v3", VERIFIER_PATH)
assert SPEC and SPEC.loader
V3 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = V3
SPEC.loader.exec_module(V3)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def file_sha(path: Path) -> str:
    return lifecycle_contract.sha256_file(path, reject_symlink=True)


class Fixture:
    def __init__(self, parent: Path, *, sensitivity: bool = False, parallel_samples: int = 2):
        self.parent = parent
        self.sensitivity = sensitivity
        self.parallel_samples = parallel_samples
        self.run_parent = parent / "runs"
        self.run_parent.mkdir()
        self.contract_root = parent / "contract"
        self.manifest, _ = contract_fixture.make_fixture(self.contract_root)
        measured_bed = self.contract_root / "measured_h1437_cn.bed"
        measured_bed.write_text("chr1\t0\t100\tneutral\n", encoding="utf-8")
        h1437 = next(item for item in self.manifest["samples"] if item["sample"] == "H1437")
        h1437["copy_number"] = {
            "availability": "measured",
            "source": "reviewed_cross_component_fixture",
            "semantics": "explicit segments; unlisted positions are neutral",
            "coordinate_system": "0_based_half_open",
            "unlisted_position_semantics": "neutral",
            "allowed_states": ["gain", "loss", "loh", "neutral"],
            "overlap_policy": "forbid",
            "reason": None,
            "cn_bed": input_contract.artifact(measured_bed),
            "cn_int_gain": None,
            "cn_int_loss": None,
            "integration_json": None,
        }
        if self.sensitivity:
            self.manifest["analysis_contract"]["task_type"] = "backbone_sensitivity"
            self.manifest["analysis_contract"]["tree_input_contract"] = V3.SENSITIVITY_TREE_INPUT_CONTRACT
            for item in self.manifest["samples"]:
                item["somatic"]["backbone_role"] = "clairs_filter_pass_sensitivity"
                item["somatic"]["tree_vcf"] = json.loads(
                    json.dumps(item["somatic"]["caller_pass_baseline_vcf"])
                )
        self.manifest_path = self.contract_root / "layered_input_manifest_v3.json"
        contract_fixture.write_json(self.manifest_path, self.manifest)
        self.lock = input_contract.validate_and_build_lock(
            self.manifest_path,
            contract_fixture.SCHEMA_PATH,
            validator_path=Path(input_contract.__file__).resolve(),
        )
        self.lifecycle = lifecycle_contract.RunLifecycle(
            self.run_parent, "synthetic_layered_v3_run", install_traps=False
        )
        self._build_lifecycle()
        self.sample_metadata = {item["sample"]: item for item in self.lock["samples"]}
        self._build_outputs()
        self.lifecycle.begin_verifying()
        self.run_root = self.lifecycle.root
        self.lock_path = self.run_root / V3.RUN_LOCK_FILENAME
        self.receipt_path = self.run_root / V3.RUN_RECEIPT_FILENAME
        self.state_path = self.run_root / V3.RUN_STATE_FILENAME

    def close(self) -> None:
        if self.lifecycle.state not in lifecycle_contract.TERMINAL_STATES:
            self.lifecycle.fail("E_TEST_END", "synthetic verifier test ended", grace_seconds=0.01)
        self.lifecycle.close()
        if getattr(self, "global_lock_fd", None) is not None:
            fcntl.flock(self.global_lock_fd, fcntl.LOCK_UN)
            os.close(self.global_lock_fd)
            self.global_lock_fd = None

    def _build_lifecycle(self) -> None:
        sources = self.parent / "lifecycle_sources"
        sources.mkdir()
        runner = sources / "runner.py"
        consumer = sources / "consumer.py"
        runner.write_text("# synthetic runner\n", encoding="utf-8")
        consumer.write_text("# synthetic consumer\n", encoding="utf-8")
        runner.chmod(0o750)
        consumer.chmod(0o640)
        self.lifecycle.begin_preflight()
        global_lock = self.run_parent / ".layered_chr1_22_7dataset_full.lock"
        self.global_lock_fd = os.open(global_lock, os.O_RDWR | os.O_CREAT, 0o600)
        fcntl.flock(self.global_lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        os.ftruncate(self.global_lock_fd, 0)
        os.write(
            self.global_lock_fd,
            (
                f"run_id={self.lifecycle.run_id} pid={os.getpid()} "
                f"pid_start_time={self.lifecycle.launcher_pid_start_time}\n"
            ).encode(),
        )
        os.fsync(self.global_lock_fd)
        process_path = self.lifecycle.root / "process_observation.json"
        lifecycle_contract.atomic_write_json(
            process_path,
            {
                "schema_name": "intersubmod.layered_process_observation",
                "schema_version": "1.0.0",
                "created_at_utc": lifecycle_contract.utc_now(),
                "observer_pid": os.getpid(),
                "observer_pid_start_time": self.lifecycle.launcher_pid_start_time,
                "global_scope_lock": {
                    "path": str(global_lock),
                    "held": True,
                    "run_id": self.lifecycle.run_id,
                },
                "match_policy": ["run_layered_v3.py", "run_layered_7samples_newbb.sh"],
                "processes_inspected": 0,
                "conflicts": [],
                "pass": True,
            },
            no_overwrite=True,
        )
        self.lifecycle.write_frozen_lock(self.lock)
        self.lifecycle.build_source_bundle(
            runner,
            Path(input_contract.__file__).resolve(),
            VERIFIER_PATH,
            [consumer],
        )
        self.lifecycle.capture_environment_lock(
            sm_allowlist=set(),
            environment={"LC_ALL": "C", "TZ": "UTC", "PYTHONHASHSEED": "0"},
            distributions=(),
            tools=(),
            storage_path=self.run_parent,
        )
        lifecycle_contract.atomic_write_bytes(
            self.lifecycle.root / "input_manifest.snapshot.json",
            self.manifest_path.read_bytes(),
            mode=0o444,
            no_overwrite=True,
        )
        self.lifecycle.publish_ready(
            {
                "mode": "full",
                "task_type": "backbone_sensitivity" if self.sensitivity else "comprehensive_validation",
                "analysis_params": {
                    "scope": "chr1-22",
                    "contigs": list(V3.AUTOSOMES),
                    "dataset_count": 7,
                    "biological_sample_count": 6,
                    "parallel_samples": self.parallel_samples,
                    "parallel_parts_per_sample": 1,
                    "VERIFY_EVERY": 1,
                    "ANALYSIS_TREE_CAP": 0,
                    "DISPLAY_TREE_CAP": 32,
                    "MINREAD": 3,
                    "MAX_SNV": 8,
                    "TIER_R": 50000,
                    "MAPQ_MIN": 20,
                    "BASEQ_MIN": 0,
                },
                "source_manifest_snapshot": {
                    "path": "input_manifest.snapshot.json",
                    "sha256": file_sha(self.manifest_path),
                },
                "process_observation": {
                    "path": "process_observation.json",
                    "sha256": file_sha(process_path),
                    "conflict_count": 0,
                },
                "validator_execution": "bundled_core_via_runner_adapter",
                "worker_input_authority": "frozen_input_lock.json",
                "tree_backbone_role": (
                    "clairs_filter_pass_sensitivity"
                    if self.sensitivity
                    else "longphase_s_recalibrated_filter_pass"
                ),
                "ledger_roles": [
                    "caller_raw_vcf",
                    "longphase_input_vcf",
                    "caller_pass_baseline_vcf",
                    "longphase_recalibrated_all_vcf",
                    "longphase_recalibrated_pass_vcf",
                    "tree_vcf",
                ],
            }
        )
        self.lifecycle.begin_running()

    def provenance(self, metadata: dict[str, Any]) -> dict[str, str]:
        return {
            "frozen_lock_sha256": file_sha(self.lifecycle.root / V3.RUN_LOCK_FILENAME),
            "launch_receipt_sha256": file_sha(self.lifecycle.root / V3.RUN_RECEIPT_FILENAME),
            "environment_lock_sha256": self.lifecycle.environment_lock_sha256,
            "source_bundle_manifest_sha256": self.lifecycle.source_bundle_manifest_sha256,
            "source_bundle_content_sha256": self.lifecycle.source_bundle_content_sha256,
            "input_set_sha256": V3.input_set_digest(metadata),
        }

    def somatic_roles(self, metadata: dict[str, Any]) -> dict[str, str]:
        somatic = metadata["somatic"]
        return {
            "longphase_input_role": "normalized_clairs_raw_all",
            "longphase_input_vcf_sha256": somatic["longphase_input_vcf"]["identity"]["sha256"],
            "caller_pass_baseline_role": (
                "clairs_filter_pass_selected_sensitivity_tree"
                if self.sensitivity
                else "clairs_filter_pass_sensitivity_only"
            ),
            "caller_pass_baseline_vcf_sha256": somatic["caller_pass_baseline_vcf"]["identity"]["sha256"],
            "tree_input_contract": (
                V3.SENSITIVITY_TREE_INPUT_CONTRACT
                if self.sensitivity
                else V3.CANONICAL_TREE_INPUT_CONTRACT
            ),
            "tree_backbone_role": (
                "clairs_filter_pass_sensitivity"
                if self.sensitivity
                else "longphase_s_recalibrated_filter_pass"
            ),
            "tree_vcf_sha256": somatic["tree_vcf"]["identity"]["sha256"],
            "ledger_role": "clairs_raw",
            "caller_raw_vcf_sha256": somatic["caller_raw_vcf"]["identity"]["sha256"],
            "longphase_recalibrated_all_vcf_sha256": somatic[
                "longphase_recalibrated_all_vcf"
            ]["identity"]["sha256"],
            "longphase_recalibrated_pass_vcf_sha256": somatic[
                "longphase_recalibrated_pass_vcf"
            ]["identity"]["sha256"],
        }

    @staticmethod
    def part_funnel() -> dict[str, Any]:
        # Deliberately non-zero unsupported branch: the valid formula is
        # 5 = singleton(1) + cap(0) + read_unsupported(2) + retained(2).
        return {
            "n_sSNV_scope_input": 5,
            "n_positional_singleton": 1,
            "n_multilocus_pre_cap_groups": 2,
            "n_multilocus_pre_cap_sSNV": 4,
            "n_groups_retained": 1,
            "n_groups_read_unsupported": 1,
            "n_sSNV_retained": 2,
            "n_sSNV_read_unsupported": 2,
            "n_groups_capped_by_MAX_SNV": 0,
            "n_sSNV_cap_excluded": 0,
            "n_sSNV_accounted": 5,
            "check_scope_conservation": True,
        }

    @staticmethod
    def aggregate_funnel() -> dict[str, Any]:
        result = {
            key: value * 5 if isinstance(value, int) and not isinstance(value, bool) else value
            for key, value in Fixture.part_funnel().items()
        }
        result["check_scope_conservation"] = True
        return result

    def _build_outputs(self) -> None:
        (self.lifecycle.root / "samples").mkdir()
        (self.lifecycle.root / "verification").mkdir()
        for sample in sorted(V3.EXPECTED_DATASETS):
            metadata = self.sample_metadata[sample]
            sample_root = self.lifecycle.root / "samples" / sample
            sample_root.mkdir()
            provenance = self.provenance(metadata)
            roles = self.somatic_roles(metadata)
            groups = []
            details = []
            regions = []
            tag_subject = metadata["read_tags"]["subject_binding"]
            for part in range(1, 6):
                chrom = V3.PART_CHROMOSOMES[part][0]
                start, end = part * 100, part * 100 + 10
                group = {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "positions": [start, end],
                    "n_sSNV": 2,
                }
                groups.append(group)
                write_json(
                    sample_root / V3.output_name(sample, f"mlhp_part_{part}"),
                    {
                        "schema_version": "2.0",
                        "sample": sample,
                        "part": part,
                        "chromosomes": list(V3.PART_CHROMOSOMES[part]),
                        "provenance": provenance,
                        "somatic_roles": roles,
                        "input_funnel": self.part_funnel(),
                        "groups": [group],
                        "read_tag_census": {
                            "identity_schema": "coordinate_join_v1",
                            "sidecar_sha256": tag_subject["sidecar_sha256"],
                            "sidecar_index_sha256": tag_subject["sidecar_index_sha256"],
                            "alignment_payload_identity_sha256": tag_subject[
                                "alignment_payload_storage_identity_sha256"
                            ],
                            "sidecar_missing": 0,
                            "sidecar_conflicts": 0,
                            "alignment_identity_allele_conflicts": 0,
                            "sidecar_duplicates": tag_subject["duplicate_count"],
                            "duplicate_identity_policy": tag_subject[
                                "duplicate_identity_policy"
                            ],
                            "sidecar_extra": 0,
                            "sidecar_malformed": 0,
                            "sidecar_exact_matches": 2,
                            "alignment_group_exposures": 2,
                        },
                    },
                )
                details.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "family": "1",
                        "unit_role": "primary_mutation_lineage",
                        "is_primary_lineage": True,
                        "reference_only": False,
                        "capped": False,
                        "verification_skipped": False,
                        "verification_status": "full_pass",
                        "analysis_candidate_set_complete": True,
                        "analysis_trees_generated": 1,
                        "n_trees": 1,
                        "n_distinct_shapes_exact": 1,
                        "trees": [{"n_hidden": 0}],
                    }
                )
                regions.append(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "hp_multiplicity": 1,
                        "lineages": [
                            {"family": "1", "is_primary_lineage": True, "reference_only": False}
                        ],
                    }
                )
            l1 = {
                "n_units_total_including_unphased": 5,
                "n_primary_lineage_units": 5,
                "n_reference_only_controls": 0,
                "n_verification_fail": 0,
                "n_eligible_skipped_V4V5": 0,
                "all_eligible_V1V7_pass": True,
            }
            write_json(
                sample_root / V3.output_name(sample, "layered_reconstruction"),
                {
                    "schema_version": "2.0",
                    "sample": sample,
                    "provenance": provenance,
                    "input_funnel": self.aggregate_funnel(),
                    "n_detail_units": 5,
                    "detail": details,
                    "L1_ssnv_algorithm": l1,
                },
            )
            write_json(
                sample_root / V3.output_name(sample, "layered_region_view"),
                {
                    "schema_version": "2.0",
                    "sample": sample,
                    "provenance": provenance,
                    "copy_number_contract": metadata["copy_number"],
                    "census": {"n_regions": 5, "funnel": self.aggregate_funnel(), "L1": l1},
                    "regions": regions,
                },
            )
            ledger = sample_root / V3.output_name(sample, "site_ledger")
            with gzip.open(ledger, "wt", encoding="utf-8", newline="") as handle:
                handle.write("sample\tchrom\tpos\n")
                for number in range(25):
                    handle.write(f"{sample}\tchr1\t{number + 1}\n")
            write_json(
                sample_root / V3.output_name(sample, "site_ledger_summary"),
                {
                    "schema_version": "2.0",
                    "sample": sample,
                    "provenance": provenance,
                    "longphase_input_contract": "clairs_raw_all",
                    "tree_contract": (
                        "clairs_PASS_input" if self.sensitivity else "longphase_recalibrated_PASS"
                    ),
                    "raw_clairs_records": 25,
                    "longphase_input_records": 25,
                    "longphase_recalibrated_records": 25,
                    "filter_transition_counts": {"LowQual->PASS": 5, "PASS->PASS": 20},
                    "duplicate_record_key_excess": {
                        "raw_clairs": 0,
                        "longphase_input": 0,
                        "longphase_recalibrated": 0,
                        "tree_input": 0,
                    },
                    "branch_counts": {
                        "positional_singleton": 5,
                        "read_unsupported": 10,
                        "retained": 10,
                    },
                    "checks": {
                        "all_raw_records_written": True,
                        "autosomal_snv_conservation": True,
                        "all_mlhp_retained_groups_joined": True,
                    },
                    "pass": True,
                },
            )
            outputs = []
            for role in V3.SCIENTIFIC_ROLES:
                path = sample_root / V3.output_name(sample, role)
                outputs.append(
                    {
                        "role": role,
                        "path": path.name,
                        "size_bytes": path.stat().st_size,
                        "sha256": file_sha(path),
                    }
                )
            write_json(
                sample_root / "output_manifest.json",
                {
                    "schema_name": V3.OUTPUT_MANIFEST_SCHEMA_NAME,
                    "schema_version": V3.OUTPUT_MANIFEST_SCHEMA_VERSION,
                    "sample": sample,
                    "biological_id": metadata["biological_id"],
                    "run_id": self.lifecycle.run_id,
                    **provenance,
                    "somatic_roles": roles,
                    "copy_number_contract": metadata["copy_number"],
                    "outputs": outputs,
                },
            )
        manifest_index = []
        for sample in sorted(V3.EXPECTED_DATASETS):
            path = self.lifecycle.root / "samples" / sample / "output_manifest.json"
            manifest_index.append(
                {
                    "sample": sample,
                    "path": str(path.relative_to(self.lifecycle.root)),
                    "sha256": file_sha(path),
                }
            )
        write_json(
            self.lifecycle.root / "output_manifests.json",
            {"dataset_count": 7, "manifests": manifest_index},
        )

    def output(self, name: str) -> Path:
        return self.run_root / "verification" / f"{name}.json"

    def run(self, name: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(VERIFIER_PATH),
                "--run-root",
                str(self.run_root),
                "--frozen-lock",
                str(self.lock_path),
                "--launch-receipt",
                str(self.receipt_path),
                "--output",
                str(self.output(name)),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def summary(self, name: str) -> dict[str, Any]:
        return json.loads(self.output(name).read_text(encoding="utf-8"))

    def rebind_output(self, sample: str, role: str) -> None:
        sample_root = self.run_root / "samples" / sample
        manifest_path = sample_root / "output_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        target = sample_root / V3.output_name(sample, role)
        for item in manifest["outputs"]:
            if item["role"] == role:
                item["size_bytes"] = target.stat().st_size
                item["sha256"] = file_sha(target)
                break
        write_json(manifest_path, manifest)
        index_path = self.run_root / "output_manifests.json"
        index = json.loads(index_path.read_text(encoding="utf-8"))
        for item in index["manifests"]:
            if item["sample"] == sample:
                item["sha256"] = file_sha(manifest_path)
                break
        write_json(index_path, index)


class FunnelSerializationTests(unittest.TestCase):
    def test_empty_funnel_serializes_all_verifier_fields_as_zero(self) -> None:
        funnel = combinations_contract.new_funnel()
        self.assertEqual(set(funnel), set(combinations_contract.FUNNEL_FIELDS))
        self.assertTrue(all(funnel[field] == 0 for field in combinations_contract.FUNNEL_FIELDS))


class StorageIdentityPathTests(unittest.TestCase):
    def test_symlinked_bam_index_keeps_frozen_logical_path(self) -> None:
        with tempfile.TemporaryDirectory(prefix="layered-v3-storage-path-") as temporary:
            root = Path(temporary)
            run_root = root / "run"
            run_root.mkdir()
            lock_path = run_root / V3.RUN_LOCK_FILENAME
            receipt_path = run_root / V3.RUN_RECEIPT_FILENAME
            lock_path.write_text("{}\n", encoding="utf-8")
            receipt_path.write_text("{}\n", encoding="utf-8")

            target_root = root / "targets"
            logical_root = root / "logical"
            target_root.mkdir()
            logical_root.mkdir()
            target_bam = target_root / "sample.bam"
            target_index = target_root / "sample.bam.bai"
            target_bam.write_bytes(b"synthetic-bam-payload")
            target_index.write_bytes(b"synthetic-index-payload")
            logical_bam = logical_root / "sample.bam"
            logical_index = logical_root / "sample.bam.bai"
            logical_bam.symlink_to(target_bam)
            logical_index.symlink_to(target_index)

            expected = input_contract.storage_identity_v1(logical_bam, logical_index)
            verifier = V3.Verifier(run_root, lock_path, receipt_path)
            _resolved, identity_sha = verifier._verify_storage_identity(
                {"path": str(logical_bam), "storage_identity_v1": expected}, "SYNTHETIC"
            )

            self.assertEqual(identity_sha, expected["identity_sha256"])
            self.assertTrue(verifier.checks[-1].passed)
            self.assertEqual(verifier.checks[-1].observed["index"]["path"], str(logical_index))


class VerifyLayeredV3Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="layered-v3-verifier-")
        self.fixture = Fixture(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.fixture.close()
        self.temporary.cleanup()

    def assert_error(self, name: str, code: str) -> dict[str, Any]:
        result = self.fixture.run(name)
        self.assertEqual(result.returncode, 7, msg=f"stdout={result.stdout}\nstderr={result.stderr}")
        summary = self.fixture.summary(name)
        self.assertFalse(summary["all_pass"])
        self.assertIn(code, summary["error_codes"])
        self.assertFalse((self.fixture.run_root / "_SUCCESS").exists())
        return summary

    def test_lifecycle_contract_happy_path_and_runner_only_success(self) -> None:
        result = self.fixture.run("valid")
        self.assertEqual(result.returncode, 0, msg=f"stdout={result.stdout}\nstderr={result.stderr}")
        summary = self.fixture.summary("valid")
        self.assertTrue(summary["all_pass"])
        self.assertEqual(summary["dataset_count"], 7)
        self.assertEqual(summary["biological_sample_count"], 6)
        self.assertEqual(summary["n_pass"], 7)
        self.assertTrue(all(item["metrics"]["n_sSNV_scope_input"] == 25 for item in summary["samples"]))
        h1437_manifest = json.loads(
            (self.fixture.run_root / "samples/H1437/output_manifest.json").read_text()
        )
        self.assertEqual(
            h1437_manifest["copy_number_contract"]["source"],
            "reviewed_cross_component_fixture",
        )
        unavailable_manifest = json.loads(
            (self.fixture.run_root / "samples/COLO829/output_manifest.json").read_text()
        )
        self.assertEqual(
            unavailable_manifest["copy_number_contract"]["reason"],
            "No reviewed measured CN source is available",
        )
        self.assertFalse((self.fixture.run_root / "_SUCCESS").exists())
        marker = self.fixture.lifecycle.succeed(self.fixture.output("valid"), {"dataset_count": 7})
        self.assertTrue(marker.is_file())
        success = json.loads(marker.read_text(encoding="utf-8"))
        self.assertEqual(success["verification_sha256"], file_sha(self.fixture.output("valid")))

    def test_explicit_clairs_pass_sensitivity_is_independently_verified(self) -> None:
        root = Path(self.temporary.name) / "sensitivity"
        root.mkdir()
        sensitivity = Fixture(root, sensitivity=True)
        try:
            result = sensitivity.run("sensitivity")
            self.assertEqual(result.returncode, 0, msg=f"stdout={result.stdout}\nstderr={result.stderr}")
            summary = sensitivity.summary("sensitivity")
            self.assertTrue(summary["all_pass"])
            self.assertEqual(summary["n_pass"], 7)
            manifest = json.loads(
                (sensitivity.run_root / "samples/HCC1395/output_manifest.json").read_text()
            )
            roles = manifest["somatic_roles"]
            self.assertEqual(roles["tree_input_contract"], V3.SENSITIVITY_TREE_INPUT_CONTRACT)
            self.assertEqual(roles["tree_backbone_role"], "clairs_filter_pass_sensitivity")
            self.assertEqual(
                roles["tree_vcf_sha256"], roles["caller_pass_baseline_vcf_sha256"]
            )
            lock = json.loads(sensitivity.lock_path.read_text())
            somatic = next(item for item in lock["samples"] if item["sample"] == "HCC1395")["somatic"]
            self.assertNotEqual(
                somatic["tree_vcf"]["path"], somatic["longphase_recalibrated_pass_vcf"]["path"]
            )
        finally:
            sensitivity.close()

    def test_parallel_sample_profile_accepts_four_and_rejects_five(self) -> None:
        root = Path(self.temporary.name)
        four_root = root / "parallel-four"
        four_root.mkdir()
        four = Fixture(four_root, parallel_samples=4)
        try:
            result = four.run("parallel-four")
            self.assertEqual(result.returncode, 0, msg=f"stdout={result.stdout}\nstderr={result.stderr}")
            self.assertTrue(four.summary("parallel-four")["all_pass"])
        finally:
            four.close()

        five_root = root / "parallel-five"
        five_root.mkdir()
        five = Fixture(five_root, parallel_samples=5)
        try:
            result = five.run("parallel-five")
            self.assertEqual(result.returncode, 7, msg=f"stdout={result.stdout}\nstderr={result.stderr}")
            summary = five.summary("parallel-five")
            self.assertFalse(summary["all_pass"])
            self.assertIn("E_STATE_INVALID", summary["error_codes"])
        finally:
            five.close()

    def test_empty_sample_set_is_not_vacuous_pass(self) -> None:
        lock = json.loads(self.fixture.lock_path.read_text())
        lock["dataset_count"] = 0
        lock["biological_sample_count"] = 0
        lock["samples"] = []
        write_json(self.fixture.lock_path, lock)
        self.assert_error("empty", "E_DATASET_SET_MISMATCH")

    def test_duplicate_sample_is_rejected(self) -> None:
        lock = json.loads(self.fixture.lock_path.read_text())
        lock["samples"][-1] = lock["samples"][0]
        write_json(self.fixture.lock_path, lock)
        self.assert_error("duplicate", "E_SAMPLE_DUPLICATE")

    def test_extra_sample_directory_is_rejected(self) -> None:
        (self.fixture.run_root / "samples" / "EXTRA").mkdir()
        self.assert_error("extra-directory", "E_DATASET_SET_MISMATCH")

    def test_post_input_storage_identity_mutation_is_rejected(self) -> None:
        bam = Path(self.fixture.sample_metadata["COLO829"]["alignment_payload"]["path"])
        bam.write_bytes(bam.read_bytes() + b"MUTATION")
        self.assert_error("input-mutation", "E_POST_INPUT_IDENTITY")

    def test_source_bundle_mutation_is_rejected(self) -> None:
        runner = self.fixture.run_root / "source_bundle/files/core/runner.py"
        runner.write_bytes(runner.read_bytes() + b"MUTATION")
        self.assert_error("source-mutation", "E_SOURCE_BUNDLE_MISMATCH")

    def test_environment_mutation_is_rejected(self) -> None:
        environment = self.fixture.run_root / V3.ENVIRONMENT_LOCK_FILENAME
        environment.write_text(environment.read_text() + " ", encoding="utf-8")
        self.assert_error("environment-mutation", "E_ENVIRONMENT_MISMATCH")

    def test_output_digest_mutation_is_rejected(self) -> None:
        sample, role = "H1437", "layered_reconstruction"
        path = self.fixture.run_root / "samples" / sample / V3.output_name(sample, role)
        path.write_text(path.read_text() + " ", encoding="utf-8")
        self.assert_error("output-mutation", "E_OUTPUT_HASH_MISMATCH")

    def test_part_conservation_with_unsupported_branch_is_recomputed(self) -> None:
        sample, role = "H2009", "mlhp_part_1"
        path = self.fixture.run_root / "samples" / sample / V3.output_name(sample, role)
        document = json.loads(path.read_text())
        document["input_funnel"]["n_sSNV_read_unsupported"] = 0
        write_json(path, document)
        self.fixture.rebind_output(sample, role)
        self.assert_error("part-conservation", "E_PART_CONSERVATION")

    def test_layered_boolean_cannot_override_recomputation(self) -> None:
        sample, role = "HCC1937", "layered_reconstruction"
        path = self.fixture.run_root / "samples" / sample / V3.output_name(sample, role)
        document = json.loads(path.read_text())
        document["L1_ssnv_algorithm"]["all_eligible_V1V7_pass"] = False
        write_json(path, document)
        self.fixture.rebind_output(sample, role)
        self.assert_error("layered-invariant", "E_LAYERED_INVARIANT")

    def test_region_multiplicity_is_recomputed(self) -> None:
        sample, role = "HCC1954", "layered_region_view"
        path = self.fixture.run_root / "samples" / sample / V3.output_name(sample, role)
        document = json.loads(path.read_text())
        document["regions"][0]["hp_multiplicity"] = 2
        write_json(path, document)
        self.fixture.rebind_output(sample, role)
        self.assert_error("region-invariant", "E_REGION_INVARIANT")

    def test_site_ledger_rows_are_streamed_and_recomputed(self) -> None:
        sample, role = "HCC1395", "site_ledger_summary"
        path = self.fixture.run_root / "samples" / sample / V3.output_name(sample, role)
        document = json.loads(path.read_text())
        document["branch_counts"]["retained"] = 9
        write_json(path, document)
        self.fixture.rebind_output(sample, role)
        self.assert_error("ledger-invariant", "E_SITE_LEDGER_INVARIANT")

    def test_part_provenance_drift_is_rejected_after_output_rehash(self) -> None:
        sample, role = "HCC1395_DORADO", "mlhp_part_2"
        path = self.fixture.run_root / "samples" / sample / V3.output_name(sample, role)
        document = json.loads(path.read_text())
        document["provenance"]["input_set_sha256"] = "0" * 64
        write_json(path, document)
        self.fixture.rebind_output(sample, role)
        self.assert_error("part-provenance", "E_SAMPLE_OUTPUT_BINDING")

    def test_false_payload_v2_claim_is_rejected(self) -> None:
        lock = json.loads(self.fixture.lock_path.read_text())
        lock["analysis_contract"]["sidecar_identity_schema"] = "alignment_identity_v2"
        lock["analysis_contract"]["sidecar_assurance"] = "global_alignment_payload_identity"
        write_json(self.fixture.lock_path, lock)
        self.assert_error("false-v2", "E_SIDECAR_CONTRACT_UNSUPPORTED")

    def test_tree_backbone_cannot_be_swapped_with_longphase_input(self) -> None:
        sample = "HCC1395"
        sample_root = self.fixture.run_root / "samples" / sample
        manifest_path = sample_root / "output_manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["somatic_roles"]["tree_vcf_sha256"] = manifest["somatic_roles"][
            "longphase_input_vcf_sha256"
        ]
        for part in range(1, 6):
            path = sample_root / V3.output_name(sample, f"mlhp_part_{part}")
            document = json.loads(path.read_text())
            document["somatic_roles"] = manifest["somatic_roles"]
            write_json(path, document)
            for output in manifest["outputs"]:
                if output["role"] == f"mlhp_part_{part}":
                    output["size_bytes"] = path.stat().st_size
                    output["sha256"] = file_sha(path)
        write_json(manifest_path, manifest)
        index_path = self.fixture.run_root / "output_manifests.json"
        index = json.loads(index_path.read_text())
        for item in index["manifests"]:
            if item["sample"] == sample:
                item["sha256"] = file_sha(manifest_path)
        write_json(index_path, index)
        self.assert_error("tree-role-swap", "E_SAMPLE_OUTPUT_BINDING")

    def test_precreated_success_marker_is_rejected_and_not_removed(self) -> None:
        marker = self.fixture.run_root / "_SUCCESS"
        marker.write_text("premature\n", encoding="utf-8")
        result = self.fixture.run("premature-success")
        self.assertEqual(result.returncode, 7)
        self.assertIn("E_STATE_INVALID", self.fixture.summary("premature-success")["error_codes"])
        self.assertTrue(marker.is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)
