#!/usr/bin/env python3
"""Tiny, non-production integration tests for ``run_layered_v3.py``.

All BAM/VCF/sidecar/scientific programs live below ``TemporaryDirectory``.
The tests never read production BAMs, never launch LongPhase, and never create
a production output root.  They exercise the real RunLifecycle/source-bundle
integration with deterministic stub validator, scientific modules, and
verifier bytes.
"""

from __future__ import annotations

import contextlib
import copy
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import textwrap
import unittest
from typing import Any

import pysam

import run_layered_v3 as runner


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def full_artifact(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "identity": {"policy": "full_sha256", "size_bytes": path.stat().st_size, "sha256": sha256(path)},
    }


def indexed_artifact(path: Path) -> dict[str, Any]:
    return {**full_artifact(path), "index": full_artifact(Path(f"{path}.tbi"))}


def write_vcf(path: Path, filter_value: str) -> None:
    plain = path.with_suffix("")
    plain.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t10\t.\tA\tC\t60\t{filter_value}\t.\n",
        encoding="utf-8",
    )
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="vcf", force=True)


VALIDATOR_STUB = r'''#!/usr/bin/env python3
import hashlib, json, os
from pathlib import Path

class ContractError(Exception):
    def __init__(self, code, message, exit_code=2, stage="schema"):
        super().__init__(message); self.code=code; self.exit_code=exit_code; self.stage=stage

def set_method_receipt_repo_root(root):
    global METHOD_RECEIPT_REPO_ROOT
    METHOD_RECEIPT_REPO_ROOT = root

def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def atomic_write_json(path, value, mode=0o644):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temp=path.parent/("."+path.name+".partial."+str(os.getpid()))
    temp.write_text(json.dumps(value, sort_keys=True, indent=2)+"\n", encoding="utf-8")
    os.chmod(temp, mode); os.replace(temp, path)

def validate_and_build_lock(manifest_path, schema_path, validator_path=None):
    document=json.loads(Path(manifest_path).read_text())
    if document.get("mode") == "reject":
        raise ContractError("E_SCHEMA_EMPTY", "fixture rejection")
    lock=document["lock_template"]
    lock["source_manifest"]={
        "path":str(Path(manifest_path).resolve()),
        "byte_sha256":digest(manifest_path),
        "canonical_sha256":hashlib.sha256(canonical(document)).hexdigest(),
    }
    lock["validator"]={
        "path":str(Path(validator_path).resolve()), "sha256":digest(validator_path),
        "schema_path":str(Path(schema_path).resolve()), "schema_sha256":digest(schema_path),
    }
    payload=dict(lock); payload.pop("lock_sha256", None)
    lock["lock_sha256"]=hashlib.sha256(canonical(payload)).hexdigest()
    return lock

def atomic_publish_frozen_lock(path, lock, lock_schema_path):
    if Path(path).exists(): raise ContractError("E_OUTPUT_EXISTS", "exists", exit_code=7)
    atomic_write_json(path, lock, mode=0o444)

def failure_document(error, manifest, frozen):
    return {"valid_lock_written":False,"error":{"code":error.code,"message":str(error)}}
'''


MLHP_STUB = r'''#!/usr/bin/env python3
import json, sys
chromosomes=sys.argv[1].split(",")
out={
 "schema_version":"2.0", "source_token":"frozen-v1",
 "input_funnel":{
  "n_sSNV_scope_input":0,"n_positional_singleton":0,"n_multilocus_pre_cap_groups":0,
  "n_multilocus_pre_cap_sSNV":0,"n_groups_retained":0,"n_sSNV_retained":0,
  "n_groups_capped_by_MAX_SNV":0,"n_sSNV_cap_excluded":0,
  "n_groups_read_unsupported":0,"n_sSNV_read_unsupported":0,"n_sSNV_accounted":0,
  "check_scope_conservation":True},
 "read_tag_census":{"sidecar_missing":0,"sidecar_conflicts":0,
  "alignment_identity_allele_conflicts":0,
  "sidecar_exact_matches":1,"alignment_group_exposures":1},
 "groups":[]}
open(sys.argv[2],"w").write(json.dumps(out))
'''


LAYERED_STUB = r'''#!/usr/bin/env python3
import json, os
open(os.environ["SM_OUT"],"w").write(json.dumps({
 "schema_version":"2.0","source_token":"frozen-v1","input_funnel":{},
 "L1_ssnv_algorithm":{"all_eligible_V1V7_pass":True},"detail":[],"n_detail_units":0}))
'''


REGION_STUB = r'''#!/usr/bin/env python3
import json, os
open(os.environ["SM_OUT"],"w").write(json.dumps({
 "schema_version":"2.0","sample":os.environ["SM_SAMPLE"],"source_token":"frozen-v1",
 "census":{"funnel":{},"n_regions":0},"regions":[]}))
'''


LEDGER_STUB = r'''#!/usr/bin/env python3
import argparse, gzip, json
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument("--sample"); p.add_argument("--caller-raw-vcf")
p.add_argument("--longphase-input-vcf"); p.add_argument("--tree-input-vcf"); p.add_argument("--tree-contract")
p.add_argument("--longphase-input-contract")
p.add_argument("--recalibrated-vcf"); p.add_argument("--mlhp-glob")
p.add_argument("--output-tsv-gz",type=Path); p.add_argument("--output-summary",type=Path); a=p.parse_args()
with gzip.open(a.output_tsv_gz,"wt") as h: h.write("sample\tchrom\tpos\n"+a.sample+"\tchr1\t10\n")
Path(str(a.output_tsv_gz)+".tbi").write_bytes(b"tiny-index")
a.output_summary.write_text(json.dumps({"schema_version":"2.0","sample":a.sample,"pass":True,
 "checks":{"fixture":True},"branch_counts":{"retained":1},"raw_clairs_records":1}))
raise SystemExit(0)
'''


EXIT_SEMANTICS_STUB = r'''#!/usr/bin/env python3
import sys

mode = sys.argv[1]
if mode == "none":
    raise SystemExit()
if mode == "zero":
    raise SystemExit(0)
if mode == "nonzero":
    raise SystemExit(7)
raise SystemExit("explicit failure message")
'''


VERIFIER_STUB = r'''#!/usr/bin/env python3
import argparse, json
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument("--run-root",type=Path); p.add_argument("--frozen-lock")
p.add_argument("--launch-receipt"); p.add_argument("--output",type=Path); a=p.parse_args()
manifests=list((a.run_root/"samples").glob("*/output_manifest.json"))
if len(manifests)!=7: raise SystemExit(7)
a.output.write_text(json.dumps({"all_pass":True,"dataset_count":7,"n_pass":7})+"\n")
(a.run_root/"stub_verifier_called.txt").write_text("bundle verifier called\n")
'''


class LayeredV3RunnerTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix="layered-v3-runner-")
        self.root = Path(self.temp.name)
        self.old_environment = dict(os.environ)
        os.environ.update({"LC_ALL": "C", "TZ": "UTC", "PYTHONHASHSEED": "0"})
        for key in list(os.environ):
            if key.startswith("SM_"):
                del os.environ[key]
        self.sources = self.root / "sources"
        self.sources.mkdir()
        self.inventory = self._make_inventory()
        self.lock_template = self._make_lock_template()

    def tearDown(self) -> None:
        os.environ.clear()
        os.environ.update(self.old_environment)
        self.temp.cleanup()

    def _source(self, name: str, payload: str) -> Path:
        path = self.sources / name
        path.write_text(textwrap.dedent(payload), encoding="utf-8")
        return path

    def test_bundled_executor_preserves_system_exit_semantics(self) -> None:
        script = self._source("exit_semantics.py", EXIT_SEMANTICS_STUB)
        bundle = runner.L.build_source_bundle(
            self.root / "exit-semantics-bundle",
            self.inventory.runner,
            self.inventory.validator,
            self.inventory.verifier,
            (script,),
        )
        self.assertEqual(runner._exec_bundled_script(bundle.manifest_path, script.name, ["none"]), 0)
        self.assertEqual(runner._exec_bundled_script(bundle.manifest_path, script.name, ["zero"]), 0)
        self.assertEqual(runner._exec_bundled_script(bundle.manifest_path, script.name, ["nonzero"]), 7)
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            code = runner._exec_bundled_script(bundle.manifest_path, script.name, ["message"])
        self.assertEqual(code, 1)
        self.assertEqual(stderr.getvalue(), "explicit failure message\n")

    def _make_inventory(self) -> runner.SourceInventory:
        validator = self._source("validate_layered_v3_inputs.py", VALIDATOR_STUB)
        verifier = self._source("verify_layered_v3.py", VERIFIER_STUB)
        scientific = (
            self._source("sm_linkage_genomewide.py", "# frozen dependency\n"),
            self._source("sm_multilocus_combinations.py", MLHP_STUB),
            self._source("tree_enumeration_solver.py", "# frozen dependency\n"),
            self._source("layered_tree_reconstruction.py", LAYERED_STUB),
            self._source("build_region_view.py", REGION_STUB),
            self._source("build_ssnv_site_ledger.py", LEDGER_STUB),
        )
        schemas = tuple(
            self._source(name, json.dumps({"fixture": name}) + "\n")
            for name in runner.SCHEMA_BASENAMES
        )
        return runner.SourceInventory(
            runner=Path(runner.__file__).resolve(),
            lifecycle=(runner.REPO_ROOT / "scripts/layered_v3_lifecycle.py").resolve(),
            validator=validator,
            verifier=verifier,
            scientific=scientific,
            schemas=schemas,
            distributions=(),
            tools=(),
        )

    def _make_lock_template(self) -> dict[str, Any]:
        artifacts = self.root / "artifacts"
        artifacts.mkdir()
        bam = artifacts / "tumor.bam"
        bam.write_bytes(b"tiny-bam")
        raw = artifacts / "caller_raw.vcf.gz"
        longphase_input = artifacts / "longphase_input.vcf.gz"
        caller_pass_baseline = artifacts / "caller_pass_baseline.vcf.gz"
        recalibrated_all = artifacts / "recalibrated_all.vcf.gz"
        tree = artifacts / "tree_pass.vcf.gz"
        for path, filter_value in (
            (raw, "PASS"),
            (longphase_input, "LowQual"),
            (caller_pass_baseline, "PASS"),
            (recalibrated_all, "LowQual"),
            (tree, "PASS"),
        ):
            write_vcf(path, filter_value)
        sidecar = artifacts / "tags.tsv.gz"
        with __import__("gzip").open(sidecar, "wt") as handle:
            handle.write("#CHROM\tSTART0\tEND0\tQNAME\tFLAG\tMAPQ\tCIGAR_B2\tHP\tPS\n")
            handle.write("chr1\t0\t20\tr1\t0\t60\t1111111111111111\t1\t10\n")
        sidecar_index = artifacts / "tags.tsv.gz.tbi"
        sidecar_index.write_bytes(b"tiny-sidecar-index")
        validation = artifacts / "sidecar_validation.json"
        producer_receipt = artifacts / "producer_capture_receipt_v3.json"
        write_json(validation, {"pass": True})
        write_json(producer_receipt, {"pass": True})
        somatic = {
            "caller_raw_vcf": indexed_artifact(raw),
            "longphase_input_vcf": indexed_artifact(longphase_input),
            "caller_pass_baseline_vcf": indexed_artifact(caller_pass_baseline),
            "longphase_recalibrated_all_vcf": indexed_artifact(recalibrated_all),
            "longphase_recalibrated_pass_vcf": indexed_artifact(tree),
            "tree_vcf": indexed_artifact(tree),
        }
        samples = []
        for sample, biological_id in runner.EXPECTED_BINDING.items():
            samples.append(
                {
                    "sample": sample,
                    "biological_id": biological_id,
                    "pass": True,
                    "alignment_payload": {
                        "path": str(bam),
                        "storage_identity_v1": {"identity_sha256": "a" * 64},
                    },
                    "somatic": copy.deepcopy(somatic),
                    "read_tags": {
                        "sidecar": full_artifact(sidecar),
                        "index": full_artifact(sidecar_index),
                        "validation": full_artifact(validation),
                        "producer_capture_receipt": full_artifact(producer_receipt),
                        "duplicate_identity_policy": "collapse_redundant_rows_with_identical_HP_PS",
                        "subject_binding": {
                            "sample": sample,
                            "duplicate_identity_policy": "collapse_redundant_rows_with_identical_HP_PS",
                            "longphase_recalibrated_pass_vcf_sha256": somatic[
                                "longphase_recalibrated_pass_vcf"
                            ]["identity"]["sha256"],
                            "longphase_input_vcf_sha256": somatic["longphase_input_vcf"]["identity"]["sha256"],
                            "caller_pass_baseline_vcf_sha256": somatic["caller_pass_baseline_vcf"]["identity"]["sha256"],
                            "mapped_alignment_count": 1,
                            "identity_unique_count": 1,
                            "duplicate_count": 0,
                            "conflict_count": 0,
                        },
                        "producer_policy": {},
                        "producer_evidence": {},
                    },
                    "copy_number": {
                        "availability": "unavailable",
                        "source": "unavailable",
                        "semantics": "missing; never interpreted neutral",
                        "coordinate_system": None,
                        "unlisted_position_semantics": "unavailable",
                        "allowed_states": [],
                        "overlap_policy": "not_applicable",
                        "reason": "No reviewed measured CN source is available",
                        "cn_bed": None,
                        "cn_int_gain": None,
                        "cn_int_loss": None,
                        "integration_json": None,
                    },
                }
            )
        template = {
            "$schema": "InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/schemas/layered_input_lock_v1.schema.json",
            "schema_name": "intersubmod.layered_input_lock",
            "schema_version": "1.0.0",
            "lock_id": "tiny.frozen",
            "created_at_utc": "2026-07-11T00:00:00Z",
            "source_manifest": {},
            "validator": {},
            "dataset_count": 7,
            "biological_sample_count": 6,
            "analysis_contract": {
                "task_type": "comprehensive_validation",
                "scope": {"name": "whole_autosomes_chr1_22", "contigs": list(runner.AUTOSOMES)},
                "longphase_input_contract": "normalized_ClairS_raw_all",
                "tree_input_contract": "longphase_s_recalibrated_FILTER_PASS",
                "duplicate_identity_policy": "collapse_redundant_rows_with_identical_HP_PS",
            },
            "production_summary": {
                "expected_dataset_count": 7,
                "completed_dataset_count": 7,
                "passed_dataset_count": 7,
                "all_pass": True,
                "datasets": [{"sample": sample, "pass": True} for sample in runner.EXPECTED_BINDING],
            },
            "all_pass": True,
            "samples": samples,
            "lock_sha256": "0" * 64,
        }
        return template

    def _manifest(self, name: str, *, mode: str = "pass", template: dict[str, Any] | None = None) -> Path:
        path = self.root / f"{name}.json"
        write_json(path, {"mode": mode, "lock_template": template or copy.deepcopy(self.lock_template)})
        return path

    def _config(self, manifest: Path, run_id: str, *, preflight_only: bool = False) -> runner.RunConfig:
        return runner.RunConfig(
            manifest=manifest,
            run_parent=self.root / "runs",
            run_id=run_id,
            parallel_samples=2,
            heartbeat_interval=0.05,
            preflight_only=preflight_only,
            min_logical_cpus=1,
            min_available_memory_gib=0.0,
            min_free_disk_gib=0.0,
            max_load_per_cpu=1000000.0,
            max_iowait_percent=100.0,
            resource_sample_seconds=0.01,
            max_nfs_read_mbps=1000000000.0,
            nfs_mountpoint=Path("/big8_disk"),
        )

    def _run_pipeline(
        self,
        config: runner.RunConfig,
        *,
        after_bundle_hook: Any = None,
        process_snapshot_provider: Any = None,
    ) -> Path:
        return runner.run_pipeline(
            config,
            self.inventory,
            after_bundle_hook=after_bundle_hook,
            process_snapshot_provider=process_snapshot_provider or (lambda: []),
        )

    def assert_no_formal_root_or_worker(self, run_id: str) -> None:
        run_parent = self.root / "runs"
        self.assertFalse((run_parent / run_id).exists())
        staging = list(run_parent.glob(f".{run_id}.staging.*"))
        self.assertEqual(len(staging), 1)
        self.assertFalse((staging[0] / "samples").exists())

    def test_reviewed_production_contract_hashes_are_still_frozen(self) -> None:
        inventory = runner.production_inventory()
        runner._validate_source_inventory(inventory)
        self.assertEqual(dict(inventory.expected_sha256), runner.FROZEN_CONTRACT_SHA256)

    def test_production_resource_defaults_and_allowlist_are_frozen(self) -> None:
        config = runner.RunConfig(
            manifest=self.root / "unused.json",
            run_parent=self.root / "unused-runs",
            run_id="production-defaults",
        )
        observed = {
            "min_logical_cpus": config.min_logical_cpus,
            "min_available_memory_gib": config.min_available_memory_gib,
            "min_free_disk_gib": config.min_free_disk_gib,
            "max_load_per_cpu": config.max_load_per_cpu,
            "max_iowait_percent": config.max_iowait_percent,
            "resource_sample_seconds": config.resource_sample_seconds,
            "max_nfs_read_mbps": config.max_nfs_read_mbps,
            "nfs_mountpoint": str(config.nfs_mountpoint),
        }
        self.assertEqual(observed, runner.PRODUCTION_RESOURCE_POLICY)
        self.assertTrue(
            {
                "SM_MIN_AVAILABLE_MEMORY_GIB",
                "SM_MIN_FREE_DISK_GIB",
                "SM_MAX_IOWAIT_PERCENT",
                "SM_RESOURCE_SAMPLE_SECONDS",
                "SM_MAX_NFS_READ_MBPS",
                "SM_NFS_MOUNTPOINT",
            }.issubset(runner.SM_ALLOWLIST)
        )
        self.assertTrue(
            {
                "longphase-s",
                "run_longphase_raw_all_production_sidecars.sh",
                "capture_longphase_tagged_bam_sidecar.py",
                "validate_streamed_longphase_sidecar.py",
            }.issubset(runner.CONFLICTING_PROCESS_BASENAMES)
        )
        runner._validate_production_resource_policy(config, runner.production_inventory())
        drifted = runner.RunConfig(**{**config.__dict__, "resource_sample_seconds": 0.01})
        with self.assertRaises(runner.ContractError) as caught:
            runner._validate_production_resource_policy(drifted, runner.production_inventory())
        self.assertEqual(caught.exception.code, "E_RESOURCE_POLICY_DRIFT")

    def test_completion_supervisor_exports_deterministic_env_and_reuses_manifest_pair(self) -> None:
        source = Path(__file__).resolve().with_name("complete_raw_all_layered_v3_validation.sh")
        text = source.read_text(encoding="utf-8")
        self.assertIn("export LC_ALL=C TZ=UTC PYTHONHASHSEED=0", text)
        self.assertIn("USE_EXISTING_MANIFESTS=1", text)
        self.assertIn('INPUT_PINS+=("$CANONICAL_MANIFEST" "$SENSITIVITY_MANIFEST")', text)
        self.assertIn("E_MANIFEST_PAIR", text)
        self.assertIn("--parallel-samples 4", text)

    def test_parallel_sample_limit_accepts_four_and_rejects_five(self) -> None:
        manifest = self._manifest("parallel-limit")
        config = self._config(manifest, "parallel-limit")
        four = runner.RunConfig(**{**config.__dict__, "parallel_samples": 4})
        runner._validate_config(four)
        five = runner.RunConfig(**{**config.__dict__, "parallel_samples": 5})
        with self.assertRaises(runner.ContractError) as caught:
            runner._validate_config(five)
        self.assertEqual(caught.exception.code, "E_PARAM_RANGE")

    def test_bundled_validator_receives_repo_root_from_manifest_ancestor(self) -> None:
        source = Path(runner.__file__).read_text(encoding="utf-8")
        self.assertIn('parent.name == "InterSubMod"', source)
        self.assertIn("validator.set_method_receipt_repo_root(repo_roots[0] if repo_roots else None)", source)

    def test_invalid_manifest_has_no_formal_root_or_worker(self) -> None:
        manifest = self._manifest("invalid", mode="reject")
        with self.assertRaises(runner.ContractError) as caught:
            self._run_pipeline(self._config(manifest, "invalid-manifest"))
        self.assertEqual(caught.exception.code, "E_SCHEMA_EMPTY")
        self.assert_no_formal_root_or_worker("invalid-manifest")

    def test_resource_gate_rejects_before_bundle_or_worker(self) -> None:
        manifest = self._manifest("resource-gate")
        config = self._config(manifest, "resource-gate")
        config = runner.RunConfig(
            **{
                **config.__dict__,
                "min_logical_cpus": (os.cpu_count() or 1) + 1,
            }
        )
        with self.assertRaises(runner.ContractError) as caught:
            self._run_pipeline(config)
        self.assertEqual(caught.exception.code, "E_RESOURCE_GATE")
        self.assert_no_formal_root_or_worker("resource-gate")
        staging = list((self.root / "runs").glob(".resource-gate.staging.*"))[0]
        evidence = json.loads((staging / "process_observation.json").read_text())
        self.assertFalse((staging / "source_bundle").exists())
        self.assertFalse(evidence["resources"]["checks"]["logical_cpus"])
        self.assertFalse(evidence["pass"])

    def test_missing_nfs_counter_writes_failure_evidence_and_rejects(self) -> None:
        manifest = self._manifest("missing-nfs")
        config = self._config(manifest, "missing-nfs")
        config = runner.RunConfig(**{**config.__dict__, "nfs_mountpoint": Path("/missing-nfs")})
        with self.assertRaises(runner.ContractError) as caught:
            self._run_pipeline(config)
        self.assertEqual(caught.exception.code, "E_NFS_BASELINE_UNAVAILABLE")
        self.assert_no_formal_root_or_worker("missing-nfs")
        staging = list((self.root / "runs").glob(".missing-nfs.staging.*"))[0]
        evidence_path = staging / "process_observation.json"
        evidence = json.loads(evidence_path.read_text())
        self.assertFalse((staging / "source_bundle").exists())
        self.assertFalse(evidence["resources"]["checks"]["evidence_available"])
        self.assertEqual(
            evidence["resources"]["error"]["code"], "E_NFS_BASELINE_UNAVAILABLE"
        )
        self.assertIn(sha256(evidence_path), caught.exception.message)

    def test_nfs_read_rate_gate_rejects_before_bundle(self) -> None:
        manifest = self._manifest("nfs-rate-gate")
        config = self._config(manifest, "nfs-rate-gate")
        config = runner.RunConfig(**{**config.__dict__, "max_nfs_read_mbps": 0.0})
        with self.assertRaises(runner.ContractError) as caught:
            self._run_pipeline(config)
        self.assertEqual(caught.exception.code, "E_RESOURCE_GATE")
        self.assert_no_formal_root_or_worker("nfs-rate-gate")
        staging = list((self.root / "runs").glob(".nfs-rate-gate.staging.*"))[0]
        evidence = json.loads((staging / "process_observation.json").read_text())
        self.assertFalse((staging / "source_bundle").exists())
        self.assertFalse(evidence["resources"]["checks"]["nfs_read_mbps"])
        baseline = evidence["resources"]["observed"]["nfs_read_baseline"]
        self.assertEqual(baseline["threshold_max_mbps_decimal"], 0.0)
        self.assertGreaterEqual(baseline["read_bytes_recv_delta"], 0)

    def test_longphase_process_conflict_rejects_before_bundle(self) -> None:
        manifest = self._manifest("longphase-conflict")
        provider = lambda: [
            {
                "pid": 424242,
                "pid_start_time": "12345",
                "argv": ["/frozen/bin/longphase-s", "--tumor-bam", "/tiny/tumor.bam"],
            }
        ]
        with self.assertRaises(runner.ContractError) as caught:
            self._run_pipeline(
                self._config(manifest, "longphase-conflict"),
                process_snapshot_provider=provider,
            )
        self.assertEqual(caught.exception.code, "E_CONFLICTING_FULL_RUN")
        self.assert_no_formal_root_or_worker("longphase-conflict")
        staging = list((self.root / "runs").glob(".longphase-conflict.staging.*"))[0]
        evidence = json.loads((staging / "process_observation.json").read_text())
        self.assertFalse((staging / "source_bundle").exists())
        self.assertEqual(evidence["conflicts"][0]["matched_basenames"], ["longphase-s"])
        self.assertEqual(
            evidence["conflicts"][0]["observation_points"],
            ["baseline_end", "baseline_start"],
        )
        self.assertFalse(evidence["pass"])

    def test_conflict_starting_during_baseline_is_detected_at_end(self) -> None:
        manifest = self._manifest("baseline-toctou")
        snapshots = iter(
            [
                [],
                [
                    {
                        "pid": 434343,
                        "pid_start_time": "54321",
                        "argv": [
                            "/bin/bash",
                            "/frozen/run_longphase_raw_all_production_sidecars.sh",
                        ],
                    }
                ],
            ]
        )

        def provider() -> Any:
            return next(snapshots)

        with self.assertRaises(runner.ContractError) as caught:
            self._run_pipeline(
                self._config(manifest, "baseline-toctou"),
                process_snapshot_provider=provider,
            )
        self.assertEqual(caught.exception.code, "E_CONFLICTING_FULL_RUN")
        self.assert_no_formal_root_or_worker("baseline-toctou")
        staging = list((self.root / "runs").glob(".baseline-toctou.staging.*"))[0]
        evidence = json.loads((staging / "process_observation.json").read_text())
        self.assertEqual(evidence["process_observations"]["start"]["conflicts"], [])
        self.assertEqual(
            evidence["process_observations"]["end"]["conflicts"][0]["matched_basenames"],
            ["run_longphase_raw_all_production_sidecars.sh"],
        )
        self.assertEqual(evidence["conflicts"][0]["observation_points"], ["baseline_end"])

    def test_empty_and_extra_dataset_sets_are_rejected_before_publish(self) -> None:
        empty = copy.deepcopy(self.lock_template)
        empty["samples"] = []
        empty["dataset_count"] = 0
        with self.assertRaises(runner.ContractError) as caught_empty:
            self._run_pipeline(self._config(self._manifest("empty", template=empty), "empty-set"))
        self.assertEqual(caught_empty.exception.code, "E_DATASET_SET_MISMATCH")
        self.assert_no_formal_root_or_worker("empty-set")

        extra = copy.deepcopy(self.lock_template)
        eighth = copy.deepcopy(extra["samples"][0])
        eighth["sample"] = "EXTRA"
        eighth["biological_id"] = "EXTRA"
        eighth["read_tags"]["subject_binding"]["sample"] = "EXTRA"
        extra["samples"].append(eighth)
        extra["dataset_count"] = 8
        with self.assertRaises(runner.ContractError) as caught_extra:
            self._run_pipeline(self._config(self._manifest("extra", template=extra), "extra-set"))
        self.assertEqual(caught_extra.exception.code, "E_DATASET_SET_MISMATCH")
        self.assert_no_formal_root_or_worker("extra-set")

    def test_preflight_only_is_terminal_and_starts_zero_scientific_workers(self) -> None:
        manifest = self._manifest("preflight")
        root = self._run_pipeline(self._config(manifest, "happy-preflight", preflight_only=True))
        report = json.loads((root / "preflight_verification.json").read_text())
        marker = json.loads((root / "_SUCCESS").read_text())
        self.assertTrue(report["pass"])
        self.assertEqual(report["scientific_workers_started"], 0)
        self.assertEqual(marker["extra"]["mode"], "preflight_only")
        self.assertFalse((root / "samples").exists())

    def test_explicit_clairs_sensitivity_preflight_keeps_longphase_tag_subject(self) -> None:
        template = copy.deepcopy(self.lock_template)
        template["analysis_contract"]["task_type"] = "backbone_sensitivity"
        template["analysis_contract"]["tree_input_contract"] = runner.SENSITIVITY_TREE_INPUT_CONTRACT
        for item in template["samples"]:
            somatic = item["somatic"]
            somatic["tree_vcf"] = copy.deepcopy(somatic["caller_pass_baseline_vcf"])
            self.assertEqual(
                item["read_tags"]["subject_binding"]["longphase_recalibrated_pass_vcf_sha256"],
                somatic["longphase_recalibrated_pass_vcf"]["identity"]["sha256"],
            )
        manifest = self._manifest("sensitivity-preflight", template=template)
        root = self._run_pipeline(
            self._config(manifest, "sensitivity-preflight", preflight_only=True)
        )
        receipt = json.loads((root / "launch_receipt.json").read_text())
        frozen = json.loads((root / "frozen_input_lock.json").read_text())
        self.assertEqual(receipt["extra"]["task_type"], "backbone_sensitivity")
        self.assertEqual(
            receipt["extra"]["tree_backbone_role"], "clairs_filter_pass_sensitivity"
        )
        self.assertTrue(all(
            item["somatic"]["tree_vcf"] == item["somatic"]["caller_pass_baseline_vcf"]
            for item in frozen["samples"]
        ))

    def test_fixed_scope_lock_rejects_a_different_run_id(self) -> None:
        run_parent = self.root / "runs"
        with runner.GlobalScopeLock(run_parent, "first-full-run"):
            with self.assertRaises(runner.ContractError) as caught:
                self._run_pipeline(
                    self._config(self._manifest("second-lock"), "different-run-id")
                )
        self.assertEqual(caught.exception.code, "E_GLOBAL_RUN_LOCKED")
        self.assertFalse((run_parent / "different-run-id").exists())
        self.assertFalse(list(run_parent.glob(".different-run-id.staging.*")))

    def test_source_mutation_after_bundle_does_not_change_worker_bytes(self) -> None:
        full_template = copy.deepcopy(self.lock_template)
        cn_bed = self.root / "artifacts/measured_cn.bed"
        cn_bed.write_text("chr1\t0\t100\tneutral\n", encoding="utf-8")
        measured = {
            "availability": "measured",
            "source": "reviewed_tiny_fixture",
            "semantics": "explicit segments; unlisted positions are neutral",
            "coordinate_system": "0_based_half_open",
            "unlisted_position_semantics": "neutral",
            "allowed_states": ["gain", "loss", "loh", "neutral"],
            "overlap_policy": "forbid",
            "reason": None,
            "cn_bed": full_artifact(cn_bed),
            "cn_int_gain": None,
            "cn_int_loss": None,
            "integration_json": None,
        }
        next(item for item in full_template["samples"] if item["sample"] == "H1437")[
            "copy_number"
        ] = measured
        manifest = self._manifest("full", template=full_template)
        mutable = self.sources / "sm_multilocus_combinations.py"

        def mutate(_bundle_manifest: Path) -> None:
            mutable.write_text("raise SystemExit(99)\n", encoding="utf-8")

        root = self._run_pipeline(
            self._config(manifest, "happy-full"),
            after_bundle_hook=mutate,
        )
        self.assertTrue((root / "_SUCCESS").is_file())
        self.assertTrue((root / "stub_verifier_called.txt").is_file())
        observation = json.loads((root / "process_observation.json").read_text())
        self.assertTrue(observation["pass"])
        self.assertEqual(observation["conflicts"], [])
        self.assertTrue(observation["global_scope_lock"]["held"])
        self.assertTrue(observation["resources"]["checks"]["nfs_read_mbps"])
        nfs = observation["resources"]["observed"]["nfs_read_baseline"]
        self.assertEqual(nfs["mountpoint"], "/big8_disk")
        self.assertGreaterEqual(nfs["read_bytes_recv_delta"], 0)
        self.assertGreater(nfs["sample_seconds"], 0)
        for sample in runner.EXPECTED_BINDING:
            part = json.loads((root / "samples" / sample / "mlhp_part_1.json").read_text())
            self.assertEqual(part["source_token"], "frozen-v1")
            output_manifest = json.loads(
                (root / "samples" / sample / "output_manifest.json").read_text()
            )
            region = json.loads(
                (root / "samples" / sample / f"layered_region_view_{sample}.json").read_text()
            )
            self.assertEqual(output_manifest["copy_number_contract"], region["copy_number_contract"])
            if sample == "H1437":
                self.assertEqual(output_manifest["copy_number_contract"]["source"], "reviewed_tiny_fixture")
            else:
                self.assertEqual(
                    output_manifest["copy_number_contract"]["reason"],
                    "No reviewed measured CN source is available",
                )
        receipt = json.loads((root / "launch_receipt.json").read_text())
        self.assertEqual(receipt["extra"]["analysis_params"]["parallel_samples"], 2)
        self.assertEqual(receipt["extra"]["analysis_params"]["contigs"], list(runner.AUTOSOMES))
        policy = receipt["extra"]["analysis_params"]["resource_policy"]
        self.assertEqual(policy["sample_seconds"], 0.01)
        self.assertEqual(policy["nfs_mountpoint"], "/big8_disk")


if __name__ == "__main__":
    unittest.main(verbosity=2)
