#!/usr/bin/env python3
"""Synthetic red-team tests for the immutable M2 release contract.

All fixtures are explicitly non-evidentiary.  They nevertheless use the full
7 datasets x 6 direct roles receipt/manifest shape so a structurally plausible
forgery is tested, rather than a toy count-only receipt.
"""

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import shutil
import stat
import sys
import tempfile
import unittest
from unittest import mock


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import freeze_m2_release_contract as freezer  # noqa: E402


SYNTHETIC_RUNTIME = {
    "python": {"executable": "/synthetic/python", "version": "0.0", "implementation": "synthetic"},
    "packages": {"numpy": "synthetic", "scipy": "synthetic", "pysam": "synthetic"},
    "samtools": {
        "path": "/synthetic/samtools", "version_line": "synthetic",
        "htslib_version_line": "synthetic",
    },
    "platform": "synthetic-test-platform",
}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_path(path: pathlib.Path) -> str:
    return sha256_bytes(path.read_bytes())


def full_identity(path: pathlib.Path) -> dict:
    payload = path.read_bytes()
    return {"policy": "full_sha256", "size_bytes": len(payload), "sha256": sha256_bytes(payload)}


def storage_identity(bam: pathlib.Path, bai: pathlib.Path) -> dict:
    logical = bam.lstat()
    resolved = bam.resolve(strict=True)
    observed = resolved.stat()
    payload = resolved.read_bytes()
    length = min(1024 * 1024, len(payload))
    offsets = (0, max(0, (len(payload) - length) // 2), max(0, len(payload) - length))
    value = {
        "policy": "storage_identity_v1",
        "assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
        "is_full_content_hash": False,
        "requested_path": str(bam),
        "realpath": str(resolved),
        "logical_is_symlink": False,
        "logical_size_bytes": logical.st_size,
        "logical_mtime_ns": logical.st_mtime_ns,
        "st_dev": observed.st_dev,
        "st_ino": observed.st_ino,
        "size_bytes": observed.st_size,
        "mtime_ns": observed.st_mtime_ns,
        "ctime_ns": observed.st_ctime_ns,
        "chunk_size_bytes": 1024 * 1024,
        "chunks": [
            {
                "label": label, "offset": offset, "length": length,
                "sha256": sha256_bytes(payload[offset:offset + length]),
            }
            for label, offset in zip(("first", "middle", "last"), offsets)
        ],
        "index": {"path": str(bai), "identity": full_identity(bai)},
    }
    value["identity_sha256"] = freezer.canonical_json_sha256(value)
    return value


def write_authenticated_json(path: pathlib.Path, document: dict) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    value = dict(document)
    value["receipt_integrity"] = {
        "scheme": "external_sha256_sidecar_v1",
        "sidecar_name": f"{path.name}.sha256",
        "covers": path.name,
    }
    payload = json.dumps(value, ensure_ascii=False, allow_nan=False, indent=2).encode() + b"\n"
    path.write_bytes(payload)
    digest = sha256_bytes(payload)
    path.with_name(f"{path.name}.sha256").write_text(f"{digest}  {path.name}\n", encoding="ascii")
    return digest


def make_pre_document(
    manifest_path: pathlib.Path,
    manifest_document: dict,
    schema_path: pathlib.Path,
    verifier_path: pathlib.Path,
    verifier_sha256: str,
    authority: freezer._Authority,
    *,
    covers: str = "synthetic_fresh.json",
) -> dict:
    derived, inventory = freezer._derive_manifest_role_contract(manifest_document)
    audit_artifacts = [{**row, "observed": row["expected"], "match": True} for row in derived]
    snapshot = {
        "manifest_sha256": authority.expected_manifest_sha256,
        "schema_sha256": authority.expected_schema_sha256,
        "datasets": list(freezer.DATASETS),
        "artifacts": [
            {
                "dataset": row["dataset"], "role": row["role"], "policy": row["policy"],
                "path": row["path"], "observed": row["expected"],
            }
            for row in derived
        ],
    }
    return {
        "schema_name": freezer.PRE_SCHEMA_NAME,
        "schema_version": freezer.PRE_SCHEMA_VERSION,
        "task_type": "B_COMPREHENSIVE_VALIDATION",
        "mode": "PRE",
        "authority_mode": authority.expected_pre_authority_mode,
        "validation_evidence_eligible": authority.validation_evidence_eligible,
        "authority": {
            "canonical_manifest_sha256": freezer.CANONICAL_MANIFEST_SHA256,
            "canonical_schema_sha256": freezer.CANONICAL_SCHEMA_SHA256,
            "canonical_schema_path": str(schema_path.absolute()),
            "selected_authority_is_canonical": authority.validation_evidence_eligible,
            "test_only_override": not authority.validation_evidence_eligible,
        },
        "scope": {
            "technical_datasets": 7, "biological_samples": 6, "chromosomes": "chr1-chr22",
            "tasks": 154, "datasets": list(freezer.DATASETS), "direct_input_artifacts": 42,
        },
        "manifest": {
            "path": str(manifest_path.absolute()),
            "sha256": authority.expected_manifest_sha256,
            "expected_sha256": authority.expected_manifest_sha256,
        },
        "canonical_schema": {
            "path": str(schema_path.absolute()), "sha256": authority.expected_schema_sha256,
            "expected_sha256": authority.expected_schema_sha256,
        },
        "verifier": {"path": str(verifier_path.absolute()), "sha256": verifier_sha256},
        "assurance": dict(freezer.PRE_ASSURANCE),
        "input_identity_snapshot": snapshot,
        "input_identity_snapshot_sha256": freezer.canonical_json_sha256(snapshot),
        "artifact_audit": {
            "artifacts": audit_artifacts, "role_inventory": inventory,
            "n_artifacts": 42, "n_unique_resolved_files": 42,
            "n_storage_identity_v1": 7, "n_full_sha256": 35,
            "n_sampled_bam_chunks": 21, "n_mismatches": 0,
        },
        "compare_to": None,
        "checks": {name: True for name in freezer.PRE_CHECK_NAMES},
        "all_pass": True,
        "receipt_integrity": {
            "scheme": "external_sha256_sidecar_v1", "sidecar_name": f"{covers}.sha256",
            "covers": covers,
        },
    }


def make_writable_for_cleanup(root: pathlib.Path) -> None:
    if not root.exists():
        return
    for path in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        try:
            if path.is_dir() and not path.is_symlink():
                path.chmod(0o755)
            elif not path.is_symlink():
                path.chmod(0o644)
        except FileNotFoundError:
            pass
    if not root.is_symlink():
        root.chmod(0o755)


class SyntheticReleaseFixture(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._cleanup)
        self.base = pathlib.Path(self._temporary.name)
        self.repo = self.base / "synthetic_repo"
        self.repo.mkdir()
        for index, (role, relative) in enumerate(freezer.SOURCE_ALLOWLIST):
            path = self.repo / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"# synthetic source {index} {role}\nVALUE = {index!r}\n", encoding="utf-8")
        self.schema_relative = pathlib.Path("synthetic_authority/schema.json")
        self.schema = self.repo / self.schema_relative
        self.schema.parent.mkdir(parents=True, exist_ok=True)
        self.schema.write_text('{"synthetic_schema":true}\n', encoding="utf-8")

        self.inputs = self.base / "inputs"
        samples = []
        for dataset_index, dataset in enumerate(freezer.DATASETS):
            folder = self.inputs / dataset
            folder.mkdir(parents=True)
            paths = {
                "bam": folder / "alignment.bam", "bai": folder / "alignment.bam.bai",
                "tree": folder / "tree.vcf.gz", "tree_index": folder / "tree.vcf.gz.csi",
                "sidecar": folder / "read_tags.tsv.gz", "sidecar_index": folder / "read_tags.tsv.gz.tbi",
            }
            for role_index, path in enumerate(paths.values()):
                path.write_bytes(f"synthetic-{dataset_index}-{role_index}-{path.name}\n".encode())
            samples.append({
                "sample": dataset,
                "biological_id": "HCC1395" if dataset == "HCC1395_DORADO" else dataset,
                "platform": "ONT_DORADO" if dataset == "HCC1395_DORADO" else "ONT",
                "replicate_role": "platform_replica" if dataset == "HCC1395_DORADO" else "canonical",
                "alignment_payload": {
                    "path": str(paths["bam"]), "index_path": str(paths["bai"]),
                    "identity_policy": "storage_identity_v1", "embedded_hp_ps_policy": "ignore",
                    "storage_identity_v1": storage_identity(paths["bam"], paths["bai"]),
                },
                "somatic": {
                    "backbone_role": "longphase_s_recalibrated_filter_pass",
                    "tree_vcf": {
                        "path": str(paths["tree"]), "identity": full_identity(paths["tree"]),
                        "index": {"path": str(paths["tree_index"]), "identity": full_identity(paths["tree_index"])},
                    },
                },
                "read_tags": {
                    "authority": "external_coordinate_sidecar", "identity_schema": "coordinate_join_v1",
                    "fallback_policy": "forbidden",
                    "sidecar": {"path": str(paths["sidecar"]), "identity": full_identity(paths["sidecar"])},
                    "index": {"path": str(paths["sidecar_index"]), "identity": full_identity(paths["sidecar_index"])},
                },
            })
        self.manifest_document = {"samples": samples, "synthetic_test_only": True}
        self.manifest = self.base / "synthetic_manifest.json"
        self.manifest.write_text(json.dumps(self.manifest_document, indent=2) + "\n", encoding="utf-8")
        self.authority = freezer._synthetic_test_authority(
            sha256_path(self.manifest), self.schema_relative, sha256_path(self.schema)
        )
        self.input_verifier_source = self.repo / dict(freezer.SOURCE_ALLOWLIST)["input_identity_verifier"]
        self.pre = self.base / "synthetic_pre_identity.json"
        pre_document = make_pre_document(
            self.manifest, self.manifest_document, self.schema, self.input_verifier_source,
            sha256_path(self.input_verifier_source), self.authority, covers=self.pre.name,
        )
        write_authenticated_json(self.pre, pre_document)
        self.contract = self.base / "synthetic_contract"
        self.fresh_calls: list[dict] = []

    def _cleanup(self) -> None:
        make_writable_for_cleanup(self.base)
        self._temporary.cleanup()

    def fresh_factory(self, **kwargs: object) -> dict:
        self.fresh_calls.append(dict(kwargs))
        return make_pre_document(
            pathlib.Path(kwargs["manifest_path"]), dict(kwargs["manifest_document"]),
            pathlib.Path(kwargs["schema_path"]), pathlib.Path(kwargs["verifier_path"]),
            str(kwargs["verifier_sha256"]), kwargs["authority"],
        )

    def freeze(self) -> dict:
        return freezer.freeze_release_contract(
            self.manifest, self.pre, self.contract, _repo_root=self.repo,
            _authority=self.authority, _runtime=SYNTHETIC_RUNTIME,
            _fresh_pre_factory=self.fresh_factory,
        )

    def verify(self, *, runtime: dict = SYNTHETIC_RUNTIME, factory=None) -> dict:
        return freezer.verify_release_contract(
            self.contract, _test_authority=self.authority, _runtime=runtime,
            _fresh_pre_factory=self.fresh_factory if factory is None else factory,
        )

    def rewrite_run_manifest(self, mutate) -> None:
        path = self.contract / "m2_run_manifest.json"
        sidecar = path.with_name(f"{path.name}.sha256")
        self.contract.chmod(0o755)
        path.chmod(0o644)
        sidecar.chmod(0o644)
        document = json.loads(path.read_text(encoding="utf-8"))
        mutate(document)
        write_authenticated_json(path, document)
        path.chmod(0o444)
        sidecar.chmod(0o444)
        self.contract.chmod(0o555)


class PositiveContractTests(SyntheticReleaseFixture):
    def test_full_synthetic_freeze_and_verify(self) -> None:
        summary = self.freeze()
        self.assertEqual(len(self.fresh_calls), 1)
        self.assertEqual(pathlib.Path(self.fresh_calls[0]["verifier_path"]), self.input_verifier_source)
        self.assertEqual(summary["n_snapshot_files"], 11)
        self.assertFalse(summary["validation_evidence_eligible"])
        path = self.contract / "m2_run_manifest.json"
        document = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(document["entrypoints"]["canonical_manifest_copy"], "input_contract/canonical_manifest.json")
        self.assertEqual(len(document["source_snapshot"]["entries"]), 11)
        self.assertEqual(document["producer"]["role"], "release_contract_freezer")
        self.assertEqual(document["runtime_semantic_sha256"], freezer.canonical_json_sha256(SYNTHETIC_RUNTIME))
        self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o444)
        for entry in document["source_snapshot"]["entries"]:
            copy = self.contract / entry["snapshot"]["path"]
            self.assertEqual(copy.stat().st_nlink, 1)
            self.assertEqual(entry["source"]["sha256"], entry["snapshot"]["sha256"])
        verification = self.verify()
        self.assertEqual(len(self.fresh_calls), 2)
        self.assertEqual(
            pathlib.Path(self.fresh_calls[1]["verifier_path"]),
            self.contract / document["source_snapshot"]["entrypoints"]["input_identity_verifier"],
        )
        self.assertTrue(verification["all_pass"])
        self.assertEqual(verification["snapshot"]["n_files"], 11)
        self.assertTrue(verification["fresh_input_identity_verification"]["exactly_equals_supplied_pre_snapshot"])

    def test_freeze_uses_latest_source_bytes(self) -> None:
        role, relative = freezer.SOURCE_ALLOWLIST[2]
        source = self.repo / relative
        source.write_text("# latest synthetic bytes\nVERSION = 2\n", encoding="utf-8")
        expected = sha256_path(source)
        self.freeze()
        document = json.loads((self.contract / "m2_run_manifest.json").read_text())
        entry = next(row for row in document["source_snapshot"]["entries"] if row["role"] == role)
        self.assertEqual(entry["source"]["sha256"], expected)

    def test_exclusive_verification_receipt(self) -> None:
        self.freeze()
        output = self.base / "verification.json"
        digest = freezer._write_authenticated_exclusive(output, self.verify())
        self.assertEqual(digest, sha256_path(output))
        self.assertEqual(output.stat().st_nlink, 1)
        self.assertEqual(output.with_name(f"{output.name}.sha256").stat().st_nlink, 1)
        with self.assertRaisesRegex(freezer.ReleaseContractError, "overwrite"):
            freezer._write_authenticated_exclusive(output, self.verify())

    def test_second_receipt_publish_failure_archives_owned_half_pair(self) -> None:
        self.freeze()
        output = self.base / "verification.json"
        sidecar = output.with_name(f"{output.name}.sha256")
        original = freezer._rename_noreplace

        def fail_second_publish(source: pathlib.Path, destination: pathlib.Path) -> None:
            if destination == sidecar:
                raise freezer.ReleaseContractError("injected second publish failure")
            original(source, destination)

        with mock.patch.object(freezer, "_rename_noreplace", side_effect=fail_second_publish):
            with self.assertRaisesRegex(freezer.ReleaseContractError, "second publish failure"):
                freezer._write_authenticated_exclusive(output, self.verify())

        self.assertFalse(output.exists(), "half-published authority receipt must be archived")
        self.assertFalse(sidecar.exists(), "failed authority sidecar must not appear")
        archives = list(self.base.glob(".failed-publication.*"))
        self.assertEqual(len(archives), 1)
        archived_files = [path for path in archives[0].iterdir() if path.is_file()]
        self.assertEqual(len(archived_files), 2)
        self.assertTrue(all(path.stat().st_nlink == 1 for path in archived_files))
        self.assertTrue(all(stat.S_IMODE(path.stat().st_mode) == 0o444 for path in archived_files))

    def test_freeze_publish_failure_archives_complete_staging_tree(self) -> None:
        original = freezer._rename_noreplace

        def fail_contract_publish(source: pathlib.Path, destination: pathlib.Path) -> None:
            if destination == self.contract:
                raise freezer.ReleaseContractError("injected freeze publish failure")
            original(source, destination)

        with mock.patch.object(freezer, "_rename_noreplace", side_effect=fail_contract_publish):
            with self.assertRaisesRegex(freezer.ReleaseContractError, "freeze publish failure"):
                self.freeze()

        self.assertFalse(self.contract.exists())
        archives = list(self.base.glob(".failed-staging.*"))
        self.assertEqual(len(archives), 1)
        archived_manifest = archives[0] / "m2_run_manifest.json"
        self.assertTrue(archived_manifest.is_file())
        self.assertTrue(archived_manifest.with_name(f"{archived_manifest.name}.sha256").is_file())
        self.assertEqual(archived_manifest.stat().st_nlink, 1)


class PreAndFreshVerifierRedTeamTests(SyntheticReleaseFixture):
    def test_forged_structurally_complete_pre_is_rejected(self) -> None:
        document = json.loads(self.pre.read_text())
        document["artifact_audit"]["artifacts"][0]["observed"]["st_ino"] += 1
        document["input_identity_snapshot"]["artifacts"][0]["observed"]["st_ino"] += 1
        document["input_identity_snapshot_sha256"] = freezer.canonical_json_sha256(document["input_identity_snapshot"])
        write_authenticated_json(self.pre, document)
        with self.assertRaisesRegex(freezer.ReleaseContractError, "manifest-derived"):
            self.freeze()

    def test_wrong_or_extra_pre_check_is_rejected(self) -> None:
        document = json.loads(self.pre.read_text())
        document["checks"]["forged_check"] = True
        write_authenticated_json(self.pre, document)
        with self.assertRaisesRegex(freezer.ReleaseContractError, "checks schema mismatch"):
            self.freeze()

    def test_synthetic_freeze_cannot_skip_fresh_verifier(self) -> None:
        with self.assertRaisesRegex(freezer.ReleaseContractError, "must inject.*fresh"):
            freezer.freeze_release_contract(
                self.manifest, self.pre, self.contract, _repo_root=self.repo,
                _authority=self.authority, _runtime=SYNTHETIC_RUNTIME,
            )

    def test_pre_sidecar_tamper_is_rejected(self) -> None:
        self.pre.with_name(f"{self.pre.name}.sha256").write_text(f"{'0' * 64}  {self.pre.name}\n")
        with self.assertRaisesRegex(freezer.ReleaseContractError, "sidecar mismatch"):
            self.freeze()

    def test_immutable_manifest_and_pre_copies_are_used_after_origins_removed(self) -> None:
        self.freeze()
        self.manifest.rename(self.manifest.with_name(f"{self.manifest.name}.removed"))
        self.pre.rename(self.pre.with_name(f"{self.pre.name}.removed"))
        pre_sidecar = self.pre.with_name(f"{self.pre.name}.sha256")
        pre_sidecar.rename(pre_sidecar.with_name(f"{pre_sidecar.name}.removed"))
        self.assertTrue(self.verify()["all_pass"])


class SourcePathAndHardlinkRedTeamTests(SyntheticReleaseFixture):
    def test_publish_boundary_source_mutation_is_rejected(self) -> None:
        _, relative = freezer.SOURCE_ALLOWLIST[0]
        target = self.repo / relative
        original = freezer._stable_regular_file
        mutated = False

        def adversarial(path, label):
            nonlocal mutated
            if label == "publish-boundary source/extractor" and not mutated:
                target.write_text("# mutated after initial snapshot\n", encoding="utf-8")
                mutated = True
            return original(path, label)

        with mock.patch.object(freezer, "_stable_regular_file", side_effect=adversarial):
            with self.assertRaisesRegex(freezer.ReleaseContractError, "drifted before publish"):
                self.freeze()
        self.assertTrue(mutated)
        self.assertFalse(self.contract.exists())

    def test_source_file_symlink_is_rejected(self) -> None:
        _, relative = freezer.SOURCE_ALLOWLIST[0]
        source = self.repo / relative
        target = self.base / "outside.py"
        target.write_text("# outside\n")
        source.rename(source.with_name(f"{source.name}.original"))
        source.symlink_to(target)
        with self.assertRaisesRegex(freezer.ReleaseContractError, "regular|symlink"):
            self.freeze()

    def test_repo_intermediate_parent_symlink_escape_is_rejected(self) -> None:
        scripts = self.repo / freezer.SCRIPTS_RELATIVE
        outside = self.base / "outside_scripts"
        shutil.move(str(scripts), str(outside))
        scripts.symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(freezer.ReleaseContractError, "parent component"):
            self.freeze()

    def test_external_hardlink_to_source_is_rejected(self) -> None:
        _, relative = freezer.SOURCE_ALLOWLIST[0]
        os.link(self.repo / relative, self.base / "external_source_hardlink")
        with self.assertRaisesRegex(freezer.ReleaseContractError, "hardlink count"):
            self.freeze()

    def test_external_hardlink_to_snapshot_is_rejected(self) -> None:
        self.freeze()
        document = json.loads((self.contract / "m2_run_manifest.json").read_text())
        copy = self.contract / document["source_snapshot"]["entries"][0]["snapshot"]["path"]
        os.link(copy, self.base / "external_snapshot_hardlink")
        with self.assertRaisesRegex(freezer.ReleaseContractError, "hardlink count|st_nlink"):
            self.verify()


class RunManifestRewriteRedTeamTests(SyntheticReleaseFixture):
    def test_runtime_rewrite_and_recomputed_digests_is_rejected(self) -> None:
        self.freeze()
        def mutate(document: dict) -> None:
            document["runtime"]["platform"] = "forged-runtime"
            document["runtime_semantic_sha256"] = freezer.canonical_json_sha256(document["runtime"])
        self.rewrite_run_manifest(mutate)
        with self.assertRaisesRegex(freezer.ReleaseContractError, "current runtime differs"):
            self.verify()

    def test_chromosome_names_rewrite_is_rejected(self) -> None:
        self.freeze()
        self.rewrite_run_manifest(lambda document: document["scope"].update({"chromosome_names": ["chr1"] * 22}))
        with self.assertRaisesRegex(freezer.ReleaseContractError, "exact scope"):
            self.verify()

    def test_producer_rewrite_is_rejected(self) -> None:
        self.freeze()
        self.rewrite_run_manifest(lambda document: document["producer"].update({"role": "impostor"}))
        with self.assertRaisesRegex(freezer.ReleaseContractError, "producer"):
            self.verify()

    def test_fresh_verifier_binding_rewrite_is_rejected(self) -> None:
        self.freeze()
        self.rewrite_run_manifest(
            lambda document: document["fresh_input_identity_verification"].update(
                {"verifier_path": "/forged/verifier.py"}
            )
        )
        with self.assertRaisesRegex(freezer.ReleaseContractError, "not bound"):
            self.verify()

    def test_production_verify_rejects_synthetic_contract(self) -> None:
        self.freeze()
        with self.assertRaisesRegex(freezer.ReleaseContractError, "synthetic.*forbidden"):
            freezer.verify_release_contract(self.contract)


if __name__ == "__main__":
    unittest.main()
