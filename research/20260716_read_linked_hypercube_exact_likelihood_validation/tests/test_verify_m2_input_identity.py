from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import stat
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "verify_m2_input_identity.py"
)
CANONICAL_MANIFEST_TEMPLATE = (
    Path(__file__).resolve().parents[3]
    / "research"
    / "20260710_layered_reconstruction_v2"
    / "data"
    / "layered_input_manifest_v3_raw_all_lps_pass_v2.json"
)
SPEC = importlib.util.spec_from_file_location("verify_m2_input_identity", SCRIPT)
assert SPEC and SPEC.loader
identity = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(identity)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def full(path: Path) -> dict:
    return {"policy": "full_sha256", "size_bytes": path.stat().st_size, "sha256": sha(path)}


def sampled(path: Path, label: str, offset: int, length: int) -> dict:
    payload = path.read_bytes()[offset : offset + length]
    return {
        "label": label,
        "offset": offset,
        "length": length,
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def storage(bam: Path, bai: Path) -> dict:
    logical = bam.lstat()
    resolved = bam.resolve(strict=True)
    target = resolved.stat()
    chunk_size = 1024 * 1024
    length = min(chunk_size, target.st_size)
    middle = max(0, (target.st_size - length) // 2)
    last = max(0, target.st_size - length)
    value = {
        "policy": "storage_identity_v1",
        "assurance": "metadata_plus_sampled_chunks_not_full_content_hash",
        "is_full_content_hash": False,
        "requested_path": str(bam),
        "realpath": str(resolved),
        "logical_is_symlink": stat.S_ISLNK(logical.st_mode),
        "logical_size_bytes": logical.st_size,
        "logical_mtime_ns": logical.st_mtime_ns,
        "st_dev": target.st_dev,
        "st_ino": target.st_ino,
        "size_bytes": target.st_size,
        "mtime_ns": target.st_mtime_ns,
        "ctime_ns": target.st_ctime_ns,
        "chunk_size_bytes": chunk_size,
        "chunks": [
            sampled(resolved, "first", 0, length),
            sampled(resolved, "middle", middle, length),
            sampled(resolved, "last", last, length),
        ],
        "index": {"path": str(bai), "identity": full(bai)},
    }
    value["identity_sha256"] = identity.canonical_sha256(value)
    return value


class Fixture:
    def __init__(self, root: Path):
        self.root = root
        self.samples: dict[str, dict[str, Path]] = {}
        document = json.loads(CANONICAL_MANIFEST_TEMPLATE.read_text())
        document["manifest_id"] = f"synthetic_m2_identity_{root.name}"
        document["created_at_utc"] = "2026-07-16T00:00:00Z"
        by_sample = {row["sample"]: row for row in document["samples"]}
        for index, dataset in enumerate(identity.DATASETS):
            directory = root / dataset
            directory.mkdir(parents=True)
            bam = directory / f"{dataset}.bam"
            bai = directory / f"{dataset}.bam.bai"
            vcf = directory / f"{dataset}.vcf.gz"
            csi = directory / f"{dataset}.vcf.gz.csi"
            tags = directory / f"{dataset}.read_tags.tsv.gz"
            tbi = directory / f"{dataset}.read_tags.tsv.gz.tbi"
            bam.write_bytes(bytes((position + index) % 251 for position in range(8192)))
            bai.write_bytes(f"bai-{dataset}".encode())
            vcf.write_bytes(f"vcf-{dataset}".encode())
            csi.write_bytes(f"csi-{dataset}".encode())
            tags.write_bytes(f"tags-{dataset}".encode())
            tbi.write_bytes(f"tbi-{dataset}".encode())
            self.samples[dataset] = {
                "bam": bam, "bai": bai, "vcf": vcf, "csi": csi, "tags": tags, "tbi": tbi
            }
            row = by_sample[dataset]
            row["alignment_payload"] = {
                "path": str(bam),
                "index_path": str(bai),
                "embedded_hp_ps_policy": "ignore",
                "identity_policy": "storage_identity_v1",
                "storage_identity_v1": storage(bam, bai),
            }
            row["somatic"]["tree_vcf"] = {
                "path": str(vcf),
                "identity": full(vcf),
                "index": {"path": str(csi), "identity": full(csi)},
            }
            row["read_tags"]["sidecar"] = {"path": str(tags), "identity": full(tags)}
            row["read_tags"]["index"] = {"path": str(tbi), "identity": full(tbi)}
            row["read_tags"]["subject_binding"]["sidecar_sha256"] = full(tags)["sha256"]
            row["read_tags"]["subject_binding"]["sidecar_index_sha256"] = full(tbi)["sha256"]
            row["read_tags"]["subject_binding"]["alignment_payload_storage_identity_sha256"] = (
                row["alignment_payload"]["storage_identity_v1"]["identity_sha256"]
            )
        self.manifest = root / "input_manifest.snapshot.json"
        self.manifest.write_text(json.dumps(document, indent=2) + "\n")


class IdentityVerificationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def verify_fixture(self, compare_to: Path | None = None) -> dict:
        return identity.verify(
            self.fixture.manifest,
            compare_to,
            schema_path=identity.CANONICAL_SCHEMA_PATH,
            expected_manifest_sha256=sha(self.fixture.manifest),
            expected_schema_sha256=identity.CANONICAL_SCHEMA_SHA256,
            _test_authority_override=True,
        )

    def run_cli(self, output: Path, compare_to: Path | None = None) -> int:
        arguments = [
            "--manifest", str(self.fixture.manifest),
            "--output", str(output),
        ]
        if compare_to is not None:
            arguments.extend(("--compare-to", str(compare_to)))
        return identity.main(arguments)

    def test_positive_scope_counts_and_bam_is_never_full_hashed(self) -> None:
        original = identity._sha256_fd
        hashed_inodes = []

        def recording(descriptor: int) -> str:
            metadata = os.fstat(descriptor)
            hashed_inodes.append((metadata.st_dev, metadata.st_ino))
            return original(descriptor)

        with mock.patch.object(identity, "_sha256_fd", side_effect=recording):
            receipt = self.verify_fixture()
        self.assertTrue(receipt["all_pass"])
        audit = receipt["artifact_audit"]
        self.assertEqual(audit["n_artifacts"], 42)
        self.assertEqual(audit["n_unique_resolved_files"], 42)
        self.assertEqual(audit["n_storage_identity_v1"], 7)
        self.assertEqual(audit["n_full_sha256"], 35)
        self.assertEqual(audit["n_sampled_bam_chunks"], 21)
        bam_inodes = {
            (paths["bam"].stat().st_dev, paths["bam"].stat().st_ino)
            for paths in self.fixture.samples.values()
        }
        self.assertTrue(bam_inodes.isdisjoint(hashed_inodes))
        self.assertFalse(receipt["assurance"]["temporal_immutability_proven"])

    def test_pre_post_test_receipts_are_authenticated_and_equal_but_not_evidence(self) -> None:
        pre = self.root / "pre.json"
        post = self.root / "post.json"
        identity.write_receipt(pre, self.verify_fixture())
        identity.write_receipt(post, self.verify_fixture(pre))
        pre_doc = json.loads(pre.read_text())
        post_doc = json.loads(post.read_text())
        self.assertEqual(pre_doc["mode"], "PRE")
        self.assertEqual(post_doc["mode"], "POST_COMPARE")
        self.assertTrue(post_doc["compare_to"]["exact_snapshot_equal"])
        self.assertFalse(pre_doc["validation_evidence_eligible"])
        self.assertEqual(pre_doc["authority_mode"], "TEST_ONLY_UNFROZEN")
        self.assertEqual(pre_doc["input_identity_snapshot"], post_doc["input_identity_snapshot"])
        self.assertEqual(pre.with_name("pre.json.sha256").read_text().split()[0], sha(pre))
        self.assertEqual(post.with_name("post.json.sha256").read_text().split()[0], sha(post))

    def test_full_sha_artifact_tamper_is_detected(self) -> None:
        target = self.fixture.samples[identity.DATASETS[0]]["vcf"]
        target.write_bytes(b"X" + target.read_bytes()[1:])
        receipt = self.verify_fixture()
        self.assertFalse(receipt["all_pass"])
        mismatch = [
            row for row in receipt["artifact_audit"]["artifacts"]
            if row["dataset"] == identity.DATASETS[0] and row["role"] == "tree_vcf"
        ]
        self.assertEqual(len(mismatch), 1)
        self.assertFalse(mismatch[0]["match"])

    def test_bam_chunk_tamper_is_detected_even_at_constant_size(self) -> None:
        target = self.fixture.samples[identity.DATASETS[1]]["bam"]
        before = target.stat()
        payload = bytearray(target.read_bytes())
        payload[10] ^= 1
        target.write_bytes(payload)
        os.utime(target, ns=(before.st_atime_ns, before.st_mtime_ns))
        receipt = self.verify_fixture()
        self.assertFalse(receipt["all_pass"])
        row = next(
            row for row in receipt["artifact_audit"]["artifacts"]
            if row["dataset"] == identity.DATASETS[1] and row["role"] == "alignment_bam"
        )
        self.assertFalse(row["match"])
        self.assertNotEqual(
            row["expected"]["chunks"][0]["sha256"], row["observed"]["chunks"][0]["sha256"]
        )

    def test_bam_metadata_tamper_is_detected(self) -> None:
        target = self.fixture.samples[identity.DATASETS[2]]["bam"]
        before = target.stat()
        os.utime(target, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
        receipt = self.verify_fixture()
        self.assertFalse(receipt["all_pass"])
        row = next(
            row for row in receipt["artifact_audit"]["artifacts"]
            if row["dataset"] == identity.DATASETS[2] and row["role"] == "alignment_bam"
        )
        self.assertNotEqual(row["expected"]["mtime_ns"], row["observed"]["mtime_ns"])

    def test_manifest_scope_tamper_fails_closed(self) -> None:
        document = json.loads(self.fixture.manifest.read_text())
        document["samples"].pop()
        document["dataset_count"] = 6
        self.fixture.manifest.write_text(json.dumps(document))
        with self.assertRaises(identity.IdentityError):
            self.verify_fixture()

    def test_schema_and_integer_types_are_fail_closed(self) -> None:
        document = json.loads(self.fixture.manifest.read_text())
        document["schema_name"] = "synthetic.canonical"
        self.fixture.manifest.write_text(json.dumps(document))
        with self.assertRaises(identity.IdentityError):
            self.verify_fixture()

        self.fixture = Fixture(self.root / "second")
        document = json.loads(self.fixture.manifest.read_text())
        document["dataset_count"] = 7.0
        self.fixture.manifest.write_text(json.dumps(document))
        with self.assertRaises(identity.IdentityError):
            self.verify_fixture()

    def test_full_identity_float_and_zero_length_chunk_are_rejected(self) -> None:
        document = json.loads(self.fixture.manifest.read_text())
        document["samples"][0]["read_tags"]["sidecar"]["identity"]["size_bytes"] = float(
            document["samples"][0]["read_tags"]["sidecar"]["identity"]["size_bytes"]
        )
        self.fixture.manifest.write_text(json.dumps(document))
        with self.assertRaises(identity.IdentityError):
            self.verify_fixture()

        self.fixture = Fixture(self.root / "second")
        document = json.loads(self.fixture.manifest.read_text())
        storage_block = document["samples"][0]["alignment_payload"]["storage_identity_v1"]
        storage_block["chunk_size_bytes"] = 0
        for chunk in storage_block["chunks"]:
            chunk["length"] = 0
            chunk["offset"] = 0
            chunk["label"] = "same"
            chunk["sha256"] = hashlib.sha256(b"").hexdigest()
        without_digest = dict(storage_block)
        without_digest.pop("identity_sha256")
        storage_block["identity_sha256"] = identity.canonical_sha256(without_digest)
        self.fixture.manifest.write_text(json.dumps(document))
        with self.assertRaises(identity.IdentityError):
            self.verify_fixture()

    def test_role_alias_to_bam_is_rejected_before_any_full_bam_hash(self) -> None:
        document = json.loads(self.fixture.manifest.read_text())
        sample = document["samples"][0]
        bam = Path(sample["alignment_payload"]["path"])
        sample["somatic"]["tree_vcf"]["path"] = str(bam)
        sample["somatic"]["tree_vcf"]["identity"] = full(bam)
        self.fixture.manifest.write_text(json.dumps(document))
        original = identity._sha256_fd
        bam_inode = (bam.stat().st_dev, bam.stat().st_ino)
        hashed_inodes = []

        def rejecting_bam(descriptor: int) -> str:
            metadata = os.fstat(descriptor)
            current = (metadata.st_dev, metadata.st_ino)
            hashed_inodes.append(current)
            if current == bam_inode:
                raise AssertionError("verifier attempted a forbidden full BAM hash")
            return original(descriptor)

        with mock.patch.object(identity, "_sha256_fd", side_effect=rejecting_bam):
            with self.assertRaises(identity.IdentityError):
                self.verify_fixture()
        self.assertNotIn(bam_inode, hashed_inodes)

    def test_symlink_swap_after_inventory_is_rejected_before_bam_read(self) -> None:
        document = json.loads(self.fixture.manifest.read_text())
        sample = next(row for row in document["samples"] if row["sample"] == identity.DATASETS[0])
        tree_path = Path(sample["somatic"]["tree_vcf"]["path"])
        tree_target = tree_path.with_name(f"{tree_path.name}.payload")
        tree_path.rename(tree_target)
        tree_path.symlink_to(tree_target)
        sample["somatic"]["tree_vcf"]["identity"] = full(tree_target)
        self.fixture.manifest.write_text(json.dumps(document))
        bam = self.fixture.samples[identity.DATASETS[0]]["bam"]
        bam_inode = (bam.stat().st_dev, bam.stat().st_ino)
        original_inventory = identity._build_role_inventory
        original_hash = identity._sha256_fd
        hashed_inodes = []

        def swap_after_inventory(by_dataset: dict) -> list[dict]:
            inventory = original_inventory(by_dataset)
            tree_path.unlink()
            tree_path.symlink_to(bam)
            return inventory

        def recording_hash(descriptor: int) -> str:
            metadata = os.fstat(descriptor)
            current = (metadata.st_dev, metadata.st_ino)
            hashed_inodes.append(current)
            if current == bam_inode:
                raise AssertionError("TOCTOU caused a forbidden full BAM hash")
            return original_hash(descriptor)

        with mock.patch.object(identity, "_build_role_inventory", side_effect=swap_after_inventory):
            with mock.patch.object(identity, "_sha256_fd", side_effect=recording_hash):
                with self.assertRaises(identity.IdentityError):
                    self.verify_fixture()
        self.assertNotIn(bam_inode, hashed_inodes)

    def test_symlink_swap_after_verified_open_is_rejected_at_boundary(self) -> None:
        document = json.loads(self.fixture.manifest.read_text())
        sample = next(row for row in document["samples"] if row["sample"] == identity.DATASETS[0])
        tree_path = Path(sample["somatic"]["tree_vcf"]["path"])
        tree_target = tree_path.with_name(f"{tree_path.name}.payload")
        replacement = tree_path.with_name(f"{tree_path.name}.replacement")
        tree_path.rename(tree_target)
        replacement.write_bytes(b"replacement-target")
        tree_path.symlink_to(tree_target)
        sample["somatic"]["tree_vcf"]["identity"] = full(tree_target)
        self.fixture.manifest.write_text(json.dumps(document))
        target_inode = (tree_target.stat().st_dev, tree_target.stat().st_ino)
        original_hash = identity._sha256_fd
        swapped = False

        def swap_during_hash(descriptor: int) -> str:
            nonlocal swapped
            metadata = os.fstat(descriptor)
            if not swapped and (metadata.st_dev, metadata.st_ino) == target_inode:
                tree_path.unlink()
                tree_path.symlink_to(replacement)
                swapped = True
            return original_hash(descriptor)

        with mock.patch.object(identity, "_sha256_fd", side_effect=swap_during_hash):
            with self.assertRaises(identity.IdentityError):
                self.verify_fixture()
        self.assertTrue(swapped)

    def test_complete_schema_rejects_missing_required_and_unknown_top_level(self) -> None:
        document = json.loads(self.fixture.manifest.read_text())
        document.pop("manifest_id")
        self.fixture.manifest.write_text(json.dumps(document))
        with self.assertRaises(identity.IdentityError):
            self.verify_fixture()

        self.fixture = Fixture(self.root / "second")
        document = json.loads(self.fixture.manifest.read_text())
        document["unexpected_top_level"] = True
        self.fixture.manifest.write_text(json.dumps(document))
        with self.assertRaises(identity.IdentityError):
            self.verify_fixture()

    def test_manifest_and_schema_sha_locks_are_fail_closed(self) -> None:
        with self.assertRaises(identity.IdentityError):
            identity.verify(
                self.fixture.manifest,
                schema_path=identity.CANONICAL_SCHEMA_PATH,
                expected_manifest_sha256="0" * 64,
                expected_schema_sha256=identity.CANONICAL_SCHEMA_SHA256,
                _test_authority_override=True,
            )
        with self.assertRaises(identity.IdentityError):
            identity.verify(
                self.fixture.manifest,
                schema_path=identity.CANONICAL_SCHEMA_PATH,
                expected_manifest_sha256=sha(self.fixture.manifest),
                expected_schema_sha256="0" * 64,
                _test_authority_override=True,
            )

    def test_compare_detects_manifest_drift_even_if_input_paths_do_not_change(self) -> None:
        pre = self.root / "pre.json"
        identity.write_receipt(pre, self.verify_fixture())
        document = json.loads(self.fixture.manifest.read_text())
        document["manifest_id"] = "synthetic_m2_identity_changed"
        self.fixture.manifest.write_text(json.dumps(document))
        receipt = self.verify_fixture(pre)
        self.assertFalse(receipt["all_pass"])
        self.assertFalse(receipt["compare_to"]["exact_snapshot_equal"])

    def test_compare_receipt_or_sidecar_tamper_is_rejected(self) -> None:
        pre = self.root / "pre.json"
        identity.write_receipt(pre, self.verify_fixture())
        pre.write_text(pre.read_text() + " ")
        with self.assertRaises(identity.IdentityError):
            self.verify_fixture(pre)

    def test_compare_receipt_requires_exact_boolean_authority_field(self) -> None:
        pre = self.root / "pre.json"
        identity.write_receipt(pre, self.verify_fixture())
        document = json.loads(pre.read_text())
        document["validation_evidence_eligible"] = ""
        pre.write_text(json.dumps(document))
        pre.with_name("pre.json.sha256").write_text(
            f"{sha(pre)}  {pre.name}\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(identity.IdentityError, "is not a boolean"):
            self.verify_fixture(pre)

    def test_full_identity_cache_never_skips_rehash_or_boundary_checks(self) -> None:
        path = self.fixture.samples[identity.DATASETS[0]]["bai"]
        dev_ino = (path.stat().st_dev, path.stat().st_ino)
        cache = {}
        original = identity._sha256_fd
        with mock.patch.object(identity, "_sha256_fd", wraps=original) as hasher:
            first = identity.full_identity(
                path,
                expected_dev_ino=dev_ino,
                forbidden_dev_inos=frozenset(),
                cache=cache,
            )
            second = identity.full_identity(
                path,
                expected_dev_ino=dev_ino,
                forbidden_dev_inos=frozenset(),
                cache=cache,
            )
        self.assertEqual(first, second)
        self.assertEqual(hasher.call_count, 2)

    def test_output_overwrite_is_rejected(self) -> None:
        output = self.root / "pre.json"
        receipt = self.verify_fixture()
        identity.write_receipt(output, receipt)
        with self.assertRaises(identity.IdentityError):
            identity.write_receipt(output, receipt)

    def test_production_cli_rejects_synthetic_authority_and_override_flags(self) -> None:
        output = self.root / "production_attempt.json"
        self.assertEqual(self.run_cli(output), 1)
        receipt = json.loads(output.read_text())
        self.assertFalse(receipt["all_pass"])
        self.assertIn("canonical v5 manifest SHA-256 mismatch", receipt["failure"])

        with self.assertRaises(SystemExit):
            identity.parse_args([
                "--manifest", str(self.fixture.manifest),
                "--output", str(self.root / "override.json"),
                "--expected-manifest-sha256", sha(self.fixture.manifest),
            ])


if __name__ == "__main__":
    unittest.main()
