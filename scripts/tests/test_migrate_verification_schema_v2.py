#!/usr/bin/env python3
"""Tests for the byte-preserving Verification schema v2 migration.

All CLI tests use TemporaryDirectory.  They never create or update the
authoritative output/schema_migrations/20260714_verification_v2 directory.
"""

from __future__ import annotations

import csv
import errno
import hashlib
import importlib.util
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[2]
MIGRATOR_PATH = REPO / "scripts" / "migration" / "migrate_verification_schema_v2.py"

MODULE_SPEC = importlib.util.spec_from_file_location("verification_schema_v2_migrator", MIGRATOR_PATH)
if MODULE_SPEC is None or MODULE_SPEC.loader is None:
    raise RuntimeError("cannot load migration module")
MIGRATOR = importlib.util.module_from_spec(MODULE_SPEC)
sys.modules[MODULE_SPEC.name] = MIGRATOR
MODULE_SPEC.loader.exec_module(MIGRATOR)


BASE_HEADERS = (
    "RegionID",
    "Chr",
    "Pos",
    "Ref",
    "Alt",
    "VerificationClass",
    "VerificationClass_Legacy",
    "Potential_LOH",
    "LOH_Subtype",
    "Significant",
    "Note",
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def tsv_rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class MigrationFixture:
    def __init__(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="verification-v2-migration-test-")
        self.root = Path(self.temporary.name)
        self.input_root = self.root / "input"
        self.input_root.mkdir()
        self.manifest = self.root / "manifest.tsv"

    def close(self):
        self.temporary.cleanup()

    def add_input(self, filename: str, data: bytes) -> Path:
        path = self.input_root / filename
        path.write_bytes(data)
        return path

    def write_manifest(self, entries, include_counts=False):
        headers = ["file", "rows", "sha256"]
        if include_counts:
            headers.extend(MIGRATOR.MANIFEST_COUNT_HEADERS)
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for entry in entries:
            writer.writerow({name: entry.get(name, "") for name in headers})
        with self.manifest.open("w", encoding="utf-8", newline="") as handle:
            handle.write(stream.getvalue())

    def manifest_entry(self, filename: str, data: bytes, **extra):
        entry = {
            "file": filename,
            "rows": MIGRATOR.count_csv_data_rows(data),
            "sha256": digest(data),
        }
        entry.update(extra)
        return entry

    def run(self, output: Path, extra_args=None, include_absence_flag=True):
        command = [
            sys.executable,
            str(MIGRATOR_PATH),
            "--manifest",
            str(self.manifest),
            "--input-root",
            str(self.input_root),
            "--output-dir",
            str(output),
        ]
        if include_absence_flag:
            command.append("--require-output-dir-absent")
        if extra_args:
            command.extend(extra_args)
        return subprocess.run(command, cwd=REPO, capture_output=True, text=True)


def simple_csv(
    current="Strong",
    legacy="Strong",
    potential_loh="true",
    loh_subtype="LOH_Strong",
    significant="true",
    headers=BASE_HEADERS,
):
    values = {
        "RegionID": "1",
        "Chr": "chr1",
        "Pos": "101",
        "Ref": "A",
        "Alt": "T",
        "VerificationClass": current,
        "VerificationClass_Legacy": legacy,
        "Potential_LOH": potential_loh,
        "LOH_Subtype": loh_subtype,
        "Significant": significant,
        "Note": "plain",
    }
    header = ",".join(headers).encode("ascii") + b"\n"
    row = ",".join(values.get(name, "") for name in headers).encode("ascii") + b"\n"
    return header + row


class RawCsvParserTests(unittest.TestCase):
    def test_parser_preserves_crlf_quotes_commas_and_embedded_newline(self):
        data = (
            b'a,b,c\r\n'
            b'1,"quoted,comma","said ""hello"""\r\n'
            b'2,"line1\r\nline2",tail'
        )
        records = list(MIGRATOR.iter_csv_records(data))
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0].newline, b"\r\n")
        self.assertEqual(records[1].fields[1], b'"quoted,comma"')
        self.assertEqual(records[1].fields[2], b'"said ""hello"""')
        self.assertEqual(MIGRATOR.decode_csv_token(records[1].fields[2]), 'said "hello"')
        self.assertEqual(records[2].fields[1], b'"line1\r\nline2"')
        self.assertEqual(records[2].newline, b"")
        self.assertEqual(MIGRATOR.count_csv_data_rows(data), 2)

    def test_parser_rejects_unterminated_quote(self):
        with self.assertRaises(MIGRATOR.CsvFormatError):
            list(MIGRATOR.iter_csv_records(b'a,b\n1,"broken\n'))

    def test_post_write_audit_rejects_each_tampered_appended_contract_field(self):
        source = simple_csv()
        source_records = list(MIGRATOR.iter_csv_records(source))
        names = MIGRATOR.header_names(source_records[0])
        indices = {name: index for index, name in enumerate(names)}
        migrated = b"".join(
            MIGRATOR._migrated_record(record, indices, is_header=(index == 0))
            for index, record in enumerate(source_records)
        )
        self.assertEqual(MIGRATOR.audit_raw_preservation(source, migrated), 1)

        migrated_records = list(MIGRATOR.iter_csv_records(migrated))
        base_count = len(source_records[0].fields)
        tamper_cases = {
            0: b"999",
            2: b"false",
            3: b"false",
            4: b"true",
            5: b"false",
            6: b"WRONG_PATH",
            7: b"LIVE",
        }
        for appended_index, bad_token in tamper_cases.items():
            with self.subTest(field=MIGRATOR.APPENDED_HEADERS[appended_index]):
                fields = list(migrated_records[1].fields)
                fields[base_count + appended_index] = bad_token
                tampered = (
                    b",".join(migrated_records[0].fields)
                    + migrated_records[0].newline
                    + b",".join(fields)
                    + migrated_records[1].newline
                )
                with self.assertRaisesRegex(
                    MIGRATOR.MigrationError,
                    MIGRATOR.APPENDED_HEADERS[appended_index],
                ):
                    MIGRATOR.audit_raw_preservation(source, tampered)

        significant_fields = list(migrated_records[1].fields)
        significant_fields[indices["Significant"]] = b"false"
        significant_tamper = (
            b",".join(migrated_records[0].fields)
            + migrated_records[0].newline
            + b",".join(significant_fields)
            + migrated_records[1].newline
        )
        with self.assertRaisesRegex(
            MIGRATOR.MigrationError,
            "Significant stable-key invariant",
        ):
            MIGRATOR.audit_raw_preservation(source, significant_tamper)


class TruthTableTests(unittest.TestCase):
    def test_authoritative_v1_to_v2_truth_table(self):
        cases = (
            ("Strong", "Strong", "Strong_Bidirectional", "BIDIRECTIONAL", "true", "true"),
            ("Strong", "Subclone", "ClusterFirstOnly", "CLUSTER_FIRST_ONLY", "false", "true"),
            ("LOH-Structure", "Weak", "LOH-Structure", "LOH_STRUCTURE", "true", "false"),
            ("MultiGroupNoLabel", "Noise", "MultiGroupNoLabel", "WITHIN_HP_MULTIGROUP", "false", "false"),
            ("LabelShift", "Weak", "LabelShift", "LABEL_SHIFT", "true", "false"),
            ("PermanovaLocation", "Noise", "PermanovaLocation", "PERMANOVA_LOCATION", "false", "false"),
            ("StructureNoLabel", "Weak", "StructureNoLabel", "HP_AUC_STRUCTURE_NO_LABEL", "true", "false"),
            ("DispersionStructure", "Noise", "DispersionStructure", "DISPERSION_STRUCTURE", "false", "false"),
            ("Noise_Uniform", "Weak", "Noise_Uniform", "NOISE_UNIFORM", "true", "false"),
            ("Noise_Chaotic", "Noise", "Noise_Chaotic", "NOISE_CHAOTIC", "false", "false"),
            ("Noise_Uncorrelated", "Noise", "Noise_Uncorrelated", "NOISE_UNCORRELATED", "false", "false"),
        )
        for current, legacy, v2, path, label, cluster in cases:
            with self.subTest(current=current, legacy=legacy):
                decision = MIGRATOR.classify_offline(current, legacy)
                self.assertEqual(decision.verification_class_v2, v2)
                self.assertEqual(decision.evidence_path, path)
                self.assertEqual(decision.label_first_support, label)
                self.assertEqual(decision.cluster_first_support, cluster)

    def test_incompatible_and_unknown_pairs_fail(self):
        cases = (
            ("Strong", "Weak"),
            ("Strong", "Noise"),
            ("Noise_Uniform", "Strong"),
            ("Noise_Uniform", "Subclone"),
            ("Unknown", "Noise"),
            ("Strong", "Unknown"),
        )
        for current, legacy in cases:
            with self.subTest(current=current, legacy=legacy):
                with self.assertRaises(MIGRATOR.RowConflict):
                    MIGRATOR.classify_offline(current, legacy)

    def test_frozen_weighted_mapping_logic(self):
        frozen_weights = {
            ("Strong", "Strong"): 59910,
            ("Strong", "Subclone"): 6931,
        }
        mapped = Counter()
        paths = Counter()
        for pair, count in frozen_weights.items():
            decision = MIGRATOR.classify_offline(*pair)
            mapped[decision.verification_class_v2] += count
            paths[decision.evidence_path] += count
        self.assertEqual(mapped["Strong_Bidirectional"], 59910)
        self.assertEqual(mapped["ClusterFirstOnly"], 6931)
        self.assertEqual(sum(mapped.values()), 66841)
        self.assertEqual(paths["BIDIRECTIONAL"] + paths["CLUSTER_FIRST_ONLY"], 66841)


class MigrationCliTests(unittest.TestCase):
    def setUp(self):
        self.fixture = MigrationFixture()

    def tearDown(self):
        self.fixture.close()

    def test_valid_complex_csv_preserves_raw_tokens_and_reports_provenance(self):
        header = b",".join(name.encode("ascii") for name in BASE_HEADERS) + b"\r\n"
        rows = (
            b'1,chr1,10,A,T,"Strong",Strong,true,"LOH_Strong",true,"quoted,comma"\r\n'
            b'2,chr1,20,C,G,Strong,Subclone,true,LOH_Subclone,true,"said ""hello"""\r\n'
            b'3,chr2,30,G,A,Noise_Uniform,Noise,false,None,false,"line1\r\nline2"\r\n'
            b'4,chr2,40,T,C,MultiGroupNoLabel,Weak,true,LOH_Weak,true,plain'
        )
        data = header + rows
        filename = "complex.csv"
        self.fixture.add_input(filename, data)
        entry = self.fixture.manifest_entry(
            filename,
            data,
            input_final_strong=2,
            legacy_strong_to_final_strong=1,
            legacy_subclone_to_final_strong=1,
            legacy_strong_or_subclone_exceptions=0,
        )
        self.fixture.write_manifest([entry], include_counts=True)
        output = self.fixture.root / "migrated"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertIn("status=VALID", process.stdout)
        migrated_path = output / filename
        self.assertTrue(migrated_path.is_file())
        migrated = migrated_path.read_bytes()

        source_records = list(MIGRATOR.iter_csv_records(data))
        output_records = list(MIGRATOR.iter_csv_records(migrated))
        self.assertEqual(len(source_records), len(output_records))
        base_count = len(BASE_HEADERS)
        self.assertEqual(output_records[0].fields[:base_count], source_records[0].fields)
        self.assertEqual(
            [MIGRATOR.decode_csv_token(token) for token in output_records[0].fields[base_count:]],
            list(MIGRATOR.APPENDED_HEADERS),
        )
        source_names = MIGRATOR.header_names(source_records[0])
        indices = {name: index for index, name in enumerate(source_names)}
        expected_classes = (
            "Strong_Bidirectional",
            "ClusterFirstOnly",
            "Noise_Uniform",
            "MultiGroupNoLabel",
        )
        for row_index, (source, target, expected_class) in enumerate(
            zip(source_records[1:], output_records[1:], expected_classes), start=2
        ):
            with self.subTest(row=row_index):
                self.assertEqual(source.newline, target.newline)
                for field_index, token in enumerate(source.fields):
                    if field_index != indices["VerificationClass"]:
                        self.assertEqual(target.fields[field_index], token)
                self.assertEqual(
                    MIGRATOR.decode_csv_token(target.fields[indices["VerificationClass"]]),
                    expected_class,
                )
                appended = target.fields[base_count:]
                self.assertEqual(appended[1], source.fields[indices["VerificationClass"]])
                self.assertEqual(appended[4], b"NA")
                self.assertEqual(appended[5], b"NA")
                self.assertEqual(appended[7], b"LEGACY_CLASS")
                self.assertEqual(appended[8], source.fields[indices["LOH_Subtype"]])

        self.assertTrue(output_records[1].fields[indices["VerificationClass"]].startswith(b'"'))
        self.assertEqual(output_records[-1].newline, b"")
        self.assertIn(b'"line1\r\nline2"', migrated)

        summary = json.loads((output / "migration_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["status"], "VALID")
        self.assertEqual(summary["input_rows"], 4)
        self.assertEqual(summary["output_rows"], 4)
        self.assertEqual(summary["after_current_counts"]["Strong_Bidirectional"], 1)
        self.assertEqual(summary["after_current_counts"]["ClusterFirstOnly"], 1)
        self.assertIn("--require-output-dir-absent", summary["command"])
        self.assertRegex(summary["generated_at"], r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")
        for audit_flag in (
            "raw_token_preservation",
            "significant_stable_key_invariant",
            "stable_key_uniqueness",
        ):
            self.assertEqual(summary[audit_flag], "PASS")
        self.assertEqual(len(tsv_rows(output / "migrated_outputs_manifest.tsv")), 1)
        file_report = tsv_rows(output / "migration_file_report.tsv")[0]
        self.assertEqual(file_report["input_sha256"], digest(data))
        self.assertEqual(file_report["output_sha256"], digest(migrated))
        status_report = tsv_rows(output / "migration_status.tsv")[0]
        self.assertEqual(status_report["status"], "VALID")
        for audit_flag in (
            "raw_token_preservation",
            "significant_stable_key_invariant",
            "stable_key_uniqueness",
        ):
            self.assertEqual(file_report[audit_flag], "PASS")
            self.assertEqual(status_report[audit_flag], "PASS")
        self.assertEqual(tsv_rows(output / "unmapped_conflicts.tsv"), [])

    def test_conflict_writes_diagnostics_and_does_not_publish_file(self):
        data = simple_csv(current="Strong", legacy="Weak", loh_subtype="LOH_Weak")
        filename = "conflict.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "conflict-output"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 1, process.stderr)
        self.assertFalse((output / filename).exists())
        conflicts = tsv_rows(output / "unmapped_conflicts.tsv")
        self.assertEqual(len(conflicts), 1)
        self.assertIn("INCOMPATIBLE_STRONG_LEGACY_PAIR", conflicts[0]["reason"])
        status = tsv_rows(output / "migration_status.tsv")[0]
        self.assertEqual(status["status"], "FAILED")
        self.assertEqual(status["failed_files"], "1")

    def test_loh_provenance_mismatch_fails_closed(self):
        data = simple_csv(
            current="Strong", legacy="Strong", potential_loh="true", loh_subtype="LOH_Noise"
        )
        filename = "loh-mismatch.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "loh-mismatch-output"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 1, process.stderr)
        self.assertFalse((output / filename).exists())
        conflicts = tsv_rows(output / "unmapped_conflicts.tsv")
        self.assertIn("LOH_PROVENANCE_MISMATCH", conflicts[0]["reason"])

    def test_missing_legacy_is_unresolved_and_no_canonical_is_published(self):
        headers = tuple(name for name in BASE_HEADERS if name != "VerificationClass_Legacy")
        data = simple_csv(headers=headers)
        filename = "missing-legacy.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "missing-legacy-output"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 1, process.stderr)
        self.assertFalse((output / filename).exists())
        unresolved = tsv_rows(output / "unresolved_files.tsv")
        self.assertEqual(unresolved[0]["status"], "UNRESOLVED_LEGACY_PROVENANCE")
        self.assertIn("MISSING_VERIFICATION_CLASS_LEGACY", unresolved[0]["reason"])
        report = tsv_rows(output / "migration_file_report.tsv")[0]
        self.assertEqual(report["input_rows"], "1")

    def test_missing_loh_provenance_is_unresolved(self):
        headers = tuple(name for name in BASE_HEADERS if name != "LOH_Subtype")
        data = simple_csv(headers=headers)
        filename = "missing-loh.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "missing-loh-output"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 1, process.stderr)
        self.assertFalse((output / filename).exists())
        unresolved = tsv_rows(output / "unresolved_files.tsv")
        self.assertEqual(unresolved[0]["status"], "UNRESOLVED_LOH_PROVENANCE")
        self.assertIn("MISSING_LOH_FIELDS", unresolved[0]["reason"])

    def test_mixed_manifest_publishes_only_valid_file(self):
        valid = simple_csv()
        conflict = simple_csv(current="Noise_Uniform", legacy="Strong", loh_subtype="LOH_Strong")
        self.fixture.add_input("valid.csv", valid)
        self.fixture.add_input("bad.csv", conflict)
        self.fixture.write_manifest(
            [
                self.fixture.manifest_entry("valid.csv", valid),
                self.fixture.manifest_entry("bad.csv", conflict),
            ]
        )
        output = self.fixture.root / "mixed-output"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 1, process.stderr)
        self.assertTrue((output / "valid.csv").is_file())
        self.assertFalse((output / "bad.csv").exists())
        status = tsv_rows(output / "migration_status.tsv")[0]
        self.assertEqual(status["valid_files"], "1")
        self.assertEqual(status["failed_files"], "1")

    def test_manifest_sha_and_row_mismatches_fail_before_output_creation(self):
        data = simple_csv()
        filename = "input.csv"
        self.fixture.add_input(filename, data)

        cases = (
            ("sha", {"file": filename, "rows": 1, "sha256": "0" * 64}, "SHA-256 mismatch"),
            ("rows", {"file": filename, "rows": 2, "sha256": digest(data)}, "row-count mismatch"),
        )
        for name, entry, message in cases:
            with self.subTest(case=name):
                self.fixture.write_manifest([entry])
                output = self.fixture.root / ("preflight-" + name)
                process = self.fixture.run(output)
                self.assertEqual(process.returncode, 2)
                self.assertIn(message, process.stderr)
                self.assertFalse(output.exists())

    def test_manifest_mapping_count_mismatch_does_not_publish_file(self):
        data = simple_csv()
        filename = "count-mismatch.csv"
        self.fixture.add_input(filename, data)
        entry = self.fixture.manifest_entry(
            filename,
            data,
            input_final_strong=99,
            legacy_strong_to_final_strong=1,
            legacy_subclone_to_final_strong=0,
            legacy_strong_or_subclone_exceptions=0,
        )
        self.fixture.write_manifest([entry], include_counts=True)
        output = self.fixture.root / "count-mismatch-output"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 1, process.stderr)
        self.assertFalse((output / filename).exists())
        report = tsv_rows(output / "migration_file_report.tsv")[0]
        self.assertEqual(report["status"], "FAILED_MANIFEST_COUNTS")
        self.assertIn("input_final_strong", report["reason"])

    def test_existing_output_is_refused_without_merging_or_overwriting(self):
        data = simple_csv()
        filename = "input.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "already-exists"
        output.mkdir()
        sentinel = output / "sentinel.txt"
        sentinel.write_text("keep", encoding="utf-8")

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 2)
        self.assertIn("already exists", process.stderr)
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep")
        self.assertEqual(list(output.iterdir()), [sentinel])

    def test_dangling_output_symlink_is_refused_without_following_target(self):
        data = simple_csv()
        filename = "input.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "dangling-output"
        target = self.fixture.root / "must-not-be-created"
        output.symlink_to(target, target_is_directory=True)

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 2)
        self.assertIn("already exists", process.stderr)
        self.assertTrue(output.is_symlink())
        self.assertEqual(os.readlink(output), str(target))
        self.assertFalse(target.exists())

    def test_atomic_no_replace_race_preserves_destination_and_cleans_staging(self):
        data = simple_csv()
        filename = "input.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "race-output"
        real_publish = MIGRATOR.atomic_publish_noreplace
        destination_before = {}

        def publish_after_destination_appears(staging, destination):
            destination.mkdir()
            sentinel = destination / "sentinel.txt"
            sentinel.write_bytes(b"keep-original")
            metadata = os.stat(destination, follow_symlinks=False)
            destination_before.update(
                inode=metadata.st_ino,
                device=metadata.st_dev,
                sentinel=sentinel.read_bytes(),
            )
            return real_publish(staging, destination)

        with mock.patch.object(
            MIGRATOR,
            "atomic_publish_noreplace",
            side_effect=publish_after_destination_appears,
        ):
            with self.assertRaisesRegex(
                MIGRATOR.MigrationError,
                "atomic no-replace publication refused",
            ):
                MIGRATOR.run_migration(
                    self.fixture.manifest,
                    self.fixture.input_root,
                    output,
                    "synthetic atomic race",
                )

        metadata_after = os.stat(output, follow_symlinks=False)
        self.assertEqual(metadata_after.st_ino, destination_before["inode"])
        self.assertEqual(metadata_after.st_dev, destination_before["device"])
        self.assertEqual((output / "sentinel.txt").read_bytes(), destination_before["sentinel"])
        self.assertEqual(list(output.iterdir()), [output / "sentinel.txt"])
        self.assertEqual(list(self.fixture.root.glob(".race-output.staging.*")), [])

    def test_unsupported_atomic_publish_refuses_fallback_and_cleans_staging(self):
        data = simple_csv()
        filename = "input.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "unsupported-output"
        unsupported = OSError(errno.EOPNOTSUPP, os.strerror(errno.EOPNOTSUPP))

        with mock.patch.object(
            MIGRATOR,
            "_call_libc_renameat2_noreplace",
            side_effect=unsupported,
        ), mock.patch.object(MIGRATOR.os, "rename") as rename_fallback:
            with self.assertRaisesRegex(
                MIGRATOR.MigrationError,
                "atomic no-replace publication is unsupported",
            ):
                MIGRATOR.run_migration(
                    self.fixture.manifest,
                    self.fixture.input_root,
                    output,
                    "synthetic unsupported renameat2",
                )

        rename_fallback.assert_not_called()
        self.assertFalse(output.exists())
        self.assertEqual(list(self.fixture.root.glob(".unsupported-output.staging.*")), [])

    def test_injected_exception_removes_only_current_staging(self):
        data = simple_csv()
        filename = "input.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "exception-output"
        unrelated = self.fixture.root / ".exception-output.staging.unrelated"
        unrelated.mkdir()
        sentinel = unrelated / "sentinel.txt"
        sentinel.write_text("do-not-delete", encoding="utf-8")

        with mock.patch.object(
            MIGRATOR,
            "write_reports",
            side_effect=RuntimeError("injected report failure"),
        ):
            with self.assertRaisesRegex(RuntimeError, "injected report failure"):
                MIGRATOR.run_migration(
                    self.fixture.manifest,
                    self.fixture.input_root,
                    output,
                    "synthetic injected failure",
                )

        self.assertFalse(output.exists())
        self.assertTrue(unrelated.is_dir())
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "do-not-delete")
        self.assertEqual(
            list(self.fixture.root.glob(".exception-output.staging.*")),
            [unrelated],
        )

    def test_cleanup_refuses_replaced_staging_identity(self):
        prefix = ".owned-output.staging."
        staging = Path(tempfile.mkdtemp(prefix=prefix, dir=str(self.fixture.root)))
        identity = MIGRATOR.capture_staging_identity(staging, self.fixture.root, prefix)
        original = self.fixture.root / "captured-original"
        staging.rename(original)
        staging.mkdir()
        sentinel = staging / "sentinel.txt"
        sentinel.write_text("replacement", encoding="utf-8")

        with self.assertRaisesRegex(
            MIGRATOR.MigrationError,
            "directory identity changed",
        ):
            MIGRATOR.cleanup_owned_staging(identity)

        self.assertTrue(original.is_dir())
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "replacement")

    def test_duplicate_stable_key_fails_and_reports_audit_flags(self):
        single = simple_csv()
        data = single + single.split(b"\n", 1)[1]
        filename = "duplicate-key.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "duplicate-key-output"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 1, process.stderr)
        self.assertFalse((output / filename).exists())
        conflicts = tsv_rows(output / "unmapped_conflicts.tsv")
        self.assertEqual(len(conflicts), 1)
        self.assertIn("DUPLICATE_STABLE_KEY", conflicts[0]["reason"])
        report = tsv_rows(output / "migration_file_report.tsv")[0]
        self.assertEqual(report["stable_key_uniqueness"], "FAIL")
        self.assertEqual(report["raw_token_preservation"], "NOT_RUN")
        self.assertEqual(report["significant_stable_key_invariant"], "NOT_RUN")
        summary = json.loads((output / "migration_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["stable_key_uniqueness"], "FAIL")
        self.assertEqual(summary["raw_token_preservation"], "FAIL")
        self.assertEqual(summary["significant_stable_key_invariant"], "FAIL")

    def test_output_inside_input_root_is_refused_as_in_place_migration(self):
        data = simple_csv()
        filename = "input.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.input_root / "nested-output"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 2)
        self.assertIn("must not be inside input root", process.stderr)
        self.assertFalse(output.exists())

    def test_required_absence_acknowledgement_flag_is_enforced(self):
        data = simple_csv()
        filename = "input.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "no-flag-output"

        process = self.fixture.run(output, include_absence_flag=False)
        self.assertEqual(process.returncode, 2)
        self.assertIn("--require-output-dir-absent", process.stderr)
        self.assertFalse(output.exists())

    def test_preexisting_v2_fields_are_refused_as_merge(self):
        headers = BASE_HEADERS + ("VerificationSchemaVersion",)
        data = simple_csv(headers=headers)
        filename = "already-v2.csv"
        self.fixture.add_input(filename, data)
        self.fixture.write_manifest([self.fixture.manifest_entry(filename, data)])
        output = self.fixture.root / "already-v2-output"

        process = self.fixture.run(output)
        self.assertEqual(process.returncode, 1, process.stderr)
        self.assertFalse((output / filename).exists())
        unresolved = tsv_rows(output / "unresolved_files.tsv")
        self.assertEqual(unresolved[0]["status"], "UNRESOLVED_SCHEMA")
        self.assertIn("OUTPUT_FIELDS_ALREADY_PRESENT", unresolved[0]["reason"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
