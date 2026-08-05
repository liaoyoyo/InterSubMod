#!/usr/bin/env python3
"""Small synthetic tests for the independent schema-v2 migration validator."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
VALIDATOR_PATH = REPO / "scripts" / "validation" / "validate_verification_schema_v2_migration.py"
MODULE_SPEC = importlib.util.spec_from_file_location("verification_schema_v2_validator", VALIDATOR_PATH)
if MODULE_SPEC is None or MODULE_SPEC.loader is None:
    raise RuntimeError("cannot load validator module")
VALIDATOR = importlib.util.module_from_spec(MODULE_SPEC)
sys.modules[MODULE_SPEC.name] = VALIDATOR
MODULE_SPEC.loader.exec_module(VALIDATOR)


FILENAME = "sample.csv"
GENERATED_AT = "2026-07-15T00:00:00Z"
COMMAND = "python3 migrate_verification_schema_v2.py --frozen-test"
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
EXPECTED = VALIDATOR.ExpectedCorpus(
    total_files=1,
    total_rows=2,
    input_final_strong=2,
    strong_bidirectional=1,
    cluster_first_only=1,
    exceptions=0,
    hcc1395_file=FILENAME,
    hcc1395_input_strong=2,
    hcc1395_cluster_first_only=1,
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def tsv_bytes(headers, rows) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(headers), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({header: row.get(header, "") for header in headers})
    return stream.getvalue().encode("utf-8")


class ValidationFixture:
    def __init__(self, duplicate_key: bool = False):
        self.temporary = tempfile.TemporaryDirectory(prefix="verification-v2-validator-test-")
        self.root = Path(self.temporary.name)
        self.input_root = self.root / "input"
        self.output_dir = self.root / "output"
        self.input_root.mkdir()
        self.output_dir.mkdir()
        self.manifest = self.root / "manifest.tsv"
        self._write_corpus(duplicate_key=duplicate_key)

    def close(self):
        self.temporary.cleanup()

    def _write_corpus(self, duplicate_key: bool = False):
        second_region = b"1" if duplicate_key else b"2"
        second_identity = b",chr1,101,A,T," if duplicate_key else b",chr1,202,C,G,"
        header = b",".join(name.encode("ascii") for name in BASE_HEADERS) + b"\r\n"
        input_rows = (
            b'1,chr1,101,A,T,"Strong",Strong,true,LOH_Strong,true,"quoted,comma"\r\n'
            + second_region
            + second_identity
            + b"Strong,Subclone,true,LOH_Subclone,false,plain\r\n"
        )
        output_header = header[:-2] + b"," + b",".join(
            name.encode("ascii") for name in VALIDATOR.APPENDED_HEADERS
        ) + b"\r\n"
        output_rows = (
            b'1,chr1,101,A,T,"Strong_Bidirectional",Strong,true,LOH_Strong,true,"quoted,comma",'
            b'2,"Strong",true,true,NA,NA,BIDIRECTIONAL,LEGACY_CLASS,LOH_Strong\r\n'
            + second_region
            + second_identity
            + b"ClusterFirstOnly,Subclone,true,LOH_Subclone,false,plain,"
            b"2,Strong,false,true,NA,NA,CLUSTER_FIRST_ONLY,LEGACY_CLASS,LOH_Subclone\r\n"
        )
        self.input_data = header + input_rows
        self.output_data = output_header + output_rows
        (self.input_root / FILENAME).write_bytes(self.input_data)
        (self.output_dir / FILENAME).write_bytes(self.output_data)
        self._write_all_provenance()

    @staticmethod
    def counters():
        return {
            "before_current_counts": {"Strong": 2},
            "before_legacy_counts": {"Strong": 1, "Subclone": 1},
            "before_loh_counts": {"LOH_Strong": 1, "LOH_Subclone": 1},
            "after_current_counts": {"ClusterFirstOnly": 1, "Strong_Bidirectional": 1},
            "after_evidence_counts": {"BIDIRECTIONAL": 1, "CLUSTER_FIRST_ONLY": 1},
        }

    def _write_all_provenance(self):
        input_sha = sha256((self.input_root / FILENAME).read_bytes())
        output_sha = sha256((self.output_dir / FILENAME).read_bytes())
        manifest_row = {
            "file": FILENAME,
            "rows": 2,
            "sha256": input_sha,
            "input_final_strong": 2,
            "legacy_strong_to_final_strong": 1,
            "legacy_subclone_to_final_strong": 1,
            "legacy_strong_or_subclone_exceptions": 0,
        }
        self.manifest.write_bytes(tsv_bytes(VALIDATOR.MANIFEST_HEADERS, [manifest_row]))
        manifest_sha = sha256(self.manifest.read_bytes())
        counters = self.counters()
        file_row = {
            "file": FILENAME,
            "status": "VALID",
            "reason": "MAPPABLE",
            "input_rows": 2,
            "output_rows": 2,
            "input_sha256": input_sha,
            "output_sha256": output_sha,
            "schema_version": 2,
            "unmapped_count": 0,
            "raw_token_preservation": "PASS",
            "significant_stable_key_invariant": "PASS",
            "stable_key_uniqueness": "PASS",
        }
        file_row.update({field: json.dumps(value, sort_keys=True, separators=(",", ":")) for field, value in counters.items()})
        (self.output_dir / "migration_file_report.tsv").write_bytes(
            tsv_bytes(VALIDATOR.FILE_REPORT_HEADERS, [file_row])
        )
        (self.output_dir / "migrated_outputs_manifest.tsv").write_bytes(
            tsv_bytes(
                VALIDATOR.OUTPUT_MANIFEST_HEADERS,
                [{"file": FILENAME, "rows": 2, "sha256": output_sha, "schema_version": 2}],
            )
        )
        (self.output_dir / "unmapped_conflicts.tsv").write_bytes(
            tsv_bytes(VALIDATOR.CONFLICT_HEADERS, [])
        )
        (self.output_dir / "unresolved_files.tsv").write_bytes(
            tsv_bytes(VALIDATOR.UNRESOLVED_HEADERS, [])
        )
        summary = {
            "status": "VALID",
            "reason": "ALL_FILES_MIGRATED",
            "schema_version": 2,
            "generated_at": GENERATED_AT,
            "command": COMMAND,
            "manifest": str(self.manifest.resolve()),
            "manifest_sha256": manifest_sha,
            "input_root": str(self.input_root.resolve()),
            "output_dir": str(self.output_dir.resolve()),
            "total_files": 1,
            "valid_files": 1,
            "failed_files": 0,
            "input_rows": 2,
            "output_rows": 2,
            "unmapped_rows": 0,
            "raw_token_preservation": "PASS",
            "significant_stable_key_invariant": "PASS",
            "stable_key_uniqueness": "PASS",
        }
        summary.update(counters)
        (self.output_dir / "migration_summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        status_row = {
            "status": "VALID",
            "reason": "ALL_FILES_MIGRATED",
            "schema_version": 2,
            "total_files": 1,
            "valid_files": 1,
            "failed_files": 0,
            "input_rows": 2,
            "output_rows": 2,
            "unmapped_rows": 0,
            "raw_token_preservation": "PASS",
            "significant_stable_key_invariant": "PASS",
            "stable_key_uniqueness": "PASS",
            "generated_at": GENERATED_AT,
        }
        (self.output_dir / "migration_status.tsv").write_bytes(
            tsv_bytes(VALIDATOR.STATUS_HEADERS, [status_row])
        )
        (self.output_dir / "migration_command.txt").write_text(COMMAND + "\n", encoding="utf-8")

    def refresh_output_hash_provenance(self):
        output_sha = sha256((self.output_dir / FILENAME).read_bytes())
        rows = list(
            csv.DictReader(
                io.StringIO((self.output_dir / "migration_file_report.tsv").read_text(encoding="utf-8")),
                delimiter="\t",
            )
        )
        rows[0]["output_sha256"] = output_sha
        (self.output_dir / "migration_file_report.tsv").write_bytes(
            tsv_bytes(VALIDATOR.FILE_REPORT_HEADERS, rows)
        )
        (self.output_dir / "migrated_outputs_manifest.tsv").write_bytes(
            tsv_bytes(
                VALIDATOR.OUTPUT_MANIFEST_HEADERS,
                [{"file": FILENAME, "rows": 2, "sha256": output_sha, "schema_version": 2}],
            )
        )

    def validate(self):
        return VALIDATOR.validate_migration(
            self.manifest,
            self.input_root,
            self.output_dir,
            expected=EXPECTED,
        )


class IndependentValidatorTests(unittest.TestCase):
    def setUp(self):
        self.fixture = ValidationFixture()

    def tearDown(self):
        self.fixture.close()

    def test_success(self):
        watched = [self.fixture.manifest]
        watched.extend(sorted(self.fixture.input_root.iterdir()))
        watched.extend(sorted(self.fixture.output_dir.iterdir()))
        before = {path: sha256(path.read_bytes()) for path in watched}
        result = self.fixture.validate()
        self.assertEqual(result.total_files, 1)
        self.assertEqual(result.total_rows, 2)
        self.assertEqual(result.cluster_first_only, 1)
        self.assertEqual({path: sha256(path.read_bytes()) for path in watched}, before)

    def test_missing_provenance_report_fails_closed(self):
        (self.fixture.output_dir / "migration_status.tsv").unlink()
        with self.assertRaisesRegex(VALIDATOR.ValidationError, "missing"):
            self.fixture.validate()

    def test_output_hash_mismatch_fails_closed(self):
        path = self.fixture.output_dir / "migrated_outputs_manifest.tsv"
        rows = list(csv.DictReader(io.StringIO(path.read_text(encoding="utf-8")), delimiter="\t"))
        rows[0]["sha256"] = "0" * 64
        path.write_bytes(tsv_bytes(VALIDATOR.OUTPUT_MANIFEST_HEADERS, rows))
        with self.assertRaisesRegex(VALIDATOR.ValidationError, "output-manifest SHA"):
            self.fixture.validate()

    def test_reported_row_count_mismatch_fails_closed(self):
        path = self.fixture.output_dir / "migrated_outputs_manifest.tsv"
        path.write_bytes(path.read_bytes().replace(b"sample.csv\t2\t", b"sample.csv\t3\t"))
        with self.assertRaisesRegex(VALIDATOR.ValidationError, "rows"):
            self.fixture.validate()

    def test_duplicate_stable_key_fails_closed(self):
        self.fixture.close()
        self.fixture = ValidationFixture(duplicate_key=True)
        with self.assertRaisesRegex(VALIDATOR.ValidationError, "duplicate stable key"):
            self.fixture.validate()

    def test_significant_raw_token_drift_fails_closed(self):
        path = self.fixture.output_dir / FILENAME
        path.write_bytes(path.read_bytes().replace(b",false,plain,2,Strong", b",true,plain,2,Strong"))
        self.fixture.refresh_output_hash_provenance()
        with self.assertRaisesRegex(VALIDATOR.ValidationError, "Significant raw token"):
            self.fixture.validate()

    def test_existing_raw_token_drift_fails_closed(self):
        path = self.fixture.output_dir / FILENAME
        path.write_bytes(path.read_bytes().replace(b",plain,2,Strong", b",PLAIN,2,Strong"))
        self.fixture.refresh_output_hash_provenance()
        with self.assertRaisesRegex(VALIDATOR.ValidationError, "existing raw token"):
            self.fixture.validate()


if __name__ == "__main__":
    unittest.main()
