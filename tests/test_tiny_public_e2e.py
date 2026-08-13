import csv
import importlib.util
import json
import shutil
import subprocess
from argparse import Namespace
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "tiny_public"
VALIDATOR_PATH = ROOT / "scripts" / "handoff" / "validate_tiny_public_e2e.py"
SPEC = importlib.util.spec_from_file_location("tiny_public_validator", VALIDATOR_PATH)
assert SPEC and SPEC.loader
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


def test_expected_significance_schema_is_exact_v199_contract():
    schema = json.loads((FIXTURE / "expected_significance_schema.json").read_text(encoding="utf-8"))
    assert schema["column_count"] == 199
    assert len(schema["columns"]) == 199
    assert len(set(schema["columns"])) == 199
    assert schema["tree_semantics"] == VALIDATOR.TREE_SEMANTICS
    assert schema["columns"][-3:] == ["RegionStratum_ID", "RegionStratum_Label", "RegionStratum_Reason"]


def test_fixture_sources_are_small_deterministic_and_tagged():
    reference = "".join(
        line.strip() for line in (FIXTURE / "reference.fa").read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )
    records = [
        line.split("\t") for line in (FIXTURE / "tumor.sam").read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("@")
    ]
    assert len(reference) == 200
    assert len(records) == 12
    assert {record[0] for record in records} == set(VALIDATOR.EXPECTED_READ_NAMES)
    assert all(len(record[9]) == 200 and len(record[10]) == 200 for record in records)
    assert all(any(field.startswith("MM:Z:C+m?") for field in record[11:]) for record in records)
    assert all(any(field.startswith("ML:B:C,") for field in record[11:]) for record in records)


def test_csv_validator_rejects_schema_drift(tmp_path):
    schema_path = FIXTURE / "expected_significance_schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    row = [""] * 199
    values = {
        "Chr": "chrTiny", "Pos": "101", "NumReads": "12", "NumCpGs": "50",
        "NReadsValid": "12", "VerificationSchemaVersion": "2",
        "RegionStratificationSchemaVersion": "1", "RegionStratum_ID": "-1",
        "RegionStratum_Label": "Unassigned", "RegionStratum_Reason": "NOT_APPLICABLE_TUMOR_ONLY",
        "Subclone_ID": "-1",
    }
    for key, value in values.items():
        row[schema["columns"].index(key)] = value
    summary = tmp_path / "significance_summary.csv"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(schema["columns"])
        writer.writerow(row)
    checks = []
    result = VALIDATOR.validate_csv(summary, schema_path, checks)
    assert result["column_count"] == 199
    assert all(check["status"] == "PASS" for check in checks)

    bad_header = schema["columns"].copy()
    bad_header[-1] = "DriftedColumn"
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(bad_header)
        writer.writerow(row)
    checks = []
    VALIDATOR.validate_csv(summary, schema_path, checks)
    assert any(check["check_id"] == "SUMMARY_HEADER_EXACT" and check["status"] == "FAIL" for check in checks)


def test_tree_validator_locks_read_dendrogram_semantics(tmp_path):
    tree_dir = tmp_path / "variants" / "chrTiny" / "clustering"
    tree_dir.mkdir(parents=True)
    leaves = ",".join(f"{name}:0.1" for name in VALIDATOR.EXPECTED_READ_NAMES)
    (tree_dir / "tree.nwk").write_text(f"({leaves})0;\n", encoding="utf-8")
    checks = []
    result = VALIDATOR.validate_tree(tmp_path, checks)
    assert all(check["status"] == "PASS" for check in checks)
    assert result["leaf_count"] == 12
    assert result["semantics"] == "read_dendrogram_from_methylation_distance_not_cellular_lineage"
    assert "cellular_lineage" in result["not_evidence_for"]


def test_runner_writes_failure_receipt_when_core_exits_nonzero(tmp_path):
    if shutil.which("samtools") is None:
        pytest.skip("samtools is required to materialize the tiny BAM fixture")
    work_dir = tmp_path / "failed_e2e"
    completed = subprocess.run(
        [
            str(ROOT / "scripts" / "handoff" / "run_tiny_public_e2e.sh"),
            "--repo-root", str(ROOT),
            "--binary", "/bin/false",
            "--work-dir", str(work_dir),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert completed.returncode == 1
    receipt = json.loads((work_dir / "tiny_public_e2e_receipt.json").read_text(encoding="utf-8"))
    assert receipt["all_pass"] is False
    core_check = next(check for check in receipt["checks"] if check["check_id"] == "CORE_PROCESS_EXIT_ZERO")
    assert core_check == {
        "check_id": "CORE_PROCESS_EXIT_ZERO",
        "status": "FAIL",
        "observed": 1,
        "expected": 0,
    }


@pytest.mark.parametrize("forbidden_option", ["--binary", "--build-dir"])
def test_source_repository_rejects_external_build_identity_before_workdir_creation(
    tmp_path, forbidden_option
):
    work_dir = tmp_path / "must-not-exist"
    completed = subprocess.run(
        [
            str(ROOT / "scripts" / "handoff" / "run_tiny_public_e2e.sh"),
            "--source-repository", str(ROOT),
            "--revision", "HEAD",
            forbidden_option, str(tmp_path / "external"),
            "--work-dir", str(work_dir),
        ],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    assert completed.returncode == 2
    assert "forbidden" in completed.stderr
    assert not work_dir.exists()


@pytest.mark.parametrize(
    "repository",
    [
        "https://token@example.invalid/owner/repo.git",
        "https://example.invalid/owner/repo.git?token=secret",
        "https://example.invalid/owner/repo.git#credential-fragment",
    ],
)
def test_source_repository_rejects_credential_bearing_url_before_workdir_creation(
    tmp_path, repository
):
    work_dir = tmp_path / "must-not-exist"
    completed = subprocess.run(
        [
            str(ROOT / "scripts" / "handoff" / "run_tiny_public_e2e.sh"),
            "--source-repository", repository,
            "--revision", "deadbeef",
            "--work-dir", str(work_dir),
        ],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    assert completed.returncode == 2
    assert "must not contain URL credentials" in completed.stderr
    assert not work_dir.exists()


def test_receipt_repository_locator_redacts_url_credentials_and_tracking_data():
    assert VALIDATOR.public_repository_locator(
        "https://user:secret@example.invalid/owner/repo.git?token=secret#fragment"
    ) == "https://example.invalid/owner/repo.git"


def test_receipt_binds_requested_revision_resolved_commit_and_binary_identity(tmp_path):
    commit = subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
        check=True, text=True, stdout=subprocess.PIPE,
    ).stdout.strip()
    binary = tmp_path / "binary"
    binary.write_bytes(b"binary-identity")
    binary.chmod(0o755)
    args = Namespace(
        repo_root=ROOT,
        fixture_dir=tmp_path / "fixture",
        output_dir=tmp_path / "output",
        schema=FIXTURE / "expected_significance_schema.json",
        receipt=tmp_path / "receipt.json",
        binary=binary,
        run_log=None,
        command_receipt=None,
        run_exit_code=1,
        source_repository=str(ROOT),
        requested_revision="HEAD",
        resolved_commit=commit,
    )
    receipt = VALIDATOR.validate(args)
    assert receipt["source_checkout"]["requested_revision"] == "HEAD"
    assert receipt["source_checkout"]["resolved_commit"] == commit
    assert receipt["source_checkout"]["commit"] == commit
    assert receipt["binary_identity"] == {
        "path": str(binary.resolve()),
        "size_bytes": len(b"binary-identity"),
        "sha256": VALIDATOR.sha256(binary),
    }
    binding = next(
        check for check in receipt["checks"]
        if check["check_id"] == "SOURCE_REVISION_RESOLVED_COMMIT_BOUND"
    )
    assert binding["status"] == "PASS"
