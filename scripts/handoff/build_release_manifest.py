#!/usr/bin/env python3
"""Build an out-of-tree InterSubMod research-handoff manifest and SHA256SUMS.

The builder reads Git assets from an immutable commit rather than from mutable
working-tree bytes. It refuses a dirty repository and never decides that a
blocked candidate is release-ready. Generated files are release artifacts and
must stay outside the repository to avoid a self-referential commit manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any, Iterable


HANDOFF_ROOT = Path("docs/handoff/20260813_完整研究資料與軟體交接_01")
AUTHORITY_PATH = HANDOFF_ROOT / "evidence/authority_manifest.json"
DENOMINATOR_PATH = HANDOFF_ROOT / "evidence/denominator_registry.tsv"
SCHEMA_PATH = HANDOFF_ROOT / "schemas/research-handoff-release-manifest.schema.json"
AGGREGATE_RECEIPT_SCHEMA_PATH = HANDOFF_ROOT / "schemas/aggregate-release-acceptance.schema.json"
LONGLINEAGE_SAFETY_SCHEMA_PATH = (
    HANDOFF_ROOT / "schemas/longlineage-preview-safety-receipt.schema.json"
)
COMMON_GATE_SCHEMA_PATH = HANDOFF_ROOT / "schemas/gate-acceptance-receipt.schema.json"
READER_RECEIPT_PATH = HANDOFF_ROOT / "evidence/reader_acceptance_receipt.json"
READER_SCHEMA_PATH = HANDOFF_ROOT / "schemas/reader-acceptance.schema.json"
READER_PROMPT_PATH = HANDOFF_ROOT / "ai_context/READER_ACCEPTANCE_PROMPT.md"
READER_VALIDATOR_PATH = Path("scripts/handoff/validate_reader_acceptance.py")
READER_VALIDATOR_TEST_PATH = Path("tests/test_validate_reader_acceptance.py")
AUTHORITY_REPLAY_PATH = HANDOFF_ROOT / "evidence/authority_replay_receipt.json"
ARTIFACT_REGISTRY_PATH = HANDOFF_ROOT / "registries/artifact_registry.json"
FROZEN_AUTHORITY_SHA256 = "a88afc589206b7b3dab4d17c1d7b02a6cf20b125d847409649577c02abf0bfa0"
FROZEN_DENOMINATOR_SHA256 = "a41f726cbf66f22b7e95ddb44216a08bf480d12961758ebd5b1ab6e3e61db6a4"
LONGLINEAGE_CANDIDATE = "b9aaa12a11fa00606bd174dabd0f172a5d112359"
LONGLINEAGE_HISTORICAL_SAFE_FAILURE_STACK = "f60b5f3274123bdf818371600b608d52626c893e"
LONGLINEAGE_SOURCE_MANIFEST_PATH = (
    HANDOFF_ROOT / "evidence/longlineage_source_to_target_manifest.json"
)
LONGLINEAGE_PUBLIC_GATE_RECEIPT_PATH = (
    HANDOFF_ROOT / "evidence/longlineage_public_gate_receipt.json"
)
LONGLINEAGE_HISTORY_RECEIPT_PATH = (
    HANDOFF_ROOT / "evidence/longlineage_history_scan_receipt.json"
)
LONGLINEAGE_PUBLIC_GATE_SCRIPT_PATH = (
    HANDOFF_ROOT / "evidence/longlineage_check_public_preview_gate.sh"
)
LONGLINEAGE_PUBLIC_GATE_OUTPUT_PATH = (
    HANDOFF_ROOT / "evidence/longlineage_public_gate_output.txt"
)
LONGLINEAGE_HISTORY_OUTPUT_PATH = (
    HANDOFF_ROOT / "evidence/longlineage_history_scan_output.json"
)
LICENSE_CONFIRMATION_BLOCKER = "INTERSUBMOD_PROJECT_LICENSE_CONFIRMATION_REQUIRED"
REQUIRED_AGGREGATE_GATES = (
    "authority_19_of_19_replay",
    "registry_and_handoff_package",
    "repo_hygiene",
    "tracked_large_asset_policy",
    "clean_release_build",
    "ctest_complete",
    "python_tests_complete",
    "tiny_synthetic_e2e",
    "bip7_fresh_clone_acceptance",
    "bip8_fresh_clone_acceptance",
    "html_browser_qa_84",
    "docs_claim_link_validation",
    "public_claim_registry_158",
    "github_publication_commit_alignment",
    "github_branch_protection_required_ci",
    "full_history_secret_scan",
    "release_asset_sha256_roundtrip",
    "reader_ai_six_question_acceptance",
    "intersubmod_project_license_confirmation",
    "longlineage_preview_safety",
)
REQUIRED_GATE_EVIDENCE_TYPES = {
    "authority_19_of_19_replay": "GIT_RECEIPT",
    "registry_and_handoff_package": "CI_RECEIPT",
    "repo_hygiene": "GIT_RECEIPT",
    "tracked_large_asset_policy": "GIT_RECEIPT",
    "clean_release_build": "CI_RECEIPT",
    "ctest_complete": "CI_RECEIPT",
    "python_tests_complete": "CI_RECEIPT",
    "tiny_synthetic_e2e": "CI_RECEIPT",
    "bip7_fresh_clone_acceptance": "HOST_RECEIPT",
    "bip8_fresh_clone_acceptance": "HOST_RECEIPT",
    "html_browser_qa_84": "CI_RECEIPT",
    "docs_claim_link_validation": "CI_RECEIPT",
    "public_claim_registry_158": "CI_RECEIPT",
    "github_publication_commit_alignment": "GITHUB_API_RECEIPT",
    "github_branch_protection_required_ci": "GITHUB_API_RECEIPT",
    "full_history_secret_scan": "CI_RECEIPT",
    "release_asset_sha256_roundtrip": "RELEASE_ASSET_RECEIPT",
    "reader_ai_six_question_acceptance": "READER_TEST_RECEIPT",
    "intersubmod_project_license_confirmation": "LICENSE_ATTESTATION",
    "longlineage_preview_safety": "RELEASE_ASSET_RECEIPT",
}
EXPECTED_DENOMINATORS = {
    "ssnv_dataset_records": ("469849", "469849", "100.0000", "strict_linkage_ready"),
    "strict_components": ("255752", "255752", "100.0000", "strict_linkage_ready"),
    "k1_linkage_abstain": ("170131", "255752", "66.5219", "strict_linkage_ready"),
    "strict_read_linked_w": ("85621", "255752", "33.4781", "strict_linkage_ready"),
    "mutation_bearing_units": ("85941", "98955", "86.8486", "all7_summary"),
    "family_complete": ("75224", "85941", "87.5298", "all7_summary"),
    "resource_abstain": ("10717", "85941", "12.4702", "all7_summary"),
    "ranked_complete": ("71955", "75224", "95.6543", "all7_summary"),
    "one_rooted_unlabeled_topology": ("63506", "71955", "88.2579", "topology_summary"),
    "methyl_formal_units": ("1045", "1045", "100.0000", "methyl_manifest"),
    "methyl_evaluable_units": ("811", "1045", "77.6077", "methyl_manifest"),
    "methyl_robust_associations": ("3", "811", "0.3699", "methyl_manifest"),
}
REQUIRED_GIT_ASSETS = (
    Path("CITATION.cff"),
    Path("CONTRIBUTING.md"),
    Path("SECURITY.md"),
    Path("CHANGELOG.md"),
    Path("LICENSE"),
    Path("README.md"),
    Path("README.zh-TW.md"),
    Path("QUICKSTART.md"),
    Path("requirements-ci.lock"),
    Path(".github/workflows/handoff-portable-ci.yml"),
    Path("config/site-profile.example.json"),
    Path("config/site-profile.schema.json"),
    HANDOFF_ROOT / "00_INDEX.md",
    HANDOFF_ROOT / "20260813_研究結論時間與Finality_01.md",
    HANDOFF_ROOT / "20260813_軟體輸入輸出與研究流程_01.md",
    HANDOFF_ROOT / "20260813_bip7_bip8操作與驗證_01.md",
    HANDOFF_ROOT / "20260813_完整研究交接總覽_01.standalone.html",
    HANDOFF_ROOT / "ai_context/CONTEXT.md",
    READER_PROMPT_PATH,
    HANDOFF_ROOT / "evidence/EVIDENCE_MANIFEST.json",
    AUTHORITY_PATH,
    DENOMINATOR_PATH,
    AUTHORITY_REPLAY_PATH,
    HANDOFF_ROOT / "evidence/algorithm_cli_claim_crosswalk.json",
    HANDOFF_ROOT / "evidence/full_claim_registry.json",
    HANDOFF_ROOT / "evidence/large_asset_migration_receipt.json",
    HANDOFF_ROOT / "evidence/longlineage_SBOM.spdx.json",
    HANDOFF_ROOT / "evidence/longlineage_clean_foundation_receipt.json",
    HANDOFF_ROOT / "evidence/longlineage_capability_matrix.md",
    HANDOFF_ROOT / "evidence/longlineage_public_safety_receipt.json",
    HANDOFF_ROOT / "evidence/longlineage_third_party_notices.md",
    HANDOFF_ROOT / "evidence/public_claim_validation_receipt.json",
    HANDOFF_ROOT / "evidence/repo_hygiene_summary_receipt.json",
    ARTIFACT_REGISTRY_PATH,
    HANDOFF_ROOT / "registries/authority_superseded_crosswalk.json",
    HANDOFF_ROOT / "registries/claim_registry.json",
    HANDOFF_ROOT / "registries/dataset_alias_registry.json",
    HANDOFF_ROOT / "registries/dataset_registry.json",
    HANDOFF_ROOT / "registries/machine_path_registry.json",
    HANDOFF_ROOT / "registries/registry_build_receipt.json",
    HANDOFF_ROOT / "registries/run_registry.json",
    HANDOFF_ROOT / "registries/storage_root_manifest.json",
    HANDOFF_ROOT / "registries/tracked_large_asset_registry.json",
    HANDOFF_ROOT / "registries/workflow_registry.json",
    SCHEMA_PATH,
    AGGREGATE_RECEIPT_SCHEMA_PATH,
    LONGLINEAGE_SAFETY_SCHEMA_PATH,
    COMMON_GATE_SCHEMA_PATH,
    READER_SCHEMA_PATH,
    HANDOFF_ROOT / "schemas/algorithm-cli-claim-crosswalk.schema.json",
    HANDOFF_ROOT / "schemas/artifact-registry.schema.json",
    HANDOFF_ROOT / "schemas/authority-crosswalk.schema.json",
    HANDOFF_ROOT / "schemas/dataset-registry.schema.json",
    HANDOFF_ROOT / "schemas/machine-path-registry.schema.json",
    HANDOFF_ROOT / "schemas/run-registry.schema.json",
    HANDOFF_ROOT / "schemas/storage-root-manifest.schema.json",
    HANDOFF_ROOT / "schemas/tracked-large-asset-registry.schema.json",
    Path("scripts/site/bootstrap"),
    Path("scripts/site/doctor"),
    Path("scripts/site/site_profile.py"),
    Path("scripts/handoff/build_algorithm_cli_crosswalk.py"),
    Path("scripts/handoff/build_large_asset_registry.py"),
    Path("scripts/handoff/build_registries.py"),
    Path("scripts/handoff/build_tiny_public_fixture.sh"),
    Path("scripts/handoff/build_workflow_registry.py"),
    Path("scripts/handoff/build_release_manifest.py"),
    READER_VALIDATOR_PATH,
    Path("scripts/handoff/replay_authority.py"),
    Path("scripts/handoff/repo_hygiene.py"),
    Path("scripts/handoff/run_tiny_public_e2e.sh"),
    Path("scripts/handoff/sync_public_claim_evidence.py"),
    Path("scripts/handoff/validate_algorithm_cli_crosswalk.py"),
    Path("scripts/handoff/validate_handoff_package.py"),
    Path("scripts/handoff/validate_tiny_public_e2e.py"),
    Path("tests/fixtures/tiny_public/README.md"),
    Path("tests/fixtures/tiny_public/expected_significance_schema.json"),
    Path("tests/fixtures/tiny_public/reference.fa"),
    Path("tests/fixtures/tiny_public/tumor.sam"),
    Path("tests/fixtures/tiny_public/variants.vcf"),
    READER_VALIDATOR_TEST_PATH,
)
APPROVED_ONLY_GIT_ASSETS = (
    READER_RECEIPT_PATH,
    LONGLINEAGE_SOURCE_MANIFEST_PATH,
    LONGLINEAGE_PUBLIC_GATE_RECEIPT_PATH,
    LONGLINEAGE_HISTORY_RECEIPT_PATH,
    LONGLINEAGE_PUBLIC_GATE_SCRIPT_PATH,
    LONGLINEAGE_PUBLIC_GATE_OUTPUT_PATH,
    LONGLINEAGE_HISTORY_OUTPUT_PATH,
)
READER_CORE_PATHS = {
    "00_INDEX.md",
    "20260813_研究結論時間與Finality_01.md",
    "20260813_軟體輸入輸出與研究流程_01.md",
    "20260813_bip7_bip8操作與驗證_01.md",
    "ai_context/CONTEXT.md",
    "ai_context/READER_ACCEPTANCE_PROMPT.md",
    "schemas/reader-acceptance.schema.json",
    "registries/artifact_registry.json",
    "registries/authority_superseded_crosswalk.json",
    "registries/claim_registry.json",
    "registries/machine_path_registry.json",
    "evidence/authority_manifest.json",
    "evidence/denominator_registry.tsv",
    "evidence/authority_replay_receipt.json",
    "evidence/longlineage_capability_matrix.md",
}
READER_QUESTION_IDS = (
    "Q_PROJECT",
    "Q_CONCLUSION",
    "Q_FINALITY",
    "Q_SOFTWARE_ROLES",
    "Q_MACHINES",
    "Q_VERIFY_CONTINUE",
)
READER_PROMOTION_IDS = (
    "NO_CELLULAR_PROMOTION",
    "NO_ANCESTRY_PROMOTION",
    "NO_882579_ACCURACY_OR_PREVALENCE",
    "NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE",
    "NO_FEATURE_AS_MAIN",
    "NO_LOCAL_AS_LIVE_PUBLISHED",
    "NO_BIP7_AS_BIP8",
    "NO_PYTHON_ROLE_CONFLATION",
)
LONGLINEAGE_GATES = ("P3", "P4", "P5", "P7", "P8")
LONGLINEAGE_SUPPORT_PATHS = {
    "license": (
        HANDOFF_ROOT / "evidence/longlineage_third_party_notices.md",
        HANDOFF_ROOT / "evidence/longlineage_SBOM.spdx.json",
    ),
    "source_origin": (
        HANDOFF_ROOT / "evidence/longlineage_public_safety_receipt.json",
        LONGLINEAGE_SOURCE_MANIFEST_PATH,
    ),
    "dependencies": (
        HANDOFF_ROOT / "evidence/longlineage_SBOM.spdx.json",
    ),
    "public_safety": (
        HANDOFF_ROOT / "evidence/longlineage_public_safety_receipt.json",
        HANDOFF_ROOT / "evidence/longlineage_clean_foundation_receipt.json",
        LONGLINEAGE_PUBLIC_GATE_RECEIPT_PATH,
        LONGLINEAGE_HISTORY_RECEIPT_PATH,
        LONGLINEAGE_PUBLIC_GATE_SCRIPT_PATH,
        LONGLINEAGE_PUBLIC_GATE_OUTPUT_PATH,
        LONGLINEAGE_HISTORY_OUTPUT_PATH,
    ),
}
SPECIAL_GATE_IDS = {
    "authority_19_of_19_replay",
    "repo_hygiene",
    "tracked_large_asset_policy",
    "longlineage_preview_safety",
}
GIT_LIMIT_BYTES = 50 * 1024 * 1024
RELEASE_LIMIT_BYTES = 2 * 1024 * 1024 * 1024
FORBIDDEN_RELEASE_ASSET_SUFFIXES = (
    ".bam",
    ".bam.bai",
    ".cram",
    ".crai",
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".vcf",
    ".vcf.gz",
    ".bcf",
    ".csi",
    ".tbi",
    ".env",
    ".key",
    ".pem",
    ".p12",
    ".pfx",
    ".log",
    ".sqlite",
    ".sqlite3",
)
FORBIDDEN_RELEASE_ASSET_NAME_TOKENS = (
    "credential",
    "password",
    "private_key",
    "secret",
    "token",
)


class ManifestError(RuntimeError):
    """A fail-closed manifest construction error."""


def run_git(repo: Path, *args: str, text: bool = True) -> str | bytes:
    result = subprocess.run(
        ["git", *args],
        cwd=repo,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    )
    if result.returncode != 0:
        stderr = result.stderr.strip() if text else result.stderr.decode("utf-8", "replace").strip()
        raise ManifestError(f"git {' '.join(args)} failed: {stderr}")
    return result.stdout


def normalize_repository_path(raw: str | Path) -> str:
    value = str(raw).replace("\\", "/")
    path = PurePosixPath(value)
    if path.is_absolute() or not value or any(part in ("", ".", "..") for part in path.parts):
        raise ManifestError(f"repository path must be normalized and relative: {raw}")
    return path.as_posix()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_blob_info(repo: Path, commit: str, repository_path: str) -> tuple[str, str, int, str]:
    output = run_git(repo, "ls-tree", "-z", commit, "--", repository_path, text=False)
    records = [record for record in output.split(b"\0") if record]
    if len(records) != 1:
        raise ManifestError(f"expected one tracked asset at {commit}:{repository_path}")
    metadata, encoded_path = records[0].split(b"\t", 1)
    mode, object_type, object_id = metadata.decode("ascii").split()
    actual_path = encoded_path.decode("utf-8")
    if actual_path != repository_path or object_type != "blob" or mode not in ("100644", "100755"):
        raise ManifestError(
            f"asset must be a regular tracked file, not a symlink/submodule: {repository_path}"
        )
    size = int(run_git(repo, "cat-file", "-s", object_id).strip())
    if size > GIT_LIMIT_BYTES:
        raise ManifestError(f"regular Git asset exceeds 50 MiB policy: {repository_path} ({size})")

    digest = hashlib.sha256()
    process = subprocess.Popen(
        ["git", "cat-file", "blob", object_id],
        cwd=repo,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    for block in iter(lambda: process.stdout.read(8 * 1024 * 1024), b""):
        digest.update(block)
    stderr = process.stderr.read().decode("utf-8", "replace") if process.stderr else ""
    if process.wait() != 0:
        raise ManifestError(f"unable to read Git blob {object_id}: {stderr.strip()}")
    return mode, object_id, size, digest.hexdigest()


def git_file_bytes(repo: Path, commit: str, repository_path: str) -> bytes:
    normalize_repository_path(repository_path)
    return run_git(repo, "show", f"{commit}:{repository_path}", text=False)


def build_git_asset(repo: Path, commit: str, repository_path: str) -> dict[str, Any]:
    mode, object_id, size, digest = git_blob_info(repo, commit, repository_path)
    path_id = hashlib.sha256(repository_path.encode("utf-8")).hexdigest()[:16]
    return {
        "asset_id": f"git/{path_id}",
        "publication_channel": "GIT",
        "published_path": repository_path,
        "repository_path": repository_path,
        "git_blob": object_id,
        "git_mode": mode,
        "size_bytes": size,
        "sha256": digest,
    }


def parse_release_asset(spec: str) -> tuple[str, Path]:
    if "=" in spec:
        published_name, raw_path = spec.split("=", 1)
    else:
        raw_path = spec
        published_name = Path(raw_path).name
    published = PurePosixPath(published_name)
    if (
        not published_name
        or published.is_absolute()
        or len(published.parts) != 1
        or published.name in (".", "..")
    ):
        raise ManifestError(f"release asset name must be one filename: {published_name}")
    return published_name, Path(raw_path).expanduser().resolve()


def build_release_asset(repo: Path, spec: str) -> dict[str, Any]:
    published_name, source_path = parse_release_asset(spec)
    try:
        source_path.relative_to(repo)
    except ValueError:
        pass
    else:
        raise ManifestError(
            f"GitHub Release asset source must be outside the repository: {source_path}"
        )
    if source_path.is_symlink() or not source_path.is_file():
        raise ManifestError(f"release asset must be a regular file: {source_path}")
    size = source_path.stat().st_size
    if size >= RELEASE_LIMIT_BYTES:
        raise ManifestError(f"GitHub Release asset must be smaller than 2 GiB: {source_path}")
    path_id = hashlib.sha256(published_name.encode("utf-8")).hexdigest()[:16]
    return {
        "asset_id": f"release/{path_id}",
        "publication_channel": "GITHUB_RELEASE",
        "published_path": f"release-assets/{published_name}",
        "release_source_name": source_path.name,
        "size_bytes": size,
        "sha256": sha256_file(source_path),
    }


def load_artifact_registry(raw: bytes) -> dict[str, dict[str, Any]]:
    try:
        document = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ManifestError("artifact registry is not JSON") from exc
    artifacts = document if isinstance(document, list) else document.get("artifacts", [])
    if not isinstance(artifacts, list):
        raise ManifestError("artifact registry is not a list")
    by_id: dict[str, dict[str, Any]] = {}
    for artifact in artifacts:
        if not isinstance(artifact, dict) or not isinstance(artifact.get("artifact_id"), str):
            raise ManifestError("artifact registry contains a record without artifact_id")
        artifact_id = artifact["artifact_id"]
        if artifact_id in by_id:
            raise ManifestError(f"artifact registry contains duplicate ID: {artifact_id}")
        by_id[artifact_id] = artifact
    return by_id


def build_governed_release_asset(
    repo: Path,
    spec: str,
    artifact_registry: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    if "::" not in spec:
        raise ManifestError(
            "--release-asset must use ARTIFACT_ID::[PUBLISHED_NAME=]PATH syntax"
        )
    artifact_id, asset_spec = spec.split("::", 1)
    if not artifact_id or not asset_spec:
        raise ManifestError(
            "--release-asset must use ARTIFACT_ID::[PUBLISHED_NAME=]PATH syntax"
        )
    record = artifact_registry.get(artifact_id)
    if record is None:
        raise ManifestError(f"release asset has no governed artifact record: {artifact_id}")
    asset = build_release_asset(repo, asset_spec)
    published_name = PurePosixPath(asset["published_path"]).name
    lowered_name = published_name.lower()
    if any(lowered_name.endswith(suffix) for suffix in FORBIDDEN_RELEASE_ASSET_SUFFIXES):
        raise ManifestError(f"raw, runtime, or secret-like release asset is forbidden: {published_name}")
    if any(token in lowered_name for token in FORBIDDEN_RELEASE_ASSET_NAME_TOKENS):
        raise ManifestError(f"secret-like release asset name is forbidden: {published_name}")
    if record.get("finality") != "FINAL_FOR_SCOPE":
        raise ManifestError(f"release asset is not FINAL_FOR_SCOPE: {artifact_id}")
    if record.get("evidence_status") not in ("AUTHORITY", "VALIDATED_DERIVED"):
        raise ManifestError(f"release asset evidence status is not publishable: {artifact_id}")
    if record.get("availability") != "GITHUB_RELEASE":
        raise ManifestError(f"release asset registry channel is not GITHUB_RELEASE: {artifact_id}")
    license_value = record.get("license")
    if not isinstance(license_value, str) or not license_value.strip():
        raise ManifestError(f"release asset license is empty: {artifact_id}")
    lowered_license = license_value.lower()
    if any(
        token in lowered_license
        for token in (
            "unknown",
            "review required",
            "tracked separately",
            "pending_maintainer_confirmation",
        )
    ):
        raise ManifestError(f"release asset license is unresolved: {artifact_id}")
    if record.get("sha256") != asset["sha256"]:
        raise ManifestError(f"release asset SHA-256 does not match artifact registry: {artifact_id}")
    if record.get("size_bytes") != asset["size_bytes"]:
        raise ManifestError(f"release asset size does not match artifact registry: {artifact_id}")
    public_uri = record.get("public_uri")
    expected_uri = asset["published_path"]
    if not isinstance(public_uri, str) or not (
        public_uri == expected_uri or public_uri.rstrip("/").endswith(f"/{published_name}")
    ):
        raise ManifestError(f"release asset public_uri does not match published name: {artifact_id}")
    asset["registry_artifact_id"] = artifact_id
    return asset


def build_gate_evidence_asset(
    repo: Path,
    commit: str,
    spec: str,
    schema: dict[str, Any],
) -> dict[str, Any]:
    asset = build_release_asset(repo, spec)
    published_name = PurePosixPath(asset["published_path"]).name
    lowered_name = published_name.lower()
    if not lowered_name.endswith(".json"):
        raise ManifestError(f"gate evidence asset must be JSON: {published_name}")
    if asset["size_bytes"] > 10 * 1024 * 1024:
        raise ManifestError(f"gate evidence receipt exceeds 10 MiB: {published_name}")
    _, source_path = parse_release_asset(spec)
    try:
        receipt = json.loads(source_path.read_bytes())
    except json.JSONDecodeError as exc:
        raise ManifestError(f"gate evidence receipt is not JSON: {source_path}") from exc
    validate_document(receipt, schema, f"gate evidence receipt {source_path}")
    verify_common_gate_receipt(repo, receipt, commit, source_path)
    asset["gate_evidence_receipt"] = True
    asset["gate_id"] = receipt["gate_id"]
    asset["execution_log_uri"] = receipt["execution"]["log_uri"]
    asset["execution_log_sha256"] = receipt["execution"]["log_sha256"]
    return asset


def build_gate_log_asset(repo: Path, spec: str) -> dict[str, Any]:
    asset = build_release_asset(repo, spec)
    if asset["size_bytes"] == 0 or asset["size_bytes"] > 10 * 1024 * 1024:
        raise ManifestError("gate execution log must be nonempty and at most 10 MiB")
    asset["gate_execution_log"] = True
    return asset


def verify_common_gate_receipt(
    repo: Path,
    receipt: Any,
    commit: str,
    source_path: Path,
) -> None:
    if not isinstance(receipt, dict):
        raise ManifestError(f"gate evidence receipt must be an object: {source_path}")
    gate_id = receipt.get("gate_id")
    if gate_id not in REQUIRED_AGGREGATE_GATES or gate_id in SPECIAL_GATE_IDS:
        raise ManifestError(f"gate evidence receipt has invalid gate_id: {source_path}")
    if receipt.get("schema_version") != "1.0.0":
        raise ManifestError(f"gate evidence receipt schema version mismatch: {source_path}")
    if receipt.get("receipt_type") != f"{gate_id}_acceptance":
        raise ManifestError(f"gate evidence receipt_type does not match gate_id: {source_path}")
    if receipt.get("source_commit") != commit:
        raise ManifestError(f"gate evidence receipt does not bind source commit: {source_path}")
    if receipt.get("overall_pass") is not True or receipt.get("blocking_gates") != []:
        raise ManifestError(f"gate evidence receipt is not all-pass: {source_path}")
    if receipt.get("scope") != "FULL":
        raise ManifestError(f"gate evidence receipt scope is not FULL: {source_path}")
    try:
        parse_generated_at(receipt["generated_at"])
    except (KeyError, ManifestError) as exc:
        raise ManifestError(f"gate evidence receipt timestamp is invalid: {source_path}") from exc
    execution = receipt.get("execution")
    if not isinstance(execution, dict):
        raise ManifestError(f"gate evidence receipt lacks execution object: {source_path}")
    expected_kind = {
        "github_publication_commit_alignment": "GITHUB_API",
        "github_branch_protection_required_ci": "GITHUB_API",
        "reader_ai_six_question_acceptance": "READER_AGENT",
        "intersubmod_project_license_confirmation": "GITHUB_API",
    }.get(gate_id, "COMMAND")
    if execution.get("kind") != expected_kind:
        raise ManifestError(f"gate evidence execution kind is invalid: {source_path}")
    for field in ("tool", "version", "command"):
        if not isinstance(execution.get(field), str) or not execution[field].strip():
            raise ManifestError(f"gate evidence execution {field} is empty: {source_path}")
    if not re.fullmatch(r"[0-9a-f]{64}", str(execution.get("log_sha256"))):
        raise ManifestError(f"gate evidence execution log hash is invalid: {source_path}")
    if not re.fullmatch(
        r"release-assets/[A-Za-z0-9._-]+", str(execution.get("log_uri"))
    ):
        raise ManifestError(f"gate evidence execution log URI is invalid: {source_path}")
    if execution.get("exit_code") != 0:
        raise ManifestError(f"gate evidence command did not exit zero: {source_path}")
    command_lower = execution["command"].lower()
    if gate_id == "full_history_secret_scan" and not any(
        token in command_lower for token in ("--all", "full-history", "full_history")
    ):
        raise ManifestError(f"secret scan command does not cover full history: {source_path}")
    result = receipt.get("result")
    if not isinstance(result, dict) or result.get("status") != "PASS":
        raise ManifestError(f"gate evidence result is not PASS: {source_path}")
    if not isinstance(result.get("checks"), int) or result["checks"] < 1:
        raise ManifestError(f"gate evidence result has no checks: {source_path}")
    if result.get("failures") != 0:
        raise ManifestError(f"gate evidence result has failures: {source_path}")

    details = receipt.get("details")
    if not isinstance(details, dict):
        raise ManifestError(f"gate evidence details are missing: {source_path}")
    verify_gate_details(gate_id, details, commit, source_path)
    if gate_id == "reader_ai_six_question_acceptance":
        expected_reader_sha = sha256_bytes(
            git_file_bytes(repo, commit, READER_RECEIPT_PATH.as_posix())
        )
        if details.get("detailed_receipt_sha256") != expected_reader_sha:
            raise ManifestError(f"reader gate summary does not hash-bind detailed receipt: {source_path}")

    if gate_id in ("bip7_fresh_clone_acceptance", "bip8_fresh_clone_acceptance"):
        host = receipt.get("host")
        if not isinstance(host, dict) or not all(
            isinstance(host.get(key), str) and host[key]
            for key in ("hostname", "mount_type")
        ):
            raise ManifestError(f"host acceptance receipt lacks host/mount identity: {source_path}")
        tool_hashes = host.get("tool_hashes")
        if not isinstance(tool_hashes, dict) or not tool_hashes:
            raise ManifestError(f"host acceptance receipt lacks tool hashes: {source_path}")
        if not all(re.fullmatch(r"[0-9a-f]{64}", str(value)) for value in tool_hashes.values()):
            raise ManifestError(f"host acceptance receipt has invalid tool hash: {source_path}")
    elif gate_id in (
        "github_publication_commit_alignment",
        "github_branch_protection_required_ci",
    ):
        github = receipt.get("github")
        if not isinstance(github, dict):
            raise ManifestError(f"GitHub acceptance receipt lacks API evidence: {source_path}")
        if github.get("repository_url") != "https://github.com/liaoyoyo/InterSubMod":
            raise ManifestError(f"GitHub acceptance repository mismatch: {source_path}")
        if not isinstance(github.get("api_endpoint"), str) or not github["api_endpoint"]:
            raise ManifestError(f"GitHub acceptance API endpoint is empty: {source_path}")
        raw_response = github.get("raw_api_response")
        if not isinstance(raw_response, str) or sha256_bytes(raw_response.encode()) != github.get(
            "response_sha256"
        ):
            raise ManifestError(f"GitHub acceptance response hash is invalid: {source_path}")


def require_exact_counts(details: dict[str, Any], pairs: tuple[tuple[str, str], ...]) -> None:
    for expected, observed in pairs:
        if details.get(expected) != details.get(observed):
            raise ManifestError(f"gate detail counts differ: {expected} != {observed}")


def verify_gate_details(
    gate_id: str,
    details: dict[str, Any],
    commit: str,
    source_path: Path,
) -> None:
    if gate_id == "registry_and_handoff_package":
        if details.get("schema_errors") != 0 or details.get("package_errors") != 0:
            raise ManifestError(f"registry/package validation failed: {source_path}")
    elif gate_id == "clean_release_build":
        if details != {"build_type": "Release", "configure_exit": 0, "build_exit": 0}:
            raise ManifestError(f"clean Release build details are invalid: {source_path}")
    elif gate_id == "ctest_complete":
        require_exact_counts(details, (("tests_total", "tests_passed"),))
        if details.get("tests_failed") != 0 or "ctest" not in details.get(
            "discovery_command", ""
        ):
            raise ManifestError(f"dynamic CTest evidence is invalid: {source_path}")
    elif gate_id == "python_tests_complete":
        require_exact_counts(details, (("tests_total", "tests_passed"),))
        if details.get("tests_failed") != 0 or not str(details.get("python_version", "")).startswith(
            "3.10"
        ):
            raise ManifestError(f"Python 3.10 test evidence is invalid: {source_path}")
    elif gate_id == "tiny_synthetic_e2e":
        if any(details.get(key) != 0 for key in ("plan_exit", "plan_side_effects", "run_exit")):
            raise ManifestError(f"tiny E2E plan/run evidence is invalid: {source_path}")
        if details.get("schema_valid") is not True or details.get("hash_valid") is not True:
            raise ManifestError(f"tiny E2E schema/hash evidence is invalid: {source_path}")
    elif gate_id in ("bip7_fresh_clone_acceptance", "bip8_fresh_clone_acceptance"):
        exit_fields = (
            "clone_exit",
            "build_exit",
            "ctest_exit",
            "python_tests_exit",
            "synthetic_smoke_exit",
            "real_data_preflight_exit",
        )
        if any(details.get(key) != 0 for key in exit_fields):
            raise ManifestError(f"host acceptance has a nonzero step: {source_path}")
        if details.get("fresh_clone") is not True or details.get("real_data_read_only") is not True:
            raise ManifestError(f"host acceptance is not fresh/read-only: {source_path}")
    elif gate_id == "html_browser_qa_84":
        if details.get("checks_total") != 84 or details.get("checks_passed") != 84:
            raise ManifestError(f"HTML QA is not 84/84: {source_path}")
    elif gate_id == "docs_claim_link_validation":
        if any(value != 0 for value in details.values()):
            raise ManifestError(f"documentation validation has failures: {source_path}")
    elif gate_id == "public_claim_registry_158":
        if details.get("rows_total") != 158 or details.get("invalid_rows") != 0:
            raise ManifestError(f"claim inventory is not valid 158 rows: {source_path}")
        if sum(details.get(key, -1000) for key in ("confirmed", "confirmed_with_limits", "unverified")) != 158:
            raise ManifestError(f"claim inventory disposition counts do not total 158: {source_path}")
        if details.get("p0_unresolved") != 0:
            raise ManifestError(f"claim inventory retains unresolved P0 rows: {source_path}")
    elif gate_id == "github_publication_commit_alignment":
        if any(details.get(key) != commit for key in ("main_commit", "wiki_commit", "pages_commit")):
            raise ManifestError(f"GitHub publication commits are not aligned: {source_path}")
        if details.get("readme_zh_http_status") != 200:
            raise ManifestError(f"default Chinese README is not HTTP 200: {source_path}")
    elif gate_id == "github_branch_protection_required_ci":
        if details.get("branch") != "main" or details.get("protected") is not True:
            raise ManifestError(f"main branch is not protected: {source_path}")
        if details.get("enforcement") is not True or not details.get("required_ci"):
            raise ManifestError(f"required CI enforcement is absent: {source_path}")
    elif gate_id == "full_history_secret_scan":
        if details.get("scope") != "FULL_HISTORY" or details.get("findings") != 0:
            raise ManifestError(f"secret scan is not full-history/zero-finding: {source_path}")
    elif gate_id == "release_asset_sha256_roundtrip":
        expected = details.get("assets_expected")
        if expected != details.get("assets_downloaded") or expected != details.get("hash_matches"):
            raise ManifestError(f"release roundtrip asset counts differ: {source_path}")
        if expected != details.get("http_200_count"):
            raise ManifestError(f"release roundtrip has non-200 assets: {source_path}")
    elif gate_id == "reader_ai_six_question_acceptance":
        if (
            details.get("questions_total") != 6
            or details.get("questions_passed") != 6
            or details.get("promotions_total") != 8
            or details.get("promotions_passed") != 8
            or details.get("validator_exit") != 0
        ):
            raise ManifestError(f"reader acceptance summary is not 6/6 + 8/8: {source_path}")
    elif gate_id == "intersubmod_project_license_confirmation":
        verify_license_owner_approval(details, commit, source_path)


def verify_license_owner_approval(
    details: dict[str, Any], commit: str, source_path: Path
) -> None:
    if details.get("maintainer_login") != "liaoyoyo" or details.get("approval_commit") != commit:
        raise ManifestError(f"license authority/commit mismatch: {source_path}")
    raw_response = details.get("raw_api_response")
    if not isinstance(raw_response, str) or sha256_bytes(raw_response.encode()) != details.get(
        "response_sha256"
    ):
        raise ManifestError(f"license approval API response hash mismatch: {source_path}")
    try:
        response = json.loads(raw_response)
    except json.JSONDecodeError as exc:
        raise ManifestError(f"license approval API response is not JSON: {source_path}") from exc
    if not isinstance(response, dict) or any(
        response.get(key) != value
        for key, value in (
            ("user_login", "liaoyoyo"),
            ("author_association", "OWNER"),
            ("state", "APPROVED"),
            ("head_sha", commit),
        )
    ):
        raise ManifestError(f"license approval is not an owner-approved PR review: {source_path}")


def verify_aggregate_evidence_bytes(
    repo: Path,
    commit: str,
    aggregate_receipt: dict[str, Any],
    declared_assets: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    by_published_path = {asset["published_path"]: asset for asset in declared_assets}
    discovered_git_assets: dict[str, dict[str, Any]] = {}
    referenced_release_assets: set[str] = set()
    validated_gate_receipts: set[str] = set()
    for gate_id in REQUIRED_AGGREGATE_GATES:
        for evidence in aggregate_receipt["gates"][gate_id]["evidence"]:
            uri = evidence["uri"]
            expected_sha = evidence["sha256"]
            if uri.startswith("release-assets/"):
                asset = by_published_path.get(uri)
                if asset is None or asset.get("sha256") != expected_sha:
                    raise ManifestError(
                        f"aggregate gate evidence bytes are not declared: {gate_id} {uri}"
                    )
                referenced_release_assets.add(uri)
                if asset.get("gate_evidence_receipt") is True:
                    if asset.get("gate_id") != gate_id:
                        raise ManifestError(
                            f"gate evidence receipt is attached to wrong gate: {gate_id} {uri}"
                        )
                    validated_gate_receipts.add(gate_id)
                continue
            repository_path = normalize_repository_path(uri)
            asset = build_git_asset(repo, commit, repository_path)
            if asset["sha256"] != expected_sha:
                raise ManifestError(
                    f"aggregate gate evidence hash does not match Git bytes: {gate_id} {uri}"
                )
            discovered_git_assets[repository_path] = asset

    for asset in declared_assets:
        if asset.get("gate_evidence_receipt") is True:
            path = asset["published_path"]
            if path not in referenced_release_assets:
                raise ManifestError(f"unreferenced --gate-evidence-asset: {path}")
    gate_receipts = [
        asset for asset in declared_assets if asset.get("gate_evidence_receipt") is True
    ]
    referenced_logs: set[str] = set()
    for gate_asset in gate_receipts:
        gate_id = gate_asset["gate_id"]
        log_uri = gate_asset["execution_log_uri"]
        log_sha = gate_asset["execution_log_sha256"]
        log_asset = by_published_path.get(log_uri)
        if (
            log_asset is None
            or log_asset.get("gate_execution_log") is not True
            or log_asset.get("sha256") != log_sha
        ):
            raise ManifestError(
                f"gate receipt does not bind a declared raw execution log: {gate_id} {log_uri}"
            )
        if not any(
            item.get("uri") == log_uri and item.get("sha256") == log_sha
            for item in aggregate_receipt["gates"][gate_id]["evidence"]
        ):
            raise ManifestError(
                f"aggregate gate does not bind raw execution log: {gate_id} {log_uri}"
            )
        referenced_logs.add(log_uri)
    declared_logs = {
        asset["published_path"]
        for asset in declared_assets
        if asset.get("gate_execution_log") is True
    }
    if declared_logs != referenced_logs:
        raise ManifestError("unreferenced or missing --gate-log-asset")
    expected_common = set(REQUIRED_AGGREGATE_GATES) - SPECIAL_GATE_IDS
    if validated_gate_receipts != expected_common:
        missing = sorted(expected_common - validated_gate_receipts)
        raise ManifestError(
            "aggregate gates lack semantic all-pass evidence receipts: " + ", ".join(missing)
        )
    return list(discovered_git_assets.values())


def verify_reader_acceptance(
    repo: Path,
    commit: str,
    aggregate_receipt: dict[str, Any],
    schema: dict[str, Any],
) -> None:
    raw = git_file_bytes(repo, commit, READER_RECEIPT_PATH.as_posix())
    require_gate_evidence(
        aggregate_receipt,
        "reader_ai_six_question_acceptance",
        sha256_bytes(raw),
        READER_RECEIPT_PATH.as_posix(),
    )
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ManifestError("reader acceptance receipt is not JSON") from exc
    validate_document(receipt, schema, "reader acceptance receipt")
    gate_evidence = aggregate_receipt["gates"]["reader_ai_six_question_acceptance"][
        "evidence"
    ]
    receipt_evidence = [
        item
        for item in gate_evidence
        if item.get("uri") == READER_RECEIPT_PATH.as_posix()
    ]
    if len(receipt_evidence) != 1:
        raise ManifestError("reader gate must bind exactly one detailed Git receipt")
    tested_commit = receipt.get("tested_git_commit")
    if not isinstance(tested_commit, str):
        raise ManifestError("reader acceptance lacks tested commit")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", tested_commit, commit],
        cwd=repo,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if ancestor.returncode != 0:
        raise ManifestError("reader-tested commit is not an ancestor of release commit")
    manifest_rows = receipt.get("package_source_manifest")
    if not isinstance(manifest_rows, list):
        raise ManifestError("reader package source manifest is absent")
    observed: dict[str, str] = {}
    for row in manifest_rows:
        if not isinstance(row, dict) or not isinstance(row.get("path"), str):
            raise ManifestError("reader package source manifest has invalid row")
        path = normalize_repository_path(row["path"])
        if path in observed:
            raise ManifestError(f"reader package source manifest duplicates {path}")
        package_path = f"{HANDOFF_ROOT.as_posix()}/{path}"
        tested_raw = git_file_bytes(repo, tested_commit, package_path)
        release_raw = git_file_bytes(repo, commit, package_path)
        if tested_raw != release_raw or row.get("sha256") != sha256_bytes(tested_raw):
            raise ManifestError(f"reader source bytes drifted or hash mismatch: {path}")
        observed[path] = row["sha256"]
    if set(observed) != READER_CORE_PATHS:
        raise ManifestError("reader receipt must hash-bind the exact 15-file core package")
    questions = receipt.get("questions", [])
    if [item.get("question_id") for item in questions if isinstance(item, dict)] != list(
        READER_QUESTION_IDS
    ):
        raise ManifestError("reader receipt does not contain the canonical six questions")
    promotions = receipt.get("prohibited_promotions", [])
    if [item.get("check_id") for item in promotions if isinstance(item, dict)] != list(
        READER_PROMOTION_IDS
    ):
        raise ManifestError("reader receipt does not contain the canonical eight anti-promotions")
    if receipt.get("verdict") != "PASS" or any(
        item.get("status") != "PASS" for item in questions + promotions
    ):
        raise ManifestError("reader acceptance receipt is not a complete PASS")
    if receipt.get("evaluator_context") != "NO_PRIOR_CONVERSATION_PACKAGE_ONLY":
        raise ManifestError("reader was not evaluated package-only with no prior conversation")

    # The external gate receipt is independently validated later; the detailed
    # Git receipt is authoritative here and is always re-executed below.

    validator_bytes = git_file_bytes(repo, commit, READER_VALIDATOR_PATH.as_posix())
    with tempfile.TemporaryDirectory(prefix="ism-reader-validator-") as temporary:
        validator_path = Path(temporary) / "validate_reader_acceptance.py"
        validator_path.write_bytes(validator_bytes)
        execution = subprocess.run(
            [
                sys.executable,
                str(validator_path),
                "--repo-root",
                str(repo),
                "--receipt",
                str(repo / READER_RECEIPT_PATH),
            ],
            cwd=repo,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    try:
        validator_summary = json.loads(execution.stdout)
    except json.JSONDecodeError as exc:
        raise ManifestError("commit-pinned reader validator emitted invalid JSON") from exc
    if (
        execution.returncode != 0
        or validator_summary.get("verdict") != "PASS"
        or validator_summary.get("questions") != 6
        or validator_summary.get("prohibited_promotions") != 8
        or validator_summary.get("required_source_files") != 15
        or validator_summary.get("source_files") != 15
    ):
        raise ManifestError(
            "commit-pinned reader validator did not PASS 6 questions, 8 anti-promotions, and 15 files: "
            + execution.stdout[:1000]
            + execution.stderr[:1000]
        )

    combined = "\n".join(
        str(item.get("answer", "")) + " " + str(item.get("explanation", ""))
        for item in questions + promotions
    ).lower()
    required_claims = (
        "confirmed cellular subclone",
        "confirmed linear ancestry",
        "88.2579",
        "methylation",
        "cn/loh",
        "longphase",
        "longlineage",
        "intersubmod",
        "python",
        "bip7",
        "bip8",
    )
    if any(token not in combined for token in required_claims):
        raise ManifestError("reader answers omit a required role, limit, host, or claim boundary")
    forbidden = (
        r"confirmed cellular subclones?\s*(?:=|:|is|are)\s*[1-9]",
        r"confirmed linear ancestry\s*(?:=|:|is|are)\s*[1-9]",
        r"88\.2579%?\s*(?:is|=|means|represents).{0,30}(?:accuracy|prevalence)",
    )
    if any(re.search(pattern, combined) for pattern in forbidden):
        raise ManifestError("reader receipt promotes a bounded research claim")


def verify_algorithm_crosswalk(repo: Path, commit: str) -> None:
    path = HANDOFF_ROOT / "evidence/algorithm_cli_claim_crosswalk.json"
    try:
        crosswalk = json.loads(git_file_bytes(repo, commit, path.as_posix()))
    except json.JSONDecodeError as exc:
        raise ManifestError("algorithm/CLI/claim crosswalk is not JSON") from exc
    gates = crosswalk.get("gates")
    expected = {
        "CROSSWALK_COMPLETENESS",
        "PUBLICATION_ASSET_READY",
        "RELEASE_ASSET_READY",
    }
    if not isinstance(gates, dict) or set(gates) != expected:
        raise ManifestError("algorithm crosswalk gate set is not exact")
    if gates["CROSSWALK_COMPLETENESS"].get("status") != "PASS":
        raise ManifestError("algorithm crosswalk is incomplete")
    for gate_id in ("PUBLICATION_ASSET_READY", "RELEASE_ASSET_READY"):
        if gates[gate_id].get("status") != "PASS":
            raise ManifestError(f"algorithm crosswalk static asset gate is not PASS: {gate_id}")
    records = crosswalk.get("records")
    if not isinstance(records, list) or len(records) != 35:
        raise ManifestError("algorithm crosswalk must contain exactly 35 rows")
    expected_ids = {f"ALG-{index:03d}" for index in range(1, 36)}
    observed_ids = {row.get("algorithm_claim_id") for row in records if isinstance(row, dict)}
    if observed_ids != expected_ids:
        raise ManifestError("algorithm crosswalk IDs are not exactly ALG-001 through ALG-035")
    allowed_dispositions = {"CONFIRMED", "CONFIRMED_WITH_LIMITS", "UNVERIFIED"}
    if any(row.get("current_disposition") not in allowed_dispositions for row in records):
        raise ManifestError("algorithm crosswalk contains an invalid claim disposition")
    if crosswalk.get("source_state") not in (
        "COMMIT_HASH_BOUND",
        "GIT_COMMIT_BOUND",
        "COMMIT_PINNED",
    ):
        raise ManifestError("algorithm crosswalk source-state marker is invalid")


def parse_generated_at(value: str | None) -> str:
    if value is None:
        return datetime.now().astimezone().isoformat(timespec="seconds")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ManifestError(f"invalid --generated-at: {value}") from exc
    if parsed.tzinfo is None:
        raise ManifestError("--generated-at must include a timezone")
    return parsed.isoformat(timespec="seconds")


def verify_authority(authority: dict[str, Any]) -> dict[str, Any]:
    scope = authority.get("scope", {})
    datasets = scope.get("technical_datasets")
    chromosomes = scope.get("chromosomes")
    if not isinstance(datasets, list) or len(datasets) != 7:
        raise ManifestError("authority scope must contain exactly 7 technical datasets")
    if scope.get("biological_ids") != 6:
        raise ManifestError("authority scope must contain exactly 6 biological IDs")
    if chromosomes != [f"chr{index}" for index in range(1, 23)]:
        raise ManifestError("authority scope must be chr1-chr22 in order")

    denominators = authority.get("denominators", {})
    graph_shape = denominators.get("one_rooted_unlabeled_topology", {})
    if graph_shape.get("count") != 63506 or graph_shape.get("denominator") != 71955:
        raise ManifestError("authority graph-shape numerator/denominator drifted")
    if graph_shape.get("percent") != 88.2579:
        raise ManifestError("authority model-conditional graph-shape percent drifted")
    if authority.get("model_contract", {}).get("methyl_role") != "pattern-conditioned association-only sidecar":
        raise ManifestError("authority methylation role is not association-only")

    forbidden = " ".join(authority.get("claim_boundary", {}).get("forbidden", [])).lower()
    for phrase in ("confirmed cellular clone", "confirmed biological mutation ancestry"):
        if phrase not in forbidden:
            raise ManifestError(f"authority forbidden-claim boundary missing: {phrase}")
    limitations = " ".join(authority.get("known_limitations", [])).lower()
    if "copy-neutral loh" not in limitations or "not integrated" not in limitations:
        raise ManifestError("authority CN/LOH non-integration limit is missing")

    return {
        "technical_datasets": datasets,
        "biological_ids": scope["biological_ids"],
        "chromosomes": chromosomes,
    }


def verify_denominator_registry(raw: bytes) -> None:
    try:
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8")), delimiter="\t")
        rows = list(reader)
    except UnicodeDecodeError as exc:
        raise ManifestError("frozen denominator registry is not UTF-8") from exc
    required_columns = {
        "metric",
        "numerator",
        "denominator",
        "percent",
        "status",
        "authority_artifact_id",
    }
    if not reader.fieldnames or not required_columns.issubset(reader.fieldnames):
        raise ManifestError("frozen denominator registry columns drifted")
    by_metric = {row["metric"]: row for row in rows}
    if len(by_metric) != len(rows):
        raise ManifestError("frozen denominator registry contains duplicate metric IDs")
    for metric, expected in EXPECTED_DENOMINATORS.items():
        row = by_metric.get(metric)
        if row is None:
            raise ManifestError(f"frozen denominator registry is missing {metric}")
        observed = (
            row["numerator"],
            row["denominator"],
            row["percent"],
            row["authority_artifact_id"],
        )
        if observed != expected or row["status"] != "current":
            raise ManifestError(
                f"frozen denominator registry drifted for {metric}: {observed} status={row['status']}"
            )


def authority_replay_expected_rows(
    authority: dict[str, Any],
) -> set[tuple[str, str, str, str]]:
    expected: set[tuple[str, str, str, str]] = set()
    for artifact in authority.get("artifacts", []):
        expected.add(
            (
                "artifact",
                str(artifact["artifact_id"]),
                str(artifact["path"]),
                str(artifact["sha256"]),
            )
        )
    implementation = authority.get("implementation", {})
    binary = implementation.get("frozen_binary", {})
    expected.add(
        (
            "binary",
            "frozen_binary",
            str(binary["path"]),
            str(binary["sha256"]),
        )
    )
    for index, snapshot in enumerate(implementation.get("source_snapshots", []), start=1):
        expected.add(
            (
                "source_snapshot",
                f"source_snapshot_{index}",
                str(snapshot["path"]),
                str(snapshot["sha256_at_handoff"]),
            )
        )
    if len(expected) != 19:
        raise ManifestError("frozen authority does not define the exact 13+1+5 replay set")
    return expected


def verify_authority_replay(raw: bytes, authority: dict[str, Any]) -> dict[str, Any]:
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ManifestError("authority replay receipt is not JSON") from exc
    expected_tally = {"MATCH": 19, "MISSING": 0, "HASH_MISMATCH": 0}
    if receipt.get("receipt_type") != "frozen_authority_hash_replay":
        raise ManifestError("authority replay receipt type mismatch")
    if receipt.get("manifest_sha256") != FROZEN_AUTHORITY_SHA256:
        raise ManifestError("authority replay receipt does not bind the frozen authority bytes")
    if receipt.get("authority_as_of_date") != "2026-08-01":
        raise ManifestError("authority replay receipt date mismatch")
    if receipt.get("pass") is not True or receipt.get("total") != 19:
        raise ManifestError("authority replay receipt is not a 19-record PASS")
    if receipt.get("tally") != expected_tally:
        raise ManifestError("authority replay receipt tally is not 19 MATCH / 0 missing / 0 mismatch")
    results = receipt.get("results")
    if not isinstance(results, list) or len(results) != 19:
        raise ManifestError("authority replay receipt must contain exactly 19 result records")
    observed_rows: set[tuple[str, str, str, str]] = set()
    for result in results:
        if not isinstance(result, dict):
            raise ManifestError("authority replay result must be an object")
        identity = (str(result.get("record_kind")), str(result.get("artifact_id")))
        if result.get("status") != "MATCH":
            raise ManifestError(f"authority replay result is not MATCH: {identity}")
        expected_sha = result.get("expected_sha256")
        actual_sha = result.get("actual_sha256")
        if not isinstance(expected_sha, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_sha):
            raise ManifestError(f"authority replay expected hash is invalid: {identity}")
        if expected_sha != actual_sha:
            raise ManifestError(f"authority replay expected/actual hash mismatch: {identity}")
        observed_rows.add(
            (
                identity[0],
                identity[1],
                str(result.get("path")),
                expected_sha,
            )
        )
    if observed_rows != authority_replay_expected_rows(authority):
        raise ManifestError(
            "authority replay rows do not exactly match the frozen 13 artifacts + binary + 5 source snapshots"
        )
    return receipt


def verify_aggregate_gate_bindings(receipt: dict[str, Any], commit: str) -> None:
    gates = receipt.get("gates")
    if not isinstance(gates, dict) or set(gates) != set(REQUIRED_AGGREGATE_GATES):
        raise ManifestError("aggregate acceptance receipt gate set is not exact")
    for gate_id in REQUIRED_AGGREGATE_GATES:
        gate = gates[gate_id]
        if gate.get("source_commit") != commit:
            raise ManifestError(f"aggregate gate does not bind source commit: {gate_id}")
        evidence = gate.get("evidence", [])
        expected_type = REQUIRED_GATE_EVIDENCE_TYPES[gate_id]
        if not any(item.get("evidence_type") == expected_type for item in evidence):
            raise ManifestError(f"aggregate gate lacks {expected_type} evidence: {gate_id}")
        for item in evidence:
            if item.get("subject_commit") != commit:
                raise ManifestError(f"aggregate gate evidence does not bind source commit: {gate_id}")


def require_gate_evidence(
    receipt: dict[str, Any], gate_id: str, sha256: str, uri: str
) -> None:
    evidence = receipt["gates"][gate_id]["evidence"]
    if not any(item.get("sha256") == sha256 and item.get("uri") == uri for item in evidence):
        raise ManifestError(f"aggregate gate does not hash-bind required receipt: {gate_id}")


def require_git_gate_receipt(
    repo: Path,
    commit: str,
    receipt: dict[str, Any],
    gate_id: str,
    repository_path: Path,
) -> dict[str, Any]:
    raw = git_file_bytes(repo, commit, repository_path.as_posix())
    require_gate_evidence(
        receipt,
        gate_id,
        sha256_bytes(raw),
        repository_path.as_posix(),
    )
    try:
        document = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ManifestError(f"required gate receipt is not JSON: {repository_path}") from exc
    if not isinstance(document, dict):
        raise ManifestError(f"required gate receipt is not an object: {repository_path}")
    return document


def verify_repo_hygiene_receipt(receipt: dict[str, Any]) -> None:
    if receipt.get("receipt_type") != "repo_hygiene_package_summary" or receipt.get("pass") is not True:
        raise ManifestError("repo hygiene receipt is not a PASS summary")
    for key in (
        "versioned_worktree_absolute_symlink_count",
        "versioned_worktree_broken_symlink_count",
    ):
        if receipt.get(key) != 0:
            raise ManifestError(f"repo hygiene receipt has nonzero {key}")
    if receipt.get("tracked_settings_local_present") is not False:
        raise ManifestError("repo hygiene receipt still has tracked local settings")
    if receipt.get("archive", {}).get("pass") is not True:
        raise ManifestError("repo hygiene archive verification is not PASS")


def verify_repo_tree_hygiene(repo: Path, commit: str) -> None:
    raw = run_git(repo, "ls-tree", "-r", "-z", commit, text=False)
    tree_paths: set[str] = set()
    symlinks: list[tuple[str, str]] = []
    for record in (item for item in raw.split(b"\0") if item):
        metadata, encoded_path = record.split(b"\t", 1)
        mode, object_type, _ = metadata.decode("ascii").split()
        path = encoded_path.decode("utf-8")
        tree_paths.add(path)
        lowered = path.lower()
        if lowered.endswith("settings.local.json") or "/settings.local." in lowered:
            raise ManifestError(f"tracked local settings remain in release tree: {path}")
        if mode == "120000" and object_type == "blob":
            target = git_file_bytes(repo, commit, path).decode("utf-8", errors="strict")
            symlinks.append((path, target))
    symlink_targets = dict(symlinks)

    def resolve_link(origin: str) -> str:
        pending = list(PurePosixPath(symlink_targets[origin]).parts)
        resolved: list[str] = list(PurePosixPath(origin).parent.parts)
        if PurePosixPath(symlink_targets[origin]).is_absolute():
            raise ManifestError(
                f"absolute symlink remains in release tree: {origin} -> {symlink_targets[origin]}"
            )
        seen: set[tuple[str, tuple[str, ...]]] = set()
        while pending:
            part = pending.pop(0)
            if part in ("", "."):
                continue
            if part == "..":
                if not resolved:
                    raise ManifestError(f"symlink escapes release tree: {origin}")
                resolved.pop()
                continue
            resolved.append(part)
            candidate = PurePosixPath(*resolved).as_posix()
            target = symlink_targets.get(candidate)
            if target is None:
                continue
            if PurePosixPath(target).is_absolute():
                raise ManifestError(f"absolute symlink remains in release tree: {candidate} -> {target}")
            state = (candidate, tuple(pending))
            if state in seen or len(seen) > len(symlink_targets) * 4 + 16:
                raise ManifestError(f"symlink cycle remains in release tree: {origin}")
            seen.add(state)
            resolved.pop()
            pending = list(PurePosixPath(target).parts) + pending

        return PurePosixPath(*resolved).as_posix()

    for path in symlink_targets:
        resolved_path = resolve_link(path)
        if resolved_path not in tree_paths and not any(
            candidate.startswith(resolved_path.rstrip("/") + "/") for candidate in tree_paths
        ):
            raise ManifestError(f"broken symlink chain remains in release tree: {path}")


def verify_tracked_large_asset_policy(
    registry: dict[str, Any], migration_receipt: dict[str, Any], commit: str
) -> None:
    summary = registry.get("summary", {})
    if summary.get("verdict") != "PASS":
        raise ManifestError(
            "tracked-large asset policy is not PASS: " + str(summary.get("verdict"))
        )
    if summary.get("tracked_assets_over_threshold") != 0 or registry.get("assets") != []:
        raise ManifestError("tracked-large asset policy still contains unadjudicated assets")
    if migration_receipt.get("remaining_policy_status") != "PASS":
        raise ManifestError("large-asset migration remaining policy is not PASS")
    if migration_receipt.get("verdict") not in (
        "PASS",
        "GITHUB_RELEASE_ROUNDTRIP_VERIFIED",
    ):
        raise ManifestError("large-asset migration receipt is not publication-complete")


def verify_actual_tracked_large_assets(repo: Path, commit: str) -> None:
    raw = str(run_git(repo, "ls-tree", "-r", "-l", commit))
    offenders: list[str] = []
    for line in raw.splitlines():
        metadata, path = line.split("\t", 1)
        fields = metadata.split()
        if len(fields) != 4 or fields[1] != "blob" or fields[3] == "-":
            continue
        if int(fields[3]) > 1024 * 1024:
            offenders.append(f"{fields[3]}:{path}")
    if offenders:
        raise ManifestError(
            "tracked-large asset policy recomputation found files over 1 MiB: "
            + "; ".join(offenders[:10])
        )


def verify_project_license(
    repo: Path, commit: str, expected_spdx: str
) -> None:
    try:
        import yaml
    except ImportError as exc:
        raise ManifestError("PyYAML is required to validate CITATION.cff") from exc
    citation_raw = git_file_bytes(repo, commit, "CITATION.cff")
    try:
        citation = yaml.safe_load(citation_raw)
    except yaml.YAMLError as exc:
        raise ManifestError("CITATION.cff is not valid YAML") from exc
    if not isinstance(citation, dict) or citation.get("license") != expected_spdx:
        raise ManifestError(
            "CITATION.cff does not contain the maintainer-confirmed project SPDX expression"
        )

    license_text = git_file_bytes(repo, commit, "LICENSE").decode("utf-8", errors="replace")
    if "GNU GENERAL PUBLIC LICENSE" not in license_text or "Version 3, 29 June 2007" not in license_text:
        raise ManifestError("LICENSE does not contain the GNU GPL version 3 license text")

    grep = subprocess.run(
        [
            "git",
            "grep",
            "-I",
            "-n",
            "-E",
            r"SPDX-License-Identifier: GPL-3\.0-(only|or-later)",
            commit,
            "--",
        ],
        cwd=repo,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if grep.returncode not in (0, 1):
        raise ManifestError(f"cannot audit source SPDX headers: {grep.stderr.strip()}")
    conflicting_headers = [
        line for line in grep.stdout.splitlines() if expected_spdx not in line
    ]
    if conflicting_headers:
        raise ManifestError(
            "source SPDX headers conflict with maintainer-confirmed project SPDX: "
            + "; ".join(conflicting_headers[:5])
        )

    registry_raw = git_file_bytes(repo, commit, ARTIFACT_REGISTRY_PATH.as_posix())
    artifacts = load_artifact_registry(registry_raw).values()
    conflicting: set[str] = set()
    for artifact in artifacts:
        license_value = str(artifact.get("license", ""))
        for expression in ("GPL-3.0-only", "GPL-3.0-or-later"):
            if expression in license_value and expression != expected_spdx:
                conflicting.add(expression)
    if conflicting:
        raise ManifestError(
            "artifact registry conflicts with maintainer-confirmed project SPDX: "
            + ", ".join(sorted(conflicting))
        )


def external_receipt(
    repo: Path,
    commit: str,
    spec: str,
    expected_type: str,
    schema: dict[str, Any],
) -> tuple[dict[str, str], dict[str, Any], dict[str, Any]]:
    _, source_path = parse_release_asset(spec)
    asset = build_release_asset(repo, spec)
    raw = source_path.read_bytes()
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ManifestError(f"receipt is not JSON: {source_path}") from exc
    validate_document(receipt, schema, f"receipt {source_path}")
    if receipt.get("receipt_type") != expected_type:
        raise ManifestError(f"receipt type mismatch: {source_path}")
    if receipt.get("source_commit") != commit:
        raise ManifestError(f"receipt does not bind source commit: {source_path}")
    if receipt.get("overall_pass") is not True or receipt.get("blocking_gates") != []:
        raise ManifestError(f"receipt is not an all-pass receipt: {source_path}")

    record = {"published_path": asset["published_path"], "sha256": sha256_bytes(raw)}
    return record, receipt, asset


def verify_longlineage_safety_bindings(
    repo: Path,
    commit: str,
    receipt: dict[str, Any],
    longlineage_repo: Path,
) -> None:
    for gate_id, paths in LONGLINEAGE_SUPPORT_PATHS.items():
        evidence = receipt["preview_safety_gates"][gate_id]["evidence"]
        for path in paths:
            raw = git_file_bytes(repo, commit, path.as_posix())
            if not any(
                item.get("uri") == path.as_posix()
                and item.get("sha256") == sha256_bytes(raw)
                for item in evidence
            ):
                raise ManifestError(
                    f"LongLineage safety gate does not bind supporting bytes: {gate_id} {path}"
                )

    source_manifest = json.loads(
        git_file_bytes(repo, commit, LONGLINEAGE_SOURCE_MANIFEST_PATH.as_posix())
    )
    mappings = source_manifest.get("mappings")
    if not isinstance(mappings, list) or len(mappings) != 21:
        raise ManifestError("LongLineage source-to-target manifest is not exactly 21 rows")
    safety_stacked_commit = receipt.get("safety_stacked_commit")
    if (
        not re.fullmatch(r"[0-9a-f]{40}", str(safety_stacked_commit))
        or source_manifest.get("license_review_status") != "APPROVED_FOR_PUBLIC_RELEASE"
    ):
        raise ManifestError("LongLineage source manifest is not bound to the approved safety stack")
    mapping_ids: set[str] = set()
    for row in mappings:
        if not isinstance(row, dict) or not isinstance(row.get("origin_id"), str):
            raise ManifestError("LongLineage source-to-target manifest has invalid row")
        mapping_id = row["origin_id"]
        if mapping_id in mapping_ids:
            raise ManifestError(f"LongLineage source mapping ID is duplicated: {mapping_id}")
        mapping_ids.add(mapping_id)
        if row.get("source_replay_status") not in (
            "MATCHED_DECLARED_ORIGIN_COMMIT",
            "MATCHED_OTHER_COMMIT",
        ):
            raise ManifestError(f"LongLineage source mapping is unresolved: {mapping_id}")
        if row.get("license_disposition") != "APPROVED_FOR_PUBLIC_RELEASE":
            raise ManifestError(f"LongLineage mapping license is not approved: {mapping_id}")
        for key in ("source_sha256", "target_sha256"):
            if not re.fullmatch(r"[0-9a-f]{64}", str(row.get(key))):
                raise ManifestError(f"LongLineage mapping has invalid {key}: {mapping_id}")
        if not re.fullmatch(r"[0-9a-f]{40}", str(row.get("source_commit"))):
            raise ManifestError(f"LongLineage mapping lacks origin commit: {mapping_id}")
        if not isinstance(row.get("target"), str) or not row["target"]:
            raise ManifestError(f"LongLineage mapping lacks target path: {mapping_id}")
        if not isinstance(row.get("license_evidence"), str) or not row["license_evidence"]:
            raise ManifestError(f"LongLineage mapping lacks license evidence: {mapping_id}")

    gate_script = git_file_bytes(repo, commit, LONGLINEAGE_PUBLIC_GATE_SCRIPT_PATH.as_posix())
    gate_output = git_file_bytes(repo, commit, LONGLINEAGE_PUBLIC_GATE_OUTPUT_PATH.as_posix())
    gate_receipt = json.loads(
        git_file_bytes(repo, commit, LONGLINEAGE_PUBLIC_GATE_RECEIPT_PATH.as_posix())
    )
    if (
        gate_receipt.get("candidate_commit") != LONGLINEAGE_CANDIDATE
        or gate_receipt.get("safety_stacked_commit") != safety_stacked_commit
        or gate_receipt.get("repository_url") != "https://github.com/liaoyoyo/LongLineage"
        or not re.fullmatch(r"[0-9a-f]{40}", str(gate_receipt.get("safety_tree")))
        or gate_receipt.get("cross_repo_verification_method")
        not in ("GIT_SHOW_EXACT_HEAD", "GITHUB_RAW_API_EXACT_HEAD")
        or gate_receipt.get("source_manifest_sha256")
        != sha256_bytes(
            git_file_bytes(repo, commit, LONGLINEAGE_SOURCE_MANIFEST_PATH.as_posix())
        )
        or gate_receipt.get("script_sha256") != sha256_bytes(gate_script)
        or gate_receipt.get("output_sha256") != sha256_bytes(gate_output)
        or gate_receipt.get("exit_code") != 0
        or gate_receipt.get("ctest_total") != 49
        or gate_receipt.get("ctest_passed") != 49
        or gate_receipt.get("ctest_failed") != 0
        or gate_receipt.get("strict_public_preview_gate_exit") != 0
        or gate_receipt.get("production_run_exit") != 6
        or gate_receipt.get("production_run_status") != "KernelBlocked"
        or gate_receipt.get("phase_gates")
        != {gate_id: "BLOCKED" for gate_id in LONGLINEAGE_GATES}
        or gate_receipt.get("overall_pass") is not True
        or gate_receipt.get("blocking_gates") != []
    ):
        raise ManifestError("LongLineage public-preview gate script/output receipt is invalid")

    support_bindings = gate_receipt.get("support_asset_bindings")
    support_paths = {
        LONGLINEAGE_SOURCE_MANIFEST_PATH.as_posix(),
        LONGLINEAGE_PUBLIC_GATE_SCRIPT_PATH.as_posix(),
        LONGLINEAGE_PUBLIC_GATE_OUTPUT_PATH.as_posix(),
        LONGLINEAGE_HISTORY_OUTPUT_PATH.as_posix(),
        (HANDOFF_ROOT / "evidence/longlineage_SBOM.spdx.json").as_posix(),
        (HANDOFF_ROOT / "evidence/longlineage_third_party_notices.md").as_posix(),
    }
    if not isinstance(support_bindings, list) or {
        item.get("ism_path") for item in support_bindings if isinstance(item, dict)
    } != support_paths:
        raise ManifestError("LongLineage cross-repository support binding set is not exact")
    actual_ll_root = Path(
        str(run_git(longlineage_repo, "rev-parse", "--show-toplevel")).strip()
    ).resolve()
    if actual_ll_root != longlineage_repo.resolve():
        raise ManifestError("--longlineage-repo-root must be the LongLineage Git top level")
    origin_url = str(
        run_git(longlineage_repo, "config", "--get", "remote.origin.url")
    ).strip()
    if origin_url.removesuffix(".git") != "https://github.com/liaoyoyo/LongLineage":
        raise ManifestError("LongLineage repository origin URL is not the governed public repository")
    advertised_head = gate_receipt.get("remote_advertised_head")
    advertised_observed_at = gate_receipt.get("remote_advertised_observed_at")
    remote_branch = gate_receipt.get("remote_branch")
    if advertised_head != safety_stacked_commit:
        raise ManifestError("LongLineage safety head is not bound to the advertised remote head")
    if (
        not isinstance(remote_branch, str)
        or not re.fullmatch(r"[A-Za-z0-9._/-]+", remote_branch)
        or ".." in remote_branch
        or remote_branch.startswith("/")
    ):
        raise ManifestError("LongLineage remote branch name is invalid")
    try:
        parse_generated_at(advertised_observed_at)
    except ManifestError as exc:
        raise ManifestError("LongLineage remote-head observation timestamp is invalid") from exc
    try:
        ls_remote = subprocess.run(
            ["git", "ls-remote", "origin", f"refs/heads/{remote_branch}"],
            cwd=longlineage_repo,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ManifestError("LongLineage remote-head reachability query failed") from exc
    if ls_remote.returncode != 0 or not any(
        line == f"{safety_stacked_commit}\trefs/heads/{remote_branch}"
        for line in ls_remote.stdout.splitlines()
        if "\t" in line
    ):
        raise ManifestError("LongLineage safety head is not advertised by an origin branch")
    try:
        resolved_safety_commit = str(
            run_git(
                longlineage_repo,
                "rev-parse",
                "--verify",
                f"{safety_stacked_commit}^{{commit}}",
            )
        ).strip()
    except ManifestError as exc:
        raise ManifestError("LongLineage safety head is not reachable in the supplied repository") from exc
    if resolved_safety_commit != safety_stacked_commit:
        raise ManifestError("LongLineage safety head did not resolve exactly")
    actual_safety_tree = str(
        run_git(longlineage_repo, "rev-parse", f"{safety_stacked_commit}^{{tree}}")
    ).strip()
    if actual_safety_tree != gate_receipt.get("safety_tree"):
        raise ManifestError("LongLineage safety tree does not match cross-repository receipt")
    for item in support_bindings:
        if (
            item.get("repository_url") != "https://github.com/liaoyoyo/LongLineage"
            or item.get("head_commit") != safety_stacked_commit
            or not isinstance(item.get("longlineage_path"), str)
            or not item["longlineage_path"]
            or item.get("sha256")
            != sha256_bytes(git_file_bytes(repo, commit, item["ism_path"]))
        ):
            raise ManifestError("LongLineage cross-repository support asset binding is invalid")
        ll_raw = git_file_bytes(
            longlineage_repo, safety_stacked_commit, item["longlineage_path"]
        )
        ism_raw = git_file_bytes(repo, commit, item["ism_path"])
        if ll_raw != ism_raw or item["sha256"] != sha256_bytes(ll_raw):
            raise ManifestError(
                "LongLineage support copy differs from exact safety-head Git bytes: "
                + item["longlineage_path"]
            )

    history_output = git_file_bytes(repo, commit, LONGLINEAGE_HISTORY_OUTPUT_PATH.as_posix())
    history_receipt = json.loads(
        git_file_bytes(repo, commit, LONGLINEAGE_HISTORY_RECEIPT_PATH.as_posix())
    )
    if (
        history_receipt.get("base_commit") != LONGLINEAGE_CANDIDATE
        or history_receipt.get("head_commit") != safety_stacked_commit
        or history_receipt.get("scope") != "FULL_HISTORY"
        or history_receipt.get("finding_count") != 0
        or history_receipt.get("findings") != []
        or history_receipt.get("raw_output_sha256") != sha256_bytes(history_output)
        or history_receipt.get("overall_pass") is not True
    ):
        raise ManifestError("LongLineage full-history scan raw receipt is invalid")

    safety_path = HANDOFF_ROOT / "evidence/longlineage_public_safety_receipt.json"
    safety = json.loads(git_file_bytes(repo, commit, safety_path.as_posix()))
    if safety.get("candidate_commit") != LONGLINEAGE_CANDIDATE:
        raise ManifestError("LongLineage public-safety evidence candidate mismatch")
    if safety.get("safety_stacked_commit") != safety_stacked_commit:
        raise ManifestError("LongLineage public-safety summary is bound to another stack")
    if safety.get("gate_status") not in ("PASS", "PUBLIC_PREVIEW_SAFETY_PASS"):
        raise ManifestError("LongLineage public-safety evidence remains fail-closed")
    if safety.get("required_verdict") not in ("PASS", "PUBLIC_PREVIEW_SAFETY_PASS"):
        raise ManifestError("LongLineage public-safety verdict does not authorize a preview")
    history = safety.get("history_hygiene", {})
    if history.get("status") != "PASS" or history.get("finding_count") != 0:
        raise ManifestError("LongLineage history-hygiene evidence is not PASS/zero-findings")
    source_replay = safety.get("source_replay", {})
    if (
        source_replay.get("mapping_rows") != 21
        or source_replay.get("unresolved_rows") != 0
        or source_replay.get("target_digest_matches") != 21
        or source_replay.get("target_digest_mismatches") != 0
        or source_replay.get("missing_unique_hashes") not in (None, [])
    ):
        raise ManifestError("LongLineage 21-row source-origin mapping is not fully resolved")
    license_review = safety.get("license_review", {})
    if (
        license_review.get("row_dispositions_pending") != 0
        or license_review.get("row_dispositions_approved") != 21
        or license_review.get("repository_status") not in ("APPROVED", "PASS")
        or license_review.get("third_party_notice_status") not in ("APPROVED", "PASS")
        or license_review.get("sbom_status") not in ("APPROVED", "PASS")
    ):
        raise ManifestError("LongLineage license/source/dependency audit is not fully approved")

    foundation_path = HANDOFF_ROOT / "evidence/longlineage_clean_foundation_receipt.json"
    foundation = json.loads(git_file_bytes(repo, commit, foundation_path.as_posix()))
    if foundation.get("candidate_baseline") != LONGLINEAGE_CANDIDATE:
        raise ManifestError("LongLineage clean-foundation evidence candidate mismatch")
    if foundation.get("validated_stacked_commit") != LONGLINEAGE_HISTORICAL_SAFE_FAILURE_STACK:
        raise ManifestError("LongLineage historical safe-failure foundation commit drifted")
    results = foundation.get("results", {})
    if (
        results.get("configure_exit") != 0
        or results.get("build_exit") != 0
        or results.get("ctest_total") != 49
        or results.get("ctest_passed") != 49
        or results.get("ctest_failed") != 0
        or results.get("strict_public_preview_gate_exit") != 1
        or results.get("strict_public_preview_gate_interpretation") != "EXPECTED_SAFE_FAILURE"
        or results.get("production_run_expected_exit") != 6
        or results.get("production_run_expected_status") != "KernelBlocked"
    ):
        raise ManifestError("LongLineage clean-foundation fail-closed contract is not verified")

    sbom_path = HANDOFF_ROOT / "evidence/longlineage_SBOM.spdx.json"
    sbom = json.loads(git_file_bytes(repo, commit, sbom_path.as_posix()))
    packages = sbom.get("packages")
    if sbom.get("spdxVersion") != "SPDX-2.3" or not isinstance(packages, list) or not packages:
        raise ManifestError("LongLineage SBOM is not a valid SPDX-2.3 package inventory")
    for package in packages:
        concluded = package.get("licenseConcluded") if isinstance(package, dict) else None
        if concluded in (None, "", "NOASSERTION", "NONE"):
            raise ManifestError("LongLineage SBOM contains unresolved package licensing")
    mapping_packages = {
        str(package.get("SPDXID", "")).removeprefix("SPDXRef-Origin-"): package.get(
            "licenseConcluded"
        )
        for package in packages
        if isinstance(package, dict)
        and str(package.get("SPDXID", "")).startswith("SPDXRef-Origin-")
    }
    if set(mapping_packages) != mapping_ids or any(
        value in (None, "", "NOASSERTION", "NONE") for value in mapping_packages.values()
    ):
        raise ManifestError("LongLineage SBOM does not exactly crosswalk the 21 source mappings")

    notices_path = HANDOFF_ROOT / "evidence/longlineage_third_party_notices.md"
    notices = git_file_bytes(repo, commit, notices_path.as_posix()).decode("utf-8")
    normalized_notices = notices.upper()
    if "STATUS: APPROVED_FOR_PUBLIC_RELEASE" not in normalized_notices or any(
        marker in normalized_notices
        for marker in (
            "PENDING",
            "NOT APPROVAL",
            "INVENTORY_ONLY_NOT_APPROVAL",
            "NOASSERTION",
        )
    ):
        raise ManifestError("LongLineage third-party notices are not approval-ready")


def ensure_output_outside_repository(repo: Path, output: Path) -> Path:
    resolved = output.expanduser().resolve()
    try:
        resolved.relative_to(repo)
    except ValueError:
        pass
    else:
        raise ManifestError(
            f"generated manifest/checksums must stay outside Git to avoid self-reference: {resolved}"
        )
    if resolved.exists():
        raise ManifestError(f"refusing to overwrite existing output: {resolved}")
    return resolved


def validate_document(document: dict[str, Any], schema: dict[str, Any], label: str) -> None:
    try:
        from jsonschema import FormatChecker
        from jsonschema.validators import validator_for
    except ImportError as exc:
        raise ManifestError("jsonschema is required to validate the release manifest") from exc
    validator_class = validator_for(schema)
    validator_class.check_schema(schema)
    errors = sorted(
        validator_class(schema, format_checker=FormatChecker()).iter_errors(document),
        key=lambda error: list(error.absolute_path),
    )
    if errors:
        rendered = []
        for error in errors:
            location = "/".join(str(part) for part in error.absolute_path) or "<root>"
            rendered.append(f"{location}: {error.message}")
        raise ManifestError(f"{label} failed schema:\n" + "\n".join(rendered))


def write_text_exclusive(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(text)


def build_manifest(args: argparse.Namespace) -> tuple[dict[str, Any], list[str]]:
    repo = args.repo_root.expanduser().resolve()
    actual_root = Path(str(run_git(repo, "rev-parse", "--show-toplevel").strip())).resolve()
    if actual_root != repo:
        raise ManifestError(f"--repo-root must be the Git top level: {repo} != {actual_root}")
    dirty = str(run_git(repo, "status", "--porcelain", "--untracked-files=all")).strip()
    if dirty:
        sample = "\n".join(dirty.splitlines()[:8])
        raise ManifestError(f"repository is not clean; commit/stash first:\n{sample}")

    requested_revision = args.revision
    commit = str(run_git(repo, "rev-parse", "--verify", f"{requested_revision}^{{commit}}")).strip()
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ManifestError(f"resolved revision is not a full commit: {commit}")
    tree = str(run_git(repo, "rev-parse", f"{commit}^{{tree}}")).strip()
    generated_at = parse_generated_at(args.generated_at)

    authority_path = AUTHORITY_PATH.as_posix()
    denominator_path = DENOMINATOR_PATH.as_posix()
    schema_path = SCHEMA_PATH.as_posix()
    authority_raw = git_file_bytes(repo, commit, authority_path)
    denominator_raw = git_file_bytes(repo, commit, denominator_path)
    schema_raw = git_file_bytes(repo, commit, schema_path)
    aggregate_schema_raw = git_file_bytes(
        repo, commit, AGGREGATE_RECEIPT_SCHEMA_PATH.as_posix()
    )
    longlineage_schema_raw = git_file_bytes(
        repo, commit, LONGLINEAGE_SAFETY_SCHEMA_PATH.as_posix()
    )
    common_gate_schema_raw = git_file_bytes(
        repo, commit, COMMON_GATE_SCHEMA_PATH.as_posix()
    )
    reader_schema_raw = git_file_bytes(repo, commit, READER_SCHEMA_PATH.as_posix())
    authority_replay_raw = git_file_bytes(repo, commit, AUTHORITY_REPLAY_PATH.as_posix())
    artifact_registry_raw = git_file_bytes(repo, commit, ARTIFACT_REGISTRY_PATH.as_posix())
    if sha256_bytes(authority_raw) != FROZEN_AUTHORITY_SHA256:
        raise ManifestError("authority manifest bytes do not match the frozen 2026-08-01 hash")
    if sha256_bytes(denominator_raw) != FROZEN_DENOMINATOR_SHA256:
        raise ManifestError("denominator registry bytes do not match the frozen 2026-08-01 hash")
    authority = json.loads(authority_raw)
    schema = json.loads(schema_raw)
    aggregate_schema = json.loads(aggregate_schema_raw)
    longlineage_schema = json.loads(longlineage_schema_raw)
    common_gate_schema = json.loads(common_gate_schema_raw)
    reader_schema = json.loads(reader_schema_raw)
    scope = verify_authority(authority)
    verify_denominator_registry(denominator_raw)
    verify_authority_replay(authority_replay_raw, authority)
    artifact_registry = load_artifact_registry(artifact_registry_raw)

    blockers = sorted(set(args.blocker))
    if args.status == "BLOCKED":
        if not blockers:
            raise ManifestError("BLOCKED status requires at least one --blocker")
        if args.approval_receipt:
            raise ManifestError("BLOCKED status cannot carry aggregate PASS receipts")
        if args.project_license_spdx is None:
            if LICENSE_CONFIRMATION_BLOCKER not in blockers:
                raise ManifestError(
                    f"unconfirmed project license requires blocker {LICENSE_CONFIRMATION_BLOCKER}"
                )
            project_license_status = "NEEDS_MAINTAINER_CONFIRMATION"
        else:
            verify_project_license(repo, commit, args.project_license_spdx)
            project_license_status = "CONFIRMED"
        approvals: list[dict[str, str]] = []
        approval_assets: list[dict[str, Any]] = []
        aggregate_receipt = None
        release_ready = False
    else:
        if blockers:
            raise ManifestError("APPROVED_RESEARCH_HANDOFF cannot contain blockers")
        if len(args.approval_receipt) != 1:
            raise ManifestError("approved handoff requires exactly one --approval-receipt")
        if args.project_license_spdx is None:
            raise ManifestError("approved handoff requires --project-license-spdx")
        approvals = []
        approval_assets = []
        record, aggregate_receipt, asset = external_receipt(
            repo,
            commit,
            args.approval_receipt[0],
            expected_type="aggregate_release_acceptance",
            schema=aggregate_schema,
        )
        verify_aggregate_gate_bindings(aggregate_receipt, commit)
        if aggregate_receipt["project_license_spdx"] != args.project_license_spdx:
            raise ManifestError("aggregate receipt project license does not match CLI selection")
        verify_project_license(repo, commit, args.project_license_spdx)
        project_license_status = "CONFIRMED"
        approvals.append(record)
        approval_assets.append(asset)
        release_ready = True

    if args.longlineage_preview_safety_status == "BLOCKED":
        if args.longlineage_safety_receipt:
            raise ManifestError("blocked LongLineage preview safety cannot carry a PASS receipt")
        if args.status == "APPROVED_RESEARCH_HANDOFF":
            raise ManifestError("approved handoff requires LongLineage preview safety VERIFIED")
        preview_status = "NOT_READY"
        safety_receipt_record = None
        safety_receipt_assets: list[dict[str, Any]] = []
    else:
        if not args.longlineage_safety_receipt:
            raise ManifestError("verified LongLineage preview safety requires --longlineage-safety-receipt")
        if args.longlineage_repo_root is None:
            raise ManifestError("verified LongLineage preview safety requires --longlineage-repo-root")
        safety_receipt_record, safety_receipt, safety_asset = external_receipt(
            repo,
            commit,
            args.longlineage_safety_receipt,
            expected_type="longlineage_public_preview_safety",
            schema=longlineage_schema,
        )
        verify_longlineage_safety_bindings(
            repo,
            commit,
            safety_receipt,
            args.longlineage_repo_root.expanduser().resolve(),
        )
        preview_status = "RESEARCH_PREVIEW_SAFETY_VERIFIED"
        safety_receipt_assets = [safety_asset]

    if aggregate_receipt is not None:
        require_gate_evidence(
            aggregate_receipt,
            "authority_19_of_19_replay",
            sha256_bytes(authority_replay_raw),
            AUTHORITY_REPLAY_PATH.as_posix(),
        )
        assert safety_receipt_record is not None
        require_gate_evidence(
            aggregate_receipt,
            "longlineage_preview_safety",
            safety_receipt_record["sha256"],
            safety_receipt_record["published_path"],
        )
        hygiene_receipt = require_git_gate_receipt(
            repo,
            commit,
            aggregate_receipt,
            "repo_hygiene",
            HANDOFF_ROOT / "evidence/repo_hygiene_summary_receipt.json",
        )
        verify_repo_hygiene_receipt(hygiene_receipt)
        verify_repo_tree_hygiene(repo, commit)
        tracked_large_registry = require_git_gate_receipt(
            repo,
            commit,
            aggregate_receipt,
            "tracked_large_asset_policy",
            HANDOFF_ROOT / "registries/tracked_large_asset_registry.json",
        )
        migration_receipt_path = HANDOFF_ROOT / "evidence/large_asset_migration_receipt.json"
        migration_receipt_raw = git_file_bytes(repo, commit, migration_receipt_path.as_posix())
        require_gate_evidence(
            aggregate_receipt,
            "tracked_large_asset_policy",
            sha256_bytes(migration_receipt_raw),
            migration_receipt_path.as_posix(),
        )
        migration_receipt = json.loads(migration_receipt_raw)
        verify_tracked_large_asset_policy(
            tracked_large_registry,
            migration_receipt,
            commit,
        )
        verify_actual_tracked_large_assets(repo, commit)
        verify_reader_acceptance(repo, commit, aggregate_receipt, reader_schema)
        verify_algorithm_crosswalk(repo, commit)

    requested_git_assets = [normalize_repository_path(path) for path in args.git_asset]
    if len(set(requested_git_assets)) != len(requested_git_assets):
        raise ManifestError("duplicate --git-asset path")
    git_asset_paths = list(
        dict.fromkeys(
            [normalize_repository_path(path) for path in REQUIRED_GIT_ASSETS]
            + (
                [normalize_repository_path(path) for path in APPROVED_ONLY_GIT_ASSETS]
                if aggregate_receipt is not None
                else []
            )
            + requested_git_assets
        )
    )

    assets = [build_git_asset(repo, commit, path) for path in git_asset_paths]
    if (
        args.release_asset or args.gate_evidence_asset or args.gate_log_asset
    ) and aggregate_receipt is None:
        raise ManifestError(
            "--release-asset, --gate-evidence-asset and --gate-log-asset require an approved aggregate receipt"
        )
    governed_release_assets = [
        build_governed_release_asset(repo, spec, artifact_registry)
        for spec in args.release_asset
    ]
    gate_evidence_assets = [
        build_gate_evidence_asset(repo, commit, spec, common_gate_schema)
        for spec in args.gate_evidence_asset
    ]
    gate_log_assets = [build_gate_log_asset(repo, spec) for spec in args.gate_log_asset]
    if aggregate_receipt is not None:
        for asset in governed_release_assets:
            for gate_id in ("release_asset_sha256_roundtrip", "full_history_secret_scan"):
                require_gate_evidence(
                    aggregate_receipt,
                    gate_id,
                    asset["sha256"],
                    asset["published_path"],
                )
    assets.extend(governed_release_assets)
    assets.extend(approval_assets)
    assets.extend(safety_receipt_assets)
    assets.extend(gate_evidence_assets)
    assets.extend(gate_log_assets)
    published_paths = [asset["published_path"] for asset in assets]
    if len(set(published_paths)) != len(published_paths):
        raise ManifestError("duplicate published asset path")
    if aggregate_receipt is not None:
        discovered_git_assets = verify_aggregate_evidence_bytes(
            repo,
            commit,
            aggregate_receipt,
            assets,
        )
        existing_paths = {asset["published_path"] for asset in assets}
        assets.extend(
            asset
            for asset in discovered_git_assets
            if asset["published_path"] not in existing_paths
        )
    assets.sort(key=lambda asset: asset["published_path"])

    manifest = {
        "schema_version": "1.0.0",
        "manifest_type": "research_handoff_release_manifest",
        "delivery_level": "research_handoff_snapshot",
        "clinical_use": False,
        "production_use": False,
        "release_status": args.status,
        "release_ready": release_ready,
        "intended_tag": args.intended_tag,
        "generated_at": generated_at,
        "source": {
            "repository_url": args.repository_url,
            "commit": commit,
            "tree": tree,
            "requested_revision": requested_revision,
            "clean_worktree_observed": True,
        },
        "authority": {
            "manifest_path": authority_path,
            "manifest_sha256": sha256_bytes(authority_raw),
            "denominator_registry_path": denominator_path,
            "denominator_registry_sha256": sha256_bytes(denominator_raw),
            "as_of_date": authority["as_of_date"],
            "scope": scope,
            "claim_ceiling": {
                "confirmed_cellular_subclone": 0,
                "confirmed_linear_ancestry": 0,
                "methylation": "association_only",
                "copy_number_loh": "NOT_INTEGRATED",
                "model_conditional_graph_shape_percent": 88.2579,
                "interpretation": "Model-conditional graph shape; not biological accuracy or prevalence.",
            },
        },
        "project_license": {
            "status": project_license_status,
            "spdx_expression": args.project_license_spdx,
            "authority_observation": (
                "GPL-3.0 license text present; project-level SPDX choice requires "
                "maintainer confirmation."
                if project_license_status == "NEEDS_MAINTAINER_CONFIRMATION"
                else "GPL-3.0 license text present; project-level SPDX choice is "
                "maintainer-confirmed and synchronized."
            ),
        },
        "longlineage": {
            "candidate_commit": LONGLINEAGE_CANDIDATE,
            "delivery_level": "research_preview_non_production",
            "preview_status": preview_status,
            "preview_safety_status": args.longlineage_preview_safety_status,
            "safety_receipt": safety_receipt_record,
            "production_ready": False,
            "phase_gate_interpretation": (
                "Expected research-preview ceiling; not InterSubMod handoff blockers "
                "and never production approval."
            ),
            "gates": {gate: "BLOCKED" for gate in LONGLINEAGE_GATES},
        },
        "blockers": blockers,
        "approval_receipts": approvals,
        "assets": assets,
        "asset_count": len(assets),
        "checksum_sidecar": {
            "algorithm": "SHA-256",
            "filename": args.output_checksums.name,
            "manifest_filename": args.output_json.name,
            "line_count": len(assets) + 1,
            "coverage": (
                "release manifest plus declared assets; "
                "checksum sidecar self-hash intentionally excluded"
            ),
        },
    }
    validate_document(manifest, schema, "generated manifest")
    manifest_bytes = (
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    checksum_lines = [f"{sha256_bytes(manifest_bytes)}  {args.output_json.name}"]
    checksum_lines.extend(f"{asset['sha256']}  {asset['published_path']}" for asset in assets)
    return manifest, checksum_lines


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--revision", default="HEAD")
    parser.add_argument("--repository-url", default="https://github.com/liaoyoyo/InterSubMod")
    parser.add_argument("--intended-tag", default="research-handoff-2026.08.1")
    parser.add_argument("--generated-at")
    parser.add_argument("--status", choices=("BLOCKED", "APPROVED_RESEARCH_HANDOFF"), required=True)
    parser.add_argument("--blocker", action="append", default=[])
    parser.add_argument(
        "--project-license-spdx",
        choices=("GPL-3.0-only", "GPL-3.0-or-later"),
        help="Maintainer-confirmed project-level SPDX choice; omit while unresolved",
    )
    parser.add_argument(
        "--approval-receipt",
        action="append",
        default=[],
        metavar="[PUBLISHED_NAME=]PATH",
        help="External aggregate acceptance receipt bound to the source commit",
    )
    parser.add_argument(
        "--longlineage-preview-safety-status",
        choices=("BLOCKED", "VERIFIED"),
        required=True,
    )
    parser.add_argument(
        "--longlineage-safety-receipt",
        metavar="[PUBLISHED_NAME=]PATH",
        help="External license/source/dependency/public-safety PASS receipt",
    )
    parser.add_argument(
        "--longlineage-candidate",
        required=True,
        help=f"Pinned preview commit; must equal {LONGLINEAGE_CANDIDATE}",
    )
    parser.add_argument(
        "--longlineage-repo-root",
        type=Path,
        help=(
            "Read-only LongLineage Git top level used to verify exact safety-head/tree "
            "and support-asset bytes; required when preview safety is VERIFIED"
        ),
    )
    parser.add_argument("--git-asset", action="append", default=[])
    parser.add_argument(
        "--release-asset",
        action="append",
        default=[],
        metavar="ARTIFACT_ID::[PUBLISHED_NAME=]PATH",
        help=(
            "Governed external file: registry record must be FINAL_FOR_SCOPE, publishable, "
            "licensed, hash/size matched, and assigned to GITHUB_RELEASE"
        ),
    )
    parser.add_argument(
        "--gate-evidence-asset",
        action="append",
        default=[],
        metavar="[PUBLISHED_NAME=]PATH",
        help=(
            "External JSON receipt bound to the source commit and referenced by exact "
            "sha256/URI from one or more aggregate gates"
        ),
    )
    parser.add_argument(
        "--gate-log-asset",
        action="append",
        default=[],
        metavar="[PUBLISHED_NAME=]PATH",
        help=(
            "Raw execution log/report referenced by exact URI and SHA-256 from a "
            "gate receipt and the aggregate gate"
        ),
    )
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-checksums", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.longlineage_candidate != LONGLINEAGE_CANDIDATE:
            raise ManifestError(
                f"--longlineage-candidate must equal pinned commit {LONGLINEAGE_CANDIDATE}"
            )
        output_json = ensure_output_outside_repository(args.repo_root.expanduser().resolve(), args.output_json)
        output_checksums = ensure_output_outside_repository(
            args.repo_root.expanduser().resolve(), args.output_checksums
        )
        if output_json == output_checksums:
            raise ManifestError("manifest and checksum outputs must be different paths")
        print(f"[INPUT] repository={args.repo_root.expanduser().resolve()}")
        print(f"[INPUT] revision={args.revision} status={args.status}")
        manifest, checksum_lines = build_manifest(args)
        write_text_exclusive(
            output_json,
            json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        )
        write_text_exclusive(output_checksums, "\n".join(checksum_lines) + "\n")
        print(f"[OUTPUT] release_manifest={output_json}")
        print(f"[OUTPUT] checksums={output_checksums}")
        print(
            f"[RESULT] release_ready={str(manifest['release_ready']).lower()} "
            f"assets={manifest['asset_count']} blockers={len(manifest['blockers'])} "
            f"source_commit={manifest['source']['commit']}"
        )
        if manifest["assets"]:
            sample = manifest["assets"][0]
            print(
                f"[SAMPLE] {sample['publication_channel']} {sample['sha256']} "
                f"{sample['published_path']}"
            )
        return 0
    except (ManifestError, json.JSONDecodeError, OSError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
