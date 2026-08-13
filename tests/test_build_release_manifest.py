from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Callable

from jsonschema import FormatChecker
from jsonschema.validators import validator_for


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path(os.environ.get("INTERSUBMOD_TEST_SOURCE_ROOT", REPO_ROOT)).resolve()
SCRIPT = REPO_ROOT / "scripts/handoff/build_release_manifest.py"
SCHEMA = REPO_ROOT / (
    "docs/handoff/20260813_完整研究資料與軟體交接_01/"
    "schemas/research-handoff-release-manifest.schema.json"
)
HANDOFF = Path("docs/handoff/20260813_完整研究資料與軟體交接_01")
LONG_LINEAGE_COMMIT = "b9aaa12a11fa00606bd174dabd0f172a5d112359"
LONG_LINEAGE_HISTORICAL_STACK = "f60b5f3274123bdf818371600b608d52626c893e"
LONG_LINEAGE_SAFETY_STACK = "e" * 40
AGGREGATE_SCHEMA = REPO_ROOT / (HANDOFF / "schemas/aggregate-release-acceptance.schema.json")
LONGLINEAGE_SCHEMA = REPO_ROOT / (
    HANDOFF / "schemas/longlineage-preview-safety-receipt.schema.json"
)
COMMON_GATE_SCHEMA = REPO_ROOT / (HANDOFF / "schemas/gate-acceptance-receipt.schema.json")
READER_SCHEMA = REPO_ROOT / (HANDOFF / "schemas/reader-acceptance.schema.json")
AUTHORITY_REPLAY = HANDOFF / "evidence/authority_replay_receipt.json"
REQUIRED_GATES = (
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
EVIDENCE_TYPES = {
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
VERIFIED_AT = "2026-08-13T19:00:00+08:00"
PROJECT_LICENSE = "GPL-3.0-or-later"
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
    HANDOFF / "00_INDEX.md",
    HANDOFF / "20260813_研究結論時間與Finality_01.md",
    HANDOFF / "20260813_軟體輸入輸出與研究流程_01.md",
    HANDOFF / "20260813_bip7_bip8操作與驗證_01.md",
    HANDOFF / "20260813_完整研究交接總覽_01.standalone.html",
    HANDOFF / "ai_context/CONTEXT.md",
    HANDOFF / "ai_context/READER_ACCEPTANCE_PROMPT.md",
    HANDOFF / "evidence/EVIDENCE_MANIFEST.json",
    HANDOFF / "evidence/authority_manifest.json",
    HANDOFF / "evidence/denominator_registry.tsv",
    AUTHORITY_REPLAY,
    HANDOFF / "evidence/algorithm_cli_claim_crosswalk.json",
    HANDOFF / "evidence/full_claim_registry.json",
    HANDOFF / "evidence/large_asset_migration_receipt.json",
    HANDOFF / "evidence/longlineage_SBOM.spdx.json",
    HANDOFF / "evidence/longlineage_clean_foundation_receipt.json",
    HANDOFF / "evidence/longlineage_capability_matrix.md",
    HANDOFF / "evidence/longlineage_public_safety_receipt.json",
    HANDOFF / "evidence/longlineage_third_party_notices.md",
    HANDOFF / "evidence/public_claim_validation_receipt.json",
    HANDOFF / "evidence/repo_hygiene_summary_receipt.json",
    HANDOFF / "registries/artifact_registry.json",
    HANDOFF / "registries/authority_superseded_crosswalk.json",
    HANDOFF / "registries/claim_registry.json",
    HANDOFF / "registries/dataset_alias_registry.json",
    HANDOFF / "registries/dataset_registry.json",
    HANDOFF / "registries/machine_path_registry.json",
    HANDOFF / "registries/registry_build_receipt.json",
    HANDOFF / "registries/run_registry.json",
    HANDOFF / "registries/storage_root_manifest.json",
    HANDOFF / "registries/tracked_large_asset_registry.json",
    HANDOFF / "registries/workflow_registry.json",
    HANDOFF / "schemas/research-handoff-release-manifest.schema.json",
    HANDOFF / "schemas/aggregate-release-acceptance.schema.json",
    HANDOFF / "schemas/longlineage-preview-safety-receipt.schema.json",
    HANDOFF / "schemas/gate-acceptance-receipt.schema.json",
    HANDOFF / "schemas/reader-acceptance.schema.json",
    HANDOFF / "schemas/algorithm-cli-claim-crosswalk.schema.json",
    HANDOFF / "schemas/artifact-registry.schema.json",
    HANDOFF / "schemas/authority-crosswalk.schema.json",
    HANDOFF / "schemas/dataset-registry.schema.json",
    HANDOFF / "schemas/machine-path-registry.schema.json",
    HANDOFF / "schemas/run-registry.schema.json",
    HANDOFF / "schemas/storage-root-manifest.schema.json",
    HANDOFF / "schemas/tracked-large-asset-registry.schema.json",
    Path("scripts/site/bootstrap"),
    Path("scripts/site/doctor"),
    Path("scripts/site/site_profile.py"),
    Path("scripts/handoff/build_algorithm_cli_crosswalk.py"),
    Path("scripts/handoff/build_large_asset_registry.py"),
    Path("scripts/handoff/build_registries.py"),
    Path("scripts/handoff/build_tiny_public_fixture.sh"),
    Path("scripts/handoff/build_workflow_registry.py"),
    Path("scripts/handoff/build_release_manifest.py"),
    Path("scripts/handoff/validate_reader_acceptance.py"),
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
    Path("tests/test_validate_reader_acceptance.py"),
)
SYNTHETIC_REQUIRED_ASSETS = {
    HANDOFF / "evidence/reader_acceptance_receipt.json",
    HANDOFF / "evidence/longlineage_source_to_target_manifest.json",
    HANDOFF / "evidence/longlineage_public_gate_receipt.json",
    HANDOFF / "evidence/longlineage_history_scan_receipt.json",
    HANDOFF / "evidence/longlineage_check_public_preview_gate.sh",
    HANDOFF / "evidence/longlineage_public_gate_output.txt",
    HANDOFF / "evidence/longlineage_history_scan_output.json",
}
APPROVED_ONLY_GIT_ASSETS = tuple(
    sorted(SYNTHETIC_REQUIRED_ASSETS, key=lambda path: path.as_posix())
)
DYNAMIC_GATE_GIT_EVIDENCE: set[str] = set()
SPECIAL_GATE_IDS = {
    "authority_19_of_19_replay",
    "repo_hygiene",
    "tracked_large_asset_policy",
    "longlineage_preview_safety",
}
COMMON_GATE_IDS = set(REQUIRED_GATES) - SPECIAL_GATE_IDS


def run(command: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)


def initialize_fixture_repository(tmp_path: Path) -> tuple[Path, str]:
    repo = tmp_path / "source"
    archive_preferred = {
        Path("CITATION.cff"),
        Path("CONTRIBUTING.md"),
        Path("SECURITY.md"),
        Path("CHANGELOG.md"),
        HANDOFF / "schemas/research-handoff-release-manifest.schema.json",
        HANDOFF / "schemas/aggregate-release-acceptance.schema.json",
        HANDOFF / "schemas/longlineage-preview-safety-receipt.schema.json",
        Path("scripts/handoff/build_release_manifest.py"),
    }
    for destination in REQUIRED_GIT_ASSETS:
        if destination == Path("CITATION.cff") or destination in SYNTHETIC_REQUIRED_ASSETS:
            continue
        candidates = (
            (REPO_ROOT / destination, SOURCE_ROOT / destination)
            if destination in archive_preferred
            else (SOURCE_ROOT / destination, REPO_ROOT / destination)
        )
        source = next((candidate for candidate in candidates if candidate.is_file()), None)
        assert source is not None, destination
        target = repo / destination
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
    (repo / "README.md").write_text("fixture\n", encoding="utf-8")
    (repo / "CITATION.cff").write_text(
        "cff-version: 1.2.0\ntitle: Fixture\ntype: software\n"
        "authors:\n  - name: Fixture\nlicense: GPL-3.0-or-later\n",
        encoding="utf-8",
    )
    hygiene_path = repo / HANDOFF / "evidence/repo_hygiene_summary_receipt.json"
    hygiene_path.write_text(
        json.dumps(
            {
                "receipt_type": "repo_hygiene_package_summary",
                "schema_version": "1.0.0",
                "pass": True,
                "versioned_worktree_absolute_symlink_count": 0,
                "versioned_worktree_broken_symlink_count": 0,
                "tracked_settings_local_present": False,
                "archive": {"pass": True},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    tracked_path = repo / HANDOFF / "registries/tracked_large_asset_registry.json"
    tracked_registry = json.loads(tracked_path.read_text(encoding="utf-8"))
    tracked_registry["source_git_head"] = "0" * 40
    tracked_registry["summary"] = {
        "tracked_assets_over_threshold": 0,
        "total_size_bytes": 0,
        "tracked_assets_over_50_mib": 0,
        "verdict": "PASS",
    }
    tracked_registry["assets"] = []
    tracked_path.write_text(json.dumps(tracked_registry) + "\n", encoding="utf-8")
    migration_path = repo / HANDOFF / "evidence/large_asset_migration_receipt.json"
    migration = json.loads(migration_path.read_text(encoding="utf-8"))
    migration["source_git_commit"] = "0" * 40
    migration["verdict"] = "GITHUB_RELEASE_ROUNDTRIP_VERIFIED"
    migration["remaining_policy_status"] = "PASS"
    migration_path.write_text(json.dumps(migration) + "\n", encoding="utf-8")
    artifact_path = repo / HANDOFF / "registries/artifact_registry.json"
    artifacts = json.loads(artifact_path.read_text(encoding="utf-8"))
    for artifact in artifacts:
        value = str(artifact.get("license", ""))
        value = value.replace(
            "PROJECT_LICENSE_PENDING_MAINTAINER_CONFIRMATION", PROJECT_LICENSE
        )
        value = value.replace("GPL-3.0-only", PROJECT_LICENSE)
        artifact["license"] = value
    artifact_path.write_text(json.dumps(artifacts) + "\n", encoding="utf-8")
    crosswalk_path = repo / HANDOFF / "evidence/algorithm_cli_claim_crosswalk.json"
    crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
    for gate_id in (
        "CROSSWALK_COMPLETENESS",
        "PUBLICATION_ASSET_READY",
        "RELEASE_ASSET_READY",
    ):
        crosswalk["gates"][gate_id]["status"] = "PASS"
        crosswalk["gates"][gate_id]["reason"] = "Fixture gate passed with direct evidence."
    crosswalk_path.write_text(json.dumps(crosswalk) + "\n", encoding="utf-8")
    mapping_path = repo / HANDOFF / "evidence/longlineage_source_to_target_manifest.json"
    mappings = [
        {
            "origin_id": f"MAP-{index:02d}",
            "source_commit": f"{index + 1:040x}",
            "source_sha256": hashlib.sha256(f"source-{index}".encode()).hexdigest(),
            "source_replay_status": "MATCHED_DECLARED_ORIGIN_COMMIT",
            "target": f"src/mapped_{index}.cpp",
            "target_sha256": hashlib.sha256(f"target-{index}".encode()).hexdigest(),
            "license_disposition": "APPROVED_FOR_PUBLIC_RELEASE",
            "license_evidence": "OWNER_APPROVED_SOURCE_ORIGIN_AUDIT",
        }
        for index in range(21)
    ]
    mapping_path.write_text(
        json.dumps(
            {
                "schema_version": "1.2.0",
                "license_review_status": "APPROVED_FOR_PUBLIC_RELEASE",
                "mappings": mappings,
            }
        )
        + "\n",
        encoding="utf-8",
    )
    gate_script_path = repo / HANDOFF / "evidence/longlineage_check_public_preview_gate.sh"
    gate_script_path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    gate_output_path = repo / HANDOFF / "evidence/longlineage_public_gate_output.txt"
    gate_output_path.write_text("PUBLIC_PREVIEW_GATE_PASS\n", encoding="utf-8")
    history_output_path = repo / HANDOFF / "evidence/longlineage_history_scan_output.json"
    history_output_path.write_text('{"findings":[]}\n', encoding="utf-8")
    (repo / HANDOFF / "evidence/longlineage_SBOM.spdx.json").write_text(
        json.dumps(
            {
                "spdxVersion": "SPDX-2.3",
                "packages": [
                    {
                        "name": f"longlineage-{row['origin_id']}",
                        "SPDXID": f"SPDXRef-Origin-{row['origin_id']}",
                        "licenseConcluded": "GPL-3.0-only",
                    }
                    for row in mappings
                ],
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (repo / HANDOFF / "evidence/longlineage_third_party_notices.md").write_text(
        "# LongLineage third-party notices\n\nStatus: APPROVED_FOR_PUBLIC_RELEASE\n",
        encoding="utf-8",
    )
    gate_receipt_path = repo / HANDOFF / "evidence/longlineage_public_gate_receipt.json"
    gate_receipt_path.write_text(
        json.dumps(
            {
                "candidate_commit": LONG_LINEAGE_COMMIT,
                "safety_stacked_commit": LONG_LINEAGE_SAFETY_STACK,
                "safety_tree": "0" * 40,
                "repository_url": "https://github.com/liaoyoyo/LongLineage",
                "cross_repo_verification_method": "GIT_SHOW_EXACT_HEAD",
                "remote_branch": "agent/approved-safety-fixture",
                "remote_advertised_head": LONG_LINEAGE_SAFETY_STACK,
                "remote_advertised_observed_at": VERIFIED_AT,
                "support_asset_bindings": [],
                "source_manifest_sha256": file_sha256(mapping_path),
                "script_sha256": file_sha256(gate_script_path),
                "output_sha256": file_sha256(gate_output_path),
                "exit_code": 0,
                "ctest_total": 49,
                "ctest_passed": 49,
                "ctest_failed": 0,
                "strict_public_preview_gate_exit": 0,
                "production_run_exit": 6,
                "production_run_status": "KernelBlocked",
                "phase_gates": {
                    gate_id: "BLOCKED" for gate_id in ("P3", "P4", "P5", "P7", "P8")
                },
                "overall_pass": True,
                "blocking_gates": [],
            }
        )
        + "\n",
        encoding="utf-8",
    )
    history_receipt_path = repo / HANDOFF / "evidence/longlineage_history_scan_receipt.json"
    history_receipt_path.write_text(
        json.dumps(
            {
                "base_commit": LONG_LINEAGE_COMMIT,
                "head_commit": LONG_LINEAGE_SAFETY_STACK,
                "scope": "FULL_HISTORY",
                "scanner": "gitleaks",
                "scanner_version": "fixture-1.0",
                "command": "gitleaks git --log-opts --all",
                "raw_output_sha256": file_sha256(history_output_path),
                "finding_count": 0,
                "findings": [],
                "overall_pass": True,
            }
        )
        + "\n",
        encoding="utf-8",
    )
    stacked_commit = LONG_LINEAGE_SAFETY_STACK
    ll_safety_path = repo / HANDOFF / "evidence/longlineage_public_safety_receipt.json"
    ll_safety_path.write_text(
        json.dumps(
            {
                "candidate_commit": LONG_LINEAGE_COMMIT,
                "safety_stacked_commit": stacked_commit,
                "gate_status": "PUBLIC_PREVIEW_SAFETY_PASS",
                "required_verdict": "PUBLIC_PREVIEW_SAFETY_PASS",
                "history_hygiene": {"status": "PASS", "finding_count": 0},
                "source_replay": {
                    "mapping_rows": 21,
                    "unresolved_rows": 0,
                    "target_digest_matches": 21,
                    "target_digest_mismatches": 0,
                    "missing_unique_hashes": [],
                },
                "license_review": {
                    "repository_status": "APPROVED",
                    "row_dispositions_pending": 0,
                    "row_dispositions_approved": 21,
                    "third_party_notice_status": "APPROVED",
                    "sbom_status": "APPROVED",
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    ll_foundation_path = repo / HANDOFF / "evidence/longlineage_clean_foundation_receipt.json"
    ll_foundation_path.write_text(
        json.dumps(
            {
                "candidate_baseline": LONG_LINEAGE_COMMIT,
                "validated_stacked_commit": LONG_LINEAGE_HISTORICAL_STACK,
                "results": {
                    "configure_exit": 0,
                    "build_exit": 0,
                    "ctest_failed": 0,
                    "ctest_total": 49,
                    "ctest_passed": 49,
                    "strict_public_preview_gate_exit": 1,
                    "strict_public_preview_gate_interpretation": "EXPECTED_SAFE_FAILURE",
                    "production_run_expected_exit": 6,
                    "production_run_expected_status": "KernelBlocked",
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (repo / HANDOFF / "evidence/longlineage_SBOM.spdx.json").write_text(
        json.dumps(
            {
                "spdxVersion": "SPDX-2.3",
                "packages": [
                    {
                        "name": f"longlineage-{row['origin_id']}",
                        "SPDXID": f"SPDXRef-Origin-{row['origin_id']}",
                        "licenseConcluded": "GPL-3.0-only",
                    }
                    for row in mappings
                ],
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (repo / HANDOFF / "evidence/longlineage_third_party_notices.md").write_text(
        "# LongLineage third-party notices\n\nStatus: APPROVED_FOR_PUBLIC_RELEASE\n",
        encoding="utf-8",
    )
    finalize_longlineage_fixture(repo)
    assert run(["git", "init", "-q"], repo).returncode == 0
    assert run(["git", "config", "user.name", "Release manifest test"], repo).returncode == 0
    assert run(["git", "config", "user.email", "test@example.invalid"], repo).returncode == 0
    assert run(["git", "add", "."], repo).returncode == 0
    assert run(["git", "commit", "-qm", "fixture"], repo).returncode == 0
    tested_commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    write_reader_receipt(repo, tested_commit)
    assert run(["git", "add", (HANDOFF / "evidence/reader_acceptance_receipt.json").as_posix()], repo).returncode == 0
    assert run(["git", "commit", "-qm", "add reader acceptance receipt"], repo).returncode == 0
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    return repo, commit


def finalize_longlineage_fixture(repo: Path) -> tuple[Path, str]:
    ll_repo = repo.parent / "longlineage"
    path_pairs = (
        (HANDOFF / "evidence/longlineage_source_to_target_manifest.json", Path("provenance/source_to_target_manifest.json")),
        (HANDOFF / "evidence/longlineage_check_public_preview_gate.sh", Path("scripts/check_public_preview_gate.sh")),
        (HANDOFF / "evidence/longlineage_public_gate_output.txt", Path("receipts/public_gate_output.txt")),
        (HANDOFF / "evidence/longlineage_history_scan_output.json", Path("receipts/history_scan_output.json")),
        (HANDOFF / "evidence/longlineage_SBOM.spdx.json", Path("SBOM.spdx.json")),
        (HANDOFF / "evidence/longlineage_third_party_notices.md", Path("THIRD_PARTY_NOTICES.md")),
    )
    for ism_relative, ll_relative in path_pairs:
        target = ll_repo / ll_relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((repo / ism_relative).read_bytes())
    assert run(["git", "init", "-q"], ll_repo).returncode == 0
    assert run(["git", "config", "user.name", "LongLineage fixture"], ll_repo).returncode == 0
    assert run(["git", "config", "user.email", "ll@example.invalid"], ll_repo).returncode == 0
    assert run(["git", "remote", "add", "origin", "https://github.com/liaoyoyo/LongLineage.git"], ll_repo).returncode == 0
    assert run(["git", "add", "."], ll_repo).returncode == 0
    assert run(["git", "commit", "-qm", "approved safety fixture"], ll_repo).returncode == 0
    ll_commit = run(["git", "rev-parse", "HEAD"], ll_repo).stdout.strip()
    ll_tree = run(["git", "rev-parse", "HEAD^{tree}"], ll_repo).stdout.strip()
    bare_remote = repo.parent / "longlineage-origin.git"
    assert run(["git", "init", "--bare", "-q", str(bare_remote)], repo.parent).returncode == 0
    rewrite_key = "url.file://" + str(bare_remote) + ".insteadOf"
    assert run(
        ["git", "config", rewrite_key, "https://github.com/liaoyoyo/LongLineage.git"],
        ll_repo,
    ).returncode == 0
    assert run(
        ["git", "push", "-q", "origin", "HEAD:refs/heads/agent/approved-safety-fixture"],
        ll_repo,
    ).returncode == 0

    gate_path = repo / HANDOFF / "evidence/longlineage_public_gate_receipt.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    gate["safety_stacked_commit"] = ll_commit
    gate["safety_tree"] = ll_tree
    gate["remote_branch"] = "agent/approved-safety-fixture"
    gate["remote_advertised_head"] = ll_commit
    gate["remote_advertised_observed_at"] = VERIFIED_AT
    bindings = []
    for ism_relative, ll_relative in path_pairs:
        digest = file_sha256(repo / ism_relative)
        bindings.append(
            {
                "repository_url": "https://github.com/liaoyoyo/LongLineage",
                "head_commit": ll_commit,
                "longlineage_path": ll_relative.as_posix(),
                "ism_path": ism_relative.as_posix(),
                "sha256": digest,
            }
        )
    gate["support_asset_bindings"] = bindings
    gate_path.write_text(json.dumps(gate) + "\n", encoding="utf-8")

    history_path = repo / HANDOFF / "evidence/longlineage_history_scan_receipt.json"
    history = json.loads(history_path.read_text(encoding="utf-8"))
    history["head_commit"] = ll_commit
    history_path.write_text(json.dumps(history) + "\n", encoding="utf-8")
    safety_path = repo / HANDOFF / "evidence/longlineage_public_safety_receipt.json"
    safety = json.loads(safety_path.read_text(encoding="utf-8"))
    safety["safety_stacked_commit"] = ll_commit
    safety_path.write_text(json.dumps(safety) + "\n", encoding="utf-8")
    return ll_repo, ll_commit


def write_reader_receipt(repo: Path, tested_commit: str) -> None:
    core_paths = (
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
    )
    question_ids = (
        "Q_PROJECT",
        "Q_CONCLUSION",
        "Q_FINALITY",
        "Q_SOFTWARE_ROLES",
        "Q_MACHINES",
        "Q_VERIFY_CONTINUE",
    )
    promotion_ids = (
        "NO_CELLULAR_PROMOTION",
        "NO_ANCESTRY_PROMOTION",
        "NO_882579_ACCURACY_OR_PREVALENCE",
        "NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE",
        "NO_FEATURE_AS_MAIN",
        "NO_LOCAL_AS_LIVE_PUBLISHED",
        "NO_BIP7_AS_BIP8",
        "NO_PYTHON_ROLE_CONFLATION",
    )
    answers = {
        "Q_PROJECT": "This research handoff snapshot covers ONT long-read sSNV read linkage, candidate mutation-state reconstruction, and methylation association. It is not a production release or cellular-lineage caller.",
        "Q_CONCLUSION": "Confirmed cellular subclones = 0 and confirmed linear ancestry relationships = 0. Methylation is association only. In the frozen reconstruction CN/LOH is not integrated. 88.2579% is a model-conditional rooted graph-shape statistic; it is not accuracy and it is not prevalence.",
        "Q_FINALITY": "An artifact can be FINAL_FOR_SCOPE only after its evidence_status, finality, producer, scope and hash are checked, while supersedes identifies replaced material. Final is bounded only to the declared scope; it is not production-ready or whole research final.",
        "Q_SOFTWARE_ROLES": "LongPhase-S/TO produces HP/PS, phased or recalibrated VCF and tagged BAM in one chain. In a second parallel provenance chain, exact-PS/LongLineage produces candidate families and read assignments. InterSubMod produces per-region methylation, distance, read clustering and statistics. A commit-pinned Python research solver can be a science producer; a validator/publication builder or HTML presenter validates and presents validated data and does not recompute science.",
        "Q_MACHINES": "On bip7 a bounded local preflight exists, but fresh-clone validation is still pending and BLOCKED. bip8 is BLOCKED because there is no host-local receipt. Each host-specific receipt cannot substitute for or prove the other host. Run doctor, then build/test, then synthetic smoke on each machine.",
        "Q_VERIFY_CONTINUE": "Replay authority hashes: 19/19 byte SHA-256 records MATCH; this 19/19 replay does not constitute a science rerun. Validate registries for schema and unique IDs, then respect publication/release/machine gates that remain BLOCKED. For a new research cycle, run pre-decision audit and pin commit, input hash, scope and denominator.",
    }
    promotion_explanations = {
        "NO_CELLULAR_PROMOTION": "Confirmed cellular subclones = 0; candidate read groups are not cellular truth and cannot be promoted to validated cellular subclones.",
        "NO_ANCESTRY_PROMOTION": "Confirmed linear ancestry = 0; a graph shape or read dendrogram is not ancestry and cannot establish a cellular phylogeny.",
        "NO_882579_ACCURACY_OR_PREVALENCE": "88.2579% is a model-conditional rooted graph-shape statistic; it is not accuracy and it is not prevalence or biological truth.",
        "NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE": "The frozen scope contains 7 technical datasets and 6 biological IDs; they are not 7 independent biological replicates or samples.",
        "NO_FEATURE_AS_MAIN": "b9aaa12 is a research preview and non-production candidate; P3/P4/P5/P7/P8 remain BLOCKED, and a feature is not main, release, or production capability.",
        "NO_LOCAL_AS_LIVE_PUBLISHED": "Local source validation is not live publication; main, Wiki and Pages remain BLOCKED until publish and re-fetch verification is completed.",
        "NO_BIP7_AS_BIP8": "A bip7 receipt cannot prove or substitute for bip8; bip8 is BLOCKED until a host-local receipt is created on that hostname.",
        "NO_PYTHON_ROLE_CONFLATION": "A commit-pinned Python research solver can be a science producer; a validator/publication builder or HTML presenter validates and presents, and the presenter does not recompute science.",
    }
    question_evidence = {
        "Q_PROJECT": ["00_INDEX.md", "20260813_軟體輸入輸出與研究流程_01.md"],
        "Q_CONCLUSION": ["20260813_研究結論時間與Finality_01.md", "evidence/denominator_registry.tsv"],
        "Q_FINALITY": ["20260813_研究結論時間與Finality_01.md", "registries/artifact_registry.json", "registries/authority_superseded_crosswalk.json"],
        "Q_SOFTWARE_ROLES": ["20260813_軟體輸入輸出與研究流程_01.md", "evidence/longlineage_capability_matrix.md"],
        "Q_MACHINES": ["20260813_bip7_bip8操作與驗證_01.md", "registries/machine_path_registry.json"],
        "Q_VERIFY_CONTINUE": ["00_INDEX.md", "evidence/authority_replay_receipt.json", "registries/claim_registry.json"],
    }
    promotion_evidence = {
        "NO_CELLULAR_PROMOTION": ["20260813_研究結論時間與Finality_01.md"],
        "NO_ANCESTRY_PROMOTION": ["20260813_研究結論時間與Finality_01.md"],
        "NO_882579_ACCURACY_OR_PREVALENCE": ["20260813_研究結論時間與Finality_01.md", "evidence/denominator_registry.tsv"],
        "NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE": ["00_INDEX.md", "evidence/authority_manifest.json"],
        "NO_FEATURE_AS_MAIN": ["evidence/longlineage_capability_matrix.md"],
        "NO_LOCAL_AS_LIVE_PUBLISHED": ["00_INDEX.md", "registries/claim_registry.json"],
        "NO_BIP7_AS_BIP8": ["20260813_bip7_bip8操作與驗證_01.md", "registries/machine_path_registry.json"],
        "NO_PYTHON_ROLE_CONFLATION": ["20260813_軟體輸入輸出與研究流程_01.md"],
    }
    receipt = {
        "schema_name": "intersubmod.reader_acceptance_receipt",
        "schema_version": "2.0.0",
        "receipt_id": "intersubmod-reader-acceptance-20260813",
        "evaluated_at": VERIFIED_AT,
        "tested_git_commit": tested_commit,
        "evaluator_context": "NO_PRIOR_CONVERSATION_PACKAGE_ONLY",
        "evaluator_id": "fixture-independent-reader",
        "package_source_manifest": [
            {
                "path": path,
                "sha256": hashlib.sha256(
                    subprocess.run(
                        ["git", "show", f"{tested_commit}:{(HANDOFF / path).as_posix()}"],
                        cwd=repo,
                        check=True,
                        stdout=subprocess.PIPE,
                    ).stdout
                ).hexdigest(),
            }
            for path in core_paths
        ],
        "questions": [
            {
                "question_id": question_id,
                "status": "PASS",
                "answer": answers[question_id],
                "evidence_paths": question_evidence[question_id],
                "limits": ["Reader comprehension is not scientific validation or publication proof."],
            }
            for question_id in question_ids
        ],
        "prohibited_promotions": [
            {
                "check_id": check_id,
                "status": "PASS",
                "explanation": promotion_explanations[check_id],
                "evidence_paths": promotion_evidence[check_id],
            }
            for check_id in promotion_ids
        ],
        "verdict": "PASS",
        "claim_ceiling": "Reader comprehension only; not scientific validation, host acceptance, or publication proof.",
    }
    path = repo / HANDOFF / "evidence/reader_acceptance_receipt.json"
    path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")


def blocked_command(repo: Path, output_root: Path) -> list[str]:
    return [
        sys.executable,
        str(SCRIPT),
        "--repo-root",
        str(repo),
        "--status",
        "BLOCKED",
        "--blocker",
        "BIP8_DATA_PREFLIGHT_BLOCKED",
        "--blocker",
        "INTERSUBMOD_PROJECT_LICENSE_CONFIRMATION_REQUIRED",
        "--longlineage-preview-safety-status",
        "BLOCKED",
        "--longlineage-candidate",
        LONG_LINEAGE_COMMIT,
        "--git-asset",
        "README.md",
        "--generated-at",
        "2026-08-13T19:00:00+08:00",
        "--output-json",
        str(output_root / "release-manifest.json"),
        "--output-checksums",
        str(output_root / "SHA256SUMS"),
    ]


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_longlineage_safety_receipt(
    path: Path, source_commit: str, repo: Path
) -> None:
    support_paths = {
        "license": (
            HANDOFF / "evidence/longlineage_third_party_notices.md",
            HANDOFF / "evidence/longlineage_SBOM.spdx.json",
        ),
        "source_origin": (
            HANDOFF / "evidence/longlineage_public_safety_receipt.json",
            HANDOFF / "evidence/longlineage_source_to_target_manifest.json",
        ),
        "dependencies": (HANDOFF / "evidence/longlineage_SBOM.spdx.json",),
        "public_safety": (
            HANDOFF / "evidence/longlineage_public_safety_receipt.json",
            HANDOFF / "evidence/longlineage_clean_foundation_receipt.json",
            HANDOFF / "evidence/longlineage_public_gate_receipt.json",
            HANDOFF / "evidence/longlineage_history_scan_receipt.json",
            HANDOFF / "evidence/longlineage_check_public_preview_gate.sh",
            HANDOFF / "evidence/longlineage_public_gate_output.txt",
            HANDOFF / "evidence/longlineage_history_scan_output.json",
        ),
    }

    def gate(gate_id: str) -> dict[str, object]:
        return {
            "status": "PASS",
            "subject_commit": LONG_LINEAGE_COMMIT,
            "verified_at": VERIFIED_AT,
            "evidence": [
                {
                    "artifact_id": f"longlineage/{support_path.name}",
                    "subject_commit": LONG_LINEAGE_COMMIT,
                    "sha256": file_sha256(repo / support_path),
                    "uri": support_path.as_posix(),
                    "verified_at": VERIFIED_AT,
                }
                for support_path in support_paths[gate_id]
            ],
        }

    gate_receipt = json.loads(
        (repo / HANDOFF / "evidence/longlineage_public_gate_receipt.json").read_text(
            encoding="utf-8"
        )
    )
    receipt = {
        "schema_version": "1.0.0",
        "receipt_type": "longlineage_public_preview_safety",
        "source_commit": source_commit,
        "longlineage_candidate_commit": LONG_LINEAGE_COMMIT,
        "safety_stacked_commit": gate_receipt["safety_stacked_commit"],
        "generated_at": VERIFIED_AT,
        "overall_pass": True,
        "blocking_gates": [],
        "preview_safety_gates": {
            gate_id: gate(gate_id)
            for gate_id in ("license", "source_origin", "dependencies", "public_safety")
        },
        "phase_gates": {
            gate_id: "BLOCKED" for gate_id in ("P3", "P4", "P5", "P7", "P8")
        },
    }
    path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")


def write_common_gate_receipts(
    receipt_root: Path, source_commit: str, repo: Path
) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for gate_id in sorted(COMMON_GATE_IDS):
        kind = {
            "github_publication_commit_alignment": "GITHUB_API",
            "github_branch_protection_required_ci": "GITHUB_API",
            "reader_ai_six_question_acceptance": "READER_AGENT",
            "intersubmod_project_license_confirmation": "GITHUB_API",
        }.get(gate_id, "COMMAND")
        details = gate_details(gate_id, source_commit, repo)
        execution_log = f"fixture execution log for {gate_id}\n"
        log_path = receipt_root / f"gate-{gate_id}.log"
        log_path.write_text(execution_log, encoding="utf-8")
        execution_command = f"fixture-verifier --gate {gate_id}"
        if gate_id == "full_history_secret_scan":
            execution_command += " --all"
        receipt: dict[str, object] = {
            "schema_version": "1.0.0",
            "receipt_type": f"{gate_id}_acceptance",
            "gate_id": gate_id,
            "source_commit": source_commit,
            "generated_at": VERIFIED_AT,
            "scope": "FULL",
            "execution": {
                "kind": kind,
                "tool": "fixture-verifier",
                "version": "1.0.0",
                "command": execution_command,
                "log_uri": f"release-assets/{log_path.name}",
                "log_sha256": hashlib.sha256(execution_log.encode()).hexdigest(),
                "exit_code": 0,
            },
            "details": details,
            "result": {"status": "PASS", "checks": 1, "failures": 0},
            "overall_pass": True,
            "blocking_gates": [],
        }
        if gate_id in ("bip7_fresh_clone_acceptance", "bip8_fresh_clone_acceptance"):
            receipt["host"] = {
                "hostname": gate_id.split("_", 1)[0],
                "mount_type": "fixturefs",
                "tool_hashes": {"inter_sub_mod": hashlib.sha256(gate_id.encode()).hexdigest()},
            }
        elif gate_id in (
            "github_publication_commit_alignment",
            "github_branch_protection_required_ci",
        ):
            raw_api_response = json.dumps({"gate": gate_id, "commit": source_commit})
            receipt["github"] = {
                "repository_url": "https://github.com/liaoyoyo/InterSubMod",
                "api_endpoint": f"https://api.github.com/fixture/{gate_id}",
                "raw_api_response": raw_api_response,
                "response_sha256": hashlib.sha256(raw_api_response.encode()).hexdigest(),
            }
        path = receipt_root / f"gate-{gate_id}.json"
        path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        paths[gate_id] = path
    return paths


def gate_details(gate_id: str, source_commit: str, repo: Path) -> dict[str, object]:
    sha = hashlib.sha256(gate_id.encode()).hexdigest()
    details: dict[str, dict[str, object]] = {
        "registry_and_handoff_package": {"registry_records": 101, "schema_errors": 0, "package_errors": 0},
        "clean_release_build": {"build_type": "Release", "configure_exit": 0, "build_exit": 0},
        "ctest_complete": {"discovery_command": "ctest --test-dir build -N", "tests_total": 270, "tests_passed": 270, "tests_failed": 0},
        "python_tests_complete": {"python_version": "3.10.14", "requirements_lock_sha256": sha, "tests_total": 39, "tests_passed": 39, "tests_failed": 0},
        "tiny_synthetic_e2e": {"plan_exit": 0, "plan_side_effects": 0, "run_exit": 0, "schema_valid": True, "hash_valid": True},
        "bip7_fresh_clone_acceptance": {"fresh_clone": True, "clone_exit": 0, "build_exit": 0, "ctest_exit": 0, "python_tests_exit": 0, "synthetic_smoke_exit": 0, "real_data_preflight_exit": 0, "real_data_read_only": True},
        "bip8_fresh_clone_acceptance": {"fresh_clone": True, "clone_exit": 0, "build_exit": 0, "ctest_exit": 0, "python_tests_exit": 0, "synthetic_smoke_exit": 0, "real_data_preflight_exit": 0, "real_data_read_only": True},
        "html_browser_qa_84": {"checks_total": 84, "checks_passed": 84, "checks_failed": 0, "desktop": True, "mobile": True, "no_js": True, "print": True, "xml_svg_parse": True, "external_requests": 0},
        "docs_claim_link_validation": {"markdown_link_errors": 0, "source_path_errors": 0, "claim_id_errors": 0, "superseded_current_index_refs": 0},
        "public_claim_registry_158": {"rows_total": 158, "confirmed": 124, "confirmed_with_limits": 34, "unverified": 0, "invalid_rows": 0, "p0_unresolved": 0},
        "github_publication_commit_alignment": {"main_commit": source_commit, "wiki_commit": source_commit, "pages_commit": source_commit, "readme_zh_http_status": 200},
        "github_branch_protection_required_ci": {"branch": "main", "protected": True, "required_ci": ["handoff-portable-ci"], "enforcement": True},
        "full_history_secret_scan": {"scope": "FULL_HISTORY", "findings": 0, "report_sha256": sha},
        "release_asset_sha256_roundtrip": {"assets_expected": 1, "assets_downloaded": 1, "hash_matches": 1, "http_200_count": 1},
        "reader_ai_six_question_acceptance": {"detailed_receipt_path": (HANDOFF / "evidence/reader_acceptance_receipt.json").as_posix(), "detailed_receipt_sha256": file_sha256(repo / HANDOFF / "evidence/reader_acceptance_receipt.json"), "schema_path": (HANDOFF / "schemas/reader-acceptance.schema.json").as_posix(), "validator_path": "scripts/handoff/validate_reader_acceptance.py", "questions_total": 6, "questions_passed": 6, "promotions_total": 8, "promotions_passed": 8, "validator_exit": 0},
    }
    if gate_id == "intersubmod_project_license_confirmation":
        raw_response = json.dumps({"user_login": "liaoyoyo", "author_association": "OWNER", "state": "APPROVED", "head_sha": source_commit}, separators=(",", ":"), sort_keys=True)
        return {"authority_kind": "GITHUB_REPOSITORY_OWNER_APPROVED_PR", "maintainer_login": "liaoyoyo", "spdx_expression": PROJECT_LICENSE, "approval_url": "https://github.com/liaoyoyo/InterSubMod/pull/1", "approval_commit": source_commit, "api_endpoint": "https://api.github.com/repos/liaoyoyo/InterSubMod/pulls/1/reviews", "raw_api_response": raw_response, "response_sha256": hashlib.sha256(raw_response.encode()).hexdigest()}
    return details[gate_id]


def write_aggregate_receipt(
    path: Path,
    repo: Path,
    source_commit: str,
    safety_receipt: Path,
    governed_assets: tuple[tuple[Path, str], ...] = (),
) -> None:
    receipt_root = path.parent
    common_receipts = {
        gate_id: receipt_root / f"gate-{gate_id}.json" for gate_id in COMMON_GATE_IDS
    }
    gates: dict[str, object] = {}
    for gate_id in REQUIRED_GATES:
        if gate_id == "authority_19_of_19_replay":
            evidence_uri = AUTHORITY_REPLAY.as_posix()
            evidence_sha = file_sha256(repo / AUTHORITY_REPLAY)
        elif gate_id == "longlineage_preview_safety":
            evidence_uri = f"release-assets/{safety_receipt.name}"
            evidence_sha = file_sha256(safety_receipt)
        elif gate_id == "repo_hygiene":
            evidence_uri = (HANDOFF / "evidence/repo_hygiene_summary_receipt.json").as_posix()
            evidence_sha = file_sha256(repo / evidence_uri)
        elif gate_id == "tracked_large_asset_policy":
            evidence_uri = (HANDOFF / "registries/tracked_large_asset_registry.json").as_posix()
            evidence_sha = file_sha256(repo / evidence_uri)
        else:
            gate_receipt = common_receipts[gate_id]
            evidence_uri = f"release-assets/{gate_receipt.name}"
            evidence_sha = file_sha256(gate_receipt)
        gates[gate_id] = {
            "status": "PASS",
            "source_commit": source_commit,
            "verified_at": VERIFIED_AT,
            "evidence": [
                {
                    "artifact_id": f"acceptance/{gate_id}",
                    "evidence_type": EVIDENCE_TYPES[gate_id],
                    "subject_commit": source_commit,
                    "sha256": evidence_sha,
                    "uri": evidence_uri,
                    "verified_at": VERIFIED_AT,
                }
            ],
        }
        if gate_id in COMMON_GATE_IDS:
            log_path = receipt_root / f"gate-{gate_id}.log"
            gate = gates[gate_id]
            assert isinstance(gate, dict)
            evidence = gate["evidence"]
            assert isinstance(evidence, list)
            evidence.append(
                {
                    "artifact_id": f"acceptance/{gate_id}/raw-log",
                    "evidence_type": EVIDENCE_TYPES[gate_id],
                    "subject_commit": source_commit,
                    "sha256": file_sha256(log_path),
                    "uri": f"release-assets/{log_path.name}",
                    "verified_at": VERIFIED_AT,
                }
            )
        if gate_id == "reader_ai_six_question_acceptance":
            reader_uri = (HANDOFF / "evidence/reader_acceptance_receipt.json").as_posix()
            gate = gates[gate_id]
            assert isinstance(gate, dict)
            evidence = gate["evidence"]
            assert isinstance(evidence, list)
            evidence.append(
                {
                    "artifact_id": "acceptance/reader-detailed-receipt",
                    "evidence_type": EVIDENCE_TYPES[gate_id],
                    "subject_commit": source_commit,
                    "sha256": file_sha256(repo / reader_uri),
                    "uri": reader_uri,
                    "verified_at": VERIFIED_AT,
                }
            )
        if gate_id == "tracked_large_asset_policy":
            migration_uri = (HANDOFF / "evidence/large_asset_migration_receipt.json").as_posix()
            gate = gates[gate_id]
            assert isinstance(gate, dict)
            evidence = gate["evidence"]
            assert isinstance(evidence, list)
            evidence.append(
                {
                    "artifact_id": "acceptance/large_asset_migration",
                    "evidence_type": EVIDENCE_TYPES[gate_id],
                    "subject_commit": source_commit,
                    "sha256": file_sha256(repo / migration_uri),
                    "uri": migration_uri,
                    "verified_at": VERIFIED_AT,
                }
            )
    for source_path, published_name in governed_assets:
        for gate_id in ("release_asset_sha256_roundtrip", "full_history_secret_scan"):
            gate = gates[gate_id]
            assert isinstance(gate, dict)
            evidence = gate["evidence"]
            assert isinstance(evidence, list)
            evidence.append(
                {
                    "artifact_id": f"release_asset/{published_name}",
                    "evidence_type": EVIDENCE_TYPES[gate_id],
                    "subject_commit": source_commit,
                    "sha256": file_sha256(source_path),
                    "uri": f"release-assets/{published_name}",
                    "verified_at": VERIFIED_AT,
                }
            )
    receipt = {
        "schema_version": "1.0.0",
        "receipt_type": "aggregate_release_acceptance",
        "repository_url": "https://github.com/liaoyoyo/InterSubMod",
        "source_commit": source_commit,
        "generated_at": VERIFIED_AT,
        "overall_pass": True,
        "blocking_gates": [],
        "project_license_spdx": PROJECT_LICENSE,
        "gates": gates,
    }
    path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")


def approval_receipts(
    tmp_path: Path,
    repo: Path,
    commit: str,
    governed_assets: tuple[tuple[Path, str], ...] = (),
) -> tuple[Path, Path]:
    receipt_root = tmp_path / "receipts"
    receipt_root.mkdir()
    safety_receipt = receipt_root / "longlineage-preview-safety.json"
    write_longlineage_safety_receipt(safety_receipt, commit, repo)
    write_common_gate_receipts(receipt_root, commit, repo)
    aggregate_receipt = receipt_root / "aggregate-release-acceptance.json"
    write_aggregate_receipt(
        aggregate_receipt,
        repo,
        commit,
        safety_receipt,
        governed_assets=governed_assets,
    )
    return aggregate_receipt, safety_receipt


def register_release_asset(
    repo: Path,
    source_path: Path,
    artifact_id: str,
    published_name: str,
    **overrides: object,
) -> str:
    registry_path = repo / HANDOFF / "registries/artifact_registry.json"
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    record = dict(registry[0])
    record.update(
        {
            "artifact_id": artifact_id,
            "artifact_type": "governed_release_asset",
            "availability": "GITHUB_RELEASE",
            "evidence_status": "VALIDATED_DERIVED",
            "finality": "FINAL_FOR_SCOPE",
            "license": PROJECT_LICENSE,
            "public_uri": f"release-assets/{published_name}",
            "sha256": file_sha256(source_path),
            "size_bytes": source_path.stat().st_size,
        }
    )
    record.update(overrides)
    registry.append(record)
    registry_path.write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    assert run(["git", "add", registry_path.relative_to(repo).as_posix()], repo).returncode == 0
    assert run(["git", "commit", "-qm", f"register release asset {artifact_id}"], repo).returncode == 0
    tested_commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    write_reader_receipt(repo, tested_commit)
    reader_path = HANDOFF / "evidence/reader_acceptance_receipt.json"
    assert run(["git", "add", reader_path.as_posix()], repo).returncode == 0
    assert run(["git", "commit", "-qm", "refresh reader acceptance receipt"], repo).returncode == 0
    return run(["git", "rev-parse", "HEAD"], repo).stdout.strip()


def approved_command(
    repo: Path,
    output_root: Path,
    aggregate_receipt: Path,
    safety_receipt: Path,
) -> list[str]:
    command = [
        sys.executable,
        str(SCRIPT),
        "--repo-root",
        str(repo),
        "--status",
        "APPROVED_RESEARCH_HANDOFF",
        "--project-license-spdx",
        PROJECT_LICENSE,
        "--approval-receipt",
        str(aggregate_receipt),
        "--longlineage-preview-safety-status",
        "VERIFIED",
        "--longlineage-safety-receipt",
        str(safety_receipt),
        "--longlineage-candidate",
        LONG_LINEAGE_COMMIT,
        "--longlineage-repo-root",
        str(repo.parent / "longlineage"),
        "--git-asset",
        "README.md",
        "--generated-at",
        VERIFIED_AT,
        "--output-json",
        str(output_root / "release-manifest.json"),
        "--output-checksums",
        str(output_root / "SHA256SUMS"),
    ]
    for receipt in sorted(aggregate_receipt.parent.glob("gate-*.json")):
        command.extend(["--gate-evidence-asset", str(receipt)])
    for log_path in sorted(aggregate_receipt.parent.glob("gate-*.log")):
        command.extend(["--gate-log-asset", str(log_path)])
    return command


def test_blocked_manifest_binds_clean_commit_and_validates(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    output_root = tmp_path / "release-output"
    result = run(blocked_command(repo, output_root), repo)
    assert result.returncode == 0, result.stderr
    assert (
        f"[RESULT] release_ready=false assets={len(REQUIRED_GIT_ASSETS)} blockers=2"
        in result.stdout
    )

    manifest = json.loads((output_root / "release-manifest.json").read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    validator_class = validator_for(schema)
    validator_class.check_schema(schema)
    errors = list(validator_class(schema, format_checker=FormatChecker()).iter_errors(manifest))
    assert errors == []
    assert manifest["source"]["commit"] == commit
    assert manifest["release_status"] == "BLOCKED"
    assert manifest["release_ready"] is False
    assert manifest["project_license"] == {
        "status": "NEEDS_MAINTAINER_CONFIRMATION",
        "spdx_expression": None,
        "authority_observation": (
            "GPL-3.0 license text present; project-level SPDX choice requires "
            "maintainer confirmation."
        ),
    }
    assert manifest["authority"]["claim_ceiling"]["confirmed_cellular_subclone"] == 0
    assert manifest["authority"]["claim_ceiling"]["confirmed_linear_ancestry"] == 0
    assert set(manifest["longlineage"]["gates"].values()) == {"BLOCKED"}
    git_paths = {
        asset["repository_path"]
        for asset in manifest["assets"]
        if asset["publication_channel"] == "GIT"
    }
    assert git_paths == {path.as_posix() for path in REQUIRED_GIT_ASSETS}

    manifest_digest = hashlib.sha256((output_root / "release-manifest.json").read_bytes()).hexdigest()
    checksum_lines = (output_root / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
    assert checksum_lines[0] == f"{manifest_digest}  release-manifest.json"
    assert set(checksum_lines[1:]) == {
        f"{asset['sha256']}  {asset['published_path']}" for asset in manifest["assets"]
    }
    assert manifest["checksum_sidecar"]["line_count"] == len(REQUIRED_GIT_ASSETS) + 1


def test_dirty_repository_fails_before_writing_outputs(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    (repo / "untracked.txt").write_text("dirty\n", encoding="utf-8")
    output_root = tmp_path / "release-output"
    result = run(blocked_command(repo, output_root), repo)
    assert result.returncode == 2
    assert "repository is not clean" in result.stderr
    assert not output_root.exists()


def test_approved_status_requires_aggregate_receipt(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    output_root = tmp_path / "release-output"
    command = blocked_command(repo, output_root)
    command[command.index("BLOCKED")] = "APPROVED_RESEARCH_HANDOFF"
    while "--blocker" in command:
        blocker_index = command.index("--blocker")
        del command[blocker_index : blocker_index + 2]
    status_index = command.index("--longlineage-preview-safety-status") + 1
    command[status_index] = "VERIFIED"
    result = run(command, repo)
    assert result.returncode == 2
    assert "approved handoff requires exactly one --approval-receipt" in result.stderr
    assert not output_root.exists()


def test_approved_handoff_keeps_longlineage_scientific_phases_blocked(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    output_root = tmp_path / "release-output"
    command = approved_command(repo, output_root, aggregate_receipt, safety_receipt)
    result = run(command, repo)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((output_root / "release-manifest.json").read_text(encoding="utf-8"))
    assert manifest["release_ready"] is True
    assert manifest["blockers"] == []
    assert manifest["project_license"]["spdx_expression"] == PROJECT_LICENSE
    assert manifest["longlineage"]["preview_status"] == "RESEARCH_PREVIEW_SAFETY_VERIFIED"
    assert manifest["longlineage"]["preview_safety_status"] == "VERIFIED"
    assert manifest["longlineage"]["production_ready"] is False
    assert "not InterSubMod handoff blockers" in manifest["longlineage"]["phase_gate_interpretation"]
    assert manifest["longlineage"]["gates"] == {
        "P3": "BLOCKED",
        "P4": "BLOCKED",
        "P5": "BLOCKED",
        "P7": "BLOCKED",
        "P8": "BLOCKED",
    }
    assert len(manifest["approval_receipts"]) == 1
    assert manifest["longlineage"]["safety_receipt"] is not None
    assert manifest["asset_count"] == (
        len(REQUIRED_GIT_ASSETS)
        + len(APPROVED_ONLY_GIT_ASSETS)
        + len(COMMON_GATE_IDS)
        + len(COMMON_GATE_IDS)
        + 2
    )
    approved_git_paths = {
        asset["repository_path"]
        for asset in manifest["assets"]
        if asset["publication_channel"] == "GIT"
    }
    assert approved_git_paths == {
        path.as_posix() for path in REQUIRED_GIT_ASSETS + APPROVED_ONLY_GIT_ASSETS
    } | DYNAMIC_GATE_GIT_EVIDENCE

    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    validator_class = validator_for(schema)
    validator = validator_class(schema, format_checker=FormatChecker())
    assert list(validator.iter_errors(manifest)) == []
    manifest["longlineage"]["gates"]["P3"] = "VERIFIED"
    assert list(validator.iter_errors(manifest)), "schema must reject promotion of LongLineage P3"


def test_fake_minimal_aggregate_receipt_is_rejected(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    aggregate_receipt.write_text(
        json.dumps(
            {
                "receipt_type": "aggregate_release_acceptance",
                "source_commit": commit,
                "overall_pass": True,
                "blocking_gates": [],
            }
        ),
        encoding="utf-8",
    )
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "failed schema" in result.stderr


def test_aggregate_receipt_missing_required_gate_is_rejected(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    receipt = json.loads(aggregate_receipt.read_text(encoding="utf-8"))
    del receipt["gates"]["bip8_fresh_clone_acceptance"]
    aggregate_receipt.write_text(json.dumps(receipt), encoding="utf-8")
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "bip8_fresh_clone_acceptance" in result.stderr


def test_aggregate_receipt_extra_gate_is_rejected(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    receipt = json.loads(aggregate_receipt.read_text(encoding="utf-8"))
    receipt["gates"]["self_declared_release_ready"] = receipt["gates"]["clean_release_build"]
    aggregate_receipt.write_text(json.dumps(receipt), encoding="utf-8")
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "Additional properties are not allowed" in result.stderr


def test_aggregate_receipt_must_hash_bind_authority_replay(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    receipt = json.loads(aggregate_receipt.read_text(encoding="utf-8"))
    receipt["gates"]["authority_19_of_19_replay"]["evidence"][0]["sha256"] = "0" * 64
    aggregate_receipt.write_text(json.dumps(receipt), encoding="utf-8")
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "does not hash-bind required receipt" in result.stderr


def test_fake_minimal_longlineage_safety_receipt_is_rejected(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    safety_receipt.write_text(
        json.dumps(
            {
                "receipt_type": "longlineage_public_preview_safety",
                "source_commit": commit,
                "longlineage_candidate_commit": LONG_LINEAGE_COMMIT,
                "overall_pass": True,
                "blocking_gates": [],
            }
        ),
        encoding="utf-8",
    )
    write_aggregate_receipt(aggregate_receipt, repo, commit, safety_receipt)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "failed schema" in result.stderr


def test_wrong_longlineage_candidate_argument_is_rejected(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    command = blocked_command(repo, tmp_path / "release-output")
    candidate_index = command.index("--longlineage-candidate") + 1
    command[candidate_index] = "0" * 40
    result = run(command, repo)
    assert result.returncode == 2
    assert f"must equal pinned commit {LONG_LINEAGE_COMMIT}" in result.stderr


def test_wrong_longlineage_candidate_in_safety_receipt_is_rejected(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    receipt = json.loads(safety_receipt.read_text(encoding="utf-8"))
    receipt["longlineage_candidate_commit"] = "0" * 40
    safety_receipt.write_text(json.dumps(receipt), encoding="utf-8")
    write_aggregate_receipt(aggregate_receipt, repo, commit, safety_receipt)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert LONG_LINEAGE_COMMIT in result.stderr


def test_longlineage_safety_receipt_must_hash_bind_supporting_documents(
    tmp_path: Path,
) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    receipt = json.loads(safety_receipt.read_text(encoding="utf-8"))
    receipt["preview_safety_gates"]["dependencies"]["evidence"][0]["sha256"] = "0" * 64
    safety_receipt.write_text(json.dumps(receipt), encoding="utf-8")
    write_aggregate_receipt(aggregate_receipt, repo, commit, safety_receipt)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "does not bind supporting bytes: dependencies" in result.stderr


def test_longlineage_fail_closed_public_safety_evidence_prevents_verification(
    tmp_path: Path,
) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    public_safety = HANDOFF / "evidence/longlineage_public_safety_receipt.json"

    def fail_safety(raw: bytes) -> bytes:
        receipt = json.loads(raw)
        receipt["gate_status"] = "FAIL_CLOSED"
        return (json.dumps(receipt) + "\n").encode("utf-8")

    commit_tamper(repo, public_safety, fail_safety)
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    receipt_parent = tmp_path / "ll-failed"
    receipt_parent.mkdir()
    aggregate_receipt, safety_receipt = approval_receipts(receipt_parent, repo, commit)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "public-safety evidence remains fail-closed" in result.stderr


def test_approved_license_rejects_conflicting_source_spdx_header(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    source = repo / "src/conflicting.cpp"
    source.parent.mkdir()
    source.write_text(
        "// SPDX-License-Identifier: GPL-3.0-only\nint fixture = 0;\n",
        encoding="utf-8",
    )
    assert run(["git", "add", source.relative_to(repo).as_posix()], repo).returncode == 0
    assert run(["git", "commit", "-qm", "add conflicting SPDX fixture"], repo).returncode == 0
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "source SPDX headers conflict" in result.stderr


def commit_tamper(
    repo: Path, relative_path: Path, transform: Callable[[bytes], bytes]
) -> None:
    path = repo / relative_path
    path.write_bytes(transform(path.read_bytes()))
    assert run(["git", "add", relative_path.as_posix()], repo).returncode == 0
    assert run(["git", "commit", "-qm", f"tamper {relative_path.name}"], repo).returncode == 0


def test_tampered_frozen_authority_is_rejected(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    authority = HANDOFF / "evidence/authority_manifest.json"
    commit_tamper(repo, authority, lambda raw: raw + b"\n")
    result = run(blocked_command(repo, tmp_path / "release-output"), repo)
    assert result.returncode == 2
    assert "authority manifest bytes do not match" in result.stderr


def test_tampered_frozen_denominator_registry_is_rejected(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    denominator = HANDOFF / "evidence/denominator_registry.tsv"
    commit_tamper(repo, denominator, lambda raw: raw.replace(b"88.2579", b"88.2580", 1))
    result = run(blocked_command(repo, tmp_path / "release-output"), repo)
    assert result.returncode == 2
    assert "denominator registry bytes do not match" in result.stderr


def test_tampered_authority_replay_is_rejected(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)

    def fail_one_replay(raw: bytes) -> bytes:
        receipt = json.loads(raw)
        receipt["results"][0]["status"] = "HASH_MISMATCH"
        return (json.dumps(receipt) + "\n").encode("utf-8")

    commit_tamper(repo, AUTHORITY_REPLAY, fail_one_replay)
    result = run(blocked_command(repo, tmp_path / "release-output"), repo)
    assert result.returncode == 2
    assert "authority replay result is not MATCH" in result.stderr


def test_self_consistent_but_fabricated_authority_replay_rows_are_rejected(
    tmp_path: Path,
) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)

    def fabricate_identity(raw: bytes) -> bytes:
        receipt = json.loads(raw)
        first = receipt["results"][0]
        first["artifact_id"] = "fabricated_unique_artifact"
        first["path"] = "/fabricated/authority/result.json"
        first["expected_sha256"] = "f" * 64
        first["actual_sha256"] = "f" * 64
        return (json.dumps(receipt) + "\n").encode("utf-8")

    commit_tamper(repo, AUTHORITY_REPLAY, fabricate_identity)
    result = run(blocked_command(repo, tmp_path / "release-output"), repo)
    assert result.returncode == 2
    assert "do not exactly match the frozen 13 artifacts + binary + 5 source snapshots" in result.stderr


def test_semantically_empty_gate_receipt_is_rejected(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    gate_path = aggregate_receipt.parent / "gate-clean_release_build.json"
    gate_path.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "receipt_type": "clean_release_build_acceptance",
                "gate_id": "clean_release_build",
                "source_commit": commit,
                "overall_pass": True,
                "blocking_gates": [],
            }
        ),
        encoding="utf-8",
    )
    aggregate = json.loads(aggregate_receipt.read_text(encoding="utf-8"))
    aggregate["gates"]["clean_release_build"]["evidence"][0]["sha256"] = file_sha256(
        gate_path
    )
    aggregate_receipt.write_text(json.dumps(aggregate), encoding="utf-8")
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "failed schema" in result.stderr
    assert "'execution' is a required property" in result.stderr


def test_tracked_large_policy_blocker_prevents_approval(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    registry_path = HANDOFF / "registries/tracked_large_asset_registry.json"

    def block_registry(raw: bytes) -> bytes:
        registry = json.loads(raw)
        registry["summary"]["verdict"] = "LARGE_ASSET_MIGRATION_BLOCKED"
        return (json.dumps(registry) + "\n").encode("utf-8")

    commit_tamper(repo, registry_path, block_registry)
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    receipt_parent = tmp_path / "after-tamper"
    receipt_parent.mkdir()
    aggregate_receipt, safety_receipt = approval_receipts(receipt_parent, repo, commit)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "tracked-large asset policy is not PASS" in result.stderr


def test_actual_git_blob_over_one_mib_prevents_approval(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    large_path = repo / "docs/too-large.html"
    large_path.parent.mkdir(exist_ok=True)
    large_path.write_bytes(b"x" * (1024 * 1024 + 1))
    assert run(["git", "add", large_path.relative_to(repo).as_posix()], repo).returncode == 0
    assert run(["git", "commit", "-qm", "add oversized tracked fixture"], repo).returncode == 0
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    receipt_parent = tmp_path / "oversized"
    receipt_parent.mkdir()
    aggregate_receipt, safety_receipt = approval_receipts(receipt_parent, repo, commit)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "recomputation found files over 1 MiB" in result.stderr


def test_governed_external_release_asset_requires_registry_and_two_gate_bindings(
    tmp_path: Path,
) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    external_asset = tmp_path / "final-report.html"
    external_asset.write_text("<!doctype html><title>final report</title>\n", encoding="utf-8")
    artifact_id = "release.final_report"
    commit = register_release_asset(
        repo,
        external_asset,
        artifact_id,
        published_name=external_asset.name,
    )
    aggregate_receipt, safety_receipt = approval_receipts(
        tmp_path,
        repo,
        commit,
        governed_assets=((external_asset, external_asset.name),),
    )
    output_root = tmp_path / "release-output"
    command = approved_command(repo, output_root, aggregate_receipt, safety_receipt)
    command.extend(
        ["--release-asset", f"{artifact_id}::{external_asset.name}={external_asset}"]
    )
    result = run(command, repo)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((output_root / "release-manifest.json").read_text(encoding="utf-8"))
    governed = [asset for asset in manifest["assets"] if "registry_artifact_id" in asset]
    assert governed == [
        {
            "asset_id": governed[0]["asset_id"],
            "publication_channel": "GITHUB_RELEASE",
            "published_path": f"release-assets/{external_asset.name}",
            "registry_artifact_id": artifact_id,
            "release_source_name": external_asset.name,
            "sha256": file_sha256(external_asset),
            "size_bytes": external_asset.stat().st_size,
        }
    ]


def test_unregistered_external_release_asset_is_rejected(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    external_asset = tmp_path / "unregistered.html"
    external_asset.write_text("unregistered\n", encoding="utf-8")
    aggregate_receipt, safety_receipt = approval_receipts(
        tmp_path,
        repo,
        commit,
        governed_assets=((external_asset, external_asset.name),),
    )
    command = approved_command(
        repo, tmp_path / "release-output", aggregate_receipt, safety_receipt
    )
    command.extend(
        ["--release-asset", f"release.not_registered::{external_asset}"]
    )
    result = run(command, repo)
    assert result.returncode == 2
    assert "no governed artifact record" in result.stderr


def test_raw_bam_release_asset_is_rejected_even_if_registered(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    external_asset = tmp_path / "raw.bam"
    external_asset.write_bytes(b"not-a-real-bam")
    artifact_id = "release.raw_bam"
    commit = register_release_asset(repo, external_asset, artifact_id, external_asset.name)
    aggregate_receipt, safety_receipt = approval_receipts(
        tmp_path,
        repo,
        commit,
        governed_assets=((external_asset, external_asset.name),),
    )
    command = approved_command(
        repo, tmp_path / "release-output", aggregate_receipt, safety_receipt
    )
    command.extend(["--release-asset", f"{artifact_id}::{external_asset}"])
    result = run(command, repo)
    assert result.returncode == 2
    assert "raw, runtime, or secret-like release asset is forbidden" in result.stderr


def test_nonfinal_release_asset_registry_record_is_rejected(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    external_asset = tmp_path / "partial.tsv"
    external_asset.write_text("metric\tvalue\npartial\t1\n", encoding="utf-8")
    artifact_id = "release.partial_table"
    commit = register_release_asset(
        repo,
        external_asset,
        artifact_id,
        external_asset.name,
        finality="NON_FINAL",
    )
    aggregate_receipt, safety_receipt = approval_receipts(
        tmp_path,
        repo,
        commit,
        governed_assets=((external_asset, external_asset.name),),
    )
    command = approved_command(
        repo, tmp_path / "release-output", aggregate_receipt, safety_receipt
    )
    command.extend(["--release-asset", f"{artifact_id}::{external_asset}"])
    result = run(command, repo)
    assert result.returncode == 2
    assert "not FINAL_FOR_SCOPE" in result.stderr


def test_release_asset_must_be_bound_by_secret_scan_gate(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    external_asset = tmp_path / "figure-bundle.zip"
    external_asset.write_bytes(b"fixture bundle")
    artifact_id = "release.figure_bundle"
    commit = register_release_asset(repo, external_asset, artifact_id, external_asset.name)
    aggregate_receipt, safety_receipt = approval_receipts(
        tmp_path,
        repo,
        commit,
        governed_assets=((external_asset, external_asset.name),),
    )
    receipt = json.loads(aggregate_receipt.read_text(encoding="utf-8"))
    receipt["gates"]["full_history_secret_scan"]["evidence"] = receipt["gates"][
        "full_history_secret_scan"
    ]["evidence"][:1]
    aggregate_receipt.write_text(json.dumps(receipt), encoding="utf-8")
    command = approved_command(
        repo, tmp_path / "release-output", aggregate_receipt, safety_receipt
    )
    command.extend(["--release-asset", f"{artifact_id}::{external_asset}"])
    result = run(command, repo)
    assert result.returncode == 2
    assert "full_history_secret_scan" in result.stderr


def test_external_gate_receipt_must_bind_commit_and_be_referenced(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    host_receipt = tmp_path / "bip8-fresh-clone-receipt.json"
    original = aggregate_receipt.parent / "gate-bip8_fresh_clone_acceptance.json"
    host_receipt.write_bytes(original.read_bytes())
    original.unlink()
    receipt = json.loads(aggregate_receipt.read_text(encoding="utf-8"))
    evidence = receipt["gates"]["bip8_fresh_clone_acceptance"]["evidence"][0]
    evidence["sha256"] = file_sha256(host_receipt)
    evidence["uri"] = f"release-assets/{host_receipt.name}"
    aggregate_receipt.write_text(json.dumps(receipt), encoding="utf-8")
    output_root = tmp_path / "release-output"
    command = approved_command(repo, output_root, aggregate_receipt, safety_receipt)
    command.extend(["--gate-evidence-asset", str(host_receipt)])
    result = run(command, repo)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((output_root / "release-manifest.json").read_text(encoding="utf-8"))
    assert any(
        asset.get("gate_evidence_receipt") is True
        and asset["published_path"] == f"release-assets/{host_receipt.name}"
        for asset in manifest["assets"]
    )


def test_unreferenced_external_gate_receipt_is_rejected(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    receipt_path = tmp_path / "unused-gate-receipt.json"
    receipt_path.write_bytes(
        (aggregate_receipt.parent / "gate-clean_release_build.json").read_bytes()
    )
    command = approved_command(
        repo, tmp_path / "release-output", aggregate_receipt, safety_receipt
    )
    command.extend(["--gate-evidence-asset", str(receipt_path)])
    result = run(command, repo)
    assert result.returncode == 2
    assert "unreferenced --gate-evidence-asset" in result.stderr


def test_output_inside_repository_is_rejected(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    command = blocked_command(repo, repo / "release-output")
    result = run(command, repo)
    assert result.returncode == 2
    assert "must stay outside Git to avoid self-reference" in result.stderr
    assert not (repo / "release-output").exists()


def test_generic_self_reported_gate_receipt_is_rejected(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    gate_path = aggregate_receipt.parent / "gate-clean_release_build.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    gate["execution"] = {"kind": "COMMAND", "locator": "true", "exit_code": 0}
    gate["result"] = {"status": "PASS", "checks": 1, "failures": 0}
    gate_path.write_text(json.dumps(gate), encoding="utf-8")
    aggregate = json.loads(aggregate_receipt.read_text(encoding="utf-8"))
    aggregate["gates"]["clean_release_build"]["evidence"][0]["sha256"] = file_sha256(
        gate_path
    )
    aggregate_receipt.write_text(json.dumps(aggregate), encoding="utf-8")
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "failed schema" in result.stderr
    assert "'tool' is a required property" in result.stderr


def test_recursive_broken_symlink_chain_is_rejected(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    (repo / "link-a").symlink_to("link-b")
    (repo / "link-b").symlink_to("missing-target")
    assert run(["git", "add", "link-a", "link-b"], repo).returncode == 0
    assert run(["git", "commit", "-qm", "add recursive broken symlink fixture"], repo).returncode == 0
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "broken symlink chain remains" in result.stderr


def test_symlink_cycle_is_rejected(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)
    (repo / "cycle-a").symlink_to("cycle-b")
    (repo / "cycle-b").symlink_to("cycle-a")
    assert run(["git", "add", "cycle-a", "cycle-b"], repo).returncode == 0
    assert run(["git", "commit", "-qm", "add symlink cycle fixture"], repo).returncode == 0
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "symlink cycle remains" in result.stderr


def test_license_gate_rejects_arbitrary_maintainer_identity(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    gate_path = aggregate_receipt.parent / "gate-intersubmod_project_license_confirmation.json"
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    gate["details"]["maintainer_login"] = "self-appointed-maintainer"
    gate_path.write_text(json.dumps(gate), encoding="utf-8")
    aggregate = json.loads(aggregate_receipt.read_text(encoding="utf-8"))
    aggregate["gates"]["intersubmod_project_license_confirmation"]["evidence"][0][
        "sha256"
    ] = file_sha256(gate_path)
    aggregate_receipt.write_text(json.dumps(aggregate), encoding="utf-8")
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "failed schema" in result.stderr
    assert "liaoyoyo" in result.stderr


def test_reader_gate_executes_commit_pinned_semantic_validator(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)

    def generic_reader(raw: bytes) -> bytes:
        receipt = json.loads(raw)
        receipt["questions"][0]["answer"] = (
            "Everything in this package looks fine and understandable to a new reader, "
            "and I have no additional project-specific observations to record here."
        )
        return (json.dumps(receipt) + "\n").encode("utf-8")

    commit_tamper(repo, HANDOFF / "evidence/reader_acceptance_receipt.json", generic_reader)
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "commit-pinned reader validator did not PASS" in result.stderr


def test_algorithm_crosswalk_static_release_asset_gate_must_pass(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)

    def block_crosswalk(raw: bytes) -> bytes:
        crosswalk = json.loads(raw)
        crosswalk["gates"]["RELEASE_ASSET_READY"]["status"] = "BLOCKED"
        return (json.dumps(crosswalk) + "\n").encode("utf-8")

    commit_tamper(
        repo,
        HANDOFF / "evidence/algorithm_cli_claim_crosswalk.json",
        block_crosswalk,
    )
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "static asset gate is not PASS: RELEASE_ASSET_READY" in result.stderr


def test_longlineage_safety_head_must_be_advertised_by_origin(tmp_path: Path) -> None:
    repo, _ = initialize_fixture_repository(tmp_path)

    def point_to_missing_branch(raw: bytes) -> bytes:
        receipt = json.loads(raw)
        receipt["remote_branch"] = "agent/not-advertised"
        return (json.dumps(receipt) + "\n").encode("utf-8")

    commit_tamper(
        repo,
        HANDOFF / "evidence/longlineage_public_gate_receipt.json",
        point_to_missing_branch,
    )
    commit = run(["git", "rev-parse", "HEAD"], repo).stdout.strip()
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "not advertised by an origin branch" in result.stderr


def test_gate_execution_log_bytes_must_match_receipt_and_aggregate(tmp_path: Path) -> None:
    repo, commit = initialize_fixture_repository(tmp_path)
    aggregate_receipt, safety_receipt = approval_receipts(tmp_path, repo, commit)
    log_path = aggregate_receipt.parent / "gate-clean_release_build.log"
    log_path.write_text("tampered raw execution log\n", encoding="utf-8")
    result = run(
        approved_command(repo, tmp_path / "release-output", aggregate_receipt, safety_receipt),
        repo,
    )
    assert result.returncode == 2
    assert "aggregate gate evidence bytes are not declared" in result.stderr
