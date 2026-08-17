from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/handoff/validate_reader_acceptance.py"
SCHEMA = ROOT / "docs/handoff/20260813_完整研究資料與軟體交接_01/schemas/reader-acceptance.schema.json"
SPEC = importlib.util.spec_from_file_location("reader_acceptance_validator", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


ANSWERS = {
    "Q_PROJECT": (
        "This research handoff snapshot covers ONT long-read sSNV read linkage, candidate mutation-state "
        "reconstruction, and methylation association. It is not a production release or cellular-lineage caller."
    ),
    "Q_CONCLUSION": (
        "Confirmed cellular subclones = 0 and confirmed linear ancestry relationships = 0. Methylation is "
        "association only. In the frozen reconstruction CN/LOH is not integrated. 88.2579% is a "
        "model-conditional rooted graph-shape statistic; it is not accuracy and it is not prevalence."
    ),
    "Q_FINALITY": (
        "An artifact can be FINAL_FOR_SCOPE only after its evidence_status, finality, producer, scope and hash "
        "are checked, while supersedes identifies replaced material. Final is bounded only to the declared "
        "scope; it is not production-ready or whole research final."
    ),
    "Q_SOFTWARE_ROLES": (
        "LongPhase-S/TO produces HP/PS, phased or recalibrated VCF and tagged BAM in one chain. In a second "
        "parallel provenance chain, exact-PS/LongLineage produces candidate families and read assignments. "
        "InterSubMod produces per-region methylation, distance, read clustering and statistics. A commit-pinned "
        "Python research solver can be a science producer; a validator/publication builder or HTML presenter "
        "validates and presents validated data and does not recompute science."
    ),
    "Q_MACHINES": (
        "On bip7 a bounded local preflight exists, but fresh-clone validation is still pending and BLOCKED. "
        "bip8 is BLOCKED because there is no host-local receipt. Each host-specific receipt cannot substitute "
        "for or prove the other host. Run doctor, then build/test, then synthetic smoke on each machine."
    ),
    "Q_VERIFY_CONTINUE": (
        "Replay authority hashes: 19/19 byte SHA-256 records MATCH; this 19/19 replay does not constitute a "
        "science rerun. Validate registries for schema and unique IDs, then respect publication/release/machine "
        "gates that remain BLOCKED. For a new research cycle, run pre-decision audit and pin commit, input hash, "
        "scope and denominator."
    ),
}

PROMOTION_EXPLANATIONS = {
    "NO_CELLULAR_PROMOTION": (
        "Confirmed cellular subclones = 0; candidate read groups are not cellular truth and cannot be promoted "
        "to validated cellular subclones."
    ),
    "NO_ANCESTRY_PROMOTION": (
        "Confirmed linear ancestry = 0; a graph shape or read dendrogram is not ancestry and cannot establish "
        "a cellular phylogeny."
    ),
    "NO_882579_ACCURACY_OR_PREVALENCE": (
        "88.2579% is a model-conditional rooted graph-shape statistic; it is not accuracy and it is not prevalence "
        "or biological truth."
    ),
    "NO_TECHNICAL_DATASET_AS_BIOLOGICAL_REPLICATE": (
        "The frozen scope contains 7 technical datasets and 6 biological IDs; they are not 7 independent "
        "biological replicates or samples."
    ),
    "NO_FEATURE_AS_MAIN": (
        "b9aaa12 is a research preview and non-production candidate; P3/P4/P5/P7/P8 remain BLOCKED, and a "
        "feature is not main, release, or production capability."
    ),
    "NO_LOCAL_AS_LIVE_PUBLISHED": (
        "Local source validation is not live publication; main, Wiki and Pages remain BLOCKED until publish "
        "and re-fetch verification is completed."
    ),
    "NO_BIP7_AS_BIP8": (
        "A bip7 receipt cannot prove or substitute for bip8; bip8 is BLOCKED until a host-local receipt is "
        "created on that hostname."
    ),
    "NO_PYTHON_ROLE_CONFLATION": (
        "A commit-pinned Python research solver can be a science producer; a validator/publication builder or "
        "HTML presenter validates and presents, and the presenter does not recompute science."
    ),
}


def run(*args: str, cwd: Path) -> str:
    result = subprocess.run(args, cwd=cwd, check=True, capture_output=True, text=True)
    return result.stdout.strip()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_repository(tmp_path: Path, *, descendant_commit: bool = False) -> tuple[Path, Path, str]:
    repo = tmp_path / "repo"
    package = repo / MODULE.PACKAGE_RELATIVE
    for relative in MODULE.REQUIRED_MANIFEST_PATHS:
        path = package / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if relative == MODULE.SCHEMA_RELATIVE.as_posix():
            path.write_bytes(SCHEMA.read_bytes())
        else:
            path.write_text(f"fixture evidence for {relative}\n", encoding="utf-8")
    run("git", "init", "-q", cwd=repo)
    run("git", "config", "user.email", "reader-test@example.invalid", cwd=repo)
    run("git", "config", "user.name", "Reader Test", cwd=repo)
    run("git", "add", ".", cwd=repo)
    run("git", "commit", "-q", "-m", "reader source candidate", cwd=repo)
    tested_commit = run("git", "rev-parse", "HEAD", cwd=repo)
    if descendant_commit:
        (repo / "UNRELATED.md").write_text("descendant commit without package-byte changes\n", encoding="utf-8")
        run("git", "add", "UNRELATED.md", cwd=repo)
        run("git", "commit", "-q", "-m", "unrelated descendant", cwd=repo)
    return repo, package, tested_commit


def make_receipt(package: Path, tested_commit: str) -> dict:
    manifest = [
        {"path": relative, "sha256": sha256(package / relative)}
        for relative in sorted(MODULE.REQUIRED_MANIFEST_PATHS)
    ]
    questions = [
        {
            "question_id": question_id,
            "status": "PASS",
            "answer": ANSWERS[question_id],
            "evidence_paths": sorted(MODULE.QUESTION_EVIDENCE[question_id]),
            "limits": ["This answer is bounded to the commit- and hash-bound research handoff snapshot."],
        }
        for question_id in MODULE.QUESTIONS
    ]
    promotions = [
        {
            "check_id": check_id,
            "status": "PASS",
            "explanation": PROMOTION_EXPLANATIONS[check_id],
            "evidence_paths": sorted(MODULE.PROMOTION_EVIDENCE[check_id]),
        }
        for check_id in MODULE.PROMOTIONS
    ]
    return {
        "schema_name": "intersubmod.reader_acceptance_receipt",
        "schema_version": "2.0.0",
        "receipt_id": "intersubmod-reader-acceptance-20260813",
        "evaluated_at": "2026-08-13T20:15:00+08:00",
        "tested_git_commit": tested_commit,
        "evaluator_context": "NO_PRIOR_CONVERSATION_PACKAGE_ONLY",
        "evaluator_id": "fresh-test-agent",
        "package_source_manifest": manifest,
        "questions": questions,
        "prohibited_promotions": promotions,
        "verdict": "PASS",
        "claim_ceiling": "Reader comprehension only; not scientific validation, host acceptance, or publication proof.",
    }


def validate_receipt(package: Path, receipt: dict):
    receipt_path = package / "reader_acceptance_test_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return MODULE.validate(receipt_path, package, require_head=True)


def test_valid_receipt_accepts_reachable_ancestor_and_dual_hashes(tmp_path):
    _, package, tested_commit = make_repository(tmp_path, descendant_commit=True)
    errors, summary = validate_receipt(package, make_receipt(package, tested_commit))
    assert errors == []
    assert summary["verdict"] == "PASS"
    assert summary["tested_commit_relation"] == "REACHABLE_ANCESTOR"
    assert summary["source_files"] == len(MODULE.REQUIRED_MANIFEST_PATHS)
    assert summary["current_hash_matches"] == len(MODULE.REQUIRED_MANIFEST_PATHS)
    assert summary["tested_commit_hash_matches"] == len(MODULE.REQUIRED_MANIFEST_PATHS)
    assert summary["questions"] == 6
    assert summary["prohibited_promotions"] == 8


def test_jsonschema_rejects_unknown_root_field(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["unexpected"] = "must fail"
    errors, _ = validate_receipt(package, receipt)
    assert any("schema violation" in error and "Additional properties" in error for error in errors)


def test_manifest_requires_every_core_document_and_registry(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["package_source_manifest"] = [
        row for row in receipt["package_source_manifest"] if row["path"] != "registries/artifact_registry.json"
    ]
    errors, _ = validate_receipt(package, receipt)
    assert any("lacks required core paths" in error and "artifact_registry.json" in error for error in errors)


def test_manifest_hash_must_match_current_bytes(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["package_source_manifest"][0]["sha256"] = "0" * 64
    errors, _ = validate_receipt(package, receipt)
    assert any("current package hash mismatch" in error for error in errors)


def test_manifest_hash_must_also_match_tested_commit_bytes(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    changed = package / "00_INDEX.md"
    changed.write_text("uncommitted package source mutation\n", encoding="utf-8")
    receipt = make_receipt(package, tested_commit)
    errors, _ = validate_receipt(package, receipt)
    assert any("tested commit hash mismatch: 00_INDEX.md" in error for error in errors)
    assert not any("current package hash mismatch: 00_INDEX.md" in error for error in errors)


def test_every_question_evidence_path_must_be_in_manifest(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["questions"][0]["evidence_paths"].append("unbound/source.md")
    errors, _ = validate_receipt(package, receipt)
    assert any("Q_PROJECT evidence path not in package_source_manifest" in error for error in errors)


def test_every_promotion_evidence_path_must_be_in_manifest(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["prohibited_promotions"][0]["evidence_paths"].append("unbound/source.md")
    errors, _ = validate_receipt(package, receipt)
    assert any("NO_CELLULAR_PROMOTION evidence path not in package_source_manifest" in error for error in errors)


@pytest.mark.parametrize("question_id", MODULE.QUESTIONS)
def test_generic_or_semantically_empty_question_answers_fail(tmp_path, question_id):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    row = next(item for item in receipt["questions"] if item["question_id"] == question_id)
    row["answer"] = (
        "This is a sufficiently long but generic response that tells the reader to see the documentation and "
        "claims that all required boundaries have been preserved without naming any concrete evidence."
    )
    errors, _ = validate_receipt(package, receipt)
    assert any(error.startswith(f"{question_id} missing required concepts") for error in errors)


@pytest.mark.parametrize("check_id", MODULE.PROMOTIONS)
def test_generic_promotion_explanations_fail(tmp_path, check_id):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    row = next(item for item in receipt["prohibited_promotions"] if item["check_id"] == check_id)
    row["explanation"] = (
        "This generic explanation merely says the answer preserves every required boundary without naming the "
        "specific scientific, software, publication, or machine constraint being tested."
    )
    errors, _ = validate_receipt(package, receipt)
    assert any(error.startswith(f"{check_id} missing required explanation concepts") for error in errors)


def test_conclusion_forbidden_accuracy_promotion_fails(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["questions"][1]["answer"] += " 88.2579% is the biological accuracy."
    errors, _ = validate_receipt(package, receipt)
    assert any("Q_CONCLUSION contains forbidden claim" in error for error in errors)


def test_wrong_zero_claims_fail(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["questions"][1]["answer"] = receipt["questions"][1]["answer"].replace(
        "Confirmed cellular subclones = 0", "Confirmed cellular subclones = 3"
    )
    errors, _ = validate_receipt(package, receipt)
    assert any("Q_CONCLUSION missing required concepts: zero confirmed cellular subclones" in error for error in errors)
    assert any("Q_CONCLUSION contains forbidden claim" in error for error in errors)


def test_unreachable_tested_commit_fails(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["tested_git_commit"] = "f" * 40
    errors, summary = validate_receipt(package, receipt)
    assert any("not a reachable Git object" in error for error in errors)
    assert summary["verdict"] == "FAIL"


def test_required_question_evidence_cannot_be_replaced_by_other_manifest_file(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["questions"][2]["evidence_paths"] = ["00_INDEX.md", "ai_context/CONTEXT.md"]
    errors, _ = validate_receipt(package, receipt)
    assert any("Q_FINALITY lacks required evidence paths" in error for error in errors)


def test_duplicate_generic_answers_fail(tmp_path):
    _, package, tested_commit = make_repository(tmp_path)
    receipt = make_receipt(package, tested_commit)
    receipt["questions"][1]["answer"] = receipt["questions"][0]["answer"]
    errors, _ = validate_receipt(package, receipt)
    assert any("duplicate generic answers" in error for error in errors)
