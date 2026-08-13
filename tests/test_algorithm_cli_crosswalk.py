from __future__ import annotations

import copy
import importlib.util
import json
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/handoff/validate_algorithm_cli_crosswalk.py"
BUILDER_SCRIPT = ROOT / "scripts/handoff/build_algorithm_cli_crosswalk.py"
SPEC = importlib.util.spec_from_file_location("validate_algorithm_cli_crosswalk", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)
BUILDER_SPEC = importlib.util.spec_from_file_location("build_algorithm_cli_crosswalk", BUILDER_SCRIPT)
assert BUILDER_SPEC is not None and BUILDER_SPEC.loader is not None
BUILDER = importlib.util.module_from_spec(BUILDER_SPEC)
BUILDER_SPEC.loader.exec_module(BUILDER)
CROSSWALK = VALIDATOR.DEFAULT_CROSSWALK


def validate_mutation(mutator):
    data = json.loads(CROSSWALK.read_text(encoding="utf-8"))
    mutator(data)
    with tempfile.TemporaryDirectory(prefix="algorithm-cli-crosswalk-test-") as temporary:
        path = Path(temporary) / "crosswalk.json"
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        return VALIDATOR.validate_crosswalk(path, ROOT)


def test_current_crosswalk_passes():
    errors, summary = VALIDATOR.validate_crosswalk(CROSSWALK, ROOT)
    assert errors == []
    assert summary["verdict"] == "PASS"
    assert summary["records"] == 35
    assert summary["source_rows"] == 35
    data = json.loads(CROSSWALK.read_text(encoding="utf-8"))
    assert data["gates"]["PUBLICATION_ASSET_READY"]["status"] == "PASS"
    assert data["gates"]["RELEASE_ASSET_READY"]["status"] == "PASS"
    assert "aggregate repository release approval" in data["gates"]["RELEASE_ASSET_READY"]["reason"]


def test_builder_is_deterministic_and_matches_tracked_candidate():
    current = json.loads(CROSSWALK.read_text(encoding="utf-8"))
    assert BUILDER.build() == current


def test_duplicate_or_missing_algorithm_row_fails_closed():
    def mutate(data):
        data["records"][1] = copy.deepcopy(data["records"][0])

    errors, summary = validate_mutation(mutate)
    assert summary["verdict"] == "FAIL"
    assert any("ALG-001..ALG-035 exactly once" in error for error in errors)


def test_source_row_hash_drift_fails_closed():
    def mutate(data):
        data["records"][0]["source_row_sha256"] = "0" * 64

    errors, _ = validate_mutation(mutate)
    assert any("ALG-001: source_row_sha256 mismatch" in error for error in errors)


def test_alg023_cannot_be_promoted_without_receipt():
    def mutate(data):
        row = next(item for item in data["records"] if item["algorithm_claim_id"] == "ALG-023")
        row["current_disposition"] = "CONFIRMED"
        row["evidence"]["status"] = "VERIFIED"
        row["release_gate"] = "SOURCE_READY"

    errors, _ = validate_mutation(mutate)
    assert any("ALG-023 invariant" in error for error in errors)


def test_alg025_representation_boundary_and_guard_are_required():
    def mutate(data):
        row = next(item for item in data["records"] if item["algorithm_claim_id"] == "ALG-025")
        row["bounded_current_claim"] = "The sidecar is the only interface."
        row["guard"]["status"] = "MISSING"
        row["guard"]["guard_ids"] = []

    errors, _ = validate_mutation(mutate)
    assert any("ALG-025 invariant" in error for error in errors)


def test_alg005_default_and_explicit_selection_boundary_is_required():
    def mutate(data):
        row = next(item for item in data["records"] if item["algorithm_claim_id"] == "ALG-005")
        row["bounded_current_claim"] = "Every request is forced to NHD."

    errors, _ = validate_mutation(mutate)
    assert any("ALG-005 invariant" in error for error in errors)


def test_hash_bound_source_snapshot_drift_fails_closed():
    def mutate(data):
        data["source_snapshot"][0]["sha256"] = "f" * 64

    errors, _ = validate_mutation(mutate)
    assert any("source snapshot sha256 mismatch" in error for error in errors)


def test_false_source_state_fails_closed():
    def mutate(data):
        data["source_state"] = (
            "WORKTREE_HASH_BOUND_PENDING_COMMIT"
            if data["source_state"] == "COMMIT_PINNED"
            else "COMMIT_PINNED"
        )

    errors, _ = validate_mutation(mutate)
    assert any("source_state" in error for error in errors)


def test_source_snapshot_tree_hash_drift_fails_closed():
    def mutate(data):
        data["source_snapshot_tree_sha256"] = "a" * 64

    errors, _ = validate_mutation(mutate)
    assert any("source_snapshot_tree_sha256 mismatch" in error for error in errors)


def test_related_claim_id_must_join_the_158_row_inventory():
    def mutate(data):
        data["records"][0]["related_claim_ids"] = ["C999"]

    errors, _ = validate_mutation(mutate)
    assert any("absent from the 158-row inventory" in error for error in errors)
