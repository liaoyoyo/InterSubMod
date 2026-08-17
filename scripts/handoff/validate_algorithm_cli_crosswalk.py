#!/usr/bin/env python3
"""Fail-closed validation for the 35-row algorithm/CLI claim crosswalk."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "docs/handoff/20260813_完整研究資料與軟體交接_01"
DEFAULT_CROSSWALK = PACKAGE / "evidence/algorithm_cli_claim_crosswalk.json"
DEFAULT_MATRIX = ROOT / "research/20260812_intersubmod_github_public_docs_full_validation/algorithm_cli_matrix.tsv"
EXPECTED_IDS = [f"ALG-{index:03d}" for index in range(1, 36)]
EXPECTED_COLUMNS = [
    "claim_id",
    "public_claim",
    "occurrence",
    "source_symbol",
    "runtime_or_test",
    "verdict",
    "minimum_rewrite",
    "caveat",
]
TOP_LEVEL_KEYS = {
    "schema_name",
    "schema_version",
    "crosswalk_id",
    "assessed_at",
    "assessed_git_commit",
    "source_state",
    "source_matrix",
    "claim_inventory",
    "guard_registry",
    "source_snapshot",
    "source_snapshot_tree_sha256",
    "disposition_semantics",
    "counts",
    "gates",
    "records",
}
RECORD_KEYS = {
    "algorithm_claim_id",
    "source_row_sha256",
    "audit_verdict",
    "original_public_claim",
    "bounded_current_claim",
    "current_disposition",
    "related_claim_ids",
    "source_correction",
    "guard",
    "evidence",
    "known_limits",
    "release_gate",
    "next_action",
}
DISPOSITIONS = {"CONFIRMED", "CONFIRMED_WITH_LIMITS", "UNVERIFIED"}
CORRECTION_STATUSES = {"NO_CHANGE_REQUIRED", "LOCAL_SOURCE_CORRECTED", "PARTIAL", "BLOCKED", "NOT_APPLICABLE"}
GUARD_STATUSES = {"PASS", "PARTIAL", "MISSING", "NOT_REQUIRED"}
EVIDENCE_STATUSES = {"VERIFIED", "VERIFIED_WITH_LIMITS", "SOURCE_INSPECTION_ONLY", "HISTORICAL_ONLY", "MISSING"}
AUDIT_VERDICTS = {"SUPPORTED", "SUPPORTED_WITH_LIMITS", "OVERSTATED", "STALE", "CONTRADICTED", "BRANCH_SCOPED", "UNVERIFIABLE"}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_matrix_row(row: dict[str, str]) -> bytes:
    return ("\t".join(row[column] for column in EXPECTED_COLUMNS) + "\n").encode("utf-8")


def load_matrix(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    errors: list[str] = []
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != EXPECTED_COLUMNS:
                errors.append(f"matrix columns differ: {reader.fieldnames!r}")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        return [], [f"cannot read source matrix: {error}"]
    ids = [row.get("claim_id", "") for row in rows]
    if ids != EXPECTED_IDS:
        errors.append(f"source matrix must contain ordered ALG-001..ALG-035 exactly once; got {ids!r}")
    return rows, errors


def check_exact_keys(value: Any, expected: set[str], label: str, errors: list[str]) -> bool:
    if not isinstance(value, dict):
        errors.append(f"{label} must be an object")
        return False
    keys = set(value)
    if keys != expected:
        errors.append(f"{label} keys differ; missing={sorted(expected - keys)!r} extra={sorted(keys - expected)!r}")
        return False
    return True


def load_guard_ids(path: Path) -> tuple[set[str], list[str]]:
    try:
        registry = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        return set(), [f"cannot read guard registry: {error}"]
    ids: set[str] = set()
    cross = registry.get("cross_document_guards", {})
    for group in ("document_rules", "context_rules"):
        for row in cross.get(group, []):
            if isinstance(row, dict) and isinstance(row.get("guard_id"), str):
                ids.add(row["guard_id"])
    for row in registry.get("claims", []):
        if isinstance(row, dict) and isinstance(row.get("claim_id"), str):
            ids.add(f"P0:{row['claim_id']}")
    return ids, []


def snapshot_tree_hash(bindings: list[dict[str, object]]) -> str:
    payload = json.dumps(bindings, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def commit_path_sha256(repo_root: Path, commit: str, relative: str) -> str | None:
    probe = subprocess.run(
        ["git", "-C", str(repo_root), "show", f"{commit}:{relative}"],
        capture_output=True,
    )
    if probe.returncode != 0:
        return None
    return hashlib.sha256(probe.stdout).hexdigest()


def validate_crosswalk(crosswalk_path: Path = DEFAULT_CROSSWALK, repo_root: Path = ROOT) -> tuple[list[str], dict[str, Any]]:
    errors: list[str] = []
    try:
        data = json.loads(crosswalk_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        return [f"cannot read crosswalk: {error}"], {}
    if not check_exact_keys(data, TOP_LEVEL_KEYS, "crosswalk", errors):
        return errors, {}
    if data.get("schema_name") != "intersubmod.algorithm_cli_claim_crosswalk":
        errors.append("schema_name must be intersubmod.algorithm_cli_claim_crosswalk")
    if data.get("schema_version") != "1.0.0":
        errors.append("schema_version must be 1.0.0")
    if data.get("crosswalk_id") != "intersubmod-algorithm-cli-claims-20260813":
        errors.append("unexpected crosswalk_id")
    commit = data.get("assessed_git_commit")
    if not isinstance(commit, str) or not COMMIT_RE.fullmatch(commit):
        errors.append("assessed_git_commit must be a full 40-character lowercase SHA")
    elif (repo_root / ".git").exists() or subprocess.run(
        ["git", "-C", str(repo_root), "rev-parse", "--git-dir"], capture_output=True, text=True
    ).returncode == 0:
        probe = subprocess.run(
            ["git", "-C", str(repo_root), "cat-file", "-e", f"{commit}^{{commit}}"],
            capture_output=True,
            text=True,
        )
        if probe.returncode != 0:
            errors.append(f"assessed_git_commit is not reachable in this repository: {commit}")

    matrix_binding = data.get("source_matrix")
    if not isinstance(matrix_binding, dict):
        errors.append("source_matrix must be an object")
        matrix_path = DEFAULT_MATRIX
    else:
        matrix_rel = matrix_binding.get("path")
        matrix_path = repo_root / matrix_rel if isinstance(matrix_rel, str) else DEFAULT_MATRIX
        if matrix_binding.get("rows") != 35:
            errors.append("source_matrix.rows must be 35")
        expected_hash = matrix_binding.get("sha256")
        if not isinstance(expected_hash, str) or not SHA256_RE.fullmatch(expected_hash):
            errors.append("source_matrix.sha256 is invalid")
        elif not matrix_path.is_file():
            errors.append(f"source matrix is missing: {matrix_path}")
        elif sha256_file(matrix_path) != expected_hash:
            errors.append("source matrix sha256 mismatch")
    matrix_rows, matrix_errors = load_matrix(matrix_path)
    errors.extend(matrix_errors)
    matrix_by_id = {row["claim_id"]: row for row in matrix_rows if row.get("claim_id")}

    claim_inventory_binding = data.get("claim_inventory")
    claim_inventory_ids: set[str] = set()
    if not isinstance(claim_inventory_binding, dict):
        errors.append("claim_inventory must be an object")
    else:
        inventory_rel = claim_inventory_binding.get("path")
        inventory_path = repo_root / inventory_rel if isinstance(inventory_rel, str) else repo_root / "__invalid_claim_inventory__"
        expected_hash = claim_inventory_binding.get("sha256")
        if claim_inventory_binding.get("rows") != 158:
            errors.append("claim_inventory.rows must be 158")
        if not inventory_path.is_file():
            errors.append(f"claim inventory is missing: {inventory_path}")
        elif not isinstance(expected_hash, str) or sha256_file(inventory_path) != expected_hash:
            errors.append("claim inventory sha256 mismatch")
        else:
            try:
                with inventory_path.open(encoding="utf-8", newline="") as handle:
                    inventory_rows = list(csv.DictReader(handle, delimiter="\t"))
                inventory_id_list = [row.get("claim_id", "") for row in inventory_rows]
                expected_inventory_ids = [f"C{index:03d}" for index in range(1, 159)]
                if inventory_id_list != expected_inventory_ids:
                    errors.append("claim inventory must contain ordered C001..C158 exactly once")
                claim_inventory_ids = set(inventory_id_list)
            except (OSError, UnicodeError, csv.Error) as error:
                errors.append(f"cannot parse claim inventory: {error}")

    guard_binding = data.get("guard_registry")
    guard_ids: set[str] = set()
    if not isinstance(guard_binding, dict):
        errors.append("guard_registry must be an object")
    else:
        guard_rel = guard_binding.get("path")
        guard_path = repo_root / guard_rel if isinstance(guard_rel, str) else repo_root / "__invalid_guard_registry__"
        expected_hash = guard_binding.get("sha256")
        if not guard_path.is_file():
            errors.append(f"guard registry is missing: {guard_path}")
        elif not isinstance(expected_hash, str) or sha256_file(guard_path) != expected_hash:
            errors.append("guard registry sha256 mismatch")
        else:
            guard_ids, guard_errors = load_guard_ids(guard_path)
            errors.extend(guard_errors)

    snapshot = data.get("source_snapshot")
    snapshot_paths: set[str] = set()
    if not isinstance(snapshot, list) or not snapshot:
        errors.append("source_snapshot must be a non-empty list")
    else:
        for index, binding in enumerate(snapshot):
            if not isinstance(binding, dict) or set(binding) != {"path", "sha256"}:
                errors.append(f"source_snapshot[{index}] must contain only path and sha256")
                continue
            relative = binding.get("path")
            expected_hash = binding.get("sha256")
            if not isinstance(relative, str) or not relative:
                errors.append(f"source_snapshot[{index}].path is invalid")
                continue
            if relative in snapshot_paths:
                errors.append(f"duplicate source_snapshot path: {relative}")
            snapshot_paths.add(relative)
            path = repo_root / relative
            if not path.is_file():
                errors.append(f"source snapshot path is missing: {relative}")
            elif not isinstance(expected_hash, str) or sha256_file(path) != expected_hash:
                errors.append(f"source snapshot sha256 mismatch: {relative}")
    snapshot_list = snapshot if isinstance(snapshot, list) else []
    if isinstance(snapshot, list):
        expected_tree_hash = snapshot_tree_hash(snapshot_list)
        if data.get("source_snapshot_tree_sha256") != expected_tree_hash:
            errors.append("source_snapshot_tree_sha256 mismatch")

    source_state = data.get("source_state")
    if source_state not in {"COMMIT_PINNED", "WORKTREE_HASH_BOUND_PENDING_COMMIT"}:
        errors.append("invalid source_state")
    elif isinstance(commit, str) and COMMIT_RE.fullmatch(commit):
        candidate_bindings = [data.get("source_matrix"), data.get("claim_inventory"), data.get("guard_registry"), *snapshot_list]
        commit_matches = []
        for binding in candidate_bindings:
            if not isinstance(binding, dict):
                continue
            relative = binding.get("path")
            expected_hash = binding.get("sha256")
            if isinstance(relative, str) and isinstance(expected_hash, str):
                commit_matches.append(commit_path_sha256(repo_root, commit, relative) == expected_hash)
        all_commit_pinned = bool(commit_matches) and all(commit_matches)
        if source_state == "COMMIT_PINNED" and not all_commit_pinned:
            errors.append("source_state COMMIT_PINNED but one or more source hashes differ from assessed_git_commit")
        if source_state == "WORKTREE_HASH_BOUND_PENDING_COMMIT" and all_commit_pinned:
            errors.append("source_state is pending even though all bound sources match assessed_git_commit")

    records = data.get("records")
    if not isinstance(records, list):
        errors.append("records must be a list")
        records = []
    record_ids = [row.get("algorithm_claim_id") for row in records if isinstance(row, dict)]
    if len(record_ids) != len(records):
        errors.append("records contains a non-object row")
    if record_ids != EXPECTED_IDS:
        errors.append(f"records must contain ordered ALG-001..ALG-035 exactly once; got {record_ids!r}")

    disposition_counts: Counter[str] = Counter()
    correction_counts: Counter[str] = Counter()
    guard_counts: Counter[str] = Counter()
    evidence_counts: Counter[str] = Counter()
    for index, record in enumerate(records):
        label = f"records[{index}]"
        if not check_exact_keys(record, RECORD_KEYS, label, errors):
            continue
        claim_id = record["algorithm_claim_id"]
        source_row = matrix_by_id.get(claim_id)
        if source_row is None:
            errors.append(f"{claim_id}: missing from source matrix")
        else:
            source_row_hash = hashlib.sha256(canonical_matrix_row(source_row)).hexdigest()
            if record["source_row_sha256"] != source_row_hash:
                errors.append(f"{claim_id}: source_row_sha256 mismatch")
            if record["audit_verdict"] != source_row["verdict"]:
                errors.append(f"{claim_id}: audit_verdict differs from source matrix")
            if record["original_public_claim"] != source_row["public_claim"]:
                errors.append(f"{claim_id}: original_public_claim differs from source matrix")

        disposition = record["current_disposition"]
        if disposition not in DISPOSITIONS:
            errors.append(f"{claim_id}: invalid current_disposition {disposition!r}")
        else:
            disposition_counts[disposition] += 1
        if record["audit_verdict"] not in AUDIT_VERDICTS:
            errors.append(f"{claim_id}: invalid audit_verdict")
        related = record["related_claim_ids"]
        if not isinstance(related, list) or any(not re.fullmatch(r"C[0-9]{3}", value) for value in related):
            errors.append(f"{claim_id}: related_claim_ids must contain only Cnnn identifiers")
        elif any(value not in claim_inventory_ids for value in related):
            errors.append(f"{claim_id}: related_claim_ids contains an ID absent from the 158-row inventory")

        correction = record["source_correction"]
        if not isinstance(correction, dict) or set(correction) != {"status", "target_paths", "note"}:
            errors.append(f"{claim_id}: invalid source_correction object")
        else:
            correction_status = correction["status"]
            if correction_status not in CORRECTION_STATUSES:
                errors.append(f"{claim_id}: invalid source correction status")
            else:
                correction_counts[correction_status] += 1
            targets = correction["target_paths"]
            if not isinstance(targets, list) or any(not isinstance(path, str) or not path for path in targets):
                errors.append(f"{claim_id}: source_correction.target_paths is invalid")
            else:
                for target in targets:
                    if target not in snapshot_paths:
                        errors.append(f"{claim_id}: correction target is not hash-bound in source_snapshot: {target}")

        guard = record["guard"]
        if not isinstance(guard, dict) or set(guard) != {"status", "guard_ids", "validator"}:
            errors.append(f"{claim_id}: invalid guard object")
        else:
            guard_status = guard["status"]
            if guard_status not in GUARD_STATUSES:
                errors.append(f"{claim_id}: invalid guard status")
            else:
                guard_counts[guard_status] += 1
            ids = guard["guard_ids"]
            if not isinstance(ids, list) or any(not isinstance(value, str) or not value for value in ids):
                errors.append(f"{claim_id}: guard_ids is invalid")
            elif guard_status == "PASS":
                if not ids:
                    errors.append(f"{claim_id}: PASS guard requires at least one guard_id")
                for guard_id in ids:
                    if guard_id not in guard_ids:
                        errors.append(f"{claim_id}: unknown PASS guard_id {guard_id}")
            elif guard_status in {"MISSING", "NOT_REQUIRED"} and ids:
                errors.append(f"{claim_id}: {guard_status} guard must not list guard_ids")

        evidence = record["evidence"]
        if not isinstance(evidence, dict) or set(evidence) != {"status", "references", "claim_ceiling"}:
            errors.append(f"{claim_id}: invalid evidence object")
        else:
            evidence_status = evidence["status"]
            if evidence_status not in EVIDENCE_STATUSES:
                errors.append(f"{claim_id}: invalid evidence status")
            else:
                evidence_counts[evidence_status] += 1
            references = evidence["references"]
            if not isinstance(references, list) or not references or any(not isinstance(value, str) or not value for value in references):
                errors.append(f"{claim_id}: evidence.references must be a non-empty string list")
            if not isinstance(evidence["claim_ceiling"], str) or not evidence["claim_ceiling"].strip():
                errors.append(f"{claim_id}: evidence.claim_ceiling must be non-empty")

        limits = record["known_limits"]
        if not isinstance(limits, list) or any(not isinstance(value, str) or not value for value in limits):
            errors.append(f"{claim_id}: known_limits is invalid")
        if disposition == "CONFIRMED" and evidence.get("status") != "VERIFIED":
            errors.append(f"{claim_id}: CONFIRMED requires VERIFIED evidence")
        if disposition == "CONFIRMED_WITH_LIMITS" and not limits:
            errors.append(f"{claim_id}: CONFIRMED_WITH_LIMITS requires known_limits")
        if disposition == "UNVERIFIED":
            if evidence.get("status") not in {"MISSING", "HISTORICAL_ONLY"}:
                errors.append(f"{claim_id}: UNVERIFIED requires MISSING or HISTORICAL_ONLY evidence")
            if record["release_gate"] != "BLOCKED":
                errors.append(f"{claim_id}: UNVERIFIED must remain release_gate BLOCKED")

    by_id = {record.get("algorithm_claim_id"): record for record in records if isinstance(record, dict)}
    alg023 = by_id.get("ALG-023", {})
    if alg023.get("current_disposition") != "UNVERIFIED":
        errors.append("ALG-023 invariant: 334/111 claim must remain UNVERIFIED")
    if alg023.get("release_gate") != "BLOCKED" or alg023.get("evidence", {}).get("status") != "MISSING":
        errors.append("ALG-023 invariant: missing replay receipt must remain MISSING/BLOCKED")

    alg025 = by_id.get("ALG-025", {})
    bounded_025 = alg025.get("bounded_current_claim", "")
    required_025 = ("InterSubMod reads HP/PS from BAM", "LongLineage/exact-PS use sidecars", "not a direct engine-to-engine runtime interface")
    if alg025.get("current_disposition") != "CONFIRMED_WITH_LIMITS" or any(value not in bounded_025 for value in required_025):
        errors.append("ALG-025 invariant: bounded claim must separate BAM tags, sidecars and the absent direct runtime bridge")
    if alg025.get("guard", {}).get("status") != "PASS" or "G_HPPS_REPRESENTATION_BOUNDARY" not in alg025.get("guard", {}).get("guard_ids", []):
        errors.append("ALG-025 invariant: G_HPPS_REPRESENTATION_BOUNDARY must PASS")

    alg005 = by_id.get("ALG-005", {})
    bounded_005 = alg005.get("bounded_current_claim", "")
    required_005 = ("no-argument default is NHD", "explicit NHD/BERNOULLI selections are honored")
    if alg005.get("current_disposition") != "CONFIRMED_WITH_LIMITS" or any(value not in bounded_005 for value in required_005):
        errors.append("ALG-005 invariant: bounded claim must distinguish default from explicit selections")
    if alg005.get("guard", {}).get("status") != "PASS" or "G_NHD_CLI_DEFAULT_SEMANTICS" not in alg005.get("guard", {}).get("guard_ids", []):
        errors.append("ALG-005 invariant: G_NHD_CLI_DEFAULT_SEMANTICS must PASS")

    expected_counts = {
        "records": len(records),
        "related_claim_links": sum(len(record.get("related_claim_ids", [])) for record in records if isinstance(record, dict)),
        "records_with_related_claim_ids": sum(bool(record.get("related_claim_ids")) for record in records if isinstance(record, dict)),
        "records_without_related_claim_ids": sum(not record.get("related_claim_ids") for record in records if isinstance(record, dict)),
        "by_current_disposition": dict(sorted(disposition_counts.items())),
        "by_source_correction": dict(sorted(correction_counts.items())),
        "by_guard_status": dict(sorted(guard_counts.items())),
        "by_evidence_status": dict(sorted(evidence_counts.items())),
    }
    if data.get("counts") != expected_counts:
        errors.append(f"counts differ from records: expected {expected_counts!r}")
    gates = data.get("gates")
    if not isinstance(gates, dict):
        errors.append("gates must be an object")
    else:
        if gates.get("CROSSWALK_COMPLETENESS", {}).get("status") != "PASS":
            errors.append("CROSSWALK_COMPLETENESS must be PASS")
        for gate in ("PUBLICATION_ASSET_READY", "RELEASE_ASSET_READY"):
            if gates.get(gate, {}).get("status") != "PASS":
                errors.append(f"{gate} must be PASS for the bounded crosswalk asset")

    summary = {
        "crosswalk": str(crosswalk_path),
        "source_matrix": str(matrix_path),
        "records": len(records),
        "source_rows": len(matrix_rows),
        "by_current_disposition": dict(sorted(disposition_counts.items())),
        "by_source_correction": dict(sorted(correction_counts.items())),
        "by_guard_status": dict(sorted(guard_counts.items())),
        "by_evidence_status": dict(sorted(evidence_counts.items())),
        "errors": len(errors),
        "verdict": "PASS" if not errors else "FAIL",
    }
    return errors, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosswalk", type=Path, default=DEFAULT_CROSSWALK)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    args = parser.parse_args()
    errors, summary = validate_crosswalk(args.crosswalk.resolve(), args.repo_root.resolve())
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
