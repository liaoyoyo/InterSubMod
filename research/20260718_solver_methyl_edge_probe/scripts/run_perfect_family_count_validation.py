#!/usr/bin/env python3
"""Build a non-overwriting receipt for the exact ``h*=m-d`` family counter."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import pathlib
import re
import sys
from datetime import datetime, timezone
from typing import Any, Mapping, Sequence

from perfect_family_count import (
    ABSTAIN_STATUS,
    EXACT_STATUS,
    EXACT_STATUSES,
    assert_finite_elapsed,
    count_manifest_case,
)


SCHEMA = "intersubmod.perfect_family_count_validation.v1"
MANIFEST_SCHEMA = "intersubmod.solver_stress_panel_manifest"
POINTER_SCHEMA = "intersubmod.solver_stress_panel.authority_pointer.v1"
REQUIRED_AUTHORITY_STATUS = "AUTHORITATIVE_R3"
EVIDENCE_SCHEMA = "intersubmod.solver_stress_panel.worker.v1"
EVIDENCE_BACKEND = "optimized"
EVIDENCE_OBJECTIVE_STATUS = "OPTIMAL_VALUE_CERTIFIED"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
HARD_CASE_IDS = (
    "k10_m10_a0944cdd1ac8fa9d",
    "k10_m9_e5b33e46c7c23c0f",
    "k11_m11_09b1da787e58efed",
    "k14_m10_cecad6897b192d47",
)
EXPECTED_HARD_COUNTS = {
    "k10_m10_a0944cdd1ac8fa9d": 104_640,
    "k11_m11_09b1da787e58efed": 122_281_152,
    "k14_m10_cecad6897b192d47": 27_360,
}
EXPECTED_RECURRENCE_CONTROL = "k6_m6_bd8aa7beafa3f719"


class ValidationError(RuntimeError):
    """Raised when a frozen input or exact oracle does not match."""


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _semantic_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return _sha256_bytes(payload)


def _manifest_content_sha256(document: Mapping[str, Any]) -> str:
    payload = copy.deepcopy(dict(document))
    payload.setdefault("integrity", {}).pop("manifest_content_sha256", None)
    return _semantic_sha256(payload)


def _require_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise ValidationError(f"{label} is not a lowercase SHA-256 digest")
    return value


def _verify_byte_sidecar(path: pathlib.Path) -> str:
    sidecar = path.with_name(path.name + ".sha256")
    if not sidecar.is_file():
        raise ValidationError(f"SHA-256 sidecar is missing: {sidecar}")
    fields = sidecar.read_text(encoding="utf-8").strip().split()
    if len(fields) != 2 or fields[1] != path.name:
        raise ValidationError(f"malformed SHA-256 sidecar: {sidecar}")
    expected = _require_sha256(fields[0], f"sidecar digest for {path}")
    observed = _sha256_file(path)
    if expected != observed:
        raise ValidationError(
            f"byte digest mismatch for {path}: {observed} != {expected}"
        )
    return observed


def _validate_case_record(
    case: Mapping[str, Any],
    *,
    label: str,
) -> tuple[str, str, int, int, Mapping[str, Any]]:
    if not isinstance(case, Mapping):
        raise ValidationError(f"{label} is not an object")
    case_id = case.get("case_id")
    if not isinstance(case_id, str) or not case_id:
        raise ValidationError(f"{label}.case_id must be a non-empty string")
    key_sha = _require_sha256(
        case.get("structural_key_sha256"),
        f"{label}.structural_key_sha256",
    )
    structural = case.get("structural_input")
    if not isinstance(structural, Mapping):
        raise ValidationError(f"{label}.structural_input is not an object")
    k = structural.get("k")
    if type(k) is not int or k < 1:
        raise ValidationError(f"{label}.structural_input.k must be a positive int")
    mask = structural.get("structural_alt_universe_mask")
    if type(mask) is not int or mask < 0 or mask >= (1 << k):
        raise ValidationError(
            f"{label}.structural_input.structural_alt_universe_mask "
            "must be an in-range int"
        )
    for field, allow_x in (("full_patterns", False), ("partial_patterns", True)):
        patterns = structural.get(field)
        if not isinstance(patterns, list):
            raise ValidationError(f"{label}.structural_input.{field} must be a list")
        allowed = {"A", "R", "X"} if allow_x else {"A", "R"}
        for index, pattern in enumerate(patterns):
            if not isinstance(pattern, str):
                raise ValidationError(
                    f"{label}.structural_input.{field}[{index}] must be a string"
                )
            if len(pattern) != k or not set(pattern) <= allowed:
                raise ValidationError(
                    f"{label}.structural_input.{field}[{index}] is invalid"
                )
    expected_key = _semantic_sha256(structural)
    if key_sha != expected_key:
        raise ValidationError(
            f"{label} structural key does not hash structural_input"
        )
    observed_alt_mask = 0
    for pattern in structural["full_patterns"] + structural["partial_patterns"]:
        for bit, symbol in enumerate(pattern):
            if symbol == "A":
                observed_alt_mask |= 1 << bit
    if mask != observed_alt_mask:
        raise ValidationError(
            f"{label} structural ALT universe differs from observed ALT union"
        )
    return case_id, key_sha, k, bin(mask).count("1"), structural


def _index_manifest_cases(document: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    raw_cases = document.get("cases")
    if not isinstance(raw_cases, list):
        raise ValidationError("manifest cases must be a list")
    cases: dict[str, Mapping[str, Any]] = {}
    for index, case in enumerate(raw_cases):
        case_id, _, _, _, _ = _validate_case_record(
            case,
            label=f"manifest.cases[{index}]",
        )
        if case_id in cases:
            raise ValidationError(f"duplicate manifest case_id: {case_id}")
        cases[case_id] = case
    raw_units = document.get("units")
    if not isinstance(raw_units, list):
        raise ValidationError("manifest units must be a list")
    unit_case_ids: set[str] = set()
    for index, unit in enumerate(raw_units):
        if not isinstance(unit, Mapping):
            raise ValidationError(f"manifest.units[{index}] is not an object")
        case_id = unit.get("case_id")
        if not isinstance(case_id, str):
            raise ValidationError(
                f"manifest.units[{index}].case_id must be a string"
            )
        if case_id not in cases:
            raise ValidationError(
                f"manifest.units[{index}] references unknown case_id {case_id!r}"
            )
        unit_case_ids.add(case_id)
    missing_units = sorted(set(cases) - unit_case_ids)
    if missing_units:
        raise ValidationError(
            "manifest contains case(s) with no unit: " + ", ".join(missing_units)
        )
    return cases


def _load_verified_manifest(
    path: pathlib.Path,
    *,
    required_authority_status: str | None,
) -> dict[str, Any]:
    if not path.is_file():
        raise ValidationError(f"stress-panel manifest is missing: {path}")
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema_name") != MANIFEST_SCHEMA:
        raise ValidationError("unexpected stress-panel manifest schema")
    expected = document.get("integrity", {}).get("manifest_content_sha256")
    _require_sha256(expected, "manifest semantic digest")
    observed = _manifest_content_sha256(document)
    if expected != observed:
        raise ValidationError("manifest semantic digest mismatch")
    checks = document.get("checks")
    if (
        not isinstance(checks, Mapping)
        or not checks
        or not all(value is True for value in checks.values())
    ):
        raise ValidationError("manifest frozen checks are not all true")
    _verify_byte_sidecar(path)
    if required_authority_status is not None:
        authority_status = document.get("authority", {}).get("status")
        if authority_status != required_authority_status:
            raise ValidationError(
                "manifest authority status mismatch: "
                f"{authority_status!r} != {required_authority_status!r}"
            )
    _index_manifest_cases(document)
    return document


def _load_frozen_manifest(path: pathlib.Path) -> dict[str, Any]:
    """Load a direct manifest only when it is explicitly authoritative R3."""
    return _load_verified_manifest(
        path.resolve(),
        required_authority_status=REQUIRED_AUTHORITY_STATUS,
    )


def _resolve_authority_pointer(
    pointer_path: pathlib.Path,
) -> tuple[pathlib.Path, dict[str, Any], dict[str, Any]]:
    pointer_path = pointer_path.resolve()
    if not pointer_path.is_file():
        raise ValidationError(f"authority pointer is missing: {pointer_path}")
    pointer_file_sha = _verify_byte_sidecar(pointer_path)
    pointer = json.loads(pointer_path.read_text(encoding="utf-8"))
    if pointer.get("schema") != POINTER_SCHEMA:
        raise ValidationError("unexpected authority-pointer schema")
    if pointer.get("status") != REQUIRED_AUTHORITY_STATUS:
        raise ValidationError(
            "authority pointer status is not "
            f"{REQUIRED_AUTHORITY_STATUS}: {pointer.get('status')!r}"
        )
    relative_value = pointer.get("authoritative_manifest")
    if not isinstance(relative_value, str) or not relative_value:
        raise ValidationError("authority pointer manifest path is invalid")
    relative_path = pathlib.Path(relative_value)
    if relative_path.is_absolute():
        raise ValidationError("authority pointer manifest path must be relative")
    pointer_root = pointer_path.parent.resolve()
    manifest_path = (pointer_root / relative_path).resolve()
    try:
        manifest_path.relative_to(pointer_root)
    except ValueError as error:
        raise ValidationError(
            "authority pointer manifest path escapes pointer directory"
        ) from error
    expected_file_sha = _require_sha256(
        pointer.get("authoritative_manifest_file_sha256"),
        "pointer authoritative_manifest_file_sha256",
    )
    observed_file_sha = _sha256_file(manifest_path)
    if observed_file_sha != expected_file_sha:
        raise ValidationError("authority pointer manifest byte digest mismatch")
    manifest = _load_verified_manifest(
        manifest_path,
        required_authority_status=REQUIRED_AUTHORITY_STATUS,
    )
    expected_content_sha = _require_sha256(
        pointer.get("authoritative_manifest_content_sha256"),
        "pointer authoritative_manifest_content_sha256",
    )
    if manifest["integrity"]["manifest_content_sha256"] != expected_content_sha:
        raise ValidationError("authority pointer manifest content digest mismatch")
    return manifest_path, manifest, {
        "path": str(pointer_path),
        "sha256": pointer_file_sha,
        "schema": POINTER_SCHEMA,
        "status": REQUIRED_AUTHORITY_STATUS,
    }


def _resolve_authoritative_manifest(
    *,
    authority_pointer_path: pathlib.Path | None,
    manifest_path: pathlib.Path | None,
) -> tuple[pathlib.Path, dict[str, Any], dict[str, Any]]:
    if (authority_pointer_path is None) == (manifest_path is None):
        raise ValidationError(
            "specify exactly one of authority_pointer_path or manifest_path"
        )
    if authority_pointer_path is not None:
        return _resolve_authority_pointer(authority_pointer_path)
    assert manifest_path is not None
    resolved = manifest_path.resolve()
    manifest = _load_frozen_manifest(resolved)
    return resolved, manifest, {
        "path": None,
        "sha256": None,
        "schema": None,
        "status": REQUIRED_AUTHORITY_STATUS,
    }


def _find_evidence_source_manifest(
    *,
    evidence_dir: pathlib.Path,
    expected_content_sha256: str,
    authoritative_manifest_path: pathlib.Path,
) -> tuple[pathlib.Path, dict[str, Any]]:
    """Resolve a historical evidence manifest by hash and verify its cases."""
    candidates: set[pathlib.Path] = {authoritative_manifest_path.resolve()}
    evidence_root = evidence_dir.resolve()
    for directory in (evidence_root, *evidence_root.parents):
        candidate = directory / "manifest.json"
        if candidate.is_file():
            candidates.add(candidate.resolve())

    matches: list[tuple[pathlib.Path, dict[str, Any]]] = []
    for candidate in sorted(candidates):
        try:
            raw = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if raw.get("schema_name") != MANIFEST_SCHEMA:
            continue
        if _manifest_content_sha256(raw) != expected_content_sha256:
            continue
        verified = _load_verified_manifest(
            candidate,
            required_authority_status=None,
        )
        matches.append((candidate, verified))
    if len(matches) != 1:
        raise ValidationError(
            "historical evidence source manifest is "
            f"{'missing' if not matches else 'ambiguous'} for content digest "
            f"{expected_content_sha256}"
        )
    return matches[0]


def _evidence_record(
    evidence_dir: pathlib.Path,
    case_id: str,
    *,
    expected_case: Mapping[str, Any],
    authoritative_manifest: Mapping[str, Any],
    authoritative_manifest_path: pathlib.Path,
) -> dict[str, Any] | None:
    paths = sorted(evidence_dir.glob(f"{case_id}.r*.optimized.json"))
    if not paths:
        return None
    if len(paths) != 1:
        raise ValidationError(f"ambiguous optimized evidence for {case_id}")
    path = paths[0]
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("case_id") != case_id:
        raise ValidationError(f"evidence case mismatch for {case_id}")
    (
        _,
        expected_structural_key,
        expected_k,
        expected_effective_m,
        expected_structural_input,
    ) = _validate_case_record(expected_case, label=f"authoritative case {case_id}")
    if document.get("structural_key_sha256") != expected_structural_key:
        raise ValidationError(f"evidence structural key mismatch for {case_id}")
    if type(document.get("k")) is not int or document["k"] != expected_k:
        raise ValidationError(f"evidence k mismatch for {case_id}")
    if (
        type(document.get("effective_m")) is not int
        or document["effective_m"] != expected_effective_m
    ):
        raise ValidationError(f"evidence effective_m mismatch for {case_id}")
    if document.get("backend") != EVIDENCE_BACKEND:
        raise ValidationError(f"evidence backend is not optimized for {case_id}")
    if document.get("schema") != EVIDENCE_SCHEMA:
        raise ValidationError(f"evidence schema mismatch for {case_id}")
    if document.get("objective_status") != EVIDENCE_OBJECTIVE_STATUS:
        raise ValidationError(f"evidence objective status mismatch for {case_id}")
    if document.get("objective_certified") is not True:
        raise ValidationError(f"evidence objective is not certified for {case_id}")
    if (
        type(document.get("objective")) is not int
        or document["objective"] < 0
    ):
        raise ValidationError(f"evidence objective is not a nonnegative int for {case_id}")
    if type(document.get("complete")) is not bool:
        raise ValidationError(f"evidence complete flag is not bool for {case_id}")
    if not isinstance(document.get("status"), str) or not document["status"]:
        raise ValidationError(f"evidence status is invalid for {case_id}")

    evidence_manifest_sha = _require_sha256(
        document.get("manifest_content_sha256"),
        f"evidence manifest content digest for {case_id}",
    )
    authoritative_manifest_sha = authoritative_manifest["integrity"][
        "manifest_content_sha256"
    ]
    if evidence_manifest_sha == authoritative_manifest_sha:
        source_manifest_path = authoritative_manifest_path
        source_manifest = authoritative_manifest
        equivalence = "CURRENT_AUTHORITATIVE_MANIFEST"
    else:
        source_manifest_path, source_manifest = _find_evidence_source_manifest(
            evidence_dir=evidence_dir,
            expected_content_sha256=evidence_manifest_sha,
            authoritative_manifest_path=authoritative_manifest_path,
        )
        equivalence = "VERIFIED_STRUCTURAL_EQUIVALENCE"
    source_cases = _index_manifest_cases(source_manifest)
    source_case = source_cases.get(case_id)
    if source_case is None:
        raise ValidationError(
            f"evidence source manifest lacks case {case_id}"
        )
    (
        _,
        source_structural_key,
        source_k,
        source_effective_m,
        source_structural_input,
    ) = _validate_case_record(
        source_case,
        label=f"evidence source case {case_id}",
    )
    if (
        source_structural_key != expected_structural_key
        or source_k != expected_k
        or source_effective_m != expected_effective_m
        or dict(source_structural_input) != dict(expected_structural_input)
    ):
        raise ValidationError(
            f"evidence source case is not structurally equivalent for {case_id}"
        )
    return {
        "path": str(path.resolve()),
        "sha256": _sha256_file(path),
        "status": document.get("status"),
        "schema": document["schema"],
        "backend": document["backend"],
        "structural_key_sha256": document["structural_key_sha256"],
        "k": document["k"],
        "effective_m": document["effective_m"],
        "complete": document.get("complete") is True,
        "objective_status": document["objective_status"],
        "objective_certified": True,
        "objective": document.get("objective"),
        "optimal_family_count": document.get("optimal_family_count"),
        "source_manifest": {
            "path": str(source_manifest_path.resolve()),
            "sha256": _sha256_file(source_manifest_path),
            "semantic_sha256": evidence_manifest_sha,
            "binding": equivalence,
        },
        "acceptance_role": (
            "OBJECTIVE_DIAGNOSTIC_ONLY"
            if not document.get("complete")
            else "COMPLETE_CONTROL_ORACLE"
        ),
    }


def build_receipt(
    *,
    authority_pointer_path: pathlib.Path | None = None,
    manifest_path: pathlib.Path | None = None,
    evidence_dir: pathlib.Path,
    source_path: pathlib.Path,
    runner_path: pathlib.Path,
) -> dict[str, Any]:
    manifest_path, manifest, authority = _resolve_authoritative_manifest(
        authority_pointer_path=authority_pointer_path,
        manifest_path=manifest_path,
    )
    cases = _index_manifest_cases(manifest)
    if len(cases) != 46:
        raise ValidationError(f"expected 46 structural keys, observed {len(cases)}")

    results: list[dict[str, Any]] = []
    indexed_results: dict[str, dict[str, Any]] = {}
    for case_id in sorted(cases):
        result = count_manifest_case(cases[case_id])
        assert_finite_elapsed(result)
        row = {
            "case_id": case_id,
            "structural_key_sha256": cases[case_id]["structural_key_sha256"],
            **result.to_dict(),
        }
        results.append(row)
        indexed_results[case_id] = row

    units_by_case: dict[str, list[Mapping[str, Any]]] = {}
    for unit in manifest["units"]:
        units_by_case.setdefault(unit["case_id"], []).append(unit)
    complete_case_ids = sorted(
        case_id
        for case_id, units in units_by_case.items()
        if any(unit["source_classification"] == "complete" for unit in units)
    )
    if len(complete_case_ids) != 8:
        raise ValidationError(
            f"expected 8 complete control keys, observed {len(complete_case_ids)}"
        )

    controls: list[dict[str, Any]] = []
    exact_matches = 0
    abstentions = 0
    for case_id in complete_case_ids:
        expected_counts = {
            int(unit["source_candidate_vertex_sets"])
            for unit in units_by_case[case_id]
            if unit["source_classification"] == "complete"
        }
        if len(expected_counts) != 1:
            raise ValidationError(f"non-unique complete count for {case_id}")
        expected_count = next(iter(expected_counts))
        observed = indexed_results[case_id]
        if observed["status"] in EXACT_STATUSES:
            passed = observed["perfect_family_count"] == expected_count
            outcome = "EXACT_COUNT_MATCH" if passed else "EXACT_COUNT_MISMATCH"
            exact_matches += int(passed)
        else:
            passed = (
                case_id == EXPECTED_RECURRENCE_CONTROL
                and observed["status"] == ABSTAIN_STATUS
            )
            outcome = (
                "EXPECTED_RECURRENCE_ABSTENTION"
                if passed
                else "UNEXPECTED_ABSTENTION"
            )
            abstentions += int(passed)
        controls.append(
            {
                "case_id": case_id,
                "expected_source_complete_count": expected_count,
                "observed_perfect_count": observed["perfect_family_count"],
                "counter_status": observed["status"],
                "outcome": outcome,
                "pass": passed,
                "optimized_evidence": _evidence_record(
                    evidence_dir,
                    case_id,
                    expected_case=cases[case_id],
                    authoritative_manifest=manifest,
                    authoritative_manifest_path=manifest_path,
                ),
            }
        )
    if exact_matches != 7 or abstentions != 1:
        raise ValidationError(
            f"control oracle failed: exact={exact_matches}, abstain={abstentions}"
        )
    if not all(control["pass"] for control in controls):
        raise ValidationError("one or more complete control oracles failed")

    hard_cases: list[dict[str, Any]] = []
    for case_id in HARD_CASE_IDS:
        result = indexed_results[case_id]
        expected_count = EXPECTED_HARD_COUNTS.get(case_id)
        if expected_count is None:
            passed = (
                result["status"] == ABSTAIN_STATUS
                and result["perfect_family_count"] == 0
            )
        else:
            passed = (
                result["status"] == EXACT_STATUS
                and result["perfect_family_count"] == expected_count
            )
        evidence = _evidence_record(
            evidence_dir,
            case_id,
            expected_case=cases[case_id],
            authoritative_manifest=manifest,
            authoritative_manifest_path=manifest_path,
        )
        if evidence is None or not evidence["objective_certified"]:
            raise ValidationError(f"missing objective evidence for {case_id}")
        recurrence_relation = "NOT_APPLICABLE"
        if result["status"] == ABSTAIN_STATUS:
            if evidence["objective"] != result["effective_m"] + 1:
                raise ValidationError(
                    f"recurrence control objective is not m+1 for {case_id}"
                )
            recurrence_relation = "CERTIFIED_H_EQ_M_PLUS_1"
        elif evidence["objective"] != result["objective"]:
            raise ValidationError(
                f"objective evidence differs from exact counter for {case_id}"
            )
        hard_cases.append(
            {
                "case_id": case_id,
                "counter_status": result["status"],
                "effective_m": result["effective_m"],
                "perfect_family_count": result["perfect_family_count"],
                "expected_perfect_family_count": expected_count,
                "objective_evidence": evidence,
                "recurrence_relation": recurrence_relation,
                "pass": passed,
            }
        )
    if not all(row["pass"] for row in hard_cases):
        raise ValidationError("hard-key oracle failed")

    payload: dict[str, Any] = {
        "schema_name": SCHEMA,
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": (
            "isolated exact perfect-family counting; no canonical routing, "
            "no VAF ranking, no recurrence-family enumeration"
        ),
        "verdict": "PASS_EXACT_PERFECT_COUNT_WITH_RECURRENCE_ABSTENTION",
        "ranking_allowed": False,
        "canonical_promotion_allowed": False,
        "method_contract": {
            "mandatory_state_handling": (
                "full-read states are exact state-presence constraints"
            ),
            "exact_gate": (
                "perfect_family_count > 0 proves h*=effective_m-d, where d "
                "is the distinct mandatory non-root state count; hard keys "
                "are root-only so d=0 and h*=effective_m"
            ),
            "zero_count_policy": ABSTAIN_STATUS,
            "time_complexity": "O(3^m poly(m))",
            "space_complexity": "O(2^m + m^2)",
            "family_materialized": False,
            "next_representation": (
                "DP traceback DAG or ZDD required before family-aware ranking"
            ),
        },
        "bindings": {
            "manifest": {
                "path": str(manifest_path.resolve()),
                "sha256": _sha256_file(manifest_path),
                "semantic_sha256": manifest["integrity"][
                    "manifest_content_sha256"
                ],
                "authority_status": REQUIRED_AUTHORITY_STATUS,
            },
            "authority_pointer": authority,
            "counter_source": {
                "path": str(source_path.resolve()),
                "sha256": _sha256_file(source_path),
            },
            "validation_runner": {
                "path": str(runner_path.resolve()),
                "sha256": _sha256_file(runner_path),
            },
            "evidence_dir": str(evidence_dir.resolve()),
        },
        "validation_summary": {
            "structural_keys_evaluated": len(results),
            "complete_control_keys": len(controls),
            "exact_count_matches": exact_matches,
            "expected_recurrence_abstentions": abstentions,
            "failures": 0,
        },
        "hard_cases": hard_cases,
        "complete_controls": controls,
        "all_case_results": results,
        "immutability": {
            "write_policy": "exclusive_create",
            "receipt_mode_after_write": "0444",
            "sidecar_mode_after_write": "0444",
        },
    }
    content_payload = copy.deepcopy(payload)
    payload["integrity"] = {
        "receipt_content_sha256": _semantic_sha256(content_payload)
    }
    return payload


def write_receipt(path: pathlib.Path, receipt: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = (
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    with path.open("xb") as handle:
        handle.write(encoded)
    sidecar = path.with_name(path.name + ".sha256")
    sidecar_line = f"{_sha256_bytes(encoded)}  {path.name}\n".encode("ascii")
    with sidecar.open("xb") as handle:
        handle.write(sidecar_line)
    os.chmod(path, 0o444)
    os.chmod(sidecar, 0o444)


def main() -> int:
    parser = argparse.ArgumentParser()
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument(
        "--authority-pointer",
        type=pathlib.Path,
        help="preferred: cryptographically resolve the AUTHORITATIVE_R3 manifest",
    )
    source_group.add_argument(
        "--manifest",
        type=pathlib.Path,
        help="direct manifest; authority.status must be AUTHORITATIVE_R3",
    )
    parser.add_argument("--evidence-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    source_path = pathlib.Path(__file__).with_name("perfect_family_count.py")
    runner_path = pathlib.Path(__file__)
    receipt = build_receipt(
        authority_pointer_path=args.authority_pointer,
        manifest_path=args.manifest,
        evidence_dir=args.evidence_dir,
        source_path=source_path,
        runner_path=runner_path,
    )
    write_receipt(args.output, receipt)
    print(
        json.dumps(
            {
                "status": receipt["verdict"],
                "output": str(args.output.resolve()),
                "structural_keys": receipt["validation_summary"][
                    "structural_keys_evaluated"
                ],
                "controls": receipt["validation_summary"][
                    "complete_control_keys"
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValidationError, ValueError, OSError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
