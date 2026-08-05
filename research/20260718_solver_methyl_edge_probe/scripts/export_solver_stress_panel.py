#!/usr/bin/env python3
"""Export a write-once structural stress panel from the frozen M2 diagnostics.

The three input gzip streams are known to be truncated.  This exporter consumes
only complete TSV records, verifies every source binding from the prior census,
and fails closed if any selected unit cannot be reconstructed losslessly.
It does not claim that the truncated streams represent a complete chromosome.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import pathlib
import sys
from datetime import datetime, timezone
from typing import Any, Dict, Mapping, Sequence, Tuple


REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]
PROBE_ROOT = REPO_ROOT / "research/20260718_solver_methyl_edge_probe"
DEFAULT_CENSUS = PROBE_ROOT / "results/m2_solver_route_census_receipt.json"
DEFAULT_OUTPUT_DIR = PROBE_ROOT / "results/solver_stress_panel_v1"
DEFAULT_CENSUS_SOURCE = PROBE_ROOT / "scripts/census_m2_solver_routes.py"
DEFAULT_PROBE_SOURCE = PROBE_ROOT / "scripts/solver_probe.py"
DEFAULT_OPTIMIZED_SOURCE = PROBE_ROOT / "scripts/optimized_hypercube_backend.py"
DEFAULT_RUNNER_SOURCE = PROBE_ROOT / "scripts/run_solver_stress_panel.py"

PANEL_CONFIG = {
    "schema": "intersubmod.solver_stress_panel.selection.v1",
    "dataset": "HCC1395_DORADO",
    "chrom": "chr6",
    "universe_mode": "observed_alt",
    "selection": {
        "tail": "candidate_generation_invoked=true and likelihood_fit_invoked=false",
        "control": "candidate_generation_invoked=true and likelihood_fit_invoked=true and candidate_vertex_sets>=100",
    },
    "panels": {
        "primary": {
            "minread": 1,
            "expected_units": 44,
            "expected_distinct_structural_keys": 25,
        },
        "sensitivity": {
            "minread": 2,
            "expected_units": 36,
            "expected_distinct_structural_keys": 21,
        },
    },
}


class PanelExportError(RuntimeError):
    """Raised when the bounded panel cannot be exported losslessly."""


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def semantic_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load_module(path: pathlib.Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise PanelExportError(f"cannot import source: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def manifest_content_sha256(document: Mapping[str, Any]) -> str:
    payload = copy.deepcopy(dict(document))
    payload.setdefault("integrity", {}).pop("manifest_content_sha256", None)
    return semantic_sha256(payload)


def verify_census_bindings(census: Mapping[str, Any]) -> Dict[str, Dict[str, str]]:
    if census.get("schema_name") != "intersubmod.m2_solver_route_census_receipt":
        raise PanelExportError("unexpected census schema")
    if not census.get("all_pass"):
        raise PanelExportError("census receipt is not all-pass")
    bindings: Dict[str, Dict[str, str]] = {}
    for name, record in census["source_files"].items():
        path = pathlib.Path(record["path"])
        if not path.is_file():
            raise PanelExportError(f"bound census input is missing: {path}")
        observed = sha256_file(path)
        expected = str(record["sha256"])
        if observed != expected:
            raise PanelExportError(
                f"bound census input hash mismatch for {name}: {observed} != {expected}"
            )
        bindings[name] = {
            "path": str(path.resolve()),
            "sha256": observed,
        }
    return bindings


def structural_payload(
    *,
    k: int,
    universe_mask: int,
    full_patterns: Sequence[str],
    partial_patterns: Sequence[str],
) -> Dict[str, Any]:
    return {
        "k": int(k),
        "structural_alt_universe_mask": int(universe_mask),
        "full_patterns": sorted(set(full_patterns)),
        "partial_patterns": sorted(set(partial_patterns)),
    }


def audit_selected_record_blocks(
    record_keys: Sequence[Tuple[str, ...]],
    *,
    selected_keys: set[Tuple[str, ...]],
    required_keys: set[Tuple[str, ...]],
    stream_name: str,
) -> Dict[str, Any]:
    """Prove selected records are closed blocks before a truncated EOF block.

    A truncated line-oriented stream can only leave its terminal unit block
    uncertain.  Therefore a selected unit is replay-safe only when every
    recovered record for it is in one contiguous block and that block ends
    strictly before the final recovered record.  Candidate rows are optional
    for units whose source diagnostic reports zero emitted vertex sets, so
    ``required_keys`` may be a subset of ``selected_keys``.
    """
    if not record_keys:
        raise PanelExportError(f"{stream_name} has no complete records")
    positions: Dict[Tuple[str, ...], list[int]] = {
        key: [] for key in selected_keys
    }
    for index, key in enumerate(record_keys):
        if key in positions:
            positions[key].append(index)

    terminal_key = record_keys[-1]
    missing_required = sorted(
        [list(key) for key in required_keys if not positions.get(key)]
    )
    noncontiguous = sorted(
        [
            {
                "unit_key": list(key),
                "first_record": values[0],
                "last_record": values[-1],
                "record_count": len(values),
            }
            for key, values in positions.items()
            if values and values[-1] - values[0] + 1 != len(values)
        ],
        key=lambda row: row["unit_key"],
    )
    terminal_selected = sorted(
        [list(key) for key, values in positions.items() if values and key == terminal_key]
    )
    observed_blocks = [
        {
            "unit_key": list(key),
            "first_record": values[0],
            "last_record": values[-1],
            "record_count": len(values),
        }
        for key, values in sorted(positions.items())
        if values
    ]
    all_observed_before_terminal = all(
        block["last_record"] < len(record_keys) - 1 for block in observed_blocks
    )
    optional_absent = sorted(
        [
            list(key)
            for key in selected_keys - required_keys
            if not positions.get(key)
        ]
    )
    checks = {
        "required_selected_units_have_records": not missing_required,
        "selected_records_form_single_contiguous_block": not noncontiguous,
        "selected_unit_is_not_terminal_eof_unit": not terminal_selected,
        "selected_blocks_end_strictly_before_truncation_boundary": (
            all_observed_before_terminal
        ),
    }
    return {
        "stream": stream_name,
        "complete_record_count": len(record_keys),
        "selected_unit_count": len(selected_keys),
        "required_selected_unit_count": len(required_keys),
        "observed_selected_unit_count": len(observed_blocks),
        "optional_absent_unit_count": len(optional_absent),
        "optional_absent_unit_keys": optional_absent,
        "terminal_unit_key": list(terminal_key),
        "selected_block_descriptor_sha256": semantic_sha256(observed_blocks),
        "violations": {
            "missing_required": missing_required,
            "noncontiguous": noncontiguous,
            "terminal_selected": terminal_selected,
        },
        "checks": checks,
        "all_pass": all(checks.values()),
    }


def _factorial_oracle(frozen, payload: Mapping[str, Any]) -> Dict[str, Any] | None:
    full = list(payload["full_patterns"])
    partial = list(payload["partial_patterns"])
    mandatory = {0, *(frozen.parse_full(pattern) for pattern in full)}
    groups = tuple(frozen.SymbolicPattern.from_string(pattern) for pattern in partial)
    reduction = frozen._reduce_partial_groups(
        groups,
        int(payload["structural_alt_universe_mask"]),
        mandatory,
    )
    required = mandatory | set(reduction.forced_vertices)
    nonroot = sorted(required - {0})
    if len(nonroot) != 1 or reduction.groups:
        return None
    weight = bin(int(nonroot[0])).count("1")
    return {
        "kind": "single_terminal_permutation_family",
        "terminal_weight": weight,
        "expected_optimal_vertex_sets": math.factorial(weight),
    }


def build_manifest(
    census_path: pathlib.Path,
    *,
    census_source: pathlib.Path = DEFAULT_CENSUS_SOURCE,
    probe_source: pathlib.Path = DEFAULT_PROBE_SOURCE,
    optimized_source: pathlib.Path = DEFAULT_OPTIMIZED_SOURCE,
    runner_source: pathlib.Path = DEFAULT_RUNNER_SOURCE,
    authority_status: str = "DRAFT_UNPUBLISHED",
) -> Dict[str, Any]:
    census = json.loads(census_path.read_text(encoding="utf-8"))
    input_bindings = verify_census_bindings(census)
    census_module = load_module(census_source, "_solver_panel_census")
    frozen_path = pathlib.Path(input_bindings["hypercube_source"]["path"])
    frozen = census_module.load_frozen_hypercube(frozen_path)

    runtime_path = pathlib.Path(input_bindings["runtime"]["path"])
    pattern_path = pathlib.Path(input_bindings["patterns"]["path"])
    candidate_path = pathlib.Path(input_bindings["candidates"]["path"])

    patterns_by_unit: Dict[Tuple[str, ...], Dict[str, Any]] = {}
    pattern_record_keys: list[Tuple[str, ...]] = []
    pattern_stream, pattern_status = census_module.iter_truncated_tsv_gzip(
        pattern_path
    )
    for row in pattern_stream:
        key = census_module.unit_key(row)
        pattern_record_keys.append(key)
        entry = patterns_by_unit.setdefault(
            key,
            {"k": int(row["k"]), "retained_patterns": set()},
        )
        if entry["k"] != int(row["k"]):
            raise PanelExportError(f"inconsistent k in pattern rows: {list(key)}")
        if census_module.parse_bool(row["structural_retained"]):
            entry["retained_patterns"].add(row["pattern_rax"])

    candidates_by_unit: Dict[Tuple[str, ...], set[str]] = {}
    candidate_record_keys: list[Tuple[str, ...]] = []
    candidate_stream, candidate_status = census_module.iter_truncated_tsv_gzip(
        candidate_path
    )
    for row in candidate_stream:
        key = census_module.unit_key(row)
        candidate_record_keys.append(key)
        candidates_by_unit.setdefault(key, set()).add(row["vertex_set_id"])

    runtime_rows: list[Dict[str, Any]] = []
    runtime_record_keys: list[Tuple[str, ...]] = []
    runtime_stream, runtime_status = census_module.iter_truncated_tsv_gzip(
        runtime_path
    )
    for order, row in enumerate(runtime_stream):
        key = census_module.unit_key(row)
        runtime_record_keys.append(key)
        classification = census_module.classify_runtime(row)
        runtime_rows.append(
            {
                "order": order,
                "key": key,
                "classification": classification,
                "candidate_seconds": float(
                    row["candidate_generation_elapsed_seconds"]
                ),
                "candidate_vertex_sets": len(candidates_by_unit.get(key, set())),
            }
        )

    selected_rows: Dict[str, list[Dict[str, Any]]] = {
        "primary": [],
        "sensitivity": [],
    }
    missing: list[Dict[str, Any]] = []
    units: list[Dict[str, Any]] = []
    cases: Dict[str, Dict[str, Any]] = {}
    unit_to_case: Dict[str, str] = {}

    for panel_name, panel_config in PANEL_CONFIG["panels"].items():
        minread = int(panel_config["minread"])
        rows = [
            row
            for row in runtime_rows
            if int(row["key"][-1]) == minread
            and (
                row["classification"] == "tail"
                or (
                    row["classification"] == "complete"
                    and row["candidate_vertex_sets"] >= 100
                )
            )
        ]
        selected_rows[panel_name] = rows
        if len(rows) != int(panel_config["expected_units"]):
            raise PanelExportError(
                f"{panel_name} selected-unit count mismatch: "
                f"{len(rows)} != {panel_config['expected_units']}"
            )

        for row in rows:
            key = row["key"]
            pattern_entry = patterns_by_unit.get(key)
            if pattern_entry is None:
                missing.append(
                    {"unit_key": list(key), "reason": "MISSING_PATTERN_RECORDS"}
                )
                continue
            retained = sorted(pattern_entry["retained_patterns"])
            k = int(pattern_entry["k"])
            if not retained:
                missing.append(
                    {"unit_key": list(key), "reason": "NO_RETAINED_PATTERNS"}
                )
                continue
            bad_lengths = sorted(
                pattern for pattern in retained if len(pattern) != k
            )
            if bad_lengths:
                missing.append(
                    {
                        "unit_key": list(key),
                        "reason": "PATTERN_LENGTH_MISMATCH",
                        "patterns": bad_lengths,
                    }
                )
                continue
            full = sorted(pattern for pattern in retained if "X" not in pattern)
            partial = sorted(pattern for pattern in retained if "X" in pattern)
            universe_mask = int(frozen.observed_alt_universe(full, partial))
            payload = structural_payload(
                k=k,
                universe_mask=universe_mask,
                full_patterns=full,
                partial_patterns=partial,
            )
            key_sha = semantic_sha256(payload)
            case_id = f"k{k}_m{bin(universe_mask).count('1')}_{key_sha[:16]}"
            unit_payload = {
                field: value
                for field, value in zip(census_module.UNIT_FIELDS, key)
            }
            unit_id = semantic_sha256(
                {
                    "schema": "intersubmod.m2_unit_identity.v1",
                    "fields": unit_payload,
                }
            )
            reason = (
                "33_INCOMPLETE_CANDIDATE_GENERATION"
                if row["classification"] == "tail"
                else "COMPLETE_WITH_VERTEX_FAMILY_GE_100"
            )
            unit = {
                "unit_id": unit_id,
                "panel": panel_name,
                "identity": unit_payload,
                "source_runtime_order": int(row["order"]),
                "selection_reason": reason,
                "source_classification": row["classification"],
                "source_candidate_vertex_sets": int(row["candidate_vertex_sets"]),
                "source_candidate_seconds": float(row["candidate_seconds"]),
                "k": k,
                "full_patterns": full,
                "partial_patterns": partial,
                "structural_alt_universe_mask": universe_mask,
                "structural_key_sha256": key_sha,
                "case_id": case_id,
            }
            units.append(unit)
            unit_to_case[unit_id] = case_id
            case = cases.setdefault(
                case_id,
                {
                    "case_id": case_id,
                    "structural_key_sha256": key_sha,
                    "structural_input": payload,
                    "panels": [],
                    "selection_reasons": [],
                    "unit_ids": [],
                    "factorial_oracle": _factorial_oracle(frozen, payload),
                },
            )
            if case["structural_input"] != payload:
                raise PanelExportError(f"structural hash collision: {case_id}")
            case["panels"].append(panel_name)
            case["selection_reasons"].append(reason)
            case["unit_ids"].append(unit_id)

    if missing:
        raise PanelExportError(
            "selected cases cannot be reconstructed losslessly: "
            + json.dumps(missing, ensure_ascii=False, sort_keys=True)
        )

    for case in cases.values():
        case["panels"] = sorted(set(case["panels"]))
        case["selection_reasons"] = sorted(set(case["selection_reasons"]))
        case["unit_ids"] = sorted(case["unit_ids"])

    selected_unit_keys = {tuple(row["key"]) for rows in selected_rows.values() for row in rows}
    candidate_required_keys = {
        tuple(row["key"])
        for rows in selected_rows.values()
        for row in rows
        if int(row["candidate_vertex_sets"]) > 0
    }
    block_audits = {
        "runtime": audit_selected_record_blocks(
            runtime_record_keys,
            selected_keys=selected_unit_keys,
            required_keys=selected_unit_keys,
            stream_name="runtime",
        ),
        "patterns": audit_selected_record_blocks(
            pattern_record_keys,
            selected_keys=selected_unit_keys,
            required_keys=selected_unit_keys,
            stream_name="patterns",
        ),
        "candidates": audit_selected_record_blocks(
            candidate_record_keys,
            selected_keys=selected_unit_keys,
            required_keys=candidate_required_keys,
            stream_name="candidates",
        ),
    }
    failed_block_audits = [
        name for name, audit in block_audits.items() if not audit["all_pass"]
    ]
    if failed_block_audits:
        raise PanelExportError(
            "selected-unit block audit failed: "
            + json.dumps(
                {
                    name: block_audits[name]["violations"]
                    for name in failed_block_audits
                },
                ensure_ascii=False,
                sort_keys=True,
            )
        )

    panel_case_ids: Dict[str, list[str]] = {}
    for panel_name, rows in selected_rows.items():
        unit_ids = {
            semantic_sha256(
                {
                    "schema": "intersubmod.m2_unit_identity.v1",
                    "fields": {
                        field: value
                        for field, value in zip(census_module.UNIT_FIELDS, row["key"])
                    },
                }
            )
            for row in rows
        }
        case_ids = sorted({unit_to_case[unit_id] for unit_id in unit_ids})
        expected = int(
            PANEL_CONFIG["panels"][panel_name][
                "expected_distinct_structural_keys"
            ]
        )
        if len(case_ids) != expected:
            raise PanelExportError(
                f"{panel_name} structural-key count mismatch: "
                f"{len(case_ids)} != {expected}"
            )
        panel_case_ids[panel_name] = case_ids

    source_paths = {
        "exporter": pathlib.Path(__file__),
        "census": census_source,
        "solver_probe": probe_source,
        "optimized_backend": optimized_source,
        "stress_runner": runner_source,
        "frozen_current_solver": frozen_path,
    }
    for name, path in source_paths.items():
        if not path.is_file():
            raise PanelExportError(f"required source is missing ({name}): {path}")

    document: Dict[str, Any] = {
        "schema_name": "intersubmod.solver_stress_panel_manifest",
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_COMPREHENSIVE_VALIDATION_COMPONENT",
        "scope": {
            "dataset": PANEL_CONFIG["dataset"],
            "chrom": PANEL_CONFIG["chrom"],
            "eligibility": "DIAGNOSTIC_ONLY_COMPLETE_RECORDS_FROM_TRUNCATED_FROZEN_V2",
            "partial_flag": True,
            "not_claimed": [
                "complete chromosome-wide ranking output",
                "cross-sample performance",
                "production-backend promotion",
            ],
        },
        "immutability": {
            "policy": "WRITE_ONCE_REFUSE_EXISTING_OUTPUT",
            "manifest_sidecar": "manifest.json.sha256",
        },
        "authority": {
            "status": authority_status,
            "resolution_policy": (
                "RUNNER_MUST_USE_AUTHORITY_POINTER_OR_EXPLICIT_MATCHING_STATUS"
            ),
        },
        "config": PANEL_CONFIG,
        "bindings": {
            "config_sha256": semantic_sha256(PANEL_CONFIG),
            "census_receipt": {
                "path": str(census_path.resolve()),
                "sha256": sha256_file(census_path),
            },
            "input_files": input_bindings,
            "source_files": {
                name: {
                    "path": str(path.resolve()),
                    "sha256": sha256_file(path),
                }
                for name, path in source_paths.items()
            },
        },
        "truncation_audit": {
            "runtime": runtime_status,
            "patterns": pattern_status,
            "candidates": candidate_status,
            "selected_case_reconstruction": "LOSSLESS_FOR_SELECTED_COMPLETE_RECORDS",
            "unrecoverable_selected_units": [],
            "global_prefix_completeness": "NOT_CLAIMED",
            "selected_unit_block_audits": block_audits,
            "candidate_zero_row_policy": (
                "candidate records are optional only when the runtime-derived "
                "source_candidate_vertex_sets count is zero"
            ),
        },
        "panels": {
            panel_name: {
                **PANEL_CONFIG["panels"][panel_name],
                "case_ids": panel_case_ids[panel_name],
                "observed_units": len(selected_rows[panel_name]),
                "observed_distinct_structural_keys": len(
                    panel_case_ids[panel_name]
                ),
            }
            for panel_name in ("primary", "sensitivity")
        },
        "units": sorted(
            units,
            key=lambda unit: (
                unit["panel"],
                unit["source_runtime_order"],
                unit["unit_id"],
            ),
        ),
        "cases": [
            cases[case_id]
            for case_id in sorted(cases)
        ],
        "unit_to_case": dict(sorted(unit_to_case.items())),
        "checks": {
            "census_all_pass": True,
            "all_bound_input_hashes_match": True,
            "all_selected_units_losslessly_reconstructed": True,
            "primary_44_units_25_keys": (
                len(selected_rows["primary"]) == 44
                and len(panel_case_ids["primary"]) == 25
            ),
            "sensitivity_36_units_21_keys": (
                len(selected_rows["sensitivity"]) == 36
                and len(panel_case_ids["sensitivity"]) == 21
            ),
            "selected_units_not_terminal_eof_units": all(
                audit["checks"]["selected_unit_is_not_terminal_eof_unit"]
                for audit in block_audits.values()
            ),
            "selected_records_are_single_contiguous_blocks": all(
                audit["checks"][
                    "selected_records_form_single_contiguous_block"
                ]
                for audit in block_audits.values()
            ),
            "selected_blocks_complete_before_truncation": all(
                audit["checks"][
                    "selected_blocks_end_strictly_before_truncation_boundary"
                ]
                for audit in block_audits.values()
            ),
        },
        "integrity": {
            "manifest_content_hash_schema": (
                "sha256(canonical-json(document without "
                "integrity.manifest_content_sha256))"
            ),
        },
    }
    document["integrity"]["manifest_content_sha256"] = manifest_content_sha256(
        document
    )
    if not all(document["checks"].values()):
        raise PanelExportError("manifest checks did not all pass")
    return document


def write_once_manifest(document: Mapping[str, Any], output_dir: pathlib.Path) -> None:
    manifest_path = output_dir / "manifest.json"
    sidecar_path = output_dir / "manifest.json.sha256"
    if output_dir.exists():
        raise PanelExportError(f"write-once output already exists: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)
    encoded = (
        json.dumps(
            document,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")
    with manifest_path.open("xb") as handle:
        handle.write(encoded)
    exact_sha = hashlib.sha256(encoded).hexdigest()
    with sidecar_path.open("x", encoding="utf-8") as handle:
        handle.write(f"{exact_sha}  manifest.json\n")


def summary(document: Mapping[str, Any]) -> Dict[str, Any]:
    cases = list(document["cases"])
    return {
        "schema": document["schema_name"],
        "manifest_content_sha256": document["integrity"][
            "manifest_content_sha256"
        ],
        "primary_units": document["panels"]["primary"]["observed_units"],
        "primary_distinct_keys": document["panels"]["primary"][
            "observed_distinct_structural_keys"
        ],
        "sensitivity_units": document["panels"]["sensitivity"][
            "observed_units"
        ],
        "sensitivity_distinct_keys": document["panels"]["sensitivity"][
            "observed_distinct_structural_keys"
        ],
        "union_distinct_keys": len(cases),
        "factorial_oracle_cases": sum(
            case["factorial_oracle"] is not None for case in cases
        ),
        "all_checks_pass": all(document["checks"].values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", type=pathlib.Path, default=DEFAULT_CENSUS)
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=DEFAULT_OUTPUT_DIR,
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--authority-status",
        default="DRAFT_UNPUBLISHED",
        help="Authority label embedded in the immutable manifest.",
    )
    args = parser.parse_args()

    try:
        document = build_manifest(
            args.census,
            authority_status=args.authority_status,
        )
        if not args.dry_run:
            write_once_manifest(document, args.output_dir)
    except (PanelExportError, FileNotFoundError, KeyError, ValueError) as error:
        print(
            json.dumps(
                {
                    "status": "FAIL_CLOSED",
                    "error": str(error),
                },
                ensure_ascii=False,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 2

    result = summary(document)
    result["status"] = "DRY_RUN_PASS" if args.dry_run else "EXPORTED"
    result["output_dir"] = str(args.output_dir.resolve())
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
