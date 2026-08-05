#!/usr/bin/env python3
"""Run a bounded analytic perfect-family count census on 208 frozen keys.

The input is the immutable historical-v5 cross-sample panel.  This runner uses
only ``perfect_family_count.py``; it never imports or invokes the current or
optimized explicit family enumerators.  Every compatible case must terminate
as either an exact recurrence-free family count or an explicit abstention.

Timing fields are diagnostic provenance only.  They are not a formal speed
comparison because unrelated external workloads were active when this census
was authorized.
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
import time
from collections import Counter
from datetime import datetime, timezone
from typing import Any, Dict, Mapping, Sequence


REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]
PROBE_ROOT = REPO_ROOT / "research/20260718_solver_methyl_edge_probe"
DEFAULT_MANIFEST = (
    PROBE_ROOT
    / "results/cross_sample_solver_stress_panel_v1/manifest.json"
)
DEFAULT_OUTPUT_DIR = (
    PROBE_ROOT
    / "results/cross_sample_solver_stress_panel_v1/perfect_count_census_v1"
)
COUNTER_SOURCE = PROBE_ROOT / "scripts/perfect_family_count.py"
EXPECTED_CASES = 208
EXPECTED_MANIFEST_SCHEMA = (
    "intersubmod.historical_v5_cross_sample_solver_stress_panel_manifest"
)
CLAIM_SCOPE = "HISTORICAL_V5_SOLVER_CORE_ONLY"


class PerfectCountCensusError(RuntimeError):
    """Raised when a census contract or authority binding fails closed."""


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def semantic_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_module(path: pathlib.Path, name: str):
    if not path.is_file():
        raise PerfectCountCensusError(f"required source is missing: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise PerfectCountCensusError(f"cannot import source: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def document_content_sha256(
    document: Mapping[str, Any],
    *,
    field_name: str,
) -> str:
    payload = copy.deepcopy(dict(document))
    payload.setdefault("integrity", {}).pop(field_name, None)
    return semantic_sha256(payload)


def _verify_identity(record: Mapping[str, Any], *, label: str) -> None:
    path_value = record.get("path")
    if not isinstance(path_value, str):
        raise PerfectCountCensusError(f"{label} lacks a path")
    path = pathlib.Path(path_value)
    if not path.is_file():
        raise PerfectCountCensusError(f"{label} is missing: {path}")
    observed_sha = sha256_file(path)
    if observed_sha != record.get("sha256"):
        raise PerfectCountCensusError(
            f"{label} SHA mismatch: {observed_sha} != {record.get('sha256')}"
        )
    if "size_bytes" in record and path.stat().st_size != record.get("size_bytes"):
        raise PerfectCountCensusError(f"{label} size mismatch")


def verify_manifest_sidecar(manifest_path: pathlib.Path, manifest_sha: str) -> str:
    sidecar_path = manifest_path.with_name(manifest_path.name + ".sha256")
    if not sidecar_path.is_file():
        raise PerfectCountCensusError(
            f"input manifest sidecar is missing: {sidecar_path}"
        )
    expected_line = f"{manifest_sha}  {manifest_path.name}"
    observed_line = sidecar_path.read_text(encoding="utf-8").strip()
    if observed_line != expected_line:
        raise PerfectCountCensusError(
            f"input manifest sidecar mismatch: {observed_line!r}"
        )
    return sha256_file(sidecar_path)


def structural_case_sha256(structural: Mapping[str, Any]) -> str:
    allowed = {
        "k",
        "structural_alt_universe_mask",
        "full_patterns",
        "partial_patterns",
    }
    if set(structural) != allowed:
        raise PerfectCountCensusError(
            "structural_input fields are not exactly compatible with the counter"
        )
    return semantic_sha256(
        {
            "k": structural["k"],
            "structural_alt_universe_mask": structural[
                "structural_alt_universe_mask"
            ],
            "full_patterns": structural["full_patterns"],
            "partial_patterns": structural["partial_patterns"],
        }
    )


def validate_case(case: Mapping[str, Any]) -> None:
    structural = case.get("structural_input")
    if not isinstance(structural, Mapping):
        raise PerfectCountCensusError("case lacks structural_input")
    observed_sha = structural_case_sha256(structural)
    expected_sha = case.get("structural_key_sha256")
    if observed_sha != expected_sha:
        raise PerfectCountCensusError(
            f"case structural SHA mismatch: {observed_sha} != {expected_sha}"
        )
    k = structural.get("k")
    mask = structural.get("structural_alt_universe_mask")
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise PerfectCountCensusError("case k is not a positive integer")
    if isinstance(mask, bool) or not isinstance(mask, int):
        raise PerfectCountCensusError("case structural mask is not an integer")
    if mask <= 0 or mask >= (1 << k):
        raise PerfectCountCensusError("case structural mask lies outside k")
    if case.get("k") != k:
        raise PerfectCountCensusError("case k differs from structural_input k")
    effective_m = bin(mask).count("1")
    if case.get("effective_m") != effective_m:
        raise PerfectCountCensusError(
            "case effective_m differs from structural ALT mask"
        )
    if case.get("sample") not in {"H1437", "H2009"}:
        raise PerfectCountCensusError("case sample is outside the frozen panel")
    if not isinstance(case.get("unit_identity"), Mapping):
        raise PerfectCountCensusError("case lacks unit_identity provenance")


def verify_authority_manifest(
    manifest_path: pathlib.Path,
    *,
    expected_cases: int = EXPECTED_CASES,
    verify_bound_sources: bool = True,
) -> tuple[Dict[str, Any], Dict[str, Any]]:
    if not manifest_path.is_file():
        raise PerfectCountCensusError(
            f"input authority manifest is missing: {manifest_path}"
        )
    manifest_file_sha = sha256_file(manifest_path)
    sidecar_sha = verify_manifest_sidecar(manifest_path, manifest_file_sha)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PerfectCountCensusError(
            f"cannot parse input authority manifest: {error}"
        ) from error
    if manifest.get("schema_name") != EXPECTED_MANIFEST_SCHEMA:
        raise PerfectCountCensusError("unexpected input-manifest schema")
    authority = manifest.get("authority", {})
    if authority.get("status") != "AUTHORITATIVE":
        raise PerfectCountCensusError("input manifest is not authoritative")
    if authority.get("claim_scope") != CLAIM_SCOPE:
        raise PerfectCountCensusError("input manifest claim scope mismatch")
    if not all(manifest.get("checks", {}).values()):
        raise PerfectCountCensusError("input manifest checks are not all-pass")

    claimed_content_sha = manifest.get("integrity", {}).get(
        "manifest_content_sha256"
    )
    observed_content_sha = document_content_sha256(
        manifest,
        field_name="manifest_content_sha256",
    )
    if claimed_content_sha != observed_content_sha:
        raise PerfectCountCensusError(
            f"input manifest semantic hash mismatch: "
            f"{observed_content_sha} != {claimed_content_sha}"
        )

    cases = manifest.get("cases")
    if not isinstance(cases, list) or len(cases) != expected_cases:
        raise PerfectCountCensusError(
            f"input manifest case count mismatch: "
            f"{len(cases) if isinstance(cases, list) else 'not-list'} "
            f"!= {expected_cases}"
        )
    keys: list[str] = []
    payloads: list[Mapping[str, Any]] = []
    for case in cases:
        if not isinstance(case, Mapping):
            raise PerfectCountCensusError("input manifest has a non-object case")
        validate_case(case)
        keys.append(str(case["structural_key_sha256"]))
    if len(keys) != len(set(keys)):
        raise PerfectCountCensusError("input manifest has duplicate structural keys")
    for case in sorted(cases, key=lambda item: item["structural_key_sha256"]):
        payloads.append(case["structural_input"])
    sorted_keys = sorted(keys)
    if sorted_keys != manifest.get("selection", {}).get(
        "selected_structural_keys"
    ):
        raise PerfectCountCensusError(
            "case keys differ from the frozen selected key list"
        )
    if semantic_sha256(sorted_keys) != manifest.get("selection", {}).get(
        "key_list_sha256"
    ):
        raise PerfectCountCensusError("frozen key-list digest mismatch")
    if semantic_sha256(payloads) != manifest.get("selection", {}).get(
        "payload_list_sha256"
    ):
        raise PerfectCountCensusError("frozen payload-list digest mismatch")

    verified_source_count = 0
    if verify_bound_sources:
        bindings = manifest.get("bindings", {})
        identity_groups = {
            "origin_manifest": [bindings.get("origin_manifest", {})],
            "sample_output_manifests": list(
                bindings.get("sample_output_manifests", {}).values()
            ),
            "data_sources": list(bindings.get("data_sources", {}).values()),
            "code_sources": list(bindings.get("code_sources", {}).values()),
        }
        for group_name, records in identity_groups.items():
            if not records:
                raise PerfectCountCensusError(
                    f"input manifest lacks bound {group_name}"
                )
            for index, record in enumerate(records):
                if not isinstance(record, Mapping):
                    raise PerfectCountCensusError(
                        f"invalid bound identity in {group_name}"
                    )
                _verify_identity(record, label=f"{group_name}[{index}]")
                verified_source_count += 1
        if len(bindings.get("data_sources", {})) != 12:
            raise PerfectCountCensusError(
                "input manifest no longer binds exactly 12 data sources"
            )

    binding = {
        "path": str(manifest_path.resolve()),
        "size_bytes": manifest_path.stat().st_size,
        "file_sha256": manifest_file_sha,
        "sidecar_path": str(
            manifest_path.with_name(manifest_path.name + ".sha256").resolve()
        ),
        "sidecar_sha256": sidecar_sha,
        "semantic_content_sha256": observed_content_sha,
        "verified_bound_source_identities": verified_source_count,
        "case_count": len(cases),
        "key_list_sha256": semantic_sha256(sorted_keys),
        "payload_list_sha256": semantic_sha256(payloads),
    }
    return manifest, binding


def _validate_counter_result(result: Any, counter: Any) -> None:
    if result.status not in counter.EXACT_STATUSES | {counter.ABSTAIN_STATUS}:
        raise PerfectCountCensusError(
            f"counter returned an unsupported status: {result.status}"
        )
    counter.assert_finite_elapsed(result)
    if result.ranking_allowed is not False:
        raise PerfectCountCensusError("counter unexpectedly allows ranking")
    if isinstance(result.perfect_family_count, bool) or not isinstance(
        result.perfect_family_count, int
    ):
        raise PerfectCountCensusError("counter family count is not an integer")
    if result.perfect_family_count < 0:
        raise PerfectCountCensusError("counter family count is negative")
    if result.exact:
        if result.status not in counter.EXACT_STATUSES:
            raise PerfectCountCensusError("exact counter result has wrong status")
        if result.perfect_family_count <= 0:
            raise PerfectCountCensusError("exact counter result has zero count")
        if isinstance(result.objective, bool) or not isinstance(
            result.objective, int
        ):
            raise PerfectCountCensusError(
                "exact counter result lacks an integer objective"
            )
    else:
        if result.status != counter.ABSTAIN_STATUS:
            raise PerfectCountCensusError("abstain result has wrong status")
        if result.perfect_family_count != 0 or result.objective is not None:
            raise PerfectCountCensusError(
                "abstain result claims a family count or objective"
            )


def run_cases(
    cases: Sequence[Mapping[str, Any]],
    *,
    counter: Any,
    max_m: int,
) -> list[Dict[str, Any]]:
    rows: list[Dict[str, Any]] = []
    for case_index, case in enumerate(cases):
        validate_case(case)
        result = counter.count_manifest_case(case, max_m=max_m)
        _validate_counter_result(result, counter)
        if result.k != case["k"] or result.effective_m != case["effective_m"]:
            raise PerfectCountCensusError(
                "counter result dimensions differ from the frozen case"
            )
        row = {
            "case_index": case_index,
            "structural_key_sha256": case["structural_key_sha256"],
            "sample": case["sample"],
            "unit_identity": copy.deepcopy(case["unit_identity"]),
            "selection_reason": case["selection_reason"],
            "historical_cap_class": case["historical_cap_class"],
            "reduced_q": case["reduced_q"],
            **result.to_dict(),
            "formal_speed_claim": False,
        }
        if row["ranking_allowed"] is not False:
            raise PerfectCountCensusError("row unexpectedly allows ranking")
        rows.append(row)
    return rows


def summarize_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    counter: Any,
    expected_cases: int = EXPECTED_CASES,
) -> Dict[str, Any]:
    if len(rows) != expected_cases:
        raise PerfectCountCensusError(
            f"census row count mismatch: {len(rows)} != {expected_cases}"
        )
    if len({row["structural_key_sha256"] for row in rows}) != len(rows):
        raise PerfectCountCensusError("census rows contain duplicate keys")
    if not all(row.get("ranking_allowed") is False for row in rows):
        raise PerfectCountCensusError("one or more census rows allows ranking")
    if not all(row.get("formal_speed_claim") is False for row in rows):
        raise PerfectCountCensusError(
            "one or more census rows makes a formal speed claim"
        )
    for row in rows:
        elapsed = row.get("elapsed_seconds")
        if (
            isinstance(elapsed, bool)
            or not isinstance(elapsed, (int, float))
            or elapsed < 0
            or not math.isfinite(float(elapsed))
        ):
            raise PerfectCountCensusError("row elapsed time is not finite")

    exact_rows = [row for row in rows if row["status"] in counter.EXACT_STATUSES]
    abstain_rows = [
        row for row in rows if row["status"] == counter.ABSTAIN_STATUS
    ]
    if len(exact_rows) + len(abstain_rows) != len(rows):
        raise PerfectCountCensusError("census is not exact-or-abstain total")
    maximum = max(
        exact_rows,
        key=lambda row: (
            int(row["perfect_family_count"]),
            row["structural_key_sha256"],
        ),
        default=None,
    )

    def distribution(
        source_rows: Sequence[Mapping[str, Any]],
        field: str,
    ) -> Dict[str, int]:
        return dict(
            sorted(
                Counter(str(row[field]) for row in source_rows).items()
            )
        )

    return {
        "total_cases": len(rows),
        "exact_cases": len(exact_rows),
        "abstain_cases": len(abstain_rows),
        "status_distribution": dict(
            sorted(Counter(row["status"] for row in rows).items())
        ),
        "exact_by_sample": distribution(exact_rows, "sample"),
        "abstain_distribution": {
            "by_sample": distribution(abstain_rows, "sample"),
            "by_effective_m": distribution(abstain_rows, "effective_m"),
            "by_historical_cap_class": distribution(
                abstain_rows,
                "historical_cap_class",
            ),
            "by_reduced_q": distribution(abstain_rows, "reduced_q"),
        },
        "objective_distribution_exact": distribution(exact_rows, "objective"),
        "counter_elapsed_seconds_sum": math.fsum(
            float(row["elapsed_seconds"]) for row in rows
        ),
        "counter_elapsed_seconds_max": max(
            (float(row["elapsed_seconds"]) for row in rows),
            default=0.0,
        ),
        "maximum_exact_family_count": (
            {
                "structural_key_sha256": maximum["structural_key_sha256"],
                "sample": maximum["sample"],
                "effective_m": maximum["effective_m"],
                "objective": maximum["objective"],
                "perfect_family_count": maximum["perfect_family_count"],
            }
            if maximum is not None
            else None
        ),
        "all_rows_total": len(exact_rows) + len(abstain_rows) == len(rows),
        "all_elapsed_finite": True,
        "ranking_allowed": False,
        "formal_speed_claim": False,
    }


def build_receipt(
    manifest_path: pathlib.Path = DEFAULT_MANIFEST,
    *,
    max_m: int = 20,
    expected_cases: int = EXPECTED_CASES,
    verify_bound_sources: bool = True,
) -> Dict[str, Any]:
    wall_started = time.perf_counter()
    manifest, manifest_binding = verify_authority_manifest(
        manifest_path,
        expected_cases=expected_cases,
        verify_bound_sources=verify_bound_sources,
    )
    counter = load_module(COUNTER_SOURCE, "perfect_family_count")
    rows = run_cases(
        manifest["cases"],
        counter=counter,
        max_m=max_m,
    )
    census = summarize_rows(
        rows,
        counter=counter,
        expected_cases=expected_cases,
    )
    runner_wall = time.perf_counter() - wall_started
    if runner_wall < 0 or not math.isfinite(runner_wall):
        raise PerfectCountCensusError("runner wall time is not finite")

    code_sources = {
        "counter": COUNTER_SOURCE,
        "runner": pathlib.Path(__file__),
    }
    for name, path in code_sources.items():
        if not path.is_file():
            raise PerfectCountCensusError(f"{name} source is missing: {path}")

    checks = {
        "input_manifest_file_sha256_verified": True,
        "input_manifest_semantic_content_sha256_verified": True,
        "input_manifest_sidecar_verified": True,
        "input_manifest_bound_sources_verified": (
            manifest_binding["verified_bound_source_identities"] > 0
            if verify_bound_sources
            else True
        ),
        "case_contract_compatible_208_of_208": len(rows) == expected_cases,
        "exact_or_abstain_total_208_of_208": census["all_rows_total"],
        "all_elapsed_finite": census["all_elapsed_finite"],
        "all_rows_ranking_disallowed": all(
            row["ranking_allowed"] is False for row in rows
        ),
        "formal_speed_claim_false": census["formal_speed_claim"] is False,
        "explicit_family_solvers_not_imported_or_invoked": True,
        "canonical_pipeline_not_modified_or_invoked": True,
    }
    if not all(checks.values()):
        raise PerfectCountCensusError("one or more census checks failed")

    receipt: Dict[str, Any] = {
        "schema_name": (
            "intersubmod.historical_v5_cross_sample_perfect_count_census_receipt"
        ),
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authority": {
            "status": "AUTHORITATIVE_BOUNDED_CENSUS",
            "claim_scope": CLAIM_SCOPE,
            "write_policy": "WRITE_ONCE_REFUSE_EXISTING_OUTPUT",
        },
        "task_type": "B_COMPREHENSIVE_VALIDATION_COMPONENT",
        "method": {
            "counter": "analytic recurrence-free perfect-family subset DP",
            "counter_scope": counter.SCOPE,
            "max_effective_m": int(max_m),
            "result_contract": "EXACT_OR_ABSTAIN_TOTAL",
            "explicit_family_materialization": False,
            "current_solver_invoked": False,
            "optimized_solver_invoked": False,
            "ranking_allowed": False,
        },
        "performance_context": {
            "runner_wall_seconds": runner_wall,
            "counter_elapsed_seconds_sum": census[
                "counter_elapsed_seconds_sum"
            ],
            "counter_elapsed_seconds_max": census[
                "counter_elapsed_seconds_max"
            ],
            "external_workload_present_at_authorization": True,
            "formal_speed_claim": False,
            "interpretation": (
                "diagnostic bounded census timing only; not comparable solver "
                "speedup evidence"
            ),
        },
        "bindings": {
            "input_authority_manifest": manifest_binding,
            "source_files": {
                name: {
                    "path": str(path.resolve()),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for name, path in code_sources.items()
            },
        },
        "census": census,
        "rows": rows,
        "checks": checks,
        "integrity": {
            "receipt_content_hash_schema": (
                "sha256(canonical JSON document without "
                "integrity.receipt_content_sha256)"
            ),
        },
    }
    receipt["integrity"]["receipt_content_sha256"] = (
        document_content_sha256(
            receipt,
            field_name="receipt_content_sha256",
        )
    )
    return receipt


def write_once_receipt(
    receipt: Mapping[str, Any],
    output_dir: pathlib.Path,
) -> None:
    if output_dir.exists():
        raise PerfectCountCensusError(
            f"write-once output already exists: {output_dir}"
        )
    output_dir.mkdir(parents=True, exist_ok=False)
    receipt_path = output_dir / "receipt.json"
    sidecar_path = output_dir / "receipt.json.sha256"
    encoded = (
        json.dumps(
            receipt,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")
    with receipt_path.open("xb") as handle:
        handle.write(encoded)
    receipt_sha = hashlib.sha256(encoded).hexdigest()
    with sidecar_path.open("x", encoding="utf-8") as handle:
        handle.write(f"{receipt_sha}  receipt.json\n")


def summary(
    receipt: Mapping[str, Any],
    *,
    status: str,
    output_dir: pathlib.Path,
) -> Dict[str, Any]:
    census = receipt["census"]
    return {
        "status": status,
        "claim_scope": receipt["authority"]["claim_scope"],
        "output_dir": str(output_dir.resolve()),
        "total_cases": census["total_cases"],
        "exact_cases": census["exact_cases"],
        "abstain_cases": census["abstain_cases"],
        "status_distribution": census["status_distribution"],
        "maximum_exact_family_count": census["maximum_exact_family_count"],
        "abstain_distribution": census["abstain_distribution"],
        "runner_wall_seconds": receipt["performance_context"][
            "runner_wall_seconds"
        ],
        "counter_elapsed_seconds_sum": census[
            "counter_elapsed_seconds_sum"
        ],
        "ranking_allowed": False,
        "formal_speed_claim": False,
        "all_checks_pass": all(receipt["checks"].values()),
        "receipt_content_sha256": receipt["integrity"][
            "receipt_content_sha256"
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=pathlib.Path,
        default=DEFAULT_MANIFEST,
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=DEFAULT_OUTPUT_DIR,
    )
    parser.add_argument("--max-m", type=int, default=20)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    try:
        receipt = build_receipt(
            args.manifest,
            max_m=args.max_m,
        )
        if not args.dry_run:
            write_once_receipt(receipt, args.output_dir)
    except (
        PerfectCountCensusError,
        FileNotFoundError,
        KeyError,
        TypeError,
        ValueError,
    ) as error:
        print(
            json.dumps(
                {"status": "FAIL_CLOSED", "error": str(error)},
                ensure_ascii=False,
                sort_keys=True,
            ),
            file=sys.stderr,
        )
        return 2

    print(
        json.dumps(
            summary(
                receipt,
                status="DRY_RUN_PASS" if args.dry_run else "EXPORTED",
                output_dir=args.output_dir,
            ),
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
