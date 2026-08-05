#!/usr/bin/env python3
"""Freeze a source-bound H1437/H2009 historical-v5 solver stress panel.

This exporter is deliberately limited to the isolated research solver core.
It reconstructs symbolic structural inputs from the canonical historical-v5
``subread_groups_by_hp`` summaries, selects a pre-registered deterministic
panel, and writes one immutable authority manifest.  It neither invokes a
solver nor claims equivalence to the M2 read-linked component/PS route.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import pathlib
import sys
from collections import Counter
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, Mapping, Sequence, Tuple


REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]
PROBE_ROOT = REPO_ROOT / "research/20260718_solver_methyl_edge_probe"
RUN_ROOT = pathlib.Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260713_layered_reconstruction_v3_raw_all_lps_pass_v5"
)
DEFAULT_OUTPUT_DIR = (
    PROBE_ROOT / "results/cross_sample_solver_stress_panel_v1"
)
SOLVER_PROBE_SOURCE = PROBE_ROOT / "scripts/solver_probe.py"
OPTIMIZED_BACKEND_SOURCE = (
    PROBE_ROOT / "scripts/optimized_hypercube_backend.py"
)

SAMPLES = ("H1437", "H2009")
EXPECTED_ORIGIN_MANIFEST_SHA256 = (
    "16f2ef66634e8592e32e5088d8383d94dead0ae2b0d32847f4d8843f8bdc1a45"
)
EXPECTED_OUTPUT_MANIFEST_SHA256 = {
    "H1437": (
        "5ddbb14bbc467f1908c9024840bbd95bbec4d6bcd1a3fbeb9f75e561802a4c70"
    ),
    "H2009": (
        "a1e9f1995ff51e569cf6314d07045e31841ba6202aa2650625c8452279baa2b8"
    ),
}
EXPECTED_SAMPLE_COUNTS = {
    "H1437": {
        "incomplete_units": 1773,
        "distinct_structural_keys": 1763,
        "q_gt_8_keys": 10,
        "selected_keys": 86,
    },
    "H2009": {
        "incomplete_units": 3759,
        "distinct_structural_keys": 3746,
        "q_gt_8_keys": 38,
        "selected_keys": 122,
    },
}
EXPECTED_UNION_STRUCTURAL_KEYS = 5506
EXPECTED_SELECTED_KEYS = 208
EXPECTED_Q_GT_8_KEYS = 48
EXPECTED_KEY_LIST_BYTES = 13937
EXPECTED_KEY_LIST_SHA256 = (
    "7ae4f85ba9f05f146303d3f7921e5f333f395c690cdcd1dbd84ea39037b75692"
)
EXPECTED_PAYLOAD_LIST_BYTES = 58243
EXPECTED_PAYLOAD_LIST_SHA256 = (
    "30a1a767f2bdf6c83b30c383d941748bd933207bc93cd80a4b8a49ca0acbd682"
)


class CrossSamplePanelError(RuntimeError):
    """Raised when the historical panel cannot be reconstructed exactly."""


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


def load_json(path: pathlib.Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CrossSamplePanelError(f"cannot read JSON source {path}: {error}") from error


def load_module(path: pathlib.Path, name: str):
    if not path.is_file():
        raise CrossSamplePanelError(f"required Python source is missing: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise CrossSamplePanelError(f"cannot import Python source: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_solver_modules():
    probe = load_module(SOLVER_PROBE_SOURCE, "solver_probe")
    optimized = load_module(
        OPTIMIZED_BACKEND_SOURCE,
        "optimized_hypercube_backend",
    )
    return probe, optimized


def structural_payload(
    pattern_counts: Mapping[str, Any],
    *,
    k: int,
    minread: int = 3,
) -> Dict[str, Any]:
    """Reconstruct the exact historical symbolic input used by the audit."""
    structural: list[str] = []
    for pattern, raw_count in pattern_counts.items():
        if isinstance(raw_count, bool) or not isinstance(raw_count, int):
            raise CrossSamplePanelError(
                f"pattern count must be an integer: {pattern}={raw_count!r}"
            )
        if raw_count < 0:
            raise CrossSamplePanelError(
                f"pattern count must be nonnegative: {pattern}={raw_count}"
            )
        if not isinstance(pattern, str) or len(pattern) != int(k):
            raise CrossSamplePanelError(
                f"pattern length differs from k={k}: {pattern!r}"
            )
        invalid = set(pattern) - {"R", "A", "X"}
        if invalid:
            raise CrossSamplePanelError(
                f"invalid pattern symbols in {pattern!r}: {sorted(invalid)}"
            )
        if raw_count >= minread and set(pattern) != {"X"}:
            structural.append(pattern)
    structural = sorted(set(structural))
    if not structural:
        raise CrossSamplePanelError("no structural patterns remain after filtering")

    full = sorted(pattern for pattern in structural if "X" not in pattern)
    partial = sorted(pattern for pattern in structural if "X" in pattern)
    universe_mask = 0
    for pattern in structural:
        for index, symbol in enumerate(pattern):
            if symbol == "A":
                universe_mask |= 1 << index
    if universe_mask == 0:
        raise CrossSamplePanelError("primary mutation-bearing unit has no ALT dimension")
    return {
        "k": int(k),
        "structural_alt_universe_mask": int(universe_mask),
        "full_patterns": full,
        "partial_patterns": partial,
    }


def historical_cap_class(cap_reason: Any) -> str:
    if not isinstance(cap_reason, str):
        raise CrossSamplePanelError(f"missing historical cap reason: {cap_reason!r}")
    if cap_reason.startswith("no feasible N within extra_cap=4"):
        return "EXTRA_CAP_4"
    if cap_reason.startswith("level e=3:"):
        return "LEVEL_BUDGET_E3"
    if cap_reason.startswith("level e=4:"):
        return "LEVEL_BUDGET_E4"
    raise CrossSamplePanelError(f"unrecognized historical cap reason: {cap_reason}")


def manifest_content_sha256(document: Mapping[str, Any]) -> str:
    payload = copy.deepcopy(dict(document))
    payload.setdefault("integrity", {}).pop("manifest_content_sha256", None)
    return semantic_sha256(payload)


def verify_origin(run_root: pathlib.Path) -> Dict[str, Any]:
    success_path = run_root / "_SUCCESS"
    if not success_path.is_file():
        raise CrossSamplePanelError(f"historical run lacks _SUCCESS: {success_path}")
    origin_path = run_root / "input_manifest.snapshot.json"
    if not origin_path.is_file():
        raise CrossSamplePanelError(f"origin manifest is missing: {origin_path}")
    observed_sha = sha256_file(origin_path)
    if observed_sha != EXPECTED_ORIGIN_MANIFEST_SHA256:
        raise CrossSamplePanelError(
            f"origin manifest hash mismatch: {observed_sha} != "
            f"{EXPECTED_ORIGIN_MANIFEST_SHA256}"
        )
    origin = load_json(origin_path)
    if origin.get("schema_name") != "intersubmod.layered_input_manifest":
        raise CrossSamplePanelError("unexpected origin-manifest schema")
    if origin.get("production_summary", {}).get("all_pass") is not True:
        raise CrossSamplePanelError("origin manifest is not all-pass")
    return {
        "path": str(origin_path.resolve()),
        "size_bytes": origin_path.stat().st_size,
        "sha256": observed_sha,
    }


def verify_sample_source_bindings(
    run_root: pathlib.Path,
    sample: str,
) -> Tuple[Dict[str, Any], Dict[str, Dict[str, Any]]]:
    sample_root = run_root / "samples" / sample
    output_manifest_path = sample_root / "output_manifest.json"
    if not output_manifest_path.is_file():
        raise CrossSamplePanelError(
            f"sample output manifest is missing: {output_manifest_path}"
        )
    observed_manifest_sha = sha256_file(output_manifest_path)
    expected_manifest_sha = EXPECTED_OUTPUT_MANIFEST_SHA256[sample]
    if observed_manifest_sha != expected_manifest_sha:
        raise CrossSamplePanelError(
            f"{sample} output manifest hash mismatch: "
            f"{observed_manifest_sha} != {expected_manifest_sha}"
        )
    output_manifest = load_json(output_manifest_path)
    if output_manifest.get("sample") != sample:
        raise CrossSamplePanelError(
            f"output manifest sample mismatch: {output_manifest.get('sample')} != {sample}"
        )
    if output_manifest.get("schema_name") != "intersubmod.layered_sample_output_manifest":
        raise CrossSamplePanelError(f"unexpected output-manifest schema for {sample}")

    by_role: Dict[str, Mapping[str, Any]] = {}
    for record in output_manifest.get("outputs", []):
        role = record.get("role")
        if not isinstance(role, str) or role in by_role:
            raise CrossSamplePanelError(
                f"missing or duplicate output role in {sample}: {role!r}"
            )
        by_role[role] = record

    required_roles = [
        *(f"mlhp_part_{index}" for index in range(1, 6)),
        "layered_reconstruction",
    ]
    source_bindings: Dict[str, Dict[str, Any]] = {}
    for role in required_roles:
        if role not in by_role:
            raise CrossSamplePanelError(f"{sample} manifest lacks source role {role}")
        record = by_role[role]
        relative_path = record.get("path")
        if not isinstance(relative_path, str) or pathlib.Path(relative_path).is_absolute():
            raise CrossSamplePanelError(
                f"{sample}/{role} has invalid relative path: {relative_path!r}"
            )
        path = sample_root / relative_path
        if not path.is_file():
            raise CrossSamplePanelError(f"bound source is missing: {path}")
        observed_sha = sha256_file(path)
        observed_size = path.stat().st_size
        if observed_sha != record.get("sha256"):
            raise CrossSamplePanelError(
                f"{sample}/{role} hash differs from output manifest"
            )
        if observed_size != record.get("size_bytes"):
            raise CrossSamplePanelError(
                f"{sample}/{role} size differs from output manifest"
            )
        source_bindings[role] = {
            "path": str(path.resolve()),
            "size_bytes": observed_size,
            "sha256": observed_sha,
            "output_manifest_role": role,
        }
    return (
        {
            "path": str(output_manifest_path.resolve()),
            "size_bytes": output_manifest_path.stat().st_size,
            "sha256": observed_manifest_sha,
        },
        source_bindings,
    )


def select_cases(
    rows: Sequence[Mapping[str, Any]],
) -> Tuple[list[Dict[str, Any]], Dict[str, Any]]:
    """Apply the pre-registered q-tail plus deterministic stratum rule."""
    q_tail = [dict(row) for row in rows if int(row["reduced_q"]) > 8]
    strata: Dict[Tuple[Any, ...], Dict[str, Any]] = {}
    for source_row in rows:
        row = dict(source_row)
        if int(row["reduced_q"]) > 8:
            continue
        stratum = (
            row["sample"],
            int(row["k"]),
            int(row["effective_m"]),
            int(row["reduced_q"]),
            row["historical_cap_class"],
        )
        previous = strata.get(stratum)
        if previous is None or row["structural_key_sha256"] < previous[
            "structural_key_sha256"
        ]:
            strata[stratum] = row

    selected_memberships = q_tail + list(strata.values())
    selected_memberships.sort(
        key=lambda row: (row["sample"], row["structural_key_sha256"])
    )
    selected_hashes = [
        row["structural_key_sha256"] for row in selected_memberships
    ]
    if len(selected_hashes) != len(set(selected_hashes)):
        raise CrossSamplePanelError(
            "selected structural keys overlap across sample memberships"
        )
    selection_summary = {
        "q_gt_8_memberships": len(q_tail),
        "strata": len(strata),
        "selected_memberships": len(selected_memberships),
        "selected_by_sample": dict(
            sorted(Counter(row["sample"] for row in selected_memberships).items())
        ),
    }
    return selected_memberships, selection_summary


def reconstruct_sample(
    run_root: pathlib.Path,
    sample: str,
    *,
    probe: Any,
    optimized: Any,
) -> Tuple[list[Dict[str, Any]], Dict[str, Any], Dict[str, Any]]:
    output_manifest_binding, source_bindings = verify_sample_source_bindings(
        run_root,
        sample,
    )

    groups: Dict[Tuple[str, int, int], Mapping[str, Any]] = {}
    groups_per_part: Dict[str, int] = {}
    for part in range(1, 6):
        role = f"mlhp_part_{part}"
        path = pathlib.Path(source_bindings[role]["path"])
        document = load_json(path)
        if document.get("sample") != sample or document.get("part") != part:
            raise CrossSamplePanelError(
                f"{sample}/{role} sample or part identity mismatch"
            )
        part_groups = document.get("groups")
        if not isinstance(part_groups, list):
            raise CrossSamplePanelError(f"{sample}/{role} lacks a groups list")
        if document.get("n_groups_analyzed") != len(part_groups):
            raise CrossSamplePanelError(
                f"{sample}/{role} group count differs from n_groups_analyzed"
            )
        groups_per_part[role] = len(part_groups)
        for group in part_groups:
            try:
                key = (
                    str(group["chrom"]),
                    int(group["start"]),
                    int(group["end"]),
                )
            except (KeyError, TypeError, ValueError) as error:
                raise CrossSamplePanelError(
                    f"invalid raw group identity in {sample}/{role}"
                ) from error
            if key in groups:
                raise CrossSamplePanelError(
                    f"duplicate raw group across parts: {sample}/{key}"
                )
            groups[key] = group

    layered_path = pathlib.Path(source_bindings["layered_reconstruction"]["path"])
    layered = load_json(layered_path)
    if layered.get("sample") != sample:
        raise CrossSamplePanelError(f"layered sample mismatch for {sample}")
    details = layered.get("detail")
    if not isinstance(details, list):
        raise CrossSamplePanelError(f"layered detail is not a list for {sample}")
    incomplete = [
        unit
        for unit in details
        if unit.get("is_primary_lineage") is True
        and unit.get("analysis_candidate_set_complete") is False
    ]

    expected = EXPECTED_SAMPLE_COUNTS[sample]
    if len(incomplete) != expected["incomplete_units"]:
        raise CrossSamplePanelError(
            f"{sample} incomplete-unit count mismatch: "
            f"{len(incomplete)} != {expected['incomplete_units']}"
        )

    rows: list[Dict[str, Any]] = []
    sample_seen: set[str] = set()
    unit_identities: set[Tuple[str, int, int, str]] = set()
    duplicate_structural_units = 0
    missing_pattern_joins: list[Dict[str, Any]] = []

    for detail_order, unit in enumerate(incomplete):
        try:
            group_key = (
                str(unit["chrom"]),
                int(unit["start"]),
                int(unit["end"]),
            )
            family = str(unit["family"])
            k = int(unit["n_sSNV"])
        except (KeyError, TypeError, ValueError) as error:
            raise CrossSamplePanelError(
                f"invalid incomplete-unit identity in {sample}"
            ) from error
        unit_identity = (*group_key, family)
        if unit_identity in unit_identities:
            raise CrossSamplePanelError(
                f"duplicate incomplete unit identity: {sample}/{unit_identity}"
            )
        unit_identities.add(unit_identity)
        group = groups.get(group_key)
        if group is None:
            missing_pattern_joins.append(
                {
                    "detail_order": detail_order,
                    "identity": list(unit_identity),
                    "reason": "MISSING_RAW_GROUP",
                }
            )
            continue
        pattern_counts = (
            ((group.get("subread_groups_by_hp") or {}).get(family)) or {}
        )
        if not isinstance(pattern_counts, Mapping) or not pattern_counts:
            missing_pattern_joins.append(
                {
                    "detail_order": detail_order,
                    "identity": list(unit_identity),
                    "reason": "MISSING_FAMILY_PATTERN_COUNTS",
                }
            )
            continue

        payload = structural_payload(pattern_counts, k=k)
        structural_hash = semantic_sha256(payload)
        if structural_hash in sample_seen:
            duplicate_structural_units += 1
            continue
        sample_seen.add(structural_hash)

        instance = probe.build_instance(
            payload["full_patterns"],
            payload["partial_patterns"],
            k,
            universe_mode="observed_alt",
        )
        if int(instance.universe_mask) != int(
            payload["structural_alt_universe_mask"]
        ):
            raise CrossSamplePanelError(
                f"solver universe differs from structural payload: {sample}/{unit_identity}"
            )
        prepared = optimized.prepare_bitset_problem(instance)
        cap_class = historical_cap_class(unit.get("cap_reason"))
        row = {
            "sample": sample,
            "detail_order": detail_order,
            "unit_identity": {
                "chrom": group_key[0],
                "start": group_key[1],
                "end": group_key[2],
                "family": family,
            },
            "k": k,
            "effective_m": bin(int(instance.universe_mask)).count("1"),
            "reduced_q": int(prepared.q),
            "historical_cap_class": cap_class,
            "selection_reason": None,
            "structural_key_sha256": structural_hash,
            "structural_input": payload,
        }
        rows.append(row)

    if missing_pattern_joins:
        raise CrossSamplePanelError(
            "incomplete units cannot be joined losslessly: "
            + json.dumps(
                missing_pattern_joins[:20],
                ensure_ascii=False,
                sort_keys=True,
            )
        )
    if len(rows) != expected["distinct_structural_keys"]:
        raise CrossSamplePanelError(
            f"{sample} distinct structural-key count mismatch: "
            f"{len(rows)} != {expected['distinct_structural_keys']}"
        )
    q_gt_8 = sum(int(row["reduced_q"]) > 8 for row in rows)
    if q_gt_8 != expected["q_gt_8_keys"]:
        raise CrossSamplePanelError(
            f"{sample} q>8 structural-key count mismatch: "
            f"{q_gt_8} != {expected['q_gt_8_keys']}"
        )
    sample_summary = {
        "historical_incomplete_units": len(incomplete),
        "distinct_structural_keys": len(rows),
        "duplicate_structural_units": duplicate_structural_units,
        "q_gt_8_keys": q_gt_8,
        "missing_pattern_joins": 0,
        "groups_per_part": groups_per_part,
        "raw_groups_total": len(groups),
    }
    return rows, sample_summary, {
        "output_manifest": output_manifest_binding,
        "sources": source_bindings,
    }


def build_manifest(run_root: pathlib.Path = RUN_ROOT) -> Dict[str, Any]:
    origin_binding = verify_origin(run_root)
    probe, optimized = load_solver_modules()

    all_rows: list[Dict[str, Any]] = []
    sample_summaries: Dict[str, Any] = {}
    source_manifests: Dict[str, Any] = {}
    data_sources: Dict[str, Any] = {}
    for sample in SAMPLES:
        rows, sample_summary, bindings = reconstruct_sample(
            run_root,
            sample,
            probe=probe,
            optimized=optimized,
        )
        all_rows.extend(rows)
        sample_summaries[sample] = sample_summary
        source_manifests[sample] = bindings["output_manifest"]
        for role, record in bindings["sources"].items():
            data_sources[f"{sample}/{role}"] = record

    if len(data_sources) != 12:
        raise CrossSamplePanelError(
            f"bound source-file count mismatch: {len(data_sources)} != 12"
        )
    union_hashes = {
        row["structural_key_sha256"] for row in all_rows
    }
    if len(union_hashes) != EXPECTED_UNION_STRUCTURAL_KEYS:
        raise CrossSamplePanelError(
            f"union structural-key count mismatch: "
            f"{len(union_hashes)} != {EXPECTED_UNION_STRUCTURAL_KEYS}"
        )

    selected, selection_summary = select_cases(all_rows)
    if selection_summary["selected_memberships"] != EXPECTED_SELECTED_KEYS:
        raise CrossSamplePanelError(
            f"selected-key count mismatch: "
            f"{selection_summary['selected_memberships']} != {EXPECTED_SELECTED_KEYS}"
        )
    if selection_summary["q_gt_8_memberships"] != EXPECTED_Q_GT_8_KEYS:
        raise CrossSamplePanelError(
            f"selected q>8 count mismatch: "
            f"{selection_summary['q_gt_8_memberships']} != {EXPECTED_Q_GT_8_KEYS}"
        )
    for sample in SAMPLES:
        observed = selection_summary["selected_by_sample"].get(sample, 0)
        expected = EXPECTED_SAMPLE_COUNTS[sample]["selected_keys"]
        if observed != expected:
            raise CrossSamplePanelError(
                f"{sample} selected-key count mismatch: {observed} != {expected}"
            )

    selected_by_hash = {
        row["structural_key_sha256"]: row for row in selected
    }
    selected_hashes = sorted(selected_by_hash)
    selected_payloads = [
        selected_by_hash[key]["structural_input"] for key in selected_hashes
    ]
    key_preimage = canonical_json_bytes(selected_hashes)
    payload_preimage = canonical_json_bytes(selected_payloads)
    key_digest = hashlib.sha256(key_preimage).hexdigest()
    payload_digest = hashlib.sha256(payload_preimage).hexdigest()
    if len(key_preimage) != EXPECTED_KEY_LIST_BYTES:
        raise CrossSamplePanelError(
            f"key-list preimage size mismatch: "
            f"{len(key_preimage)} != {EXPECTED_KEY_LIST_BYTES}"
        )
    if key_digest != EXPECTED_KEY_LIST_SHA256:
        raise CrossSamplePanelError(
            f"key-list digest mismatch: {key_digest} != {EXPECTED_KEY_LIST_SHA256}"
        )
    if len(payload_preimage) != EXPECTED_PAYLOAD_LIST_BYTES:
        raise CrossSamplePanelError(
            f"payload-list preimage size mismatch: "
            f"{len(payload_preimage)} != {EXPECTED_PAYLOAD_LIST_BYTES}"
        )
    if payload_digest != EXPECTED_PAYLOAD_LIST_SHA256:
        raise CrossSamplePanelError(
            f"payload-list digest mismatch: "
            f"{payload_digest} != {EXPECTED_PAYLOAD_LIST_SHA256}"
        )

    selected_cases: list[Dict[str, Any]] = []
    selected_hash_set = set(selected_hashes)
    for row in selected:
        selected_row = copy.deepcopy(row)
        selected_row["selection_reason"] = (
            "ALL_REDUCED_Q_GT_8"
            if int(row["reduced_q"]) > 8
            else "LEXICOGRAPHIC_MINIMUM_STRUCTURAL_SHA_PER_STRATUM"
        )
        selected_cases.append(selected_row)
    selected_cases.sort(
        key=lambda row: row["structural_key_sha256"]
    )
    if {
        row["structural_key_sha256"] for row in selected_cases
    } != selected_hash_set:
        raise AssertionError("selected-case materialization drifted from digest input")

    code_sources = {
        "exporter": pathlib.Path(__file__),
        "solver_probe": SOLVER_PROBE_SOURCE,
        "optimized_backend": OPTIMIZED_BACKEND_SOURCE,
    }
    for name, path in code_sources.items():
        if not path.is_file():
            raise CrossSamplePanelError(f"code binding is missing ({name}): {path}")

    checks = {
        "historical_run_success_marker_present": True,
        "origin_manifest_sha256_matches": True,
        "two_output_manifest_sha256_match": len(source_manifests) == 2,
        "twelve_source_files_bound_and_verified": len(data_sources) == 12,
        "missing_pattern_joins_zero": all(
            summary["missing_pattern_joins"] == 0
            for summary in sample_summaries.values()
        ),
        "historical_incomplete_counts_match": all(
            sample_summaries[sample]["historical_incomplete_units"]
            == EXPECTED_SAMPLE_COUNTS[sample]["incomplete_units"]
            for sample in SAMPLES
        ),
        "sample_structural_key_counts_match": all(
            sample_summaries[sample]["distinct_structural_keys"]
            == EXPECTED_SAMPLE_COUNTS[sample]["distinct_structural_keys"]
            for sample in SAMPLES
        ),
        "union_5506_keys": len(union_hashes) == EXPECTED_UNION_STRUCTURAL_KEYS,
        "all_48_q_gt_8_keys_selected": (
            selection_summary["q_gt_8_memberships"] == EXPECTED_Q_GT_8_KEYS
        ),
        "selected_208_keys": len(selected_cases) == EXPECTED_SELECTED_KEYS,
        "selected_sample_counts_86_122": (
            selection_summary["selected_by_sample"]
            == {
                sample: EXPECTED_SAMPLE_COUNTS[sample]["selected_keys"]
                for sample in SAMPLES
            }
        ),
        "key_list_digest_matches": key_digest == EXPECTED_KEY_LIST_SHA256,
        "payload_list_digest_matches": (
            payload_digest == EXPECTED_PAYLOAD_LIST_SHA256
        ),
        "solver_timing_not_invoked": True,
        "canonical_pipeline_not_modified_or_invoked": True,
    }
    if not all(checks.values()):
        raise CrossSamplePanelError("one or more authority checks failed")

    document: Dict[str, Any] = {
        "schema_name": (
            "intersubmod.historical_v5_cross_sample_solver_stress_panel_manifest"
        ),
        "schema_version": "1.0.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "authority": {
            "status": "AUTHORITATIVE",
            "claim_scope": "HISTORICAL_V5_SOLVER_CORE_ONLY",
            "write_policy": "WRITE_ONCE_REFUSE_EXISTING_OUTPUT",
        },
        "task_type": "B_COMPREHENSIVE_VALIDATION_COMPONENT",
        "research_goals": ["G3", "G4", "G5"],
        "scope": {
            "samples": list(SAMPLES),
            "historical_source": str(run_root.resolve()),
            "input_semantics": (
                "canonical-v5 subread_groups_by_hp; MINREAD>=3; "
                "historical primary incomplete units"
            ),
            "selection_rule": (
                "all reduced q>8 keys plus lexicographically minimum structural "
                "SHA per sample×k×effective_m×reduced_q×historical_cap_class"
            ),
            "not_claimed": [
                "M2 read-linked component or phase-set route equivalence",
                "cross-sample end-to-end completion rate",
                "solver runtime or speedup",
                "canonical or production backend promotion",
                "biological clone confirmation",
            ],
        },
        "bindings": {
            "origin_manifest": origin_binding,
            "sample_output_manifests": source_manifests,
            "data_source_file_count": len(data_sources),
            "data_sources": dict(sorted(data_sources.items())),
            "code_sources": {
                name: {
                    "path": str(path.resolve()),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for name, path in code_sources.items()
            },
        },
        "source_census": {
            "samples": sample_summaries,
            "distinct_structural_keys_union": len(union_hashes),
            "cross_sample_duplicate_structural_keys": (
                sum(
                    summary["distinct_structural_keys"]
                    for summary in sample_summaries.values()
                )
                - len(union_hashes)
            ),
        },
        "selection": {
            **selection_summary,
            "selected_global_distinct_keys": len(selected_hashes),
            "key_list_preimage_bytes": len(key_preimage),
            "key_list_sha256": key_digest,
            "payload_list_preimage_bytes": len(payload_preimage),
            "payload_list_sha256": payload_digest,
            "selected_structural_keys": selected_hashes,
        },
        "cases": selected_cases,
        "checks": checks,
        "integrity": {
            "manifest_content_hash_schema": (
                "sha256(canonical JSON document without "
                "integrity.manifest_content_sha256)"
            ),
        },
    }
    document["integrity"]["manifest_content_sha256"] = manifest_content_sha256(
        document
    )
    return document


def write_once_manifest(
    document: Mapping[str, Any],
    output_dir: pathlib.Path,
) -> None:
    if output_dir.exists():
        raise CrossSamplePanelError(
            f"write-once output already exists: {output_dir}"
        )
    output_dir.mkdir(parents=True, exist_ok=False)
    manifest_path = output_dir / "manifest.json"
    sidecar_path = output_dir / "manifest.json.sha256"
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
    manifest_sha = hashlib.sha256(encoded).hexdigest()
    with sidecar_path.open("x", encoding="utf-8") as handle:
        handle.write(f"{manifest_sha}  manifest.json\n")


def summary(
    document: Mapping[str, Any],
    *,
    status: str,
    output_dir: pathlib.Path,
) -> Dict[str, Any]:
    return {
        "status": status,
        "claim_scope": document["authority"]["claim_scope"],
        "output_dir": str(output_dir.resolve()),
        "source_files_verified": document["bindings"]["data_source_file_count"],
        "historical_incomplete_units": {
            sample: document["source_census"]["samples"][sample][
                "historical_incomplete_units"
            ]
            for sample in SAMPLES
        },
        "distinct_structural_keys": {
            sample: document["source_census"]["samples"][sample][
                "distinct_structural_keys"
            ]
            for sample in SAMPLES
        },
        "union_structural_keys": document["source_census"][
            "distinct_structural_keys_union"
        ],
        "q_gt_8_selected": document["selection"]["q_gt_8_memberships"],
        "selected_by_sample": document["selection"]["selected_by_sample"],
        "selected_global_keys": document["selection"][
            "selected_global_distinct_keys"
        ],
        "key_list_sha256": document["selection"]["key_list_sha256"],
        "payload_list_sha256": document["selection"]["payload_list_sha256"],
        "manifest_content_sha256": document["integrity"][
            "manifest_content_sha256"
        ],
        "all_checks_pass": all(document["checks"].values()),
        "solver_timing_invoked": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=pathlib.Path, default=RUN_ROOT)
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=DEFAULT_OUTPUT_DIR,
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    try:
        document = build_manifest(args.run_root)
        if not args.dry_run:
            write_once_manifest(document, args.output_dir)
    except (
        CrossSamplePanelError,
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
                document,
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
