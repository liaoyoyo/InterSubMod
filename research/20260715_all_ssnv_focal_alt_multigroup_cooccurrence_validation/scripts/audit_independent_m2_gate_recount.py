#!/usr/bin/env python3
"""Logic-independent M2 recount using the screen TSV and stable assignments."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import sys
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping

from scipy.stats import chi2, f, ncf, ncx2


SCHEMA_NAME = "intersubmod.independent_m2_gate_recount"
SCHEMA_VERSION = "2.0.0"
ASSIGNMENT_SCHEMA = "intersubmod.all_ssnv_stable_multigroup_read_assignment"
GATE_CONTRACT = (
    "m2-measured-axis-v4_asymmetric-positive-confound-and-observed-categorical-levels"
)
P_MAX = 0.05
PERMUTATIONS = 499
P_FLOOR = 1.0 / (PERMUTATIONS + 1)
TARGET_POWER = 0.80
POWER_ALPHA = 0.05
MIN_GROUP_N = 5
MAX_GROUPS = 10
CATEGORY_LEVEL_CEILINGS = {"hp_exact": 7, "hp_family": 5, "strand": 2}
AXIS_SPECS = (
    ("hp_exact", "categorical", "v", 0.30),
    ("hp_family", "categorical", "v", 0.30),
    ("strand", "categorical", "v", 0.30),
    ("start", "continuous", "eta2", 0.14),
    ("end", "continuous", "eta2", 0.14),
    ("length", "continuous", "eta2", 0.14),
    ("mapq", "continuous", "eta2", 0.14),
    ("cpg_called", "continuous", "eta2", 0.14),
)
INDETERMINATE_STATUSES = frozenset({"HIGH_EFFECT_P_NOT_PASS", "LOW_POWER", "MISSING"})

SiteKey = tuple[str, str, int, str, str]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def open_text(path: Path) -> Any:
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8", newline=""
    )


def parse_bool(value: Any, *, field: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no"}:
        return False
    raise ValueError(f"Invalid boolean for {field}: {value!r}")


def finite_float(value: Any, *, field: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid numeric value for {field}: {value!r}") from error
    if not math.isfinite(parsed):
        raise ValueError(f"Non-finite numeric value for {field}: {value!r}")
    return parsed


def parse_cluster_sizes(value: Any) -> dict[str, int]:
    if isinstance(value, str):
        value = json.loads(value)
    if not isinstance(value, Mapping) or not value:
        raise ValueError("cluster_sizes must be a non-empty object")
    parsed = {str(label): int(count) for label, count in value.items()}
    if any(not label or count <= 0 for label, count in parsed.items()) or len(parsed) < 2:
        raise ValueError(f"Invalid cluster_sizes: {value!r}")
    return parsed


def hp_family(tag: Any) -> str:
    value = str(tag).strip()
    if value in {"1", "HP1", "1-1", "1-2"}:
        return "HP1-side"
    if value in {"2", "HP2", "2-1", "2-2"}:
        return "HP2-side"
    if value == "3":
        return "HP3-ambiguous"
    if value == "4":
        return "HP4-both"
    return "untagged"


def _categorical_power(n: int, groups: int, levels: int, threshold: float) -> float:
    degrees = (groups - 1) * (levels - 1)
    scale = min(groups - 1, levels - 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        critical = chi2.ppf(1.0 - POWER_ALPHA, degrees)
        return float(ncx2.sf(critical, degrees, n * threshold * threshold * scale))


def _continuous_power(n: int, groups: int, threshold: float) -> float:
    if n <= groups:
        return 0.0
    between = groups - 1
    within = n - groups
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        critical = f.ppf(1.0 - POWER_ALPHA, between, within)
        noncentrality = n * threshold / (1.0 - threshold)
        return float(ncf.sf(critical, between, within, noncentrality))


@lru_cache(maxsize=None)
def minimum_n(prefix: str, kind: str, groups: int, threshold: float) -> int:
    start = max(groups + 1, groups * MIN_GROUP_N)
    for n in range(start, 10_001):
        power = (
            _categorical_power(n, groups, CATEGORY_LEVEL_CEILINGS[prefix], threshold)
            if kind == "categorical"
            else _continuous_power(n, groups, threshold)
        )
        if math.isfinite(power) and power >= TARGET_POWER:
            return n
    raise ValueError(f"Power target unattainable: {prefix} groups={groups}")


def assignment_key(record: Mapping[str, Any]) -> SiteKey:
    posthoc = record.get("posthoc")
    if not isinstance(posthoc, Mapping):
        raise ValueError("Assignment lacks posthoc metadata")
    return (
        str(record.get("dataset") or record.get("sample")),
        str(record["chrom"]),
        int(record["pos"]),
        str(posthoc["ref"]),
        str(posthoc["alt"]),
    )


def screen_key(row: Mapping[str, Any]) -> SiteKey:
    return str(row["dataset"]), str(row["chrom"]), int(row["pos"]), str(row["ref"]), str(row["alt"])


def load_assignment_proofs(path: Path) -> dict[SiteKey, dict[str, Any]]:
    proofs: dict[SiteKey, dict[str, Any]] = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            if record.get("schema_name") != ASSIGNMENT_SCHEMA or record.get("schema_version") != "1.0.0":
                raise ValueError(f"Unexpected assignment schema at line {line_number}")
            key = assignment_key(record)
            if key in proofs:
                raise ValueError(f"Duplicate assignment key: {key}")
            core_reads = record.get("core_reads")
            if not isinstance(core_reads, list) or not core_reads:
                raise ValueError(f"Assignment has no core reads: {key}")
            for read in core_reads:
                if not isinstance(read, Mapping) or any(
                    field not in read or read[field] is None
                    for field in ("label", "latest_hp", "strand")
                ):
                    raise ValueError(f"Assignment core read lacks proof fields: {key}")
            groups = Counter(str(read["label"]) for read in core_reads)
            proofs[key] = {
                "group_sizes": dict(groups),
                "category_levels": {
                    "hp_exact": len({str(read["latest_hp"]) for read in core_reads}),
                    "hp_family": len({hp_family(read["latest_hp"]) for read in core_reads}),
                    "strand": len({str(read["strand"]) for read in core_reads}),
                },
            }
    return proofs


def classify_axis(
    row: Mapping[str, Any],
    *,
    prefix: str,
    kind: str,
    effect_suffix: str,
    threshold: float,
    group_sizes: Mapping[str, int],
    category_levels: Mapping[str, int],
) -> dict[str, Any]:
    effect_field = f"{prefix}_{effect_suffix}"
    p_field = f"{prefix}_p_perm"
    n_field = f"{prefix}_n"
    aligned_field = f"{prefix}_aligned"
    n = int(row[n_field])
    if n != sum(group_sizes.values()):
        raise ValueError(f"Axis/core count drift for {prefix}: {n} != {sum(group_sizes.values())}")
    declared_aligned = parse_bool(row[aligned_field], field=aligned_field)
    effect_missing = row[effect_field] in {None, ""}
    p_missing = row[p_field] in {None, ""}
    if effect_missing != p_missing:
        raise ValueError(f"Partial axis statistic for {prefix}")
    observed_levels = category_levels[prefix] if kind == "categorical" else None
    if effect_missing:
        if declared_aligned:
            raise ValueError(f"Missing statistic declared aligned for {prefix}")
        if kind == "categorical" and observed_levels == 1:
            return {"status": "CONSTANT", "aligned": False, "sufficient": True}
        raise ValueError(
            f"Missing statistic lacks constant-axis proof for {prefix}: levels={observed_levels}"
        )
    if kind == "categorical" and observed_levels == 1:
        raise ValueError(f"Statistic exists for constant categorical axis {prefix}")
    effect = finite_float(row[effect_field], field=effect_field)
    p_value = finite_float(row[p_field], field=p_field)
    if not 0.0 <= effect <= 1.0 or not 0.0 <= p_value <= 1.0:
        raise ValueError(f"Axis statistic outside [0,1] for {prefix}")
    grid = p_value * (PERMUTATIONS + 1)
    if p_value < P_FLOOR - 1e-12 or not math.isclose(grid, round(grid), abs_tol=1e-9):
        raise ValueError(f"Permutation p-grid drift for {prefix}: {p_value}")
    groups = len(group_sizes)
    min_required = minimum_n(prefix, kind, groups, threshold)
    power = (
        _categorical_power(n, groups, CATEGORY_LEVEL_CEILINGS[prefix], threshold)
        if kind == "categorical"
        else _continuous_power(n, groups, threshold)
    )
    sufficient = min(group_sizes.values()) >= MIN_GROUP_N and n >= min_required and power >= TARGET_POWER
    aligned = effect >= threshold and p_value < P_MAX
    if aligned != declared_aligned:
        raise ValueError(f"Declared/recomputed alignment drift for {prefix}")
    if aligned:
        status = "ALIGNED"
    elif effect >= threshold:
        status = "HIGH_EFFECT_P_NOT_PASS"
    elif not sufficient:
        status = "LOW_POWER"
    else:
        status = "LOW_EFFECT"
    return {
        "status": status,
        "aligned": aligned,
        "sufficient": sufficient,
        "aligned_below_negative_power": aligned and not sufficient,
    }


def classify_site(row: Mapping[str, Any], proof: Mapping[str, Any]) -> dict[str, Any]:
    groups = parse_cluster_sizes(row["cluster_sizes"])
    if groups != proof["group_sizes"]:
        raise ValueError(f"Screen/assignment coarse group-count drift: {groups} != {proof['group_sizes']}")
    if len(groups) > MAX_GROUPS:
        return {"bucket": "not_evaluable_group_count_gt10", "axes": {}}
    levels = proof["category_levels"]
    for prefix, ceiling in CATEGORY_LEVEL_CEILINGS.items():
        if not 1 <= int(levels[prefix]) <= ceiling:
            raise ValueError(f"Observed category-level ceiling drift for {prefix}: {levels[prefix]}")
    axes = {
        prefix: classify_axis(
            row,
            prefix=prefix,
            kind=kind,
            effect_suffix=effect_suffix,
            threshold=threshold,
            group_sizes=groups,
            category_levels=levels,
        )
        for prefix, kind, effect_suffix, threshold in AXIS_SPECS
    }
    aligned = {prefix for prefix, result in axes.items() if result["aligned"]}
    indeterminate = {
        prefix for prefix, result in axes.items() if result["status"] in INDETERMINATE_STATUSES
    }
    hp_confound = parse_bool(row["hp_axis_confound"], field="hp_axis_confound")
    technical = parse_bool(row["technical_axis_confound"], field="technical_axis_confound")
    expected_hp = bool(aligned.intersection({"hp_exact", "hp_family"}))
    expected_technical = bool(aligned.difference({"hp_exact", "hp_family"}))
    if hp_confound != expected_hp or technical != expected_technical:
        raise ValueError("High-level confound flag drift")
    stable = parse_bool(row["stable_null_multigroup"], field="stable_null_multigroup")
    residual = parse_bool(row["residual_unexplained_multigroup"], field="residual_unexplained_multigroup")
    robust = parse_bool(
        row["phase_anchored_robust_epigenetic_candidate"],
        field="phase_anchored_robust_epigenetic_candidate",
    )
    if residual != (stable and not hp_confound and not technical) or (robust and not residual):
        raise ValueError("High-level residual/robust flag drift")
    if indeterminate:
        bucket = "not_evaluable_axis_indeterminate"
    elif hp_confound or technical or not residual or not robust:
        bucket = "evaluable_ineligible"
    else:
        bucket = "eligible"
    return {"bucket": bucket, "axes": axes}


def recount(site_results: Path, assignments: Path) -> dict[str, Any]:
    proofs = load_assignment_proofs(assignments)
    expected_assignment_keys = set(proofs)
    observed_stable_keys: set[SiteKey] = set()
    counts: Counter[str] = Counter()
    per_dataset: defaultdict[str, Counter[str]] = defaultdict(Counter)
    group_distribution: Counter[int] = Counter()
    constant_axis_counts: Counter[str] = Counter()
    aligned_below_power_axis_counts: Counter[str] = Counter()
    aligned_below_power_evaluable_sites = 0
    examples_group_count_gt10: list[dict[str, Any]] = []
    with open_text(site_results) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "dataset", "chrom", "pos", "ref", "alt", "stable_null_multigroup",
            "cluster_sizes", "hp_axis_confound", "technical_axis_confound",
            "residual_unexplained_multigroup", "phase_anchored_robust_epigenetic_candidate",
        }
        required.update(
            field
            for prefix, _kind, suffix, _threshold in AXIS_SPECS
            for field in (f"{prefix}_{suffix}", f"{prefix}_p_perm", f"{prefix}_n", f"{prefix}_aligned")
        )
        missing = sorted(required.difference(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"Site TSV lacks required fields: {missing}")
        for row in reader:
            counts["all_rows"] += 1
            if not parse_bool(row["stable_null_multigroup"], field="stable_null_multigroup"):
                continue
            counts["m1_stable_rows"] += 1
            key = screen_key(row)
            if key in observed_stable_keys:
                raise ValueError(f"Duplicate stable screen key: {key}")
            observed_stable_keys.add(key)
            if key not in proofs:
                raise ValueError(f"Stable screen key lacks assignment proof: {key}")
            group_count = len(parse_cluster_sizes(row["cluster_sizes"]))
            group_distribution[group_count] += 1
            result = classify_site(row, proofs[key])
            bucket = result["bucket"]
            counts[bucket] += 1
            dataset = key[0]
            per_dataset[dataset]["m1_stable_rows"] += 1
            per_dataset[dataset][bucket] += 1
            site_aligned_below_power = False
            for prefix, axis in result["axes"].items():
                if axis["status"] == "CONSTANT":
                    constant_axis_counts[prefix] += 1
                if axis.get("aligned_below_negative_power"):
                    aligned_below_power_axis_counts[prefix] += 1
                    site_aligned_below_power = True
            if site_aligned_below_power and bucket in {"eligible", "evaluable_ineligible"}:
                aligned_below_power_evaluable_sites += 1
            if group_count > MAX_GROUPS:
                examples_group_count_gt10.append(
                    {"dataset": key[0], "chrom": key[1], "pos": key[2], "ref": key[3], "alt": key[4]}
                )
    if observed_stable_keys != expected_assignment_keys:
        raise ValueError(
            "Stable screen/assignment key-set drift: "
            f"screen={len(observed_stable_keys)} assignments={len(expected_assignment_keys)} "
            f"missing={sorted(expected_assignment_keys - observed_stable_keys)[:3]} "
            f"extra={sorted(observed_stable_keys - expected_assignment_keys)[:3]}"
        )
    outcome_total = sum(
        counts[name]
        for name in (
            "eligible", "evaluable_ineligible", "not_evaluable_axis_indeterminate",
            "not_evaluable_group_count_gt10",
        )
    )
    return {
        "counts": dict(sorted(counts.items())),
        "per_dataset": {
            dataset: dict(sorted(values.items())) for dataset, values in sorted(per_dataset.items())
        },
        "group_count_distribution": {str(k): v for k, v in sorted(group_distribution.items())},
        "constant_axis_assignment_proof_counts": dict(sorted(constant_axis_counts.items())),
        "aligned_below_negative_evaluability_power_axis_counts": dict(
            sorted(aligned_below_power_axis_counts.items())
        ),
        "n_evaluable_sites_with_aligned_axis_below_negative_evaluability_power": (
            aligned_below_power_evaluable_sites
        ),
        "examples_group_count_gt10": examples_group_count_gt10,
        "conservation": {
            "stable_screen_assignment_keys_exact": observed_stable_keys == expected_assignment_keys,
            "stable_screen_key_count": len(observed_stable_keys),
            "stable_assignment_key_count": len(expected_assignment_keys),
            "m1_outcome_categories_sum": outcome_total,
            "m1_outcome_categories_equal_m1_stable_rows": outcome_total == counts["m1_stable_rows"],
            "group_distribution_sum": sum(group_distribution.values()),
            "group_distribution_equals_m1_stable_rows": (
                sum(group_distribution.values()) == counts["m1_stable_rows"]
            ),
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-results", type=Path, required=True)
    parser.add_argument("--stable-assignments", type=Path, required=True)
    parser.add_argument("--claim-contract", type=Path, required=True)
    parser.add_argument("--production-gate-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-all-rows", type=int, default=469_849)
    parser.add_argument("--expected-m1-stable-rows", type=int, default=102_842)
    parser.add_argument("--expected-eligible", type=int, default=919)
    parser.add_argument("--expected-evaluable-ineligible", type=int, default=948)
    parser.add_argument("--expected-axis-indeterminate", type=int, default=100_974)
    parser.add_argument("--expected-group-count-gt10", type=int, default=1)
    parser.add_argument("--expected-aligned-below-power-evaluable", type=int, default=121)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paths = {
        "site_results": args.site_results.expanduser().resolve(),
        "stable_assignments": args.stable_assignments.expanduser().resolve(),
        "claim_contract": args.claim_contract.expanduser().resolve(),
        "production_gate_source_reference_only": args.production_gate_source.expanduser().resolve(),
    }
    for label, path in paths.items():
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"Missing {label}: {path}")
    output = args.output.expanduser().resolve()
    if os.path.lexists(output):
        raise FileExistsError(f"Refusing to overwrite existing output: {output}")
    observed = recount(paths["site_results"], paths["stable_assignments"])
    counts = observed["counts"]
    checks = {
        "all_rows_match_expected": counts.get("all_rows") == args.expected_all_rows,
        "m1_stable_rows_match_expected": counts.get("m1_stable_rows") == args.expected_m1_stable_rows,
        "eligible_matches_expected": counts.get("eligible") == args.expected_eligible,
        "evaluable_ineligible_matches_expected": (
            counts.get("evaluable_ineligible") == args.expected_evaluable_ineligible
        ),
        "axis_indeterminate_matches_expected": (
            counts.get("not_evaluable_axis_indeterminate") == args.expected_axis_indeterminate
        ),
        "group_count_gt10_matches_expected": (
            counts.get("not_evaluable_group_count_gt10") == args.expected_group_count_gt10
        ),
        "aligned_below_power_evaluable_matches_expected": (
            observed["n_evaluable_sites_with_aligned_axis_below_negative_evaluability_power"]
            == args.expected_aligned_below_power_evaluable
        ),
        "stable_screen_assignment_keys_exact": observed["conservation"][
            "stable_screen_assignment_keys_exact"
        ],
        "m1_outcome_categories_conserved": observed["conservation"][
            "m1_outcome_categories_equal_m1_stable_rows"
        ],
        "group_distribution_conserved": observed["conservation"][
            "group_distribution_equals_m1_stable_rows"
        ],
        "group_limit_examples_reconcile": (
            len(observed["examples_group_count_gt10"])
            == counts.get("not_evaluable_group_count_gt10", 0)
        ),
    }
    payload = {
        "schema_name": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "task_type": "B_comprehensive_validation",
        "command": [str(Path(sys.executable).resolve()), *sys.argv],
        "contract": GATE_CONTRACT,
        "logic_independence": {
            "production_gate_imported": False,
            "production_gate_functions_called": False,
            "screen_effect_and_p_values_reused_as_frozen_inputs": True,
            "assignment_categories_and_coarse_group_counts_recomputed": True,
        },
        "inputs": {label: artifact(path) for label, path in paths.items()},
        "code": {"independent_recount": artifact(Path(__file__).resolve())},
        "expected": {
            "all_rows": args.expected_all_rows,
            "m1_stable_rows": args.expected_m1_stable_rows,
            "eligible": args.expected_eligible,
            "evaluable_ineligible": args.expected_evaluable_ineligible,
            "not_evaluable_axis_indeterminate": args.expected_axis_indeterminate,
            "not_evaluable_group_count_gt10": args.expected_group_count_gt10,
            "aligned_below_power_evaluable": args.expected_aligned_below_power_evaluable,
        },
        **observed,
        "checks": checks,
        "interpretation": {
            "positive_aligned_confound_can_fail_even_below_negative_evaluability_power": True,
            "non_alignment_requires_adequate_planning_power_or_constant_axis_proof": True,
            "eligible_is_not_final_G1_G2_or_subclone_discovery": True,
            "not_evaluable_is_not_a_biological_negative": True,
        },
        "pass_semantics": "logic_independent_gate_recount_integrity_not_scientific_confirmation",
        "pass": all(checks.values()),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"output": str(output), "checks": checks, "pass": payload["pass"]}))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
