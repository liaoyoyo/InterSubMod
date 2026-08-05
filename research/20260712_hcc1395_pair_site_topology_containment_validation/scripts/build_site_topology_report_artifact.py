#!/usr/bin/env python3
"""Build the canonical portable-report artifact for HCC1395 site topology.

The report consumes the machine outputs from ``build_site_topology_containment.py``
and stops at the canonical ``artifact.json`` boundary.  The final HTML must be
packaged with the Data Analytics ``deliver_portable_artifact.mjs`` renderer; this
script deliberately does not implement a second HTML or chart runtime.

The fixed analytical population is 5,720 exact-coordinate, complete-both region
pairs.  The report keeps the full read-compatible candidate space separate from
the VAF exact-score top set, retains VAF ties, and labels a missing conditional
null as a blocking access issue.  It never upgrades technical compatibility to
biological clone-tree truth.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sqlite3
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence


TOPIC = Path("research/20260712_hcc1395_pair_site_topology_containment_validation")
DEFAULT_DATA_DIR = TOPIC / "data"
REPORT_DIR = Path(
    "docs/reports/in_progress/2026/07/"
    "20260712_HCC1395逐位點拓撲子結構跨來源驗證_01"
)
REPORT_STEM = "20260712_HCC1395逐位點拓撲子結構跨來源驗證_01"
HCC1395_KB_PATH = Path("/big8_disk/liaoyoyo2001/knowledge/02_samples/hcc1395.md")
ANNOTATION_CONTEXT_PATH = Path(
    "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/"
    "hcc1395_annotation_reproducibility.tsv"
)

FIXED_DENOMINATOR = 5_720
ENDPOINTS = ("read_full", "vaf_official")
MAPPINGS = ("ordered_HP", "whole_region_HP_swap_tolerant")
PREFIX_BY_MAPPING = {
    "ordered_HP": "ordered",
    "whole_region_HP_swap_tolerant": "swap_tolerant",
}
ENDPOINT_ZH = {
    "read_full": "Read full candidate space",
    "vaf_official": "VAF exact-score top set",
    "vaf_normalized": "VAF normalized sensitivity",
}
MAPPING_ZH = {
    "ordered_HP": "Ordered HP",
    "whole_region_HP_swap_tolerant": "Global HP-swap tolerant",
    "identity": "Ordered HP",
    "global_hp_swap": "Global HP swap only",
    "swap_tolerant": "Best of ordered / one global swap",
}
OUTCOMES = (
    "exact",
    "a_proper_subset_b",
    "b_proper_subset_a",
    "overlap",
    "disjoint",
    "not_evaluable",
)
OUTCOME_ZH = {
    "exact": "候選集完全相同",
    "a_proper_subset_b": "HCC 候選集 ⊂ DORADO",
    "b_proper_subset_a": "DORADO 候選集 ⊂ HCC",
    "overlap": "僅部分候選相容",
    "disjoint": "共享位點拓撲衝突",
    "not_evaluable": "不可判定",
}
ROBUST_CATEGORIES = frozenset(
    {"exact", "a_proper_subset_b", "b_proper_subset_a"}
)
SITE_RELATION_ZH = {
    "equal": "位點集相同",
    "a_proper_subset_b": "HCC 位點集 ⊂ DORADO",
    "b_proper_subset_a": "DORADO 位點集 ⊂ HCC",
    "partial_overlap": "部分重疊（非嵌套）",
    "disjoint": "無共享位點",
}
FIXED_OUTCOME_ZH = {
    "strict_full_exact": "完整 exact",
    "A_induced_substructure": "HCC 子結構",
    "B_induced_substructure": "DORADO 子結構",
    "resolution_A_more_specific": "HCC 候選較窄",
    "resolution_B_more_specific": "DORADO 候選較窄",
    "candidate_overlap": "候選部分交集",
    "conflict": "關係衝突",
    "shared_core_only": "shared core",
    "not_evaluable": "不可判定",
}
FIXED_OUTCOME_ORDER = tuple(FIXED_OUTCOME_ZH)
SITE_PAIR_OUTCOME_ZH = {
    "determined_same": "唯一同關係",
    "ambiguous_equal": "多值集相同",
    "one_sided_contained": "單側被包含",
    "ambiguous_overlap": "部分交集",
    "conflict": "衝突",
    "not_evaluable": "不可判定",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=REPORT_DIR)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else Path.open
    with opener(path, mode="rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"expected JSON object: {path}")
    return value


def as_int(value: Any) -> int:
    if value in (None, ""):
        raise ValueError("empty integer")
    return int(float(value))


def as_float(value: Any) -> float:
    if value in (None, ""):
        raise ValueError("empty float")
    return float(value)


def as_bool(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"1", "true", "yes", "pass"}


def safe_int(value: Any) -> int | None:
    try:
        return as_int(value)
    except (TypeError, ValueError):
        return None


def safe_float(value: Any) -> float | None:
    try:
        return as_float(value)
    except (TypeError, ValueError):
        return None


def pct(value: int | float, denominator: int | float) -> float | None:
    return float(value) / float(denominator) if denominator else None


def fmt_pct(value: float | None, digits: int = 2) -> str:
    return "—" if value is None else f"{100 * value:.{digits}f}%"


def fmt_n(value: int | float | None) -> str:
    if value is None:
        return "—"
    return f"{int(value):,}" if float(value).is_integer() else f"{value:.3f}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_columns(path: Path, rows: Sequence[dict[str, Any]], required: Iterable[str]) -> None:
    if not rows:
        raise RuntimeError(f"required input is empty: {path}")
    missing = sorted(set(required) - set(rows[0]))
    if missing:
        raise RuntimeError(f"{path} is missing required columns: {missing}")


def wilson_interval(numerator: int, denominator: int, z: float = 1.959963984540054) -> tuple[float | None, float | None]:
    if denominator <= 0:
        return None, None
    proportion = numerator / denominator
    scale = 1 + z * z / denominator
    center = (proportion + z * z / (2 * denominator)) / scale
    half = z * math.sqrt(
        proportion * (1 - proportion) / denominator + z * z / (4 * denominator * denominator)
    ) / scale
    return max(0.0, center - half), min(1.0, center + half)


def chrom_key(row: dict[str, str]) -> tuple[int, int, str]:
    chrom = row.get("chrom", "").removeprefix("chr")
    chrom_index = int(chrom) if chrom.isdigit() else 10_000
    region = row.get("region", "")
    try:
        start = int(region.split(":", 1)[1].split("-", 1)[0])
    except (IndexError, ValueError):
        start = 10**18
    return chrom_index, start, region


def endpoint_fields(endpoint: str, mapping: str) -> dict[str, str]:
    prefix = f"{PREFIX_BY_MAPPING[mapping]}_{endpoint}"
    return {
        key: f"{prefix}_{key}"
        for key in (
            "status",
            "reason",
            "category",
            "candidate_tree_product_a",
            "candidate_tree_product_b",
            "topology_set_size_a",
            "topology_set_size_b",
            "intersection_size",
            "union_size",
            "jaccard",
            "coverage_a",
            "coverage_b",
            "compatible",
            "unique_exact",
            "n_informative_components",
            "n_relation_pairs",
            "selected_mapping",
        )
    }


def endpoint_stats(regions: Sequence[dict[str, str]], endpoint: str, mapping: str) -> dict[str, Any]:
    fields = endpoint_fields(endpoint, mapping)
    evaluable = [row for row in regions if row[fields["status"]] == "evaluable"]
    categories = Counter(row[fields["category"]] for row in evaluable)
    unknown = sorted(set(categories) - set(OUTCOMES))
    if unknown:
        raise RuntimeError(f"unexpected {endpoint}/{mapping} categories: {unknown}")
    robust = sum(categories[key] for key in ROBUST_CATEGORIES)
    possible = robust + categories["overlap"]
    unique_exact = sum(as_bool(row[fields["unique_exact"]]) for row in evaluable)
    return {
        "endpoint": endpoint,
        "mapping": mapping,
        "population": len(regions),
        "evaluable": len(evaluable),
        "not_evaluable": len(regions) - len(evaluable),
        "categories": categories,
        "robust": robust,
        "possible": possible,
        "unique_exact": unique_exact,
        "evaluable_coverage": pct(len(evaluable), len(regions)),
        "robust_yield": pct(robust, len(regions)),
        "robust_conditional": pct(robust, len(evaluable)),
        "possible_yield": pct(possible, len(regions)),
        "possible_conditional": pct(possible, len(evaluable)),
        "conflict_yield": pct(categories["disjoint"], len(regions)),
    }


def metrics_index(rows: Sequence[dict[str, str]]) -> dict[tuple[str, str, str, str], dict[str, str]]:
    result: dict[tuple[str, str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["endpoint"], row["comparison"], row["population"], row["metric"])
        if key in result:
            raise RuntimeError(f"duplicate metrics row: {key}")
        result[key] = row
    return result


def validate_inputs(
    regions_path: Path,
    regions: Sequence[dict[str, str]],
    metrics_path: Path,
    metrics: Sequence[dict[str, str]],
    checks_path: Path,
    checks: Sequence[dict[str, str]],
    summary: dict[str, Any],
) -> list[dict[str, Any]]:
    required_region_columns = {
        "match_id",
        "chrom",
        "region",
        "hp_count_a",
        "hp_count_b",
        "primary_universe_k_a",
        "primary_universe_k_b",
        "primary_universe_shared_k",
        "primary_universe_site_relation",
    }
    for endpoint in ENDPOINTS:
        for mapping in MAPPINGS:
            required_region_columns.update(endpoint_fields(endpoint, mapping).values())
    require_columns(regions_path, regions, required_region_columns)
    require_columns(
        metrics_path,
        metrics,
        {"endpoint", "comparison", "population", "metric", "numerator", "denominator", "percent"},
    )
    require_columns(checks_path, checks, {"check", "pass", "observed", "expected"})

    report_checks: list[dict[str, Any]] = []

    def check(check_id: str, observed: Any, expected: Any, detail: str) -> None:
        passed = observed == expected
        report_checks.append(
            {
                "check": check_id,
                "pass": "PASS" if passed else "FAIL",
                "observed": observed,
                "expected": expected,
                "detail": detail,
            }
        )
        if not passed:
            raise RuntimeError(f"report QA failed {check_id}: {observed!r} != {expected!r}")

    check("fixed_denominator", len(regions), FIXED_DENOMINATOR, "One row per exact-coordinate complete-both region.")
    check(
        "unique_match_id",
        len({row["match_id"] for row in regions}),
        FIXED_DENOMINATOR,
        "No region-pair duplication.",
    )
    check(
        "upstream_checks_pass",
        sum(as_bool(row["pass"]) for row in checks),
        len(checks),
        "All producer checks must pass before report packaging.",
    )
    check(
        "summary_scope_denominator",
        safe_int(summary.get("source_counts", {}).get("exact_coordinate_complete_both")),
        FIXED_DENOMINATOR,
        "Summary and region TSV share the fixed population.",
    )

    index = metrics_index(metrics)
    for endpoint in ENDPOINTS:
        for mapping in MAPPINGS:
            stats = endpoint_stats(regions, endpoint, mapping)
            key_base = (endpoint, mapping, "all_exact_complete_both")
            expected_metrics = {
                "population_regions": stats["population"],
                "evaluable_regions": stats["evaluable"],
                "compatible_regions": stats["possible"],
                **{
                    f"category_{category}": stats["categories"][category]
                    for category in OUTCOMES[:-1]
                },
            }
            for metric, expected in expected_metrics.items():
                key = (*key_base, metric)
                if key not in index:
                    raise RuntimeError(f"metrics input lacks required row: {key}")
                observed = as_int(index[key]["numerator"])
                check(
                    f"metric_reconcile_{endpoint}_{mapping}_{metric}",
                    observed,
                    expected,
                    "Independent report derivation reconciles to producer metrics.",
                )
            outcome_sum = stats["evaluable"] + stats["not_evaluable"]
            check(
                f"outcome_sum_{endpoint}_{mapping}",
                outcome_sum,
                FIXED_DENOMINATOR,
                "Every layer/mapping returns one terminal region outcome.",
            )
            check(
                f"lower_not_above_upper_{endpoint}_{mapping}",
                stats["robust"] <= stats["possible"],
                True,
                "Exact/nested robust compatibility cannot exceed any-overlap compatibility.",
            )
    return report_checks


def build_outcome_rows(all_stats: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for stats in all_stats:
        for outcome in OUTCOMES:
            count = (
                stats["not_evaluable"]
                if outcome == "not_evaluable"
                else stats["categories"][outcome]
            )
            rows.append(
                {
                    "endpoint": stats["endpoint"],
                    "layer": ENDPOINT_ZH[stats["endpoint"]],
                    "mapping": stats["mapping"],
                    "mapping_label": MAPPING_ZH[stats["mapping"]],
                    "layer_mapping": f"{ENDPOINT_ZH[stats['endpoint']]} · {MAPPING_ZH[stats['mapping']]}",
                    "outcome": outcome,
                    "outcome_label": OUTCOME_ZH[outcome],
                    "n": count,
                    "fixed_denominator": FIXED_DENOMINATOR,
                    "share": count / FIXED_DENOMINATOR,
                }
            )
    return rows


def build_bounds_rows(all_stats: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for stats in all_stats:
        context = {
            "endpoint": stats["endpoint"],
            "layer": ENDPOINT_ZH[stats["endpoint"]],
            "mapping": stats["mapping"],
            "mapping_label": MAPPING_ZH[stats["mapping"]],
            "layer_mapping": f"{ENDPOINT_ZH[stats['endpoint']]} · {MAPPING_ZH[stats['mapping']]}",
            "population": stats["population"],
            "evaluable": stats["evaluable"],
            "evaluable_coverage": stats["evaluable_coverage"],
        }
        rows.extend(
            [
                {
                    **context,
                    "bound": "robust_lower",
                    "bound_label": "Robust exact/nested",
                    "n": stats["robust"],
                    "all_pair_yield": stats["robust_yield"],
                    "conditional_rate": stats["robust_conditional"],
                },
                {
                    **context,
                    "bound": "any_overlap_upper",
                    "bound_label": "Any-overlap optimistic",
                    "n": stats["possible"],
                    "all_pair_yield": stats["possible_yield"],
                    "conditional_rate": stats["possible_conditional"],
                },
            ]
        )
    return rows


def build_site_inventory(regions: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    counts = Counter(row["primary_universe_site_relation"] for row in regions)
    unknown = sorted(set(counts) - set(SITE_RELATION_ZH))
    if unknown:
        raise RuntimeError(f"unexpected site-set relation values: {unknown}")
    return [
        {
            "relation": relation,
            "relation_label": SITE_RELATION_ZH[relation],
            "n": counts[relation],
            "fixed_denominator": FIXED_DENOMINATOR,
            "share": counts[relation] / FIXED_DENOMINATOR,
            "directional_containment_eligible": relation
            in {"a_proper_subset_b", "b_proper_subset_a"},
        }
        for relation in SITE_RELATION_ZH
    ]


def build_caller_site_inventory(summary: dict[str, Any]) -> list[dict[str, Any]]:
    """Build the audited caller-level site inventory from the VCF postprocess."""

    counts = summary["site_inventory"]["caller_region_relations"]
    unknown = sorted(set(counts) - set(SITE_RELATION_ZH))
    if unknown:
        raise RuntimeError(f"unexpected caller site-set relation values: {unknown}")
    rows: list[dict[str, Any]] = []
    for relation in SITE_RELATION_ZH:
        count = as_int(counts.get(relation, 0))
        rows.append(
            {
                "relation": relation,
                "relation_label": SITE_RELATION_ZH[relation],
                "n": count,
                "fixed_denominator": FIXED_DENOMINATOR,
                "share": count / FIXED_DENOMINATOR,
                "directional_containment_eligible": (
                    "yes" if relation in {"a_proper_subset_b", "b_proper_subset_a"} else "no"
                ),
            }
        )
    return rows


def build_fixed_outcome_rows(
    compatibility_rows: Sequence[dict[str, str]],
) -> list[dict[str, Any]]:
    """Normalize mutually exclusive, whole-region swap-tolerant outcomes."""

    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for row in compatibility_rows:
        if row["mapping"] != "whole_region_HP_swap_tolerant":
            continue
        endpoint = row["layer"]
        outcome = row["outcome"]
        if endpoint not in ENDPOINT_ZH or outcome not in FIXED_OUTCOME_ZH:
            raise RuntimeError(f"unexpected fixed outcome row: {endpoint}/{outcome}")
        key = (endpoint, outcome)
        if key in seen:
            raise RuntimeError(f"duplicate fixed outcome row: {key}")
        seen.add(key)
        count = as_int(row["n"])
        denominator = as_int(row["denominator"])
        if denominator != FIXED_DENOMINATOR:
            raise RuntimeError(f"fixed outcome denominator drift: {key} -> {denominator}")
        rows.append(
            {
                "endpoint": endpoint,
                "layer": ENDPOINT_ZH[endpoint],
                "mapping": row["mapping"],
                "mapping_label": "Best of ordered / one global HP swap",
                "layer_mapping": f"{ENDPOINT_ZH[endpoint]} · fixed swap-tolerant",
                "outcome": outcome,
                "outcome_label": FIXED_OUTCOME_ZH[outcome],
                "n": count,
                "fixed_denominator": denominator,
                "share": count / denominator,
            }
        )
    expected = {
        (endpoint, outcome)
        for endpoint in ENDPOINT_ZH
        for outcome in FIXED_OUTCOME_ORDER
    }
    if seen != expected:
        raise RuntimeError(
            f"fixed outcome grid incomplete: missing={sorted(expected - seen)} extra={sorted(seen - expected)}"
        )
    return rows


def summarize_fixed_outcomes(
    fixed_rows: Sequence[dict[str, Any]],
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """Return endpoint headline metrics and an auditable structural summary."""

    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in fixed_rows:
        counts[row["endpoint"]][row["outcome"]] = as_int(row["n"])
    stats: dict[str, dict[str, Any]] = {}
    summary_rows: list[dict[str, Any]] = []
    metric_specs = (
        (
            "strict_full_exact",
            "完整位點＋候選拓撲集 exact",
            "Caller 與每個 mapped HP 位點集相同，且投影候選拓撲集相同。",
        ),
        (
            "real_induced_substructure",
            "方向性 induced substructure",
            "一側 HP-specific 位點集真包含於另一側，且共享位點候選關係同向相容。",
        ),
        (
            "strong_exact_or_induced",
            "嚴格 exact 或真誘導子結構",
            "主張結構一致時使用的嚴格合計；不含只因候選不確定度不同而形成的集合嵌套。",
        ),
        (
            "resolution_nesting",
            "同位點但候選解析度不同",
            "兩側位點相同，一側候選空間較窄；是 constraint-space nesting，不是額外／缺失位點的子樹。",
        ),
        (
            "weak_shared_core_or_overlap",
            "只支持 shared core／部分候選交集",
            "有共享關係但不足以稱完整結構或方向性誘導子結構。",
        ),
        ("conflict", "共享位點關係衝突", "投影候選拓撲集在共享位點上無交集。"),
        ("not_evaluable", "不可判定", "HP 對應、共享關係、recurrence 或 VAF availability 阻擋評估。"),
    )
    for endpoint in ENDPOINT_ZH:
        c = counts[endpoint]
        evaluable = FIXED_DENOMINATOR - c["not_evaluable"]
        induced = c["A_induced_substructure"] + c["B_induced_substructure"]
        strong = c["strict_full_exact"] + induced
        resolution = c["resolution_A_more_specific"] + c["resolution_B_more_specific"]
        weak = c["shared_core_only"] + c["candidate_overlap"]
        values = {
            "strict_full_exact": c["strict_full_exact"],
            "real_induced_substructure": induced,
            "strong_exact_or_induced": strong,
            "resolution_nesting": resolution,
            "weak_shared_core_or_overlap": weak,
            "conflict": c["conflict"],
            "not_evaluable": c["not_evaluable"],
        }
        stats[endpoint] = {
            "population": FIXED_DENOMINATOR,
            "evaluable": evaluable,
            "evaluable_coverage": evaluable / FIXED_DENOMINATOR,
            "strong": strong,
            "strong_yield": strong / FIXED_DENOMINATOR,
            "strong_conditional": strong / evaluable if evaluable else None,
            "strict_full_exact": c["strict_full_exact"],
            "induced": induced,
            "resolution": resolution,
            "weak": weak,
            "conflict": c["conflict"],
            "not_evaluable": c["not_evaluable"],
        }
        for metric, label, definition in metric_specs:
            value = values[metric]
            summary_rows.append(
                {
                    "endpoint": endpoint,
                    "layer": ENDPOINT_ZH[endpoint],
                    "mapping": "whole_region_HP_swap_tolerant",
                    "mapping_label": "Best of ordered / one global HP swap",
                    "metric": metric,
                    "metric_label": label,
                    "n": value,
                    "fixed_denominator": FIXED_DENOMINATOR,
                    "share": value / FIXED_DENOMINATOR,
                    "conditional_share": value / evaluable if evaluable else None,
                    "definition": definition,
                    "exclusive": "no" if metric == "strong_exact_or_induced" else "yes",
                }
            )
    return stats, summary_rows


def build_site_pair_summary(summary: dict[str, Any]) -> list[dict[str, Any]]:
    """Aggregate per-shared-site-pair relation outcomes without relabeling ambiguity as certainty."""

    rows: list[dict[str, Any]] = []
    for endpoint, counts in summary["site_pair_evidence"].items():
        total = sum(as_int(value) for value in counts.values())
        not_evaluable = as_int(counts["not_evaluable"])
        evaluable = total - not_evaluable
        for outcome in SITE_PAIR_OUTCOME_ZH:
            count = as_int(counts.get(outcome, 0))
            rows.append(
                {
                    "endpoint": endpoint,
                    "layer": ENDPOINT_ZH[endpoint],
                    "outcome": outcome,
                    "outcome_label": SITE_PAIR_OUTCOME_ZH[outcome],
                    "n": count,
                    "total_site_pairs": total,
                    "evaluable_site_pairs": evaluable,
                    "all_pair_share": count / total if total else None,
                    "conditional_share": (
                        count / evaluable if evaluable and outcome != "not_evaluable" else None
                    ),
                }
            )
    return rows


def build_postprocess_examples(rows: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    return [
        {
            "layer": ENDPOINT_ZH[row["layer"]],
            "outcome": row["outcome"],
            "outcome_label": FIXED_OUTCOME_ZH[row["outcome"]],
            "example_rank": as_int(row["example_rank"]),
            "chrom": row["chrom"],
            "region": row["region"],
            "site_set_relation": row["caller_site_set_relation"],
            "site_set_relation_label": SITE_RELATION_ZH[row["caller_site_set_relation"]],
            "shared_k": as_int(row["caller_shared_k"]),
            "hp_pair": f"{row['hp_count_a']} vs {row['hp_count_b']}",
            "selected_mapping": row["selected_mapping"],
            "candidate_relation": row["candidate_relation"],
            "candidate_tree_product_a": safe_int(row["candidate_tree_product_a"]),
            "candidate_tree_product_b": safe_int(row["candidate_tree_product_b"]),
            "topology_set_size_a": safe_int(row["topology_set_size_a"]),
            "topology_set_size_b": safe_int(row["topology_set_size_b"]),
            "outcome_reason": row["outcome_reason"],
        }
        for row in rows
    ]


def k_bin(row: dict[str, str], _fields: dict[str, str]) -> str:
    value = as_int(row["primary_universe_shared_k"])
    if value <= 1:
        return "0–1"
    if value <= 5:
        return str(value)
    return "≥6"


def hp_bin(row: dict[str, str], _fields: dict[str, str]) -> str:
    left, right = as_int(row["hp_count_a"]), as_int(row["hp_count_b"])
    return f"{left} vs {right}"


def candidate_bin(row: dict[str, str], fields: dict[str, str]) -> str:
    left = safe_int(row[fields["candidate_tree_product_a"]])
    right = safe_int(row[fields["candidate_tree_product_b"]])
    if left is None or right is None:
        return "not available"
    value = left * right
    if value == 1:
        return "1"
    if value <= 5:
        return "2–5"
    if value <= 20:
        return "6–20"
    return ">20"


def build_complexity_rows(regions: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    specs: list[tuple[str, str, Callable[[dict[str, str], dict[str, str]], str], list[str]]] = [
        ("shared_k", "Shared-site k", k_bin, ["0–1", "2", "3", "4", "5", "≥6"]),
        ("hp_pair", "Primary HP components", hp_bin, []),
        (
            "candidate_pairs",
            "Candidate-tree pair count",
            candidate_bin,
            ["1", "2–5", "6–20", ">20", "not available"],
        ),
    ]
    rows: list[dict[str, Any]] = []
    for endpoint in ENDPOINTS:
        for mapping in MAPPINGS:
            fields = endpoint_fields(endpoint, mapping)
            for stratum_type, stratum_label, classifier, fixed_order in specs:
                grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
                for row in regions:
                    grouped[classifier(row, fields)].append(row)
                order = fixed_order or sorted(grouped)
                for order_index, value in enumerate(order):
                    selected = grouped.get(value, [])
                    if not selected:
                        continue
                    evaluable = [row for row in selected if row[fields["status"]] == "evaluable"]
                    robust = sum(row[fields["category"]] in ROBUST_CATEGORIES for row in evaluable)
                    conflict = sum(row[fields["category"]] == "disjoint" for row in evaluable)
                    possible_only = sum(row[fields["category"]] == "overlap" for row in evaluable)
                    low, high = wilson_interval(robust, len(evaluable))
                    rows.append(
                        {
                            "endpoint": endpoint,
                            "layer": ENDPOINT_ZH[endpoint],
                            "mapping": mapping,
                            "mapping_label": MAPPING_ZH[mapping],
                            "stratum_type": stratum_type,
                            "stratum_label": stratum_label,
                            "stratum_value": value,
                            "stratum_order": order_index,
                            "population": len(selected),
                            "evaluable": len(evaluable),
                            "evaluable_coverage": pct(len(evaluable), len(selected)),
                            "robust": robust,
                            "robust_yield": pct(robust, len(selected)),
                            "robust_conditional": pct(robust, len(evaluable)),
                            "ci95_low": low,
                            "ci95_high": high,
                            "possible_only": possible_only,
                            "conflict": conflict,
                        }
                    )
    return rows


def build_structural_summary(
    regions: Sequence[dict[str, str]], all_stats: Sequence[dict[str, Any]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    stats_index = {(row["endpoint"], row["mapping"]): row for row in all_stats}
    definitions = {
        "same_site_unique_exact": "位點集相同，且雙側投影後都只有一個拓撲 signature。",
        "same_site_candidate_set_exact": "位點集相同，且完整候選 topology-signature 集相同；可能仍有多顆候選。",
        "directional_induced_substructure": "一側位點集是另一側真子集，且共享位點上的候選空間屬 exact/nested。",
        "shared_core_only": "位點集只部分重疊但共享投影非衝突；不算方向性子結構。",
        "candidate_space_overlap_only": "完整候選集只有一部分交集；屬樂觀 possible-only。",
        "conflict": "在可判定的共享位點上，候選 topology-signature 集無交集。",
        "not_evaluable": "HP 對應、共享關係或 VAF top set 不足，fail closed。",
    }
    for endpoint in ENDPOINTS:
        for mapping in MAPPINGS:
            fields = endpoint_fields(endpoint, mapping)
            evaluable = [row for row in regions if row[fields["status"]] == "evaluable"]
            counts = {
                "same_site_unique_exact": sum(
                    row["primary_universe_site_relation"] == "equal"
                    and as_bool(row[fields["unique_exact"]])
                    for row in evaluable
                ),
                "same_site_candidate_set_exact": sum(
                    row["primary_universe_site_relation"] == "equal"
                    and row[fields["category"]] == "exact"
                    for row in evaluable
                ),
                "directional_induced_substructure": sum(
                    row["primary_universe_site_relation"]
                    in {"a_proper_subset_b", "b_proper_subset_a"}
                    and row[fields["category"]] in ROBUST_CATEGORIES
                    for row in evaluable
                ),
                "shared_core_only": sum(
                    row["primary_universe_site_relation"] == "partial_overlap"
                    and row[fields["category"]] != "disjoint"
                    for row in evaluable
                ),
                "candidate_space_overlap_only": sum(
                    row[fields["category"]] == "overlap" for row in evaluable
                ),
                "conflict": sum(row[fields["category"]] == "disjoint" for row in evaluable),
                "not_evaluable": stats_index[(endpoint, mapping)]["not_evaluable"],
            }
            for metric, count in counts.items():
                rows.append(
                    {
                        "endpoint": endpoint,
                        "layer": ENDPOINT_ZH[endpoint],
                        "mapping": mapping,
                        "mapping_label": MAPPING_ZH[mapping],
                        "metric": metric,
                        "metric_label": metric.replace("_", " "),
                        "n": count,
                        "fixed_denominator": FIXED_DENOMINATOR,
                        "share": count / FIXED_DENOMINATOR,
                        "definition": definitions[metric],
                        "exclusive": metric
                        in {"candidate_space_overlap_only", "conflict", "not_evaluable"},
                    }
                )
    return rows


def build_examples(regions: Sequence[dict[str, str]]) -> list[dict[str, Any]]:
    ordered = sorted(regions, key=chrom_key)
    read = endpoint_fields("read_full", "ordered_HP")
    vaf = endpoint_fields("vaf_official", "ordered_HP")

    selectors: list[tuple[str, str, Callable[[dict[str, str]], bool]]] = [
        (
            "unique_exact",
            "位點與投影拓撲都唯一且完全相同",
            lambda row: row[read["status"]] == "evaluable"
            and row["primary_universe_site_relation"] == "equal"
            and as_bool(row[read["unique_exact"]]),
        ),
        (
            "hcc_induced_substructure",
            "HCC 較小位點集與 DORADO 共享投影相容",
            lambda row: row[read["status"]] == "evaluable"
            and row["primary_universe_site_relation"] == "a_proper_subset_b"
            and row[read["category"]] in ROBUST_CATEGORIES,
        ),
        (
            "dorado_induced_substructure",
            "DORADO 較小位點集與 HCC 共享投影相容",
            lambda row: row[read["status"]] == "evaluable"
            and row["primary_universe_site_relation"] == "b_proper_subset_a"
            and row[read["category"]] in ROBUST_CATEGORIES,
        ),
        (
            "candidate_overlap_only",
            "僅一部分 read-compatible 候選相容",
            lambda row: row[read["status"]] == "evaluable"
            and row[read["category"]] == "overlap",
        ),
        (
            "hard_conflict",
            "共享位點關係無任何候選交集",
            lambda row: row[read["status"]] == "evaluable"
            and row[read["category"]] == "disjoint",
        ),
        (
            "vaf_rescue",
            "Read 層只 partial/conflict，VAF top 縮到 exact/nested",
            lambda row: row[read["status"]] == "evaluable"
            and row[vaf["status"]] == "evaluable"
            and row[read["category"]] in {"overlap", "disjoint"}
            and row[vaf["category"]] in ROBUST_CATEGORIES,
        ),
        (
            "vaf_divergence",
            "Read 層 robust，但 VAF top 在共享位點上衝突",
            lambda row: row[read["status"]] == "evaluable"
            and row[vaf["status"]] == "evaluable"
            and row[read["category"]] in ROBUST_CATEGORIES
            and row[vaf["category"]] == "disjoint",
        ),
        (
            "not_evaluable",
            "Ordered HP 下 fail-closed 不可判定",
            lambda row: row[read["status"]] != "evaluable"
            or row[vaf["status"]] != "evaluable",
        ),
    ]

    rows: list[dict[str, Any]] = []
    for example_type, label, predicate in selectors:
        selected = next((row for row in ordered if predicate(row)), None)
        if selected is None:
            continue
        rows.append(
            {
                "example_type": example_type,
                "selection_rule": label,
                "region": selected["region"],
                "chrom": selected["chrom"],
                "site_set_relation": selected["primary_universe_site_relation"],
                "site_set_relation_label": SITE_RELATION_ZH[selected["primary_universe_site_relation"]],
                "k_a": as_int(selected["primary_universe_k_a"]),
                "k_b": as_int(selected["primary_universe_k_b"]),
                "shared_k": as_int(selected["primary_universe_shared_k"]),
                "hp_pair": f"{selected['hp_count_a']} vs {selected['hp_count_b']}",
                "read_status": selected[read["status"]],
                "read_outcome": selected[read["category"]] or selected[read["reason"]],
                "read_q_a": safe_int(selected[read["topology_set_size_a"]]),
                "read_q_b": safe_int(selected[read["topology_set_size_b"]]),
                "vaf_status": selected[vaf["status"]],
                "vaf_outcome": selected[vaf["category"]] or selected[vaf["reason"]],
                "vaf_q_a": safe_int(selected[vaf["topology_set_size_a"]]),
                "vaf_q_b": safe_int(selected[vaf["topology_set_size_b"]]),
            }
        )
    return rows


def build_definitions() -> list[dict[str, str]]:
    return [
        {
            "term": "Fixed population",
            "definition": "5,720 chr1–22 exact-coordinate region pairs where both sources have complete primary topology outputs.",
            "boundary": "All-pair yield always uses 5,720; conditional rates use evaluable regions only.",
        },
        {
            "term": "Read full candidate space",
            "definition": "All minimal trees compatible with the observed full/partial same-HP read patterns.",
            "boundary": "It is a constraint set, not an absolute or biological truth tree.",
        },
        {
            "term": "VAF exact-score top set",
            "definition": "Candidates maximizing the exact Fraction same-HP read-AF ordering score; all exact ties remain.",
            "boundary": "Uses the same reads, is uncorrected for CN/purity/multiplicity, and is not an independent validation.",
        },
        {
            "term": "Shared-site projection",
            "definition": "For each common genomic mutation pair, retain forward ancestry, reverse ancestry, or parallel relation; private events are omitted.",
            "boundary": "Recurrence at a shared site fails closed; arbitrary unlabeled graph embedding is forbidden.",
        },
        {
            "term": "Directional induced substructure",
            "definition": "One complete site set is a strict subset of the other and the projected candidate spaces on shared sites are exact or nested.",
            "boundary": "Partial non-nested overlap is only a shared core, not a subtree claim.",
        },
        {
            "term": "Candidate-set exact / subset / overlap / disjoint",
            "definition": "Set relation between unique projected topology signatures from the two sources.",
            "boundary": "Subset describes uncertainty-space nesting; it does not by itself prove a biological ancestor tree.",
        },
        {
            "term": "Ordered HP",
            "definition": "HP family labels must match across sources.",
            "boundary": "Primary endpoint; an HP component-count mismatch remains not evaluable.",
        },
        {
            "term": "Global HP-swap sensitivity",
            "definition": "Compare identity mapping with one whole-region HP1↔HP2 mapping and retain the better evaluable mapping.",
            "boundary": "No per-component cherry-picking is allowed.",
        },
    ]


def build_claim_ceiling(null_available: bool) -> list[dict[str, str]]:
    null_boundary = (
        "Projected candidate-set exact/nesting can be compared with a within-region B-site-label permutation null."
        if null_available
        else "The within-region B-site-label permutation null is absent, so chance enrichment is not established."
    )
    return [
        {
            "evidence": "Read-full projected candidate relation",
            "allowed": "The two technical sources reproduce exact or nested read-compatible constraints on shared sSNVs.",
            "not_allowed": "The same true clone tree was recovered.",
        },
        {
            "evidence": "Directional site-set containment",
            "allowed": "The smaller observed site structure is compatible with the larger source after private-site contraction.",
            "not_allowed": "Both complete trees are identical.",
        },
        {
            "evidence": "VAF top-set relation",
            "allowed": "The same-HP read-AF heuristic prioritizes compatible candidate signatures.",
            "not_allowed": "VAF confirms the true ancestor direction.",
        },
        {
            "evidence": "Conditional-null gate",
            "allowed": null_boundary,
            "not_allowed": "Reuse its p-value for fixed strict-exact+true-induced, or claim biological accuracy without external truth.",
        },
        {
            "evidence": "Same-cell-line cross-source pair",
            "allowed": "Partial HCC1395 technical reproducibility under the frozen pipeline.",
            "not_allowed": "Independent biological replication or generalization to patients/cell lines.",
        },
    ]


def normalize_null(path: Path) -> list[dict[str, Any]]:
    rows = read_tsv(path)
    required = {
        "endpoint",
        "comparison",
        "population",
        "k_bin",
        "metric",
        "denominator",
        "observed_n",
        "observed_rate",
        "null_mean_rate",
        "null_q025_rate",
        "null_q975_rate",
        "observed_minus_null_mean",
        "empirical_p_ge",
        "block_bootstrap_excess_ci95_low",
        "block_bootstrap_excess_ci95_high",
        "loco_positive_chromosomes",
        "loco_chromosomes",
        "permutations",
        "null_method",
        "vaf_condition",
        "k_definition",
    }
    require_columns(path, rows, required)
    normalized: list[dict[str, Any]] = []
    for row in rows:
        if row["k_bin"] != "all":
            continue
        normalized.append(
            {
                "endpoint": row["endpoint"],
                "layer": ENDPOINT_ZH.get(row["endpoint"], row["endpoint"]),
                "comparison": row["comparison"],
                "mapping_label": MAPPING_ZH.get(row["comparison"], row["comparison"]),
                "population": row["population"],
                "metric": row["metric"],
                "denominator": as_int(row["denominator"]),
                "observed_n": as_int(row["observed_n"]),
                "observed": as_float(row["observed_rate"]),
                "null_mean": as_float(row["null_mean_rate"]),
                "null_q025": as_float(row["null_q025_rate"]),
                "null_q975": as_float(row["null_q975_rate"]),
                "delta": as_float(row["observed_minus_null_mean"]),
                "p_empirical": as_float(row["empirical_p_ge"]),
                "block_ci_low": as_float(row["block_bootstrap_excess_ci95_low"]),
                "block_ci_high": as_float(row["block_bootstrap_excess_ci95_high"]),
                "loco_positive_chromosomes": as_int(row["loco_positive_chromosomes"]),
                "loco_chromosomes": as_int(row["loco_chromosomes"]),
                "permutations": as_int(row["permutations"]),
                "null_method": row["null_method"],
                "vaf_condition": row["vaf_condition"],
                "k_definition": row["k_definition"],
            }
        )
    return normalized


def source_spec(
    source_id: str,
    label: str,
    path: Path,
    description: str,
    generated_at: str,
    *,
    filters: Sequence[str] = (),
    metric_definitions: Sequence[str] = (),
) -> dict[str, Any]:
    return {
        "id": source_id,
        "label": label,
        "path": path.as_posix(),
        "query": {
            "engine": "python",
            "language": "python",
            "description": description,
            "executed_at": generated_at,
            "tables_used": [path.as_posix()],
            "filters": list(filters),
            "metric_definitions": list(metric_definitions),
        },
    }


def sql_literal(value: Any) -> str:
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (int, float)):
        if isinstance(value, float) and not math.isfinite(value):
            return "NULL"
        return repr(value)
    return "'" + str(value).replace("'", "''") + "'"


def freeze_dataset_with_sql(
    dataset_name: str, rows: list[dict[str, Any]]
) -> tuple[str, list[dict[str, Any]]]:
    """Materialize reviewed snapshot rows through a literal SQLite query."""

    if not rows:
        raise RuntimeError(f"cannot freeze empty snapshot dataset: {dataset_name}")
    fields = list(rows[0])
    if any(list(row) != fields for row in rows):
        raise RuntimeError(f"snapshot dataset has unstable columns: {dataset_name}")
    identifiers = ", ".join(f'"{field}"' for field in fields)
    values = ",\n    ".join(
        "(" + ", ".join(sql_literal(row[field]) for field in fields) + ")" for row in rows
    )
    query = (
        f'WITH "{dataset_name}" ({identifiers}) AS (\n'
        f"  VALUES\n    {values}\n"
        f')\nSELECT {identifiers} FROM "{dataset_name}"'
    )
    connection = sqlite3.connect(":memory:")
    try:
        cursor = connection.execute(query)
        frozen = [dict(zip(fields, result)) for result in cursor.fetchall()]
    finally:
        connection.close()
    if len(frozen) != len(rows):
        raise RuntimeError(f"snapshot SQL row-count mismatch: {dataset_name}")
    return query, frozen


def columns(*specs: tuple[str, str, dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"field": field, "label": label, **options} for field, label, options in specs]


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir
    output_dir = args.output_dir
    paths = {
        "regions": data_dir / "hcc1395_site_topology_regions.tsv",
        "metrics": data_dir / "hcc1395_site_topology_metrics.tsv",
        "checks": data_dir / "hcc1395_site_topology_checks.tsv",
        "summary": data_dir / "hcc1395_site_topology_summary.json",
        "signatures": data_dir / "hcc1395_site_topology_signature_sets.jsonl.gz",
        "pair_outcomes": data_dir / "hcc1395_site_topology_pair_outcomes.tsv",
        "compatibility": data_dir / "hcc1395_topology_compatibility_metrics.tsv",
        "outcome_complexity": data_dir / "hcc1395_topology_outcome_complexity_strata.tsv",
        "examples": data_dir / "hcc1395_topology_examples.tsv",
        "containment_summary": data_dir / "hcc1395_topology_containment_summary.json",
        "containment_checks": data_dir / "hcc1395_topology_containment_checks.tsv",
        "site_relation_evidence": data_dir / "hcc1395_site_relation_evidence.tsv.gz",
        "allele_identity": data_dir / "hcc1395_site_allele_identity.tsv.gz",
        "null_checks": data_dir / "hcc1395_topology_null_checks.tsv",
    }
    missing = [path for path in paths.values() if not path.exists()]
    if not HCC1395_KB_PATH.exists():
        missing.append(HCC1395_KB_PATH)
    if not ANNOTATION_CONTEXT_PATH.exists():
        missing.append(ANNOTATION_CONTEXT_PATH)
    if missing:
        raise FileNotFoundError(
            "required containment outputs are not ready:\n" + "\n".join(str(path) for path in missing)
        )

    regions = read_tsv(paths["regions"])
    metrics = read_tsv(paths["metrics"])
    checks = read_tsv(paths["checks"])
    summary = load_json(paths["summary"])
    pair_outcomes = read_tsv(paths["pair_outcomes"])
    compatibility_source = read_tsv(paths["compatibility"])
    outcome_complexity_source = read_tsv(paths["outcome_complexity"])
    example_source = read_tsv(paths["examples"])
    containment_summary = load_json(paths["containment_summary"])
    containment_checks = read_tsv(paths["containment_checks"])
    site_relation_source = read_tsv(paths["site_relation_evidence"])
    allele_identity_source = read_tsv(paths["allele_identity"])
    null_checks_source = read_tsv(paths["null_checks"])
    require_columns(
        paths["pair_outcomes"],
        pair_outcomes,
        {
            "match_id",
            "chrom",
            "region",
            "fixed_denominator",
            "caller_shared_k",
            "k_bin",
            "caller_site_set_relation",
            "allele_identity_pass",
            "read_full_selected_mapping",
            "read_full_outcome",
            "vaf_official_selected_mapping",
            "vaf_official_outcome",
        },
    )
    require_columns(
        paths["compatibility"],
        compatibility_source,
        {"layer", "mapping", "outcome", "n", "denominator", "share"},
    )
    require_columns(
        paths["outcome_complexity"],
        outcome_complexity_source,
        {"layer", "mapping", "k_bin", "outcome", "n", "denominator", "share"},
    )
    require_columns(
        paths["examples"],
        example_source,
        {
            "layer",
            "outcome",
            "example_rank",
            "chrom",
            "region",
            "caller_site_set_relation",
            "caller_shared_k",
            "selected_mapping",
            "candidate_relation",
            "outcome_reason",
        },
    )
    require_columns(
        paths["containment_checks"],
        containment_checks,
        {"check", "pass", "observed", "expected"},
    )
    require_columns(
        paths["site_relation_evidence"],
        site_relation_source,
        {
            "match_id",
            "chrom",
            "region",
            "left_position",
            "right_position",
            "read_full_outcome",
            "vaf_official_outcome",
        },
    )
    require_columns(
        paths["allele_identity"],
        allele_identity_source,
        {"match_id", "chrom", "region", "position", "allele_identity"},
    )
    require_columns(
        paths["null_checks"],
        null_checks_source,
        {"check", "status", "observed", "expected"},
    )
    annotation_source_rows = read_tsv(ANNOTATION_CONTEXT_PATH)
    require_columns(
        ANNOTATION_CONTEXT_PATH,
        annotation_source_rows,
        {
            "feature",
            "feature_label",
            "n_present",
            "present_agreement",
            "difference_present_minus_absent_pp",
            "permutation_null_q025_pp",
            "permutation_null_q975_pp",
            "permutation_p_two_sided",
        },
    )
    annotation_by_feature = {row["feature"]: row for row in annotation_source_rows}
    annotation_keys = (
        "cgc_body_any",
        "dgidb_approved_antineoplastic_body_any",
    )
    absent_annotation_keys = [key for key in annotation_keys if key not in annotation_by_feature]
    if absent_annotation_keys:
        raise RuntimeError(f"annotation context lacks required rows: {absent_annotation_keys}")
    annotation_context = [
        {
            "feature": annotation_by_feature[key]["feature_label"],
            "n_present": as_int(annotation_by_feature[key]["n_present"]),
            "present_agreement": as_float(annotation_by_feature[key]["present_agreement"]),
            "present_minus_absent_pp": as_float(
                annotation_by_feature[key]["difference_present_minus_absent_pp"]
            ),
            "null_q025_pp": as_float(annotation_by_feature[key]["permutation_null_q025_pp"]),
            "null_q975_pp": as_float(annotation_by_feature[key]["permutation_null_q975_pp"]),
            "p_two_sided": as_float(annotation_by_feature[key]["permutation_p_two_sided"]),
            "interpretation": "Context only; not topology truth and not a clinical-actionability result.",
        }
        for key in annotation_keys
    ]
    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")

    report_checks = validate_inputs(
        paths["regions"],
        regions,
        paths["metrics"],
        metrics,
        paths["checks"],
        checks,
        summary,
    )

    def post_check(check_id: str, observed: Any, expected: Any, detail: str) -> None:
        passed = observed == expected
        report_checks.append(
            {
                "check": check_id,
                "pass": "PASS" if passed else "FAIL",
                "observed": observed,
                "expected": expected,
                "detail": detail,
            }
        )
        if not passed:
            raise RuntimeError(
                f"report postprocess QA failed {check_id}: {observed!r} != {expected!r}"
            )

    post_check(
        "fixed_pair_outcome_denominator",
        len(pair_outcomes),
        FIXED_DENOMINATOR,
        "Corrected fixed mutually exclusive region outcomes.",
    )
    post_check(
        "containment_checks_pass",
        sum(as_bool(row["pass"]) for row in containment_checks),
        len(containment_checks),
        "Postprocess site inventory, allele identity, outcome, and evidence QA.",
    )
    post_check(
        "shared_allele_identity_rows",
        len(allele_identity_source),
        15_713,
        "One ledger row per shared caller sSNV locus.",
    )
    post_check(
        "shared_allele_identity_pass",
        sum(as_bool(row["allele_identity"]) for row in allele_identity_source),
        len(allele_identity_source),
        "CHROM/POS/REF/ALT identity on every shared locus.",
    )
    post_check(
        "site_pair_evidence_rows",
        len(site_relation_source),
        as_int(containment_summary["evidence_rows"]),
        "Shared-site pair relation evidence is complete.",
    )
    post_check(
        "null_checks_pass",
        sum(as_bool(row["status"]) for row in null_checks_source),
        len(null_checks_source),
        "Within-region label-permutation, block-bootstrap, and LOCO QA.",
    )
    # The core region TSV is retained for candidate-space diagnostics, while all
    # report headline claims use the corrected fixed whole-region outcome ledger.
    all_stats = [
        endpoint_stats(regions, endpoint, mapping)
        for endpoint in ENDPOINTS
        for mapping in MAPPINGS
    ]
    bounds_rows = build_bounds_rows(all_stats)
    complexity_rows = build_complexity_rows(regions)
    complexity_k = [
        row
        for row in complexity_rows
        if row["mapping"] == "ordered_HP"
        and row["stratum_type"] == "shared_k"
        and row["stratum_value"] != "0–1"
    ]
    outcome_rows = build_fixed_outcome_rows(compatibility_source)
    fixed_stats, structural_summary = summarize_fixed_outcomes(outcome_rows)
    outcome_chart_rows: list[dict[str, Any]] = []
    for endpoint in ENDPOINT_ZH:
        endpoint_stats_fixed = fixed_stats[endpoint]
        grouped = (
            ("strong", "嚴格一致", endpoint_stats_fixed["strong"]),
            (
                "unresolved_compatible",
                "未決相容",
                endpoint_stats_fixed["resolution"] + endpoint_stats_fixed["weak"],
            ),
            ("conflict", "衝突", endpoint_stats_fixed["conflict"]),
            ("not_evaluable", "不可判定", endpoint_stats_fixed["not_evaluable"]),
        )
        for outcome, label, count in grouped:
            outcome_chart_rows.append(
                {
                    "endpoint": endpoint,
                    "layer_mapping": f"{ENDPOINT_ZH[endpoint]} · fixed swap-tolerant",
                    "outcome": outcome,
                    "outcome_label": label,
                    "n": count,
                    "share": count / FIXED_DENOMINATOR,
                }
            )
    read_primary = fixed_stats["read_full"]
    vaf_primary = fixed_stats["vaf_official"]
    site_inventory = build_caller_site_inventory(containment_summary)
    site_pair_summary = build_site_pair_summary(containment_summary)
    site_pair_chart_rows: list[dict[str, Any]] = []
    site_pair_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in site_pair_summary:
        site_pair_counts[row["endpoint"]][row["outcome"]] = as_int(row["n"])
    for endpoint in ENDPOINT_ZH:
        c = site_pair_counts[endpoint]
        grouped = (
            ("determined_same", "確定一致", c["determined_same"]),
            (
                "unresolved_compatible",
                "未決相容",
                c["ambiguous_equal"] + c["one_sided_contained"] + c["ambiguous_overlap"],
            ),
            ("conflict", "衝突", c["conflict"]),
            ("not_evaluable", "不可判定", c["not_evaluable"]),
        )
        for outcome, label, count in grouped:
            site_pair_chart_rows.append(
                {
                    "endpoint": endpoint,
                    "layer": ENDPOINT_ZH[endpoint],
                    "outcome": outcome,
                    "outcome_label": label,
                    "n": count,
                    "share": count / 8_096,
                }
            )
    examples = build_postprocess_examples(example_source)
    definitions = build_definitions()

    null_candidates = [
        data_dir / "hcc1395_topology_null.tsv",
        data_dir / "hcc1395_site_topology_null.tsv",
    ]
    null_path = next((path for path in null_candidates if path.exists()), None)
    null_rows = normalize_null(null_path) if null_path else []
    null_available = bool(null_rows)
    claim_ceiling = build_claim_ceiling(null_available)

    def null_row(endpoint: str, metric: str) -> dict[str, Any]:
        selected = [
            row
            for row in null_rows
            if row["endpoint"] == endpoint
            and row["comparison"] == "swap_tolerant"
            and row["metric"] == metric
        ]
        if len(selected) != 1:
            raise RuntimeError(
                f"expected one swap-tolerant null row for {endpoint}/{metric}, got {len(selected)}"
            )
        return selected[0]

    null_headline = {
        (endpoint, metric): null_row(endpoint, metric)
        for endpoint in ("read_full", "vaf_official")
        for metric in ("exact", "robust_nested")
    } if null_available else {}

    inventory_index = {row["relation"]: row for row in site_inventory}
    equal_sites = inventory_index["equal"]["n"]
    nested_sites = (
        inventory_index["a_proper_subset_b"]["n"]
        + inventory_index["b_proper_subset_a"]["n"]
    )
    directional_read = read_primary["induced"]
    directional_vaf = vaf_primary["induced"]

    site_pair_index = {
        (row["endpoint"], row["outcome"]): row for row in site_pair_summary
    }
    read_pair_determined = site_pair_index[("read_full", "determined_same")]
    vaf_pair_determined = site_pair_index[("vaf_official", "determined_same")]
    read_pair_conflict = site_pair_index[("read_full", "conflict")]
    vaf_pair_conflict = site_pair_index[("vaf_official", "conflict")]

    headline = [
        {
            "fixed_denominator": FIXED_DENOMINATOR,
            "read_evaluable_coverage": read_primary["evaluable_coverage"],
            "read_strong_n": read_primary["strong"],
            "read_strong_yield": read_primary["strong_yield"],
            "read_strong_conditional": read_primary["strong_conditional"],
            "read_pair_determined_n": read_pair_determined["n"],
            "read_pair_determined_conditional": read_pair_determined["conditional_share"],
            "read_pair_conflict_n": read_pair_conflict["n"],
            "vaf_evaluable_coverage": vaf_primary["evaluable_coverage"],
            "vaf_strong_n": vaf_primary["strong"],
            "vaf_strong_yield": vaf_primary["strong_yield"],
            "vaf_strong_conditional": vaf_primary["strong_conditional"],
            "vaf_pair_determined_n": vaf_pair_determined["n"],
            "vaf_pair_determined_conditional": vaf_pair_determined["conditional_share"],
            "vaf_pair_conflict_n": vaf_pair_conflict["n"],
            "site_equal_share": equal_sites / FIXED_DENOMINATOR,
            "site_nested_share": nested_sites / FIXED_DENOMINATOR,
            "shared_allele_identity": as_int(
                containment_summary["site_inventory"]["caller_shared_allele_identity"]
            ),
            "read_exact_null_excess": (
                null_headline[("read_full", "exact")]["delta"] if null_available else None
            ),
            "read_nested_null_excess": (
                null_headline[("read_full", "robust_nested")]["delta"]
                if null_available
                else None
            ),
            "vaf_exact_null_excess": (
                null_headline[("vaf_official", "exact")]["delta"]
                if null_available
                else None
            ),
            "vaf_nested_null_excess": (
                null_headline[("vaf_official", "robust_nested")]["delta"]
                if null_available
                else None
            ),
        }
    ]

    post_check(
        "read_fixed_strong_count",
        read_primary["strong"],
        1_599,
        "Read strict-full-exact plus real induced-substructure outcomes.",
    )
    post_check(
        "vaf_fixed_strong_count",
        vaf_primary["strong"],
        1_790,
        "Official VAF strict-full-exact plus real induced-substructure outcomes.",
    )
    post_check(
        "site_pair_evidence_total_per_layer",
        sorted({row["total_site_pairs"] for row in site_pair_summary}),
        [8_096],
        "All evidence layers use the same shared-site-pair ledger.",
    )

    checks_dataset = [
        {
            "check": row["check"],
            "pass": "PASS" if as_bool(row["pass"]) else "FAIL",
            "observed": row["observed"],
            "expected": row["expected"],
            "detail": "Producer QA",
        }
        for row in checks
    ] + [
        {
            "check": row["check"],
            "pass": "PASS" if as_bool(row["pass"]) else "FAIL",
            "observed": row["observed"],
            "expected": row["expected"],
            "detail": "Postprocess QA",
        }
        for row in containment_checks
    ] + [
        {
            "check": row["check"],
            "pass": "PASS" if as_bool(row["status"]) else "FAIL",
            "observed": row["observed"],
            "expected": row["expected"],
            "detail": "Site-label null QA",
        }
        for row in null_checks_source
    ] + report_checks

    datasets: dict[str, list[dict[str, Any]]] = {
        "headline": headline,
        "outcomes": outcome_rows,
        "outcome_chart": outcome_chart_rows,
        "candidate_bounds": bounds_rows,
        "site_inventory": site_inventory,
        "site_pair_summary": site_pair_summary,
        "site_pair_chart": site_pair_chart_rows,
        "complexity": complexity_rows,
        "complexity_k": complexity_k,
        "structural_summary": structural_summary,
        "examples": examples,
        "definitions": definitions,
        "claim_ceiling": claim_ceiling,
        "annotation_context": annotation_context,
        "checks": checks_dataset,
    }
    if null_available:
        datasets["conditional_null"] = null_rows

    filters = [
        "autosomes chr1-22",
        "HCC1395 versus HCC1395_DORADO",
        "exact-coordinate complete-both region pairs",
        "fixed population n=5,720",
    ]
    generated_metric_definitions = [
        "evaluable coverage = evaluable regions / 5,720",
        "robust compatibility = exact candidate set + either proper-subset candidate-set relation",
        "all-pair robust yield = robust compatible regions / 5,720",
        "conditional robust rate = robust compatible regions / evaluable regions",
        "any-overlap optimistic rate additionally counts partial candidate-set overlap",
    ]
    sources = [
        source_spec(
            "site_regions",
            "HCC1395 site-topology region outcomes",
            paths["regions"],
            "Full 5,720-row site-aware read-full and VAF-top candidate-set comparison.",
            generated_at,
            filters=filters,
            metric_definitions=generated_metric_definitions,
        ),
        source_spec(
            "site_metrics",
            "HCC1395 site-topology producer metrics",
            paths["metrics"],
            "Producer-level reconciled endpoint, mapping, population, and outcome counts.",
            generated_at,
            filters=filters,
            metric_definitions=generated_metric_definitions,
        ),
        source_spec(
            "site_checks",
            "HCC1395 site-topology validation checks",
            paths["checks"],
            "Producer conservation and candidate re-enumeration checks.",
            generated_at,
        ),
        source_spec(
            "fixed_outcomes",
            "Corrected fixed region outcomes",
            paths["pair_outcomes"],
            (
                "One mutually exclusive whole-region outcome per exact-coordinate pair, "
                "using each layer's precomputed identity-or-one-global-HP-swap mapping."
            ),
            generated_at,
            filters=filters,
        ),
        source_spec(
            "site_pair_evidence",
            "Shared-site-pair relation evidence",
            paths["site_relation_evidence"],
            (
                "One row per shared genomic-position pair; read-full selected mapping is frozen "
                "and reused for VAF display to prevent per-site mapping cherry-picking."
            ),
            generated_at,
            filters=filters,
        ),
        source_spec(
            "allele_identity",
            "Shared caller sSNV allele identity ledger",
            paths["allele_identity"],
            "CHROM/POS/REF/ALT identity audit for every caller-shared sSNV locus.",
            generated_at,
            filters=filters,
        ),
        source_spec(
            "method_contract",
            "Site-topology pre-decision audit",
            TOPIC / "pre-decision-audit.md",
            "Frozen operational definitions, falsifiers, matched-null gate, and claim ceiling.",
            generated_at,
        ),
        source_spec(
            "hcc1395_kb",
            "HCC1395 sample provenance KB",
            Path("02_samples/hcc1395.md"),
            (
                "Canonical local Knowledge Base entry: 5 kHz material was sequenced by the Clair team; "
                "the DORADO public dataset was released with Google DeepSomatic and its raw HCC1395 sequencing originated at NYGC."
            ),
            generated_at,
        ),
        source_spec(
            "annotation_context",
            "HCC1395 cancer-gene/drug annotation matched-null context",
            ANNOTATION_CONTEXT_PATH,
            (
                "Previously validated coarse-topology agreement stratified by COSMIC CGC and "
                "DGIdb approved-antineoplastic gene-body annotations with a chromosome and region-length matched null."
            ),
            generated_at,
            filters=filters,
            metric_definitions=[
                "present-minus-absent difference is measured in percentage points",
                "permutation p is two-sided under the chromosome plus region-length-decile conditional null",
            ],
        ),
    ]
    if null_available and null_path is not None:
        sources.append(
            source_spec(
                "conditional_null",
                "Within-region shared-site-label permutation null",
                null_path,
                (
                    "Permutes dataset-B shared-site labels within each observed mapped HP component; "
                    "preserves region, k, candidate-set size, candidate multiplicity, and HP mapping. "
                    "Chromosome blocks are used for uncertainty, not for the permutation itself."
                ),
                generated_at,
                filters=filters,
            )
        )

    dataset_source_inputs: dict[str, list[Path]] = {
        "headline": [paths["pair_outcomes"], paths["site_relation_evidence"]],
        "outcomes": [paths["compatibility"]],
        "outcome_chart": [paths["compatibility"]],
        "candidate_bounds": [paths["regions"], paths["metrics"]],
        "site_inventory": [paths["containment_summary"], paths["allele_identity"]],
        "site_pair_summary": [paths["containment_summary"], paths["site_relation_evidence"]],
        "site_pair_chart": [paths["containment_summary"], paths["site_relation_evidence"]],
        "complexity": [paths["regions"]],
        "complexity_k": [paths["regions"]],
        "structural_summary": [paths["compatibility"]],
        "examples": [paths["examples"]],
        "definitions": [TOPIC / "pre-decision-audit.md"],
        "claim_ceiling": [TOPIC / "pre-decision-audit.md"],
        "annotation_context": [ANNOTATION_CONTEXT_PATH],
        "checks": [
            paths["checks"],
            paths["containment_checks"],
            paths["null_checks"],
            paths["regions"],
            paths["metrics"],
        ],
    }
    if null_available and null_path is not None:
        dataset_source_inputs["conditional_null"] = [null_path]

    frozen_datasets: dict[str, list[dict[str, Any]]] = {}
    for dataset_id, rows in datasets.items():
        query, frozen = freeze_dataset_with_sql(dataset_id, rows)
        frozen_datasets[dataset_id] = frozen
        sources.append(
            {
                "id": f"dataset_{dataset_id}",
                "label": f"Frozen report dataset: {dataset_id}",
                "path": (
                    TOPIC / "scripts/build_site_topology_report_artifact.py"
                ).as_posix(),
                "query": {
                    "engine": "sqlite",
                    "language": "sql",
                    "sql": query,
                    "description": (
                        f"Literal reviewed snapshot rows for the {dataset_id} report dataset; "
                        "materialized and row-count checked in SQLite by the report builder."
                    ),
                    "executed_at": generated_at,
                    "tables_used": [
                        path.as_posix() for path in dataset_source_inputs[dataset_id]
                    ],
                    "filters": filters,
                    "metric_definitions": generated_metric_definitions,
                },
            }
        )
    datasets = frozen_datasets

    cards = [
        {
            "id": "population_card",
            "description": "Exact-coordinate complete-both HCC1395 technical region pairs.",
            "dataset": "headline",
            "sourceId": "fixed_outcomes",
            "metrics": [
                {"label": "Fixed denominator", "field": "fixed_denominator", "format": "number"},
                {"label": "Equal site sets", "field": "site_equal_share", "format": "percent"},
                {"label": "Shared alleles identical", "field": "shared_allele_identity", "format": "number"},
            ],
        },
        {
            "id": "read_card",
            "description": "Full read-pattern candidate space; fixed best of ordered or one global HP swap.",
            "dataset": "headline",
            "sourceId": "fixed_outcomes",
            "metrics": [
                {"label": "Strict exact / real induced", "field": "read_strong_yield", "format": "percent"},
                {"label": "Evaluable coverage", "field": "read_evaluable_coverage", "format": "percent"},
                {"label": "Conditional strict", "field": "read_strong_conditional", "format": "percent"},
            ],
        },
        {
            "id": "vaf_card",
            "description": "Exact-score VAF top set with ties retained; same fixed whole-region mapping rule.",
            "dataset": "headline",
            "sourceId": "fixed_outcomes",
            "metrics": [
                {"label": "Strict exact / real induced", "field": "vaf_strong_yield", "format": "percent"},
                {"label": "Evaluable coverage", "field": "vaf_evaluable_coverage", "format": "percent"},
                {"label": "Conditional strict", "field": "vaf_strong_conditional", "format": "percent"},
            ],
        },
    ]

    charts = [
        {
            "id": "outcome_composition_chart",
            "title": "Read-full 與 VAF-top 的固定 whole-region outcome",
            "subtitle": "每個 evidence layer 使用事先計算的 identity 或一次 global HP swap；每列 n=5,720。",
            "type": "horizontalStackedBar100",
            "dataset": "outcome_chart",
            "sourceId": "fixed_outcomes",
            "intent": "composition",
            "question": "每個證據層在完整 exact、真誘導子結構、解析度嵌套、shared core、conflict 與不可判定之間如何分配？",
            "rationale": "A 100% composition view preserves the fixed denominator and exposes coverage loss instead of conditioning it away.",
            "comparisonContext": {
                "denominator": "5,720 exact-coordinate complete-both regions per layer",
                "grain": "region pair",
                "normalization": "within layer to 100%",
                "unit": "share",
            },
            "encodings": {
                "x": {"field": "layer_mapping", "type": "nominal", "label": "Evidence layer / HP mapping"},
                "y": {"field": "n", "type": "quantitative", "label": "Regions"},
                "color": {"field": "outcome_label", "type": "nominal", "label": "Candidate-set relation"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True, "title": "關係類別"},
        },
        {
            "id": "site_inventory_chart",
            "title": "區域主要 sSNV 位點集對應",
            "subtitle": "HCC1395 對 HCC1395_DORADO；全 5,720 區域。",
            "type": "horizontalBar",
            "dataset": "site_inventory",
            "sourceId": "allele_identity",
            "intent": "comparison",
            "question": "有多少區域具備相同、方向性嵌套或只部分重疊的位點集？",
            "rationale": "Horizontal bars make the dominant equal-site population and the much smaller containment-eligible population explicit.",
            "comparisonContext": {
                "denominator": "5,720 exact-coordinate complete-both regions",
                "grain": "region pair",
                "unit": "regions",
            },
            "encodings": {
                "x": {"field": "relation_label", "type": "nominal", "label": "Site-set relation"},
                "y": {"field": "n", "type": "quantitative", "label": "Regions"},
            },
            "valueFormat": "number",
            "palette": {"kind": "sequential"},
            "labels": {"values": "all"},
        },
        {
            "id": "site_pair_outcome_chart",
            "title": "每一對共享 sSNV 的關係支持結果",
            "subtitle": "Read-full selected mapping 先固定，VAF 在同一 mapping 下顯示；ambiguous equal 不等於 determined same。",
            "type": "horizontalStackedBar100",
            "dataset": "site_pair_chart",
            "sourceId": "site_pair_evidence",
            "intent": "composition",
            "question": "8,096 個共享位點對中，有多少是雙側唯一且同一關係、不確定集合嵌套、衝突或不可判定？",
            "rationale": "The composition separates resolved agreement from candidate ambiguity instead of merging both into one optimistic compatibility rate.",
            "comparisonContext": {
                "denominator": "8,096 shared-site pairs per evidence layer",
                "grain": "shared genomic-position pair under a fixed region mapping",
                "normalization": "within evidence layer to 100%",
                "unit": "share",
            },
            "encodings": {
                "x": {"field": "layer", "type": "nominal", "label": "Evidence layer"},
                "y": {"field": "n", "type": "quantitative", "label": "Shared-site pairs"},
                "color": {"field": "outcome_label", "type": "nominal", "label": "Relation support"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True, "title": "逐位點對 outcome"},
        },
        {
            "id": "candidate_bounds_chart",
            "title": "Robust 與 any-overlap 的候選不確定邊界",
            "subtitle": "All-pair yield；固定分母 n=5,720，不將不可判定區域排除。",
            "type": "bar",
            "dataset": "candidate_bounds",
            "sourceId": "site_regions",
            "intent": "comparison",
            "question": "只要有一部分候選交集時，樂觀一致率會比 robust exact/nested 高多少？",
            "rationale": "Paired bars show the multiplicity-sensitive optimistic bound beside the stricter complete candidate-space relation.",
            "comparisonContext": {
                "denominator": "5,720 exact-coordinate complete-both regions",
                "grain": "region pair",
                "unit": "all-pair yield",
            },
            "encodings": {
                "x": {"field": "layer_mapping", "type": "nominal", "label": "Evidence layer / HP mapping"},
                "y": {"field": "all_pair_yield", "type": "quantitative", "label": "All-pair yield", "format": "percent"},
                "color": {"field": "bound_label", "type": "nominal", "label": "Compatibility definition"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True, "title": "邊界"},
        },
        {
            "id": "complexity_k_chart",
            "title": "Shared-site k 分層的 robust all-pair yield",
            "subtitle": "Ordered HP primary；每個 k stratum 的區域數為分母。",
            "type": "bar",
            "dataset": "complexity_k",
            "sourceId": "site_regions",
            "intent": "comparison",
            "question": "結構相容是否只由低 k 簡單區域驅動？",
            "rationale": "Grouped bars compare read-full and VAF-top yields across increasing shared-site complexity without mixing denominators.",
            "comparisonContext": {
                "denominator": "all regions within each shared-site k stratum",
                "grain": "region pair",
                "unit": "robust all-pair yield",
            },
            "encodings": {
                "x": {"field": "stratum_value", "type": "ordinal", "label": "Shared-site k"},
                "y": {"field": "robust_yield", "type": "quantitative", "label": "Robust yield", "format": "percent"},
                "color": {"field": "layer", "type": "nominal", "label": "Evidence layer"},
            },
            "valueFormat": "percent",
            "palette": {"kind": "categorical"},
            "legend": {"position": "bottom", "interactive": True, "title": "證據層"},
        },
    ]

    tables = [
        {
            "id": "outcomes_table",
            "title": "固定分母、互斥 whole-region outcome 完整數據",
            "subtitle": "Read full / official VAF / normalized VAF sensitivity；每層 9 個 outcome 回加為 5,720。",
            "dataset": "outcomes",
            "sourceId": "fixed_outcomes",
            "density": "dense",
            "defaultSort": {"field": "n", "direction": "desc"},
            "columns": columns(
                ("layer", "Evidence layer", {"type": "text"}),
                ("mapping_label", "HP mapping", {"type": "text"}),
                ("outcome_label", "Candidate-set relation", {"type": "text"}),
                ("n", "Regions", {"format": "number"}),
                ("fixed_denominator", "Fixed denominator", {"format": "number"}),
                ("share", "All-pair share", {"format": "percent"}),
            ),
        },
        {
            "id": "structural_summary_table",
            "title": "一致拓撲與方向性子結構的可審核計數",
            "subtitle": "strong exact-or-induced 是組合列，其餘主要列為互斥計數；不把 resolution nesting 誤稱為真子樹。",
            "dataset": "structural_summary",
            "sourceId": "fixed_outcomes",
            "density": "dense",
            "defaultSort": {"field": "n", "direction": "desc"},
            "columns": columns(
                ("layer", "Evidence layer", {"type": "text"}),
                ("mapping_label", "HP mapping", {"type": "text"}),
                ("metric_label", "Structural metric", {"type": "text"}),
                ("n", "Regions", {"format": "number"}),
                ("share", "Share of 5,720", {"format": "percent"}),
                ("conditional_share", "Share of evaluable", {"format": "percent"}),
                ("exclusive", "Exclusive metric", {"type": "text"}),
                ("definition", "Definition / boundary", {"type": "text"}),
            ),
        },
        {
            "id": "site_inventory_table",
            "title": "位點集嵌套適用母群",
            "subtitle": "只有真子集關係可直接稱方向性 induced substructure。",
            "dataset": "site_inventory",
            "sourceId": "allele_identity",
            "density": "spacious",
            "defaultSort": {"field": "n", "direction": "desc"},
            "columns": columns(
                ("relation_label", "Site-set relation", {"type": "text"}),
                ("n", "Regions", {"format": "number"}),
                ("share", "Share", {"format": "percent"}),
                ("directional_containment_eligible", "Directional containment eligible", {"type": "text"}),
            ),
        },
        {
            "id": "site_pair_summary_table",
            "title": "逐共享位點對關係支持完整計數",
            "subtitle": "determined same 必須是雙側 singleton 且同一 F/R/P 關係；ambiguous equal 只表示候選關係集合相同。",
            "dataset": "site_pair_summary",
            "sourceId": "site_pair_evidence",
            "density": "dense",
            "defaultSort": {"field": "n", "direction": "desc"},
            "columns": columns(
                ("layer", "Evidence layer", {"type": "text"}),
                ("outcome_label", "Site-pair outcome", {"type": "text"}),
                ("n", "Site pairs", {"format": "number"}),
                ("total_site_pairs", "All site pairs", {"format": "number"}),
                ("evaluable_site_pairs", "Evaluable", {"format": "number"}),
                ("all_pair_share", "Share of 8,096", {"format": "percent"}),
                ("conditional_share", "Share of evaluable", {"format": "percent"}),
            ),
        },
        {
            "id": "bounds_table",
            "title": "Candidate multiplicity 下界與樂觀上界",
            "subtitle": "Robust 不計 overlap；optimistic 額外計入任一 projected-signature 交集。",
            "dataset": "candidate_bounds",
            "sourceId": "site_regions",
            "density": "dense",
            "defaultSort": {"field": "all_pair_yield", "direction": "desc"},
            "columns": columns(
                ("layer", "Evidence layer", {"type": "text"}),
                ("mapping_label", "HP mapping", {"type": "text"}),
                ("bound_label", "Definition", {"type": "text"}),
                ("n", "Compatible", {"format": "number"}),
                ("evaluable", "Evaluable", {"format": "number"}),
                ("evaluable_coverage", "Evaluable coverage", {"format": "percent"}),
                ("conditional_rate", "Conditional rate", {"format": "percent"}),
                ("all_pair_yield", "All-pair yield", {"format": "percent"}),
            ),
        },
        {
            "id": "complexity_table",
            "title": "Shared k、HP components 與 candidate multiplicity 分層",
            "subtitle": "Wilson 95% interval 是區域級描述區間；尚未取代 chromosome-block uncertainty。",
            "dataset": "complexity",
            "sourceId": "site_regions",
            "density": "dense",
            "defaultSort": {"field": "population", "direction": "desc"},
            "columns": columns(
                ("layer", "Evidence layer", {"type": "text"}),
                ("mapping_label", "HP mapping", {"type": "text"}),
                ("stratum_label", "Stratum", {"type": "text"}),
                ("stratum_value", "Level", {"type": "text"}),
                ("population", "Population", {"format": "number"}),
                ("evaluable", "Evaluable", {"format": "number"}),
                ("robust", "Robust", {"format": "number"}),
                ("robust_yield", "All-pair yield", {"format": "percent"}),
                ("robust_conditional", "Conditional robust", {"format": "percent"}),
                ("ci95_low", "Wilson low", {"format": "percent"}),
                ("ci95_high", "Wilson high", {"format": "percent"}),
                ("possible_only", "Possible only", {"format": "number"}),
                ("conflict", "Conflict", {"format": "number"}),
            ),
        },
        {
            "id": "examples_table",
            "title": "代表區域（每層每類固定取前兩例）",
            "subtitle": "來自 corrected fixed-outcome ledger，可回查 mapping、candidate relation 與候選數。",
            "dataset": "examples",
            "sourceId": "fixed_outcomes",
            "density": "dense",
            "defaultSort": {"field": "chrom", "direction": "asc"},
            "columns": columns(
                ("layer", "Evidence layer", {"type": "text"}),
                ("outcome_label", "Fixed outcome", {"type": "text"}),
                ("example_rank", "Rank", {"format": "number"}),
                ("chrom", "Chromosome", {"type": "text"}),
                ("region", "Region", {"type": "text"}),
                ("site_set_relation_label", "Site relation", {"type": "text"}),
                ("shared_k", "Shared k", {"format": "number"}),
                ("hp_pair", "HP count", {"type": "text"}),
                ("selected_mapping", "Selected mapping", {"type": "text"}),
                ("candidate_relation", "Candidate relation", {"type": "text"}),
                ("topology_set_size_a", "HCC topology-set Q", {"format": "number"}),
                ("topology_set_size_b", "DORADO topology-set Q", {"format": "number"}),
                ("outcome_reason", "Reason", {"type": "text"}),
            ),
        },
        {
            "id": "definitions_table",
            "title": "比對單位、指標與邊界",
            "subtitle": "先定義 shared-site projection，再解讀數字。",
            "dataset": "definitions",
            "sourceId": "method_contract",
            "density": "spacious",
            "defaultSort": {"field": "term", "direction": "asc"},
            "columns": columns(
                ("term", "Term", {"type": "text"}),
                ("definition", "Operational definition", {"type": "text"}),
                ("boundary", "Interpretation boundary", {"type": "text"}),
            ),
        },
        {
            "id": "claim_table",
            "title": "證據可支持與不可支持的說法",
            "subtitle": "本報告上限是 same-cell-line cross-source partial technical reproducibility。",
            "dataset": "claim_ceiling",
            "sourceId": "method_contract",
            "density": "spacious",
            "defaultSort": {"field": "evidence", "direction": "asc"},
            "columns": columns(
                ("evidence", "Evidence", {"type": "text"}),
                ("allowed", "Allowed claim", {"type": "text"}),
                ("not_allowed", "Not supported", {"type": "text"}),
            ),
        },
        {
            "id": "checks_table",
            "title": "Producer 與報告端 QA checks",
            "subtitle": "包含 5,720 恒等、outcome 回加、metrics reconciliation 與 lower≤upper。",
            "dataset": "checks",
            "sourceId": "site_checks",
            "density": "dense",
            "defaultSort": {"field": "pass", "direction": "asc"},
            "columns": columns(
                ("check", "Check", {"type": "text"}),
                ("pass", "Status", {"type": "text"}),
                ("observed", "Observed", {"type": "text"}),
                ("expected", "Expected", {"type": "text"}),
                ("detail", "Scope", {"type": "text"}),
            ),
        },
        {
            "id": "annotation_context_table",
            "title": "Cancer-gene / drug annotation context",
            "subtitle": "Previously validated matched-null results; context only, not a topology truth set.",
            "dataset": "annotation_context",
            "sourceId": "annotation_context",
            "density": "spacious",
            "defaultSort": {"field": "p_two_sided", "direction": "asc"},
            "columns": columns(
                ("feature", "Annotation stratum", {"type": "text"}),
                ("n_present", "Present n", {"format": "number"}),
                ("present_agreement", "Coarse agreement", {"format": "percent"}),
                ("present_minus_absent_pp", "Present - absent (pp)", {"format": "number"}),
                ("null_q025_pp", "Null q2.5 (pp)", {"format": "number"}),
                ("null_q975_pp", "Null q97.5 (pp)", {"format": "number"}),
                ("p_two_sided", "Matched-null p", {"format": "number"}),
                ("interpretation", "Boundary", {"type": "text"}),
            ),
        },
    ]

    if null_available:
        tables.append(
            {
                "id": "null_table",
                "title": "Within-region shared-site-label permutation null",
                "subtitle": "B 端共享位點標籤在 mapped HP component 內置換；chromosome block 只用於 uncertainty。",
                "dataset": "conditional_null",
                "sourceId": "conditional_null",
                "density": "dense",
                "defaultSort": {"field": "delta", "direction": "desc"},
                "columns": columns(
                    ("layer", "Evidence layer", {"type": "text"}),
                    ("mapping_label", "HP mapping", {"type": "text"}),
                    ("metric", "Metric", {"type": "text"}),
                    ("observed", "Observed", {"format": "percent"}),
                    ("null_mean", "Null mean", {"format": "percent"}),
                    ("null_q025", "Null q2.5", {"format": "percent"}),
                    ("null_q975", "Null q97.5", {"format": "percent"}),
                    ("delta", "Excess", {"format": "percent", "movement": True}),
                    ("p_empirical", "Empirical p", {"format": "number"}),
                    ("permutations", "Permutations", {"format": "number"}),
                    ("block_ci_low", "Block CI low", {"format": "percent"}),
                    ("block_ci_high", "Block CI high", {"format": "percent"}),
                    ("loco_positive_chromosomes", "LOCO positive chr", {"format": "number"}),
                    ("loco_chromosomes", "LOCO chr", {"format": "number"}),
                    ("null_method", "Permutation method", {"type": "text"}),
                    ("vaf_condition", "VAF condition", {"type": "text"}),
                ),
            }
        )

    technical_summary = (
        "## 技術總結：支持部分跨來源結構再現；不等於已證明方法準確\n\n"
        f"- **固定母群與位點 QA：** {FIXED_DENOMINATOR:,} 個 chr1–22 exact-coordinate complete-both regions；"
        f"caller 位點集 {equal_sites:,}（{fmt_pct(equal_sites / FIXED_DENOMINATOR)}）完全相同、"
        f"{nested_sites:,}（{fmt_pct(nested_sites / FIXED_DENOMINATOR)}）為單方真子集；"
        f"15,713/15,713 個共享 sSNV 的 CHROM/POS/REF/ALT 相同。\n"
        f"- **Read 直接支持的完整候選空間：** 可判定 {read_primary['evaluable']:,} 區"
        f"（{fmt_pct(read_primary['evaluable_coverage'])}）；strict-full-exact 或真正 induced substructure 為 "
        f"**{read_primary['strong']:,}/{FIXED_DENOMINATOR:,} = {fmt_pct(read_primary['strong_yield'])}**，"
        f"在可判定區為 {fmt_pct(read_primary['strong_conditional'])}；其中真正 induced 為 {directional_read:,} 區。\n"
        f"- **VAF 推測的 exact-score top set：** 可判定 {vaf_primary['evaluable']:,} 區"
        f"（{fmt_pct(vaf_primary['evaluable_coverage'])}）；strict-full-exact 或真正 induced 為 "
        f"**{vaf_primary['strong']:,}/{FIXED_DENOMINATOR:,} = {fmt_pct(vaf_primary['strong_yield'])}**，"
        f"在可判定區為 {fmt_pct(vaf_primary['strong_conditional'])}；其中真正 induced 為 {directional_vaf:,} 區。"
        "VAF 與候選樹來自相關 read evidence，並非獨立驗證。\n"
        f"- **逐共享位點對：** Read 層 determined-same 為 {read_pair_determined['n']:,}/"
        f"{read_pair_determined['evaluable_site_pairs']:,} = {fmt_pct(read_pair_determined['conditional_share'])}；"
        f"VAF 層為 {vaf_pair_determined['n']:,}/{vaf_pair_determined['evaluable_site_pairs']:,} = "
        f"{fmt_pct(vaf_pair_determined['conditional_share'])}。ambiguous-equal 不併入 determined-same。\n"
        + (
            f"- **Site-label null：** projected candidate-set exact 的 observed-minus-null 為 "
            f"Read **{100 * null_headline[('read_full', 'exact')]['delta']:.2f} pp**、"
            f"VAF **{100 * null_headline[('vaf_official', 'exact')]['delta']:.2f} pp**；"
            "兩者 empirical p=1/5001、block CI>0、LOCO 22/22。這個 p 值不直接檢驗上面的 fixed strict-exact+true-induced 合計。"
            if null_available
            else "- **Null gate 未完成：** 尚無 within-region site-label permutation baseline。"
        )
    )

    freeze_body = (
        "**PARTIAL · ENGINEERING VALIDATION ONLY — historical layered-v2 snapshot.** "
        "同一 HCC1395 cell line 的兩個 technical sequencing/processing sources，不是獨立 biological replicates；"
        "clean-v3 release closeout 仍無效，且沒有 external biological truth。"
    )
    verdict_body = (
        "### 目前判定\n\n"
        "**PARTIAL technical reproducibility; SCIENTIFIC NO-GO for proof of accuracy/effectiveness.** "
        "逐位點比對比單看粗分類比例更貼近「結構是否重現」；"
        "但 fixed structural endpoint 尚無直接 null，且在 depth-matched split-half 與 single-cell/synthetic truth 之前，"
        "證據只能支持這一對 technical datasets 中的部分結構現象再現。"
    )

    blocks: list[dict[str, Any]] = [
        {
            "id": "title",
            "type": "markdown",
            "body": "# HCC1395 逐位點拓撲子結構跨來源驗證",
        },
        {"id": "partial_ribbon", "type": "markdown", "body": freeze_body},
        {
            "id": "technical_summary",
            "type": "markdown",
            "body": technical_summary,
            "sourceId": "fixed_outcomes",
        },
        {"id": "verdict", "type": "markdown", "body": verdict_body},
        {
            "id": "dataset_provenance",
            "type": "markdown",
            "body": (
                "### 兩個 HCC1395 是 same-cell-line cross-source technical datasets\n\n"
                "Local HCC1395 KB 記載：5 kHz 資料由 Clair 團隊定序；DORADO 資料為 Google DeepSomatic "
                "公開來源，HCC1395 原始定序取自 NYGC。因此最安全的口徑是「同一 cell line 的跨來源技術資料」；"
                "aliquot、passaging 與 molecule-level independence 尚不明，不稱兩個獨立生物複本。"
            ),
            "sourceId": "hcc1395_kb",
        },
        {
            "id": "headline_metrics",
            "type": "metric-strip",
            "cardIds": ["population_card", "read_card", "vaf_card"],
        },
        {
            "id": "outcome_heading",
            "type": "markdown",
            "body": (
                "## 固定 whole-region outcome：真子結構與候選解析度差異分開\n\n"
                "Read-full 與 VAF-top 都以 shared-site F/R/P relation signature 比對。"
                "主結論只把 strict-full-exact 與真正 induced substructure 合計；"
                "相同位點但一側候選空間較窄者另列 resolution nesting，不冒充真子樹。"
            ),
            "sourceId": "fixed_outcomes",
        },
        {"id": "outcome_chart", "type": "chart", "chartId": "outcome_composition_chart"},
        {"id": "outcomes_table", "type": "table", "tableId": "outcomes_table"},
        {
            "id": "structural_count_heading",
            "type": "markdown",
            "body": (
                "### 「一致拓撲」與「可接受子結構」必須分開計數\n\n"
                "位點集完全相同時，檢查完整 projected candidate-set exact；"
                "只有一側位點集是另一側真子集時，才使用方向性 induced-substructure 說法。"
                "strong exact-or-induced 是組合列，不得再與其 constituent rows 相加。"
            ),
            "sourceId": "fixed_outcomes",
        },
        {"id": "structural_summary", "type": "table", "tableId": "structural_summary_table"},
        {
            "id": "site_inventory_heading",
            "type": "markdown",
            "body": (
                f"## 位點集嵌套規則只直接適用於 {nested_sites:,}/{FIXED_DENOMINATOR:,} 區\n\n"
                f"Caller VCF inventory 中，{equal_sites:,} 區的位點集完全相同；"
                f"{nested_sites:,} 區是方向性真子集。"
                "部分重疊位點集即使共享投影一致，也只能稱 shared core，不稱子樹。"
            ),
            "sourceId": "allele_identity",
        },
        {"id": "site_inventory_chart", "type": "chart", "chartId": "site_inventory_chart"},
        {"id": "site_inventory_table", "type": "table", "tableId": "site_inventory_table"},
        {
            "id": "site_pair_heading",
            "type": "markdown",
            "body": (
                "## 每一對共享 sSNV 都以其支持的 F/R/P 關係集合比對\n\n"
                "F 表示較小 genomic position 為祖先、R 表示反向祖先、P 表示平行；"
                "只有雙側都只剩同一個關係才是 determined same。"
                "關係支持集相等或單側被包含仍保留為 ambiguity/containment，不升格成已確定拓撲。"
            ),
            "sourceId": "site_pair_evidence",
        },
        {"id": "site_pair_chart", "type": "chart", "chartId": "site_pair_outcome_chart"},
        {"id": "site_pair_table", "type": "table", "tableId": "site_pair_summary_table"},
        {
            "id": "bounds_heading",
            "type": "markdown",
            "body": (
                "## 補充診斷：candidate multiplicity 的 robust 下界與 any-overlap 上界\n\n"
                "候選數越多，「存在一對可對上」越容易偶合。"
                "此 ordered/global-swap candidate-space 診斷不是上方 fixed structural headline；"
                "any-overlap 額外計入 partial overlap，只作樂觀 sensitivity。"
            ),
            "sourceId": "site_regions",
        },
        {"id": "bounds_chart", "type": "chart", "chartId": "candidate_bounds_chart"},
        {"id": "bounds_table", "type": "table", "tableId": "bounds_table"},
        {
            "id": "complexity_heading",
            "type": "markdown",
            "body": (
                "## 補充診斷：複雜度分層是否只由簡單區驅動\n\n"
                "報告分開 shared-site k、HP component pair 與 candidate-tree pair count。"
                "All-pair yield 保留每個 stratum 的不可判定成本；conditional robust 則只用可判定區。"
                "Wilson 區間未校正區域相關性，不可代替 chromosome-block bootstrap。"
            ),
            "sourceId": "site_regions",
        },
        {"id": "complexity_chart", "type": "chart", "chartId": "complexity_k_chart"},
        {"id": "complexity_table", "type": "table", "tableId": "complexity_table"},
    ]

    if null_available:
        blocks.extend(
            [
                {
                    "id": "null_heading",
                    "type": "markdown",
                    "body": (
                        "## Within-region site-label null：一致訊號超越標籤任意配對\n\n"
                        "在同一 region、同一 mapped HP component 內置換 B 端共享位點標籤；"
                        "保留 k、Q-set size、candidate multiplicity 與 HP mapping。"
                        "Chromosome block bootstrap 與 LOCO 只評估不確定度；VAF top set 凍結且沒有重新 scoring。"
                        "此 null 檢驗 projected candidate-set relation，不是 fixed strict-exact+true-induced 合計的直接 p 值。"
                    ),
                    "sourceId": "conditional_null",
                },
                {"id": "null_table", "type": "table", "tableId": "null_table"},
            ]
        )
    else:
        blocks.append(
            {
                "id": "null_pending",
                "type": "markdown",
                "body": (
                    "## Site-label null 尚未產出：這是不能宣稱方法已驗證的主要關卡\n\n"
                    "需保留 region、shared-site k、HP components 與 candidate-size 的 conditional null。"
                    "在此之前，報告只描述 observed compatibility，不提供 chance-corrected 證明。"
                ),
            }
        )

    blocks.extend(
        [
            {
                "id": "examples_heading",
                "type": "markdown",
                "body": (
                    "## 代表區域涵蓋 fixed outcome 各類別\n\n"
                    "每個 evidence layer／非空 outcome 依固定順序取前兩筆。"
                    "它們用來解釋 classifier，不是對總體的抽樣估計。"
                ),
                "sourceId": "fixed_outcomes",
            },
            {"id": "examples_table", "type": "table", "tableId": "examples_table"},
            {
                "id": "scope_heading",
                "type": "markdown",
                "body": (
                    "## 範圍與指標定義：比的是相同 genomic sSNV 上的關係支持\n\n"
                    "每個 shared sSNV pair 都編碼為正向祖先、反向祖先或 parallel；private events 投影時移除。"
                    "隱藏／private vertices 可在 shared-site projection 中收縮，但所有共享 mutation label、可達方向與 HP mapping 必須保留。"
                ),
                "sourceId": "method_contract",
            },
            {"id": "definitions_table", "type": "table", "tableId": "definitions_table"},
            {
                "id": "method_heading",
                "type": "markdown",
                "body": (
                    "## 方法：Read full 與 VAF top 使用相同投影合約，但證據強度不同\n\n"
                    "1. 重枚舉每個 region×HP 的完整 minimal-tree candidates，並與 frozen n_T/digest 對齊。\n"
                    "2. 在 HP-specific genomic universe 交集上，將每顆樹轉成所有共享位點對的 F/R/P relation vector。\n"
                    "3. 對投影後唯一 signature 集計算 exact、A⊂B、B⊂A、overlap 或 disjoint；同時保留原始 candidate-tree 數。\n"
                    "4. Fixed region outcome 在每個 evidence layer 預先選 identity 或整個 region 一次 HP1↔HP2 swap；不允許逐 component cherry-pick。\n"
                    "5. Site-pair ledger 先固定 Read-selected region mapping，再以同一 mapping 顯示 Read/VAF，避免逐位點選最好結果。\n"
                    "6. VAF 層使用 exact Fraction-score argmax，保留所有 ties；normalized score 只作 comparison-count sensitivity，不取代 official endpoint。\n"
                    "7. Null 在同 region、mapped HP component 內置換 B-site labels；VAF top set 固定、不重新 scoring。"
                ),
                "sourceId": "method_contract",
            },
            {
                "id": "limitations_heading",
                "type": "markdown",
                "body": (
                    "## 限制與穩健性：技術再現不等於生物準確\n\n"
                    "- 兩個 dataset 來自同一 HCC1395 cell line，可測 technical reproducibility，但不是 independent biological replication。\n"
                    "- Read candidate 與 VAF top 都來自相關的 read evidence；VAF 不是 orthogonal confirmation。\n"
                    "- Historical layered-v2 是 engineering snapshot；clean-v3 release closeout 未作為本報告主輸入。\n"
                    "- Raw read-AF 未校正 purity、CN 與 mutation multiplicity；softmax/argmax 不是 calibrated posterior。\n"
                    "- Official 與 normalized VAF top set 在可比較 ranked units 中僅 HCC 71.47%、DORADO 63.73% 相同，顯示排序對分母定義敏感。\n"
                    "- Candidate-set subset 描述不確定空間嵌套；只有 genomic site-set 也嵌套時，才稱 directional induced substructure。\n"
                    "- Site-label null 支持 projected candidate relation 高於任意標籤基準，但不是 fixed strict-exact+true-induced 的直接 null。"
                ),
            },
            {
                "id": "claim_heading",
                "type": "markdown",
                "body": (
                    "## Claim ceiling：可說「部分結構現象可在兩技術來源重現」，不可說「方法已證明準確有效」\n\n"
                    "最強安全說法仍受三個關卡限制：fixed structural endpoint 的直接 null、clean-v3 freeze 與 external biological truth。"
                    "除非三者補齊，不將技術相容詞寫成 clone-tree accuracy。"
                ),
                "sourceId": "method_contract",
            },
            {"id": "claim_table", "type": "table", "tableId": "claim_table"},
            {
                "id": "annotation_context_heading",
                "type": "markdown",
                "body": (
                    "## 癌症基因與藥物註記沒有提供額外的 topology validation\n\n"
                    "先前粗拓撲 exact-coordinate complete-both 配對的一致為 3,969/5,720 = **69.39%**，Cohen's κ=**0.4973**；"
                    "這與本報告的逐位點／候選集 endpoint 不同，不能互換。"
                    "其 chromosome + region-length-decile matched-null 結果顯示："
                    "COSMIC CGC gene-body stratum 的 present-vs-absent difference 為 **3.15 pp**（p=**0.5855**）；"
                    "DGIdb approved-antineoplastic stratum 為 **3.11 pp**（p=**0.4295**）。"
                    "兩者都不顯著，因此只能作 biological context / face validity，"
                    "不能證明逐位點拓撲正確，也不是用藥證據。"
                ),
                "sourceId": "annotation_context",
            },
            {
                "id": "annotation_context_table",
                "type": "table",
                "tableId": "annotation_context_table",
            },
            {
                "id": "validation_heading",
                "type": "markdown",
                "body": (
                    "## 驗證：三層 checks、重跑與逐位點 allele identity 全部通過\n\n"
                    "Core producer 18/18、postprocess 23/23、null 25/25 checks PASS；"
                    "5,720 fixed outcomes 各層回加守恆、15,713/15,713 共享 allele identity，"
                    "核心與 null 同參數重跑均 deterministic。"
                ),
                "sourceId": "site_checks",
            },
            {"id": "checks_table", "type": "table", "tableId": "checks_table"},
            {
                "id": "next_steps",
                "type": "markdown",
                "body": (
                    "## 下一步：將描述性技術再現升級為可驗證科學證據\n\n"
                    "1. 對 fixed strict-full-exact+true-induced endpoint 建立直接、locked 的 conditional null。\n"
                    "2. 保留現有 chromosome-block bootstrap/LOCO，另完成 depth-matched within-dataset split-half。\n"
                    "3. 在 clean layered-v3 closeout 後全量重算，不以 historical-v2 做 release claim。\n"
                    "4. 以 single-cell DNA、colony truth 或 synthetic spike-in 分開「可重現」與「正確」。"
                ),
            },
            {
                "id": "further_questions",
                "type": "markdown",
                "body": (
                    "## 尚待回答的問題\n\n"
                    "- 移除 low-k、Direct-heavy 與 large-candidate strata 後，robust excess 是否仍存在？\n"
                    "- HP component mismatch 是 biological LOH/split 現象，還是 caller/basecaller 差異？\n"
                    "- 加入 CN/purity/multiplicity 校正後，VAF top-set 是否仍比 read-full 更集中且不增加 conflict？"
                ),
            },
        ]
    )

    access_issues: list[dict[str, Any]] = []
    if not null_available:
        access_issues.append(
            {
                "id": "conditional_null_pending",
                "scope": "report",
                "dataset": "conditional_null",
                "message": (
                    "Complexity-matched conditional null is not present in the analysis output; "
                    "chance-corrected enrichment and proof-of-effectiveness remain unavailable."
                ),
            }
        )

    # Native cards/charts/tables must resolve to the literal SQL snapshot that
    # actually materialized their reviewed rows. Markdown blocks keep their
    # upstream file/method provenance because they do not expose widget data.
    for card in cards:
        card["sourceId"] = f"dataset_{card['dataset']}"
    for chart in charts:
        chart["sourceId"] = f"dataset_{chart['dataset']}"
    for table in tables:
        table["sourceId"] = f"dataset_{table['dataset']}"

    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": "HCC1395 逐位點拓撲子結構跨來源驗證",
            "description": (
                "5,720 exact-coordinate complete-both HCC1395 technical region pairs; "
                "read-full versus VAF-top shared-site topology equality, nesting, overlap, conflict, complexity strata, and claim ceiling."
            ),
            "generatedAt": generated_at,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": sources,
            "blocks": blocks,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "partial",
            "datasets": datasets,
            "accessIssues": access_issues,
        },
        "sources": sources,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    artifact_path = output_dir / "artifact.json"
    artifact_path.write_text(
        json.dumps(artifact, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    source_notes = {
        "generated_at": generated_at,
        "task_type": "B comprehensive validation",
        "audience": "technical",
        "delivery_mode": "portable HTML via canonical artifact",
        "fixed_denominator": FIXED_DENOMINATOR,
        "null_status": "available" if null_available else "pending",
        "input_files": [
            {"path": path.as_posix(), "sha256": sha256(path)} for path in paths.values()
        ]
        + [
            {"path": ANNOTATION_CONTEXT_PATH.as_posix(), "sha256": sha256(ANNOTATION_CONTEXT_PATH)},
            {"path": HCC1395_KB_PATH.as_posix(), "sha256": sha256(HCC1395_KB_PATH)},
        ],
        "chart_map": [
            {
                "chart": chart["id"],
                "question": chart["question"],
                "type": chart["type"],
                "dataset": chart["dataset"],
                "claim": chart["comparisonContext"],
            }
            for chart in charts
        ],
        "omissions": [
            "Conditional-null chart omitted because the required matched-null TSV was not produced."
        ]
        if not null_available
        else [],
        "claim_ceiling": "PARTIAL technical reproducibility; SCIENTIFIC NO-GO for accuracy/effectiveness.",
    }
    notes_path = output_dir / "source_notes.json"
    notes_path.write_text(
        json.dumps(source_notes, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"INPUT data dir: {data_dir.resolve()}")
    print(f"INPUT regions: {paths['regions'].resolve()}")
    print(f"INPUT metrics: {paths['metrics'].resolve()}")
    print(f"INPUT checks: {paths['checks'].resolve()}")
    print(f"INPUT summary: {paths['summary'].resolve()}")
    print(f"OUTPUT artifact: {artifact_path.resolve()}")
    print(f"OUTPUT source notes: {notes_path.resolve()}")
    print(
        "RESULT "
        f"regions={len(regions)} report_checks={len(report_checks)} "
        f"charts={len(charts)} tables={len(tables)} null={'available' if null_available else 'pending'}"
    )


if __name__ == "__main__":
    main()
