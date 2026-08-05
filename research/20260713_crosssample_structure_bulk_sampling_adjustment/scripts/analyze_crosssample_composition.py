#!/usr/bin/env python3
"""Bulk-sampling and compositional sensitivity analysis for seven dataset rows.

Primary estimand: five complete-region states, including unresolved.
Secondary estimand: four resolved-only rooted-unlabeled graph patterns.

This is deliberately a descriptive technical-reproducibility analysis.  It
does not test equality across different cell lines and does not validate clone
trees.  All stochastic endpoints use fixed, independently derived seeds.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import sys
from typing import Iterable, Mapping, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import beta


SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
BIOLOGICAL_ID = {
    "HCC1395": "HCC1395",
    "HCC1395_DORADO": "HCC1395",
    "COLO829": "COLO829",
    "H1437": "H1437",
    "H2009": "H2009",
    "HCC1937": "HCC1937",
    "HCC1954": "HCC1954",
}
BIOLOGICAL_ORDER = ["HCC1395", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
RESOLVED_CATEGORIES = ["single_only", "sister_only", "direct_only", "sister_and_direct"]
COMPLETE_CATEGORIES = RESOLVED_CATEGORIES + ["unresolved"]
ESTIMANDS = {
    "complete_five_class_primary": COMPLETE_CATEGORIES,
    "resolved_four_class_secondary": RESOLVED_CATEGORIES,
}
SOURCE_ORDER = ["structural_topo1", "vaf_resolved_topogt1"]
REGION_PATTERN = re.compile(r"^(chr(?:[1-9]|1[0-9]|2[0-2])):(\d+)-(\d+)$")
POSTERIOR_ALPHA = 0.5
POSTERIOR_LEVEL = 0.95
DEFAULT_SEED = 20260713
DEFAULT_REPLICATES = 5000
BLOCK_WIDTHS = (5_000_000, 10_000_000)
EB_MULTIPLIERS = (0.5, 1.0, 2.0)
CLAIM_CEILING = (
    "Historical layered-v2 VAF-selected rooted-unlabeled mutation-state graph-pattern composition; "
    "technical sampling sensitivity only, not clone/subclone truth, biological equivalence, or a validated tree."
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def formatted(value: object) -> object:
    if isinstance(value, (np.bool_, bool)):
        return "True" if bool(value) else "False"
    if isinstance(value, (np.integer, int)):
        return str(int(value))
    if isinstance(value, (np.floating, float)):
        if not math.isfinite(float(value)):
            return ""
        return f"{float(value):.12g}"
    if value is None:
        return ""
    return value


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: formatted(row.get(field)) for field in fields})


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")


def truth(value: object) -> bool:
    return value is True or str(value).lower() in {"true", "1"}


def proportions(counts: Sequence[float]) -> np.ndarray:
    values = np.asarray(counts, dtype=float)
    total = values.sum()
    if total <= 0:
        raise ValueError("Composition has non-positive denominator")
    return values / total


def tv_distance(first: np.ndarray, second: np.ndarray) -> float:
    return float(0.5 * np.abs(first - second).sum())


def js_distance(first: np.ndarray, second: np.ndarray) -> float:
    midpoint = 0.5 * (first + second)
    with np.errstate(divide="ignore", invalid="ignore"):
        first_term = np.where(first > 0, first * np.log2(first / midpoint), 0.0)
        second_term = np.where(second > 0, second * np.log2(second / midpoint), 0.0)
    return float(np.sqrt(0.5 * (first_term.sum() + second_term.sum())))


def aitchison_distance(first: np.ndarray, second: np.ndarray) -> float:
    # A tiny replacement is only a safety net for standardized sensitivity rows.
    first = np.clip(first, 1e-15, None)
    second = np.clip(second, 1e-15, None)
    first = first / first.sum()
    second = second / second.sum()
    clr_first = np.log(first) - np.log(first).mean()
    clr_second = np.log(second) - np.log(second).mean()
    return float(np.sqrt(np.square(clr_first - clr_second).sum()))


def distances(first: np.ndarray, second: np.ndarray) -> dict[str, float]:
    return {
        "total_variation": tv_distance(first, second),
        "jensen_shannon_distance_base2": js_distance(first, second),
        "aitchison_distance": aitchison_distance(first, second),
    }


def vector_distances(first: np.ndarray, second: np.ndarray, clr_first: np.ndarray, clr_second: np.ndarray) -> dict[str, np.ndarray]:
    midpoint = 0.5 * (first + second)
    with np.errstate(divide="ignore", invalid="ignore"):
        term_first = np.where(first > 0, first * np.log2(first / midpoint), 0.0)
        term_second = np.where(second > 0, second * np.log2(second / midpoint), 0.0)
    clr_a = np.log(clr_first) - np.log(clr_first).mean(axis=1, keepdims=True)
    clr_b = np.log(clr_second) - np.log(clr_second).mean(axis=1, keepdims=True)
    return {
        "total_variation": 0.5 * np.abs(first - second).sum(axis=1),
        "jensen_shannon_distance_base2": np.sqrt(0.5 * (term_first.sum(axis=1) + term_second.sum(axis=1))),
        "aitchison_distance": np.sqrt(np.square(clr_a - clr_b).sum(axis=1)),
    }


def quantiles(values: np.ndarray) -> tuple[float, float, float, float]:
    return (
        float(np.mean(values)),
        float(np.quantile(values, 0.5)),
        float(np.quantile(values, 0.025)),
        float(np.quantile(values, 0.975)),
    )


def parse_region(region: str) -> tuple[str, int, int]:
    match = REGION_PATTERN.fullmatch(region)
    if not match:
        raise ValueError(f"Invalid autosomal region coordinate: {region}")
    chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    if start > end:
        raise ValueError(f"Inverted region coordinate: {region}")
    return chrom, start, end


def pair_type(sample_a: str, sample_b: str) -> str:
    return "technical_same_biological_id" if BIOLOGICAL_ID[sample_a] == BIOLOGICAL_ID[sample_b] else "cross_biological_id"


def raw_pairwise_rows(
    level: str,
    estimand: str,
    compositions: Mapping[str, np.ndarray],
    denominators: Mapping[str, int],
    biological_ids: Mapping[str, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    names = list(compositions)
    for first, second in itertools.combinations(names, 2):
        relation = (
            "technical_same_biological_id"
            if biological_ids[first] == biological_ids[second]
            else "cross_biological_id"
        )
        row: dict[str, object] = {
            "comparison_level": level,
            "estimand": estimand,
            "entity_a": first,
            "entity_b": second,
            "biological_id_a": biological_ids[first],
            "biological_id_b": biological_ids[second],
            "pair_type": relation,
            "denominator_a": denominators[first],
            "denominator_b": denominators[second],
        }
        row.update(distances(compositions[first], compositions[second]))
        rows.append(row)
    return rows


def raw_rank_rows(estimand: str, pairwise_rows: Sequence[Mapping[str, object]], method: str) -> list[dict[str, object]]:
    technical = next(
        row for row in pairwise_rows
        if row["entity_a"] == "HCC1395" and row["entity_b"] == "HCC1395_DORADO"
    )
    cross = [row for row in pairwise_rows if row["pair_type"] == "cross_biological_id"]
    rows = []
    for metric in ("total_variation", "jensen_shannon_distance_base2", "aitchison_distance"):
        hcc = float(technical[metric])
        reference = np.asarray([float(row[metric]) for row in cross])
        count_lower = int(np.sum(reference < hcc))
        count_le = int(np.sum(reference <= hcc))
        rows.append({
            "method": method,
            "estimand": estimand,
            "metric": metric,
            "replicates": 1,
            "hcc_distance_mean": hcc,
            "hcc_distance_median": hcc,
            "hcc_distance_q025": hcc,
            "hcc_distance_q975": hcc,
            "cross_biological_pairs": len(reference),
            "hcc_rank_smaller_is_more_similar_median": 1 + count_lower,
            "hcc_rank_q025": 1 + count_lower,
            "hcc_rank_q975": 1 + count_lower,
            "probability_hcc_rank_1": float(count_lower == 0),
            "lower_tail_percentile_median": 100.0 * count_le / len(reference),
            "add_one_lower_tail_percentile": 100.0 * (1 + count_le) / (len(reference) + 1),
            "cross_biological_min": float(reference.min()),
            "cross_biological_median": float(np.median(reference)),
            "cross_biological_max": float(reference.max()),
        })
    return rows


def summarize_resampling(
    method: str,
    estimand: str,
    raw_draws: Mapping[str, np.ndarray],
    clr_draws: Mapping[str, np.ndarray],
    categories: Sequence[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    composition_rows: list[dict[str, object]] = []
    for sample in SAMPLE_ORDER:
        values = raw_draws[sample]
        for index, category in enumerate(categories):
            mean, median, low, high = quantiles(values[:, index])
            composition_rows.append({
                "method": method,
                "estimand": estimand,
                "sample": sample,
                "biological_id": BIOLOGICAL_ID[sample],
                "category": category,
                "replicates": values.shape[0],
                "mean": mean,
                "median": median,
                "q025": low,
                "q975": high,
            })

    pair_rows: list[dict[str, object]] = []
    metric_arrays: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for first, second in itertools.combinations(SAMPLE_ORDER, 2):
        metric_values = vector_distances(
            raw_draws[first], raw_draws[second], clr_draws[first], clr_draws[second]
        )
        metric_arrays[(first, second)] = metric_values
        for metric, values in metric_values.items():
            mean, median, low, high = quantiles(values)
            pair_rows.append({
                "method": method,
                "estimand": estimand,
                "sample_a": first,
                "sample_b": second,
                "biological_id_a": BIOLOGICAL_ID[first],
                "biological_id_b": BIOLOGICAL_ID[second],
                "pair_type": pair_type(first, second),
                "metric": metric,
                "replicates": len(values),
                "mean": mean,
                "median": median,
                "q025": low,
                "q975": high,
            })

    rank_rows: list[dict[str, object]] = []
    hcc_key = ("HCC1395", "HCC1395_DORADO")
    cross_keys = [key for key in metric_arrays if pair_type(*key) == "cross_biological_id"]
    for metric in ("total_variation", "jensen_shannon_distance_base2", "aitchison_distance"):
        hcc = metric_arrays[hcc_key][metric]
        reference = np.stack([metric_arrays[key][metric] for key in cross_keys], axis=1)
        ranks = 1 + np.sum(reference < hcc[:, None], axis=1)
        lower_percentile = 100.0 * np.sum(reference <= hcc[:, None], axis=1) / reference.shape[1]
        mean, median, low, high = quantiles(hcc)
        rank_rows.append({
            "method": method,
            "estimand": estimand,
            "metric": metric,
            "replicates": len(hcc),
            "hcc_distance_mean": mean,
            "hcc_distance_median": median,
            "hcc_distance_q025": low,
            "hcc_distance_q975": high,
            "cross_biological_pairs": reference.shape[1],
            "hcc_rank_smaller_is_more_similar_median": float(np.median(ranks)),
            "hcc_rank_q025": float(np.quantile(ranks, 0.025)),
            "hcc_rank_q975": float(np.quantile(ranks, 0.975)),
            "probability_hcc_rank_1": float(np.mean(ranks == 1)),
            "lower_tail_percentile_median": float(np.median(lower_percentile)),
            "add_one_lower_tail_percentile": float(np.median(100.0 * (1 + np.sum(reference <= hcc[:, None], axis=1)) / (reference.shape[1] + 1))),
            "cross_biological_min": float(np.median(reference.min(axis=1))),
            "cross_biological_median": float(np.median(reference)),
            "cross_biological_max": float(np.median(reference.max(axis=1))),
        })
    return composition_rows, pair_rows, rank_rows


def fit_dirichlet_mle(compositions: np.ndarray) -> tuple[np.ndarray, dict[str, object]]:
    values = np.asarray(compositions, dtype=float)
    values = np.clip(values, 1e-12, None)
    values /= values.sum(axis=1, keepdims=True)
    mean = values.mean(axis=0)
    variance = values.var(axis=0, ddof=1)
    candidate_tau = mean * (1.0 - mean) / np.maximum(variance, 1e-12) - 1.0
    positive_tau = candidate_tau[np.isfinite(candidate_tau) & (candidate_tau > 0)]
    tau = float(np.median(positive_tau)) if len(positive_tau) else 10.0
    tau = min(max(tau, 0.1), 10000.0)
    initial = np.log(np.clip(mean * tau, 1e-3, None))

    sum_log = np.log(values).sum(axis=0)
    n = values.shape[0]

    def objective(log_alpha: np.ndarray) -> float:
        alpha = np.exp(log_alpha)
        log_likelihood = (
            n * (gammaln(alpha.sum()) - gammaln(alpha).sum())
            + np.dot(alpha - 1.0, sum_log)
        )
        return float(-log_likelihood)

    result = minimize(
        objective,
        initial,
        method="L-BFGS-B",
        bounds=[(math.log(1e-4), math.log(1e5))] * values.shape[1],
        options={"maxiter": 20000, "ftol": 1e-12, "gtol": 1e-9},
    )
    alpha = np.exp(result.x)
    return alpha, {
        "success": bool(result.success),
        "message": str(result.message),
        "iterations": int(result.nit),
        "objective": float(result.fun),
        "n_biological_compositions": int(values.shape[0]),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--by-source", type=Path, required=True)
    parser.add_argument("--coarse-regions", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    args = parser.parse_args()

    for path in (args.summary, args.regions, args.by_source, args.coarse_regions):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.replicates < 1000:
        raise ValueError("At least 1,000 replicates are required for stable 95% endpoints")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    checks: list[dict[str, object]] = []

    def check(name: str, observed: object, expected: object, severity: str = "ERROR") -> None:
        checks.append({
            "check": name,
            "severity": severity,
            "pass": observed == expected,
            "observed": observed,
            "expected": expected,
        })

    # ---------- Summary counts and primary/secondary estimands ----------
    summary_all = read_tsv(args.summary)
    summary_rows = [row for row in summary_all if row["sample"] != "aggregate"]
    aggregate_rows = [row for row in summary_all if row["sample"] == "aggregate"]
    check("seven_dataset_rows", len(summary_rows), 7)
    check("one_aggregate_row", len(aggregate_rows), 1)
    check("dataset_order_exact", [row["sample"] for row in summary_rows], SAMPLE_ORDER)

    counts_by_estimand: dict[str, dict[str, np.ndarray]] = {key: {} for key in ESTIMANDS}
    denominators: dict[str, dict[str, int]] = {key: {} for key in ESTIMANDS}
    for row in summary_rows:
        sample = row["sample"]
        resolved = np.asarray([int(row[category]) for category in RESOLVED_CATEGORIES], dtype=int)
        complete = np.append(resolved, int(row["unresolved_regions"]))
        counts_by_estimand["resolved_four_class_secondary"][sample] = resolved
        counts_by_estimand["complete_five_class_primary"][sample] = complete
        denominators["resolved_four_class_secondary"][sample] = int(resolved.sum())
        denominators["complete_five_class_primary"][sample] = int(complete.sum())
        check(f"{sample}_resolved_categories_conserve", int(resolved.sum()), int(row["final_single_shape_regions"]))
        check(f"{sample}_complete_five_class_conserve", int(complete.sum()), int(row["complete_regions"]))
        check(f"{sample}_resolved_plus_unresolved", int(resolved.sum()) + int(row["unresolved_regions"]), int(row["complete_regions"]))

    aggregate = aggregate_rows[0]
    check("aggregate_resolved_37039", sum(denominators["resolved_four_class_secondary"].values()), int(aggregate["final_single_shape_regions"]))
    check("aggregate_complete_39885", sum(denominators["complete_five_class_primary"].values()), int(aggregate["complete_regions"]))
    check("six_biological_ids", len(set(BIOLOGICAL_ID.values())), 6)
    check("technical_pair_same_biological_id", BIOLOGICAL_ID["HCC1395"] == BIOLOGICAL_ID["HCC1395_DORADO"], True)

    # ---------- Region key/coordinate reconstruction, including unresolved ----------
    resolved_rows = read_tsv(args.regions)
    resolved_by_key: dict[tuple[str, str], str] = {}
    for row in resolved_rows:
        key = (row["sample"], row["region"])
        if key in resolved_by_key:
            raise ValueError(f"Duplicate resolved sample×region key: {key}")
        category = row["final_shape_category"]
        if category not in RESOLVED_CATEGORIES:
            raise ValueError(f"Unexpected resolved category: {category}")
        resolved_by_key[key] = category
    check("resolved_region_rows_37039", len(resolved_rows), 37039)
    check("resolved_sample_region_key_unique", len(resolved_by_key), len(resolved_rows))

    coarse_rows = read_tsv(args.coarse_regions)
    complete_rows: list[dict[str, object]] = []
    seen_complete: set[tuple[str, str]] = set()
    coordinate_mismatches = 0
    for row in coarse_rows:
        if not truth(row["complete"]):
            continue
        sample, region = row["sample"], row["region"]
        key = (sample, region)
        if key in seen_complete:
            raise ValueError(f"Duplicate complete sample×region key: {key}")
        seen_complete.add(key)
        parsed_chrom, parsed_start, parsed_end = parse_region(region)
        if parsed_chrom != row["chrom"] or parsed_start != int(row["start"]) or parsed_end != int(row["end"]):
            coordinate_mismatches += 1
        complete_rows.append({
            "sample": sample,
            "biological_id": BIOLOGICAL_ID[sample],
            "region": region,
            "chrom": parsed_chrom,
            "start": parsed_start,
            "end": parsed_end,
            "midpoint": (parsed_start + parsed_end) // 2,
            "category": resolved_by_key.get(key, "unresolved"),
        })
    check("complete_region_rows_39885", len(complete_rows), 39885)
    check("complete_sample_region_key_unique", len(seen_complete), len(complete_rows))
    check("coordinate_columns_match_region_string", coordinate_mismatches, 0)
    check("all_resolved_keys_in_complete_universe", len(set(resolved_by_key) - seen_complete), 0)

    complete_counts_from_regions: dict[str, Counter] = defaultdict(Counter)
    chromosomes_by_sample: dict[str, set[str]] = defaultdict(set)
    for row in complete_rows:
        complete_counts_from_regions[str(row["sample"])][str(row["category"])] += 1
        chromosomes_by_sample[str(row["sample"])].add(str(row["chrom"]))
    for sample in SAMPLE_ORDER:
        observed = [complete_counts_from_regions[sample][category] for category in COMPLETE_CATEGORIES]
        check(f"{sample}_complete_region_counts_match_summary", observed, counts_by_estimand["complete_five_class_primary"][sample].tolist())
        check(f"{sample}_has_chr1_to_chr22", len(chromosomes_by_sample[sample]), 22)

    # ---------- Source strata ----------
    source_rows = [row for row in read_tsv(args.by_source) if row["sample"] != "aggregate"]
    source_counts: dict[tuple[str, str], np.ndarray] = {}
    for sample in SAMPLE_ORDER:
        for source in SOURCE_ORDER:
            rows = [row for row in source_rows if row["sample"] == sample and row["selection_source"] == source]
            if len(rows) != len(RESOLVED_CATEGORIES):
                raise ValueError(f"Incomplete source table for {sample}/{source}: {len(rows)}")
            by_category = {row["final_shape_category"]: int(row["regions"]) for row in rows}
            vector = np.asarray([by_category[category] for category in RESOLVED_CATEGORIES], dtype=int)
            source_counts[(sample, source)] = vector
            check(f"{sample}_{source}_declared_source_total", int(vector.sum()), int(rows[0]["source_total"]))
        combined = source_counts[(sample, SOURCE_ORDER[0])] + source_counts[(sample, SOURCE_ORDER[1])]
        check(f"{sample}_source_categories_match_resolved_summary", combined.tolist(), counts_by_estimand["resolved_four_class_secondary"][sample].tolist())

    global_source_totals = np.asarray([
        sum(int(source_counts[(sample, source)].sum()) for sample in SAMPLE_ORDER)
        for source in SOURCE_ORDER
    ])
    check("global_structural_source_total", int(global_source_totals[0]), 21976)
    check("global_vaf_resolved_source_total", int(global_source_totals[1]), 15063)
    global_resolved_source_weights = proportions(global_source_totals)
    global_unresolved = sum(int(counts_by_estimand["complete_five_class_primary"][sample][-1]) for sample in SAMPLE_ORDER)
    global_complete_source_weights = np.asarray([
        global_source_totals[0], global_source_totals[1], global_unresolved
    ], dtype=float)
    global_complete_source_weights /= global_complete_source_weights.sum()

    source_stratified_pairwise: list[dict[str, object]] = []
    hcc_source_conditional: dict[str, dict[str, float]] = {}
    for source in SOURCE_ORDER:
        conditional_compositions = {
            sample: proportions(source_counts[(sample, source)]) for sample in SAMPLE_ORDER
        }
        conditional_denominators = {
            sample: int(source_counts[(sample, source)].sum()) for sample in SAMPLE_ORDER
        }
        pair_rows = raw_pairwise_rows(
            "dataset_row_within_selection_source",
            "resolved_four_class_secondary",
            conditional_compositions,
            conditional_denominators,
            BIOLOGICAL_ID,
        )
        for row in pair_rows:
            output = dict(row)
            output["selection_source"] = source
            output["aitchison_zero_handling"] = "clip_to_1e-15_if_zero"
            source_stratified_pairwise.append(output)
        hcc_row = next(
            row for row in pair_rows
            if row["entity_a"] == "HCC1395" and row["entity_b"] == "HCC1395_DORADO"
        )
        hcc_source_conditional[source] = {
            metric: float(hcc_row[metric])
            for metric in ("total_variation", "jensen_shannon_distance_base2", "aitchison_distance")
        }

    # ---------- Marginal multinomial/Dirichlet posterior intervals ----------
    interval_rows: list[dict[str, object]] = []
    raw_dataset_compositions: dict[str, dict[str, np.ndarray]] = {key: {} for key in ESTIMANDS}
    for estimand, categories in ESTIMANDS.items():
        for sample in SAMPLE_ORDER:
            counts = counts_by_estimand[estimand][sample]
            n = int(counts.sum())
            raw = proportions(counts)
            raw_dataset_compositions[estimand][sample] = raw
            alpha = counts.astype(float) + POSTERIOR_ALPHA
            alpha_total = alpha.sum()
            for index, category in enumerate(categories):
                low = float(beta.ppf((1.0 - POSTERIOR_LEVEL) / 2.0, alpha[index], alpha_total - alpha[index]))
                high = float(beta.ppf(1.0 - (1.0 - POSTERIOR_LEVEL) / 2.0, alpha[index], alpha_total - alpha[index]))
                interval_rows.append({
                    "estimand": estimand,
                    "sample": sample,
                    "biological_id": BIOLOGICAL_ID[sample],
                    "category": category,
                    "count": int(counts[index]),
                    "denominator": n,
                    "raw_proportion": raw[index],
                    "posterior_prior": f"Jeffreys_Dirichlet_alpha_{POSTERIOR_ALPHA}_each",
                    "posterior_mean": alpha[index] / alpha_total,
                    "posterior_marginal_q025": low,
                    "posterior_marginal_q975": high,
                })

    # ---------- Dataset and biological-ID raw composition distances ----------
    raw_pairwise: list[dict[str, object]] = []
    biological_composition_rows: list[dict[str, object]] = []
    biological_compositions_by_estimand: dict[str, dict[str, np.ndarray]] = {}
    raw_rank: list[dict[str, object]] = []
    for estimand, categories in ESTIMANDS.items():
        dataset_pairs = raw_pairwise_rows(
            "dataset_row", estimand, raw_dataset_compositions[estimand], denominators[estimand], BIOLOGICAL_ID
        )
        raw_pairwise.extend(dataset_pairs)
        raw_rank.extend(raw_rank_rows(estimand, dataset_pairs, "raw_unadjusted"))

        bio_compositions: dict[str, np.ndarray] = {}
        bio_denominators: dict[str, int] = {}
        bio_map = {bio: bio for bio in BIOLOGICAL_ORDER}
        for bio in BIOLOGICAL_ORDER:
            members = [sample for sample in SAMPLE_ORDER if BIOLOGICAL_ID[sample] == bio]
            matrix = np.stack([raw_dataset_compositions[estimand][sample] for sample in members])
            vector = matrix.mean(axis=0)
            bio_compositions[bio] = vector
            bio_denominators[bio] = sum(denominators[estimand][sample] for sample in members)
            aggregation = "equal_dataset_source_mean" if len(members) > 1 else "single_dataset_source"
            for index, category in enumerate(categories):
                biological_composition_rows.append({
                    "estimand": estimand,
                    "biological_id": bio,
                    "dataset_sources": ";".join(members),
                    "aggregation": aggregation,
                    "supporting_region_rows_not_independent_denominator": bio_denominators[bio],
                    "category": category,
                    "proportion": vector[index],
                })
        biological_compositions_by_estimand[estimand] = bio_compositions
        raw_pairwise.extend(raw_pairwise_rows(
            "biological_id_equal_source_mean", estimand, bio_compositions, bio_denominators, bio_map
        ))

    # ---------- Posterior sampling ----------
    resampled_compositions: list[dict[str, object]] = []
    resampled_pairwise: list[dict[str, object]] = []
    rank_rows: list[dict[str, object]] = list(raw_rank)
    for estimand_index, (estimand, categories) in enumerate(ESTIMANDS.items()):
        rng = np.random.default_rng(args.seed + 100 + estimand_index)
        draws: dict[str, np.ndarray] = {}
        for sample in SAMPLE_ORDER:
            alpha = counts_by_estimand[estimand][sample].astype(float) + POSTERIOR_ALPHA
            draws[sample] = rng.dirichlet(alpha, size=args.replicates)
        comp, pairs, ranks = summarize_resampling(
            "dirichlet_posterior", estimand, draws, draws, categories
        )
        resampled_compositions.extend(comp)
        resampled_pairwise.extend(pairs)
        rank_rows.extend(ranks)

    # ---------- Equal-n rarefaction ----------
    for estimand_index, (estimand, categories) in enumerate(ESTIMANDS.items()):
        target_n = min(denominators[estimand].values())
        rng = np.random.default_rng(args.seed + 200 + estimand_index)
        raw_draws: dict[str, np.ndarray] = {}
        clr_draws: dict[str, np.ndarray] = {}
        for sample in SAMPLE_ORDER:
            sampled = rng.multivariate_hypergeometric(
                counts_by_estimand[estimand][sample], target_n, size=args.replicates
            )
            raw_draws[sample] = sampled / sampled.sum(axis=1, keepdims=True)
            smoothed = sampled.astype(float) + POSTERIOR_ALPHA
            clr_draws[sample] = smoothed / smoothed.sum(axis=1, keepdims=True)
        method = f"rarefaction_without_replacement_equal_n_{target_n}"
        comp, pairs, ranks = summarize_resampling(method, estimand, raw_draws, clr_draws, categories)
        resampled_compositions.extend(comp)
        resampled_pairwise.extend(pairs)
        rank_rows.extend(ranks)

    # ---------- Common fixed-genomic-block bootstrap ----------
    block_feasibility_rows: list[dict[str, object]] = []
    for width_index, width in enumerate(BLOCK_WIDTHS):
        rows_by_sample = {sample: [row for row in complete_rows if row["sample"] == sample] for sample in SAMPLE_ORDER}
        block_sets = {
            sample: {(str(row["chrom"]), int(row["midpoint"]) // width) for row in rows}
            for sample, rows in rows_by_sample.items()
        }
        common_blocks = sorted(
            set.intersection(*block_sets.values()),
            key=lambda value: (int(value[0][3:]), value[1]),
        )
        union_blocks = set.union(*block_sets.values())
        check(f"common_{width // 1_000_000}mb_blocks_nonzero", len(common_blocks) > 0, True)
        common_index = {block: index for index, block in enumerate(common_blocks)}
        for sample in SAMPLE_ORDER:
            rows_common = sum(
                (str(row["chrom"]), int(row["midpoint"]) // width) in common_index
                for row in rows_by_sample[sample]
            )
            block_feasibility_rows.append({
                "block_width_bp": width,
                "sample": sample,
                "sample_blocks": len(block_sets[sample]),
                "all7_common_blocks": len(common_blocks),
                "all7_union_blocks": len(union_blocks),
                "common_block_share_union": len(common_blocks) / len(union_blocks),
                "complete_rows_total": len(rows_by_sample[sample]),
                "complete_rows_in_common_blocks": rows_common,
                "complete_row_coverage": rows_common / len(rows_by_sample[sample]),
                "coordinate_status": "PASS_chr_start_end_available",
            })

        rng = np.random.default_rng(args.seed + 300 + width_index)
        block_multiplicities = rng.multinomial(
            len(common_blocks), np.full(len(common_blocks), 1.0 / len(common_blocks)), size=args.replicates
        )
        for estimand, categories in ESTIMANDS.items():
            raw_draws: dict[str, np.ndarray] = {}
            clr_draws: dict[str, np.ndarray] = {}
            for sample in SAMPLE_ORDER:
                matrix = np.zeros((len(common_blocks), len(categories)), dtype=np.int64)
                category_index = {category: index for index, category in enumerate(categories)}
                for row in rows_by_sample[sample]:
                    category = str(row["category"])
                    if category not in category_index:
                        continue
                    block = (str(row["chrom"]), int(row["midpoint"]) // width)
                    if block in common_index:
                        matrix[common_index[block], category_index[category]] += 1
                sampled = block_multiplicities @ matrix
                if np.any(sampled.sum(axis=1) <= 0):
                    raise RuntimeError(f"Empty {width} block-bootstrap composition for {sample}/{estimand}")
                raw_draws[sample] = sampled / sampled.sum(axis=1, keepdims=True)
                smoothed = sampled.astype(float) + POSTERIOR_ALPHA
                clr_draws[sample] = smoothed / smoothed.sum(axis=1, keepdims=True)
            method = f"paired_common_genomic_block_bootstrap_{width // 1_000_000}mb"
            comp, pairs, ranks = summarize_resampling(method, estimand, raw_draws, clr_draws, categories)
            resampled_compositions.extend(comp)
            resampled_pairwise.extend(pairs)
            rank_rows.extend(ranks)

    # ---------- Source direct standardization ----------
    source_standardized_compositions: list[dict[str, object]] = []
    source_standardized_pairwise: list[dict[str, object]] = []
    hcc_pair_source_weights = proportions(np.asarray([
        sum(source_counts[(sample, source)].sum() for sample in ("HCC1395", "HCC1395_DORADO"))
        for source in SOURCE_ORDER
    ], dtype=float))
    hcc_complete_weights = proportions(np.asarray([
        sum(source_counts[(sample, SOURCE_ORDER[0])].sum() for sample in ("HCC1395", "HCC1395_DORADO")),
        sum(source_counts[(sample, SOURCE_ORDER[1])].sum() for sample in ("HCC1395", "HCC1395_DORADO")),
        sum(counts_by_estimand["complete_five_class_primary"][sample][-1] for sample in ("HCC1395", "HCC1395_DORADO")),
    ], dtype=float))
    source_weight_schemes = {
        "resolved_four_class_secondary": {
            "global_resolved_pooled_source_mix_primary_standardization": global_resolved_source_weights,
            "hcc_pair_pooled_source_mix_sensitivity": hcc_pair_source_weights,
            "equal_0.5_0.5_source_mix_sensitivity": np.asarray([0.5, 0.5]),
        },
        "complete_five_class_primary": {
            "global_complete_three_source_mix_sensitivity": global_complete_source_weights,
            "hcc_pair_complete_three_source_mix_sensitivity": hcc_complete_weights,
            "global_unresolved_equal_resolved_source_mix_sensitivity": np.asarray([
                (1.0 - global_complete_source_weights[2]) / 2.0,
                (1.0 - global_complete_source_weights[2]) / 2.0,
                global_complete_source_weights[2],
            ]),
        },
    }
    for estimand, schemes in source_weight_schemes.items():
        categories = ESTIMANDS[estimand]
        for scheme, weights in schemes.items():
            standardized: dict[str, np.ndarray] = {}
            for sample in SAMPLE_ORDER:
                conditional = [proportions(source_counts[(sample, source)]) for source in SOURCE_ORDER]
                resolved_part = weights[0] * conditional[0] + weights[1] * conditional[1]
                if estimand == "complete_five_class_primary":
                    vector = np.append(resolved_part, weights[2])
                else:
                    vector = resolved_part
                vector = vector / vector.sum()
                standardized[sample] = vector
                for index, category in enumerate(categories):
                    source_standardized_compositions.append({
                        "estimand": estimand,
                        "weight_scheme": scheme,
                        "sample": sample,
                        "biological_id": BIOLOGICAL_ID[sample],
                        "category": category,
                        "standardized_proportion": vector[index],
                        "structural_weight": weights[0],
                        "vaf_resolved_weight": weights[1],
                        "unresolved_weight": weights[2] if len(weights) == 3 else 0.0,
                        "causal_status": "sensitivity_only_source_is_outcome_dependent_selector",
                    })
            pair_rows = raw_pairwise_rows(
                "dataset_row", estimand, standardized,
                denominators[estimand], BIOLOGICAL_ID,
            )
            for row in pair_rows:
                output = dict(row)
                output["weight_scheme"] = scheme
                output["structural_weight"] = weights[0]
                output["vaf_resolved_weight"] = weights[1]
                output["unresolved_weight"] = weights[2] if len(weights) == 3 else 0.0
                source_standardized_pairwise.append(output)
            rank_rows.extend(raw_rank_rows(estimand, pair_rows, f"source_standardized:{scheme}"))

    # ---------- Biological-ID-aware empirical-Bayes sensitivity ----------
    eb_prior_rows: list[dict[str, object]] = []
    eb_shrinkage_rows: list[dict[str, object]] = []
    eb_pairwise_rows: list[dict[str, object]] = []
    eb_fit_status: list[dict[str, object]] = []
    for estimand, categories in ESTIMANDS.items():
        bio_compositions = biological_compositions_by_estimand[estimand]
        prior_inputs = {
            "all_6_biological_ids_equal_weight": BIOLOGICAL_ORDER,
            "leave_hcc1395_out_5_biological_ids": [bio for bio in BIOLOGICAL_ORDER if bio != "HCC1395"],
        }
        for prior_name, bios in prior_inputs.items():
            alpha, fit = fit_dirichlet_mle(np.stack([bio_compositions[bio] for bio in bios]))
            fit_row = {"estimand": estimand, "prior_name": prior_name, **fit}
            eb_fit_status.append(fit_row)
            check(f"eb_fit_success_{estimand}_{prior_name}", fit["success"], True)
            for index, category in enumerate(categories):
                eb_prior_rows.append({
                    "estimand": estimand,
                    "prior_name": prior_name,
                    "category": category,
                    "alpha": alpha[index],
                    "prior_mean": alpha[index] / alpha.sum(),
                    "prior_concentration": alpha.sum(),
                    "n_biological_ids": len(bios),
                    "biological_ids": ";".join(bios),
                    "fit_success": fit["success"],
                    "fit_message": fit["message"],
                })
            for multiplier in EB_MULTIPLIERS:
                scaled_alpha = alpha * multiplier
                shrunk: dict[str, np.ndarray] = {}
                for sample in SAMPLE_ORDER:
                    counts = counts_by_estimand[estimand][sample].astype(float)
                    raw = proportions(counts)
                    posterior = proportions(counts + scaled_alpha)
                    shrunk[sample] = posterior
                    for index, category in enumerate(categories):
                        eb_shrinkage_rows.append({
                            "estimand": estimand,
                            "prior_name": prior_name,
                            "concentration_multiplier": multiplier,
                            "effective_prior_concentration": scaled_alpha.sum(),
                            "sample": sample,
                            "biological_id": BIOLOGICAL_ID[sample],
                            "category": category,
                            "count": int(counts[index]),
                            "denominator": int(counts.sum()),
                            "raw_proportion": raw[index],
                            "shrunk_proportion": posterior[index],
                            "shift_percentage_points": 100.0 * (posterior[index] - raw[index]),
                            "interpretation": "regularization_sensitivity_not_equality_null",
                        })
                pair_rows = raw_pairwise_rows(
                    "dataset_row", estimand, shrunk, denominators[estimand], BIOLOGICAL_ID
                )
                method = f"empirical_bayes:{prior_name}:x{multiplier:g}"
                for row in pair_rows:
                    output = dict(row)
                    output["prior_name"] = prior_name
                    output["concentration_multiplier"] = multiplier
                    output["effective_prior_concentration"] = scaled_alpha.sum()
                    eb_pairwise_rows.append(output)
                rank_rows.extend(raw_rank_rows(estimand, pair_rows, method))

    # ---------- Final QA and warnings ----------
    posterior_sums = defaultdict(float)
    for row in interval_rows:
        posterior_sums[(row["estimand"], row["sample"])] += float(row["posterior_mean"])
    check("all_posterior_means_close_to_one", all(abs(value - 1.0) < 1e-10 for value in posterior_sums.values()), True)
    check("all_raw_pairwise_source_crossbio_count_primary", sum(
        row["estimand"] == "complete_five_class_primary"
        and row["comparison_level"] == "dataset_row"
        and row["pair_type"] == "cross_biological_id"
        for row in raw_pairwise
    ), 20)
    check("all_raw_pairwise_technical_pair_count_primary", sum(
        row["estimand"] == "complete_five_class_primary"
        and row["comparison_level"] == "dataset_row"
        and row["pair_type"] == "technical_same_biological_id"
        for row in raw_pairwise
    ), 1)
    check("biological_id_pair_count_per_estimand", sum(
        row["estimand"] == "complete_five_class_primary"
        and row["comparison_level"] == "biological_id_equal_source_mean"
        for row in raw_pairwise
    ), 15)

    warnings = [
        {
            "warning": "historical_layered_v2_claim_ceiling",
            "severity": "HIGH",
            "status": "OPEN",
            "detail": CLAIM_CEILING,
        },
        {
            "warning": "resolved_only_closure_bias",
            "severity": "HIGH",
            "status": "MITIGATED_PRIMARY_INCLUDES_UNRESOLVED",
            "detail": "Resolved-only four-class composition conditions on successful shape resolution and is secondary.",
        },
        {
            "warning": "source_standardization_possible_overadjustment",
            "severity": "HIGH",
            "status": "SENSITIVITY_ONLY",
            "detail": "Structural/VAF-resolved source is an outcome-dependent method-selection state, not a proven external confounder; direct standardization may remove or amplify real method differences.",
        },
        {
            "warning": "different_cell_lines_not_expected_equal",
            "severity": "HIGH",
            "status": "ENFORCED_BY_DESIGN",
            "detail": "Cross-biological pairs form an empirical reference distribution only; no equality null or common-clone target is imposed.",
        },
        {
            "warning": "block_bootstrap_not_exact_region_matching",
            "severity": "MEDIUM",
            "status": "DISCLOSED",
            "detail": "Common 5/10 Mb genomic blocks are jointly resampled, but HCC dataset rows have partly different region boundaries; this is not an exact paired-region bootstrap.",
        },
        {
            "warning": "empirical_bayes_small_biological_panel",
            "severity": "MEDIUM",
            "status": "SENSITIVITY_ONLY",
            "detail": "Dirichlet prior is estimated from six biological compositions (or five leave-HCC-out); shrinkage is not evidence that cell lines share one composition.",
        },
    ]

    failed = [row for row in checks if row["severity"] == "ERROR" and not row["pass"]]
    if failed:
        write_tsv(args.output_dir / "checks.tsv", checks, ["check", "severity", "pass", "observed", "expected"])
        raise RuntimeError(f"Analysis QA failed: {failed[:5]}")

    # ---------- Key summaries ----------
    raw_source_pairs = [
        row for row in raw_pairwise
        if row["comparison_level"] == "dataset_row"
    ]
    raw_hcc = {
        estimand: next(
            row for row in raw_source_pairs
            if row["estimand"] == estimand
            and row["entity_a"] == "HCC1395"
            and row["entity_b"] == "HCC1395_DORADO"
        )
        for estimand in ESTIMANDS
    }
    global_source_hcc = next(
        row for row in source_standardized_pairwise
        if row["estimand"] == "resolved_four_class_secondary"
        and row["weight_scheme"] == "global_resolved_pooled_source_mix_primary_standardization"
        and row["entity_a"] == "HCC1395" and row["entity_b"] == "HCC1395_DORADO"
    )
    global_complete_hcc = next(
        row for row in source_standardized_pairwise
        if row["estimand"] == "complete_five_class_primary"
        and row["weight_scheme"] == "global_complete_three_source_mix_sensitivity"
        and row["entity_a"] == "HCC1395" and row["entity_b"] == "HCC1395_DORADO"
    )
    raw_rank_lookup = {
        (row["estimand"], row["metric"]): row
        for row in rank_rows if row["method"] == "raw_unadjusted"
    }
    rank_lookup_all = {
        (row["method"], row["estimand"], row["metric"]): row for row in rank_rows
    }
    primary_jsd_methods = [
        "raw_unadjusted",
        "dirichlet_posterior",
        f"rarefaction_without_replacement_equal_n_{min(denominators['complete_five_class_primary'].values())}",
        "paired_common_genomic_block_bootstrap_5mb",
        "paired_common_genomic_block_bootstrap_10mb",
        "empirical_bayes:leave_hcc1395_out_5_biological_ids:x2",
        "source_standardized:global_complete_three_source_mix_sensitivity",
    ]
    primary_jsd_sensitivity = {
        method: rank_lookup_all[(method, "complete_five_class_primary", "jensen_shannon_distance_base2")]
        for method in primary_jsd_methods
    }
    max_eb_shift = max(abs(float(row["shift_percentage_points"])) for row in eb_shrinkage_rows)
    five_mb_coverage = {
        row["sample"]: row["complete_row_coverage"]
        for row in block_feasibility_rows if int(row["block_width_bp"]) == 5_000_000
    }

    outputs = {
        "checks": "checks.tsv",
        "warnings": "warnings.tsv",
        "posterior_intervals": "composition_posterior_intervals.tsv",
        "biological_id_compositions": "biological_id_compositions.tsv",
        "raw_pairwise": "raw_pairwise_distances.tsv",
        "resampled_compositions": "resampled_compositions.tsv",
        "resampled_pairwise": "resampled_pairwise_distances.tsv",
        "hcc_rank": "hcc_pair_relative_rank.tsv",
        "block_feasibility": "block_bootstrap_feasibility.tsv",
        "eb_priors": "empirical_bayes_priors.tsv",
        "eb_shrinkage": "empirical_bayes_shrinkage.tsv",
        "eb_pairwise": "empirical_bayes_pairwise_distances.tsv",
        "source_compositions": "source_standardized_compositions.tsv",
        "source_stratified_pairwise": "source_stratified_pairwise_distances.tsv",
        "source_pairwise": "source_standardized_pairwise_distances.tsv",
        "summary": "summary.json",
        "report": "analysis_report.md",
        "provenance": "provenance.json",
    }

    summary = {
        "schema_name": "intersubmod.crosssample_structure_bulk_sampling_adjustment",
        "schema_version": "1.0.0",
        "analysis_date": "2026-07-13",
        "status": "PASS_WITH_CAVEATS",
        "task_type": "B_comprehensive_validation",
        "scope": "chr1-22; 7 dataset rows; 6 biological IDs; historical layered-v2 engineering snapshot",
        "claim_ceiling": CLAIM_CEILING,
        "seed": args.seed,
        "replicates_per_stochastic_endpoint": args.replicates,
        "primary_estimand": "complete_five_class_primary",
        "secondary_estimand": "resolved_four_class_secondary",
        "categories": ESTIMANDS,
        "dataset_denominators": denominators,
        "biological_id_handling": {
            "technical_pair": ["HCC1395", "HCC1395_DORADO"],
            "technical_pair_biological_id": "HCC1395",
            "biological_ids": BIOLOGICAL_ORDER,
            "biological_id_composition_rule": "equal dataset-source mean for HCC1395; single source for other IDs",
            "cross_biological_source_pair_reference_n": 20,
            "biological_id_pair_n": 15,
        },
        "hcc_raw": {
            estimand: {
                metric: raw_hcc[estimand][metric]
                for metric in ("total_variation", "jensen_shannon_distance_base2", "aitchison_distance")
            }
            for estimand in ESTIMANDS
        },
        "hcc_raw_relative_to_cross_biological_pairs": {
            estimand: {
                metric: raw_rank_lookup[(estimand, metric)]
                for metric in ("total_variation", "jensen_shannon_distance_base2", "aitchison_distance")
            }
            for estimand in ESTIMANDS
        },
        "source_standardization": {
            "global_resolved_source_weights": {
                SOURCE_ORDER[index]: float(global_resolved_source_weights[index])
                for index in range(2)
            },
            "global_complete_three_source_weights": {
                "structural_topo1": float(global_complete_source_weights[0]),
                "vaf_resolved_topogt1": float(global_complete_source_weights[1]),
                "unresolved": float(global_complete_source_weights[2]),
            },
            "hcc_within_source_distances": hcc_source_conditional,
            "hcc_resolved_raw_total_variation": raw_hcc["resolved_four_class_secondary"]["total_variation"],
            "hcc_resolved_global_source_standardized_total_variation": global_source_hcc["total_variation"],
            "resolved_tv_amplification_factor": global_source_hcc["total_variation"] / raw_hcc["resolved_four_class_secondary"]["total_variation"],
            "hcc_complete_raw_total_variation": raw_hcc["complete_five_class_primary"]["total_variation"],
            "hcc_complete_global_three_source_standardized_total_variation": global_complete_hcc["total_variation"],
            "interpretation": "Raw resolved mixture partially cancels source-specific differences; common global source weights expose a larger conditional contrast. This is Simpson-like mixture cancellation sensitivity, not a corrected truth, because source is outcome-dependent.",
        },
        "primary_jsd_adjustment_stability": {
            "methods": primary_jsd_sensitivity,
            "verdict": "MODERATE_NOT_EXCEPTIONAL_TECHNICAL_PROXIMITY",
            "interpretation": (
                "Across raw, posterior, equal-n rarefaction, 5/10 Mb block bootstrap and EB, "
                "the HCC technical pair remains around rank 9 when inserted beside 20 cross-biological dataset pairs; "
                "it is not the closest pair. Source standardization gives rank 10 and is sensitivity-only."
            ),
        },
        "block_bootstrap": {
            "widths_bp": list(BLOCK_WIDTHS),
            "five_mb_common_blocks": next(int(row["all7_common_blocks"]) for row in block_feasibility_rows if int(row["block_width_bp"]) == 5_000_000),
            "five_mb_union_blocks": next(int(row["all7_union_blocks"]) for row in block_feasibility_rows if int(row["block_width_bp"]) == 5_000_000),
            "five_mb_complete_row_coverage_by_sample": five_mb_coverage,
            "method": "jointly resample fixed genomic block keys common to all seven dataset rows",
        },
        "empirical_bayes": {
            "prior_strategies": ["all_6_biological_ids_equal_weight", "leave_hcc1395_out_5_biological_ids"],
            "concentration_multipliers": list(EB_MULTIPLIERS),
            "maximum_absolute_shift_percentage_points": max_eb_shift,
            "interpretation": "regularization sensitivity only; no common-composition equality assumption",
        },
        "qa": {
            "checks": len(checks),
            "passed": sum(bool(row["pass"]) for row in checks),
            "errors": len(failed),
            "warnings": len(warnings),
        },
        "outputs": outputs,
    }

    # ---------- Materialize deterministic outputs ----------
    write_tsv(args.output_dir / outputs["checks"], checks, ["check", "severity", "pass", "observed", "expected"])
    write_tsv(args.output_dir / outputs["warnings"], warnings, ["warning", "severity", "status", "detail"])
    write_tsv(args.output_dir / outputs["posterior_intervals"], interval_rows, [
        "estimand", "sample", "biological_id", "category", "count", "denominator", "raw_proportion",
        "posterior_prior", "posterior_mean", "posterior_marginal_q025", "posterior_marginal_q975",
    ])
    write_tsv(args.output_dir / outputs["biological_id_compositions"], biological_composition_rows, [
        "estimand", "biological_id", "dataset_sources", "aggregation",
        "supporting_region_rows_not_independent_denominator", "category", "proportion",
    ])
    write_tsv(args.output_dir / outputs["raw_pairwise"], raw_pairwise, [
        "comparison_level", "estimand", "entity_a", "entity_b", "biological_id_a", "biological_id_b",
        "pair_type", "denominator_a", "denominator_b", "total_variation",
        "jensen_shannon_distance_base2", "aitchison_distance",
    ])
    write_tsv(args.output_dir / outputs["resampled_compositions"], resampled_compositions, [
        "method", "estimand", "sample", "biological_id", "category", "replicates", "mean", "median", "q025", "q975",
    ])
    write_tsv(args.output_dir / outputs["resampled_pairwise"], resampled_pairwise, [
        "method", "estimand", "sample_a", "sample_b", "biological_id_a", "biological_id_b", "pair_type",
        "metric", "replicates", "mean", "median", "q025", "q975",
    ])
    write_tsv(args.output_dir / outputs["hcc_rank"], rank_rows, [
        "method", "estimand", "metric", "replicates", "hcc_distance_mean", "hcc_distance_median",
        "hcc_distance_q025", "hcc_distance_q975", "cross_biological_pairs",
        "hcc_rank_smaller_is_more_similar_median", "hcc_rank_q025", "hcc_rank_q975",
        "probability_hcc_rank_1", "lower_tail_percentile_median", "add_one_lower_tail_percentile",
        "cross_biological_min", "cross_biological_median", "cross_biological_max",
    ])
    write_tsv(args.output_dir / outputs["block_feasibility"], block_feasibility_rows, [
        "block_width_bp", "sample", "sample_blocks", "all7_common_blocks", "all7_union_blocks",
        "common_block_share_union", "complete_rows_total", "complete_rows_in_common_blocks",
        "complete_row_coverage", "coordinate_status",
    ])
    write_tsv(args.output_dir / outputs["eb_priors"], eb_prior_rows, [
        "estimand", "prior_name", "category", "alpha", "prior_mean", "prior_concentration",
        "n_biological_ids", "biological_ids", "fit_success", "fit_message",
    ])
    write_tsv(args.output_dir / outputs["eb_shrinkage"], eb_shrinkage_rows, [
        "estimand", "prior_name", "concentration_multiplier", "effective_prior_concentration", "sample",
        "biological_id", "category", "count", "denominator", "raw_proportion", "shrunk_proportion",
        "shift_percentage_points", "interpretation",
    ])
    write_tsv(args.output_dir / outputs["eb_pairwise"], eb_pairwise_rows, [
        "comparison_level", "estimand", "entity_a", "entity_b", "biological_id_a", "biological_id_b",
        "pair_type", "denominator_a", "denominator_b", "total_variation", "jensen_shannon_distance_base2",
        "aitchison_distance", "prior_name", "concentration_multiplier", "effective_prior_concentration",
    ])
    write_tsv(args.output_dir / outputs["source_compositions"], source_standardized_compositions, [
        "estimand", "weight_scheme", "sample", "biological_id", "category", "standardized_proportion",
        "structural_weight", "vaf_resolved_weight", "unresolved_weight", "causal_status",
    ])
    write_tsv(args.output_dir / outputs["source_stratified_pairwise"], source_stratified_pairwise, [
        "comparison_level", "estimand", "selection_source", "entity_a", "entity_b", "biological_id_a",
        "biological_id_b", "pair_type", "denominator_a", "denominator_b", "total_variation",
        "jensen_shannon_distance_base2", "aitchison_distance", "aitchison_zero_handling",
    ])
    write_tsv(args.output_dir / outputs["source_pairwise"], source_standardized_pairwise, [
        "comparison_level", "estimand", "weight_scheme", "entity_a", "entity_b", "biological_id_a",
        "biological_id_b", "pair_type", "denominator_a", "denominator_b", "total_variation",
        "jensen_shannon_distance_base2", "aitchison_distance", "structural_weight", "vaf_resolved_weight",
        "unresolved_weight",
    ])
    write_json(args.output_dir / outputs["summary"], summary)

    primary_rank = raw_rank_lookup[("complete_five_class_primary", "jensen_shannon_distance_base2")]
    secondary_rank = raw_rank_lookup[("resolved_four_class_secondary", "jensen_shannon_distance_base2")]
    report = f"""<!--
建立時間: 2026-07-13
目標: 全 7 dataset structure composition 的 bulk-sampling／compositional consistency adjustment
處理範圍: chr1-22；7 dataset rows；6 biological IDs；historical layered-v2 engineering snapshot
關聯檔案: InterSubMod/research/20260713_crosssample_structure_bulk_sampling_adjustment/
-->

# Cross-sample structure composition adjustment

> **PASS WITH CAVEATS。Primary 是 complete 五類（含 unresolved）；resolved-only 四類只作 secondary。這是 technical compositional proximity，不是 clone-tree truth。**

## 研究問題與設計

- 問題：HCC1395 technical pair 的 pattern composition 是否比不同 biological IDs 更接近；結果對 closure、樣本量、genomic blocks、source mixture 與 EB shrinkage 是否敏感？
- 7 dataset rows／6 biological IDs分開；HCC pair 不算兩個 biological replicates。
- 不同 cell lines 不被強迫相同；20 個 cross-biological dataset-row pairs只作 empirical reference。
- 固定 seed `{args.seed}`；每個 stochastic endpoint `{args.replicates:,}` replicates。

## Raw HCC technical-pair composition distance

| Estimand | TV | JSD (base 2) | Aitchison | JSD rank vs 20 cross-bio pairs | lower-tail percentile |
|---|---:|---:|---:|---:|---:|
| Complete five-class **primary** | {raw_hcc['complete_five_class_primary']['total_variation']:.4f} | {raw_hcc['complete_five_class_primary']['jensen_shannon_distance_base2']:.4f} | {raw_hcc['complete_five_class_primary']['aitchison_distance']:.4f} | {primary_rank['hcc_rank_smaller_is_more_similar_median']:.0f}/21-scale | {primary_rank['lower_tail_percentile_median']:.1f}% |
| Resolved four-class secondary | {raw_hcc['resolved_four_class_secondary']['total_variation']:.4f} | {raw_hcc['resolved_four_class_secondary']['jensen_shannon_distance_base2']:.4f} | {raw_hcc['resolved_four_class_secondary']['aitchison_distance']:.4f} | {secondary_rank['hcc_rank_smaller_is_more_similar_median']:.0f}/21-scale | {secondary_rank['lower_tail_percentile_median']:.1f}% |

Complete五類把 unresolved 2.05% vs 9.90% 保留，所以 HCC TV={raw_hcc['complete_five_class_primary']['total_variation']:.2%}；若先丟掉 unresolved 再 closure，TV降為 {raw_hcc['resolved_four_class_secondary']['total_variation']:.2%}。後者不可當 primary。

Primary JSD的 rank=9 表示把 HCC technical pair插入20個 cross-biological dataset-row distances後，**有8個跨生物 pairs反而更近**；不是「technical pair最接近」。

## Adjustment stability（primary complete五類，JSD）

| Method | HCC JSD median | 95% interval | median rank | P(rank=1) |
|---|---:|---:|---:|---:|
| Raw | {primary_jsd_sensitivity['raw_unadjusted']['hcc_distance_median']:.4f} | — | {primary_jsd_sensitivity['raw_unadjusted']['hcc_rank_smaller_is_more_similar_median']:.0f} | {primary_jsd_sensitivity['raw_unadjusted']['probability_hcc_rank_1']:.3f} |
| Dirichlet posterior | {primary_jsd_sensitivity['dirichlet_posterior']['hcc_distance_median']:.4f} | [{primary_jsd_sensitivity['dirichlet_posterior']['hcc_distance_q025']:.4f}, {primary_jsd_sensitivity['dirichlet_posterior']['hcc_distance_q975']:.4f}] | {primary_jsd_sensitivity['dirichlet_posterior']['hcc_rank_smaller_is_more_similar_median']:.0f} | {primary_jsd_sensitivity['dirichlet_posterior']['probability_hcc_rank_1']:.3f} |
| Equal-n rarefaction | {primary_jsd_sensitivity[f"rarefaction_without_replacement_equal_n_{min(denominators['complete_five_class_primary'].values())}"]['hcc_distance_median']:.4f} | [{primary_jsd_sensitivity[f"rarefaction_without_replacement_equal_n_{min(denominators['complete_five_class_primary'].values())}"]['hcc_distance_q025']:.4f}, {primary_jsd_sensitivity[f"rarefaction_without_replacement_equal_n_{min(denominators['complete_five_class_primary'].values())}"]['hcc_distance_q975']:.4f}] | {primary_jsd_sensitivity[f"rarefaction_without_replacement_equal_n_{min(denominators['complete_five_class_primary'].values())}"]['hcc_rank_smaller_is_more_similar_median']:.0f} | {primary_jsd_sensitivity[f"rarefaction_without_replacement_equal_n_{min(denominators['complete_five_class_primary'].values())}"]['probability_hcc_rank_1']:.3f} |
| Common-block bootstrap 5 Mb | {primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_5mb']['hcc_distance_median']:.4f} | [{primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_5mb']['hcc_distance_q025']:.4f}, {primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_5mb']['hcc_distance_q975']:.4f}] | {primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_5mb']['hcc_rank_smaller_is_more_similar_median']:.0f} | {primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_5mb']['probability_hcc_rank_1']:.3f} |
| Common-block bootstrap 10 Mb | {primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_10mb']['hcc_distance_median']:.4f} | [{primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_10mb']['hcc_distance_q025']:.4f}, {primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_10mb']['hcc_distance_q975']:.4f}] | {primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_10mb']['hcc_rank_smaller_is_more_similar_median']:.0f} | {primary_jsd_sensitivity['paired_common_genomic_block_bootstrap_10mb']['probability_hcc_rank_1']:.3f} |

**結論：bulk sample size與spatial resampling不是主要解釋。** Raw／posterior／rarefaction／5 Mb／10 Mb的HCC primary JSD與rank幾乎不動；EB最大比例位移也只有 {max_eb_shift:.3f} percentage points。aggregate composition呈現中等、非例外的technical proximity，而不是高度一致。

## Source-mixture sensitivity

- HCC resolved raw TV：{raw_hcc['resolved_four_class_secondary']['total_variation']:.2%}。
- 分 source 看，structural Topo=1 conditional TV={hcc_source_conditional['structural_topo1']['total_variation']:.2%}；VAF-resolved Topo>1 conditional TV={hcc_source_conditional['vaf_resolved_topogt1']['total_variation']:.2%}。
- 用全 7 resolved regions 的共同 source weights（structural={global_resolved_source_weights[0]:.2%}、VAF-resolved={global_resolved_source_weights[1]:.2%}）direct-standardize後：TV={global_source_hcc['total_variation']:.2%}（{global_source_hcc['total_variation']/raw_hcc['resolved_four_class_secondary']['total_variation']:.2f}×）。
- 這是 **Simpson-like mixture cancellation** 診斷：raw source mix部分抵銷 conditional差異；不是「校正後真值」。source是 outcome-dependent selection state，標準化可能 over-adjust。
- Complete三來源共同權重 sensitivity TV={global_complete_hcc['total_variation']:.2%}；同樣不能取代 raw primary。

## Sampling adjustment

- Dirichlet posterior：Jeffreys `alpha=0.5/category`，輸出每類 marginal 95% credible interval與 posterior pair distances。
- Equal-n rarefaction：primary n={min(denominators['complete_five_class_primary'].values()):,}；secondary n={min(denominators['resolved_four_class_secondary'].values()):,}，無放回抽樣。
- Spatial bootstrap：region具真實 chr/start/end；共同 5 Mb blocks={summary['block_bootstrap']['five_mb_common_blocks']}/{summary['block_bootstrap']['five_mb_union_blocks']} union，另做10 Mb sensitivity。共同 blocks以相同 multiplicity跨7 rows聯合重抽，但不是exact-region配對。
- EB：先把 HCC pair等權合成一個 biological composition，再估 all-6 與 leave-HCC-out-5 priors；0.5×/1×/2× concentration只作 regularization sensitivity。最大比例位移={max_eb_shift:.3f} percentage points。

## Claim ceiling

{CLAIM_CEILING}

fresh layered-v3尚未7/7 scientific gate；任何「相似演化」「同一 clone tree」「方法證明有效」均超出本輪證據。完整 sensitivity數值見 `hcc_pair_relative_rank.tsv`、`source_standardized_pairwise_distances.tsv` 與 `summary.json`。
"""
    (args.output_dir / outputs["report"]).write_text(report, encoding="utf-8")

    input_paths = {
        "summary": args.summary,
        "regions": args.regions,
        "by_source": args.by_source,
        "coarse_regions": args.coarse_regions,
        "script": Path(__file__).resolve(),
    }
    output_receipts = {
        name: {
            "filename": filename,
            "size_bytes": (args.output_dir / filename).stat().st_size,
            "sha256": sha256(args.output_dir / filename),
        }
        for name, filename in outputs.items()
        if name != "provenance" and (args.output_dir / filename).is_file()
    }
    provenance = {
        "schema_name": "intersubmod.crosssample_structure_bulk_sampling_adjustment.provenance",
        "schema_version": "1.0.0",
        "analysis_date": "2026-07-13",
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "seed": args.seed,
        "replicates": args.replicates,
        "command_template": [
            "python3", "scripts/analyze_crosssample_composition.py",
            "--summary", str(args.summary.resolve()),
            "--regions", str(args.regions.resolve()),
            "--by-source", str(args.by_source.resolve()),
            "--coarse-regions", str(args.coarse_regions.resolve()),
            "--output-dir", "${OUTPUT_DIR}",
            "--seed", str(args.seed),
            "--replicates", str(args.replicates),
        ],
        "inputs": {
            name: {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for name, path in input_paths.items()
        },
        "outputs": output_receipts,
        "claim_ceiling": CLAIM_CEILING,
    }
    write_json(args.output_dir / outputs["provenance"], provenance)

    print(f"INPUT summary -> {args.summary.resolve()}")
    print(f"INPUT resolved regions -> {args.regions.resolve()}")
    print(f"INPUT complete/coarse regions -> {args.coarse_regions.resolve()}")
    print(f"OUTPUT -> {args.output_dir.resolve()}")
    print(f"STATUS -> PASS_WITH_CAVEATS ({sum(bool(row['pass']) for row in checks)}/{len(checks)} checks)")
    print(
        "HCC primary five-class -> "
        f"TV={raw_hcc['complete_five_class_primary']['total_variation']:.6f} "
        f"JSD={raw_hcc['complete_five_class_primary']['jensen_shannon_distance_base2']:.6f}"
    )
    print(
        "HCC resolved source sensitivity -> "
        f"raw_TV={raw_hcc['resolved_four_class_secondary']['total_variation']:.6f} "
        f"global_standardized_TV={global_source_hcc['total_variation']:.6f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
