#!/usr/bin/env python3
"""Summarize PyClone-VI conditional clustering and HCC1395 reproducibility."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Mapping, MutableMapping, NamedTuple, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr, spearmanr


CLONAL_THRESHOLD = 0.90


class ResultRow(NamedTuple):
    mutation_id: str
    sample_id: str
    cluster_id: int
    cellular_prevalence: float
    cellular_prevalence_std: float
    assignment_prob: float


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_results(path: Path) -> List[ResultRow]:
    rows: List[ResultRow] = []
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(ResultRow(
                row["mutation_id"], row["sample_id"], int(row["cluster_id"]),
                float(row["cellular_prevalence"]), float(row["cellular_prevalence_std"]),
                float(row["cluster_assignment_prob"]),
            ))
    return rows


def comb2(value: int) -> float:
    return value * (value - 1) / 2


def label_metrics(labels_a: Sequence[int], labels_b: Sequence[int]) -> Mapping[str, float]:
    if len(labels_a) != len(labels_b) or not labels_a:
        return {"ari": math.nan, "nmi": math.nan, "hungarian_agreement": math.nan}
    a_values = sorted(set(labels_a))
    b_values = sorted(set(labels_b))
    a_index = {value: index for index, value in enumerate(a_values)}
    b_index = {value: index for index, value in enumerate(b_values)}
    table = np.zeros((len(a_values), len(b_values)), dtype=np.int64)
    for first, second in zip(labels_a, labels_b):
        table[a_index[first], b_index[second]] += 1
    row_sums = table.sum(axis=1)
    column_sums = table.sum(axis=0)
    n = int(table.sum())
    sum_cells = float(sum(comb2(int(value)) for value in table.ravel()))
    sum_rows = float(sum(comb2(int(value)) for value in row_sums))
    sum_columns = float(sum(comb2(int(value)) for value in column_sums))
    total_pairs = comb2(n)
    expected = sum_rows * sum_columns / total_pairs if total_pairs else 0.0
    maximum = 0.5 * (sum_rows + sum_columns)
    ari = (sum_cells - expected) / (maximum - expected) if maximum != expected else 1.0
    probabilities = table / n
    p_a = row_sums / n
    p_b = column_sums / n
    mi = 0.0
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            if probabilities[i, j] > 0:
                mi += probabilities[i, j] * math.log(probabilities[i, j] / (p_a[i] * p_b[j]))
    h_a = -sum(value * math.log(value) for value in p_a if value > 0)
    h_b = -sum(value * math.log(value) for value in p_b if value > 0)
    nmi = 2 * mi / (h_a + h_b) if h_a + h_b else 1.0
    row_match, column_match = linear_sum_assignment(-table)
    matched = int(table[row_match, column_match].sum())
    mapping_b_to_a = {b_values[j]: a_values[i] for i, j in zip(row_match, column_match)}
    return {
        "ari": float(ari),
        "nmi": float(nmi),
        "hungarian_agreement": matched / n,
        "clusters_a": len(a_values),
        "clusters_b": len(b_values),
        "hungarian_mapping_b_to_a": mapping_b_to_a,
    }


def safe_correlation(first: Sequence[float], second: Sequence[float], method: str) -> float:
    if len(first) < 2 or len(set(first)) < 2 or len(set(second)) < 2:
        return math.nan
    result = pearsonr(first, second) if method == "pearson" else spearmanr(first, second)
    return float(result.statistic)


def js_divergence(first: Sequence[float], second: Sequence[float]) -> float:
    p = np.asarray(first, dtype=float)
    q = np.asarray(second, dtype=float)
    if p.sum() <= 0 or q.sum() <= 0:
        return math.nan
    p = p / p.sum()
    q = q / q.sum()
    midpoint = 0.5 * (p + q)
    with np.errstate(divide="ignore", invalid="ignore"):
        first_terms = np.where(p > 0, p * np.log2(p / midpoint), 0.0)
        second_terms = np.where(q > 0, q * np.log2(q / midpoint), 0.0)
    return float(0.5 * (first_terms.sum() + second_terms.sum()))


def by_sample(rows: Sequence[ResultRow]) -> Mapping[str, Mapping[str, ResultRow]]:
    result: Dict[str, Dict[str, ResultRow]] = defaultdict(dict)
    for row in rows:
        if row.mutation_id in result[row.sample_id]:
            raise ValueError(f"duplicate result row: {row.sample_id} {row.mutation_id}")
        result[row.sample_id][row.mutation_id] = row
    return result


def summarize_fit(bundle_id: str, rows: Sequence[ResultRow]) -> Tuple[Mapping[str, object], List[Mapping[str, object]]]:
    sample_map = by_sample(rows)
    all_probabilities = [row.assignment_prob for row in rows]
    cluster_profiles: List[Mapping[str, object]] = []
    per_sample: Dict[str, Mapping[str, object]] = {}
    for sample, mutation_rows in sample_map.items():
        cluster_counts = Counter(row.cluster_id for row in mutation_rows.values())
        cluster_cp = {row.cluster_id: row.cellular_prevalence for row in mutation_rows.values()}
        cluster_std = {row.cluster_id: row.cellular_prevalence_std for row in mutation_rows.values()}
        clonal_count = sum(
            count for cluster, count in cluster_counts.items() if cluster_cp[cluster] >= CLONAL_THRESHOLD
        )
        mutation_count = len(mutation_rows)
        clonal_probabilities = [
            row.assignment_prob
            for row in mutation_rows.values()
            if row.cellular_prevalence >= CLONAL_THRESHOLD
        ]
        subclonal_probabilities = [
            row.assignment_prob
            for row in mutation_rows.values()
            if row.cellular_prevalence < CLONAL_THRESHOLD
        ]
        per_sample[sample] = {
            "mutations": mutation_count,
            "clusters": len(cluster_counts),
            "clonal_mutations": clonal_count,
            "subclonal_mutations": mutation_count - clonal_count,
            "clonal_fraction": clonal_count / mutation_count,
            "subclonal_fraction": (mutation_count - clonal_count) / mutation_count,
            "clonal_assignment_probability_mean": (
                float(np.mean(clonal_probabilities)) if clonal_probabilities else math.nan
            ),
            "subclonal_assignment_probability_mean": (
                float(np.mean(subclonal_probabilities)) if subclonal_probabilities else math.nan
            ),
            "subclonal_assignment_probability_median": (
                float(np.median(subclonal_probabilities)) if subclonal_probabilities else math.nan
            ),
            "subclonal_assignment_probability_lt_0_8_fraction": (
                sum(value < 0.8 for value in subclonal_probabilities) / len(subclonal_probabilities)
                if subclonal_probabilities else math.nan
            ),
        }
        for cluster in sorted(cluster_counts):
            cluster_probabilities = [
                row.assignment_prob for row in mutation_rows.values() if row.cluster_id == cluster
            ]
            cluster_profiles.append({
                "bundle_id": bundle_id,
                "sample_id": sample,
                "cluster_id": cluster,
                "mutation_count": cluster_counts[cluster],
                "mutation_fraction": cluster_counts[cluster] / mutation_count,
                "cellular_prevalence": cluster_cp[cluster],
                "cellular_prevalence_std": cluster_std[cluster],
                "assignment_probability_mean": float(np.mean(cluster_probabilities)),
                "assignment_probability_median": float(np.median(cluster_probabilities)),
                "assignment_probability_lt_0_8_fraction": (
                    sum(value < 0.8 for value in cluster_probabilities) / len(cluster_probabilities)
                ),
                "clone_class": "clonal" if cluster_cp[cluster] >= CLONAL_THRESHOLD else "subclonal",
            })
    mutation_ids = {row.mutation_id for row in rows}
    summary = {
        "bundle_id": bundle_id,
        "samples": sorted(sample_map),
        "mutations": len(mutation_ids),
        "rows": len(rows),
        "clusters": len({row.cluster_id for row in rows}),
        "assignment_probability_mean": float(np.mean(all_probabilities)),
        "assignment_probability_median": float(np.median(all_probabilities)),
        "assignment_probability_lt_0_8_fraction": sum(value < 0.8 for value in all_probabilities) / len(all_probabilities),
        "clonal_threshold": CLONAL_THRESHOLD,
        "per_sample": per_sample,
    }
    return summary, cluster_profiles


def compare_joint(rows: Sequence[ResultRow]) -> Mapping[str, object]:
    sample_map = by_sample(rows)
    samples = sorted(sample_map)
    if len(samples) != 2:
        raise ValueError(f"Joint comparison expects two samples, got {samples}")
    first, second = samples
    mutations = sorted(set(sample_map[first]) & set(sample_map[second]))
    if len(mutations) != len(sample_map[first]) or len(mutations) != len(sample_map[second]):
        raise ValueError("Joint result is not a complete shared mutation matrix")
    first_rows = [sample_map[first][mutation] for mutation in mutations]
    second_rows = [sample_map[second][mutation] for mutation in mutations]
    cluster_counts = Counter(row.cluster_id for row in first_rows)
    first_cp = {row.cluster_id: row.cellular_prevalence for row in first_rows}
    second_cp = {row.cluster_id: row.cellular_prevalence for row in second_rows}
    weighted_abs_delta = sum(
        cluster_counts[cluster] * abs(first_cp[cluster] - second_cp[cluster]) for cluster in cluster_counts
    ) / len(mutations)
    clusters = sorted(cluster_counts)
    abundance_first = [cluster_counts[cluster] * first_cp[cluster] for cluster in clusters]
    abundance_second = [cluster_counts[cluster] * second_cp[cluster] for cluster in clusters]
    clonal_first = [row.cellular_prevalence >= CLONAL_THRESHOLD for row in first_rows]
    clonal_second = [row.cellular_prevalence >= CLONAL_THRESHOLD for row in second_rows]
    both_clonal = sum(a and b for a, b in zip(clonal_first, clonal_second))
    either_clonal = sum(a or b for a, b in zip(clonal_first, clonal_second))
    return {
        "samples": samples,
        "mutations": len(mutations),
        "clusters": len(clusters),
        "weighted_mean_absolute_cluster_prevalence_delta": weighted_abs_delta,
        "mutation_level_prevalence_pearson": safe_correlation(
            [row.cellular_prevalence for row in first_rows],
            [row.cellular_prevalence for row in second_rows], "pearson",
        ),
        "mutation_level_prevalence_spearman": safe_correlation(
            [row.cellular_prevalence for row in first_rows],
            [row.cellular_prevalence for row in second_rows], "spearman",
        ),
        "mutation_weighted_cluster_cp_profile_js_divergence_bits": js_divergence(
            abundance_first, abundance_second
        ),
        "clonal_subclonal_state_agreement": sum(a == b for a, b in zip(clonal_first, clonal_second)) / len(mutations),
        "clonal_mutation_jaccard": both_clonal / either_clonal if either_clonal else 1.0,
        "assignment_note": (
            "Joint-fit cluster labels are shared by model construction; label agreement is not an independent reproducibility metric. "
            "The mutation-weighted cluster CP profile is a descriptive profile, not a direct estimate of clone abundance."
        ),
    }


def compare_separate(rows_a: Sequence[ResultRow], rows_b: Sequence[ResultRow]) -> Mapping[str, object]:
    map_a = by_sample(rows_a)
    map_b = by_sample(rows_b)
    if len(map_a) != 1 or len(map_b) != 1:
        raise ValueError("Separate-fit comparison expects one sample per fit")
    sample_a = next(iter(map_a))
    sample_b = next(iter(map_b))
    a = map_a[sample_a]
    b = map_b[sample_b]
    mutations = sorted(set(a) & set(b))
    labels_a = [a[mutation].cluster_id for mutation in mutations]
    labels_b = [b[mutation].cluster_id for mutation in mutations]
    metrics = dict(label_metrics(labels_a, labels_b))
    mapping = metrics.pop("hungarian_mapping_b_to_a")
    counts_a = Counter(labels_a)
    counts_b_mapped: MutableMapping[object, int] = Counter()
    for label, count in Counter(labels_b).items():
        counts_b_mapped[mapping.get(label, f"unmatched_b_{label}")] += count
    axes = sorted(set(counts_a) | set(counts_b_mapped), key=str)
    prevalence_a = [a[mutation].cellular_prevalence for mutation in mutations]
    prevalence_b = [b[mutation].cellular_prevalence for mutation in mutations]
    clonal_a = [value >= CLONAL_THRESHOLD for value in prevalence_a]
    clonal_b = [value >= CLONAL_THRESHOLD for value in prevalence_b]
    both_clonal = sum(first and second for first, second in zip(clonal_a, clonal_b))
    either_clonal = sum(first or second for first, second in zip(clonal_a, clonal_b))
    both_subclonal = sum(not first and not second for first, second in zip(clonal_a, clonal_b))
    either_subclonal = sum(not first or not second for first, second in zip(clonal_a, clonal_b))
    a_only_subclonal = sum(not first and second for first, second in zip(clonal_a, clonal_b))
    b_only_subclonal = sum(first and not second for first, second in zip(clonal_a, clonal_b))
    observed_state_agreement = (both_clonal + both_subclonal) / len(mutations)
    a_subclonal_fraction = (both_subclonal + a_only_subclonal) / len(mutations)
    b_subclonal_fraction = (both_subclonal + b_only_subclonal) / len(mutations)
    expected_state_agreement = (
        a_subclonal_fraction * b_subclonal_fraction
        + (1 - a_subclonal_fraction) * (1 - b_subclonal_fraction)
    )
    state_kappa = (
        (observed_state_agreement - expected_state_agreement) / (1 - expected_state_agreement)
        if expected_state_agreement < 1 else 1.0
    )
    mcc_denominator = math.sqrt(
        (both_subclonal + a_only_subclonal)
        * (both_subclonal + b_only_subclonal)
        * (both_clonal + a_only_subclonal)
        * (both_clonal + b_only_subclonal)
    )
    state_mcc = (
        (both_subclonal * both_clonal - a_only_subclonal * b_only_subclonal) / mcc_denominator
        if mcc_denominator else math.nan
    )
    both_subclonal_ids = [
        mutation for mutation, first, second in zip(mutations, clonal_a, clonal_b)
        if not first and not second
    ]
    both_subclonal_cluster_metrics = (
        label_metrics(
            [a[mutation].cluster_id for mutation in both_subclonal_ids],
            [b[mutation].cluster_id for mutation in both_subclonal_ids],
        )
        if both_subclonal_ids else {}
    )
    metrics.update({
        "samples": [sample_a, sample_b],
        "mutations_intersection": len(mutations),
        "mutations_a": len(a),
        "mutations_b": len(b),
        "mutation_universe_identical": set(a) == set(b),
        "hungarian_mapping_b_to_a": {str(key): value for key, value in mapping.items()},
        "cluster_mutation_count_profile_js_divergence_bits": js_divergence(
            [counts_a.get(axis, 0) for axis in axes],
            [counts_b_mapped.get(axis, 0) for axis in axes],
        ),
        "mutation_level_prevalence_mae": float(np.mean(np.abs(np.asarray(prevalence_a) - np.asarray(prevalence_b)))),
        "mutation_level_prevalence_pearson": safe_correlation(prevalence_a, prevalence_b, "pearson"),
        "mutation_level_prevalence_spearman": safe_correlation(prevalence_a, prevalence_b, "spearman"),
        "clonal_subclonal_state_agreement": observed_state_agreement,
        "clonal_subclonal_state_kappa": state_kappa,
        "clonal_subclonal_state_mcc": state_mcc,
        "clonal_mutation_jaccard": both_clonal / either_clonal if either_clonal else 1.0,
        "subclonal_mutation_jaccard": both_subclonal / either_subclonal if either_subclonal else 1.0,
        "subclonal_state_f1": (
            2 * both_subclonal / (2 * both_subclonal + a_only_subclonal + b_only_subclonal)
            if both_subclonal or a_only_subclonal or b_only_subclonal else 1.0
        ),
        "subclonal_contingency": {
            "both_subclonal": both_subclonal,
            "a_only_subclonal": a_only_subclonal,
            "b_only_subclonal": b_only_subclonal,
            "both_clonal": both_clonal,
            "either_subclonal": either_subclonal,
        },
        "both_subclonal_cluster_metrics": both_subclonal_cluster_metrics,
        "imbalance_note": (
            "Overall label agreement is dominated by the clonal majority. Subclonal Jaccard, "
            "binary-state kappa/MCC, and the both-subclonal cluster metrics are the primary "
            "minor-clone reproducibility checks."
        ),
    })
    return metrics


def compare_sensitivity(rows_main: Sequence[ResultRow], rows_sensitivity: Sequence[ResultRow]) -> Mapping[str, object]:
    main = by_sample(rows_main)
    sensitivity = by_sample(rows_sensitivity)
    common_samples = sorted(set(main) & set(sensitivity))
    sample_metrics = {}
    for sample in common_samples:
        mutations = sorted(set(main[sample]) & set(sensitivity[sample]))
        labels_main = [main[sample][mutation].cluster_id for mutation in mutations]
        labels_sensitivity = [sensitivity[sample][mutation].cluster_id for mutation in mutations]
        metrics = dict(label_metrics(labels_main, labels_sensitivity))
        metrics.pop("hungarian_mapping_b_to_a")
        cp_main = [main[sample][mutation].cellular_prevalence for mutation in mutations]
        cp_sensitivity = [sensitivity[sample][mutation].cellular_prevalence for mutation in mutations]
        main_subclonal = [value < CLONAL_THRESHOLD for value in cp_main]
        sensitivity_subclonal = [value < CLONAL_THRESHOLD for value in cp_sensitivity]
        both_subclonal = sum(
            first and second for first, second in zip(main_subclonal, sensitivity_subclonal)
        )
        either_subclonal = sum(
            first or second for first, second in zip(main_subclonal, sensitivity_subclonal)
        )
        observed_state_agreement = sum(
            first == second for first, second in zip(main_subclonal, sensitivity_subclonal)
        ) / len(mutations)
        main_subclonal_fraction = sum(main_subclonal) / len(mutations)
        sensitivity_subclonal_fraction = sum(sensitivity_subclonal) / len(mutations)
        expected_state_agreement = (
            main_subclonal_fraction * sensitivity_subclonal_fraction
            + (1 - main_subclonal_fraction) * (1 - sensitivity_subclonal_fraction)
        )
        metrics.update({
            "mutations_intersection": len(mutations),
            "mutations_main": len(main[sample]),
            "mutations_sensitivity": len(sensitivity[sample]),
            "prevalence_mae": float(np.mean(np.abs(np.asarray(cp_main) - np.asarray(cp_sensitivity)))),
            "prevalence_spearman": safe_correlation(cp_main, cp_sensitivity, "spearman"),
            "main_subclonal_mutations_intersection": sum(main_subclonal),
            "sensitivity_subclonal_mutations": sum(sensitivity_subclonal),
            "both_subclonal_mutations": both_subclonal,
            "either_subclonal_mutations": either_subclonal,
            "subclonal_mutation_jaccard": (
                both_subclonal / either_subclonal if either_subclonal else 1.0
            ),
            "clonal_subclonal_state_agreement": observed_state_agreement,
            "clonal_subclonal_state_kappa": (
                (observed_state_agreement - expected_state_agreement) / (1 - expected_state_agreement)
                if expected_state_agreement < 1 else 1.0
            ),
        })
        sample_metrics[sample] = metrics
    return {"samples": sample_metrics}


def compare_joint_to_separate(
    rows_joint: Sequence[ResultRow], rows_separate: Sequence[ResultRow], sample: str
) -> Mapping[str, object]:
    """Measure modeling-mode sensitivity for one sample on the shared mutation universe."""
    joint_map = by_sample(rows_joint)
    separate_map = by_sample(rows_separate)
    if sample not in joint_map or list(separate_map) != [sample]:
        raise ValueError(f"Joint/separate comparison sample mismatch for {sample}")
    joint = joint_map[sample]
    separate = separate_map[sample]
    mutations = sorted(set(joint) & set(separate))
    metrics = dict(label_metrics(
        [joint[mutation].cluster_id for mutation in mutations],
        [separate[mutation].cluster_id for mutation in mutations],
    ))
    metrics.pop("hungarian_mapping_b_to_a")
    joint_cp = [joint[mutation].cellular_prevalence for mutation in mutations]
    separate_cp = [separate[mutation].cellular_prevalence for mutation in mutations]
    joint_subclonal = [value < CLONAL_THRESHOLD for value in joint_cp]
    separate_subclonal = [value < CLONAL_THRESHOLD for value in separate_cp]
    both_subclonal = sum(first and second for first, second in zip(joint_subclonal, separate_subclonal))
    either_subclonal = sum(first or second for first, second in zip(joint_subclonal, separate_subclonal))
    metrics.update({
        "sample": sample,
        "mutations_intersection": len(mutations),
        "joint_subclonal_mutations": sum(joint_subclonal),
        "separate_subclonal_mutations": sum(separate_subclonal),
        "both_subclonal_mutations": both_subclonal,
        "either_subclonal_mutations": either_subclonal,
        "subclonal_mutation_jaccard": both_subclonal / either_subclonal if either_subclonal else 1.0,
        "cellular_prevalence_mae": float(
            np.mean(np.abs(np.asarray(joint_cp) - np.asarray(separate_cp)))
        ),
        "cellular_prevalence_spearman": safe_correlation(joint_cp, separate_cp, "spearman"),
        "interpretation_note": (
            "Joint fitting shares cluster structure across technical sources, whereas separate fits do not. "
            "Their subclonal-set overlap is a modeling-mode sensitivity result, not biological truth."
        ),
    })
    return metrics


def parse_args() -> argparse.Namespace:
    topic = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=topic / "runs" / "pyclone_vi")
    parser.add_argument("--output-dir", type=Path, default=topic / "analysis" / "pyclone_vi")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    status_paths = sorted(args.run_dir.glob("*/status.json"))
    pass_jobs = {}
    rows_by_bundle = {}
    for status_path in status_paths:
        status = json.loads(status_path.read_text())
        if status.get("status") != "PASS":
            continue
        bundle_id = status["bundle_id"]
        results_path = Path(status["outputs"]["results"])
        pass_jobs[bundle_id] = status
        rows_by_bundle[bundle_id] = read_results(results_path)
    if not rows_by_bundle:
        raise ValueError(f"No PASS PyClone-VI jobs under {args.run_dir}")

    fit_summaries = {}
    cluster_profiles: List[Mapping[str, object]] = []
    for bundle_id, rows in rows_by_bundle.items():
        summary, profiles = summarize_fit(bundle_id, rows)
        fit_summaries[bundle_id] = summary
        cluster_profiles.extend(profiles)

    joint_comparisons = {}
    for bundle_id, rows in rows_by_bundle.items():
        if "hcc1395_pair" in bundle_id and "joint" in bundle_id:
            joint_comparisons[bundle_id] = compare_joint(rows)

    separate_comparisons = {}
    for suffix in ("main", "near_integer"):
        first_id = f"hcc1395_pair_primary_separate_HCC1395_{suffix}"
        second_id = f"hcc1395_pair_primary_separate_HCC1395_DORADO_{suffix}"
        if first_id in rows_by_bundle and second_id in rows_by_bundle:
            separate_comparisons[suffix] = compare_separate(rows_by_bundle[first_id], rows_by_bundle[second_id])

    model_mode_comparisons = {}
    joint_id = "hcc1395_pair_primary_joint_main"
    if joint_id in rows_by_bundle:
        for sample in ("HCC1395", "HCC1395_DORADO"):
            separate_id = f"hcc1395_pair_primary_separate_{sample}_main"
            if separate_id in rows_by_bundle:
                model_mode_comparisons[sample] = compare_joint_to_separate(
                    rows_by_bundle[joint_id], rows_by_bundle[separate_id], sample
                )

    sensitivity_pairs = [
        ("hcc1395_pair_primary_joint_main", "hcc1395_pair_primary_joint_near_integer", "hcc1395_primary_joint_main_vs_near_integer"),
        ("hcc1395_pair_primary_separate_HCC1395_main", "hcc1395_pair_primary_separate_HCC1395_near_integer", "HCC1395_separate_main_vs_near_integer"),
        ("hcc1395_pair_primary_separate_HCC1395_DORADO_main", "hcc1395_pair_primary_separate_HCC1395_DORADO_near_integer", "HCC1395_DORADO_separate_main_vs_near_integer"),
        ("hcc1395_pair_pass_union_joint_main", "hcc1395_pair_pass_union_joint_near_integer", "hcc1395_pass_union_joint_main_vs_near_integer"),
        ("H1437_individual_main", "H1437_individual_near_integer", "H1437_main_vs_near_integer"),
        ("H2009_individual_main", "H2009_individual_near_integer", "H2009_main_vs_near_integer"),
        ("HCC1954_individual_main", "HCC1954_individual_near_integer", "HCC1954_main_vs_near_integer"),
        ("hcc1395_pair_primary_joint_main", "hcc1395_pair_pass_union_joint_main", "hcc1395_primary_vs_pass_union_main"),
    ]
    sensitivity = {}
    for main_id, sensitivity_id, label in sensitivity_pairs:
        if main_id in rows_by_bundle and sensitivity_id in rows_by_bundle:
            sensitivity[label] = compare_sensitivity(rows_by_bundle[main_id], rows_by_bundle[sensitivity_id])

    output = {
        "schema_name": "intersubmod.pyclone_vi_conditional_clustering_analysis",
        "schema_version": "1.0.0",
        "claim_ceiling": (
            "Conditional external mutation clustering under caller DP/AF, supplied SAVANA CN/purity, and selected truth-like mutation universes. "
            "PyClone-VI does not infer a phylogenetic tree and is not independent ground truth for InterSubMod."
        ),
        "clonal_threshold": CLONAL_THRESHOLD,
        "fit_summaries": fit_summaries,
        "hcc1395_joint_comparisons": joint_comparisons,
        "hcc1395_separate_fit_comparisons": separate_comparisons,
        "hcc1395_model_mode_comparisons": model_mode_comparisons,
        "sensitivity_comparisons": sensitivity,
        "pass_bundle_count": len(rows_by_bundle),
        "pending_or_failed_bundles": sorted(
            path.parent.name for path in status_paths if json.loads(path.read_text()).get("status") != "PASS"
        ),
        "source_receipts": {
            bundle_id: {
                "results_path": status["outputs"]["results"],
                "results_sha256": status["outputs"]["results_sha256"],
                "status_path": str((args.run_dir / bundle_id / "status.json").resolve()),
            }
            for bundle_id, status in pass_jobs.items()
        },
    }
    json_path = args.output_dir / "analysis_summary.json"
    json_path.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=True) + "\n")

    fit_path = args.output_dir / "fit_summaries.tsv"
    with fit_path.open("w") as handle:
        handle.write(
            "bundle_id\tsamples\tmutations\tclusters\tassignment_prob_mean\tassignment_prob_median\t"
            "assignment_prob_lt_0.8_fraction\n"
        )
        for bundle_id, summary in sorted(fit_summaries.items()):
            handle.write(
                f"{bundle_id}\t{','.join(summary['samples'])}\t{summary['mutations']}\t{summary['clusters']}\t"
                f"{summary['assignment_probability_mean']:.8g}\t{summary['assignment_probability_median']:.8g}\t"
                f"{summary['assignment_probability_lt_0_8_fraction']:.8g}\n"
            )
    cluster_path = args.output_dir / "cluster_profiles.tsv"
    with cluster_path.open("w") as handle:
        columns = [
            "bundle_id", "sample_id", "cluster_id", "mutation_count", "mutation_fraction",
            "cellular_prevalence", "cellular_prevalence_std", "assignment_probability_mean",
            "assignment_probability_median", "assignment_probability_lt_0_8_fraction", "clone_class",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        writer.writerows(cluster_profiles)

    joint_path = args.output_dir / "hcc1395_joint_comparisons.tsv"
    with joint_path.open("w") as handle:
        columns = [
            "bundle_id", "mutations", "clusters", "weighted_mean_abs_cp_delta", "cp_pearson",
            "cp_spearman", "mutation_weighted_cluster_cp_profile_jsd_bits", "clonal_state_agreement",
            "clonal_jaccard",
        ]
        handle.write("\t".join(columns) + "\n")
        for bundle_id, values in sorted(joint_comparisons.items()):
            handle.write("\t".join(map(str, [
                bundle_id, values["mutations"], values["clusters"],
                values["weighted_mean_absolute_cluster_prevalence_delta"],
                values["mutation_level_prevalence_pearson"], values["mutation_level_prevalence_spearman"],
                values["mutation_weighted_cluster_cp_profile_js_divergence_bits"],
                values["clonal_subclonal_state_agreement"],
                values["clonal_mutation_jaccard"],
            ])) + "\n")

    separate_path = args.output_dir / "hcc1395_separate_fit_comparisons.tsv"
    with separate_path.open("w") as handle:
        columns = [
            "scenario", "mutations_intersection", "clusters_a", "clusters_b", "ari", "nmi",
            "hungarian_agreement", "cluster_mutation_count_profile_jsd_bits", "cp_mae", "cp_pearson",
            "cp_spearman",
            "clonal_state_agreement", "clonal_jaccard",
        ]
        handle.write("\t".join(columns) + "\n")
        for scenario, values in sorted(separate_comparisons.items()):
            handle.write("\t".join(map(str, [
                scenario, values["mutations_intersection"], values["clusters_a"], values["clusters_b"],
                values["ari"], values["nmi"], values["hungarian_agreement"],
                values["cluster_mutation_count_profile_js_divergence_bits"], values["mutation_level_prevalence_mae"],
                values["mutation_level_prevalence_pearson"], values["mutation_level_prevalence_spearman"],
                values["clonal_subclonal_state_agreement"], values["clonal_mutation_jaccard"],
            ])) + "\n")

    sensitivity_path = args.output_dir / "sensitivity_comparisons.tsv"
    with sensitivity_path.open("w") as handle:
        columns = [
            "comparison", "sample_id", "mutations_intersection", "mutations_main", "mutations_sensitivity",
            "clusters_main", "clusters_sensitivity", "ari", "nmi", "hungarian_agreement",
            "prevalence_mae", "prevalence_spearman",
        ]
        handle.write("\t".join(columns) + "\n")
        for comparison, comparison_values in sorted(sensitivity.items()):
            for sample, values in sorted(comparison_values["samples"].items()):
                handle.write("\t".join(map(str, [
                    comparison, sample, values["mutations_intersection"], values["mutations_main"],
                    values["mutations_sensitivity"], values["clusters_a"], values["clusters_b"],
                    values["ari"], values["nmi"], values["hungarian_agreement"],
                    values["prevalence_mae"], values["prevalence_spearman"],
                ])) + "\n")

    topic = Path(__file__).resolve().parents[1]
    input_qa_paths = sorted((topic / "data" / "pyclone_inputs").glob("*.qa.json"))
    status_by_bundle = {path.parent.name: json.loads(path.read_text()) for path in status_paths}
    readiness_path = args.output_dir / "all_ready_summary.tsv"
    readiness_rows = []
    with readiness_path.open("w") as handle:
        columns = [
            "bundle_id", "status", "mutations", "samples", "clusters", "num_restarts", "seed", "results_path",
        ]
        handle.write("\t".join(columns) + "\n")
        for qa_path in input_qa_paths:
            qa = json.loads(qa_path.read_text())
            bundle_id = qa["bundle_id"]
            status = status_by_bundle.get(bundle_id)
            if status is None:
                values = [bundle_id, "PENDING", qa["counters"]["included_mutations"], ",".join(qa["samples"]), "", "", "", ""]
            else:
                result_shape = status.get("result_shape", {})
                values = [
                    bundle_id, status.get("status", "UNKNOWN"), status["input"]["shape"]["mutations"],
                    ",".join(status["input"]["shape"]["samples"]), result_shape.get("clusters", ""),
                    status.get("parameters", {}).get("num_restarts", ""), status.get("parameters", {}).get("seed", ""),
                    status.get("outputs", {}).get("results", ""),
                ]
            readiness_rows.append(values)
            handle.write("\t".join(map(str, values)) + "\n")
    output["all_ready"] = bool(readiness_rows) and all(row[1] == "PASS" for row in readiness_rows)
    output["expected_bundle_count"] = len(readiness_rows)
    output["readiness_tsv"] = str(readiness_path.resolve())
    json_path.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=True) + "\n")

    config = json.loads((topic / "config" / "pyclone_validation_config.json").read_text())
    blocked_path = args.output_dir / "blocked_samples.tsv"
    with blocked_path.open("w") as handle:
        handle.write("sample\tstatus\treason\n")
        for sample, reason in sorted(config["blocked_samples"].items()):
            handle.write(f"{sample}\tBLOCKED\t{reason}\n")
    provenance = {
        "analyzer": str(Path(__file__).resolve()),
        "analyzer_sha256": sha256_file(Path(__file__).resolve()),
        "run_dir": str(args.run_dir.resolve()),
        "analysis_summary": str(json_path.resolve()),
        "fit_summaries": str(fit_path.resolve()),
        "cluster_profiles": str(cluster_path.resolve()),
        "hcc1395_joint_comparisons": str(joint_path.resolve()),
        "hcc1395_separate_fit_comparisons": str(separate_path.resolve()),
        "sensitivity_comparisons": str(sensitivity_path.resolve()),
        "all_ready_summary": str(readiness_path.resolve()),
        "blocked_samples": str(blocked_path.resolve()),
    }
    (args.output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"Analyzed {len(rows_by_bundle)} PASS PyClone-VI bundles")
    print(f"Analysis summary: {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
