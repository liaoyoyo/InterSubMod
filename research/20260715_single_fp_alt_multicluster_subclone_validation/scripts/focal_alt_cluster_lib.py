#!/usr/bin/env python3
"""Shared focal-ALT clustering and confound diagnostics for the validation round."""

from __future__ import annotations

import csv
import hashlib
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score, silhouette_score


MIN_SIZE = 3
SEP_MIN = 1.3
RNULL = 40
MODAL_CONFIDENCE = 0.7
HIDDEN_HETEROGENEITY_FRACTION = 0.30


def parse_float(value: str) -> float:
    value = value.strip()
    return np.nan if value in {"", "NA", "nan", "NaN"} else float(value)


def load_matrix(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return [], np.empty((0, 0), dtype=float)
    identifiers: list[str] = []
    values: list[list[float]] = []
    for row in rows[1:]:
        if not row:
            continue
        identifiers.append(row[0])
        values.append([parse_float(value) for value in row[1:]])
    return identifiers, np.asarray(values, dtype=float)


def load_reads(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return {row["read_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def is_tumor(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def peel_complete(distance: np.ndarray) -> list[int]:
    indices = list(range(distance.shape[0]))
    while indices:
        sub = distance[np.ix_(indices, indices)]
        bad = (sub < 0) | np.isnan(sub)
        np.fill_diagonal(bad, False)
        if not bad.any():
            return indices
        indices.remove(indices[int(np.argmax(bad.sum(axis=1)))])
        if len(indices) < 2 * MIN_SIZE:
            return indices
    return indices


def link_tree(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    symmetric = distance.copy()
    np.fill_diagonal(symmetric, 0.0)
    symmetric = np.maximum(symmetric, symmetric.T)
    if np.any(~np.isfinite(symmetric)) or np.any(symmetric < 0):
        raise ValueError("Distance matrix is not finite and complete")
    return linkage(squareform(symmetric, checks=False), method="average"), symmetric


def forced_silhouette_split(distance: np.ndarray) -> dict[str, Any]:
    n_reads = distance.shape[0]
    if n_reads < 2 * MIN_SIZE:
        return {"k": 0, "silhouette": None, "labels": None}
    try:
        tree, symmetric = link_tree(distance)
    except (ValueError, FloatingPointError):
        return {"k": 0, "silhouette": None, "labels": None}
    best: tuple[int, np.ndarray, float] | None = None
    for groups in range(2, min(6, n_reads // MIN_SIZE) + 1):
        labels = fcluster(tree, groups, criterion="maxclust")
        sizes = Counter(labels)
        if len(sizes) < groups or min(sizes.values()) < MIN_SIZE:
            continue
        try:
            score = float(silhouette_score(symmetric, labels, metric="precomputed"))
        except ValueError:
            continue
        if best is None or score > best[2]:
            best = (groups, labels, score)
    if best is None:
        return {"k": 0, "silhouette": None, "labels": None}
    return {"k": best[0], "silhouette": best[2], "labels": best[1].astype(int).tolist()}


def bernoulli_distance(methylation: np.ndarray) -> np.ndarray:
    n_reads = methylation.shape[0]
    observed = (~np.isnan(methylation)).astype(float)
    filled = np.where(np.isnan(methylation), 0.0, methylation)
    weights = 2.0 * np.abs(filled - 0.5) * observed
    weighted_probability = weights * filled
    sum_weights = weights @ weights.T
    numerator = (
        weighted_probability @ weights.T
        + weights @ weighted_probability.T
        - 2.0 * (weighted_probability @ weighted_probability.T)
    )
    common = observed @ observed.T
    with np.errstate(divide="ignore", invalid="ignore"):
        distance = numerator / sum_weights
    distance[(common < MIN_SIZE) | (sum_weights < 1e-9)] = -1.0
    np.fill_diagonal(distance, 0.0)
    if distance.shape != (n_reads, n_reads):
        raise RuntimeError("Invalid BERNOULLI matrix shape")
    return distance


def tree_maps(distance: np.ndarray) -> tuple[np.ndarray, dict[int, list[int]], dict[int, tuple[int, int]], int]:
    tree, symmetric = link_tree(distance)
    n_reads = distance.shape[0]
    descendants: dict[int, list[int]] = {index: [index] for index in range(n_reads)}
    children: dict[int, tuple[int, int]] = {}
    for row_index, row in enumerate(tree):
        first, second = int(row[0]), int(row[1])
        node = n_reads + row_index
        descendants[node] = descendants[first] + descendants[second]
        children[node] = (first, second)
    return symmetric, descendants, children, 2 * n_reads - 2


def between_within_ratio(
    distance: np.ndarray, first: np.ndarray | list[int], second: np.ndarray | list[int]
) -> float | None:
    first = np.asarray(first, dtype=int)
    second = np.asarray(second, dtype=int)
    between = distance[np.ix_(first, second)]
    between = between[between >= 0]
    within_parts: list[np.ndarray] = []
    for group in (first, second):
        if len(group) >= 2:
            part = distance[np.ix_(group, group)][np.triu_indices(len(group), 1)]
            within_parts.append(part[part >= 0])
    within = np.concatenate(within_parts) if within_parts else np.asarray([], dtype=float)
    if not between.size or not within.size or float(within.mean()) <= 1e-6:
        return None
    return float(between.mean()) / float(within.mean())


def permute_methylation(
    methylation: np.ndarray, rng: np.random.Generator, mode: str = "column"
) -> np.ndarray:
    permuted = methylation.copy()
    if mode == "column":
        for column_index in range(permuted.shape[1]):
            valid = np.where(~np.isnan(permuted[:, column_index]))[0]
            if valid.size > 1:
                permuted[valid, column_index] = permuted[rng.permutation(valid), column_index]
        return permuted
    if mode == "row_circular":
        if permuted.shape[1] <= 1:
            return permuted
        for row_index in range(permuted.shape[0]):
            # Uniformly sample the complete circular-shift group, including identity.
            shift = int(rng.integers(0, permuted.shape[1]))
            permuted[row_index] = np.roll(permuted[row_index], shift)
        return permuted
    raise ValueError(f"Unsupported methylation null mode: {mode}")


def null_distribution(
    methylation: np.ndarray,
    leaves: list[int],
    rng: np.random.Generator,
    percentile: float,
    replicates: int,
    mode: str = "column",
) -> list[float]:
    subset = methylation[np.asarray(leaves, dtype=int)]
    null_ratios: list[float] = []
    for _ in range(replicates):
        permuted = permute_methylation(subset, rng, mode=mode)
        null_distance = bernoulli_distance(permuted)
        try:
            symmetric, descendants, children, root = tree_maps(null_distance)
        except (ValueError, FloatingPointError):
            continue
        if root not in children:
            continue
        first, second = children[root]
        ratio = between_within_ratio(symmetric, descendants[first], descendants[second])
        if ratio is not None:
            null_ratios.append(ratio)
    return null_ratios


def phylo_label(
    distance: np.ndarray,
    methylation: np.ndarray,
    rng: np.random.Generator,
    null_pct: float = 95.0,
    rnull: int = RNULL,
    null_mode: str = "column",
    empirical_alpha: float | None = None,
    min_valid_null_fraction: float = 0.8,
    trace: list[dict[str, Any]] | None = None,
) -> list[str]:
    symmetric, descendants, children, root = tree_maps(distance)
    n_reads = distance.shape[0]
    labels: list[str | None] = [None] * n_reads

    def descend(node: int) -> tuple[int | None, list[int]]:
        current = node
        quarantined: list[int] = []
        while current in children:
            first, second = children[current]
            first_size, second_size = len(descendants[first]), len(descendants[second])
            if min(first_size, second_size) >= MIN_SIZE:
                return current, quarantined
            small, large = (first, second) if first_size < second_size else (second, first)
            quarantined.extend(descendants[small])
            current = large
        return None, quarantined

    def passes(node: int) -> bool:
        first, second = children[node]
        observed_ratio = between_within_ratio(symmetric, descendants[first], descendants[second])
        if observed_ratio is None or observed_ratio < SEP_MIN:
            if trace is not None:
                trace.append(
                    {
                        "n_node": len(descendants[node]),
                        "child_sizes": [len(descendants[first]), len(descendants[second])],
                        "observed_between_within": observed_ratio,
                        "null_percentile": null_pct,
                        "null_threshold": None,
                        "empirical_p": None,
                        "n_valid_null": 0,
                        "passed": False,
                        "failure": "below_sep_min_or_undefined",
                    }
                )
            return False
        null_ratios = null_distribution(
            methylation, descendants[node], rng, null_pct, rnull, mode=null_mode
        )
        minimum_valid = max(1, int(np.ceil(rnull * min_valid_null_fraction)))
        if len(null_ratios) < minimum_valid:
            if trace is not None:
                trace.append(
                    {
                        "n_node": len(descendants[node]),
                        "child_sizes": [len(descendants[first]), len(descendants[second])],
                        "observed_between_within": observed_ratio,
                        "null_percentile": null_pct,
                        "null_threshold": None,
                        "empirical_p": None,
                        "n_valid_null": len(null_ratios),
                        "minimum_valid_null": minimum_valid,
                        "null_mode": null_mode,
                        "passed": False,
                        "failure": "insufficient_valid_null",
                    }
                )
            return False
        threshold = float(np.percentile(null_ratios, null_pct))
        empirical_p = (1 + sum(value >= observed_ratio for value in null_ratios)) / (len(null_ratios) + 1)
        passed_threshold = observed_ratio > threshold
        passed_empirical = empirical_alpha is None or empirical_p <= empirical_alpha
        passed = passed_threshold and passed_empirical
        if trace is not None:
            trace.append(
                {
                    "n_node": len(descendants[node]),
                    "child_sizes": [len(descendants[first]), len(descendants[second])],
                    "observed_between_within": observed_ratio,
                    "null_percentile": null_pct,
                    "null_threshold": threshold,
                    "empirical_p": empirical_p,
                    "n_valid_null": len(null_ratios),
                    "null_mode": null_mode,
                    "passed": passed,
                    "failure": (
                        None
                        if passed
                        else "not_above_null_threshold"
                        if not passed_threshold
                        else "empirical_p_above_alpha"
                    ),
                }
            )
        return passed

    def recurse(node: int, label: str) -> None:
        leaves = descendants[node]
        if len(leaves) < 2 * MIN_SIZE:
            for index in leaves:
                labels[index] = label
            return
        balanced_node, quarantined = descend(node)
        if balanced_node is not None and passes(balanced_node):
            for index in quarantined:
                labels[index] = "outlier"
            first, second = children[balanced_node]
            large, small = (
                (first, second)
                if len(descendants[first]) >= len(descendants[second])
                else (second, first)
            )
            recurse(large, label + "-1")
            recurse(small, label + "-2")
        else:
            for index in leaves:
                labels[index] = label

    recurse(root, "1")
    concrete = [label or "outlier" for label in labels]
    sizes = Counter(label for label in concrete if label != "outlier")
    too_small = {label for label, size in sizes.items() if size < MIN_SIZE}
    concrete = ["outlier" if label in too_small else label for label in concrete]
    if concrete.count("outlier") >= MIN_SIZE:
        concrete = ["other" if label == "outlier" else label for label in concrete]
    return concrete


def number_of_groups(labels: Iterable[str]) -> int:
    return len({label for label in labels if label not in {"outlier", "other"}})


def analyze_phylo(
    distance: np.ndarray,
    methylation: np.ndarray,
    base_seed: int = 20260622,
    seeds: int = 10,
    rnull: int = RNULL,
    null_mode: str = "column",
    empirical_alpha: float | None = None,
    min_valid_null_fraction: float = 0.8,
) -> dict[str, Any]:
    coarse_runs: list[list[str]] = []
    coarse_counts: list[int] = []
    coarse_traces: list[list[dict[str, Any]]] = []
    for seed_index in range(seeds):
        run_trace: list[dict[str, Any]] = []
        labels = phylo_label(
            distance,
            methylation,
            np.random.default_rng(base_seed + seed_index * 101),
            null_pct=95.0,
            rnull=rnull,
            null_mode=null_mode,
            empirical_alpha=empirical_alpha,
            min_valid_null_fraction=min_valid_null_fraction,
            trace=run_trace,
        )
        coarse_runs.append(labels)
        coarse_counts.append(number_of_groups(labels))
        coarse_traces.append(run_trace)
    frequencies = Counter(coarse_counts)
    modal_groups, modal_count = frequencies.most_common(1)[0]
    representative_index = next(index for index, count in enumerate(coarse_counts) if count == modal_groups)
    representative = coarse_runs[representative_index]
    fine_trace: list[dict[str, Any]] = []
    fine_labels = phylo_label(
        distance,
        methylation,
        np.random.default_rng(base_seed),
        null_pct=90.0,
        rnull=rnull,
        null_mode=null_mode,
        empirical_alpha=empirical_alpha,
        min_valid_null_fraction=min_valid_null_fraction,
        trace=fine_trace,
    )
    modal_fraction = modal_count / seeds
    n_other = representative.count("other")
    modal_indices = [index for index, count in enumerate(coarse_counts) if count == modal_groups]
    assignment_aris = [
        float(adjusted_rand_score(coarse_runs[first], coarse_runs[second]))
        for offset, first in enumerate(modal_indices)
        for second in modal_indices[offset + 1 :]
    ]
    if assignment_aris:
        assignment_ari_median = float(np.median(assignment_aris))
        assignment_ari_min = min(assignment_aris)
    else:
        assignment_ari_median = 1.0
        assignment_ari_min = 1.0
    return {
        "coarse_ng": modal_groups,
        "modal_fraction": modal_fraction,
        "fine_ng": number_of_groups(fine_labels),
        "n_other": n_other,
        "n_outlier": representative.count("outlier"),
        "unstable": modal_fraction < MODAL_CONFIDENCE,
        "ng_min": min(coarse_counts),
        "ng_max": max(coarse_counts),
        "seed_group_counts": coarse_counts,
        "modal_assignment_pair_count": len(assignment_aris),
        "modal_assignment_ari_median": assignment_ari_median,
        "modal_assignment_ari_min": assignment_ari_min,
        "modal_assignment_all_pairs_ari_ge_0_8": assignment_ari_min >= 0.8,
        "hidden_heterogeneity": n_other / len(representative) > HIDDEN_HETEROGENEITY_FRACTION,
        "coarse_labels": representative,
        "fine_labels": fine_labels,
        "coarse_split_trace": coarse_traces[representative_index],
        "fine_split_trace": fine_trace,
    }


def cramer_v(table: np.ndarray) -> float:
    if table.shape[0] < 2 or table.shape[1] < 2 or table.sum() <= 0:
        return 0.0
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / table.sum()
    valid = expected > 0
    chi_square = float((((table - expected) ** 2 / np.where(valid, expected, 1.0)) * valid).sum())
    denominator = table.sum() * min(table.shape[0] - 1, table.shape[1] - 1)
    return float(np.sqrt(chi_square / denominator)) if denominator > 0 else 0.0


def categorical_permutation_association(
    values: list[str], labels: list[str], seed: int, permutations: int = 499
) -> dict[str, Any]:
    pairs = [(value, label) for value, label in zip(values, labels) if value not in {"", ".", "NA"}]
    categories = sorted({value for value, _ in pairs})
    groups = sorted({label for _, label in pairs})
    if len(categories) < 2 or len(groups) < 2:
        return {"v": None, "p_perm": None, "n": len(pairs), "aligned": False}
    category_index = {value: index for index, value in enumerate(categories)}
    group_index = {value: index for index, value in enumerate(groups)}

    def make_table(group_values: list[str]) -> np.ndarray:
        table = np.zeros((len(categories), len(groups)), dtype=int)
        for (category, _), group in zip(pairs, group_values):
            table[category_index[category], group_index[group]] += 1
        return table

    group_values = [label for _, label in pairs]
    observed = cramer_v(make_table(group_values))
    rng = np.random.default_rng(seed)
    exceed = 0
    array = np.asarray(group_values, dtype=object)
    for _ in range(permutations):
        if cramer_v(make_table(rng.permutation(array).tolist())) >= observed - 1e-12:
            exceed += 1
    p_value = (exceed + 1) / (permutations + 1)
    return {
        "v": observed,
        "p_perm": p_value,
        "n": len(pairs),
        "aligned": observed >= 0.30 and p_value < 0.05,
    }


def eta_squared(values: np.ndarray, labels: list[str]) -> float:
    finite = np.isfinite(values)
    values = values[finite]
    labels_array = np.asarray(labels, dtype=object)[finite]
    if values.size < 2 or len(set(labels_array.tolist())) < 2:
        return 0.0
    grand_mean = float(values.mean())
    total = float(np.sum((values - grand_mean) ** 2))
    if total <= 1e-12:
        return 0.0
    between = 0.0
    for label in set(labels_array.tolist()):
        group = values[labels_array == label]
        between += len(group) * (float(group.mean()) - grand_mean) ** 2
    return between / total


def continuous_permutation_association(
    values: list[float], labels: list[str], seed: int, permutations: int = 499
) -> dict[str, Any]:
    array = np.asarray(values, dtype=float)
    observed = eta_squared(array, labels)
    if np.isfinite(array).sum() < 2 or len(set(labels)) < 2:
        return {"eta2": None, "p_perm": None, "n": int(np.isfinite(array).sum()), "aligned": False}
    rng = np.random.default_rng(seed)
    exceed = 0
    label_array = np.asarray(labels, dtype=object)
    for _ in range(permutations):
        if eta_squared(array, rng.permutation(label_array).tolist()) >= observed - 1e-12:
            exceed += 1
    p_value = (exceed + 1) / (permutations + 1)
    return {
        "eta2": observed,
        "p_perm": p_value,
        "n": int(np.isfinite(array).sum()),
        "aligned": observed >= 0.14 and p_value < 0.05,
    }


def hp_family(tag: str) -> str:
    tag = str(tag).strip()
    if tag in {"1", "HP1", "1-1", "1-2"}:
        return "HP1-side"
    if tag in {"2", "HP2", "2-1", "2-2"}:
        return "HP2-side"
    if tag == "3":
        return "HP3-ambiguous"
    if tag == "4":
        return "HP4-both"
    return "untagged"


def stable_seed(sample: str, chrom: str, position: int, offset: int = 0) -> int:
    payload = f"{sample}|{chrom}|{position}|{offset}".encode()
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "big") % (2**32)
