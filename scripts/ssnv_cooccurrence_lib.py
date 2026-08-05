#!/usr/bin/env python3
"""Read-level sSNV co-occurrence and fixed-margin inference primitives."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pysam
from scipy.stats import beta, binomtest, chi2_contingency, fisher_exact


DEFAULT_EXACT_STATE_SPACE_CEILING = 250_000
FOUR_STATE_ERROR_CEILING = 0.02
FOUR_STATE_CONFIDENCE = 0.95
FOUR_STATE_MIN_SPLIT_COUNT = 3
FOUR_STATE_RELATION_FAMILY_SIZE = 3
FOUR_STATE_MULTIPLICITY_METHOD = "bonferroni_three_relation_models"


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str


def call_snv(alignment: pysam.AlignedSegment, variant: Variant, min_base_quality: int = 20) -> str:
    """Return R/A/O/X for one biallelic SNV on an already selected alignment."""
    target = variant.pos - 1
    query_position: int | None = None
    for query_pos, reference_pos in alignment.get_aligned_pairs(matches_only=False):
        if reference_pos == target:
            query_position = query_pos
            break
        if reference_pos is not None and reference_pos > target:
            break
    if query_position is None:
        return "X"
    sequence = alignment.query_sequence or ""
    qualities = alignment.query_qualities
    if query_position >= len(sequence) or qualities is None or query_position >= len(qualities):
        return "X"
    if int(qualities[query_position]) < min_base_quality:
        return "X"
    base = sequence[query_position].upper()
    if base == variant.ref.upper():
        return "R"
    if base == variant.alt.upper():
        return "A"
    return "O"


def call_snvs(
    alignment: pysam.AlignedSegment,
    variants: Sequence[Variant],
    min_base_quality: int = 20,
) -> dict[int, str]:
    """Call many markers with one CIGAR traversal; result is keyed by one-based position."""
    targets = {variant.pos - 1: variant for variant in variants}
    query_by_reference: dict[int, int] = {}
    for query_pos, reference_pos in alignment.get_aligned_pairs(matches_only=True):
        if reference_pos in targets and query_pos is not None:
            query_by_reference[int(reference_pos)] = int(query_pos)
    sequence = alignment.query_sequence or ""
    qualities = alignment.query_qualities
    calls: dict[int, str] = {}
    for target, variant in targets.items():
        query_position = query_by_reference.get(target)
        if (
            query_position is None
            or qualities is None
            or query_position >= len(sequence)
            or query_position >= len(qualities)
            or int(qualities[query_position]) < min_base_quality
        ):
            calls[variant.pos] = "X"
            continue
        base = sequence[query_position].upper()
        calls[variant.pos] = "R" if base == variant.ref.upper() else "A" if base == variant.alt.upper() else "O"
    return calls


def cramer_v(table: np.ndarray) -> float:
    if table.ndim != 2 or min(table.shape) < 2 or table.sum() <= 0:
        return 0.0
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / table.sum()
    valid = expected > 0
    chi_square = float((((table - expected) ** 2 / np.where(valid, expected, 1.0)) * valid).sum())
    denominator = float(table.sum() * min(table.shape[0] - 1, table.shape[1] - 1))
    return math.sqrt(chi_square / denominator) if denominator > 0 else 0.0


def contingency_table(
    labels: Sequence[str], categories: Sequence[str]
) -> tuple[list[str], list[str], np.ndarray]:
    if len(labels) != len(categories):
        raise ValueError("labels/categories length mismatch")
    groups = sorted(set(labels))
    states = sorted(set(categories))
    table = np.zeros((len(groups), len(states)), dtype=int)
    group_index = {value: index for index, value in enumerate(groups)}
    state_index = {value: index for index, value in enumerate(states)}
    for label, category in zip(labels, categories):
        table[group_index[label], state_index[category]] += 1
    return groups, states, table


def group_allele_table(labels: Sequence[str], calls: Sequence[str]) -> tuple[list[str], np.ndarray]:
    pairs = [(label, call) for label, call in zip(labels, calls) if call in {"R", "A"}]
    if not pairs:
        return [], np.zeros((0, 2), dtype=int)
    groups = sorted({label for label, _ in pairs})
    index = {group: offset for offset, group in enumerate(groups)}
    table = np.zeros((len(groups), 2), dtype=int)
    for label, call in pairs:
        table[index[label], 0 if call == "R" else 1] += 1
    return groups, table


def _direction_summary(groups: Sequence[str], table: np.ndarray) -> dict[str, object]:
    row_totals = table.sum(axis=1)
    if table.shape[0] < 2 or np.any(row_totals <= 0):
        return {
            "dominant_group": None,
            "direction_contrast": None,
            "direction": "NOT_IDENTIFIABLE_DEGENERATE_TABLE",
        }
    largest = int(row_totals.max())
    dominant_indices = np.flatnonzero(row_totals == largest)
    if len(dominant_indices) != 1:
        return {
            "dominant_group": None,
            "direction_contrast": None,
            "direction": "NOT_IDENTIFIABLE_TIED_DOMINANT_GROUP",
        }
    dominant = int(dominant_indices[0])
    other_total = int(row_totals.sum() - row_totals[dominant])
    if other_total <= 0:
        return {
            "dominant_group": None,
            "direction_contrast": None,
            "direction": "NOT_IDENTIFIABLE_DEGENERATE_TABLE",
        }
    dominant_alt = float(table[dominant, 1] / row_totals[dominant])
    other_alt = float((table[:, 1].sum() - table[dominant, 1]) / other_total)
    contrast = dominant_alt - other_alt
    if contrast > 1e-12:
        direction = "PARTNER_ALT_ENRICHED_IN_DOMINANT_METHYL_GROUP"
    elif contrast < -1e-12:
        direction = "PARTNER_ALT_DEPLETED_IN_DOMINANT_METHYL_GROUP"
    else:
        direction = "NO_DIRECTIONAL_CONTRAST"
    return {
        "dominant_group": groups[dominant],
        "direction_contrast": contrast,
        "direction": direction,
    }


def methyl_group_marker_association(
    labels: Sequence[str],
    calls: Sequence[str],
    min_total: int = 10,
    min_group: int = 3,
    min_allele: int = 3,
) -> dict[str, object]:
    """Build a Kx2 table; analytic p is descriptive and never a formal gate."""
    groups, table = group_allele_table(labels, calls)
    informative = int(table.sum())
    result: dict[str, object] = {
        "testable": False,
        "reason": "insufficient_joint_information",
        "n_informative": informative,
        "groups": groups,
        "table": table.tolist(),
        "cramers_v": None,
        "delta_alt_fraction": None,
        "p_analytic": None,
        "analytic_test": None,
        **_direction_summary(groups, table),
    }
    if len(groups) < 2 or informative < min_total:
        return result
    if int(table.sum(axis=1).min()) < min_group:
        result["reason"] = "methyl_group_below_minimum"
        return result
    if int(table.sum(axis=0).min()) < min_allele:
        result["reason"] = "partner_allele_below_minimum"
        return result
    alt_fractions = table[:, 1] / table.sum(axis=1)
    if table.shape == (2, 2):
        p_value = float(fisher_exact(table)[1])
        test = "fisher_exact_2x2_descriptive_only"
    else:
        p_value = float(chi2_contingency(table, correction=False)[1])
        test = "pearson_chi_square_descriptive_only"
    result.update(
        {
            "testable": True,
            "reason": "ok",
            "cramers_v": cramer_v(table),
            "delta_alt_fraction": float(alt_fractions.max() - alt_fractions.min()),
            "p_analytic": p_value,
            "analytic_test": test,
            **_direction_summary(groups, table),
        }
    )
    return result


def _log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


class _StateSpaceLimitExceeded(RuntimeError):
    pass


def fisher_freeman_halton_kx2(
    table: Sequence[Sequence[int]] | np.ndarray,
    state_space_ceiling: int = DEFAULT_EXACT_STATE_SPACE_CEILING,
) -> dict[str, object]:
    """Probability-ordered, two-sided fixed-margins exact test for a Kx2 table.

    Each state is the vector of ALT counts per methyl group. Row totals and the
    total ALT count are fixed, so state weights are proportional to the product
    of ``choose(row_total, alt_count)`` terms.
    """
    observed = np.asarray(table)
    if observed.ndim != 2 or observed.shape[1] != 2:
        raise ValueError("Fisher-Freeman-Halton Kx2 test requires a two-column table")
    if state_space_ceiling < 1:
        raise ValueError("state_space_ceiling must be positive")
    if not np.issubdtype(observed.dtype, np.integer):
        if not np.all(np.equal(observed, np.floor(observed))):
            raise ValueError("contingency table must contain integer counts")
        observed = observed.astype(np.int64)
    else:
        observed = observed.astype(np.int64, copy=False)
    if np.any(observed < 0):
        raise ValueError("contingency table cannot contain negative counts")
    observed = observed[observed.sum(axis=1) > 0]
    column_totals = observed.sum(axis=0) if observed.size else np.asarray([0, 0])
    base = {
        "method": "fisher_freeman_halton_fixed_margins_probability_ordering",
        "p_value": None,
        "identifiable": False,
        "state_space_status": "NOT_IDENTIFIABLE_DEGENERATE_TABLE",
        "state_space_size": None,
        "state_space_lower_bound": None,
        "state_space_ceiling": state_space_ceiling,
        "fallback": "fail_closed_not_identifiable",
        "observed_table_probability": None,
    }
    if observed.shape[0] < 2 or np.any(column_totals <= 0):
        return base

    row_totals = observed.sum(axis=1).astype(int).tolist()
    observed_alt = observed[:, 1].astype(int).tolist()
    target_alt = int(column_totals[1])
    log_choose_by_row = [
        [_log_choose(row_total, alt_count) for alt_count in range(row_total + 1)]
        for row_total in row_totals
    ]
    observed_log_weight = math.fsum(
        log_choose_by_row[index][alt_count]
        for index, alt_count in enumerate(observed_alt)
    )
    suffix_capacity = [0] * (len(row_totals) + 1)
    for index in range(len(row_totals) - 1, -1, -1):
        suffix_capacity[index] = suffix_capacity[index + 1] + row_totals[index]

    state_count = 0
    total_log_weight = -math.inf
    tail_log_weight = -math.inf

    def visit(log_weight: float) -> None:
        nonlocal state_count, total_log_weight, tail_log_weight
        state_count += 1
        if state_count > state_space_ceiling:
            raise _StateSpaceLimitExceeded
        total_log_weight = float(np.logaddexp(total_log_weight, log_weight))
        if log_weight <= observed_log_weight + 1e-12:
            tail_log_weight = float(np.logaddexp(tail_log_weight, log_weight))

    def enumerate_states(index: int, remaining_alt: int, log_weight: float) -> None:
        row_total = row_totals[index]
        if index == len(row_totals) - 1:
            if 0 <= remaining_alt <= row_total:
                visit(log_weight + log_choose_by_row[index][remaining_alt])
            return
        remaining_capacity = suffix_capacity[index + 1]
        lower = max(0, remaining_alt - remaining_capacity)
        upper = min(row_total, remaining_alt)
        for alt_count in range(lower, upper + 1):
            enumerate_states(
                index + 1,
                remaining_alt - alt_count,
                log_weight + log_choose_by_row[index][alt_count],
            )

    try:
        enumerate_states(0, target_alt, 0.0)
    except _StateSpaceLimitExceeded:
        return {
            **base,
            "state_space_status": "NOT_IDENTIFIABLE_STATE_SPACE_LIMIT",
            "state_space_lower_bound": state_space_ceiling + 1,
        }
    if state_count == 0 or not math.isfinite(total_log_weight):
        return base
    p_value = min(1.0, math.exp(tail_log_weight - total_log_weight))
    observed_probability = min(1.0, math.exp(observed_log_weight - total_log_weight))
    return {
        **base,
        "p_value": p_value,
        "identifiable": True,
        "state_space_status": "EXACT_ENUMERATED",
        "state_space_size": state_count,
        "observed_table_probability": observed_probability,
    }


def categorical_association_table(
    labels: Sequence[str],
    categories: Sequence[str],
    min_total: int = 10,
    min_group: int = 3,
    min_category: int = 3,
) -> dict[str, object]:
    groups, states, table = contingency_table(labels, categories)
    result: dict[str, object] = {
        "testable": False,
        "reason": "insufficient_joint_information",
        "groups": groups,
        "categories": states,
        "table": table.tolist(),
        "n": int(table.sum()),
        "cramers_v": None,
    }
    if len(groups) < 2 or len(states) < 2 or int(table.sum()) < min_total:
        return result
    if int(table.sum(axis=1).min()) < min_group:
        result["reason"] = "group_below_minimum"
        return result
    if int(table.sum(axis=0).min()) < min_category:
        result["reason"] = "category_below_minimum"
        return result
    result.update({"testable": True, "reason": "ok", "cramers_v": cramer_v(table)})
    return result


def conditional_label_permutation(
    labels: Sequence[str],
    categories: Sequence[str],
    strata: Sequence[str],
    seed: int,
    permutations: int = 999,
) -> dict[str, object]:
    """Permute methyl labels within fixed strata using Cramer's V as statistic."""
    if not (len(labels) == len(categories) == len(strata)):
        raise ValueError("labels/categories/strata length mismatch")
    if permutations < 1:
        raise ValueError("permutations must be positive")
    groups, states, observed_table = contingency_table(labels, categories)
    base: dict[str, object] = {
        "p_conditional_perm": None,
        "permutations": 0,
        "permutable": False,
        "status": "NOT_IDENTIFIABLE_DEGENERATE_TABLE",
        "n_strata": len(set(strata)),
        "n_exchangeable_strata": 0,
        "n_informative_exchangeable_strata": 0,
    }
    if len(groups) < 2 or len(states) < 2:
        return base
    labels_array = np.asarray(labels, dtype=object)
    indices_by_stratum: dict[str, list[int]] = defaultdict(list)
    for index, stratum in enumerate(strata):
        indices_by_stratum[stratum].append(index)
    exchangeable = [
        indices
        for indices in indices_by_stratum.values()
        if len(indices) >= 2 and len({labels_array[index] for index in indices}) >= 2
    ]
    informative = [
        indices
        for indices in exchangeable
        if len({categories[index] for index in indices}) >= 2
    ]
    base["n_exchangeable_strata"] = len(exchangeable)
    base["n_informative_exchangeable_strata"] = len(informative)
    if not informative:
        return {**base, "status": "NOT_IDENTIFIABLE_NONEXCHANGEABLE"}

    group_index = {value: index for index, value in enumerate(groups)}
    state_index = {value: index for index, value in enumerate(states)}

    def table_for(candidate_labels: np.ndarray) -> np.ndarray:
        table = np.zeros((len(groups), len(states)), dtype=int)
        for label, category in zip(candidate_labels, categories):
            table[group_index[str(label)], state_index[category]] += 1
        return table

    observed = cramer_v(observed_table)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(permutations):
        shuffled = labels_array.copy()
        for indices in exchangeable:
            shuffled[indices] = rng.permutation(shuffled[indices])
        if cramer_v(table_for(shuffled)) >= observed - 1e-12:
            exceed += 1
    return {
        **base,
        "p_conditional_perm": (exceed + 1) / (permutations + 1),
        "permutations": permutations,
        "permutable": True,
        "status": "PERMUTABLE",
    }


def _binomial_violation_gate(
    violations: int,
    denominator: int,
    error_ceiling: float,
    confidence: float,
) -> dict[str, object]:
    result: dict[str, object] = {
        "violations": violations,
        "denominator": denominator,
        "rate": violations / denominator if denominator else None,
        "p_exact_greater": None,
        "upper_bound": None,
        "threshold": error_ceiling,
        "confidence": confidence,
        "status": "NOT_IDENTIFIABLE_NO_DENOMINATOR",
    }
    if denominator <= 0:
        return result
    p_value = float(binomtest(violations, denominator, error_ceiling, alternative="greater").pvalue)
    upper = (
        1.0
        if violations == denominator
        else float(beta.ppf(confidence, violations + 1, denominator - violations))
    )
    alpha = 1.0 - confidence
    if upper <= error_ceiling + 1e-15:
        status = "COMPATIBLE_WITH_FIXED_ERROR_CEILING"
    elif p_value <= alpha + 1e-15:
        status = "VIOLATES_FIXED_ERROR_CEILING"
    else:
        status = "NOT_IDENTIFIABLE_LOW_PRECISION"
    return {
        **result,
        "p_exact_greater": p_value,
        "upper_bound": upper,
        "status": status,
    }


def four_state_summary(
    focal_calls: Sequence[str],
    partner_calls: Sequence[str],
    error_ceiling: float = FOUR_STATE_ERROR_CEILING,
    confidence: float = FOUR_STATE_CONFIDENCE,
) -> dict[str, object]:
    if len(focal_calls) != len(partner_calls):
        raise ValueError("focal_calls/partner_calls length mismatch")
    if not 0.0 < error_ceiling < 1.0 or not 0.0 < confidence < 1.0:
        raise ValueError("invalid fixed-error model parameters")
    per_relation_confidence = 1.0 - (1.0 - confidence) / FOUR_STATE_RELATION_FAMILY_SIZE
    counts: Counter[str] = Counter()
    for focal, partner in zip(focal_calls, partner_calls):
        if focal not in {"R", "A"}:
            continue
        if partner in {"R", "A"}:
            counts[f"{focal}{partner}"] += 1
        elif partner == "O":
            counts["O"] += 1
        else:
            counts["X"] += 1
    rr, ar, ra, aa = (counts[key] for key in ("RR", "AR", "RA", "AA"))
    called = rr + ar + ra + aa
    focal_ref = rr + ra
    focal_alt = ar + aa
    partner_ref = rr + ar
    partner_alt = ra + aa
    focal_gate = _binomial_violation_gate(
        ra, partner_alt, error_ceiling, per_relation_confidence
    )
    partner_gate = _binomial_violation_gate(
        ar, focal_alt, error_ceiling, per_relation_confidence
    )
    branching_denominator = ar + ra + aa
    branching_gate = _binomial_violation_gate(
        aa, branching_denominator, error_ceiling, per_relation_confidence
    )
    complete = called >= 10 and min(focal_ref, focal_alt, partner_ref, partner_alt) >= 3
    focal_split = ar >= FOUR_STATE_MIN_SPLIT_COUNT and aa >= FOUR_STATE_MIN_SPLIT_COUNT
    partner_split = ra >= FOUR_STATE_MIN_SPLIT_COUNT and aa >= FOUR_STATE_MIN_SPLIT_COUNT
    branching_split = ar >= FOUR_STATE_MIN_SPLIT_COUNT and ra >= FOUR_STATE_MIN_SPLIT_COUNT
    focal_compatible = complete and focal_split and focal_gate["status"] == "COMPATIBLE_WITH_FIXED_ERROR_CEILING"
    partner_compatible = (
        complete and partner_split and partner_gate["status"] == "COMPATIBLE_WITH_FIXED_ERROR_CEILING"
    )
    branching_compatible = (
        complete and branching_split and branching_gate["status"] == "COMPATIBLE_WITH_FIXED_ERROR_CEILING"
    )
    relevant_statuses = []
    if focal_split:
        relevant_statuses.append(str(focal_gate["status"]))
    if partner_split:
        relevant_statuses.append(str(partner_gate["status"]))
    if branching_split:
        relevant_statuses.append(str(branching_gate["status"]))
    compatible_models = []
    if focal_compatible:
        compatible_models.append("FOCAL_ANCESTOR")
    if partner_compatible:
        compatible_models.append("PARTNER_ANCESTOR")
    if branching_compatible:
        compatible_models.append("BRANCHING")
    if not complete:
        relation = "NOT_IDENTIFIABLE_INSUFFICIENT_FOUR_STATE_DEPTH"
    elif len(compatible_models) > 1:
        relation = "MULTIPLE_MUTATION_ORDER_MODELS_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    elif compatible_models == ["FOCAL_ANCESTOR"]:
        relation = "FOCAL_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    elif compatible_models == ["PARTNER_ANCESTOR"]:
        relation = "PARTNER_ANCESTOR_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    elif compatible_models == ["BRANCHING"]:
        relation = "BRANCHING_COMPATIBLE_UNDER_FIXED_ERROR_MODEL"
    elif not relevant_statuses or any(status.startswith("NOT_IDENTIFIABLE") for status in relevant_statuses):
        relation = "NOT_IDENTIFIABLE_FIXED_ERROR_CEILING"
    else:
        relation = "INCOMPATIBLE_OR_COMPLEX_UNDER_FIXED_ERROR_MODEL"
    minimum_zero_depth = math.ceil(
        math.log(1.0 - per_relation_confidence) / math.log(1.0 - error_ceiling)
    )
    return {
        "n_joint": called,
        "n_called_depth": called,
        "state_counts": {key: counts[key] for key in ("RR", "AR", "RA", "AA")},
        "state_counts_with_noncallable": {
            key: counts[key] for key in ("RR", "AR", "RA", "AA", "O", "X")
        },
        "n_focal_ref": focal_ref,
        "n_focal_alt": focal_alt,
        "n_partner_ref": partner_ref,
        "n_partner_alt": partner_alt,
        "focal_ancestor_violation_rate": focal_gate["rate"],
        "partner_ancestor_violation_rate": partner_gate["rate"],
        "focal_ancestor_violation": focal_gate,
        "partner_ancestor_violation": partner_gate,
        "branching_violation": branching_gate,
        "error_ceiling": error_ceiling,
        "confidence": per_relation_confidence,
        "familywise_confidence": confidence,
        "relation_family_size": FOUR_STATE_RELATION_FAMILY_SIZE,
        "multiplicity_method": FOUR_STATE_MULTIPLICITY_METHOD,
        "minimum_zero_violation_depth": minimum_zero_depth,
        "complete_four_state_testable": complete,
        "relation_compatibility": relation,
        "compatible_relation_models": compatible_models,
        "n_compatible_relation_models": len(compatible_models),
        "claim_guardrail": (
            "One or more mutation-order models may be compatible under the Bonferroni-adjusted "
            "familywise fixed-error model; compatibility is not proof of cellular ancestry or an "
            "identifiable unique order."
        ),
    }


def benjamini_hochberg(p_values: Iterable[float | None]) -> list[float | None]:
    values = list(p_values)
    valid = [(index, float(value)) for index, value in enumerate(values) if value is not None and math.isfinite(value)]
    valid.sort(key=lambda pair: pair[1])
    adjusted: list[float | None] = [None] * len(values)
    running = 1.0
    total = len(valid)
    for rank_index in range(total - 1, -1, -1):
        original_index, p_value = valid[rank_index]
        rank = rank_index + 1
        running = min(running, p_value * total / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def benjamini_yekutieli(p_values: Iterable[float | None]) -> list[float | None]:
    values = list(p_values)
    n_valid = sum(value is not None and math.isfinite(float(value)) for value in values)
    if n_valid == 0:
        return [None] * len(values)
    harmonic = math.fsum(1.0 / rank for rank in range(1, n_valid + 1))
    return [
        None if q_value is None else min(1.0, float(q_value) * harmonic)
        for q_value in benjamini_hochberg(values)
    ]


def spatially_separated_positions(
    positions: Iterable[int], minimum_separation: int = 20
) -> list[int]:
    if minimum_separation < 1:
        raise ValueError("minimum_separation must be positive")
    selected: list[int] = []
    for position in sorted(set(int(value) for value in positions)):
        if not selected or position - selected[-1] >= minimum_separation:
            selected.append(position)
    return selected


def count_spatially_separated_positions(positions: Iterable[int], minimum_separation: int = 20) -> int:
    return len(spatially_separated_positions(positions, minimum_separation))


def count_independent_positions(positions: Iterable[int], minimum_separation: int = 20) -> int:
    """Deprecated alias; spacing does not establish statistical independence."""
    return count_spatially_separated_positions(positions, minimum_separation)
