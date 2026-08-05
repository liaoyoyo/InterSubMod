#!/usr/bin/env python3
"""Fail-closed M2 measured-axis gate shared by producer and final integrator."""

from __future__ import annotations

import math
import json
import warnings
from functools import lru_cache
from typing import Any, Mapping

from scipy.stats import chi2, f, ncf, ncx2


GATE_CONTRACT = (
    "m2-measured-axis-v4_asymmetric-positive-confound-and-observed-categorical-levels"
)
P_MAX = 0.05
PERMUTATIONS = 499
PERMUTATION_P_FLOOR = 1.0 / (PERMUTATIONS + 1)
TARGET_POWER = 0.80
POWER_ALPHA = 0.05
MIN_GROUP_N = 5
MAX_METHYL_GROUPS = 10

# The source-locked screen was planned with these level ceilings. Keep that
# planning model unchanged; separately require observed levels from the stable
# assignment to prove that a missing categorical statistic means a constant axis.
CATEGORICAL_LEVEL_CEILINGS = {
    "hp_exact": 7,
    "hp_family": 5,
    "strand": 2,
}
CATEGORICAL_LEVEL_KEYS = tuple(CATEGORICAL_LEVEL_CEILINGS)

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

HIGH_LEVEL_FIELDS = (
    "stable_null_multigroup",
    "modal_assignment_ari_min",
    "cluster_sizes",
    "hp_axis_confound",
    "technical_axis_confound",
    "residual_unexplained_multigroup",
    "phase_anchored_robust_epigenetic_candidate",
)

AXIS_SOURCE_FIELDS = tuple(
    field
    for prefix, _kind, effect_suffix, _threshold in AXIS_SPECS
    for field in (
        f"{prefix}_{effect_suffix}",
        f"{prefix}_p_perm",
        f"{prefix}_n",
        f"{prefix}_aligned",
    )
)

INDETERMINATE_AXIS_STATUS = "INDETERMINATE_EFFECT_ABOVE_THRESHOLD_WITHOUT_P_LT_0_05"
LOW_POWER_AXIS_STATUS = "INDETERMINATE_INSUFFICIENT_INFORMATION_FOR_EFFECT_THRESHOLD"
MISSING_STATISTIC_AXIS_STATUS = "INDETERMINATE_AXIS_STATISTIC_MISSING"
CONSTANT_AXIS_STATUS = "NOT_ALIGNED_AXIS_HAS_NO_OBSERVED_VARIATION"
LOW_EFFECT_AXIS_STATUS = "NOT_ALIGNED_EFFECT_BELOW_THRESHOLD_WITH_ADEQUATE_POWER"
ALIGNED_AXIS_STATUS = "ALIGNED_EFFECT_AND_PERMUTATION_P_PASS"


class M2GateError(ValueError):
    """Raised when declared source flags disagree with recomputed axis evidence."""


def parse_bool(value: Any, *, field: str) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no"}:
        return False
    raise M2GateError(f"Invalid boolean for {field}: {value!r}")


def finite_float(value: Any, *, field: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise M2GateError(f"Invalid numeric value for {field}: {value!r}") from error
    if not math.isfinite(parsed):
        raise M2GateError(f"Non-finite numeric value for {field}: {value!r}")
    return parsed


def parse_cluster_sizes(
    value: Any, *, max_groups: int | None = MAX_METHYL_GROUPS
) -> dict[str, int]:
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError as error:
            raise M2GateError(f"Invalid cluster_sizes JSON: {value!r}") from error
    if not isinstance(value, Mapping) or not value:
        raise M2GateError("cluster_sizes must be a non-empty object")
    parsed: dict[str, int] = {}
    for label, count_value in value.items():
        label_text = str(label).strip()
        try:
            count = int(count_value)
        except (TypeError, ValueError) as error:
            raise M2GateError(f"Invalid cluster size for {label_text!r}: {count_value!r}") from error
        if not label_text or count <= 0 or str(count_value).strip() != str(count):
            raise M2GateError(f"Invalid cluster size for {label_text!r}: {count_value!r}")
        parsed[label_text] = count
    if len(parsed) < 2 or (max_groups is not None and len(parsed) > max_groups):
        supported = f"2-{max_groups}" if max_groups is not None else ">=2"
        raise M2GateError(
            f"M2 requires {supported} observed methyl groups; found {len(parsed)}"
        )
    return parsed


def parse_categorical_level_counts(value: Any) -> dict[str, int]:
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError as error:
            raise M2GateError(f"Invalid categorical-level-count JSON: {value!r}") from error
    if not isinstance(value, Mapping):
        raise M2GateError("categorical level counts must be an object")
    missing = sorted(set(CATEGORICAL_LEVEL_KEYS).difference(value))
    extra = sorted(set(value).difference(CATEGORICAL_LEVEL_KEYS))
    if missing or extra:
        raise M2GateError(
            "Categorical-level-count keys drift: "
            f"missing={missing} extra={extra}"
        )
    parsed: dict[str, int] = {}
    for prefix in CATEGORICAL_LEVEL_KEYS:
        raw = value[prefix]
        try:
            count = int(raw)
        except (TypeError, ValueError) as error:
            raise M2GateError(
                f"Invalid observed categorical level count for {prefix}: {raw!r}"
            ) from error
        if count < 1 or str(raw).strip() != str(count):
            raise M2GateError(
                f"Invalid observed categorical level count for {prefix}: {raw!r}"
            )
        ceiling = CATEGORICAL_LEVEL_CEILINGS[prefix]
        if count > ceiling:
            raise M2GateError(
                f"Observed categorical levels exceed declared ceiling for {prefix}: "
                f"observed={count} ceiling={ceiling}"
            )
        parsed[prefix] = count
    return parsed


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
    degrees_between = groups - 1
    degrees_within = n - groups
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        critical = f.ppf(1.0 - POWER_ALPHA, degrees_between, degrees_within)
        noncentrality = n * threshold / (1.0 - threshold)
        return float(ncf.sf(critical, degrees_between, degrees_within, noncentrality))


@lru_cache(maxsize=None)
def minimum_n_for_target_power(
    prefix: str,
    kind: str,
    groups: int,
    threshold: float,
) -> int:
    if not 2 <= groups <= MAX_METHYL_GROUPS:
        raise M2GateError(f"Invalid methyl-group count for power calculation: {groups}")
    start = max(groups + 1, groups * MIN_GROUP_N)
    for n in range(start, 10_001):
        if kind == "categorical":
            levels = CATEGORICAL_LEVEL_CEILINGS[prefix]
            power = _categorical_power(n, groups, levels, threshold)
        elif kind == "continuous":
            power = _continuous_power(n, groups, threshold)
        else:
            raise M2GateError(f"Unknown M2 axis kind: {kind!r}")
        if math.isfinite(power) and power >= TARGET_POWER:
            return n
    raise M2GateError(
        f"Unable to attain target power for {prefix}: groups={groups} threshold={threshold}"
    )


def classify_axis(
    row: Mapping[str, Any],
    *,
    prefix: str,
    kind: str,
    effect_suffix: str,
    threshold: float,
    cluster_sizes: Mapping[str, int],
    categorical_level_counts: Mapping[str, int],
) -> dict[str, Any]:
    effect_field = f"{prefix}_{effect_suffix}"
    p_field = f"{prefix}_p_perm"
    n_field = f"{prefix}_n"
    aligned_field = f"{prefix}_aligned"
    missing = [
        field
        for field in (effect_field, p_field, n_field, aligned_field)
        if field not in row
    ]
    if missing:
        raise M2GateError(f"Missing M2 axis fields for {prefix}: {','.join(missing)}")
    try:
        n = int(row[n_field])
    except (TypeError, ValueError) as error:
        raise M2GateError(f"Invalid count for {n_field}: {row[n_field]!r}") from error
    if n < 0:
        raise M2GateError(f"Negative count for {n_field}: {n}")
    cluster_total = sum(cluster_sizes.values())
    if n != cluster_total:
        raise M2GateError(
            f"M2 axis count/core-cluster drift for {prefix}: n={n} cluster_total={cluster_total}"
        )
    groups = len(cluster_sizes)
    min_group_n = min(cluster_sizes.values())
    observed_levels = categorical_level_counts[prefix] if kind == "categorical" else None
    planning_levels = CATEGORICAL_LEVEL_CEILINGS[prefix] if kind == "categorical" else None
    declared_aligned = parse_bool(row[aligned_field], field=aligned_field)
    effect_missing = row[effect_field] in {None, ""}
    p_missing = row[p_field] in {None, ""}
    if effect_missing != p_missing:
        raise M2GateError(
            f"Partial M2 axis evidence for {prefix}: {effect_field}={row[effect_field]!r}, "
            f"{p_field}={row[p_field]!r}"
        )
    if effect_missing:
        if declared_aligned:
            raise M2GateError(f"Null axis evidence cannot be aligned for {prefix}")
        if kind != "categorical":
            return {
                "kind": kind,
                "effect_metric": effect_suffix,
                "effect": None,
                "effect_threshold": threshold,
                "p_perm": None,
                "p_threshold_exclusive": P_MAX,
                "permutations": PERMUTATIONS,
                "permutation_p_floor": PERMUTATION_P_FLOOR,
                "n": n,
                "n_groups": groups,
                "min_group_n": min_group_n,
                "declared_aligned": False,
                "recomputed_aligned": False,
                "power_target": TARGET_POWER,
                "power_at_effect_threshold": None,
                "minimum_n_for_target_power": None,
                "observed_category_levels": None,
                "planning_category_levels": None,
                "adequate_information_for_non_alignment": False,
                "positive_alignment_overrides_power": False,
                "status": MISSING_STATISTIC_AXIS_STATUS,
            }
        if observed_levels != 1:
            raise M2GateError(
                f"Missing categorical statistic is not explained by a constant axis for {prefix}: "
                f"observed_levels={observed_levels}"
            )
        return {
            "kind": kind,
            "effect_metric": effect_suffix,
            "effect": None,
            "effect_threshold": threshold,
            "p_perm": None,
            "p_threshold_exclusive": P_MAX,
            "permutations": PERMUTATIONS,
            "permutation_p_floor": PERMUTATION_P_FLOOR,
            "n": n,
            "n_groups": groups,
            "min_group_n": min_group_n,
            "declared_aligned": False,
            "recomputed_aligned": False,
            "power_target": TARGET_POWER,
            "power_at_effect_threshold": None,
            "minimum_n_for_target_power": None,
            "observed_category_levels": observed_levels,
            "planning_category_levels": planning_levels,
            "adequate_information_for_non_alignment": True,
            "positive_alignment_overrides_power": False,
            "status": CONSTANT_AXIS_STATUS,
        }
    if kind == "categorical" and observed_levels is not None and observed_levels < 2:
        raise M2GateError(
            f"Categorical statistic exists for a constant observed axis {prefix}: "
            f"observed_levels={observed_levels}"
        )
    effect = finite_float(row[effect_field], field=effect_field)
    p_value = finite_float(row[p_field], field=p_field)
    if not 0.0 <= effect <= 1.0:
        raise M2GateError(f"Effect outside [0,1] for {effect_field}: {effect}")
    if not 0.0 <= p_value <= 1.0:
        raise M2GateError(f"P value outside [0,1] for {p_field}: {p_value}")
    if p_value < PERMUTATION_P_FLOOR - 1e-12:
        raise M2GateError(
            f"Permutation p value below add-one floor for {p_field}: {p_value}"
        )
    grid_value = p_value * (PERMUTATIONS + 1)
    if not math.isclose(grid_value, round(grid_value), abs_tol=1e-9):
        raise M2GateError(
            f"Permutation p value is not on the {PERMUTATIONS}-permutation add-one grid "
            f"for {p_field}: {p_value}"
        )
    minimum_n = minimum_n_for_target_power(prefix, kind, groups, threshold)
    power = (
        _categorical_power(
            n,
            groups,
            CATEGORICAL_LEVEL_CEILINGS[prefix],
            threshold,
        )
        if kind == "categorical"
        else _continuous_power(n, groups, threshold)
    )
    sufficiently_informative = (
        min_group_n >= MIN_GROUP_N and n >= minimum_n and power >= TARGET_POWER
    )
    recomputed_aligned = effect >= threshold and p_value < P_MAX
    if declared_aligned != recomputed_aligned:
        raise M2GateError(
            f"Declared/recomputed alignment drift for {prefix}: "
            f"declared={declared_aligned} recomputed={recomputed_aligned}"
        )
    if recomputed_aligned:
        status = ALIGNED_AXIS_STATUS
    elif effect >= threshold:
        status = INDETERMINATE_AXIS_STATUS
    elif not sufficiently_informative:
        status = LOW_POWER_AXIS_STATUS
    else:
        status = LOW_EFFECT_AXIS_STATUS
    return {
        "kind": kind,
        "effect_metric": effect_suffix,
        "effect": effect,
        "effect_threshold": threshold,
        "p_perm": p_value,
        "p_threshold_exclusive": P_MAX,
        "permutations": PERMUTATIONS,
        "permutation_p_floor": PERMUTATION_P_FLOOR,
        "n": n,
        "n_groups": groups,
        "min_group_n": min_group_n,
        "declared_aligned": declared_aligned,
        "recomputed_aligned": recomputed_aligned,
        "power_target": TARGET_POWER,
        "power_at_effect_threshold": power,
        "minimum_n_for_target_power": minimum_n,
        "observed_category_levels": observed_levels,
        "planning_category_levels": planning_levels,
        "adequate_information_for_non_alignment": sufficiently_informative,
        "positive_alignment_overrides_power": (
            recomputed_aligned and not sufficiently_informative
        ),
        "status": status,
    }


def evaluate_m2_screen(
    row: Mapping[str, Any],
    *,
    categorical_level_counts: Mapping[str, Any] | str | None = None,
) -> dict[str, Any]:
    missing_high_level = [
        field
        for field in HIGH_LEVEL_FIELDS
        if row.get(field) is None or row.get(field) == ""
    ]
    missing_axis = [field for field in AXIS_SOURCE_FIELDS if field not in row]
    if missing_high_level or missing_axis:
        missing = sorted(set(missing_high_level + missing_axis))
        return {
            "contract": GATE_CONTRACT,
            "eligible": False,
            "evaluable": False,
            "status": "NOT_EVALUABLE_M2_SCREEN_FIELDS_MISSING:" + ",".join(missing),
            "axis_statuses": {},
            "indeterminate_axes": [],
            "low_power_axes": [],
            "aligned_axes": [],
            "constant_axes": [],
            "aligned_below_negative_evaluability_power_axes": [],
        }
    cluster_sizes = parse_cluster_sizes(row["cluster_sizes"], max_groups=None)
    if categorical_level_counts is None:
        categorical_level_counts = row.get("m2_categorical_level_counts")
    observed_categorical_levels = (
        None
        if categorical_level_counts is None or categorical_level_counts == ""
        else parse_categorical_level_counts(categorical_level_counts)
    )
    if len(cluster_sizes) > MAX_METHYL_GROUPS:
        return {
            "contract": GATE_CONTRACT,
            "eligible": False,
            "evaluable": False,
            "status": (
                "NOT_EVALUABLE_M2_GROUP_COUNT_EXCEEDS_PLANNING_MODEL_MAXIMUM:"
                f"observed={len(cluster_sizes)}:maximum={MAX_METHYL_GROUPS}"
            ),
            "axis_statuses": {},
            "indeterminate_axes": [],
            "low_power_axes": [],
            "aligned_axes": [],
            "constant_axes": [],
            "aligned_below_negative_evaluability_power_axes": [],
            "categorical_level_counts": observed_categorical_levels,
        }
    if observed_categorical_levels is None:
        return {
            "contract": GATE_CONTRACT,
            "eligible": False,
            "evaluable": False,
            "status": "NOT_EVALUABLE_M2_CATEGORICAL_LEVEL_COUNTS_MISSING",
            "axis_statuses": {},
            "indeterminate_axes": [],
            "low_power_axes": [],
            "aligned_axes": [],
            "constant_axes": [],
            "aligned_below_negative_evaluability_power_axes": [],
        }
    axes = {
        prefix: classify_axis(
            row,
            prefix=prefix,
            kind=kind,
            effect_suffix=effect_suffix,
            threshold=threshold,
            cluster_sizes=cluster_sizes,
            categorical_level_counts=observed_categorical_levels,
        )
        for prefix, kind, effect_suffix, threshold in AXIS_SPECS
    }
    aligned_axes = sorted(
        prefix
        for prefix, evidence in axes.items()
        if evidence["status"] == ALIGNED_AXIS_STATUS
    )
    indeterminate_axes = sorted(
        prefix
        for prefix, evidence in axes.items()
        if evidence["status"] in {
            INDETERMINATE_AXIS_STATUS,
            LOW_POWER_AXIS_STATUS,
            MISSING_STATISTIC_AXIS_STATUS,
        }
    )
    low_power_axes = sorted(
        prefix
        for prefix, evidence in axes.items()
        if evidence["status"] == LOW_POWER_AXIS_STATUS
    )
    constant_axes = sorted(
        prefix
        for prefix, evidence in axes.items()
        if evidence["status"] == CONSTANT_AXIS_STATUS
    )
    aligned_below_power_axes = sorted(
        prefix
        for prefix, evidence in axes.items()
        if evidence["status"] == ALIGNED_AXIS_STATUS
        and evidence["positive_alignment_overrides_power"]
    )
    stable = parse_bool(row["stable_null_multigroup"], field="stable_null_multigroup")
    ari = finite_float(row["modal_assignment_ari_min"], field="modal_assignment_ari_min")
    hp_confound = parse_bool(row["hp_axis_confound"], field="hp_axis_confound")
    technical_confound = parse_bool(
        row["technical_axis_confound"], field="technical_axis_confound"
    )
    residual = parse_bool(
        row["residual_unexplained_multigroup"], field="residual_unexplained_multigroup"
    )
    phase_robust = parse_bool(
        row["phase_anchored_robust_epigenetic_candidate"],
        field="phase_anchored_robust_epigenetic_candidate",
    )
    expected_hp_confound = any(axis in aligned_axes for axis in ("hp_exact", "hp_family"))
    expected_technical_confound = any(
        axis in aligned_axes
        for axis in ("strand", "start", "end", "length", "mapq", "cpg_called")
    )
    if hp_confound != expected_hp_confound:
        raise M2GateError(
            f"hp_axis_confound drift: declared={hp_confound} expected={expected_hp_confound}"
        )
    if technical_confound != expected_technical_confound:
        raise M2GateError(
            "technical_axis_confound drift: "
            f"declared={technical_confound} expected={expected_technical_confound}"
        )
    expected_residual = stable and not hp_confound and not technical_confound
    if residual != expected_residual:
        raise M2GateError(
            f"residual_unexplained_multigroup drift: declared={residual} "
            f"expected={expected_residual}"
        )
    if phase_robust and not residual:
        raise M2GateError("Phase-anchored robust flag requires residual measured-axis structure")

    base = {
        "contract": GATE_CONTRACT,
        "axis_statuses": axes,
        "indeterminate_axes": indeterminate_axes,
        "low_power_axes": low_power_axes,
        "aligned_axes": aligned_axes,
        "constant_axes": constant_axes,
        "categorical_level_counts": observed_categorical_levels,
        "aligned_below_negative_evaluability_power_axes": aligned_below_power_axes,
    }
    if not stable:
        return {
            **base,
            "eligible": False,
            "evaluable": False,
            "status": "NOT_RUN_M2_WITHOUT_M1_STABLE_MULTIGROUP",
        }
    if ari < 0.8:
        raise M2GateError(f"Stable M1 row has modal ARI below 0.8: {ari}")
    if indeterminate_axes:
        return {
            **base,
            "eligible": False,
            "evaluable": False,
            "status": (
                "NOT_EVALUABLE_M2_AXIS_INDETERMINATE:"
                + ",".join(
                    f"{axis}={axes[axis]['status']}" for axis in indeterminate_axes
                )
            ),
        }
    if hp_confound:
        return {
            **base,
            "eligible": False,
            "evaluable": True,
            "status": "INELIGIBLE_SCREEN_HP_CONFOUND",
        }
    if technical_confound:
        return {
            **base,
            "eligible": False,
            "evaluable": True,
            "status": "INELIGIBLE_SCREEN_TECHNICAL_CONFOUND",
        }
    if not residual:
        return {
            **base,
            "eligible": False,
            "evaluable": True,
            "status": "INELIGIBLE_NOT_RESIDUAL_UNEXPLAINED",
        }
    if not phase_robust:
        return {
            **base,
            "eligible": False,
            "evaluable": True,
            "status": "INELIGIBLE_NOT_M2_PHASE_ANCHORED_ROBUST",
        }
    return {
        **base,
        "eligible": True,
        "evaluable": True,
        "status": "ELIGIBLE_M2_RESIDUAL_UNEXPLAINED_AND_AXES_DETERMINATE",
    }
