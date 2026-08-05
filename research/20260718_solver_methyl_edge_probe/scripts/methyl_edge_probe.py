#!/usr/bin/env python3
"""Synthetic methylation edge tie-breaker with explicit abstention.

The score is deliberately separate from the sSNV state likelihood.  It can only
state which parent edge is more compatible with methylation similarity under the
declared model; it does not prove ancestry.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import random
import statistics
from dataclasses import asdict, dataclass
from typing import Mapping, Sequence


@dataclass(frozen=True)
class MethylRead:
    read_id: str
    stratum: str
    probabilities: Mapping[int, float]


@dataclass
class EdgePreference:
    status: str
    delta_m: float | None
    distance_10_11: float | None
    distance_01_11: float | None
    bootstrap_ci95: tuple[float, float] | None
    bootstrap_valid_replicates: int
    bootstrap_direction_fraction: float | None
    n_reads_10: int
    n_reads_01: int
    n_reads_11: int
    stratum: str | None
    min_common_cpg: int
    delta_min: float
    interpretation: str

    def to_dict(self) -> dict:
        return asdict(self)


def confidence_weight(probability: float) -> float:
    if not 0.0 <= probability <= 1.0:
        raise ValueError("methylation probability must lie in [0,1]")
    return 2.0 * abs(float(probability) - 0.5)


def expected_bernoulli_mismatch(left: float, right: float) -> float:
    return float(left) + float(right) - 2.0 * float(left) * float(right)


def read_distance(
    left: MethylRead,
    right: MethylRead,
    *,
    min_common_cpg: int,
) -> float | None:
    common = sorted(set(left.probabilities) & set(right.probabilities))
    if len(common) < min_common_cpg:
        return None
    numerator = 0.0
    denominator = 0.0
    for cpg in common:
        p_left = float(left.probabilities[cpg])
        p_right = float(right.probabilities[cpg])
        weight = confidence_weight(p_left) * confidence_weight(p_right)
        numerator += weight * expected_bernoulli_mismatch(p_left, p_right)
        denominator += weight
    if denominator <= 1e-15:
        return None
    return numerator / denominator


def state_distance(
    left: Sequence[MethylRead],
    right: Sequence[MethylRead],
    *,
    min_common_cpg: int,
) -> float | None:
    distances = [
        distance
        for left_read in left
        for right_read in right
        if (
            distance := read_distance(
                left_read,
                right_read,
                min_common_cpg=min_common_cpg,
            )
        )
        is not None
    ]
    return None if not distances else float(statistics.median(distances))


def edge_delta(
    reads_10: Sequence[MethylRead],
    reads_01: Sequence[MethylRead],
    reads_11: Sequence[MethylRead],
    *,
    min_common_cpg: int,
) -> tuple[float, float, float] | None:
    distance_10_11 = state_distance(
        reads_10,
        reads_11,
        min_common_cpg=min_common_cpg,
    )
    distance_01_11 = state_distance(
        reads_01,
        reads_11,
        min_common_cpg=min_common_cpg,
    )
    if distance_10_11 is None or distance_01_11 is None:
        return None
    return (
        float(distance_01_11 - distance_10_11),
        float(distance_10_11),
        float(distance_01_11),
    )


def quantile(values: Sequence[float], probability: float) -> float:
    if not values:
        raise ValueError("quantile requires values")
    if not 0.0 <= probability <= 1.0:
        raise ValueError("invalid quantile probability")
    ordered = sorted(float(value) for value in values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def bootstrap_delta(
    reads_10: Sequence[MethylRead],
    reads_01: Sequence[MethylRead],
    reads_11: Sequence[MethylRead],
    *,
    min_common_cpg: int,
    replicates: int,
    seed: int,
) -> list[float]:
    rng = random.Random(seed)
    result: list[float] = []
    for _ in range(replicates):
        sampled_10 = [rng.choice(reads_10) for _ in reads_10]
        sampled_01 = [rng.choice(reads_01) for _ in reads_01]
        sampled_11 = [rng.choice(reads_11) for _ in reads_11]
        current = edge_delta(
            sampled_10,
            sampled_01,
            sampled_11,
            min_common_cpg=min_common_cpg,
        )
        if current is not None:
            result.append(current[0])
    return result


def evaluate_edge_preference(
    reads_10: Sequence[MethylRead],
    reads_01: Sequence[MethylRead],
    reads_11: Sequence[MethylRead],
    *,
    min_reads_per_state: int = 20,
    min_common_cpg: int = 10,
    delta_min: float = 0.05,
    bootstrap_replicates: int = 500,
    seed: int = 20260718,
) -> EdgePreference:
    sizes = (len(reads_10), len(reads_01), len(reads_11))
    if min(sizes, default=0) < min_reads_per_state:
        return EdgePreference(
            "ABSTAIN_INSUFFICIENT_STATE_READS",
            None,
            None,
            None,
            None,
            0,
            None,
            *sizes,
            None,
            min_common_cpg,
            delta_min,
            "At least one exact state lacks the preregistered molecule count.",
        )
    strata = {
        read.stratum
        for reads in (reads_10, reads_01, reads_11)
        for read in reads
    }
    if len(strata) != 1:
        return EdgePreference(
            "ABSTAIN_MIXED_HP_PS_STRATA",
            None,
            None,
            None,
            None,
            0,
            None,
            *sizes,
            None,
            min_common_cpg,
            delta_min,
            "State profiles do not belong to one HP-family × PS stratum.",
        )
    point = edge_delta(
        reads_10,
        reads_01,
        reads_11,
        min_common_cpg=min_common_cpg,
    )
    stratum = next(iter(strata))
    if point is None:
        return EdgePreference(
            "ABSTAIN_INSUFFICIENT_COMMON_CPG_OR_CONFIDENCE",
            None,
            None,
            None,
            None,
            0,
            None,
            *sizes,
            stratum,
            min_common_cpg,
            delta_min,
            "Pairwise profiles lack enough shared informative CpGs.",
        )
    delta_m, distance_10_11, distance_01_11 = point
    bootstrap = bootstrap_delta(
        reads_10,
        reads_01,
        reads_11,
        min_common_cpg=min_common_cpg,
        replicates=bootstrap_replicates,
        seed=seed,
    )
    if len(bootstrap) < max(1, math.ceil(0.9 * bootstrap_replicates)):
        return EdgePreference(
            "ABSTAIN_BOOTSTRAP_INVALID",
            delta_m,
            distance_10_11,
            distance_01_11,
            None,
            len(bootstrap),
            None,
            *sizes,
            stratum,
            min_common_cpg,
            delta_min,
            "Too many resamples lacked a valid common-CpG distance.",
        )
    ci = (quantile(bootstrap, 0.025), quantile(bootstrap, 0.975))
    if delta_m > 0:
        direction_fraction = sum(value > 0 for value in bootstrap) / len(bootstrap)
    elif delta_m < 0:
        direction_fraction = sum(value < 0 for value in bootstrap) / len(bootstrap)
    else:
        direction_fraction = sum(value == 0 for value in bootstrap) / len(bootstrap)
    if (
        abs(delta_m) < delta_min
        or ci[0] <= 0.0 <= ci[1]
        or direction_fraction < 0.9
    ):
        status = "TIE"
        interpretation = "Methylation evidence does not stably separate the two parent edges."
    elif delta_m > 0:
        status = "MODEL_FAVORED_10_TO_11"
        interpretation = (
            "11 is methylation-closer to 10 than to 01 under the declared model; "
            "this is not ancestry proof."
        )
    else:
        status = "MODEL_FAVORED_01_TO_11"
        interpretation = (
            "11 is methylation-closer to 01 than to 10 under the declared model; "
            "this is not ancestry proof."
        )
    return EdgePreference(
        status,
        delta_m,
        distance_10_11,
        distance_01_11,
        ci,
        len(bootstrap),
        direction_fraction,
        *sizes,
        stratum,
        min_common_cpg,
        delta_min,
        interpretation,
    )


def make_reads(
    state: str,
    profile: Sequence[float],
    *,
    n_reads: int = 12,
    stratum: str = "HP1|PS100",
    cpg_offset: int = 1000,
) -> list[MethylRead]:
    probabilities = {
        cpg_offset + index: float(value)
        for index, value in enumerate(profile)
    }
    return [
        MethylRead(f"{state}_{index}", stratum, probabilities)
        for index in range(n_reads)
    ]


def synthetic_cases() -> list[dict]:
    low = [0.1] * 12
    high = [0.9] * 12
    alternating = [0.1 if index % 2 == 0 else 0.9 for index in range(12)]
    definitions = [
        ("positive_10", low, high, [0.15] * 12, "MODEL_FAVORED_10_TO_11"),
        ("positive_01", low, high, [0.85] * 12, "MODEL_FAVORED_01_TO_11"),
        ("exact_tie", low, high, alternating, "TIE"),
        ("same_profile_tie", high, high, high, "TIE"),
    ]
    rows: list[dict] = []
    for case_id, profile_10, profile_01, profile_11, expected in definitions:
        result = evaluate_edge_preference(
            make_reads("10", profile_10),
            make_reads("01", profile_01),
            make_reads("11", profile_11),
            min_reads_per_state=8,
            min_common_cpg=5,
            bootstrap_replicates=200,
        )
        row = {
            "case_id": case_id,
            "expected_status": expected,
            "result": result.to_dict(),
            "pass": result.status == expected,
        }
        rows.append(row)

    uninformative = evaluate_edge_preference(
        make_reads("10", [0.5] * 12),
        make_reads("01", [0.5] * 12),
        make_reads("11", [0.5] * 12),
        min_reads_per_state=8,
        min_common_cpg=5,
        bootstrap_replicates=50,
    )
    rows.append(
        {
            "case_id": "uninformative_probability",
            "expected_status": "ABSTAIN_INSUFFICIENT_COMMON_CPG_OR_CONFIDENCE",
            "result": uninformative.to_dict(),
            "pass": (
                uninformative.status
                == "ABSTAIN_INSUFFICIENT_COMMON_CPG_OR_CONFIDENCE"
            ),
        }
    )

    missing = evaluate_edge_preference(
        make_reads("10", low, cpg_offset=1000),
        make_reads("01", high, cpg_offset=1000),
        make_reads("11", low, cpg_offset=2000),
        min_reads_per_state=8,
        min_common_cpg=5,
        bootstrap_replicates=50,
    )
    rows.append(
        {
            "case_id": "missing_common_cpg",
            "expected_status": "ABSTAIN_INSUFFICIENT_COMMON_CPG_OR_CONFIDENCE",
            "result": missing.to_dict(),
            "pass": (
                missing.status
                == "ABSTAIN_INSUFFICIENT_COMMON_CPG_OR_CONFIDENCE"
            ),
        }
    )

    mixed = evaluate_edge_preference(
        make_reads("10", low, stratum="HP1|PS100"),
        make_reads("01", high, stratum="HP2|PS100"),
        make_reads("11", low, stratum="HP1|PS100"),
        min_reads_per_state=8,
        min_common_cpg=5,
        bootstrap_replicates=50,
    )
    rows.append(
        {
            "case_id": "mixed_hp_ps_confound",
            "expected_status": "ABSTAIN_MIXED_HP_PS_STRATA",
            "result": mixed.to_dict(),
            "pass": mixed.status == "ABSTAIN_MIXED_HP_PS_STRATA",
        }
    )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = synthetic_cases()
    receipt = {
        "schema": "intersubmod.solver_methyl_edge_probe.methyl.v1",
        "scope": "SYNTHETIC_ONLY",
        "formula": {
            "read_distance": (
                "weighted mean of p_i + p_j - 2*p_i*p_j over common CpGs; "
                "weight=(2|p_i-0.5|)(2|p_j-0.5|)"
            ),
            "state_distance": "median pairwise read distance",
            "delta_m": "D_M(01,11)-D_M(10,11)",
        },
        "cases": rows,
        "all_pass": all(row["pass"] for row in rows),
        "claim_ceiling": (
            "Synthetic recovery/abstention under the declared methyl-similarity "
            "model; not true ancestry, clone/subclone, or edge accuracy."
        ),
        "real_data_floor": {
            "exact_reads_per_state_min": 20,
            "preferred_reads_per_state": 30,
            "common_cpg_min": 10,
            "same_hp_family_and_ps": True,
            "cn_loh_gate_required": True,
            "held_out_or_cross_fit_required": True,
        },
    }
    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "all_pass": receipt["all_pass"],
                "cases": [
                    {
                        "case_id": row["case_id"],
                        "status": row["result"]["status"],
                        "delta_m": row["result"]["delta_m"],
                    }
                    for row in rows
                ],
            },
            ensure_ascii=False,
        )
    )
    return 0 if receipt["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
