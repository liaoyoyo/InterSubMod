#!/usr/bin/env python3
"""Build a fixed-denominator red-team cause-decomposition audit.

The analysis joins the frozen 5,720-region topology population to exact shared
alleles, paired caller VAF/depth records, the exact-coordinate coarse topology
comparison, and the currently available shared SAVANA CN context.  All VAF and
depth associations are diagnostic because they are derived from the same reads
used by the topology workflow; they are not orthogonal validation evidence.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm


SEED = 20260712
BOOTSTRAP_REPLICATES = 500
STRICT_OUTCOMES = {
    "strict_full_exact",
    "A_induced_substructure",
    "B_induced_substructure",
}


def parse_args() -> argparse.Namespace:
    script = Path(__file__).resolve()
    default_repo = script.parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=default_repo)
    parser.add_argument("--output-dir", type=Path, default=script.parent)
    parser.add_argument("--bootstrap-replicates", type=int, default=BOOTSTRAP_REPLICATES)
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args()


def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().map({"true": True, "false": False}).fillna(False)


def outcome_group(value: Any) -> str:
    text = str(value)
    if text == "not_evaluable":
        return "not_evaluable"
    if text == "strict_full_exact":
        return "strict_full_exact"
    if text in {"A_induced_substructure", "B_induced_substructure"}:
        return "true_induced_substructure"
    if text.startswith("resolution_"):
        return "resolution_difference"
    if text == "shared_core_only":
        return "shared_core_only"
    if text == "candidate_overlap":
        return "candidate_overlap"
    if text == "conflict":
        return "conflict"
    return text


def multiplicity_bin(value: Any, vaf: bool = False) -> str:
    if pd.isna(value):
        return "unavailable"
    number = float(value)
    if number <= 1:
        return "1"
    if vaf:
        return ">=2"
    if number <= 5:
        return "2-5"
    if number <= 20:
        return "6-20"
    return ">20"


def exact_allele_key(value: str) -> Tuple[str, int, str, str]:
    records = json.loads(value)
    if len(records) != 1 or len(records[0]) != 4:
        raise ValueError(f"Expected one exact allele, observed: {value}")
    chrom, pos, ref, alt = records[0]
    return str(chrom), int(pos), str(ref), str(alt)


def mode_or_nan(series: pd.Series) -> Any:
    clean = series.dropna()
    if clean.empty:
        return np.nan
    modes = clean.mode(dropna=True)
    return modes.iloc[0] if not modes.empty else clean.iloc[0]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def write_tsv_gz(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )


def add_raw_region_metrics(
    regions: pd.DataFrame,
    alleles: pd.DataFrame,
    raw: pd.DataFrame,
    source: str,
    slug: str,
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    source_rows = raw.loc[raw["source"] == source].copy()
    if source_rows["key"].duplicated().any():
        raise ValueError(f"Duplicate exact allele keys in raw source {source}")

    required = [
        "key",
        "af_hcc1395",
        "dp_hcc1395",
        "alt_count_hcc1395",
        "af_dorado",
        "dp_dorado",
        "alt_count_dorado",
    ]
    joined = alleles.merge(source_rows[required], on="key", how="left", validate="one_to_one")
    joined["raw_hit"] = joined["dp_hcc1395"].notna() & joined["dp_dorado"].notna()
    joined["pair_min_dp"] = joined[["dp_hcc1395", "dp_dorado"]].min(axis=1)
    joined["abs_dp_delta"] = (joined["dp_dorado"] - joined["dp_hcc1395"]).abs()
    joined["abs_log2_dp_ratio"] = np.abs(np.log2(joined["dp_dorado"] / joined["dp_hcc1395"]))
    joined["vaf_delta"] = joined["af_dorado"] - joined["af_hcc1395"]
    joined["abs_vaf_delta"] = joined["vaf_delta"].abs()
    pooled = (
        joined["alt_count_hcc1395"] + joined["alt_count_dorado"]
    ) / (joined["dp_hcc1395"] + joined["dp_dorado"])
    joined["expected_binomial_variance"] = pooled * (1.0 - pooled) * (
        1.0 / joined["dp_hcc1395"] + 1.0 / joined["dp_dorado"]
    )
    finite = (
        joined["raw_hit"]
        & np.isfinite(joined["expected_binomial_variance"])
        & (joined["expected_binomial_variance"] > 0)
    )
    joined["binomial_compatible"] = np.nan
    joined.loc[finite, "binomial_compatible"] = (
        joined.loc[finite, "abs_vaf_delta"]
        <= 1.96 * np.sqrt(joined.loc[finite, "expected_binomial_variance"])
    ).astype(float)
    joined["observed_delta_sq"] = joined["vaf_delta"] ** 2

    grouped = joined.groupby("match_id", sort=False)
    metrics = grouped.agg(
        n_expected_sites=("key", "size"),
        n_sites=("raw_hit", "sum"),
        mean_dp_hcc=("dp_hcc1395", "mean"),
        min_dp_hcc=("dp_hcc1395", "min"),
        mean_dp_dorado=("dp_dorado", "mean"),
        min_dp_dorado=("dp_dorado", "min"),
        mean_pair_min_dp=("pair_min_dp", "mean"),
        min_pair_dp=("pair_min_dp", "min"),
        mean_abs_dp_delta=("abs_dp_delta", "mean"),
        mean_abs_log2_dp_ratio=("abs_log2_dp_ratio", "mean"),
        mean_abs_vaf_delta=("abs_vaf_delta", "mean"),
        mean_signed_vaf_delta=("vaf_delta", "mean"),
        binomial_testable_n=("binomial_compatible", "count"),
        binomial_compatible_fraction=("binomial_compatible", "mean"),
        observed_delta_sq_sum=("observed_delta_sq", "sum"),
        expected_binomial_variance_sum=("expected_binomial_variance", "sum"),
    ).reset_index()
    metrics["coverage_fraction"] = metrics["n_sites"] / metrics["n_expected_sites"]
    metrics["complete"] = metrics["n_sites"] == metrics["n_expected_sites"]
    metrics["binomial_noise_rmse_ratio"] = np.sqrt(
        metrics["observed_delta_sq_sum"] / metrics["expected_binomial_variance_sum"]
    )
    metrics["binomial_all_sites_compatible"] = (
        (metrics["binomial_testable_n"] > 0)
        & (metrics["binomial_compatible_fraction"] == 1.0)
    )
    metrics = metrics.drop(columns=["observed_delta_sq_sum", "expected_binomial_variance_sum"])
    metrics = metrics.rename(
        columns={column: f"raw_{slug}_{column}" for column in metrics.columns if column != "match_id"}
    )
    regions = regions.merge(metrics, on="match_id", how="left", validate="one_to_one")

    n_hits = int(joined["raw_hit"].sum())
    complete_regions = int(metrics[f"raw_{slug}_complete"].sum())
    partial_regions = int(
        (
            (metrics[f"raw_{slug}_n_sites"] > 0)
            & (~metrics[f"raw_{slug}_complete"])
        ).sum()
    )
    zero_regions = int((metrics[f"raw_{slug}_n_sites"] == 0).sum())
    return regions, {
        "source": source,
        "allele_hits": n_hits,
        "allele_expected": int(len(joined)),
        "complete_regions": complete_regions,
        "partial_regions": partial_regions,
        "zero_regions": zero_regions,
    }


def add_cn_context(
    regions: pd.DataFrame,
    alleles: pd.DataFrame,
    metadata: pd.DataFrame,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    metadata = metadata.copy()
    metadata["key"] = (
        metadata["chrom"].astype(str)
        + ":"
        + metadata["pos"].astype(int).astype(str)
        + ":"
        + metadata["ref"].astype(str)
        + ":"
        + metadata["alt"].astype(str)
    )
    if metadata["key"].duplicated().any():
        raise ValueError("Duplicate exact allele keys in PyClone site metadata")
    fields = [
        "key",
        "segment_id",
        "total_cn_raw",
        "minor_cn_raw",
        "total_cn_discrete",
        "major_cn",
        "minor_cn",
        "near_integer",
    ]
    joined = alleles.merge(metadata[fields], on="key", how="left", validate="one_to_one")
    joined["cn_hit"] = joined["segment_id"].notna()
    joined["loh"] = np.where(joined["cn_hit"], joined["minor_cn"] == 0, np.nan)

    def summarize(group: pd.DataFrame) -> pd.Series:
        hit = group.loc[group["cn_hit"]].copy()
        n_expected = int(len(group))
        n_sites = int(len(hit))
        coverage = n_sites / n_expected if n_expected else np.nan
        return pd.Series(
            {
                "cn_n_expected_sites": n_expected,
                "cn_n_sites": n_sites,
                "cn_coverage_fraction": coverage,
                "cn_complete": n_sites == n_expected,
                "cn_segment_count": int(hit["segment_id"].nunique()) if n_sites else 0,
                "cn_total_cn_mode": mode_or_nan(hit["total_cn_discrete"]),
                "cn_major_cn_mode": mode_or_nan(hit["major_cn"]),
                "cn_minor_cn_mode": mode_or_nan(hit["minor_cn"]),
                "cn_loh_fraction": float(hit["loh"].mean()) if n_sites else np.nan,
                "cn_near_integer_fraction": float(hit["near_integer"].mean()) if n_sites else np.nan,
            }
        )

    cn = joined.groupby("match_id", sort=False).apply(summarize).reset_index()
    cn["cn_status"] = "partial_coverage"
    cn.loc[cn["cn_n_sites"] == 0, "cn_status"] = "unavailable"
    full = cn["cn_complete"]
    cn.loc[full & (cn["cn_loh_fraction"] == 0), "cn_status"] = "full_non_loh"
    cn.loc[full & (cn["cn_loh_fraction"] == 1), "cn_status"] = "full_all_loh"
    cn.loc[full & cn["cn_loh_fraction"].between(0, 1, inclusive="neither"), "cn_status"] = "full_mixed_loh"
    cn["cn_loh_full_factor"] = "unavailable_or_partial"
    cn.loc[full & (cn["cn_loh_fraction"] == 0), "cn_loh_full_factor"] = "no_loh"
    cn.loc[full & (cn["cn_loh_fraction"] > 0), "cn_loh_full_factor"] = "any_loh"
    regions = regions.merge(cn, on="match_id", how="left", validate="one_to_one")
    return regions, {
        "allele_hits": int(joined["cn_hit"].sum()),
        "allele_expected": int(len(joined)),
        "complete_regions": int(cn["cn_complete"].sum()),
        "partial_regions": int(((cn["cn_n_sites"] > 0) & (~cn["cn_complete"])).sum()),
        "zero_regions": int((cn["cn_n_sites"] == 0).sum()),
    }


def quartile_labels(series: pd.Series, prefix: str) -> Tuple[pd.Series, List[float]]:
    labels = [f"Q1_{prefix}_low", f"Q2_{prefix}", f"Q3_{prefix}", f"Q4_{prefix}_high"]
    values, bins = pd.qcut(series, q=4, labels=labels, retbins=True, duplicates="raise")
    return values.astype(str), [float(value) for value in bins]


def endpoint_definitions() -> Dict[str, Dict[str, str]]:
    return {
        "coarse_agree_fixed": {"y": "coarse_agree", "eligible": "all_eligible", "scope": "fixed_5720"},
        "read_evaluable_fixed": {"y": "read_evaluable", "eligible": "all_eligible", "scope": "fixed_5720"},
        "read_full_exact_fixed": {"y": "read_full_exact", "eligible": "all_eligible", "scope": "fixed_5720"},
        "read_exact_or_induced_fixed": {"y": "read_strict", "eligible": "all_eligible", "scope": "fixed_5720"},
        "read_full_exact_given_evaluable": {
            "y": "read_full_exact",
            "eligible": "read_evaluable",
            "scope": "read_evaluable_only",
        },
        "read_exact_or_induced_given_evaluable": {
            "y": "read_strict",
            "eligible": "read_evaluable",
            "scope": "read_evaluable_only",
        },
        "read_conflict_given_evaluable": {
            "y": "read_conflict",
            "eligible": "read_evaluable",
            "scope": "read_evaluable_only",
        },
        "vaf_evaluable_fixed": {"y": "vaf_evaluable", "eligible": "all_eligible", "scope": "fixed_5720"},
        "vaf_full_exact_fixed": {"y": "vaf_full_exact", "eligible": "all_eligible", "scope": "fixed_5720"},
        "vaf_exact_or_induced_fixed": {"y": "vaf_strict", "eligible": "all_eligible", "scope": "fixed_5720"},
        "vaf_full_exact_given_evaluable": {
            "y": "vaf_full_exact",
            "eligible": "vaf_evaluable",
            "scope": "vaf_evaluable_only",
        },
        "vaf_exact_or_induced_given_evaluable": {
            "y": "vaf_strict",
            "eligible": "vaf_evaluable",
            "scope": "vaf_evaluable_only",
        },
        "vaf_conflict_given_evaluable": {
            "y": "vaf_conflict",
            "eligible": "vaf_evaluable",
            "scope": "vaf_evaluable_only",
        },
    }


def build_factor_summary(
    regions: pd.DataFrame,
    factors: Mapping[str, str],
    endpoints: Mapping[str, Mapping[str, str]],
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for factor_name, column in factors.items():
        for level, level_rows in regions.groupby(column, dropna=False, sort=True):
            level_text = "NA" if pd.isna(level) else str(level)
            for endpoint, definition in endpoints.items():
                eligible = bool_series(level_rows[definition["eligible"]])
                denominator = int(eligible.sum())
                numerator = int(bool_series(level_rows.loc[eligible, definition["y"]]).sum())
                rows.append(
                    {
                        "factor": factor_name,
                        "factor_column": column,
                        "level": level_text,
                        "factor_level_n": int(len(level_rows)),
                        "endpoint": endpoint,
                        "endpoint_scope": definition["scope"],
                        "denominator": denominator,
                        "numerator": numerator,
                        "rate": numerator / denominator if denominator else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def build_numeric_outcome_summary(regions: pd.DataFrame) -> pd.DataFrame:
    features = [
        "raw_baseline_mean_pair_min_dp",
        "raw_baseline_min_pair_dp",
        "raw_baseline_mean_abs_dp_delta",
        "raw_baseline_mean_abs_log2_dp_ratio",
        "raw_baseline_mean_abs_vaf_delta",
        "raw_baseline_binomial_compatible_fraction",
        "raw_baseline_binomial_noise_rmse_ratio",
        "raw_latest_coverage_fraction",
        "read_candidate_topology_max",
        "vaf_candidate_topology_max",
        "cn_coverage_fraction",
        "cn_loh_fraction",
        "caller_extra_site_count",
    ]
    groupings = {
        "read": "read_outcome_group",
        "vaf": "vaf_outcome_group",
        "coarse": "coarse_outcome_group",
    }
    rows: List[Dict[str, Any]] = []
    for layer, group_column in groupings.items():
        for outcome, group in regions.groupby(group_column, sort=True):
            for feature in features:
                values = pd.to_numeric(group[feature], errors="coerce").dropna()
                rows.append(
                    {
                        "layer": layer,
                        "outcome_group": outcome,
                        "feature": feature,
                        "group_n": int(len(group)),
                        "nonmissing_n": int(len(values)),
                        "mean": float(values.mean()) if len(values) else np.nan,
                        "q25": float(values.quantile(0.25)) if len(values) else np.nan,
                        "median": float(values.median()) if len(values) else np.nan,
                        "q75": float(values.quantile(0.75)) if len(values) else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def odds_ratio_haldane(a: int, b: int, c: int, d: int) -> float:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def bh_adjust(p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite_indices = np.where(np.isfinite(values))[0]
    if not len(finite_indices):
        return adjusted
    order = finite_indices[np.argsort(values[finite_indices])]
    ranked = values[order]
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    adjusted[order] = np.minimum(q, 1.0)
    return adjusted


def build_contrast_effects(
    regions: pd.DataFrame,
    contrasts: Sequence[Mapping[str, Any]],
    endpoints: Mapping[str, Mapping[str, str]],
    bootstrap_replicates: int,
    seed: int,
) -> pd.DataFrame:
    blocks = sorted(regions["block_1mb"].unique())
    block_index = {block: index for index, block in enumerate(blocks)}
    rng = np.random.default_rng(seed)
    weights = rng.multinomial(
        len(blocks),
        np.repeat(1.0 / len(blocks), len(blocks)),
        size=bootstrap_replicates,
    )

    rows: List[Dict[str, Any]] = []
    bootstrap_inputs: List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
    for contrast in contrasts:
        endpoint_names = contrast.get("endpoints", list(endpoints))
        for endpoint in endpoint_names:
            definition = endpoints[endpoint]
            factor = regions[contrast["column"]].astype(str)
            eligible = bool_series(regions[definition["eligible"]])
            mask_a = eligible & (factor == contrast["level_a"])
            mask_b = eligible & (factor == contrast["level_b"])
            y = bool_series(regions[definition["y"]]).astype(int)
            denominator_a = int(mask_a.sum())
            denominator_b = int(mask_b.sum())
            numerator_a = int(y.loc[mask_a].sum())
            numerator_b = int(y.loc[mask_b].sum())
            rate_a = numerator_a / denominator_a if denominator_a else np.nan
            rate_b = numerator_b / denominator_b if denominator_b else np.nan

            table = [
                [numerator_a, denominator_a - numerator_a],
                [numerator_b, denominator_b - numerator_b],
            ]
            if denominator_a and denominator_b:
                fisher_p = float(stats.fisher_exact(table, alternative="two-sided").pvalue)
                odds_ratio = odds_ratio_haldane(*table[0], *table[1])
            else:
                fisher_p = np.nan
                odds_ratio = np.nan

            arrays: List[np.ndarray] = []
            for mask, kind in ((mask_a, "num"), (mask_a, "den"), (mask_b, "num"), (mask_b, "den")):
                subset = regions.loc[mask, ["block_1mb"]].copy()
                subset["value"] = y.loc[mask].to_numpy() if kind == "num" else 1
                grouped = subset.groupby("block_1mb")["value"].sum()
                array = np.zeros(len(blocks), dtype=float)
                for block, value in grouped.items():
                    array[block_index[block]] = float(value)
                arrays.append(array)
            bootstrap_inputs.append((arrays[0], arrays[1], arrays[2], arrays[3]))
            rows.append(
                {
                    "contrast": contrast["name"],
                    "factor_column": contrast["column"],
                    "level_a": contrast["level_a"],
                    "level_b": contrast["level_b"],
                    "endpoint": endpoint,
                    "endpoint_scope": definition["scope"],
                    "denominator_a": denominator_a,
                    "numerator_a": numerator_a,
                    "rate_a": rate_a,
                    "denominator_b": denominator_b,
                    "numerator_b": numerator_b,
                    "rate_b": rate_b,
                    "rate_difference_a_minus_b": rate_a - rate_b,
                    "odds_ratio_haldane": odds_ratio,
                    "naive_fisher_p": fisher_p,
                    "bootstrap_replicates": bootstrap_replicates,
                    "bootstrap_unit": "chromosome_plus_1mb_block",
                }
            )

    for row, arrays in zip(rows, bootstrap_inputs):
        num_a, den_a, num_b, den_b = arrays
        boot_num_a = weights @ num_a
        boot_den_a = weights @ den_a
        boot_num_b = weights @ num_b
        boot_den_b = weights @ den_b
        valid = (boot_den_a > 0) & (boot_den_b > 0)
        differences = boot_num_a[valid] / boot_den_a[valid] - boot_num_b[valid] / boot_den_b[valid]
        row["bootstrap_valid_replicates"] = int(valid.sum())
        row["rate_difference_ci95_low"] = float(np.quantile(differences, 0.025)) if len(differences) else np.nan
        row["rate_difference_ci95_high"] = float(np.quantile(differences, 0.975)) if len(differences) else np.nan

    frame = pd.DataFrame(rows)
    frame["bh_q_exploratory"] = bh_adjust(frame["naive_fisher_p"].to_numpy(dtype=float))
    frame["inference_warning"] = (
        "same-read diagnostic; Fisher/BH assumes independent regions and is exploratory; "
        "block CI captures genomic resampling only"
    )
    return frame


def build_multivariable_models(regions: pd.DataFrame) -> pd.DataFrame:
    """Fit diagnostic clustered logistic models, never accuracy models."""
    frame = regions.copy()
    frame["x_log2_mean_pair_min_depth"] = np.log2(frame["raw_baseline_mean_pair_min_dp"])
    frame["x_depth_imbalance_per_0_1"] = frame["raw_baseline_mean_abs_log2_dp_ratio"] / 0.1
    frame["x_abs_vaf_delta_per_0_05"] = frame["raw_baseline_mean_abs_vaf_delta"] / 0.05
    frame["x_k3_vs_k2"] = (frame["k_bin"] == "k=3").astype(int)
    frame["x_k4_vs_k2"] = (frame["k_bin"] == "k=4").astype(int)
    frame["x_kge5_vs_k2"] = (frame["k_bin"] == "k>=5").astype(int)
    frame["x_site_set_mismatch"] = (frame["site_set_status"] == "mismatch").astype(int)
    frame["x_cn_full_any_loh"] = (frame["cn_loh_full_factor"] == "any_loh").astype(int)
    frame["x_cn_partial_or_unavailable"] = (
        frame["cn_loh_full_factor"] == "unavailable_or_partial"
    ).astype(int)
    frame["x_latest_partial"] = (frame["latest_raw_status"] == "partial").astype(int)
    frame["x_hp_count_mismatch"] = (frame["hp_count_status"] == "mismatch").astype(int)
    frame["x_read_candidate_log2"] = np.log2(frame["read_candidate_topology_max"])
    frame["x_vaf_candidate_log2"] = np.log2(frame["vaf_candidate_topology_max"])

    base_predictors = [
        "x_log2_mean_pair_min_depth",
        "x_depth_imbalance_per_0_1",
        "x_abs_vaf_delta_per_0_05",
        "x_k3_vs_k2",
        "x_k4_vs_k2",
        "x_kge5_vs_k2",
        "x_site_set_mismatch",
        "x_cn_full_any_loh",
        "x_cn_partial_or_unavailable",
        "x_latest_partial",
    ]
    all_rows = pd.Series(True, index=frame.index)
    hp_matched = frame["hp_count_status"] == "matched"
    specs = [
        {
            "model": "coarse_agree_context_fixed",
            "outcome": "coarse_agree",
            "mask": all_rows,
            "scope": "fixed_5720",
            "role": "measurement_context",
            "predictors": base_predictors + ["x_hp_count_mismatch"],
        },
        {
            "model": "read_evaluable_context_hp_matched",
            "outcome": "read_evaluable",
            "mask": hp_matched,
            "scope": "hp_count_matched_only",
            "role": "measurement_context",
            "predictors": base_predictors,
        },
        {
            "model": "vaf_evaluable_context_hp_matched",
            "outcome": "vaf_evaluable",
            "mask": hp_matched,
            "scope": "hp_count_matched_only",
            "role": "measurement_context",
            "predictors": base_predictors,
        },
        {
            "model": "read_strict_context_given_evaluable",
            "outcome": "read_strict",
            "mask": frame["read_evaluable"],
            "scope": "read_evaluable_only",
            "role": "measurement_context",
            "predictors": base_predictors,
        },
        {
            "model": "read_strict_plus_candidate_given_evaluable",
            "outcome": "read_strict",
            "mask": frame["read_evaluable"],
            "scope": "read_evaluable_only",
            "role": "adds_structural_identifiability_descriptor",
            "predictors": base_predictors + ["x_read_candidate_log2"],
        },
        {
            "model": "vaf_strict_context_given_evaluable",
            "outcome": "vaf_strict",
            "mask": frame["vaf_evaluable"],
            "scope": "vaf_evaluable_only",
            "role": "measurement_context",
            "predictors": base_predictors,
        },
        {
            "model": "vaf_strict_plus_candidate_given_evaluable",
            "outcome": "vaf_strict",
            "mask": frame["vaf_evaluable"],
            "scope": "vaf_evaluable_only",
            "role": "adds_structural_identifiability_descriptor",
            "predictors": base_predictors + ["x_vaf_candidate_log2"],
        },
        {
            "model": "vaf_conflict_context_given_evaluable",
            "outcome": "vaf_conflict",
            "mask": frame["vaf_evaluable"],
            "scope": "vaf_evaluable_only",
            "role": "measurement_context",
            "predictors": base_predictors,
        },
        {
            "model": "vaf_conflict_plus_candidate_given_evaluable",
            "outcome": "vaf_conflict",
            "mask": frame["vaf_evaluable"],
            "scope": "vaf_evaluable_only",
            "role": "adds_structural_identifiability_descriptor",
            "predictors": base_predictors + ["x_vaf_candidate_log2"],
        },
    ]
    labels = {
        "const": "intercept",
        "x_log2_mean_pair_min_depth": "per_2x_mean_pair_min_depth",
        "x_depth_imbalance_per_0_1": "per_0.1_mean_abs_log2_depth_ratio",
        "x_abs_vaf_delta_per_0_05": "per_0.05_mean_abs_vaf_delta",
        "x_k3_vs_k2": "k3_vs_k2",
        "x_k4_vs_k2": "k4_vs_k2",
        "x_kge5_vs_k2": "k_ge_5_vs_k2",
        "x_site_set_mismatch": "site_set_mismatch_vs_equal",
        "x_cn_full_any_loh": "full_cn_any_loh_vs_full_no_loh",
        "x_cn_partial_or_unavailable": "cn_partial_or_unavailable_vs_full_no_loh",
        "x_latest_partial": "latest_partial_vs_complete_site_coverage",
        "x_hp_count_mismatch": "hp_count_mismatch_vs_matched",
        "x_read_candidate_log2": "per_2x_read_candidate_topology_max",
        "x_vaf_candidate_log2": "per_2x_vaf_candidate_topology_max",
    }

    rows: List[Dict[str, Any]] = []
    for spec in specs:
        columns = [spec["outcome"], "block_1mb"] + spec["predictors"]
        subset = frame.loc[bool_series(spec["mask"]), columns].dropna().copy()
        design = sm.add_constant(subset[spec["predictors"]], has_constant="add")
        model = sm.GLM(
            subset[spec["outcome"]].astype(int),
            design,
            family=sm.families.Binomial(),
        )
        result = model.fit(
            cov_type="cluster",
            cov_kwds={"groups": subset["block_1mb"]},
        )
        confidence = result.conf_int(alpha=0.05)
        for predictor in result.params.index:
            coefficient = float(result.params[predictor])
            low = float(confidence.loc[predictor, 0])
            high = float(confidence.loc[predictor, 1])
            rows.append(
                {
                    "model": spec["model"],
                    "model_role": spec["role"],
                    "outcome": spec["outcome"],
                    "scope": spec["scope"],
                    "n": int(len(subset)),
                    "positive_n": int(subset[spec["outcome"]].sum()),
                    "block_clusters": int(subset["block_1mb"].nunique()),
                    "converged": bool(result.converged),
                    "deviance_fraction_explained": float(
                        1.0 - result.deviance / result.null_deviance
                    ),
                    "predictor": predictor,
                    "predictor_interpretation": labels[predictor],
                    "coefficient": coefficient,
                    "cluster_robust_se": float(result.bse[predictor]),
                    "odds_ratio": float(np.exp(coefficient)),
                    "odds_ratio_ci95_low": float(np.exp(low)),
                    "odds_ratio_ci95_high": float(np.exp(high)),
                    "cluster_robust_p": float(result.pvalues[predictor]),
                    "warning": (
                        "diagnostic model, not causal or accuracy evidence; same-read VAF/depth circularity, "
                        "conditional selection, shared CN, and residual genomic dependence remain"
                    ),
                }
            )

    output = pd.DataFrame(rows)
    non_intercept = output["predictor"] != "const"
    output["bh_q_across_all_nonintercept_coefficients"] = np.nan
    output.loc[non_intercept, "bh_q_across_all_nonintercept_coefficients"] = bh_adjust(
        output.loc[non_intercept, "cluster_robust_p"].to_numpy(dtype=float)
    )
    return output


def add_check(rows: List[Dict[str, Any]], name: str, observed: Any, expected: Any) -> None:
    rows.append(
        {
            "check": name,
            "observed": observed,
            "expected": expected,
            "pass": str(observed) == str(expected),
        }
    )


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> None:
    args = parse_args()
    repo = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    topology_dir = repo / "research/20260712_hcc1395_pair_site_topology_containment_validation/data"
    validation_dir = repo / "research/20260712_vaf_ccf_subclone_crosssoftware_validation"
    paths = {
        "outcomes": topology_dir / "hcc1395_site_topology_pair_outcomes.tsv",
        "alleles": topology_dir / "hcc1395_site_allele_identity.tsv.gz",
        "raw": validation_dir / "results/raw_vaf_validation_v1/data/hcc1395_pair_exact_shared_records.tsv.gz",
        "coarse": repo
        / "research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/hcc1395_pair_matches.tsv",
        "cn": validation_dir / "data/pyclone_inputs/hcc1395_pair_primary_joint_main.site_metadata.tsv.gz",
    }
    for name, path in paths.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing {name}: {path}")

    regions = pd.read_csv(paths["outcomes"], sep="\t")
    if len(regions) != 5720 or regions["match_id"].nunique() != 5720:
        raise ValueError("The fixed topology population is not 5,720 unique regions")
    if regions["fixed_denominator"].nunique() != 1 or int(regions["fixed_denominator"].iloc[0]) != 5720:
        raise ValueError("Unexpected fixed denominator")

    alleles = pd.read_csv(paths["alleles"], sep="\t")
    parsed = alleles["alleles_a"].map(exact_allele_key)
    alleles[["parsed_chrom", "parsed_pos", "ref", "alt"]] = pd.DataFrame(parsed.tolist(), index=alleles.index)
    if not bool_series(alleles["allele_identity"]).all():
        raise ValueError("Non-identical alleles found in fixed shared-site mapping")
    if not (
        (alleles["parsed_chrom"] == alleles["chrom"])
        & (alleles["parsed_pos"] == alleles["position"])
    ).all():
        raise ValueError("Allele JSON and coordinate columns disagree")
    alleles["key"] = (
        alleles["chrom"].astype(str)
        + ":"
        + alleles["position"].astype(int).astype(str)
        + ":"
        + alleles["ref"].astype(str)
        + ":"
        + alleles["alt"].astype(str)
    )
    if alleles["key"].duplicated().any() or len(alleles) != 15713:
        raise ValueError("Expected 15,713 unique exact shared alleles")
    allele_counts = alleles.groupby("match_id").size().rename("allele_identity_site_count").reset_index()
    regions = regions.merge(allele_counts, on="match_id", how="left", validate="one_to_one")

    raw = pd.read_csv(paths["raw"], sep="\t")
    regions, baseline_join = add_raw_region_metrics(
        regions, alleles, raw, "baseline_caller_pass", "baseline"
    )
    regions, latest_join = add_raw_region_metrics(regions, alleles, raw, "latest_lps_pass", "latest")

    coarse = pd.read_csv(paths["coarse"], sep="\t")
    coarse = coarse.loc[
        (coarse["scenario"] == "exact_coordinate")
        & bool_series(coarse["complete_a"])
        & bool_series(coarse["complete_b"]),
        ["match_id", "category_a", "category_b", "category_agree"],
    ].copy()
    if len(coarse) != 5720 or coarse["match_id"].nunique() != 5720:
        raise ValueError("Expected 5,720 exact-coordinate complete-both coarse matches")
    coarse = coarse.rename(
        columns={
            "category_a": "coarse_category_hcc1395",
            "category_b": "coarse_category_dorado",
            "category_agree": "coarse_agree",
        }
    )
    regions = regions.merge(coarse, on="match_id", how="left", validate="one_to_one")
    regions["coarse_agree"] = bool_series(regions["coarse_agree"])

    metadata = pd.read_csv(paths["cn"], sep="\t")
    regions, cn_join = add_cn_context(regions, alleles, metadata)

    regions["read_evaluable"] = bool_series(regions["read_full_evaluable"])
    regions["vaf_evaluable"] = bool_series(regions["vaf_official_evaluable"])
    regions["read_full_exact"] = regions["read_full_outcome"] == "strict_full_exact"
    regions["vaf_full_exact"] = regions["vaf_official_outcome"] == "strict_full_exact"
    regions["read_strict"] = regions["read_full_outcome"].isin(STRICT_OUTCOMES)
    regions["vaf_strict"] = regions["vaf_official_outcome"].isin(STRICT_OUTCOMES)
    regions["read_conflict"] = regions["read_full_outcome"] == "conflict"
    regions["vaf_conflict"] = regions["vaf_official_outcome"] == "conflict"
    regions["all_eligible"] = True
    regions["read_outcome_group"] = regions["read_full_outcome"].map(outcome_group)
    regions["vaf_outcome_group"] = regions["vaf_official_outcome"].map(outcome_group)
    regions["coarse_outcome_group"] = np.where(regions["coarse_agree"], "agree", "different")
    regions["site_set_status"] = np.where(
        regions["caller_site_set_relation"] == "equal", "equal", "mismatch"
    )
    regions["hp_count_status"] = np.where(
        bool_series(regions["hp_count_mismatch"]), "mismatch", "matched"
    )
    regions["caller_extra_site_count"] = (
        regions["caller_k_a"] + regions["caller_k_b"] - 2 * regions["caller_shared_k"]
    )
    regions["read_candidate_topology_max"] = regions[
        ["read_full_topology_set_size_a", "read_full_topology_set_size_b"]
    ].max(axis=1)
    regions["vaf_candidate_topology_max"] = regions[
        ["vaf_official_topology_set_size_a", "vaf_official_topology_set_size_b"]
    ].max(axis=1)
    regions["read_candidate_multiplicity_bin"] = regions["read_candidate_topology_max"].map(
        multiplicity_bin
    )
    regions["vaf_candidate_multiplicity_bin"] = regions["vaf_candidate_topology_max"].map(
        lambda value: multiplicity_bin(value, vaf=True)
    )
    regions["baseline_binomial_status"] = np.select(
        [
            regions["raw_baseline_binomial_testable_n"] == 0,
            regions["raw_baseline_binomial_all_sites_compatible"],
        ],
        ["not_testable", "all_testable_sites_compatible"],
        default="one_or_more_testable_incompatible",
    )
    regions["latest_raw_status"] = np.where(
        regions["raw_latest_complete"], "complete", "partial"
    )
    regions["region_start"] = regions["region"].str.extract(r":(\d+)-", expand=False).astype(int)
    regions["block_1mb"] = regions["chrom"].astype(str) + ":" + (
        regions["region_start"] // 1_000_000
    ).astype(str)

    regions["baseline_depth_quartile"], depth_bins = quartile_labels(
        regions["raw_baseline_mean_pair_min_dp"], "depth"
    )
    regions["baseline_depth_imbalance_quartile"], imbalance_bins = quartile_labels(
        regions["raw_baseline_mean_abs_log2_dp_ratio"], "depth_imbalance"
    )
    regions["baseline_vaf_delta_quartile"], vaf_delta_bins = quartile_labels(
        regions["raw_baseline_mean_abs_vaf_delta"], "vaf_delta"
    )

    endpoints = endpoint_definitions()
    factors = {
        "shared_site_complexity_k": "k_bin",
        "caller_site_set_relation": "site_set_status",
        "hp_count_match": "hp_count_status",
        "baseline_mean_pair_min_depth": "baseline_depth_quartile",
        "baseline_depth_imbalance": "baseline_depth_imbalance_quartile",
        "baseline_mean_abs_vaf_delta": "baseline_vaf_delta_quartile",
        "baseline_binomial_compatibility": "baseline_binomial_status",
        "latest_release_site_coverage": "latest_raw_status",
        "read_candidate_topology_multiplicity": "read_candidate_multiplicity_bin",
        "vaf_candidate_topology_multiplicity": "vaf_candidate_multiplicity_bin",
        "shared_cn_loh_context": "cn_status",
    }
    factor_summary = build_factor_summary(regions, factors, endpoints)
    numeric_summary = build_numeric_outcome_summary(regions)

    general_endpoints = list(endpoints)
    contrasts: List[Dict[str, Any]] = [
        {
            "name": "highest_vs_lowest_baseline_depth_quartile",
            "column": "baseline_depth_quartile",
            "level_a": "Q4_depth_high",
            "level_b": "Q1_depth_low",
            "endpoints": general_endpoints,
        },
        {
            "name": "highest_vs_lowest_baseline_depth_imbalance_quartile",
            "column": "baseline_depth_imbalance_quartile",
            "level_a": "Q4_depth_imbalance_high",
            "level_b": "Q1_depth_imbalance_low",
            "endpoints": general_endpoints,
        },
        {
            "name": "highest_vs_lowest_baseline_abs_vaf_delta_quartile",
            "column": "baseline_vaf_delta_quartile",
            "level_a": "Q4_vaf_delta_high",
            "level_b": "Q1_vaf_delta_low",
            "endpoints": general_endpoints,
        },
        {
            "name": "k_ge_5_vs_k_2",
            "column": "k_bin",
            "level_a": "k>=5",
            "level_b": "k=2",
            "endpoints": general_endpoints,
        },
        {
            "name": "site_set_mismatch_vs_equal",
            "column": "site_set_status",
            "level_a": "mismatch",
            "level_b": "equal",
            "endpoints": general_endpoints,
        },
        {
            "name": "hp_count_mismatch_vs_matched",
            "column": "hp_count_status",
            "level_a": "mismatch",
            "level_b": "matched",
            "endpoints": general_endpoints,
        },
        {
            "name": "one_or_more_testable_binomial_incompatible_vs_all_testable_compatible",
            "column": "baseline_binomial_status",
            "level_a": "one_or_more_testable_incompatible",
            "level_b": "all_testable_sites_compatible",
            "endpoints": general_endpoints,
        },
        {
            "name": "latest_partial_vs_complete_site_coverage",
            "column": "latest_raw_status",
            "level_a": "partial",
            "level_b": "complete",
            "endpoints": general_endpoints,
        },
        {
            "name": "full_cn_any_loh_vs_no_loh",
            "column": "cn_loh_full_factor",
            "level_a": "any_loh",
            "level_b": "no_loh",
            "endpoints": general_endpoints,
        },
        {
            "name": "read_candidate_topology_gt20_vs_1",
            "column": "read_candidate_multiplicity_bin",
            "level_a": ">20",
            "level_b": "1",
            "endpoints": [
                "read_full_exact_given_evaluable",
                "read_exact_or_induced_given_evaluable",
                "read_conflict_given_evaluable",
            ],
        },
        {
            "name": "vaf_candidate_topology_ge2_vs_1",
            "column": "vaf_candidate_multiplicity_bin",
            "level_a": ">=2",
            "level_b": "1",
            "endpoints": [
                "vaf_full_exact_given_evaluable",
                "vaf_exact_or_induced_given_evaluable",
                "vaf_conflict_given_evaluable",
            ],
        },
    ]
    contrast_effects = build_contrast_effects(
        regions,
        contrasts,
        endpoints,
        bootstrap_replicates=args.bootstrap_replicates,
        seed=args.seed,
    )
    multivariable_models = build_multivariable_models(regions)

    checks: List[Dict[str, Any]] = []
    add_check(checks, "fixed_region_rows", len(regions), 5720)
    add_check(checks, "fixed_unique_match_ids", regions["match_id"].nunique(), 5720)
    add_check(checks, "exact_shared_alleles", len(alleles), 15713)
    add_check(checks, "allele_identity_pass", int(bool_series(alleles["allele_identity"]).sum()), 15713)
    add_check(checks, "baseline_raw_allele_hits", baseline_join["allele_hits"], 15713)
    add_check(checks, "baseline_raw_complete_regions", baseline_join["complete_regions"], 5720)
    add_check(checks, "latest_raw_allele_hits", latest_join["allele_hits"], 15636)
    add_check(checks, "latest_raw_complete_regions", latest_join["complete_regions"], 5647)
    add_check(checks, "coarse_agree_regions", int(regions["coarse_agree"].sum()), 3969)
    add_check(checks, "read_evaluable_regions", int(regions["read_evaluable"].sum()), 4038)
    add_check(checks, "read_strict_exact_or_induced_regions", int(regions["read_strict"].sum()), 1599)
    add_check(checks, "vaf_evaluable_regions", int(regions["vaf_evaluable"].sum()), 3860)
    add_check(checks, "vaf_strict_exact_or_induced_regions", int(regions["vaf_strict"].sum()), 1790)
    add_check(checks, "vaf_conflict_regions", int(regions["vaf_conflict"].sum()), 1234)
    add_check(checks, "cn_allele_hits", cn_join["allele_hits"], 14369)
    add_check(checks, "cn_complete_regions", cn_join["complete_regions"], 4542)
    checks_frame = pd.DataFrame(checks)
    if not bool_series(checks_frame["pass"]).all():
        failed = checks_frame.loc[~bool_series(checks_frame["pass"])]
        raise AssertionError(f"Audit checks failed:\n{failed.to_string(index=False)}")

    region_path = output_dir / "region_cause_decomposition.tsv.gz"
    factor_path = output_dir / "factor_strata_summary.tsv"
    numeric_path = output_dir / "outcome_numeric_feature_summary.tsv"
    contrast_path = output_dir / "contrast_effects.tsv"
    models_path = output_dir / "multivariable_logistic_models.tsv"
    checks_path = output_dir / "checks.tsv"
    write_tsv_gz(regions.sort_values(["chrom", "region_start", "match_id"]), region_path)
    write_tsv(factor_summary, factor_path)
    write_tsv(numeric_summary, numeric_path)
    write_tsv(contrast_effects, contrast_path)
    write_tsv(multivariable_models, models_path)
    write_tsv(checks_frame, checks_path)

    output_paths = [region_path, factor_path, numeric_path, contrast_path, models_path, checks_path]
    summary = {
        "schema_version": "1.0",
        "task_type": "B_comprehensive_validation_red_team_audit",
        "population": {
            "regions": int(len(regions)),
            "exact_shared_alleles": int(len(alleles)),
            "site_set_equal_regions": int((regions["site_set_status"] == "equal").sum()),
            "site_set_mismatch_regions": int((regions["site_set_status"] == "mismatch").sum()),
            "hp_count_mismatch_regions": int((regions["hp_count_status"] == "mismatch").sum()),
        },
        "endpoints": {
            "coarse_agree": int(regions["coarse_agree"].sum()),
            "read_evaluable": int(regions["read_evaluable"].sum()),
            "read_full_exact": int(regions["read_full_exact"].sum()),
            "read_exact_or_induced": int(regions["read_strict"].sum()),
            "read_conflict": int(regions["read_conflict"].sum()),
            "vaf_evaluable": int(regions["vaf_evaluable"].sum()),
            "vaf_full_exact": int(regions["vaf_full_exact"].sum()),
            "vaf_exact_or_induced": int(regions["vaf_strict"].sum()),
            "vaf_conflict": int(regions["vaf_conflict"].sum()),
        },
        "joins": {
            "baseline_raw": baseline_join,
            "latest_raw": latest_join,
            "shared_savana_cn": cn_join,
        },
        "quartile_boundaries": {
            "baseline_mean_pair_min_depth": depth_bins,
            "baseline_mean_abs_log2_depth_ratio": imbalance_bins,
            "baseline_mean_abs_vaf_delta": vaf_delta_bins,
        },
        "bootstrap": {
            "seed": args.seed,
            "replicates": args.bootstrap_replicates,
            "unit": "chromosome_plus_1mb_block",
        },
        "diagnostic_multivariable_models": {
            "models": int(multivariable_models["model"].nunique()),
            "coefficient_rows": int(len(multivariable_models)),
            "covariance": "cluster_robust_by_chromosome_plus_1mb_block",
            "multiple_testing": "BH across all non-intercept coefficients",
        },
        "claim_ceiling": {
            "same_read_circularity": (
                "Depth, paired VAF, VAF-selected topology, and read topology are derived from overlapping reads; "
                "associations diagnose measurement behavior and do not independently validate a true tree."
            ),
            "selection": (
                "The estimand is conditional on 5,720 exact-coordinate complete-both regions; missing or failed regions "
                "are outside the fixed population."
            ),
            "cn": (
                "CN/LOH is the shared HCC1395 SAVANA context used for both technical sources and is available only "
                "for the primary PyClone selection; it is not source-specific corroboration."
            ),
            "release": (
                "Baseline raw records exactly cover the frozen topology alleles; latest LPS records are a cross-release "
                "sensitivity join with 77 missing alleles across 73 regions."
            ),
            "truth": (
                "No single-cell, multi-region, synthetic, or orthogonal exact-tree truth is available for these regions."
            ),
        },
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in paths.items()},
        "outputs": {path.name: {"path": str(path), "sha256": sha256(path)} for path in output_paths},
        "checks_passed": int(bool_series(checks_frame["pass"]).sum()),
        "checks_total": int(len(checks_frame)),
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(json_ready(summary), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(f"repo_root={repo}")
    print(f"output_dir={output_dir}")
    print(f"regions={len(regions)} alleles={len(alleles)} checks={len(checks_frame)}/{len(checks_frame)} PASS")
    print(
        "baseline_hits={}/{} latest_hits={}/{} cn_hits={}/{}".format(
            baseline_join["allele_hits"],
            baseline_join["allele_expected"],
            latest_join["allele_hits"],
            latest_join["allele_expected"],
            cn_join["allele_hits"],
            cn_join["allele_expected"],
        )
    )
    print(
        "coarse_agree={} read_strict={} vaf_strict={} vaf_conflict={}".format(
            int(regions["coarse_agree"].sum()),
            int(regions["read_strict"].sum()),
            int(regions["vaf_strict"].sum()),
            int(regions["vaf_conflict"].sum()),
        )
    )
    print(f"summary={summary_path}")


if __name__ == "__main__":
    main()
