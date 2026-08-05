#!/usr/bin/env python3
"""Whole-autosome raw-VAF technical reproducibility analysis.

This program intentionally reconstructs ALT counts from FORMAT/DP and FORMAT/AF:

    alt_count = floor(DP * AF + 0.5)
    ref_count = DP - alt_count

FORMAT/AD is retained only as an ALT-count diagnostic.  AD[0] is never used as
the reference count because LongPhase-S recalibrated records can carry AD[0]=0
while DP still represents total depth.
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import platform
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.spatial.distance import jensenshannon


AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
REGIONS = ",".join(AUTOSOMES)
SOURCE_FIELDS = {
    "latest_lps_pass": ("somatic", "tree_vcf", "path"),
    "baseline_caller_pass": ("somatic", "caller_pass_baseline_vcf", "path"),
}
TRUTH_PATHS = {
    "HCC1395": "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz",
    "HCC1395_DORADO": "/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz",
    "COLO829": "/big8_disk/data/COLO829/NYGC/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz",
    "H1437": "/big8_disk/data/H1437/orthogonal-tools-benchmark/H1437_orthogonal-tools-benchmark_somatic-only.vcf.gz",
    "H2009": "/big8_disk/data/H2009/orthogonal-tools-benchmark/H2009_orthogonal-tools-benchmark_somatic-only.vcf.gz",
    "HCC1937": "/big8_disk/data/HCC1937/orthogonal-tools-benchmark/HCC1937_orthogonal-tools-benchmark_somatic-only.vcf.gz",
    "HCC1954": "/big8_disk/data/HCC1954/orthogonal-tools-benchmark/HCC1954_orthogonal-tools-benchmark_somatic-only.vcf.gz",
}
TRUTH_LABELS = {
    "HCC1395": "SEQC2 v1.2.1 HC sSNV",
    "HCC1395_DORADO": "SEQC2 v1.2.1 HC sSNV",
    "COLO829": "NYGC consensus",
    "H1437": "orthogonal-tools benchmark",
    "H2009": "orthogonal-tools benchmark",
    "HCC1937": "orthogonal-tools benchmark",
    "HCC1954": "orthogonal-tools benchmark",
}
SAMPLE_ORDER = [
    "HCC1395",
    "HCC1395_DORADO",
    "COLO829",
    "H1437",
    "H2009",
    "HCC1937",
    "HCC1954",
]
SOURCE_ORDER = ["latest_lps_pass", "baseline_caller_pass"]
VAF_BINS = [0.0, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.0000001]
VAF_BAND_LABELS = ["[0,.05)", "[.05,.10)", "[.10,.20)", "[.20,.40)", "[.40,.60)", "[.60,.80)", "[.80,1]"]
COLORS = {
    "HCC1395": "#16697A",
    "HCC1395_DORADO": "#E07A5F",
    "COLO829": "#6A4C93",
    "H1437": "#2A9D8F",
    "H2009": "#E9C46A",
    "HCC1937": "#457B9D",
    "HCC1954": "#C44536",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--bootstrap-reps", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260712)
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def nested_get(record: Mapping[str, Any], keys: Sequence[str]) -> Any:
    value: Any = record
    for key in keys:
        value = value[key]
    return value


def run_checked(command: Sequence[str]) -> str:
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    return completed.stdout


def vcf_sample_names(vcf: Path) -> List[str]:
    output = run_checked(["bcftools", "query", "-l", str(vcf)])
    return [line.strip() for line in output.splitlines() if line.strip()]


def parse_float(value: str) -> float:
    if value in {"", "."}:
        return float("nan")
    return float(value.split(",")[0])


def parse_int(value: str) -> float:
    if value in {"", "."}:
        return float("nan")
    return float(int(value))


def parse_ad_alt(value: str) -> float:
    """Return only the reported ALT depth; never inspect or use AD[0]."""
    if value in {"", "."}:
        return float("nan")
    parts = value.split(",")
    if len(parts) < 2 or parts[-1] in {"", "."}:
        return float("nan")
    return float(int(parts[-1]))


def read_vaf_vcf(vcf: Path, sample: str, source: str) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    names = vcf_sample_names(vcf)
    if len(names) != 1:
        raise RuntimeError(f"Expected exactly one sample column in {vcf}; found {names}")

    expression = 'TYPE="snp" && N_ALT=1 && FILTER="PASS"'
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%DP\\t%AF\\t%AD]\\n"
    command = [
        "bcftools",
        "query",
        "-r",
        REGIONS,
        "-i",
        expression,
        "-f",
        fmt,
        str(vcf),
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert process.stdout is not None
    rows: List[Dict[str, Any]] = []
    queried_rows = 0
    invalid_missing = 0
    invalid_range = 0
    for line in process.stdout:
        queried_rows += 1
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 8:
            process.kill()
            raise RuntimeError(f"Unexpected bcftools query row ({len(fields)} fields): {line[:300]}")
        chrom, pos, ref, alt, filt, dp_text, af_text, ad_text = fields
        dp_value = parse_int(dp_text)
        af_value = parse_float(af_text)
        if not np.isfinite(dp_value) or not np.isfinite(af_value):
            invalid_missing += 1
            continue
        dp = int(dp_value)
        af = float(af_value)
        if dp <= 0 or af < 0.0 or af > 1.0:
            invalid_range += 1
            continue
        # Half-up rounding implements the declared round(DP*AF) count contract.
        alt_count = int(math.floor(dp * af + 0.5))
        ref_count = dp - alt_count
        if alt_count < 0 or ref_count < 0 or alt_count + ref_count != dp:
            raise RuntimeError(f"Count reconstruction invariant failed: {sample} {chrom}:{pos} {dp=} {af=}")
        ad_alt = parse_ad_alt(ad_text)
        rows.append(
            {
                "sample": sample,
                "source": source,
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "filter": filt,
                "dp": dp,
                "af": af,
                "alt_count": alt_count,
                "ref_count": ref_count,
                "ad_raw": ad_text,
                "ad_alt_reported": ad_alt,
                "alt_minus_ad_alt": (alt_count - ad_alt) if np.isfinite(ad_alt) else float("nan"),
            }
        )
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"bcftools query failed ({return_code}) for {vcf}: {stderr}")

    frame = pd.DataFrame(rows)
    if frame.empty:
        raise RuntimeError(f"No valid PASS biallelic autosomal SNVs found in {vcf}")
    frame["key"] = (
        frame["chrom"].astype(str)
        + ":"
        + frame["pos"].astype(str)
        + ":"
        + frame["ref"].astype(str)
        + ":"
        + frame["alt"].astype(str)
    )
    duplicate_mask = frame.duplicated("key", keep=False)
    duplicate_rows = int(duplicate_mask.sum())
    duplicate_keys = int(frame.loc[duplicate_mask, "key"].nunique())
    if duplicate_keys:
        conflict = (
            frame.loc[duplicate_mask]
            .groupby("key", sort=False)[["dp", "af", "alt_count", "ref_count"]]
            .nunique(dropna=False)
            .gt(1)
            .any(axis=1)
        )
        if bool(conflict.any()):
            examples = conflict[conflict].index[:5].tolist()
            raise RuntimeError(f"Conflicting duplicate exact keys in {vcf}: {examples}")
        frame = frame.drop_duplicates("key", keep="first").copy()

    af_reconstructed = frame["alt_count"] / frame["dp"]
    # FORMAT/AF is rounded in the VCF, so equality is validated at count scale.
    count_reconstruction_ok = (
        frame["alt_count"] == np.floor(frame["dp"] * frame["af"] + 0.5).astype(int)
    ) & (frame["ref_count"] == frame["dp"] - frame["alt_count"])
    ad_comparable = frame["ad_alt_reported"].notna()
    stats_record = {
        "sample": sample,
        "source": source,
        "path": str(vcf),
        "bcftools_filter": expression,
        "bcftools_sample_column": names[0],
        "queried_rows": queried_rows,
        "valid_unique_rows": int(len(frame)),
        "invalid_missing_dp_or_af": invalid_missing,
        "invalid_dp_or_af_range": invalid_range,
        "duplicate_rows": duplicate_rows,
        "duplicate_keys": duplicate_keys,
        "count_reconstruction_ok_n": int(count_reconstruction_ok.sum()),
        "count_reconstruction_ok_fraction": float(count_reconstruction_ok.mean()),
        "af_quantization_max_abs_error": float(np.max(np.abs(af_reconstructed - frame["af"]))),
        "ad_alt_comparable_n": int(ad_comparable.sum()),
        "ad_alt_exact_n": int((frame.loc[ad_comparable, "alt_minus_ad_alt"] == 0).sum()),
        "ad_alt_exact_fraction": float((frame.loc[ad_comparable, "alt_minus_ad_alt"] == 0).mean()) if ad_comparable.any() else None,
        "reference_count_definition": "DP - round_half_up(DP*AF); FORMAT/AD[0] never used",
    }
    return frame.sort_values(["chrom", "pos", "ref", "alt"]).reset_index(drop=True), stats_record


def read_truth_keys(vcf: Path) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    expression = 'TYPE="snp" && N_ALT=1'
    command = [
        "bcftools",
        "query",
        "-r",
        REGIONS,
        "-i",
        expression,
        "-f",
        "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\n",
        str(vcf),
    ]
    output = run_checked(command)
    rows = []
    for line in output.splitlines():
        chrom, pos, ref, alt, filt = line.split("\t")
        rows.append((chrom, int(pos), ref, alt, filt))
    frame = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "filter"])
    frame["key"] = frame["chrom"] + ":" + frame["pos"].astype(str) + ":" + frame["ref"] + ":" + frame["alt"]
    duplicate_keys = int(frame.duplicated("key", keep=False).sum())
    frame = frame.drop_duplicates("key").sort_values(["chrom", "pos", "ref", "alt"]).reset_index(drop=True)
    return frame, {
        "path": str(vcf),
        "bcftools_filter": expression,
        "biallelic_autosomal_snv_count": int(len(frame)),
        "duplicate_rows_removed": duplicate_keys,
    }


def distribution_summary(frame: pd.DataFrame, scope: str) -> Dict[str, Any]:
    return {
        "scope": scope,
        "n": int(len(frame)),
        "vaf_mean": float(frame["af"].mean()) if len(frame) else float("nan"),
        "vaf_std": float(frame["af"].std(ddof=1)) if len(frame) > 1 else float("nan"),
        "vaf_q05": float(frame["af"].quantile(0.05)) if len(frame) else float("nan"),
        "vaf_q25": float(frame["af"].quantile(0.25)) if len(frame) else float("nan"),
        "vaf_median": float(frame["af"].median()) if len(frame) else float("nan"),
        "vaf_q75": float(frame["af"].quantile(0.75)) if len(frame) else float("nan"),
        "vaf_q95": float(frame["af"].quantile(0.95)) if len(frame) else float("nan"),
        "dp_median": float(frame["dp"].median()) if len(frame) else float("nan"),
        "dp_q05": float(frame["dp"].quantile(0.05)) if len(frame) else float("nan"),
        "dp_q95": float(frame["dp"].quantile(0.95)) if len(frame) else float("nan"),
    }


def histogram_rows(frame: pd.DataFrame, sample: str, source: str, scope: str) -> List[Dict[str, Any]]:
    edges = np.linspace(0.0, 1.0, 51)
    counts, _ = np.histogram(frame["af"].to_numpy(dtype=float), bins=edges)
    total = int(counts.sum())
    result = []
    for index, count in enumerate(counts):
        result.append(
            {
                "sample": sample,
                "source": source,
                "scope": scope,
                "bin_index": index,
                "bin_left": float(edges[index]),
                "bin_right": float(edges[index + 1]),
                "count": int(count),
                "fraction": float(count / total) if total else float("nan"),
            }
        )
    return result


def vaf_band_rows(frame: pd.DataFrame, sample: str, source: str, scope: str) -> List[Dict[str, Any]]:
    categories = pd.cut(frame["af"], bins=VAF_BINS, labels=VAF_BAND_LABELS, right=False, include_lowest=True)
    counts = categories.value_counts(sort=False).reindex(VAF_BAND_LABELS, fill_value=0)
    total = int(counts.sum())
    return [
        {
            "sample": sample,
            "source": source,
            "scope": scope,
            "band": band,
            "count": int(counts.loc[band]),
            "fraction": float(counts.loc[band] / total) if total else float("nan"),
        }
        for band in VAF_BAND_LABELS
    ]


def concordance_correlation(x: np.ndarray, y: np.ndarray) -> float:
    mean_x = float(np.mean(x))
    mean_y = float(np.mean(y))
    var_x = float(np.mean((x - mean_x) ** 2))
    var_y = float(np.mean((y - mean_y) ** 2))
    covariance = float(np.mean((x - mean_x) * (y - mean_y)))
    denominator = var_x + var_y + (mean_x - mean_y) ** 2
    return float(2.0 * covariance / denominator) if denominator > 0 else float("nan")


def paired_metric_values(union: pd.DataFrame) -> Dict[str, float]:
    present_a = union["present_hcc1395"].to_numpy(dtype=bool)
    present_b = union["present_dorado"].to_numpy(dtype=bool)
    shared_mask = present_a & present_b
    n_union = int(len(union))
    n_shared = int(shared_mask.sum())
    values: Dict[str, float] = {
        "callset_jaccard": float(n_shared / n_union) if n_union else float("nan"),
    }
    if n_shared < 2:
        return values
    shared = union.loc[shared_mask]
    x = shared["af_hcc1395"].to_numpy(dtype=float)
    y = shared["af_dorado"].to_numpy(dtype=float)
    delta = y - x
    values.update(
        {
            "vaf_pearson": float(np.corrcoef(x, y)[0, 1]),
            "vaf_spearman": float(stats.spearmanr(x, y).statistic),
            "vaf_ccc": concordance_correlation(x, y),
            "vaf_mae": float(np.mean(np.abs(delta))),
            "vaf_median_abs_delta": float(np.median(np.abs(delta))),
            "vaf_mean_signed_delta_dorado_minus_hcc1395": float(np.mean(delta)),
            "vaf_median_signed_delta_dorado_minus_hcc1395": float(np.median(delta)),
            "vaf_within_0.02": float(np.mean(np.abs(delta) <= 0.02)),
            "vaf_within_0.05": float(np.mean(np.abs(delta) <= 0.05)),
            "vaf_within_0.10": float(np.mean(np.abs(delta) <= 0.10)),
            "vaf_ks_statistic": float(stats.ks_2samp(x, y, alternative="two-sided", method="auto").statistic),
            "vaf_wasserstein": float(stats.wasserstein_distance(x, y)),
        }
    )
    hist_x, _ = np.histogram(x, bins=np.linspace(0.0, 1.0, 51))
    hist_y, _ = np.histogram(y, bins=np.linspace(0.0, 1.0, 51))
    values["vaf_js_divergence_50bin_nats"] = float(
        jensenshannon(hist_x.astype(float) + 1e-12, hist_y.astype(float) + 1e-12, base=np.e) ** 2
    )

    alt_a = shared["alt_count_hcc1395"].to_numpy(dtype=float)
    alt_b = shared["alt_count_dorado"].to_numpy(dtype=float)
    dp_a = shared["dp_hcc1395"].to_numpy(dtype=float)
    dp_b = shared["dp_dorado"].to_numpy(dtype=float)
    pooled_p = (alt_a + alt_b) / (dp_a + dp_b)
    expected_variance = pooled_p * (1.0 - pooled_p) * (1.0 / dp_a + 1.0 / dp_b)
    finite = np.isfinite(expected_variance) & (expected_variance > 0)
    if finite.any():
        observed_mse = float(np.mean(delta[finite] ** 2))
        expected_mse = float(np.mean(expected_variance[finite]))
        values["binomial_noise_variance_ratio"] = observed_mse / expected_mse
        values["binomial_noise_rmse_ratio"] = math.sqrt(observed_mse / expected_mse)
        values["binomial_95pct_compatible_fraction"] = float(
            np.mean(np.abs(delta[finite]) <= 1.96 * np.sqrt(expected_variance[finite]))
        )
    else:
        values["binomial_noise_variance_ratio"] = float("nan")
        values["binomial_noise_rmse_ratio"] = float("nan")
        values["binomial_95pct_compatible_fraction"] = float("nan")
    return values


def build_pair_union(frame_a: pd.DataFrame, frame_b: pd.DataFrame) -> pd.DataFrame:
    columns = ["key", "chrom", "pos", "ref", "alt", "af", "dp", "alt_count", "ref_count"]
    a = frame_a[columns].rename(
        columns={
            "chrom": "chrom_hcc1395",
            "pos": "pos_hcc1395",
            "ref": "ref_hcc1395",
            "alt": "alt_hcc1395",
            "af": "af_hcc1395",
            "dp": "dp_hcc1395",
            "alt_count": "alt_count_hcc1395",
            "ref_count": "ref_count_hcc1395",
        }
    )
    b = frame_b[columns].rename(
        columns={
            "chrom": "chrom_dorado",
            "pos": "pos_dorado",
            "ref": "ref_dorado",
            "alt": "alt_dorado",
            "af": "af_dorado",
            "dp": "dp_dorado",
            "alt_count": "alt_count_dorado",
            "ref_count": "ref_count_dorado",
        }
    )
    union = a.merge(b, on="key", how="outer", validate="one_to_one")
    union["present_hcc1395"] = union["af_hcc1395"].notna()
    union["present_dorado"] = union["af_dorado"].notna()
    union["chrom"] = union["chrom_hcc1395"].fillna(union["chrom_dorado"])
    union["pos"] = union["pos_hcc1395"].fillna(union["pos_dorado"]).astype(int)
    union["ref"] = union["ref_hcc1395"].fillna(union["ref_dorado"])
    union["alt"] = union["alt_hcc1395"].fillna(union["alt_dorado"])
    union["block_1mb"] = union["chrom"] + ":" + ((union["pos"] - 1) // 1_000_000).astype(str)
    return union.sort_values(["chrom", "pos", "ref", "alt"]).reset_index(drop=True)


def block_bootstrap_metrics(
    union: pd.DataFrame, reps: int, seed: int
) -> Tuple[Dict[str, Tuple[float, float, int]], List[Dict[str, float]]]:
    block_codes, block_names = pd.factorize(union["block_1mb"], sort=True)
    block_count = len(block_names)
    row_indices = np.arange(len(union), dtype=int)
    rng = np.random.default_rng(seed)
    distributions: MutableMapping[str, List[float]] = {}
    flat_rows: List[Dict[str, float]] = []
    for replicate in range(reps):
        sampled_block_counts = rng.multinomial(block_count, np.full(block_count, 1.0 / block_count))
        row_weights = sampled_block_counts[block_codes]
        sampled_indices = np.repeat(row_indices, row_weights)
        sampled = union.iloc[sampled_indices]
        values = paired_metric_values(sampled)
        row: Dict[str, float] = {"bootstrap_replicate": float(replicate)}
        for metric, value in values.items():
            row[metric] = value
            if np.isfinite(value):
                distributions.setdefault(metric, []).append(float(value))
        flat_rows.append(row)
    intervals: Dict[str, Tuple[float, float, int]] = {}
    for metric, values in distributions.items():
        array = np.asarray(values, dtype=float)
        intervals[metric] = (float(np.quantile(array, 0.025)), float(np.quantile(array, 0.975)), int(len(array)))
    return intervals, flat_rows


def save_figure(fig: plt.Figure, output_base: Path) -> None:
    fig.savefig(output_base.with_suffix(".png"), dpi=180, bbox_inches="tight", metadata={"Software": "InterSubMod raw VAF validation"})
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def plot_truth_distributions(histogram: pd.DataFrame, output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=True)
    for axis, source in zip(axes, SOURCE_ORDER):
        subset = histogram[(histogram["source"] == source) & (histogram["scope"] == "truth_confirmed")]
        for sample in SAMPLE_ORDER:
            values = subset[subset["sample"] == sample].sort_values("bin_index")
            if values.empty:
                continue
            centers = (values["bin_left"].to_numpy() + values["bin_right"].to_numpy()) / 2.0
            density = values["fraction"].to_numpy() / (values["bin_right"].to_numpy() - values["bin_left"].to_numpy())
            n = int(values["count"].sum())
            axis.plot(centers, density, lw=1.8, color=COLORS[sample], label=f"{sample} (n={n:,})")
        axis.set_title("Latest LongPhase-S PASS" if source == "latest_lps_pass" else "Caller-PASS baseline sensitivity")
        axis.set_xlabel("Raw caller VAF")
        axis.grid(alpha=0.22, lw=0.7)
    axes[0].set_ylabel("Density (0.02-wide bins)")
    axes[1].legend(fontsize=8, frameon=False, loc="upper right")
    fig.suptitle("Truth-confirmed autosomal biallelic SNV raw-VAF distributions\nTechnical profiles only; not purity/CN-corrected clone fractions", fontsize=13)
    fig.tight_layout()
    save_figure(fig, output_dir / "all_sample_truth_confirmed_vaf_distribution")


def plot_band_composition(bands: pd.DataFrame, output_dir: Path) -> None:
    palette = ["#264653", "#2A9D8F", "#8AB17D", "#E9C46A", "#F4A261", "#E76F51", "#8C2F39"]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
    for axis, source in zip(axes, SOURCE_ORDER):
        subset = bands[(bands["source"] == source) & (bands["scope"] == "truth_confirmed")]
        bottoms = np.zeros(len(SAMPLE_ORDER))
        for band, color in zip(VAF_BAND_LABELS, palette):
            values = []
            for sample in SAMPLE_ORDER:
                row = subset[(subset["sample"] == sample) & (subset["band"] == band)]
                values.append(float(row["fraction"].iloc[0]) if len(row) else 0.0)
            axis.bar(SAMPLE_ORDER, values, bottom=bottoms, color=color, label=band, width=0.78)
            bottoms += np.asarray(values)
        axis.set_title("Latest LongPhase-S PASS" if source == "latest_lps_pass" else "Caller-PASS baseline sensitivity")
        axis.tick_params(axis="x", rotation=35, labelsize=8)
        axis.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Fraction of truth-confirmed calls")
    axes[1].legend(title="Raw VAF band", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=8)
    fig.suptitle("Raw-VAF band composition by dataset\nBands are descriptive and must not be interpreted as clone labels", fontsize=13)
    fig.tight_layout()
    save_figure(fig, output_dir / "all_sample_truth_confirmed_vaf_band_composition")


def plot_hcc_hexbin(pair_unions: Mapping[str, pd.DataFrame], output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
    for axis, source in zip(axes, SOURCE_ORDER):
        union = pair_unions[source]
        shared = union[union["present_hcc1395"] & union["present_dorado"]]
        x = shared["af_hcc1395"].to_numpy(dtype=float)
        y = shared["af_dorado"].to_numpy(dtype=float)
        artist = axis.hexbin(x, y, gridsize=45, extent=(0, 1, 0, 1), bins="log", mincnt=1, cmap="viridis")
        axis.plot([0, 1], [0, 1], ls="--", color="#D1495B", lw=1.3)
        pearson = float(np.corrcoef(x, y)[0, 1])
        ccc = concordance_correlation(x, y)
        axis.text(0.03, 0.97, f"n={len(shared):,}\nr={pearson:.3f}\nCCC={ccc:.3f}", transform=axis.transAxes, va="top", fontsize=9,
                  bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})
        axis.set_title("Latest LongPhase-S PASS" if source == "latest_lps_pass" else "Caller-PASS baseline sensitivity")
        axis.set_xlabel("HCC1395 raw VAF")
        axis.grid(alpha=0.15)
        fig.colorbar(artist, ax=axis, label="log10 bin count")
    axes[0].set_ylabel("HCC1395_DORADO raw VAF")
    fig.suptitle("Exact shared-allele raw-VAF concordance\nSame cell line, different sequencing/basecalling source", fontsize=13)
    fig.tight_layout()
    save_figure(fig, output_dir / "hcc1395_pair_exact_shared_vaf_hexbin")


def plot_bland_altman(pair_unions: Mapping[str, pd.DataFrame], output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
    for axis, source in zip(axes, SOURCE_ORDER):
        union = pair_unions[source]
        shared = union[union["present_hcc1395"] & union["present_dorado"]]
        x = shared["af_hcc1395"].to_numpy(dtype=float)
        y = shared["af_dorado"].to_numpy(dtype=float)
        mean = (x + y) / 2.0
        delta = y - x
        artist = axis.hexbin(mean, delta, gridsize=50, extent=(0, 1, -1, 1), bins="log", mincnt=1, cmap="magma")
        bias = float(np.mean(delta))
        sd = float(np.std(delta, ddof=1))
        axis.axhline(0, color="#666666", lw=0.9)
        axis.axhline(bias, color="#0077B6", lw=1.4, label=f"bias={bias:.3f}")
        axis.axhline(bias + 1.96 * sd, color="#D1495B", ls="--", lw=1.1)
        axis.axhline(bias - 1.96 * sd, color="#D1495B", ls="--", lw=1.1, label=f"±1.96 SD ({1.96*sd:.3f})")
        axis.set_title("Latest LongPhase-S PASS" if source == "latest_lps_pass" else "Caller-PASS baseline sensitivity")
        axis.set_xlabel("Mean raw VAF")
        axis.legend(frameon=False, fontsize=8, loc="lower left")
        axis.grid(alpha=0.15)
        fig.colorbar(artist, ax=axis, label="log10 bin count")
    axes[0].set_ylabel("DORADO − HCC1395 raw VAF")
    fig.suptitle("Bland–Altman view of exact shared alleles\nLimits quantify technical disagreement, not biological clone divergence", fontsize=13)
    fig.tight_layout()
    save_figure(fig, output_dir / "hcc1395_pair_exact_shared_vaf_bland_altman")


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", lineterminator="\n")


def write_tsv_gz(frame: pd.DataFrame, path: Path) -> None:
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            frame.to_csv(compressed, sep="\t", index=False, na_rep="NA", lineterminator="\n")


def main() -> None:
    args = parse_args()
    manifest_path = args.manifest.resolve()
    output_dir = args.output_dir.resolve()
    data_dir = output_dir / "data"
    figure_dir = output_dir / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.hashsalt": str(args.seed)})

    with manifest_path.open() as handle:
        manifest = json.load(handle)
    samples = {entry["sample"]: entry for entry in manifest["samples"]}
    if list(samples) != SAMPLE_ORDER:
        missing = sorted(set(SAMPLE_ORDER) - set(samples))
        extra = sorted(set(samples) - set(SAMPLE_ORDER))
        if missing or extra:
            raise RuntimeError(f"Manifest sample mismatch: missing={missing}, extra={extra}")

    input_provenance: Dict[str, Any] = {
        "manifest": {"path": str(manifest_path), "sha256": sha256_file(manifest_path)},
        "caller_vcfs": {},
        "truth_vcfs": {},
    }
    truth_frames: Dict[str, pd.DataFrame] = {}
    truth_stats: Dict[str, Dict[str, Any]] = {}
    truth_path_cache: Dict[str, Tuple[pd.DataFrame, Dict[str, Any], str]] = {}
    for sample in SAMPLE_ORDER:
        path = Path(TRUTH_PATHS[sample])
        if not path.exists():
            raise FileNotFoundError(f"Truth VCF unavailable for {sample}: {path}")
        cache_key = str(path.resolve())
        if cache_key not in truth_path_cache:
            truth_frame, truth_stat = read_truth_keys(path)
            truth_path_cache[cache_key] = (truth_frame, truth_stat, sha256_file(path))
        truth_frame, truth_stat, truth_hash = truth_path_cache[cache_key]
        truth_frames[sample] = truth_frame
        truth_stats[sample] = dict(truth_stat, sample=sample, truth_label=TRUTH_LABELS[sample])
        input_provenance["truth_vcfs"][sample] = {
            "path": str(path),
            "sha256": truth_hash,
            "label": TRUTH_LABELS[sample],
            "biallelic_autosomal_snv_count": int(len(truth_frame)),
        }

    frames: Dict[Tuple[str, str], pd.DataFrame] = {}
    query_qa: List[Dict[str, Any]] = []
    all_records: List[pd.DataFrame] = []
    summary_rows: List[Dict[str, Any]] = []
    histogram_output: List[Dict[str, Any]] = []
    band_output: List[Dict[str, Any]] = []
    truth_intersections: List[Dict[str, Any]] = []
    for source in SOURCE_ORDER:
        input_provenance["caller_vcfs"][source] = {}
        for sample in SAMPLE_ORDER:
            manifest_entry = samples[sample]
            vcf_path = Path(nested_get(manifest_entry, SOURCE_FIELDS[source])).resolve()
            if not vcf_path.exists():
                raise FileNotFoundError(f"Caller VCF unavailable for {sample}/{source}: {vcf_path}")
            frame, qa = read_vaf_vcf(vcf_path, sample, source)
            truth_keys = set(truth_frames[sample]["key"])
            frame["truth_confirmed"] = frame["key"].isin(truth_keys)
            frame["truth_label"] = TRUTH_LABELS[sample]
            frames[(source, sample)] = frame
            all_records.append(frame)
            query_qa.append(qa)
            input_provenance["caller_vcfs"][source][sample] = {
                "path": str(vcf_path),
                "sha256": sha256_file(vcf_path),
                "manifest_declared_sha256": nested_get(manifest_entry, SOURCE_FIELDS[source][:-1] + ("identity", "sha256")),
            }
            for scope, scoped_frame in (
                ("all_caller", frame),
                ("truth_confirmed", frame[frame["truth_confirmed"]]),
            ):
                summary_rows.append(dict(sample=sample, source=source, **distribution_summary(scoped_frame, scope)))
                histogram_output.extend(histogram_rows(scoped_frame, sample, source, scope))
                band_output.extend(vaf_band_rows(scoped_frame, sample, source, scope))
            intersection = int(frame["truth_confirmed"].sum())
            truth_count = int(len(truth_frames[sample]))
            truth_intersections.append(
                {
                    "sample": sample,
                    "source": source,
                    "truth_label": TRUTH_LABELS[sample],
                    "caller_biallelic_autosomal_snv_n": int(len(frame)),
                    "truth_biallelic_autosomal_snv_n": truth_count,
                    "exact_truth_intersection_n": intersection,
                    "exact_intersection_over_caller": float(intersection / len(frame)),
                    "exact_intersection_over_truth": float(intersection / truth_count) if truth_count else float("nan"),
                }
            )

    records = pd.concat(all_records, ignore_index=True).sort_values(["source", "sample", "chrom", "pos", "ref", "alt"])
    summary = pd.DataFrame(summary_rows).sort_values(["source", "sample", "scope"])
    histogram = pd.DataFrame(histogram_output).sort_values(["source", "sample", "scope", "bin_index"])
    bands = pd.DataFrame(band_output).sort_values(["source", "sample", "scope", "band"])
    truth_intersection_frame = pd.DataFrame(truth_intersections).sort_values(["source", "sample"])
    query_qa_frame = pd.DataFrame(query_qa).sort_values(["source", "sample"])
    truth_stats_frame = pd.DataFrame(truth_stats.values()).sort_values("sample")

    pair_unions: Dict[str, pd.DataFrame] = {}
    pair_metric_rows: List[Dict[str, Any]] = []
    pair_shared_rows: List[pd.DataFrame] = []
    bootstrap_rows: List[pd.DataFrame] = []
    shared_truth_keys = set(truth_frames["HCC1395"]["key"])
    for source_index, source in enumerate(SOURCE_ORDER):
        union = build_pair_union(frames[(source, "HCC1395")], frames[(source, "HCC1395_DORADO")])
        union["truth_confirmed"] = union["key"].isin(shared_truth_keys)
        pair_unions[source] = union
        shared = union[union["present_hcc1395"] & union["present_dorado"]].copy()
        shared["source"] = source
        shared["vaf_delta_dorado_minus_hcc1395"] = shared["af_dorado"] - shared["af_hcc1395"]
        shared["vaf_mean"] = (shared["af_dorado"] + shared["af_hcc1395"]) / 2.0
        pair_shared_rows.append(shared)
        for scope_index, (scope, scoped_union) in enumerate(
            [
                ("exact_caller_union", union),
                ("exact_caller_union_intersect_seqc2", union[union["truth_confirmed"]].copy()),
            ]
        ):
            values = paired_metric_values(scoped_union)
            intervals, boot = block_bootstrap_metrics(
                scoped_union,
                reps=args.bootstrap_reps,
                seed=args.seed + source_index * 100_000 + scope_index * 10_000,
            )
            bootstrap_frame = pd.DataFrame(boot)
            bootstrap_frame.insert(0, "scope", scope)
            bootstrap_frame.insert(0, "source", source)
            bootstrap_rows.append(bootstrap_frame)
            n_shared = int((scoped_union["present_hcc1395"] & scoped_union["present_dorado"]).sum())
            for metric, value in values.items():
                ci_low, ci_high, valid_reps = intervals.get(metric, (float("nan"), float("nan"), 0))
                pair_metric_rows.append(
                    {
                        "source": source,
                        "scope": scope,
                        "metric": metric,
                        "value": value,
                        "ci95_low_1mb_block_bootstrap": ci_low,
                        "ci95_high_1mb_block_bootstrap": ci_high,
                        "bootstrap_valid_reps": valid_reps,
                        "bootstrap_requested_reps": args.bootstrap_reps,
                        "n_union": int(len(scoped_union)),
                        "n_shared": n_shared,
                        "n_hcc1395_only": int((scoped_union["present_hcc1395"] & ~scoped_union["present_dorado"]).sum()),
                        "n_dorado_only": int((~scoped_union["present_hcc1395"] & scoped_union["present_dorado"]).sum()),
                        "n_1mb_blocks": int(scoped_union["block_1mb"].nunique()),
                    }
                )

    pair_metrics = pd.DataFrame(pair_metric_rows).sort_values(["source", "scope", "metric"])
    pair_shared = pd.concat(pair_shared_rows, ignore_index=True).sort_values(["source", "chrom", "pos", "ref", "alt"])
    bootstrap_frame = pd.concat(bootstrap_rows, ignore_index=True).sort_values(["source", "scope", "bootstrap_replicate"])

    write_tsv_gz(records, data_dir / "raw_vaf_records.tsv.gz")
    write_tsv(summary, data_dir / "raw_vaf_summary.tsv")
    write_tsv(histogram, data_dir / "raw_vaf_histogram_0.02.tsv")
    write_tsv(bands, data_dir / "raw_vaf_band_composition.tsv")
    write_tsv(truth_intersection_frame, data_dir / "truth_intersection_stats.tsv")
    write_tsv(query_qa_frame, data_dir / "bcftools_query_count_reconstruction_qa.tsv")
    write_tsv(truth_stats_frame, data_dir / "truth_vcf_inventory.tsv")
    write_tsv(pair_metrics, data_dir / "hcc1395_pair_metrics.tsv")
    write_tsv_gz(pair_shared, data_dir / "hcc1395_pair_exact_shared_records.tsv.gz")
    write_tsv_gz(bootstrap_frame, data_dir / "hcc1395_pair_1mb_block_bootstrap_replicates.tsv.gz")

    plot_truth_distributions(histogram, figure_dir)
    plot_band_composition(bands, figure_dir)
    plot_hcc_hexbin(pair_unions, figure_dir)
    plot_bland_altman(pair_unions, figure_dir)

    bcftools_version = run_checked(["bcftools", "--version"]).splitlines()[:2]
    output_hashes = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "provenance_and_claim_ceiling.json":
            output_hashes[str(path.relative_to(output_dir))] = {
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
    provenance = {
        "schema": "intersubmod.raw_vaf_technical_reproducibility.v1",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "working_directory": os.getcwd(),
        "scope": {
            "task_type": "comprehensive_validation",
            "contigs": AUTOSOMES,
            "datasets": SAMPLE_ORDER,
            "sources": SOURCE_ORDER,
            "variant_contract": "FILTER=PASS, biallelic SNV, exact CHROM:POS:REF:ALT",
        },
        "count_contract": {
            "alt_count": "round_half_up(DP * AF) = floor(DP * AF + 0.5)",
            "ref_count": "DP - alt_count",
            "ad_policy": "FORMAT/AD retained only for ALT-depth QA; AD[0] is never read or used",
        },
        "hcc_pair_contract": {
            "biological_relationship": "same HCC1395 cell line; different sequencing/basecalling/processing datasets",
            "join": "exact CHROM:POS:REF:ALT",
            "truth_subset": "shared exact alleles intersected with the same SEQC2 v1.2.1 HC sSNV key set",
            "bootstrap": {
                "unit": "1 Mb genomic block (chromosome + floor((POS-1)/1,000,000))",
                "method": "nonparametric block resampling with replacement",
                "repetitions": args.bootstrap_reps,
                "seed": args.seed,
                "interval": "2.5th and 97.5th percentiles",
            },
            "binomial_noise_ratio": "observed paired delta MSE divided by mean pooled-binomial sampling variance; RMSE ratio is its square root",
        },
        "claim_ceiling": {
            "level": "raw_vaf_technical_reproducibility_only",
            "allowed": [
                "describe raw caller-VAF distributions",
                "quantify exact-locus technical concordance between HCC1395 and HCC1395_DORADO",
                "identify whether observed paired differences exceed a simple independent-binomial sampling model",
            ],
            "not_allowed": [
                "infer clone or subclone truth directly from raw VAF",
                "claim a unique or accurate tumor phylogeny",
                "treat the two HCC1395 datasets as independent biological replicates",
                "attribute discordance to biological evolution without CN, purity, multiplicity, and orthogonal evidence",
            ],
            "confounders": [
                "copy number and LOH",
                "tumor purity",
                "mutant allele multiplicity",
                "coverage and caller/filter selection",
                "basecalling and alignment differences",
                "truth-set scope differences across cell lines",
            ],
        },
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "bcftools": bcftools_version,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "inputs": input_provenance,
        "qa": {
            "query_rows": query_qa,
            "truth_inventory": list(truth_stats.values()),
        },
        "outputs": output_hashes,
    }
    provenance_path = output_dir / "provenance_and_claim_ceiling.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, ensure_ascii=False, allow_nan=False) + "\n")

    latest_truth = truth_intersection_frame[truth_intersection_frame["source"] == "latest_lps_pass"]
    latest_metrics = pair_metrics[(pair_metrics["source"] == "latest_lps_pass") & (pair_metrics["scope"] == "exact_caller_union")]
    metric_lookup = dict(zip(latest_metrics["metric"], latest_metrics["value"]))
    print(f"PASS datasets={len(SAMPLE_ORDER)} sources={len(SOURCE_ORDER)} valid_records={len(records):,}")
    print(
        "PASS count_reconstruction="
        f"{int(query_qa_frame['count_reconstruction_ok_n'].sum()):,}/"
        f"{int(query_qa_frame['valid_unique_rows'].sum()):,}"
    )
    print(
        "LATEST truth intersections="
        + ", ".join(
            f"{row.sample}:{int(row.exact_truth_intersection_n):,}/{int(row.caller_biallelic_autosomal_snv_n):,}"
            for row in latest_truth.itertuples()
        )
    )
    print(
        "LATEST HCC pair exact metrics "
        f"Jaccard={metric_lookup['callset_jaccard']:.4f} "
        f"Pearson={metric_lookup['vaf_pearson']:.4f} "
        f"Spearman={metric_lookup['vaf_spearman']:.4f} "
        f"CCC={metric_lookup['vaf_ccc']:.4f} "
        f"MAE={metric_lookup['vaf_mae']:.4f}"
    )
    print(f"OUTPUT {output_dir}")


if __name__ == "__main__":
    main()
