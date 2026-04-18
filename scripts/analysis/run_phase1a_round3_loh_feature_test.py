#!/usr/bin/env python3
"""Phase 1A Round 3: LOH feature incremental test.

4-way comparison on the existing paired-multibio sample637 shard:
  1. context_only       — baseline (no methyl, no LOH)
  2. methyl+context     — Round 2 best (methyl features added)
  3. loh+context        — LOH features added (Potential_LOH, LOH_Subtype)
  4. methyl+loh+context — full feature set

Key hypothesis: LOH features will specifically rescue H2009's negative
delta F1 (-0.0039 in Round 2), since 68% of H2009 FP regions are in LOH.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

from research_common import ensure_dir, write_json, write_tsv_rows
from run_phase1a_read_classifier_benchmark import (
    CATEGORICAL_FEATURES_CONTEXT,
    NUMERIC_FEATURES_CONTEXT,
    NUMERIC_FEATURES_CONTEXT_ONLY,
    build_logistic_pipeline,
    compute_metrics,
    load_benchmark_rows,
    pick_discovery_holdout,
)


DEFAULT_INPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_manifest_shard_export_paired_multibio_sample637/"
    "phase1_shard_read_training_table.tsv"
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260328_phase1a_round3_loh_feature_v1"
)

# LOH categorical features already present in shard
LOH_CATEGORICAL_FEATURES = ["Potential_LOH", "LOH_Subtype"]

# 4 model definitions: (name, numeric_features, categorical_features)
MODEL_CONFIGS: List[Tuple[str, List[str], List[str]]] = [
    ("context_only", list(NUMERIC_FEATURES_CONTEXT_ONLY), list(CATEGORICAL_FEATURES_CONTEXT)),
    ("methyl+context", list(NUMERIC_FEATURES_CONTEXT), list(CATEGORICAL_FEATURES_CONTEXT)),
    ("loh+context", list(NUMERIC_FEATURES_CONTEXT_ONLY), list(CATEGORICAL_FEATURES_CONTEXT) + LOH_CATEGORICAL_FEATURES),
    ("methyl+loh+context", list(NUMERIC_FEATURES_CONTEXT), list(CATEGORICAL_FEATURES_CONTEXT) + LOH_CATEGORICAL_FEATURES),
]

RESULT_FIELDS = [
    "model_name",
    "evaluation_split",
    "train_rows",
    "test_rows",
    "train_regions",
    "test_regions",
    "F1",
    "precision",
    "recall",
    "accuracy",
    "balanced_accuracy",
    "delta_F1_vs_context_only",
    "bootstrap_delta_f1_mean",
    "bootstrap_delta_f1_ci_low",
    "bootstrap_delta_f1_ci_high",
    "bootstrap_delta_f1_positive_fraction",
    "incremental_supported",
    "mcnemar_pvalue_vs_context",
]

DATASET_FIELDS = [
    "model_name",
    "evaluation_split",
    "dataset_id",
    "dataset_label",
    "sample",
    "harmonization_group",
    "rows_total",
    "regions_total",
    "loh_fp_region_pct",
    "F1",
    "delta_F1_vs_context_only",
    "bootstrap_delta_f1_ci_low",
    "bootstrap_delta_f1_ci_high",
    "bootstrap_delta_f1_positive_fraction",
    "incremental_supported",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tsv", default=str(DEFAULT_INPUT))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--test-size", type=float, default=0.2)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--bootstrap-reps", type=int, default=500)
    return parser.parse_args()


def prepare_loh_features(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure LOH categorical features are string-encoded for OneHotEncoder."""
    df = df.copy()
    # Potential_LOH: convert to string "True"/"False" for consistent encoding
    df["Potential_LOH"] = df["Potential_LOH"].astype(str)
    # LOH_Subtype: fill NaN, ensure string
    df["LOH_Subtype"] = df["LOH_Subtype"].fillna("None").astype(str)
    return df


def bootstrap_delta(
    test_df: pd.DataFrame,
    baseline_pred: np.ndarray,
    model_pred: np.ndarray,
    reps: int,
    seed: int,
) -> Dict[str, float]:
    y_true = test_df["label"].to_numpy(dtype=int)
    region_keys = test_df["region_key"].astype(str).to_numpy()
    unique_regions = pd.unique(region_keys)
    region_to_indices = {
        rk: np.flatnonzero(region_keys == rk) for rk in unique_regions
    }
    rng = np.random.default_rng(seed)
    delta_f1s: List[float] = []
    for _ in range(reps):
        sampled = rng.choice(unique_regions, size=len(unique_regions), replace=True)
        idx = np.concatenate([region_to_indices[rk] for rk in sampled])
        f1_base = compute_metrics(y_true[idx], baseline_pred[idx])["F1"]
        f1_model = compute_metrics(y_true[idx], model_pred[idx])["F1"]
        delta_f1s.append(float(f1_model) - float(f1_base))
    arr = np.asarray(delta_f1s)
    return {
        "bootstrap_delta_f1_mean": float(np.mean(arr)),
        "bootstrap_delta_f1_ci_low": float(np.quantile(arr, 0.025)),
        "bootstrap_delta_f1_ci_high": float(np.quantile(arr, 0.975)),
        "bootstrap_delta_f1_positive_fraction": float(np.mean(arr > 0)),
    }


def mcnemar_pvalue(correct_a: np.ndarray, correct_b: np.ndarray) -> float:
    a_only = int(np.logical_and(correct_a, ~correct_b).sum())
    b_only = int(np.logical_and(~correct_a, correct_b).sum())
    total = a_only + b_only
    if total == 0:
        return 1.0
    chi2 = ((abs(a_only - b_only) - 1.0) ** 2) / float(total)
    return float(math.erfc(math.sqrt(chi2 / 2.0)))


def loh_fp_region_pct(region_df: pd.DataFrame) -> float:
    """Fraction of FP regions that have Potential_LOH=True."""
    fp_regions = region_df[region_df["truth_status"] == "FP"]
    if fp_regions.empty:
        return 0.0
    return float((fp_regions["Potential_LOH"].astype(str) == "True").mean())


def run_all_models(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    random_state: int,
    bootstrap_reps: int,
    evaluation_split: str,
) -> Tuple[List[Dict], List[Dict], Dict[str, np.ndarray]]:
    """Fit all 4 models and return result rows, dataset rows, and predictions dict."""
    predictions: Dict[str, np.ndarray] = {}

    for model_name, numeric_feats, categorical_feats in MODEL_CONFIGS:
        features = numeric_feats + categorical_feats
        pipeline = build_logistic_pipeline(numeric_feats, categorical_feats)
        pipeline.fit(train_df[features], train_df["label"])
        predictions[model_name] = pipeline.predict(test_df[features]).astype(int)

    baseline_pred = predictions["context_only"]
    y_true = test_df["label"].to_numpy(dtype=int)
    baseline_correct = baseline_pred == y_true

    result_rows: List[Dict] = []
    for model_name, _, _ in MODEL_CONFIGS:
        pred = predictions[model_name]
        metrics = compute_metrics(y_true, pred)
        baseline_f1 = compute_metrics(y_true, baseline_pred)["F1"]
        delta_f1 = float(metrics["F1"]) - float(baseline_f1)

        boot = bootstrap_delta(
            test_df=test_df,
            baseline_pred=baseline_pred,
            model_pred=pred,
            reps=bootstrap_reps,
            seed=random_state,
        )
        pred_correct = pred == y_true
        pvalue = mcnemar_pvalue(baseline_correct, pred_correct)

        result_rows.append({
            "model_name": model_name,
            "evaluation_split": evaluation_split,
            "train_rows": int(len(train_df)),
            "test_rows": int(len(test_df)),
            "train_regions": int(train_df["region_key"].nunique()),
            "test_regions": int(test_df["region_key"].nunique()),
            "F1": float(metrics["F1"]),
            "precision": float(metrics["precision"]),
            "recall": float(metrics["recall"]),
            "accuracy": float(metrics["accuracy"]),
            "balanced_accuracy": float(metrics["balanced_accuracy"]),
            "delta_F1_vs_context_only": delta_f1,
            **boot,
            "incremental_supported": bool(
                not math.isnan(boot["bootstrap_delta_f1_ci_low"])
                and boot["bootstrap_delta_f1_ci_low"] > 0
            ),
            "mcnemar_pvalue_vs_context": pvalue,
        })

    # Per-dataset breakdown
    dataset_rows: List[Dict] = []
    region_df = test_df.drop_duplicates("region_key")

    for model_name, _, _ in MODEL_CONFIGS:
        pred = predictions[model_name]
        working = test_df.copy()
        working["pred"] = pred
        baseline_working = test_df.copy()
        baseline_working["pred"] = baseline_pred

        grouped = working.groupby(["dataset_id", "dataset_label", "sample", "harmonization_group"], sort=True, dropna=False)
        for (dataset_id, dataset_label, sample, harmonization_group), sub in grouped:
            sub_pred = sub["pred"].to_numpy(dtype=int)
            sub_base = baseline_working.loc[sub.index, "pred"].to_numpy(dtype=int)
            sub_true = sub["label"].to_numpy(dtype=int)
            metrics = compute_metrics(sub_true, sub_pred)
            base_f1 = compute_metrics(sub_true, sub_base)["F1"]
            delta = float(metrics["F1"]) - float(base_f1)

            sub_regions = region_df[region_df["dataset_id"] == dataset_id]
            fp_loh_pct = loh_fp_region_pct(sub_regions)

            boot = bootstrap_delta(
                test_df=sub,
                baseline_pred=sub_base,
                model_pred=sub_pred,
                reps=bootstrap_reps,
                seed=random_state,
            )
            dataset_rows.append({
                "model_name": model_name,
                "evaluation_split": evaluation_split,
                "dataset_id": dataset_id,
                "dataset_label": dataset_label,
                "sample": sample,
                "harmonization_group": harmonization_group,
                "rows_total": int(len(sub)),
                "regions_total": int(sub["region_key"].nunique()),
                "loh_fp_region_pct": round(fp_loh_pct, 4),
                "F1": float(metrics["F1"]),
                "delta_F1_vs_context_only": delta,
                **{k: v for k, v in boot.items() if "ci" in k or "positive" in k},
                "incremental_supported": bool(
                    not math.isnan(boot["bootstrap_delta_f1_ci_low"])
                    and boot["bootstrap_delta_f1_ci_low"] > 0
                ),
            })

    return result_rows, dataset_rows, predictions


def make_overall_f1_figure(result_rows: List[Dict], output_dir: Path) -> None:
    """Bar chart: F1 per model per evaluation split."""
    df = pd.DataFrame(result_rows)
    splits = df["evaluation_split"].unique()
    models = [m for m, _, _ in MODEL_CONFIGS]
    colors = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759"]

    fig, axes = plt.subplots(1, len(splits), figsize=(8 * len(splits), 6))
    if len(splits) == 1:
        axes = [axes]

    for ax, split in zip(axes, splits):
        sub = df[df["evaluation_split"] == split]
        f1s = [sub[sub["model_name"] == m]["F1"].values[0] if len(sub[sub["model_name"] == m]) > 0 else 0.0 for m in models]
        deltas = [sub[sub["model_name"] == m]["delta_F1_vs_context_only"].values[0] if len(sub[sub["model_name"] == m]) > 0 else 0.0 for m in models]
        bars = ax.bar(models, f1s, color=colors, edgecolor="black", linewidth=0.7)
        for bar, f1, delta in zip(bars, f1s, deltas):
            sign = "+" if delta > 0 else ""
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.002,
                    f"{f1:.4f}\n({sign}{delta:.4f})", ha="center", va="bottom", fontsize=8.5)
        ax.set_ylim(max(0.0, min(f1s) - 0.05), min(1.0, max(f1s) + 0.04))
        ax.set_title(f"{split}", fontsize=12, fontweight="bold")
        ax.set_ylabel("F1 Score")
        ax.set_xticklabels(models, rotation=20, ha="right")
        ax.axhline(f1s[0], color="gray", linestyle="--", linewidth=0.8, label="baseline (context_only)")

    fig.suptitle("Phase 1A Round 3: 4-Model F1 Comparison\n(delta vs context_only in parentheses)", fontsize=13)
    fig.tight_layout()
    fig.savefig(output_dir / "fig01_overall_f1.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[round3] wrote fig01_overall_f1.png")


def make_per_dataset_delta_figure(dataset_rows: List[Dict], output_dir: Path) -> None:
    """Per-dataset delta F1 vs context_only for all 3 non-baseline models."""
    df = pd.DataFrame(dataset_rows)
    val_df = df[df["evaluation_split"] == "external_validation"].copy()
    if val_df.empty:
        return

    models_to_plot = ["methyl+context", "loh+context", "methyl+loh+context"]
    colors = {"methyl+context": "#f28e2b", "loh+context": "#59a14f", "methyl+loh+context": "#e15759"}

    samples = sorted(val_df["sample"].unique())
    x = np.arange(len(samples))
    width = 0.25

    fig, ax = plt.subplots(figsize=(12, 6))
    for i, model in enumerate(models_to_plot):
        mdf = val_df[val_df["model_name"] == model].copy()
        mdf_indexed = mdf.set_index("sample")
        deltas = [float(mdf_indexed.loc[s, "delta_F1_vs_context_only"]) if s in mdf_indexed.index else 0.0 for s in samples]
        ci_lows = [float(mdf_indexed.loc[s, "bootstrap_delta_f1_ci_low"]) if s in mdf_indexed.index else 0.0 for s in samples]
        ci_highs = [float(mdf_indexed.loc[s, "bootstrap_delta_f1_ci_high"]) if s in mdf_indexed.index else 0.0 for s in samples]
        err_low = [d - cl for d, cl in zip(deltas, ci_lows)]
        err_high = [ch - d for d, ch in zip(deltas, ci_highs)]
        ax.bar(x + (i - 1) * width, deltas,
               width=width * 0.9, color=colors[model], edgecolor="black",
               linewidth=0.5, label=model)
        ax.errorbar(x + (i - 1) * width, deltas,
                    yerr=[err_low, err_high], fmt="none",
                    color="black", capsize=3, linewidth=0.8)

    # Annotate LOH FP pct for each sample
    loh_df = val_df[val_df["model_name"] == "loh+context"].set_index("sample")
    for j, s in enumerate(samples):
        if s in loh_df.index:
            pct = float(loh_df.loc[s, "loh_fp_region_pct"])
            ax.text(x[j], ax.get_ylim()[0] - 0.005, f"LOH-FP:{pct:.0%}", ha="center", va="top", fontsize=7.5, color="#666666")

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=20, ha="right")
    ax.set_ylabel("delta F1 vs context_only")
    ax.set_title("Per-Dataset Delta F1 (external validation)\n[LOH-FP%=fraction of FP regions with LOH]", fontsize=12)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(output_dir / "fig02_per_dataset_delta_f1.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[round3] wrote fig02_per_dataset_delta_f1.png")


def make_loh_distribution_figure(test_df: pd.DataFrame, output_dir: Path) -> None:
    """LOH region distribution by sample and truth_status."""
    rdf = test_df.drop_duplicates("region_key").copy()
    rdf["Potential_LOH_str"] = rdf["Potential_LOH"].astype(str)

    samples = sorted(rdf["sample"].unique())
    fig, axes = plt.subplots(1, len(samples), figsize=(3 * len(samples), 5), sharey=False)
    if len(samples) == 1:
        axes = [axes]

    colors = {"TP": "#4e79a7", "FP": "#e15759"}
    for ax, sample in zip(axes, samples):
        sdf = rdf[rdf["sample"] == sample]
        for status in ["TP", "FP"]:
            sub = sdf[sdf["truth_status"] == status]
            if sub.empty:
                continue
            loh_true = int((sub["Potential_LOH_str"] == "True").sum())
            total = len(sub)
            loh_false = total - loh_true
            ax.bar([f"{status}\nLOH=N", f"{status}\nLOH=Y"],
                   [loh_false, loh_true],
                   color=[colors[status] + "88", colors[status]],
                   edgecolor="black", linewidth=0.5)
            ax.text(1, loh_true + 0.2, f"{loh_true}/{total}\n({loh_true/total:.0%})",
                    ha="center", va="bottom", fontsize=7.5)
        ax.set_title(sample, fontsize=9, fontweight="bold")
        ax.set_xlabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", fontsize=7)

    fig.suptitle("LOH Region Distribution per Sample (region-level)\nblue=TP, red=FP", fontsize=11)
    fig.tight_layout()
    fig.savefig(output_dir / "fig03_loh_region_distribution.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[round3] wrote fig03_loh_region_distribution.png")


def make_bucket_figure(
    test_df: pd.DataFrame,
    predictions: Dict[str, np.ndarray],
    output_dir: Path,
    focus_samples: List[str] = ["H2009", "HCC1937", "HCC1954"],
) -> None:
    """Error rate per bucket for focus samples: context_only vs loh+context vs methyl+loh+context."""
    working = test_df.copy()
    for model_name, pred in predictions.items():
        working[f"pred_{model_name}"] = pred
    working["y_true"] = working["label"].to_numpy(dtype=int)

    fig_rows: List[Dict] = []
    for sample in focus_samples:
        sub = working[working["sample"] == sample]
        if sub.empty:
            continue
        for bucket_name, bucket_mask in [
            ("TP+PassedGating", (sub["truth_status"] == "TP") & (sub["PassedGating"].astype(str) == "True")),
            ("TP+NotGated", (sub["truth_status"] == "TP") & (sub["PassedGating"].astype(str) == "False")),
            ("FP+LOH=Y", (sub["truth_status"] == "FP") & (sub["Potential_LOH"].astype(str) == "True")),
            ("FP+LOH=N", (sub["truth_status"] == "FP") & (sub["Potential_LOH"].astype(str) == "False")),
        ]:
            bucket_sub = sub[bucket_mask]
            if len(bucket_sub) == 0:
                continue
            for model_name in ["context_only", "loh+context", "methyl+loh+context"]:
                if f"pred_{model_name}" not in bucket_sub.columns:
                    continue
                pred = bucket_sub[f"pred_{model_name}"].to_numpy(dtype=int)
                y = bucket_sub["y_true"].to_numpy(dtype=int)
                err_rate = float(np.mean(pred != y))
                fig_rows.append({
                    "sample": sample,
                    "bucket": bucket_name,
                    "model": model_name,
                    "error_rate": err_rate,
                    "n": len(bucket_sub),
                })

    if not fig_rows:
        return

    fig_df = pd.DataFrame(fig_rows)
    bucket_order = ["TP+PassedGating", "TP+NotGated", "FP+LOH=Y", "FP+LOH=N"]
    model_colors = {"context_only": "#4e79a7", "loh+context": "#59a14f", "methyl+loh+context": "#e15759"}
    model_labels = list(model_colors.keys())
    bar_width = 0.25

    fig, axes = plt.subplots(1, len(focus_samples), figsize=(7 * len(focus_samples), 6), sharey=False)
    if len(focus_samples) == 1:
        axes = [axes]

    for ax, sample in zip(axes, focus_samples):
        sub_df = fig_df[fig_df["sample"] == sample]
        if sub_df.empty:
            ax.set_title(f"{sample}\n(no data)")
            continue
        x = np.arange(len(bucket_order))
        for i, model in enumerate(model_labels):
            mdf = sub_df[sub_df["model"] == model].set_index("bucket")
            vals = [float(mdf.loc[b, "error_rate"]) if b in mdf.index else 0.0 for b in bucket_order]
            ax.bar(x + (i - 1) * bar_width, vals, width=bar_width * 0.9,
                   color=model_colors[model], edgecolor="black", linewidth=0.5,
                   label=model if sample == focus_samples[0] else "")
        ax.set_xticks(x)
        ax.set_xticklabels(bucket_order, rotation=25, ha="right", fontsize=8.5)
        ax.set_ylabel("Error Rate")
        ax.set_title(f"{sample}", fontsize=11, fontweight="bold")
        ax.set_ylim(0, 1.0)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", fontsize=9)
    fig.suptitle("Bucket Error Rate: context_only vs loh+context vs methyl+loh+context\n(focus: H2009, HCC1937, HCC1954)", fontsize=11)
    fig.tight_layout()
    fig.savefig(output_dir / "fig04_bucket_error_rates.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("[round3] wrote fig04_bucket_error_rates.png")


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_tsv).resolve()
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    print(f"[round3] loading {input_path}")
    df = load_benchmark_rows(input_path)
    df = prepare_loh_features(df)

    discovery_df = df[df["dataset_role"] == "discovery"].copy()
    validation_df = df[df["dataset_role"] == "validation"].copy()

    print(f"[round3] discovery rows={len(discovery_df)} | validation rows={len(validation_df)}")

    discovery_train_df, discovery_holdout_df = pick_discovery_holdout(
        discovery_df, args.test_size, args.random_state
    )

    all_result_rows: List[Dict] = []
    all_dataset_rows: List[Dict] = []
    all_predictions: Dict[str, np.ndarray] = {}

    # Discovery holdout
    print("[round3] running discovery holdout...")
    r_rows, d_rows, preds = run_all_models(
        discovery_train_df, discovery_holdout_df,
        args.random_state, args.bootstrap_reps, "discovery_holdout"
    )
    all_result_rows.extend(r_rows)
    all_dataset_rows.extend(d_rows)
    all_predictions["discovery_holdout"] = preds

    # External validation
    if not validation_df.empty:
        print("[round3] running external validation...")
        r_rows, d_rows, preds = run_all_models(
            discovery_df, validation_df,
            args.random_state, args.bootstrap_reps, "external_validation"
        )
        all_result_rows.extend(r_rows)
        all_dataset_rows.extend(d_rows)
        all_predictions["external_validation"] = preds

    # Write TSVs
    write_tsv_rows(output_dir / "round3_model_comparison.tsv", RESULT_FIELDS, all_result_rows)
    write_tsv_rows(output_dir / "round3_dataset_comparison.tsv", DATASET_FIELDS, all_dataset_rows)
    write_json(output_dir / "run_context.json", {
        "task": "Phase 1A Round 3 LOH feature test",
        "input_tsv": str(input_path),
        "random_state": args.random_state,
        "test_size": args.test_size,
        "bootstrap_reps": args.bootstrap_reps,
        "models": [m for m, _, _ in MODEL_CONFIGS],
        "loh_features": LOH_CATEGORICAL_FEATURES,
    })

    # Figures
    print("[round3] generating figures...")
    make_overall_f1_figure(all_result_rows, output_dir)
    make_per_dataset_delta_figure(all_dataset_rows, output_dir)

    # Use validation data for LOH distribution and bucket figures
    eval_data = validation_df if not validation_df.empty else discovery_holdout_df
    make_loh_distribution_figure(eval_data, output_dir)
    if not validation_df.empty:
        make_bucket_figure(validation_df, all_predictions["external_validation"], output_dir)

    # Print summary
    print("\n=== Phase 1A Round 3 Summary ===")
    result_df = pd.DataFrame(all_result_rows)
    for split in result_df["evaluation_split"].unique():
        sub = result_df[result_df["evaluation_split"] == split]
        print(f"\n{split}:")
        for _, row in sub.iterrows():
            supported = "✓" if row["incremental_supported"] else " "
            sign = "+" if row["delta_F1_vs_context_only"] > 0 else ""
            print(f"  {row['model_name']:25s}  F1={row['F1']:.4f}  "
                  f"delta={sign}{row['delta_F1_vs_context_only']:.4f}  "
                  f"CI=[{row['bootstrap_delta_f1_ci_low']:+.4f}, {row['bootstrap_delta_f1_ci_high']:+.4f}]  "
                  f"{supported}")

    print(f"\n[round3] outputs written to {output_dir}")


if __name__ == "__main__":
    main()
