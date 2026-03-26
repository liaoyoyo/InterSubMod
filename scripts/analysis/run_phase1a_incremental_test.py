#!/usr/bin/env python3
"""Quantify Phase 1A methyl incremental gain over context-only baselines."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np
import pandas as pd

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
    "20260325_phase1_manifest_shard_export_sample200/phase1_shard_read_training_table.tsv"
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_incremental_test_sample200_v1"
)


COMPARISON_FIELDS = [
    "evaluation_split",
    "random_state",
    "bootstrap_reps",
    "train_rows",
    "test_rows",
    "train_regions",
    "test_regions",
    "context_only_F1",
    "methyl_context_F1",
    "delta_F1",
    "context_only_accuracy",
    "methyl_context_accuracy",
    "delta_accuracy",
    "context_only_balanced_accuracy",
    "methyl_context_balanced_accuracy",
    "delta_balanced_accuracy",
    "context_only_correct_methyl_context_wrong",
    "context_only_wrong_methyl_context_correct",
    "discordant_total",
    "mcnemar_pvalue",
    "bootstrap_delta_f1_mean",
    "bootstrap_delta_f1_ci_low",
    "bootstrap_delta_f1_ci_high",
    "bootstrap_delta_f1_positive_fraction",
    "bootstrap_delta_accuracy_mean",
    "bootstrap_delta_accuracy_ci_low",
    "bootstrap_delta_accuracy_ci_high",
    "bootstrap_delta_accuracy_positive_fraction",
    "incremental_supported",
]


DATASET_FIELDS = [
    "evaluation_split",
    "random_state",
    "dataset_id",
    "dataset_label",
    "harmonization_group",
    "rows_total",
    "regions_total",
    "context_only_F1",
    "methyl_context_F1",
    "delta_F1",
    "context_only_accuracy",
    "methyl_context_accuracy",
    "delta_accuracy",
    "context_only_correct_methyl_context_wrong",
    "context_only_wrong_methyl_context_correct",
    "discordant_total",
    "mcnemar_pvalue",
    "bootstrap_delta_f1_mean",
    "bootstrap_delta_f1_ci_low",
    "bootstrap_delta_f1_ci_high",
    "bootstrap_delta_f1_positive_fraction",
    "incremental_supported",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tsv", default=str(DEFAULT_INPUT), help="Read-level Phase 1A shard TSV.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="Output directory.")
    parser.add_argument("--test-size", type=float, default=0.2, help="Discovery region holdout fraction.")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed for discovery holdout.")
    parser.add_argument(
        "--bootstrap-reps",
        type=int,
        default=500,
        help="Region-bootstrap repetitions for incremental CI estimates.",
    )
    return parser.parse_args()


def fit_context_models(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    context_only_features = list(NUMERIC_FEATURES_CONTEXT_ONLY) + list(CATEGORICAL_FEATURES_CONTEXT)
    context_features = list(NUMERIC_FEATURES_CONTEXT) + list(CATEGORICAL_FEATURES_CONTEXT)

    logistic_context_only = build_logistic_pipeline(NUMERIC_FEATURES_CONTEXT_ONLY, CATEGORICAL_FEATURES_CONTEXT)
    logistic_context_only.fit(train_df[context_only_features], train_df["label"])
    context_only_pred = logistic_context_only.predict(test_df[context_only_features]).astype(int)

    logistic_methyl_context = build_logistic_pipeline(NUMERIC_FEATURES_CONTEXT, CATEGORICAL_FEATURES_CONTEXT)
    logistic_methyl_context.fit(train_df[context_features], train_df["label"])
    methyl_context_pred = logistic_methyl_context.predict(test_df[context_features]).astype(int)

    return context_only_pred, methyl_context_pred


def mcnemar_pvalue(correct_a: np.ndarray, correct_b: np.ndarray) -> tuple[int, int, int, float]:
    a_only = int(np.logical_and(correct_a, np.logical_not(correct_b)).sum())
    b_only = int(np.logical_and(np.logical_not(correct_a), correct_b).sum())
    discordant_total = a_only + b_only
    if discordant_total == 0:
        return a_only, b_only, discordant_total, 1.0

    chi2 = ((abs(a_only - b_only) - 1.0) ** 2) / float(discordant_total)
    pvalue = math.erfc(math.sqrt(chi2 / 2.0))
    return a_only, b_only, discordant_total, float(pvalue)


def percentile(values: Sequence[float], q: float) -> float:
    if not values:
        return float("nan")
    return float(np.quantile(np.asarray(values, dtype=float), q))


def bootstrap_delta_metrics(
    test_df: pd.DataFrame,
    context_only_pred: np.ndarray,
    methyl_context_pred: np.ndarray,
    reps: int,
    seed: int,
) -> Dict[str, float]:
    if reps <= 0 or test_df.empty:
        return {
            "bootstrap_delta_f1_mean": float("nan"),
            "bootstrap_delta_f1_ci_low": float("nan"),
            "bootstrap_delta_f1_ci_high": float("nan"),
            "bootstrap_delta_f1_positive_fraction": float("nan"),
            "bootstrap_delta_accuracy_mean": float("nan"),
            "bootstrap_delta_accuracy_ci_low": float("nan"),
            "bootstrap_delta_accuracy_ci_high": float("nan"),
            "bootstrap_delta_accuracy_positive_fraction": float("nan"),
        }

    y_true = test_df["label"].to_numpy(dtype=int)
    region_keys = test_df["region_key"].astype(str).to_numpy()
    unique_regions = pd.unique(region_keys)
    region_to_indices = {
        region_key: np.flatnonzero(region_keys == region_key) for region_key in unique_regions
    }
    rng = np.random.default_rng(seed)

    delta_f1_values: List[float] = []
    delta_accuracy_values: List[float] = []
    for _ in range(reps):
        sampled_regions = rng.choice(unique_regions, size=len(unique_regions), replace=True)
        sampled_indices = np.concatenate([region_to_indices[region_key] for region_key in sampled_regions])
        context_metrics = compute_metrics(y_true[sampled_indices], context_only_pred[sampled_indices])
        methyl_metrics = compute_metrics(y_true[sampled_indices], methyl_context_pred[sampled_indices])
        delta_f1_values.append(float(methyl_metrics["F1"]) - float(context_metrics["F1"]))
        delta_accuracy_values.append(float(methyl_metrics["accuracy"]) - float(context_metrics["accuracy"]))

    return {
        "bootstrap_delta_f1_mean": float(np.mean(delta_f1_values)),
        "bootstrap_delta_f1_ci_low": percentile(delta_f1_values, 0.025),
        "bootstrap_delta_f1_ci_high": percentile(delta_f1_values, 0.975),
        "bootstrap_delta_f1_positive_fraction": float(np.mean(np.asarray(delta_f1_values) > 0)),
        "bootstrap_delta_accuracy_mean": float(np.mean(delta_accuracy_values)),
        "bootstrap_delta_accuracy_ci_low": percentile(delta_accuracy_values, 0.025),
        "bootstrap_delta_accuracy_ci_high": percentile(delta_accuracy_values, 0.975),
        "bootstrap_delta_accuracy_positive_fraction": float(np.mean(np.asarray(delta_accuracy_values) > 0)),
    }


def comparison_row(
    evaluation_split: str,
    random_state: int,
    bootstrap_reps: int,
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    context_only_pred: np.ndarray,
    methyl_context_pred: np.ndarray,
) -> Dict[str, object]:
    context_metrics = compute_metrics(test_df["label"], context_only_pred)
    methyl_metrics = compute_metrics(test_df["label"], methyl_context_pred)
    context_correct = context_only_pred == test_df["label"].to_numpy(dtype=int)
    methyl_correct = methyl_context_pred == test_df["label"].to_numpy(dtype=int)
    a_only, b_only, discordant_total, pvalue = mcnemar_pvalue(context_correct, methyl_correct)
    bootstrap = bootstrap_delta_metrics(
        test_df=test_df,
        context_only_pred=context_only_pred,
        methyl_context_pred=methyl_context_pred,
        reps=bootstrap_reps,
        seed=random_state,
    )

    return {
        "evaluation_split": evaluation_split,
        "random_state": random_state,
        "bootstrap_reps": bootstrap_reps,
        "train_rows": int(len(train_df.index)),
        "test_rows": int(len(test_df.index)),
        "train_regions": int(train_df["region_key"].nunique()),
        "test_regions": int(test_df["region_key"].nunique()),
        "context_only_F1": float(context_metrics["F1"]),
        "methyl_context_F1": float(methyl_metrics["F1"]),
        "delta_F1": float(methyl_metrics["F1"]) - float(context_metrics["F1"]),
        "context_only_accuracy": float(context_metrics["accuracy"]),
        "methyl_context_accuracy": float(methyl_metrics["accuracy"]),
        "delta_accuracy": float(methyl_metrics["accuracy"]) - float(context_metrics["accuracy"]),
        "context_only_balanced_accuracy": float(context_metrics["balanced_accuracy"]),
        "methyl_context_balanced_accuracy": float(methyl_metrics["balanced_accuracy"]),
        "delta_balanced_accuracy": float(methyl_metrics["balanced_accuracy"]) - float(context_metrics["balanced_accuracy"]),
        "context_only_correct_methyl_context_wrong": a_only,
        "context_only_wrong_methyl_context_correct": b_only,
        "discordant_total": discordant_total,
        "mcnemar_pvalue": pvalue,
        **bootstrap,
        "incremental_supported": bool(
            not math.isnan(float(bootstrap["bootstrap_delta_f1_ci_low"]))
            and float(bootstrap["bootstrap_delta_f1_ci_low"]) > 0
        ),
    }


def comparison_rows_by_dataset(
    evaluation_split: str,
    random_state: int,
    bootstrap_reps: int,
    test_df: pd.DataFrame,
    context_only_pred: np.ndarray,
    methyl_context_pred: np.ndarray,
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    working = test_df.copy()
    working["context_only_pred"] = context_only_pred
    working["methyl_context_pred"] = methyl_context_pred

    grouped = working.groupby(["dataset_id", "dataset_label", "harmonization_group"], sort=True, dropna=False)
    for (dataset_id, dataset_label, harmonization_group), sub in grouped:
        context_pred = sub["context_only_pred"].to_numpy(dtype=int)
        methyl_pred = sub["methyl_context_pred"].to_numpy(dtype=int)
        context_metrics = compute_metrics(sub["label"], context_pred)
        methyl_metrics = compute_metrics(sub["label"], methyl_pred)
        context_correct = context_pred == sub["label"].to_numpy(dtype=int)
        methyl_correct = methyl_pred == sub["label"].to_numpy(dtype=int)
        a_only, b_only, discordant_total, pvalue = mcnemar_pvalue(context_correct, methyl_correct)
        bootstrap = bootstrap_delta_metrics(
            test_df=sub,
            context_only_pred=context_pred,
            methyl_context_pred=methyl_pred,
            reps=bootstrap_reps,
            seed=random_state,
        )
        rows.append(
            {
                "evaluation_split": evaluation_split,
                "random_state": random_state,
                "dataset_id": dataset_id,
                "dataset_label": dataset_label,
                "harmonization_group": harmonization_group,
                "rows_total": int(len(sub.index)),
                "regions_total": int(sub["region_key"].nunique()),
                "context_only_F1": float(context_metrics["F1"]),
                "methyl_context_F1": float(methyl_metrics["F1"]),
                "delta_F1": float(methyl_metrics["F1"]) - float(context_metrics["F1"]),
                "context_only_accuracy": float(context_metrics["accuracy"]),
                "methyl_context_accuracy": float(methyl_metrics["accuracy"]),
                "delta_accuracy": float(methyl_metrics["accuracy"]) - float(context_metrics["accuracy"]),
                "context_only_correct_methyl_context_wrong": a_only,
                "context_only_wrong_methyl_context_correct": b_only,
                "discordant_total": discordant_total,
                "mcnemar_pvalue": pvalue,
                "bootstrap_delta_f1_mean": bootstrap["bootstrap_delta_f1_mean"],
                "bootstrap_delta_f1_ci_low": bootstrap["bootstrap_delta_f1_ci_low"],
                "bootstrap_delta_f1_ci_high": bootstrap["bootstrap_delta_f1_ci_high"],
                "bootstrap_delta_f1_positive_fraction": bootstrap["bootstrap_delta_f1_positive_fraction"],
                "incremental_supported": bool(
                    not math.isnan(float(bootstrap["bootstrap_delta_f1_ci_low"]))
                    and float(bootstrap["bootstrap_delta_f1_ci_low"]) > 0
                ),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_tsv).resolve()
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    df = load_benchmark_rows(input_path)
    discovery_df = df[df["dataset_role"] == "discovery"].copy()
    validation_df = df[df["dataset_role"] == "validation"].copy()
    discovery_train_df, discovery_holdout_df = pick_discovery_holdout(discovery_df, args.test_size, args.random_state)

    overall_rows: List[Dict[str, object]] = []
    dataset_rows: List[Dict[str, object]] = []

    discovery_context_pred, discovery_methyl_pred = fit_context_models(discovery_train_df, discovery_holdout_df)
    overall_rows.append(
        comparison_row(
            evaluation_split="discovery_holdout",
            random_state=args.random_state,
            bootstrap_reps=args.bootstrap_reps,
            train_df=discovery_train_df,
            test_df=discovery_holdout_df,
            context_only_pred=discovery_context_pred,
            methyl_context_pred=discovery_methyl_pred,
        )
    )
    dataset_rows.extend(
        comparison_rows_by_dataset(
            evaluation_split="discovery_holdout",
            random_state=args.random_state,
            bootstrap_reps=args.bootstrap_reps,
            test_df=discovery_holdout_df,
            context_only_pred=discovery_context_pred,
            methyl_context_pred=discovery_methyl_pred,
        )
    )

    if not validation_df.empty:
        validation_context_pred, validation_methyl_pred = fit_context_models(discovery_df, validation_df)
        overall_rows.append(
            comparison_row(
                evaluation_split="external_validation",
                random_state=args.random_state,
                bootstrap_reps=args.bootstrap_reps,
                train_df=discovery_df,
                test_df=validation_df,
                context_only_pred=validation_context_pred,
                methyl_context_pred=validation_methyl_pred,
            )
        )
        dataset_rows.extend(
            comparison_rows_by_dataset(
                evaluation_split="external_validation",
                random_state=args.random_state,
                bootstrap_reps=args.bootstrap_reps,
                test_df=validation_df,
                context_only_pred=validation_context_pred,
                methyl_context_pred=validation_methyl_pred,
            )
        )

    write_tsv_rows(output_dir / "incremental_comparison.tsv", COMPARISON_FIELDS, overall_rows)
    write_tsv_rows(output_dir / "incremental_dataset_comparison.tsv", DATASET_FIELDS, dataset_rows)
    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1A incremental test",
            "input_tsv": str(input_path),
            "random_state": args.random_state,
            "test_size": args.test_size,
            "bootstrap_reps": args.bootstrap_reps,
            "output_dir": str(output_dir),
        },
    )

    notes = [
        "# Phase 1A Incremental Test",
        "",
        f"- input_tsv: `{input_path}`",
        f"- random_state: `{args.random_state}`",
        f"- test_size: `{args.test_size}`",
        f"- bootstrap_reps: `{args.bootstrap_reps}`",
        "",
        "## Outputs",
        "",
        "- `incremental_comparison.tsv`",
        "- `incremental_dataset_comparison.tsv`",
        "- `run_context.json`",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"[phase1a-incremental] wrote {output_dir / 'incremental_comparison.tsv'}")
    print(f"[phase1a-incremental] wrote {output_dir / 'incremental_dataset_comparison.tsv'}")


if __name__ == "__main__":
    main()
