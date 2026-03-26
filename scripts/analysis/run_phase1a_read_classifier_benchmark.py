#!/usr/bin/env python3
"""Run a Phase 1A read-level ALT/REF benchmark on sharded training-table exports."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, precision_score, recall_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from research_common import ensure_dir, write_json, write_tsv_rows


DEFAULT_INPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1_manifest_shard_export_sample80/phase1_shard_read_training_table.tsv"
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_read_classifier_benchmark_sample80_v1"
)


NUMERIC_FEATURES_METHYL = [
    "mapq",
    "num_cpg_observed",
    "methyl_na_fraction",
    "methyl_mean",
    "methyl_std",
    "methyl_median",
    "methyl_high_fraction",
    "methyl_low_fraction",
]


NUMERIC_FEATURES_CONTEXT = NUMERIC_FEATURES_METHYL + [
    "PairwiseMedianDist",
    "AlleleDelta",
    "CramersV",
]


NUMERIC_FEATURES_CONTEXT_ONLY = [
    "mapq",
    "PairwiseMedianDist",
    "AlleleDelta",
    "CramersV",
]


CATEGORICAL_FEATURES_CONTEXT = [
    "hp",
    "mode",
    "VerificationClass",
    "PassedGating",
]


METRIC_FIELDS = [
    "model_name",
    "evaluation_split",
    "train_rows",
    "test_rows",
    "train_regions",
    "test_regions",
    "truth_total",
    "calls_total",
    "TP",
    "FP",
    "FN",
    "TN",
    "precision",
    "recall",
    "F1",
    "accuracy",
    "balanced_accuracy",
    "positive_label",
    "input_tsv",
]


GROUP_METRIC_FIELDS = [
    "model_name",
    "evaluation_split",
    "harmonization_group",
    "rows_total",
    "truth_total",
    "calls_total",
    "TP",
    "FP",
    "FN",
    "TN",
    "precision",
    "recall",
    "F1",
    "accuracy",
    "balanced_accuracy",
]


DATASET_METRIC_FIELDS = [
    "model_name",
    "evaluation_split",
    "dataset_id",
    "dataset_label",
    "harmonization_group",
    "rows_total",
    "truth_total",
    "calls_total",
    "TP",
    "FP",
    "FN",
    "TN",
    "precision",
    "recall",
    "F1",
    "accuracy",
    "balanced_accuracy",
]


COEF_FIELDS = [
    "model_name",
    "feature_name",
    "coefficient",
    "abs_coefficient",
]


SUMMARY_FIELDS = [
    "input_tsv",
    "discovery_rows",
    "validation_rows",
    "discovery_regions",
    "validation_regions",
    "models_evaluated",
]


PREDICTION_FIELDS = [
    "model_name",
    "evaluation_split",
    "dataset_id",
    "dataset_label",
    "dataset_role",
    "harmonization_group",
    "region_key",
    "read_id",
    "phase1a_read_label",
    "predicted_label",
    "is_error",
    "truth_status",
    "VerificationClass",
    "PassedGating",
]


@dataclass(frozen=True)
class ThresholdRule:
    threshold: float
    alt_if_leq: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tsv", default=str(DEFAULT_INPUT), help="Read-level Phase 1A shard TSV.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="Output directory.")
    parser.add_argument("--test-size", type=float, default=0.2, help="Discovery region holdout fraction.")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed for group split.")
    return parser.parse_args()


def load_benchmark_rows(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df[df["phase1a_read_label"].isin(["ALT", "REF"])].copy()
    df["label"] = (df["phase1a_read_label"] == "ALT").astype(int)

    numeric_columns = sorted(set(NUMERIC_FEATURES_CONTEXT))
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def pick_discovery_holdout(discovery_df: pd.DataFrame, test_size: float, random_state: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_idx, test_idx = next(splitter.split(discovery_df, groups=discovery_df["region_key"]))
    return discovery_df.iloc[train_idx].copy(), discovery_df.iloc[test_idx].copy()


def compute_metrics(y_true: Sequence[int], y_pred: Sequence[int]) -> Dict[str, object]:
    y_true_arr = np.asarray(y_true, dtype=int)
    y_pred_arr = np.asarray(y_pred, dtype=int)
    tp = int(((y_true_arr == 1) & (y_pred_arr == 1)).sum())
    fp = int(((y_true_arr == 0) & (y_pred_arr == 1)).sum())
    fn = int(((y_true_arr == 1) & (y_pred_arr == 0)).sum())
    tn = int(((y_true_arr == 0) & (y_pred_arr == 0)).sum())
    return {
        "truth_total": int(len(y_true_arr)),
        "calls_total": int(len(y_pred_arr)),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "precision": float(precision_score(y_true_arr, y_pred_arr, zero_division=0)),
        "recall": float(recall_score(y_true_arr, y_pred_arr, zero_division=0)),
        "F1": float(f1_score(y_true_arr, y_pred_arr, zero_division=0)),
        "accuracy": float(accuracy_score(y_true_arr, y_pred_arr)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true_arr, y_pred_arr)),
    }


def evaluate_predictions(
    model_name: str,
    evaluation_split: str,
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    predictions: Sequence[int],
    input_tsv: str,
) -> Dict[str, object]:
    metrics = compute_metrics(test_df["label"], predictions)
    metrics.update(
        {
            "model_name": model_name,
            "evaluation_split": evaluation_split,
            "train_rows": int(len(train_df.index)),
            "test_rows": int(len(test_df.index)),
            "train_regions": int(train_df["region_key"].nunique()),
            "test_regions": int(test_df["region_key"].nunique()),
            "positive_label": "ALT",
            "input_tsv": input_tsv,
        }
    )
    return metrics


def evaluate_by_group(
    model_name: str,
    evaluation_split: str,
    test_df: pd.DataFrame,
    predictions: Sequence[int],
) -> List[Dict[str, object]]:
    result_rows: List[Dict[str, object]] = []
    working = test_df.copy()
    working["prediction"] = np.asarray(predictions, dtype=int)
    for group_name, sub in working.groupby("harmonization_group", sort=True):
        metrics = compute_metrics(sub["label"], sub["prediction"])
        metrics.update(
            {
                "model_name": model_name,
                "evaluation_split": evaluation_split,
                "harmonization_group": group_name,
                "rows_total": int(len(sub.index)),
            }
        )
        result_rows.append(metrics)
    return result_rows


def evaluate_by_dataset(
    model_name: str,
    evaluation_split: str,
    test_df: pd.DataFrame,
    predictions: Sequence[int],
) -> List[Dict[str, object]]:
    result_rows: List[Dict[str, object]] = []
    working = test_df.copy()
    working["prediction"] = np.asarray(predictions, dtype=int)
    grouped = working.groupby(["dataset_id", "dataset_label", "harmonization_group"], sort=True, dropna=False)
    for (dataset_id, dataset_label, harmonization_group), sub in grouped:
        metrics = compute_metrics(sub["label"], sub["prediction"])
        metrics.update(
            {
                "model_name": model_name,
                "evaluation_split": evaluation_split,
                "dataset_id": dataset_id,
                "dataset_label": dataset_label,
                "harmonization_group": harmonization_group,
                "rows_total": int(len(sub.index)),
            }
        )
        result_rows.append(metrics)
    return result_rows


def build_prediction_rows(
    model_name: str,
    evaluation_split: str,
    test_df: pd.DataFrame,
    predictions: Sequence[int],
) -> List[Dict[str, object]]:
    label_map = {0: "REF", 1: "ALT"}
    rows: List[Dict[str, object]] = []
    for (_, row), pred in zip(test_df.iterrows(), predictions):
        predicted_label = label_map[int(pred)]
        truth_label = str(row["phase1a_read_label"])
        rows.append(
            {
                "model_name": model_name,
                "evaluation_split": evaluation_split,
                "dataset_id": row["dataset_id"],
                "dataset_label": row["dataset_label"],
                "dataset_role": row["dataset_role"],
                "harmonization_group": row["harmonization_group"],
                "region_key": row["region_key"],
                "read_id": row["read_id"],
                "phase1a_read_label": truth_label,
                "predicted_label": predicted_label,
                "is_error": predicted_label != truth_label,
                "truth_status": row["truth_status"],
                "VerificationClass": row["VerificationClass"],
                "PassedGating": row["PassedGating"],
            }
        )
    return rows


def majority_ref_predict(rows: int) -> np.ndarray:
    return np.zeros(rows, dtype=int)


def fit_methyl_threshold(train_df: pd.DataFrame) -> ThresholdRule:
    observed = train_df["methyl_mean"].dropna().sort_values().unique()
    if len(observed) == 0:
        return ThresholdRule(threshold=0.5, alt_if_leq=True)

    candidates = sorted(set([float(observed[0]), float(observed[-1]), 0.5] + [float(value) for value in np.quantile(observed, np.linspace(0.05, 0.95, 19))]))
    best_score = -math.inf
    best_rule = ThresholdRule(threshold=0.5, alt_if_leq=True)

    y_true = train_df["label"].to_numpy(dtype=int)
    methyl = train_df["methyl_mean"].to_numpy(dtype=float)
    for threshold in candidates:
        for alt_if_leq in [True, False]:
            if alt_if_leq:
                y_pred = (methyl <= threshold).astype(int)
            else:
                y_pred = (methyl >= threshold).astype(int)
            score = f1_score(y_true, y_pred, zero_division=0)
            if score > best_score:
                best_score = score
                best_rule = ThresholdRule(threshold=float(threshold), alt_if_leq=alt_if_leq)
    return best_rule


def apply_methyl_threshold(df: pd.DataFrame, rule: ThresholdRule) -> np.ndarray:
    methyl = df["methyl_mean"].to_numpy(dtype=float)
    if rule.alt_if_leq:
        return (methyl <= rule.threshold).astype(int)
    return (methyl >= rule.threshold).astype(int)


def build_logistic_pipeline(
    numeric_features: Sequence[str],
    categorical_features: Sequence[str],
) -> Pipeline:
    transformer_steps = [
        (
            "num",
            Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                ]
            ),
            list(numeric_features),
        )
    ]
    if categorical_features:
        transformer_steps.append(
            (
                "cat",
                Pipeline(
                    [
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("encoder", OneHotEncoder(handle_unknown="ignore")),
                    ]
                ),
                list(categorical_features),
            )
        )

    preprocessor = ColumnTransformer(transformer_steps)
    return Pipeline(
        [
            ("preprocessor", preprocessor),
            ("classifier", LogisticRegression(max_iter=2000, class_weight="balanced")),
        ]
    )


def extract_coefficients(model_name: str, pipeline: Pipeline) -> List[Dict[str, object]]:
    preprocessor: ColumnTransformer = pipeline.named_steps["preprocessor"]
    classifier: LogisticRegression = pipeline.named_steps["classifier"]
    feature_names = list(preprocessor.get_feature_names_out())
    coefficients = classifier.coef_[0]
    rows = [
        {
            "model_name": model_name,
            "feature_name": feature_name,
            "coefficient": float(coef),
            "abs_coefficient": float(abs(coef)),
        }
        for feature_name, coef in zip(feature_names, coefficients)
    ]
    rows.sort(key=lambda row: row["abs_coefficient"], reverse=True)
    return rows


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_tsv).resolve()
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    df = load_benchmark_rows(input_path)
    discovery_df = df[df["dataset_role"] == "discovery"].copy()
    validation_df = df[df["dataset_role"] == "validation"].copy()
    discovery_train_df, discovery_holdout_df = pick_discovery_holdout(discovery_df, args.test_size, args.random_state)

    metrics_rows: List[Dict[str, object]] = []
    group_rows: List[Dict[str, object]] = []
    dataset_rows: List[Dict[str, object]] = []
    coef_rows: List[Dict[str, object]] = []
    prediction_rows: List[Dict[str, object]] = []

    evaluations = [
        ("discovery_holdout", discovery_train_df, discovery_holdout_df),
        ("external_validation", discovery_df, validation_df),
    ]

    for evaluation_split, train_df, test_df in evaluations:
        majority_pred = majority_ref_predict(len(test_df.index))
        metrics_rows.append(
            evaluate_predictions("majority_ref", evaluation_split, train_df, test_df, majority_pred, str(input_path))
        )
        group_rows.extend(evaluate_by_group("majority_ref", evaluation_split, test_df, majority_pred))
        dataset_rows.extend(evaluate_by_dataset("majority_ref", evaluation_split, test_df, majority_pred))
        prediction_rows.extend(build_prediction_rows("majority_ref", evaluation_split, test_df, majority_pred))

        threshold_rule = fit_methyl_threshold(train_df)
        threshold_pred = apply_methyl_threshold(test_df, threshold_rule)
        threshold_metrics = evaluate_predictions(
            "methyl_mean_threshold", evaluation_split, train_df, test_df, threshold_pred, str(input_path)
        )
        threshold_metrics["threshold"] = threshold_rule.threshold
        threshold_metrics["alt_if_leq"] = threshold_rule.alt_if_leq
        metrics_rows.append(threshold_metrics)
        group_rows.extend(evaluate_by_group("methyl_mean_threshold", evaluation_split, test_df, threshold_pred))
        dataset_rows.extend(evaluate_by_dataset("methyl_mean_threshold", evaluation_split, test_df, threshold_pred))
        prediction_rows.extend(build_prediction_rows("methyl_mean_threshold", evaluation_split, test_df, threshold_pred))

        logistic_methyl = build_logistic_pipeline(NUMERIC_FEATURES_METHYL, [])
        logistic_methyl.fit(train_df[NUMERIC_FEATURES_METHYL], train_df["label"])
        logistic_methyl_pred = logistic_methyl.predict(test_df[NUMERIC_FEATURES_METHYL])
        metrics_rows.append(
            evaluate_predictions(
                "logistic_methyl_only",
                evaluation_split,
                train_df,
                test_df,
                logistic_methyl_pred,
                str(input_path),
            )
        )
        group_rows.extend(evaluate_by_group("logistic_methyl_only", evaluation_split, test_df, logistic_methyl_pred))
        dataset_rows.extend(evaluate_by_dataset("logistic_methyl_only", evaluation_split, test_df, logistic_methyl_pred))
        prediction_rows.extend(build_prediction_rows("logistic_methyl_only", evaluation_split, test_df, logistic_methyl_pred))
        coef_rows.extend(extract_coefficients("logistic_methyl_only__" + evaluation_split, logistic_methyl)[:20])

        logistic_context_only = build_logistic_pipeline(NUMERIC_FEATURES_CONTEXT_ONLY, CATEGORICAL_FEATURES_CONTEXT)
        context_only_features = list(NUMERIC_FEATURES_CONTEXT_ONLY) + list(CATEGORICAL_FEATURES_CONTEXT)
        logistic_context_only.fit(train_df[context_only_features], train_df["label"])
        logistic_context_only_pred = logistic_context_only.predict(test_df[context_only_features])
        metrics_rows.append(
            evaluate_predictions(
                "logistic_context_only",
                evaluation_split,
                train_df,
                test_df,
                logistic_context_only_pred,
                str(input_path),
            )
        )
        group_rows.extend(evaluate_by_group("logistic_context_only", evaluation_split, test_df, logistic_context_only_pred))
        dataset_rows.extend(evaluate_by_dataset("logistic_context_only", evaluation_split, test_df, logistic_context_only_pred))
        prediction_rows.extend(build_prediction_rows("logistic_context_only", evaluation_split, test_df, logistic_context_only_pred))
        coef_rows.extend(extract_coefficients("logistic_context_only__" + evaluation_split, logistic_context_only)[:30])

        logistic_context = build_logistic_pipeline(NUMERIC_FEATURES_CONTEXT, CATEGORICAL_FEATURES_CONTEXT)
        context_features = list(NUMERIC_FEATURES_CONTEXT) + list(CATEGORICAL_FEATURES_CONTEXT)
        logistic_context.fit(train_df[context_features], train_df["label"])
        logistic_context_pred = logistic_context.predict(test_df[context_features])
        metrics_rows.append(
            evaluate_predictions(
                "logistic_methyl_context",
                evaluation_split,
                train_df,
                test_df,
                logistic_context_pred,
                str(input_path),
            )
        )
        group_rows.extend(evaluate_by_group("logistic_methyl_context", evaluation_split, test_df, logistic_context_pred))
        dataset_rows.extend(evaluate_by_dataset("logistic_methyl_context", evaluation_split, test_df, logistic_context_pred))
        prediction_rows.extend(build_prediction_rows("logistic_methyl_context", evaluation_split, test_df, logistic_context_pred))
        coef_rows.extend(extract_coefficients("logistic_methyl_context__" + evaluation_split, logistic_context)[:30])

    summary_rows = [
        {
            "input_tsv": str(input_path),
            "discovery_rows": int(len(discovery_df.index)),
            "validation_rows": int(len(validation_df.index)),
            "discovery_regions": int(discovery_df["region_key"].nunique()),
            "validation_regions": int(validation_df["region_key"].nunique()),
            "models_evaluated": 5,
        }
    ]

    write_tsv_rows(output_dir / "benchmark_metrics.tsv", METRIC_FIELDS + ["threshold", "alt_if_leq"], metrics_rows)
    write_tsv_rows(output_dir / "benchmark_group_metrics.tsv", GROUP_METRIC_FIELDS, group_rows)
    write_tsv_rows(output_dir / "benchmark_dataset_metrics.tsv", DATASET_METRIC_FIELDS, dataset_rows)
    write_tsv_rows(output_dir / "benchmark_top_coefficients.tsv", COEF_FIELDS, coef_rows)
    write_tsv_rows(output_dir / "benchmark_predictions.tsv", PREDICTION_FIELDS, prediction_rows)
    write_tsv_rows(output_dir / "benchmark_summary.tsv", SUMMARY_FIELDS, summary_rows)
    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1A read classifier benchmark",
            "input_tsv": str(input_path),
            "test_size": args.test_size,
            "random_state": args.random_state,
            "output_dir": str(output_dir),
        },
    )

    notes = [
        "# Phase 1A Read Classifier Benchmark",
        "",
        f"- input_tsv: `{input_path}`",
        f"- test_size: `{args.test_size}`",
        f"- random_state: `{args.random_state}`",
        "",
        "## Outputs",
        "",
        "- `benchmark_metrics.tsv`",
        "- `benchmark_group_metrics.tsv`",
        "- `benchmark_dataset_metrics.tsv`",
        "- `benchmark_top_coefficients.tsv`",
        "- `benchmark_predictions.tsv`",
        "- `benchmark_summary.tsv`",
        "- `run_context.json`",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"[phase1a-benchmark] wrote {output_dir / 'benchmark_metrics.tsv'}")
    print(f"[phase1a-benchmark] wrote {output_dir / 'benchmark_group_metrics.tsv'}")
    print(f"[phase1a-benchmark] wrote {output_dir / 'benchmark_dataset_metrics.tsv'}")
    print(f"[phase1a-benchmark] wrote {output_dir / 'benchmark_top_coefficients.tsv'}")


if __name__ == "__main__":
    main()
