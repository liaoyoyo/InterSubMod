#!/usr/bin/env python3
"""Compare Phase 1A model error rates across datasets and verification buckets."""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from typing import Dict, List

import pandas as pd

from research_common import ensure_dir, write_json, write_tsv_rows

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    UNKNOWN_CURRENT_CLASS,
    V1_CURRENT_CLASSES,
    SchemaContractError,
    read_evidence,
    select_current_view,
)


DEFAULT_INPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_read_classifier_benchmark_paired_multibio_sample637_v1/benchmark_predictions.tsv"
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_model_error_compare_v1"
)


PAIR_KEYS = [
    "evaluation_split",
    "dataset_id",
    "dataset_label",
    "dataset_role",
    "harmonization_group",
    "region_key",
    "read_id",
    "phase1a_read_label",
    "truth_status",
    "VerificationClass",
    "PassedGating",
    "VerificationSchemaStatus",
    "VerificationSchemaVersion",
    "VerificationClass_V1_Deprecated",
    "VerificationClass_Legacy",
    "LabelFirstSupport",
    "ClusterFirstSupport",
    "WithinHPSupport",
    "DispersionWarning",
    "EvidencePath",
    "EvidenceDerivation",
]


BUCKET_FIELDS = [
    "evaluation_split",
    "dataset_id",
    "dataset_label",
    "harmonization_group",
    "truth_status",
    "VerificationClass",
    "PassedGating",
    "VerificationSchemaStatus",
    "VerificationSchemaVersion",
    "VerificationClass_V1_Deprecated",
    "VerificationClass_Legacy",
    "LabelFirstSupport",
    "ClusterFirstSupport",
    "WithinHPSupport",
    "DispersionWarning",
    "EvidencePath",
    "EvidenceDerivation",
    "rows_total",
    "model_a",
    "model_a_error_count",
    "model_a_error_rate",
    "model_b",
    "model_b_error_count",
    "model_b_error_rate",
    "delta_error_count",
    "delta_error_rate",
    "improved_rows",
    "worsened_rows",
]


DATASET_FIELDS = [
    "evaluation_split",
    "dataset_id",
    "dataset_label",
    "harmonization_group",
    "VerificationSchemaStatus",
    "VerificationSchemaVersion",
    "rows_total",
    "model_a",
    "model_a_error_count",
    "model_a_error_rate",
    "model_b",
    "model_b_error_count",
    "model_b_error_rate",
    "delta_error_count",
    "delta_error_rate",
    "improved_rows",
    "worsened_rows",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions-tsv", default=str(DEFAULT_INPUT), help="Input benchmark_predictions.tsv path.")
    parser.add_argument("--evaluation-split", default="external_validation", help="Filter by evaluation split.")
    parser.add_argument("--model-a", default="logistic_context_only", help="Baseline model name.")
    parser.add_argument("--model-b", default="logistic_methyl_context", help="Comparison model name.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="Output directory.")
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help="Explicitly authorize a historical unversioned v1 current taxonomy; evidence remains unavailable.",
    )
    return parser.parse_args()


PROVENANCE_FIELDS = [
    "VerificationClass_V1_Deprecated",
    "VerificationClass_Legacy",
    "LabelFirstSupport",
    "ClusterFirstSupport",
    "WithinHPSupport",
    "DispersionWarning",
    "EvidencePath",
    "EvidenceDerivation",
]


def attach_current_taxonomy(
    df: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> tuple[pd.DataFrame, Dict[str, object]]:
    """Validate current taxonomy identity and preserve v2 evidence provenance."""
    result = df.copy()
    if "VerificationSchemaVersion" in result.columns:
        view = select_current_view(result)
        missing = [field for field in PROVENANCE_FIELDS if field not in result.columns]
        if missing:
            raise SchemaContractError(
                "Phase 1A predictions: missing v2 provenance fields: " + ", ".join(missing)
            )
        read_evidence(result)
        result["VerificationClass"] = view.values
        result["VerificationSchemaStatus"] = view.schema_status
        metadata = view.metadata()
        metadata["evidence_fields"] = list(PROVENANCE_FIELDS)
        return result, metadata

    if not allow_unversioned_v1:
        raise SchemaContractError(
            "Phase 1A predictions are unversioned; --allow-unversioned-v1 is required"
        )
    if "VerificationClass" not in result.columns:
        raise SchemaContractError("Phase 1A predictions: VerificationClass is missing")

    raw = result["VerificationClass"].astype("string")
    valid = raw.isin(V1_CURRENT_CLASSES)
    unknown_counts = {
        str(key): int(value)
        for key, value in raw[~valid].fillna("<MISSING>").value_counts().sort_index().items()
    }
    result["VerificationClass"] = raw.where(valid, UNKNOWN_CURRENT_CLASS)
    result["VerificationSchemaVersion"] = ""
    result["VerificationSchemaStatus"] = "UNVERSIONED_V1"
    for field in PROVENANCE_FIELDS:
        result[field] = ""
    message = (
        "UNVERSIONED: Phase 1A current taxonomy accepted only under "
        "--allow-unversioned-v1; v2 evidence provenance is unavailable"
    )
    warnings.warn(message, UserWarning, stacklevel=2)
    return result, {
        "selection_field": "VerificationClass",
        "schema_status": "UNVERSIONED_V1",
        "categories": list(V1_CURRENT_CLASSES) + [UNKNOWN_CURRENT_CLASS],
        "unknown_counts": unknown_counts,
        "warnings": [message],
        "evidence_fields": [],
    }


def build_pairwise_table(df: pd.DataFrame, model_a: str, model_b: str) -> pd.DataFrame:
    filtered = df[
        (df["model_name"].isin([model_a, model_b]))
    ].copy()

    pivot = filtered.pivot_table(
        index=PAIR_KEYS,
        columns="model_name",
        values="is_error",
        aggfunc="first",
    ).reset_index()
    if model_a not in pivot.columns or model_b not in pivot.columns:
        return pivot.iloc[0:0].copy()

    pivot["model_a_error"] = pivot[model_a].astype(bool)
    pivot["model_b_error"] = pivot[model_b].astype(bool)
    pivot["improved"] = pivot["model_a_error"] & (~pivot["model_b_error"])
    pivot["worsened"] = (~pivot["model_a_error"]) & pivot["model_b_error"]
    return pivot


def summarize_group(
    sub: pd.DataFrame,
    model_a: str,
    model_b: str,
) -> Dict[str, object]:
    rows_total = int(len(sub.index))
    model_a_error_count = int(sub["model_a_error"].sum())
    model_b_error_count = int(sub["model_b_error"].sum())
    return {
        "rows_total": rows_total,
        "model_a": model_a,
        "model_a_error_count": model_a_error_count,
        "model_a_error_rate": float(model_a_error_count / rows_total) if rows_total else 0.0,
        "model_b": model_b,
        "model_b_error_count": model_b_error_count,
        "model_b_error_rate": float(model_b_error_count / rows_total) if rows_total else 0.0,
        "delta_error_count": model_b_error_count - model_a_error_count,
        "delta_error_rate": float((model_b_error_count - model_a_error_count) / rows_total) if rows_total else 0.0,
        "improved_rows": int(sub["improved"].sum()),
        "worsened_rows": int(sub["worsened"].sum()),
    }


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.predictions_tsv, sep="\t", low_memory=False, keep_default_na=False)
    try:
        df, verification_metadata = attach_current_taxonomy(
            df,
            allow_unversioned_v1=args.allow_unversioned_v1,
        )
    except SchemaContractError as exc:
        raise SystemExit(f"[phase1a-compare][schema-contract] {exc}") from exc
    df = df[df["evaluation_split"] == args.evaluation_split].copy()
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    pairwise = build_pairwise_table(df, args.model_a, args.model_b)

    bucket_rows: List[Dict[str, object]] = []
    if not pairwise.empty:
        grouped = pairwise.groupby(
            [
                "evaluation_split",
                "dataset_id",
                "dataset_label",
                "harmonization_group",
                "truth_status",
                "VerificationClass",
                "PassedGating",
                "VerificationSchemaStatus",
                "VerificationSchemaVersion",
                "VerificationClass_V1_Deprecated",
                "VerificationClass_Legacy",
                "LabelFirstSupport",
                "ClusterFirstSupport",
                "WithinHPSupport",
                "DispersionWarning",
                "EvidencePath",
                "EvidenceDerivation",
            ],
            sort=True,
            dropna=False,
        )
        for keys, sub in grouped:
            row = {field: value for field, value in zip(
                [
                    "evaluation_split",
                    "dataset_id",
                    "dataset_label",
                    "harmonization_group",
                    "truth_status",
                    "VerificationClass",
                    "PassedGating",
                    "VerificationSchemaStatus",
                    "VerificationSchemaVersion",
                    "VerificationClass_V1_Deprecated",
                    "VerificationClass_Legacy",
                    "LabelFirstSupport",
                    "ClusterFirstSupport",
                    "WithinHPSupport",
                    "DispersionWarning",
                    "EvidencePath",
                    "EvidenceDerivation",
                ],
                keys,
            )}
            row.update(summarize_group(sub, args.model_a, args.model_b))
            bucket_rows.append(row)

    dataset_rows: List[Dict[str, object]] = []
    if not pairwise.empty:
        grouped = pairwise.groupby(
            [
                "evaluation_split",
                "dataset_id",
                "dataset_label",
                "harmonization_group",
                "VerificationSchemaStatus",
                "VerificationSchemaVersion",
            ],
            sort=True,
            dropna=False,
        )
        for keys, sub in grouped:
            row = {field: value for field, value in zip(
                [
                    "evaluation_split",
                    "dataset_id",
                    "dataset_label",
                    "harmonization_group",
                    "VerificationSchemaStatus",
                    "VerificationSchemaVersion",
                ],
                keys,
            )}
            row.update(summarize_group(sub, args.model_a, args.model_b))
            dataset_rows.append(row)

    write_tsv_rows(output_dir / "model_error_bucket_summary.tsv", BUCKET_FIELDS, bucket_rows)
    write_tsv_rows(output_dir / "model_error_dataset_summary.tsv", DATASET_FIELDS, dataset_rows)
    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1A model error comparison",
            "predictions_tsv": args.predictions_tsv,
            "evaluation_split": args.evaluation_split,
            "model_a": args.model_a,
            "model_b": args.model_b,
            "verification_taxonomy": verification_metadata,
            "output_dir": str(output_dir),
        },
    )

    notes = [
        "# Phase 1A Model Error Comparison",
        "",
        f"- predictions_tsv: `{args.predictions_tsv}`",
        f"- evaluation_split: `{args.evaluation_split}`",
        f"- model_a: `{args.model_a}`",
        f"- model_b: `{args.model_b}`",
        "",
        "## Outputs",
        "",
        "- `model_error_bucket_summary.tsv`",
        "- `model_error_dataset_summary.tsv`",
        "- `run_context.json`",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"[phase1a-compare] wrote {output_dir / 'model_error_bucket_summary.tsv'}")
    print(f"[phase1a-compare] wrote {output_dir / 'model_error_dataset_summary.tsv'}")


if __name__ == "__main__":
    main()
