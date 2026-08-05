#!/usr/bin/env python3
"""Summarize Phase 1A prediction errors from benchmark prediction tables."""

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
    "20260325_phase1a_read_classifier_benchmark_sample200_v1/benchmark_predictions.tsv"
)

DEFAULT_OUTPUT = Path(
    "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/"
    "20260325_phase1a_to_pure_failure_diagnosis_v1"
)


SUMMARY_FIELDS = [
    "model_name",
    "evaluation_split",
    "harmonization_group",
    "truth_status",
    "VerificationClass",
    "PassedGating",
    "phase1a_read_label",
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
    "error_count",
    "error_rate",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions-tsv", default=str(DEFAULT_INPUT), help="Input benchmark_predictions.tsv path.")
    parser.add_argument("--model-name", default="logistic_methyl_context", help="Filter by model_name.")
    parser.add_argument("--evaluation-split", default="external_validation", help="Filter by evaluation split.")
    parser.add_argument("--harmonization-group", default="ONT_Dorado|to-pure", help="Filter by harmonization group.")
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


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.predictions_tsv, sep="\t", low_memory=False, keep_default_na=False)
    try:
        df, verification_metadata = attach_current_taxonomy(
            df,
            allow_unversioned_v1=args.allow_unversioned_v1,
        )
    except SchemaContractError as exc:
        raise SystemExit(f"[phase1a-failure][schema-contract] {exc}") from exc
    output_dir = ensure_dir(Path(args.output_dir).resolve())

    filtered = df[
        (df["model_name"] == args.model_name)
        & (df["evaluation_split"] == args.evaluation_split)
        & (df["harmonization_group"] == args.harmonization_group)
    ].copy()

    group_columns = [
        "model_name",
        "evaluation_split",
        "harmonization_group",
        "truth_status",
        "VerificationClass",
        "PassedGating",
        "phase1a_read_label",
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

    rows: List[Dict[str, object]] = []
    if not filtered.empty:
        grouped = filtered.groupby(group_columns, dropna=False, sort=True)
        for keys, sub in grouped:
            row = {column: value for column, value in zip(group_columns, keys)}
            error_count = int(sub["is_error"].astype(bool).sum())
            row.update(
                {
                    "rows_total": int(len(sub.index)),
                    "error_count": error_count,
                    "error_rate": float(error_count / len(sub.index)),
                }
            )
            rows.append(row)

    total_row = {
        "model_name": args.model_name,
        "evaluation_split": args.evaluation_split,
        "harmonization_group": args.harmonization_group,
        "VerificationSchemaStatus": verification_metadata["schema_status"],
        "VerificationSchemaVersion": 2 if verification_metadata["schema_status"] == "V2" else "",
        "rows_total": int(len(filtered.index)),
        "error_count": int(filtered["is_error"].astype(bool).sum()) if not filtered.empty else 0,
        "error_rate": float(filtered["is_error"].astype(bool).mean()) if not filtered.empty else 0.0,
    }

    write_tsv_rows(output_dir / "failure_breakdown.tsv", SUMMARY_FIELDS, rows)
    write_tsv_rows(
        output_dir / "failure_summary.tsv",
        [
            "model_name",
            "evaluation_split",
            "harmonization_group",
            "VerificationSchemaStatus",
            "VerificationSchemaVersion",
            "rows_total",
            "error_count",
            "error_rate",
        ],
        [total_row],
    )
    write_json(
        output_dir / "run_context.json",
        {
            "task": "Phase 1A prediction failure summary",
            "predictions_tsv": args.predictions_tsv,
            "model_name": args.model_name,
            "evaluation_split": args.evaluation_split,
            "harmonization_group": args.harmonization_group,
            "verification_taxonomy": verification_metadata,
            "output_dir": str(output_dir),
        },
    )

    notes = [
        "# Phase 1A Failure Diagnosis",
        "",
        f"- predictions_tsv: `{args.predictions_tsv}`",
        f"- model_name: `{args.model_name}`",
        f"- evaluation_split: `{args.evaluation_split}`",
        f"- harmonization_group: `{args.harmonization_group}`",
        "",
        "## Outputs",
        "",
        "- `failure_breakdown.tsv`",
        "- `failure_summary.tsv`",
        "- `run_context.json`",
    ]
    (output_dir / "round_summary.md").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"[phase1a-failure] wrote {output_dir / 'failure_breakdown.tsv'}")
    print(f"[phase1a-failure] wrote {output_dir / 'failure_summary.tsv'}")


if __name__ == "__main__":
    main()
