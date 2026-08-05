#!/usr/bin/env python3
"""Extract explicit label-first metrics from significance_summary.csv."""

import argparse
import csv
import os
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    SchemaContractError,
    read_evidence,
    select_current_view,
)


OUTPUT_COLUMNS = [
    "RegionKey",
    "RegionID",
    "NumReads",
    "PairwiseMedianDist",
    "ClusterGlobalP",
    "ClusterCramersV",
    "ClusterPassedGating",
    "ClusterPermanovaF",
    "ClusterPermanovaP",
    "ClusterPermanovaValid",
    "ClusterDispersionP",
    "ClusterDispersionWarn",
    "HPMergedDelta",
    "HPMergedP",
    "HPMergedSig",
    "HPFineF",
    "HPFineP",
    "HPFineSig",
    "HPFineNGroups",
    "AlleleDelta",
    "AlleleP",
    "AlleleSig",
    "LabelHPPermanovaF",
    "LabelHPPermanovaP",
    "LabelHPPermanovaValid",
    "LabelHPDispersionP",
    "LabelHPDispersionWarn",
    "LabelAllelePermanovaF",
    "LabelAllelePermanovaP",
    "LabelAllelePermanovaValid",
    "LabelAlleleDispersionP",
    "LabelAlleleDispersionWarn",
    "DominantLabel",
    "VerificationClass",
    "Significant",
    "SuggestFilter",
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


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-csv", required=True, help="Input significance_summary.csv")
    parser.add_argument(
        "--output-tsv",
        default=None,
        help="Output TSV path (default: same directory/label_first_metrics.tsv)",
    )
    return parser.parse_args()


def validate_label_first_schema(rows):
    """Require schema-v2 evidence; label-first support is never inferred from class names."""
    df = pd.DataFrame(rows)
    current = select_current_view(df)
    evidence = read_evidence(df)
    required = [
        "Chr",
        "Pos",
        "Ref",
        "Alt",
        "VerificationClass_V1_Deprecated",
        "VerificationClass_Legacy",
    ]
    missing = [field for field in required if field not in df.columns]
    if missing:
        raise SchemaContractError(
            "label-first metrics: missing required input fields: " + ", ".join(missing)
        )
    df["VerificationClass"] = current.values
    # ``read_evidence`` is a validator here.  Preserve the serialized source
    # values verbatim so the derived TSV remains a true pass-through artifact
    # (not a pandas-specific re-serialization of nullable booleans).
    del evidence
    return df.to_dict(orient="records"), current.metadata()


def main():
    args = parse_args()
    output_tsv = args.output_tsv
    if output_tsv is None:
        output_tsv = os.path.join(os.path.dirname(args.summary_csv), "label_first_metrics.tsv")

    with open(args.summary_csv, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    try:
        rows, _verification_metadata = validate_label_first_schema(rows)
    except SchemaContractError as exc:
        raise SystemExit(f"[extract_label_first_metrics][schema-contract] {exc}") from exc

    os.makedirs(os.path.dirname(output_tsv) or ".", exist_ok=True)
    with open(output_tsv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()

        for row in rows:
            region_key = f"{row['Chr']}:{row['Pos']}:{row['Ref']}:{row['Alt']}"
            writer.writerow(
                {
                    "RegionKey": region_key,
                    "RegionID": row.get("RegionID", ""),
                    "NumReads": row.get("NumReads", ""),
                    "PairwiseMedianDist": row.get("PairwiseMedianDist", ""),
                    "ClusterGlobalP": row.get("GlobalP", ""),
                    "ClusterCramersV": row.get("CramersV", ""),
                    "ClusterPassedGating": row.get("PassedGating", ""),
                    "ClusterPermanovaF": row.get("ClusterPermanovaF", ""),
                    "ClusterPermanovaP": row.get("ClusterPermanovaP", ""),
                    "ClusterPermanovaValid": row.get("ClusterPermanovaValid", ""),
                    "ClusterDispersionP": row.get("ClusterDispersionP", ""),
                    "ClusterDispersionWarn": row.get("ClusterDispersionWarn", ""),
                    "HPMergedDelta": row.get("HPMergedDelta", ""),
                    "HPMergedP": row.get("HPMergedP", ""),
                    "HPMergedSig": row.get("HPMergedSig", ""),
                    "HPFineF": row.get("HPFineF", ""),
                    "HPFineP": row.get("HPFineP", ""),
                    "HPFineSig": row.get("HPFineSig", ""),
                    "HPFineNGroups": row.get("HPFineNGroups", ""),
                    "AlleleDelta": row.get("AlleleDelta", ""),
                    "AlleleP": row.get("AlleleP", ""),
                    "AlleleSig": row.get("AlleleSig", ""),
                    "LabelHPPermanovaF": row.get("LabelHPPermanovaF", ""),
                    "LabelHPPermanovaP": row.get("LabelHPPermanovaP", ""),
                    "LabelHPPermanovaValid": row.get("LabelHPPermanovaValid", ""),
                    "LabelHPDispersionP": row.get("LabelHPDispersionP", ""),
                    "LabelHPDispersionWarn": row.get("LabelHPDispersionWarn", ""),
                    "LabelAllelePermanovaF": row.get("LabelAllelePermanovaF", ""),
                    "LabelAllelePermanovaP": row.get("LabelAllelePermanovaP", ""),
                    "LabelAllelePermanovaValid": row.get("LabelAllelePermanovaValid", ""),
                    "LabelAlleleDispersionP": row.get("LabelAlleleDispersionP", ""),
                    "LabelAlleleDispersionWarn": row.get("LabelAlleleDispersionWarn", ""),
                    "DominantLabel": row.get("DominantLabel", ""),
                    "VerificationClass": row.get("VerificationClass", ""),
                    "Significant": row.get("Significant", ""),
                    "SuggestFilter": row.get("SuggestFilter", ""),
                    "VerificationSchemaVersion": row.get("VerificationSchemaVersion", ""),
                    "VerificationClass_V1_Deprecated": row.get("VerificationClass_V1_Deprecated", ""),
                    "VerificationClass_Legacy": row.get("VerificationClass_Legacy", ""),
                    "LabelFirstSupport": row.get("LabelFirstSupport", ""),
                    "ClusterFirstSupport": row.get("ClusterFirstSupport", ""),
                    "WithinHPSupport": row.get("WithinHPSupport", ""),
                    "DispersionWarning": row.get("DispersionWarning", ""),
                    "EvidencePath": row.get("EvidencePath", ""),
                    "EvidenceDerivation": row.get("EvidenceDerivation", ""),
                }
            )

    print(f"[extract_label_first_metrics] Wrote {len(rows)} rows to {output_tsv}")


if __name__ == "__main__":
    main()
