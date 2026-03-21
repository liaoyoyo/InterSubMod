#!/usr/bin/env python3
"""Extract explicit label-first metrics from significance_summary.csv."""

import argparse
import csv
import os


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


def main():
    args = parse_args()
    output_tsv = args.output_tsv
    if output_tsv is None:
        output_tsv = os.path.join(os.path.dirname(args.summary_csv), "label_first_metrics.tsv")

    with open(args.summary_csv, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

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
                }
            )

    print(f"[extract_label_first_metrics] Wrote {len(rows)} rows to {output_tsv}")


if __name__ == "__main__":
    main()
