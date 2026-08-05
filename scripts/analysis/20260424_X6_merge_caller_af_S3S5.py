#!/usr/bin/env python3
"""
X6 — Merge caller AF (from VCF FORMAT:AF field) into ISM TSV,
     then rerun Thread B S3/S5 filter validation across 6 TO samples.

Fixes X5 v2 limitation: `AlleleDelta` ≠ caller AF; Thread B v2 report used VCF caller AF.
"""
from __future__ import annotations
import argparse
import gzip
import json
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (
    SchemaContractError,
    VERIFICATION_PROVENANCE_COLUMNS,
    select_loh_legacy,
)

X1_BATCH = Path("/big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/x1_archive_to_rerun")
OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data")

SAMPLES = {
    "HCC1395": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot",
    "HCC1395_DORADO": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot",
    "H1437": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h1437_to_pilot_fastresume",
    "H2009": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_h2009_to_pilot_fastresume",
    "HCC1937": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1937_to_pilot_fastresume",
    "HCC1954": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/20260318_hcc1954_to_pilot_fastresume",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize historical unversioned LOH_Subtype input; all "
            "available verification provenance columns are retained through the AF merge"
        ),
    )
    return parser.parse_args()


def read_summary_with_loh_provenance(
    path: Path,
    allow_unversioned_v1: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """Read required features while retaining all available schema provenance."""
    base_columns = [
        "Chr", "Pos", "Ref", "Alt", "HPFineNGroups",
        "Coverage_Multiple", "Diploid_Coverage_Used",
    ]
    available = set(pd.read_csv(path, nrows=0).columns)
    provenance_columns = [
        column for column in VERIFICATION_PROVENANCE_COLUMNS if column in available
    ]
    frame = pd.read_csv(
        path,
        low_memory=False,
        usecols=base_columns + provenance_columns,
    )
    view = select_loh_legacy(frame, allow_unversioned_v1=allow_unversioned_v1)
    frame["_loh_subtype_legacy"] = view.values
    metadata = {
        "input_path": str(path),
        "selection_field": view.field,
        "schema_status": view.schema_status,
        "allow_unversioned_v1": allow_unversioned_v1,
        "provenance_columns_passthrough": provenance_columns,
        "warnings": list(view.warning_messages),
        "row_count": len(frame),
    }
    return frame, metadata


def parse_vcf_af(vcf_path: str) -> pd.DataFrame:
    """Parse VCF FORMAT:AF to DataFrame (Chr, Pos, Ref, Alt, caller_AF)."""
    rows = []
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, _, ref, alt, _, _, _, fmt, sample = parts[:10]
            fmt_fields = fmt.split(":")
            try:
                af_idx = fmt_fields.index("AF")
            except ValueError:
                continue
            sample_vals = sample.split(":")
            try:
                af = float(sample_vals[af_idx])
            except (ValueError, IndexError):
                continue
            rows.append((chrom, int(pos), ref, alt, af))
    return pd.DataFrame(rows, columns=["Chr", "Pos", "Ref", "Alt", "caller_AF"])


def load_and_merge(
    sample: str,
    archive: str,
    allow_unversioned_v1: bool = False,
):
    tp_tsv = X1_BATCH / f"{sample}_TO_tp" / "significance_summary.csv"
    fp_tsv = X1_BATCH / f"{sample}_TO_fp" / "significance_summary.csv"
    tp_vcf = f"{archive}/step02_benchmark_clairs_to/filtered_snv_tp.vcf.gz"
    fp_vcf = f"{archive}/step02_benchmark_clairs_to/filtered_snv_fp.vcf.gz"

    dfs = []
    schema_rows = []
    for tsv_path, vcf_path, lbl in [(tp_tsv, tp_vcf, 1), (fp_tsv, fp_vcf, 0)]:
        tsv, schema_metadata = read_summary_with_loh_provenance(
            tsv_path,
            allow_unversioned_v1=allow_unversioned_v1,
        )
        schema_metadata["truth_label"] = "TP" if lbl == 1 else "FP"
        vcf_af = parse_vcf_af(vcf_path)
        merged = tsv.merge(vcf_af, on=["Chr", "Pos", "Ref", "Alt"], how="left")
        missing_provenance = set(schema_metadata["provenance_columns_passthrough"]) - set(
            merged.columns
        )
        if missing_provenance:
            raise SchemaContractError(
                "X6 caller-AF merge dropped provenance columns: "
                + ", ".join(sorted(missing_provenance))
            )
        merged["tp_label"] = lbl
        dfs.append(merged)
        schema_rows.append(schema_metadata)
    schema_modes = {(row["selection_field"], row["schema_status"]) for row in schema_rows}
    if len(schema_modes) != 1:
        raise SchemaContractError(
            f"X6 TP/FP inputs use mixed LOH schema modes: {sorted(schema_modes)}"
        )
    df = pd.concat(dfs, ignore_index=True)
    missing = df["caller_AF"].isna().sum()
    return df, missing, schema_rows


def af_class(af):
    if pd.isna(af): return "NA"
    if af < 0.1 or af > 0.9: return "Extreme"
    if 0.4 <= af <= 0.6: return "Near-half"
    return "Intermediate"


def cn_tier(cm):
    if pd.isna(cm): return "NA"
    if cm < 0.65: return "T0"
    if cm < 0.99: return "T1"
    if cm < 1.33: return "T2"
    if cm < 1.82: return "T3"
    return "T4"


def classify(df):
    df = df.copy()
    df["AF_class"] = df["caller_AF"].apply(af_class)
    df["CN_tier"] = df["Coverage_Multiple"].apply(cn_tier)
    df["LOH_class"] = df["_loh_subtype_legacy"].astype(str)
    df["is_S1"] = (df["LOH_class"] == "LOH_Strong") & (df["AF_class"] == "Extreme")
    df["is_S2"] = df["LOH_class"].isin(["LOH_Weak", "LOH_Subclone"])
    df["is_S3"] = (df["LOH_class"] == "None") & (df["AF_class"] == "Near-half") & df["CN_tier"].isin(["T1", "T2"])
    df["is_S4"] = (df["LOH_class"] == "None") & (df["AF_class"] == "Extreme")
    df["is_S5"] = (df["is_S1"] | df["is_S2"] | df["is_S3"]) & (~df["is_S4"])
    df["is_S6"] = df["is_S1"] & (df["HPFineNGroups"].fillna(0) >= 3)
    df["is_S7"] = df["is_S5"] & (df["HPFineNGroups"].fillna(0) >= 3)
    return df


def scheme_stats(df, flag):
    sub = df[df[flag]]
    n = len(sub); tp = int(sub["tp_label"].sum()); fp = n - tp
    total_tp = int((df["tp_label"] == 1).sum()); total_fp = int((df["tp_label"] == 0).sum())
    return {
        "n": n, "tp": tp, "fp": fp,
        "tp_rate": (tp / n) if n else None,
        "fp_reduction": (1 - fp / total_fp) if total_fp else None,
        "tp_recall": (tp / total_tp) if total_tp else None,
    }


def main():
    args = parse_args()
    all_rows = []
    loh_schema_provenance = []
    for sample, archive in SAMPLES.items():
        print(f"=== {sample} ===")
        df, missing, schema_rows = load_and_merge(
            sample,
            archive,
            allow_unversioned_v1=args.allow_unversioned_v1,
        )
        for schema_row in schema_rows:
            loh_schema_provenance.append({"sample": sample, **schema_row})
        dc = df["Diploid_Coverage_Used"].iloc[0]
        print(f"  total={len(df):,}  DC={dc}  AF missing={missing}")
        print(f"  caller_AF stats: {df['caller_AF'].describe()[['count','mean','50%']].to_dict()}")
        df = classify(df)
        af_dist = df["AF_class"].value_counts().to_dict()
        print(f"  AF_class: {af_dist}")

        entry = {"sample": sample, "diploid_cov": float(dc), "rows": len(df),
                 "caller_af_missing": int(missing),
                 "total_tp": int((df["tp_label"] == 1).sum()),
                 "total_fp": int((df["tp_label"] == 0).sum())}
        for sc in ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]:
            s = scheme_stats(df, f"is_{sc}")
            entry[sc] = s
            mark = "⭐" if s["tp_rate"] and s["tp_rate"] >= 0.90 and s["n"] >= 20 else " "
            print(f"  {mark}{sc}  n={s['n']:>6} tp={s['tp_rate']} fp_red={s['fp_reduction']}")
        all_rows.append(entry)
        print()

    # Cross-sample S3 stability
    print("=== Cross-sample stability (S3 / S5) ===")
    s3_ok = sum(1 for r in all_rows if r["S3"]["tp_rate"] and r["S3"]["tp_rate"] >= 0.85 and r["S3"]["n"] >= 20)
    s5_ok = sum(1 for r in all_rows if r["S5"]["tp_rate"] and r["S5"]["tp_rate"] >= 0.85 and r["S5"]["n"] >= 50)
    print(f"  S3 TP≥0.85 & n≥20 : {s3_ok}/6 samples")
    print(f"  S5 TP≥0.85 & n≥50 : {s5_ok}/6 samples")

    # Wilcoxon S3 vs baseline
    wilc = None
    pairs = [(r["S3"]["tp_rate"], r["total_tp"] / (r["total_tp"] + r["total_fp"]))
             for r in all_rows if r["S3"]["tp_rate"] is not None and r["S3"]["n"] >= 20]
    if len(pairs) >= 3:
        try:
            w = stats.wilcoxon([p[0] for p in pairs], [p[1] for p in pairs],
                               alternative="greater", method="exact", zero_method="wilcox")
            wilc = {"W": float(w.statistic), "p": float(w.pvalue), "n": len(pairs)}
            print(f"  Wilcoxon S3 > baseline: W={w.statistic} p={w.pvalue:.4g} n={len(pairs)}")
        except Exception as e:
            print(f"  Wilcoxon failed: {e}")

    out = OUT_DIR / "X6_caller_AF_S3S5.json"
    with open(out, "w") as f:
        json.dump({
            "analysis": "X6_caller_AF_merged_S3S5",
            "date": "2026-04-24",
            "af_source": "VCF FORMAT:AF (ClairS-TO filtered)",
            "loh_schema_provenance": loh_schema_provenance,
            "per_sample": all_rows,
            "cross_sample_stability": {"s3_pass_4of6": s3_ok, "s5_pass_4of6": s5_ok},
            "s3_vs_baseline_wilcoxon": wilc,
        }, f, indent=2)
    print(f"\nJSON → {out}")


if __name__ == "__main__":
    main()
