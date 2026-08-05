#!/usr/bin/env python3
"""
X5 — Cross-sample obs18 re-computation + Thread B S3/S5 filter validation
     using KDE-corrected master from X1 Archive TO re-run.

Prerequisite: X1 batch completed (6 samples × TP+FP = 12 ISM runs).
Input:  /big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/x1_archive_to_rerun/
            {SAMPLE}_TO_{tp,fp}/significance_summary.csv
Output: research/tpfp_loh_af_kde_discrimination/data/X5_crosssample_v2.{tsv,json}
        research/tpfp_loh_af_kde_discrimination/figures/x5/ (optional)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (
    SchemaContractError,
    select_loh_legacy,
)

X1_BATCH = Path("/big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/x1_archive_to_rerun")
OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data")

SAMPLES = ["HCC1395", "HCC1395_DORADO", "H1437", "H2009", "HCC1937", "HCC1954"]
BUCKET_COLS = ["HPFineN_HP1", "HPFineN_HP1S", "HPFineN_HP2", "HPFineN_HP2S"]
NEEDED = ["RegionID", "Chr", "Pos", "AlleleDelta", "HPFineNGroups",
          "Potential_LOH", "LOH_Subtype_LegacyVC", "LOH_Subtype",
          "VerificationSchemaVersion", "Coverage_Multiple", "Diploid_Coverage_Used",
          "Coverage_Category", "NumReads"] + BUCKET_COLS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize historical inputs that only have deprecated LOH_Subtype. "
            "Versioned inputs must carry LOH_Subtype_LegacyVC and an exact LOH_Subtype alias."
        ),
    )
    return parser.parse_args()


def categorize_combo(r1, r1s, r2, r2s) -> str:
    if r1 and r1s and not r2 and not r2s: return "same_HP1"
    if r2 and r2s and not r1 and not r1s: return "same_HP2"
    if r1 and r2s and not r1s and not r2: return "cross_het"
    if r1s and r2 and not r1 and not r2s: return "cross_het_inv"
    return "other"


def af_class(af):
    if pd.isna(af): return "NA"
    a = max(0, min(1, af))
    if a < 0.1 or a > 0.9: return "Extreme"
    if 0.4 <= a <= 0.6: return "Near-half"
    return "Intermediate"


def load_sample(sample: str, allow_unversioned_v1: bool = False):
    tp_csv = X1_BATCH / f"{sample}_TO_tp" / "significance_summary.csv"
    fp_csv = X1_BATCH / f"{sample}_TO_fp" / "significance_summary.csv"
    if not tp_csv.exists() or not fp_csv.exists():
        print(f"  [!] {sample}: missing {tp_csv.name if not tp_csv.exists() else fp_csv.name}")
        return None
    dfs = []
    loh_contracts = []
    for path, lbl in [(tp_csv, 1), (fp_csv, 0)]:
        try:
            # use intersect of wanted cols (some may miss if older run)
            header = pd.read_csv(path, nrows=0).columns
            cols = [c for c in NEEDED if c in header]
            df = pd.read_csv(path, low_memory=False, usecols=cols)
            loh_view = select_loh_legacy(
                df,
                allow_unversioned_v1=allow_unversioned_v1,
            )
            df["LOH_Subtype_LegacyVC"] = loh_view.values
            loh_contracts.append({
                "input_path": str(path),
                "canonical_analysis_field": "LOH_Subtype_LegacyVC",
                "allow_unversioned_v1": allow_unversioned_v1,
                **loh_view.metadata(),
            })
            df["tp_label"] = lbl
            df["sample"] = sample
            dfs.append(df)
        except SchemaContractError:
            raise
        except Exception as e:
            print(f"  [!] {sample}/{path.name}: {e}")
            return None
    return pd.concat(dfs, ignore_index=True), loh_contracts


def obs18_split(df: pd.DataFrame) -> pd.DataFrame:
    ng2 = df[df["HPFineNGroups"] == 2].copy()
    for c in BUCKET_COLS:
        ng2[c] = ng2[c].fillna(0).astype(int)
    ng2["combo"] = ng2.apply(
        lambda r: categorize_combo(int(r["HPFineN_HP1"] > 0), int(r["HPFineN_HP1S"] > 0),
                                   int(r["HPFineN_HP2"] > 0), int(r["HPFineN_HP2S"] > 0)), axis=1)
    ng2["loh_side"] = ng2["Potential_LOH"].apply(
        lambda x: "Inner" if str(x).lower() in ("true", "1") else "Outer")
    ng2["AF_class"] = ng2["AlleleDelta"].apply(af_class)
    ext = ng2[ng2["AF_class"] == "Extreme"].copy()
    g = (ext.groupby(["sample", "loh_side", "combo"])
         .agg(n=("tp_label", "size"), n_tp=("tp_label", "sum")).reset_index())
    g["tp_rate"] = g["n_tp"] / g["n"]
    return g


def thread_b_s3_s5(df: pd.DataFrame) -> Dict:
    """Thread B S3 (Diploid Het) + S5 (combo) filter on per-sample data."""
    diploid = df[df["Coverage_Category"].isin(["Normal"])].copy()
    has_hetaf = (diploid["AlleleDelta"].between(0.4, 0.6, inclusive="both"))
    nonloh = (diploid["Potential_LOH"].astype(str).str.lower().isin(["false", "0"]))

    # S3 = Diploid Het
    s3 = diploid[has_hetaf & nonloh]
    s3_tp = s3["tp_label"].sum() if "tp_label" in s3.columns else 0
    s3_n = len(s3)
    s3_rate = s3_tp / s3_n if s3_n else None

    # S5 = S3 ∩ NG>=4 (combo)
    if "HPFineNGroups" in s3.columns:
        s5 = s3[s3["HPFineNGroups"] >= 4]
        s5_tp = s5["tp_label"].sum()
        s5_n = len(s5)
        s5_rate = s5_tp / s5_n if s5_n else None
    else:
        s5_tp = s5_n = 0; s5_rate = None

    total_tp = df["tp_label"].sum() if "tp_label" in df.columns else 0
    total_fp = (df["tp_label"] == 0).sum() if "tp_label" in df.columns else 0

    s3_fp = s3_n - s3_tp
    fp_reduction_s3 = 1 - (s3_fp / total_fp) if total_fp else None
    s5_fp = s5_n - s5_tp
    fp_reduction_s5 = 1 - (s5_fp / total_fp) if total_fp else None

    return {
        "total_tp": int(total_tp), "total_fp": int(total_fp),
        "S3_n": int(s3_n), "S3_tp_rate": s3_rate, "S3_fp_reduction": fp_reduction_s3,
        "S5_n": int(s5_n), "S5_tp_rate": s5_rate, "S5_fp_reduction": fp_reduction_s5,
    }


def main():
    args = parse_args()
    all_obs18 = []
    per_sample_summary = {}
    loh_schema_contract = {}

    for sample in SAMPLES:
        print(f"=== {sample} ===")
        loaded = load_sample(sample, allow_unversioned_v1=args.allow_unversioned_v1)
        if loaded is None: continue
        df, sample_loh_contracts = loaded
        loh_schema_contract[sample] = sample_loh_contracts
        statuses = sorted({str(item["schema_status"]) for item in sample_loh_contracts})
        print(f"  LOH schema status={','.join(statuses)}")
        dc = df["Diploid_Coverage_Used"].iloc[0] if "Diploid_Coverage_Used" in df.columns else "NA"
        print(f"  rows={len(df):,}  Diploid_Coverage_Used={dc}")

        o = obs18_split(df)
        all_obs18.append(o)

        s = thread_b_s3_s5(df)
        s["sample"] = sample
        s["diploid_cov"] = dc
        per_sample_summary[sample] = s
        print(f"  S3 n={s['S3_n']} TP rate={s['S3_tp_rate']} FP↓={s['S3_fp_reduction']}")
        print(f"  S5 n={s['S5_n']} TP rate={s['S5_tp_rate']} FP↓={s['S5_fp_reduction']}")
        print()

    if all_obs18:
        out = pd.concat(all_obs18, ignore_index=True)
        out_tsv = OUT_DIR / "X5_crosssample_obs18.tsv"
        out.to_csv(out_tsv, sep="\t", index=False)
        print(f"obs18 cross-sample → {out_tsv}")

        # Inner same_HP1 vs Outer cross_het gap per sample (replicate B1 with new data)
        pivot_gap = []
        for s in SAMPLES:
            inner = out[(out["sample"] == s) & (out["loh_side"] == "Inner") &
                        (out["combo"] == "same_HP1")]
            outer = out[(out["sample"] == s) & (out["loh_side"] == "Outer") &
                        (out["combo"] == "cross_het")]
            i_rate = inner["tp_rate"].iloc[0] if not inner.empty else None
            o_rate = outer["tp_rate"].iloc[0] if not outer.empty else None
            gap = (i_rate - o_rate) if (i_rate is not None and o_rate is not None) else None
            pivot_gap.append({"sample": s, "inner_same_HP1": i_rate,
                              "outer_cross_het": o_rate, "gap": gap})
        pivot_df = pd.DataFrame(pivot_gap)
        print("\n=== Cross-sample gap (post-KDE master) ===")
        print(pivot_df.to_string(index=False))

        gaps = [p["gap"] for p in pivot_gap if p["gap"] is not None]
        wilcoxon_result = None
        if len(gaps) >= 2:
            try:
                w = stats.wilcoxon(gaps, alternative="greater", method="exact", zero_method="wilcox")
                wilcoxon_result = {"W": float(w.statistic), "p": float(w.pvalue), "n": len(gaps)}
                print(f"\nWilcoxon W={w.statistic} p={w.pvalue:.4g} (n={len(gaps)})")
            except Exception as e:
                print(f"Wilcoxon failed: {e}")

    out_json = OUT_DIR / "X5_crosssample_summary.json"
    with open(out_json, "w") as f:
        json.dump({
            "analysis": "X5_crosssample_obs18_S3S5",
            "date": "2026-04-24",
            "source_batch": str(X1_BATCH),
            "loh_schema_contract": loh_schema_contract,
            "per_sample": per_sample_summary,
            "obs18_pivot": pivot_gap if 'pivot_gap' in locals() else [],
            "wilcoxon_replicate_B1": wilcoxon_result if 'wilcoxon_result' in locals() else None,
        }, f, indent=2)
    print(f"\nSummary JSON → {out_json}")


if __name__ == "__main__":
    main()
