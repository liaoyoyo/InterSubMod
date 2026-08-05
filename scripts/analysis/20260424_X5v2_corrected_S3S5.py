#!/usr/bin/env python3
"""
X5 v2 — Corrected S3/S5 definition per Thread B v2 report:
  S3 Diploid Het = LOH_Subtype_LegacyVC=None AND AF_class=Near-half AND CN_tier in (T1, T2)
  S5 Combo = (S1 ∨ S2 ∨ S3) \\ S4

CN tiers (post-KDE, SEQC2 empirical):
  T0: CovM < 0.65 (CN Loss)
  T1: 0.65 ≤ CovM < 0.99 (CN=1-2)
  T2: 0.99 ≤ CovM < 1.33 (CN=2, diploid)
  T3: 1.33 ≤ CovM < 1.82 (CN=3)
  T4: CovM ≥ 1.82 (CN≥4)

AF classes: Extreme (<0.1 or >0.9), Near-half (0.4-0.6), Intermediate (else)

S1 LOH Strong + Extreme AF
S2 LOH Weak/Subclone + any AF
S3 None + Near-half + T1/T2
S4 None + Extreme AF (ambiguous)
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import select_loh_legacy

X1_BATCH = Path("/big7_disk/liaoyoyo2001/big7_disk_output/kde_smoke_test/x1_archive_to_rerun")
OUT_DIR = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tpfp_loh_af_kde_discrimination/data")
SAMPLES = ["HCC1395", "HCC1395_DORADO", "H1437", "H2009", "HCC1937", "HCC1954"]
LOH_SCHEMA_COLUMNS = [
    "VerificationSchemaVersion",
    "LOH_Subtype_LegacyVC",
    "LOH_Subtype",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-unversioned-v1",
        action="store_true",
        help=(
            "Explicitly authorize historical unversioned LOH_Subtype input; "
            "versioned input requires canonical LOH_Subtype_LegacyVC plus exact alias"
        ),
    )
    return parser.parse_args()


def attach_loh_legacy_view(
    frame: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> pd.DataFrame:
    view = select_loh_legacy(frame, allow_unversioned_v1=allow_unversioned_v1)
    prepared = frame.copy()
    prepared["_loh_subtype_legacy"] = view.values
    prepared.attrs["loh_schema_contract"] = {
        "selection_field": view.field,
        "schema_status": view.schema_status,
        "allow_unversioned_v1": allow_unversioned_v1,
        "warnings": list(view.warning_messages),
    }
    return prepared


def read_sample_summary(path: Path) -> pd.DataFrame:
    """Read analysis fields plus every available LOH schema identity column."""
    base_columns = [
        "RegionID", "AlleleDelta", "HPFineNGroups", "Coverage_Multiple",
        "Diploid_Coverage_Used",
    ]
    available = set(pd.read_csv(path, nrows=0).columns)
    usecols = base_columns + [column for column in LOH_SCHEMA_COLUMNS if column in available]
    return pd.read_csv(path, low_memory=False, usecols=usecols)


def af_class(af):
    if pd.isna(af): return "NA"
    a = max(0, min(1, af))
    if a < 0.1 or a > 0.9: return "Extreme"
    if 0.4 <= a <= 0.6: return "Near-half"
    return "Intermediate"


def cn_tier(cm):
    if pd.isna(cm): return "NA"
    if cm < 0.65: return "T0"
    if cm < 0.99: return "T1"
    if cm < 1.33: return "T2"
    if cm < 1.82: return "T3"
    return "T4"


def load_sample(sample: str, allow_unversioned_v1: bool = False):
    tp = X1_BATCH / f"{sample}_TO_tp" / "significance_summary.csv"
    fp = X1_BATCH / f"{sample}_TO_fp" / "significance_summary.csv"
    dfs = []
    for path, lbl in [(tp, 1), (fp, 0)]:
        df = read_sample_summary(path)
        df["tp_label"] = lbl
        df["sample"] = sample
        dfs.append(df)
    return attach_loh_legacy_view(
        pd.concat(dfs, ignore_index=True),
        allow_unversioned_v1=allow_unversioned_v1,
    )


def classify(df):
    df = df.copy()
    df["AF_class"] = df["AlleleDelta"].apply(af_class)
    df["CN_tier"] = df["Coverage_Multiple"].apply(cn_tier)
    df["LOH_class"] = df["_loh_subtype_legacy"].astype(str)

    # S1 LOH Strong + Extreme
    df["is_S1"] = (df["LOH_class"] == "LOH_Strong") & (df["AF_class"] == "Extreme")
    # S2 LOH Weak/Subclone (any AF)
    df["is_S2"] = df["LOH_class"].isin(["LOH_Weak", "LOH_Subclone"])
    # S3 Diploid Het
    df["is_S3"] = (df["LOH_class"] == "None") & (df["AF_class"] == "Near-half") & df["CN_tier"].isin(["T1", "T2"])
    # S4 None + Extreme (ambiguous)
    df["is_S4"] = (df["LOH_class"] == "None") & (df["AF_class"] == "Extreme")
    # S5 (S1 ∨ S2 ∨ S3) \ S4
    df["is_S5"] = (df["is_S1"] | df["is_S2"] | df["is_S3"]) & (~df["is_S4"])
    # S6 S1 + NG≥3
    df["is_S6"] = df["is_S1"] & (df["HPFineNGroups"].fillna(0) >= 3)
    # S7 S5 + NG≥3
    df["is_S7"] = df["is_S5"] & (df["HPFineNGroups"].fillna(0) >= 3)
    return df


def stats_on(df_sample, flag_col):
    subset = df_sample[df_sample[flag_col]]
    n = len(subset)
    tp = subset["tp_label"].sum()
    fp = n - tp
    tp_rate = tp / n if n else None
    total_fp = (df_sample["tp_label"] == 0).sum()
    fp_reduction = 1 - (fp / total_fp) if total_fp else None
    total_tp = (df_sample["tp_label"] == 1).sum()
    tp_recall = tp / total_tp if total_tp else None
    return {"n": int(n), "tp": int(tp), "fp": int(fp),
            "tp_rate": tp_rate, "fp_reduction": fp_reduction, "tp_recall": tp_recall}


def main():
    args = parse_args()
    all_rows = []
    for sample in SAMPLES:
        print(f"=== {sample} ===")
        df = load_sample(sample, allow_unversioned_v1=args.allow_unversioned_v1)
        loh_schema_contract = dict(df.attrs["loh_schema_contract"])
        df = classify(df)
        dc = df["Diploid_Coverage_Used"].iloc[0]

        entry = {"sample": sample, "diploid_cov": float(dc),
                 "loh_schema_contract": loh_schema_contract,
                 "total_tp": int((df["tp_label"] == 1).sum()),
                 "total_fp": int((df["tp_label"] == 0).sum())}
        for scheme in ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]:
            s = stats_on(df, f"is_{scheme}")
            entry[scheme] = s
            mark = "⭐" if s["tp_rate"] and s["tp_rate"] >= 0.90 and s["n"] >= 20 else " "
            print(f"  {mark}{scheme}  n={s['n']:>6} tp_rate={s['tp_rate']} fp_red={s['fp_reduction']} tp_rec={s['tp_recall']}")
        all_rows.append(entry)
        print()

    # cross-sample consistency on S3
    print("=== Cross-sample S3 Diploid Het ===")
    s3_stable = sum(1 for r in all_rows if r["S3"]["tp_rate"] and r["S3"]["tp_rate"] >= 0.85 and r["S3"]["n"] >= 20)
    print(f"  S3 TP rate ≥0.85 & n≥20: {s3_stable}/6 samples")

    print("\n=== Cross-sample S5 Combo ===")
    s5_stable = sum(1 for r in all_rows if r["S5"]["tp_rate"] and r["S5"]["tp_rate"] >= 0.85 and r["S5"]["n"] >= 50)
    print(f"  S5 TP rate ≥0.85 & n≥50: {s5_stable}/6 samples")

    # S3 TP rate list for stat test
    s3_tp_rates = [r["S3"]["tp_rate"] for r in all_rows if r["S3"]["tp_rate"] is not None and r["S3"]["n"] >= 20]
    baseline_tp = [r["total_tp"] / (r["total_tp"] + r["total_fp"]) for r in all_rows]

    # Wilcoxon: S3 TP rate vs overall baseline
    wilcoxon_res = None
    if len(s3_tp_rates) >= 3:
        # pair S3 TP with sample baseline
        pairs = [(r["S3"]["tp_rate"], r["total_tp"] / (r["total_tp"] + r["total_fp"]))
                 for r in all_rows if r["S3"]["tp_rate"] is not None and r["S3"]["n"] >= 20]
        s3_rates = [p[0] for p in pairs]
        baseline_rates = [p[1] for p in pairs]
        try:
            w = stats.wilcoxon(s3_rates, baseline_rates, alternative="greater", method="exact", zero_method="wilcox")
            wilcoxon_res = {"W": float(w.statistic), "p": float(w.pvalue), "n": len(pairs)}
            print(f"\nS3 TP rate vs baseline Wilcoxon W={w.statistic} p={w.pvalue:.4g} n={len(pairs)}")
        except Exception as e:
            print(f"Wilcoxon failed: {e}")

    out = OUT_DIR / "X5v2_corrected_S3S5.json"
    with open(out, "w") as f:
        json.dump({
            "analysis": "X5v2_corrected_S3S5",
            "date": "2026-04-24",
            "scheme_definitions": {
                "S1": "LOH_Subtype_LegacyVC=LOH_Strong AND AF=Extreme",
                "S2": "LOH_Subtype_LegacyVC in (LOH_Weak, LOH_Subclone)",
                "S3": "LOH_Subtype_LegacyVC=None AND AF=Near-half AND CN_tier in (T1,T2)",
                "S4": "LOH_Subtype_LegacyVC=None AND AF=Extreme (ambiguous)",
                "S5": "(S1∨S2∨S3) \\\\ S4",
                "S6": "S1 + NG≥3", "S7": "S5 + NG≥3",
            },
            "per_sample": all_rows,
            "cross_sample": {
                "s3_stable_4of6": s3_stable,
                "s5_stable_4of6": s5_stable,
                "s3_vs_baseline_wilcoxon": wilcoxon_res,
            },
        }, f, indent=2)
    print(f"\nJSON → {out}")


if __name__ == "__main__":
    main()
