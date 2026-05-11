#!/usr/bin/env python3
"""
Step 6: E1-E5 deep validation for v2 TP/FP discrimination report.

E1: SEQC2-vs-scheme TP:FP ratio fold-improvement (per sample-mode)
E2: Per-sample x scheme decomposition with Wilson CI
E3: Statistical power check (n thresholds)
E4: TP vs FP feature distribution comparisons (Mann-Whitney + Cliff's delta)
E5: Sub-scheme proposals based on E4 findings

Limitations:
  - Only HCC1395 has TO mode (other samples: paired_full only)
  - No FN data -> no true F1; report precision + TP:FP fold-improvement
  - paired_full per-sample FP is very low for H1437/H2009/HCC1954
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).parent))
from utils_io import DATA_DIR, SAMPLE_ORDER, CN_STRATEGIES, assign_cn_tier


MASTER_TSV = DATA_DIR / "merged_7samples_paired_full_plus_hcc1395_to.tsv.gz"

SCHEME_DEFS = [
    ("S1_LOH_Strong_Extreme", lambda d: (d["LOH_Subtype"] == "LOH_Strong") & (d["AF_class"] == "Extreme")),
    ("S2_Subclonal_LOH_Inter",
        lambda d: (d["LOH_Subtype"].isin(["LOH_Subclone", "LOH_Weak"])) & (d["AF_class"] == "Intermediate")),
    ("S3_Diploid_Het",
        lambda d: (d["LOH_Subtype"] == "None") & (d["AF_class"] == "Near-half")
                  & (d["cn_tier_F"].astype(str).isin(["T1", "T2"]))),
    ("S4_NonLOH_Extreme", lambda d: (d["LOH_Subtype"] == "None") & (d["AF_class"] == "Extreme")),
]


def add_combo_masks(df: pd.DataFrame) -> dict:
    """Returns dict scheme_name -> boolean Series."""
    m = {name: fn(df) for name, fn in SCHEME_DEFS}
    m["S5_Combo"] = (m["S1_LOH_Strong_Extreme"] | m["S2_Subclonal_LOH_Inter"] | m["S3_Diploid_Het"]) & ~m["S4_NonLOH_Extreme"]
    m["S6_S1_NG>=3"] = m["S1_LOH_Strong_Extreme"] & (df["HPFineNGroups"] >= 3)
    m["S7_Combo_NG>=3"] = m["S5_Combo"] & (df["HPFineNGroups"] >= 3)
    return m


def wilson_ci(k: int, n: int, z: float = 1.96) -> tuple:
    if n == 0:
        return (np.nan, np.nan, np.nan)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z / denom) * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return p, max(0.0, center - half), min(1.0, center + half)


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Cliff's delta = P(X>Y) - P(X<Y). Uses rank-based approach for speed."""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan
    # rank-based: more positive = x > y more often
    combined = np.concatenate([x, y])
    ranks = pd.Series(combined).rank().values
    rx = ranks[:nx].sum()
    # Mann-Whitney U equivalent
    U = rx - nx * (nx + 1) / 2.0
    return 2 * U / (nx * ny) - 1


def tp_fp_ratio(n_tp: int, n_fp: int) -> float:
    if n_fp == 0:
        return float("inf") if n_tp > 0 else np.nan
    return n_tp / n_fp


# ----------------------------------------------------------------------------
def e1_baseline_vs_scheme(df: pd.DataFrame, out_tsv: Path) -> pd.DataFrame:
    """Per sample-mode: baseline ratio vs each scheme ratio + fold-improvement."""
    rows = []
    for (samp, mode), g in df.groupby(["sample", "mode"], observed=True):
        baseline_tp = int(g["tp_label"].sum())
        baseline_fp = int((g["tp_label"] == 0).sum())
        baseline_ratio = tp_fp_ratio(baseline_tp, baseline_fp)
        baseline_prec = baseline_tp / len(g) if len(g) > 0 else np.nan

        masks = add_combo_masks(g)
        for scheme, mask in masks.items():
            sub = g[mask]
            n = len(sub); n_tp = int(sub["tp_label"].sum()); n_fp = n - n_tp
            scheme_ratio = tp_fp_ratio(n_tp, n_fp)
            scheme_prec = n_tp / n if n > 0 else np.nan
            fold = (scheme_ratio / baseline_ratio) if (not np.isinf(baseline_ratio) and baseline_ratio > 0 and not np.isnan(scheme_ratio)) else np.nan
            p, lo, hi = wilson_ci(n_tp, n)
            rows.append(dict(
                sample=samp, mode=mode,
                baseline_n=len(g), baseline_tp=baseline_tp, baseline_fp=baseline_fp,
                baseline_ratio=baseline_ratio, baseline_precision=baseline_prec,
                scheme=scheme, n=n, n_tp=n_tp, n_fp=n_fp,
                scheme_ratio=scheme_ratio, scheme_precision=scheme_prec,
                ci_lo=lo, ci_hi=hi, ci_width=hi - lo if not (np.isnan(lo) or np.isnan(hi)) else np.nan,
                fold_improvement=fold,
                tp_recall=n_tp / baseline_tp if baseline_tp > 0 else np.nan,
                fp_recall=n_fp / baseline_fp if baseline_fp > 0 else np.nan,
            ))
    out = pd.DataFrame(rows)
    out.to_csv(out_tsv, sep="\t", index=False)
    return out


# ----------------------------------------------------------------------------
def e3_power_tags(e1_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate each row with statistical power tag."""
    def tag(n):
        if n < 20:
            return "NO_USE"
        elif n < 50:
            return "WEAK"
        elif n < 200:
            return "OK"
        else:
            return "STRONG"
    out = e1_df.copy()
    out["power_tag"] = out["n"].apply(tag)
    return out


# ----------------------------------------------------------------------------
def e4_feature_diffs(df: pd.DataFrame, out_tsv: Path) -> pd.DataFrame:
    """For each scheme, compare TP vs FP distributions of key features."""
    features = ["NumReads", "NumCpGs", "AlleleDelta", "CovM_used", "HPFineNGroups"]
    rows = []

    for sample_mode, g in df.groupby(["sample", "mode"], observed=True):
        samp, mode = sample_mode
        masks = add_combo_masks(g)
        for scheme, mask in masks.items():
            sub = g[mask]
            tp = sub[sub["tp_label"] == 1]
            fp = sub[sub["tp_label"] == 0]
            if len(tp) < 10 or len(fp) < 10:
                continue
            for feat in features:
                if feat not in sub.columns:
                    continue
                tp_vals = tp[feat].dropna().values
                fp_vals = fp[feat].dropna().values
                if len(tp_vals) < 10 or len(fp_vals) < 10:
                    continue
                try:
                    u, p = mannwhitneyu(tp_vals, fp_vals, alternative="two-sided")
                except ValueError:
                    u, p = np.nan, np.nan
                d = cliffs_delta(tp_vals, fp_vals)
                rows.append(dict(
                    sample=samp, mode=mode, scheme=scheme, feature=feat,
                    n_tp=len(tp_vals), n_fp=len(fp_vals),
                    median_tp=float(np.median(tp_vals)),
                    median_fp=float(np.median(fp_vals)),
                    q25_tp=float(np.quantile(tp_vals, 0.25)),
                    q75_tp=float(np.quantile(tp_vals, 0.75)),
                    q25_fp=float(np.quantile(fp_vals, 0.25)),
                    q75_fp=float(np.quantile(fp_vals, 0.75)),
                    mannwhitney_u=float(u) if not np.isnan(u) else np.nan,
                    mannwhitney_p=float(p) if not np.isnan(p) else np.nan,
                    cliffs_delta=float(d) if not np.isnan(d) else np.nan,
                ))
    out = pd.DataFrame(rows)
    out.to_csv(out_tsv, sep="\t", index=False)
    return out


# ----------------------------------------------------------------------------
def e5_subscheme_proposal(df: pd.DataFrame, feature_diffs: pd.DataFrame, out_tsv: Path) -> pd.DataFrame:
    """Apply candidate sub-schemes based on E4 findings."""
    rows = []

    for sample_mode, g in df.groupby(["sample", "mode"], observed=True):
        samp, mode = sample_mode
        if samp == "HCC1395" and mode == "to_pileup":
            label = "HCC1395_TO"
        else:
            label = f"{samp}_{mode}"
        masks = add_combo_masks(g)
        m_s4 = masks["S4_NonLOH_Extreme"]
        m_s5 = masks["S5_Combo"]
        m_s1 = masks["S1_LOH_Strong_Extreme"]

        def reg(name, mask, description):
            sub = g[mask]
            n = len(sub); n_tp = int(sub["tp_label"].sum()); n_fp = n - n_tp
            baseline_tp = int(g["tp_label"].sum()); baseline_fp = int((g["tp_label"]==0).sum())
            baseline_ratio = tp_fp_ratio(baseline_tp, baseline_fp)
            scheme_ratio = tp_fp_ratio(n_tp, n_fp)
            fold = (scheme_ratio/baseline_ratio) if (not np.isinf(baseline_ratio) and baseline_ratio>0 and not np.isnan(scheme_ratio)) else np.nan
            p, lo, hi = wilson_ci(n_tp, n)
            rows.append(dict(
                sample=samp, mode=mode, scheme_id=name, description=description,
                n=n, n_tp=n_tp, n_fp=n_fp,
                tp_rate=p, ci_lo=lo, ci_hi=hi,
                tp_fp_ratio=scheme_ratio,
                baseline_ratio=baseline_ratio, fold_improvement=fold,
                tp_recall=n_tp/baseline_tp if baseline_tp>0 else np.nan,
            ))

        nr = g["NumReads"]; ng = g["HPFineNGroups"]; cov = g["CovM_used"]

        # S4 refinements
        reg("S4a_S4+NR>=80", m_s4 & (nr >= 80), "S4 + NumReads>=80 (high coverage)")
        reg("S4b_S4+CovM_in_diploid", m_s4 & (cov >= 0.7) & (cov <= 1.4),
            "S4 + CovM in diploid range")
        reg("S4c_S4+NG=4", m_s4 & (ng == 4), "S4 + NG=4 only (subclone marker)")
        reg("S4d_S4+NR>=80+NG>=3", m_s4 & (nr >= 80) & (ng >= 3),
            "S4 combined refinements")

        # S5 refinements
        reg("S5a_S5+NG=4", m_s5 & (ng == 4), "S5 + NG=4 (strictest union)")
        reg("S5b_S5+NR>=80", m_s5 & (nr >= 80), "S5 + NumReads>=80")

        # Strictest
        reg("Sx_S1orS3+NG=4+NR>=80",
            (m_s1 | masks["S3_Diploid_Het"]) & (ng == 4) & (nr >= 80),
            "LOH_Strong+Extreme OR DiploidHet + NG=4 + NR>=80")

    out = pd.DataFrame(rows)
    out.to_csv(out_tsv, sep="\t", index=False)
    return out


# ----------------------------------------------------------------------------
def main() -> None:
    print(f"[step6] loading master")
    df = pd.read_csv(MASTER_TSV, sep="\t", compression="gzip", low_memory=False)
    df["LOH_Subtype"] = df["LOH_Subtype"].fillna("None")
    df["cn_tier_F"] = assign_cn_tier(df["CovM_used"], CN_STRATEGIES["F_SEQC2_grounded"])
    print(f"[step6] shape: {df.shape}")

    # E1 + E3
    print("\n[E1] computing baseline vs scheme TP:FP ratios...")
    e1 = e1_baseline_vs_scheme(df, DATA_DIR / "tpfp_baseline_ratio.tsv")
    e1 = e3_power_tags(e1)
    e1.to_csv(DATA_DIR / "tpfp_per_sample_scheme_full.tsv", sep="\t", index=False)
    print(f"[E1+E3] wrote tpfp_baseline_ratio.tsv + tpfp_per_sample_scheme_full.tsv: {e1.shape}")

    # Print highlight table for HCC1395 TO and COLO829 paired_full
    print("\n[E1 highlight] HCC1395 TO mode:")
    sub = e1[(e1["sample"] == "HCC1395") & (e1["mode"] == "to_pileup")]
    print(sub[["scheme","n","n_tp","n_fp","scheme_precision","scheme_ratio","fold_improvement","tp_recall","power_tag"]].to_string(index=False))

    print("\n[E1 highlight] COLO829 paired_full:")
    sub = e1[(e1["sample"] == "COLO829") & (e1["mode"] == "paired_full")]
    print(sub[["scheme","n","n_tp","n_fp","scheme_precision","scheme_ratio","fold_improvement","tp_recall","power_tag"]].to_string(index=False))

    # E3 summary
    print("\n[E3] power tag distribution:")
    print(e1.groupby(["power_tag"]).size().to_string())

    # E4
    print("\n[E4] computing TP vs FP feature diffs...")
    e4 = e4_feature_diffs(df, DATA_DIR / "tpfp_feature_diffs.tsv")
    print(f"[E4] wrote tpfp_feature_diffs.tsv: {e4.shape}")

    # Print top feature gaps (HCC1395 TO, by |cliffs_delta|)
    print("\n[E4 highlight] HCC1395 TO feature gaps (top |delta|):")
    sub = e4[(e4["sample"] == "HCC1395") & (e4["mode"] == "to_pileup")].copy()
    sub["abs_delta"] = sub["cliffs_delta"].abs()
    print(sub.sort_values("abs_delta", ascending=False).head(15)[
        ["scheme","feature","n_tp","n_fp","median_tp","median_fp","mannwhitney_p","cliffs_delta"]
    ].to_string(index=False))

    # E5
    print("\n[E5] running sub-scheme proposals...")
    e5 = e5_subscheme_proposal(df, e4, DATA_DIR / "tpfp_subschemes.tsv")
    print(f"[E5] wrote tpfp_subschemes.tsv: {e5.shape}")

    print("\n[E5 highlight] HCC1395 TO sub-schemes:")
    sub = e5[(e5["sample"] == "HCC1395") & (e5["mode"] == "to_pileup")]
    print(sub[["scheme_id","n","n_tp","n_fp","tp_rate","tp_fp_ratio","fold_improvement","tp_recall"]].to_string(index=False))

    print("\n[E5 highlight] COLO829 paired_full sub-schemes:")
    sub = e5[(e5["sample"] == "COLO829") & (e5["mode"] == "paired_full")]
    print(sub[["scheme_id","n","n_tp","n_fp","tp_rate","tp_fp_ratio","fold_improvement","tp_recall"]].to_string(index=False))

    print("\n[step6] DONE")


if __name__ == "__main__":
    main()
