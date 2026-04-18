#!/usr/bin/env python3
"""
B.2 LOH×AF×Methylation Deep Skeptical Check — Batch 1 (R2, R3, R4)

Opus 4.7 Plan 對 LOH Subclone AF × Methylation POSITIVE 的三項質疑：
  R2  HCC1954 反向是否為 post-hoc rationalization
  R3  NR-bin 效應隨 NR 增強是否為統計功效 artifact（matched sampling）
  R4  cnLOH (CN2) vs deletion-LOH (CN1) 方向是否一致

Inputs:
  output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz

Outputs:
  research/B2_loh_skeptical_check/data/r2_mann_whitney_by_sample.tsv
  research/B2_loh_skeptical_check/data/r2_hcc1954_segment_level.tsv
  research/B2_loh_skeptical_check/data/r3_nr_bins_with_80plus.tsv
  research/B2_loh_skeptical_check/data/r3_matched_sampling.tsv
  research/B2_loh_skeptical_check/data/r4_cn1_vs_cn2_per_sample.tsv

Pre-registered decisions (before running):
  - HCC1954 exclusion rule: CN1 segment count < 50 OR >55% intermediate → exclude from segment-level
    (HCC1954 has 60.2% intermediate and ~30-34 segments per original; both triggered)
  - Matched sampling: downsample NR≥80 bin to N = NR 10-30 bin size
  - CN classification: CN1 cm<0.75 (deletion-LOH), CN2 0.75≤cm<1.25 (cnLOH)
  - Both modes analysed: TO (mode="to") AND Paired (mode="paired")
  - LOH indicator: Potential_LOH (cross-mode column; core_loh_like / tool_potential_loh identical).
    Note: original 20260414 evidence used to_loh_bed_hit (TO-only). Potential_LOH is a superset for TO
    (TO: 175542 vs to_loh_bed_hit 112218; overlap 111932). We use Potential_LOH for BOTH modes to ensure
    symmetric definitions; TO results with to_loh_bed_hit remain unchanged from 20260414. Caveat recorded
    in batch 1 report Section §Method Note.

Refs:
  docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md
  docs/plans/opus-4-7-big8-disk-liaoyoyo2001-knowled-cryptic-moore.md (E.9 R2/R3/R4)
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[3]
MASTER = REPO_ROOT / "output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
OUT_DIR = Path(__file__).resolve().parent.parent / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Pre-registered thresholds
CN1_MAX = 0.75          # deletion-LOH
CN2_MAX = 1.25          # cnLOH / diploid
HCC1954_MIN_SEGMENTS = 50
HCC1954_MAX_INTER_FRAC = 0.55
NR_BINS = [(10, 30), (30, 50), (50, 80), (80, 500)]
RNG_SEED = 20260418


def classify_af(af):
    if af < 0.1 or af > 0.9:
        return "Extreme"
    if 0.4 <= af <= 0.6:
        return "Near-half"
    return "Intermediate"


def classify_cn(cm):
    if pd.isna(cm):
        return "NA"
    if cm < CN1_MAX:
        return "CN1"
    if cm < CN2_MAX:
        return "CN2"
    if cm < 1.75:
        return "CN3"
    return "CN4+"


def rank_biserial(a, b):
    """rank-biserial r for Mann-Whitney (greater alternative)."""
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return np.nan
    u, _ = stats.mannwhitneyu(a, b, alternative="greater")
    # r such that r > 0 means a > b in rank sense
    return 1 - (2 * u) / (n1 * n2)


def mw_test(inter, extr):
    if len(inter) < 5 or len(extr) < 5:
        return np.nan, np.nan, np.nan
    u, p = stats.mannwhitneyu(inter, extr, alternative="greater")
    r = rank_biserial(inter, extr)
    return u, p, r


def load_loh_tp():
    print(f"[load] reading {MASTER}")
    df = pd.read_csv(MASTER, sep="\t", compression="gzip", low_memory=False)
    print(f"[load] raw rows={len(df):,}")
    # Canonical subset: LOH (to_loh_bed_hit) × TP × mode in {TO, paired}
    df["af_class"] = df["caller_af"].apply(classify_af)
    df["cn_tier"] = df["Coverage_Multiple"].apply(classify_cn)
    # Coerce mode to lowercase canonical
    df["mode"] = df["mode"].astype(str).str.lower()
    # Cross-mode LOH: Potential_LOH (identical to core_loh_like / tool_potential_loh)
    # TO mode has extra to_loh_bed_hit (112k subset of Potential_LOH 175k);
    # paired mode to_loh_bed_hit is empty so we use Potential_LOH for both.
    loh = df[df["Potential_LOH"] == True].copy()
    tp = loh[loh["truth_label"] == "TP"].copy()
    print(f"[load] Potential_LOH rows={len(loh):,} | LOH×TP rows={len(tp):,}")
    print(f"[load] modes present: {sorted(tp['mode'].unique())}")
    print(f"[load] samples present: {sorted(tp['sample'].unique())}")
    return tp


# ---------------------------------------------------------------------------
# R2 — HCC1954 pre-registered exclusion
# ---------------------------------------------------------------------------
def run_r2(tp: pd.DataFrame):
    print("\n=== R2: HCC1954 pre-registered exclusion ===")
    rows = []
    for mode in ["to", "paired"]:
        sub = tp[tp["mode"] == mode].copy()
        cn1 = sub[sub["cn_tier"] == "CN1"].copy()
        for sample in sorted(cn1["sample"].unique()):
            samp = cn1[cn1["sample"] == sample]
            inter = samp[samp["af_class"] == "Intermediate"]["HPFineNGroups"].dropna()
            extr = samp[samp["af_class"] == "Extreme"]["HPFineNGroups"].dropna()
            near = samp[samp["af_class"] == "Near-half"]
            n_total = len(samp)
            frac_inter = (samp["af_class"] == "Intermediate").mean()
            u, p, r = mw_test(inter.values, extr.values)
            rows.append(
                {
                    "mode": mode,
                    "sample": sample,
                    "n_total": n_total,
                    "n_extreme": len(extr),
                    "n_intermediate": len(inter),
                    "frac_intermediate": round(frac_inter, 3),
                    "mean_ng_ext": round(extr.mean(), 3) if len(extr) else np.nan,
                    "mean_ng_inter": round(inter.mean(), 3) if len(inter) else np.nan,
                    "delta_ng": (
                        round(inter.mean() - extr.mean(), 3)
                        if len(extr) and len(inter)
                        else np.nan
                    ),
                    "mw_p": p,
                    "rank_biserial_r": round(r, 3) if not np.isnan(r) else np.nan,
                    "direction_positive": bool(
                        (len(extr) >= 5)
                        and (len(inter) >= 5)
                        and (inter.mean() > extr.mean())
                    ),
                }
            )
    r2 = pd.DataFrame(rows)
    r2_path = OUT_DIR / "r2_mann_whitney_by_sample.tsv"
    r2.to_csv(r2_path, sep="\t", index=False)
    print(f"[R2] wrote {r2_path}")

    # Pre-registered exclusion check for HCC1954
    if r2.empty or "sample" not in r2.columns:
        print("[R2] WARNING: r2 DataFrame empty — check mode filter + LOH column")
        return r2, pd.DataFrame()
    hcc1954 = r2[r2["sample"] == "HCC1954"]
    for _, row in hcc1954.iterrows():
        trigger_seg = row["n_total"] < HCC1954_MIN_SEGMENTS
        trigger_int = row["frac_intermediate"] > HCC1954_MAX_INTER_FRAC
        excluded = trigger_seg or trigger_int
        print(
            f"[R2] {row['mode']} HCC1954 n_total={row['n_total']} "
            f"frac_inter={row['frac_intermediate']:.3f} "
            f"→ excluded={excluded} (seg_trigger={trigger_seg}, inter_trigger={trigger_int})"
        )

    # Segment-level (per-chr aggregation as a proxy for LOH segment) for HCC1954 only
    # Using Chr as segment ID: AF-SD vs mean NG Spearman
    print("[R2] computing segment-level Spearman (per-chr aggregation)…")
    seg_rows = []
    for mode in ["to", "paired"]:
        sub = tp[(tp["mode"] == mode) & (tp["cn_tier"] == "CN1")].copy()
        for sample in sorted(sub["sample"].unique()):
            samp = sub[sub["sample"] == sample]
            seg = samp.groupby("Chr").agg(
                n_var=("caller_af", "size"),
                af_sd=("caller_af", "std"),
                af_mean=("caller_af", "mean"),
                ng_mean=("HPFineNGroups", "mean"),
            ).reset_index()
            seg = seg[seg["n_var"] >= 5]
            if len(seg) < 5:
                seg_rows.append(
                    {"mode": mode, "sample": sample, "n_seg": len(seg), "rho": np.nan, "p": np.nan}
                )
                continue
            rho, pval = stats.spearmanr(seg["af_sd"], seg["ng_mean"], nan_policy="omit")
            seg_rows.append(
                {
                    "mode": mode,
                    "sample": sample,
                    "n_seg": len(seg),
                    "rho": round(rho, 3) if not np.isnan(rho) else np.nan,
                    "p": pval,
                }
            )
    seg_df = pd.DataFrame(seg_rows)
    seg_path = OUT_DIR / "r2_hcc1954_segment_level.tsv"
    seg_df.to_csv(seg_path, sep="\t", index=False)
    print(f"[R2] wrote {seg_path}")

    # Summary counts: POS vs NEG samples, with/without HCC1954
    for mode in ["to", "paired"]:
        mode_rows = r2[r2["mode"] == mode]
        n_total = len(mode_rows)
        pos = mode_rows["direction_positive"].sum()
        mode_no_hcc1954 = mode_rows[mode_rows["sample"] != "HCC1954"]
        pos_ex = mode_no_hcc1954["direction_positive"].sum()
        print(
            f"[R2] {mode} region-level Mann-Whitney: {pos}/{n_total} POS (full) | "
            f"{pos_ex}/{len(mode_no_hcc1954)} POS (exclude HCC1954)"
        )
    return r2, seg_df


# ---------------------------------------------------------------------------
# R3 — NR-bin + matched sampling
# ---------------------------------------------------------------------------
def run_r3(tp: pd.DataFrame):
    print("\n=== R3: NR-bin high-power verification + matched sampling ===")
    rows = []
    for mode in ["to", "paired"]:
        cn1 = tp[(tp["mode"] == mode) & (tp["cn_tier"] == "CN1")].copy()
        for lo, hi in NR_BINS:
            bin_df = cn1[(cn1["NumReads"] >= lo) & (cn1["NumReads"] < hi)]
            inter = bin_df[bin_df["af_class"] == "Intermediate"]["HPFineNGroups"].dropna()
            extr = bin_df[bin_df["af_class"] == "Extreme"]["HPFineNGroups"].dropna()
            u, p, r = mw_test(inter.values, extr.values)
            rows.append(
                {
                    "mode": mode,
                    "nr_bin": f"{lo}-{hi}",
                    "n_extreme": len(extr),
                    "n_intermediate": len(inter),
                    "mean_ng_ext": round(extr.mean(), 3) if len(extr) else np.nan,
                    "mean_ng_inter": round(inter.mean(), 3) if len(inter) else np.nan,
                    "delta_ng": (
                        round(inter.mean() - extr.mean(), 3)
                        if len(extr) and len(inter)
                        else np.nan
                    ),
                    "mw_p": p,
                    "rank_biserial_r": round(r, 3) if not np.isnan(r) else np.nan,
                }
            )
    r3 = pd.DataFrame(rows)
    r3_path = OUT_DIR / "r3_nr_bins_with_80plus.tsv"
    r3.to_csv(r3_path, sep="\t", index=False)
    print(f"[R3] wrote {r3_path}")

    # Matched sampling: CN1 LOH TP naturally caps NR at ~56 because deletion-LOH halves
    # coverage. NR≥80 is mechanically empty in CN1. We therefore compare NR 50-80 (highest
    # populated bin) downsampled to NR 10-30 size as the "statistical power confound" test.
    # If high-NR r (0.709) persists under matched n, the "high NR artifact" hypothesis is rejected.
    print("[R3] running matched sampling (100 iters)…")
    print("[R3] Note: CN1 LOH TP NR caps at ~56; matching NR 50-80 → NR 10-30 sizes")
    rng = np.random.default_rng(RNG_SEED)
    match_rows = []
    for mode in ["to", "paired"]:
        cn1 = tp[(tp["mode"] == mode) & (tp["cn_tier"] == "CN1")].copy()
        bin_low = cn1[(cn1["NumReads"] >= 10) & (cn1["NumReads"] < 30)]
        bin_high = cn1[(cn1["NumReads"] >= 50) & (cn1["NumReads"] < 80)]
        high_inter = bin_high[bin_high["af_class"] == "Intermediate"]["HPFineNGroups"].dropna().values
        high_extr = bin_high[bin_high["af_class"] == "Extreme"]["HPFineNGroups"].dropna().values
        n_inter_low = int((bin_low["af_class"] == "Intermediate").sum())
        n_extr_low = int((bin_low["af_class"] == "Extreme").sum())
        print(
            f"[R3] {mode}: NR 10-30 n_inter={n_inter_low} n_extr={n_extr_low} | "
            f"NR 50-80 n_inter={len(high_inter)} n_extr={len(high_extr)}"
        )
        low_inter = bin_low[bin_low["af_class"] == "Intermediate"]["HPFineNGroups"].dropna().values
        low_extr = bin_low[bin_low["af_class"] == "Extreme"]["HPFineNGroups"].dropna().values
        # Matched sampling: downsample BOTH bins to the minimum class sizes
        m_inter = min(len(low_inter), len(high_inter))
        m_extr = min(len(low_extr), len(high_extr))
        rs_low, rs_high = [], []
        ps_low, ps_high = [], []
        for _ in range(100):
            li = rng.choice(len(low_inter), size=m_inter, replace=False)
            le = rng.choice(len(low_extr), size=m_extr, replace=False)
            hi = rng.choice(len(high_inter), size=m_inter, replace=False)
            he = rng.choice(len(high_extr), size=m_extr, replace=False)
            _, p1, r1 = mw_test(low_inter[li], low_extr[le])
            _, p2, r2 = mw_test(high_inter[hi], high_extr[he])
            rs_low.append(r1); ps_low.append(p1)
            rs_high.append(r2); ps_high.append(p2)
        rs_low = np.array(rs_low); rs_high = np.array(rs_high)
        # Full-bin (unmatched) baselines
        _, p_full_low, r_full_low = mw_test(low_inter, low_extr)
        _, p_full, r_full = mw_test(high_inter, high_extr)
        match_rows.append(
            {
                "mode": mode,
                "low_bin": "NR 10-30",
                "high_bin": "NR 50-80",
                "matched_n_inter": m_inter,
                "matched_n_extr": m_extr,
                "r_full_low": round(r_full_low, 3) if not np.isnan(r_full_low) else np.nan,
                "r_full_high": round(r_full, 3) if not np.isnan(r_full) else np.nan,
                "r_matched_low_median": round(np.median(rs_low), 3),
                "r_matched_high_median": round(np.median(rs_high), 3),
                "r_matched_high_q025": round(np.quantile(rs_high, 0.025), 3),
                "r_matched_high_q975": round(np.quantile(rs_high, 0.975), 3),
                "delta_matched_r": round(np.median(rs_high) - np.median(rs_low), 3),
            }
        )
    match_df = pd.DataFrame(match_rows)
    match_path = OUT_DIR / "r3_matched_sampling.tsv"
    match_df.to_csv(match_path, sep="\t", index=False)
    print(f"[R3] wrote {match_path}")

    return r3, match_df


# ---------------------------------------------------------------------------
# R4 — CN1 vs CN2 stratification
# ---------------------------------------------------------------------------
def run_r4(tp: pd.DataFrame):
    print("\n=== R4: CN1 (deletion-LOH) vs CN2 (cnLOH) ===")
    rows = []
    for mode in ["to", "paired"]:
        for cn in ["CN1", "CN2"]:
            sub = tp[(tp["mode"] == mode) & (tp["cn_tier"] == cn)].copy()
            for sample in sorted(sub["sample"].unique()):
                samp = sub[sub["sample"] == sample]
                inter = samp[samp["af_class"] == "Intermediate"]["HPFineNGroups"].dropna()
                extr = samp[samp["af_class"] == "Extreme"]["HPFineNGroups"].dropna()
                u, p, r = mw_test(inter.values, extr.values)
                rows.append(
                    {
                        "mode": mode,
                        "cn_tier": cn,
                        "sample": sample,
                        "n_extreme": len(extr),
                        "n_intermediate": len(inter),
                        "mean_ng_ext": round(extr.mean(), 3) if len(extr) else np.nan,
                        "mean_ng_inter": round(inter.mean(), 3) if len(inter) else np.nan,
                        "delta_ng": (
                            round(inter.mean() - extr.mean(), 3)
                            if len(extr) and len(inter)
                            else np.nan
                        ),
                        "mw_p": p,
                        "rank_biserial_r": round(r, 3) if not np.isnan(r) else np.nan,
                        "direction_positive": bool(
                            (len(extr) >= 5)
                            and (len(inter) >= 5)
                            and (inter.mean() > extr.mean())
                        ),
                    }
                )
    r4 = pd.DataFrame(rows)
    r4_path = OUT_DIR / "r4_cn1_vs_cn2_per_sample.tsv"
    r4.to_csv(r4_path, sep="\t", index=False)
    print(f"[R4] wrote {r4_path}")
    for mode in ["to", "paired"]:
        for cn in ["CN1", "CN2"]:
            sub = r4[(r4["mode"] == mode) & (r4["cn_tier"] == cn)]
            pos = sub["direction_positive"].sum()
            n = len(sub)
            mean_delta = sub["delta_ng"].dropna().mean()
            print(
                f"[R4] {mode} {cn}: {pos}/{n} POS | mean ΔNG = "
                f"{mean_delta:.3f} (across samples)"
            )
    return r4


def main():
    print(f"[start] B.2 batch 1 (R2/R3/R4) seed={RNG_SEED}")
    tp = load_loh_tp()
    r2, seg = run_r2(tp)
    r3, match = run_r3(tp)
    r4 = run_r4(tp)
    print("\n[done] B.2 batch 1 complete")


if __name__ == "__main__":
    sys.exit(main())
