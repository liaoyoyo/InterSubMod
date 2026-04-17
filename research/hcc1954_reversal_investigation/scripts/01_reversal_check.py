#!/usr/bin/env python3
"""
B.2-1 HCC1954 反向排除標準調查
================================
Question: HCC1954 step3 ρ=-0.30 是真反向還是小樣本雜訊？
          Step2 Δ=+0.69 POS 與 Step3 ρ=-0.30 NEG 如何調和？
          移除 HCC1954 後「7/7 一致」是否變成「6/7」？

輸入：
  - research/loh_subclone_af/data/step3_segment_statistics.tsv (TO segment-level)
  - research/loh_subclone_af_paired/data/step3_paired_segment_statistics.tsv (Paired segment-level)
  - all_region_rows.tsv.gz (for region-level reconciliation)

輸出：
  - data/rho_bootstrap_ci.tsv
  - data/sensitivity_leave_one_out.tsv
  - data/hcc1954_region_level_rho.tsv
  - figures/01_bootstrap_ci_forest.png
  - figures/02_hcc1954_vs_others_segment_dist.png
"""

import pandas as pd
import numpy as np
from scipy import stats as scipy_stats
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# ── Paths ──
ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
MASTER = "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"
TO_SEG = ROOT / "research/loh_subclone_af/data/step3_segment_statistics.tsv"
PAIRED_SEG = ROOT / "research/loh_subclone_af_paired/data/step3_paired_segment_statistics.tsv"
OUTDIR = ROOT / "research/hcc1954_reversal_investigation"
DATA = OUTDIR / "data"
FIG = OUTDIR / "figures"
DATA.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]

RNG = np.random.default_rng(42)


def bootstrap_rho_ci(x, y, n_boot=2000, seed=42):
    """Bootstrap 95% CI for Spearman rho."""
    rng = np.random.default_rng(seed)
    n = len(x)
    if n < 3:
        return np.nan, np.nan, np.nan
    rhos = []
    idx = np.arange(n)
    for _ in range(n_boot):
        s = rng.choice(idx, size=n, replace=True)
        xs, ys = x[s], y[s]
        if np.std(xs) == 0 or np.std(ys) == 0:
            continue
        r, _ = scipy_stats.spearmanr(xs, ys)
        if not np.isnan(r):
            rhos.append(r)
    if not rhos:
        return np.nan, np.nan, np.nan
    rhos = np.array(rhos)
    r_obs, _ = scipy_stats.spearmanr(x, y)
    lo, hi = np.percentile(rhos, [2.5, 97.5])
    return r_obs, lo, hi


def analyze_segment_level(seg_df, mode_name):
    """Per-sample bootstrap CI — 重現原 step3 設定：CN1 only, ρ(af_sd, mean_ngroups)."""
    rows = []
    for sample in SAMPLE_ORDER:
        d = seg_df[(seg_df["sample"] == sample) & (seg_df["dominant_cn"] == "CN1")]
        n = len(d)
        if n < 3:
            rows.append(dict(mode=mode_name, sample=sample, n_cn1_segments=n,
                             rho=np.nan, ci_low=np.nan, ci_high=np.nan,
                             p_value=np.nan, ci_contains_zero=np.nan))
            continue
        r, p = scipy_stats.spearmanr(d["af_sd"], d["mean_ngroups"])
        r_obs, lo, hi = bootstrap_rho_ci(d["af_sd"].values, d["mean_ngroups"].values, n_boot=2000)
        rows.append(dict(mode=mode_name, sample=sample, n_cn1_segments=n,
                         rho=r_obs, ci_low=lo, ci_high=hi,
                         p_value=p, ci_contains_zero=(lo <= 0 <= hi)))
    return pd.DataFrame(rows)


def leave_one_out_meta(df_rho):
    """Sensitivity: remove each sample, compute pooled meta-rho (Fisher z)."""
    def fisher_z(r):
        r = np.clip(r, -0.9999, 0.9999)
        return 0.5 * np.log((1 + r) / (1 - r))

    rows = []
    for mode in df_rho["mode"].unique():
        sub = df_rho[df_rho["mode"] == mode].dropna(subset=["rho"])
        # full
        zs = [fisher_z(r) for r in sub["rho"]]
        ws = [n - 3 for n in sub["n_cn1_segments"]]  # inverse variance ~ n-3
        mean_z_full = np.sum(np.array(zs) * np.array(ws)) / np.sum(ws)
        mean_r_full = (np.exp(2 * mean_z_full) - 1) / (np.exp(2 * mean_z_full) + 1)
        rows.append(dict(mode=mode, excluded="NONE",
                         n_samples=len(sub), meta_rho=mean_r_full,
                         min_sample_rho=sub["rho"].min(), max_sample_rho=sub["rho"].max()))
        for i, s in enumerate(sub["sample"]):
            zs_loo = [fisher_z(r) for j, r in enumerate(sub["rho"]) if sub["sample"].iloc[j] != s]
            ws_loo = [n - 3 for j, n in enumerate(sub["n_cn1_segments"]) if sub["sample"].iloc[j] != s]
            mean_z = np.sum(np.array(zs_loo) * np.array(ws_loo)) / np.sum(ws_loo)
            mean_r = (np.exp(2 * mean_z) - 1) / (np.exp(2 * mean_z) + 1)
            sub_loo = sub[sub["sample"] != s]
            rows.append(dict(mode=mode, excluded=s,
                             n_samples=len(sub_loo), meta_rho=mean_r,
                             min_sample_rho=sub_loo["rho"].min(),
                             max_sample_rho=sub_loo["rho"].max()))
    return pd.DataFrame(rows)


def analyze_af_mean_all_cn(seg_df, mode_name):
    """對照分析：所有 CN、ρ(af_mean, mean_ngroups)。檢查是否有不同方向。"""
    rows = []
    for sample in SAMPLE_ORDER:
        d = seg_df[seg_df["sample"] == sample]
        n = len(d)
        if n < 3:
            continue
        r, p = scipy_stats.spearmanr(d["af_mean"], d["mean_ngroups"])
        rows.append(dict(mode=mode_name, sample=sample, n_all_cn=n,
                         rho_af_mean=r, p_af_mean=p))
    return pd.DataFrame(rows)


def region_level_rho(sample_names=("HCC1954",)):
    """Reconcile Step2 vs Step3: compute per-sample region-level ρ(caller_af, HPFineNGroups)
    in LOH regions, then compare with segment-level ρ."""
    print("Loading master dataset for region-level reconciliation...")
    cols = ["sample", "mode", "truth_label", "caller_af", "HPFineNGroups", "to_loh_bed_hit"]
    df = pd.read_csv(MASTER, sep="\t", usecols=cols, low_memory=False)

    rows = []
    for sample in SAMPLE_ORDER:
        for mode in ["to", "paired"]:
            sub = df[(df["sample"] == sample) &
                     (df["mode"] == mode) &
                     (df["truth_label"] == "TP") &
                     (df["to_loh_bed_hit"].astype(str).str.lower() == "true")]
            n = len(sub)
            if n < 3 or sub["caller_af"].std() == 0:
                rows.append(dict(mode=mode, sample=sample, n_regions=n,
                                 rho_region=np.nan, p_region=np.nan))
                continue
            r, p = scipy_stats.spearmanr(sub["caller_af"], sub["HPFineNGroups"])
            rows.append(dict(mode=mode, sample=sample, n_regions=n,
                             rho_region=r, p_region=p))
    return pd.DataFrame(rows)


def plot_forest(df_rho, out_path):
    """Forest plot of ρ per sample with 95% CI, separated by mode."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, mode in zip(axes, ["TO", "Paired"]):
        sub = df_rho[df_rho["mode"] == mode].dropna(subset=["rho"]).copy()
        sub = sub.set_index("sample").reindex(SAMPLE_ORDER).reset_index()
        ypos = np.arange(len(sub))
        colors = ["red" if s == "HCC1954" else "steelblue" for s in sub["sample"]]
        ax.axvline(0, color="black", ls=":", lw=1)
        for i, row in sub.iterrows():
            if pd.isna(row["rho"]):
                continue
            ax.errorbar(row["rho"], i, xerr=[[row["rho"] - row["ci_low"]],
                                             [row["ci_high"] - row["rho"]]],
                        fmt="o", color=colors[i], ecolor=colors[i], capsize=4, ms=8)
            ax.text(row["rho"], i + 0.25,
                    f"ρ={row['rho']:.2f} (n_CN1={row['n_cn1_segments']})",
                    ha="center", fontsize=8)
        ax.set_yticks(ypos)
        ax.set_yticklabels(sub["sample"])
        ax.set_xlabel("Spearman ρ (segment-level af_mean ↔ mean_ngroups)")
        ax.set_title(f"{mode} mode — HCC1954 red")
        ax.set_xlim(-1, 1)
        ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def plot_segment_size_dist(to_seg_df, out_path):
    """HCC1954 vs others: segment count and variant-per-segment distribution."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    counts = to_seg_df.groupby("sample").size().reindex(SAMPLE_ORDER)
    colors = ["red" if s == "HCC1954" else "steelblue" for s in counts.index]
    axes[0].bar(counts.index, counts.values, color=colors)
    axes[0].set_ylabel("n segments (TP ≥2 variants)")
    axes[0].set_title("Segment count per sample (TO)")
    axes[0].tick_params(axis="x", rotation=45)

    for s in SAMPLE_ORDER:
        d = to_seg_df[to_seg_df["sample"] == s]["n_variants"]
        c = "red" if s == "HCC1954" else "gray"
        axes[1].hist(np.log10(d.clip(lower=1)), bins=25, alpha=0.5, label=s, color=c)
    axes[1].set_xlabel("log10(n_variants per segment)")
    axes[1].set_ylabel("count")
    axes[1].set_title("Variants-per-segment distribution")
    axes[1].legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def main():
    print("=== B.2-1 HCC1954 Reversal Investigation ===\n")

    print("[S1] Load segment-level stats")
    to_seg = pd.read_csv(TO_SEG, sep="\t")
    paired_seg = pd.read_csv(PAIRED_SEG, sep="\t")
    print(f"  TO segments: {len(to_seg)}")
    print(f"  Paired segments: {len(paired_seg)}")

    print("\n[S2] Bootstrap 95% CI for per-sample ρ (CN1 only, af_sd vs mean_ngroups, B=2000)")
    rho_to = analyze_segment_level(to_seg, "TO")
    rho_paired = analyze_segment_level(paired_seg, "Paired")
    df_rho = pd.concat([rho_to, rho_paired], ignore_index=True)
    df_rho.to_csv(DATA / "rho_bootstrap_ci.tsv", sep="\t", index=False)
    print(df_rho.to_string(index=False))

    print("\n[S2b] 對照分析：所有 CN + af_mean vs mean_ngroups（檢查不同角度）")
    df_alt_to = analyze_af_mean_all_cn(to_seg, "TO")
    df_alt_paired = analyze_af_mean_all_cn(paired_seg, "Paired")
    df_alt = pd.concat([df_alt_to, df_alt_paired], ignore_index=True)
    df_alt.to_csv(DATA / "rho_af_mean_all_cn.tsv", sep="\t", index=False)
    print(df_alt.to_string(index=False))

    print("\n[S3] HCC1954 specific:")
    hcc = df_rho[df_rho["sample"] == "HCC1954"]
    for _, r in hcc.iterrows():
        verdict = "NOISE (CI含0)" if r["ci_contains_zero"] else "TRUE REVERSAL (CI不含0)"
        print(f"  {r['mode']}: ρ={r['rho']:.3f} 95%CI=[{r['ci_low']:.3f}, {r['ci_high']:.3f}] n={r['n_cn1_segments']} → {verdict}")

    print("\n[S4] Leave-one-out sensitivity (meta-rho via Fisher z, weighted n-3)")
    df_loo = leave_one_out_meta(df_rho)
    df_loo.to_csv(DATA / "sensitivity_leave_one_out.tsv", sep="\t", index=False)
    print(df_loo.to_string(index=False))

    print("\n[S5] Region-level ρ for HCC1954 (Step2 vs Step3 reconciliation)")
    df_reg = region_level_rho()
    df_reg.to_csv(DATA / "hcc1954_region_level_rho.tsv", sep="\t", index=False)
    print(df_reg.to_string(index=False))

    print("\n[S6] Plot forest + segment distribution")
    plot_forest(df_rho, FIG / "01_bootstrap_ci_forest.png")
    plot_segment_size_dist(to_seg, FIG / "02_hcc1954_vs_others_segment_dist.png")

    # Final verdict
    print("\n=== FINAL VERDICT (B.2-1) ===")
    hcc_to = df_rho[(df_rho["sample"] == "HCC1954") & (df_rho["mode"] == "TO")].iloc[0]
    hcc_paired = df_rho[(df_rho["sample"] == "HCC1954") & (df_rho["mode"] == "Paired")].iloc[0]

    if hcc_to["ci_contains_zero"] and hcc_paired["ci_contains_zero"]:
        print("HCC1954 95% CI 含 0 於 TO + Paired → **反向為統計噪音**（小樣本）")
        print("→ 結論 11 實為「7/7 方向一致、1/7 HCC1954 功效不足（n<35）」")
    else:
        print(f"HCC1954 TO CI=[{hcc_to['ci_low']:.3f}, {hcc_to['ci_high']:.3f}]")
        print(f"HCC1954 Paired CI=[{hcc_paired['ci_low']:.3f}, {hcc_paired['ci_high']:.3f}]")
        print("→ 需 pre-registered exclusion criterion（基於 biology 而非 rho 本身）")

    # LOO check
    to_loo = df_loo[(df_loo["mode"] == "TO")]
    full_rho = to_loo[to_loo["excluded"] == "NONE"]["meta_rho"].iloc[0]
    hcc_excluded = to_loo[to_loo["excluded"] == "HCC1954"]["meta_rho"].iloc[0]
    print(f"\nMeta-ρ (TO): full 7 = {full_rho:.3f}, minus HCC1954 = {hcc_excluded:.3f}")
    delta = hcc_excluded - full_rho
    print(f"→ 排除 HCC1954 meta-ρ 變化 Δ={delta:+.3f}")
    if abs(delta) < 0.05:
        print("→ 影響極小，HCC1954 包含與否不改變結論方向")
    elif delta > 0:
        print("→ 排除 HCC1954 後 meta-ρ 增加 → HCC1954 確實壓低整體訊號")
    else:
        print("→ 反向：排除 HCC1954 後 meta-ρ 下降（HCC1954 並非拖累）")


if __name__ == "__main__":
    main()
