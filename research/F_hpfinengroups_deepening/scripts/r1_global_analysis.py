#!/usr/bin/env python3
"""R1-Global analysis — HCC1395 TO arm on full 40,239 variants.

Scope:
  - per-feature AUC on global regions (TP vs FP) — raw + bootstrap 95% CI
  - Residualized AUC on NumReads + Coverage_Multiple (no caller_af column in ISM output;
    use NumReads as NR proxy and Coverage_Multiple as CN/depth normalization)
  - Compare vs F pilot subset (CL-025a ⭐3) for overfit assessment
  - Write step8_r1_global_to_arm.md with reasoning chain

Outputs:
  research/F_hpfinengroups_deepening/observations/step8_r1_global_to_arm.md
  research/F_hpfinengroups_deepening/figures/step8_global_feature_auc.png
  research/F_hpfinengroups_deepening/figures/step8_global_vs_subset_overfit.png
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LinearRegression
from scipy.stats import mannwhitneyu

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
GLOBAL_DIR = ROOT / "output/hcc1395_normal_pilot_global/TO"
SUBSET_DIR = ROOT / "output/hcc1395_normal_pilot/TO"
DATA = ROOT / "research/F_hpfinengroups_deepening/data"
OBS = ROOT / "research/F_hpfinengroups_deepening/observations"
FIG = ROOT / "research/F_hpfinengroups_deepening/figures"
RNG = np.random.default_rng(2026)

FEATURES = [
    "SampleASM_Delta", "HP_Residual_Delta", "HP_Signed_Residual",
    "Combined_HP_Signed_Delta", "Tumor_HP_Delta", "Normal_HP_Delta",
    "NormalBaseline_Mean", "NormalBaseline_Coverage",
    "Fisher_Frac_Sig", "Entropy_Imbalance", "Epipoly_Delta",
    "PairwiseMedianDist", "HPFineF", "HPFineNGroups", "CramersV",
]


def load_data(sig_path, truth_tsv):
    """Load significance summary + merge truth labels."""
    sig = pd.read_csv(sig_path, low_memory=False)
    truth = pd.read_csv(truth_tsv, sep="\t")
    truth.columns = [c.strip() for c in truth.columns]
    sig["key"] = sig["Chr"].astype(str) + ":" + sig["Pos"].astype(str)
    truth["key"] = truth["Chr"].astype(str) + ":" + truth["Pos"].astype(str)
    merged = sig.merge(truth[["key", "truth_label"]], on="key", how="left")
    merged["is_tp"] = (merged["truth_label"].astype(str) == "TP").astype(int)
    return merged


def auc_ci(y_true, score, n_boot=500, seed=2026):
    """AUC + 95% bootstrap CI. Flip if inverse direction (AUC<0.5)."""
    score = np.asarray(score, dtype=float)
    y_true = np.asarray(y_true, dtype=int)
    mask = ~np.isnan(score) & ~np.isnan(y_true)
    y_true = y_true[mask]
    score = score[mask]
    if len(np.unique(y_true)) < 2 or len(y_true) < 20:
        return np.nan, np.nan, np.nan, len(y_true)
    auc = roc_auc_score(y_true, score)
    if auc < 0.5:
        score = -score
        auc = 1 - auc
    rng = np.random.default_rng(seed)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], score[idx]))
    if not aucs:
        return auc, np.nan, np.nan, len(y_true)
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return auc, lo, hi, len(y_true)


def residualize(df, feature, confounds):
    """OLS residualize feature on confounds (dropna on all columns). Returns residual series."""
    cols = [feature] + confounds
    sub = df.dropna(subset=cols).copy()
    if len(sub) < 50:
        return None, sub.index
    X = sub[confounds].values
    y = sub[feature].values
    lr = LinearRegression().fit(X, y)
    resid = y - lr.predict(X)
    return pd.Series(resid, index=sub.index), sub.index


def compute_all(df, label="global", confounds=("NumReads", "Coverage_Multiple")):
    """Per-feature raw + residualized AUC with CI."""
    confounds = list(confounds)
    rows = []
    for feat in FEATURES:
        if feat not in df.columns:
            continue
        auc_raw, lo_raw, hi_raw, n_raw = auc_ci(df["is_tp"].values, df[feat].values)
        resid, idx = residualize(df, feat, confounds)
        if resid is None:
            auc_res, lo_res, hi_res, n_res = np.nan, np.nan, np.nan, 0
        else:
            auc_res, lo_res, hi_res, n_res = auc_ci(df.loc[idx, "is_tp"].values, resid.values)
        rows.append({
            "feature": feat,
            "auc_raw": auc_raw, "ci_raw_lo": lo_raw, "ci_raw_hi": hi_raw, "n_raw": n_raw,
            "auc_res": auc_res, "ci_res_lo": lo_res, "ci_res_hi": hi_res, "n_res": n_res,
        })
    return pd.DataFrame(rows).sort_values("auc_raw", ascending=False, na_position="last")


def compare_pilot_vs_global(global_df, subset_df):
    """Compare per-feature AUC between global and F pilot subset (overfit assessment)."""
    rows = []
    for feat in FEATURES:
        if feat not in global_df.columns or feat not in subset_df.columns:
            continue
        auc_g, _, _, n_g = auc_ci(global_df["is_tp"].values, global_df[feat].values, n_boot=200)
        auc_s, _, _, n_s = auc_ci(subset_df["is_tp"].values, subset_df[feat].values, n_boot=200)
        # Check F pilot canonical filter subset in global (NR>=80 + NG>=4 + NonLOH)
        # Note: F pilot also had AF<0.4 filter, but AF not in ISM output; NG>=4 proxies it
        nonloh = (global_df.get("Potential_LOH", "False").astype(str).str.lower() == "false")
        filt = global_df[(global_df["NumReads"] >= 80) & (global_df["HPFineNGroups"] >= 4) & nonloh]
        if len(filt) >= 50 and filt["is_tp"].nunique() >= 2:
            auc_gf, _, _, n_gf = auc_ci(filt["is_tp"].values, filt[feat].values, n_boot=200)
        else:
            auc_gf, n_gf = np.nan, 0
        rows.append({
            "feature": feat,
            "auc_global": auc_g, "n_global": n_g,
            "auc_subset_old": auc_s, "n_subset_old": n_s,
            "auc_global_fpilot_filter": auc_gf, "n_global_fpilot": n_gf,
        })
    return pd.DataFrame(rows)


def plot_feature_auc(rdf, out_path, title):
    """Horizontal bar plot of raw + residualized AUC with CI error bars."""
    fig, ax = plt.subplots(figsize=(9, 6))
    y = np.arange(len(rdf))
    ax.barh(y - 0.2, rdf["auc_raw"], height=0.38, label="Raw", color="#3465a4", alpha=0.85)
    ax.errorbar(rdf["auc_raw"], y - 0.2,
                xerr=[rdf["auc_raw"] - rdf["ci_raw_lo"], rdf["ci_raw_hi"] - rdf["auc_raw"]],
                fmt='none', ecolor='#204a87', capsize=2)
    ax.barh(y + 0.2, rdf["auc_res"], height=0.38, label="Residualized", color="#cc0000", alpha=0.85)
    ax.errorbar(rdf["auc_res"], y + 0.2,
                xerr=[rdf["auc_res"] - rdf["ci_res_lo"], rdf["ci_res_hi"] - rdf["auc_res"]],
                fmt='none', ecolor='#a40000', capsize=2)
    ax.set_yticks(y)
    ax.set_yticklabels(rdf["feature"])
    ax.axvline(0.5, color="gray", ls="--", lw=1, alpha=0.5)
    ax.axvline(0.58, color="orange", ls=":", lw=1, alpha=0.7, label="CL-008 ceiling 0.58")
    ax.axvline(0.65, color="green", ls=":", lw=1, alpha=0.7, label="Promising 0.65")
    ax.set_xlabel("AUC (TP vs FP)")
    ax.set_xlim(0.4, 0.85)
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=8)
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[+] Saved {out_path}")


def plot_global_vs_subset(cmp_df, out_path):
    fig, ax = plt.subplots(figsize=(9, 6))
    y = np.arange(len(cmp_df))
    w = 0.27
    ax.barh(y - w, cmp_df["auc_subset_old"], height=w, label="F pilot subset (n≈983)", color="#6e8d8d")
    ax.barh(y, cmp_df["auc_global_fpilot_filter"], height=w, label="Global ∩ F pilot filter", color="#3465a4")
    ax.barh(y + w, cmp_df["auc_global"], height=w, label="Global full (40,239)", color="#cc0000")
    ax.set_yticks(y)
    ax.set_yticklabels(cmp_df["feature"])
    ax.axvline(0.5, color="gray", ls="--", lw=1, alpha=0.5)
    ax.axvline(0.58, color="orange", ls=":", lw=1, alpha=0.7)
    ax.set_xlabel("AUC (TP vs FP)")
    ax.set_xlim(0.4, 0.85)
    ax.set_title("Global vs F pilot subset — overfit assessment")
    ax.legend(loc="lower right", fontsize=8)
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[+] Saved {out_path}")


def write_report(global_df, global_rdf, subset_rdf, cmp_df, fig_auc, fig_cmp, out_md):
    tp = int(global_df["is_tp"].sum())
    fp = len(global_df) - tp

    # Identify promising features (residualized AUC ≥ 0.60 with CI lower ≥ 0.55)
    promising = global_rdf[
        (global_rdf["auc_res"] >= 0.60) & (global_rdf["ci_res_lo"] >= 0.55)
    ].copy()

    lines = [
        "# Step 8 — R1-Global HCC1395 Normal BAM Phase 2 TO arm",
        "",
        "**Date**: 2026-04-21",
        "**Scope**: Full 40,239 TO variants (28,396 TP + 11,843 FP) × inter_sub_mod Phase 2 pipeline × full 144GB Normal BAM random access",
        "**Purpose**: 驗證 F pilot CL-025a（SampleASM_Delta ⭐3 characterization-only）在全域 region 是否保持；評估 paired arm 放棄後 TO arm 是否有獨立的 F1-filter 潛力。",
        "",
        "## 1. 規模",
        "",
        f"- Global regions processed: {len(global_df)} (vs F pilot subset 801)",
        f"- TP={tp}, FP={fp}, TP_rate={tp/(tp+fp):.4f} (vs F pilot subset 0.9700)",
        f"- FP={fp} is {fp/12:.0f}x F pilot subset's 12 FP — bootstrap CI now statistically meaningful",
        "",
        "## 2. 全域 per-feature AUC（raw + residualized on NumReads+Coverage_Multiple）",
        "",
        f"![]({fig_auc.relative_to(fig_auc.parent.parent)})",
        "",
        "| Feature | Raw AUC [95% CI] | n | Residualized AUC [95% CI] | n |",
        "|---------|------------------|---|---------------------------|---|",
    ]
    for _, r in global_rdf.iterrows():
        raw = f"{r['auc_raw']:.3f} [{r['ci_raw_lo']:.3f},{r['ci_raw_hi']:.3f}]" if pd.notna(r['auc_raw']) else "—"
        res = f"{r['auc_res']:.3f} [{r['ci_res_lo']:.3f},{r['ci_res_hi']:.3f}]" if pd.notna(r['auc_res']) else "—"
        lines.append(f"| {r['feature']} | {raw} | {int(r['n_raw'])} | {res} | {int(r['n_res']) if pd.notna(r['n_res']) else 0} |")

    lines += [
        "",
        "### 2.1 Promising features（residualized AUC ≥0.60 且 CI 下界 ≥0.55）",
        "",
    ]
    if len(promising) > 0:
        lines.append("| Feature | Residualized AUC [95% CI] |")
        lines.append("|---------|---------------------------|")
        for _, r in promising.iterrows():
            res = f"{r['auc_res']:.3f} [{r['ci_res_lo']:.3f},{r['ci_res_hi']:.3f}]"
            lines.append(f"| {r['feature']} | {res} |")
    else:
        lines.append("**無特徵通過標準** — 全域 Phase 2 TO arm 無 F1-filter 候選。")

    lines += [
        "",
        "## 3. 與 F pilot subset 對照（overfit 評估）",
        "",
        f"![]({fig_cmp.relative_to(fig_cmp.parent.parent)})",
        "",
        "| Feature | F pilot subset old (n=801) | Global ∩ F pilot filter | Global full (40,239) |",
        "|---------|----------------------------|-------------------------|----------------------|",
    ]
    for _, r in cmp_df.iterrows():
        s = f"{r['auc_subset_old']:.3f} (n={int(r['n_subset_old'])})" if pd.notna(r['auc_subset_old']) else "—"
        gf = f"{r['auc_global_fpilot_filter']:.3f} (n={int(r['n_global_fpilot'])})" if pd.notna(r['auc_global_fpilot_filter']) else "—"
        g = f"{r['auc_global']:.3f} (n={int(r['n_global'])})" if pd.notna(r['auc_global']) else "—"
        lines.append(f"| {r['feature']} | {s} | {gf} | {g} |")

    lines += [
        "",
        "**Overfit heuristic**：若 `auc_subset_old` >> `auc_global_fpilot_filter` >> `auc_global` 遞減 → F pilot 子集內的 AUC 為 pre-selection overfit；若 `auc_global_fpilot_filter` ≈ `auc_subset_old` → filter 真有 enrichment 效果；若 `auc_global` ≈ `auc_subset_old` → 特徵本身 robust。",
        "",
        "## 4. 推論與結論",
        "",
        "### 4.1 reasoning chain",
        "",
        "1. **Premise**：F pilot subset Step 7 顯示 paired `Fisher_Frac_Sig` 在 941-region 子集內 AUC=0.726，但殘差化後 CI 跨隨機 + TP 飽和（99.5%），paired F1-filter 方向 2026-04-21 放棄。",
        "2. **Question**：TO arm `SampleASM_Delta` 在同 subset 殘差化 0.610（CL-025a ⭐3）是否為真實訊號？",
        "3. **Test**：全域 40,239 region 不經 F pilot filter 直接評估；若 residualized AUC ≥0.60 且 CI 下界 ≥0.55 → 真訊號；若顯著衰減 → subset artifact。",
        "4. **Result**：見 §2.1 與 §3 表。",
        "5. **Interpretation**：",
    ]
    if len(promising) > 0:
        lines.append("   - **有通過標準的特徵** → TO arm 有 F1-filter 候選，下一步 R12 跨樣本驗證。")
    else:
        lines.append("   - **無特徵通過標準** → TO arm Phase 2 SampleASM_Delta 在全域 region 無 F1-filter 價值；CL-025a 降級為 characterization-only（或 concluded NEGATIVE 若連 characterization 都無）。")

    lines += [
        "",
        "### 4.2 對 Registry 的影響",
        "",
        "- **CL-025a**：原 ⭐3 (characterization-only, under confound scrutiny)",
    ]
    sampleasm_row = global_rdf[global_rdf["feature"] == "SampleASM_Delta"]
    if not sampleasm_row.empty:
        r = sampleasm_row.iloc[0]
        if pd.notna(r["auc_res"]) and r["auc_res"] >= 0.60 and r["ci_res_lo"] >= 0.55:
            lines.append(f"  → **升級為 ⭐4**（global residualized AUC={r['auc_res']:.3f} CI [{r['ci_res_lo']:.3f},{r['ci_res_hi']:.3f}]，全域穩定）")
        elif pd.notna(r["auc_res"]) and r["auc_res"] >= 0.55:
            lines.append(f"  → **維持 ⭐3**（global residualized AUC={r['auc_res']:.3f} 接近但未達 0.60 閾值，保留 characterization 角色）")
        else:
            lines.append(f"  → **降級為 ⭐2 concluded NEGATIVE**（global residualized AUC={r['auc_res']:.3f} 不足支撐）")
    lines += [
        "- **CL-008**（Beyond-AUC ≤0.58 ceiling）：視是否有特徵突破決定",
        "",
        "## 5. 下一步候選",
        "",
        "1. 若 §2.1 有通過特徵 → R12 跨 6 樣本擴展驗證",
        "2. 若無通過特徵 → 批次 3 轉向其他方向（R8 Per-CpG ASM residualize、R12 non-HCC1395 Phase 2A）",
        "",
        "## 6. Provenance",
        "",
        "- Binary: `build/bin/inter_sub_mod`",
        f"- Global VCF: `output/hcc1395_normal_pilot_global/HCC1395_TO_global.vcf.gz` (40,239 variants)",
        "- Tumor BAM: `tumor_tagged.bam` (TO arm, 2026-03-07 step03_longphase_to)",
        "- Normal BAM: `HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam` (full 144GB, BAM index random access)",
        "- LOH bed: `tumor_phased_LOH.bed` (1,100 regions)",
        "- CLI: `-w 5000 -j 16 --distance-metric BERNOULLI`",
        "- Pipeline script: `r1_run_phase2_global_to.sh`",
        "- Analysis script: `r1_global_analysis.py`",
        "",
        "## 7. Registry links",
        "",
        "- CL-025a (SampleASM_Delta characterization)",
        "- CL-025c (Cross-mode concordance 2026-04-21)",
        "- CL-008 (Beyond-AUC ceiling)",
    ]
    out_md.write_text("\n".join(lines))
    print(f"[+] Wrote {out_md}")


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    OBS.mkdir(parents=True, exist_ok=True)

    global_sig = GLOBAL_DIR / "significance_summary.csv"
    subset_sig = SUBSET_DIR / "significance_summary.csv"
    truth_global = DATA / "r1_hcc1395_filter_snvs_to_global.tsv"
    truth_subset = DATA / "r1_hcc1395_filter_snvs_to.tsv"

    for p in [global_sig, subset_sig, truth_global, truth_subset]:
        if not p.exists():
            print(f"[!] Missing: {p}")
            return

    print(f"Loading global: {global_sig}")
    global_df = load_data(global_sig, truth_global)
    print(f"  n={len(global_df)}, TP={int(global_df['is_tp'].sum())}, FP={len(global_df) - int(global_df['is_tp'].sum())}")

    print(f"Loading F pilot subset: {subset_sig}")
    subset_df = load_data(subset_sig, truth_subset)
    print(f"  n={len(subset_df)}, TP={int(subset_df['is_tp'].sum())}, FP={len(subset_df) - int(subset_df['is_tp'].sum())}")

    print("\nComputing global feature AUC (raw + residualized)...")
    global_rdf = compute_all(global_df, "global")
    print(global_rdf[["feature", "auc_raw", "auc_res", "ci_res_lo", "ci_res_hi", "n_raw"]].to_string(index=False))

    print("\nComputing F pilot subset feature AUC...")
    subset_rdf = compute_all(subset_df, "subset")

    print("\nComputing global vs subset comparison...")
    cmp_df = compare_pilot_vs_global(global_df, subset_df)
    print(cmp_df.to_string(index=False))

    fig_auc = FIG / "step8_global_feature_auc.png"
    fig_cmp = FIG / "step8_global_vs_subset_overfit.png"
    plot_feature_auc(global_rdf, fig_auc, "R1-Global HCC1395 TO arm — Phase 2 per-feature AUC (N=40,239)")
    plot_global_vs_subset(cmp_df, fig_cmp)

    global_rdf.to_csv(ROOT / "research/F_hpfinengroups_deepening/data/step8_global_feature_auc.tsv", sep="\t", index=False)
    cmp_df.to_csv(ROOT / "research/F_hpfinengroups_deepening/data/step8_global_vs_subset.tsv", sep="\t", index=False)

    out_md = OBS / "step8_r1_global_to_arm.md"
    write_report(global_df, global_rdf, subset_rdf, cmp_df, fig_auc, fig_cmp, out_md)

    print("\n[+] R1-Global analysis complete.")


if __name__ == "__main__":
    main()
