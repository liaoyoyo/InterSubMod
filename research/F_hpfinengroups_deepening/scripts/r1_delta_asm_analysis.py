#!/usr/bin/env python3
"""R1 Δ_ASM analysis — TO + paired arm independent packaging.

TO arm  : Primary goal = FP removal / F1 improvement
Paired arm: Primary goal = subclone characterization (TO 無法分辨的 region)

Inputs:
  output/hcc1395_normal_pilot/TO/significance_summary.csv   (801 regions)
  output/hcc1395_normal_pilot/paired/significance_summary.csv (983 regions)
  data/r1_hcc1395_filter_snvs_to.tsv   (TP/FP truth per region)
  data/r1_hcc1395_filter_snvs_paired.tsv

Outputs:
  research/F_hpfinengroups_deepening/observations/step7_hcc1395_normal_TO_pilot.md
  research/F_hpfinengroups_deepening/observations/step7_hcc1395_normal_paired_pilot.md
"""
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod")
PILOT = ROOT / "output/hcc1395_normal_pilot"
DATA = ROOT / "research/F_hpfinengroups_deepening/data"
OBS = ROOT / "research/F_hpfinengroups_deepening/observations"


def load_arm(mode):
    sig = pd.read_csv(PILOT / mode / "significance_summary.csv", low_memory=False)
    truth = pd.read_csv(DATA / f"r1_hcc1395_filter_snvs_{'to' if mode == 'TO' else mode}.tsv", sep="\t")
    truth.columns = [c.strip() for c in truth.columns]

    # Merge on Chr+Pos
    sig["key"] = sig["Chr"].astype(str) + ":" + sig["Pos"].astype(str)
    truth["key"] = truth["Chr"].astype(str) + ":" + truth["Pos"].astype(str)
    merged = sig.merge(truth[["key", "truth_label"]], on="key", how="left")
    return merged


def feature_auc(df, feature, label_col="is_tp"):
    sub = df[df[feature].notna() & df[label_col].notna()].copy()
    if sub[label_col].nunique() < 2 or len(sub) < 20:
        return None, len(sub)
    try:
        val = sub[feature].astype(float)
        auc = roc_auc_score(sub[label_col], val)
        # Flip if inversed
        return max(auc, 1 - auc), len(sub)
    except Exception:
        return None, len(sub)


def analyze_arm(mode):
    df = load_arm(mode)
    print(f"\n=== {mode} arm ===")
    print(f"Total regions: {len(df)}")
    print(f"Regions with truth: {df['truth_label'].notna().sum()}")

    df["is_tp"] = (df["truth_label"].astype(str) == "TP").astype(int)
    tp = int(df["is_tp"].sum())
    fp = len(df) - tp
    print(f"TP={tp}, FP={fp}, TP_rate={tp/(tp+fp):.4f}")

    features = [
        "SampleASM_Delta", "HP_Residual_Delta", "HP_Signed_Residual",
        "Combined_HP_Signed_Delta", "Tumor_HP_Delta", "Normal_HP_Delta",
        "NormalBaseline_Mean", "NormalBaseline_Coverage",
        "Fisher_Frac_Sig", "Entropy_Imbalance", "Epipoly_Delta",
        "PairwiseMedianDist", "HPFineF", "HPFineNGroups", "CramersV",
    ]
    results = []
    for feat in features:
        if feat not in df.columns:
            continue
        auc, n = feature_auc(df, feat)
        results.append({"feature": feat, "auc": auc, "n": n})
    rdf = pd.DataFrame(results).sort_values("auc", ascending=False, na_position="last")
    print("\n[Full subset AUC]")
    print(rdf.to_string(index=False))

    # Apply F pilot canonical filter intersect: NR≥80 + AF<0.4 (NumReads proxy, caller_af 不在 ISM output)
    # 這裡用 HPFineNGroups≥4 作為已 F pilot 語義上一致的 proxy (NG≥4 subset)
    nr_col = "NumReads"
    filt = df[(df[nr_col] >= 80) & (df["HPFineNGroups"] >= 4)].copy()
    print(f"\n[F pilot filter subset (NR≥80, NG≥4)] n={len(filt)} (TP={int(filt['is_tp'].sum())}, FP={len(filt)-int(filt['is_tp'].sum())})")
    if len(filt) >= 20 and filt["is_tp"].nunique() >= 2:
        fresults = []
        for feat in features:
            if feat not in filt.columns:
                continue
            auc, n = feature_auc(filt, feat)
            fresults.append({"feature": feat, "auc": auc, "n": n})
        frdf = pd.DataFrame(fresults).sort_values("auc", ascending=False, na_position="last")
        print(frdf.to_string(index=False))
    else:
        frdf = pd.DataFrame()
        print("  [!] Insufficient data for F pilot subset analysis")

    # Sample ASM positive vs negative count
    if "SampleASM_Sig" in df.columns:
        sig_pos = df[df["SampleASM_Sig"].astype(str).str.lower() == "true"]
        print(f"\n[SampleASM_Sig=true]: {len(sig_pos)} regions; TP rate={sig_pos['is_tp'].mean():.4f} (vs all {df['is_tp'].mean():.4f})")
    if "HP_Residual_Sig" in df.columns:
        res_pos = df[df["HP_Residual_Sig"].astype(str).str.lower() == "true"]
        print(f"[HP_Residual_Sig=true]: {len(res_pos)} regions; TP rate={res_pos['is_tp'].mean():.4f}")

    return df, rdf, frdf


def write_report(mode, df, rdf, frdf, out_md):
    tp = int(df["is_tp"].sum())
    fp = len(df) - tp
    lines = [
        f"# Step 7 — HCC1395 Normal BAM pilot, {mode} arm",
        "",
        "**Date**: 2026-04-19",
        f"**Scope**: F pilot filter subset (NG=4+AF<0.4+NR≥80 NonLOH) → HCC1395 {mode}-mode subset VCF → inter_sub_mod Phase 2 pipeline → Δ_ASM 分析",
        "**Pipeline**: `r1_run_phase2_both_arms.sh` (parallel TO + paired)；binary: `build/bin/inter_sub_mod`；Normal BAM: HCC1395BL region-subset 1.8GB",
        "",
        f"## 1. 樣本規模",
        "",
        f"- Input subset VCF variants: {801 if mode == 'TO' else 983}",
        f"- Regions with truth_label merged: {df['truth_label'].notna().sum()}",
        f"- TP={tp}, FP={fp}, TP_rate={tp/(tp+fp):.4f}",
        "",
        f"## 2. 全子集 per-feature AUC（TP vs FP）",
        "",
        "| Feature | AUC | n |",
        "|---------|-----|---|",
    ]
    for _, r in rdf.iterrows():
        auc = f"{r['auc']:.4f}" if pd.notna(r['auc']) else "—"
        lines.append(f"| {r['feature']} | {auc} | {r['n']} |")
    lines += [
        "",
        "## 3. F pilot 子集（NR≥80 + NG≥4）per-feature AUC",
        "",
    ]
    if len(frdf) > 0:
        lines += [
            "| Feature | AUC | n |",
            "|---------|-----|---|",
        ]
        for _, r in frdf.iterrows():
            auc = f"{r['auc']:.4f}" if pd.notna(r['auc']) else "—"
            lines.append(f"| {r['feature']} | {auc} | {r['n']} |")
    else:
        lines.append("_Insufficient filter subset data_")

    # Add SampleASM & Residual per label
    lines += [
        "",
        "## 4. Sample ASM / Residual 顯著性",
        "",
    ]
    if "SampleASM_Sig" in df.columns:
        sig_pos = df[df["SampleASM_Sig"].astype(str).str.lower() == "true"]
        lines.append(f"- `SampleASM_Sig=true`：{len(sig_pos)} regions；TP rate={sig_pos['is_tp'].mean():.4f}（vs overall {df['is_tp'].mean():.4f}）")
    if "HP_Residual_Sig" in df.columns:
        res_pos = df[df["HP_Residual_Sig"].astype(str).str.lower() == "true"]
        lines.append(f"- `HP_Residual_Sig=true`：{len(res_pos)} regions；TP rate={res_pos['is_tp'].mean():.4f}")

    lines += [
        "",
        "## 5. 小結",
        "",
    ]
    if mode == "TO":
        lines += [
            "**TO arm 主目標**：Δ_ASM / HP residual 作為 post-filter feature 對 F1 的實際效應。",
            "",
            "若全子集最高 AUC ≥0.60（residualize on NR/AF 後仍 ≥0.60）→ TO arm 有 filter 潛力；反之 → characterization-only。",
            "",
            "**下一步（若成功）**：擴展到其他 5 in-scope 樣本（R12 batch 3）；驗證跨樣本 Δ_ASM AUC 至少 3/6 保持 ≥0.60。",
        ]
    else:
        lines += [
            "**Paired arm 主目標**：Sample ASM characterization，不追求 F1；看 ΔNG / ΔASM 殘差訊號穩定性。",
            "",
            "Paired mode 已飽和（CL-017），本 arm 主要驗證 Sample ASM 在 paired 下是否能獨立於 F1 提供 characterization 訊號。",
            "",
            "**下一步**：將 paired Δ_ASM 與 TO Δ_ASM 在相同 region 子集比較（cross-mode concordance）。",
        ]
    lines += [
        "",
        "## 6. Provenance",
        "",
        f"- Binary: `build/bin/inter_sub_mod`",
        f"- Normal BAM: `/big7_disk/liaoyoyo2001/InterSubMod/output/hcc1395_normal_pilot/HCC1395_normal_subset.bam` (1.8GB, region-subset from 136GB source)",
        f"- Tumor BAM ({mode}): see `r1_run_phase2_both_arms.sh`",
        f"- REF: `/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta`",
        f"- CLI: `-w 5000 -j 16 --distance-metric BERNOULLI`",
        f"- Output dir: `output/hcc1395_normal_pilot/{mode}/`",
        f"- Subset VCF script: `r1_build_subset_vcf_bcftools.sh`（bcftools 1.13）",
        f"- Analysis script: `r1_delta_asm_analysis.py`",
        "",
        "## 7. Registry 連結",
        "",
        "- 對應 CL-013, CL-014（Phase 2 code + B/C/D validation）；R1 batch 2 執行",
        "- 5-axis matrix: `docs/reports/research_landscape/data/10_Registry_5axis_matrix.tsv`",
    ]
    out_md.write_text("\n".join(lines))
    print(f"[+] Report written: {out_md}")


def main():
    for mode in ["TO", "paired"]:
        df, rdf, frdf = analyze_arm(mode)
        out_md = OBS / f"step7_hcc1395_normal_{mode.lower()}_pilot.md"
        out_md.parent.mkdir(parents=True, exist_ok=True)
        write_report(mode, df, rdf, frdf, out_md)


if __name__ == "__main__":
    main()
