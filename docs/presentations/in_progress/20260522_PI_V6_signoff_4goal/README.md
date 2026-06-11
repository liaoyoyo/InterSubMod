# PI V6 Signoff — 4 Goal Hybrid HTML Report

> **建立日期**: 2026-05-21
> **目標受眾**: PI (1-on-1 / Lab meeting)
> **狀態**: in_progress (waiting for Phase A 補測完成)

## 結構

```
20260522_PI_V6_signoff_4goal/
├── index.html                                   ← 頂層 hub + executive summary
├── target1_V6_vs_baseline_tagbam.html           ← 目標 1: BAM 層 HP 譜系比較
├── target2_ISM_V6_vs_baseline.html              ← 目標 2: ISM region 層比較
├── target3_5features_distributions.html         ← 目標 3: 5 features 獨立+組合
├── target4_TO_pure_F1_methodology.html          ← 目標 4: TO 純樣本 F1 方法論
├── assets/
│   ├── figures/{target1,target2,target3,target4}/
│   ├── data/{target1,target3,target4}/
│   ├── js/   (Plotly + DataTables + Alpine.js inline)
│   └── css/  (MVP.css + custom)
└── README.md                                     ← 本檔
```

## 4 個目標摘要

| Target | 內容 | A0 audit verdict | A 補測 |
|---|---|---|---|
| **1** | V6 vs baseline tag.bam HP 6 值分佈 + SP1/2/3 + priority bug 17.3:1 | 既有素材 PASS minor (Doc 1 §8.2 row-mislabel, 17to1 §4 ~95% framing) | A1 (HP 6 值) + A2 (ALT-only ratio) |
| **2** | ISM region 層 HCC1395 4-way + 5 樣本 caveat | 既有素材 PASS 100% (20260519 110/110) | 重用既有 |
| **3** | 5 features individual + 9-cell heatmap + cross-sample direction | (建構中) | A3 (5 features 全 AUC + 9-cell + Spearman ρ) |
| **4** | TO 純樣本 F1 4-algorithm LOSO + 方法論 sanity check | (建構中) | A4 (LR/DT/RF/XGB LOSO) + A4-Ext (其他算法盤點) + A4-F1Audit (F1 step audit) |

## 引用 ground rules

詳見 `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A0_existing_artifacts_verification.md` §4 表格。重點：
- **避免引用**: 20260514 priority bug eng .md §8.2 (row-mislabel)、20260514 17to1 HTML §4 "~95% GT2" framing、4-Layer "Phase N50 ~20,000+" (用 verified 8,109)
- **直接引用**: PI signoff email (22/22 verified), phase2 trust framework HTML (LOSO 全 verified), 20260519 V6 vs baseline HTML (110/110), 20260515 V6 TPFP characterization HTML (14/14)
- master TSV col count = **112** (不是 116)

## JS 套件 inline (offline standalone)

- Plotly.js 2.27.0 (互動 boxplot / heatmap / scatter)
- DataTables 1.13.7 (raw data TSV table)
- Alpine.js 3.13 (tab / collapse / expand)
- MVP.css (base style)

## 4 級敘述分層 (每子 HTML 通用)

- L0 ⭐⭐⭐ Headline (≤ 50 字, 永遠顯示)
- L1 ⭐⭐ Top 3 Findings (3-5 條 + tier badge + 1 圖)
- L2 ⭐ Evidence Cards (默認 collapsed, click expand)
- L3 Raw Data (DataTables, hidden, toggle)
