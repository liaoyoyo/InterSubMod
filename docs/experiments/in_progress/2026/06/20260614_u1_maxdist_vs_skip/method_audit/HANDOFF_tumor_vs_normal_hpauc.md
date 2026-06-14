---
title: "HANDOFF — tumor HP-AUC < normal HP-AUC 檢核觀察報告生成指南（給接手 AI session）"
date: 2026-06-14
type: handoff / observation-guide
audience: 另一 session AI（獨立接手生成檢核觀察報告）
binary_baseline: develop（含 SKIP 預設 + HP-AUC 標準欄, commit 185db26）
data_sources: wg_hpauc.json, hp_auc.json, /tmp/wg_hpauc/significance_summary.csv（會清，可重生）
---

# 接手任務：生成「tumor HP-AUC < normal HP-AUC」檢核觀察報告

> 你（接手的 AI）的目標：產出一份**檢核觀察報告**，驗證並描述「tumor 的 germline-HP 距離區分度
> 系統性低於 normal」這個現象，**並嚴格排除 confound**，最後誠實判定它是不是 somatic 訊號。
> **這是一個 L2 線索，不是定論**——你的報告要幫忙把它推向 L1 或證偽。

## 1. 背景（必讀，否則會誤判）

**HP-AUC 是什麼**：`HP-AUC = P( dist(不同HP的read對) > dist(同HP的read對) )`。
- ground truth = `HP tag`（longphase 給的 **germline haplotype** 標籤，**不是** somatic subclone）。
- 1.0=距離完美區分 HP，0.5=與 HP 無關，−1=undefined（該 subset 沒有「同HP+異HP」兩種對，常見於 tumor 單一 HP）。
- 演算法：`RegionProcessor::compute_hp_auc`（rank-based，NaN/SKIP invalid pair 自動跳過）。
- **已是 ISM 標準輸出欄**：`HP_AUC_Normal` / `HP_AUC_Tumor` / `HP_AUC_All`（significance_summary.csv）。

**核心現象（本 session 全基因組 30490 region 觀察，← wg_hpauc.json）**：
- HP_AUC_Normal: median **0.788**, 76.4% ≥0.7, 僅 1.4% anti → 距離對 germline-HP 強有效（**地基已確認**）。
- **同 both-HP region（n=20449）配對比較**：normal HP-AUC median **0.785** vs tumor **0.641**；
  normal strong(≥0.7) 75.3% vs tumor 39.6%；normal anti 1.7% vs tumor 10.6%。
- **若純 germline，tumor 應 ≈ normal（同 haplotype）；tumor 明顯偏低 = 線索**。

**三個競爭假說（你的報告要逐一處理）**：
- **H_somatic**：tumor 的 germline-HP 甲基結構被 **somatic 甲基改變打亂** → tumor HP-AUC 降低。
- **H_cnv**：tumor CNV/LOH 造成 HP read 不平衡 + noise → HP-AUC 降低（**confound，非 somatic**）。
- **H_coverage/quality**：tumor read 深度/品質差異 → HP-AUC 降低（**confound**）。

## 2. 資料來源 + 重生方式（summary 在 /tmp 會清）

```bash
# 全基因組 SKIP run（重生 significance_summary.csv，~12min）
BIN=<develop build>/build/bin/inter_sub_mod   # 含 SKIP 預設 + HP-AUC 欄
$BIN -t .../HCC1395/tumor.bam -n .../normal.bam -r .../hg38.fa \
     -v .../filtered_snv_tp.vcf.gz -w 5000 -j 16 \
     --distance-metric BERNOULLI --no-output-distance-matrix -o /tmp/wg_hpauc
# (SKIP 已是預設, 不必加 --nan-distance-strategy)
```
關鍵欄位：`HP_AUC_Normal/Tumor/All`, `Tumor_HP1/HP2`, `Normal_HP1/HP2`, `Coverage_Category`,
`Coverage_Multiple`, `Potential_LOH`, `LOH_Subtype`, `NumReads`, `NTumorReads`, `NNormalReads`,
`Significant`, `CramersV`。

## 3. 報告必涵蓋的觀察指標

### 3a 總體分布（已有 `wg_hpauc_analyze.py`，env U1_SKIP 跑即得）
- HP_AUC_Normal/Tumor/All 的 median / mean / frac(≥0.7,0.6-0.7,0.5-0.6,<0.5) / undefined%。

### 3b 配對 delta（核心）— per-region `normal − tumor`
- 只取 both-HP region（`Tumor_HP1>0 & Tumor_HP2>0` 且兩個 AUC 都 ≥0）。
- 算 `delta = HP_AUC_Normal − HP_AUC_Tumor` 分布（median, %delta>0, %delta>0.2）。
- delta>0 普遍 = tumor 系統性低於 normal（量化現象強度）。

### 3c 🔴 confound 分層（**最重要，決定是不是 somatic**）
把 delta 按以下分層，看 tumor<normal 是否在「乾淨」層仍成立：
- **CNV**：`Coverage_Category` ∈ {Normal, Elevated, CNV_Gain, CNV_Loss, High_Copy, Low}。
  → 若 tumor<normal 只在 CNV 層、在 `Normal` 層消失 → **H_cnv（confound），非 somatic**。
  → 若在 `Normal`(diploid) 層 tumor 仍 < normal → 排除 CNV，支持 H_somatic。
- **LOH**：`Potential_LOH==True` 的 region 單獨看（LOH 本身打亂 HP）。
- **coverage**：tumor read 數 vs normal read 數（`NTumorReads` vs `NNormalReads`）；
  控制 read 數後 delta 是否仍在（避免「tumor read 少 → AUC 估計差」假象）。
- **HP 平衡**：`Tumor_HP1/HP2` 比例（極不平衡 → AUC 不穩）。

### 3d 與顯著性/CramersV 的關係
- tumor HP-AUC 低的 region，是否更/更不顯著（Significant/CramersV）？

## 4. 案例選取（怎麼挑代表案例畫圖）

挑「明顯 tumor << normal 且 confound 乾淨」的代表 region：
- 條件：`HP_AUC_Normal ≥0.8` 且 `HP_AUC_Tumor ≤0.5` 且 both-HP 且 `Coverage_Category==Normal`
  且 `Potential_LOH==False`（= diploid 乾淨層仍 tumor 崩 → 最強 somatic 候選）。
- **對照案例**：再挑一個 `CNV_Gain` 或 `LOH` 的 tumor<normal（展示 confound 版本，對比）。
- 每染色體各挑 1-2 個，避免 chr1 偏差（見樹見林）。

**畫圖**（範本 `plot_cases.py`，需重跑該 region 開 distance matrix + methylation 輸出）：
- distance heatmap：read 按 (is_tumor, HP) 排序 → 看 tumor block 是否比 normal block 鬆散/混亂。
- methylation heatmap：tumor reads vs normal reads 分層 → 看 tumor 的 HP 甲基模式是否被打亂。
- 建小 VCF（選定 SNV）重跑（不加 `--no-output-distance-matrix`），參考 hp_auc.py 的 case 載入邏輯。

## 5. 報告誠實鐵則（必標，否則會過度宣稱）
- **HP = germline 真值**：HP-AUC 高/低都是「對 germline-HP 結構」的度量，**不是 somatic subclone 直接證據**。
- **normal>tumor 是 L2 線索**：只有在 **3c confound 分層後**（diploid 乾淨層仍 tumor<normal）才可升 L1-somatic；
  否則維持「可能 CNV/LOH/coverage confound」。
- **單樣本 HCC1395 paired**：跨樣本（COLO829 等）未驗 → 結論標 partial。
- 每個數字標來源（significance_summary 欄 / wg_hpauc.json），grep 得到才寫（[[feedback_no_fabricated_numbers_in_reports]]）。

## 6. 現有產物（接手可直接用 / 參考）
- `wg_hpauc_analyze.py` + `wg_hpauc.json` — 全基因組分布 + 部分分層（已含 tumor_both_hp 配對）。
- `hp_auc.py` + `hp_auc.json` — case 級 HP-AUC（含 C++/python 交叉驗證、tumor-only=null 範例）。
- `plot_cases.py` + `figures/` — distance/methylation 熱圖範本（4 case）。
- `tumor_landscape.py` + `tumor_landscape.json` — tumor 單-HP 32% / LOH 8.5% / CNV 分布（confound 基數）。
- 上游脈絡：`../README.md`（U1）、`../gap2/README.md`（germline 主導 + hp_residual 缺陷）、
  `InterSubMod/docs/plans/20260614_ism_method_validation_plan_01.md`（P1-P5 計劃）。

## 7. 交付物（你的報告應產出）
1. `tumor_vs_normal_hpauc_report.md`（L0 結論 + 3a/3b/3c 表 + 3-5 案例圖 + caveat）。
2. confound 分層判定：tumor<normal 在 diploid 乾淨層**是否仍成立** → H_somatic vs H_cnv 結論。
3. 可重現 script（延續 env U1_SKIP 慣例）。
