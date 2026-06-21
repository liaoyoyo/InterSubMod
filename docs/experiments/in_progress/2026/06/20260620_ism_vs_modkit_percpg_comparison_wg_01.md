---
title: ISM joint/cluster vs modkit-equiv per-CpG marginal — 全基因組對照（HCC1395）
date: 2026-06-20
status: in_progress
tier: 3
sample: HCC1395（單樣本 ⭐2-3 characterization）
scope: 全基因組 22 autosome TP SNV，N=30,490（100% full）
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/percpg_compare_summary.json
build_commit: e57a7ab
observation_standard: true
---

# ISM(read 層 joint/cluster)vs modkit-equiv(per-CpG marginal)— 全基因組對照

> **核心問題**：modkit 式逐-CpG 率差能涵蓋多少 ISM 抓到的訊號？ISM 的 read 層 joint 多抓多少？
> **方法**：modkit = `per-CpG-by-label`（Fisher，= modkit `dmr` marginal 檢定本質；binary 此 checkout 缺，equiv 忠實重現）；ISM = `joint PERMANOVA`（read×read BERNOULLI 距離）+ 無監督 cluster。BH-FDR：per-CpG within-locus、joint across-loci per axis。
> **HTML**：`InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/20260620_ism_vs_modkit_percpg_comparison_01.standalone.html`

## 0. modkit 實際能力（為何只能用 equiv，不能直接驗證）

modkit `pileup --phased` 依 BAM **HP tag 分區** → 每區 bedMethyl（per-CpG 分數）；`dmr pair -a -b` 做兩區**逐位點 pairwise 率差**（beta-binomial）。本質 = **per-CpG marginal + 兩區 pairwise**。
- **hp 軸（HP1 vs HP2）**：原生分區 ✓。
- **carrier 軸（germline{1,2} vs somatic{1-1,2-1}）**：是 HP 值的**合併**，須先 retag/合併 bedMethyl，非原生。
- **allele 軸（REF vs ALT）**：**不是 HP tag**，須先給 read 加自訂 allele tag（= ISM 內部 alt_support）。
- modkit **只做 per-CpG marginal pairwise，無 read 層 joint / 無監督 clustering**。
→ 為何不能用 modkit 驗證：(1) per-CpG 部分是同一 marginal 計算（近循環）；(2) carrier/allele 軸要先用我們標籤 retag（非獨立）；(3) ISM 的 joint/clustering modkit 無對應功能。equiv（Fisher）涵蓋 modkit 唯一能貢獻那塊；真 binary 只在 beta-binomial 細節 + runtime 多一點，不改結論。

## 1. 4-cell Venn（K=1：≥1 顯著 CpG = modkit 偵測到；N=30,490）

| 軸 | 可比 n | both（modkit也抓） | ISM-only（joint顯著,per-CpG貧乏） | modkit-only（per-CpG有,joint不顯著） | neither |
|---|---:|---:|---:|---:|---:|
| hp | 16,510 | 12,201 | 1,152（6.98%） | 611（3.70%） | 2,546 |
| carrier | 24,191 | 16,451 | 2,020（8.35%） | 1,228（5.08%） | 4,492 |
| allele | 24,695 | 18,773 | 1,564（6.33%） | 1,085（4.39%） | 3,273 |
| **整體 any-軸** | 30,490 | **19,637** | **1,593（5.22%）** | **878（2.88%）** | 8,382 |

ISM-only / modkit-only = **1.81×**。

## 2. 各 ISM state 中「modkit per-CpG 一無所獲（=0）」比例

| ISM state | 位點數 | per-CpG 顯著數 中位 | per-CpG=0（modkit 找不到） |
|---|---:|---:|---:|
| ⑤ 對齊真結構 | 6,111 | 17 | 10（0.16%） |
| ④ 可分未對齊 | 13,390 | 4 | 4,473（**33.41%**） |
| ③ 切不出但 joint 顯著 | 5,779 | 2 | 902（**15.61%**） |
| ② 無訊號 | 4,797 | 0 | 4,177（87.08%） |
| ① 不可驗證 | 413 | 0 | 413（100%） |

## 3. 解讀（誠實，論文口徑）

1. **modkit 式 per-CpG 涵蓋多數** — both=19,637 是有訊號位點絕大部分；強 ASM 焦點訊號（BRCA2 型）modkit 也抓得到（⑤ per-CpG=0 僅 0.16%、中位 17 顯著 CpG）。**per-CpG 定位不是 ISM 創新。**
2. **ISM read 層 joint 淨多抓 5.22%（1,593 位點）** > modkit-only 2.88%（878）≈ 1.8×。這些是「分散弱訊號跨多 CpG」（單點各自不顯著、聯合距離抓得到）。
3. **最清楚的 ISM 優勢證據：③（切不出但 joint 顯著）15.61% 位點 per-CpG 完全=0** → ISM 說「有結構」、modkit 逐點一無所獲。
4. ⚠ **④ 未對齊 33.41% per-CpG=0** → 三分之一「可分未對齊」切群無 per-CpG 支持 = 雜訊嫌疑（呼應 [[project_decisionflow_5state_classification_wg]] ④ 桶異質）。
5. **此 Venn 只比 supervised per-CpG vs joint**；ISM 無監督 clustering（⑤/④ 把 read 切群、發現未知 epiallele 群）modkit 結構上無對應 → 另一層 ISM-only 能力，未計入此表。

**一句話定位**：modkit/DSS 式逐-CpG 率差涵蓋多數已知 ASM 焦點；ISM 賣點 = **read-level joint（多抓 ~5% 分散弱訊號，尤其 ③ 中 16% modkit 全漏）+ 無監督 epiallele 群發現**，不是 per-CpG 定位本身。

## 3b. S3 per-CpG=0 dual-panel 肉眼確認（2026-06-20，7 案抽樣 + 用戶判讀）

從 567 個「S3 + per-CpG=0」（location-clean 472 / dispersion 95）抽 7 案做 dual-panel（HTML `20260620_s3_percpg0_joint_observation_01.standalone.html`），依 joint-sig 軸排序：
- **location-clean（4 案，dispP≥.05；如 chr17:72044253 allele F=5.084 dispP=0.91）**：甲基近均勻、無單一 CpG 分群（=per-CpG=0 對），但 read×read 距離依顯著軸有**微弱塊** = diffuse 水平差，modkit marginal 漏。
- **dispersion（3 案，dispP<.05；如 chr7:116324245 allele F=7.23 dispP=0.01）**：標籤群**離散度**差（+可能 location）。

**🔧 框架修正（用戶判讀後）**：用戶肉眼判 7 案**皆 label-consistent 真結構**（外部標籤一致對齊距離塊＝非循環證據）。據此修正先前「dispersion=modkit 對」的過度簡化 → **dispersion 案是 ISM 偵測到「異質性結構」（一群 read 甲基更雜），modkit per-CpG 結構上看不到此維度，故仍是 ISM-only 價值**；惟須區分生物異質（subclone 混 epiallele）vs 技術 spread（覆蓋/品質）。

**強化口徑**：ISM 偵測 read-level label-associated 結構含三種 modkit 皆做不到者 — (1) diffuse location 差、(2) heterogeneity/dispersion 差、(3) 無監督群發現。🔴 仍守邊界：訊號弱（F 2.4-7.2 ≪ BRCA2 15.4）、單樣本 ⭐2-3、「有結構」≠ subclone；per-CpG 定位本身 modkit 夠用，ISM 賣點是結構偵測非定位。

## 4. 驗證表

| 數字 | 值 | 來源 key | 重算 | L |
|---|---|---|---|---|
| 合併位點 | 30,490 | summary.n_merged | percpg∩decisionflow | L1 |
| both（any-軸） | 19,637 | overall_K1.both | joint-sig & percpg-sig | L1 |
| ISM-only | 1,593（5.22%） | overall_K1.ism_only | joint-sig & ~percpg-sig | L1 |
| modkit-only | 878（2.88%） | overall_K1.mod_only | percpg-sig & ~joint-sig | L1 |
| carrier ISM-only | 2,020（8.35%） | venn_K1.carrier.ism_only | per-軸 | L1 |
| S3 per-CpG=0 | 902（15.61%） | by_state.S3.percpg_zero | max-axis n_sig==0 | L1 |
| S4 per-CpG=0 | 4,473（33.41%） | by_state.S4.percpg_zero | 同上 | L1 |
| S5 per-CpG=0 | 10（0.16%） | by_state.S5.percpg_zero | 同上 | L1 |

## 5. 限制

- **modkit-equiv（Fisher）非真 modkit binary**：本質同 marginal 檢定；beta-binomial vs Fisher 的尾部細節 + runtime 未測（binary 此 checkout 缺）。
- **單樣本 HCC1395 ⭐2-3**；joint/per-CpG 皆 cis-ASM characterization 非 subclone。
- **K=1 門檻**（≥1 顯著 CpG = modkit 偵測）；K=3 更嚴（HTML 內附，ISM-only 比例略升）。
- compute 因共享機器 CPU 競爭拖至 128min（結果不受影響，err=0）。

## 6. 重生
```bash
cd InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot
python3 scripts/percpg_wg.py        # ~所有位點逐CpG三軸Fisher+BH → percpg_records.json
python3 scripts/percpg_compare.py   # merge decisionflow joint → percpg_compare_summary.json
python3 scripts/build_cmp_report.py # → 對照 HTML
```
