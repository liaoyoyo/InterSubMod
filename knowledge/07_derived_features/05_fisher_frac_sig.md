---
id: ism-kb-07-derived-features-fisher-frac-sig
name: "Fisher_Frac_Sig"
description: "F1-filter 方向已放棄的衍生信號；Fisher_Frac_Sig CI 跨隨機 + F pilot TP 99.5% 飽和；2026-04-21 characterization-only 決策。"
status: active
last_verified: 2026-04-22
content_nature: historical-note
doc_type: reference
verified_scope: "Fisher_Frac_Sig abandonment decision 2026-04-21"
related_ids:
  - ism-kb-07-derived-features-index
  - ism-kb-09-conclusions-concluded-negative
tags: [features, fisher, f1-filter, abandoned, historical]
canonical_paths: [07_derived_features/05_fisher_frac_sig.md]
alias_paths: []
---

# Fisher_Frac_Sig

- 一句結論：🔴 F1-filter 方向已放棄（2026-04-21 決策）；Fisher_Frac_Sig CI 跨隨機 + F pilot TP 99.5% 飽和；characterization-only
- 適用對象：理解已放棄方向、避免重新調查
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  grep -l "Fisher_Frac" /big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/ -r | head
  ```

---

## 定義

**Fisher_Frac_Sig**：Fisher exact test 顯著比例（per-region 或 per-sample 聚合）

**計算**：對多重檢定調整後，顯著位點占總位點的比例

---

## 🔴 已放棄決策（2026-04-21）

**決策**：放棄 F1-filter 方向；僅保留 characterization 用途

**理由**：
1. **Fisher_Frac_Sig CI 跨隨機**：Confidence interval 涵蓋隨機值
2. **F pilot TP 99.5% 飽和**：TP rate 極高但無法選擇性過濾 FP
3. **characterization-only 須全域 region**：不適合 variant-specific filter

**MEMORY**：`project_paired_f1_filter_abandoned`

---

## 為何放棄

### 問題 1：CI 跨隨機
- 隨機 baseline 的 Fisher_Frac_Sig 與觀測值重疊
- 信號強度 < 隨機擾動

### 問題 2：飽和
- F pilot 中 TP rate 達 99.5%
- 無 gradient 可用來區分 TP / FP

### 問題 3：scope mismatch
- 計算在全域 region 才有意義
- 對單一 variant 的 ISM filter 無 selectivity

---

## 歷史用途

- **舊版 ISM filter 候選**：曾作為 F1 提升候選
- **2026-04-21 前**：嘗試多個 threshold 變體
- **2026-04-21 後**：僅保留在 `significance_summary.csv` 供 characterization

---

## 不要做的事

- ❌ 基於 Fisher_Frac_Sig 設計新的 F1-filter
- ❌ 使用 Fisher_Frac_Sig 作 variant-specific 過濾
- ❌ 與其他 sig metric 並用計算 composite score（已證混淆）

---

## 仍可做的事

- ✅ Characterization：跨樣本 signal stability 分析
- ✅ Debug：診斷 pipeline 不穩定
- ✅ Reference：供新 feature 設計時比較 baseline

---

## 相關

- 索引：[00_index.md](00_index.md)
- NEGATIVE 目錄：[../09_conclusions/03_concluded_negative.md](../09_conclusions/03_concluded_negative.md)
- CURRENT_FOCUS：[../../docs/CURRENT_FOCUS.md](../../docs/CURRENT_FOCUS.md)
