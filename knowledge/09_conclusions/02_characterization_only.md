---
id: ism-kb-09-conclusions-characterization-only
name: "Characterization Only Findings"
description: "能描述但不能做 variant filter 的發現：HPFineNGroups filter NEGATIVE、Zone-Aware QS NEGATIVE、F1-filter 方向放棄。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: reference
verified_scope: "characterization-only findings against MEMORY"
related_ids:
  - ism-kb-09-conclusions-index
  - ism-kb-07-derived-features-hpfinengroups
  - ism-kb-07-derived-features-zone-aware-framework
  - ism-kb-09-conclusions-concluded-negative
tags: [conclusions, characterization, index]
canonical_paths: [09_conclusions/02_characterization_only.md]
alias_paths: []
---

# Characterization Only Findings

- 一句結論：**能描述真實現象，但用作 variant filter 會產生負 F1 增益**；這是 ISM 最重要的 caveat 類結論
- 適用對象：避免把 characterization finding 當 filter 誤用
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  grep -l "characterization" /bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/*.md | head
  ```

---

## Characterization Only 清單

### 🟡 C1：HPFineNGroups Subclone Marker（filter 方向 NEGATIVE）
- **描述**：N≥4 region TP rate 高於 baseline（characterization POSITIVE）
- **但**：用作 `HPFineNGroups <= N` 作 FP filter → 負 F1 增益
- **原因**：subclone 結構在 TP 中也有，filter 傷害 recall
- **KB 詳細**：[../07_derived_features/01_hpfinengroups.md](../07_derived_features/01_hpfinengroups.md)

### 🟡 C2：Zone-Aware Framework（QS 調整 NEGATIVE）
- **描述**：Zone TP rate 差異真實（H1 POSITIVE）
- **但**：QS simulation ΔF1 無顯著改善（H3 NEGATIVE）
- **結論**：zones 是真實，但 QS-based filter 無效
- **KB 詳細**：[../07_derived_features/03_zone_aware_framework.md](../07_derived_features/03_zone_aware_framework.md)

### 🟡 C3：F1-filter 方向放棄（2026-04-21 決策）
- **描述**：Fisher_Frac_Sig 等 F1-filter 候選全部 NEGATIVE
- **理由**：
  1. Fisher_Frac_Sig CI 跨隨機
  2. F pilot TP 99.5% 飽和
  3. characterization-only 須全域 region
- **KB 詳細**：[../07_derived_features/05_fisher_frac_sig.md](../07_derived_features/05_fisher_frac_sig.md)
- **MEMORY**：`project_paired_f1_filter_abandoned`

### 🟡 C4：Read-Level Germline FP（CONDITIONAL NO-GO）
- **描述**：LOSO AUC 0.721（有判別力）
- **但**：FP removal=0%（無實際過濾效果）
- **場景**：低純度樣本可能有潛力
- **MEMORY**：`project_read_level_germline_fp`

### 🟡 C5：Per-CpG ASM 與 Epiallele 指標
- **描述**：2026-04-15 characterization positive
- **但**：無法作 variant filter
- **狀態**：特徵化用途

---

## 核心原則

對所有 characterization-only finding：

✅ **可以做**：
- 論文 characterization results
- 描述生物學現象
- 作為分析特徵（不是 filter）

❌ **不可做**：
- 設計 variant filter 基於此特徵
- 宣稱可提升 F1
- 排除 variant（會傷害 recall）

---

## 為何 characterization ≠ filter?

### 根本原因
ISM 特徵反映**真實生物學結構**（subclone、zone、epigenetic heterogeneity），但：
- TP 與 FP **都可能**有這些結構
- Filter 需要 **TP/FP 差異性**，而 ISM 特徵主要反映 **universal biology**

### 實務例子
- HPFineNGroups ≥4：高 TP rate（89-92%）
- 但 ~10% region 是 FP → filter 會誤刪 TP（若 FP 分布相似）
- 實測：用此 filter ΔF1 為負

---

## 相關

- Positive：[01_positive_findings.md](01_positive_findings.md)
- NEGATIVE：[03_concluded_negative.md](03_concluded_negative.md)
- Landscape `03_ISM分析價值界定.md`：[../../docs/reports/research_landscape/03_ISM分析價值界定.md](../../docs/reports/research_landscape/03_ISM分析價值界定.md)
