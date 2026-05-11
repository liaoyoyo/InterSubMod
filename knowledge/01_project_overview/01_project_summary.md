---
id: ism-kb-01-project-overview-project-summary
name: "ISM Project Summary"
description: "InterSubMod 專案定義：C++17 long-read epigenetic 分析工具；核心價值 read-level epigenetic characterization；論文定位 2026-04 修正。"
status: active
last_verified: 2026-04-22
content_nature: reference-summary
doc_type: explanation
verified_scope: "project positioning against .claude/CLAUDE.md and README_PROJECT_SUMMARY.md"
related_ids:
  - ism-kb-01-project-overview-index
  - ism-kb-01-project-overview-five-research-goals
  - ism-kb-01-project-overview-breakthrough-strategy
tags: [overview, project, summary, positioning]
canonical_paths: [01_project_overview/01_project_summary.md]
alias_paths: []
---

# ISM Project Summary

- 一句結論：InterSubMod 是用 ONT 長讀+甲基化偵測腫瘤亞克隆結構的 C++17 工具；核心價值 = read-level epigenetic characterization（非 variant filter）
- 適用對象：新進者、論文 abstract 撰寫者、跨專案協作
- 可直接執行命令（驗證日期：2026-04-22）：`cat /big7_disk/liaoyoyo2001/InterSubMod/README_PROJECT_SUMMARY.md`

---

## 專案定義

**InterSubMod (Inter-Subclonal Methylation Analysis)**：生物資訊工具，用於 ONT 長讀測序數據偵測腫瘤樣本中的亞克隆結構。整合：
- 甲基化模式（Methylation patterns）
- 體細胞變異（Somatic SNVs）
- 單倍體型（Haplotypes）

分析 epigenetic heterogeneity。

---

## 核心技術特點

| 特點 | 說明 |
|------|------|
| 高效能 C++17 核心 | 結合 OpenMP 平行；單 region 處理 <300ms |
| 精確甲基化解析 | BAM MM/ML tag，精確定位 CpG |
| 多元距離度量 | NHD / L1 / L2 / BERNOULLI / JACCARD / CORR |
| 統計顯著性 | PERMANOVA / Fisher / Cramér's V |
| 自動化視覺化 | 距離熱圖 + 分群熱圖 |

---

## 當前定位（2026-04 確認）

> **ISM 的核心價值在 read-level epigenetic characterization，而非 variant filter。**

**Phase 1A F1 優化已鎖定**：
- paired-pure ΔF1 = +0.0112（95% CI [+0.0044, +0.0188]）
- TO 模式甲基化增益為負（-0.0206）

**當前重心**：Phase 2 生物學特徵化方向（Normal Methylation Reference + CN/Purity-aware correction）

---

## 論文定位

- 目標：**投稿 bioinformatics / epigenetics / cancer genomics** journal
- 重點：read-level epigenetic heterogeneity framework
- 2026-04-18 決策：**分析完整性與對領域貢獻 > 投稿速度**（MEMORY: `feedback_paper_positioning_de_prioritized`）

---

## 核心依賴

| 依賴 | 用途 |
|------|------|
| HTSlib | BAM 處理 |
| OpenMP | 平行運算 |
| Eigen3 | 線性代數 |
| GoogleTest | 單元測試 |
| jemalloc 5.3.0 | 記憶體分配 |
| Python3 + matplotlib/seaborn/scipy/pandas | 視覺化 |

---

## 開發重點

1. **Phase 2 Normal Methylation Reference (方向 A+D)**：當前最高優先
2. **甲基化解析精確度**：MM/ML + CIGAR 校正
3. **距離計算效能**：OpenMP + 稀疏矩陣
4. **統計檢驗**：PERMANOVA + Cramér's V
5. **視覺化品質**：HP tag / Allele 熱圖

---

## 相關

- 五大目標：[02_five_research_goals.md](02_five_research_goals.md)
- 突破策略：[04_breakthrough_strategy.md](04_breakthrough_strategy.md)
- 當前 phase：[03_current_phase.md](03_current_phase.md)
- README：[../../README.md](../../README.md)
- README_PROJECT_SUMMARY：[../../README_PROJECT_SUMMARY.md](../../README_PROJECT_SUMMARY.md)
