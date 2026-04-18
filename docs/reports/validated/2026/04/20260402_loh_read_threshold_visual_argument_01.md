<!--
建立時間: 2026-04-02 14:00
目標: 視覺化論證 ISM LOH 判定缺少最低 read 門檻的影響與根因分析
處理範圍: 6 張分析圖表 + 數據表 + 雙問題分解結論
關聯檔案:
  - scripts/analysis/build_loh_read_threshold_analysis.py
  - src/core/RegionProcessor.cpp (is_potential_loh, compute_hp_ratio)
  - docs/reports/validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md
-->

# LOH 判定無最低 Read 門檻 — 視覺化論證報告

> **結論預覽**：ISM `is_potential_loh()` 確實缺少最低 read count 門檻，但此缺陷僅影響 ~3% LOH。97% 的過量 LOH 源自 LongPhase-TO self-phasing scaffold 汙染，在所有 read count tier（含 200+ reads）穩定存在。

---

## 1. 程式碼缺陷確認

```cpp
// src/core/RegionProcessor.cpp
double compute_hp_ratio(int hp1_family_n, int hp2_family_n) {
    double total = static_cast<double>(hp1_family_n + hp2_family_n) + 0.002;
    return (static_cast<double>(hp1_family_n) + 0.001) / total;
}

bool is_potential_loh(double hp_ratio) {
    return (hp_ratio < 0.1) || (hp_ratio > 0.9);  // ← NO read count check
}
```

**問題**：只檢查 HP ratio 閾值（<0.1 或 >0.9），**完全沒有最低 read 數量要求**。當 effective_hp=1 時，HP_Ratio 必然為 0 或 1 → 100% 被判為 LOH。

---

## 2. Fig01 — LOH Rate vs Effective HP Reads（核心論證）

![Fig01: LOH Rate vs Effective HP Reads](figures/20260402_loh_read_threshold_figures/fig01_loh_rate_vs_effective_hp.png)

三條線揭示完整問題結構：

| 線條 | 含義 | 關鍵觀察 |
|------|------|---------|
| 綠色虛線（Binomial null） | 隨機 HP 分配下的期望 LOH rate | ehp=1 時 100%，ehp≥10 趨近 0% |
| 紅線（TO observed） | 實際 TO LOH rate | **全區間遠高於 null**，ehp=200+ 仍 27.0% |
| 藍線（Paired observed） | Paired LOH rate | 高 read 區（200+）僅 1.9% |

**下方柱狀圖**：TO excess over null 在所有 tier 穩定 +36-64% → 不是低 read 噪音，是 **self-phasing scaffold 系統性汙染**。

---

## 3. Fig02 — HP_Ratio 散佈圖（8 個 ehp tier）

![Fig02: HP_Ratio scatter by ehp tier](figures/20260402_loh_read_threshold_figures/fig02_hp_ratio_scatter_by_ehp_tier.png)

8 面板從 ehp=1 → ehp≥200 展示 HP_Ratio 分佈演變：

- **ehp=1**：所有點鎖死在 0 或 1（1 條 read → 必然全在一個 HP）
- **ehp=2-4**：大量 0/1 極端值，幾乎無中間值
- **ehp≥50**：散佈展開，但 LOH 區（紅色背景）仍有大量點
- **ehp≥200**：read 充足，Paired 下 LOH 幾乎消失（1.9%），TO 仍 27.0%

---

## 4. Fig03 — LOH Excess Heatmap（7 samples x 8 tiers）

![Fig03: LOH excess heatmap](figures/20260402_loh_read_threshold_figures/fig03_loh_excess_heatmap_sample_ehp.png)

跨樣本一致性驗證：

- **全部 7 個樣本**在所有 ehp tier 都呈現深紅色（高 LOH rate）
- 藍色數字（頂軸）= Paired 參考值，普遍低於 TO
- ehp 200+ tier：所有樣本 TO LOH 19-44%，Paired 僅 0-6%
- **一致性 7/7** → 非個別樣本特異性

---

## 5. Fig04 — Cumulative LOH Contribution

![Fig04: Cumulative LOH contribution](figures/20260402_loh_read_threshold_figures/fig04_cumulative_loh_contribution.png)

Lorenz 曲線風格展示低 read 的 LOH 佔比：

| 門檻 | 佔 TO regions | 佔 TO LOH | 解讀 |
|------|-------------|----------|------|
| ehp < 10 | 1.7% | 2.7% | 略不成比例，但**影響極小** |
| ehp < 30 | 11.7% | 15.7% | 中等不成比例 |
| ehp < 50 | 28.2% | 37.4% | 含大量 self-phasing LOH |

**關鍵**：低 read LOH 僅佔 ~3%。加 read 門檻只解決冰山一角，**97% LOH 來自中高 read 區的 self-phasing**。

---

## 6. Fig05 — 高 Read 區 TO vs Paired（決定性論證）

![Fig05: High-read LOH TO vs Paired](figures/20260402_loh_read_threshold_figures/fig05_high_read_loh_to_vs_paired.png)

三面板展示 read 充足區域的 LOH rate：

| ehp Tier | TO TP | TO FP | Paired TP | Paired FP |
|----------|-------|-------|-----------|-----------|
| 50-99 | 44.4% | 30.2% | 33.8% | 54.7% |
| 100-199 | 26.4% | 29.5% | 11.9% | 61.0% |
| **200+** | **27.1%** | **26.7%** | **1.9%** | **0.0%** |

**ehp≥200 是決定性證據**：
- Read 完全充足（200+ reads），Binomial noise = 0
- Paired 幾乎無 LOH（1.9%）
- TO 仍有 27% LOH
- **唯一解釋：TO phasing scaffold 品質問題**

---

## 7. Fig06 — 雙問題分解（最終結論）

![Fig06: Two-problem decomposition](../../../../../output/synthesis/observation_workspaces/20260402_loh_read_threshold_analysis/figures/fig06_two_problem_decomposition.png)

### 左圖 — Problem A: Low-Read Noise

- TO observed vs Binomial null 並列
- ehp=1 兩者相同（100%）→ 數學必然，非 bug
- ehp 2-4 起差距出現（86% vs 25%）→ **excess = self-phasing 造成**
- 所有 tier 的 excess：+61%, +67%, +64%, +49%, +55%, +40%, +27%, +27%

### 右圖 — Problem B: Self-Phasing Scaffold

- TO vs Paired 並列
- **所有 tier 都有 TO >> Paired**
- ehp 200+：TO 27% vs Paired 2%（**14.0x**）
- TO/Paired 倍率隨 read count 增加反而**擴大**（非收斂）

---

## 8. 量化總結

| ehp Tier | TO N | TO LOH% | Paired LOH% | Binomial Null% | TO Excess |
|----------|------|---------|-------------|----------------|-----------|
| 0 | 1,070 | 0.0% | 0.0% | 0.0% | 0.0% |
| 1 | 547 | 100.0% | 100.0% | 100.0% | 0.0% |
| 2-4 | 1,861 | 86.1% | 74.2% | 25.0% | +61.1% |
| 5-9 | 3,703 | 68.2% | 43.2% | 1.6% | +66.6% |
| 10-19 | 15,571 | 64.5% | 49.3% | 0.1% | +64.4% |
| 20-29 | 26,357 | 48.8% | 33.8% | 0.0% | +48.8% |
| 30-49 | 69,432 | 54.8% | 49.0% | 0.0% | +54.8% |
| 50-99 | 218,719 | 40.0% | 33.9% | 0.0% | +40.0% |
| 100-199 | 74,511 | 27.1% | 12.0% | 0.0% | +27.1% |
| 200+ | 7,921 | 27.0% | 1.9% | 0.0% | +27.0% |

---

## 9. 結論：雙問題定性

| 問題 | 描述 | 影響範圍 | 修復難度 | 優先級 |
|------|------|---------|---------|--------|
| **Problem A** | `is_potential_loh()` 無 read 門檻 | ~2.7% LOH（ehp<10） | 簡單（加 `ehp >= 10`） | 低 |
| **Problem B** | Self-phasing scaffold 汙染 | ~97% 過量 LOH | 困難（LongPhase-TO 架構） | **高** |

### 行動建議

1. **Quick Fix（Problem A）**：在 `is_potential_loh()` 加入 `effective_hp >= 10` 下限 → 移除 ~3% trivial LOH
2. **Root Cause（Problem B）**：從 phasing scaffold 入手改善 TO haplotagging 品質 → 這是唯一能根治 97% 過量 LOH 的路徑
3. **已確認**：即使加了 read 門檻，ehp≥200 的 TO LOH 仍為 Paired 的 14 倍 → read count 不是根因

---

*分析腳本*: `scripts/analysis/build_loh_read_threshold_analysis.py`
*數據輸出*: `big7_disk_output/synthesis/observation_workspaces/20260402_loh_read_threshold_analysis/`
*前置報告*: `20260402_longphase_to_vs_s_causal_chain_report_01.md`
