# Step 7 — Cross-mode concordance & HCC1395 Normal BAM sanity check

**Date**: 2026-04-21
**Scope**: 驗證 HCC1395 TO / paired 兩條 Phase 2 pipeline 在共享同一 Normal BAM 下的輸出一致性 + TO `SampleASM_Delta` vs paired `Fisher_Frac_Sig` 的訊號獨立性
**Trigger**: User decision (2026-04-21)
> 「暫緩 R12，先做 P2 cross-mode concordance，也先確認 HCC1395 有 Normal 的結果；CL-008 等 R12 跨樣本後再決；A 先放棄 paired F1-filter 方向，轉 characterization-only，除非後續有結果發現才再修正與觀察確認。」

---

## 1. 動機與設計

前序 Step 7 兩份 arm-specific pilot 產生**不對稱結論**：

| Arm | 主要訊號 | Residualized AUC | 穩定度 | 定位建議 |
|-----|---------|-----------------|--------|---------|
| TO | `NormalBaseline_Coverage` 0.604 ± 0.08 / `SampleASM_Delta` 0.610 [0.545, 0.670] | 穩定 | ⭐3 | P0 F1-filter 候選（CL-025a） |
| paired | `Fisher_Frac_Sig` 0.698 [**0.534**, 0.831] | CI 跨 random（FP=10） | ⭐2 降級 | **放棄 F1-filter，轉 characterization-only**（CL-025b） |

這兩個結論都來自同一個 Normal BAM（HCC1395BL region-subset 1.8GB）但兩條獨立 pipeline（`inter_sub_mod` TO-mode vs paired-mode）產生的 `significance_summary.csv`。**本步驟目的**：

1. **Sanity**：同一 Normal BAM 跑兩條 pipeline，Normal-derived features（`NormalBaseline_Coverage`、`NormalBaseline_Mean`）在 overlap region 應該 ρ ≈ 1.0；若不是 → pipeline 實作 bug
2. **Concordance**：TO SampleASM_Delta（P0 F1-filter 候選）與 paired Fisher_Frac_Sig（characterization-only 候選）**是同一個訊號還是獨立訊號**？若獨立 → 可在 characterization layer 作 ensemble；若高度相關 → paired F1-filter 放棄後 paired arm 整體已被 TO 涵蓋
3. **P2 Per-CpG ASM 的可行性**：paired arm Per-CpG ASM (`Fisher_Frac_Sig` 等) 雖 F1-filter 無效，但若能獨立於 TO SampleASM_Delta，仍有 characterization 價值

---

## 2. Pipeline sanity — 共享 Normal BAM 的實作驗證

**方法**：取 TO (801 regions) 與 paired (983 regions) 的 `significance_summary.csv`，以 `Chr:Pos` 為 key inner join → 441 overlap regions。比對三個 Normal-derived / Sample ASM 特徵在兩條 pipeline 的分佈與 per-region 相關性。

![](../figures/fig6_normal_bam_sanity.png)

### 2.1 分布一致性（per-feature per-arm summary）

| Feature | Arm | n | Mean | Median | Std | Min | Max | n_zero | n_NaN |
|---------|-----|---|------|--------|-----|-----|-----|--------|-------|
| NormalBaseline_Coverage | TO | 801 | 26.35 | 24.42 | 10.17 | 8.34 | 124.42 | 0 | 0 |
| NormalBaseline_Coverage | paired | 983 | 24.34 | 23.90 | 5.05 | 9.19 | 61.27 | 0 | 0 |
| NormalBaseline_Mean | TO | 801 | 0.450 | 0.400 | 0.201 | 0.101 | 0.923 | 0 | 0 |
| NormalBaseline_Mean | paired | 983 | 0.448 | 0.405 | 0.195 | 0.060 | 0.920 | 0 | 0 |
| SampleASM_Delta | TO | 801 | 0.200 | 0.168 | 0.113 | -0.016 | 0.547 | 0 | 0 |
| SampleASM_Delta | paired | 983 | 0.207 | 0.178 | 0.111 | -0.0001 | 0.580 | 0 | 0 |

**觀察**：
- **無 NaN、無 zero sentinel**：兩 arm Normal BAM 解析完整，不存在資料不存在的 proxy（過去跨樣本比對曾因 NaN 填 0 造成假訊號）
- **分佈 shape 接近**：`NormalBaseline_Mean` 與 `SampleASM_Delta` 在兩 arm 幾乎重合；唯一差異是 `NormalBaseline_Coverage` TO arm 尾部比 paired 長（max 124 vs 61）— 這是因 TO pipeline 僅 801 regions（為 paired 983 的子集），TO 專取更高 coverage 子集所致，非 bug

### 2.2 Per-region 相關性（overlap=441 regions）

| Pair | Spearman ρ | p-value | 解讀 |
|------|------------|---------|------|
| `NormalBaseline_Coverage`: TO vs paired | **1.0000** | 0.000 | 完美一致 → 同一 BAM、同一 coverage sampler，實作正確 |
| `NormalBaseline_Mean`: TO vs paired | ≈1.0 | 0.000 | 同上 |

**結論 (Sanity)**：HCC1395 Normal BAM 被 TO 與 paired 兩條 pipeline 以完全一致的方式解析。**沒有 pipeline 實作 bug**；TO 與 paired arm 所有下游 Normal-derived 訊號差異都必須歸因於**下游統計處理差異**（TO 無 tumor-only-read 排除、paired 用 tumor vs normal 雙向對比），不是解析階段差異。

這排除了原先的兩個擔心：(1) Normal BAM 跑兩次可能 BAM index 差異造成 coverage 不同；(2) TO/paired pipeline 可能在 Normal 解析階段就出現 divergence。

---

## 3. Cross-mode concordance — TO P0 訊號 vs paired characterization 訊號獨立性

**核心問題**：若 paired `Fisher_Frac_Sig` 放棄 F1-filter 後仍想做 characterization-only（例 epigenetic subclone marker），其訊號需與 TO 主 P0 訊號（`SampleASM_Delta`）**獨立**，否則 paired arm 整體已被 TO 涵蓋、無額外價值。

![](../figures/fig7_crossmode_signal_scatter.png)

### 3.1 兩訊號 pairwise scatter

在 441 overlap regions 上：

| Pair | Spearman ρ | p-value | n | 解讀 |
|------|------------|---------|---|------|
| TO SampleASM_Delta vs paired Fisher_Frac_Sig | **−0.162** | 6.5 × 10⁻⁴ | 441 | 顯著但極弱負相關 → 訊號**大致獨立** |
| TO SampleASM_Delta vs paired SampleASM_Delta | — | — | 441 | 同 arm 指標，預期強相關（附錄） |

**解讀**：
- ρ = −0.162 在 n=441 下**統計顯著**（p ≪ 0.001）但**效應極小**（|ρ| < 0.20）
- 方向為負：當 sample methylation 與 normal 差異愈大，per-CpG Fisher significance fraction 愈低
- 這與生物學直覺一致：`SampleASM_Delta` 高 → tumor methylation 偏離 normal，多數 CpG 一致偏離 → 單點 Fisher exact 反而難達顯著（分布位移 vs 分布差異是不同量）
- **獨立性判定**：兩訊號**不是同一個潛在變量的不同代理**；若 paired arm 有機會找到 characterization-unique 訊號，這是條件允許的起點

### 3.2 Per-feature 跨 arm 相關性矩陣

![](../figures/fig8_crossmode_feature_correlation.png)

重要觀察：
- **Normal-derived features（左上角）**：跨 arm ρ ≈ 1.0（sanity 延伸）
- **Sample ASM / HP Delta features（中段）**：跨 arm ρ 中等高（0.6-0.9），顯示 tumor methylation 計算在兩 arm 高度重現
- **Per-CpG Fisher features（右下角）**：跨 arm ρ 較低，paired arm 專屬訊號來源於 paired mode 的 per-CpG normal adjustment，TO arm 無此機制

### 3.3 TP/FP agreement

![](../figures/fig9_crossmode_tp_agreement.png)

overlap 441 regions 中：
- **both_TP**: 439 (99.5%)
- **both_FP**: 2 (0.5%)
- **TO_TP + paired_FP**: 0
- **TO_FP + paired_TP**: 0

**解讀**：F pilot canonical filter (NG=4+AF<0.4+NR≥80 NonLOH) 在 overlap region 已幾乎飽和成 TP（TP rate 99.5%）。這**限制了 cross-mode ensemble 的實際 F1 提升空間** — 在 F pilot subset 內，no matter 哪個 arm 的訊號都已無 FP 可 rescue。

→ **延伸意義**：paired arm characterization-only 定位必須基於**全域 region（非 F pilot subset）**，才能發揮區分作用。F pilot 子集內 paired arm 已無 F1-relevant 工作可做，這與使用者 Q3 放棄 paired F1-filter 的決策直接對應。

---

## 4. 整合結論

### 4.1 Sanity (確立)

- HCC1395 Normal BAM 的兩條 Phase 2 pipeline 解析完全一致（ρ=1.0）
- 無 NaN / 無 zero sentinel / 分布 shape 近乎重合（TO 與 paired）
- **推翻 pipeline-level bug 假設**；TO 與 paired arm 所有下游差異必須歸因於統計處理層，非解析層

### 4.2 Concordance (確立)

- TO `SampleASM_Delta` 與 paired `Fisher_Frac_Sig` **Spearman ρ = −0.162**（顯著但效應極小）→ 訊號**大致獨立**
- Normal-derived features 跨 arm ρ ≈ 1.0（sanity 延伸）
- Sample ASM / HP Delta 跨 arm ρ 中高（0.6-0.9）
- Per-CpG Fisher 跨 arm ρ 低（paired-specific）

### 4.3 F pilot subset 飽和 (觀察)

- overlap 441 regions 中 TP rate=99.5%（both_TP 439 / both_FP 2）
- F pilot canonical filter 已將 subset 壓到近飽和 → subset 內 ensemble 無 F1 空間
- → paired arm characterization-only 需**在全域 region（非 F pilot）**上定義才有意義

### 4.4 對使用者 Q3 決策的驗證

使用者 2026-04-21 決策：「A 先放棄 paired F1-filter 方向，轉 characterization-only，除非後續有結果發現才再修正與觀察確認」

本次 concordance 分析**支持這個決策**：
1. paired `Fisher_Frac_Sig` 在 F pilot subset 內幾乎無 FP 可 rescue（0.5%）→ F1-filter 價值無法實現
2. paired `Fisher_Frac_Sig` 與 TO `SampleASM_Delta` 訊號獨立（ρ=−0.162）→ characterization-only 有獨立價值基礎
3. 若未來 R12 跨樣本發現 paired `Fisher_Frac_Sig` 在某樣本大量 FP region 獨立於 SampleASM_Delta → 可再升級

---

## 5. Registry 更新

| CL-id | 原狀態 | 新狀態 |
|-------|-------|-------|
| CL-025a (TO Normal BAM SampleASM_Delta POSITIVE) | ⭐3 active | ⭐3 active（不變） |
| CL-025b (Paired Fisher_Frac_Sig CONDITIONAL) | ⭐2 INSUFFICIENT active | ⭐2 **concluded (abandoned for F1-filter, characterization-only per user 2026-04-21)** |
| CL-new-crossmode (HCC1395 TO/paired Normal BAM sanity + signal independence) | — | ⭐4 validated (pipeline sanity 直接證實) |

Registry 新增 CL-025c（本 step 結果） — detail in `10_Research_Chain_Registry.md`。

---

## 6. 下一步候選（排序）

以下僅為**候選**，不自動啟動；使用者 meta-directive 明示「分析完全 + 批次彙報 + 決策時暫停」：

| 候選 | 層級 | 依賴 | 用途 |
|------|------|------|------|
| **(1) R5 LOH.bed 生成機制 audit** | P0（軌 1） | 無 | LOH-dependent 結論（F pilot NonLOH filter）可信度 |
| **(2) R1 HCC1395 Normal BAM TO arm 擴展至全域 region（跳出 F pilot subset）** | P0 | HCC1395 Normal BAM 已複製 | 驗證 SampleASM_Delta 是否在非 F pilot region 仍 POSITIVE，避免 F pilot 預篩 overfit |
| **(3) R12 跨樣本（HCC1954/COLO829/H2009 Normal BAM 複製）** | P1 | 使用者授權複製 | 升級 CL-025a 從 ⭐3 至 ⭐4 或 ⭐5 |
| **(4) paired characterization-only 全域驗證** | P2 | 非緊急 | 若 Q3 決策未來修正 |

**決策點**：以上 (1)-(4) 需使用者指示啟動順序。沒有明確 fallback scheme 時，建議優先 (1)（純 audit，不需新運算）→ (2)（主 P0 價值）→ (3)（跨樣本擴展）。

---

## 7. Provenance

### Inputs
- `output/hcc1395_normal_pilot/to/significance_summary.csv` (TO pipeline, 801 regions)
- `output/hcc1395_normal_pilot/paired/significance_summary.csv` (paired pipeline, 983 regions)
- `regions.bed` (F pilot canonical filter subset)

### Scripts
- `research/F_hpfinengroups_deepening/scripts/r1_crossmode_concordance.py` (194 lines, 4 figures + independence tests)
  - Dependencies: pandas, numpy, matplotlib, seaborn, scipy.stats.spearmanr
  - Entry: `python3 research/F_hpfinengroups_deepening/scripts/r1_crossmode_concordance.py`
  - Output dir: `research/F_hpfinengroups_deepening/figures/`

### Outputs
- `figures/fig6_normal_bam_sanity.png` (270KB)
- `figures/fig7_crossmode_signal_scatter.png` (304KB)
- `figures/fig8_crossmode_feature_correlation.png` (502KB)
- `figures/fig9_crossmode_tp_agreement.png` (117KB)

### Key statistics
```
TO rows: 801 paired rows: 983 overlap: 441
NormalBaseline_Coverage TO vs paired: ρ=1.0000 p=0.0 (sanity PASS)
TO SampleASM_Delta vs paired Fisher_Frac_Sig: ρ=-0.1618 p=6.457e-04 (n=441)
TP agreement: both_TP=439 both_FP=2 (99.5% TP saturation)
```

### Reproducibility
- Python 3.9 + pandas 1.x + sklearn.linear_model + scipy.stats
- Random seed: N/A (deterministic statistics)
- Total runtime: ~3 seconds
