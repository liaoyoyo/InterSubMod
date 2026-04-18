# LOH Subclone AF × Methylation: Paired Mode 延伸驗證

> **狀態**：executing
> **建立日期**：2026-04-15
> **專案目錄**：`research/loh_subclone_af_paired/`
> **前置研究**：`research/loh_subclone_af/`（TO mode, POSITIVE）

## 背景與動機

TO mode 研究已確認 LOH 區域 intermediate AF 與 HPFineNGroups 的正向關聯（7/7 p<10⁻³⁹）。Paired mode 使用 matched normal BAM 做 somatic calling，AF 和 TP/FP 分類更準確。需驗證此生物學信號在更準確的數據中是否同樣存在。

## 假說

| ID | 假說 | Positive 標準 | TO baseline |
|----|------|-------------|-------------|
| H1p | Paired LOH intermediate AF → 更高 NGroups | ≥5/7 MW p<0.05 | 7/7 p<10⁻³⁹ |
| H2p | Paired 效應量 ≥ TO | median \|r\| ≥ TO median | TO median \|r\|≈0.58 |
| H3p | Paired segment AF-SD ∝ NGroups | ≥5/7 ρ>0 | 6/7 ρ>0 |
| H4p | 跨模式一致性 | ≥5/7 效應方向一致 | — |

## 方法

### 數據來源

- Master dataset: `all_region_rows.tsv.gz`（748K rows, paired = 328,699 rows）
- LOH.bed: 7 files from TO pipeline（用於標註 paired variants 的 LOH 區域）

### 分析步驟

```
Step 1: Paired LOH AF 分布 baseline → 驗證: 圖表+統計表產出
Step 2: NGroups × AF 交叉分析 + NR 控制 → 驗證: ≥5/7 MW p<0.05
Step 3: Segment 空間分析 → 驗證: ≥5/7 ρ>0
Step 4: Paired vs TO 直接比較 → 驗證: 效應方向跨模式一致
```

## 已知風險

1. **Paired FP 極少**：部分樣本 FP<30，FP 分析 power 不足
2. **LOH.bed 來自 TO pipeline**：用 TO 的 LOH 區域標註 paired variants（LOH 為基因體性質，非 pipeline 性質）
3. **NGroups 非獨立**：paired/TO 分析同一 tumor BAM reads
