# LOH 區域 Intermediate AF 與亞克隆偵測研究計劃

> **狀態**：completed
> **建立日期**：2026-04-14
> **專案目錄**：`research/loh_subclone_af/`
> **結論**：POSITIVE — 雙重證據鏈確認

## 背景與動機

觀察到 LOH 區域中部分 variants 的 allele frequency (AF) 偏離 clonal LOH 的期望值。在 purity = 1.0 的純腫瘤 cell line 中，deletion LOH (CN=1) 期望 AF = 0 或 1.0，cnLOH (CN=2) 期望 AF = 0、0.5 或 1.0。Intermediate AF (0.1-0.4 / 0.6-0.9) 在這些情境下不應出現，除非存在 subclonal LOH 事件。

文獻支持：TITAN、Battenberg、FACETS 均以 intermediate BAF 偵測 subclonal CNA。CAMDAC 證明 clonal LOH → ASM 消失、subclonal LOH → ASM 部分保留。但目前無工具結合 AF + methylation 在長讀長上偵測 subclonal LOH。

## 假說

### H1：Deletion LOH (CN=1) 的 intermediate AF = subclonal LOH

**前提條件**：Cell line purity = 1.0；Coverage_Multiple < 0.75 作為 CN≈1 proxy

**已知 Confound**：NumReads（HPFineNGroups 與 coverage 正相關）

**驗證標準**：
- **Positive**：Intermediate NGroups > Extreme NGroups (Mann-Whitney p<0.05, ≥5/7 samples)
- **Negative**：兩組 NGroups 無差異

### H3：多位點連鎖 + 甲基化 = 更強的 subclone 證據

**前提條件**：LOH.bed segments 可用

**驗證標準**：
- **Positive**：Segment AF-SD vs NGroups positive Spearman ρ (p<0.05, ≥5/7)
- **Negative**：無相關或負相關

### H4：HPFineNGroups 反映亞克隆複雜度

**驗證標準**：
- **Positive**：NR-bin controlled 後效應增強或不變
- **Negative**：效應消失

## 方法

### 數據來源

| 數據集 | 路徑 | 描述 |
|--------|------|------|
| Master dataset | `all_region_rows.tsv.gz` | 748K rows, 7 samples × 2 modes |
| LOH.bed | 7 files (see step3 script) | LongPhase-TO generated LOH segments |

### 分析步驟

```
Step 1: LOH AF 分布 baseline → 驗證: 圖表+統計表產出, intermediate AF 比例可計算
Step 2: Intermediate AF × NGroups 交叉分析 → 驗證: 7/7 MW test p<0.05, NR-controlled
Step 3: LOH segment 空間分析 → 驗證: Spearman ρ>0 in ≥5/7 samples
```

## 已知風險

1. **NumReads confound**：控制方式 = NR-bin 分層分析
2. **Coverage_Multiple 非精確 CN**：僅為 proxy，r=0.831 vs truth CN
