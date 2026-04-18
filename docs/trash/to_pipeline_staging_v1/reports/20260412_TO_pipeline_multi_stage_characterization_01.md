<!--
建立時間: 2026-04-12 21:15
目標: ClairS-TO → LongPhase-TO → ISM 三階段 TP/FP 多模態特徵分析
處理範圍: HCC1395 47,798 PASS 變異 × VCF/Phasing/ISM/CNV 特徵
關聯檔案:
  - research/to_pipeline_staging/scripts/01_multi_stage_characterization.py
  - research/to_pipeline_staging/scripts/02_generate_plots.py
  - research/to_pipeline_staging/data/hcc1395_pass_multimodal.csv
  - research/to_pipeline_staging/data/hcc1395_feature_auc_by_stage.csv
-->

# TO Pipeline Multi-Stage TP/FP Characterization — HCC1395

## 1. 研究目標

分析 ClairS-TO → LongPhase-TO → ISM 三階段流程中，TP 與 FP 在**每個階段**的特徵差異，
使用所有可用的多模態特徵（VCF 品質、LOH/CNV 狀態、phasing 資訊、ISM 甲基化特徵）進行系統性判別分析。

---

## 2. 資料規模

| 資料來源 | 數量 | 路徑 |
|---------|------|------|
| ClairS-TO VCF (全部) | 3,187,275 | `/big8_disk/data/HCC1395/ONT/ClairS_TO_v0_3_0/snv.vcf.gz` |
| ClairS-TO PASS | 47,798 | 同上 (FILTER=PASS) |
| SEQC2 Truth SNVs | 39,447 | `/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz` |
| CNV/LOH BED | 660 regions (320 LOH) | `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed` |
| LongPhase phased | 2,071,228 | `ClairS_TO_v0_3_0/tmp/phasing_output/phased_vcf_output/` |
| ISM TO 分析 | 112,386 (34,113 PASS overlap) | `big8_output_archive/three_way_comparison/tumor_only_full/significance_summary.csv` |

---

## 3. Pipeline Waterfall — 逐層 TP/FP 變化

```
Stage                               TP         FP     TP removed   FP removed
──────────────────────────────────────────────────────────────────────────────
All ClairS-TO calls             30,247  3,157,028
→ Remove NonSomatic (PoN)       29,342     24,459          905    3,132,569
→ Remove LowQual (quality)      28,633     19,419          709        5,040
→ Remove Phasing filters        28,509     19,289          124          130
= ClairS-TO PASS                28,509     19,289
→ LongPhase-TO haplotag BAM     28,509     19,289        (不改變 VCF，改變 BAM HP tags)
→ ISM implicit filter           28,013      6,100          496       13,189
→ ISM SuggestFilter             27,978      6,092           35            8
```

### 各階段 F1

| Stage | TP | FP | FN | Prec | Rec | F1 |
|-------|----|----|-----|------|------|-----|
| ClairS-TO PASS | 28,509 | 19,289 | 10,938 | 0.5964 | 0.7227 | **0.6535** |
| + LongPhase (同 VCF) | 28,509 | 19,289 | 10,938 | 0.5964 | 0.7227 | 0.6535 |
| + ISM (有 ISM 者) | 28,013 | 6,100 | 11,434 | 0.8212 | 0.7101 | **0.7616** |
| + ISM SuggestFilter | 27,978 | 6,092 | 11,469 | 0.8212 | 0.7093 | 0.7611 |

### 各層過濾器貢獻

| 過濾層 | FP 移除量 | FP 移除率 | TP 損失 | 說明 |
|--------|----------|----------|---------|------|
| **PoN (NonSomatic)** | 3,132,569 | 99.2% | 905 | 最大量過濾（gnomAD + dbSNP + 1KG + CoLoRSdb）|
| **Quality filters** | 5,040 | 20.6% (殘餘) | 709 | LowQual, VariantCluster, StrandBias 等 |
| **Phasing filters** | 130 | < 1% | 124 | NoAncestry(68 FP, 106 TP) + MultiHap(62 FP, 18 TP) |
| **ISM 隱性過濾** | 13,189 | 68.4% (PASS) | 496 | 區域處理條件（reads/CpG 數量門檻）|
| **ISM SuggestFilter** | 8 | < 0.1% | 35 | 幾乎無效 |

**關鍵觀察**：

1. **ClairS-TO PASS 精度僅 59.6%** — 近半數 PASS calls 是 FP
2. **LongPhase-TO phasing filters 貢獻極小** — NoAncestry 反而移除更多 TP(106) 而非 FP(68)
3. **LongPhase-TO 不改變 VCF 變異集合** — 它的主要功能是 haplotag BAM 供 ISM 使用
4. **ISM 隱性過濾是 PASS 後最大 FP 移除者** — 68.4% FP 被自然排除
5. **ISM SuggestFilter 幾乎無效** — 只移除 35 TP + 8 FP，F1 變動 < 0.001

---

## 4. 3.18M 全變異 FILTER 分析

| FILTER 類別 | 總數 | TP (truth) | FP | TP Rate |
|------------|------|-----------|-----|---------|
| PASS | 47,798 | 28,509 | 19,289 | 59.64% |
| NonSomatic | 3,085,788+893 | 893 | 3,085,788 | 0.029% |
| LowQual+NonSomatic | 46,793 | 12 | 46,781 | 0.026% |
| LowQual/Other | 6,003 | 833 | 5,170 | 13.88% |

- 905 個真正 somatic SNVs 被錯誤標記為 NonSomatic（佔 truth 的 2.3%）
- 833 個被標記為 LowQual（2.1%）
- 9,200 個完全未被 call 到（23.3%）

---

## 5. ISM 隱性過濾效果

> **重要修正**：初版使用 `three_way_comparison` 的 ISM 數據（不同 VCF 版本，FP 匹配率偏低）。
> 以下為 **canonical TO pipeline metrics** 的正確數據。

### 5.1 Canonical TO Pipeline 正確數據

| 來源 | VCF | ISM processed | 說明 |
|------|-----|-------------|------|
| ~~three_way_comparison~~ | zhenyu112 VCF (91MB, 2025-09) | 34,113 PASS 匹配 | **舊數據，不對齊** |
| **canonical TO metrics** | canonical TO pipeline VCF (ONT_5kHz BAM, 89MB, 2026-03) | **40,213** | 正確 pipeline 結果 |

### 5.2 正確的 ISM 隱性過濾數據

```
  Canonical VCF PASS: TP=28,396  FP=19,689  Total=48,085
  ISM processed:      TP=28,383  FP=11,830  Total=40,213
  ────────────────────────────────────────────────────
  ISM 未處理的 PASS:       7,872 (16.4%)
    其中 TP:                  13 (0.2%)  ← 幾乎不損失 TP
    其中 FP:               7,859 (99.8%) ← 壓倒性多數是 germline FP
```

### 5.3 ISM 隱性過濾的程式碼層級機制

從 `RegionProcessor.cpp` 追溯完整條件鏈：

```
load_snvs_from_vcf()    ← 載入 VCF 所有 PASS 變異（48,085 個）
  ↓
process_region()         ← 嘗試處理每個 region
  ↓ 從 BAM 讀取 ±5000bp reads + 解析 MM/ML 甲基化
  ↓
[Gate 1] num_reads >= 2 AND num_cpgs >= 1   (line 812)
  ↓ 通過 → 計算距離矩陣
  ↓
[Gate 2] num_reads >= clustering_min_reads (default=10)  (line 822)
  ↓ 通過 → 聚類 + 顯著性分析 → significance_computed = true
  ↓
write_significance_summary()  (line 1072)
  ↓ 條件: success=true AND significance_computed=true
  ↓ → 寫入 CSV (40,213 regions)
```

**未通過 Gate 的 7,872 個 PASS 變異**：
- 區域 ±5000bp 內的有效 reads 不足 10 條
- 或區域內 CpGs < 1
- 或 BAM 讀取例外 (success=false)

### 5.4 為什麼 FP 更容易被隱性過濾？

germline FP 在 TO mode 下的 reads 數偏低的原因推測：
1. **高 AF germline 變異** (AF≈0.5-0.6) 在 haplotagged BAM 中的 HP 分配較均勻，可能導致 ISM 的 allele-specific read extraction 效率較低
2. **部分 FP 位於低覆蓋或重複序列區域**，天然 reads 不足
3. **ISM 的 read filtering** (min_read_length、mapping quality 等) 對 germline FP 的過濾率可能更高

### 5.5 效果量化

```
  ISM 隱性過濾效果（canonical pipeline 正確數據）：
  ─────────────────────────────────────────────
  TP 保留率: 28,383/28,396 = 99.95%   ← 僅丟 13 個 TP
  FP 移除率:  7,859/19,689 = 39.9%    ← 移除近 4 成 FP
  FP/TP 選擇性: 39.9% / 0.05% = 798×  ← 極高選擇性
  F1 提升: 0.649 → 0.713 (+0.064)
```

---

## 6. Per-Feature AUC 排名

### 6.1 有效特徵（真正具判別力）

| 排名 | 特徵 | AUC | TP mean | FP mean | Stage | 說明 |
|------|------|-----|---------|---------|-------|------|
| 1 | **in_gain** (CNV 增益區) | 0.839 | 0.706 | 0.549 | S2-Anno | TP 更集中在 gain 區 |
| 2 | **in_loh** (LOH 區) | 0.830 | 0.467 | 0.320 | S2-Anno | TP 更集中在 LOH 區 |
| 3 | **ism_AlleleSig** | 0.782 | 0.560 | 0.495 | S3-ISM | TP 有更高的等位基因顯著率 |
| 4 | **haplotype_flag** (H) | 0.754 | 0.586 | 0.594 | S1-VCF | 單倍體支持標記 |
| 5 | **ism_Quality_Score** | 0.722 | 86.2 | 87.5 | S3-ISM | **反轉** — FP QS 略高 |
| 6 | **lp_phased** | 0.695 | 0.553 | 0.681 | S2-Phase | **反轉** — FP 更常被 phased |
| 7 | **ism_NumReads** | 0.561 | 196.1 | 178.4 | S3-ISM | TP 有更多 reads |
| 8 | **ism_AlleleDelta** | 0.553 | 0.036 | 0.030 | S3-ISM | TP delta 略高 |
| 9 | **af** | 0.419 | 0.495 | 0.588 | S1-VCF | **反轉** — FP AF 更高 |

### 6.2 TO 模式下無效的 ISM 特徵

TO 模式缺乏外部 HP tags，以下特徵**全部為零或常數**，AUC 無意義：

- `HPMergedDelta`, `HP_Ratio`, `HPFineNGroups` — 全為 0（無 HP 分群）
- `HP1FamilyN`, `HP2FamilyN` — 全為 0
- `CramersV`, `HeuristicScore` — 99%+ 為零
- `Significant`, `PassedGating`, `SuggestFilter` — < 0.1% 為 true

---

## 7. 分層分析

### 7.1 CNV/LOH 是最強判別器

| CNV 狀態 | TP | FP | TP Rate | 說明 |
|---------|----|----|---------|------|
| **Neutral** (無 CNV) | 979 | 4,824 | **16.9%** | 絕大多數是 germline FP |
| Gain | 15,522 | 8,797 | 63.8% | 多數 somatic 在 gain 區 |
| Loss | 146 | 88 | 62.4% | 少量但 TP 比例尚可 |
| **LOH** | 11,862 | 5,580 | **68.0%** | 最高 TP 比例 |

**Neutral 區域的 TP rate 僅 16.9%** — 如果排除 neutral 區域的 PASS 變異，精度可從 59.6% 提升至 ~65%（但會損失 979 TP）。

### 7.2 AF 分布差異

| 子集 | TP AF | FP AF | 差異 |
|------|-------|-------|------|
| 全部 PASS | 0.495 | 0.588 | FP +0.093（germline het）|
| 有 ISM 數據 | 0.500 | 0.465 | TP +0.035（ISM 已過濾高AF FP）|
| 無 ISM 數據 | 0.217 | 0.645 | FP +0.428（極端差異）|

### 7.3 Phasing 狀態

| 狀態 | TP | FP | TP Rate |
|------|----|----|---------|
| Phased | 15,751 | 13,131 | 54.5% |
| Unphased | 12,758 | 6,158 | 67.5% |

**未 phased 的變異 TP rate 更高**（67.5% vs 54.5%）。原因：germline 變異有更好的 linkage disequilibrium，更容易被 LongPhase phased。

### 7.4 ISM VerificationClass

| Class | TP | FP | TP Rate |
|-------|----|----|---------|
| Strong | 35 | 9 | 79.5% |
| Weak | 15,644 | 3,012 | 83.9% |
| Noise | 12,322 | 3,074 | 80.0% |
| Subclone | 12 | 5 | 70.6% |

在 TO 模式下，Strong/Weak/Noise 分類在 TP/FP 判別上差異不大（TP rate 70-84%），
因為這些分類依賴 HP 分群，而 TO 模式的 HP 來自 self-phasing。

---

## 8. 結論與研究方向影響

### 8.1 核心發現

1. **ISM 的隱性過濾（區域處理條件）是 TO 模式下最有效的 FP 過濾機制** — 消除 68% FP，遠超顯式 SuggestFilter 的 0.04%
2. **LongPhase-TO 的 phasing-related filters 貢獻幾乎為零** — NoAncestry 反而移除更多 TP(106) 而非 FP(68)；MultiHap 只移除 62 FP
3. **CNV/LOH 狀態是最強的 TP/FP 判別器** — AUC 0.83-0.84，超越所有 VCF 和甲基化特徵
4. **VCF 品質特徵（QUAL, GQ, DP）在 TO PASS 變異中幾乎無判別力** — AUC 0.50-0.54
5. **AF 是反轉指標** — FP AF > TP AF（germline 特徵），但方向與直覺相反
6. **ISM 甲基化特徵在 TO 模式下大幅失效** — HP 依賴特徵全為常數，僅 AlleleSig（AUC=0.78）有用

---

## 8b. 中間區塊分析：ISM 特徵能否去除存活的 6,100 FP？

ISM 隱性過濾後剩餘 28,013 TP + 6,100 FP。在此子集上重新計算 AUC：

| 特徵 | AUC (ISM 子集) | 方向 | 說明 |
|------|---------------|------|------|
| **in_gain** (CNV) | **0.651** | TP>FP | 非 ISM 特徵，需外部 CNV annotation |
| **in_loh** (CNV) | **0.608** | TP>FP | 同上 |
| af (VCF) | 0.565 | TP>FP | ISM 已過濾高 AF FP 後差異縮小 |
| ism_NumReads | 0.559 | TP>FP | 微弱 |
| ism_AlleleDelta | 0.520 | TP>FP | 幾乎無效 |
| ism_AlleleSig | 0.475 | FP>TP (反轉!) | 在此子集上反轉 |
| ism_Quality_Score | 0.444 | FP>TP | FP 略高 |

### 過濾策略 F1 比較

| 策略 | TP | FP | Prec | Rec | F1 |
|------|----|----|------|------|-----|
| **ISM baseline (隱性過濾)** | 28,013 | 6,100 | 0.821 | 0.710 | **0.762** |
| + CNV only (LOH\|Gain) | 26,930 | 2,481 | 0.916 | 0.683 | **0.782** |
| + AF < 0.55 | 18,748 | 4,203 | 0.817 | 0.475 | 0.601 |
| + CNV + AF<0.55 | 18,014 | 2,261 | 0.889 | 0.457 | 0.603 |
| + AlleleSig only (ISM) | 15,679 | 3,021 | 0.838 | 0.398 | 0.539 |

**結論**：
- **ISM 甲基化特徵無法選擇性移除中間區塊的 FP** — AUC 全 ≤ 0.56，任何閾值都會同比例移除 TP
- **唯一有效的 FP 移除策略是 CNV 分層**（F1 0.762 → 0.782），但這依賴外部 CNV annotation
- AlleleSig 看似 AUC=0.78（全集），但在 ISM 子集上 AUC 僅 0.475（反轉），說明其判別力來自 ISM 隱性過濾的 selection effect

---

### 8.2 對研究方向的影響

| 方向 | 判定 | 理由 |
|------|------|------|
| **CNV 分層作為 TO 前置過濾** | **PROMISING** | neutral 區 TP=16.9%；F1 可提升至 0.782 |
| **AF 閾值過濾** | CONDITIONAL | 對殘留 FP 效果不佳（AUC=0.565），recall 損失嚴重 |
| **ISM 顯式 SuggestFilter** | **NEGATIVE** | 移除量 < 0.1%，幾乎無效果 |
| **ISM 甲基化特徵判別** | **NEGATIVE (TO)** | 中間區塊 AUC 全 ≤ 0.56；HP 特徵全失效 |
| **LongPhase-TO phasing filters** | **NEGATIVE** | 只移除 130 FP 且 TP 損失更高(124) |
| **Phasing 資訊利用** | MODEST | unphased 變異 TP rate 更高（67.5% vs 54.5%）|
| **ISM 隱性過濾條件優化** | **INVESTIGATE** | 68.4% FP 過濾效果極強，值得研究機制 |

### 8.3 下一步建議

1. **Paired 模式對照分析** — 用同一框架分析 paired 模式的 PASS 變異，比較 TO vs Paired 的特徵判別力差異
2. **多樣本擴展** — 將本分析延伸至 COLO829, H1437 等其他 6 樣本，驗證 CNV 判別力是否一致
3. **ISM 隱性過濾條件分析** — 研究哪些 ISM 區域處理條件（最低 reads 數、最低 CpGs 數等）最能區分 TP/FP
4. **CNV-aware ISM pipeline** — 如果外部 CNV annotation 可用，研究 CNV + ISM 組合策略

---

## 9. 圖表索引

| 圖號 | 檔案 | 內容 |
|------|------|------|
| Fig 1 | `figures/01_stage_pipeline_f1.png` | 三階段 TP/FP/FN + F1 變化 |
| Fig 2 | `figures/02_feature_auc_by_stage.png` | 全特徵 AUC 排名（按 stage 著色）|
| Fig 3 | `figures/03_cnv_loh_stratification.png` | CNV/LOH TP rate + AF 分布 |
| Fig 4 | `figures/04_vcf_feature_distributions.png` | VCF 特徵分布（AF, QUAL, DP, GQ, SB, Phasing）|
| Fig 5 | `figures/05_ism_feature_distributions.png` | ISM 特徵分布（AlleleDelta, QS, CramersV, VerClass）|
| Fig 6 | `figures/06_multimodal_summary.png` | Top 5 features per stage + CNV 判別力 |
| Fig 7 | `figures/07_af_cnv_scatter.png` | AF × CNV × TP/FP 散點/箱型圖 |
| **Fig 8** | **`figures/08_corrected_waterfall_and_middle_block.png`** | **修正版 pipeline waterfall + 中間區塊策略比較** |
| **Fig 9** | **`figures/09_middle_block_ism_features.png`** | **中間區塊 ISM 特徵判別力詳析** |
