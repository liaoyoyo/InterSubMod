<!--
建立時間: 2026-04-03 19:00
目標: 驗證 PON-only phased VCF 經 haplotag 後對 ISM HP_Ratio 和 F1 的影響
處理範圍: Haplotag (baseline vs PON-only) → ISM 重跑 → HP_Ratio/F1/LOH 三方比較 (Paired vs TO-Baseline vs TO-PON-Only)
關聯檔案:
  - research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md (Step 1: Phasing 驗證)
  - /big7_disk/liaoyoyo2001/longphase-to-mod/output/ (全部實驗結果)
  - /big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/20260119_all-with-w5000_2/ (Paired ISM baseline)
-->

# PON-Only Haplotag + ISM 驗證報告

## 摘要

本報告為 PON-only phasing 驗證的第二步：將 phased VCF 經 haplotag 產生 HP-tagged BAM，再由 ISM 分析 HP_Ratio、VerificationClass 和 F1。比較三種模式：Paired (LongPhase-S + normal)、TO-Baseline (原始 LongPhase-TO)、TO-PON-Only (修改後)。

### 核心結論

| 發現 | 判定 |
|------|------|
| PON-only phasing VCF 層面 | **正確** — self-phasing 消除、N50 翻倍 |
| LOH.bed | **不受影響** — Jaccard=1.0000 |
| **Haplotag HP tag 分配** | **有缺陷** — somatic 位點所有 reads 被統一標記為 HP:i:21 |
| ISM HP_Ratio | **惡化** — extreme rate 從 60% 升至 99.9% |
| ISM F1 (SuggestFilter) | **無顯著差異** — 0.0965 → 0.0979 |
| ISM-only LOH excess | **惡化** — 從 15.4% 升至 54.8% |

**根因**：LongPhase-TO haplotag 在處理 PON-only phased VCF 的 somatic 位點時（GT=`0|0`, GT2=`.|1`），將所有 overlapping reads 統一標記為 `HP:i:21`（somatic HP2），而非正確地根據各 read 的 germline variant 分配 HP1/HP2。

---

## 1. 實驗設計

### 流程

```
Step 1 (已完成): LongPhase-TO phase
  → baseline/tumor_phased.vcf  (原始算法)
  → pononly/tumor_phased.vcf   (PON-only phasing)

Step 2 (本報告): LongPhase-TO haplotag
  → baseline/tumor_tagged.bam  (baseline phased VCF → haplotag)
  → pononly/tumor_tagged.bam   (PON-only phased VCF → haplotag)

Step 3 (本報告): ISM 分析
  → 6 組 ISM 結果: 3 modes × (TP + FP)
  → 比較 HP_Ratio, F1, LOH excess
```

### 比較架構

| 模式 | Phased VCF 來源 | Haplotagged BAM | ISM 結果 |
|------|-----------------|-----------------|----------|
| Paired | LongPhase-S + normal BAM | `/big8_disk/.../HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam` | 既有 (20260119) |
| TO-Baseline | LongPhase-TO 原始算法 | `output/baseline/tumor_tagged.bam` | 本次新跑 |
| TO-PON-Only | LongPhase-TO PON-only | `output/pononly/tumor_tagged.bam` | 本次新跑 |

### ISM 共用參數

- VCF: 30,490 TP + 4,822 FP (Paired pileup mode TP/FP)
- Normal BAM: `HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam`
- Distance: BERNOULLI, Window: ±5000bp, Threads: 48

---

## 2. Haplotag 結果

| 指標 | TO-Baseline | TO-PON-Only | 變化 |
|------|-------------|-------------|------|
| 執行時間 | 39m12s | **29m41s** | **1.32x 快** |
| Total alignments | 40,859,727 | 40,859,727 | 相同 |
| Tagged alignments | 19,571,246 | **20,128,638** | **+2.9%** |
| Tag rate | 47.9% | **49.3%** | +1.4pp |

PON-only haplotag 更快、tagged 更多 reads — 因為 phasing 品質更高 (N50 翻倍)。

---

## 3. HP_Ratio 分佈

### 3.1 總覽

| Dataset | N | Median | Mean | Extreme (>0.9/<0.1) | Extreme% |
|---------|---|--------|------|---------------------|----------|
| Paired-TP | 30,476 | **0.5000** | 0.5051 | 14,272 | 46.8% |
| TO-Baseline-TP | 30,476 | 0.8358 | 0.6746 | 18,285 | 60.0% |
| TO-PONOnly-TP | 30,476 | 0.0000 | 0.4432 | **30,452** | **99.9%** |
| Paired-FP | 4,822 | 0.5405 | 0.5209 | 2,885 | 59.8% |
| TO-Baseline-FP | 4,822 | 0.9555 | 0.6895 | 3,267 | 67.8% |
| TO-PONOnly-FP | 4,822 | 0.9999 | 0.5559 | **4,755** | **98.6%** |

### 3.2 HP Family 分佈（TP variants）

| 指標 | Paired | TO-Baseline | TO-PON-Only |
|------|--------|-------------|-------------|
| HP1 total | 975,993 (44.9%) | 1,266,344 (**58.2%**) | 780,055 (35.8%) |
| HP2 total | 970,993 (44.6%) | 735,008 (33.8%) | 1,248,256 (**57.4%**) |
| HP3 total | 59,443 (2.7%) | 10,008 (0.5%) | 5,804 (0.3%) |
| HP0 total | 169,684 (7.8%) | 164,753 (7.6%) | 141,998 (**6.5%**) |
| HP1 per-variant median | **30** | 39 | **0** |
| HP2 per-variant median | **30** | 12 | **41** |
| Balanced (0.3-0.7) | **36.9%** | 26.8% | **0.04%** |
| HP_Ratio=0.0 | 14.9% | 6.4% | **55.1%** |
| HP_Ratio=1.0 | 16.0% | 35.1% | **43.6%** |

### 3.3 解讀

- **Paired**：HP1 ≈ HP2 平衡（germline phasing 正確分割兩個 haplotype）
- **TO-Baseline**：HP1 dominant（58.2%）— self-phasing 將 somatic reads 偏向 HP1（17.3:1 bias 的下游效果）
- **TO-PON-Only**：**二元分佈** — 每個 variant 的 reads 要嘛全 HP1 要嘛全 HP2，幾乎無混合。HP2 整體 dominant（57.4%）

---

## 4. 根因分析：Haplotag HP Tag 分配缺陷

### 4.1 單位點驗證（chr1:117740955，非 LOH 區域）

**Phased VCF genotype：**

| | GT | GT2 | PS |
|-|----|----|-----|
| Baseline | `1\|0` | `./.` | 117272953 |
| PON-only | `0\|0` | `.\|1` | 117272953 |

**BAM HP tag 分佈（130 reads at position）：**

| HP Tag | 意義 | Paired | TO-Baseline | TO-PON-Only |
|--------|------|--------|-------------|-------------|
| HP:Z:1 / HP:i:1 | Germline HP1 | 1 | 25 | **0** |
| HP:Z:1-1 / HP:i:11 | Somatic HP1 | 57 | 33 | **0** |
| HP:Z:2 / HP:i:2 | Germline HP2 | 72 | 72 | **0** |
| HP:i:21 | Somatic HP2 | 0 | 0 | **129** |
| None | Untagged | 0 | 0 | 1 |
| → HP_Ratio | | **0.45** | **0.45** | **0.00** |

### 4.2 機制

1. **Baseline** phased VCF：variant GT=`1|0` → haplotag 視為 phased heterozygous → ALT reads → HP:i:11 (somatic HP1)，REF reads → HP:i:2 (germline HP2)，部分 reads → HP:i:1 (germline HP1)
2. **PON-only** phased VCF：variant GT=`0|0`, GT2=`.|1` → haplotag 視為 homozygous ref（primary GT），somatic info 在 GT2 → haplotag 將**所有** overlapping reads 標記為 `HP:i:21`（somatic HP2）

### 4.3 問題根因

LongPhase-TO haplotag 在處理 GT=`0|0` + GT2=`.|1` 的位點時，不根據各 read 的 germline variant 分配 HP1/HP2，而是統一將所有 reads 標記為 somatic haplotype（HP:i:21）。

**預期行為**：
- 攜帶 ALT allele 的 reads → HP:i:21（somatic HP2）
- 攜帶 REF allele 的 reads → HP:i:1 或 HP:i:2（根據 surrounding germline variants 決定）

**實際行為**：
- **所有** reads → HP:i:21（somatic HP2）

### 4.4 6,485 個非 LOH 平衡位點驗證

選取 Paired 模式下 HP_Ratio 0.4-0.6（平衡）且不在 LOH.bed 中的 6,485 個 TP variants：
- **TO-PON-Only HP_Ratio**: median=0.0000, mean=0.0303
- **Extreme rate**: 6,483/6,485 (**100.0%**)

確認此非單一位點問題，而是 **系統性缺陷**。

---

## 5. ISM 指標比較

### 5.1 VerificationClass 分佈（TP）

| Class | Paired | TO-Baseline | TO-PON-Only |
|-------|--------|-------------|-------------|
| Strong | 7,279 (23.9%) | 7,139 (23.4%) | 6,353 (**20.8%**) |
| Weak | 11,335 (37.2%) | 10,868 (35.7%) | 10,045 (**33.0%**) |
| Subclone | 1,428 (4.7%) | 1,462 (4.8%) | 1,442 (4.7%) |
| Noise | 10,434 (34.2%) | 11,007 (36.1%) | 12,636 (**41.5%**) |

TO-PON-Only 的 Noise 佔比增加（+5.4pp），Strong 和 Weak 減少。

### 5.2 Quality Score

| Dataset | QS Median | QS Mean |
|---------|-----------|---------|
| Paired-TP | **85.0** | **78.8** |
| TO-Baseline-TP | 75.0 | 75.4 |
| TO-PON-Only-TP | 75.0 | **65.1** |

### 5.3 ISM Filter F1（SuggestFilter）

| 模式 | FP caught | TP lost | Precision | Recall | **F1** |
|------|-----------|---------|-----------|--------|--------|
| Paired | 235 | 113 | 0.6753 | 0.0487 | **0.0909** |
| TO-Baseline | 250 | 112 | 0.6906 | 0.0518 | **0.0965** |
| TO-PON-Only | 252 | 73 | **0.7754** | 0.0523 | **0.0979** |

F1 差異極小（<0.01），三者基本一致。值得注意的是 TO-PON-Only 的 TP lost 最低（73 vs 112/113），precision 最高（0.78）。

### 5.4 ISM-only LOH Excess

| 模式 | ISM LOH=true | In LOH.bed | NOT in bed (excess) | Excess rate | Concordance |
|------|-------------|------------|---------------------|-------------|-------------|
| Paired | 14,272 | 13,086 | 1,186 | **3.9%** | **91.7%** |
| TO-Baseline | 18,285 | 13,606 | 4,679 | 15.4% | 74.4% |
| TO-PON-Only | 30,452 | 13,736 | 16,716 | **54.8%** | 45.1% |

Excess 從 Paired 3.9% → TO-Baseline 15.4% → TO-PON-Only **54.8%**。但這並非 phasing 惡化，而是 haplotag HP tag 分配缺陷導致 ISM 錯判。

---

## 6. LOH.bed 結果確認

| 指標 | Baseline | PON-Only |
|------|----------|----------|
| Region 數量 | 1,094 | 1,094 |
| 總覆蓋量 | 1,632,234,970 bp | 1,632,234,970 bp |
| Jaccard (相互) | — | **1.0000** |

**LOH.bed 完全不受 PON-only phasing 或 haplotag 影響。**

---

## 7. 結論

### 已確認

1. **PON-only phasing 在 VCF 層面完全正確**：self-phasing 消除（17.3:1 → 0）、N50 翻倍、phased rate +23.6pp
2. **LOH.bed 不受任何影響**：1,094 regions 完全一致
3. **Haplotag 有 HP tag 分配缺陷**：somatic 位點所有 reads 統一標記為 somatic haplotype，而非正確分割
4. **ISM F1 基本不變**：SuggestFilter F1 差異 <0.01
5. **ISM HP_Ratio 和 LOH excess 因 haplotag 缺陷而惡化**：非 phasing 本身問題

### Paired vs TO 的根本差異

| 維度 | Paired | TO |
|------|--------|-----|
| Germline phasing 來源 | Normal BAM (高品質) | PON filtering (adequate) |
| Somatic variant 在 phasing 中 | 不參與 | Baseline: 污染; PON-only: 排除 |
| Haplotag 行為 | Germline HP1/HP2 + somatic HP1-1 | PON-only: **全部 somatic tag** |
| HP_Ratio 語義 | 兩個 haplotype 的 read 比例 | 無法正確反映 |

### 修復方向

| 方案 | 複雜度 | 效果 |
|------|--------|------|
| A: 修改 haplotag 代碼 | 中 | 使 somatic 位點 reads 按 germline variants 正確分割 HP |
| B: ISM 端忽略 somatic HP tags | 低 | ISM 只使用 HP:i:1/HP:i:2，忽略 HP:i:11/HP:i:21/HP:i:33 |
| C: Post-process BAM | 中 | 重新分配 somatic 位點的 HP tags |
| D: ISM 重新定義 TO 模式 HP_Ratio | 中 | 不依賴 HP tag，改用 allele-based assignment |

**建議優先方案 B**：在 ISM ReadParser 中，遇到 TO 模式 haplotagged BAM 時，只計算 HP:i:1 和 HP:i:2 reads（忽略 somatic tags），最為簡單且不需修改 LongPhase-TO。

---

## 8. 檔案清單

| 檔案 | 說明 |
|------|------|
| `output/baseline/tumor_tagged.bam` | Baseline haplotagged BAM (260G) |
| `output/baseline/tumor_tagged.bam.bai` | Index |
| `output/pononly/tumor_tagged.bam` | PON-only haplotagged BAM |
| `output/pononly/tumor_tagged.bam.bai` | Index |
| `output/ism_baseline_tp/significance_summary.csv` | TO-Baseline ISM TP 結果 |
| `output/ism_baseline_fp/significance_summary.csv` | TO-Baseline ISM FP 結果 |
| `output/ism_pononly_tp/significance_summary.csv` | TO-PON-Only ISM TP 結果 |
| `output/ism_pononly_fp/significance_summary.csv` | TO-PON-Only ISM FP 結果 |
| `compare_haplotag_ism.py` | 比較分析腳本 |
| `run_haplotag_and_ism.sh` | 執行腳本 |

所有結果目錄：`/big7_disk/liaoyoyo2001/longphase-to-mod/output/`
