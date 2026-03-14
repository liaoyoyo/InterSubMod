<!--
建立時間: 2026-01-13 18:50
更新時間: 2026-03-12 00:40
目標: 甲基化顯著性與 F1 最佳化研究總覽
處理範圍: 探索甲基化特徵與 SNV TP/FP 分類的關係
關聯檔案:
  - docs/plans/2026/01/20260113_甲基化F1研究計劃_01.md
  - data/README.md
-->

# 甲基化顯著性與 F1 最佳化研究

## 研究動機

InterSubMod 透過甲基化異質性分析偵測體細胞變異 (SNV)，但目前：
- **Cramér's V AUC 僅 0.5194**：現有甲基化指標無法有效區分 TP/FP
- **TP 顯著率僅 6.11%**：93.9% 的真陽性位點沒有顯著的甲基化訊號
- **直接過濾會惡化 F1**：從 0.8155 降至 0.0899

本研究探索是否存在 **組合特徵** 或 **條件過濾策略** 可以改善分類效果。

---

## 基準數據

| 指標 | 明恩 SNV |
|:---|---:|
| TP | 30,490 |
| FP | 4,842 |
| FN | 8,960 |
| Precision | 0.8630 |
| Recall | 0.7729 |
| **F1-score** | **0.8155** |
| 理論上限 | 0.8736 |

### F1 優化的取捨關係

- **嚴格方向**（降 FP）：每移除 1 FP，最多可誤刪 0.69 TP → 比值 > 1.45 有利
- **寬鬆方向**（降 FN）：每救回 1 FN，最多可誤抓 1.449 FP → 比值 < 1.45 有利

**判斷公式**：`FP 移除數 / TP 移除數 > 1.45` → 策略有效

---

## 現行目錄邊界

```
methylation_f1_optimization/
├── README.md
├── 2026/01/*.md                       # 正式研究文檔
└── assets/2026/01/                    # 文檔直接引用的圖表資產
```

### 邊界說明

1. `docs/research/methylation_f1_optimization/`
   - 只保留研究文檔與最小必要圖表資產
2. `scripts/analysis/legacy/methylation_f1_optimization/`
   - 承接歷史研究分析腳本
3. `docs/archive/2026/03/20260312_docs_round2_pending_review_01/`
   - 承接原本混在 `docs/research/` 內的空骨架與舊工作區結構，保留待審

## 研究文件

1. [docs/research/methylation_f1_optimization/2026/01/20260115_phase1_report_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_f1_optimization/2026/01/20260115_phase1_report_01.md)
2. [docs/research/methylation_f1_optimization/2026/01/20260115_phase2_report_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_f1_optimization/2026/01/20260115_phase2_report_01.md)
3. [docs/research/methylation_f1_optimization/2026/01/20260115_phase3_verification_report_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_f1_optimization/2026/01/20260115_phase3_verification_report_01.md)
4. [docs/research/methylation_f1_optimization/2026/01/20260115_phase4_methylation_combination_report_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_f1_optimization/2026/01/20260115_phase4_methylation_combination_report_01.md)
5. [docs/research/methylation_f1_optimization/2026/01/20260115_final_conclusion_report_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_f1_optimization/2026/01/20260115_final_conclusion_report_01.md)

## 支援腳本

1. [scripts/analysis/legacy/methylation_f1_optimization/phase1_analysis.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/legacy/methylation_f1_optimization/phase1_analysis.py)
2. [scripts/analysis/legacy/methylation_f1_optimization/phase2_analysis.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/legacy/methylation_f1_optimization/phase2_analysis.py)
3. [scripts/analysis/legacy/methylation_f1_optimization/phase3_verification.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/legacy/methylation_f1_optimization/phase3_verification.py)
4. [scripts/analysis/legacy/methylation_f1_optimization/phase4_methylation_combination.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/legacy/methylation_f1_optimization/phase4_methylation_combination.py)

## 圖表資產

1. `docs/research/methylation_f1_optimization/assets/2026/01/phase1_plots/`
2. `docs/research/methylation_f1_optimization/assets/2026/01/phase2_plots/`
3. `docs/research/methylation_f1_optimization/assets/2026/01/phase3_plots/`
4. `docs/research/methylation_f1_optimization/assets/2026/01/phase4_plots/`
5. `docs/research/methylation_f1_optimization/assets/2026/01/longphase_s_qual_af_analysis.png`

---

## 可觀察特徵 (15 個面向)

### A. VCF 原始特徵

| 編號 | 特徵 | 說明 | 預期效果 |
|:---:|:---|:---|:---|
| F01 | QUAL | ClairS 品質分數 | TP 應較高 |
| F02 | AF | Variant Allele Frequency | 可能有特定範圍 |
| F03 | DP | Read Depth | 高深度更可靠 |
| F04 | H flag | Haplotype-specific variant | 與 phasing 相關 |
| F05 | Strand Bias | FAU/RAU 比值 | FP 可能有更強 bias |

### B. 甲基化分析特徵

| 編號 | 特徵 | 說明 | 預期效果 |
|:---:|:---|:---|:---|
| F06 | CramersV | 效應量 (0-1) | 現有 AUC=0.52 |
| F07 | LabelDelta | 標籤距離差異 | 已知有微弱效果 |
| F08 | GlobalP | 全域 p-value | 低 P 可能更可信 |
| F09 | VerificationClass | Strong/Subclone/Weak/Noise | Strong 應有更高 TP 率 |
| F10 | DominantLabel | hp/allele/none | 可能與變異類型相關 |

### C. 衍生與交互特徵

| 編號 | 特徵 | 說明 | 預期效果 |
|:---:|:---|:---|:---|
| F11 | NumReads | 甲基化分析深度 | 低深度不可靠 |
| F12 | NumCpGs | CpG 位點數量 | 少 CpG 增加噪音 |
| F13 | Regional Clustering | 50kb 窗口聚集度 | 高聚集可能是 FP |
| F14 | QUAL × CramersV | 品質與效應量交互 | 尋找最佳組合閾值 |
| F15 | AF × LabelDelta | VAF 與標籤差異交互 | 可能發現隱藏規律 |

---

## 判斷標準

### 特徵有效性

| AUC 範圍 | 判定 | 後續行動 |
|:---|:---|:---|
| > 0.6 | 有效 | 納入組合模型 |
| 0.55 - 0.6 | 弱有效 | 與其他特徵組合測試 |
| < 0.55 | 無效 | 記錄並跳過 |

### 過濾策略有效性

滿足以下 **全部** 條件：
1. **FP/TP 移除比 > 1.45**
2. **新 F1 > 0.8155** (baseline)
3. **策略可解釋**（非過擬合）

---

## 實驗優先順序

### Phase 1: 快速驗證

| 順序 | 實驗 | 工具 |
|:---:|:---|:---|
| 1 | VCF QUAL 分布分析 | bcftools + Python |
| 2 | VCF AF 分布分析 | bcftools + Python |
| 3 | Strand Bias 分析 | bcftools + Python |
| 4 | Regional Clustering | Python |
| 5 | 特徵組合 AUC | sklearn |

### Phase 2: 深度分析

| 順序 | 實驗 | 工具 |
|:---:|:---|:---|
| 6 | VCF × 甲基化交互 | Python |
| 7 | 多閾值 Grid Search | Python |
| 8 | 機器學習特徵選擇 | sklearn |
| 9 | samtools 抽樣驗證 | samtools |
| 10 | SEQC2 比對 | bcftools |

### Phase 3: 程式碼實驗

需開 git 分支，修改 C++ 程式碼

---

## 常用命令

### 提取 VCF 特徵
```bash
bcftools query -f '%CHROM\t%POS\t%QUAL\t[%AF]\t[%DP]\n' \
  data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz > tp_features.tsv
```

### 觀察特定區域
```bash
samtools view -h data/bam/HCC1395/tumor.bam chr1:877000-878000 | less
```

### Python 分析範例
```python
import pandas as pd
from sklearn.metrics import roc_auc_score

tp = pd.read_csv('output/.../filtered_snv_tp/significance_summary.csv')
fp = pd.read_csv('output/.../filtered_snv_fp/significance_summary.csv')

combined = pd.concat([tp.assign(label=1), fp.assign(label=0)])
auc = roc_auc_score(combined['label'], combined['CramersV'])
print(f"CramersV AUC: {auc:.4f}")
```

---

## 關鍵資源

| 資源 | 位置 |
|:---|:---|
| TP VCF | `data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz` |
| FP VCF | `data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz` |
| SEQC2 金標準 | `data/answer/SEQC2/` |
| TP 結果 | `output/.../filtered_snv_tp/significance_summary.csv` |
| FP 結果 | `output/.../filtered_snv_fp/significance_summary.csv` |
| 分析腳本 | `tools/compare_vcf_results.py` |
| F1 優化腳本 | `tools/f1_optimization_analysis.py` |

---

## 聯絡與記錄

- 計劃文件：`docs/plans/2026/01/20260113_甲基化F1研究計劃_01.md`
- 數據說明：`data/README.md`
- AI 對話記錄：`docs/provenance/ai_sessions/2026/01/`
