<!--
建立時間: 2026-01-13 19:10
目標: VCF 欄位特徵分析記錄
處理範圍: QUAL, AF, DP, Strand Bias 等 VCF 欄位
-->

# VCF 欄位分析

## 待分析欄位

| 欄位 | 說明 | 提取命令 | 狀態 |
|:---|:---|:---|:---|
| QUAL | ClairS 品質分數 | `bcftools query -f '%QUAL\n'` | 待分析 |
| AF | Variant Allele Frequency | `bcftools query -f '[%AF]\n'` | 待分析 |
| DP | Read Depth | `bcftools query -f '[%DP]\n'` | 待分析 |
| H flag | Haplotype-specific | `bcftools view -H \| grep ';H;'` | 待分析 |
| Strand Bias | FAU/RAU 等 | 需自訂腳本 | 待分析 |

## 分析結果

### QUAL 分析

（待填寫）

### AF 分析

（待填寫）

### DP 分析

（待填寫）

### Strand Bias 分析

（待填寫）
