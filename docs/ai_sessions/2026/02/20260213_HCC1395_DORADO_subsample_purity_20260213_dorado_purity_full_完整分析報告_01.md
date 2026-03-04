<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# HCC1395_DORADO Subsample Purity（run-tag: 20260213_dorado_purity_full）完整分析報告
- 產生時間：2026-02-13
- 執行腳本：`scripts/analysis/run_purity_and_standard_verification.sh`
- Run-tag：`20260213_dorado_purity_full`（獨立新目錄，不覆蓋舊結果）
- 樣本來源：`/big8_disk/data/HCC1395/ONT_Dorado/subsample/`
- Truth VCF：`/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
- Truth BED：`/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed`
- Truth SNV（BED 內）總數：**39447**

## 0. 知識依據（Knowledge）
- 根據 `Knowledge/02_samples/subsample_purity.md`：`t00_n25` 為 0% tumor baseline，主 purity 驗證可排除。
- 根據 `Knowledge/03_file_formats/bam_format.md`：MM/ML 為甲基化必要標籤，缺失會影響下游甲基化分析。
- 根據 `Knowledge/06_workflows/methylation_analysis.md`：InterSubMod 需依賴 BAM 中的修飾訊號（MM/ML）。

## 1. 執行覆蓋與流程完成性
### 1.1 Subsample 覆蓋
`ONT_Dorado/subsample` 實際存在：`t00_n25`, `t7_n29`, `t19_n29`, `t30_n20`, `t40_n10`, `t50_n00`。

本次設定：
- `t00_n25`：**跳過**（zero-tumor，符合需求）
- 其餘 5 組：**全部完成**

### 1.2 每組步驟完成狀態
來源：`purity_status.tsv`

| Subsample | 推估 purity | Step01 LongPhase-S | Step02 InterSubMod | Step03 Filter | cleanup |
|---|---:|---|---|---|---|
| t7_n29 | 19.4% | OK | OK | OK | REMOVED |
| t19_n29 | 39.6% | OK | OK | OK | REMOVED |
| t30_n20 | 60.0% | OK | OK | OK | REMOVED |
| t40_n10 | 80.0% | OK | OK | OK | REMOVED |
| t50_n00 | 100.0% | OK | OK | OK | REMOVED |

### 1.3 空間清理驗證
- 所有 purity 執行後 `HCC1395_DORADO_tagged.bam` 皆已刪除。
- `find ... -name 'HCC1395_DORADO_tagged.bam'` 回傳空結果。

## 2. GQ 解析問題與流程修正驗證
### 2.1 問題現象
`pileup_filter.vcf` header 定義為 `GQ Type=Integer`，但內容有 `GQ=.` / 浮點狀態，直接給 `bcftools` 會報錯。

實測（不修正）
```bash
bcftools view -f PASS pileup_filter.vcf -Ov -o /tmp/no_patch.vcf
# [E::vcf_parse_format_fill5] Invalid character '.' in 'GQ' FORMAT field ...
```

### 2.2 已納入流程的修正
`scripts/pipeline/steps/01_longphase_s.sh` 會先把 header 的：
- `ID=GQ,Number=1,Type=Integer`
改為：
- `ID=GQ,Number=1,Type=Float`
再送進 `bcftools view -f PASS`。

實測（修正後）
- `rc=0`，可正常輸出 PASS VCF。

## 3. MM/ML 是否因混合流程遺失：驗證方法與結果
### 3.1 驗證方法
1. 檢查來源 BAM（tumor/normal）是否有 MM/ML。
2. 以 `samtools view` 抽樣建立小型 tumor+normal 混和 BAM。
3. 比較混和前後 MM/ML 計數是否維持。

### 3.2 結果（DORADO）
- `HCC1395_DORADO tumor`（前 2000 reads）：MM=2000, ML=2000
- `HCC1395_DORADO normal`（前 2000 reads）：MM=2000, ML=1999
- 小型混和 BAM（前 2000 reads）：MM=2000, ML=2000

結論：**混和操作本身不會必然造成 MM/ML 消失**。若後續遇到 MM/ML=0，應優先檢查來源 BAM 或上游流程是否已失去標籤。

## 4. 各 purity 三階段指標（TP/FP/FN/P/R/F1）
來源：
- `purity_status.tsv`
- `purity_stage_metrics.tsv`

| Subsample | Purity | MM/ML | ClairS TP | ClairS FP | ClairS FN | ClairS P | ClairS R | ClairS F1 | LongPhase TP | LongPhase FP | LongPhase FN | LongPhase P | LongPhase R | LongPhase F1 | Final TP | Final FP | Final FN | Final P | Final R | Final F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| t7_n29 | 19.4 | 1000/1000 | 12090 | 176 | 27357 | 0.9857 | 0.3065 | 0.4676 | 12042 | 45 | 27405 | 0.9963 | 0.3053 | 0.4674 | 10173 | 42 | 29274 | 0.9959 | 0.2579 | 0.4097 |
| t19_n29 | 39.6 | 1000/1000 | 24476 | 527 | 14971 | 0.9789 | 0.6205 | 0.7595 | 24447 | 136 | 15000 | 0.9945 | 0.6197 | 0.7636 | 23886 | 132 | 15561 | 0.9945 | 0.6055 | 0.7527 |
| t30_n20 | 60.0 | 1000/1000 | 28611 | 1036 | 10836 | 0.9651 | 0.7253 | 0.8282 | 28546 | 255 | 10901 | 0.9911 | 0.7237 | 0.8366 | 28499 | 252 | 10948 | 0.9912 | 0.7225 | 0.8358 |
| t40_n10 | 80.0 | 1000/1000 | 30044 | 2023 | 9403 | 0.9369 | 0.7616 | 0.8402 | 29902 | 489 | 9545 | 0.9839 | 0.7580 | 0.8563 | 29890 | 479 | 9557 | 0.9842 | 0.7577 | 0.8563 |
| t50_n00 | 100.0 | 1000/1000 | 30642 | 3510 | 8805 | 0.8972 | 0.7768 | 0.8327 | 30409 | 774 | 9038 | 0.9752 | 0.7709 | 0.8611 | 30394 | 754 | 9053 | 0.9758 | 0.7705 | 0.8611 |

## 5. F1 變化與特殊現象
### 5.1 ClairS → LongPhase-S
- 平均 `ΔF1 = +0.0114`
- 5 組中 4 組提升，1 組（t7_n29）微幅下降 `-0.0002`。
- 最高提升出現在 `t50_n00`：`+0.0284`。

### 5.2 LongPhase-S → InterSubMod + Methylation Filter（目前預設標準）
- 平均 `ΔF1 = -0.0139`
- t40_n10、t50_n00 為持平（四捨五入 0.0000）
- t30_n20、t19_n29 輕微下降
- **t7_n29 明顯下降（-0.0577）**：TP 被移除 1869、FP 僅移除 3，Recall 下降明顯。

### 5.3 現象判讀
1. 在低 purity（19.4%）下，預設 methylation filter 過於嚴格，對 TP 傷害遠大於 FP 改善。
2. 在中高 purity（60%~100%）下，filter 影響趨近中性（F1 接近不變）。
3. 因 Dorado subsample 具備 MM/ML，這次結果可視為「可有效驗證甲基化濾波行為」；結論是：
   - **目前單一固定門檻不適合所有 purity**。
   - 需要 purity-aware（依 purity 調整）或更保守的條件。

## 6. 結論
1. `HCC1395_DORADO` 的非 zero-tumor purity 已全部完整跑完，流程與輸出正確，且未覆蓋舊結果。
2. `GQ` 造成 `bcftools` 解析失敗的問題已在 pipeline 內處理（Integer→Float patch）。
3. 各 purity 的 `TP/FP/FN/Precision/Recall/F1` 三階段數據已完整驗證。
4. 目前預設 Methylation Filter 在此批 purity 的整體效果 **未帶來穩定 F1 提升**，尤其低 purity 會明顯傷害 Recall。
5. 下一步擴充到其他樣本前，建議先導入 purity-aware filter 或先鎖定中高 purity（>=60%）條件測試。

## 7. 跨樣本 MM/ML 快速核對（前 1000 reads）
| 樣本 | MM 計數 | ML 計數 | 判讀 |
|---|---:|---:|---|
| HCC1395/ONT | 0 | 0 | 無甲基標籤 |
| HCC1395/ONT_Dorado | 1000 | 1000 | 有甲基標籤 |
| COLO829/ONT_PAO | 1000 | 999 | 有甲基標籤 |
| H1437/ONT | 1000 | 1000 | 有甲基標籤 |
| H2009/ONT | 1000 | 1000 | 有甲基標籤 |
| HCC1937/ONT | 768 | 768 | 有甲基標籤（比例較低） |
| HCC1954/ONT | 1000 | 1000 | 有甲基標籤 |

## 8. 主要輸出位置
- Run log：`/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395_DORADO/purity_runs/20260213_dorado_purity_full/run.log`
- Purity 狀態：`/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395_DORADO/purity_runs/20260213_dorado_purity_full/purity_status.tsv`
- 階段指標：`/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395_DORADO/purity_runs/20260213_dorado_purity_full/purity_stage_metrics.tsv`
- 各 purity 輸出：`/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395_DORADO/purity_t*_20260213_dorado_purity_full/`
