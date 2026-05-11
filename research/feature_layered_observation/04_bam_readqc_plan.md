# 04 — BAM Read-Level QC 抽取計畫 (Group G4)

**建立日期**：2026-04-23
**作者**：AI session (Phase B2)
**對應觀察組別**：G4 — BAM read QC

## 1. 目的

從每個樣本的 `tumor_tagged.bam` 抽取 per-region read-level QC 統計，供 G4
觀察使用：檢查 MapQ、NM、softclip、strand bias、read length 等是否在 TP vs
FP region 之間存在系統性差異。

## 2. 特徵清單

| 欄位 | 定義 | 單位 |
|------|------|------|
| `n_reads` | 視窗內 primary+mapped+非 dup 的 read 數 | count |
| `MapQ_mean` | MapQ 平均 | 0-60 |
| `MapQ_median` | MapQ 中位數 | 0-60 |
| `MapQ_p10` | MapQ 10 百分位 | 0-60 |
| `MapQ_p90` | MapQ 90 百分位 | 0-60 |
| `LowMQ20_Frac` | MapQ<20 reads 比例 | 0-1 |
| `NM_mean` | NM tag 平均（mismatch 數）| count |
| `Softclip_Frac` | 有 softclip 的 read 比例 | 0-1 |
| `Strand_Bias` | |forward - reverse| / n_reads | 0-1 |
| `Read_Length_mean` | read.query_length 平均 | bp |

Window = 5000 bp（Pos ± 2500），與 ISM pipeline region window 一致。

排除條件（與 primary-read analysis 對齊）：`is_secondary`、`is_supplementary`、
`is_duplicate`、`is_unmapped`。

## 3. BAM 路徑對應表 (14 筆)

| Sample | Mode | BAM |
|---|---|---|
| HCC1395 | paired_full | `output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam` |
| HCC1395_DORADO | paired_full | `output/canonical/HCC1395_DORADO/paired_full/20260315_HCC1395_DORADO_paired_full_full_complete_matrix/longphase_s/HCC1395_DORADO_tagged.bam` |
| H1437 | paired_full | `output/canonical/H1437/paired_full/20260315_H1437_paired_full_full_complete_matrix/longphase_s/H1437_tagged.bam` |
| H2009 | paired_full | `output/canonical/H2009/paired_full/20260315_H2009_paired_full_full_complete_matrix/longphase_s/H2009_tagged.bam` |
| HCC1937 | paired_full | `output/canonical/HCC1937/paired_full/20260315_HCC1937_paired_full_full_complete_matrix/longphase_s/HCC1937_tagged.bam` |
| HCC1954 | paired_full | `output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix/longphase_s/HCC1954_tagged.bam` |
| COLO829 | paired_full | `output/canonical/COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix/longphase_s/COLO829_tagged.bam` |
| HCC1395 | to_pileup | `output/20260307_hcc1395_to_pilot_1/step03_longphase_to/tumor_tagged.bam` |
| HCC1395_DORADO | to_pileup | `output/synthesis/research_rounds/archive/202603_early_pilots/20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_tagged.bam` |
| H1437 | to_pileup | `.../20260318_h1437_to_pilot_fastresume/step03_longphase_to/tumor_tagged.bam` |
| H2009 | to_pileup | `.../20260318_h2009_to_pilot_fastresume/step03_longphase_to/tumor_tagged.bam` |
| HCC1937 | to_pileup | `.../20260318_hcc1937_to_pilot_fastresume/step03_longphase_to/tumor_tagged.bam` |
| HCC1954 | to_pileup | `.../20260318_hcc1954_to_pilot_fastresume/step03_longphase_to/tumor_tagged.bam` |
| COLO829 | to_pileup | `output/synthesis/research_rounds/20260423_colo829_to_pilot/step03_longphase_to/tumor_tagged.bam` |

全部 14 筆 BAM + BAI 已確認存在。

**注意**：`master_extended.tsv.gz` 目前僅涵蓋 13 組合（COLO829 to_pileup 的
master rows 尚未併入）。該組合的輸出會是 empty（0 regions），不影響 QC
抽取邏輯；後續 master 更新後可用 `--force` 重跑該一筆。

## 4. 樣本 × 變異數

| Sample | paired_full | to_pileup | 合計 |
|---|---:|---:|---:|
| COLO829 | 37,458 | 0 (master pending) | 37,458 |
| H1437 | 67,476 | 58,915 | 126,391 |
| H2009 | 132,995 | 137,695 | 270,690 |
| HCC1395_DORADO | 30,129 | 40,428 | 70,557 |
| HCC1395 | 30,381 | 40,115 | 70,496 |
| HCC1937 | 12,588 | 24,655 | 37,243 |
| HCC1954 | 17,938 | 67,286 | 85,224 |
| **總計** | **328,965** | **369,094** | **698,059** |

## 5. 並行策略

- 單 BAM：`pysam.AlignmentFile(..., threads=4)`（解壓 4 執行緒）
- 跨 BAM：`xargs -P 4`（4 BAMs 同時處理）
- 全域 concurrent threads 上限：**16**（節制 I/O，避免 NVMe/BAM 解碼飽和）
- 實測 warm cache：HCC1395 paired chr19 500 regions in 92 s ≈ 5.4 reg/s（60× 深覆蓋；ONT ultra-long）
- 預估跨 14 組合 total：≈ 698k / (5 reg/s × 4 parallel) ≈ 9.7 hr 並行 wall-time；實際依樣本覆蓋差異

## 6. 輸出規格

路徑：`research/feature_layered_observation/data/bam_readqc/bam_readqc_{sample}_{mode}.tsv`

欄位（17）：`RegionID, sample, mode, Chr, Pos, win_start, win_end, n_reads, MapQ_mean, MapQ_median, MapQ_p10, MapQ_p90, LowMQ20_Frac, NM_mean, Softclip_Frac, Strand_Bias, Read_Length_mean`

Atomic write 規範：寫入 `.tmp`，成功後 `rename` 為正式檔；script 可 resume
（有 output 且 non-empty 則 skip）。

## 7. 執行方式

### Smoke test (已通過)

```bash
python3 research/feature_layered_observation/scripts/bam_readqc_per_region.py \
    --sample HCC1395 --mode paired_full --chr chr19 --limit 20
```

結果：`bam_readqc_HCC1395_paired_full_chr19.tsv` 20 rows、17 cols。
MapQ_mean=59.6-60、Softclip_Frac≈0.86-0.95、Strand_Bias≈0.03-0.27、
Read_Length_mean≈15-23 kb。

### 全量背景

```bash
nohup bash research/feature_layered_observation/scripts/run_bam_readqc_all.sh \
    > research/feature_layered_observation/data/bam_readqc_nohup.out 2>&1 &
```

## 8. 驗證標準（完成後）

1. 14 個 TSV 均存在、含 header + >0 data rows（COLO829 TO 例外：0 rows 可接受）
2. 每個 TSV 的 row count = master_extended 對應 (sample, mode) rows 數
3. 各欄位無全 NaN；`n_reads>0` rate 應 >95%（variant 位置通常有 coverage）
4. Join back to master_extended 的 RegionID 無 mismatch（以 RegionID 作 key）

## 9. 預期用途

G4 觀察：
- 比對 TP vs FP regions 的 MapQ/NM/Softclip 分佈差異
- 檢查 LowMQ20_Frac 是否與 germline FP 相關
- Strand_Bias 異常 region 是否 enrich 於特定 sample/chr
- Read_Length_mean 極短 region 是否 QC 偏誤 hotspot
