<!--
建立時間: 2026-06-22
類型: 工具執行 + 三方交叉比對結果
狀態: in_progress
主軸: Subclonal reconstruction — CN/SV 整合層
data_sources: /big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/wakhan_out/2.77_0.92_0.8/, .../savana_wgs/cna_normalhet/, /big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed
build_branch: research/subclonal-reconstruction-202606
-->

# Wakhan HCC1395 haplotype-specific CN + 三方交叉比對（Wakhan ↔ SAVANA ↔ SEQC2）

> **版本** v1.0 (2026-06-22) · 單樣本 **⭐3** · 證據 tier：**L1**=第一手磁碟數據重現
> 腳本 `scripts/analysis/compare_wakhan_savana_seqc2.py`

## 0. 一句話結論

**Wakhan 在 HCC1395 跑通且與 SAVANA 獨立收斂**：Wakhan rank-1 解 **purity 0.92 / ploidy 2.77** ≈ SAVANA `cna_normalhet` **0.96 / 2.79**（兩個獨立 ONT CN caller 自動選到同一 near-triploid 解）；LOH 三方一致（Wakhan↔SEQC2 Jaccard **0.819**，Wakhan↔SAVANA 0.832，三方一致 88.2%）。Wakhan 額外提供 **HP1/HP2 分離 integer CN + SV-linked 段邊界**（SAVANA 無）。

## 1. 執行（job bhmc9i1nt）

| 步驟 | 耗時 | 結果 |
|---|---|---|
| Severus（tumor-normal, `--phasing-vcf germline_phased`）| 134 min | somatic SV **1,789**（巧合接近 SEQC2 Talsania truth 1,777）/ all SV 22,588 |
| Wakhan `all`（`--normal-phased-vcf germline_phased --breakpoints severus_somatic.vcf --purity-range 0.9-1.0`）| 108 min | 2 解，rank-1=2.77/0.92 |

- 輸入：tumor `HCC1395.bam`(292GB) + ref GRCh38_no_alt + **`germline_phased_merged.vcf.gz`**（Clair3 on HCC1395BL，normal-het phased，避 SAVANA -vg 汙染）+ Severus somatic SV。
- Wakhan 2 解（purity-ploidy 識別曲線）：**rank-1 `2.77_0.92`（conf 0.8）** / rank-2 `1.42_0.96`（conf 0.51，near-diploid 替代支，= SAVANA -vg 版誤入的解）。**Wakhan 自動選對 near-triploid 解**。
- 輸出：`wakhan_out/2.77_0.92_0.8/bed_output/`：`*_segments_HP_1.bed` / `HP_2.bed`（711 段，欄 chr/start/end/coverage/copynumber_state/confidence/**svs_breakpoints_ids**）+ `loh_regions.bed`(129) + subclonal segments + per-chr CN HTML。

## 2. 三方交叉比對（10kb bin, chr1-22, 共同覆蓋 2,801 Mb，L1）

### purity / ploidy — 跨工具收斂 ✅
| 工具 | purity | ploidy |
|---|---|---|
| Wakhan rank-1 | 0.92 | 2.77 |
| SAVANA cna_normalhet | 0.96 | 2.79 |

兩獨立工具（SAVANA normal-het + Wakhan phased-germline）**收斂到同一 near-triploid 解** → 強化此解的可信度（純細胞系 purity~1.0、hyper-triploid 符 HCC1395 karyotype）。

### LOH 三方 Jaccard ✅
| 對 | LOH 覆蓋率 | Jaccard |
|---|---|---|
| Wakhan vs SEQC2 | 58.5% vs 53.3% | **0.819** |
| SAVANA vs SEQC2 | 55.2% vs 53.3% | **0.962** |
| Wakhan vs SAVANA | — | **0.832** |
| 三方一致（all-LOH 或 all-not）| — | **88.2%** |

→ Wakhan LOH 與 truth 高度一致（0.819）但**略 over-call**（58.5%），不如 SAVANA 緊（0.962）。三方 88.2% 一致。

### total CN — 中度一致
- Wakhan(HP1+HP2) vs SAVANA copyNumber：**Spearman 0.699 / Pearson 0.546**（n=280,075 bin）。
- 中位數 Wakhan 3.00 ≈ SAVANA 2.86（皆 ~triploid baseline）。
- ⚠ SAVANA 高 CN 達 19.6（chr8 大擴增），**Wakhan 封頂 ~9** → 高擴增區壓縮是 Pearson 偏低主因；Spearman(rank) 0.699 較穩。

## 3. 兩工具定位（互補）

| | SAVANA cna_normalhet | Wakhan |
|---|---|---|
| LOH vs SEQC2 | **0.962**（最緊）| 0.819（略 over-call）|
| allele CN | total + minor | **HP1/HP2 分離 integer CN**（不可替代）|
| SV 整合 | 可選 -bp | **段邊界連 Severus SV ID** |
| subclonal CN | （未啟用）| `*_subclonal_segments_HP_*.bed` |
| purity/ploidy | 0.96/2.79 | 0.92/2.77（**獨立佐證**）|

→ **SAVANA = LOH/total CN 最準的固定參考**；**Wakhan = haplotype-resolved + SV-linked CN 的加值層**。兩者交叉確認 near-triploid 解。

## 4. caveat / tier
- 單樣本 ⭐3；跨樣本需 COLO829。
- Wakhan LOH 略 over-call（0.819 < SAVANA 0.962）；high-CN 封頂 9（vs SAVANA 19.6）。
- Wakhan rank-2（1.42/0.96 near-diploid）= 識別曲線替代支，**勿用**（conf 0.51 < rank-1 0.8）。
- SEQC2 = segment-level 定性 truth；LOH 為穩健軸。
- Severus somatic SV 1,789 未獨立 vs Talsania truth 逐一比對（數量巧合，需正式 concordance 才可宣稱）。

## 5. Provenance
- 腳本 `scripts/analysis/compare_wakhan_savana_seqc2.py`；pipeline log `cnv_sv_work/logs/wakhan_pipeline.log`。
- L1 數字皆 `wakhan_out/2.77_0.92_0.8/bed_output/` + `cna_normalhet/` + SEQC2 truth 第一手 bin 比對。
