<!--
建立時間: 2026-06-22
目標: Q1 論文 motivation 量化基石 — HCC1395 sSNV 間距 × ONT 讀長覆蓋
處理範圍: 單樣本 HCC1395（paired pileup TP sSNV + tumor ONT BAM），全基因組（autosomes）
關聯檔案:
  - scripts/analysis/ssnv_spacing_coverage.py
  - state/analysis_outputs/ssnv_spacing_coverage.json
-->
---
title: Q1 — HCC1395 sSNV 間距 × ONT 讀長覆蓋（論文 motivation 量化基石）
date: 2026-06-22
status: in_progress
task: T-Q1-COVERAGE
audience: 論文作者 + 跨 agent
partial_flag: "PARTIAL — 單樣本（HCC1395）；跨樣本一致性未做"
data_sources:
  - state/analysis_outputs/ssnv_spacing_coverage.json
  - /big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/filtered_snv_tp.vcf.gz
  - /big8_disk/liaoyoyo2001/data/vcf/ClairS_ssrs/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/pileup/filtered_snv_fp.vcf.gz
  - /big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam
build_branch: feat/q1-coverage
related_memory:
  - project_subclone_snv_difficulty_methylation_framework
---

# Q1 — sSNV 間距 × ONT 讀長覆蓋

> **一句話結論**：HCC1395 全基因組 30,490 個 somatic sSNV 平均每 ~94 kb 才一個（中位相鄰間距 51.4 kb），而 ONT tumor read 的基因組跨度中位僅 4.9 kb；**單一 read 能同時連到 ≥2 個 sSNV 的比例只有 ~0.7–1.0%**，94% 的 read 一個 sSNV 都碰不到。這就是「只用 sSNV 在 read 層直接重建 subclone 困難」的量化基石，也正當化把甲基化當每-read 的補位訊號。

> ⚠ **PARTIAL flag**：本報告為**單樣本（HCC1395）**結果。跨 6 樣本一致性未做（屬後續 WS-B 任務）。所有數字來源見 §6 溯源表。

---

## 1. 資料與口徑（provenance）

| 項目 | 值 / 路徑 | 來源 |
|---|---|---|
| sSNV 集（主） | TP = ClairS pileup 對 SEQC2 truth 確認的 somatic 集 | `filtered_snv_tp.vcf.gz` |
| sSNV 集（對照） | FP = ClairS pileup 落在 truth 外的 call | `filtered_snv_fp.vcf.gz` |
| tumor ONT BAM | `HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam`（283 GB） | symlink `data/bam/HCC1395/tumor.bam` |
| pipeline 對應 | `scripts/run_vcf_all_snv.sh:46-48`（paired 模式 canonical 路徑） | run script |
| 染色體範圍 | autosomes chr1–22；**0 個 non-autosome、0 個 non-PASS** | JSON `ssnv_totals` |
| 讀長定義 | **reference-consumed span**（`reference_end − reference_start`）= read 覆蓋的基因組區間（決定能碰幾個 sSNV 座標）；query/SEQ 長另記 | JSON `scope` |
| 基因組大小分母 | autosome 長度總和（hg38，取自 BAM header）= 2,875,001,522 bp | JSON `genome` |

> **口徑選擇理由**：論文 motivation 問的是「真 somatic sSNV 有多稀疏、read 能否連 ≥2 個」。TP 集（= ISM 消費為 somatic 的集合）是正確分母；FP 僅作上下文對照。讀長用 reference-span 而非 SEQ-length，因為「跨幾個 sSNV」由基因組區間決定（read 含 indel 時 SEQ 長與基因組跨度不同）。

---

## 2. (a) sSNV 總數

| 集 | 總記錄 | PASS | autosome 用於分析 | non-autosome |
|---|---|---|---|---|
| **TP（主）** | **30,490** | 30,490（100%） | **30,490** | 0 |
| FP（對照） | 4,842 | — | 4,842 | — |

- TP = **30,490**，與交接包「~30,490 sSNV」精確吻合；全部 PASS、全部落 autosome（無 chrX/Y/alt contig）。

## 3. (b) 全基因組相鄰 sSNV 間距分布（TP）

| 統計 | 值 (bp) |
|---|---|
| 相鄰 pair 數 | 30,468（= 30,490 − 22 染色體首位點） |
| **中位間距** | **51,425** |
| 平均間距 | 87,809 |
| 均勻假設下期望間距（genome/N） | 94,293 |
| p10 / p25 / p75 / p90 | 5,466 / 18,983 / 112,361 / 202,543 |
| min / max | 1 / 18,567,723 |
| 平均密度 | **10.6 sSNV / Mb**（= 1.06 / 100 kb） |
| 相鄰兩 sSNV 落入同一 ±5 kb window（間距 ≤10 kb） | **15.4%** |

> **解讀**：中位間距 51.4 kb < 均勻期望 94.3 kb → sSNV 呈**過度離散群聚**（Poisson over-dispersed），少數熱點區密集、大片區域稀疏（最大空隙 18.6 Mb）。只有 15.4% 的相鄰 sSNV 對近到可被同一個 ±5 kb 甲基 window 同時涵蓋。

## 4. (c) ONT tumor read 長度分布（取樣 200,000 reads，primary mapped）

| 統計 | reference-span (bp) | query-len (bp) |
|---|---|---|
| 中位 | **4,898** | 5,016 |
| 平均 | 9,625 | 9,841 |
| p10 / p25 / p75 | 700 / 1,555 / 13,864 | — |
| p90 / p95 / p99 | 25,877 / 33,985 / 49,251 | — |

> **解讀**：讀長右偏（中位 ~4.9 kb，但長尾到 p99 ~49 kb）。reference-span 與 query-len 中位接近 → indel 對基因組跨度影響不大，read 主要由純比對段構成。

## 5. (d) 單讀可跨 ≥2 個 sSNV 的比例

**兩法交叉驗證**（密度 = 30,490 / 2.875 Gb = 1.06×10⁻⁵ sSNV/bp）：

| 法 | P(≥1 sSNV) | **P(≥2 sSNV)** | 備註 |
|---|---|---|---|
| **解析（Poisson，按讀長分布加權）** | 9.1% | **0.99%** | λ = density × ref_span，對 200k read 跨度分布平均 |
| **經驗（50,000 reads 實際與 sSNV 座標 intersect）** | 5.8% | **0.74%** | 直接 bisect read 區間內 sSNV 數 |

**經驗 sSNV-per-read 直方圖**（50,000 reads）：

| read 上 sSNV 數 | read 數 | 佔比 |
|---|---|---|
| 0 | 47,081 | **94.2%** |
| 1 | 2,549 | 5.1% |
| 2 | 322 | 0.64% |
| 3 | 45 | 0.09% |
| 4 | 3 | 0.006% |

- 單一 read 上最多只看到 **4** 個 sSNV。
- **94.2% 的 read 一個 somatic sSNV 都碰不到** → 絕大多數 read 對「以 sSNV 連結 subclone lineage」零貢獻。
- 經驗 P(≥2)=0.74% 略低於 Poisson 0.99%：解析法假設密度均勻，但實際 sSNV 群聚 + 經驗樣本取 chr1-start 起 genome-ordered read，兩者口徑差異所致；**兩法同數量級（~1%）**。

> **核心 motivation**：即使 ONT read 中位 ~5 kb、長尾到數十 kb，受限於 somatic sSNV 的稀疏密度（10.6/Mb），**單一 read 幾乎無法直接 phase 兩個 somatic 變異**（~99% 的 read 跨 0 或 1 個 sSNV）。bulk 上「以 read 同時連結 ≥2 sSNV 重建 lineage」在物理上極度受限——這正是需要正交、每-read 都有的訊號（甲基化）來補位的量化理由。

## 6. (e) ±5000 bp window 對讀長分布的涵蓋率

ISM 每個 sSNV 取 ±5 kb（= 10 kb 寬）甲基 window。對照 read 跨度分布：

| 指標 | 值 |
|---|---|
| read ref-span ≥ 全 window（10 kb） | **33.2%** |
| read ref-span ≥ 半 window（5 kb，可從一側完整觸及 SNV） | **49.5%** |
| 中位 read 涵蓋的完整 window 數（median span / 10 kb） | 0.49 |

> **解讀**：約 1/3 的 read 跨度足以涵蓋完整 ±5 kb window，約一半至少達半 window。中位 read（4.9 kb）約涵蓋半個 window → ISM 的 ±5 kb 甲基 window 設定與實際讀長分布相稱：多數 read 能貢獻 SNV 鄰近的甲基 context（即便該 read 只碰一個 sSNV）。

---

## 7. 對論文 motivation 的意涵（敘述紅線內）

1. **「只用 sSNV 難」已量化**：30,490 sSNV、中位間距 51 kb、密度 10.6/Mb，使 ~99% 的 ONT read 無法同時連結 ≥2 個 somatic sSNV（empirical 0.74% / Poisson 0.99%）。bulk read 層 sSNV co-occurrence 物理上稀疏。
2. **正當化甲基補位（不過度宣稱）**：甲基化是每-read 都有的密集正交訊號，可在 sSNV 連不到的 ~94% read 上提供額外結構資訊。**此處僅論「sSNV 稀疏 → 需要補位訊號」，不宣稱甲基化驅動偵測 / 當 filter / genome-wide tree**（守 cognition 文件 §5 紅線 1）。
3. **±5 kb window 設定合理**：與實際讀長分布相稱（1/3 read 跨完整 window，半數跨半 window）。

---

## 8. 限制與後續

- **單樣本**：HCC1395 only。跨 6 樣本一致性（H1437/H2009/HCC1937/HCC1954/COLO829）未做。
- **sSNV 集口徑**：用 ClairS pileup TP（vs SEQC2）；若改用其他 truth 定義或含 FP，密度會變（FP 另有 4,842 個，未納入主分析）。
- **取樣**：讀長/經驗 intersect 取 genome-ordered 前 N read（chr1-start 起），非隨機跨基因組均勻取樣；對 read 長度分布的偏差預期小（mapped primary read），但 P(≥2) 經驗值對取樣口徑略敏感。
- **後續**：(i) 跨樣本重跑同 script；(ii) 若要更嚴謹的 P(≥2)，可改全 BAM 掃描而非前 50k read。

---

## 9. (溯源) 每個數字 → 檔案

所有數字取自 `state/analysis_outputs/ssnv_spacing_coverage.json`（本輪 `scripts/analysis/ssnv_spacing_coverage.py` 跑出、exit 0、Read 回確認）：

| 數字 | JSON key |
|---|---|
| sSNV 30,490 / PASS / autosome | `ssnv_totals.tp_n_total_records` / `tp_n_pass` / `tp_n_autosome_used` |
| FP 4,842 | `ssnv_totals.fp_n_total_records` |
| genome 2,875,001,522 / 密度 10.6/Mb | `genome.autosome_total_bp` / `spacing.mean_density_ssnv_per_mb` |
| 中位間距 51,425 / 平均 87,809 / 期望 94,293 | `spacing.gap_median_bp` / `gap_mean_bp` / `expected_uniform_spacing_bp` |
| 間距 ≤10kb 15.4% | `spacing.frac_gaps_le_window_2x` |
| 讀長中位 4,898 / 平均 9,625 / p99 49,251 | `read_length.ref_span_median_bp` / `ref_span_mean_bp` / `ref_span_percentiles_bp.p99` |
| P(≥2) Poisson 0.99% / 經驗 0.74% | `span_analytic_poisson.P_ge2_ssnv_read_weighted` / `span_empirical.P_ge2_ssnv` |
| 94.2% read 0 sSNV / hist / max 4 | `span_empirical.frac_reads_0_ssnv` / `ssnv_per_read_hist` / `max_ssnv_on_single_read` |
| window 33.2% / 49.5% / 0.49 | `window_coverage.frac_reads_span_ge_full_window` / `frac_reads_span_ge_half_window` / `median_windows_covered_per_read` |
