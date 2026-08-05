<!--
建立時間: 2026-06-20
類型: 工具可行性執行結果 (A-pilot execution result)
狀態: in_progress（全基因組固定參考層執行中）
主軸: Subclonal reconstruction — CN/SV 輔助整合層
plan: InterSubMod/docs/plans/20260620_ont_cnv_sv_subclone_verification_feasibility_01.md
data_sources: /big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_chr2_chr8/cna/HCC1395_segmented_absolute_copy_number.tsv, /big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_chr2_chr8/cna/HCC1395_fitted_purity_ploidy.tsv, /big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed
-->

# SAVANA ONT CNV/SV 可行性確認 — HCC1395 chr2+chr8 vs SEQC2

> **版本** v1.0 (2026-06-20) · **裁決：可行性 CONFIRMED（GREEN）** · 全基因組固定參考層執行中
> 證據 tier：**L1**=第一手磁碟數據重現重算 ｜ L2=工具文件 ｜ L3=推論

---

## 0. 一句話結論

**SAVANA 1.3.7 在我們的 ONT HCC1395（R10 5kHz, Dorado）跑得通，且 CN/LOH 與 SEQC2 truth 高度一致**（chr8 LOH Jaccard **0.977**、chr2 **0.920**），purity 估計 **0.96**（純細胞系期望 ~1.0，命中）。工具可用性確認 → 啟動全基因組「固定 CN 參考層」供後續 k_ISM vs k_CN 整合分析。

---

## 1. 環境與輸入（L1）

| 項目 | 值 |
|---|---|
| 工具 | SAVANA 1.3.7（bioconda, env `cnvtools`）；severus 1.7；Wakhan 0.4.3（GitHub，待移獨立 env）|
| 安裝 | conda classic solver 卡死 → 改 `--solver=libmamba`（已裝於 base）；SAVANA/Wakhan **不在 PyPI**（bioconda/GitHub only），pip 同名 severus 為冒牌 7.5kB 套件勿用 |
| Tumor BAM | `/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam`（292 GB）|
| Normal BAM | 同目錄 `HCC1395BL.bam`（149 GB）|
| Reference | `/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta`（chr 長度與 BAM @SQ 一致，已驗）|
| het SNP 源 | SAVANA 內建 `hg38_g1000_biallelic_AF0.25-0.75.vcf.gz`（`-vg 1000g_hg38`；無 chr 前綴，SAVANA `allele_counter.py` 偵測 ref 自動補 `chr`）|
| 資源 | 48 cores / 503 GB RAM；TMPDIR→bip7（4.1T，避開 big7 98%）|

---

## 2. chr2+chr8 執行結果（L1）

| 步驟 | 耗時 | 輸出 |
|---|---|---|
| `savana run`（SV）| **14 min** | 12,886 breakpoints（chr2 7874 / chr8 5012）；`HCC1395.sv_breakpoints_read_support.tsv`（14.6MB，3 欄 `VARIANT_ID｜TUMOUR_SUPPORTING_READS｜NORMAL_SUPPORTING_READS` = read-level SV 標籤）|
| `savana cna`（CN）| **49 min** | `*_segmented_absolute_copy_number.tsv`（111 段）；`*_fitted_purity_ploidy.tsv` |

- **purity = 0.96, ploidy = 2.68**（rank 1, distance 0.147）。purity 對純細胞系命中；ploidy 2.68 與 HCC1395 hyper-diploid 相符（僅 2 染色體擬合，全基因組會微調）。
- CN segment 欄位：`chromosome start end segment_id bin_count ... copyNumber minorAlleleCopyNumber meanBAF no_hetSNPs`。chr8 total CN 達 19.3（大段擴增）。
- ⚠ 此 cna 僅用 chr2+chr8 擬合 purity/ploidy（不權威，全基因組才正規）；但 LOH/CN 結構已與 truth 高度吻合（見 §3）。allele counting 是 cna 主要時間瓶頸、purity/ploidy 擬合單緒。

---

## 3. 驗證 vs SEQC2 truth（L1，10kb bin，限 SAVANA covered 區）

> 第一手比對 `HCC1395_segmented_absolute_copy_number.tsv` ↔ `/big8_disk/.../SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed`

| 指標 | chr8 | chr2 |
|---|---|---|
| LOH：SAVANA(minorCN<0.5) vs SEQC2 loh | **98.4% vs 96.3%** | 54.5% vs 50.8% |
| **LOH Jaccard(SAVANA,SEQC2)** | **0.977** | **0.920** |
| LOH per-bin 一致率 | 97.8% | 95.6% |
| gain：SAVANA(CN>2.5) vs SEQC2 gain | 66.4% vs 64.3% | 53.2% vs 57.7% |
| 非-LOH 段（minor≥0.5）| 3/50 | 32/61 |

- **chr8 大段 LOH 如預期重現**（SEQC2 96% ↔ SAVANA 98.4%，Jaccard 0.977）。
- chr2 約半 LOH（與 truth 50.8% 吻合），證明非全染色體假性 LOH（早期前幾段 minor≈0 的疑慮已排除）。
- gain 比例兩染色體皆接近（ploidy-2.68 aneuploid baseline 下 gain/loss 閾值用 2 為基準僅供方向參考；LOH 指標最穩健且最一致）。

---

## 4. 整合框架：k_ISM vs k_CN（用戶提案 — 非循環的 CN 輔助用法）

> **定位修正**：先前 refute 的是「CN↔甲基對齊→證明 subclone」（循環，LOH-unmask 製造對齊）。本框架**不同**：把全基因組 CN 當**固定觀察參考層 / 期望基線**，比較「ISM 切的群數 k_ISM」與「該區 CN 狀態數 k_CN」。**不宣稱對齊=subclone，只當比較尺規與分流依據**（同專案 excess-over-null 邏輯，非循環）。

| 關係 | 解讀 | subclone 意涵 |
|---|---|---|
| **k_ISM < k_CN** | 甲基未解析所有 CN 預測狀態 | 甲基訊號**受 CN 上界所限** → 支持論文 A-framing「甲基化有界（bounded characterize）」|
| **k_ISM ≈ k_CN** | 甲基群數對上 CN 狀態數 | 最可能 **CN-mirror**（W2/I2）→ 不加 subclone 資訊但內部一致 |
| **k_ISM > k_CN** | 🎯 甲基解析超過 CN 能解釋 | **candidate 甲基-subclone 訊號**（excess-over-CN）→ 才進一步 normal cis-control + multi-sSNV CCF 查真偽 |

**價值**：用 CN 把全基因組各區**分流成「不用再問」（k≤k_CN）vs「值得深究」（k>k_CN）**，指導甲基分群該怎麼詮釋。**補既有死角**：無監督「抽象幾群」無乾淨解（[[project_subcluster_cluster_count_determination]]），CN 給 k 一個外部錨。

**4 caveat**：①k_ISM>k_CN ≠ 證明 subclone（candidate，也可能 ASM/技術/read-內甲基相關）②k_CN 定義須講死（建議 allele/haplotype-specific 狀態數）③CN 是估計非真值（purity/ploidy 依賴）= 正是「非解答輔助資訊」④k_ISM 取 a-priori-conditioned 群數（reclassify-v2/SubcloneDbeta @ `feat/summary-nreadsvalid`），非不可靠的無監督 k。

---

## 5. 全基因組固定參考層（已完成 + het-site bug 診斷 + 正解，2026-06-21）

> 全基因組過程揭露一個重要方法學 bug 並修正，最終得到**可信已驗證的固定 CN 參考層**。

### 5a. 第一版 `-vg 1000g`（job btr4kbdrw, cna 345min+SV 100min）→ 🔴 over-call LOH
- auto-fit **purity 0.76 / ploidy 1.83**；SV genome-wide **139,334 breakpoints**。
- 全基因組 LOH **SAVANA 86% vs SEQC2 53%, Jaccard 僅 0.616**（chr7 95% vs 16%、chr21 100% vs 2%）；chr8 仍 0.973（本就真 LOH）。
- 診斷根因二重（L1）：① `--min_cellularity 0.90` re-fit（job bb2a2yj3k）SAVANA **主動拒絕** 0.96/2.79（log: `"NOT an acceptable solution, main CN change step=2"`）堅持 0.76/1.83 → 約束無效；② **`-vg 1000g` 未用 normal 過濾** → allele_counts 混入 germline-homozygous 位點（實例 `chr1:770988 A=0/G=54 BAF=1.0`）→ BAF 膨脹既 over-call LOH 又把 purity 拉低。

### 5b. 正解 `-v` matched-normal het SNP（job bpci11g00, 320min）→ ✅ 決定性成功
- het-site 源 = `germline_phased_merged.vcf.gz`（Clair3 對 **HCC1395BL（normal）** 跑，**1,698,116 筆 100% het**，genome-wide）。
- 結果：**purity 0.96 / ploidy 2.79**（細胞系 purity 命中 + near-triploid 符合 karyotype）。

| 指標（genome-wide pooled, chr1-22, 2801Mb）| `-vg` 汙染版 | **`-v` 正解版** |
|---|---|---|
| purity / ploidy | 0.76 / 1.83 | **0.96 / 2.79** |
| **LOH Jaccard** | 0.616 | **0.962** |
| LOH per-bin 一致 | 67.0% | **97.9%** |
| gain Jaccard | 0.383 | **0.814** |
| 代表 per-chrom LOH Jaccard | chr7 0.171 | **chr7 0.998 / chr1 0.964 / chr8 0.987 / chr5 0.992 / chr12 0.996** |

- **canonical 固定 CN 參考層** = `/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna_normalhet/HCC1395_segmented_absolute_copy_number.tsv`（713 段，欄 `copyNumber`+`minorAlleleCopyNumber`+`meanBAF`）。
- 🔴 **方法定論**：SAVANA tumor-normal **必用 `-v` matched-normal het VCF，不可用 `-vg 1000g`**（後者未過濾 germline-hom 位點 → 汙染）。教訓：固定 CN 參考層必須先 vs SEQC2 驗證才能用於 k 比較（over-call LOH 會把 k_CN 灌成 1）。

---

## 6. 下一步

1. 全基因組 cna 完成 → 取 segment-level total/minor/haplotype CN。
2. 定位 ISM 每區 **k_ISM**（cluster count）來源欄位（significance_summary.csv 或 reclassify-v2 輸出；須對齊 `feat/summary-nreadsvalid` binary）。
3. 全基因組 **k_ISM vs k_CN** 比較 → 每區 regime 標記（<, ≈, >）→ 交甲基分群詮釋。
4. （並行軸）SV `read_support.tsv` 的 tumour-supporting reads → 「攜帶 SV/不攜帶」二分群 → 餵 ISM 甲基 PERMANOVA（非循環 anchor）。
5. Wakhan 移獨立 env，跑 haplotype-specific CN（HP1/HP2 分離）補 SAVANA。

---

## 7. Provenance / caveat

- **L1（第一手重現）**：purity 0.96/ploidy 2.68、SV 12886、chr8 LOH 98.4% vs 96.3% Jaccard 0.977、chr2 54.5% vs 50.8% Jaccard 0.920 — 皆本機 SAVANA 輸出 + SEQC2 truth `awk`/python 直算。
- **caveat**：①chr2+chr8 purity/ploidy 非全基因組擬合（不權威，全基因組進行中）②gain/loss 方向受 aneuploid baseline 影響，LOH 指標最穩健 ③Wakhan 降了 cnvtools env 的 numpy/scipy（SAVANA 仍正常，但 Wakhan 應移獨立 env）④SEQC2 為 segment-level 定性 truth（imprecise breakpoint），只驗方向/區段非 per-base integer CN（見 plan doc）。
