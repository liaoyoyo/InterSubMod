<!--
建立時間: 2026-06-21
類型: 資產清單 + 敘述索引（navigation hub / asset manifest）
狀態: in_progress（全基因組固定參考層執行中）
主軸: Subclonal reconstruction — CN/SV 輔助整合層
build_branch: research/subclonal-reconstruction-202606
data_sources: /big7_disk/liaoyoyo2001/cnv_sv_work/ (ls 盤點 2026-06-21), /big8_disk/data/HCC1395/, /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools
-->

# CNV/SV 工作資產清單 + 敘述索引（HCC1395, SAVANA/Wakhan/Severus）

> **這是導航中樞**：把本主題所有重要檔案（repo 文件 + repo 外工作產出 + env + memory）串成一張表，附敘述脈絡。檔案分散在 repo / `big7_disk_output` 外的 `cnv_sv_work/` / `big8_disk` 唯讀資料 / `bip7` conda env，故需此索引。
> 路徑前綴：repo 內檔 `InterSubMod/...`；repo 外工作檔給絕對路徑。

---

## 0. 敘述脈絡（一頁讀懂這條線在做什麼）

**問題**（用戶）：ONT 上最好的 CNV 偵測軟體是哪個？能否定位每位置 CN 數量？能否用 CN（+SV 斷點）當標籤，去確認 ISM 甲基 read-clustering 的 subclone 分群、與標籤關係，並用 SEQC2 CNV 解答驗證？

**研究裁決**（workflow wf_5497277c-620, 10 agents）：
- 工具 = **SAVANA（務實 PRIMARY）+ Wakhan（haplotype-specific 互補）+ Severus（SV）**。
- per-segment integer/allele/haplotype CN ✅；per-read CN / read→subclone 指派 = **無任何工具**（segment-level 鴻溝）。
- 🔴 「**CN 分群↔甲基分群對齊→證明 subclone**」= **REFUTED**（循環：W2/I2 LOH-unmask 製造對齊）。
- SEQC2 只能驗 **gain/loss/loh 方向**，不能驗 per-base/allele/read-level integer CN。

**用戶 reframing（採納，非循環）**：CN 不當 validator，而當**固定觀察參考層**比較 **k_ISM vs k_CN** → k<k_CN=甲基有界 / k≈=CN-mirror / **k>k_CN=candidate subclone**。用 CN 分流各區甲基詮釋，補無監督 k 無外部錨的死角。

**執行驗證**（2026-06-20）：SAVANA 在 ONT HCC1395 跑通，chr2+chr8 vs SEQC2 → **chr8 LOH Jaccard 0.977 / chr2 0.920, purity 0.96** = 可行性 CONFIRMED → 啟動**全基因組固定 CN 參考層**（過夜，執行中）。

---

## 1. Repo 文件索引（敘述 + 結論）

| 檔案（InterSubMod/...）| 角色 | tier |
|---|---|---|
| `docs/plans/20260620_ont_cnv_sv_subclone_verification_feasibility_01.md` | **研究裁決 plan**：工具選型 + per-locus CN + confound + 7 子問題 GO/PROBE/NO-GO + SEQC2 路徑 | L1/L2/L3 標註 |
| `docs/experiments/in_progress/2026/06/20260620_savana_hcc1395_cnv_sv_feasibility_result_01.md` | **執行結果**：chr2+chr8 vs SEQC2 驗證表 + k_ISM vs k_CN 整合框架 | L1 數據 |
| `docs/experiments/in_progress/2026/06/20260620_cnv_sv_work_assets_manifest_01.md` | **本檔**（資產清單 + 敘述索引）| — |
| `docs/experiments/INDEX.md`（最新區頂條目）| 實驗 SoT 索引條目 | — |
| memory `project_ont_cnv_sv_subclone_verification_feasibility` | 跨 session 記憶（裁決 + 執行 + k 框架）| — |

---

## 2. 工作產出索引（repo 外 `/big7_disk/liaoyoyo2001/cnv_sv_work/`）

### 2a. chr2+chr8 pilot（已完成，可行性證據）
| 檔案 | 大小 | 內容 |
|---|---|---|
| `HCC1395/savana_chr2_chr8/run/HCC1395.sv_breakpoints.vcf` | 10 MB | 12,886 SV breakpoints（chr2 7874/chr8 5012；raw 未分類 somatic/germline）|
| `HCC1395/savana_chr2_chr8/run/HCC1395.sv_breakpoints_read_support.tsv` | 14.6 MB | **read-level SV 標籤**（`VARIANT_ID｜TUMOUR_SUPPORTING_READS｜NORMAL_SUPPORTING_READS`）|
| `HCC1395/savana_chr2_chr8/run/HCC1395.sv_breakpoints.bedpe` | 718 KB | SV bedpe |
| `HCC1395/savana_chr2_chr8/cna/HCC1395_segmented_absolute_copy_number.tsv` | 9.7 KB | **integer/allele CN per segment**（111 段；`copyNumber`+`minorAlleleCopyNumber`+`meanBAF`）|
| `HCC1395/savana_chr2_chr8/cna/HCC1395_fitted_purity_ploidy.tsv` | 46 B | **purity 0.96 / ploidy 2.68**（rank1）|
| `HCC1395/savana_chr2_chr8/cna/HCC1395_allele_counts_hetSNPs.bed` | 20 MB | het SNP allele counts |
| `HCC1395/savana_chr2_chr8/cna/HCC1395_read_counts_mnorm_log2r_segmented.tsv` | 4.2 MB | log2 ratio 分段 |

### 2b. 全基因組固定參考層（✅ 完成 2026-06-21）
| 路徑 | 狀態 |
|---|---|
| **`HCC1395/savana_wgs/cna_normalhet/HCC1395_segmented_absolute_copy_number.tsv`** | ✅ **canonical 固定 CN 參考層**（`-v` normal-het；purity 0.96/ploidy 2.79；713 段；全基因組 LOH Jaccard vs SEQC2 **0.962**）|
| `HCC1395/savana_wgs/run/HCC1395.sv_breakpoints*` | ✅ 全基因組 SV **139,334 breakpoints** + read_support.tsv（163MB）|
| ~~`HCC1395/savana_wgs/cna/`~~（`-vg 1000g`）| 🔴 **棄用**（het-site 汙染 over-call LOH，Jaccard 僅 0.616；診斷見結果 doc §5a）|
| ~~`HCC1395/savana_wgs/cna_purity1/`~~ | 🔴 棄用（cellularity 約束無效，SAVANA 仍回 0.76/1.83）|

### 2b-2. Wakhan haplotype-specific CN + Severus SV（✅ 完成 2026-06-22，job bhmc9i1nt）
| 路徑 | 狀態 |
|---|---|
| `HCC1395/severus_out/somatic_SVs/severus_somatic.vcf` | ✅ Severus somatic SV **1,789**（+ all_SVs 22,588；`--phasing-vcf germline_phased`）|
| **`HCC1395/wakhan_out/2.77_0.92_0.8/bed_output/`** | ✅ **Wakhan rank-1**（purity 0.92/ploidy 2.77 ≈ SAVANA 0.96/2.79）：`*_segments_HP_1.bed`/`HP_2.bed`(711 段 HP-specific integer CN + svs_breakpoints_ids)+`loh_regions.bed`(129)+subclonal |
| `HCC1395/wakhan_out/1.42_0.96_0.51/` | 🔴 rank-2（near-diploid 替代支，conf 0.51，勿用）|
| 三方比對 | Wakhan↔SEQC2 LOH Jaccard 0.819 / Wakhan↔SAVANA 0.832 / total CN Spearman 0.699；報告 `20260622_wakhan_hcc1395_3way_cn_comparison_01.md` · script `scripts/analysis/compare_wakhan_savana_seqc2.py` |

### 2c. 分析腳本（可重現）
| 檔案 | 用途 |
|---|---|
| `/big7_disk/liaoyoyo2001/cnv_sv_work/compare_savana_seqc2_loh.py` | **SAVANA LOH/CN vs SEQC2 concordance**（產生 chr8 0.977 / chr2 0.920 那組數字；可套全基因組/COLO829）|
| `HCC1395/contigs_chr2_chr8.txt` | SAVANA `--contigs` 染色體限定檔 |

### 2d. 日誌（`cnv_sv_work/logs/`）
`install_libmamba.log`（成功安裝）· `savana_run_chr2_chr8.log`（14min）· `savana_cna_chr2_chr8.log`（49min）· `savana_wgs.log`（全基因組，執行中）· `install_pip.log`/`install_cnvtools.log`（失敗的早期嘗試紀錄）

---

## 3. 環境與工具

| 項目 | 值 |
|---|---|
| conda env | `cnvtools` @ `/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools`（bioconda, `--solver=libmamba`）|
| SAVANA | 1.3.7（subcommands: run/cna/to/classify）|
| Severus | 1.7（somatic SV, `--output-read-ids`）|
| Wakhan | 0.4.3 @ `cnv_sv_work/Wakhan/`（GitHub clone+pip；🔴 降了 env numpy/scipy，**待移獨立 env**；entry_points.txt 在 egg-info）|
| 🔴 安裝陷阱 | SAVANA/Wakhan **不在 PyPI**（bioconda/GitHub only）；pip `severus` 是冒牌 7.5KB 套件；conda classic solver 卡死需 `--solver=libmamba` |

---

## 4. 輸入資料（`/big8_disk` 唯讀）

| 項目 | 路徑 | 大小 |
|---|---|---|
| Tumor BAM | `/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam` | 292 GB |
| Normal BAM | 同目錄 `HCC1395BL.bam` | 149 GB |
| Reference | `/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta` | chr-prefixed（與 BAM @SQ 一致）|
| SEQC2 CNV truth | `/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed` | 660 段定性 gain/loss/loh |
| 內部 CN proxy | `big7_disk_output/canonical/HCC1395/paired_full/20260420_*/intersubmod_tp/significance_summary.csv` | Coverage_Multiple/Potential_LOH/LOH_Bed_Overlap（29,754 位點）|

---

## 5. 執行紀錄（jobs + timing）

| job | 範圍 | 耗時 | 結果 |
|---|---|---|---|
| `savana run` | chr2+chr8 SV | 14 min | 12,886 breakpoints ✓ |
| `savana cna` | chr2+chr8 CN | 49 min | purity 0.96 / chr8 LOH Jaccard 0.977 ✓ |
| `btr4kbdrw`（全基因組）| cna→SV | 執行中（allele-counting ~2h+）| ⏳ 過夜 |

---

## 6. 進行中 + 下一步

1. ⏳ **全基因組固定 CN 參考層完成**（job btr4kbdrw）→ 取 segment-level total/minor/haplotype CN。
2. **定位 ISM 每區 k_ISM 來源欄位**（須對齊 `feat/summary-nreadsvalid` 的 a-priori cluster count，非無監督 k）。
3. **全基因組 k_ISM vs k_CN 比較** → 每區 regime 標記（<, ≈, >）。
4. SV `read_support.tsv` → 「攜帶 SV/不攜帶」二分群 → 餵 ISM 甲基 PERMANOVA（非循環 anchor）。
5. Wakhan 移獨立 env → haplotype-specific CN（HP1/HP2）補 SAVANA。

---

## 7. Provenance / caveat

- **L1（第一手）**：purity 0.96、chr8 LOH 98.4% vs 96.3% Jaccard 0.977、chr2 54.5% vs 50.8% Jaccard 0.920、SV 12,886 — 皆本機 SAVANA 輸出 + SEQC2 truth 經 `compare_savana_seqc2_loh.py` / awk 直算。
- **caveat**：chr2+chr8 purity/ploidy 非全基因組擬合（全基因組進行中）；gain/loss 方向受 aneuploid baseline 影響（LOH 最穩健）；Wakhan 降 env 套件待隔離；SEQC2 = segment-level 定性 truth；「CN↔甲基對齊=subclone」REFUTED（k 比較框架是不同的非循環用法）；單樣本封頂 ⭐3。
