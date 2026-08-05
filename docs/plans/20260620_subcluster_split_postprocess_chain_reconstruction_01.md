# Subcluster「切割總帳」後處理鏈 provenance 重建（HCC1395 → 可參數化）

> **建立**：2026-06-20 · **scope**：read-only provenance 偵察（不改檔/不切 branch）
> **目標物**：`docs/methodology/_assets/20260618_subcluster_pilot/method_comparison.json`（含 `fisher_v.cansplit=5997`=A、`venn_clean.C=8833`=C、`Jaccard=0.118`）
> **腳本所在 branch**：`feat/summary-nreadsvalid`（當前 working tree=`research/subclonal-reconstruction-202606`，scripts 子目錄不在當前 branch；全程用 `git show feat/summary-nreadsvalid:<path>` 讀回）
> **binary provenance**（README）：`feat/summary-nreadsvalid @ 5c39051`；tumor=ClairS_pileup_v040、normal=5khz_simplex（分析端 filter 至 tumor）
> **所有路徑 hardcode 到 worktree** `WT=/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra`（套到他樣本必參數化，見 §5）

---

## 0. 重要釐清 — 「兩個 compute 根」+ 「k-sweep 是分岔支鏈」

`method_comparison.json` 由 **`method_comparison_data.py`** 產生，它讀 **2 個獨立的 ISM-binary compute 根**：

| 欄位 | 含義 | 來源檔 | 真正 compute 腳本 |
|------|------|--------|------------------|
| **A** = `fisher_v.cansplit` (5997) | UPGMA best_k≥2（切得出≥2 平衡群） | `records_wg2.json`（gitignored）| **`wg_contingency.py`**（跑 ISM binary 22 chr） |
| **C** = `venn_clean.C` (8833) | 任一軸 PERMANOVA p<0.05 且非 dispersion-warn | `output/_wg_bdcprime_verify/significance_summary.csv` | **ISM C++ 本體**（pre-existing 外部輸出，無 Python 腳本產生它）|
| `fourgroup` | Q2 2組vs4組 PERMANOVA | `permanova_clean_4group.json` | `permanova_clean_and_4group.py`（讀 `_thresh_cal_2122/` distance matrices）|

> ⚠ **關鍵更正 vs 任務初始假設**：`records_wg2.json` 的真正產生者是 **`wg_contingency.py`**（line 109/111），**不是** `ksweep_kprofile.py`、**不是** `ksweep_wg_records.json` 的別名。
> `ksweep_wg.py` / `ksweep_wg_v2.py` / `ksweep_kprofile.py` / `ksweep_analyze.py` 是**另一條獨立支鏈**（silhouette-k vs alignment-k 對齊式 k-selection 分析），輸出 `ksweep_wg_summary.json` / `kprofile_summary.json`，**完全不餵入 `method_comparison.json`**。任務 prompt 把「method_comparison_data.py 讀 records_wg2.json 的 best_k」與「ksweep_wg_records.json」混為一鏈，實測為兩鏈。

### records_wg2.json 與 ksweep_wg_records.json 的 schema 差異（佐證兩鏈不同）
- `records_wg2.json`（wg_contingency）：每筆 = `{chrom, pos, all:{n,n_complete,best_k,best_sil,contingency,labels}, within:{1-1,2-1,1,2}}` — best_k 是「all-read UPGMA 掃 k=2..min(6,n//3) 取 max silhouette」的單一最佳 k。
- `ksweep_wg_records.json`（ksweep_wg）：每筆 = `{chrom, pos, n, n_k, per_k:[{k,sil,V_hp,p_hp,e_hp,V_carrier,...,V_allele,...}]}` — 保留每個 k 的 sil + 三軸 Cramér's V（為後續 align-k 分析）。

---

## 1. 切割總帳 scripts 目錄全清單（33 .py）

`git ls-tree feat/summary-nreadsvalid --name-only docs/methodology/_assets/20260618_subcluster_pilot/scripts/`

**compute（跑 ISM binary）**：`wg_contingency.py`、`ksweep_wg.py`、`ksweep_wg_v2.py`、`ksweep_pilot.py`、`threshold_calibration.py`、`subcluster_pilot.py`、`subcluster_cases.py`、`precondition_coverage.py`(部分)
**零-compute 後處理（讀 records/sig）**：`method_comparison_data.py`、`permanova_clean_and_4group.py`、`ksweep_kprofile.py`、`ksweep_analyze.py`、`classify_contingency.py`、`classify_nosignal.py`、`cantsplit_validation.py`、`section8_data.py`、`tumor_vs_paired.py`、`sil_lowthresh.py`、`threshold_justify.py`、`wg_contingency.py`(同上)、`select_cantsplit_examples.py`、`select_kprofile_examples.py`、`ksweep_fdr_dump.py`、`ksweep_analyze.py`
**HTML/圖 builder**：`build_method_comparison.py`、`build_correspondence_html.py`、`build_spectrum_html.py`、`build_validation_html.py`、`build_precondition.py`、`build_cantsplit_validation.py`、`build_fdr_workstation.py`、`build_ksweep_explainer.py`、`build_kprofile_explainer.py`、`plot_cantsplit_heatmaps.py`、`plot_fdr_heatmaps.py`、`plot_kprofile_heatmaps.py`

---

## 2. 完整可執行鏈（ISM significance_summary.csv + distance matrices → method_comparison.json）

### 步驟 0（前置 — 外部 ISM 輸出，非本鏈腳本）
ISM C++ 本體已先跑出 **`output/_wg_bdcprime_verify/significance_summary.csv`**（全 WG paired，含 `LabelHPPermanovaP/Valid/DispersionWarn`、`LabelAllelePermanovaP/Valid/DispersionWarn`、`GermlineAsmDbeta` 等欄）。
- **無任何 Python 腳本產生此檔**（grep `_bdcprime_verify.*-o` / `subprocess.*bdcprime` 全空）→ 它是 pre-existing 外部依賴。套到他樣本必須先各自跑 ISM 產生對應 significance_summary.csv。

### 步驟 1 — wg_contingency.py（產 A 的根；跑 ISM binary）
- **腳本**：`scripts/wg_contingency.py`
- **輸入**：ISM `BIN`（per chr 即時跑）；in: TUMOR/NORMAL BAM + REF + `VCFDIR/filtered_snv_tp_chrN.vcf.gz`
- **每 chr 即時產**：`output/_wgc_tmp_chrN/**/distance/BERNOULLI/matrix.csv` + `.../reads/reads.tsv`（跑完即 `rmtree` 刪除，disk-safe）
- **輸出**：`records_wg2.json`（gitignored，~18MB）
- **演算法**：per region → tumor-only reads（`is_tumor` true）× `hp∈{1,1-1,2,2-1}` → peel NaN/負值 → UPGMA(average) 掃 k=2..min(6,n//3) → 取 max silhouette 為 `best_k/best_sil` + cluster×label 列聯表；另對 4 標籤各自做 within-label 子分群。
- **命令**：`python3 scripts/wg_contingency.py`
- **耗時**（README）：~31min

### 步驟 2 — threshold_calibration.py（產 Q2 4-group 用的 distance matrices；跑 ISM binary，僅 chr21+chr22）
- **腳本**：`scripts/threshold_calibration.py`
- **輸出 distance matrices**：`output/_thresh_cal_2122/**/distance/BERNOULLI/matrix.csv`（**保留不刪**，供步驟 3 的 permanova_clean 讀）+ `threshold_calibration.json`
- **命令**：`python3 scripts/threshold_calibration.py`
- 註：此步是 chr21+chr22 pilot（非全 WG），只為 Q2 的 Python-PERMANOVA 2組vs4組示範提供 matrix 來源。

### 步驟 3 — permanova_clean_and_4group.py（產 fourgroup；零 ISM，重算 PERMANOVA）
- **腳本**：`scripts/permanova_clean_and_4group.py`
- **輸入**：`records_wg2.json`（取 cansplit）+ `SIG`(=`_wg_bdcprime_verify/significance_summary.csv`，取 clean) + `output/_thresh_cal_2122/**/matrix.csv`（Q2 Python-PERMANOVA）
- **輸出**：`permanova_clean_4group.json`（其中 `Q2` 被步驟 4 取用）
- **命令**：`python3 scripts/permanova_clean_and_4group.py`

### 步驟 4 — method_comparison_data.py（最終彙整 → method_comparison.json）
- **腳本**：`scripts/method_comparison_data.py`
- **輸入**：
  - `records_wg2.json`（A = `r['all']['best_k']>=2` → cansplit=5997）
  - `SIG`=`_wg_bdcprime_verify/significance_summary.csv`（permanova_raw / dispersion / venn / Δβ-by-class / C=clean）
  - `permanova_clean_4group.json` 的 `["Q2"]`（fourgroup）
- **輸出**：`method_comparison.json`（N=30490；A=5997；C=8833；AnC=1560；Jaccard=0.118；A_only=4437；C_only=7273）
- **命令**：`python3 scripts/method_comparison_data.py`

### 步驟 5（可選）— build_method_comparison.py（HTML，§13-A 注入）
- **輸入**：`method_comparison.json` + `figs_cantsplit/*.png` + `figs_fdr/*.png`（base64 嵌入）
- **輸出**：`20260620_method_comparison_fisher_v_vs_permanova_01.standalone.html`
- **命令**：`python3 scripts/build_method_comparison.py`

### 鏈摘要圖（純文字）
```
[ISM C++ 本體, 外部] ─► output/_wg_bdcprime_verify/significance_summary.csv ──┐
                                                                              │
wg_contingency.py ─(跑ISM 22chr,即時distance+reads)─► records_wg2.json ───────┤
                                                                              ├─► permanova_clean_and_4group.py ─► permanova_clean_4group.json(Q2)
threshold_calibration.py ─(跑ISM chr21+22)─► _thresh_cal_2122/**/matrix.csv ──┘                                            │
                                                                                                                          ▼
records_wg2.json + significance_summary.csv + permanova_clean_4group.json[Q2] ─► method_comparison_data.py ─► method_comparison.json ─► build_method_comparison.py ─► HTML
```

---

## 2b. k-sweep 支鏈（旁系，**不**餵 method_comparison.json — 列出供完整性）

| 步驟 | 腳本 | 輸入 | 輸出 | 命令 |
|------|------|------|------|------|
| ks1 | `ksweep_wg.py`（跑 ISM 22chr, tumor-only）| BAM+REF+VCF | `ksweep_wg_records.json`(gitignored) | `python3 scripts/ksweep_wg.py` |
| ks2 | `ksweep_wg_v2.py`（跑 ISM 22chr, tumor+merged dual-mode）| BAM+REF+VCF | `ksweep_records_tumor.json` + `ksweep_records_merged.json`(gitignored) | `python3 scripts/ksweep_wg_v2.py` |
| ks3 | `ksweep_analyze.py`（零 compute, gate+FDR）| `ksweep_wg_records.json` | `ksweep_wg_summary.json` | `python3 scripts/ksweep_analyze.py [DELTA] [FDR] [INP] [OUTP]` |
| ks4 | `ksweep_kprofile.py`（零 compute, margin+三態）| `ksweep_wg_records.json`(tumor) + `ksweep_records_merged.json` | `kprofile_summary.json` + `kprofile_loci_{tumor,merged}.json` | `python3 scripts/ksweep_kprofile.py` |

> 註：`ksweep_kprofile.py` 內部自算 `best_k=max(per_k by sil)`，與 `wg_contingency` 的 best_k 是**不同 compute 路徑**（ksweep best_k 來自 per_k silhouette；wg_contingency best_k 來自獨立 cluster_bestk）。兩者不可互換餵 method_comparison。

---

## 3. records_wg2.json 產生者（明確結論）

**`wg_contingency.py`**（`scripts/wg_contingency.py:109` 逐 chr dump + `:111` 最終 dump）。
- 它**不是** kprofile 的輸出（kprofile 讀 ksweep_wg_records.json，輸出 kprofile_*.json）。
- 它**不是** ksweep_wg_records.json 的別名（schema 不同，見 §0）。
- 它**是另一獨立腳本**：直接跑 ISM binary 22 chr，per region 算 all-read UPGMA best_k + 列聯表 + within-label 子分群。
- gitignored（README：~18MB 中間檔，從 binary 5c39051 重生）。

---

## 4. compute 量級估計

| 腳本 | 跑 ISM？ | per-region 讀 distance CSV？ | scope | 量級 |
|------|:---:|:---:|------|------|
| `wg_contingency.py` | ✅ 22 chr | ✅（即時產→讀→刪） | 全 WG 30,490 region | **~31min**（README 實測）≈ 與 1 次全 WG ISM scan 同級 + Python UPGMA/silhouette overhead；I/O 重（每 chr 寫滿 distance/reads 再 glob 讀再 rmtree）|
| `ksweep_wg.py` | ✅ 22 chr | ✅ 同上 | 全 WG | 同 ~ISM scan 級（多掃每 k 的三軸 Cramér's V，CPU 略高於 wg_contingency）|
| `ksweep_wg_v2.py` | ✅ 22 chr | ✅ 同上 | 全 WG ×2 read-set | 略高於 ks1（同一 binary pass 但 tumor+merged 雙分群）|
| `threshold_calibration.py` | ✅ 2 chr | ✅ 保留不刪 | chr21+chr22 pilot | 遠小於全 WG（~2/22 region + B=30 permutation/locus）|
| `method_comparison_data.py` | ❌ | ❌（讀 json+csv）| — | 秒級（純讀 records_wg2.json + 一次 CSV pass）|
| `permanova_clean_and_4group.py` | ❌ | ✅（讀 _thresh_cal matrix）| chr21+22 matrices | Q1 秒級；Q2 對 chr21+22 region 各跑 2組+4組 B=199 permutation PERMANOVA → 分鐘級（受 region 數×199 影響）|

**結論**：compute 重心是「跑 ISM binary」的 3 個全 WG 腳本（`wg_contingency` / `ksweep_wg` / `ksweep_wg_v2`），各 ~ISM 全 WG scan 同級（~30min）+ Python 後處理 overhead。I/O 模式皆為 per-chr 即時產 distance matrix → glob 讀 → 立即 rmtree（disk-safe 但反覆讀寫）。最終彙整步驟（`method_comparison_data.py`）秒級、零 ISM。

---

## 5. Hardcode 路徑清單（套到他樣本必參數化的點）

> 所有 compute/後處理腳本均 hardcode 到 worktree `ism-review-infra` + HCC1395 BAM/VCF/REF + `_wg_bdcprime_verify`。下表為 method_comparison.json 鏈 + 其支撐腳本的關鍵 hardcode（檔:行:變數）。

### 5.1 worktree / asset 根（每個腳本都有）
| 檔 | 行 | 變數 | 值 |
|----|----|------|----|
| `wg_contingency.py` | 12 | `WT` | `/big7_disk/.../worktrees/ism-review-infra` |
| `wg_contingency.py` | 13 | `BIN` | `{WT}/build/bin/inter_sub_mod` |
| `wg_contingency.py` | 18 | `ASSET` | `{WT}/docs/methodology/_assets/20260618_subcluster_pilot` |
| `method_comparison_data.py` | 5 | `WT` | 同上 |
| `method_comparison_data.py` | 6 | `A` | `{WT}/docs/.../20260618_subcluster_pilot` |
| `permanova_clean_and_4group.py` | 5 | `WT` | 同上 |
| `permanova_clean_and_4group.py` | 6 | `A` | 同上 |
| `ksweep_wg.py` / `ksweep_wg_v2.py` | 14 / 13 | `WT`/`BIN` | 同上 |
| `ksweep_kprofile.py` | (頂) | `A` | hardcode 全路徑（無 WT 變數，整串寫死）|
| `ksweep_analyze.py` | (頂) | `ASSET` | hardcode 全路徑 |
| `threshold_calibration.py` | 15-16 | `WT`/`OUT`/`A` | `OUT={WT}/output/_thresh_cal_2122` |

### 5.2 樣本專屬 BAM/VCF/REF（HCC1395 寫死 — 換樣本必改）
| 檔 | 行 | 變數 | 值 |
|----|----|------|----|
| `wg_contingency.py` | 14 | `TUMOR` | `/big8_disk/.../data/bam/HCC1395/tumor.bam` |
| `wg_contingency.py` | 15 | `NORMAL` | `/big8_disk/.../data/bam/HCC1395/normal.bam` |
| `wg_contingency.py` | 16 | `REF` | `/big8_disk/.../data/ref/hg38.fa` |
| `wg_contingency.py` | 17 | `VCFDIR` | `/big8_disk/.../data/vcf/HCC1395/pileup` |
| `wg_contingency.py` | 70 | (inline) | `filtered_snv_tp_{chrom}.vcf.gz` 命名 pattern |
| `ksweep_wg.py` / `ksweep_wg_v2.py` | 15-18 / 14-17 | `TUMOR/NORMAL/REF/VCFDIR` | 同 HCC1395（完全重複）|
| `threshold_calibration.py` | (內) | BAM/VCF/REF | 同 HCC1395 |

### 5.3 外部 ISM 輸出（pre-existing，整鏈共用 — 換樣本必各自重產）
| 檔 | 行 | 變數 | 值 |
|----|----|------|----|
| `method_comparison_data.py` | 8 | `SIG` | `{WT}/output/_wg_bdcprime_verify/significance_summary.csv` |
| `permanova_clean_and_4group.py` | 8 | `SIG` | 同上 |
| `precondition_coverage.py` | 30 | `SIG` | 同上 |
| `cantsplit_validation.py` | 9 / `classify_contingency.py` 9 / `classify_nosignal.py` 8 / `section8_data.py` 9 / `sil_lowthresh.py` 11 / `threshold_justify.py` 11 / `tumor_vs_paired.py` 8 / `select_cantsplit_examples.py` 7 | (各) | `SIG`/`SUM` 皆指向 `_wg_bdcprime_verify/significance_summary.csv`（HCC1395 paired，無 sample 變數）|

### 5.4 中間檔名（records_wg2.json 等 — 跨樣本需加 sample 前綴避免覆蓋）
| 檔 | 行 | 變數/檔名 |
|----|----|----------|
| `wg_contingency.py` | 109,111 | `{ASSET}/records_wg2.json`（寫死，無 sample 區分）|
| `method_comparison_data.py` | 9 | 讀 `{A}/records_wg2.json` |
| `permanova_clean_and_4group.py` | 9 | 讀 `{A}/records_wg2.json` |
| `permanova_clean_and_4group.py` | 78 | glob `_thresh_cal_2122/**/matrix.csv`（dir 名寫死）|

### 5.5 其他寫死常數（演算法參數，跨樣本通常保留但須知曉）
- `wg_contingency.py:21` `MIN_SZ=3`；`:22` `LABMAP`（HP 標籤映射）
- `wg_contingency.py` UPGMA k 範圍：`range(2,min(6,n//MIN_SZ)+1)`（best_k 掃 k=2..5）
- ISM 參數寫死：`-w 5000`、`-j 16`、`--distance-metric BERNOULLI`、`--nan-distance-strategy SKIP`（每個 compute 腳本的 subprocess call）
- `os.environ["TMPDIR"]="/big7_disk/liaoyoyo2001/tmp"`（每個 compute 腳本）
- best_k≥2 門檻（A 定義）：`method_comparison_data.py:15` `r['all']['best_k']>=2`
- clean 定義：`method_comparison_data.py` p<0.05 & `not DispersionWarn`（per HP/Allele 軸）

---

## 6. 參數化 checklist（套用他樣本最小改動集）
1. `WT` / `ASSET` / `A`：改為目標 worktree 或新 asset dir（建議每樣本獨立 asset dir）。
2. `TUMOR/NORMAL/REF/VCFDIR`：指向目標樣本 BAM/REF/VCF（3 個 compute 腳本各自有一份，須同步改）。
3. **先各自跑 ISM 產 `<sample>/significance_summary.csv`**（_wg_bdcprime_verify 是 HCC1395 pre-existing，他樣本無），再把所有 `SIG=` 指向該樣本檔。
4. `records_wg2.json` 等中間檔名加 sample 前綴（否則跨樣本互相覆蓋）。
5. VCF 命名 pattern `filtered_snv_tp_{chrom}.vcf.gz` 須與目標樣本 VCF 命名一致。
6. 演算法常數（MIN_SZ / k 範圍 / -w 5000 / BERNOULLI / SKIP）通常保持與 HCC1395 一致以保可比性。

---

## 7. UNKNOWN / blocker
- **`_wg_bdcprime_verify/significance_summary.csv` 的確切產生命令**：UNKNOWN（無 Python 腳本產生它；README 只說 binary=5c39051、tumor=ClairS_pileup_v040 / normal=5khz_simplex；ISM 本體呼叫的完整 flag 組合未在此 asset dir 記錄）。套他樣本前需從 ISM run script（`scripts/run_*.sh` 或 ism-review-infra worktree 的 run log）回溯該 verify-run 的完整命令。
- **records_wg2.json / significance_summary.csv 的實體檔**：gitignored / 在 worktree 內，未在當前 branch；本次偵察未實際 Read 其磁碟內容（只讀 committed `method_comparison.json` 確認下游數字）。重跑前須確認 worktree `ism-review-infra` 的 `output/_wg_bdcprime_verify/` 仍存在或可重產。
- **clean 定義的兩腳本微差**：`method_comparison_data.py`（C=8833）的 clean 用「(hsig & !hw) OR (asig & !aw)」逐軸；`permanova_clean_and_4group.py` 的 clean 用「hp_clean OR al_clean」同義 — 兩者數字應一致，但若日後一邊改門檻會 drift，套用時以 `method_comparison_data.py` 為 SoT。
