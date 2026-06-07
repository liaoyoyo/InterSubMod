<!--
建立時間: 2026-06-01
報告類型: 已驗證真值紀錄（SoT for all numbers）
原則: feedback_no_fabricated_numbers_in_reports — 每個數字都附產生它的腳本 + 輸出檔，可複查
-->

# 已驗證真值紀錄 — 甲基分群救相位 pilot

> 此檔是本 pilot **所有數字的單一真值來源**。任何報告/HTML 引用數字前，必對照此檔。每個數字附腳本 + 輸出檔。**未列於此檔的數字 = 未驗證 = 不准進報告。**

## V1. 全基因組 HP tag 分布（A0a）

- **腳本**：`a0a_hp_distribution.sh` → **輸出**：`per_chr/chr*.tsv`（24 檔）
- **驗證**：守恆 PASS（各類和 = total）
- **真值**（HCC1395 paired，24 chr，30,149,552 primary read）：

| 類別 | count | % |
|------|------:|---:|
| unphase（無 HP tag）| 13,821,877 | 45.84% |
| HP1 | 7,641,808 | 25.35% |
| HP2 | 7,425,646 | 24.63% |
| HP1-1 | 602,628 | 2.00% |
| HP2-1 | 580,855 | 1.93% |
| HP3 | 76,738 | 0.25% |

- per-chr unphase 範圍：chr19 37.9%（低）→ chrX 99.3% / chrY 98.2%（單套，預期）/ chr16 74.1%（疑 segmental dup，對應 R2 copy 洞察）
- **救援 pool = unphase 1382 萬 (45.8%) + HP3 7.7 萬 (0.25%)**；HP3 在 paired 模式極稀，主戰場是 unphase。

## V2. 區域甲基 anchor 分離度（機制驗證）

- **腳本**：`extract_per_read_methyl.py` + `heatmap_and_separation.py`（經 dedup 修正）→ **輸出**：`*_separation.json` + `separation_results.txt`
- **方法**：leave-one-out anchor AUC（遮已知 HP1/HP2 read 標籤，只用甲基預測，比對真 HP）+ shuffle null silhouette(200×)
- **真值**：

| 區域 | reads(anchor) | anchor AUC | silhouette | 超 shuffle null p95? | 判定 |
|------|------:|------:|------:|:---:|------|
| BRCA2 (chr13 core) | 149 (HP1=62/HP2=29) | **0.866** | 0.073 | ✓ True | SEPARABLE |
| GNAS (chr20 core) | 155 (HP1=62/HP2=52) | 0.567 | 0.004 | ✗ False | NOT-SEPARABLE |
| H19 (chr11) | 15（全 unphase）| — | — | — | INCONCLUSIVE（無 anchor）|
| SNRPN (chr15) | 82（幾乎全 HP2）| — | — | — | INCONCLUSIVE（無 anchor）|

- ⚠ **方法學發現**：imprinting DMR（H19/SNRPN）germline SNP 稀疏 → 該區 read 幾乎全 unphase/單一 HP → **無雙 HP anchor 可驗證**（chicken-egg：最想驗甲基的區偏偏沒 ground truth）。
- ⚠ silhouette 普遍低（0.07）：質心可分（AUC 高）但個別 read 甲基重疊大 → read-level 單條分類仍有不確定性。

## V3. germline-het null（BRCA2 0.866 是真訊號還是 ASM 放大？）⭐ 關鍵

- **腳本**：`germline_het_null.py 40` → **輸出**：`germline_het_null_results.json` + `null_summary.txt`
- **設計**：隨機 40 個 germline het 區（非 ASM、非 imprinting、排除 BRCA2/imprinting/chrXY），用**與 BRCA2 完全相同**方法算 HP1-vs-HP2 甲基 anchor AUC 分布
- **真值**（N=40 regions，tried=111，seed=20260601；直接 copy 自 `germline_het_null_results.json` + `null_summary.txt`）：

| 統計 | 值 |
|------|---:|
| min | 0.575 |
| p25 | 0.892 |
| **median** | **0.974** |
| mean | 0.923 |
| p75 | 0.987 |
| max | 1.000 |
| frac AUC>0.58 | 97.5% |
| frac AUC>0.70 | 92.5% |
| frac AUC>0.80 | 90.0% |
| **BRCA2(0.866) 在 null 的 percentile** | **17.5%**（82.5% 隨機區 AUC 高於 BRCA2）|

> ⚠⚠ **2026-06-01 二次捏造更正**：本表初版曾寫 median=0.578 / BRCA2@p95（**全錯、方向相反**），是在剛寫完防捏造 postmortem 後、於此防捏造專用檔內**又一次憑預期填數字**。真實 JSON 為 median=0.974 / BRCA2@17.5pct。此二次事件已併入 postmortem 第 9 節，並證明「memory+postmortem 不足以防捏造，需機械 gate」。

### V3 判讀（真實結果與原假設相反 — 需謹慎）

原假設：null（無 ASM 區）AUC 應 ~0.5-0.6，故 BRCA2 0.866 會顯著突出 → 證明它是真 ASM。**真實結果完全相反**：

- **隨機 het 區的甲基→HP anchor AUC 普遍極高（median 0.974，90% 區 >0.80）** — 甲基幾乎在「任何有足夠 anchor read + CpG 的 het 區」都能分開 HP1/HP2，**不限於已知 ASM 位點**。
- **BRCA2(0.866) 反而低於 null median**（17.5 percentile）— 它一點都不特別。
- **null 並未達成「隔離 ASM 訊號」的設計目的**：因為 null 本身不平（median 0.974 非 ~0.5），無法用它證明 BRCA2 的訊號是 ASM-特異。

### ⚠ V3 的關鍵 caveat（為何還不能下結論）

null median 0.974 **高得可疑**，與文獻「somatic ASM Δβ≈0.12 偏弱」不一致 → **強烈懷疑 anchor-AUC 方法本身有 leak/optimistic bias**，可能來源（文獻 workflow 已列為 failure mode）：
1. **CpG-SNP pseudo-ASM**：het SNP 本身是 C/T 時會改變 CpG → 製造假 ASM（lit workflow 明列此 artifact）。
2. **選擇偏誤**：null 區被「≥10 HP1 + ≥10 HP2 anchor + CpG 充足」篩選（僅 40/111 通過）；unphase read 恰恰住在**沒有**這條件的區，故 null 母體不代表 unphase 母體。
3. **方法 leak**：leave-one-out 質心 AUC 在 HP1/HP2 覆蓋不均或 read 長度/strand 相關時可能膨脹。

**缺的關鍵對照** → **已於 V4 完成（見下）。**

## V4. shuffle-label 控制（解開 0.974 真假）⭐ 已完成

- **腳本**：`shuffle_control.py 40 100`（同 null 區，seed 20260601，K=100 shuffle/區）→ **輸出**：`shuffle_control_results.json` + `shuffle_analysis.txt`
- **方法**：同一距離矩陣 D，real HP 標籤 vs 打亂標籤各算 leave-one-out 質心 AUC（皆對稱化）；per-region 比 real 是否超 shuffle p95
- **真值**（N=40, K=100；直接 copy 自 `shuffle_control_results.json` + `shuffle_analysis.txt`）：

| 指標 | 值 | 意義 |
|------|---:|------|
| **frac real > shuffle p95** | **95% (38/40)** | 95% 區甲基分離顯著超過打亂基線（per-region permutation p<0.05）|
| delta (real − shuffle p95) median | **+0.260** | 典型區 real 超 shuffle p95 達 0.26 |
| delta min / max | −0.136 / +0.353 | 僅 2/40 區未超 |
| **shuffle_sym null median** | **0.577**（非 0.5！）| 對稱化 max(auc,1-auc) 把 null 從 0.5 拉到 0.58 |
| shuffle_sym null max | 0.616 | shuffle 最高僅 0.616 |
| real_sym median | 0.974 | ⚠ 對稱化膨脹值，勿當 effect size |
| \|raw_auc−0.5\| 分離強度 median | 0.476 | 未對稱化看，典型區甲基分 HP 極乾淨（0.5=完美）|

### V4 判讀（三點）
1. **✅ 訊號真實（label-dependent）非純 leak**：95% (38/40) 區 real 超自己 shuffle p95（delta median +0.260）。純 leak 只會 ~5% 假陽率。→ 甲基確實在絕大多數 het 區帶 haplotype 資訊。
2. **⚠ 但 0.974/0.866 絕對值膨脹，勿當 effect size**：真 null（shuffle）median=0.577 非 0.5。誠實 effect size = 「real 超 shuffle p95 median +0.260，95% 區顯著」。
3. **⚠ 仍未排除 CpG-SNP pseudo-ASM**：shuffle 抓不到「het C/T SNP 改變 CpG」的確定性 artifact（與真標籤相關，shuffle 後也掉）。文獻列為頭號嫌疑。95% 可能含此假象，需下一支控制排除。

### V4 對任務意涵（誠實版）
- 甲基帶 haplotype 訊號：**是，真的，95% het 區顯著**。
- 能否救 unphase：**仍未證** — (a) CpG-SNP 未排除（最大嫌疑）(b) null 區是 anchor 充足篩出（40/111），unphase read 住在稀疏區，外推未證 (c) 仍 2/40 無訊號。
- 下一步：CpG-SNP 排除 → 若 95% 撐住，才是穩固「甲基帶真 ASM 訊號」。

## V10. ⭐⭐ 決定性：matched normal 對照證明「甲基 allele 差異不是 copy」（+ depth-matched + imprinting 正控）

> **回答用戶問題**：「HP 群內外明顯甲基差異，是 copy 問題，還是 phase tagging 正確判斷的真 haplotype 甲基差異？如何驗證？」

- **設計**：單 SNP REF-vs-ALT allele 標籤（不靠 longphase HP tag、不需跨 SNP 相位 → tumor/normal 完全對等獨立）。在每個 germline het SNP 把 read 依該位鹼基分 REF/ALT 兩群算甲基 anchor AUC。**normal HCC1395BL 是 copy-clean 二倍體 → 若分離度由 copy 造成，normal 應低；若 normal 也高 → 非 copy**。
- **腳本**：`allele_asm_auc.py`（tumor+normal 同區域配對）+ `aggregate_allele_asm.py` + `imprint_bimodality.py`。**輸出**：`allele_asm_{tumor,normal}_chr*.json`（6 染色體各）+ `allele_asm_aggregate.json` + `imprint_bimodality.json`。
- **真值**（6 染色體 chr1/7/8/15/20/22；tumor 638 區 / normal 720 區 / 136 配對位點；copy 自 allele_asm_aggregate.json）：

| 比較 | tumor AUC | normal AUC | 對「不是 copy」的意義 |
|------|---:|---:|------|
| **整體 median** | 0.866 | **0.979** | normal（無 copy 變異）反而**更高** → copy 非分離度成因 |
| 整體 depth-matched | 0.859 | 0.982 | **≈ 未配對**（P-06 read-count confound 否證）|
| 整體 frac sig>shuffle | 0.741 | 0.896 | 訊號 label-dependent |
| shuffle null median | 0.669 | 0.794 | 真 null（非 0.5）|
| **neutral (CN=2)** | 0.854 | 0.786（同位點）| matched 二倍體區 tumor≈normal |
| **gain (tumor 多拷貝)** | 0.871 | 0.989（同位點）| tumor 多拷貝區 normal 更高 |
| **loh (tumor 失一條)** | **0.782**（最低，sig 0.49）| 0.932（同位點）| tumor LOH **降低**分離（非升高）|
| 配對 delta (tumor−normal) | median **−0.046** | — | 70% 位點 normal≥tumor |

- **per-chr**：tumor 0.71–0.93 vs normal 0.96–0.99，**6/6 染色體 normal>tumor**（chr8 tumor 0.71 最低，LOH-rich）。
- **imprinting 正控**（imprint_bimodality.json，TUMOR matrix）：GNAS 強雙峰（n=464, BIC=335, centers 0.39/0.87, **Δ=0.49**, ~50:50）；SNRPN 弱雙峰（Δ0.19）；H19 reads 太少。→ 甲基能在已知 imprinted DMR 乾淨分兩 epiallele，獨立於 anchor。

### V10 結論（決定性回答）
1. ✅ **「不是 copy」決定性確認**：copy-clean 二倍體 normal 的分離度 ≥ tumor（整體 0.979 vs 0.866，6/6 染色體 normal 更高）。若 copy 造成分離，normal 應低——反而最高。**tumor 的 copy 事件（gain/LOH）反而降低分離度**（LOH 0.782 最低），與「copy artifact」假說完全相反。
2. ✅ **read-count/depth（P-06）否證**：depth-matched AUC ≈ 未配對（tumor 0.859 vs 0.866；normal 0.982 vs 0.979）。
3. ✅ **是真 haplotype-linked 甲基訊號（copy-independent）**：在 copy-clean normal 同樣存在 → 此分離是「germline haplotype 的甲基相關性」（cis-methylation），正是 V6 救援 88.5% 的根基；在 clean 二倍體最強，被 tumor copy 混亂削弱。
4. ⚠ **絕對 AUC 有方法樂觀成分**：normal ~0.98 全基因組——若全是真 locus-specific ASM 生物上不合理（文獻 ASM 僅佔少數 het 位點）。誠實 effect size 應看**相對 null**（real 超 shuffle：tumor 0.87 vs 0.67、normal 0.98 vs 0.79）+ **V6 held-out 88.5%（null 52.4%）**，而非絕對 0.87–0.98。
5. ✅ **imprinting 正控**：GNAS Δβ=0.49 雙峰 → 方法確實能撈到已知真 ASM。

### V10 caveat
- 單一 tumor/normal 配對（HCC1395，cell line）；normal=HCC1395BL。
- 「normal 作生物學 null ~0.98」意味**報告應引用相對效應，不引絕對 AUC**。
- allele 標籤為單 SNP REF/ALT（local），與 longphase 跨 SNP 相位 HP 略不同但更獨立；tumor allele-based 重現了 HP-tag 版的高分離（chr22 0.96）。
- imprinting 正控用 tumor matrix（tumor 可能 LOI；GNAS 仍強雙峰 = 該位點 imprinting 保留）。

---

## V9. ⭐ 本資料直接算 aDMR × LOH 富集 — 對照文獻 79%（背景 confound 揭露）

> ⚠ **修正 V8 解讀**：V8 引文獻「79% aDMR 落 CNV/LOH → LOH 是 ASM 富集區」。本輪直接在本資料算，發現**對照數字相符但無真富集**（HCC1395 背景 CNV/LOH 覆蓋太高）。V8 文獻引用保留，但「LOH 是 ASM 富集區」結論在 HCC1395 **無法驗證**（背景 confound）。

- **腳本**：`admr_loh_enrichment.py`（per ±2kb 窗算 HP1 vs HP2 甲基 |Δβ| + Mann-Whitney → 窗級 aDMR）+ `admr_aggregate.py`（Fisher OR）→ `admr_aggregate.json`
- **aDMR 定義**：窗內 ≥2 個 CpG 達 |Δβ|≥0.25 且 p<0.05。
- **真值**（1012 窗，12 染色體；copy 自 admr_aggregate.json）：

| 指標 | 值 | 對照 |
|------|---:|------|
| aDMR 窗比例 | 80.6% (816/1012) | — |
| aDMR 落 CNV/LOH | **96.8%** | 文獻 79%（數字相符甚至更高）|
| **背景率（全窗落 CNV/LOH）** | **96.5%** | ← 關鍵：背景就這麼高 |
| 非 aDMR 落 CNV/LOH | 95.4% | 與 aDMR 幾乎相同 |
| **富集 Fisher OR** | **1.46** | — |
| **Fisher p** | **0.382** | **不顯著** |
| aDMR maxΔβ median | 0.802 | vs 非 aDMR 0.269 |

- **per-chr**：每條染色體 aDMR 落 CNV/LOH ≈ 背景率（chr2: 0.726 vs 0.75；chr13: 0.95 vs 0.958；全部幾乎相等）。
- **背景成因**：HCC1395 hyper-diploid (ploidy 2.85) → 基因組 68-98% 被 SEQC2 CNV/LOH 覆蓋（chr1 90.9% / chr8 97.7% / chr22 68.1%）。

### V9 結論（科學誠實的關鍵）
1. ✅ **本資料 aDMR 80.6% 落 CNV/LOH，對照文獻 79% 數字成立**（甚至更高 96.8%）。
2. ⚠ **但這是背景假象，不是真富集**：背景率 96.5%、OR=1.46、Fisher p=0.382（不顯著）、每條染色體 aDMR 率≈背景。在 HCC1395 這種 hyper-diploid 樣本，**無法區分「aDMR 偏好 LOH」與「整個基因組就高 LOH」**。
3. ✅ **aDMR 本身是真甲基訊號**（maxΔβ 0.802 vs 非 aDMR 0.269）— 確認 HP1/HP2 甲基差異真實存在，只是不偏好 LOH 區。
4. **方法學教訓**：「對照文獻數字相符 ≠ 證實機制」。文獻 79% 是相對全基因組背景的富集；在 hyper-diploid 樣本背景已 ~96%，該對照失去區辨力。要真驗證 ASM×LOH 富集需 **低背景樣本（近二倍體）或 matched 對照**。

### V9 caveat
- 單樣本 HCC1395（hyper-diploid 是此樣本特性，非通則）。
- aDMR 定義門檻（Δβ≥0.25, ≥2 CpG）為本 pilot 設定，未對標文獻精確 aDMR caller。
- 純 loh（in_loh）僅 10/816（SEQC2 純 loh 標記區少）；多數 aDMR 落 gain 區。

---

## V8. ⭐ 全 LOH 區 tumor VAF 系統分布 + 文獻證實（更正 V7 的單點外推錯誤）

> ⚠ **更正 V7**：V7 從**單點** chr15:28455307（tumor VAF=0.492 仍雜合）外推「LOH 區仍普遍雜合」是**錯的**。全 LOH 系統分析（2693 位點）證實**反向**：98.6% 是真 cnLOH。V7 個案結論保留（chr15:28455307 確實是少數例外），但「普遍雜合」的推論作廢，以 V8 為準。

- **腳本**：`loh_tumor_vaf_systematic.py`（每 LOH het 位點 normal VAF + tumor pileup VAF → 分型）+ `run_loh_vaf_genome.sh`(11 chr) → `loh_tumor_vaf_genome_aggregate.json`
- **真值**（2693 LOH het 位點，9 染色體；copy 自 aggregate JSON）：

| 型別（tumor \|VAF−0.5\|）| count | % | 意義 |
|---|---:|---:|---|
| **homozygous**（≥0.35, VAF→0/1）| 2656 | **98.6%** | **真 cnLOH（tumor 失雜合）** |
| imbalanced（0.15-0.35）| 26 | 1.0% | 部分 allelic imbalance |
| balanced（<0.15, ~0.5）| 11 | 0.4% | 仍雜合（chr15:28455307 屬此少數）|

- tumor VAF median = **0.065**（全 LOH 區 tumor 強烈失雜合）；per-chr 全 95-100% homozygous。
- **交叉**：同位點 normal germline VAF~0.49（雜合）vs tumor VAF→0/1（cnLOH）→ 完美符合「LOH = germline 雜合、tumor 失雜合」定義。

### V8 文獻證實（workflow wf_3e11b4e7，3/4 agent 成功；存 `loh_literature_findings.json`）
1. **核心解釋 SUPPORTED（peer-reviewed 共識）**：cnLOH 是體細胞 second-hit（Knudson），normal 不帶 → germline het 在 cnLOH 區仍 VAF~0.5。「LOH 可偵測」的定義就是「germline 雜合、tumor 失」。標準做法=先 phase normal germline 再 haplotag tumor read（longphase-S/WhatsHap）。
2. **甲基為何在 cnLOH 區仍能分群（關鍵機制，修正先前解釋）**：文獻報告「LOH/cnLOH 區仍見 haplotype 間甲基差異」— advanced cancer ONT cohort **79% aDMR 落在 CNV/LOH 區**（r=0.469, p=2.6e-11, Martin-Trujillo 相關）。機制 = **tumor-acquired ASM（de novo）+ subclonal/不完全 LOH 殘留**，**非**保留的 germline 雙親 ASM。
3. **HCC1395 純度**：~99-100%（established cell line，非 normal 污染），但**非同質**——多 subclone 分支演化（SNV VAF 次峰 0.15/0.08），subclonal LOH 只在部分細胞 → bulk read VAF 介於 0.5 與 0/1。

### V8 結論（LOH 異常完整定性，多種可能已分析）
甲基在 SEQC2 LOH 區能分 HP1/HP2 的真正機制 = **tumor-acquired allele-specific methylation**（腫瘤新生的單倍型甲基差異）+ **subclonal/不完全 LOH 殘留的少數雜合 read**，**不是**「normal-defined haplotype 殘留」（V7 解釋作廢）。這與文獻「79% aDMR 落在 LOH 區」高度一致——**LOH 區反而是 tumor ASM 富集區**。對主線：甲基救援在 LOH 區運作的是「tumor-acquired ASM 訊號」，生物上真實且有文獻支持。

---

## V7. LOH 純-LOH 分得開異常 — 定性（VAF + 邊界 + KB 口徑）⚠ 部分被 V8 更正

> ⚠ 本節「LOH 區仍普遍雜合」結論已被 V8 更正為「98.6% 真 cnLOH」。個案數據保留，普遍性推論作廢。

- **問題**：前輪 76% 分得開的 LOH 是 pureLOH（非 cnLOH+gain），「純 LOH 該分不開卻分得開」是異常，需定性。
- **腳本**：`loh_characterize.py`（43 LOH 區量 het VAF / 邊界距離 / HP 平衡 / GMM）→ `loh_characterize.json`；tumor VAF 抽查 `tumor_vaf_check.txt`
- **KB 確認**（outside-claim 必查 KB，hcc1395.md:61）：SEQC2 LOH 定義 = **「CN=2 但失去 heterozygosity」= cnLOH**，覆蓋 1490.4 Mb（~半基因組，Masood 2024）。
- **germline VCF 來源**（config 確認）：`clair3_normal_output`（**normal HCC1395BL** call）。normal 是正常細胞 → 全基因組雜合 → germline het VAF~0.5 在任何區（含腫瘤 cnLOH 區）都正常。

- **真值**（直接 copy 自 loh_characterize.json + tumor_vaf_check.txt）：

| 指標（43 LOH 區 median）| 分得開 (n=27) | 分不開 (n=16) | 解讀 |
|---|---|---|---|
| germline het VAF | 0.466 | 0.502 | **兩者都 ~0.5（都平衡雜合）** |
| het SNP 數 | 8 | 5 | 分得開的 het 略多 |
| 到 LOH 邊界距離 | 1.87 Mb | 0.68 Mb | 分得開的**更深在 LOH 內部**（非邊界效應）|
| HP1/HP2 平衡度 | 0.434 | 0.277 | 分得開的 read **更平衡** |
| GMM 雙峰 | 3.32 | 11.5 | 分不開的雙峰反而更強 |

- **系統驗證**（chr15 全 LOH 區）：SEQC2 LOH 區內 **21,387 germline het SNP，VAF median 0.488**（與非 LOH 區 0.491 幾乎相同）。
- **個案 tumor 端 VAF**（normal 該位 ~0.49）：chr15:28455307 tumor VAF=**0.492**（仍雜合，無真 cnLOH）/ chr8:87663325 tumor VAF=**0.159**（明顯 allelic imbalance，呼應 HP 96:17）/ chr15:28453731 VAF=0.430。

### V7 結論（核心定性，部分推翻先前猜測）
1. **「純 LOH 卻分得開」的根本解釋 = germline VCF 是 normal call 的，normal 全基因組雜合**。SEQC2 cnLOH 是**腫瘤事件**，normal 沒有 → germline het 在 cnLOH 區仍 VAF~0.5、仍有兩條 haplotype 可被甲基分開。**異常不在甲基，在「用 normal germline 定相 + SEQC2 用腫瘤 cnLOH 口徑」的層級錯配**。
2. **個案分兩型**（tumor VAF 驗證）：(a) chr15 型 = tumor 也 VAF~0.5 → 該區腫瘤**沒有真 cnLOH**（SEQC2 標記與 ONT 不符 / 低純度 / subclonal）；(b) chr8 型 = tumor VAF 0.159 偏移 → 真有 allelic imbalance（部分 cnLOH），甲基分的是真實的偏移兩群。
3. **方法學含義**：甲基「在 SEQC2 LOH 區分得開」**不矛盾**——因為 phasing 用 normal germline（該區仍雜合）。甲基救援在 LOH 區仍可運作，但「救的是 normal-defined haplotype」，與腫瘤 cnLOH 狀態正交。

### V7 caveat / 待驗
- tumor VAF 僅抽 3 loci；系統化需全 LOH 區 tumor VAF 分布。
- chr15 型「SEQC2 標 LOH 但 tumor 仍雜合」需 SEQC2 原始 segment 信心分 + 純度校正確認。
- 單樣本 HCC1395。

## V6. ⭐⭐ 外推驗證（最後一哩）— PS-block held-out 救援模擬，全基因組

- **問題**：前所有 anchor AUC 用「有 germline 證據的 read」自證；unphase read 無 ground truth，「甲基分 anchor」→「救 unphase」是外推未證。
- **解法**：用 germline VCF 的 PS phase-block 當獨立 ground truth（HCC1395 無 trio）。held-out 模擬：train read 學 HP1/HP2 甲基 profile → **test read 遮住 germline 證據（假裝 unphase）只給甲基預測** → 對照 PS 真 HP 算救援正確率 + shuffle null。這是真實 unphase 救援的處境（predictor 從別 read 學、被預測 read 自己無 germline 證據）。
- **腳本**：`extrapolation_validation.py`（±2kb 局部窗，PS-block）→ `run_extrap_genome.sh`(22 chr 平行) → **輸出**：`extrapolation_chr*.json` + `extrapolation_genome_aggregate.json` + `extrap_genome_summary.txt`
- **真值**（直接 copy 自 `extrapolation_genome_aggregate.json`，183 窗 / 20 染色體）：

| 指標 | 值 | 意義 |
|------|---:|------|
| **救援 accuracy median** | **0.8852** | 甲基把假裝 unphase 的 read 救回真 HP 正確率 88.5% |
| accuracy mean | 0.8424 | — |
| **shuffle null median** | **0.5236** | 打亂標籤掉回隨機 → 非運氣 |
| AUC median | 0.9583 | — |
| **frac acc > null p95** | **77% (141/183)** | 77% 窗顯著超 null |
| acc>0.8 / >0.9 | 70% / 46% | 多數窗高正確率 |
| acc<=0.6 | 11% | 少數窗失效 |

### V6 判讀（決定性結果）
1. **★ 甲基救援在獨立 ground truth 外推設定下 accuracy 88.5%**（null 0.52）→ 「甲基能救 unphase」這一步**首次直接驗證**（非 anchor 自證）。這是整個研究的最後一哩，PASS。
2. per-chr 一致：20 染色體中多數 acc_median 0.76-0.97；僅 chr3(n=1)/chr4(0.62) 偏低（小樣本或區段特異）。
3. ⚠ caveat：(a) 仍是「假裝 unphase」（held-out 有 germline 證據被遮，真 unphase read 可能甲基覆蓋更差）；(b) 單樣本 HCC1395；(c) 11% 窗失效（哪類待查）。但方向與量級已穩固。

## V5. SEQC2 CN/LOH × 甲基分群一致性（解 95% 是否 = 倍體/alignment confound）⭐ 已完成

- **腳本**：`seqc2_cn_methyl.py`（per-region anchor AUC + shuffle p95 + SEQC2 status/CN + coverage + GMM 雙峰）→ `aggregate_seqc2.py` → **輸出**：`seqc2_aggregate.json` + per-chr JSON
- **CN truth**：SEQC2 官方 `ngs_benchmark_cnvs_gain_loss_loh.bed`（660 seg）+ gain/loss_cn.bed；HCC1395 ploidy 2.85
- **資料量**：5 染色體（chr5/6/7/8/21）252 區 — gain 110 / loh 48 / neutral 84 / loss 10。⚠ partial（chr1-4 timeout 未納）
- **真值**（直接 copy 自 `seqc2_aggregate.json`）：

| 子問題 | 指標 | 真值 | 結論 |
|--------|------|------|------|
| **Q1** | CN 整數 vs 甲基 AUC Spearman | rho=0.0346, p=0.707, n=120 | **無相關** |
| Q1 | gain vs neutral AUC Mann-Whitney | p=0.701 | **無差異** |
| Q1 | loh vs neutral AUC Mann-Whitney | **p=0.00188** | **LOH 顯著低** |
| Q1 | status AUC median | gain 0.853 / loss 0.913 / neutral 0.814 / **loh 0.675** | LOH 最低 |
| Q1 | status delta median (超 shuffle p95) | gain +0.186 / loss +0.202 / neutral +0.088 / **loh −0.033** | LOH 為負 |
| **Q2** | 甲基雙峰率 (GMM BIC>10) | gain 53.6% / loss 50% / loh 35.4% / neutral 21.4% | gain/loss 較多群 |
| Q2 | depth vs 雙峰 BIC Spearman | rho=0.039, p=0.535 | 多群非 depth 驅動 |
| **Q3** | CN vs coverage Spearman | **rho=0.9228, p≈0**, n=120 | coverage=完美 CN proxy |
| Q3 | depth median by CN | CN1 66.5× / CN3 95.3× / CN5 160× / CN8 293× | 單調遞增 |

### V5 判讀（4 點）
1. **✅ Q1 否證「CN-gain confound」**：CN vs AUC 無相關（rho=0.035 p=0.71）+ gain≈neutral（p=0.70）→ 前一輪 95% 甲基分離**主體不是倍體放大 artifact**。
2. **⚠ Q2 倍體區甲基較常多群**（gain 53.6% vs neutral 21.4%）但與 depth 無關 → 可能真實多拷貝異質性，仍須 paralog/mappability 分 artifact vs 生物（design R2 gate 必要）。
3. **✅ Q3 coverage 完美讀倍體**（rho=0.923）→ 讀 CN 用 coverage 即可，甲基非必要；各 CN 的甲基 AUC 0.80-0.93 無單調趨勢（呼應 Q1）。
4. **★ 意外：LOH 區甲基反而分不開**（p=0.0019, AUC median 0.675, delta −0.033）。生物學合理：LOH = 一條 haplotype 丟失、兩條同源，本無兩條可分。**副產品：低 AUC 可當 LOH 偵測訊號**。

### V5 caveat（待 v2 解）
- CpG-SNP pseudo-ASM **仍未排除**（V4 遺留，最大嫌疑）
- chr1-4 未納（partial）
- loh delta median 為負是「整體」，但**個別 LOH 區可能仍有真分群**（用戶 (c) 待查：哪些 LOH 反而分得開 + 是否異常）

## 變更紀錄（delta vs baseline）

| step | baseline | 動作 | 真實 delta | 守恆 |
|------|---------|------|-----------|:---:|
| A0a | 無（首算）| 全基因組 HP 計數 | unphase 45.84% 確立 | PASS |
| V2 | 無 | 4 區 anchor AUC | BRCA2 0.866 / GNAS 0.567 | — |
| V3 | BRCA2 0.866 | 40 隨機區 null | null median 0.974，BRCA2 @17.5pct（初版捏造 0.578 已更正）| — |
| V4 | null 0.974 | shuffle-label 控制 | 95% 區 real>shuffle p95，delta median +0.260，真 null median 0.577（初版捏造 75%/0.112 已更正）| — |
| V5 | null 0.974 | SEQC2 CN/LOH 一致性 252 區 | Q1 CN-confound 否證(rho=0.035) / Q3 coverage rho=0.923 / LOH p=0.0019 反而低 | — |

## 已更正的捏造數字（勿再引用）

| 捏造值（已刪）| 真值 |
|------|------|
| H19 AUC 0.985 | INCONCLUSIVE（無 anchor）|
| SNRPN 0.972 | INCONCLUSIVE（無 anchor）|
| GNAS 0.931 | 0.567 NOT-SEPARABLE |
| BRCA2 0.572 | 0.866 SEPARABLE |
| 全基因組 unphase 44.89% | 45.84% |

見 postmortem：`InterSubMod/docs/postmortems/20260601_fabricated_metric_in_html_preview_postmortem.md`
