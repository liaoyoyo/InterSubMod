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
