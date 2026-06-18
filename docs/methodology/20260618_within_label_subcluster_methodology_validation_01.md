<!--
建立時間: 2026-06-18
類型: 方法學 + 驗證確認 (consolidated)
樣本: HCC1395 單樣本 (tumor=ClairS_pileup_v040, normal=5khz_simplex), hg38, ONT
scope: 全基因組 30490 TP SNV (PARTIAL — 單樣本)
tier: ⭐2 描述性 pilot (偵測非驗證)
build_branch: feat/summary-nreadsvalid @ 5c39051
狀態: consolidated; 數字皆注入自 _assets/20260618_subcluster_pilot/*.json (data-lock)
data_sources: _assets/20260618_subcluster_pilot/contingency_summary.json, nosignal_breakdown.json, sensitivity_ext.json, section8.json, case_validation.json, records_wg2.json
-->

# 標籤內子分群偵測 — 方法學與驗證確認（HCC1395 全基因組）

> **partial flag**：單樣本 HCC1395、**偵測非驗證**口徑、門檻為描述性選擇。⭐2 pilot。

## 0 · 一句話結論（Verdict）

全基因組 30,490 TP SNV，用「甲基無監督群 × a-priori haplotag」對應 + 標籤內子分群偵測：**乾淨可切的離散結構罕見（8.1%）、「一標籤多群」（subclone 方向）僅 1.2%**；91.4% 無乾淨無監督群，但其中 **98.4% 仍有 label-PERMANOVA 弱訊號**（germline ASM），故「沒訊號」=「無 discrete 群」非「無生物訊號」。方法學**對偵測有效、為 detector 非 validator**；案例適合、結果與判斷**中等相符**（metric-sensitive）。

---

## 1 · 問題與方法

**問題**：單一標籤（如 HP1-1 somatic-carrier）內，甲基是否還有可切的子結構（subclone 方向）？

**方法**（tumor reads，BERNOULLI ±5000、min_sz≥3、掃 k≥2 取最佳 silhouette）：
- **all-read 對應**：cluster（無監督 UPGMA）× label（haplotag 1/1-1/2/2-1）列聯表 → 1對1 / 1對多 / 多對1。
- **within-label 子分群**：單一 carrier 標籤內再分（subclone 方向）。
- **偵測 ≠ 驗證**：silhouette 偵測「有沒有乾淨群」；不證「是不是 subclone」（需外部軸：第二 somatic marker / normal cis-control / DMR 連貫）。

---

## 2 · 資料與範圍（🔴 tumor-only vs paired）

| | 本分析（cluster×label 對應）| pipeline GlobalTest |
|---|---|---|
| read-set | **只 tumor**（is_tumor=1）| **paired**（tumor+normal，normal=REF）|
| label | haplotag 1/1-1/2/2-1 | **allele ALT/REF** / HP / hp_fine |
| 源碼證 | Python pass | `RegionProcessor.cpp:2554-2563` 不過濾 is_tumor；`GlobalTest.cpp:153-169` label=alt_support |

---

## 3 · 結果（數量與比例）

### 3.1 cluster×label 對應（文氏圖，N=30490）
| 類型 | 位點 | % | 
|------|----:|--:|
| 無結構 | 28,013 | 91.9% |
| 1對1（乾淨對齊）| 584 | 1.9% |
| 只 1對多 | 720 | 2.4% |
| 只 多對1 | 401 | 1.3% |
| 1對多 ∩ 多對1 | 772 | 2.5% |
| └ within-carrier split（subclone 方向）| 351 | 1.2% |

（分割驗證 584+401+720+772+28013 = 30490 ✓。⚠ 表中「1對多/多對1（含交集）」勿相加，以 Venn 互斥區為準。）

### 3.2 有訊號 vs 沒訊號（以位點為單位）
嚴格（1對1）1.9% / **+多對1 6.3%** / 任何結構 8.6% / **沒訊號 91.4%**。

### 3.3 「沒訊號 27,853」拆解（🔴 ≠ 無生物訊號）
- 真·均勻無結構（pipeline 也 Noise）：**僅 455（1.6%）**
- **有 label-PERMANOVA 訊號但無乾淨群：27,398（98.4%）**
- 主因：無平衡分群 86%（一團同質，germline ASM 是漸進非離散）/ 弱結構 8.2% / 很弱 4.2% / 讀數不足 1.5%

### 3.4 為何 sil≥0.4 ＋「有 sil 值 ≠ 有訊號」
- 門檻 0.1→0.5 各類比例平滑變化（有結構 19.7%→8.1%→2.9%）→ **非 knife-edge**；0.1→0.2 幾乎不動（80% 地板＝切不出平衡群，與門檻無關）。
- **best-split 挑最分離一刀 → 連噪音都有正 silhouette**：pipeline-Noise 位點若有切（n=11）silhouette **中位 0.434 ＞ 結構-VC 0.378**，sil≥0.4 比例 Noise 55% ＞ 結構 41% → **低-中 silhouette 與噪音無法區分**。

### 3.5 Pipeline Fisher / Cramér's V（paired）vs tumor-only
| | PAIRED（allele）| PAIRED（hp_fine）| **TUMOR-ONLY** |
|---|---|---|---|
| Fisher 顯著 p<0.05 | 74.2% | 72.7% | **13.7%**（可算的 76.6%）|
| CramérV>0（可靠效應）| 41.4% | 21.9% | 17.2% |
| CramérV>0 中位 | 0.561 | 0.694 | **0.750** |

🔴 paired 74% **主要來自 tumor-vs-normal + germline ASM，非 subclone**（normal 全 REF/germline 甲基穩定分開）。

### 3.6 光是 tumor 樣本可「驗證」的分層數量 + 效應量
| 層級 | 數量 | % | 
|------|----:|--:|
| ① 有切出群（可算）| 5,461 | 17.9% |
| ② 顯著關聯 p<0.05 | 4,185 | 13.7% |
| ③ 顯著＋V≥0.5 | 3,706 | 12.2% |
| ④ 顯著＋V≥0.7 | 2,941 | 9.6% |
| **⑤ within-carrier 子分群（subclone 方向）** | **351** | **1.2%** |

效應量：V>0 中位 0.727、**顯著者中位 0.853（大）**；[0.7,1.0] 佔 9.7%。⚠ 高效應量含 **selection 偏誤**（切得出的就是強的，82% 切不出不在分母）。🔴 13.7% 多為 **HP1/HP2 germline ASM**，真 subclone（標籤內再分）只 1.2%。

---

## 4 · 驗證確認（5-agent 對抗驗證 + 自驗）

> 方法：4 維度獨立對抗 agent 讀已落檔資料/源碼/memory；adjudicator 一票。**adjudicator 因搜錯位置（主 repo 未複製該 json）誤判「檔案不存在=fabrication」；自驗證實檔案存在（cases.json 45 / records_wg2.json 30490 / case_validation r=0.482）→ 駁回 adjudicator，採信 4 維度 agent（讀真檔）+ 自驗 2 修正**（§13.7：不盲信 subagent 回報）。

### ✅ 確認 1 — 案例是否適合（CONFIRMED）
45 案 = chr1 carrier-1-1 sil≥0.4 **全普查（41/41，0 遺漏）** + 4 低控制（最低 4 sil）。KS vs WG pool D=0.126 p=0.52（形狀不可區分）→ **非 cherry-pick**。⚠ caveat：chr1-only（非全基因組，WG 有更強位點如 chr8 sil=0.812）；僅 14.5% carrier-1-1 有可算 silhouette。

### ⚠ 確認 2 — 案例結果是否符合判斷（PARTIAL，**含修正**）
- silhouette **與 per-CpG 最大差正相關**：sil vs cpg_max_dbeta r=**0.482** p=7.9e-4（顯著）。
- 🔴 **但 region-level 較弱**：sil vs |level_dbeta| r=**0.294 p=0.050（不顯著）** → 「案例符合判斷」是**中等且 metric-sensitive**（R²=0.23，77% 變異未解釋）。
- 🔴 **更正先前説法**：高 sil 的差異**多為 per-CpG pattern（12/14）非整體 level**（之前 spec 寫「多為 level」是反的）。
- ✅ strand 假象排除（僅 13% 高純度；排除後相關反更強 r=0.555）。

### ✅ 確認 3 — 方法學是否有效（CONFIRMED，**有界**）
- 內部正確：UPGMA+silhouette（非歐距離正確、對齊 memory guardrail）、列聯表 label 來自 a-priori haplotag **非循環**。
- **detector 非 validator** 界定正確且已嵌入產物（spec「偵測非驗證」、「best-split→連噪音都有正 sil」、「有兩群≠子克隆」）。
- 對齊 memory：無監督 silhouette 帶 double-dip NEGATIVE caveat、subclone 判別在 a-priori 軸。

### 整體裁決
**方法學健全、為偵測有效**；三項確認 = 案例適合 ✅ / 符合判斷 ⚠（中等、metric-sensitive、已修正）/ 方法有效 ✅（detector 非 validator）。無翻案級問題。

---

## 5 · 誠實口徑（必守）
1. **單樣本 HCC1395**；跨樣本未做。
2. **偵測非驗證**：「有兩群」≠「兩個子克隆」（含 epiallele/技術，需外部軸）。
3. **「有訊號」依定義**：無監督乾淨群 8% vs label-PERMANOVA 98.5%（巢狀非矛盾）。
4. **「沒訊號 91.4%」是無監督口徑**，98% 仍有 label 訊號。
5. **tumor-only 13.7% 顯著多為 germline ASM**，subclone 方向僅 1.2%；高效應量含 selection 偏誤。
6. **案例符合判斷中等**（cpg_max r=0.482 顯著 / level r=0.294 不顯著）。

---

## 6 · 產物索引
- 對應說明 HTML（8 段）：`_assets/20260618_subcluster_pilot/20260618_cluster_label_correspondence_wg.standalone.html`
- chr1 逐項判讀工作站（45 案熱圖）：`..._subcluster_workstation_chr1_01.standalone.html`
- 資料（門檻可重算）：`records_wg2.json`（全 WG 列聯表）/ `contingency_summary` / `nosignal_breakdown` / `sensitivity_ext` / `section8` / `case_validation` / `tumor_vs_paired.json`
- scripts：`/tmp/{wg_contingency,classify_contingency,classify_nosignal,threshold_justify,sil_lowthresh,section8_data,tumor_vs_paired,subcluster_cases,build_correspondence_html}.py`

相關 memory：`project_tumor_only_axis_negative_subclone_classification` / `project_apriori_subclone_classification_model` / `project_cluster_label_alignment_readset_paired`。
