<!--
建立時間: 2026-06-18（分析跑至 2026-06-19）
類型: 觀察標準紀錄（verification-table-as-observation-standard）
狀態: in_progress · 單樣本 HCC1395 · ⭐2 偵測非驗證
-->
---
title: WG k-sweep 對齊式 k-selection + 切不出原因 — 觀察標準紀錄
sample: HCC1395 (tumor reads, BERNOULLI ±5000, all 22 chr)
tier: 2
scope: 全基因組單樣本（partial: single-sample, cross-sample 未做）
build_commit: 7789f4a
branch: feat/summary-nreadsvalid
binary_commit: 5c39051
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/ksweep_wg_summary.json,docs/methodology/_assets/20260618_subcluster_pilot/ksweep_wg_records.json,docs/methodology/_assets/20260618_subcluster_pilot/split_accounting.json,docs/methodology/_assets/20260618_subcluster_pilot/cantsplit_reasons.json,docs/methodology/_assets/20260618_subcluster_pilot/cantsplit_apriori_rescue.json,docs/methodology/_assets/20260618_subcluster_pilot/ksweep_summary_tumor.json,docs/methodology/_assets/20260618_subcluster_pilot/ksweep_summary_merged.json,docs/methodology/_assets/20260618_subcluster_pilot/kprofile_summary.json
memory: project_subcluster_cluster_count_determination, feedback_verification_table_per_data_answer
---

# WG k-sweep 對齊式 k-selection + 切不出原因 — 觀察標準

> **本檔職責**：依 `feedback_verification_table_per_data_answer` 標準，把本輪觀察的 headline 數字 + 完整驗證表落檔成「觀察標準」，供後續 drift 比對。所有數字 L1（從原始 `records_wg2.json` / WG records 重算或交叉驗證）。
> **partial flag**：單樣本 HCC1395 ⭐2 偵測非驗證；cross-sample 未做。

---

## 總結論（觀察 A-D 整合 + dual-panel 視覺確認）

| # | 結論 | headline 數據（L1）| 對應觀察 |
|---|---|---|---|
| ① | **silhouette 系統性切太少** | ~20% 可分群對齊最佳 k≠silhouette k（86-93% 更高 k）；FDR 硬 union **69**（HP18/CAR57/ALL40）| A |
| ② | **「夠好 vs 次好」取決於判準** | sil margin 0.047(唯一2.9%) vs **align 0.188(唯一35.2%, 4×)** → alignment 是更好判準 | D |
| ③ | **多解析度真階層** | multi-resolution **31.6%**；不同 k 對齊不同軸（chr14:97275848 k2→allele/k3→HP/k4→carrier）| D |
| ④ | **切不出 ≠ 沒讀數/沒訊號** | 24493(80.3%) 中 98.3% 讀數充足但同質；multi-label 92.4% 有 a-priori 訊號、真 null 7.6% | B |
| ⑤ | **subclone 用 tumor-only；merged=cis-control** | merged HP 軸 FDR rate ×6.4=germline-cis 層；CARRIER/ALLELE 幾乎不動 | C |

**🔴 dual-panel 視覺確認（2026-06-19~20）**：所有觀察熱圖改 SoT `ism_heatmap_std.mpl_dual_panel`（左甲基 read×CpG + 右 read×read 距離）。右側距離圖**對角塊明顯** → 把 FDR/multi-resolution 切法從「統計顯著」升級「**視覺確認是真結構非錯覺**」（chr8:131074695 3 對角塊=align-k=3 真；chr14 多層塊=真階層）。

**🔴 紅線（不變）**：multi-resolution ≠ subclone（跨軸多含 HP/allele=germline-cis，需 normal cis-control）；ambiguous ≠ 無結構（k-歧義）；單樣本 ⭐2；margin 門檻(0.05/0.15)是約定。

---

## 觀察 A — 全基因組 k-sweep：silhouette-k vs 對齊-a-priori-軸-k

**問題**：現行 pipeline 用 silhouette 選 best_k；改用「切群 vs a-priori 標籤(HP/CARRIER/ALLELE)的 CramérV 對齊」選 k，有多少 loci 的 silhouette-k 其實不是最對齊的 k？（加顯著性 gate + 多重比較校正）

**方法**：22 chr 重跑 ISM binary → 每 region 掃 k=2..min(5,n//3) → 每 k 記 silhouette + 三軸 CramérV/chi2-p/Cochran 期望格。gate = ΔV>0.1 + min 期望格≥5；校正 = within-locus Bonferroni×(該軸掃的 k 數) + genome-wide BH-FDR(q<0.05)。非循環（標籤源自 BAM HP tag / variant base，獨立於甲基距離；C++ 源碼層已驗 `ReadParser.cpp:126/165`）。

> **🔴 三軸精確定義（源碼校正 2026-06-19，回應 CARRIER 命名質疑）**：三軸用**同一個 BAM `HP` tag + anchor 鹼基**、僅切法不同 —
> - **ALLELE** = read 在 **anchor 位點**鹼基 ALT/REF（`determine_alt_support` `ReadParser.cpp:165,262`；**anchor 專一**）
> - **HP(family)** = 親代染色體 `{1,1-1}`=HP1 / `{2,2-1}`=HP2（`get_hp_family_label` `FisherExact.cpp:415`）
> - **CARRIER**（保留名+加註）= germline 相位塊 `{1,2}` vs somatic 相位塊 `{1-1,2-1}`（**合併 HP**）。⚠ `1-1/2-1` = longphase「**somatic phase block**」（被 somatic 變異定相的 read，`ReadParser.cpp:154`），**非 anchor ALT**（那是 ALLELE）；帶源碼循環依賴警告（`germline_hp_only` 會降為 0）；且合併 HP，比 pipeline 逐-HP `SubcloneDbeta`（`NGerm`/`NCarrier`）粗。「carrier」沿用 pipeline 用語但僅在逐-HP 才精確。→ **subclone 解讀優先 ALLELE(anchor) + 原生 SubcloneDbeta(逐-HP)**。

### 驗證表 V1（每數字溯源 + 重算）

| # | 數字 | 值 | 來源:key | 重算/加總/交叉 | L | 狀態 |
|---|---|---|---|---|---|---|
| 1 | WG 可分群 region | **5997** | `ksweep_wg_records.json` len | **= `split_accounting.cansplit` 5997**（獨立第二次全跑相符） | L1 | ✓ |
| 2 | silhouette-opt-k == best_k | 5956/5997 (99.3%) | join `records_wg2.json` | 41 差=tie-break；同演算法 99.3% 一致 | L1 | ✓ |
| 3 | HP comparable | 4110 | `HP_axis.n_comparable` | eq 3250 + ne 860 = 4110 ✓ | L1 | ✓ |
| 4 | HP sil-k≠對齊-k raw | 860 (20.9%) | `HP_axis.n_silK_ne_alignK_raw` | 更高 k 740 (86.0%) | L1 | ✓ |
| 5 | **HP 真不是最對齊 k (FDR)** | **18 (0.44%)** | `HP_axis.n_sil_not_most_aligned_FDR` | Bonf 18 = FDR 18 | L1 | ✓ |
| 6 | CARRIER comparable | 5376 | `CARRIER_axis.n_comparable` | eq 4333 + ne 1043 = 5376 ✓ | L1 | ✓ |
| 7 | CARRIER sil-k≠對齊-k raw | 1043 (19.4%) | `CARRIER_axis.n_silK_ne_alignK_raw` | 更高 k 968 (92.8%) | L1 | ✓ |
| 8 | **CARRIER 真不是最對齊 k (FDR)** | **57 (1.06%)** | `CARRIER_axis.n_sil_not_most_aligned_FDR` | Bonf 61 → FDR 57 | L1 | ✓ |
| 9 | ALLELE comparable | 5446 | `ALLELE_axis.n_comparable` | eq 4349 + ne 1097 = 5446 ✓ | L1 | ✓ |
| 10 | ALLELE sil-k≠對齊-k raw | 1097 (20.1%) | `ALLELE_axis.n_silK_ne_alignK_raw` | 更高 k 955 (87.1%) | L1 | ✓ |
| 11 | **ALLELE 真不是最對齊 k (FDR)** | **40 (0.73%)** | `ALLELE_axis.n_sil_not_most_aligned_FDR` | Bonf 40 = FDR 40 | L1 | ✓ |

### 結論 A
- **raw**：~20% 可分群 loci 的對齊最佳 k ≠ silhouette 最佳 k，且 86–93% 是**更高 k** → silhouette **系統性 under-split**（相對生物軸保守）。
- **嚴格 gate + 雙重校正後硬數字**：HP 18 / **CARRIER 57** / ALLELE 40（0.4–1.1%）= 高信心「silhouette 切太少、某更高 k 顯著更對齊」位點。CARRIER（subclone 方向）57 個最值得追。
- 強例：chr8:131074695 (n=140) CARRIER：sil k=2 V=0.293 → k=3 V=0.874 (p_bonf=1.8e-23)。
- **caveat**：對齊式 k 是「少數位點(0.4–1.1%)精煉」，非全面取代 silhouette；單樣本 ⭐2。這批 = plan「silhouette/balance gate 切太嚴漏掉的 minor 子結構」WG+FDR 量化。

---

## 觀察 B — 切不出（24493）原因分析

**問題**：為何 80.3% loci 切不出 ≥2 群？是讀數不足、無訊號、還是別的？

### 驗證表 V2

| # | 數字 | 值 | 來源:key | 重算/加總/交叉 | L | 狀態 |
|---|---|---|---|---|---|---|
| 1 | 切不出 total | 24493 | `cantsplit_reasons.cant_total` | A+B+C=412+1+24080 ✓；= N−cansplit | L1 | ✓ |
| 2 | A 讀數不足 (n<6) | 412 | `A_insuff_n_lt6` | = `split_accounting.insuff` | L1 | ✓ |
| 3 | B peel-fail (NaN) | 1 | `B_peel_fail_NaN` | 幾乎不存在 | L1 | ✓ |
| 4 | **C 讀數夠但甲基同質** | 24080 | `C_cohesive_homogeneous` | **median n=59 reads**（非讀數問題） | L1 | ✓ |
| 5 | C single-label | 2201 | `C_single_label` | single+multi=2201+21879=24080 ✓ | L1 | ✓ |
| 6 | C multi-label | 21879 | `C_multi_label` | — | L1 | ✓ |
| 7 | multi-label 有 a-priori 訊號 | 20222 (92.4%) | `cantsplit_apriori_rescue.C_multi_ANY_apriori_sig` | matched 21879/21879 | L1 | ✓ |
| 8 | 其中 SubcloneDbeta（真 subclone 軸） | 3090 (14.1%) | `C_multi_subclone_dbeta_sig` | — | L1 | ✓ |
| 9 | 真 null（所有 a-priori 軸） | 1657 (7.6%) | `C_multi_truly_null` | 20222+1657=21879 ✓ | L1 | ✓ |

### 結論 B
- **切不出 ≠ 沒讀數**：98.3% 是「讀數充足（median 59）但甲基同質」，非覆蓋不足。
- **切不出 ≠ 沒訊號**：multi-label 92.4% 有 a-priori 訊號（PERMANOVA/Δβ 測得到），真 null 只 7.6%；但真 subclone 軸只 14.1%。
- **深層原因**：clustering 需「離散雙峰 gap」，a-priori 檢定只要「平均差」。切不出多為 mean-shift/漸層/少數 CpG，無乾淨 gap → silhouette 過不了門檻 → clustering 正確回 1 群。
- **與觀察 A 的分界**：k-sweep（A）精煉「已可分群的 5997」；切不出的 24493 任何 k 都形不成群，只能靠 a-priori Δβ（B）。

---

## 觀察 C — tumor-only vs merged(tumor+normal) 對比

**問題**：把 normal reads 合併進分群+標籤（= paired），結果與 tumor-only 差多少？意義？
**方法**：dual-mode（`ksweep_wg_v2.py`，同一次 binary pass 算兩種 read-set）。tumor-only=`is_tumor=1`（重現觀察 A）；merged=tumor+normal（normal 全 germline `{1,2}`/REF、無 somatic 標籤；allele 排除 normal 的 UNKNOWN）。各自 analyze → 分開 summary。

### 驗證表 V3（tumor vs merged，每數字溯源）

| # | 數字 | tumor | merged | 來源:key | 檢查 | L | 狀態 |
|---|---|---|---|---|---|---|---|
| 0 | 交叉驗證 tumor(v2) FDR | 18/57/40 | — | summary_tumor vs ksweep_wg_summary | **= 原值 ✓**（v2 腳本正確）| L1 | ✓ |
| 1 | 可分群 region | 5997 | **18864** | records 檔 len | merged 3.1×（normal 補足 ≥6 read）| L1 | ✓ |
| 2 | HP comparable | 4110 | 18771 | `HP_axis.n_comparable` | — | L1 | ✓ |
| 3 | **HP FDR硬 (rate)** | 18 (0.44%) | **529 (2.82%)** | `HP_axis.n_sil_not_most_aligned_FDR` | **rate ×6.4** | L1 | ✓ |
| 4 | CARRIER FDR硬 (rate) | 57 (1.06%) | 271 (1.44%) | `CARRIER_axis...FDR` | rate ×1.36 | L1 | ✓ |
| 5 | ALLELE FDR硬 (rate) | 40 (0.73%) | 173 (0.92%) | `ALLELE_axis...FDR` | rate ×1.26 | L1 | ✓ |

### 結論 C
- **merged 可分群 3×（5997→18864）**：normal 補足 read 過門檻，**非新 subclone**。
- **FDR rate 變化才是真訊號**：**HP 軸暴增 ×6.4**（normal 是乾淨 germline HP1/HP2 → 強化 germline 單倍型 ASM = **cis-ASM/imprinting 層**）；CARRIER ×1.36、ALLELE ×1.26 幾乎不動（normal 無 somatic，per-locus 不幫忙；多出的絕對數來自 3× 大分母）。
- **與 69 tumor-FDR 位點預覽不矛盾**：那 69 個加 normal 後對齊 V 下降（~80%）= tumor-specific 結構被稀釋；全量 merged 點亮的是**另一批 HP/germline 軸位點**（normal 自帶的 germline ASM）。
- **🔴 定案**：**subclone 偵測用 tumor-only**（CARRIER/ALLELE 在 tumor 乾淨，merged 對這兩軸 per-locus 無益）；**merged = cis-control / germline 層揭露**（HP 軸 ×6.4 對應「≥2 群≠subclone、85% 可能 cis-ASM」紅線）。

---

## 觀察 D — 切法品質 k-profile（夠好 vs 次好 + 唯一性 + 多解析度）

**問題**：k-sweep 切法是「夠好(唯一最佳)」還是「次好(近平手)」？每位點切法是否唯一？k=2 切成 k=3 是否也有意義？
**方法**：`ksweep_kprofile.py` 讀 per_k 曲線算 margin(best-k − 2nd-best-k；silhouette + alignment-V 兩種)；分類 single-k-forced(len=1) + k-choice 三態。代表位點重跑 binary 畫 k=2/3/4 三聯熱圖。Phase 1-2 零 compute。

### 驗證表 V4（tumor，每數字溯源）

| # | 數字 | 值 | 來源:key | 重算/加總 | L | 狀態 |
|---|---|---|---|---|---|---|
| 1 | 可分群 | 5997 | `tumor.N` | =cansplit | L1 | ✓ |
| 2 | single-k-forced(只能k=2) | 3782 (63.1%) | `single_k_forced` | single+choice=N ✓ | L1 | ✓ |
| 3 | k-choice 子集 | 2215 | `k_choice_n` | — | L1 | ✓ |
| 4 | sil-margin mean/唯一(≥.15) | 0.047 / 64 (2.9%) | `sil_margin` | — | L1 | ✓ |
| 5 | **align-margin mean/唯一(≥.15)** | **0.188 / 737 (35.2%)** | `align_margin` | **mean 4× sil、唯一 11×** | L1 | ✓ |
| 6 | align>sil margin | 57.3% (1198/2091) | `align_gt_sil` | — | L1 | ✓ |
| 7 | 三態 multi-resolution | 700 (31.6%) | `three_state` | 700+103+1412=2215 ✓ | L1 | ✓ |
| 8 | 三態 confident-unique | 103 (4.7%) | | | L1 | ✓ |
| 9 | 三態 ambiguous-near-tie | 1412 (63.7%) | | | L1 | ✓ |

### 結論 D
- **「夠好 vs 次好」取決於判準**：silhouette mean margin 0.047（唯一 2.9%）= 多數次好/近平手；**alignment mean margin 0.188（唯一 35.2%，4×）= 多數變夠好/唯一**。→ **「最新版更好」的量化證明**：alignment 把 silhouette 的近平手消歧義（57% align>sil）。merged 同向（align 0.229 vs sil 0.043，67% align>sil）→ 結論穩健。
- **三態**（k-choice 2215）：multi-resolution 31.6% / confident-unique 4.7% / ambiguous 63.7%；另 63.1% single-k-forced（n 太小無從選 k）。
- **multi-resolution 是真階層**：代表位點不同 k 對齊不同軸（如 chr14:97275848 k2→allele(.96)/k3→HP(1.0)/k4→carrier(.82)）→ **k=2 切成 k=3 確實也有意義**（捕捉不同生物結構）。
- **🔴 紅線**：multi-resolution ≠ subclone（跨軸多含 HP/allele=germline-cis，需 normal cis-control）；ambiguous ≠ 無結構（是 k-歧義）；單樣本 ⭐2；margin 門檻(0.05/0.15)是約定。
- 產物：`kprofile_summary.json` + `kprofile_loci_{tumor,merged}.json` + `20260618_kprofile_explainer_01.standalone.html`（11 代表熱圖）+ `figs_kprofile/`(29)。

---

## Provenance & 重現
- binary: `5c39051`（tumor=ClairS_pileup_v040 / normal=5khz_simplex；分析端 filter 至 tumor）
- compute: `scripts/ksweep_wg.py`（22 chr disk-safe，1675s）→ `ksweep_wg_records.json`
- 後分析: `scripts/ksweep_analyze.py 0.10 0.05` → `ksweep_wg_summary.json`
- 觀察 C dual-mode: `scripts/ksweep_wg_v2.py`（2046s）→ `ksweep_records_{tumor,merged}.json` → `ksweep_analyze.py … {records}.json {summary}.json` → `ksweep_summary_{tumor,merged}.json`（tumor 重現原值=交叉驗證）
- 觀察 D k-profile: `scripts/ksweep_kprofile.py`（零 compute）→ `kprofile_{summary,loci_tumor,loci_merged}.json`；`select_kprofile_examples.py`→mini-VCF→binary→`plot_kprofile_heatmaps.py`→`figs_kprofile/`(29)；`build_kprofile_explainer.py`→explainer HTML
- 切不出: 全量讀 `records_wg2.json` → `cantsplit_reasons.json` + join `output/_wg_bdcprime_verify/significance_summary.csv` → `cantsplit_apriori_rescue.json`
- 全部腳本/json 在 `docs/methodology/_assets/20260618_subcluster_pilot/`

## 下一步候選
1. 把 CARRIER 57 個 FDR-顯著位點逐一觀察（熱圖 + k=2 vs align-k 對照），確認是真子結構非 artifact。
2. 對齊式 k-selection 整合進 within-HP 偵測（plan D 系列）。
3. cross-sample 重現（COLO829/H1437）。
