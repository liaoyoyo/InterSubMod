<!--
建立時間: 2026-06-27
類型: L9 單 sSNV 甲基外推 + HP 排解 cis-ASM + 多→單 sSNV 校正（數據驗證；clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3）
data_sources: data/sm_single_locus_methyl.json, data/sm_single_locus_methyl_perlocus.tsv, data/sm_methyl_genetic_concordance.json, data/sm_methyl_sufficiency_audit.json
-->

# L9 — 單 sSNV 甲基外推 / HP 排解 cis-ASM / 多→單 sSNV 校正（數據驗證）

> **緣起（使用者連續提問）**：(1) 甲基資訊不夠還是無？(2) 能否 characterize 已確認 clone？(3) 能否外推到單 sSNV？(4) 單 sSNV ALT 內甲基多群=subclone 嗎？(5) cis-ASM 能否用 HP 排解？(6) LOH 單一單倍型不就無 ASM？(7) 能否用多-sSNV 真值校正單-sSNV 門檻？
> 本層用**三個新 BAM 分析**逐一數據驗證。資料：單位點 `sm_single_locus_methyl`（1,610 單 sSNV）、校正 `sm_methyl_genetic_concordance`（740 多-sSNV）、cis-CN `sm_methyl_sufficiency_audit`。

## §0 一句話總答

🔴 **甲基有真實但有界的資訊，能 characterize 已由基因型確認的 clone；但不能反向偵測 subclone — 單 sSNV 外推與「多→單校正」皆不成立。** 三個數據鐵證（下）。

## §1 資訊「不夠」非「零」（Q1）

- 甲基差異**存在且 power-gated**：多-sSNV 高功率區（popB_n≥20）recover 達 **54.9%**（`sm_methyl_sufficiency_audit.json`）。非雜訊。
- 但**有界**：整體 powered-rate 僅 **10.68%**，覆蓋限制使絕大多數區達不到高功率。
→ **不夠，非零。**

## §2 能 characterize 已確認 clone（Q2）— ✅

genotype-anchored（用 sSNV 定好的群）後，6.6% 區的甲基**確實區分**遺傳群 + chr17 per-CpG 可歸因。**對已確認 clone 做表觀刻畫＝有效**；**獨立偵測新 clone＝0**（見 L8）。

## §3 單 sSNV 外推（Q3）— 方法可行，但最不可解釋

**單位點 REF-vs-ALT 甲基 MWU+FDR（1,610 單 sSNV，n_links=0）**（`sm_single_locus_methyl_perlocus.tsv`）：

| 指標 | 值 |
|---|---|
| 可測（ALT/REF 各≥3 read）| **1,267 / 1,610** |
| **ASM+（≥1 sig CpG）** | **405（32.0%）** ← **5× 多-sSNV 的 6.6%** |
| LOH | 310/985 = **31.5%** |
| neutral | 95/282 = **33.7%** |
| max Δβ median | **0.965**（all-or-nothing）|

🔴 **32% 看似高，實為「最不可解釋」**：單一位點**無連鎖錨**，那 32% 混了 ① germline cis-ASM（neutral）② subclone/somatic-cis（LOH）③ **C>T 突變直接破壞 CpG 的序列假象**（單位點特有：ALT read 該位點不再是 CpG → 假甲基差異）。**「單 sSNV 有 ASM」不能推論 subclone。** 門檻低（易抓差異）反而使陽性最不可信。

## §4 ALT 多群 = subclone？（Q4）— ❌ read×read 分群在此資料 data-starved

校正分析（`sm_methyl_genetic_concordance.json`）試圖在 740 多-sSNV 區用 **read×read 甲基距離 + PERMANOVA** 測「甲基能否 recover 遺傳 partition」：

| 結果 | 值 |
|---|---|
| 有遺傳 partition（≥2 geno 群）| 740 |
| **read×read PERMANOVA 可測** | 🔴 **僅 1 / 754** |
| sparse_shared_cpg（read 對共享 CpG <5）| **734** |
| recover（p<0.05）| **0** |

🔑 **關鍵發現**：跨多 sSNV 的 read 本就少（需單分子跨完所有 sSNV），其兩兩 CpG 重疊太稀疏 → **無監督 read×read 甲基分群在 99.9% 區根本算不出來**。所以「ALT 內甲基多群」這個做法**在此資料上不可行**（非「分不出 subclone」，是「分群本身 data-starved」）。穩健的 per-CpG 測試（不需 read 兩兩重疊）才得到 6.6%（L8）。

## §5 cis-ASM 能否用 HP 排解（Q5）— ❌ 此資料不可行

- **hp_control_eval = 0/740（多-sSNV）+ 0/1,267（單-sSNV）** —— HP1-vs-HP2 對照**從未可評估一次**。
- 原因（接 Q6）：corroborated 區以 **LOH 為主＝同質＝無 germline 雜合 → longphase 無法 phasing → BAM 根本沒有 HP tag**。沒有 HP 就無從用 HP 排解。
→ **HP 排解 cis-ASM 在 HCC1395 此資料不可行**（黃金備援仍是 matched-normal）。

## §6 LOH 單一單倍型 → 無 germline ASM（Q6）— ✅ 使用者正確（已落 commit dac85ea）

- **49 corroborated 中 41（84%）是 LOH** → germline 等位 ASM 結構上不可能（只剩一個等位）。
- 故 LOH 的 REF-vs-ALT 甲基差異 **非 germline cis-ASM**，是 subclone/somatic-cis 候選（仍需排 double-dip + somatic-cis）；neutral 8 才可能 germline cis。
- LOH 的 hp_control_eval=0 是「對照 moot」非「沒能排除」。

## §7 多→單 sSNV 校正（Q7）— ❌ 無可轉移門檻

要從多-sSNV 真值校正一個甲基門檻轉移到單 sSNV，需要 (a) 多-sSNV 區甲基能 recover 已知遺傳 subclone，且 (b) 該 recover 強度有清楚門檻分開「有/無 subclone」。兩者皆**不成立**：
1. **read×read 分群 data-starved**（1/754）→ 沒有可校正的 read-clustering 統計量。
2. **穩健 per-CpG recover 由覆蓋驅動**（power-gated，6.6%→54.9% 隨覆蓋），**非乾淨的甲基強度門檻** → 即使轉移也是在轉移「覆蓋度」非「subclone 訊號」。
→ **無法用多-sSNV 校正單-sSNV 的 subclone 偵測門檻。**

## §8 數據驗證總表（每數字可 grep）

| 問題 | 數據 | 答案 |
|---|---|---|
| 資訊夠嗎 | high-power recover 54.9% / powered 10.68% | 有界非零 |
| characterize 已確認 clone | 6.6% recover + chr17 歸因 | ✅ 可 |
| 單 sSNV 外推 | 單位點 ASM 32%（但無錨+CpG假象）| 方法可、判讀不可 |
| ALT 多群=subclone | read×read 1/754 可測 | ❌ data-starved |
| HP 排解 cis-ASM | hp_eval 0/740+0/1267 | ❌ LOH 無 HP |
| LOH 無 germline ASM | 41/49 LOH | ✅ 正確 |
| 多→單校正門檻 | 上述 2 條皆敗 | ❌ 無門檻 |

**結論**：甲基在 clone/subclone 是**「對已確認結構的有界表觀註解」**，不是**偵測器**；單 sSNV / read-clustering / HP-排解 / 跨尺度校正在 HCC1395 單樣本資料上**均不足以獨立判定 subclone**。要前進需 **matched-normal cis-control**（B 軸）+ 更深覆蓋。
