<!--
建立時間: 2026-07-06
類型: 觀察方法 + pilot 結果紀錄 (完整驗證進行中)
狀態: in_progress（pilot 完成;6 樣本完整版 permutation+FDR+cis-control 執行中）
主軸: Subclonal reconstruction — 甲基是否與 sSNV 多點分群(lineage)軸有關
data_sources: /big7_disk/liaoyoyo2001/cnv_sv_work/methyl_lineage/methyl_vs_ssnv_lineage_*.json, scripts/analysis/methyl_vs_ssnv_lineage.py
build_branch: research/subclonal-reconstruction-202606
關聯 memory: project_ont_cnv_sv_subclone_verification_feasibility
-->

# 甲基 vs sSNV 多點分群(lineage)軸 — 觀察方法 + pilot 紀錄

> **版本** v0.3 (2026-07-06) · **狀態：CN-aware 雙峰率完成、裁決修正定案** · tier：CN 率=**L1**、甲基 residual=**L3**
> **一句話（修正）**：群內甲基雙峰真實存在(~15%)，**雙峰率跨 CN 平坦(H2009 gain 0.170≈neutral 0.163=1.04x)→ 是 CN-獨立的真實現象、非 CN-gain 假象**（3b 的「77-86% 在 gain」是 base-rate 假象已修正）；但屬 within-lineage 隱藏子結構 L3、bulk 不能證 subclone → 甲基=有界輔助（有 CN-獨立內容）。難建樹區 Q5 輔助候選 7 個。

## 0. 為何做這個（修正之前的盲點）

之前甲基 vs 遺傳軸的驗證**只用 HP（單倍型）+ 單 sSNV ALT/REF**當標籤（`perread_table` 欄 `hp`/`alt_support`），**從未用 sSNV 多點共現分群(population/lineage) + 樹結構**當標籤。而 sSNV 共現 lineage 才是**非循環的 subclone 骨幹**。本分析補這個軸：**同一多-sSNV 樹區內，不同 sSNV-lineage 的 read 甲基是否不同 / 甲基非監督分群是否對齊 lineage**。

## 1. 觀察方法（read-only，最快最安全）

> tagged BAM **同時含** sSNV allele（read 序列）+ 5mC（MM/ML tag）→ **不需外部數據、不改 BAM**。

`scripts/analysis/methyl_vs_ssnv_lineage.py`（read-only pysam pileup tagged BAM）：
- **樹區來源**：`ml_part_*.json` 的 group（n_sSNV≥2 且 n_populations≥2 且 n_pop_with_ALT≥2、span≤100kb 去 Mb-artifact）。ref/alt 從 `sm_linkage_genomewide.json::census`。
- **per read**（跨區）：① 各 sSNV 位點 allele(R/A) → 多點 genotype 字串 = **lineage/population**；② 區內 5mC 甲基均值（ML qual>128=methylated）。
- **per region 指標**：
  - `delta_beta` = population 間甲基均值最大兩兩差（pops ≥4 read）
  - `perm_p` = **permutation null**（打亂 population 標籤 1000×，Δβ 經驗 p）→ **取代 pilot 用的 n-驅動 MWU 假顯著**
  - `single_sSNV_max_db` = 各單 sSNV 的 allele-Δβ(ALT vs REF) 最大 = **per-locus cis-ASM proxy**；`lineage_beyond_cis = delta_beta − single_sSNV_max_db > 0` → 訊號**超越單點 cis**（多-lineage/subclone 層）
  - `ari_2grp` = 甲基中位分割 vs top-2 population 的 ARI = **對齊** proxy
  - `perm_q_fdr` = BH FDR 校正
- **cis-control 邏輯**：`lineage_beyond_cis` 判「是否只是單點 cis-ASM」；完整 cis 天花板仍需 matched-normal DMR 對照（⭐4）。

## 2. Pilot 結果（H2009, 10 唯一樹區, L3 初步）

> 🔴 pilot 用**未校正的 MWU**（完整版已改 permutation）。

| 指標 | 值 |
|---|---|
| median Δβ（lineage 間甲基差）| **0.114** |
| MWU_p<0.05（未 FDR）| 2/10 |
| 最強 | chr22:16022148（3-sSNV,110read,5 群）Δβ=**0.199**，`ARR` 群甲基 0.324 vs 其他 ~0.14 |

**pilot 初判（L3，不可當定論）**：**部分 sSNV-lineage 之間確有甲基差異**（比 CN 軸更直接的訊號；chr22 有真 Δβ 0.2），但 **modest + 不一致**（2/10、median 0.114），且部分「顯著」是大 n 驅動（如 chr13 Δβ 0.048/p=0.007 = 假訊號）。

## 3. 🔴 完整版轉向 — 用並行 session 成熟的 methyl_auxiliary + 我的 SAVANA CN（2026-07-06 更新）

> §restraint：pilot 的「Δβ between populations」在**重造**並行 session 已有的 `methyl_auxiliary_annotation.py`（更成熟：per-genotype-group GMM 雙峰 + HP/CN 排除）。ml_part 被並行 session 重整移除 → pilot approach 廢；改用既有管線 + **我的 SAVANA CN 補 `n_residual_candidate_cnclean`**（原全 0 因 cn=unknown）。

### 3a. methyl_auxiliary residual-candidate（既有，L3）
每 sSNV-genotype 群測甲基雙峰(GMM BIC + 峰差≥0.2 + 子組≥4)，扣 HP + CN(gain) 後仍雙峰 = residual-candidate。三樣本：

| 樣本 | tested | residual-candidate | (原) cnclean |
|---|---|---|---|
| H1437 | 4523 | 216 (4.8%) | 0（cn=unknown）|
| H2009 | 3428 | 589 (17.2%) | 0 |
| HCC1954 | 2234 | 294 (13.2%) | 0 |

### 3b. 🎯 套 SAVANA SM_CNBED（我的貢獻，L1 CN 定量）— 揭露 CN-gain confound
腳本 scratchpad `residual_cnclean_with_savana.py`（read-only）；輸出 `cnv_sv_work/methyl_lineage/residual_cnclean_savana.json`。

| 樣本 | residual CN 分布 | **gain(CN-confound)** | cn-clean(≠gain) | 純 neutral | clean Δβ median |
|---|---|---|---|---|---|
| H1437 | gain 166 / loh 49 / neutral 1 | **166 (77%)** | 50 | 1 | 0.258 |
| H2009 | gain 178 / loh 260 / neutral 132 / loss 19 | 178 (30%) | 411 | 132 | 0.250 |
| HCC1954 | gain 253 / neutral 21 / loh 19 / loss 1 | **253 (86%)** | 41 | 21 | 0.334 |

## 3c. 🔴 CN-aware 雙峰**率** by CN state（因果證明，2026-07-06 修正 3b 讀法）

> 腳本 `scripts/analysis/methyl_bimodality_cn_rate.py`（8-chunk 平行,GMM 同 methyl_auxiliary,cn 從 SM_CNBED）。**率 = 每 CN 狀態下「雙峰群數/測試群數」**（分母=所有測試群,3b 的 count 缺分母）。

| 樣本 | gain 率 | neutral 率 | gain/neutral | 判讀 |
|---|---|---|---|---|
| H1437 | 0.049 (167/3399) | 0.020 (1/49) | 2.41x | ⚠ neutral n=49 太小不可靠 |
| **H2009** | 0.170 (194/1142) | 0.163 (**132/808**) | **1.04x** | ✅ 最有 power → **幾乎相等** |
| HCC1954 | 0.154 (300/1949) | 0.123 (22/179) | 1.25x | 略高 |

🔴 **修正 3b**：3b 的「77-86% residual 在 gain」是 **base-rate 假象**（gain 區佔基因組大比例→大多數群本來就在 gain,分母效應）,**非 CN-gain 造成更多雙峰**。**用率的正確因果結論：CN-gain 未顯著提高群內甲基雙峰率**（最有 power 的 H2009=1.04x 平坦）→ **雙峰率 ~15% 與 CN 狀態無關**。

## 4. 裁決（修正版，L1 CN 率 + L3 甲基）

- **群內甲基雙峰確實存在**（H2009/HCC1954 ~15% 率、H1437 ~5%）。
- 🔴 **率跨 CN 平坦（H2009 gain 0.170≈neutral 0.163）→ 雙峰是 CN-獨立的真實現象、非 CN-gain 假象** — 這比 3b 的 count 讀法更正確(count 被 base-rate 誤導)。→ **甲基確有超越 CN 的獨立內容**（比先前更正面）。
- **但仍是「群內隱藏子結構」(同一 sSNV lineage 內)非「lineage 之間」** → **L3 候選**,可能 cis-ASM/HP-殘留/技術,**bulk 不能證 subclone**。
- **總結**：甲基有 **CN-獨立的群內雙峰內容(真實)**,但屬**有界輔助 L3**(within-lineage 隱藏子結構,非確認 subclone)→ 與主軸一致。

## 3e. 🔴 前判斷分類（determinacy）× 甲基可處理 — 完整交叉統計（07-08）

> 腳本 `methyl_bimodality_cn_rate.py`(加 det_ct 記錄)+`methyl_cnrate_merge.py`；3 CN 樣本 pooled 20,450 可測區全覆蓋(含 E_subcube 9460 gap#1 救回區)。每區取其 genotype 群最佳甲基狀態。

| 前判斷分類 | 總區 | residual(甲基可助L3) | 比例 | cn_gain | hp | none | Q5 |
|---|---|---|---|---|---|---|---|
| A_determined | 4965 | 204 | 4.1% | 341 | 13 | 4407 | 0 |
| B_pairwise | 2197 | 83 | 3.8% | 63 | 2 | 2049 | 0 |
| E_subcube(gap#1) | 9460 | 2 | 0.0% | 4 | 2 | 9452 | 0 |
| A_ambiguous | 167 | 8 | 4.8% | 3 | 0 | 156 | 0 |
| C_underdetermined | 850 | 28 | 3.3% | 9 | 2 | 811 | 7 |
| recurrence_required | 672 | 25 | 3.7% | 29 | 1 | 617 | 2 |
| incompatible | 332 | 1 | 0.3% | 0 | 0 | 331 | 0 |
| other | 1807 | 84 | 4.6% | 61 | 5 | 1657 | 0 |
| **ALL** | **20450** | **435** | **2.1%** | 510 | 25 | 19480 | **9** |

**🔴 決定性結論**：① 甲基可助全體只 **2.1%**(95.3% 無訊號)；② **難建樹區(2021)甲基可助只 62=3.1%、Q5 只 9** → 甲基對「定不出樹」區幾乎沒幫助(比例還低於已定樹 4.1%)；③ **E_subcube(gap#1 救回 9460=46%)甲基訊號幾近零(0.02%)** — partial-read 群太稀疏(多 <MINGRP=12 不可測)。residual 集中已定樹/非難建區 → 甲基=已知結構弱佐證、非未知結構解算器。檔 `cnv_sv_work/methyl_lineage/determinacy_x_methyl_complete.json`。

## 3d. Q5 輔助旗標 — 難建樹區 ALT-cluster 甲基雙峰候選

在 `C_underdetermined`(sSNV 定不出樹) 區 + ALT-cluster + cn-clean(neutral) + 扣 HP → **7 個候選**(H2009 6 + HCC1954 1),Δβ 0.20-0.30。例 H2009 chr10:31944510 AA 34read@0.75 vs 11read@0.46。→ **可作難建樹區的 L3 甲基輔助標記**(候選隱藏 subclone/branching),但稀少+需 normal-cis/single-cell 確認。觀察 HTML `cnv_sv_work/methyl_lineage/`。

## 5. Provenance / caveat（更新）

- 3b 用 SAVANA SM_CNBED（H1437/H2009/HCC1954 = usable purity；COLO829/HCC1937 mis-fit 排除）。
- residual-candidate = 並行 session `methyl_auxiliary_annotation.py`（GMM 雙峰保守判準）；L3「候選隱藏次結構」非確認（其自身 interpretation：全域 ρ≈0.18）。
- cn-clean = cn≠gain（match methyl_aux 的 cn_flag=(cn=='gain')）；純 neutral 才是最嚴；LOH 仍有 LOH-unmask confound。
- 完整 subclone vs cis-ASM 天花板仍需 matched-normal DMR + single-cell（⭐4）。
- 🔴 pilot 的 `methyl_vs_ssnv_lineage.py`（Δβ approach）保留為方法備查，但**主結論用 3b 的 methyl_auxiliary+CN**。

## 4. Provenance / caveat

- pilot = L3（未 FDR/未 cis-control/單樣本）；完整版加 permutation+FDR+cis-proxy+6 樣本。
- 甲基萃取：ML qual>128 二值化（未用連續機率）；每 read 需 ≥3 CpG。
- lineage = 多-sSNV genotype 字串（未接 topology 樹的祖先-後代方向；完整版可再接 tree lineage 位置）。
- 🔴 **現階段不宣稱「甲基與 sSNV-lineage 很有關係」** — 待完整版 permutation+FDR+cis 確認。
