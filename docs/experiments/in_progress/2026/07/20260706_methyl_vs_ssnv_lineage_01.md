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

> **版本** v0.2 (2026-07-06) · **狀態：轉向用既有 methyl_auxiliary + SAVANA CN、3 樣本裁決完成** · tier：CN=**L1**、甲基 residual=**L3**
> **一句話**：甲基群內雙峰(超越 HP)確實存在(5–17%)，**但我的 SAVANA CN 揭露 2/3 樣本 77–86% 是 CN-gain(擴增)假象**；扣掉後 cn-clean 真訊號稀少(純 neutral H1437 僅 1)、L3 候選 → 甲基-lineage 關聯**重度 CN-confound**、bulk 不能證 subclone。與主軸一致。

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

## 4. 裁決（L1 CN + L3 甲基）

- **甲基群內雙峰(超越 HP)確實存在**（5–17% genotype 群）→ 甲基與 sSNV-lineage 軸**有關聯**。
- 🔴 **但 2/3 樣本 77–86% 落 CN-gain 區 = 擴增 multiplicity 假象（就是 CNV）** — 之前 cn=unknown 看不到，**我的 SAVANA CN 才揭露**。H2009 較少(30%)因其 gain 區少。
- **扣 gain 後 cn-clean residual 50/411/41（Δβ 強 0.25–0.33）= 超越 CN 的候選**，但**多在 LOH(另有 confound)、純 neutral 極少（H1437 僅 1）** → **L3 候選、非確認 subclone**。
- **總結**：甲基-lineage 關聯部分成立但**重度 CN-gain confound**；cn-clean 真訊號稀少且 L3 → 與主軸「甲基=有界輔助、bulk 不能證 subclone」一致，現有 **CN 定量證據**。

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
