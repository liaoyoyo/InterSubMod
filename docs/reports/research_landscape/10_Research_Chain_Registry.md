---
title: Research Chain Registry（研究結論鏈註冊表）
description: 14 個月研究所有核心結論的 provenance chain + reasoning graph + R-id 連結；配合 5-axis 450-cell 組合矩陣 TSV 使用
type: report
date: 2026-04-19
status: active
tags: [registry, reasoning, provenance, stability, audit]
related: ["06_結論穩定性審查.md", "09_Part_B.md", "data/10_Registry_5axis_matrix.tsv"]
---

# Research Chain Registry — 研究結論鏈註冊表

**目的**：將 InterSubMod 14 個月研究中所有核心結論以 provenance chain 方式註冊：每條結論含前提／方法／數據佐證／推論鏈／反駁可能／R-id 追蹤連結，讓任何新方向可快速對照過去結論與質疑狀態。

**結構**：
- **5 軸組合矩陣**：TSV 在 `data/10_Registry_5axis_matrix.tsv`，240 cells（mode × zone × LOH × AF band × CN tier），標示每 cell 的 TP/FP 密度、研究狀態、效應方向。
- **結論條目**：本文件 25 個 entries，分為六個層次（Zone / LOH / AF / CN / Methylation / Heterogeneity / Meta）。

**對照資料源**：
- Master dataset：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz`（748,391 rows）
- `06_結論穩定性審查.md`（stability ratings）
- `09_Part_B.md`（Part B 質疑整合）
- `research/F_hpfinengroups_deepening/`（F pilot + Step 6B + 7 Normal BAM pilot）

---

## 統計速覽（5 軸矩陣）

| 研究狀態 | cell 數 |
|---------|--------|
| no_data（組合不存在） | 133 |
| unexplored | 54 |
| pilot_subclone（Z2 探索） | 18 |
| validated_B2_batch1（LOH × Extreme） | 14 |
| underpowered（<20 rows） | 9 |
| validated_main_signal（Z1 NonLOH） | 6 |
| concluded_Z3_pilot（Z3 LOH Extreme） | 6 |
| **Total** | **240** |

| 效應方向 | cell 數 |
|---------|--------|
| POS_high_TP (TP rate ≥0.8) | 66 |
| POS_mid (0.6-0.8) | 18 |
| neutral (0.4-0.6) | 4 |
| NEG_low_TP (<0.4) | 4 |
| NA（cell 太小） | 148 |

**Top-density cells**（n_rows >10K）見 `data/10_Registry_5axis_matrix.tsv`；高優先 unexplored cells 詳列於 §附錄 B。

---

## 結論註冊表（25 entries）

### 核心結論（14 條，對照 06_結論穩定性審查.md）

#### CL-001 Self-phasing 17.3:1 bias 真實存在

- **claim**：LongPhase-TO 在 BAM tagging 階段產生 17.3:1 haplotype assignment bias，造成 62% ISM HP_Ratio LOH artifact。
- **mode**：TO
- **layer**：Meta（基礎建構）
- **data_source**：V5 Somatic Fallback 全量重跑（`output/canonical/HCC1395/V5/...` 278GB tagged BAM）；HP tag histogram 顯示 HP1:HP2 ≈ 17:1 偏移
- **method**：HP tag distribution per sample + V3-fixed / V5 對照實驗；F1=0.7153 vs 0.7154 validated
- **reasoning_chain**：LongPhase-TO 於 haplotype 未決時預設 fallback → self-phasing 循環依賴 → HP assignment 嚴重偏斜 → 下游 HP_Ratio 特徵 62% region 觀察到虛假 LOH-like pattern
- **evidence_rating**：⭐5
- **R-id_links**：-（歷史累積，非單一 R）
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-002 Paired-pure delta F1=+0.0112（ISM FP filter 上限）

- **claim**：ClairS Paired mode 下，ISM 作為 secondary filter 的 F1 delta 上限為 +0.0112（HCC1395 SEQC2 validated）。
- **mode**：paired
- **layer**：Meta（F1 ceiling）
- **data_source**：`docs/experiments/validated/2026/01/`（Phase 1A F1 優化）
- **method**：60+ 特徵 × 748K regions RF/XGBoost；paired F1 = 0.9650 → 0.9762 (+0.0112)
- **reasoning_chain**：Paired mode 已有 ClairS 內部 filter → 甲基化加持僅能篩出極少量殘留 FP → 上限小但統計顯著
- **evidence_rating**：⭐5
- **R-id_links**：-
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-003 ISM TO FP filter NO-GO

- **claim**：60+ ISM 特徵 × 748K regions 無法提供 TO mode 可用的 FP filter（RF ceiling 0.69-0.77，甲基化全 AUC<0.58）。
- **mode**：TO
- **layer**：Meta
- **data_source**：Phase 1A F1 優化 + O1-O10 系統觀察 82 張圖表
- **method**：完整 RF/XGBoost ablation；HP-free dual path AUC=0.564；Fine-Pairwise 6 距離全無效
- **reasoning_chain**：特徵空間耗盡（`project_beyond_auc_exhaustion_confirmed.md`）→ TO mode 無法用 ISM 為 production filter
- **evidence_rating**：⭐5
- **R-id_links**：-
- **status**：concluded
- **last_reviewed**：2026-04-19

#### CL-004 HP tag 整數相容性修正

- **claim**：BAM HP tag 從 str→int 修正後，全樣本重跑 downstream 特徵皆有效。
- **mode**：both
- **layer**：Meta
- **data_source**：`project_hp_integer_tag_fix.md`
- **method**：全 7 樣本 × 3 模式 canonical baseline 重跑（19 runs）
- **reasoning_chain**：原實作 HP tag parse 類型不匹配 → 10-15% read 被丟 → 下游 ISM 特徵功效降低；修正後重跑取得當前 baseline
- **evidence_rating**：⭐4
- **R-id_links**：-
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-005 V5 Somatic Fallback（AMB 17.5→8.0%, F1=0.7154）

- **claim**：V5 somatic fallback 機制將 ambiguous read 比例從 17.5% 降至 8.0%，clean blocks 95%。
- **mode**：TO
- **layer**：Meta
- **data_source**：`project_v5_somatic_fallback_verification.md`；SEQC2 HCC1395 F1=0.7154
- **method**：V3-fixed vs V5 對照；AMB% 分布 + F1 per-sample
- **reasoning_chain**：V3 fixed 未處理 somatic-only SNV fallback → AMB 高；V5 加 somatic-only tag 回填 → AMB 大降且 F1 不退步
- **evidence_rating**：⭐4（外推性待其他樣本）
- **R-id_links**：-
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-006 LOH.bed 定義驅動 AUC=0.70 artifact

- **claim**：LOH.bed 是 "somatic-only flavored het site" 的定義，region-level LOH 特徵的 AUC=0.70 表現是 LOH.bed 定義 circular，非真實獨立訊號。
- **mode**：both
- **layer**：LOH
- **data_source**：`07_LOH_CN_AF_研究總整理.md`
- **method**：LOH.bed 生成邏輯 audit + AUC 來源分解
- **reasoning_chain**：LOH.bed 本身利用 AF + phased GT ratio → region-level LOH 特徵當然與 AF / GT 高相關 → 0.70 AUC 為 self-reference artifact
- **evidence_rating**：⭐4
- **R-id_links**：[R5]
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-007 Germline ASM >> somatic passenger SNV (3-6×)

- **claim**：Germline haplotype-associated methylation（ASM）在 somatic passenger SNV region 強度為 baseline 的 3-6 倍。
- **mode**：both
- **layer**：Methylation
- **data_source**：`project_snv_methylation_association.md`
- **method**：ASM 32-66% per-sample；per-read epigenetic context 分析
- **reasoning_chain**：Somatic SNV 周邊保留大量 germline methylation 結構 → ASM 高並非腫瘤特有，而是 germline pattern persistence
- **evidence_rating**：⭐4
- **R-id_links**：-
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-008 Beyond-AUC 特徵空間耗盡（甲基化 AUC ≤0.58）

- **claim**：純 methylation 特徵的 per-region AUC 上限 ≤0.58，特徵空間對 TO FP filter 已耗盡。
- **mode**：TO
- **layer**：Methylation
- **data_source**：`project_beyond_auc_exhaustion_confirmed.md`；2026-04-21 R1-Global Phase 2 confirmation `step8_r1_global_to_arm.md`
- **method**：60+ methylation 變異 per-region AUC 全 ≤0.58（包含 distance / entropy / pairwise）；2026-04-21 擴展至 Phase 2 Normal BAM 整合（15 features × 40,237 region）均 ≤0.533
- **reasoning_chain**：TO mode 甲基化訊號維度已全面探索 → 無新特徵組合突破 ≤0.58 上限；Phase 2 Normal BAM + Sample ASM 整合後 ceiling 仍未突破（最高 raw AUC 0.532 `Epipoly_Delta`，最高 residualized 0.533 `Epipoly_Delta`/`HPFineNGroups`），CL-025a subset 0.64-0.69 已確認為 pre-selection overfit
- **evidence_rating**：⭐4 → **⭐5**（Phase 2 R1-Global 40,237 region 驗證後升級）
- **R-id_links**：[R1, R1-Global, CL-025a, CL-025b]
- **status**：concluded (strengthened)
- **last_reviewed**：2026-04-21

#### CL-009 62% ISM HP_Ratio LOH artifact

- **claim**：62% region 的 ISM HP_Ratio ≡LOH 觀察為 self-phasing bias 驅動的虛假 LOH-like pattern（非 LOH.bed region-level artifact）。
- **mode**：TO
- **layer**：LOH
- **data_source**：V5 tagged BAM HP distribution + HP_Ratio histogram
- **method**：HP_Ratio per-region 計算 + self-phasing vs PON-only HP 分布差
- **reasoning_chain**：Self-phasing 17.3:1 偏斜 → HP_Ratio 極值區域放大 → 62% region 出現 LOH-like；但 LOH.bed 本身不受影響（R5 已驗證 Jaccard=1.0000）
- **evidence_rating**：⭐3 → ⭐4 ⬆️（R5 已排除 LOH.bed 汙染假設）
- **R-id_links**：[R5]
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-010 HPFineNGroups somatic heterogeneity marker POSITIVE

- **claim**：HPFineNGroups ≥4 + AF<0.4 + NR≥80 NonLOH filter 下 TP rate=0.9281（5/7 ≥0.85；6/7 in-scope 後 AF<0.4 stratified Cohen's h 4/7 medium+）。
- **mode**：both
- **layer**：Heterogeneity
- **data_source**：`research/F_hpfinengroups_deepening/observations/step3_filter_mechanism_cross_sample.md`；748K rows filtered
- **method**：per-sample Cohen's h (AF<0.4 stratified) + chr-shuffle Z=43.5 + CN tiers 0.90-0.94
- **reasoning_chain**：HPFineNGroups 高 → 更多 fine-resolution haplotype 分組 → 更可能為 TP（真實 heterogeneity）；NR≥80 高功效支持；AF<0.4 排除 germline contamination；NonLOH 排除 LOH penalty
- **evidence_rating**：⭐4（F pilot Step 3-6B 已通過 4 層質疑驗證）
- **R-id_links**：[Part B.1-1/2/3/4]
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-011 LOH Subclone AF × Methylation 雙證據鏈 POSITIVE

- **claim**：TO ΔNG=+0.705 / Paired ΔNG=+0.787 (Extreme vs Intermediate LOH AF)；region-level 7/7 directional (含 HCC1954)；pre-reg 排除 frac_intermediate>0.55 後 6/6 robust POSITIVE。
- **mode**：both
- **layer**：LOH / AF / Methylation（交叉層）
- **data_source**：`docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md`；B.2 批次 1（R2/R3/R4）驗證報告
- **method**：Mann-Whitney U per-sample + NR-bin stratified + cnLOH vs deletion-LOH 分層 + HCC1954 pre-reg 排除
- **reasoning_chain**：LOH AF extreme → subclone dominant；LOH AF intermediate → subclone mixture → 同一 subclone 內 methylation heterogeneity 更高 → ΔNG (inter vs intra) 正向增加
- **evidence_rating**：⭐3 → ⭐4 ⬆️（B.2 batch 1 三項質疑 + Z3 pilot mechanism + R5 LOH.bed 信任度 三獨立驗證）
- **R-id_links**：[R2, R3, R4, R5, Z3]
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-012 Per-CpG ASM CONDITIONAL POSITIVE

- **claim**：24 per-CpG ASM metrics AUC 全 ≤0.54 但 fisher_frac_sig diff=0.125 邊緣，characterization 有效、filter 無效。
- **mode**：both
- **layer**：Methylation
- **data_source**：`docs/experiments/in_progress/2026/04/20260415_PerCpG_ASM_Epiallele_Metrics_01.md`
- **method**：per-CpG Fisher exact + NME/Entropy/Epipolymorphism 24 變異統計
- **reasoning_chain**：Fisher_frac_sig 邊緣顯著但 AUC 低 → 可能是 characterization 用（per-region profiling）而非 filter 用；residualized on NR+AF+CovM 後待驗證（R8 pending）
- **evidence_rating**：⭐3
- **R-id_links**：[R8]
- **status**：active（R8 batch 3 待驗證）
- **last_reviewed**：2026-04-19

#### CL-013 Phase 2 (A+D) 程式碼完工（173/173 tests）

- **claim**：Phase 2 Normal BAM + Sample ASM + LOH BED annotation + 4-group subclone + cross-region 程式碼 2026-04-13 完工，173/173 單元測試通過。
- **mode**：both
- **layer**：Meta（infra）
- **data_source**：`docs/experiments/validated/2026/04/20260413_Phase2_integration_tests.md`
- **method**：GoogleTest 173 cases 涵蓋 dual-BAM read parser / LOH bed overlap / 4-group subclone / cross-region aggregation
- **reasoning_chain**：所有 Phase 2 單元路徑有測試覆蓋 → 程式碼正確性驗證；但生物學聲稱仍需跨樣本實驗驗證
- **evidence_rating**：⭐4（程式碼）；⭐3（生物學聲稱）
- **R-id_links**：[R1, R12]
- **status**：validated（code）；active（biology）
- **last_reviewed**：2026-04-19

#### CL-014 Phase B/C/D 三層驗證（HCC1395 single sample）

- **claim**：HCC1395 單樣本 Sample ASM 97.3% / LOH concordance 94.1% / 4-group subclone 成功識別。
- **mode**：both
- **layer**：Methylation / LOH / Heterogeneity
- **data_source**：`docs/experiments/validated/2026/04/20260414_Phase_BCD_single_sample.md`
- **method**：region-level Sample ASM calculation + LOH.bed overlap + 4-group CN1/CN2 × LOH/NonLOH subclone
- **reasoning_chain**：HCC1395 單樣本驗證 → Phase 2 pipeline 生物學可行；但**未跨樣本驗證**（R12 batch 3 pending）
- **evidence_rating**：⭐3（外推性限制）
- **R-id_links**：[R1, R12]
- **status**：active
- **last_reviewed**：2026-04-19

### F pilot 新結論（5 條）

#### CL-015 HPFineNGroups 非單調性（NG=2 germline het 污染）

- **claim**：NG=2 TP rate 0.643 為所有 NG bin 最低，係 germline heterozygous SNV 污染；NG≥4 才是純 somatic heterogeneity。
- **mode**：both
- **layer**：Heterogeneity
- **data_source**：F pilot Step 1-2；748K rows 全掃
- **method**：per-sample × per-NG TP rate curve
- **reasoning_chain**：NG=2 最常見於 heterozygous germline SNV (HP1 vs HP2 明顯兩群)；NG=3-4 則需額外 fine-level group → 更可能為真實 subclone
- **evidence_rating**：⭐4
- **R-id_links**：-
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-016 HCC1954 AF confound（pseudo-tetraploid 特例）

- **claim**：HCC1954 在 HPFineNGroups filter 下方向與其他 6 樣本相反（AF<0.2 Cohen's h=+0.775 large），因 HER2+ pseudo-tetraploid 架構 AF 分布扭曲。
- **mode**：both
- **layer**：AF / CN
- **data_source**：F pilot Step 4；HCC1954 per-sample Cohen's h
- **method**：AF<0.4 stratified per-sample 驗證 + chr-level FP enrichment 追蹤（chr5/8/17 85%）
- **reasoning_chain**：HCC1954 HER2+ 腫瘤 → CNV arm-level 重排 → AF 分布偏向 intermediate (0.1-0.4 / 0.6-0.9) → 傳統 AF<0.4 filter 意義轉變 → 與其他 6 樣本方向相反但 mechanism 已解
- **evidence_rating**：⭐4（兩獨立路徑驗證）
- **R-id_links**：[R2, Z3]
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-017 Paired mode 飽和（ΔF1 無提升空間）

- **claim**：Paired mode 在 ClairS 內部 filter 後已飽和，F pilot filter 加持 ΔF1 未顯著提升；但 ΔNG characterization 訊號仍強（+0.787）。
- **mode**：paired
- **layer**：Meta
- **data_source**：F pilot Step 2-3 paired baseline
- **method**：paired baseline F1 vs F pilot filter 應用後 F1
- **reasoning_chain**：ClairS Paired 已有 normal 對照 filter → FP 多為 edge case → F pilot filter 無餘地；但 ΔNG subclone characterization 訊號獨立於 F1
- **evidence_rating**：⭐4
- **R-id_links**：-
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-018 FP 染色體層級富集（HCC1954 chr5/8/17 70-85%）

- **claim**：HCC1954 Z3 FP 大宗集中於 chr5/8/17（85%），非 spatial hotspot 而是 arm-level CNV 架構驅動。
- **mode**：TO
- **layer**：CN / LOH
- **data_source**：F pilot Step 4；Z3 pilot §3-4
- **method**：chr-shuffle permutation test + FP per-chr enrichment
- **reasoning_chain**：HCC1954 chr8p loss + chr17p TP53 LOH + chr17q HER2 amp + chr5 複合重排 → FP 大宗落在這些 arm → blacklist whole chr5/8/17 ∩ Z3 HCC1954-only ΔF1=+0.0065
- **evidence_rating**：⭐4
- **R-id_links**：[Z3]
- **status**：concluded（HCC1954-specific mechanism）
- **last_reviewed**：2026-04-19

#### CL-019 COLO829 permanent out-of-scope（ONT R10 無 methylation）

- **claim**：COLO829 ONT R10 dataset 無 methylation basecall，F pilot 所有 methylation-based 結論排除此樣本。
- **mode**：both
- **layer**：Meta
- **data_source**：COLO829 BAM metadata audit
- **method**：MM/ML tag check + basecaller version
- **reasoning_chain**：R10 flowcell 早期 pipeline 未啟用 5mC tagging → 無 Methylation_freq → F pilot in-scope 改為 6/7 樣本
- **evidence_rating**：⭐5
- **R-id_links**：-
- **status**：concluded
- **last_reviewed**：2026-04-19

### B.2 批次 1 結論（3 條）

#### CL-020 R2 HCC1954 反向可解（pre-reg frac_intermediate>0.55 → 6/6 + 7/7 directional）

- **claim**：B.2 HCC1954 TO frac_intermediate=0.727 / Paired=0.624，其他 6 樣本 0.118-0.383；pre-reg 排除後 6/6 robust POSITIVE，未排除 7/7 directional 一致。
- **mode**：both
- **layer**：AF / LOH
- **data_source**：`docs/experiments/in_progress/2026/04/20260418_B2_LOH_Subclone_Deep_Skeptical_Check_01.md` §R2
- **method**：per-sample AF band distribution 分析 + Mann-Whitney 兩層聲稱
- **reasoning_chain**：HCC1954 CNV 架構 → AF intermediate 異常富集 → frac_intermediate>0.55 pre-reg 規則非事後合理化；CL-018 Z3 pilot 獨立復現此 mechanism
- **evidence_rating**：⭐4
- **R-id_links**：[R2, Z3, CL-016, CL-018]
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-021 R3 NR-bin 高功效非 artifact（ρ=0.483→0.709 真實）

- **claim**：B.2 LOH×AF Spearman ρ 隨 NR-bin 增強（10-30→50-80：0.483→0.709）非「高 NR 功效 artifact」，NR matched sampling 下效應保持。
- **mode**：both
- **layer**：AF / Methylation
- **data_source**：`20260418_B2_LOH_Subclone_Deep_Skeptical_Check_01.md` §R3；NR≥80 單獨 bin 驗證
- **method**：matched sampling NR≥80 down-sample 至 10-30 size → 重算 ρ
- **reasoning_chain**：若為功效 artifact → matched sampling 後效應消失；實測 matched ρ 仍在 NR 10-30 預測範圍內 → 真實效應
- **evidence_rating**：⭐4
- **R-id_links**：[R3]
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-022 R4 LOH 類型分層（cnLOH + deletion-LOH 方向一致）

- **claim**：cnLOH (CovM 1.3-2.3) 與 deletion-LOH (CovM<1.3) 的 ΔNG 方向一致（非混合掩蓋），兩層 LOH 皆支持 LOH×AF×Methylation 證據鏈。
- **mode**：both
- **layer**：CN / LOH / Methylation
- **data_source**：`20260418_B2_LOH_Subclone_Deep_Skeptical_Check_01.md` §R4
- **method**：Coverage_Multiple 分層 cnLOH vs deletion-LOH → 各層 Mann-Whitney 比較
- **reasoning_chain**：若 cnLOH 與 deletion-LOH 反向 → 原結論為 LOH 類型不均衡 artifact；實測兩層方向一致 → 原結論穩固（CovM 非精確 CN 的 caveat 納入 CL-023）
- **evidence_rating**：⭐3（CovM 非精確 CN proxy 為 known limitation）
- **R-id_links**：[R4, CL-023]
- **status**：validated
- **last_reviewed**：2026-04-19

### Z3 pilot + 批次 2 新結論（3 條 + 批次 2 R1 pilot 2 條）

#### CL-025a R1 TO arm Normal BAM 特徵 POSITIVE（HCC1395 single sample）

- **claim**：Phase 2 Normal BAM 整合後，TO arm `NormalBaseline_Coverage` AUC=0.687 / `SampleASM_Delta` AUC=0.643（F pilot subset, HCC1395）— 首次出現 Normal BAM 驅動的 filter-potential 訊號。
- **mode**：TO
- **layer**：Methylation / Meta
- **data_source**：`research/F_hpfinengroups_deepening/observations/step7_hcc1395_normal_to_pilot.md`；F pilot filter subset (743 regions; TP=665, FP=78)
- **method**：ISM binary `build/bin/inter_sub_mod` + Normal BAM region-subset；per-feature AUC (TP vs FP)
- **reasoning_chain**：原 TO-only pipeline 無 Normal BAM 維度；Phase 2 加 normal coverage + tumor/normal ASM delta → TP 集中 region 有 tumor-specific coverage anomaly 與 ASM 差異 → AUC 0.64-0.69；**須跨樣本驗證**
- **evidence_rating**：⭐3 → **⭐2 concluded NEGATIVE-for-F1-filter**（2026-04-21 R1-Global 推翻）
- **R-id_links**：[R1, R1-Global, R12, CL-008]
- **status**：**concluded (NEGATIVE for F1-filter; characterization-only retained but no POSITIVE claim)** — 2026-04-21 R1-Global 40,237 region 全域驗證後判定
- **last_reviewed**：2026-04-21
- **residualization_result_subset**：`NormalBaseline_Coverage` 0.686→0.604；`SampleASM_Delta` 0.645→0.610 (F pilot subset, n=801, FP=81)
- **R1-Global result (2026-04-21)**：`SampleASM_Delta` residualized AUC **0.527 [0.520, 0.533]** (n=40,237, TP=28,396, FP=11,841)；`NormalBaseline_Coverage` **0.512 [0.507, 0.519]**；**所有 15 個 ISM 特徵 global CI 上界均 <0.54**，CL-008 Beyond-AUC ≤0.58 ceiling 在 Phase 2 下未被突破。subset AUC 0.64-0.69 為 pre-registered selection (NG=4+AF<0.4+NR≥80+NonLOH) 的 TP/FP 分離 artifact，無法 transfer 到未選過的 region。
- **finding_doc**：`research/F_hpfinengroups_deepening/observations/step8_r1_global_to_arm.md`

#### CL-025b R1 Paired arm Per-CpG ASM Fisher CONDITIONAL POSITIVE（CL-008 challenger）

- **claim**：Paired arm `Fisher_Frac_Sig` AUC=0.726（F pilot subset, HCC1395）— **首次突破 0.70 閾值**，挑戰 CL-008 Beyond-AUC ≤0.58 結論。
- **mode**：paired
- **layer**：Methylation
- **data_source**：`research/F_hpfinengroups_deepening/observations/step7_hcc1395_normal_paired_pilot.md`；F pilot filter subset (942 regions; TP=932, FP=10)
- **method**：inter_sub_mod Per-CpG ASM Fisher exact + Phase 2 normal baseline；per-feature AUC
- **reasoning_chain**：過去 Phase 1A 未測試 Per-CpG ASM Fisher statistics → CL-008 ≤0.58 結論未涵蓋此特徵；paired mode 下 normal baseline 有效去除 germline ASM 背景，暴露 tumor-specific per-CpG 訊號；但 FP=10 小樣本警訊
- **evidence_rating**：⭐3 → **⭐2 INSUFFICIENT**（2026-04-19 殘差化後 bootstrap CI 下界 0.534 跨入隨機；FP=10 不足）
- **R-id_links**：[R1, R8, R12, CL-008, CL-012, CL-025c]
- **status**：**concluded (abandoned for F1-filter, characterization-only per user 2026-04-21)** — paired F1-filter 方向放棄；characterization-only 保留但須全域 region 驗證（F pilot subset 內 TP 99.5% 飽和，無 F1 空間）
- **last_reviewed**：2026-04-21
- **caveat**：殘差化 Fisher_Frac_Sig ~ Fisher_N_Tested + NumReads 後 AUC 0.698 (CI 0.534-0.831)；raw 0.736 (CI 0.554-0.890)；confound control 下 CI 下界仍跨入隨機 → **不能拒絕 null**
- **revision_trigger**：若未來 R12 跨樣本發現某樣本大量 FP region 且 paired Fisher_Frac_Sig 獨立於 TO SampleASM_Delta → 可重新評估 F1-filter 定位

#### CL-025c Cross-mode concordance — HCC1395 TO/paired sanity + signal independence ⭐4

- **claim**：HCC1395 TO 與 paired Phase 2 pipeline 共享同一 Normal BAM 解析完全一致（`NormalBaseline_Coverage` ρ=1.0, n=441）；TO `SampleASM_Delta` 與 paired `Fisher_Frac_Sig` 訊號大致獨立（ρ=−0.162, p=6.5e-4, n=441）；overlap region TP 近飽和（439/441=99.5%）。
- **mode**：both
- **layer**：Methylation / Meta
- **data_source**：`research/F_hpfinengroups_deepening/observations/step7_crossmode_concordance.md`；TO 801 + paired 983 → overlap 441 regions
- **method**：inner join `Chr:Pos` → per-feature cross-arm Spearman ρ + TP/FP agreement cross-tabulation
- **reasoning_chain**：若 TO/paired Normal 解析不一致 → pipeline bug；實測 ρ=1.0 → bug 排除。若 TO SampleASM_Delta 與 paired Fisher_Frac_Sig 高度相關 → paired arm 無額外 characterization 價值；實測 ρ=−0.162（顯著但極弱）→ paired arm characterization-only 有獨立價值基礎，但須全域 region 驗證。
- **evidence_rating**：⭐4（pipeline sanity 直接證實；獨立性結論為 n=441 穩定）
- **R-id_links**：[R1, CL-025a, CL-025b, CL-008]
- **status**：validated
- **last_reviewed**：2026-04-21



#### CL-023 Coverage_Multiple 非獨立 CN proxy（Opus 4.7 plan 假設 4 FALSIFIED）

- **claim**：Coverage_Multiple 作為 CN 代理受 GC bias / mappability / purity 影響，與 NGroups 共同受 "local complexity" 驅動；per-sample z-score normalize 後 HCC1954 Z3 chr5/8/17 CovM z>2 仍顯著（訊號真實）。
- **mode**：both
- **layer**：CN / Meta
- **data_source**：R6 z-score normalization report；Master dataset CovM 分布
- **method**：per-sample non-LOH region CovM mean/std → LOH region z-score → 檢查 chr5/8/17 偏高是否保留
- **reasoning_chain**：若 CovM 非 CN proxy → normalize 後 HCC1954 特定 chr 偏高應消失；實測仍顯著 → 訊號真實但 CovM 不能等同精確 CN；建議未來整合 Delly/Manta/sequenza
- **evidence_rating**：⭐4
- **R-id_links**：[R6, B.2-2]
- **status**：validated
- **last_reviewed**：2026-04-19

#### CL-024 Z3 amplicon blacklist CONDITIONAL-NEGATIVE-for-canonical ⭐5

- **claim**：Zone 3 amplicon blacklist 僅 HCC1954-only 有效（ΔF1=+0.0065），其他 5/6 樣本 collateral damage（mean ΔF1=−0.0044）；**不納入 canonical production filter**，僅作 characterization。
- **mode**：TO
- **layer**：Zone / LOH / AF
- **data_source**：Z3 pilot status=CONDITIONAL（20260419）；per-sample × strategy ΔF1 matrix
- **method**：S2 whole-chr5/8/17 ∩ Z3 blacklist per-sample ΔF1；S4 ceiling 驗證
- **reasoning_chain**：HCC1954 特殊 CNV 架構 → blacklist 命中其 FP 富集；其他樣本無此 chr-level CNV → blacklist 誤傷 TP → 跨樣本不穩
- **evidence_rating**：⭐5 NEGATIVE-for-canonical
- **R-id_links**：[Z3, CL-018]
- **status**：concluded
- **last_reviewed**：2026-04-19

#### CL-025 LOH.bed 自由於 self-phasing 汙染（Jaccard=1.0000）

- **claim**：LongPhase-TO --loh flag 使用 phased VCF region-level GT ratio（非 BAM HP tag）→ PON-only vs self-phasing 產生的 LOH.bed Jaccard=1.0000，F1=96.2% vs SEQC2。
- **mode**：TO
- **layer**：LOH / Meta
- **data_source**：R5 audit report（`20260419_LOH_bed_generation_audit_01.md`）
- **method**：PON-only vs self-phasing 兩版 haplotag run → LOH.bed diff + Jaccard 計算
- **reasoning_chain**：Self-phasing 17.3:1 bias 在 BAM tagging 層；但 LOH.bed 由 phased VCF 計算（GT ratio）→ 不依賴 BAM HP → 兩版 LOH.bed 完全一致 → 所有 LOH-dependent 結論（CL-006, CL-011）解除 haplotag 汙染疑慮
- **evidence_rating**：⭐5
- **R-id_links**：[R5, CL-006, CL-009, CL-011]
- **status**：validated
- **last_reviewed**：2026-04-19

---

## 附錄 A — R-id 交叉索引表

| R-id | 狀態 | 最新結論連結 |
|------|------|-------------|
| R1 | in-progress（batch 2 pipeline 完成，Δ_ASM 分析中） | CL-013, CL-014 → step7 pilot findings |
| R2 | validated | CL-020 |
| R3 | validated | CL-021 |
| R4 | validated | CL-022 |
| R5 | validated | CL-025 |
| R6 | validated | CL-023 |
| R7 | pending (P1 batch 3) | - |
| R8 | pending (P1 batch 3) | CL-012 |
| R9 | pending (step 6B follow-up) | - |
| R10 | pending (step 6B follow-up) | - |
| R11 | pending (step 6B follow-up) | - |
| R12 | pending (batch 3 cross-sample) | CL-013, CL-014 |
| R13 | concluded | CL-019 |
| R14 | deferred (P3 by user) | - |
| R15 | external dependency | - |
| R16 | deferred | - |
| R17 | pending (HCC1954 Normal BAM; batch 3 候選) | CL-016 |

## 附錄 B — 5 軸矩陣高優先 unexplored cells（建議探索）

前 5 大 n_rows unexplored cells（來自 `data/10_Registry_5axis_matrix.tsv`）：

| mode | zone | loh_status | af_band | cn_tier | n_rows | TP_rate | 建議研究方向 |
|------|------|-----------|---------|---------|--------|---------|------------|
| to | Z5 | NonLOH | half | CN1 | 85,725 | 0.663 | Z5 大宗 FP 殘留區；中 TP rate 值得探索 filter 補強 |
| to | Z5 | NonLOH | intermediate | CN1 | 80,268 | 0.562 | TP rate 低；TO FP filter 最大機會區 |
| paired | Z5 | NonLOH | half | CN1 | 60,254 | 0.989 | Paired 高 TP；characterization-only |
| paired | Z5 | NonLOH | intermediate | CN1 | 54,090 | 0.985 | 同上 |
| paired | Z4 | NonLOH | intermediate | CN1 | 12,866 | 0.996 | Homogeneous subclone 高 TP 驗證 |

**Z5 NonLOH（TO）為下一波 filter 探索最高潛力區**；Z5 之外的組合多數落在 no_data（no-existence）或 underpowered。

## 附錄 C — 穩定度更新流程

- 每條 CL entry 更新時同步 `06_結論穩定性審查.md`
- R-id 完成後：
  1. 更新該 R-id 對應 CL entry 的 `status` + `last_reviewed`
  2. 若提升穩定度 → 同步 06_結論穩定性審查 + MEMORY.md 相關 project note
  3. 若反駁原結論 → retracted + 新增 retraction note
- 新發現 → 先寫 finding（`docs/experiments/in_progress/`）→ 穩定後新增 CL-NNN（下一可用編號為 CL-026）

## 附錄 D — 使用方式

**研究新方向時先查**：
1. Part I：5 軸矩陣 `data/10_Registry_5axis_matrix.tsv` → 找目標 cell 的 research_status + effect_direction
2. Part II：本 Registry → 對照 cell 所在 zone/layer → 相關 CL 編號
3. Part III：`06_結論穩定性審查.md` → 深度 reasoning 與 caveat
4. Part IV：R-id 對應 finding 文件 → 具體數據與 script

**新發現記入**：
- 先寫 finding → 審閱後新增 CL entry + 更新 5 軸矩陣（cell 切換 `unexplored` → `validated/active/concluded`）

**死路標記**：
- NO-GO / NEGATIVE / concluded 方向 **必須** 保留 CL entry 避免 re-investigation（MEMORY 已覆蓋但 Registry 作 on-demand reference）
