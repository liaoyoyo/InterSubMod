---
title: 甲基差異分析 — 3 個待驗證科學提問（交接其他 session AI）
date: 2026-06-08
type: research-question handoff brief
status: open / for-handoff
self_contained: true
related_ledger: 20260604_ISM_complete_TPFPFN_existence_cis, 20260607_meaningful_ASM_filter_iteration, 20260608_methylation_difference_pipeline_capstone
related_memory: project_ism_complete_tpfpfn_existence_cis, project_zar1l_brca2_asm_verification, project_phase2_cycle1_global_fp_filter
---

# 交接 brief：甲基差異位點分析的 3 個待驗證問題

> 給接手的 AI：本文件**自包含**。先讀「背景/已知結論」再看「3 問」。**務必先 query KB + 既有 ledger/memory，不要重做已 concluded 的方向**（見「禁區」）。

## 0. 一句話背景

我們用 unmodified ISM（`build/bin/inter_sub_mod`）跑全 6 樣本 × TP/FP/FN（~28 萬位點）找 allele-specific methylation（ASM）位點，並比較兩種「甲基差異」定義：**ISM CramersV（分群式）** vs **Δβ=HPMergedDelta（平均偏移式）**。核心發現：ASM 真實存在但**不是可用的 TP/FP filter**（弱 + 跨樣本不一致 + Δβ branch FP-enriched）。

## 1. 已知結論（必須尊重，勿推翻/重做）

| 結論 | 證據 | tier |
|------|------|------|
| ASM 廣泛存在但 **TP/FP significant 率 3.95% vs 1.07%**（subhap-matched 3.77 vs 1.09 ~3.5×）**但低 sensitivity(~4%)+COLO829 TP≈FP+FN≈TP → 非 usable filter** | ISM 完整掃描 6 樣本 | ⭐4 |
| **高 Δβ → FP 聚集** 由「無分群(OR 8.6)+LOH(OR 4.1)+假性小子群(OR 5.8)」造成，**非低覆蓋(OR 0.87)**；三者齊 TP 0.86% vs **FP 7.97%（9×）** | HCC1395 條件→FP Fisher | ⭐4 |
| Δβ branch 加入後 TP/FP 判別比 **3.7→~1.0 崩潰**（聯集捕捉多但非 TP-specific）| iteration I0-I3 | ⭐4 |
| **甲基→F1 filter DEAD**（HCC1395 in-dist ΔF1=+0.022 經 LOSO 證實 100% circularity，held-out −0.00004，HCC1954 transfer −0.377）；methylation 5th-rank vestigial | LOSO | ⭐2 L4 |
| strong-ASM **FP-enriched OR=0.194**（極端 Δβ 在 LOH/低覆蓋 regression-to-extreme）；ASM anti-discriminative | 全基因組 | concluded |

**🚫 禁區（已 concluded NEGATIVE，勿重開 without C1新數據/C2新方法/C3新前置）**：
- 甲基化當 TP/FP filter 提升 F1（DEAD ⭐2 L4，LOSO circularity）
- strong-ASM 判別 TP/FP（anti-discriminative）

## 2. 現有資料/工具（接手可直接用）

- ISM significance_summary（117 欄/位點）：`/big7_disk/liaoyoyo2001/ism_existence_scan/<sample>_{tp,fp,fn}/significance_summary.csv`（6 樣本 × 3 類）。含 CramersV、HPMergedDelta(Δβ)、4 檢定 p、Potential_LOH、HPFineN_*、NumReads、Significant、NormalBaseline_*、HP_Residual_*（HCC1395 cis）。
- 聚合 JSON：`research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/`{`ism_aggregate.json`, `iteration_history.json`, `condition_fp_consensus.json`}。
- caller-SEQC2 benchmark：各樣本 `.../longphase_s/benchmark_comparison.tsv`（3 階段 F1）。HCC1395 = ClairS 0.8443 / LongPhase-S 0.8522 / InterSubMod 0.8532。
- 樣本 FP 數（限制）：HCC1395 627 / COLO829 2232 / 其餘 8-195（FP 太少不可單獨統計）。
- 非甲基特徵可用：caller_af、coverage、LOH、CN（HCC1395 有 SEQC2）。**TO-pure 既有：caller_af 單獨 AUC≈0.654 已超越全 ISM 甲基特徵**（memory project_to_pure）。

## 3. 三個待驗證提問

### Q1：甲基是否必要，還是不需要甲基也能有效分析/區分？
**精確問法**：在判別 somatic 真實性（TP vs FP）或找「有意義差異位點」時，**甲基特徵（CramersV/Δβ）相對於非甲基特徵（caller_af + coverage + LOH + CN）有沒有 incremental 判別力**？
- 做法建議：建兩個模型 — (a) 只非甲基特徵 (b) 非甲基 + 甲基 → 比 AUC/F1 的 ΔAUC（用 within-group OLS + LOSO 防 circularity，見 `/auc-confound-guard`）。
- 預期方向（依已知）：甲基**很可能不必要**（caller_af 已 AUC 0.654，甲基 anti-discriminative）。但需**正式量化 incremental value** 確認。
- **答案標準**：ΔAUC（甲基 over 非甲基）的 LOSO 跨樣本均值 + CI；若 ≤0.01（Cohen small）→ 甲基非必要。

### Q2：這些「有意義甲基差異」位點的狀況與數量，是否足夠提升 F1？
**精確問法**：給定 meaningful-ASM 位點集（如 consensus 628 / 各 condition 組合），**有沒有任何甲基-based rule 能提升 caller F1（現 0.8522）**？最大可達 ΔF1 是多少？
- 做法建議：對每個候選 rule（移除 high-Δβ+LOH+無分群=FP-prone 的位點 OR rescue 某 FN），算 ΔF1，**必跑 LOSO**（禁 in-distribution，已知 circularity）。
- 預期方向（依已知）：**很可能不足**（filter DEAD ⭐2，TP/FP 比 1.0）。但「移除 FP-prone（高Δβ+LOH+無分群，9× FP）」這個**新角度**或許能小幅降 FP → 需量化是否 LOSO-positive。
- **答案標準**：LOSO held-out ΔF1 + Wilcoxon 跨樣本方向一致性；若 ≥+0.005 且 ≥4/6 樣本同向 → 有效。

### Q3：在其他樣本也有效嗎？
**精確問法**：Q1/Q2 的結論 + 「條件→FP 機制（無分群+LOH→FP）」是否**跨樣本復現**？
- 做法建議：條件→FP 在 **COLO829**（唯一另一 FP 夠多樣本 2232）重跑 Fisher；meaningful 位點數 + TP-specificity 在 6 樣本一致性檢查（`/multi-sample-consistency`）。
- 限制：其他 4 樣本 FP 太少（8-195）無法單獨做 FP 統計 → 只能看 meaningful 位點數 + 跨樣本方向。
- **答案標準**：COLO829 條件→FP OR 方向是否與 HCC1395 一致；meaningful 率跨 6 樣本 CV。

## 4. 交接注意

- 全程 §13：數字先寫檔→讀回→才寫報告；LOSO 防 circularity（in-dist 數字無效）。
- 升 tier 前對齊 `/scientific-rigor` §2-§7 + `/auc-confound-guard` 3-gate。
- 源頭 significance_summary.csv **不可改**（baseline 凍結）；所有新分析 post-hoc 另存。
- 完整 pipeline 說明見 `InterSubMod/docs/experiments/in_progress/2026/06/20260608_methylation_difference_pipeline_capstone_01.standalone.html`。
