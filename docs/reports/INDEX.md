<!--
建立時間: 2026-04-11 19:00
目標: docs/reports/ 層級報告索引 — 依研究主題分類，涵蓋所有報告檔案
處理範圍: finalized/, validated/, research_landscape/, 根目錄報告
關聯檔案:
  - docs/experiments/INDEX.md
  - docs/CURRENT_FOCUS.md
  - docs/README.md
-->

<!-- STALE-REDIRECT-BANNER (scripts/stale_redirect_banner.sh) -->
> ⚠ **此檔為 2026-05-05 前後快照，可能已過時** — 現役主軸/狀態以 `InterSubMod/docs/CURRENT_FOCUS.md` 為準（主軸已於 2026-06-11 pivot 至 Subclonal reconstruction（取代 G6））。本檔僅供歷史對照，勿據此判斷現況。


# InterSubMod 報告索引

> **最後更新**: 2026-04-11 | **報告總數**: ~110 份 | **涵蓋期間**: 2026-01 ~ 2026-04

本索引依研究主題分類，覆蓋 `docs/reports/` 下所有報告。多檔案報告以其 `00_INDEX.md` 代表。

**圖例**: 📋 finalized | 📊 validated | 🗺️ research_landscape | 📌 root-level

---

## 目錄

1. [策略總覽與階段報告](#1-策略總覽與階段報告)
2. [研究全景文件系統](#2-研究全景文件系統)
3. [週報](#3-週報)
4. [甲基化顯著性分析](#4-甲基化顯著性分析)
5. [大規模驗證與過濾策略](#5-大規模驗證與過濾策略)
6. [多樣本基線與 Confounding 診斷](#6-多樣本基線與-confounding-診斷)
7. [TO Rescue 與甲基化特徵探索](#7-to-rescue-與甲基化特徵探索)
8. [Subclone 觀察與 FP 機制](#8-subclone-觀察與-fp-機制)
9. [TO FP 來源分解](#9-to-fp-來源分解)
10. [LOH Evidence Panel](#10-loh-evidence-panel)
11. [系統性觀察 O1-O10](#11-系統性觀察-o1-o10)
12. [Self-Phasing 因果鏈與修正](#12-self-phasing-因果鏈與修正)
13. [TO 特徵深度研究與關閉](#13-to-特徵深度研究與關閉)
14. [專案重整與基礎建設](#14-專案重整與基礎建設)

---

## 1. 策略總覽與階段報告

| 日期 | 報告 | 摘要 |
|------|------|------|
| 04-08 | 📌 [InterSubMod Stage Report v1](20260408_InterSubMod_Stage_Report_v1.md) | 2025-11~2026-04 全方向階段性報告：系統性否定 + 正面成果 + 可執行後續 |
| 04-08 | 📌 [Actionable Code Tasks](20260408_Actionable_Code_Tasks.md) | 不依賴 PON-only 重跑、現在可執行的 C++ / Python 改進清單 |
| 04-03 | 📌 [研究佈局與策略架構總覽](20260403_研究佈局與策略架構總覽_01.md) | O1-O15 系統性觀察 + G1-G7 + 因果鏈 + QS + Phase 1A ML 完整推論鏈 |
| 04-04 | 📌 [研究思路驗證與推論鏈整理](20260404_研究思路驗證與推論鏈整理_01.md) | 驗證使用者對 Paired→TO 動機、FP、LOH、HP fix、PON-only 理解是否正確 |
| 03-21 | 📊 [研究進度全貌總覽](validated/2026/03/20260321_研究進度全貌總覽_01.md) | 2025-11~2026-03-17 全方向進度路線圖導覽 |
| 03-21 | 📊 [目前研究進度一頁式摘要](validated/2026/03/20260321_目前研究進度一頁式摘要_01.md) | 快速交代已完成 / 未完成 / 下一輪建議 |

---

## 2. 研究全景文件系統

> 多檔案文件，詳見 [research_landscape/00_INDEX.md](research_landscape/00_INDEX.md)

| # | 章節 | 核心問題 |
|---|------|---------|
| 01 | [TO FP 問題全貌](research_landscape/01_TO_FP問題全貌.md) | TO 30.6% FP 從何而來？PON 過濾了什麼？ |
| 02 | [Self-Phasing 根因](research_landscape/02_Self_Phasing根因.md) | TO HP tag 不可信的完整因果鏈 |
| 03 | [ISM 分析價值界定](research_landscape/03_ISM分析價值界定.md) | ISM 能做什麼？HP 依賴 vs 非依賴分類 |
| 04 | [暫停判定與重評估](research_landscape/04_暫停判定與重評估.md) | 哪些結論需修正後重測？ |
| 05 | [證據鏈總覽](research_landscape/05_證據鏈總覽.md) | 8 條完整證據鏈 |
| 06 | [結論穩定性審查](research_landscape/06_結論穩定性審查.md) | 14 個結論的穩定度評分 |

---

## 3. 週報

| 日期 | 報告 | 涵蓋期間 | 主題 |
|------|------|---------|------|
| 04-06 | 📊 [研究週報](validated/2026/04/20260406_研究週報_20260331_20260406_LOH雙定義與特徵探索全面關閉_01.md) | 03-31~04-06 | LOH 雙定義 + O11-O15 + G1-G7 NO-GO + 因果鏈 |
| 03-30 | 📊 [研究週報](validated/2026/03/20260330_研究週報_20260325_20260330_HP_bug_fix與LOH_evidence_panel_01.md) | 03-25~03-30 | HP bug fix + LOH Evidence Panel Rounds 1-4 |
| 03-24 | 📊 [研究週報](validated/2026/03/20260324_研究週報_20260318_20260324_方法學審查與TO收尾整合_01.md) | 03-18~03-24 | 方法學審查 + TO FP provenance closeout |
| 03-11 | 📊 [研究主線週報 (phase2 整合版)](validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md) | 03-05~03-11 | Paired/TO 主線 + phase 2 + annotation 回接 |
| 03-11 | 📊 [週報驗證分析報告](validated/2026/03/20260311_週報驗證分析報告_01.md) | — | 上述週報的數據驗證與邏輯審查 |
| 03-10 | 📊 [研究主線週報](validated/2026/03/20260310_研究主線週報_20260305_20260310_01.md) | 03-05~03-10 | DORADO 三路 + rescue + LongPhase-TO pilot |
| 03-12 | 📊 [AI 自動化 F1 研究主線與 PPT 藍圖](validated/2026/03/20260312_AI自動化F1研究主線與PPT藍圖整合報告_01.md) | — | AI 研究主線結論 + 16 頁 PPT 頁序設計 |

---

## 4. 甲基化顯著性分析

| 日期 | 報告 | 摘要 |
|------|------|------|
| 01-16 | 📊 [甲基化顯著性分析研究報告](validated/2026/01/20260116_甲基化顯著性分析研究報告_01.md) | 甲基化關聯顯著性的方法、限制與後續方向 |
| 01-18 | 📊 [標籤甲基化關聯顯著性分析方法](validated/2026/01/20260118_標籤甲基化關聯顯著性分析方法_01.md) | 技術架構：ReadParser → LabelTest → SignificanceAnalyzer |
| 01-18 | 📋 [Methylation Significance Report](finalized/2026/01/20260118_methylation_significance_report_01.md) | 方法與限制 + 文獻回顧 + HCC1395 輸出分析 |
| 01-19 | 📊 [多層驗證策略驗證報告](validated/2026/01/20260119_多層驗證策略驗證報告_01.md) | 多層驗證策略驗證結果 |
| 01-19 | 📊 [甲基化顯著性多層驗證實作報告](validated/2026/01/20260119_甲基化顯著性多層驗證實作報告_01.md) | LOH 子類型分類、Coverage 驗證、Quality Score 計算 |
| 01-19 | 📊 [甲基化顯著性深入診斷與多層驗證策略](validated/2026/01/20260119_甲基化顯著性深入診斷與多層驗證策略報告_01.md) | 深入診斷與策略報告 |

---

## 5. 大規模驗證與過濾策略

| 日期 | 報告 | 摘要 |
|------|------|------|
| 02-11 | 📊 [大規模驗證流程規格](validated/2026/02/20260211_大規模驗證流程規格_01.md) | 完整腳本規格、執行模式、資料管理策略 |
| 02-11 | 📊 [甲基化過濾策略綜合分析報告](validated/2026/02/20260211_甲基化過濾策略綜合分析報告_01.md) | 過濾參數意義、驗證數據、結論與改進 |
| 02-16 | 📊 [跨樣本跨純度 TP/FP/F1 綜合分析](validated/2026/02/20260216_跨樣本跨純度TP_FP_F1綜合分析報告_01.md) | 7 純樣本 + 2 組 purity + 過濾器分解 + FN 回收 |
| 02-17 | 📊 [Purity-Aware Analysis Report](validated/2026/02/20260217_purity_aware_analysis_report_01.md) | 不同 purity 下 TP/FP 差異 + QUAL 邊界 FN 回收 |

---

## 6. 多樣本基線與 Confounding 診斷

| 日期 | 報告 | 摘要 |
|------|------|------|
| 03-01 | 📊 [甲基化欄位對 TPFP 與 subclone 驗證比較](validated/2026/03/20260301_甲基化欄位對TPFP與subclone驗證比較_01.md) | 甲基化欄位在 TP/FP 與 subclone 下的比較 |
| 03-02 | 📊 [Subsample Methylation Confounding Analysis](validated/2026/03/20260302_subsample_methylation_confounding_analysis_01.md) | Subsample 驗證表現不佳的根本原因 + 可行方案 |
| 03-03 | 📊 [PMC12424993 InterSubMod 完整比較](validated/2026/03/20260303_PMC12424993_InterSubMod完整比較與相似研究整理_01.md) | 指定論文 vs InterSubMod 比較 + 相似研究檢索 |
| 03-05 | 📊 [全樣本放寬標準拆解與 F1 比較](validated/2026/03/20260305_全樣本放寬標準拆解與TPFPF1比較報告_01.md) | 7 純樣本 grid search + Shapley + Ablation |
| 03-05 | 📊 [多樣本自動化測試總覽](validated/2026/03/20260305_多樣本自動化測試總覽_01.md) | 7 樣本基本可用性確認 + chr19 距離分布基準線 |
| 03-05 | 📊 [提升幅度偏小根因與 LongPhase-TO 驗證計畫](validated/2026/03/20260305_提升幅度偏小根因與LongPhaseTO驗證計畫_01.md) | 全樣本提升幅度小的根因 + TO 純樣本可行性 |

---

## 7. TO Rescue 與甲基化特徵探索

| 日期 | 報告 | 摘要 |
|------|------|------|
| 03-08 | 📊 [HCC1395 5kHz TO Borderline Rescue 特徵證據鏈](validated/2026/03/20260308_HCC1395_5kHz_TO_borderline_rescue特徵證據鏈整理_01.md) | Rescue 關鍵特徵、規則定義、數據結果的證據鏈 |
| 03-09 | 📊 [甲基 Rescue 是否穩定有效的跨樣本判讀](validated/2026/03/20260309_甲基rescue是否穩定有效的跨樣本判讀_01.md) | 跨 HCC1395 5kHz/DORADO × paired/TO 判讀 |
| 03-11 | 📊 [HCC1395 四象限甲基 Rescue 整合報告](validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md) | 四象限 benchmark 分層 + 甲基特徵矩陣 + TO diagnostics |
| 03-12 | 📊 [Pool B FN Integration Closeout](validated/2026/03/20260312_PoolB_FN_integration_closeout_01.md) | Pool B FN caller-side rescue 正式關閉 |
| 03-12 | 📊 [TO Snapshot Scope Same-Scope Control](validated/2026/03/20260312_TO_snapshot_scope_same_scope_control_01.md) | BAM subset 是否扭曲 read-level diagnostics |
| 03-12 | 📊 [TO 雙模型可得性 Closeout](validated/2026/03/20260312_TO雙模型可得性closeout_01.md) | TO pileup vs full model 可得性問題正式關閉 |
| 03-14 | 📊 [75x 樣本放寬 GQ 與 Downstream 過濾可行性](validated/2026/03/20260314_75x樣本放寬GQ與downstream過濾可行性分析_01.md) | 放寬 GQ 依賴後段 ISM/LOH/CNV 壓回 FP 的可行性 |
| 03-14 | 📊 [GQ 救回 TP 與 FP 差異機制分析](validated/2026/03/20260314_GQ救回TP與FP差異機制分析_01.md) | GQ rescue 能救 TP 而不放回太多 FP 的機制拆解 |
| 03-15 | 📊 [四個 Universe 共同特徵交叉驗證表](validated/2026/03/20260315_四個universe共同特徵交叉驗證表_01.md) | Final kept / old full / boundary / TO 四 universe 交叉驗證 |
| 03-15 | 📊 [甲基候選研究框架與觀察資產盤點](validated/2026/03/20260315_甲基候選研究框架與觀察資產盤點_01.md) | 甲基研究資產盤點 + 候選類型 + 可驗證假設 |
| 03-17 | 📊 [跨樣本甲基特徵 TP/FP 分離觀察報告](validated/2026/03/20260317_跨樣本甲基特徵TP_FP分離觀察報告_01.md) | 7 樣本 / 16 runs Quality_Score 箱形圖有效性 |

---

## 8. Subclone 觀察與 FP 機制

| 日期 | 報告 | 摘要 |
|------|------|------|
| 03-15 | 📊 [Subclone Map Pilot](validated/2026/03/20260315_subclone_map_pilot_01.md) | 1,412 Subclone 位點的染色體分布與 Strong Count 統計 |
| 03-15 | 📊 [甲基化辨識位點分析 Pilot](validated/2026/03/20260315_甲基化辨識位點分析_pilot_01.md) | VerificationClass 分布與關鍵特徵中位數比較 (30,367 位點) |
| 03-15 | 📊 [甲基特異位點亞克隆整合分析 Pilot](validated/2026/03/20260315_甲基特異位點亞克隆整合分析_pilot_01.md) | Phase 1-3 甲基化特徵與亞克隆地圖建構 pilot |
| 03-16 | 📊 [FP 問題深度分析：SEQC2 INDEL 相鄰與 MNP 機制](validated/2026/03/20260316_FP問題深度分析_SEQC2_INDEL相鄰與MNP機制_01.md) | FP-B1 (INDEL 鄰近型) + FP-B2 (MNP 拆分型) 全面量化 |
| 03-16 | 📊 [Phase 4 甲基亞克隆視覺觀察與 FP 機制分析](validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md) | 5 TP + 4 FP 案例深度觀察與推論 |
| 03-17 | 📊 [甲基分類方法有效性觀察與改進重點](validated/2026/03/20260317_甲基分類方法有效性觀察與InterSubMod改進重點_01.md) | 9 位點系統評估 ISM 甲基分類指標 TP/FP 辨識效能 |

---

## 9. TO FP 來源分解

| 日期 | 報告 | 摘要 |
|------|------|------|
| 03-22 | 📊 [TO TP/FP Germline LOH Characterization](validated/2026/03/20260322_TO_TP_FP_germline_LOH_characterization_01.md) | TO 全量 TP 畫像 + FP Germline/LOH/Artifact 分類 |
| 03-22 | 📊 [TO FP Provenance Methylation Analysis](validated/2026/03/20260322_TO_FP_provenance_methylation_analysis_01.md) | TO FP 五層分類 + 甲基化特徵 + F1 可行性 |
| 03-22 | 📊 [TO FP 來源分解與 Paired 對照分析](validated/2026/03/20260322_TO_FP來源分解與paired對照分析_01.md) | Raw FP 來源分解 + paired oracle + TO-only 規則可行性 |
| 03-22 | 📊 [TO FP 來源分解摘要](validated/2026/03/20260322_TO_FP來源分解摘要_01.md) | 一頁式摘要：provenance + oracle + 規則驗證結論 |
| 03-22 | 📊 [Cross-Sample TO LOH FP Enrichment](validated/2026/03/20260322_cross_sample_TO_LOH_FP_enrichment_01.md) | 7 個 TO 樣本 LOH-like FP 富集與甲基化 AF 梯度 |
| 03-22 | 📊 [Cross-Sample TO ISM Gradient Analysis](validated/2026/03/20260322_cross_sample_TO_ISM_gradient_analysis_01.md) | 7 TO 樣本 ISM 甲基化特徵 AF 梯度驗證 |
| 03-22 | 📊 [TO Rescue Rules Autoresearch Report](validated/2026/03/20260322_TO_rescue_rules_autoresearch_report_01.md) | 12 假設 × 8 輪自動化研究循環結果 |
| 03-23 | 📊 [TO Residual FP Deep Dive](validated/2026/03/20260323_TO_residual_FP_deep_dive_01.md) | Raw_absent 細分 + cross-platform recurrence + paired persistent |
| 03-24 | 📊 [Paired Persistent Final FP 深度追蹤](validated/2026/03/20260324_paired_persistent_final_fp_深度追蹤_01.md) | 87+77 持續性 FP 的機制與 ISM 改善可能性 |

---

## 10. LOH Evidence Panel

### 核心報告

| 日期 | 報告 | 摘要 |
|------|------|------|
| 03-27 | 📊 [LOH Round 1 Cross-Sample Audit](validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md) | 7 paired + 7 TO 的 LOH/HP imbalance/same-locus 比較 |
| 03-27 | 📊 [LOH Round 2 Support HP0 Analysis](validated/2026/03/20260327_LOH_round2_support_hp0_analysis_01.md) | Effective_hp support 分層 + HP0 來源分析 (含 HP fix 勘誤) |
| 03-27 | 📊 [LOH Round 3 Methyl HP0 Filter](validated/2026/03/20260327_LOH_round3_methyl_hp0_filter_01.md) | HP0 Filter 驗證 + LOH×Methylation + Tier 閾值敏感度 |
| 03-28 | 📊 [LOH Evidence Panel Final Report](validated/2026/03/20260328_LOH_evidence_panel_final_report_01.md) | 四輪彙整 + 基線 F1 修正 + Phase 1 Feature 清單 |
| 04-04 | 📊 [LOH Evidence Panel Post TO HP Fix Final](validated/2026/04/20260404_LOH_evidence_panel_post_TO_HP_fix_final_report_01.md) | HP fix 後重寫的 LOH 最終結論（複寫 03-28 版） |

### 補充報告

| 日期 | 報告 | 摘要 |
|------|------|------|
| 03-29 | 📊 [LOH Round 1 v2](validated/2026/03/20260329_LOH_round1_cross_sample_audit_v2_01.md) | Fig01-06 觀察補充 + 事實修正 + 下輪建議 |
| 03-30 | 📊 [LOH Round 2 PS Export & TO Block Audit](validated/2026/03/20260330_LOH_round2_ps_export_and_to_block_audit_01.md) | PS export 邊界 + PS-block linkage 問題 |
| 03-30 | 📊 [TO LOH Enrichment Post HP Fix](validated/2026/03/20260330_TO_LOH_enrichment_post_hp_fix_01.md) | HP fix 後 TO LOH enrichment 重新評估 (7 樣本) |
| 03-30 | 📊 [Post HP Fix TO LOH Investigation](validated/2026/03/20260330_post_hp_fix_to_loh_investigation_01.md) | Enrichment 方向差異原因 + TP rescue 潛力 + Tier 重塑 |
| 04-01 | 📊 [LOH Enrichment Paired/TO Corrected Analysis](validated/2026/04/20260401_LOH_enrichment_paired_to_corrected_analysis_01.md) | 修正版 LOH enrichment — Paired vs TO 完整解讀 |

### HP Integer Tag Fix 影響

| 日期 | 報告 | 摘要 |
|------|------|------|
| 03-30 | 📊 [TO HP Integer Tag Fix 影響評估](validated/2026/03/20260330_TO_HP_integer_tag_fix影響評估與修正建議_01.md) | 7 TO 樣本影響範圍 + 重跑與修正文案建議 |
| 03-30 | 📊 [TO HP Fix 重跑驗證與 Validated 文件修正分級](validated/2026/03/20260330_TO_HP_fix重跑驗證與validated文件修正分級_01.md) | 需重跑項目新輸出驗證 + validated 報告修正優先級 |

### LOH 雙定義交叉分析 (多檔案)

> 10 章節 + 166 張圖。詳見 [00_INDEX.md](validated/2026/04/20260406_LOH雙定義交叉分析報告/00_INDEX.md)

核心判定：J9 LOH 不可作為 variant filter（確定）| J13 多特徵組合 AUC=0.577 不可行（確定）| J16 AlleleDelta 是唯一 confound-free 信號但 AUC=0.556 不足（確定）

### LOH 週報審閱 (多檔案)

> 10 份主題文件 + 3 份審查 + 1 份修正紀錄。詳見 [00_INDEX.md](validated/2026/04/20260401_LOH_weekly_review/00_INDEX.md)

---

## 11. 系統性觀察 O1-O10

| 日期 | 報告 | 摘要 |
|------|------|------|
| 04-01 | 📊 [系統性觀察 O1-O10 完整報告](validated/2026/04/20260401_comprehensive_observation_O1_O10_report_01.md) | 748K rows × 116 cols, 9 主題, 82 張圖表 |
| 04-01 | 📊 [系統性觀察 O1-O10 Summary](validated/2026/04/20260401_systematic_observation_O1_O10_summary_01.md) | Level A 結論整合摘要 |
| 04-02 | 📊 [LOH Read Threshold Visual Argument](validated/2026/04/20260402_loh_read_threshold_visual_argument_01.md) | LOH 判定缺少最低 read 門檻的影響視覺化論證 |

---

## 12. Self-Phasing 因果鏈與修正

| 日期 | 報告 | 摘要 |
|------|------|------|
| 04-02 | 📊 [LongPhase-TO vs S 因果鏈報告](validated/2026/04/20260402_longphase_to_vs_s_causal_chain_report_01.md) | 異質性變異保留策略差異 → Phasing/HP/LOH 不一致性的完整因果鏈 |
| 04-02 | 📊 [Purity-Dependent Self-Phasing Validation](validated/2026/04/20260402_purity_dependent_self_phasing_validation_01.md) | HCC1395 purity 梯度驗證 self-phasing circular dependency |
| 04-02 | 📊 [Read-Level Germline FP Research Report](validated/2026/04/20260402_read_level_germline_fp_research_report_01.md) | Site→Read-level 系統性探索 + 文獻 + pilot + 全量實驗 |
| 04-29 | 📊 [longphase-to-mod V5 Somatic Fallback PI 審核](validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md) | V5 4-commits 漸進修補；6 大判決齊備；Pass 1 only 條件下測得 |
| **05-05** | **🚨📊 [V5 Data Provenance Audit](validated/2026/05/20260505_self_phasing_V5_data_provenance_audit_01.md)** | **CRITICAL — V5 是 5-commits 不是 4；PI 報告數據為 Pass 1 only；4-30/5-01 已產 Pass 2 觸發數據待對比** |

---

## 13. TO 特徵深度研究與關閉

| 日期 | 報告 | 摘要 |
|------|------|------|
| 04-06 | 📊 [肉眼檢視推理鏈與 TP/FP 可區分性分析](validated/2026/04/20260406_肉眼檢視推理鏈與TP_FP可區分性分析_01.md) | 40 站點肉眼檢視 + 全量 7 樣本定量 + R1-R5 特徵設計局限 |
| 04-08 | 📊 [TO LOH 內外 ISM 特徵區分力完整推論鏈](validated/2026/04/20260408_TO_LOH內外ISM特徵區分力完整推論鏈報告_01.md) | TO LOH 內/外 ISM 特徵區分力的完整證據鏈 |
| 04-08 | 📊 [TO LOH 額外研究：遮罩與反轉分析](validated/2026/04/20260408_TO_LOH額外研究_遮罩與反轉分析_01.md) | Masking + reversal 補充分析 |

---

## 14. 專案重整與基礎建設

### 2026-02 專案盤點 (finalized)

| 日期 | 報告 | 摘要 |
|------|------|------|
| 02-28 | 📋 [整體報告索引與執行摘要](finalized/2026/02/20260228_InterSubMod整體報告索引與執行摘要_01.md) | 專案整體索引 (legacy) |
| 02-28 | 📋 [目標可達成性與風險評估](finalized/2026/02/20260228_InterSubMod目標可達成性與風險評估_01.md) | 目標可行性與風險評估 (legacy) |
| 02-28 | 📋 [知識庫與最新文獻對照分析](finalized/2026/02/20260228_InterSubMod知識庫與最新文獻對照分析_01.md) | 知識庫 vs 最新文獻比較 (legacy) |
| 02-28 | 📋 [研究現況與數據盤點報告](finalized/2026/02/20260228_InterSubMod研究現況與數據盤點報告_01.md) | 數據盤點 (legacy) |
| 02-28 | 📋 [專案架構內容資料結構審查與整理決策清單](finalized/2026/02/20260228_專案架構內容資料結構審查與整理決策清單_01.md) | 架構審查 (legacy) |
| 02-28 | 📋 [文件重整路徑映射索引](finalized/2026/02/20260228_文件重整路徑映射索引_01.md) | 路徑映射 (legacy) |
| 03-01 | 📋 [專案重整執行結果報告](finalized/2026/02/20260301_專案重整執行結果報告_01.md) | 重整執行結果 (legacy) |

### 2026-03 文件結構重整 (finalized + validated)

| 日期 | 報告 | 摘要 |
|------|------|------|
| 03-01 | 📋 [文件驗證基線與修復計畫](finalized/2026/03/20260301_文件驗證基線與修復計畫_01.md) | docs 重整後驗證基線 + 修復優先順序 |
| 03-03 | 📋 [docs 與專案檔案全量盤點報告](finalized/2026/03/20260303_docs與專案檔案全量盤點報告_01.md) | 全專案與 docs 的檔案盤點分類 |
| 03-04 | 📋 [scripts/tools/output 流程重整與空間遷移](finalized/2026/03/20260304_scripts_tools_output流程重整與空間遷移報告_01.md) | scripts/tools 整理 + output 空間遷移 |
| 03-01 | 📊 [AI Agent 可操作性盤點與改善方案](validated/2026/03/20260301_AI_Agent可操作性盤點與改善方案_01.md) | AI Agent 操作可行性盤點 |
| 03-01 | 📊 [output 專案外儲存策略利弊分析](validated/2026/03/20260301_output專案外儲存策略利弊分析_01.md) | 外部儲存策略比較 |
| 03-11 | 📊 [docs 架構診斷與重整規劃](validated/2026/03/20260311_docs架構診斷與重整規劃_01.md) | 結構混亂來源 + 可分批執行的重整方案 |
| 03-11 | 📊 [docs Round 1 結構重整驗證報告](validated/2026/03/20260311_docs_round1結構重整驗證報告_01.md) | Round 1 搬移內容 + 驗證結果 |
| 03-12 | 📊 [docs Round 2 結構重整驗證報告](validated/2026/03/20260312_docs_round2結構重整驗證報告_01.md) | Round 2 版本治理補強 + 路徑修正 |
| 03-12 | 📊 [docs Round 3 查詢入口與導航補強](validated/2026/03/20260312_docs_round3查詢入口與導航補強報告_01.md) | 查詢入口補強，不再依賴檔名記憶 |
| 03-12 | 📊 [presentations 一版一資料夾治理報告](validated/2026/03/20260312_presentations一版一資料夾治理報告_01.md) | presentations 輸出治理方式收斂 |

### 附屬資產

| 路徑 | 說明 |
|------|------|
| [validated/2026/03/assets/](validated/2026/03/assets/) | 03 月報告共用資產（TSV、IGV sessions、圖表，45 檔） |
