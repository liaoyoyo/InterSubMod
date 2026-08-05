<!--
建立時間: 2026-07-15T04:55:00+08:00
目標: 提供外部 Claude Code 對單一 FP focal-ALT 甲基多群報告的 read-only adversarial audit prompt
處理範圍: 主報告、machine-readable report dataset、frozen FP/TP/strict/topology/legacy outputs
關聯檔案: InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/20260715_單一FP_sSNV_ALT_read甲基多群與subclone假說全量驗證_01.md
-->

# External Claude Code read-only audit

你是外部、adversarial scientific reviewer。請在 `/big7_disk/liaoyoyo2001/InterSubMod` 內只讀審查，禁止修改任何檔案。

**輸出硬性要求**：這是第二次執行；第一次只在 stdout 留下一段摘要，完整 6 段審查沒有被保留。本次不得進入文件寫入/plan handoff 流程，不得說「上述內容已交付」，不得要求使用者授權寫檔。你必須把完整 6 段審查逐字放在本次 final stdout；若篇幅不足，優先保留 Critical findings 與 recalculation table，不能只給一句話總結。

## 主張問題

判斷「truth-labeled FP 單一 sSNV 的 focal-ALT reads 經 InterSubMod 無監督分成多個 methylation groups，是否可視為 high-probability subclone，或因可能 linear evolution 而支持 linear topology」。

## 必讀輸入

1. `research/20260715_single_fp_alt_multicluster_subclone_validation/20260715_單一FP_sSNV_ALT_read甲基多群與subclone假說全量驗證_01.md`
2. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/report_dataset.json`
3. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/per_sample_metrics.tsv`
4. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/strict_followup_candidates.tsv`
5. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/latest_full_v3_frozen/latest_summary.json`
6. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/latest_full_v3_frozen/latest_topology_context_summary.json`
7. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/fp_matched_tp_comparison_v1/fp_matched_tp_comparison_summary.json`
8. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/fp_matched_tp_robustness_v1/fp_matched_tp_robustness_summary.json`
9. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/ref_background_v2/ref_background_summary.json`
10. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/strict_survival_v1/strict_null_assignment_summary.json`
11. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/latest_canonical_compatibility_audit.json`
12. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/latest_input_preflight.json`
13. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/latest_tagged_subset_materialization.retry1.json`
14. `research/20260715_single_fp_alt_multicluster_subclone_validation/results/latest_tagged_subset_materialization.retry2.json`
15. `research/20260715_single_fp_alt_multicluster_subclone_validation/reviews/20260715_multiagent_consensus.md`
16. `/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md`
17. `/big8_disk/liaoyoyo2001/knowledge/05_tools/intersubmod.md`

## 強制審查項目

1. 確認主分析是否真的只用同一次最新 LongPhase-S recalibrated `FILTER=PASS` sSNV。
2. 確認 focal ALT 是否由 `alt_support=ALT`，而不是誤用 HP tag。
3. 確認新 BAM 是否真的沒有覆寫原始/canonical BAM，MM/ML 與 sidecar coverage 數字是否正確。
4. 重新計算至少 8 個關鍵數字：FP/evaluable/stable/residual/high-threshold/strict 比例、matched FP-TP high-threshold 差與 p、HCC leave-out、HP 1-2/2-2 數、topology Fisher、legacy distance identity。
5. 檢查 forced silhouette 與 stable 是否被錯畫成巢狀 funnel。
6. 檢查 matched TP 配對、residual imbalance、HCC1395 驅動與 6-biological-sample reweighting 的敘述是否充分。
7. 檢查 strict 10 是否被誤寫為 FDR/PPV、independent sites 或 confirmed clones。
8. 檢查 `phase_anchored_robust_epigenetic_candidate` 是否被過度解釋；依 KB 驗證 HP 1-1/2-1/1-2/2-2 語義。
9. 檢查一個 focal sSNV 是否可能識別 linear vs branching；區分 same-pipeline regional context 與 orthogonal confirmation。
10. 檢查 normal AD/AF、tumor REF methylation 與 matched-normal methylation 三者是否被混淆。
11. 檢查 prior InterSubMod figures/data 的使用界線；legacy distance/null semantics 是否容許直接比較 prevalence。
12. 檢查所有表格、figure caption、案例敘述與 source JSON/TSV 是否一致；列出任何無法由來源支持的句子。
13. 檢查有無遺漏會翻轉 verdict 的替代解釋或 confounder。
14. 評估主結論「只支持 strict epigenetic follow-up candidate，不支持 high-probability subclone 或 linear evolution」是否可信。

## 輸出格式

請用繁體中文，依序輸出：

1. `Overall verdict`: VALID / VALID WITH CORRECTIONS / NOT VALID。
2. `Critical findings`: 依嚴重度排序，每項附 source 路徑、實際重算數字與建議修正文字。
3. `Independent recalculation table`: 至少 8 個數字，列 report 值、source 重算值、PASS/FAIL。
4. `Claim boundary audit`: read-level heterogeneity / subclone / topology 各自可支持到哪一層。
5. `Missing evidence`: 若要升級 subclone-supported / linear topology，還缺什麼。
6. `Final wording`: 提供一段可直接放入終版的精確結論。

若沒有 critical error，必須明確寫「未發現會翻轉 verdict 的數據錯誤」。不要因為主報告寫得保守就略過數字重算。
