<!--
建立時間: 2026-06-01
目標: 清理稽核 workflow Area [3/6] — research/ 已結案/DEAD 專案磁碟稽核（唯讀）
處理範圍: 26 個 research/ 子目錄；逐一判定 trust_tier / conclusion_status / cleanup_verdict / reclaimable_bytes
agent: cleanup-audit subagent (read-only; 絕不刪/移)
-->

# Area 3 — research/ 已結案 / DEAD 專案清理稽核

> **唯讀稽核**。任何刪除/搬移皆為 Hard Gate，由主 agent 在用戶確認後執行。
> **最高原則**：結案 NEGATIVE ≠ 自動可刪。SAFE_DELETE 只給「重跑可重生中間檔 / 明確重複副本 / 純 transient」。撐起已發佈結論/論文段落的原始證據 → 至少 ARCHIVE。不確定 → NEEDS_USER_DECISION。

## 關鍵發現（橫切所有專案）

1. **所有大型 data/ TSV 皆 UNTRACKED（已 .gitignore）**：實測 `git check-ignore` 證實 `research/*/figures/` 與 `research/*/data/` 子目錄被 .gitignore 覆蓋。大體積（>5MB）`*_master_augmented.tsv` / `bam_readqc_*.tsv` / `step1_verdict_vs_seqc2.tsv` / `*_master.tsv.gz` 全是「從 canonical master + 腳本可重生」的中間特徵檔，**結論不依賴它們**（結論存在小體積 tracked 報告/manifest/summary TSV 內）。→ 這是本區主要可回收空間來源。
2. **結論證據 = 小體積 tracked 檔**：每個結案專案的 verdict 都寫在 `manifest.yaml` / `00_INDEX.md` / `*_report.md` / `summary.tsv`（KB / 數百 KB 級），這些一律 **KEEP**（撐 INDEX / MEMORY / evidence_ledger）。
3. **G6 論文 methods 證據鏈**：`methyl_augmented_filter_phase2`（Phase2 LOSO NEGATIVE）被 docs/ 引用 **37 次**（全區最高），是 G6 論文「ISM characterization + LR sample-level negative」§3 的 methods 證據 → 整個專案 **不可 SAFE_DELETE**；只有其 78MB cycle2 中間特徵檔可 ARCHIVE。

---

## 完整 artifact 表

| path | 大小 | bytes | purpose | trust_tier | conclusion_status | cleanup_verdict | reclaimable_bytes | referenced_by | rationale |
|------|------|-------|---------|-----------|-------------------|-----------------|-------------------|---------------|-----------|
| `research/methyl_augmented_filter_phase2/` (全) | 87M | 90770393 | 甲基當 FP filter Phase2：TP rescue + cross-sample + ⭐4 升級路徑。Predecessor v6 step5 pilot ⭐3 marginal | DEPRECATED (⭐2 L4 DOWNGRADED, filter direction FAILED) | CONCLUDED_NEGATIVE | NEEDS_USER_DECISION | — | docs/ 37 處 | **G6 論文 §3 methods 證據鏈（LOSO NEGATIVE）**。整體不可刪；報告/manifest/summary KEEP，只 cycle2 大中間檔可 ARCHIVE（下行） |
| └ `…/cycle2/data/` (master_augmented TSV) | 77M | 80832855 | cycle2 per-sample `*_master_augmented.tsv`（H2009 40MB/H1437 20MB…）= 特徵注入中間檔 | DEPRECATED | TRANSIENT (regenerable) | ARCHIVE | 80832855 | none(直接) | 從 canonical master + augment 腳本可重生；wilcoxon_summary/delta_f1 小結論檔已分離保留。保守 ARCHIVE（撐 LOSO 結論的原始計算輸入，不純 transient） |
| `research/z3_internal_feature_exploration/` | 11M | 11347737 | Z3 (TO LOH+extreme AF+NG≤1) 內二階區分特徵探索 | SUPERSEDED (concluded NEGATIVE main) | CONCLUDED_NEGATIVE | NEEDS_USER_DECISION | — | docs/ 5 處 | 00_INDEX 結論 NEGATIVE（無跨樣本特徵）+ HCC1954 amplicon CONDITIONAL。報告/圖 KEEP；data/ 10.5MB 中間 TSV 可 ARCHIVE |
| └ `…/data/` | 11M | 10504460 | z3_feature_auc / af_germline_band / mechanism TSV | SUPERSEDED | TRANSIENT (regenerable) | ARCHIVE | 10504460 | none | 從 master 可重生；結論在 00_INDEX + report |
| `research/clairs_to_verdict_pilot/` (全) | 285M | 298034050 | ClairS-TO Verdict 模組 characterization（HCC1395 t20_n30 purity0.40）| SUPERSEDED (concluded NEGATIVE on F1) | CONCLUDED_NEGATIVE | NEEDS_USER_DECISION | — | docs/ 2 處 | Verdict_Germline 全落 LowQual → production FP removal=0。00_INDEX + confusion_matrix.tsv KEEP；298MB 單檔中間 TSV → ARCHIVE/SAFE_DELETE |
| └ `…/data/step1_verdict_vs_seqc2.tsv` | 284M | 297834159 | 單一 298MB per-variant Verdict×SEQC2 join 中間表 | SUPERSEDED | TRANSIENT (regenerable) | NEEDS_USER_DECISION | 297828293 | none | UNTRACKED、可由 VCF+SEQC2 重生；結論已蒸餾進 step1_confusion_matrix.tsv(1.2KB)。**最大單一可回收**；偏 SAFE_DELETE 但結案 pilot 原始 join 表保守標 NEEDS_USER_DECISION |
| `research/fp_provenance/` (全) | 58M | 59753338 | TO FP 來源溯源（caller/PON/longphase/postprocess 分層）| PRE-FIX (2026-03 5kHz TO, 早於多次 fix) | CONCLUDED_NEGATIVE | NEEDS_USER_DECISION | — | docs/ 5 處 | 結論「TO FP 無法 ISM 過濾；FN rescue 才是方向」已落 MEMORY concluded。報告/小 TSV KEEP；57MB master.tsv.gz ARCHIVE |
| └ `…/20260322_hcc1395_5khz_to/hcc1395_to_fp_provenance_master.tsv.gz` | 57M | 59700609 | 2.7M-row FP provenance master（gz）| PRE-FIX | TRANSIENT (regenerable) | ARCHIVE | 59700609 | none | UNTRACKED；pre-fix 口徑（早於 HP-fix/KDE-fix），重生需舊 binary→保守 ARCHIVE 非 SAFE_DELETE |
| `research/zone_aware_validation/` | 1.2M | 1133548 | Zone-Aware Confidence Framework H1/H3 驗證 | SUPERSEDED (characterization only) | CONCLUDED_NEGATIVE | KEEP | 0 | docs/ 4 處 | 小體積、報告+圖；Zone TP rate 差異真實但 QS 調整 NEGATIVE。撐 concept doc，全 KEEP |
| `research/vaf_alleledelta_filter_study/` (全) | 9.4M | 9742194 | VAF/AlleleDelta + CramersV 過濾策略研究（最早 2026-04-05）| SUPERSEDED | CONCLUDED_NEGATIVE | NEEDS_USER_DECISION | — | docs/ **0 處** | 8 份 .md 報告 + cramersv sites。docs/ 零引用（早期孤兒）但報告即原始證據。報告 KEEP；6.2MB data/ ARCHIVE |
| └ `…/data/` | 6.0M | 6207432 | CramersV/VAF 中間 TSV + 子集 | SUPERSEDED | TRANSIENT | ARCHIVE | 6207432 | none | 可重生中間檔；結論在 8 份 report |
| `research/hcc1954_reversal_investigation/` | 56K | 36393 | HCC1954 step3 ρ=−0.30 真反向 vs 小樣本雜訊（bootstrap）| SUPERSEDED (REJECTED=noise) | CONCLUDED_NEGATIVE | KEEP | 0 | docs/ 3 處 | 極小；manifest verdict=REJECTED(noise)。撐「5/7 明確+2/7 功效不足」修正口徑，KEEP |
| `research/F_hpfinengroups_deepening/` | 3.5M | 3387788 | HPFineNGroups 多維 confound 深度質疑（AF/NR/spatial/CN）| PRE-FIX→降級 (⭐4→⭐3) | CONCLUDED_NEGATIVE (pipeline-dependent) | KEEP | 0 | docs/ 10 處 | verdict POSITIVE_REFINED 但週報降級 pipeline-dependent。docs 10 引用 + 機制重詮釋為 phasing signature（撐 LIVE 主軸 phasing），KEEP |
| `research/hpfinengroups_saturation_check/` | 96K | 74034 | HPFineNGroups 89.1% TP rate 是訊號/NR飽和/confound（B.1）| PRE-FIX | CONCLUDED_NEGATIVE (saturation disproved) | KEEP | 0 | docs/ 2 處 | 極小；manifest 含 verdict + residualized AUC。撐 F deepening 結論鏈，KEEP |
| `research/P4_second_hit_order_pilot/` | 252K | 247032 | ISM 目標3 二次打擊順序 region-level 可區分性 | NA (pilot) | CONCLUDED_NEGATIVE | KEEP | 0 | docs/ (pilot) | manifest verdict=CONDITIONAL_NEGATIVE（summary stat 無法區分 order；需 per-read epigenotype）。小、KEEP |
| `research/P3_window_aggregation_pilot/` | 260K | 254112 | gene/window aggregation 能否突破 region-level AUC≤0.58 | NA (pilot, gating) | CONCLUDED_NEGATIVE | KEEP | 0 | docs/ (pilot) | manifest verdict=NEGATIVE（提升全來自空間 auto-correlation artifact）。撐 spatial-autocorr confound memory，小、KEEP |
| `research/coverage_multiple_validation/` | 1.7M | 1659258 | Coverage_Multiple 作 CN proxy 可靠性 | PRE-FIX (master stale, expected_cov=75 artifact) | LIVE? (status=blocked, conclusion=null) | NEEDS_USER_DECISION | — | docs/ 4 處 | ⚠ manifest status=**blocked** 非 concluded（待 master rerun + step4 oracle CN）。雖列「暫緩」但未結案；保守 NEEDS_USER_DECISION（見 anomaly A1） |
| `research/ng_kde_rescaling/` (全) | 35M | 35932303 | NG×KDE CN-tier rescaling + TP/FP discrimination（8 step）| PRE-FIX→ | UNKNOWN (無 README/verdict 檔) | NEEDS_USER_DECISION | — | docs/ **13 處** | ⚠ 無 README/manifest/00_，僅 scripts + step logs；docs 13 引用（高）。purpose 由 script 名推得。結論狀態不明 → NEEDS_USER_DECISION（見 anomaly A2）。35MB data/ 可 ARCHIVE |
| └ `…/data/` | 35M | 35646879 | step0 master + CN-tier/TPFP 中間 TSV | PRE-FIX | TRANSIENT (regenerable) | ARCHIVE | 35646879 | none | 可重生；但專案結論狀態未確認 → data ARCHIVE 而非 delete |
| `research/partB_effect_size_cn_stratification/` | 44K | 34685 | Part B 合併：B.1-3 effect size + B.2-5 cnLOH + B.2-2 CovM bimodal | PRE-FIX | CONCLUDED_NEGATIVE (MIXED/REJECTED) | KEEP | 0 | docs/ (partB) | 極小、純結論 manifest+TSV；COLO829 負向 d=−0.17 重要記錄。KEEP |
| `research/partB_high_nr_validation/` | 44K | 31353 | B.2-3 高 NR bin + NR-matched sampling 驗 LOH×AF×methyl | PRE-FIX | CONCLUDED_POSITIVE (artifact 假設 REJECTED) | KEEP | 0 | docs/ (partB) | 極小；verdict 證實效應非 NR artifact（撐 LOH×AF POSITIVE 保留結論）。KEEP |
| `research/seqc2_cnv_stratification/` (全) | 31M | 31808276 | SEQC2 CNV truth × zone/CN 分層 AUC（HCC1395 CN-tier 校準）| PRE-FIX | CONCLUDED_NEGATIVE | NEEDS_USER_DECISION | — | docs/ 10 處 | ⚠ 無 README/manifest，僅 data+figures（7 圖 6.4MB + 15 TSV）。docs 10 引用（高，撐 CovM/CN-zone 結論）。圖+小 TSV KEEP；25MB data ARCHIVE（見 anomaly A2） |
| └ `…/data/` | 25M | 25216358 | annotated CNV + zone/interaction AUC 中間 TSV | PRE-FIX | TRANSIENT (mostly regenerable) | ARCHIVE | 25216358 | none | annotated_hcc1395_cnv 等可重生；zone_auc_matrix 等小結論表保守保留評估 |
| `research/independent_analyses_20260411/` | 20K | 13470 | H2009 root-cause + PON removal rate 獨立診斷（txt/tsv）| PRE-FIX | CONCLUDED_NEGATIVE (diagnostic) | KEEP | 0 | docs/ (—) | 極小 3 檔；H2009 negative root cause 診斷紀錄。KEEP |
| `research/kde_fix_validation/` (全) | 6.6M | 6789560 | KDE expected_coverage C++ fix 後下游 4 大影響量化 | CURRENT (post-fix validation) | CONCLUDED_POSITIVE | NEEDS_USER_DECISION | — | docs/ 3 處 | README status=complete；驗證 KDE fix 下游（H-CN1/Z3/COLO829）。撐 KDE-fix 結論。報告/step summary KEEP；outputs 6.7MB 多為小 TSV，部分圖可 ARCHIVE |
| └ `…/outputs/` | 6.5M | 6714440 | step1-4 per-cn metrics + fig + impacted_conclusions.md | CURRENT | mixed (圖 regenerable, .md 結論) | NEEDS_USER_DECISION | — | none | impacted_conclusions.md(9 筆結論) 必 KEEP；圖可 ARCHIVE。混合 → 細分需用戶定 |
| `research/loh_cn_af_verification/` | 756K | 759772 | 07_LOH_CN_AF 總整理中 6 關鍵結論的數據驗證 | PRE-FIX | CONCLUDED (verification) | KEEP | 0 | docs/ 3 處 | 小；撐 research_landscape report §07。報告+圖 KEEP |
| `research/feature_layered_observation/` (全) | 440M | 459657807 | G1-G10 標準化 TP/FP 特徵分層觀察（Step1-6 方法學封裝）+ BAM readQC | PRE-FIX→部分 LIVE 方法學 | mixed (方法學 LIVE, data TRANSIENT) | NEEDS_USER_DECISION | — | docs/ 7 處 | ⚠ /feature-layered-observation **是 LIVE P3 skill**；此目錄是其首次大規模執行產物。00_main + scripts + 01-04 plan KEEP（撐 skill 方法學）。406MB data + 33MB figures 是中間檔 |
| └ `…/data/` (G4-G10 + bam_readqc) | 406M | 425081074 | per-group G4..G10 cell_delta/auc TSV + bam_readqc per-sample（H2009 14MB…）| PRE-FIX | TRANSIENT (regenerable) | ARCHIVE | 425081074 | none | UNTRACKED 大中間特徵檔；G*_auc_table/cell_delta 小結論表已在 data/ 頂層分離。保守 ARCHIVE（量大，但撐方法學示範） |
| └ `…/figures/` | 33M | 33800444 | G1-G10 分層觀察圖 | PRE-FIX | TRANSIENT (regenerable) | ARCHIVE | 33800444 | none | 可由腳本+data 重生；.gitignored |
| `research/literature_validation/` | 1.7M | 1679212 | 60+ 篇癌症甲基化文獻假說 vs ISM 實證（L1-L4）| CURRENT (literature) | CONCLUDED_NEGATIVE (L1/L2app/L4) | KEEP | 0 | docs/ 4 處 | README+4 報告+圖；L1 directional ASM NEGATIVE 等撐 ASM characterization 方向。小、KEEP |
| `research/external_tools/` | 12K | 9985 | Wakhan+SAVANA external CN/SV pilot **plan only**（pending approval）| CURRENT (plan) | LIVE (plan, not executed) | KEEP | 0 | docs/ (—) | 純 plan .md，pending user approval；關聯 LIVE phasing 主軸 CN 驗證。極小、KEEP |
| `research/igv_sessions/` | 280K | 237140 | 跨 baseline/V3F/V5/V6/paired IGV audit session XML（統一管理）| CURRENT (infra) | LIVE (audit infra) | KEEP | 0 | docs/ (—) | README 明訂「固定路徑、不要分散」；V5/V6 比對 audit 基礎設施。KEEP（含 _archive/ 已標 BROKEN 的 1 檔，屬已歸檔，留置） |
| `research/implementation_notes_skill_spec/` | 16K | 13655 | /implementation-notes skill 自身 living doc（dogfooding）| CURRENT | LIVE (skill spec) | KEEP | 0 | docs/ (—) | 極小；對應 LIVE skill 的設計紀錄。KEEP |
| `research/pre_decision_audit_skill_spec/` | 16K | 12955 | /pre-decision-audit skill 自身 meta dogfood（verdict_GO）| CURRENT | LIVE (skill spec) | KEEP | 0 | docs/ (—) | 極小；對應 LIVE skill。KEEP |
| `research/to_pipeline_staging/` (全) | 34M | 35357682 | TO pipeline 三階段 TP/FP 多模態特徵（v2 校正版；v1 已移 docs/trash）| PRE-FIX (2026-04 canonical) | CONCLUDED (characterization) | NEEDS_USER_DECISION | — | docs/ 3 處 | README 詳述 v1→v2 校正、canonical 來源。報告+圖+scripts KEEP；33MB data/ ARCHIVE |
| └ `…/data/` | 33M | 34115955 | 三階段 ISM significance_summary + 中間 TSV | PRE-FIX | TRANSIENT (regenerable) | ARCHIVE | 34115955 | none | UNTRACKED；canonical pipeline 可重生。結論在 reports/ + README |

---

## 彙總（bytes）

| verdict | 說明 | bytes |
|---------|------|-------|
| KEEP | 結論證據 / LIVE infra / skill spec / 小體積報告（含上述大專案的非-data 部分粗估 ~12MB） | ~12,000,000 |
| ARCHIVE | 可重生大中間檔（撐已發佈結論的原始計算輸入；保守不刪）| 711,233,910 |
| SAFE_DELETE | （本區無純 transient 候選達標；最接近的 clairs 298MB 因屬結案 pilot 原始 join 表，降為 NEEDS_USER_DECISION）| 0 |
| NEEDS_USER_DECISION | 結論狀態不明 / blocked / 混合 / 最大單檔 | 304,548,599 |

> 註：feature_layered_observation 因 /feature-layered-observation 為 LIVE skill，其 406MB+33MB data/figures 雖列 ARCHIVE，建議用戶確認是否保留首次示範執行產物。

## Anomalies（與先驗清單衝突 / 需注意）

- **A1 `coverage_multiple_validation` 非 concluded**：先驗清單列「降級/暫緩」，但 manifest `status: blocked`、`conclusion: null`，hypotheses verdict 仍 null（待 master rerun + oracle CN）。不是結案 → 不應視為可清理；保守 NEEDS_USER_DECISION。
- **A2 兩個高引用專案無 README/manifest/verdict 檔**：`ng_kde_rescaling`（docs 13 引用）與 `seqc2_cnv_stratification`（docs 10 引用）只有 scripts+data+figures，**無任何結論文件**。conclusion_status 無法從目錄自證（purpose 由 script 名推得）。高 docs 引用 → 撐外部結論，**勿輕動**；建議補 README 或由用戶確認結案狀態。
- **A3 `methyl_augmented_filter_phase2` 是全區最高引用（37 次）且為 G6 論文 methods**：雖 ⭐2 L4 DOWNGRADED，整體不可 SAFE_DELETE；僅 cycle2 大中間檔可 ARCHIVE。與先驗清單一致（先驗即標「至少 ARCHIVE」）。
- **A4 `vaf_alleledelta_filter_study` docs/ 零引用**：8 份報告但無人引用（早期 2026-04-05 孤兒）。報告本身是原始證據仍 KEEP，但可能候選 memory-consolidation 歸檔；非磁碟清理範疇。
- **A5 所有 data/figures 已 .gitignore**：可回收空間不影響 git 歷史；ARCHIVE = 移到冷儲存即可，無 git 牽連。
- **A6 manifest status drift**：多個 manifest 仍寫 `status: in_progress`/`initiated`（如 methyl `initiated`、clairs `in_progress`、hpfinengroups_saturation `in_progress`），但實際 verdict 欄位已填 REJECTED/結案。狀態欄與 verdict 欄不同步 — 判定時以 verdict 欄為準。
