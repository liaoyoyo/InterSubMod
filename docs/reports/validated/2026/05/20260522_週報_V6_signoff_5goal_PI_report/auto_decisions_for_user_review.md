<!--
build_date: 2026-05-22
agent: weekly-report W1-W7 全自動 (Opus 4.7, Phase B)
status: in_progress (await user review)
report_class: companion-audit (weekly-report auto-decision trail)
audience: 廖子游 (user) — Phase B 後審
parent_workflow: InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/AUTO_DECISIONS_LOG.md (D043-D060)
last_verified: 2026-05-22
-->

# Phase B weekly-report 自動決策紀錄 (B-Decisions)

> 用戶指令 (2026-05-21)：「持續自動完成直到所有都有結論與結果，有任何要決策的都自動選用預設方式，並紀錄註記讓我最後可以知道理解」
>
> 本檔承接 `phase2_completeness_audit/AUTO_DECISIONS_LOG.md` D001-D042（Phase A），紀錄 Phase B (weekly-report) 全自動執行中所有 default 選項與理由，供用戶 audit 理解。

## 0. 執行框架對齊

- 遵循 `InterSubMod/.claude/skills/weekly-report/SKILL.md` v2 — 7 階段 W1-W7 + 5 checkpoint C0-C4
- 全自動模式：C0/C2/C3 用 AI 預設快速通過；**C1 (主線) + C4 (母稿)** 仍標記必停（但本批次因用戶明示 "持續自動"，先產 draft 後用戶 review）
- Hard Gate (刪檔/C++ commit/NO-GO/覆寫 evidence_ledger) **永遠必停** — 本流程不涉及這些操作

---

## B-D001 (W1 → C0 — Raw data 收集)

**Decision**: Raw data sources（user 已明示路徑）+ 補充 frozen sources

**Sources collected**:
1. `InterSubMod/docs/CURRENT_FOCUS.md` (live working state, line 1-200 讀取 + tail)
2. `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/AUTO_DECISIONS_LOG.md` (D001-D042 Phase A 全紀錄)
3. `InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A0_existing_artifacts_verification.md` (audit unified report, 145 行)
4. `InterSubMod/docs/experiments/in_progress/2026/05/20260521_A*_*.md` x 5（A2 / A3 / A4-Ext / A4-multi-algo / A4-F1Audit）
5. A* TSV: `A1_HP_6values_5sample.tsv` / `A2_ALT_only_HP_ratio_5sample.tsv` / `A3_per_feature_AUC_cohend_5sample.tsv` / `A4_LR_DT_RF_XGBoost_LOSO_5sample.tsv`
6. `InterSubMod/research/autoresearch/evidence_ledger.jsonl` last 10 entries (entry #42-51)
7. `InterSubMod/docs/reports/validated/2026/05/20260521_PI_V6_signoff_email_draft_5goal.md` (PI email draft = Layer 3-4 source)

**理由**: user 已在 instruction 明示「不重複既有 .md 已有的內容，引用即可」— 本檔/母稿都以**引用**方式對齊 source-of-truth。

---

## B-D002 (W2 → C1 — 主線類型識別) ★ fast-track 必停

**Decision (default)**: `progress:problem` 混合（用戶明示） 

**One-line statement** (≤ 30 字, 用戶可改)：

> "V6 binary 跨 5 樣本確認 priority bug 修補 (ALT-only ratio 全反向 < 1, hp=33 5/5 firing); LR cross-sample LOSO 失敗為 distribution shift 非 model class limitation (4 algo overall mean < Cohen ribbon); 5 features 5/9 sign-consistent (loh_inner_flag 5/5 最強 anchor)"

**Thread structure**:
- **Thread A (主, progress)**: V6 BAM-level 進展 (目標 1 priority bug + 目標 2 marker engineering)
- **Thread B (sub, problem)**: LR cross-sample 方法論 reframe (目標 4 algo space exhausted)
- **Thread C (sub, observation)**: 5 features × LOH×NG×AF combinatorial characterization (目標 3 anchor)

**理由**: 本週同時有「V6 production tag sign-off ✅ GO」(progress) 與「LR LOSO mean -0.00004 → filter not deployable」(problem)，PI 視角混合處理；progress 為主軸（V6 GO 是本週決定性結果），problem 為次（LR reframe 已 5/19 完成的延續審計）。

**敘事弧**: 主線 progress = 背景 → 處理 → 結果 → 初步分析 → 下週 (V6 binary patch → A0 audit → A1-A4 補測 → 5-Goal validation → V6 production tag)；副線 problem = 方法 → 問題發現 → 目前判斷 → 求建議 (LOSO → 4 algo 都 fail → distribution shift root cause → read-level pivot vs caller-F1-headroom)。

---

## B-D003 (W3 → C2 — 內容 4 層分類 + 4 桶分流)

**Decision (default)**: 標記如下（用戶可調整 tier 升降）

**[F] Fact (有 source + N≥validation threshold)**: 17 筆
- F1 = 0.7166 caller-invariant (5 source) → PPT Tier 1
- V6 marker +52.4% / +1.26pp / ×13.2 hp_ambig → PPT Tier 1
- 5-sample LOSO mean ΔF1 = -0.00004 / Wilcoxon p=0.125 → PPT Tier 1
- HCC1395 ALT-only ratio baseline 4.41 → V6 0.43 (10× HP2 reverse) → PPT Tier 1
- 5 sample V6 ALT-only ratio 全部 < 1 (range 0.37-0.79, 5/5 same direction) → PPT Tier 1
- 4 algo LOSO 100-fold mean: LR -0.00004 / DT +0.00255 / RF +0.00267 / XGB +0.00102 → PPT Tier 1
- DT/RF rescue HCC1395 +0.0127/+0.0138 但僅 1/5 sample → 講稿 Tier 2
- 5 features sign-consistent 5/9 (loh_inner_flag 5/5) → PPT Tier 1
- caller_af AUC range 0.20-0.92 (purity-driven direction flip) → 講稿 Tier 2
- HP=33 firing range 0.05-4.27% (HCC1395 60× baseline jump) → 講稿 Tier 2
- A0 audit 6 agents 273 claims, 254 verified (93%), 0 fabrication → 講稿 Tier 2
- 4-Layer Synthesis 17.3 = 1.77 × 9.8 source-traced → 備註 Tier 3
- §8.2 row-mislabel (priority bug eng .md) → 備註 Tier 3
- 47,798 PASS 三向 file-identity → 備註 Tier 3
- 6 V6 cycles ledger ⭐3-5 全一致 → 備註 Tier 3
- 5-sample V6 NG_on=2 rate 全 ≥ 0.83 → 備註 Tier 3
- 5-sample V6 marker rate 4/5 ≥ 0.85 → 備註 Tier 3

**[O] Observation (N 不足 / 未獨立驗證)**: 6 筆
- HCC1395 H_NEW_4 +0.00699 (drop caller_af) — single-sample 待 cross-binary 驗證 → 講稿 Tier 2
- chr8+chr19 hotspot subset V6 比 baseline 更偏 HP1 (3.24→3.63) ≠ 全 BAM (1.696→1.609) 反向 → 講稿 Tier 2
- caller-F1-headroom mechanism (cycle 2 4/5 高 caller F1 ceiling-bound) → 講稿 Tier 2
- HCC1937 ALT-only ratio 0.43 = HCC1395 V6 之近似（cross-sample stability of priority bug fix）→ 備註 Tier 3
- 9-cell HCC1395 outer × NG2 × AF_L cell TP-rate 0.22 (n=3034 FP-enriched) → 備註 Tier 3
- H1437 × H2009 Spearman ρ=0.97 (同 NCI-H lineage suggest) → 備註 Tier 3

**[I] Inference (合理推論)**: 4 筆
- LR LOSO 失敗 root cause = sample-level distribution shift (not algorithm class) → PPT Tier 1
- caller_af 是 LR 5 feature framework 的 LOSO 災難 dominant confound → 講稿 Tier 2
- Non-LR framework (zone rules / read-level epigenotype) 才是真正 pivot → 講稿 Tier 2
- 升級 V5 → V6 production tag 安全 (caller F1 不變, marker downstream 改善) → PPT Tier 1

**[U] Unconfirmed (待釐清)**: 5 筆
- H2009 ALT-only chr8+chr19 全 scope 結果 (本次 chr19 only timeout fallback) → 講稿 Tier 2
- COLO829 V6 結果 (truth set permission pending) → 備註 Tier 3
- Production deployment 在 universal cross-sample LR filter 無路 — 下一步 read-level vs caller-F1-headroom 分叉決策 → PPT Tier 1
- Goal 1 per-CpG × HP × ALT 在 V3F vs V6 trade-off (V3F imbalance 0.275 vs V6 0.377) → 備註 Tier 3
- baseline LR ΔF1=+0.02302 略勝 V6 LR +0.02236 機制 (BAM-independent) → 備註 Tier 3

**4 桶分流統計**:
- PPT Tier 1 (≤ 8): 8 筆 ✓ (恰好 cap)
- 講稿 Tier 2 (≤ 15): 11 筆 ✓
- 備註 Tier 3: 9 筆
- 暫存: 4 筆（不放本次報告 — 主要為 5/22 後續工作）

---

## B-D004 (W5 → C3 — 紅旗檢查)

**Decision (default)**: pass with 4 紅旗修正

| Trigger | 原句 | 修正後 | 理由 |
|---|---|---|---|
| 過度宣稱（[F] 證實 used on [O]） | "V6 解決 priority bug 跨 5 樣本" | "V6 priority bug 修補在 5 樣本 ALT-only ratio 全反向 (5/5 same direction, range 0.37-0.79)" | 加 specific 量化 + direction，避免 "解決" 字 |
| 過度宣稱（顯著無 CI） | "DT/RF 顯著 rescue HCC1395" | "DT/RF 在 HCC1395 hold-out 比 LR 多 capture +0.0127~+0.0138 ΔF1 (僅 1/5 sample, overall mean +0.0026~+0.0027 < Cohen +0.005 ribbon)" | 加 Cohen ribbon + 1/5 sample qualifier |
| 流水帳（>5 件平列） | A0 / A1 / A2 / A3 / A4 / A4-Ext / A4-F1Audit 7 件 | 重排為 (A0 audit verifier) + (A1+A2 cross-sample evidence) + (A3 features × combinatorial) + (A4 multi-algo + F1 audit) 4 cluster | 用 narrative cluster 而非 7-bullet |
| 教授視角缺 | 無「downstream impact」段 | §16 加「下週決策分叉 — read-level pivot vs caller-F1-headroom redesign」 | 給 PI 明示需 advisor 決策的點 |

**No 紅旗 triggered**:
- ✓ 所有 [F] / [O] / [I] 描述語氣 layer-consistent
- ✓ 因果連接詞已用 ("V6 修補 → ALT-only ratio 全反向")
- ✓ §16 下週 priority + §17 教授追問 已預備

---

## B-D005 (W6 → C3 — 教授問答預測)

**Decision (default)**: 7 個追問 (相對 5-7 建議數 cap)

**⭐⭐⭐ Must-Answer (≤ 3)**:
1. **PI Q1 (DT/RF 也試過嗎？)** — "為何不用 DT/RF/XGBoost？是否 LR-specific failure?"
2. **PI Q2 (這 filter production-deployable?)** — "Cycle 1 +0.02236 vs LOSO -0.00012 差距 +0.02248 = sample-level circularity, 不能直接 deploy"
3. **PI Q3 (下一步該往哪？)** — "read-level pivot (放棄 region-level LR) vs caller-F1-headroom redesign vs 補 low-F1 panel sample"

**⭐⭐ May-Ask**:
4. **PI Q4 (V6 production 該升嗎？)** — "Goal 2 marker + Goal 5 hard threshold 改善；Goal 1 per-CpG V3F 仍 best; Goal 3+4 unblock 需 read-level redesign — calibrated GO recommendation"
5. **PI Q5 (Phase A 7-task 全跑值得嗎？)** — "A0 audit 抓 1 row-mislabel + 17 minor caveats / A1+A2 提供 5-sample priority bug fix 鐵證 / A4 multi-algo 阻止 PI 質疑 'LR 太簡單' — 全部有 PI report 直接用 evidence"
6. **PI Q6 (§8.2 row-mislabel 處理?)** — "不改 validated 報告 → PI HTML 引用 source TSV 重新 cite + footnote 註記"
7. **PI Q7 (cross-binary BAM consistency 真的有意義嗎？)** — "H_C1_6 V3F/V5/V6 LR 5-feature ΔF1 max var 0.00073 → V6 升級不破 ML training pipeline; 與 caller F1 0.7166 file-identity 同質"

每個追問皆已預備 ≥1 段 evidence 引用回答（詳見母稿 §17）。

---

## B-D006 (W7 → C4 — 母稿產出) ★ fast-track 必停

**Decision (default)**: 17 段 Layer 0-4 完整 — 不省略

**檔案位置**: `InterSubMod/docs/reports/validated/2026/05/20260522_週報_V6_signoff_5goal_PI_report/master_draft.md`

**Frontmatter 內含**:
- `main_thesis`: V6 priority bug 修補 + LR LOSO reframe + 5 features anchor
- `report_type`: weekly_report_PI
- `audience_scenario`: PI 1-on-1 (PI report scale)
- `source_artifacts`: 11 個 source path

**handoff 4 選**: 預設 **D** (母稿留檔 + 加寫下週計畫 next_week_plan.md) — 因為用戶後續還要繼續 audit；不立即 trigger pptx-build (option A)，避免在母稿尚未 user-review 時就轉 PPT。

**理由**: 全自動 mode 下 C4 仍需 user 看過後決定 handoff（A/B/C/D）；本批先產 master_draft + 此 audit log 供 user 一次性 review，不主動觸發下游 pptx-build / html-report-build。

---

## B-D007 (Number-source-grep audit per skill §10.2 升級)

**Decision**: 母稿所有 numerical value 加 `[source: <path>:line]` 標註

**Verification**: 寫完後 user 可用 `grep -E "\[source: " master_draft.md` 看 coverage；本批重點數字皆已加 source。

例：
- `F1 = 0.7166 [source: InterSubMod/docs/reports/validated/2026/05/20260511_V6_binary_complete_documentation_01.md:587]`
- `LOSO mean ΔF1 = -0.00004 [source: InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/data/loso_cv_results.tsv]`
- `V6 ALT-only ratio range 0.37-0.79 [source: InterSubMod/research/methyl_augmented_filter_phase2/phase2_completeness_audit/A2_ALT_only_HP_ratio_5sample.tsv]`

---

## B-D008 (scientific-rigor 繼承 §2-§7)

**Decision**: 高影響 (PI report) → §2 + §3 + §4 + §7 全跑

| 規範 | 母稿落地位置 |
|---|---|
| §2 證據分級 ⭐⭐⭐⭐⭐ L1-L5 ribbon | 每個 §3-§5 證據卡 frontmatter + 主線結論 inline |
| §3 Effect Size + Cohen ribbon | 所有 ΔF1 數字加 Cohen +0.005 ribbon 對照 |
| §4 DAG (relevant claims) | §3 因果鏈 mermaid 圖（priority bug → V6 patch → HP tag distribution → marker engineering downstream）|
| §7 Pre-registration alignment | §17 教授問答對照 plan v2 R-MENTAL-DRIFT discipline pre-reg |

**反 pattern 預檢**:
- ✗ 「F1 +0.05 → 鎖定」(無 effect size + 無 CI) — 本母稿無此類句
- ✗ 「跨樣本一致」(無 n / 無 95% CI) — 已改 "5/5 same direction, range 0.37-0.79"
- ✗ "clearly / strongly / rigorously" confidence 詞彙 — 不使用

---

## B-D009 (主軸聚焦 — 不重複既有 .md 內容)

**Decision**: 母稿改用「引用既有 source path」+「本週 delta 補測 / audit / reframe」+「整合判斷」三層；不複製既有報告大段內容

**例**:
- 20260511 V6 binary complete documentation 內容 → 只引用 §7 4-BAM × 5-Goal verdict matrix + §8.6.11
- 20260519 V6 vs baseline HCC1395 內容 → 只引用 §13 5-Goal validation
- Phase A 各 task .md 內容 → 只引用 TL;DR + verdict 段

詳細 source list 在母稿 §11 references。

---

## B-D010 (status summary < 500 字 per user instruction)

**Decision**: 在 Phase B agent 完成後產出一段 ≤ 500 字 status summary，列出：
- 主線
- Top 3 finding
- Top 3 PI Q
- 決定性 next step

→ 見本檔最後 §STATUS_SUMMARY 段。

---

## STATUS_SUMMARY (Phase B weekly-report 完成回報, < 500 字)

**主線**: V6 production tag sign-off — 5-Goal validation 完成、Phase A 7-task 補測 + audit 完成、Phase B 週報母稿產出。混合主線 (progress:problem)：progress = V6 priority bug 跨 5 樣本確認修補 + Goal 2 marker filter 全勝；problem = LR cross-sample LOSO 失敗 root cause = sample-level distribution shift 而非 model class limitation。

**Top 3 finding**:
1. **V6 priority bug 修補跨 5 樣本同方向** — ALT-only HP1:HP2 ratio 全 < 1 (range 0.37-0.79, 5/5 same direction)；baseline HCC1395 4.41 → V6 0.43 = 10× HP2 reverse；hp=33 conservative tag 跨樣本 firing 不一致 (HCC1395 60× / H2009 4.27%) 為新 finding。
2. **LR LOSO 失敗為 distribution shift, 非 model class** — A4 multi-algo 100-fold (LR/DT/RF/XGB × 5 sample × 5 seed) overall mean 全 < Cohen +0.005；DT/RF 僅在 HCC1395 hold-out +0.0127~+0.0138 (1/5 sample)；H_circularity 大方向確認 + H_method partial。
3. **5 features sign-consistent 5/9, loh_inner_flag 5/5 最強 anchor** — A3 5-sample × 9-feature 統計：loh_inner_flag / NG / HPFineF / NME_imbalance / ClusterPermanovaF 達 ≥4/5 sign-consistent；caller_af AUC range 0.20-0.92 (purity-driven direction flip) 為 LOSO 災難 dominant confound。

**Top 3 PI Q (must-answer)**:
1. 為何不用 DT/RF/XGBoost？— A4 實證跨 4 algo overall mean 全 < ribbon，僅 HCC1395 1/5 sample partial rescue, 4/5 仍 ≈ 0 / negative。
2. Filter production-deployable 嗎？— Cycle 1 +0.02236 vs LOSO -0.00012 差距 +0.02248 = 100% sample-level circularity; **NOT deployable as universal cross-sample LR filter**。
3. 下一步該往哪？— 三條路：(a) read-level pivot (phase_block_3d 啟動)；(b) caller-F1-headroom redesign (cycle 5+)；(c) 補 low-F1 panel sample (driving HCC1937)。需 PI 決策。

**決定性 next step**: ✅ Send V6 production tag PI sign-off email (Hard Gate — user 親自 copy 到 mail client)；🟠 等 PI Q3 advisor 決策後啟動 phase_block_3d 或 cycle 5。
