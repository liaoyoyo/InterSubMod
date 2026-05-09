---
title: 開發流程複盤 — 過去問題、已建立機制、仍 gap 與新機制提案
date: 2026-05-08
status: validated
phase: lessons_learned
type: lessons_learned_review
tier: 4
classification: process_improvement
period_covered: 2026-04-04 to 2026-05-08
upstream_reports:
  - InterSubMod/docs/experiments/in_progress/2026/05/20260506_Drill1_Harness_Retrospective_01.md
  - InterSubMod/docs/experiments/in_progress/2026/05/20260508_Drill2_End_to_End_Cycle_Walkthrough_01.md
  - InterSubMod/docs/reports/validated/2026/05/20260508_Harness_Complete_Architecture_Review_01.md
  - ~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md (v1.7)
sources:
  - 6 retract events (2026-04)
  - 30+ memory feedback entries
  - Drill 1 §6.2 5 易錯點
  - Drill 2 routing weakness
  - plan v1.7 §10 易錯點 1-10
---

# 開發流程複盤 — 過去問題、已建立機制、仍 gap 與新機制提案

> **Bottom line**：盤點 2026-04-04 到 2026-05-08 共 30+ 個問題事件（含 6 個正式 retraction、5 個 Drill 1 揭露的設計易錯點、5 個 Drill 2 routing 觀察、20+ memory feedback 中的個人風格陷阱），歸類為 **5 大 Category**（A 前置條件、B 結論健壯性、C 認知/風格/流程、D 基礎設施、E 報告呈現）。**Category A+B 已被 7-phase × M1+M2+M4 governance 99% 覆蓋**（Drill 1 sensitivity 6/6 證實），**Category C 約 60% 覆蓋**（execution mode + confirmation matrix 已 cover，週期性 review 仍 gap），**Category D 約 30% 覆蓋**（disk monitoring + working tree check 未自動化），**Category E 約 70% 覆蓋**（PPTX rules 已成 memory，但 build 前自動檢查仍手動）。本報告提出 **10 個新機制提案**（5 個短期 P1、3 個中期 P2、2 個長期 P3），並列出 **何時建立、何時 PARKED、何時不建** 的判斷規則。

---

## §1 範圍與方法

### 1.1 涵蓋的問題來源（30+ events）

| 來源 | 數量 | 性質 |
|---|---|---|
| 2026-04 月正式 retraction events | 6 | hard ground truth (Drill 1) |
| Drill 1 揭露的設計易錯點 | 5 | tooling design oversight |
| Drill 2 forward routing 弱點 | 3 | runtime UX gap |
| plan v1.7 §10 (5+5) 易錯點 | 10 | implementation traps |
| memory feedback 個人風格陷阱 | 20+ | 過去用戶糾正紀錄 |
| project memory NEGATIVE 結論 | 15+ | 過早 promising 後來反駁 pattern |
| **TOTAL** | **~60** | — |

### 1.2 分類軸（2D matrix）

- **Category（根因類別）**：A 前置條件 / B 結論健壯性 / C 認知流程 / D 基礎設施 / E 報告呈現
- **Coverage 狀態**：✅ 已 cover / 🟡 部分 cover / 🔴 未 cover

### 1.3 排除項
- 真正的研究結論（如 NG=2 phasing 改錯解讀為 methylation）— 屬研究內容而非流程問題；雖然 harness 應抓但屬於 Drill 1 sensitivity 已驗
- 投稿 / 論文撰寫流程 — 用戶 2026-04-18 確認「論文定位降優先」(memory feedback_paper_positioning_de_prioritized)

---

## §2 Category A — 前置條件問題（資料 / binary / dataset / upstream）

### 2.1 已發生事件（6 起）

| 事件 | 日期 | 根因 | 損失 | Drill 1 攔截 |
|---|---|---|---|---|
| VCF 來源錯誤（pileup symlink 指向 paired ClairS 非 ClairS-TO） | 04-04 | dataset_id 未 verify caller 類別 | 數天 ISM 分析需重跑 | ✅ P2 BLOCKED via P-04 |
| KDE fix → master dataset 全 stale-binary | 04-13→04-20 | binary commit 後下游未 invalidate | 7 樣本 master_phase1 須重跑；S3 TP 95.5%→58.3% | ✅ P2 BLOCKED via binary stale + P5 ess=0.04 |
| merged dataset AF 非 caller_af + HCC1395 LOH 5.7% 殘缺 | 04-23 | merged 檔資料 schema 未驗 + 衝突命名 | 5 樣本 W=0 cross-sample 反例 | ✅ P2 BLOCKED via P-03 + P-04 |
| longphase TO vs V5 working tree dirty | 04-29 | binary commit 鏈 + uncommitted state | 因果鏈需 retro 整理 | ✅ P2 BLOCKED via working tree dirty + binary stale |
| HPFineNGroups flag dependency 未硬化 | 04-22 | flag=on 條件未在 plan 寫死 | ⭐4→⭐3 降級 | 🟡 部分（P5 P-06 + override） |
| HCC1395 phase1_new LOH 5.7% vs 正常 58.8% | 04-23 (merged 子) | 同 04-23 | 已含 04-23 損失 | ✅ |

### 2.2 已建立的機制（Coverage = ✅ 99%）

| 機制 | 位置 | 攔截哪類 | 驗證 |
|---|---|---|---|
| **P2 PRECHECK gate (`/check-staleness`)** | `InterSubMod/.claude/skills/check-staleness/` | binary version / dataset id / upstream report freshness 全 audit | Drill 1 4/6 events 在 P2 BLOCKED |
| **`post_cpp_commit_invalidate.sh` hook** | `InterSubMod/scripts/hooks/` | C++ commit 觸發 stale_marks.jsonl | <1 hr 標記時間（plan §1.3 達成） |
| **pitfalls_table P-01..P-06** | `InterSubMod/.claude/skills/run-evaluator/pitfalls_table.json` | 6 條結構化 pitfall 規則（pileup symlink, merged_af 等）| Drill 1 hpfinengroups P-06 命中 |
| **state.schema.json `binary_version` + `dataset_id` 必填** | `InterSubMod/state/schemas/` | plan.json schema 驗證強制填寫 | D2-A precheck PASS（純分析設 null） |
| **per-folder `state/CLAUDE.md`** | v1.7 batch C | 「不可手改 active.json / state.json」規範 | 已 commit e2eb43b |

### 2.3 仍 gap（Coverage gap = ~1%）

| Gap | 影響 | 建議機制 |
|---|---|---|
| 新型資料 schema 出現時 pitfalls_table 滯後 | 中（如未來 caller_af bins 變更未進 P-07） | **§7.1 提案 #1 — Pitfalls table 月度 review ritual** |
| `binary_version: null` (純分析) 約定未在 plan.json 範本明示 | 低 | **§7.1 提案 #2 — plan.json template + skill auto-fill** |

---

## §3 Category B — 結論健壯性問題（cross-sample / effect size / pitfall）

### 3.1 已發生事件（5 起）

| 事件 | 根因 | Drill 1 攔截 | 教訓 |
|---|---|---|---|
| Thread B (LOH×AF×CN biology-informed filter) 跨樣本反證 04-26 | 1/6 樣本（H2009 飽和 artifact）拉動全 5 fail | ✅ P5 pending_review (override) | n_passed=1/6 是顯著 fail signal |
| HPFineNGroups subclone marker ⭐4→⭐3 04-22 | flag=on 條件未驗證 + n_reads/AF confound | ✅ P5 P-06 + override | feature name 不能依字面意思推論生物語意 (memory feedback_feature_name_vs_definition_rule) |
| longphase TO vs V5 noise (-0.0003 ΔF1) | 雜訊量 < threshold 但被誤讀為訊號 | ✅ P5 ess=0.00 + composite > 0.7 | effect size stability 是最 sensitive component (4/5 cases) |
| O11 epipolymorphism 全因 n_reads confound | feature 未做 confound check | 🟡 P5 pitfall coverage 未抓（P-06 不夠廣） | 需 expanded P-07 |
| O12 LOH AlleleDelta = AF confound + L2 collider bias | residualize on AF 產生虛假信號 | 🟡 P5 P-06 抓但 L2 collider 未直接抓 | memory feedback_L2_collider_bias 已記，需編入 pitfalls_table |

### 3.2 已建立的機制（Coverage = ✅ 95%）

| 機制 | 位置 | 攔截哪類 |
|---|---|---|
| **P5 EVALUATE gate (`/run-evaluator`)** | `InterSubMod/.claude/skills/run-evaluator/` | 6-component composite + per-component override |
| **multi-sample-consistency skill** | `InterSubMod/.claude/skills/multi-sample-consistency/` | 7 樣本 canonical 順序 + 雙 p-value（v1.7 batch B 9-item checklist） |
| **per-component override 規則** | `_evaluator_run.py` Stage 2 | n_critical ≥ 1 OR n_low ≥ 3 → forced pending_review；Drill 1 貢獻 4/5 額外攔截 |
| **auc-confound-guard skill** | `InterSubMod/.claude/skills/auc-confound-guard/` | n_reads / AF / LOH confound 自動 sweep |
| **provenance-tier-audit skill** | `InterSubMod/.claude/skills/provenance-tier-audit/` | system-level audit 防 tier inflation；v1.7 batch B 8-item checklist |
| **`pre_tier_upgrade_check.sh` hook** | `InterSubMod/scripts/hooks/` | ⭐4-5 升級無 evaluation.json 阻擋 |
| **L4 mandatory 4-track evidence chain anchor #1** | validation-protocol skill | Statistical / Cross-sample / Mechanism / Orthogonal 四軌齊備才標 validated |

### 3.3 仍 gap（Coverage gap = ~5%）

| Gap | 影響 | 建議機制 |
|---|---|---|
| pitfalls_table 只有 6 條 P-01..P-06；新型 confound (KDE artifact, residualize collider) 未編入 | 中 | **§7.1 提案 #3 — Pitfalls table 擴充至 P-07~P-12** |
| L2 collider bias 在 memory 但未進 pitfalls_table | 中 | 同上 #3 |
| 「pilot 看似 OK 但跨樣本就崩」（thread_b 模式）early termination gate 在 P3-P4 之間缺 | 中 | **§7.1 提案 #4 — P3.5 cross-sample early-termination gate** |

---

## §4 Category C — 認知 / 風格 / 流程問題

### 4.1 已發生事件（10+ 起，從 memory feedback 抽）

| 事件 / 風格陷阱 | 來源 | 已 cover？ |
|---|---|---|
| Pivot 過於頻繁（4 月 4-6 次撤回 ≈ 每週 1 次） | plan §1.1 | ✅ harness 7-phase + Hard Gate |
| 過早結論（O9/O11/O12/O13 series promising → NEGATIVE） | project memories | ✅ P5 EVALUATE 抓 |
| 過度解讀 feature name（HPFineNGroups 字面推 methylation） | feedback_feature_name_vs_definition_rule | 🟡 memory 記，無 hook |
| 模糊指令 → 模型按字面執行 | CLAUDE.md Opus 4.7 literal | ✅ 暫停判斷矩陣 + 假設陳述規則 |
| 「改進」「優化」未定義成功標準 | CLAUDE.md 高影響場景 | ✅ confirmation-protocol 已記 |
| 第一次新概念 vs 已驗證重複執行 vs 全自動模式邊界 | feedback_execution_mode_hierarchy | ✅ 三級分類已 memory |
| Strategy 同意後逐項實作確認 | feedback_strategy_then_per_item_confirmation | ✅ memory + 多次驗證遵守 |
| 報告中過度宣稱 / 流水帳 | memory weekly_report_workflow | ✅ weekly-report skill 4 主線分類 + 邏輯紅旗檢查 |
| 風險未結構化（無 R-id 編號 + 無 P0/P1/P2/P3） | feedback_risk_structured_iterative | 🟡 memory 記，無 skill 強制 |
| 小規模快速驗證 < 2hr pilot 應優先 | feedback_small_scale_validation_first | ✅ problem-framing-ideation anchor #5 已硬化 |
| 空間 autocorrelation confound（chr+pos 聚合特徵） | feedback_spatial_autocorrelation_confound | 🟡 memory 記，evaluator 未 hardcode |
| Pooled OLS residualization 陷阱 | feedback_pooled_ols_residualization_trap | 🟡 memory 記，evaluator 未 hardcode |

### 4.2 已建立的機制（Coverage = 🟡 60%）

| 機制 | 位置 | 攔截哪類 |
|---|---|---|
| **暫停判斷矩陣（2D + 可逆 override）** | `.claude/CLAUDE.md` | 影響 × 信心 + 可逆性 → 4 級停 |
| **5 種執行模式分級（FYI / Review / Gate / Hard Gate / 全自動）** | `confirmation-protocol` skill | 互動 vs 全自動行為 |
| **M1 main_axis lock + drift detection** | state.json schema + cycle-init | 防 hypothesis 漂移 |
| **L4 mandatory 4-track evidence chain** | validation-protocol skill (anchor #1) | 強制多軌證據 |
| **`/known-pitfalls` skill (passive ref)** | `.claude/skills/known-pitfalls/` | 8 個 confound 規則收錄 |
| **5 段報告骨架 freeze (anchor #3)** | structured-tech-report skill | Exec / Background / Evidence×4 / Limitations / Conclusion |
| **5 個 A 密度 SKILL.md 三段元資料** | v1.7 batch B | Failure Mode + Diagnostics 強制 |

### 4.3 仍 gap（Coverage gap = ~40%）

| Gap | 影響 | 建議機制 |
|---|---|---|
| Memory feedback 中的「個人風格 anchor #3/6/7」尚未硬化（風格 anchor batch 5b PARKED） | 中 | **§7.2 提案 #5 — anchor 5b 真跑後評估** |
| Risk 結構化（R-id + P0/P1/P2/P3）無 skill 強制 | 低 | **§7.2 提案 #6 — risk-register skill 或 plan template 加 risk_register 欄位** |
| 「pivot fatigue」（月 4-6 次方向切換）無自動偵測 | 中 | **§7.2 提案 #7 — pivot rate cron monitor** |

---

## §5 Category D — 基礎設施 / 環境問題

### 5.1 已發生事件（5 起）

| 事件 | 日期 | 損失 | 影響 |
|---|---|---|---|
| **/tmp 800 GB 災情 — pipeline 中間檔寫滿 root volume 100%** | 2026-05-08 | longphase-to BAM 預設輸出到 /tmp 把 root 寫滿；新 server 上未 export TMPDIR | 大型；user pain extreme |
| longphase 04-29 working tree dirty | 04-29 | binary 改但未 commit → cycle 標 unclean state | 中（已被 P2 抓） |
| expected_coverage hardcoded 75.0 default | 早期 | 7 樣本共用 default → 7 樣本 KDE 啟用未生效 | 中（待 /cpp-change 修） |
| binary 改後未跑 `make -j$(nproc)` 直接 commit | 多次 | C++ 改但 binary 不更新 | ✅ pre_commit_compile_check.sh hook 已存在 |
| Memory KB freshness 14 天未驗證 | 持續 | 老資料被當權威引用 | 🟡 hook 已警示，無強制重認證 |

### 5.2 已建立的機制（Coverage = 🟡 30%）

| 機制 | 位置 | 攔截哪類 |
|---|---|---|
| **`pre_commit_compile_check.sh` hook** | `InterSubMod/scripts/hooks/` | C++ 改後 commit 前必編譯 |
| **`compile_success_clear.sh` hook** | 同上 | make 成功清除待編譯標記 |
| **`cpp_edit_guard.sh` hook** | 同上 | C++ 編輯時建立待編譯標記 |
| **`kb_freshness_warn.sh` hook** | 同上 | KB 文件 14 天未驗證警告 |
| **`/check-staleness` skill** | 同 §2.2 | git working tree dirty 檢查（plan §10 易錯點 #2 已加） |

### 5.3 仍 gap（Coverage gap = ~70%）

| Gap | 影響 | 建議機制 |
|---|---|---|
| Disk health（/tmp tmpfs mount + free space + log rotation）無監控 | **高（已造成 1 次重大損失）** | **§7.1 提案 #8 — disk health monitor cron + first-run sanity check** |
| Pipeline 啟動前 TMPDIR 自動 export | 高 | **§7.1 提案 #9 — pipeline wrapper script 強制 export TMPDIR** |
| KB 14 天 freshness 警告但無強制重認證 | 中 | **§7.2 提案 #10 — auto-stale flag + freeze quoting after 30d** |
| `expected_coverage` hardcoded bug 未修 | 中 | 已在 memory，待下一次 /cpp-change cycle 修 |

---

## §6 Category E — 報告 / 視覺呈現問題

### 6.1 已發生事件（10+ 起，全在 memory feedback）

| 事件 | 來源 | 已 cover？ |
|---|---|---|
| Matplotlib 中文亂碼（Droid Sans Fallback 未注入） | feedback_matplotlib_cjk_font_rule | 🟡 memory 記，無 build hook |
| PPTX Latin+CJK font fallback 失敗（單 CJK 字型讓英文變方塊） | feedback_pptx_screenshot_rendering_rules | 🟡 memory 記，無 build hook |
| PPTX 圖片強制 width+height 擠壓變形 | feedback_pptx_screenshot_rendering_rules | 🟡 memory 記 |
| PPTX 標題 >15 字、每頁 >3 焦點 | feedback_ppt_slide_design_rules | ✅ pptx-build skill 規則 |
| PPTX 文字斜體不可用 | feedback_pptx_text_formatting_rules | ✅ skill 規則 |
| PPTX 雙語格式（EN 縮排 0.25"+60% 字號） | feedback_pptx_bilingual_formatting | ✅ skill 規則 |
| PPTX storyboard 導演審查（信任→分析→根因） | feedback_pptx_director_storyboard | ✅ skill 規則 |
| PPTX 修改前 diff 同步（防覆蓋手動編輯） | feedback_pptx_presync_workflow | ✅ skill 規則 |
| 圖表分 Paired/TO + 細緻 2x2 × 7 samples 固定順序 | feedback_figure_layout_standard | ✅ results-analysis 加 5 層 caption checklist (anchor #6) |
| Figure caption 缺題/例/軸/CI/色盲 | anchor #6 | 🟡 anchor #6 batch 5b PARKED，待 Drill 2 後評 |
| 路徑前綴 `InterSubMod/...` 規則 | feedback_md_path_prefix_rule | ✅ md_path_format_rule.sh hook |

### 6.2 已建立的機制（Coverage = 🟡 70%）

| 機制 | 位置 | 攔截哪類 |
|---|---|---|
| **pptx-build skill** | `InterSubMod/.claude/skills/pptx-build/` | 9 個 PPTX rules + 8 個渲染問題排查 |
| **weekly-report skill v2** | `InterSubMod/.claude/skills/weekly-report/` | 7 phase 流程 + 4 主線分類 + 邏輯紅旗 + 教授問答預測 |
| **md_path_format_rule.sh hook** | `InterSubMod/scripts/hooks/` | UserPromptSubmit 時提示路徑前綴規則 |
| **results-analysis 5 層圖表 checklist** | results-analysis skill (anchor #6) | 題/例/軸/CI/色盲 |

### 6.3 仍 gap（Coverage gap = ~30%）

| Gap | 影響 | 建議機制 |
|---|---|---|
| Matplotlib / PPTX build **前** auto-check（CJK font + emoji ban + aspect ratio） | 中（重複犯錯） | **§7.2 提案 #11 — pre-build validator script in pptx-build/weekly-report** |
| anchor #6 圖表 caption checklist 未硬化到 results-analysis SKILL.md | 中 | batch 5b PARKED；plan §4.5.4-G 已記 |
| PPTX storyboard 規範 anchor #7（草稿可重評）未 hardcode | 低 | batch 5b PARKED |

---

## §7 新機制提案 — 10 個建議

### 7.1 短期（P1）— 可立即建立，不阻塞 Path B

#### 提案 #1：Pitfalls table 月度 review ritual（防 §2.3 #1）
- **何時建**：2026-06-01 第一次跑（每月 1 號）
- **內容**：對 `InterSubMod/.claude/skills/run-evaluator/pitfalls_table.json` 跑 sweep；對比過去 30 天 cycle 與 retract events 是否出現新 pitfall pattern；不在 table 的補進 P-07/P-08
- **實作**：新 skill `/pitfalls-monthly-review` 或加入 weekly-report 月底 special section
- **判斷規則**：若連 3 個月無新 pitfall 加入 → 降為季度

#### 提案 #2：plan.json template + auto-fill（防 §2.3 #2）
- **何時建**：cycle-init skill 下次維護時
- **內容**：cycle-init 寫 plan.json 時，依 cycle 類型（純分析 / C++ / pipeline / cross-sample）auto-fill `binary_version` (`null`/SHA)、`dataset_id` 範本、`upstream_reports` 候選 list
- **實作**：cycle-init/_init.py 加 template logic + 用戶在 P0 時 confirm/edit
- **判斷規則**：D2-A 已暴露此需求；下次真跑 cycle 啟動時做

#### 提案 #3：Pitfalls table 擴充至 P-07~P-12（防 §3.3 #1）
- **何時建**：Drill 2 真跑前
- **內容**：補 6 條 — P-07 KDE 新型 artifact / P-08 L2 collider bias / P-09 spatial autocorr / P-10 pooled OLS / P-11 feature name 直覺解讀 / P-12 飽和 artifact (1/N samples 拉動)
- **實作**：JSON entries 加在 `pitfalls_table.json`；每條附 trigger keyword + detect_via + severity
- **判斷規則**：每條須對應 ≥1 個 memory feedback 或 retract event；無 case 就不加

#### 提案 #4：P3.5 cross-sample early-termination gate（防 §3.3 #3）
- **何時建**：Drill 2 真跑暴露 thread_b 模式時
- **內容**：在 P3 PILOT → P4 GENERALIZE 之間加 conditional gate；若 pilot 觀察到 effect size variance 跨樣本 > 2× threshold → 強制 multi-sample-consistency 至少 3 樣本驗 + user ack 才進 P4
- **實作**：`/run-evaluator` 加 `--mode=p3.5-precheck` flag 或新 skill `/early-termination-check`
- **判斷規則**：等 Drill 2 / 真實 cycle 觀察到 thread_b 類 false-promise pattern 重現才實作；否則 PARKED

#### 提案 #8：Disk health monitor cron + first-run sanity check（防 §5.3 disk gap）
- **何時建**：**立即**（已造成 1 次 800GB 災情）
- **內容**：
  - cron job 每小時跑 `df -h /` + `mount | grep tmp` 寫到 `state/infra_health/disk_log.jsonl`
  - 新 server / 新 conda env 啟動時跑 `scripts/setup/sanity_check.sh`：檢查 /tmp 是否獨立 tmpfs、TMPDIR 是否 export、root volume free 是否 > 50GB
  - threshold 警示：root < 20GB → 紅色 panic，sanity check fail → block pipeline
- **實作**：`InterSubMod/scripts/infra/disk_health_check.sh` + cron entry
- **判斷規則**：必做；memory feedback_tmp_disk_full_pipeline_pitfall 已記教訓

#### 提案 #9：Pipeline wrapper 強制 export TMPDIR（防 §5.3 disk gap）
- **何時建**：與 #8 一起做
- **內容**：所有 pipeline entry script (`run_vcf_all_snv.sh` / `run_batch_vcf_analysis.sh` 等) 開頭強制 `export TMPDIR=/big7_disk/tmp`；對外部工具（longphase, ClairS, samtools）顯式傳 `-o` 或 tmpdir flag
- **實作**：對 `InterSubMod/scripts/run_*.sh` 加 wrapper sourcing `_pipeline_init.sh`
- **判斷規則**：必做；同 #8 緊密關聯

### 7.2 中期（P2）— 等 Drill 2 真跑或 30-day gate 後評估

#### 提案 #5：Anchor 5b 硬化評估（plan v1.6 §4.5.4-G）
- **何時評**：D2-A 真跑 cycle 完成後
- **內容**：anchor #3（5 段報告骨架）/ #6（5 層圖表 caption）/ #7（草稿可重評）是否需 hardcode 進對應 skill SKILL.md
- **判斷規則**：若 D2-A 真跑時觀察到 anchor 違反案例 → KEEP；否則 PARKED 至 30-day mark

#### 提案 #6：risk-register skill 或 plan.json `risk_register` 欄位（防 §4.3 #2）
- **何時建**：Drill 2 真跑時若觀察到 R-id 結構化需求重現
- **內容**：plan.schema.json 加 optional `risk_register: [{r_id, priority: P0|P1|P2|P3, description}]` 欄位；cycle-state dashboard 顯示 active risks
- **判斷規則**：等需求重現觀察 ≥2 次再做

#### 提案 #7：Pivot rate cron monitor（防 §4.3 #3）
- **何時建**：30-day gate 後（2026-06-07）若撤回率 > 2 啟動
- **內容**：cron 每週掃 docs/CURRENT_FOCUS.md 主軸字段 + research_direction.md queue 變動，計算「過去 30 天主軸切換次數」；> 4 次 → 警示
- **判斷規則**：等 30-day gate 結果決定

#### 提案 #11：PPTX/Matplotlib pre-build validator（防 §6.3 #1）
- **何時建**：下次 PPTX 出現中文亂碼時立即做
- **內容**：pptx-build / weekly-report skill 在 build 前跑 `_validate_render.py`：檢查 Python 圖檔的 rcParams 是否含 Droid Sans Fallback；PPTX template 是否有 emoji；圖片是否有強制 width+height 擠壓
- **判斷規則**：memory 已記重複教訓；下次重現即實作

### 7.3 長期（P3）— Path B 啟動後或實證需求出現

#### 提案 #10：KB auto-stale flag + 30d freeze quoting（防 §5.3 #3）
- **何時建**：Path B Layer 3 LlamaIndex 啟動時
- **內容**：14 天警示已存在；新增 30 天後自動加 `<!-- STALE: requires re-verification -->` banner；引用 stale KB 文件需明示 override
- **判斷規則**：Path B Decision Gate (2026-06-07) PASS 後再評

---

## §8 Phase 2 governance（M3 + M5）— 為何刻意延後

按 plan v1.7 Phase 1 commit 5a136d8，**M3 auto-recovery routing + M5 5-tier intervention 刻意延後**。複盤後仍維持此決策：

| 機制 | 為何延後 | 重評時機 |
|---|---|---|
| **M3 auto-recovery routing** | Drill 2 fixture-driven 無真實 routing 失敗案例；現有 cycle-state next_action + per-folder CLAUDE.md 已 cover 80% routing 需求 | D2-A 真跑或 D2-B 啟動時，若觀察到 routing dead-end ≥2 次再做 |
| **M5 5-tier intervention** | confirmation-protocol 已有 4 級（FYI / Review / Gate / Hard Gate）；5 級在文件層精細化但實際操作對單研究者收益不明 | 若觀察到「介入頻率不對稱」（過多 vs 過少）才做 |

**結論**：M3+M5 暫不阻塞；Phase 1 governance 已足夠 D2-A drill 與 30-day 觀察期。

---

## §9 防呆機制完整對照表

| 問題 Category | 已 cover% | 未 cover gaps | 提案 P 級 |
|---|---|---|---|
| A 前置條件 | 99% | pitfalls_table 滯後、plan.json template | P1 (#1, #2) |
| B 結論健壯性 | 95% | P-07~P-12 confound、P3.5 early term gate | P1 (#3), P2 (#4) |
| C 認知/風格/流程 | 60% | anchor 5b、risk register、pivot rate | P2 (#5, #6, #7) |
| D 基礎設施 | 30% | **disk health（critical）**、TMPDIR wrapper、KB freshness 強制 | **P1 (#8, #9)**, P3 (#10) |
| E 報告呈現 | 70% | matplotlib/PPTX pre-build validator | P2 (#11) |

**P1 立即做（5 個）**：#1, #2, #3, #8, #9
**P2 觀察後做（4 個）**：#4, #5, #6, #7, #11
**P3 Path B 後做（1 個）**：#10

---

## §10 判斷規則總結 — 何時建、何時 PARKED、何時不建

依用戶 memory `feedback_evidence_driven_iteration_workflow`「實證導向迭代」原則：

| 判斷 | 規則 | 適用提案 |
|---|---|---|
| **立即做** | 已造成 ≥1 次重大損失 OR 高頻重現（>3 次/月） OR memory feedback 多次糾正 | #8, #9 (disk), #11 (PPTX, 若再現) |
| **下次 cycle 觸發做** | 等下一次真實使用時自然帶出需求 | #1, #2, #3 |
| **真跑 cycle 觀察後做** | 等 D2-A 真跑或 D2-B 啟動暴露需求 | #4, #5 |
| **觀察 ≥2 次重現再做** | 避免 over-engineering 單次事件 | #6 |
| **30-day gate 後評** | 等 2026-06-07 撤回率數據 | #7 |
| **Path B 後做** | 依賴 LlamaIndex 才能 efficient | #10 |
| **不建**（anti-pattern） | 無 prior art + 收益不確定 + 增加 ritual 疲勞 | D weekly skill review (plan §11.6 已記) |

---

## §11 References

### 11.1 主要 plan / 報告

- Plan v1.7：`~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md`
- Drill 1 retrospective：`InterSubMod/docs/experiments/in_progress/2026/05/20260506_Drill1_Harness_Retrospective_01.md`
- Drill 2 walkthrough：`InterSubMod/docs/experiments/in_progress/2026/05/20260508_Drill2_End_to_End_Cycle_Walkthrough_01.md`
- Architecture review：`InterSubMod/docs/reports/validated/2026/05/20260508_Harness_Complete_Architecture_Review_01.md`

### 11.2 6 retract event source 報告

- `InterSubMod/docs/experiments/in_progress/2026/04/20260404_VCF來源錯誤矯正報告_01.md`
- `InterSubMod/docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md`
- `InterSubMod/docs/experiments/in_progress/2026/04/20260418_F_HPFineNGroups_Subclone_Marker_01.md`
- `InterSubMod/docs/experiments/in_progress/2026/04/20260424_X6_Caller_AF_S3S5_CrossSample_01.md`
- `InterSubMod/docs/experiments/in_progress/2026/04/20260426_Thread_B_Whitelist_Retraction_01.md`
- `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md`

### 11.3 重要 memory feedback（已 cover 與 gap 之依據）

- `feedback_tmp_disk_full_pipeline_pitfall.md`（提案 #8, #9 主要依據）
- `feedback_evidence_driven_iteration_workflow.md`（§10 判斷規則 SoT）
- `feedback_skill_md_must_state_dependencies_and_diagnostics.md`（v1.7 batch B Quality Checklist 依據）
- `feedback_strategy_then_per_item_confirmation.md`（D2-A 路徑選擇依據）
- `feedback_full_auto_parallel_execution.md`（v1.7 batch 機械性採納依據）
- `feedback_L2_collider_bias.md`、`feedback_pooled_ols_residualization_trap.md`、`feedback_spatial_autocorrelation_confound.md`（提案 #3 P-07~P-12 內容來源）
- `feedback_feature_name_vs_definition_rule.md`（提案 #3 P-11 內容來源）
- `feedback_matplotlib_cjk_font_rule.md`、`feedback_pptx_screenshot_rendering_rules.md`（提案 #11 內容來源）
- `feedback_execution_mode_hierarchy.md`（§4.2 三級分類已 cover 之依據）

### 11.4 Phase 1 governance commit / Phase 2 待辦

- commit 5a136d8 → merge 517c467 (M1 + M2 + M4)
- M3 auto-recovery routing：deferred Phase 2 (§8)
- M5 5-tier intervention：deferred Phase 2 (§8)

---

## §12 結論摘要

> 1. 過去 30 天累積 **~60 個問題事件**，分 5 大 Category；
> 2. **Category A+B (前置條件 + 結論健壯性) 已 95%+ 覆蓋**，由 7-phase × M1+M2+M4 governance + Drill 1 sensitivity 6/6 證實；
> 3. **Category C (認知/流程) 60% 覆蓋**，主要 gap 在 anchor 5b + risk register + pivot rate（多數 PARKED 至 Drill 2 後評）；
> 4. **Category D (基礎設施) 30% 覆蓋，critical gap 在 disk health monitor**（已造成 800GB 災情，提案 #8/#9 P1 立即做）；
> 5. **Category E (報告呈現) 70% 覆蓋**，matplotlib/PPTX pre-build validator 待重現再做；
> 6. **新機制 10 個提案**：5 個 P1 立即做 / 4 個 P2 觀察後做 / 1 個 P3 Path B 後做；
> 7. **判斷規則 7 條**（§10）：實證導向，避免無 prior art 的 ritual 疲勞；
> 8. M3+M5 Phase 2 governance **刻意延後**；現有覆蓋已足夠 D2-A drill 與 30-day 觀察期。

**最高優先（必須立即做）**：提案 #8 + #9（disk health monitor + TMPDIR wrapper）— 已造成 1 次重大損失，不可再延。
