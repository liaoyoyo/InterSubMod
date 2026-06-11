# Memory Index

## Skills（新機 baseline 17 個風格 skills）

- 已在 `.claude/skills/` 拷貝；用 `/help` 查可用列表

## Active Tasks

- _TBD：新專案任務記錄_

## Pending Items

- _TBD：新專案 pending 記錄_

## Methodology & Preferences（41 個風格 feedback memory）

### 暫停判斷 / 確認協議
- [執行模式分級：第一次處理需確認 / 已驗證可推進 / 全自動為終態](feedback_execution_mode_hierarchy.md) — 🟠 新概念彙報後等確認 / 🟡 重複執行一行告知 / 🟢 全自動 fan-out
- [Strategy 同意後逐項實作確認](feedback_strategy_then_per_item_confirmation.md) — plan 決 strategy / 實作每變動 diff preview + 等 ack
- [全自動 + multi-agent 並行執行偏好](feedback_full_auto_parallel_execution.md) — 盤點後執行階段並行 subagent + 吃滿機器
- [實證導向迭代：batch 中途暫停由驗證指導後續](feedback_evidence_driven_iteration_workflow.md) — 接受 batch 不完整風險，跳 Drill 驗證後再決定
- [緊急 handoff 任務 3-step workflow](feedback_task_first_then_doc_then_plan.md) — Step 1 task type + anchor / Step 2 等文件 / Step 3 規劃

### Task type / Scope
- [觀察 scope default = 全範圍 + 見樹也見林 + 6 類 task type](feedback_observation_scope_default_comprehensive.md) — 啟動 prompt keyword 主動 recall
- [風險結構化 + 批次彙報](feedback_risk_structured_iterative.md) — R-id 固定編號 + P0/P1/P2/P3 排序
- [小規模快速驗證優先](feedback_small_scale_validation_first.md) — 新方向先 <2h pilot 確認潛力

### 方法論嚴謹
- [Cynefin 域分類為暫停判斷前置層](feedback_cynefin_domain_classification.md) — 新方向前先 5 域分類因果性質
- [Productive Failure reopen 3 條件 C1/C2/C3](feedback_productive_failure_reopen_threshold.md) — NEGATIVE/NO-GO 結論 reopen 門檻
- [既有 artifact 必先驗證才能引用](feedback_existing_artifacts_must_verify.md) — 既有 HTML/MD/TSV 不可直接信任
- [outside-claim 必先查 KB](feedback_outside_claim_must_query_kb.md) — 外部 tool/論文 claim 必查 KB
- [researcher claim 必須實測驗證才能定論](feedback_researcher_claim_needs_empirical_verification.md) — researcher 推測屬 L3，移除前必須升 L1
- [L2 Collider Bias 警告](feedback_L2_collider_bias.md) — residualize on AF 產生虛假信號
- [Pooled OLS Residualization 陷阱](feedback_pooled_ols_residualization_trap.md) — 必須 within-group OLS
- [空間 auto-correlation confound 警告](feedback_spatial_autocorrelation_confound.md) — chr+pos 聚合特徵必須 mid-TP-rate window 驗證

### PPT / 簡報設計
- [簡報與報告設計準則 12 條 canonical](feedback_design_principles_canonical.md) — Tufte + CRAP + Assertion-Evidence + Zen + Nature + WCAG
- [PPT 極簡 visual-first 偏好](feedback_ppt_minimal_visual_first.md) — slide 重點圖 + 標題 + ≤50 字
- [PI 報告 Hybrid HTML 結構與 JS 套件選型](feedback_pi_report_html_stack.md) — 1 index + N 子 HTML / Plotly + DataTables + Alpine + MVP.css
- [PPTX 截圖字型 + 圖片比例強制規則](feedback_pptx_screenshot_rendering_rules.md) — Latin+CJK fallback + fit-within 等比
- [HTML preview PPT 字級規範](feedback_pptx_html_preview_font_sizing.md) — `.slide-canvas` 16:9 Body ≥16px
- [PPTX 生成完整流程規範](feedback_pptx_generation_workflow.md) — aspect ratio + 自適應 + 截圖驗證
- [PPTX 視覺優先設計哲學](feedback_pptx_visual_first_philosophy.md) — 物件圖示為主、文字為輔
- [PPTX 雙語格式偏好](feedback_pptx_bilingual_formatting.md) — EN 縮排 0.25"+60% 字號
- [PPTX 導演審查要求](feedback_pptx_director_storyboard.md) — Storyboard 先信任→分析→根因
- [PPT Slide 設計規則](feedback_ppt_slide_design_rules.md) — 標題≤15字；每頁2-3焦點
- [PPTX 文字格式規則](feedback_pptx_text_formatting_rules.md) — 禁止斜體；標題可粗體
- [PPTX 術語解釋與首次出現規則](feedback_pptx_term_explanation_rule.md) — 首次出現必配圖；圖示 60-70%
- [PPTX 修改前 diff 同步](feedback_pptx_presync_workflow.md) — build script 前必先 diff 偵測手動編輯
- [self_phasing_synthesis_PI 已廢棄 build_html.py](feedback_self_phasing_PI_no_build_html_py.md) — 直接 Edit HTML，禁止 build script 覆蓋手動 commit

### 報告 / 寫作工作流
- [週報工作流程偏好](feedback_weekly_report_workflow.md) — 受眾=PI；7 Phase；AI 整理 + 用戶確認
- [.md 路徑必須以 project 名稱開頭](feedback_md_path_prefix_rule.md) — UserPromptSubmit hook 提醒
- [論文定位降優先](feedback_paper_positioning_de_prioritized.md) — 分析完整性 > 投稿速度
- [SKILL.md 必說明依賴鏈與診斷](feedback_skill_md_must_state_dependencies_and_diagnostics.md) — Phase&Chain / Dependencies / Failure Mode 3 段
- [研究專案子目錄結構](feedback_project_subfolder_structure.md) — 00_ 索引 + StepN_ 前綴

### 視覺化
- [Matplotlib 圖片 CJK 字型必注入規則](feedback_matplotlib_cjk_font_rule.md) — rcParams 必加 Droid Sans Fallback
- [肉眼檢視文件需求](feedback_visual_inspection_requirement.md) — per-region 視覺化 + IGV 截圖
- [圖表排版規範](feedback_figure_layout_standard.md) — 分 Paired/TO；2x2 × N samples 固定順序

### 基礎設施
- [/tmp 800 GB 災情 — pipeline 中間檔不可寫 root /tmp](feedback_tmp_disk_full_pipeline_pitfall.md) — TMPDIR 設大碟

### 通用化後仍有價值的特定案例
- [Feature name 直覺解讀陷阱 — 必查原始碼](feedback_feature_name_vs_definition_rule.md)
- [Merged 合成檔 AF/LOH 資料陷阱](feedback_merged_dataset_af_and_loh_pitfalls.md)
- [研究方向三項核心修正](feedback_research_direction_corrections_20260405.md)

---

## Concluded (prevent re-investigation)

- _TBD：新專案 concluded 結論記錄_

---

> **MEMORY.md 維護建議**：
> 1. 保持 <200 行（context 自動載入限制）
> 2. 每條 entry 一行 ~150 字以下
> 3. 完成的 active task → 搬到 Concluded
> 4. 新 feedback memory 必同步寫到對應分群
> 5. 每月跑一次 `/memory-consolidation` 整理過時 entry
