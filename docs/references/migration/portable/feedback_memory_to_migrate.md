# 風格 Feedback Memory 遷移清單（41 個）

> 用戶 home memory 目錄: `/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory/`
> 新機目錄: `~/.claude/projects/<encoded-new-project-path>/memory/`
> （路徑編碼規則：`/` → `-`，前綴加 `-`）

---

## 拷貝指令（一次性）

```bash
OLD_MEM=/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory
NEW_MEM=$(ls -d ~/.claude/projects/*<new-project-name>*/memory/ 2>/dev/null | head -1)
[ -z "$NEW_MEM" ] && { echo "ERROR: 新機 memory 目錄不存在；先用 claude 開新 session 觸發目錄建立"; exit 1; }

# 拷貝 41 個 feedback memory
for f in \
  feedback_observation_scope_default_comprehensive.md \
  feedback_task_first_then_doc_then_plan.md \
  feedback_existing_artifacts_must_verify.md \
  feedback_pi_report_html_stack.md \
  feedback_design_principles_canonical.md \
  feedback_outside_claim_must_query_kb.md \
  feedback_tmp_disk_full_pipeline_pitfall.md \
  feedback_md_path_prefix_rule.md \
  feedback_full_auto_parallel_execution.md \
  feedback_execution_mode_hierarchy.md \
  feedback_strategy_then_per_item_confirmation.md \
  feedback_skill_md_must_state_dependencies_and_diagnostics.md \
  feedback_evidence_driven_iteration_workflow.md \
  feedback_merged_dataset_af_and_loh_pitfalls.md \
  feedback_feature_name_vs_definition_rule.md \
  feedback_paper_positioning_de_prioritized.md \
  feedback_risk_structured_iterative.md \
  feedback_small_scale_validation_first.md \
  feedback_spatial_autocorrelation_confound.md \
  feedback_research_direction_corrections_20260405.md \
  feedback_figure_layout_standard.md \
  feedback_L2_collider_bias.md \
  feedback_pooled_ols_residualization_trap.md \
  feedback_visual_inspection_requirement.md \
  feedback_project_subfolder_structure.md \
  feedback_researcher_claim_needs_empirical_verification.md \
  feedback_cynefin_domain_classification.md \
  feedback_productive_failure_reopen_threshold.md \
  feedback_matplotlib_cjk_font_rule.md \
  feedback_pptx_screenshot_rendering_rules.md \
  feedback_pptx_html_preview_font_sizing.md \
  feedback_self_phasing_PI_no_build_html_py.md \
  feedback_pptx_generation_workflow.md \
  feedback_pptx_visual_first_philosophy.md \
  feedback_pptx_bilingual_formatting.md \
  feedback_pptx_director_storyboard.md \
  feedback_ppt_slide_design_rules.md \
  feedback_pptx_text_formatting_rules.md \
  feedback_pptx_term_explanation_rule.md \
  feedback_weekly_report_workflow.md \
  feedback_pptx_presync_workflow.md \
  feedback_ppt_minimal_visual_first.md \
; do
  if [ -f "$OLD_MEM/$f" ]; then
    cp "$OLD_MEM/$f" "$NEW_MEM/"
    echo "  cp $f"
  else
    echo "  MISSING: $f"
  fi
done

echo "完成；接下來建立 MEMORY.md 索引（用 portable/MEMORY_template.md 起手）"
```

---

## 分群索引（用 MEMORY.md 引用）

### 暫停判斷 / 確認協議（5）
- `feedback_execution_mode_hierarchy` — 第一次新概念彙報後等確認 / 重複執行一行告知 / 全自動 fan-out
- `feedback_strategy_then_per_item_confirmation` — strategy 同意後逐項實作確認
- `feedback_full_auto_parallel_execution` — 盤點後執行階段，並行 subagent 吃滿機器
- `feedback_evidence_driven_iteration_workflow` — batch 中途暫停由驗證指導後續
- `feedback_task_first_then_doc_then_plan` — task 排好 → 等文件 → 依文件規劃 3-step

### Task type / Scope（3）
- `feedback_observation_scope_default_comprehensive` — 觀察 scope default 全範圍 + 見樹也見林 + 6 類 task type
- `feedback_risk_structured_iterative` — R-id 固定編號 + P0/P1/P2/P3 排序 + 批次彙報
- `feedback_small_scale_validation_first` — 新方向先 <2h pilot 確認潛力

### 方法論嚴謹（8）
- `feedback_cynefin_domain_classification` — Cynefin 域分類為暫停判斷前置層
- `feedback_productive_failure_reopen_threshold` — NEGATIVE / NO-GO reopen C1/C2/C3 3 條件
- `feedback_existing_artifacts_must_verify` — 既有 artifact 必先驗證才能引用
- `feedback_outside_claim_must_query_kb` — outside-claim 必先查 KB
- `feedback_researcher_claim_needs_empirical_verification` — researcher claim 必須實測驗證才能定論
- `feedback_L2_collider_bias` — residualize on AF 產生虛假信號；必須 L3 AF-bin 交叉驗證
- `feedback_pooled_ols_residualization_trap` — pooled OLS 殘差保留分組信息；必須 within-group OLS
- `feedback_spatial_autocorrelation_confound` — chr+pos 聚合特徵 AUC 必須 mid-TP-rate window 驗證

### PPT / 簡報設計（11）
- `feedback_design_principles_canonical` — Tufte + CRAP + Assertion-Evidence + Zen + Nature + WCAG 12 條 canonical
- `feedback_ppt_minimal_visual_first` — slide 重點圖 + 標題 + ≤50 字
- `feedback_pi_report_html_stack` — 1 index + N 子 HTML hybrid；Plotly + DataTables + Alpine + MVP.css
- `feedback_pptx_screenshot_rendering_rules` — Latin+CJK 混合 font fallback + fit-within 等比
- `feedback_pptx_html_preview_font_sizing` — `.slide-canvas` PPT 16:9 字級 Body ≥16px
- `feedback_pptx_generation_workflow` — aspect ratio + 自適應佈局 + 截圖驗證
- `feedback_pptx_visual_first_philosophy` — 物件圖示為主、文字為輔
- `feedback_pptx_bilingual_formatting` — EN 縮排 0.25"+60% 字號
- `feedback_pptx_director_storyboard` — Storyboard 審查；先信任→分析→根因
- `feedback_ppt_slide_design_rules` — 標題≤15字；每頁2-3焦點；結論句非數據
- `feedback_pptx_text_formatting_rules` — 禁止斜體；適合標題可粗體
- `feedback_pptx_term_explanation_rule` — 首次出現必配圖；圖示 60-70%
- `feedback_pptx_presync_workflow` — 修改 build_pptx.py 前必先 diff 偵測手動編輯
- `feedback_self_phasing_PI_no_build_html_py` — 直接 Edit HTML，禁止 build_html.py（事故教訓）

### 報告 / 寫作工作流（5）
- `feedback_weekly_report_workflow` — 受眾=PI；7 Phase 流程；AI 整理草稿用戶確認
- `feedback_md_path_prefix_rule` — .md 路徑必須以 `<project>/` 開頭
- `feedback_paper_positioning_de_prioritized` — 分析完整性 > 投稿速度
- `feedback_skill_md_must_state_dependencies_and_diagnostics` — SKILL.md 必含 Phase&Chain / Dependencies / Failure Mode
- `feedback_project_subfolder_structure` — 多步驟研究建專案子資料夾，含 00_ 索引 + StepN_ 前綴

### 視覺化（3）
- `feedback_matplotlib_cjk_font_rule` — Python 生圖 rcParams 必加 Droid Sans Fallback
- `feedback_visual_inspection_requirement` — per-region 視覺化 + IGV 截圖供人眼驗證
- `feedback_figure_layout_standard` — 分 Paired/TO；細緻觀察 2x2 × N samples 固定順序

### 基礎設施（1）
- `feedback_tmp_disk_full_pipeline_pitfall` — pipeline 中間檔不可寫 root /tmp；TMPDIR 設大碟

### 其他特定但通用化後仍有價值（3）
- `feedback_research_direction_corrections_20260405` — 研究方向三項核心修正範例
- `feedback_merged_dataset_af_and_loh_pitfalls` — 合成檔欄位陷阱範例
- `feedback_feature_name_vs_definition_rule` — Feature name 直覺解讀陷阱

---

## 清洗建議

部分 feedback memory 引用 InterSubMod 特定事件（如「5/24 HKU handoff」「Cycle 4 LR 失敗」「HCC1395 +0.057 fabrication」）。

**新機新任務後，可選兩條路**：
1. **保留歷史 context**：feedback 內具體事件不刪，作為「lesson came from」歷史紀錄
2. **重述為通用 lesson**：用 `/memory-consolidation` skill 跑一次清洗，把特定事件改為「past incident」抽象表述

**判斷標準**：若 feedback 規則本身與事件無關（如 design_principles_canonical），保留原樣即可；若高度依賴具體事件理解（如 self_phasing_PI_no_build_html_py），用戶可酌情改寫。

---

## 不遷移的 memory（不在此清單）

- `project_*` 系列 38 個（active research / pending / concluded）— 純 ISM 研究 artifact
- `reference_*` 系列 19 個（ISM source code 行號等）— 純 ISM 引用

如新機要做不同任務（如 web app 開發），這些都不適用。
