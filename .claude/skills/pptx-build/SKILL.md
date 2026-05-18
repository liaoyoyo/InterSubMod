---
name: pptx-build
description: InterSubMod PPT 製作 sub-skill — 從母稿（weekly-report 產出）或從零開始產 PPTX。3 層確認 outline → section → slide + 6 報告模板識別 + §20 主軸聚焦 6 階段過濾 + Tier 1/2/3 三層分流（slide / speaker note / oral-optional）+ Claude Vision 10-checkpoint review。支援 --from-draft 介面從 weekly-report master_draft.md 自動讀取 main thesis / report_type / audience，跳過 P1 部分項目。觸發：「簡報」「PPT」「PPTX」「deck」「教授報告」「教授級簡報」「投影片」「presentation」「製作簡報」。對外 PPT 結論／升 tier ⭐4-5 / PI 簡報前必繼承 InterSubMod/.claude/skills/scientific-rigor/SKILL.md §2-§7（slide title 禁 confidence 詞彙、數字 slide 必含 Cohen/CI ribbon、因果 claim 配對 DAG 圖）。SKIP WHEN 週報母稿尚未產出（先用 weekly-report W1-W7）、僅是場景識別／路徑分流（用 myPPT 入口）、PDF/markdown report 即可不需 PPTX、純 figure 重做（用 image-gen + image-vision-check）、PPTX 已存在僅微調文字（直接 Edit 即可）。
allowed-tools: Read, Write, Edit, Bash, Glob, Grep, AskUserQuestion, Agent
user-invocable: true
---

# pptx-build Skill

PPT 製作子流程，承擔原 myPPT skill 的 P1-P5（outline → section → slide → visual review → speaker script）內容。**不負責 raw data 收集 / 母稿撰寫**（這是 weekly-report skill 的責任）。

> **2026-05-05 從 myPPT 拆出**：原 myPPT 6-stage pipeline 拆分為 (a) weekly-report W1-W7 母稿生成 + (b) pptx-build P1-P5 PPTX 製作 + (c) myPPT 輕量總入口（場景識別 → 委派）。詳見 `InterSubMod/.claude/skills/weekly-report/SKILL.md`（前身 weekly-master-draft 已 merge）。

---

## 執行模式

- **互動模式**（預設）：每個 P1-P5 checkpoint 暫停確認
- **全自動**：保 P1（main thesis 鎖定）+ P5（speaker script）必停；P2/P3/P4 30 秒倒數 fast-track

## 觸發路徑

### 路徑 A：從 weekly-report 母稿（handoff A 觸發）
```
weekly-report C4 母稿確認 → handoff A 選 → /pptx-build --from-draft <path>
```
自動讀 frontmatter（report_type / main_statement / audience / suggested_pptx_template），跳過 P1 main thesis 鎖定 + 報告類型識別。

### 路徑 B：從零開始（直接觸發）
用戶說「我要做簡報/教授報告」→ 走完整 P1-P5。

---

## 5 階段流程（P1-P5）+ C1-C5 Checkpoint

```
[U: 我要做簡報] 或 [from-draft]
    ↓
P1. Audience & Goal             → C1 main thesis ≤ 30 字 + audience tier ★ 必停
    ↓
P2. Outline                     → C2 5-7 段 + 每段 thesis sentence
    ↓
P3. Section batch (×N)          → C3 每 section 5-7 slide + focal point ≤ 20 字
    ↓
P4. Slide build (×24)           → C4 三層分流 + Vision 10-check 逐張
    ↓
P5. Speaker script              → C5 字數 ↔ 時長 + Tier 3 拆遷 ★ 必停
    ↓
[Output: PPTX + speaker script + audit]
```

---

## 6 報告模板識別（§11 標準）

| 觸發 keyword | 模板 | narrative skeleton |
|------------|------|--------------------|
| 修補/修復/優化/fix | `templates/improvement_report.md` | Before → Problem → Cause → Solution → Verify → Impact |
| A vs B / 對照 / benchmark | `templates/comparison_report.md` | Setup → A → B → Side-by-side → Verdict |
| 月會/PI/KPI/quarterly | `templates/executive_summary.md` | TL;DR → Top 3 → Risks → Asks |
| 結果/實驗/統計 | `templates/data_showcase.md` | Hypothesis → Data → Stats → Caveats |
| 教學/解釋/方法論 | `templates/concept_walkthrough.md` | Why → Define → Mechanism → Examples → Boundary |
| 教授/同儕審查/thesis | `templates/academic_defense.md` | Background → Question → Method → Result → Discussion → Future |

混合用例：取主敘事用模板（如 v3 = improvement 主 + academic 深度）。

---

## 強制 Checkpoint（C1-C5）

- **C1（P1 結束）**：Main Thesis ≤ 30 字 + Audience tier ★ fast-track 必停
- **C2（P2 結束）**：5-7 段 outline + 每段 thesis sentence
- **C3（每 section ×N）**：5-7 slide 標題 + focal point ≤ 20 字
- **C4（每 slide ×24）**：三層分流（Tier 1 slide / Tier 2 note / Tier 3 oral）+ Vision 10-check
- **C5（P5 結束）**：字數 vs 預設時長 + Tier 3 拆遷建議 ★ 必停

詳細 prompt 模板 → `prompts/{outline,section,slide,focal_point_audit,tier_evaluation,visual_review,from_draft_loader}_confirm.md`

---

## §20 主軸聚焦 6 階段過濾

詳見 `playbook.md` §20。每張 slide 設計前強制 6 階段：

- **A. Main Thesis 鎖定**（≤ 30 字，整份簡報唯一）
- **B. Slide Focal Point 鎖定**（≤ 20 字，每張唯一）
- **C. Tier S/A/B/C/D 評分**（每筆素材分流前強制評分後處置）
- **D. Definition / Prerequisite / Body / Conclusion 比例**（Body ≥ 60%、Def+Prereq ≤ 25%）
- **E. 6 問 audit**（slide 進 build 前）
- **F. 雜訊紅旗清單**

---

## Tier 三層分流規則

| Tier | 位置 | 範圍 | 字數 |
|:-:|------|------|------|
| 1 | on slide | 直接構成 focal point | ≤ 6 elements，中文 ≤ 60 字、英文 ≤ 30 word |
| 2 | speaker note | 必講細節（含 caveat / source）| 75-90 sec/slide ≈ 中 300-360 字 |
| 3 | oral-optional | 依現場時間決定 | speaker note 末尾 `[ORAL-OPTIONAL]` 標 |

---

## Visual Review 10-Checkpoint

每張 slide build 後跑 Claude Vision 10-check（`prompts/visual_review.md`）：
1. Heading = thesis sentence？
2. 1 idea/slide？
3. ≤ 6 elements？
4. 文字密度（中 ≤ 60 / 英 ≤ 30）？
5. Distracted takeaway 可抓主旨？
6. 視覺對比 / 字型 / 字級？
7. 圖表 vs 文字 ≥ 60% 視覺？
8. 引文 / 數據來源？
9. 1 minute timing check？
10. 技術災難備援？

---

## 紅旗清單（看到立即削減）

**Cognitive overload**：
- 標題是 generic label（「結果」「分析」「討論」）→ 改 assertion sentence
- 一張 slide > 6 element → 拆兩張
- 中文 > 60 字 / 英文 > 30 word → 削減
- 一張 slide 多 idea → 拆兩張

**雜訊用語**（speaker note）：
- 「順便提一下」/「附帶說明」/「另外」（非主線）→ 削減
- speaker note > 360 字 → 拆 Tier 3

**比例**：
- Definition + Prerequisite > 25% → 合併
- Body slide < 60% → 重結構
- Conclusion < 15% → 補

---

## --from-draft 介面

當用 `/pptx-build --from-draft <path>` 觸發：

1. 讀 master_draft.md frontmatter：
   - `main_statement` → 鎖定 main thesis（跳過 P1 main thesis 鎖定）
   - `report_type` → 載入對應 template（improvement/comparison/...）（跳過 P1 報告類型識別）
   - `audience` / `target_duration_min` → 計算 slide 數
   - `suggested_pptx_template` → 預設模板
2. 從母稿萃取：
   - Layer 0.1 + §1 + §2 → cover slide
   - Layer 2 Thread × N → 1-2 slide per Thread
   - Layer 3 整合 → 1 integration slide
   - Layer 4 §16 + §17 → future + Q&A backup slides
3. 仍走 P2-P5（outline / section / slide / speaker checkpoint）

詳見 `prompts/from_draft_loader.md`。

---

## Style Library 調用

`style_library/colors/palette.yaml`（6 色 token）+ `style_library/objects/*.yaml`（15 物件）+ `style_library/layouts/*.yaml`（12 版型）+ `style_library/examples/case_studies.md`。

每個 object 含 `tier_recommendation`；每個 layout 含 `focal_point_zone`。

詳見 `playbook.md` §13/§14/§15/§16/§17。

---

## ppt_toolkit 引用

Reusable Python helpers：`InterSubMod/tools/ppt_toolkit/`（不在 skill 內，跨 skill 共用）

```python
from ppt_toolkit import (
    add_text_with_fallback,    # Latin/CJK font fallback
    fit_image_within,           # 等比 fit + fallback
    set_speaker_note,           # 字數上下界檢查
    assertion_title,            # 強制 thesis-style 標題
    tier_aware_speaker_note,    # 自動拆 Tier 2/3
    claude_vision_review,       # 10-checkpoint 表格
)
from ppt_toolkit.style_library_loader import load_object, load_layout
```

API 索引 → `tools/README.md`

---

## 嚴謹度繼承（/scientific-rigor）

PPTX / HTML slide 對外發佈權重最高，必嚴格繼承 `InterSubMod/.claude/skills/scientific-rigor/SKILL.md`：

- **§2 證據分級**: 每張 claim slide 必標 ⭐⭐⭐⭐⭐ ribbon（slide 角落 metadata）
- **§3 Effect Size ribbon**: 任何「+X」「-Y」數字 slide 必含 Cohen / CI 註腳；slide title 禁止 confidence 詞彙
- **§4 DAG**: TP/FP/AUC 因果 claim slide 配對 DAG 子圖（必要時 1 slide 1 DAG）
- **§7 Pre-registration**: PI report PPTX 結尾段必含 pre-reg vs actual 對照表
- **§10.1 Feynman**: slide 結論句必通過「12 歲能懂」測試（無 jargon stack）

**最小可用子集**:
- PI report / 對外論文 PPTX: §2 + §3 + §4 + §7 全跑（強制 §0.5 最高層級）
- 內部 weekly review PPTX: §2 + §3
- 草稿 / 學習材料: §2

**Slide 紅旗**（看到立即削減）:
- 1 張 slide 多個未標 evidence tier 的 claim
- F1 / AUC 數字無 ribbon
- 對外「更好」「最佳」「顯著」無數據佐證

## 與其他 skill 關聯

- **上游觸發** → `weekly-report` C4 後 handoff A 自動觸發 `/pptx-build --from-draft <path>`
- **平行不重疊** → `myPPT`（總入口場景識別 → 委派 weekly-report 或 pptx-build）
- **規範來源** → `confirmation-protocol`（C1-C5 對應 Hard Gate / Gate / Review）/ `doc-standards`（PPTX 輸出命名）

## 輸出位置

```
{master_draft_dir}/pptx_build/
├── 00_storyboard.md           # P3.5 Storyboard 三幕
├── 01_full_narrative.md       # P2 outline
├── 02_slide_outline.md        # P3 section list
├── 03_slide_layout_script.md  # P4 + P5
├── build_pptx.py              # 生成腳本
├── output.pptx                # 最終 PPTX
└── wireframes/                # 截圖驗證
```

或 `InterSubMod/docs/presentations/validated/YYYY/MM/<title>/`（直接觸發無 master_draft 時）。

## 注意事項

1. **不直接產母稿**（那是 weekly-report 責任）
2. **Vision 10-check 不可跳過**（每張 slide 強制）
3. **每張 slide 進 build 前必填 focal point ≤ 20 字**
4. **C1 / C5 fast-track 必停**（main thesis / speaker script timing）

---

## 用戶修正分類 + 個人風格累積（v2.4 強制機制）★

任何用戶提出修正/改善建議時，**強制 AskUserQuestion 分類**：
- 通用必要 → 寫入 `style_library/personal_style_log.md`，後續所有 PPT 自動套用
- 本次特定 → 只本次修，不污染通用規則
- 不確定 → AI 用 5 維度評估後建議
- 兩者皆是 → 標 [PROVISIONAL]，≥3 次升級通用

詳見 `prompts/feedback_classification.md` + `style_library/personal_style_log.md`。

每張 slide build 前 **自動讀取 personal_style_log active 規則**，加進 §20.E / visual_review / multi_agent_review 各層檢核。
