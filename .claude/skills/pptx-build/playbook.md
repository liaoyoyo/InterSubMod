# pptx-build playbook (§1-§24 主方法論)

> 從原 myPPT 設計的 §1-§24 落地。SKILL.md 是觸發殼，本檔是細則。
> 含 PLOS Ten Simple Rules + Assertion-Evidence + Duarte 三段敘事 + Tier 三層分流 + §20 主軸聚焦 6 階段過濾。

## §1 Audience & Goal（前置必填）

| 欄位 | 說明 |
|------|------|
| 受眾 | 教授 / PI / 同儕 / 業界 |
| 知識前提 | 假設聽眾已知 / 不知道哪些 |
| 預設時長（min）| 通常 20-30 |
| 聽眾目標 transformation | What is → New bliss |
| 可接受 cognitive load | low/medium/high |
| **Main Thesis** | **≤ 30 字（整份簡報唯一論點）** |
| Audience tier | T1 (PI/教授) / T2 (同儕) / T3 (業界) |

**fast-track 必停 checkpoint**。從 weekly-report 母稿觸發時，main_statement 自動鎖入此欄。

## §2 Outline 階段（C2）

- Duarte 三段對齊（What is → What could be → New bliss）
- 大段數 5-7 段
- 每段 assertion sentence 一句
- 字數 ↔ 時長預估表（中 400 字/min；slide 數 ≈ time min）
- **必停**：用戶批 outline 才能進 §3

## §3 Section 階段（C3 ×N）

- 每 section 列 5-7 slide 標題（thesis sentence form）
- 每張 slide 核心 assertion + 主視覺類型（schematic / chart / IGV / table）
- **focal point ≤ 20 字**標註（每張 slide 唯一焦點）
- **必停**：用戶批一個 section 後才進下一個（不一次給全部）

## §3.5 Storyboard（從原 myPPT/Phase 3.5 搬入）

Three-Act 故事弧：定錨 → 信任 → 核心 → 根因+轉向 → 行動。
產 `00_storyboard.md`：每頁 storyboard（核心訊息 + 觀眾心理 + 佈局草稿 + 偏離檢查）。

## §4 Slide 階段（C4 逐張）

每張 slide 三層分流：
- **Tier 1 (slide)**：assertion 標題 + 主視覺 + ≤ 6 element 文字
- **Tier 2 (speaker note)**：必講細節（含 caveat）；75-90 sec/slide ≈ 中 300-360 字
- **Tier 3 (oral-optional)**：依現場時間決定；speaker note 末尾 `[ORAL-OPTIONAL]` 標

每張 slide 走完整 build → wireframe → Vision 10-check → 用戶確認 → 下一張。

## §5 Visual Review 10-Checkpoint（PLOS + AE + cognitive load）

1. Heading = thesis sentence vs generic label
2. 1 idea/slide
3. ≤ 6 elements
4. 文字密度（中 ≤ 60 / 英 ≤ 30 word per slide，不含標題）
5. Distracted takeaway（無口述能抓主旨）
6. 視覺對比（顏色 / 字型 / 字級）
7. 圖表 vs 文字 ≥ 60% 視覺
8. 引文 / 數據來源
9. 1 minute timing check（speaker note 字數推估）
10. 技術災難備援（PDF / 圖片 fallback / 無動畫依賴）

詳見 `prompts/visual_review.md`。

## §6 Speaker Script 階段（C5）

- 每張 slide Tier 2 speaker note 字數 → 預估 sec
- 加總 vs 預設時長
- 超時：標 `[ORAL-OPTIONAL]` 拆 Tier 3
- **fast-track 必停**

## §7 Build 自動化基礎設施

從 `tools/ppt_toolkit/`（reusable Python，不在 skill 內）：

```python
from ppt_toolkit import (
    add_text_with_fallback,    # Latin/CJK font fallback
    fit_image_within,           # 等比 fit + 缺失圖 fallback
    set_speaker_note,           # 上下界字數檢查
    assertion_title,            # 強制 thesis-style 標題
    tier_aware_speaker_note,    # 自動拆 Tier 2/3 + 計算時長
    claude_vision_review,       # 10-checkpoint 表格輸出
)
from ppt_toolkit.style_library_loader import load_object, load_layout
```

screenshot_all.py 跑 wireframe + 結構驗證。

## §8 案例對照（v3 24-slide audit 示範）

`audits/v3_self_phasing_24slide_audit.md`（384 行）。識別 v3 為 12-slide（不是 24，v2 才是 24）；main thesis 推斷為「Self-Phasing 17.3:1 由 4-commit 修補到 1:1，解鎖 ISM 五大目標」（29 字）；6-Q audit 通過率 84.0%；PLOS+AE 通過率 87.5%。

## §9 反例（從 v1/v2/v3 提取）

5-8 個具體案例：
- 「purity 0.95 觸發鎖一開始放 speaker note 而非 slide → 教授可能漏關鍵 root cause」
- 「Rule 5（1 min/slide）全 12 張 FAIL → 系統性問題」
- 「S5/S6 焦點不足 → 可棄用」
- 「S11/S20 多 idea → 拆兩張」
- 「standard label 標題（"結果"/"分析"/"討論"）→ 改 assertion sentence」
- 「speaker note 含「順便提一下」「附帶說明」→ 削減」
- 「主視覺以外另有第二張不相關圖 → 拆/棄」
- 「中文 > 60 字 / 英文 > 30 word per slide → 削減」

## §10 互動 Protocol 模板

| Phase | Prompt 檔案 | 必停 (fast-track) |
|-------|------------|------------------|
| P1 → C1 | `prompts/outline_confirm.md`（含 main thesis 鎖定）| ✅ |
| P2 → C2 | `prompts/outline_confirm.md` | ❌ |
| P3 → C3 | `prompts/section_confirm.md` | ❌ |
| §20 階段 E | `prompts/focal_point_audit.md`（6 問）| ❌ |
| §20 階段 C | `prompts/tier_evaluation.md`（S/A/B/C/D）| ❌ |
| P4 → C4 | `prompts/slide_confirm.md` + `prompts/visual_review.md` | ❌ |
| P5 → C5 | `prompts/visual_review.md`（time-budget alarm）| ✅ |
| from-draft | `prompts/from_draft_loader.md` | (auto) |

## §11 6 模板識別（已實現於 templates/）

詳見 `templates/{improvement_report,comparison_report,executive_summary,data_showcase,concept_walkthrough,academic_defense}.md`。觸發 keyword 對應 SKILL.md 表。

## §12 Style Library 架構

`style_library/`（已落地 30 檔）：
- `colors/palette.yaml` — 6 色 token
- `objects/*.yaml` — 14 物件（caveat_red_strip / insight_green_card / method_blue_box / ...）
- `layouts/*.yaml` — 12 版型（before_after_split / 4quadrant_matrix / ...）
- `examples/case_studies.md` — v3 7 案例對照

每個 object 含 `tier_recommendation`；每個 layout 含 `focal_point_zone`。

## §13 物件 YAML schema

詳見 `style_library/objects/caveat_red_strip.yaml` 為範例。標準欄位：
- name, category, purpose
- visual (background, border, font)
- content_template, required_fields, optional_fields
- position_rules, alignment
- example_good, example_bad
- when_to_use, when_not_to_use
- **tier_recommendation** (S/A/B/C/D，§23)

## §14 顏色 Palette Token

詳見 `style_library/colors/palette.yaml`。6 色 + 用途定義 + WCAG AA + colorblind safe。max 3 colors per slide。

## §15 通用版型分類

每模板對應常見頁面類型，詳見 `templates/*.md` 與 `style_library/examples/case_studies.md`。

## §16 物件對齊規則

- 標題：左對齊，距上 0.4 in，距左 0.5 in
- 主視覺中心：水平 center；垂直 45-55% canvas
- 多欄寬度：2 欄 5.8 in、3 欄 3.8 in、4 欄 2.8 in；統一 0.3 in gap
- caveat strip：對齊被警示元素左邊緣，距下 0.1 in
- citation footnote：右下角，9pt #757575

## §17 Build 調用 API

```python
from ppt_toolkit.style_library import load_object, load_layout

caveat = load_object("caveat_red_strip")
caveat.render(slide, x=1.0, y=5.5, content="⚠ 來源報告寫 38%, 實為 16.5%")

layout = load_layout("before_after_split")
layout.populate(slide, {
    "title": "Baseline tag voting vs V5 三層投票",
    "left.heading": "Baseline (somatic-first)",
    ...
})
```

## §18 物件 / 版型新增 SOP

1. 確認需求（哪個 template / 解決什麼版面問題）
2. 找 example（v1/v2/v3 是否有原型）
3. 寫 YAML（依 §13 schema）
4. 寫例子（good + bad 各 2-3）
5. 跑 wireframe 渲染
6. 進 case_studies.md
7. PR review checklist

## §19 v3 案例對照

`style_library/examples/case_studies.md`（251 行）— v3 7 案例：
- S1 amber 動機 strip → motivation_amber_strip
- S11 紅色根框 + 3 葉樹 → root_cause_tree_with_trigger
- S13 雙欄程式碼 diff → before_after_split + code_diff_box
- S15 主圖 + 數據 caveat → data_main_with_caveat
- S20 業界家族樹 → family_tree_with_2x2
- S22 五大目標卡片 → 5card_grid
- S24 Take-home + Next + Q&A → tldr_3section（待補）

## §20 主軸聚焦與雜訊過濾（核心方法論）★

**6 階段強制過濾**：

### §20.A Main Thesis 鎖定（整份簡報層）
- 整份簡報唯一 main thesis（≤ 30 字）
- 所有 section thesis 必能 derive from main thesis
- 不能直接服務 main thesis 的內容 → backup / oral / 棄用
- §1 Audience & Goal 強制鎖定

### §20.B 每張 Slide Focal Point 鎖定
- 每張 slide 唯一 focal point（≤ 20 字）
- 必須 derive from section thesis
- 必須出現在標題（assertion-evidence headline）
- §3 + §4 強制標註

### §20.C Tier S/A/B/C/D 評分（每筆素材）

| Tier | 名稱 | 評分標準 | 處置 |
|:-:|------|---------|-----|
| S | ESSENTIAL | 直接構成 focal point；無此 slide 失意義 | slide 必現（Tier 1）|
| A | SUPPORTING | 強化 focal point；可不在 slide | speaker note 必講（Tier 2）|
| B | CONTEXTUAL | 提供背景；非當前 focal point 必需 | oral-optional（Tier 3）|
| C | TANGENTIAL | 與 focal point 相關但非必要 | backup / Q&A 預備 |
| D | NOISE | 與 focal point 無直接邏輯關係 | **棄用，不進任何文檔** |

**禁止「先放著看看」— 必先評分後處置**。

### §20.D 比例分流

| 類別 | 比例上限 |
|------|--------|
| Definition slide | ≤ 10% |
| Prerequisite slide | ≤ 15% |
| **Body slide** | **≥ 60%** |
| Conclusion slide | ≥ 15% |

24-slide 範例：Def 2 / Prereq 3 / Body 15 / Conclusion 4 = 24。
**Def + Prereq > 25% → 強制刪減或合併**。

### §20.E 6 問 audit（slide 進 build 前）

詳見 `prompts/focal_point_audit.md`：

1. 此 slide 的 focal point 一句話是？（≤ 20 字，否則拆兩張或棄用）
2. focal point 服務哪個 section thesis？
3. slide 上每個元素都直接支撐 focal point 嗎？
4. 移除哪些元素不影響理解？（→ Tier B/C 或棄用）
5. 哪些細節是 oral-optional？（標 `[ORAL-OPTIONAL]`）
6. **此 slide 移除整體簡報會少什麼？**（若答「沒少什麼」→ 立即棄用）

### §20.F 雜訊紅旗清單

slide / speaker note 出現以下視為紅旗：
- 「順便提一下」「附帶說明」「另外」（非主線連接）
- 「值得注意的是」（非當前 focal point）
- 「供參考」「這部分如果有時間可以再細講」（speaker note → **直接棄用該段**）
- 主視覺以外另有第二張不相關圖
- 多於 3 個 bullet point（拆 slide 或棄用）
- 標題 generic label（「結果」「分析」「討論」）

## §21 §20 與既有章節整合

| 既有章節 | §20 整合點 |
|---------|---------|
| §1 Audience & Goal | 加 Main Thesis 鎖定必填 |
| §2 Outline checkpoint | 加 section thesis derive from main thesis 檢查 |
| §3 Section checkpoint | 加每 slide focal point 草擬步驟 |
| §4 Slide checkpoint | 加 §20.E 6 問 + §20.F 紅旗檢查 |
| §5 Visual review 10-checkpoint | 第 2 條升級為「1 focal point + §20.B 通過」|
| §10 互動 protocol | 加 focal point 確認 + Tier 評分 prompt |

## §22 v3 24-slide audit（示範）

`audits/v3_self_phasing_24slide_audit.md` 為 §20 audit 示範。
- main thesis: 「Self-Phasing 17.3:1 由 4-commit 修補到 1:1，解鎖 ISM 五大目標」（29 字）
- 6-Q audit 通過率 84.0%
- PLOS+AE 通過率 87.5%
- 5 個雜訊紅旗實例
- 字數削減目標 21%（13,847 → 10,950）

## §23 物件層 Tier / Focal Point 整合

```yaml
# style_library/objects/caveat_red_strip.yaml
tier_recommendation: A  # SUPPORTING
when_to_promote_to_S: caveat 直接構成 main thesis 一部分
when_to_demote_to_B: caveat 屬背景補充而非當前 focal point 必需

# style_library/layouts/before_after_split.yaml
focal_point_zone:
  position: "delta_strip (bottom center)"
  description: "Δ metric 作為視覺與認知焦點"
  rule: "若 delta 不是 focal point，請改用其他版型"
```

## §24 PPTX 三件套輸出

從 weekly-report PPTX_PROTOCOL.md 搬入（DEPRECATED 後已搬到本檔）：

`{output_dir}/`:
- `00_storyboard.md` — Three-Act 三幕故事弧
- `01_full_narrative.md` — P2 outline 全敘事
- `02_slide_outline.md` — P3 section list 含 thesis sentence
- `03_slide_layout_and_script.md` — P4 slide 三層分流 + P5 speaker script
- `build_pptx.py` — 生成腳本
- `output.pptx` — 最終 PPTX
- `wireframes/*.png` — 截圖驗證

設計常數（從 PPTX_PROTOCOL.md 搬入）：
- 色彩：dark=#1E2A44, bg=#F7F3EC, accent=#A85540, 綠=#009E73, 紅=#D55E00
- 字體：標題 Arial Bold ≥ 32pt, 內文 Arial ≥ 14pt, 下限 9pt
- 版面：視覺 55-65%, 留白 20-30%, 文字 ≤ 15%, 每頁 ≤ 4 bullet
- 雙語：中文主 + 英文副（60% 字號 + 縮排 0.25"）
- 7 樣本固定順序：HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829
- PPTX 圖片：fit-within + centered

### Phase 5 多 Agent 驗證（PPTX 6 Agent，從 weekly-report Phase 5B 搬入）

Wave 1（結構+視覺）：
- Agent-T（字體）+ Agent-C（色彩）+ Agent-L（佈局）+ Agent-B（雙語）

Wave 2（內容+整合）：
- Agent-S（Speaker Notes）+ Agent-D（數據準確性）

每波平行執行（`feature-dev:code-reviewer` subagent_type），波次間依序修正。

## §25 success criteria（驗收標準）

| 指標 | 通過標準 |
|------|---------|
| Skill 自動載入 | 觸發 keyword 自動載入 |
| 三層確認強制性 | 跳過 outline / section / slide checkpoint 會被擋下 |
| 6 模板識別 | 觸發詞自動載入對應 template |
| Tier 分流落實 | speaker note 出現 `[ORAL-OPTIONAL]`；slide 文字 ≤ 60 字 |
| Focal point 標註 | 每張 slide 有 ≤ 20 字 focal point 註記 |
| 雜訊紅旗 | §20.F 紅旗清單在 prompts 中可自動掃描 |
| Style library 完整 | 14 objects + 12 layouts，每個有 good/bad + when_to_use |
| ppt_toolkit 可用 | `from ppt_toolkit import load_object` 可運作 |
| --from-draft 介面 | 自動讀 frontmatter，跳過 P1 main thesis 鎖定 |


---

## §24.1 6 Agent 多代理驗證詳細規範（v2.3 補強）

playbook §24 列了 6 Agent 名稱，**詳細 context + 步驟 + Vision 整合 → `prompts/multi_agent_review.md`**。

### 三層保險機制（每張 slide）

| 層 | 工具 | 何時跑 | 通過標準 |
|:-:|------|--------|---------|
| 1 | §20.E 6 問 self-audit | build 前 | 6 問全有答 + focal point ≤ 20 字 |
| 2 | visual_review 10-check（單 Agent + Vision PNG）| build 後 | ≥ 8/10 PASS |
| 3 | multi_agent_review Wave 1（4 Agent 並行）| build + 截圖後 | ≥ 16/20 PASS |
| 4 | multi_agent_review Wave 2（2 Agent 並行）| 整份完成後 | 數據錯誤=0 + caveat 缺=0 + 雜訊=0 |

### 6 Agent 並行調度

```python
# Wave 1: 每張 slide 4 Agent 並行（feature-dev:code-reviewer subagent_type）
for slide_n in range(N):
    parallel_results = parallel_spawn([
        Agent(type='feature-dev:code-reviewer', ctx='Agent-T 字體', png=f'slide_{slide_n}.png'),
        Agent(type='feature-dev:code-reviewer', ctx='Agent-C 色彩', png=f'slide_{slide_n}.png'),
        Agent(type='feature-dev:code-reviewer', ctx='Agent-L 佈局', png=f'slide_{slide_n}.png'),
        Agent(type='feature-dev:code-reviewer', ctx='Agent-B 雙語', png=f'slide_{slide_n}.png'),
    ])

# Wave 2: 整份完成後 2 Agent 並行
final_results = parallel_spawn([
    Agent(type='feature-dev:code-reviewer', ctx='Agent-S Speaker Notes', files=['03_slide_layout_script.md']),
    Agent(type='feature-dev:code-reviewer', ctx='Agent-D Data Accuracy', files=['01_*', '02_*', '03_*', master_draft]),
])
```

每個 Agent 獨立 context，避免污染（關鍵設計）。

### Agent 觸發語（在 multi_agent_review.md 詳列）

| Agent | 必看 | 重點檢查 |
|-------|------|--------|
| T 字體 | PNG + thesis_title_bar.yaml | 字級、Latin/CJK fallback、雙語縮排 |
| C 色彩 | PNG + palette.yaml | 6 色 token、WCAG AA、colorblind |
| L 佈局 | PNG + {used_layout}.yaml | 對齊、邊距、focal_point_zone、≤6 element |
| B 雙語 | PNG | 中文 ≤ 60 字、英文 ≤ 30 word、60% 字級 |
| S 講稿 | speaker_script.md + master_draft | 75-90 sec/slide、[ORAL-OPTIONAL] |
| D 數據 | 全 PPTX 檔案 + source_artifacts | 數字一致、caveat 完整、無過度宣稱 |

### 修正循環

每 Wave 後：
- PASS → 進下階段
- FAIL → 主 conversation 修正後 re-spawn 該 Agent re-review

詳見 `prompts/multi_agent_review.md` §「修正循環」表。
