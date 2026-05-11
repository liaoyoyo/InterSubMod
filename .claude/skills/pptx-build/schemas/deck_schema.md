---
title: deck.yaml schema 規格說明 (v1.0.0)
type: pptx_schema_spec
date: 2026-05-08
status: pending_review
---

# deck.yaml schema 規格 (v1.0.0)

> 此檔案搭配 `deck_schema.yaml`（範例）使用。本檔解釋每個欄位語意、必填/選填、與既有規則（R-G1/G2/G3 active）的對應。
>
> **2026-05-08 Stage 2.1 寫入。等用戶 review ack 後才動 Stage 2.3 (parser/renderer 重構)**。

---

## 1. 設計原則（用戶 5/08 已 ack 的 5 個決策）

| 決策 | 對應 schema 設計 |
|---|---|
| Q1=(c) YAML frontmatter + Markdown body 混用 | `slides[].body` / `tier2_note` 用 `\|` literal block scalar 接 Markdown |
| Q2=(a) 只 PPTX renderer | `output.pptx` 為主輸出；HTML/Reveal.js 不在此 schema |
| Q3=(c) SCQA 兩層 | deck-level `narrative_framework` + slide-level `tier`/`focal_point` |
| Q4=(b) Stage 3 暫不排 | schema 不含 reveal_js / canva_target 欄位 |
| Q5=(c) Claude Design 遷移 | `handoff.claude_design_export_target` 預留 |

額外原則：
- **schema_version**：每次 break-change 升 minor；新增 optional 欄位升 patch
- **R-G1/G2/G3 active 規則**：內嵌於 `validation_rules` section，由 deck_renderer.py 強制驗證
- **保留向下相容**：v1.0.x 內 schema 改動不破壞既有 deck.yaml

---

## 2. 7 個 sections 完整定義

### Section 1: `deck` — Deck-level metadata

| 欄位 | 必填 | 型別 | 說明 |
|------|:-:|------|------|
| `title` | ✓ | str | 簡報主標題 |
| `subtitle` | – | str | 副標題（可空）|
| `audience` | ✓ | enum | PI / advisor / committee / lab / public |
| `audience_familiarity` | – | list[str] | 受眾熟悉/不熟悉領域標明，影響 R-G2 術語密度 |
| `target_minutes` | ✓ | int | 預設時長（min）|
| `total_slides_estimate` | – | int | 估計 slide 數 |
| `date` | ✓ | str | 報告日期 ISO 格式 YYYY-MM-DD |
| `period` | – | str | 涵蓋期間 |
| `main_thesis` | ✓ | str | ≤ 30 字（C1 必停強制）|
| `report_type` | ✓ | enum | 對應 `templates/{type}.md` 6 模板之一 |
| `glossary` | – | list[{term, definition}] | Mini-glossary，每條 ≤ 30 字 |

### Section 2: `narrative_framework` — Deck-level SCQA/Pyramid

| 欄位 | 必填 | 說明 |
|------|:-:|------|
| `type` | ✓ | SCQA / Pyramid / IMRAD / Custom |
| `situation` | (SCQA 必)| Situation |
| `complication` | (SCQA 必)| Complication |
| `question` | (SCQA 必)| Question |
| `answer` | (SCQA 必)| Answer |
| `pyramid_top` | – | Pyramid 頂層答案 |
| `pyramid_supporting` | – | list[str] 3 supporting points |

選一即可，executive_summary 模板必填 `pyramid_*`。

### Section 3: `source_artifacts` — 全 deck source-of-truth refs

每筆：
| 欄位 | 必填 | 說明 |
|------|:-:|------|
| `path` | ✓ | InterSubMod/... 相對路徑（依 doc-standards）|
| `role` | ✓ | critical / supporting / context |

### Section 4: `slides` — 內容主體

每張 slide：
| 欄位 | 必填 | 說明 |
|------|:-:|------|
| `id` | ✓ | S01, S02, ... 唯一識別 |
| `section` | ✓ | section 名稱（對應 P3 section batch）|
| `layout` | ✓ | 對應 `style_library/layouts/{layout}.yaml` |
| `focal_point` | ✓ | ≤ 20 字（§20.B 強制）|
| `tier` | ✓ | S/A/B/C/D（§20.C；D 應已棄用不會在此）|
| `title` | ✓ | assertion-evidence headline ≤ 30 字 |
| `title_color` / `title_bg_color` | – | 用 palette token 名 |
| `body` | ✓ | Markdown，slide 上內容（中 ≤ 60 字 / 英 ≤ 30 word）|
| `tier2_note` | ✓ | speaker note Markdown，300-360 字 |
| `tier3_oral` | – | 純 oral 段，可選 |
| `evidence` | ✓ | list[evidence_entry]，R-G1 強制 |
| `figure` | – | str 圖片相對路徑或 null |
| `footer` | – | citation footnote（≤ 9pt 渲染）|

`evidence_entry` 結構：
| 欄位 | 必填 | 說明 |
|------|:-:|------|
| `claim` | ✓ | 此 evidence 支撐的具體 claim |
| `source_path` | ✓ | 來源 .md/.cpp/.tsv 路徑（`InterSubMod/...`）|
| `source_section` | ✓ | 段落/行號（如「§3.1 line 53」「PhasingProcess.cpp:142-220」）|
| `scope` | ✓ | **R-G1 強制**：whole genome / threshold_compare 全 BAM / 15-site cherry / multi-sample 等 |
| `verified_date` | ✓ | grep 確認日期 ISO 格式 |
| `confidence` | – | high / medium / low |
| `note` | – | 補充說明（如全 BAM 預期值）|

### Section 5: `output` — 輸出檔案配置

| 欄位 | 必填 | 預設 |
|------|:-:|------|
| `pptx` | ✓ | "output.pptx" |
| `speaker_script` | – | "speaker_script.md" |
| `audit` | – | "audit.md" |
| `figures_dir` | – | "figures/" |
| `wireframes_dir` | – | "wireframes/" |

### Section 6: `validation_rules` — R-G1/G2/G3 active 規則內嵌

| 欄位 | 預設 | 對應規則 |
|------|---|---|
| `metric_source_required` | true | R-G1 |
| `scope_label_required` | true | R-G1 |
| `unfamiliar_terms_per_slide` | 3 | R-G2 |
| `glossary_required_for_excess` | true | R-G2 |
| `wave1_agent_n_enabled` | true | R-G3 |
| `body_chinese_max_per_slide` | 60 | playbook §6 |
| `body_english_max_per_slide` | 30 | playbook §6 |
| `speaker_note_chinese_min` | 300 | playbook §6 |
| `speaker_note_chinese_max` | 360 | playbook §6 |
| `speaker_note_chinese_max_total` | 10000 | 25 min × 400 字/min |

deck_renderer.py 在 build PPTX 時會讀此 section 強制驗證每張 slide。違反規則應 raise ValidationError 或 warn。

### Section 7: `handoff` — 與 weekly-report / Claude Design 銜接

| 欄位 | 必填 | 說明 |
|------|:-:|------|
| `source_master_draft` | – | weekly-report 上游母稿路徑 |
| `weekly_report_handoff` | – | A/B/C/D（4 選 handoff）|
| `next_week_plan_path` | – | handoff=D 時填 |
| `claude_design_export_target` | – | Stage 4 export bundle 目標檔名 |

---

## 3. Markdown body 解析規則（YAML literal block scalar 內的 Markdown）

`slides[].body` 與 `tier2_note` 使用 YAML `|` 保留換行：

```yaml
body: |
  ## Caveat
  [!] 4/29 PI 報告引用的所有 V5 數值都是 Pass 1 only 條件
```

deck_parser.py 解析規則：
1. **`## Section` heading**：對應 layout zone 名（`data_main_with_caveat` 含 zone：caveat / 因果鋪墊 / Evidence / 關鍵數字）
2. **`[!] xxx`**：caveat strip（紅）
3. **`[OK] xxx`** / **`[X] xxx`**：tag 樣式（不渲染為 emoji，純 ASCII，R-G1 強制）
4. **`- item`**：bullet point
5. **table（`| ... |`）**：渲染為 PPTX table
6. **`*footnote*`**：表格 footer（小字 9pt 灰）
7. **`[ORAL-OPTIONAL] xxx`**：tier2_note 中的 Tier 3 段落標記

不支援的 Markdown 元素（會 warn）：
- 圖片 inline `![](path)` — 改用 slides[].figure 欄位
- HTML embed
- 巢狀 list（≥ 3 層）

---

## 4. 與既有 build_pptx.py 對照

| 既有 build_pptx.py | deck.yaml 對應 |
|---|---|
| 800-1500 行 inline Python 寫每張 slide | `slides[].body` Markdown 描述內容 |
| `add_assertion_title(s, "...")` | `slides[].title` |
| `add_caveat_strip(s, ..., "[!] ...")` | `body` 中 `## Caveat` section + `[!] ...` |
| `add_method_box / add_neutral_box / add_insight_card` | `body` 中對應 `## Method` / `## Evidence` / `## Insight` heading |
| `set_speaker_note(s, "...")` | `slides[].tier2_note` |
| `add_focal_marker(s, "...")` | `slides[].focal_point` |
| `fit_image_within(s, FIGURES_DIR / "...")` | `slides[].figure` |
| `add_footer_citation(s, "...")` | `slides[].footer` |

deck_renderer.py 負責把 deck.yaml 對應到既有 ppt_toolkit primitives，無需重寫渲染邏輯。

---

## 5. R-G1/G2/G3 active 規則整合

### R-G1 — Metric source verification

每個 `evidence_entry` 必須有：
- `source_path` ✓
- `source_section` ✓
- `scope` ✓（whole genome / cherry-pick / etc.）
- `verified_date` ✓

deck_parser.py 在 load 時驗證；缺項 raise ValidationError。

### R-G2 — PI 不熟術語 ≤ 3 / slide

deck_renderer.py 解析 `slides[].body`，比對 `deck.audience_familiarity` 與 `deck.glossary`，計算每張 slide 的「PI 不熟術語」數。
- ≤ 3：PASS
- 4-5：warn + 建議補 footnote
- ≥ 6：FAIL（強制要求 glossary 或 footnote）

### R-G3 — Wave 1 Agent-N（數據驗證）

`validation_rules.wave1_agent_n_enabled: true` 啟用後，Wave 1 review 自動產生 Agent-N prompt，要求對每張 slide `evidence` 列每筆做 grep 驗證。

---

## 6. 範例 deck.yaml 完整結構

見 `deck_schema.yaml`（同目錄）— 含 deck metadata + narrative_framework + 2 範例 slides（S04 + S07）+ output + validation_rules + handoff。

實際 deck（試用第 2 次 18-slide）將在 Stage 2.5 落地於 `InterSubMod/docs/reports/validated/2026/05/20260505_週報_self_phasing_v2.1試用第2次/pptx_build/deck.yaml`。

---

## 7. Schema 升級政策

| 升級類型 | 例子 | 版號規則 |
|---|---|---|
| 新增 optional 欄位 | 加 `slides[].animation` | patch (1.0.x) |
| 修改既有欄位語意 | `tier` 從 enum S/A/B/C/D 改 1-5 數字 | minor (1.x.0) |
| 移除欄位 / 結構性破壞 | 把 `slides` 改為 `sections[].slides[]` | major (2.0.0) |

deck_parser.py 開頭應有 `assert deck['schema_version'].startswith('1.')` 守門，未來 v2 schema 需另寫 v2 parser。

---

## 8. 待用戶 review 的 confirm 點（C2=a 暫停）

請對 Section 1-7 schema 確認：

1. **deck-level metadata** 欄位是否完整？是否需要加什麼欄位？（例如：tags / version / collaborators）
2. **narrative_framework** SCQA/Pyramid 雙層必要嗎？或只保留 SCQA 即可（Pyramid 用 executive_summary 模板自帶）
3. **slides[].body Markdown** 解析規則 7 條（## section / [!] / [OK] / table / etc.）是否充分？有沒有你常用但漏掉的格式
4. **evidence_entry** 5 必填欄位（claim / source_path / source_section / scope / verified_date）是否合理？scope 是否要做 enum 限制（whole_genome / cherry_picked / multi_sample / single_sample）
5. **validation_rules** 內嵌規則是合理（內嵌不同 deck 可微調）還是應 hard-code 在 renderer（避免被 deck 作者繞過）？
6. **glossary** 結構（list of `{term, definition}`）vs 其他格式（dict / 單行 yaml inline）
7. **schema_version "1.0.0"** 起點是否合適？

我建議答案：**1=可加 collaborators / 2=保留雙層（彈性）/ 3=充分 / 4=合理且 scope 用 enum / 5=內嵌（每 deck 可微調）/ 6=list of dict / 7=合適**

請對 1-7 給答覆（「同意全部」或標出特定點），我即動 Stage 2.3 (deck_parser.py + deck_renderer.py 重構) + Stage 2.5 (試用第 2 次 18 slide 重做為 deck.yaml)。
