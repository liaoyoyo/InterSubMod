# Design Principles & Pre-publish Checklist

> 12 規則 + 5 秒測試。整合 Tufte / Garr Reynolds / Nancy Duarte / Robin Williams CRAP / Assertion-Evidence / Nature figure guidelines / Nielsen Norman / WCAG。
> Skill prompts 在生成 HTML/slides 前**應該先 Read 本檔**。

---

## A. 訊息層（What to communicate）

### Rule 1 — 5 秒測試（NN/g + Berkeley dataviz checklist）
- **Test**: 給同事看 5 秒，能否說出 takeaway？答不出 → hierarchy 或 focal point 不夠。
- **Apply to InterSubMod**: 每張 slide / 每個報告 section header / 每張 SVG figure 都要過 5 秒測試。
- **Tool in skill**: `stat-grid` 4-6 個焦點數字、`verdict-banner` 一句話結論、`conclusion-arrow` 段落結尾。

### Rule 2 — 3 秒法則 / Glance media (Nancy Duarte)
- **Source**: Duarte《slide:ology》— slide 是 billboard，不是 document
- **Test**: 標題 + 主視覺 3 秒內傳達主訊息（不需讀內文）
- **Apply**: slide title ≤ 15 字 (中) / 80 字元 (英)；不要把段落塞進 slide
- **Anti-pattern**: "Method" 這種空標題；改為 "Method: Cross-Sample Validation Removes Cell-Line Bias"

### Rule 3 — Assertion-Evidence 結構（Penn State, NSF-backed）
- **Source**: Michael Alley《The Craft of Scientific Presentations》
- **Rule**: Slide 標題 = 一句完整斷言（聲明結論）；slide 內容 = 視覺證據支持
- **Effect**: 認知負擔降低（Kuwait Univ 2025 實驗 108 學生 lower cognitive load）
- **Apply to slide_prompt.md**:
  ```
  ❌ "結果"
  ✅ "V6 修正 100% chr19 priority bug victims (752/752 reads)"
  ```
- **Apply to standalone_prompt.md**: `<h2>` 標題本身就是該段結論句，不只是 topic 標籤。

---

## B. 視覺層（How to design）

### Rule 4 — Tufte data-ink ratio（最大化 data-ink，移除 chartjunk）
- **Source**: 《The Visual Display of Quantitative Information》(1983, 2001 expanded)
- **Rule**: data-ink / total-ink 比例越高越好；刪掉所有不傳達數據的 ink
- **Apply**:
  - 圖表去掉 gridline、3D 效果、漸層背景、不必要 legend
  - Bar chart 直接寫數字在 bar 旁邊（取代 y-axis tick）
  - 表格去掉所有 inner border（只留 header underline）— "Plain table" 風格
- **InterSubMod hit**: 既有 `metric-table` 有 inner border — **未來可選用 Tufte plain 風格**（standalone mode）

### Rule 5 — CRAP 四原則（Robin Williams 1993）
| Principle | Rule | Apply in skill |
|-----------|------|----------------|
| **Contrast** | 強調用大小/字型/顏色差異 | `.stat-box .number` 32px 對比 14px label；`.callout-critical` 紅框對比白底 |
| **Repetition** | 視覺元素一致重複建立 cohesion | 8 section colors 在整份報告 / deck 一致；`stat-box` 用同一個 class |
| **Alignment** | 物件靠 grid 對齊，不浮散 | `.layout` grid-template-columns: 220px 1fr；`.stat-grid` repeat(auto-fit, minmax(140px, 1fr)) |
| **Proximity** | 相關物件靠近，無關物件用空白分隔 | `.errata-card` 內 summary + body 緊湊；卡片之間 `margin: var(--sp-4)` |

**Self-check before output**: 開瀏覽器看 — 視線是否被引導到重點？相關元素是否「群聚」？

### Rule 6 — 留白 / Whitespace (Garr Reynolds Presentation Zen)
- **Rule**: 留白不是「沒設計」— 留白本身就是設計
- **Apply**:
  - Slide 不要把畫面填滿；canvas padding 28px 40px 預留邊界
  - Report 段落間距 ≥ var(--sp-5) (24px)
  - Section heading 上 margin var(--sp-7) (48px)
- **Anti-pattern**: footer-glossary 擠到 canvas 邊緣；text 撞到 border

### Rule 7 — Hierarchy（視覺層級）
- **Source**: Berkeley BPM dataviz checklist + NN/g
- **Rule**: 最重要的元素應有最強視覺權重（大小 / 顏色 / 位置 / 留白）
- **Apply in standalone**:
  1. Verdict banner（最頂、gradient）
  2. Stat-grid（焦點數字，醒目色）
  3. Core errata cards（open by default、紅框）
  4. Minor errata cards（collapsed、無強調）
  5. Footer references（最小、灰色）

---

## C. 色彩 + 可及性

### Rule 8 — 1-2 primary + 1 accent（dataviz best practice 2026）
- **Rule**: 整份報告只用 1-2 個主色 + 1 個強調色；色彩濫用 = 視覺噪音
- **Apply in skill**: design_tokens.css 已限制到 4 個 semantic color + 1 accent
- **InterSubMod 偏好**: Anthropic palette accent #D97757；semantic = critical/warning/success/info
- **Anti-pattern**: rainbow chart with 8+ colors（用 viridis sequential 或 small multiples 取代）

### Rule 9 — Colorblind-safe（WCAG + Nature 2026 expectation）
- **Source**: Nature 2025-2026 figure guidelines + WCAG 2.1
- **Rule**: 不靠**單一色彩**編碼意義；同時用 pattern / shape / label / position
- **Apply**:
  - row-bug / row-safe / row-warn **加上文字 verdict 欄**（✓/✗/⚠）
  - Stat-box `.number` 紅綠對比 → 加文字 label "bug" / "fixed"
  - SVG bar chart 顏色不同 → 加 numeric label 在 bar 右側
- **Test**: 用 Coblis simulator 看色盲版本是否仍能讀懂
- **Palette 建議**: viridis / cividis / blue-orange（替代 red-green）

### Rule 10 — WCAG contrast ratio ≥ 4.5:1 (AA) / 7:1 (AAA)
- **Apply**: 內文 text 對 bg 至少 4.5:1
  - `--c-text: #141413` on `--c-bg: #FAF9F5` → 15.8:1 ✅
  - `--c-text-soft: #6B6B66` on `--c-bg: #FAF9F5` → 5.2:1 ✅
  - `--c-accent: #D97757` on `--c-bg: #FAF9F5` → 3.4:1 ⚠（限定大字 / 非內文用）
- **Test**: Chrome DevTools Lighthouse accessibility audit

---

## D. 結構 + 一致性

### Rule 11 — Vector first（Nature / Cell figure guideline）
- **Rule**: 圖表優先用 SVG（vector）而非 PNG raster
- **Why**: 列印無鋸齒、可編輯、可 zoom、檔案小
- **Apply in skill**: `components/svg_*.html` 三類（flow / compare bar / icon）已提供
- **PNG only when**: 真實數據 scatter >50 點、IGV 截圖、照片

### Rule 12 — Pre-publish checklist（必過 6 項）

在 Write 輸出 HTML 前自審：

```
☐ 1. Decision intent clear — first-time viewer 看了知道該做什麼決定？
☐ 2. 5 秒測試 — 主訊息 5 秒內讀得到？
☐ 3. Chart-to-task fit — line 趨勢 / bar 比較 / scatter 相關 用對了？
☐ 4. CRAP — Contrast / Repetition / Alignment / Proximity 都做到？
☐ 5. Colorblind-safe — 不靠色彩單獨編碼？有 pattern/label 加強？
☐ 6. Hierarchy — 最重要元素有最強視覺權重？
```

任一條 ✗ → 重新生成 / 補強，不要 ship。

---

## 對應 skill 元素的 Quick Reference

| 業界原則 | InterSubMod 對應元素 | Source |
|---------|---------------------|--------|
| 5 秒測試 | `.stat-grid` (3-6 焦點) + `.verdict-banner` | NN/g, Berkeley BPM |
| 3 秒 glance | Slide title ≤ 15 中字 | Duarte |
| Assertion-Evidence | Slide title = 結論句；body = 視覺證據 | Alley, Penn State |
| Data-ink ratio | SVG 無 gridline / 表格少邊框 | Tufte |
| CRAP-Contrast | `.stat-box .number` 32px vs label 14px | R. Williams |
| CRAP-Repetition | 8 section colors 一致使用 | R. Williams |
| CRAP-Alignment | `.layout` grid + `.stat-grid` auto-fit | R. Williams |
| CRAP-Proximity | `.errata-card` summary + body 緊湊 | R. Williams |
| 留白 | canvas padding + section margin | G. Reynolds |
| Hierarchy | verdict > stat-grid > core cards > minor cards | NN/g |
| 1-2 primary + 1 accent | 4 semantic + 1 accent palette | dataviz BP 2026 |
| Colorblind-safe | row-class + text verdict 雙重編碼 | Nature 2026, WCAG |
| WCAG ≥ 4.5:1 | tokens 內文配色 15.8:1 通過 | WCAG 2.1 |
| Vector first | `components/svg_*.html` | Nature, Cell |

---

## 6-Taboo audit 與 Design Principles 的關係

`SKILL.md §6-Taboo Audit` 是「禁忌清單」（負面表列）。
`references/design_principles.md`（本檔）是「設計準則」（正面表列）。

**Workflow**：
1. 生成 HTML 草稿時 → 套用設計準則（Rule 1-12，正面）
2. Write 輸出前 → 6-Taboo audit grep（負面 reject）
3. 兩關都過 → ship

兩者互補：「正面準則告訴 LLM 怎麼做對」、「負面禁忌告訴 LLM 不要做錯」。

---

## Rule 13 — 場景對應量化標準（2026-05-20 Issue #4 升級）

> 不同 PI 報告場景對應不同 slide 數 / chars / rows / visual ratio 上限。生成前必對齊 frontmatter `audience_scenario`。

### 場景量化標準表

| 場景 | 時長 | Slide 數 | 純文字 chars/slide | 含表 chars/slide | Table rows | Main visual 佔比 | Reading test |
|------|:----:|:-------:|:------------------:|:----------------:|:----------:|:----------------:|:------------:|
| **PI 1-on-1（週報）** | 5-10 min | 3-6 | ≤ 150 | ≤ 250 | ≤ 5 | 60-75% | 3-second |
| **Lab meeting** | 15-30 min | 8-15 | ≤ 250 | ≤ 350 | ≤ 8 | 50-65% | 5-second |
| **Conference talk** | 45-60 min | 20-40 | ≤ 350 | ≤ 450 | ≤ 12 | 40-60% | 10-second |
| **Lab informal pitch** | 5 min | 3-5 | ≤ 100 | ≤ 200 | ≤ 4 | 70-80% | 3-second |

### 配套：speaker note 機制對齊

- **PI 1-on-1 / Lab meeting**：speaker note 用 `<details>` 折疊（`components/speaker_note_details.html`）
- **Conference**：speaker note 用 visible block（`components/speaker_note.html`）
- **Lab informal**：可省略 speaker note，slide 已含全部資訊

### 落地 audit

每張 slide 生成後 self-audit：
```
chars_visible = count_chars(slide content excl <details>)
table_rows = count(<tr>) per slide
if chars_visible > scenario_limit OR table_rows > scenario_row_limit:
    flag for compression / split / move to speaker note
```

### Cross-reference
- Memory: `reference-pi-scenario-quantitative-standards.md`
- Audit report: `InterSubMod/docs/solutions/optimization/2026/05/20260520_pptx_workflow_skill_audit_01.md` Issue #4
- Companion: [[feedback-ppt-minimal-visual-first]]

---

## Rule 14 — Interactive Dashboard / Web-UI Patterns（2026-06-03 新增；補 Rule 1-13 print/figure 導向之互動 HTML 缺口）

> **適用**：互動 HTML 儀表板 / 工作板（`state/focus_board.html`、goal-landscape workboard、harness_health dashboard）—— 非單向 print/slide，而是**可探索、可收納、可追蹤變化**的介面。Rule 1-13（Tufte/CRAP/Nature）仍全適用；本條補「互動層」。
> **總綱（Shneiderman's Mantra）**：**Overview first, zoom & filter, details on demand.** 先全貌 → 再下鑽。

### 14.1 漸進揭露（Progressive Disclosure）
- 預設只露**當下決策必需**；其餘收摺（native `<details>` 零 JS，或 `.collapsed` class）。
- 分層揭露對應既有 L0-L3：**L0 一眼焦點 → L1 主軸 → L2 關聯/譜系 → L3 細節（預設收）**。
- 反模式：一次攤開所有資訊（牆）= 認知過載（NN/g cognitive load）。

### 14.2 資訊密度層級 + WIP-limit
- 高資訊頁用**密度梯度**：hero（最大）> 卡 > 列 > 收摺細節；非均一密度。
- 焦點槽 **WIP-limit ≤2**（逼聚焦，借 Kanban / Disco Thought Cabinet）。
- 提供**密度切換**（寬鬆/緊湊 toggle）讓用戶自選資訊密度（Refactoring UI: "design for density choice"）。

### 14.3 可掃描性（Scannability, F-pattern）
- section 標題用**問句**（「現在做什麼?」「結果怎麼累積?」「卡在哪?」）= 資訊氣味（information scent）。
- 視覺權重梯度引導眼睛：accent hero → section header → body。
- 重要狀態用**雙重編碼**（色 + icon/文字，呼應 Rule 9 colorblind）：✓ 綠 / ⚠ 黃 / ✗ 紅 + 文字。

### 14.4 可操作性 affordance + 即時回饋
- 可點擊元素必有 affordance：`cursor:pointer` + `:hover` 變化 + chevron（▾/▸ 旋轉）標收合態。
- 狀態回饋即時：toggle 後立刻反映；無「點了沒反應」。
- 反模式：看起來可點但不可點 / 可點但無 hover 提示。

### 14.5 Responsive grid + overflow 安全
- grid 用 `repeat(auto-fit, minmax(...))` 或 `@media` 斷點；窄螢幕降單欄。
- **overflow 鐵則**：flex/grid 子項加 `min-width:0`（經典 flexbox 溢出根因）+ 長文字 `overflow-wrap:anywhere; word-break:break-word`。長 id/title 放 `title=` tooltip + short label 顯示。

### 14.6 Detail-on-demand（下鑽）
- 摘要 + 點開細節（`<details>` / tooltip / 連結到 source 檔）。
- 數字/結論可追溯：hover/點擊 → 來源（呼應 §13 provenance；如 lineage node 附 ledger entry id）。

### 14.7 狀態與「變化/進退」回饋（temporal）
- 進度用 **stepper（P0-P6 階段點，current 高亮）** + **roll-up %**（goal← cycle phase），非假甘特日期（Rule：研究無確定日期勿假精確）。
- **「自上次檢視變化」diff**（localStorage 存上次快照 → 比對 → 標 ▲推進 / ▼退 / +新）= 直接服務「任務進退理解」。⚠ localStorage 僅 UI 偏好/快照，非 SoT（§13 + [[feedback_self_phasing_PI_no_build_html_py]]）。

### 14.8 互動層 restraint（與 Rule 1-13 一致的克制）
- **inline vanilla JS 優先**（零 CDN / 零 framework / 零 build）→ 維持 offline 單檔、可 grep、可版控。需重互動（live filter / sortable table）才考慮 Alpine（已在 PI stack，但引 CDN/vendor 依賴）；**永不** React/Vue/webpack（破壞單檔特性）。
- 數字一律 **by-construction 注入**（derive script 讀真檔 render，§13 Layer A），非手打。
- 遊戲化只取「服務聚焦」的部分（HUD pin / 進度條 / 安靜 checkmark）；**剝離** celebration 特效 / 集點 / 排行榜 / RNG / 隱藏成就（§13.7 未驗證滿足感 + Goodhart）。

### Dashboard 專屬 pre-publish 檢查（疊加 Rule 12 的 6 項）
- [ ] **5 秒儀表板測試**：一眼能答「現在該做什麼 + 下一步」？
- [ ] 預設收摺態合理（細節收、焦點露）？
- [ ] 所有可點元素有 affordance + hover？
- [ ] 窄螢幕不溢出（min-width:0 + overflow-wrap）？
- [ ] 每個數字可追溯真檔（derive，非手打）？

### 對應 focus_board.html 元素
| Rule 14 | focus_board 實作 |
|---------|------------------|
| 14.1 漸進揭露 | L0-L3 + `.sec` collapsible + `<details>` lineage |
| 14.2 密度 + WIP | hero>卡>列梯度 + pinned≤2 + ⇲緊湊 toggle |
| 14.3 可掃描 | 問句 section header + accent hero + 雙重編碼 badge |
| 14.4 affordance | `.sec-h` cursor:pointer + hover + chevron 旋轉 |
| 14.5 responsive | grid auto + `min-width:0` + `overflow-wrap` |
| 14.7 進退 | phase stepper + goal roll-up % + since-last-view diff ✓（focus_board sincebar localStorage 已實作 2026-06-03）|
| 14.8 restraint | inline vanilla JS 零 CDN + Layer A derive |

---

## Sources

- Edward Tufte — *The Visual Display of Quantitative Information* (1983, 2001)
- Garr Reynolds — *Presentation Zen* + *Presentation Zen Design*
- Nancy Duarte — *slide:ology* + *Resonate*
- Robin Williams — *The Non-Designer's Design Book* (CRAP 1993)
- Michael Alley — *The Craft of Scientific Presentations* (Assertion-Evidence)
- Penn State Assertion-Evidence approach (NSF-backed)
- Nature / Cell / Science figure guidelines (2025-2026)
- WCAG 2.1 contrast guidelines
- Nielsen Norman Group dataviz heuristics
- Berkeley BPM dataviz checklist (PDF)
- Depict Data Studio dataviz checklist
- dataviz best practices 2026 (Techment, Julius AI, Omni Analytics)
- **（Rule 14 互動層）** Ben Shneiderman — "Overview first, zoom & filter, details on demand" mantra
- **（Rule 14）** Refactoring UI — Adam Wathan & Steve Schoger（density, hierarchy, affordance）
- **（Rule 14）** NN/g — Progressive Disclosure + Cognitive Load heuristics
- **（Rule 14）** GOV.UK Design System / Material Design / Apple HIG（interaction affordance, responsive）
