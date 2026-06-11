<!--
建立時間: 2026-05-20
agent: main session (Task #4 skill workflow audit)
status: draft (awaiting user review + fix execution)
report_class: skill-workflow-audit
audience: 廖子游 (self) · skill maintainer
scope: 2026-05-20 weekly PI deck session (Version A through H, 8 versions iteration)
trigger: Task #4 pending from PPT iteration session
verdict: 10 issues identified · 4 P0 · 4 P1 · 2 P2 · estimated 4-6 hr total fix
last_verified: 2026-05-20
-->

# PPT Workflow Skill Audit — 2026-05-20 Session Retrospective

> **Scope**: 2026-05-20 weekly PI deck session 暴露的 `myPPT → weekly-report → html-report-build → pptx-build` 工作流問題。本次 session 從 myPPT 場景識別開始，產出 8 個版本 HTML (A-H)，過程中暴露 10 個系統性 issue 值得內建到 skill。

> **核心發現**：用戶**完全 skip weekly-report**（無 master_draft.md）直接從 myPPT 跳 html-report-build；場景識別未問「PI 1-on-1 vs Lab meeting」導致先 G (lab) 才補 H (1-on-1)；critical data fabrication (+0.057) 跨 4 版本傳遞未自審；figure 路徑 off-by-one 跨 4 版本未偵測。

---

## 0. TL;DR — 10 個 issue × 4 skill files

| # | Issue | Severity | 影響 skill | Fix effort |
|---|------|:--------:|------------|:----------:|
| **1** | Figure 相對路徑 5 vs 6 levels up off-by-one bug | 🔴 P0 | html-report-build | 30 min |
| **2** | 數字 fabrication 跨版本傳遞（T=3 +0.057） | 🔴 P0 | html-report-build, weekly-report | 1 hr |
| **3** | 場景識別缺「報告場景」（PI 1-on-1 / Lab meeting / Conference） | 🔴 P0 | myPPT | 20 min |
| **4** | PI 適合量化標準未內建（≤150 / ≤250 / ≤350 chars 按場景） | 🔴 P0 | html-report-build, pptx-build | 45 min |
| **5** | Speaker note `<details>` 折疊未預設機制 | 🟠 P1 | html-report-build | 30 min |
| **6** | 多版差異化（A-D 並排）缺 differentiation guidance | 🟠 P1 | html-report-build | 30 min |
| **7** | Self-audit number-source-grep 強制機制缺 | 🟠 P1 | html-report-build, structured-tech-report | 1 hr |
| **8** | weekly-report skip 路徑未明確（LLM-direct fallback） | 🟠 P1 | myPPT, weekly-report | 30 min |
| **9** | 路徑前綴 `InterSubMod/...` 規則未 self-audit | 🟡 P2 | html-report-build, weekly-report | 15 min |
| **10** | Evaluator polish notes 後續處理 loop 未串接 | 🟡 P2 | html-report-build | 20 min |

**Total**: ~ 4.5-6 hr fix。建議 P0 4 條優先（≈ 2.5 hr），P1 4 條後續（≈ 2.5 hr），P2 2 條隨手做。

---

## 1. Critical issues (P0 · 必修)

### Issue 1: Figure 相對路徑 off-by-one bug

**Symptom**: 本次 session A-D 4 個版本 HTML 全部圖片 broken
- 我用 `../../../../../research/...`（5 levels up）
- 正確應為 `../../../../../../research/...`（6 levels up）
- HTML 在 `docs/presentations/in_progress/2026/05/weekly_PI_20260520_4versions/` (6 levels deep)，回到 root 需 6 `../`

**Root cause**: LLM 計算 relative path 時 off-by-one；無 build step 驗證。

**Affected files**:
- `InterSubMod/.claude/skills/html-report-build/SKILL.md` § Failure Modes
- `InterSubMod/.claude/skills/html-report-build/prompts/slide_prompt.md`

**Fix proposal**:
1. **加 path verification step to slide_prompt.md**：
   ```
   Before Write: for each <img src="..">, count "../" tokens.
   Expected = depth_of(target_file) - depth_of(source_dir) + 1
   If mismatch → recompute or escalate
   ```
2. **SKILL.md Failure Modes table 加新 row**：
   ```
   | Symptom | Cause | Fix |
   |---------|-------|-----|
   | <img> broken in browser | Wrong number of ../ in relative path | grep + count depth; LLM 易 off-by-one |
   ```
3. **Post-Write hook** (optional): bash 腳本 grep all `src="..."` paths, verify file existence with `ls`.

**Effort**: 30 min（SKILL.md edit + prompt edit）

---

### Issue 2: 數字 fabrication 跨版本傳遞

**Symptom**:
- Session 初期我寫 Version C V6 sign-off matrix 含 "ΔF1 +0.030 (T=2) · +0.057 (T=3) strict dominance"
- 真實 T=3 ΔF1 = **+0.1605**（source `InterSubMod/docs/experiments/in_progress/2026/05/20260519_V6_vs_baseline_HCC1395_TPFP_comparison_01.md` §11.3 line 531）
- C 還含整個 fabricated F1 sweep table（baseline T=3=0.6541 應為 0.4308 / V6 T=3=0.7106 應為 0.5913）
- 錯誤跨 C → D → E 三版未自審，evaluator 第二輪審查才揪出

**Root cause**:
- LLM 生成大量 deck 時內插「合理範圍」數字而非嚴格引用 source
- Self-audit 缺強制 number-source-grep
- 多版差異化迭代複製錯誤

**Affected files**:
- `InterSubMod/.claude/skills/html-report-build/SKILL.md` § 嚴謹度繼承
- `InterSubMod/.claude/skills/html-report-build/prompts/slide_prompt.md`
- `InterSubMod/.claude/skills/scientific-rigor/SKILL.md` § §8.4 provenance audit

**Fix proposal**:
1. **強制 number-source-grep audit step**（before Write 每個 slide）：
   ```
   For each numerical value mentioned in slide:
     - Tag with [source: <md_path>:line] inline comment
     - grep source .md to confirm value exists verbatim or rounded ≤ 2 sf
     - If not found → FAIL slide, ask user OR use only ranges
   ```
2. **多版迭代時 source citation 鎖**：
   - 第一版產出後鎖 source set（list of .md + line ranges）
   - 後續版本只能 reuse 同 source set 的數字，不可新增
3. **6-Taboo Audit 加第 7 條**：「number without source citation → 紅旗」

**Effort**: 1 hr（slide_prompt 大改 + scientific-rigor §8.4 補強）

---

### Issue 3: 場景識別缺「報告場景」

**Symptom**:
- 本次 myPPT 場景識別只問「週報 / 已有母稿 / 從零 / 完整 pipeline」4 option
- 沒問「PI 1-on-1（5-10 min） vs Lab meeting（15-30 min） vs Conference（45-60 min）」
- 結果先產 G (lab meeting 6 slide 偏 dense) → 用戶才指出要 PI 1-on-1 → 補產 H

**Root cause**: myPPT 場景識別維度僅「母稿 vs 從零」單軸；缺「報告時長/受眾」第二軸。

**Affected files**:
- `InterSubMod/.claude/skills/myPPT/SKILL.md` § 場景識別 4 options
- `InterSubMod/.claude/skills/pptx-build/SKILL.md` § P1 audience tier

**Fix proposal**:
1. **myPPT 場景識別改成 2-axis matrix**：
   ```
   Axis 1: 來源 (週報母稿 / 已有母稿 / 從零)
   Axis 2: 場景 (PI 1-on-1 / Lab meeting / Conference / Lab informal)
   ```
   AskUserQuestion 拆兩個 question (each 4 options)
2. **每場景對應預設參數**：
   ```yaml
   PI 1-on-1 (5-10 min):
     slide_count: 3-6
     chars_per_slide_max: 150 純文字 / 250 含表
     table_rows_max: 5
     hero_visual_ratio: 60-75%
     reading_test: 3-second

   Lab meeting (15-30 min):
     slide_count: 8-15
     chars_per_slide_max: 250 純文字 / 350 含表
     table_rows_max: 8
     hero_visual_ratio: 50-65%
     reading_test: 5-second

   Conference (45-60 min):
     slide_count: 20-40
     chars_per_slide_max: 350 純文字 / 450 含表
     table_rows_max: 12
     hero_visual_ratio: 40-60%
     reading_test: 10-second
   ```

**Effort**: 20 min（myPPT SKILL.md + pptx-build P1 補）

---

### Issue 4: PI 適合量化標準未內建

**Symptom**:
- A cover slide 過密 (V6 + 2 pillar cards + statement 一頁)
- D Version content-adaptive 過密
- G S5 「5 目標 + 還能驗證什麼」雙 table 16 rows 過密
- 用戶反覆指出「資訊太多」「字太雜」「資訊量不夠精簡」

**Root cause**: skill 內無 quantitative 標準對照表；LLM 依直覺生成。

**Affected files**:
- `InterSubMod/.claude/skills/html-report-build/references/design_principles.md`
- `InterSubMod/.claude/skills/pptx-build/SKILL.md`

**Fix proposal**:
1. **design_principles.md 加新 §13 場景量化標準表**（如 Issue 3 提案的 YAML 表）
2. **slide_prompt.md / pptx-build P4 (slide build) 加 per-slide audit**：
   ```
   After each slide draft:
     - count chars (excl speaker note)
     - count table rows
     - if exceed scenario threshold → flag for compression
   ```
3. **6-Taboo Audit 加第 8 條**：「scenario chars/rows exceeded → 紅旗」

**Effort**: 45 min（design_principles 加 §13 + 兩 skill prompts 補 audit step）

---

## 2. Important issues (P1 · 建議修)

### Issue 5: Speaker note `<details>` 折疊未預設機制

**Symptom**:
- G / H 才採用「`<details>` 折疊 speaker note」設計
- A-F 沒有，導致 slide 上 narrative + caveat + provenance 全擠進可見區
- 用戶最後才明示「詳細放講稿」

**Root cause**: html-report-build slide mode 預設不含 `<details>` block；`components/speaker_note.html` 是 visible 不是 collapsed。

**Affected files**:
- `InterSubMod/.claude/skills/html-report-build/components/speaker_note.html`
- `InterSubMod/.claude/skills/html-report-build/prompts/slide_prompt.md`
- 新建 `components/speaker_note_details.html`

**Fix proposal**:
1. 新建 `components/speaker_note_details.html` 含 `<details>` 折疊 markup + CSS
2. slide_prompt.md Rule 加：「PI 1-on-1 + Lab meeting 場景預設用 details；Conference 場景用 visible note」
3. 連結 Issue 3 場景識別：場景決定哪種 speaker note 機制

**Effort**: 30 min（新 component + prompt update）

---

### Issue 6: 多版差異化缺 differentiation guidance

**Symptom**:
- A-D 4 版本「並排」由用戶要求生成，但很多 slide 內容/結構重複
- C/D/E 都有同樣 fabricated +0.057 issue（複製錯誤）
- 沒有明確的「每版 differentiation axis」guidance

**Root cause**: html-report-build 沒有 multi-version output 機制；LLM 自由發揮導致 quality variance。

**Affected files**:
- `InterSubMod/.claude/skills/html-report-build/SKILL.md` § Quick Usage
- 新建 `prompts/multi_version_differentiation.md`

**Fix proposal**:
1. 新建 `prompts/multi_version_differentiation.md` 含：
   ```
   When user requests N versions for comparison:
     1. Define N differentiation axes (e.g. minimal / standard / dense / hybrid)
     2. Each version pick ONE axis, NOT N=N independent designs
     3. Common source-of-truth: shared metadata block at top
     4. After all versions: 1-page comparison index.html
   ```
2. SKILL.md § Quick Usage 加「Multi-version output」段落

**Effort**: 30 min（新 prompt + SKILL.md 補一段）

---

### Issue 7: Self-audit number-source-grep 強制機制缺

**Symptom**: 同 Issue 2 root cause 的另一面 — 沒有 self-audit gate。Evaluator 第二輪才揪出 +0.057。

**Root cause**: scientific-rigor SKILL.md §8.4 provenance 對 .md 報告 enforced，但對 .html slide 沒延伸。

**Affected files**:
- `InterSubMod/.claude/skills/scientific-rigor/SKILL.md` §8.4
- `InterSubMod/.claude/skills/html-report-build/SKILL.md` § 嚴謹度繼承

**Fix proposal**:
1. scientific-rigor §8.4 改寫：「**任何 derivative artifact（HTML / PPTX / email draft）must inherit provenance audit gate**」
2. html-report-build § 嚴謹度繼承 加：「Slide mode 數字 ≥ 3 個 → 必跑 number-source-grep self-audit」

**Effort**: 1 hr（兩 skill 嚴謹度章節對齊 + audit step 落地）

---

### Issue 8: weekly-report skip 路徑未明確

**Symptom**:
- 本次 session 用戶說「整理確認這週的重點」→ myPPT 場景識別 → 直接跳 html-report-build
- 完全 skip weekly-report W1-W7 母稿生成
- 沒有 master_draft.md，所有 narrative 直接 LLM-direct 內插
- 後續多版迭代沒 source-of-truth 母稿

**Root cause**: myPPT 場景識別當用戶說「整理」時應推 weekly-report 路徑；但 LLM 走捷徑直接到 html-report-build slide mode。

**Affected files**:
- `InterSubMod/.claude/skills/myPPT/SKILL.md` § 委派邏輯
- `InterSubMod/.claude/skills/weekly-report/SKILL.md` § 觸發

**Fix proposal**:
1. myPPT 場景識別觸發語明確化：
   ```
   「整理這週的重點」「給 PI 看的進度」「研究進度」「週報」「lab meeting」
     → 必經 weekly-report W1-W7（不可 skip）
   「已有母稿」「直接產 PPT」「我要做簡報」
     → 走 pptx-build / html-report-build slide mode
   ```
2. weekly-report SKILL.md 加 「**Skip 條件**」：「除非用戶明示『不需要母稿』，否則本 skill 不可被略過」
3. 連結 Issue 6：master_draft.md 為多版差異化的 source-of-truth

**Effort**: 30 min（兩 skill SKILL.md 補章節）

---

## 3. Minor issues (P2 · 隨手做)

### Issue 9: 路徑前綴 `InterSubMod/...` 規則未 self-audit

**Symptom**: 早期版本 footer / source 引用有時用 `docs/...`、有時 `/big7_disk/.../InterSubMod/docs/...`，未統一。

**Root cause**: UserPromptSubmit hook 已 enforce 但 LLM 在生成 HTML 內容時未自審。

**Fix proposal**: slide_prompt.md 加：「**all .md path references in footer/source must start with `InterSubMod/`**」

**Effort**: 15 min

---

### Issue 10: Evaluator polish notes 後續處理 loop 未串接

**Symptom**: Evaluator audit 出 3 個 polish notes (S5 title 數字、S4 noise floor 註、S1 hedge)，用戶選 ship as-is，但下次重做同類 deck 時這些 polish lessons 不會自動套用。

**Root cause**: evaluator 輸出未持久化進 known-pitfalls 或 skill 改進。

**Fix proposal**: evaluator agent prompt 加：「**若 polish notes 被用戶選 ship as-is，將 note 追加到 known-pitfalls .md 與 skill MEMORY.md，作為下次預防**」

**Effort**: 20 min

---

## 4. Recommended action plan

### Phase 1 · P0 quick wins（2.5 hr 可一次做完）

執行順序：
1. **Issue 3 myPPT 場景識別 2-axis 改造**（20 min）→ 阻擋未來「先做完才發現場景錯」
2. **Issue 4 PI 適合量化標準 design_principles §13**（45 min）→ 場景對應字數/row 表內建
3. **Issue 1 path off-by-one prompt update**（30 min）→ 阻擋 figure broken
4. **Issue 2 number-source-grep self-audit**（1 hr）→ 阻擋 fabrication

### Phase 2 · P1 後續（2.5 hr）

5. **Issue 5 speaker note details component**（30 min）
6. **Issue 6 multi-version differentiation prompt**（30 min）
7. **Issue 7 嚴謹度跨 skill 對齊**（1 hr）
8. **Issue 8 weekly-report skip 路徑明確化**（30 min）

### Phase 3 · P2 polish（35 min）

9. **Issue 9 path prefix self-audit**（15 min）
10. **Issue 10 evaluator polish loop**（20 min）

---

## 5. Top 3 quick wins（若只能做 3 條）

| 排序 | Issue | 為何優先 | Effort |
|:-:|------|---------|:------:|
| 1 | **#3 場景識別 2-axis** | 阻擋「先做完才補另一版本」 = 直接省下 1-2 hr 重做 | 20 min |
| 2 | **#1 path off-by-one** | 每次 figure broken = 用戶看不到圖 = 重做 | 30 min |
| 3 | **#2 number-source-grep** | fabrication 影響 trust 嚴重，evaluator catch 在第 2 輪太晚 | 1 hr |

**Total Top 3**: ~ 1 hr 50 min，可阻擋本次 session 重複出現的 80% 問題。

---

## 6. 觀察 metadata

- **Session 長度**: 約 4-5 hr (artist 估)
- **產出版本數**: 8 個 HTML (A/B/C/D/E v1/E v2/F/G/H) + 1 index
- **Evaluator 跑次數**: 2 (E v2 / G)
- **發現 critical bug 階段**: evaluator 第 1 輪 catch sign-off email range 但不 specific T=3；evaluator 第 2 輪 + 主 agent grep catch +0.057 fabrication
- **重做率**: A → E 為 differentiation iteration（非重做），E v2 → F → G → H 為 narrative pivot 重做 ≈ 4 次

---

## 7. Cross-reference

### 相關 known-pitfalls (應更新)
- `InterSubMod/.claude/skills/known-pitfalls/SKILL.md` → 應加新 entry：「number fabrication in HTML slide derived from .md report」
- `InterSubMod/.claude/skills/known-pitfalls/SKILL.md` → 應加：「relative path off-by-one in HTML output」

### 相關 memory entries (應更新)
- `feedback-ppt-minimal-visual-first.md` ✅ 本 session 已建立
- 應加：`reference-pi-scenario-quantitative-standards.md`（PI 1-on-1 / Lab meeting / Conference 量化標準表）

### 相關 hook 機會
- Post-Write hook on `*.html` in `docs/presentations/`: grep `<img src=>` paths + verify file exists
- Post-Write hook on `*.html`: grep numerical values + sample compare to source .md

---

## 8. Decision points for 廖子游

1. **是否執行 Phase 1 P0 4 條 quick fix？**（2.5 hr，阻擋 80% 重複問題）
2. **是否優先做 Top 3 only？**（~ 1 hr 50 min）
3. **memory 是否加 PI scenario quantitative standards reference**？
4. **是否將本 audit report 與 audit lessons 加進 known-pitfalls？**

---

**End of PPT Workflow Skill Audit**

Generated 2026-05-20 by main session Task #4. Sources:
- Session artifacts: `InterSubMod/docs/presentations/in_progress/2026/05/weekly_PI_20260520_4versions/` (8 HTML versions + index)
- Skill files audited: `InterSubMod/.claude/skills/{myPPT,weekly-report,html-report-build,pptx-build}/SKILL.md`
- Evaluator audit results: 2 rounds (E v2 PASS L1 + G PASS L1)
- Memory: `feedback-ppt-minimal-visual-first.md` (2026-05-20 created)
