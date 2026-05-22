---
name: narrative-frame
description: |
  全域敘述框架挑選 + 套用 + 自審 — 50+ 業界 framework catalog (SCQA / Pyramid / MECE /
  BLUF / STAR / CAR / PREP / OREO / ABT / Pixar Spine / 三幕劇 / Hero's Journey /
  Duarte Sparkline / Golden Circle / Monroe / Assertion-Evidence / 5W1H / 5 Whys /
  Fishbone / A3 / DMAIC / PDCA / OODA / Diátaxis / Feynman / Bloom / KWL / ELI5 /
  SBI / DESC / GROW / Pendleton / AIDA / FAB / PAS / BAB / WRAP / Pre-Mortem /
  Cynefin / ADR / DACI / DECIDE / TED / 3-Act Pitch / Resonate / PechaKucha / PEEL...).
  6-step workflow: N1 場景識別 (5W) → N2 框架推薦 (主 + 1-2 備 + 解釋為什麼) →
  N3 source 萃取 (多文件 [F/O/I/U] tier 標記) → N4 mapping 表 (重點 ↔ section) →
  N5 套用 + 補齊 (產出 structured narrative) → N6 缺漏標 + 5 秒測試 + 用戶 ack.
  取代既有 7 報告類 skill (structured-tech-report / weekly-report / pptx-build /
  results-report / conclude-research / report / myPPT) 固定範本之上層動態挑選機制.
  與 /pre-decision-audit (decision 層) 正交 — 本 skill 是 expression/communication 層.
  USE WHEN: 用戶説「整理 / 報告 / 説明 / 彙報 / 總結 / pitch / 答辯 / 簡報 / 教 / 解釋 /
  寫 / 整合」、輸出 ≥200 字且跨 ≥2 概念、多文件統整、複雜分析回覆、AI 回覆減少理解負擔.
  SKIP WHEN: 純 code edit / typo / build / commit / single-line answer / factual lookup /
  簡單問答 / 已在 active framework 中.
allowed-tools: Read, Write, Edit, Glob, Grep, Bash, Agent, AskUserQuestion
user-invocable: true
paths:
  - "**/*.md"
  - "docs/reports/**/*.md"
  - "docs/presentations/**/*.md"
  - "research/**/*.md"
tags: ["narrative", "framework", "catalog", "expression", "communication", "scqa", "pyramid", "star", "pixar", "abt", "golden-circle", "a3"]
---

# narrative-frame — 全域敘述框架庫

> **取代既有 7 報告類 skill 固定範本之上層動態挑選機制**。50+ 業界 framework catalog；場景自適配；對話層級啟用；減少用戶理解負擔。

> **路徑硬性規則**：本 skill 所有 .md 路徑列給用戶時必以 `InterSubMod/...` 前綴。

---

## §1 Phase & Chain Position

**Cross-cutting**：覆蓋所有「整理 / 報告 / 説明」AI 回覆 — 跨 P0-P6 phase、跨對話模式、跨輸出形式（口頭 / .md / slide / HTML / paper）。

### 啟用層級（3 Tier — 對應 AGENTS.md §15 框架聲明維度）

| Tier | 觸發條件 | 行為 |
|------|---------|------|
| **Tier 1** 簡單問答 | factual lookup / single-line / yes-no | **skip** framework（無 overhead） |
| **Tier 2** 中度整理 | 200-500 字 / 跨 2-3 概念 | 回覆**首行**聲明 framework（如「用 PREP：」），結構化但不跑完整 N1-N6 |
| **Tier 3** 複雜回覆 | ≥500 字 / 跨 ≥3 概念 / 多文件 / 結論性報告 | 完整跑 **N1-N6** + structured output + source mapping |

### 與「假説驗證三層樓」關係

- **decision / hypothesis** 層：`/pre-decision-audit`（pre）→ `/implementation-notes`（process）→ `/run-evaluator`（post）
- **expression / communication** 層：`/narrative-frame`（本 skill）
- **正交不衝突**：可同時用（pre-decision-audit 決定要不要做 + narrative-frame 決定怎麼講）

---

## §2 Dependencies

### Uses（本 skill 引用）

- `/scientific-rigor` §2-§7 — L1-L5 evidence tier 借用為本 skill 的 [F]/[O]/[I]/[U] 對應
- `/pre-decision-audit` §0 Cynefin — 場景分類 front-gate（Complex 域強制 PROBE，禁套 best-practice）
- `/research-context-loader` Tier 1/2/3 — N3 source 載入
- `/doc-standards` — filename + tier ribbon 規範
- `/known-pitfalls` P-01-P-14 — N6 自審 cross-check
- `references/framework_catalog.md` — 50+ framework SoT（單一 source of truth）
- `references/scenario_to_framework.md` — 5W → 推薦框架對照
- `references/framework_business_sources.md` — URL + ISBN + 一句引用
- `references/design_principles.md` — 5 秒測試 + 6-item pre-publish（借用 html-report-build）

### Used by（誰呼叫本 skill）

- `/structured-tech-report`（thin wrapper）— 預設 `apply A3+ADR+Postmortem-hybrid`（即原 13 段）
- `/weekly-report`（thin wrapper）— 預設 `apply Multi-Thread-Narrative`（即原 17 段 4 主線）
- `/pptx-build`（thin wrapper）— 場景識別後 forward
- `/results-report`（thin wrapper）— 預設 `apply Data-Showcase`
- `/conclude-research`（thin wrapper）— 預設 `apply Verdict-Pyramid`
- `/report`（thin wrapper）— 預設 `apply AI-Session-Companion`
- `/myPPT`（thin wrapper）— 場景路由總入口
- **AI 對話層**：UserPromptSubmit hook `narrative_frame_advisor.sh` 偵測 keyword 推薦觸發
- `/html-report-build` standalone mode — 可在 L0/L1 引用 framework summary

### Reads

- 用戶 raw request（當輪 prompt）
- N 份 source `.md` / `.tsv` / `.csv`（用戶提供或專案內 grep）
- `MEMORY.md` — 過往 framework usage / preference
- `InterSubMod/docs/CURRENT_FOCUS.md` — live context
- Active cycle state: `InterSubMod/state/cycles/<cycle_id>/`

### Writes

- Inline 對話回覆（首行 framework 聲明 + 結構化內容）
- `.md` 報告（含 frontmatter `framework:` + section 對應 source citation）
- 缺漏 ⚠ gap 標記
- （可選）`state/narrative_frame/<topic>_framework_choice.json` — 後續 audit 用

---

## §3 6-step Workflow + 3 Checkpoint

```
[Input: 用戶 raw request + N 份 source docs]
    ↓
N1. 場景識別（5W: Who / Why / What / When / How）
    ↓ C1 確認場景分類 ★ 必停（互動模式）/ fast-track 30s 倒數（全自動）
N2. 框架推薦（從 catalog 挑主 + 1-2 備 + 解釋為什麼選 + 為什麼不選其他）
    ↓ C2 確認 framework ★ 必停（用戶可 override 換框架）
N3. Source 萃取（呼叫 narrative-organizer agent if N≥3 docs）
    每筆素材標 [F]/[O]/[I]/[U] tier + source file:line
N4. Mapping 表（重點 ↔ 框架 section 對應）
N5. 套用 + 補齊（產出 structured narrative）
    自審：每 section 有對應 source？無 → 標 ⚠ gap
N6. Gap 標 + 5 秒測試 + 6-item pre-publish + 用戶確認
    ↓ C3 ack ★ 必停
[Output: framework 標題 + 結構化 narrative + source mapping + gap 標]
```

### 執行模式（與 confirmation-protocol 對齊）

- **互動模式**（預設）：C1 / C2 / C3 暫停等用戶確認
- **全自動**（`auto` / `全自動`）：保 **C2 framework 鎖定**必停；C1 / C3 用 AI 預設標籤 30s 倒數通過

---

## §4 N1 — 5W 場景識別（C1 確認）

對用戶 raw request 連問 5 維：

| 維度 | 選項 | 影響 |
|------|------|------|
| **Who** 受眾 | PI / 教授 / 同儕 / 大眾 / 自己 / 學生 / 投資人 / 評審 | 影響 framework formality + 術語層 |
| **Why** 目的 | 説服 / 解釋 / 報告進度 / 探索 / 紀錄 / 答辯 / 教學 / 比較 / 決策 | 影響主框架類別 §1-§10 |
| **What** 內容類型 | 新發現 / 進度 / 概念 / 比較 / 故事 / 決策 / 案例 / 數據 / 流程 | 影響 framework specifics |
| **When** 時長 | 30s（pitch）/ 5min（lab） / 18min（TED）/ 60min（lecture）/ 紙本（無限） | 影響細節粒度 + 段數 |
| **How** 輸出形式 | 口頭 / slide / .md / HTML / paper / email / Slack / 對話 | 影響 framework rendering |

### 場景特例（直接路由 — 不需 5 維全跑）

| 觸發語 | 路由 |
|--------|------|
| 「整理本週 / 週報 / lab meeting」 | weekly-report wrapper → `apply Multi-Thread-Narrative` |
| 「整理 commit / fix / 修補」 | structured-tech-report wrapper → `apply A3+ADR+Postmortem-hybrid` |
| 「整理實驗結果」 | results-report wrapper → `apply Data-Showcase` |
| 「整理 PI 簡報 / pitch」 | pptx-build wrapper → `apply Audience-Scenario + framework` |
| 「整理 AI session」 | report wrapper → `apply AI-Session-Companion` |
| 「整理研究收尾 / 假説 verdict」 | conclude-research wrapper → `apply Verdict-Pyramid` |

---

## §5 N2 — 框架推薦邏輯（C2 確認）

讀 `references/scenario_to_framework.md` 對照表，挑 **主框架 + 1-2 備框架**，輸出格式：

```
[narrative-frame N2]
主框架: <name>（<業界源>）
理由: <why this framework fits 場景 5W>
備框架 1: <name>（<差異點 — 何時切換到這個更好>）
備框架 2: <name>（同上）

不選擇:
- <name>: <為什麼不適合>
```

### 框架推薦規則（heuristic + Cynefin）

| 場景特徵 | 推薦類 | 例 |
|---------|--------|----|
| Who=PI/教授 + Why=報告進度 + When=5min | §1 結論先行 + §2 經驗敘事 | SCQA + STAR per item |
| Who=同儕 + Why=説服 + When=18min | §4 故事敘事 | Pixar Spine / Duarte Sparkline + ABT |
| Who=自己 + Why=分析 + What=決策 | §9 決策評估 + §5 問題分析 | WRAP + Pre-Mortem |
| Why=答辯 + When=30s | §3 短時口頭 | PREP / ABT |
| Why=教學 + What=概念 | §6 教學解釋 | Feynman / ELI5 / Diátaxis |
| Why=回饋 + Who=同事 | §7 回饋對話 | SBI / Radical Candor |
| What=bug postmortem | §5 問題分析 | A3 + 5 Whys + Inverted Pyramid |
| Why=説服募款 | §8 行銷説服 + §4 故事敘事 | Monroe + Hero's Journey |
| How=slide deck（任何） | §10 簡報架構 為外殼 + §1-§9 為內裏 | Assertion-Evidence + (主題 framework) |

### Cynefin Domain Gate（front-gate）

- **Clear**: 可直套 best-practice framework（SCQA / STAR / PREP）
- **Complicated**: 框架可選，但建議 §5 problem analysis（5 Whys / Fishbone）
- **Complex**: **強制 probe-first**（不可套 deterministic framework）；改用 §6 教學（Feynman 反問）或 §9 WRAP（reality-test）
- **Chaotic**: 先穩定再敘述；不跑 framework

借用 `/pre-decision-audit §0` 邏輯，不重複定義。

---

## §6 N3-N5 — 萃取、Mapping、套用

### N3 Source 萃取

| Source 數量 | 處理方式 |
|------------|---------|
| 1 份 | 直接 inline Read + grep 重點 |
| 2 份 | 順序 Read + 標 [F/O/I/U] tier |
| **≥3 份** | **啟動 `narrative-organizer` agent 並行萃取** + cross-file 主題聚類 |

每筆素材必標：
- **Tier**：[F]act / [O]bservation / [I]nference / [U]nconfirmed（借用 weekly-report W3 4 層分類）
- **Source**：`InterSubMod/...:line` 完整 path
- **Confidence**：L1-L5 ribbon（借用 /scientific-rigor §2）

### N4 Mapping 表

```markdown
| Framework section | 重點 | Source | Tier |
|-------------------|------|--------|------|
| <e.g., SCQA Situation> | <key fact> | `path:line` | [F] L3 |
| <SCQA Complication> | ... | ... | ... |
| <SCQA Question> | ... | ... | ... |
| <SCQA Answer> | ... | ... | ... |
```

**自審強制**：每個 framework section **必有 ≥1 對應重點**；無 → 標 ⚠ gap 進 N6。

### N5 套用 + 補齊

讀 `templates/<framework>_skeleton.md` → 填入 N4 mapping 內容 → 輸出 structured narrative。

**輸出 mode**（用戶可選）：
- `--mode=outline` — markdown outline
- `--mode=script` — 口頭講稿（含 timing 標記）
- `--mode=slide` — slide outline（接 pptx-build）
- `--mode=companion` — .md 報告附加 framework summary section
- `--mode=inline`（預設 Tier 2）— 對話 inline 回覆首行聲明 framework

---

## §7 N6 — 缺漏標 + 自審 + 用戶確認（C3）

### 自審 checklist（執行前必跑）

1. **每 framework section 都有對應 source？** 無 → ⚠ gap
2. **重要素材沒漏？** 跑 N3 萃取 list vs N4 mapping 對照
3. **5 秒測試**（NN/g + Berkeley dataviz）：給同事看 5 秒能否説出 takeaway？
4. **3 秒法則**（Nancy Duarte）：標題 + 主視覺 3 秒內傳達主訊息？
5. **Assertion-Evidence**：每段標題 = 結論句（不只是 topic 標籤）？
6. **6-item pre-publish**（design_principles.md Rule 12）：data-ink / CRAP / 留白 / hierarchy / 色 / colorblind 全 ✓？
7. **過度宣稱 check**（borrow weekly-report W5）：[O]/[I] 用了「證實 / 確認 / 解決」→ 改謙詞
8. **路徑前綴**：所有 .md ref 用 `InterSubMod/...`？

### 缺漏標格式

```markdown
## ⚠ Gap（N6 自審）

- **SCQA Complication 段缺 source**: 推測但無 evidence；補 file:line 或降為 [I]
- **stat 5 秒測試 fail**: 首段太密；建議拆 stat-grid 4 個焦點數字
- **過度宣稱**: §X 用「證實」但 tier=[O]；改「初步觀察到」
```

### C3 用戶確認形式

```
[narrative-frame C3 確認]
1. 套用 <framework> 完成；產出 X 字 / Y sections
2. ⚠ Gap: N 項（見上）
3. 5 秒測試: PASS / FAIL（理由）
4. 是否接受？(y/n/edit/換框架)
```

---

## §8 3 sub-commands

### `pick <scenario>`

純框架挑選，不套用 — 對既有 .md / 草稿 / 對話脈絡跑 N1-N2，輸出推薦結果即停。

**用例**：用戶有報告草稿但不確定 framework；用 pick 拿 3 候選對比。

### `apply <framework> <source>`

跳過 N1-N2，直接套用指定 framework — 給 advanced user 已知要用 SCQA / Pyramid / Pixar Spine 等時加速。

**用例**：thin wrapper skill 內部呼叫（如 `/structured-tech-report` → `apply A3+ADR+Postmortem-hybrid <commit>`）。

### `audit <output>`

對已產出回覆 / .md 做框架審查 — 跑 N6 自審 checklist 標 gap，不重產出。

**用例**：既有 .md 報告 retrofit 檢查；PR 前 final check。

---

## §9 Framework Catalog 速覽（10 大類）

> **完整定義 50+ entries 見 `references/framework_catalog.md`**。本表為入口索引。

| 類 | 主要 framework | 適合場景 |
|----|-------------|---------|
| **§1 結論先行 / 商務報告** | SCQA / Pyramid / MECE / BLUF / Inverted Pyramid | exec summary, PI report, news lead |
| **§2 經驗敘事 / 行為案例** | STAR / CAR / PAR / SOAR / STARR | 面試, 履歷, performance review |
| **§3 短時口頭答辯** | PREP / OREO / ABT / Elevator Pitch / Rule of Three | Q&A, 答辯, 30s pitch |
| **§4 故事敘事** | Pixar Story Spine / 三幕劇 / 英雄之旅 / Freytag / Duarte Sparkline / Golden Circle / Monroe / Assertion-Evidence | TED, 科普, keynote, paper abstract |
| **§5 問題分析因果** | 5W1H / 5 Whys / Fishbone / A3 / DMAIC / PDCA / OODA | postmortem, RCA, 改善 |
| **§6 教學解釋** | Feynman / Diátaxis / Bloom / KWL / ELI5 / CPP | tutorial, doc, 教案 |
| **§7 回饋對話 Coaching** | SBI / DESC / GROW / Pendleton / Sandwich / Radical Candor | 1-on-1, mentoring |
| **§8 行銷説服** | AIDA / AIDCA / FAB / PAS / 4P / BAB | landing page, copywriting |
| **§9 決策評估** | WRAP / Pre-Mortem / Cynefin / ADR / DACI / DECIDE / Eisenhower | 決策, risk audit |
| **§10 簡報架構** | TED / 3-Act / Resonate / Presentation Zen / Lessig / Takahashi / PechaKucha / PEEL | slide deck 外殼 |

業界源頭：Barbara Minto / Randy Olson / Joseph Campbell / Syd Field / Nancy Duarte / Simon Sinek / Michael Alley / Toyota / Six Sigma / John Boyd / Daniele Procida / Bloom / CCL / Heath / Klein / Snowden / Kenn Adams / Alan Monroe — 詳見 `references/framework_business_sources.md`。

---

## §10 Failure Modes & Diagnostics

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| 框架選錯 用戶體感不適 | N2 推薦 heuristic 過簡 | 強制 N2 出 1-2 備選 + C2 用戶 ack；用戶説「換 X」立刻跳 N5 重套 |
| Section 缺 source | N4 mapping 漏 | N6 自審強制 ⚠ gap；用戶可選補 source 或降 tier |
| 過度宣稱（[O]/[I] 用「證實」）| AI 預設高 confidence | N6 過度宣稱 check（borrow weekly-report W5 紅旗） |
| Hook noisy 每次推薦 framework | UserPromptSubmit 過敏感 | 200 字 / 2 概念門檻；用戶可任何時 `不用框架` override |
| Catalog 認知過載（50+ 太多） | 用戶不熟 catalog | `pick` 模式只挑 3 候選；scenario_to_framework 快查不必讀全 catalog |
| 與 /pre-decision-audit 混淆 | 兩者角色重疊 | pre-decision-audit = decision 層；narrative-frame = expression 層；可同時用 |
| Multi-doc 萃取慢 (>15min) | source ≥10 份 / 大檔 ≥1000 行 | narrative-organizer agent 預估 timeout；建議用戶縮 source 或分批 |
| 用戶不確定 framework 名 | 中英文 alias 不足 | scenario_to_framework 用「想做 X」自然語對照（不必背 framework 名） |
| Thin wrapper 既有 skill regression | 預設 framework entry 對應錯 | V3 smoke test：跑既有 skill → 對比舊輸出結構（不可 regression） |
| Tier 2 inline 首行聲明過長 | 用戶看到 framework 名陌生 | 中文括號註：「用 SCQA（McKinsey 結論先行）：」 |

---

## §11 嚴謹度繼承（/scientific-rigor）

任何套用 framework 的結論性回覆（Tier 3）必繼承 `/scientific-rigor §2-§7`：

- **§2 證據分級**：N3 萃取每筆素材必標 L1-L5 tier（[F]=L1 fully verified / [O]=L2-L3 / [I]=L4 / [U]=L5）
- **§3 Effect Size**：含數字 claim 必配 Cohen ribbon + CI；不只列點數
- **§4 DAG**：因果敘事段必標 confounder / collider / mediator（特別在 §5 問題分析類 framework）
- **§7 Pre-registration**：N2 框架選擇 + N5 預期 narrative section 寫入 audit log，避免事後 HARKing
- **§8.4 Provenance**：每 framework 套用後輸出含 commit hash + cycle_id + 生成時間 footer

**最小可用子集**:
- Tier 3 結論性報告 framework: §2 + §3 + §4 + §7 + §8.4 全跑
- Tier 2 中度整理 framework: §2 + §3
- Tier 1 簡單問答: skip framework，無嚴謹度層

---

## §12 與 既有 7 skill thin wrapper 對應

**取代策略**：保 skill 名 + INDEX 引用不壞；內部 body 改 thin wrapper 引述對應 framework entry。

| 既有 skill | thin wrapper 預設 framework |
|-----------|-------------------------|
| `/structured-tech-report` | `apply A3+ADR+Postmortem-hybrid`（原 13 段：TL;DR / 背景 / 流程 / 問題 / 5 Whys / ADR / 修改 / 比較 / Step→Verify / 影響 / 風險 / 後續 / 結論）|
| `/weekly-report` | `apply Multi-Thread-Narrative`（原 17 段 4 主線：進展 / 問題 / 求協助 / 探索 + Layer 0-4） |
| `/pptx-build` | 場景 audience 識別 → `apply <framework>` 動態 — 既有 6 模板（improvement / comparison / executive / data_showcase / concept / academic）保留為 framework entry |
| `/results-report` | `apply Data-Showcase`（Hypothesis → Data → Stats → Caveats） |
| `/conclude-research` | `apply Verdict-Pyramid`（Pyramid + WRAP falsifier observable） |
| `/report` | `apply AI-Session-Companion`（timeline + key decisions + provenance） |
| `/myPPT` | 場景識別後 forward 到本 skill N1（總入口） |

**重要**：thin wrapper 內仍可加 skill 特有 hook（如 structured-tech-report 涉 .cpp 改動的 `/cpp-change` 觸發），但 narrative 骨架改用 catalog entry。

---

## §13 對話層級啟用機制

### UserPromptSubmit hook `scripts/hooks/narrative_frame_advisor.sh`

偵測 keyword（中英對照）：
- 中：「整理 / 報告 / 説明 / 彙報 / 總結 / 介紹 / 講解 / 解釋 / pitch / 答辯 / 簡報 / 教 / 寫 / 整合」
- 英：「explain / summarize / report / pitch / present / teach / walk through / integrate / outline」

注入 advisory：
```
[narrative-frame]: 本回覆建議先聲明 framework（執行 /narrative-frame N1）
                   參考 InterSubMod/.claude/skills/narrative-frame/references/scenario_to_framework.md
                   或用 Tier 2 inline 模式快速套（首行聲明 framework）
```

**不強制觸發**；只提醒。用戶可：
- (a) 跑完整 N1-N6
- (b) 用 Tier 2 inline 首行聲明（如「用 PREP：」）
- (c) 説「不用框架」直接走自由格式

### AGENTS.md §15 回應分級「框架聲明」維度

- Tier 1：skip
- Tier 2：首行聲明 framework + 結構化內容
- Tier 3：完整 N1-N6 + structured output

### CLAUDE.md §11（新章）「敘述框架預設啟用」

≥200 字且跨 ≥2 概念的整理性回覆 → 預設套 framework；用戶可任何時 override。

---

## §14 範例

| 場景 | 框架 | 範例檔 |
|------|------|--------|
| PI 1-on-1 報告新發現（5min） | SCQA + STAR per item | `examples/pi_report_scqa_example.md` |
| Debug postmortem | A3 + 5 Whys + Inverted Pyramid | `examples/postmortem_a3_example.md` |
| Lab meeting 説服轉方向 | Golden Circle + Duarte Sparkline | `examples/pitch_golden_circle_example.md` |
| 教 ISM 概念給新生 | Feynman + ELI5 + CPP | `examples/teaching_feynman_example.md` |

---

## §15 對應其他 skill 的流程圖

```
                          ┌──────────────────────────┐
 用戶 raw request ───────→ │ N1: 5W 場景識別          │
 N 份 source docs          └─────────────┬────────────┘
                                         │
                          ┌──────────────▼─────────────┐
                          │ N2: 框架推薦              │ ← references/scenario_to_framework.md
                          │     (主 + 1-2 備 + 為什麼) │ ← references/framework_catalog.md
                          └──────────────┬─────────────┘
                                         │ C2 ★
                          ┌──────────────▼─────────────┐
                          │ N3: source 萃取           │ → narrative-organizer agent (≥3 docs)
                          │     [F/O/I/U] tier        │ ← /scientific-rigor §2 L1-L5
                          └──────────────┬─────────────┘
                                         │
                          ┌──────────────▼─────────────┐
                          │ N4: mapping 表             │
                          │     重點 ↔ section         │
                          └──────────────┬─────────────┘
                                         │
                          ┌──────────────▼─────────────┐
                          │ N5: 套用 + 補齊           │ ← templates/<framework>_skeleton.md
                          │     output mode 多選       │
                          └──────────────┬─────────────┘
                                         │
                          ┌──────────────▼─────────────┐
                          │ N6: gap + 5 秒測試 + ack   │ ← references/design_principles.md
                          └──────────────┬─────────────┘
                                         │ C3 ★
                          ┌──────────────▼─────────────┐
                          │ Output:                    │
                          │ - framework 聲明           │
                          │ - 結構化 narrative         │
                          │ - source mapping           │
                          │ - ⚠ gap 標                 │
                          └────────────────────────────┘
```

---

## §16 不適用場景（DO NOT USE WHEN）

- 純 code edit / typo / build / commit（無敘述需求）
- single-line factual answer（如「`X` 等於多少？」「commit hash 是什麼？」）
- 已在 active framework 中（避免遞迴）
- C++ build / test / pre-commit hook（用 `/cpp-change` 6 步驟 PDD）
- 純 git diff review（用 `/code-reviewer`）
- 純文件路徑 lookup（用 Glob / Grep）

---

## §17 Quality Checklist — 交付 narrative-frame 套用結果前自我檢查

- [ ] 5W 場景已識別（C1 確認）
- [ ] 主框架 + 1-2 備選已列出 + 為什麼 + 不選 X 理由（C2 確認）
- [ ] N3 萃取每筆素材標 [F/O/I/U] + source `path:line`
- [ ] N4 mapping 每 framework section ≥1 對應重點（無對應 → ⚠ gap）
- [ ] N5 套用 output mode 明確（outline / script / slide / companion / inline）
- [ ] N6 自審 8 項全 ✓（或 gap 列明）
- [ ] 5 秒測試 + 3 秒法則 + Assertion-Evidence + 6-item pre-publish 全 ✓
- [ ] 過度宣稱 check（[O]/[I] 不用「證實 / 確認 / 解決」）
- [ ] 路徑前綴 `InterSubMod/...`
- [ ] Tier 3 結論性報告繼承 /scientific-rigor §2 + §3 + §4 + §7 + §8.4
- [ ] 用戶 C3 ack 完成

---

## See Also

- **Spec**: 本 SKILL.md
- **Catalog (SoT)**: `InterSubMod/.claude/skills/narrative-frame/references/framework_catalog.md`
- **Scenario lookup**: `InterSubMod/.claude/skills/narrative-frame/references/scenario_to_framework.md`
- **Business sources**: `InterSubMod/.claude/skills/narrative-frame/references/framework_business_sources.md`
- **Design principles**: `InterSubMod/.claude/skills/html-report-build/references/design_principles.md`
- **Plan**: `/bip7_disk/liaoyoyo2001/.claude/plans/scientific-rigor-skill-draft-lazy-muffin.md`
- **Multi-doc sub-agent**: `InterSubMod/.claude/agents/narrative-organizer.md`
- **Advisor hook**: `InterSubMod/scripts/hooks/narrative_frame_advisor.sh`
- **Companion 假説驗證三層樓**: `/pre-decision-audit` + `/implementation-notes` + `/run-evaluator`
