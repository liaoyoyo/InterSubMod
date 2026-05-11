---
name: structured-tech-report
description: 工程／研究修改的結構化技術報告 skill — 把「問題→根因→修改→驗證→影響」流程整理為 13 段式報告，給「具備基本背景但不熟此專案」的讀者讀。整合 Toyota A3 + ADR + Google SRE Postmortem + Diátaxis 多受眾分層；繼承 InterSubMod doc-standards／v5_audit_suite 範本。USE WHEN：「寫技術報告」「整理修改」「13段報告」「self-phasing 報告」「pipeline 變更說明」「PR/commit 整理」「整理 bug fix」「整理工程改動」。與 results-report（實驗結果分析）／weekly-report（週進度多主題）／conclude-research（研究收尾假說驗證）／report（AI 對話紀錄）責任不同：本 skill 專責**單一工程／流程修改的深度敘事**。SKIP WHEN 多主題週進度（用 weekly-report）、實驗結果報告（用 results-report）、研究收尾（用 conclude-research）、AI session log（用 report）、單行修改不需深度敘事、純 build / commit。
allowed-tools: Read, Write, Edit, Glob, Grep
user-invocable: true
tags: ["report", "documentation", "engineering", "rca", "adr", "postmortem", "13-section"]
---

# Structured Tech Report（結構化技術報告）

把工程／pipeline／方法學的單一修改整理成 **13 段式** 技術報告，對齊 InterSubMod 既有最佳範本（`v5_audit_suite/07_paired_ground_truth_concordance.md`）。

> **路徑硬性規則**：本 skill 產出的所有 .md 路徑列給用戶時必以 `InterSubMod/...` 前綴（對齊專案 UserPromptSubmit hook）。

---

## 為何需要

- 既有 5 個報告類 skill（`results-report`／`weekly-report`／`conclude-research`／`doc-standards`／`report`）**沒有覆蓋**「工程／流程修改的完整敘事」。
- `docs/experiments/in_progress/` 出現孤兒文件（短、無 INDEX、無 cross-link）— 缺乏標準骨架是主要原因。
- v5_audit_suite 19 份報告證明 13 段式結構可承載真實深度議題（self-phasing），但僅憑慣例維持，沒有 skill 強制執行。

---

## 何時不要用

| 場景 | 應改用 |
|------|--------|
| 跑完實驗後寫結果分析 | `/results-report`（聚焦 metrics 解讀） |
| 多主題本週進度 | `/weekly-report` |
| 假說驗證收尾 | `/conclude-research` |
| AI 對話本身的決策紀錄 | `/report` |
| 命名／元數據規格疑問 | `/doc-standards`（被本 skill 引用） |
| 跨樣本一致性觀察 | `/feature-layered-observation`（對應特徵的 10 章節） |

---

## 執行流程（5 Step + 自檢）

### Step 1 · 定位（4 個前置問題）

對用戶連問 4 題（**答完才進 Step 2**）：

1. **受眾**：誰會讀？（PI／工程同事／外部審稿／AI session）→ 影響 §7 雙語語氣
2. **報告類別**：bug fix／pipeline 變更／方法學改進／架構決策／其他 → 影響根因章用 5 Whys 還是 ADR Decision Drivers
3. **修改範圍**：單檔／多檔／跨模組／涉外部依賴 → 影響 §3 與 §8 的詳盡度
4. **是否含 .cpp 改動**：是 → 強制呼叫 `/methodology-audit` + `/cpp-change`，並要求 §7.2 列 commit hash

> 13 段框架為**硬性預設**，不再詢問長度。

### Step 2 · 盤點（引用檔案先列出）

使用 `Grep` / `Glob` 列出將被引用的：
- `.md`（既有報告、INDEX、規範）— 路徑前綴 `InterSubMod/...`
- `.cpp` / `.hpp` / `.h`（含 file:line）
- `.py` / `.sh`（分析腳本與 pipeline）
- `.tsv` / `.csv`（數據出處）
- `.png`（圖表，含落點 `figures/{section}/...`）
- MEMORY 條目（`memory/project_*.md`）

**輸出 Step 2 草稿**到報告草稿頂端 HTML comment：

```html
<!--
build_date: YYYY-MM-DD HH:MM
agent: <claude session id 或 free text>
status: draft|validated|finalized
inputs:
  - InterSubMod/path/to/source1.md
  - InterSubMod/src/core/File.cpp:L120-L188
outputs:
  - InterSubMod/docs/reports/.../<this_report>.md
  - figures/<section>/figXX_*.png
verdict: <pending|POSITIVE|NEGATIVE|CONDITIONAL>
-->
```

### Step 3 · 填骨架

呼叫 `references/13_section_template.md`，逐段填寫。**資料不足處**標 `> ⚠ 待確認：<具體問題>`，**不停工**。

雙語規則（§7 硬性）：
- §7.1 **非工程版**：白話描述「改了什麼」，不出現 file path、function name、commit hash
- §7.2 **工程版**：精確描述「具體技術調整」，含 file:line、commit hash、API 簽章變更

### Step 4 · 自檢（25 項 checklist）

對 `references/checklist.md` 25 項逐一 ✅／⚠／❌。**底線**：
- ≥ 22 ✅ → 可進 Step 5
- 18-21 ✅ → 修補 ⚠ 項後重檢
- < 18 ✅ → 退回 Step 3

### Step 5 · 登錄

1. 報告檔案落點按「報告類別」決定：
   - bug fix／pipeline 變更（已驗證） → `InterSubMod/docs/reports/validated/YYYY/MM/`
   - 方法學改進中 → `InterSubMod/docs/experiments/in_progress/YYYY/MM/`
   - 架構決策 → `InterSubMod/docs/decisions/YYYY/MM/`
2. **追加一行**至對應 INDEX：
   - `InterSubMod/docs/experiments/INDEX.md`（experiments）
   - `InterSubMod/docs/reports/research_landscape/00_INDEX.md`（若涉 14 結論）
3. 提示用戶確認 MEMORY 是否需新增 reference entry（不自動寫，避免汙染）。

---

## 輸入模式

| 模式 | 觸發詞 | 預設動作 |
|------|--------|---------|
| **fresh** | 「寫技術報告 about X」 | 從零建檔，Step 1 全跑 |
| **from-commit** | 「整理 commit abc123 的修改」 | Step 2 自動 grep `git show abc123 --stat` 抓檔案清單 |
| **from-pr** | 「整理 PR/branch 變動」 | Step 2 跑 `git diff main...HEAD --stat` |
| **revise** | 「重整 InterSubMod/docs/.../XX.md 為 13 段」 | 讀既有 .md → 對照 13 段空缺 → 補齊 |

---

## 13 段速覽（完整模板見 `references/13_section_template.md`）

| # | 段名 | 一句目的 |
|---|-----|---------|
| 0 | TL;DR | 3-5 行給沒空讀全文的人 |
| 1 | 報告目的 | 為什麼寫這份 |
| 2 | 系統背景 | 模組原本負責什麼 |
| 3 | 原本流程 | Before |
| 4 | 問題描述 | 在什麼情境出什麼錯，量化 |
| 5 | 根本原因 | 5 Whys 或 Ishikawa |
| 6 | 修改方向 | 候選方案 + 採用理由（ADR Decision） |
| 7.1 | 修改內容（非工程版） | 白話 |
| 7.2 | 修改內容（工程版） | file:line + commit |
| 8 | 新舊比較 | 表格 |
| 9 | 驗證方式 | Step → Verify 格式（外部可觀察） |
| 10 | 影響範圍 | 受影響使用者／資料／下游 |
| 11 | 風險與限制 | 仍未解的、副作用 |
| 12 | 後續工作 | Action items |
| 13 | 結論 | 3-5 句總結 |

---

## 必繼承的既有規範

- 命名 `YYYYMMDD_主題_NN.md` → 沿用 `/doc-standards`
- HTML comment metadata → 沿用 `v5_audit_suite/07_*.md` 範式
- 圖表落點 `figures/{section}/` → 沿用 `/weekly-report`
- 涉統計驗證的 §9 → 必呼叫 `/auc-confound-guard` 三關
- 涉 .cpp 改動的 §7.2 → 必呼叫 `/methodology-audit` + `/cpp-change`

詳見 `references/linkage_map.md`。

---

## 與其他 skill 的對應圖

```
                       ┌──────────────────┐
 用戶觸發本 skill ───→ │ Step 1: 定位      │
                       └─────────┬────────┘
                                 │
                       ┌─────────▼─────────┐
                       │ Step 2: 盤點       │ ← /doc-standards
                       └─────────┬─────────┘
                                 │
                       ┌─────────▼─────────┐
                       │ Step 3: 填骨架     │ ← references/13_section_template.md
                       │ 　涉 .cpp？        │ → /methodology-audit, /cpp-change
                       │ 　涉統計？         │ → /auc-confound-guard
                       └─────────┬─────────┘
                                 │
                       ┌─────────▼─────────┐
                       │ Step 4: 自檢       │ ← references/checklist.md
                       └─────────┬─────────┘
                                 │
                       ┌─────────▼─────────┐
                       │ Step 5: 登錄       │ → INDEX.md / MEMORY.md
                       └───────────────────┘
```

---

## 範例

見 `examples/longphase_TO_vs_V5_example.md` —
此 skill 首個示範案例：`InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md`

## DO NOT USE WHEN（v1.7 batch A）

- **一般週報**（多主題進度彙整）— 用 `/weekly-report`（17 段母稿）
- **研究結論收尾**（假說最終判定）— 用 `/conclude-research`（P6 COMMIT）
- **單實驗結果分析** — 用 `/results-report`（前置 `/results-analysis`）
- **AI session 紀錄** — 用 `/report`（docs/provenance/ai_sessions/）
- **PPTX 投影片** — 用 `/pptx-build`

## Quality Checklist — 交付 13 段技術報告前自我檢查（v1.7 batch B）

- [ ] 13 段全寫齊（不可跳；缺段加 `(N/A，理由：...)` 而非靜默省略）
- [ ] §Evidence Layer 1-4 對齊 anchor #3「5 段骨架」（Exec / Background / Evidence / Limitations / Conclusion）
- [ ] §Limitations 誠實揭露範圍（已知 confound / 樣本限制 / 假設未驗）
- [ ] §Decision Log entry 格式（如有 ADR 性質：date / decision / rationale / alternatives considered）
- [ ] **每個 claim 引 evidence file path**（不可孤立陳述）
- [ ] 路徑 `InterSubMod/-` 前綴規範遵守
- [ ] reviewer-friendly：不熟此專案的人可讀懂（術語首次出現解釋）
- [ ] 4 軌證據鏈覆蓋表（anchor #1）：Statistical / Cross-sample / Mechanism / Orthogonal 各對應 evidence
