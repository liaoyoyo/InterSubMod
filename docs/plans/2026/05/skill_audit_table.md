---
title: 40 Skill Audit Table — description form + evals.json
date: 2026-05-11
status: planning
type: audit
---

# 40 Skill Audit Table

**生成時間**：2026-05-11
**目的**：改進 #3 的 fan-out 基礎 — 哪些要改 description、哪些要補 evals、優先級排序
**參考**：Anthropic best-practice article 4 (Skills) + plan `/bip7_disk/liaoyoyo2001/.claude/plans/serialized-swinging-spindle.md`

## 範圍校正

- 原任務假設 38 skills；**實際盤點 40**（`/big7_disk/liaoyoyo2001/InterSubMod/.claude/skills/`）
- `grill-me` 目錄已刪除（git status 確認）— 不計入
- 額外新增：`validation-protocol`、`weekly-report`（皆已存在但 Tier A/B/C 原列表沒蓋到）

## 統計總覽

| 指標 | 數量 | % |
|------|------|---|
| 總 skill 數 | 40 | 100% |
| Tier A (核心協議/閘門) | 22 | 55% |
| Tier B (支援/分析) | 11 | 28% |
| Tier C (工具/新增) | 7 | 18% |
| 描述式 description (需改寫) | 11 | 28% |
| 觸發式 description (已 OK) | 29 | 73% |
| 已有 evals.json | 4 | 10% |
| 缺 evals.json | 36 | 90% |

**判定原則**（description form）：
- **Directive (OK)**：含「USE WHEN」「觸發：」「觸發詞：」「SKIP WHEN」「DO NOT USE WHEN」之一，且首句以動作 / 場景識別開頭
- **Descriptive (NEEDS REWRITE)**：純名詞片段開頭、無觸發場景、或僅一句功能描述未列觸發詞

---

## Tier A audit (22 skills，最高優先)

| # | name | description form | evals | trigger keywords | 改寫建議（1 句） |
|---|------|----------------|-------|------------------|-----------------|
| 1 | auc-confound-guard | directive | NO | 「feature AUC」「discriminative」「predictive power」「TP vs FP 區分」「residualize」 | (OK — confirm SKIP WHEN added) |
| 2 | check-staleness | directive | NO | 「check staleness」「P2 precheck」「驗 precondition」「跑 pilot 前」 | (OK — already directive with state path refs) |
| 3 | conclude-research | directive | NO | 「conclude-research」「結束研究」「撰寫結論」 | Runs research project closure (validate artifacts → report skeleton → INDEX/MEMORY/ledger update) when cycle reaches verdict. USE WHEN ⭐4-5 cycle wraps, P5→P6 transition. SKIP WHEN cycle still in P1-P4. |
| 4 | confirmation-protocol | directive | YES | Hard Gate / Gate / Review / FYI、刪檔/C++ commit/NO-GO/git push | (OK — already directive, evals shipped) |
| 5 | cpp-change | directive | NO | 「開始實作 [審查文件名]」「執行方案 B」 | (OK — has DO NOT USE WHEN; consider adding USE WHEN keywords) |
| 6 | cycle-init | directive | NO | 「cycle init」「new cycle」「P0 REGISTER」「啟動新假說 cycle」 | (OK — distinct from /init-research clarified) |
| 7 | cycle-state | directive | NO | 「cycle state」「dashboard」「現在哪些 cycle 在跑」「active cycles」 | (OK — read-only label clear) |
| 8 | feature-layered-observation | directive | NO | 「觀察特徵」「feature layered」「TP FP 分層」「P3 PILOT」 | (OK — full Step 1-6 protocol described) |
| 9 | html-preview | directive | YES | 「想看看排版」「給 PI 看 preview」「html preview」「companion HTML」 | (OK — shipped with full SKIP WHEN list) |
| 10 | infra-ops | directive | YES | pipeline preflight、post-OOM、server migration、/tmp suspicions | (OK — shipped with SKIP WHEN) |
| 11 | init-research | directive | NO | 「init-research」「新專案」「多週研究」「長期研究」 | (OK — distinct from /cycle-init clarified) |
| 12 | inject-hypothesis | descriptive | NO | (none — pure description) | Injects new hypothesis into hypothesis_queue.json. USE WHEN 「注入假設」「新假設」「加入待測佇列」「from observation」「from paper」. SKIP WHEN hypothesis already registered or in active cycle. |
| 13 | methodology-audit | directive | NO | 「審查 XX 方法」「評估是否要改 XX」、修改 .cpp/.hpp 前 | (OK — already directive) |
| 14 | multi-sample-consistency | directive | NO | 「cross-sample」「multi-sample」「7 樣本驗證」「跨樣本一致」 | (OK — clear trigger list) |
| 15 | myPPT | directive | NO | 「我要做簡報」「教授報告」「週報」「PPT」「向教授彙報」 | (OK — dispatcher role explicit) |
| 16 | pivot-direction | directive | NO | 「換方向」「pivot」「我注意到 X 現象」「這方向不行」 | (OK) |
| 17 | pptx-build | directive | NO | 「簡報」「PPT」「PPTX」「deck」「教授報告」「投影片」 | (OK — should add SKIP WHEN: weekly-report 母稿尚未產出時) |
| 18 | provenance-tier-audit | directive | NO | 「週報 evidence」「tier 分佈」「audit ledger」「provenance 審計」 | (OK — system-level vs cycle-level boundary clear) |
| 19 | research-context-loader | directive | YES | 「延續研究」「研究方向」「結論驗證」「landscape」「過去結論」 | (OK — shipped with SKIP WHEN) |
| 20 | research-dashboard | directive | NO | 「研究看板」「dashboard」「研究狀態」「現在在哪」「給我 status」 | (OK — should add boundary vs cycle-state) |
| 21 | research-loop | directive | NO | 「開始研究迴圈」「research loop」「測試新假設」「P1 PLAN」 | (OK — boundary clear) |
| 22 | run-evaluator | directive | NO | 「run evaluator」「P5 evaluate」「retraction risk」「tier 升級前」 | (OK — MANDATORY label clear) |
| 23 | verification-loop | directive | NO | 「驗證程式碼」「verify」「check quality」「PR 前 review」 | (OK — boundary vs /validate command clear) |
| 24 | weekly-report | directive | NO | 「週報」「weekly report」「整理本週」「向教授報告」「lab meeting」 | (OK — has DO NOT USE WHEN) |

---

## Tier B audit (11 skills)

| # | name | description form | evals | trigger keywords | 改寫建議（1 句） |
|---|------|----------------|-------|------------------|-----------------|
| 25 | citation-verification | directive | NO | 「撰寫論文」「檢查既有引用」「產出 bibliography」 | (OK — has USE WHEN list) |
| 26 | data-audit | directive | NO | 「檢查檔案組織」「audit」、研究完成後整理 docs/ | (OK — has USE WHEN list) |
| 27 | doc-standards | directive | NO | 建立新 .md 檔案或重組文件目錄時 | (OK — has DO NOT USE WHEN) |
| 28 | known-pitfalls | directive | NO | OLS/residualization、VCF 來源、特徵設計、AUC 分析時 | Add explicit USE WHEN keywords: 「known pitfalls」「踩雷」「avoid mistake」「common bug」. Currently only mentions situational triggers. |
| 29 | memory-consolidation | directive | NO | 「整理記憶」「memory consolidation」「記憶太多」「MEMORY.md 太長」 | (OK) |
| 30 | observation-analysis | directive | NO | 「建立新觀察分析」「O-系列」「.py 分析腳本」 | (OK — has USE WHEN) |
| 31 | problem-framing-ideation | directive | NO | 「腦力激盪」「brainstorm」「研究構想」「gap analysis」「5W1H」 | (OK — boundary 3-fold explicit) |
| 32 | results-analysis | directive | NO | 「分析結果」「統計分析」「比較效能」「check significance」 | (OK) |
| 33 | results-report | directive | NO | 「寫實驗報告」「summarize results」「實驗複盤」「results report」 | (OK — prerequisites 標明) |
| 34 | review-evidence | directive | NO | 「查閱證據」「review evidence」「過去測試結果」、研究方向決策前 | (OK) |
| 35 | structured-tech-report | directive | NO | 「寫技術報告」「整理修改」「13段報告」「pipeline 變更說明」 | (OK — boundary vs results-report/weekly-report/conclude-research/report 4-fold clear) |

---

## Tier C audit (7 skills)

| # | name | description form | evals | trigger keywords | 改寫建議（1 句） |
|---|------|----------------|-------|------------------|-----------------|
| 36 | fast-learning-coach | directive | NO | 「教我 X」「想學 Y」「快速學習」「解釋 X 給我聽」 | (OK — 觸發條件相對寬鬆 already noted) |
| 37 | image-gen | directive | NO | 「需要示意圖」「畫個流程圖」「補個 icon」「生圖」 | (OK — dual-track + SKIP WHEN clear) |
| 38 | image-vision-check | directive | NO | 「圖夠不夠好」「圖達標了嗎」「vision check」「品檢」 | (OK — SKIP WHEN clear) |
| 39 | report | directive | NO | 「寫報告」「session report」、Stop hook 提醒時 | Add SKIP WHEN: 寫週報 (用 weekly-report)、寫實驗報告 (用 results-report)、寫技術報告 (用 structured-tech-report). Currently no boundary statement. |
| 40 | validation-protocol | directive | NO | 「驗證假說」「validation protocol」「怎麼驗證」 | Add SKIP WHEN: 程式碼驗證 (用 verification-loop)、研究結論判定 (用 run-evaluator). Currently very short description. |

---

## Priority queue for T3.3-T3.6

**Phase 1 (Tier A, 22 skills — fan-out 6 parallel subagents):**
- **Group A1** (2 skills): inject-hypothesis (descriptive→rewrite), conclude-research (rewrite SKIP WHEN)
- **Group A2** (4 skills): cpp-change, cycle-init, cycle-state, init-research (boundary refresh)
- **Group A3** (4 skills): research-loop, research-dashboard, run-evaluator, check-staleness (P-cycle 套件)
- **Group A4** (4 skills): feature-layered-observation, multi-sample-consistency, methodology-audit, auc-confound-guard
- **Group A5** (4 skills): myPPT, pptx-build, weekly-report, pivot-direction
- **Group A6** (4 skills): provenance-tier-audit, verification-loop, html-preview (re-confirm), infra-ops (re-confirm)
- (confirmation-protocol, research-context-loader 已 shipped — 不重複改)

**Phase 2 (Tier B, 11 skills — fan-out 3 parallel):**
- **Group B1** (4 skills): citation-verification, data-audit, doc-standards, known-pitfalls (補 USE WHEN)
- **Group B2** (4 skills): memory-consolidation, observation-analysis, problem-framing-ideation, results-analysis
- **Group B3** (3 skills): results-report, review-evidence, structured-tech-report

**Phase 3 (Tier C, 7 skills — fan-out 2 parallel):**
- **Group C1** (4 skills): fast-learning-coach, image-gen, image-vision-check, report (補 SKIP WHEN)
- **Group C2** (1 skill): validation-protocol (補 SKIP WHEN + boundary)

## Skills already shipped with directive description (no rewrite needed)

| skill | shipped via | evals.json |
|-------|------------|------------|
| html-preview | improvement #2 Phase 2 | YES |
| infra-ops | improvement #5 | YES |
| research-context-loader | improvement #1 | YES |
| confirmation-protocol | improvement #2 | YES |

T3.6 evals.json mass-creation 仍需處理另外 36 skills。

## Notes / Anomalies

- **`grill-me`** directory deleted (git status: `D .agents/skills/grill-me/SKILL.md`) — 不在 40 計數內
- **`validation-protocol`** 與 **`weekly-report`** 不在原任務 prompt 的 Tier A/B/C 列表，本表將之列入 Tier C / Tier A 補齊
- **`pptx-build` / `myPPT`** boundary：myPPT 為入口 dispatcher、pptx-build 為下游製作 — 已分離，無重疊
- **`research-dashboard` vs `cycle-state`**：dashboard = project-level、cycle-state = cycle-level — 已在 cycle-state description 註明，dashboard 可考慮加對應註記
- **`init-research` vs `cycle-init`**：已明確分離（research/ vs state/cycles/）
- **`verification-loop` vs `/validate`**：已明確分離（程式碼級 vs 研究結論級）
- **首句動作 + USE WHEN + SKIP WHEN 三件式**：以下 4 skills 是完整三件式範本可參照：
  1. confirmation-protocol（英文 directive）
  2. html-preview（英文 directive）
  3. infra-ops（英文 directive）
  4. research-context-loader（英文 directive）
- **大多數中文 directive skills 缺 SKIP WHEN / DO NOT USE WHEN**：T3.3-T3.6 改寫時應補上反向觸發條件以提升 routing precision
