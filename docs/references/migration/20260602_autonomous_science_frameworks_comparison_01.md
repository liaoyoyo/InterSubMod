---
title: 自主科研 Agent 框架 vs InterSubMod harness — 架構對照（Community Comparison Round 3）
date: 2026-06-02
type: governance / comparison
status: active
evidence: L2（官方 blog/arxiv paper）/ L3（二手文章）；本 harness 全 L1
data_sources: docs/references/migration/COMMUNITY_COMPARISON.md
---

# 自主科研 Agent 框架對照 — 整體方法架構能否更好（Round 3）

> **前情**：`COMMUNITY_COMPARISON.md` Round 1/2 已比過 **Claude Code 研究 plugin**（flonat/matsengrp/K-Dense 等）。本 Round 3 補**真正的類比軸：自主科研 agent 框架**（有已發表架構），回答「整個方法/細節能否更好架構」。
> **一句話結論**：InterSubMod 在**人機協作嚴謹度 / 負結果紀律 / 反捏造 provenance / harness 自我治理**四軸**領先**已發表系統；它們的「全自主 + tree-search」是**不同設計目標（autonomy）**，照搬反而會**摧毀**本 harness 的核心強項。可借的只有 2-3 個輕量「概念」，非重構。

---

## 1. 對照對象（5 個已發表系統）

| 系統 | 架構核心 | 驗證/反幻覺 | 人類角色 | 證據 |
|------|---------|-----------|---------|------|
| **Sakana AI Scientist v2** | Agentic **Tree Search**（hypothesis/experiment 分支+回溯）；Idea-Gen / Experiment-Manager(sandbox code-gen) / Viz-Analysis / VLM-reviewer | LLM/VLM reviewer + 形式 peer-review | **全自主**（無人介入寫出 workshop paper）| L2 (arxiv 2504.08066) |
| **Google AI co-scientist** | **6 named agent**（Generation/Reflection/Ranking/Evolution/Proximity/Meta-review）+ Supervisor planner；3 相（Generation→Debate→Evolution）| **Reflection agent = 虛擬 peer reviewer** + **idea tournament**（pairwise Elo 辯論排名）| 人機協作（scientist collaborator）| L2 (Nature 2026 / arxiv 2502.18864) |
| **FutureHouse PaperQA2 / Robin** | agentic **RAG**（query→retrieve→rerank RCS→cite）；**citation traversal**；Robin 串 Crow(lit)+exp agent | **每個 fact 附 source**（~9% error）；RAG=by-construction 反幻覺 | 人類給方向 | L2 (futurehouse.org / arxiv 2312.07559) |
| **Agent Laboratory** | 人 idea-in → 自主 report+code pipeline；可調人介入度 | report+code repo 產出 | 可變介入 | L2 (arxiv 2501.04227) |
| **Coscientist / ChatBattery** | 多模組（web/doc/code/robot）閉環到 wet-lab | **「log 所有 tool 版本+參數」reproducibility 非選項** | 閉環 + human-in-loop | L3 |

**Survey 定錨**（arxiv 2508.14111 + ClawdLab 2602.19810，L2）：現行系統**多為單 agent 或緊耦合 pipeline**，**缺**「持久研究群組 / 角色分工 / 結構化對抗辯論 / **claim 經獨立重複的累積驗證**」——這正是傳統科學知識生產方式。

---

## 2. InterSubMod 對照：領先 / 相等 / 落後

### 2.1 本 harness **領先**（已發表系統普遍缺）

| 軸 | 本 harness（L1）| 對照系統 | 判定 |
|----|------|---------|------|
| **人機協作嚴謹度** | 影響×信心矩陣 + Hard Gate + 6 類 task-type + NO-GO 必停 | AI Scientist v2 = fire-and-forget 全自主（品質受質疑）| **領先**（高風險癌症研究本就該人 gate）|
| **累積驗證 / 負結果紀律** | evidence_ledger（append-only）+ multi-sample-consistency + LOSO + Productive Failure + Concluded 區防 re-investigation | survey 明指這是現行系統**最缺**的一塊 | **領先（decisive）** |
| **反捏造 / 數字 provenance** | 3 層（fill_report by-construction / number_provenance gate / audit）+ §13.0/§13.7 | PaperQA2 做 citation grounding（~9% err）但**無數字 provenance gate** | **領先（更具體）** |
| **harness 自我治理** | harness_health 6 燈 + creation_guard + registry_sync + drift 偵測 | 無系統做「harness 自我稽核」| **領先（獨有）** |
| **restraint / 防 over-adoption** | feedback_harness_restraint + 對抗驗證採用 | 學界傾向加 agent/能力 | **領先（哲學）** |

### 2.2 **相等**（概念對應，實作不同）

| 軸 | 本 harness | 對照 | 備註 |
|----|------|------|------|
| 假說生命週期 | hypothesis_queue + inject + pivot-direction | co-scientist Generation/Evolution + Proximity clustering | 概念相等（見 §3 borrow B6）|
| 對抗審查 | 4 read-only verifier + pre-decision-audit | co-scientist Reflection + AI Scientist VLM-review | 相等但本 harness 較**散落**（見 B7）|
| generator/evaluator 分離 | run-evaluator 6-component **確定性**（非 LLM-judge）| co-scientist Ranking（LLM Elo）| 本 harness 確定性，抗 judge 漂移 |

### 2.3 **落後 / 缺**（誠實揭露）

| 軸 | 對照系統有 | 本 harness | 判定 |
|----|----------|-----------|------|
| **Agentic Tree Search**（實驗設計分支+回溯）| AI Scientist v2 核心 | 7-Phase **線性** waterfall + pivot（無系統 backtrack）| ⚠ **paradigm 差異非缺陷**：那是全自主用；人機協作下 linear+pivot 更合適。**不採**（見 §4）|
| **idea tournament**（pairwise Elo 對抗排名）| co-scientist | hypothesis_queue 用**手設 priority 整數** | 🟡 可借輕量概念（B6）|
| **citation traversal / RCS**（lit-RAG SOTA）| PaperQA2 | researcher + knowledge MCP + citation-verification | 🟡 P2 論文期（已在 Round1 G2/G3）|
| **閉環 wet-lab** | Coscientist | N/A（計算基因體，無 robotic lab）| 不適用 |

---

## 3. 可借鑑的「概念」（輕量，非重構；restraint-tagged）

| # | 借鑑 | 來源 | 處置 |
|---|------|------|------|
| **B6** | **假說 pairwise 對抗排名**（取代純手設 priority 整數）— 用 novelty×feasibility×impact 兩兩比較定序 | co-scientist idea tournament | 🟡 可選：併進 `/inject-hypothesis` / `/pivot-direction`（輕量，非全 Elo 引擎）|
| **B7** | **正名一個 Reflection / red-team 步驟**（G7 devils-advocate 當初 deferred）— cycle-init 前獨立紅隊假說 | co-scientist Reflection agent | 🟡 可選：pre-decision-audit 已部分覆蓋，可加一條明確 red-team gate |
| **B8** | **「每個 fact 附 source」延伸到文獻 claim** | PaperQA2 citation-by-construction | 🟢 與 §13 反捏造同源；citation-verification 可加 source-for-each-fact |
| **B9** | **log 所有 tool 版本+參數**（reproducibility 非選項）| Coscientist | ✅ 已有 binary_versions.jsonl + commit-hash binding；可在 pipeline-manifest 補強 |

---

## 4. 整體架構裁決（回答「是否可以更好架構整個方法」）

**裁決：不重構。本 harness 架構健全，且在 4 軸領先已發表系統。**

關鍵論證：
1. **設計目標不同 → 照搬會傷**：AI Scientist v2 / co-scientist 優化 **autonomy（最少人介入）**；InterSubMod 優化 **human-augmented rigor（單一研究者 + 高風險癌症 + restraint）**。把 harness 改成「全自主 tree-search」會**摧毀**它最強的人 gate + NO-GO + 負結果紀律。
2. **本 harness 已具備 survey 點名「現行系統最缺」的能力**：cumulative verification through independent replication（evidence_ledger + LOSO + multi-sample）。這是領先，不是落後。
3. **真正的改進方向不是「更多 agent / 更自主」，而是 §3 的 2-3 個輕量概念**（假說對抗排名 / 正名 red-team / fact-source 延伸），+ 持續修 latent bug + drift（延續 restraint ROI）。
4. **論文期才需的**（citation traversal / LaTeX / referee）已在 Round 1 G2-G4 列 P2 採購單，不提前。

> 與 [[feedback_harness_restraint_over_adoption]] 一致：類比系統越強，越凸顯「加能力 ≠ 更好」；本 harness 的 edge 在嚴謹 + 治理，不在 autonomy。

---

## 5. Lineage
- 比較：Claude Code Opus 4.8 session 2026-06-02；主迴圈 targeted survey（非 workflow）。
- 來源 WebSearch：Sakana AI-Scientist-v2 (arxiv 2504.08066)、Google AI co-scientist (arxiv 2502.18864 / Nature 2026)、FutureHouse PaperQA2 (arxiv 2312.07559) / Robin、Agent Laboratory (arxiv 2501.04227)、survey 2508.14111、ClawdLab 2602.19810。
- 證據：對照系統 L2/L3（官方 paper/blog；未跑行為測試）；本 harness L1（讀實際 SKILL.md/hook）。
- 銜接 Round 1/2：`COMMUNITY_COMPARISON.md`。
