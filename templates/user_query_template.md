<!--
建立時間: 2026-05-18
parent: InterSubMod/.claude/skills/scientific-rigor/SKILL.md (P3 audit M15 fix)
purpose: 用戶端 query template — 6 欄結構化提問，降低 AI 反問輪次（一次提問減 3-5 輪反問）
status: optional template, copy when needed
-->

# User Query Template — 6 欄結構化提問

> **使用方式**: 對複雜任務 / 跨 session 接續 / 研究方向決策，複製此 6 欄填空後一次貼給 AI。
> **不適用**: 純 build / commit / 單檔 doc / 簡單問答（這些直接問即可）。
> **預期效果**: 一次提問減少 3-5 輪反問；AI 回應直入主題不需先 5W1H 反推。

---

## Template

```
## 1. 目標 (Goal)
我要達成什麼？（一句話 + 可測量結果，對齊 SMART）
例：「驗證 V6 binary 在 6 樣本 caller F1 無 regression（|ΔF1| < 0.005）」

## 2. 已做 (What's Done)
目前進度 / 已驗證 / 已產出 artifacts？（含路徑）
例：「2026-05-15 已跑 HCC1395 5 樣本 V6 ISM，結果在 InterSubMod/research/v6_bam_tpfp_hp_loh_cn/step1_master_three_way.tsv」

## 3. 問題 (Problem Encountered)
卡在哪？預期 vs 實際差異？
例：「COLO829 truth set 0600 權限，無法擴展到 6 樣本」

## 4. 懷疑原因 (Suspected Cause)
我推測是 X 或 Y（不確定，需驗證）
例：「ONT 原廠權限政策 / 或檔案 owner 未轉移」

## 5. 限制 (Constraints)
時間 / 成本 / 不能動的部分 / Hard Gate
例：「W3 deadline 5/22；不能 sudo chmod；不能影響其他樣本」

## 6. 期望答案形式 (Expected Answer Format)
回答需要什麼結構？（決策表 / 程式碼 / Step→Verify / 短答）
例：「3 個 Path 比較表（A/B/C），含工時+風險，最後標 ⭐ 推薦」
```

---

## 簡化版（4 欄速答用）

對於中等複雜任務（不到 6 欄但比一句問題多）：

```
目標：<一句>
已做：<bullet>
問題：<bullet>
期望答案：<格式>
```

---

## 對應業界框架

| 本模板欄 | 對應方法論 |
|---------|----------|
| 1. 目標 | SMART 法則 (Specific/Measurable/...) |
| 1+5 目標+限制 | 5W2H (What/Why/Who/When/Where/How + How Much) |
| 2+3 已做+問題 | 解決問題流程 Step 1-3 (描述/影響/原因) |
| 4 懷疑原因 | 假說驅動 Hypothesis-Driven |
| 6 期望答案形式 | 30 秒法則 (要 TL;DR 還是完整分析) |

---

## 與 InterSubMod skills 的關係

- `/problem-framing-ideation`: 當用戶提問模糊時，AI 主動用 5W1H 反問補齊本模板各欄
- `/scientific-rigor §7.1 Pre-registration`: 研究方向開跑前的 H_預測 / 否證條件 / decision_threshold = 本模板 6 欄的特化（更嚴格）
- `/confirmation-protocol §1 暫停判斷矩陣`: 第 5 欄「限制」對齊 Hard Gate 觸發條件

---

## When to use vs not use

**Use（推薦）**:
- 跨 session 接續複雜任務
- 研究方向決策 / pivot
- 多選項比較（請 AI 給推薦）
- 跨多個 skill / agent 的協作任務
- 首次描述新研究方向

**Don't use（過度結構化）**:
- 「make 一下」「git status」這類純命令
- 簡單事實查詢（「X 是什麼？」）
- 對話延續（「然後呢？」「OK 繼續」）
- 已有 plan / artifact 路徑明確的 follow-up

---

## 元層說明

本模板源自 2026-05-18 P3 Methodology Audit M15 fix（Eisenhower priority P0）。原始分析見 `InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/p3_methodology_alignment.standalone.html`。
