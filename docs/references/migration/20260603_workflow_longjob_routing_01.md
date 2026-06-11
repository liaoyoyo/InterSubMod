---
title: Dynamic Workflow 壽命限制 vs 長 compute job — 路由規則 + 操作建議
date: 2026-06-03
type: governance / routing-rule
status: active
evidence: workflow wf_13833d53-153（evidence + synthesis + critique；tech agent failed StructuredOutput 但技術事實由 prompt 餵足）；critique verdict=needs_trim 已套用
data_sources:
  - .claude/CLAUDE.md
  - .claude/workflows/cross_sample_benchmark.js
---

# Workflow vs 長 compute job — 路由規則

> **裁決（校準，非誇大）**：**不是「根本不相容」**。不相容的是「把**產數字的長 compute** 塞進 workflow `agent()` step」這個特定用法——不是 workflow 本身 vs 長 job。workflow 的強項是 fan-out **已算好的結果** 做匯總；長 compute 該在主回合可見地跑。已落地 §8 +1 行 + 鐵則。

## 為何那個特定用法不相容（機制）
1. workflow `agent()` step 是**短壽 bounded 執行** + 一律 **acceptEdits 黑箱** + **無 mid-run user input/pause**。
2. 把長 compute（ISM C++ 全跑 / BAM 處理 / `run_*.sh`）放進 step → 中斷則在黑箱裡**半完成、無法 Read-back 驗證** = **直接撞 §13.0 捏造根因**（產數字與驗證無法分離）。
3. Bash 單呼叫上限 ~10min；長 job 需跨-turn 背景執行，而那是**主 agent loop** 的能力（run_in_background 跨 turn re-invoke），非 workflow step。
4. **實證**：cross-sample ASM workflow 4 樣本曾**被殺改背景補跑**（memory `project_cross_sample_asm_reproducibility`：`wf_76655c69-bea 後者4樣本被殺改背景補跑`）。

## 操作決策樹
| Q | 情況 | 用什麼 |
|---|------|--------|
| Q1 | 無長 compute | §8 既有 4 行（fan-out→workflow / 互動→主回合 / Hard Gate→主 agent / 5min→主回合）|
| Q2 | 長 compute · 單樣本/單 pipeline | 主回合 `Bash(run_in_background)` 或 `scripts/run_*.sh`；§1 用戶明示直接跑+一行告知，模型自判仍 🔴。**不入 workflow** |
| Q3 | 長 compute · 跨 ≥3 樣本同質 | **Hybrid 兩段**：① `parallel-benchmark` agent **實跑腳本**（須先 resource preflight；N 個 30-45min BAM 平行恐撞 CPU/mem → **預設改循序背景跑**）落檔 → ② workflow 只 fan-out **讀已落檔結果** 匯總 |
| Q3p | 長 compute · subset pilot（task type A）| 主回合背景跑 subset → Read 驗 → 不一定需 workflow（量小直接匯總）|
| Q4 | 長 compute · 需中間決策 / 含 Hard Gate | 全程主 agent 編排，連匯總都不入 workflow |
| Q5 | 純 fan-out 已算結果 / 文獻 cross-check / NO-GO stress-test（無 compute）| **workflow 本命**：`workflow` keyword 觸發 / `cross_sample_benchmark.js` |

## 標準 Hybrid（落地 §13.0 序列鐵則）
```
① 主回合 Bash(run_in_background) 或 scripts/run_*.sh 跑長 compute → 數字落 .json/.tsv
② Read 回真值，確認非 error/INCONCLUSIVE/未完成（§13.0 step 2）
③ 才 workflow fan-out：agent() step 只「讀已落檔結果 + 算 delta + 匯總」，不重跑 compute
```
（長 compute 與 workflow 在不同階段、不同 tool batch，與 §13.0 物理隔離一致）

## restraint — 砍掉的 over-build（critique）
- ❌ **不加** stage-2 grep PreToolUse 偵測 hook（advisory-only、workflow 本就繞 exit-2、假陽性）—— 判斷錯誤屬**路由類非 §13 數字捏造類**，純文字規則即足（§13 是因「同 session 純文字失效三次」才上機械層；此處無此實證）。
- ❌ **不加** `cross_sample_benchmark.js` 的 `args.mode='read-only-aggregate'` 結構守衛（2-模板 harness 過度工程）。
- ❌ **不改** AGENTS.md §4 五問（跨檔 pointer = scope creep；§8 規則自足）。
- ✅ 純文字規則 + 模板 1-line pointer **永久維持**，除非出現實際 in-the-wild 事故才升機械層。

## 落地
- `.claude/CLAUDE.md §8`：表 +1 行（產數字長 compute）+ 「長 job × workflow 鐵則」段。
- `.claude/workflows/cross_sample_benchmark.js`：頂部註解 +1 line pointer（指向 §8 鐵則）。
- 完整 workflow 分析：`wf_13833d53-153`（含 evidence + decision tree + critique）。
