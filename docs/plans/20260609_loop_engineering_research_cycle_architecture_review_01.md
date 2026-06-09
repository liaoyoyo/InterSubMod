---
title: Loop Engineering 對照 — InterSubMod 研究循環架構的「能否改進與穩定」評估
date: 2026-06-09
status: ANALYSIS_COMPLETE / 1-candidate-pending-user-decision
type: ADR (Architecture Decision Record) + Restraint Audit
owner: liaotzuyu000@gmail.com
framework: Verdict-Pyramid (BLUF → 證據 → 候選裁決 → 單一建議)
data_sources: >
  scripts/harness_health.py (8-light dashboard, grep-verified),
  .claude/settings.local.json (39 hooks / 0 cron),
  state/schemas/plan.schema.json (expected_effect.stop_criteria),
  research/autoresearch/research_direction.md (no-auto-execute rule),
  workflow wf_67aa09b0-29b (17 agents: CC/Codex 文件驗證 + 科研風險文獻 + 對抗裁決)
provenance_note: >
  harness 計數與機制 = grep 驗證 (L1)。industry 數字/arXiv ID = workflow subagent WebSearch
  回報 (L3，未由主 agent 逐一重抓)；核心論證不依賴精確數字，只靠質性發現 + L1 harness 事實。
---

# Loop Engineering 對照評估 — 研究循環架構

> 對應任務：「依照 Loop Engineering 文章 + 其他網路方法，確認是否可建立/改進現在的研究流程循環，並分析是否可改進與穩定整個架構。」

---

## §0 BLUF（一句結論）

**這個 harness 不需要在上面再裝一層「自動迴圈」——它本身就已經是一個迴圈**：一個*人類節奏、預註冊停止條件、外部裁判、反捏造閘控*的研究迴圈，而且**獨立收斂到了 2025-2026 文獻推薦的「無人值守科研迴圈」幾乎全部安全護欄**，同時*刻意排除了文章裡唯一危險的那一塊（自動執行心跳）*。

文章的五要素 + 記憶層對照下來：**5/6 已超前文章、1/6（自動心跳）是刻意留白的安全選擇**。7 個候選改進經對抗裁決後 **3 已覆蓋 / 3 拒絕 / 1 採納（附護欄）**。唯一值得做的是 **C7：把既有 `harness_health.py` 漂移偵測改成「變動觸發」的唯讀心跳**——純穩定化、零捏造面、需用戶批准。

---

## §1 文章的「Loop = 5 + 1」模型（精簡）

| # | 要素 | 文章定義 |
|---|------|---------|
| 1 | **Automations（心跳）** | 排程自動跑 discovery+triage，findings 進 inbox，空跑自動歸檔；`/loop`（週期）+ `/goal`（跑到可驗證停止條件，獨立小模型判 done）+ cron + GitHub Actions |
| 2 | **Worktrees** | 平行隔離，避免互相覆寫 |
| 3 | **Skills** | 把專案知識寫下來（SKILL.md），不再每次重講 |
| 4 | **Plugins/Connectors** | MCP 接真實工具：開 PR、更新 ticket、Slack |
| 5 | **Sub-agents** | maker 寫、*不同的* checker 驗（生成/評估分離）|
| 6 | **Memory** | 對話之外的 markdown/board，「repo 記得，模型會忘」 |

文章自己的警告（關鍵，與本評估高度相關）：token 成本暴衝、slop 風險、「**Verification is still on you**」、comprehension debt（理解力腐蝕）、cognitive surrender（認知投降）、「**Build the loop. Stay the engineer.**」、「兩個人同樣的 loop 可得相反結果」。

---

## §2 五+一要素 → harness 覆蓋度（逐要素盤點，evidence grep-verified）

| 要素 | 覆蓋度 | 實作證據（真實 artifact） | 缺口 |
|------|:---:|------|------|
| **Skills** | 🟢 **STRONG_EXCEEDS** | 47 skills 含 P0-P6 七階段狀態機 + 34 個含 dependency-chain header + 326 reference 檔；不只「靜態知識」而是可執行的 cycle state machine | 無 |
| **Sub-agents（maker/checker）** | 🟢 **STRONG_EXCEEDS** | 18 agents；checker 角色**展開成 5 個正交 fresh-context 唯讀驗證者**（evaluator/reviewer/methodology-reviewer/security-reviewer/reproducibility-audit），全部 tool-locked 無 Write；4-state status enum | 無（**比文章的「一個 checker」強**） |
| **Memory** | 🟢 **STRONG_EXCEEDS** | file-memory + MEMORY.md + CURRENT_FOCUS(SessionStart 注入) + active.json + cycles/ + evidence_ledger.jsonl(append-only, lineage) + hypothesis_queue + harness_health(8 燈) | 無 |
| **Worktrees** | 🟢 **STRONG_EXCEEDS** | 三層協調：12/18 agent frontmatter `isolation:worktree` + Workflow per-step `isolation:'worktree'` + 文件化 restraint（worktree 貴，僅平行衝突才用）| 無 |
| **Connectors** | 🟢 **STRONG** | knowledge/biorxiv/ensembl MCP + claude.ai(Gmail/Calendar/Drive/Notion) + context7 | 「自動開 PR/更新 ticket」的*寫側*未 wire（**刻意**：github plugin `enabled=false`）|
| **Automations（心跳）** | 🟡 **PARTIAL** | 39 hooks 全 7 event-driven；`headless-research` agent；Dynamic Workflow；`/loop`+`CronCreate` 工具*可用* | **唯一真缺口**：`CronList`=空、無 timer-driven 心跳；health/cycle-state/kb_freshness **只能手動觸發** |

---

## §3 核心洞見：harness 已經是一個迴圈，且取了「安全子集」

文章把 loop engineering 講成「在 harness 之上再蓋一層心跳」。但實測顯示，**本 harness 的 P0→P6 七階段狀態機就是 loop 的 body**，且兩個關鍵設計證明它*刻意*只採用 loop engineering 的安全子集：

1. **`/goal` 的「可驗證停止條件」schema 早已存在且更嚴謹。**
   `plan.schema.json` 的 `expected_effect` 強制 `{metric, min_threshold, direction, stop_criteria}` + `confound_sweep_plan` + `conditional_branches` + `expected_fail_mode_fallback`（OSF 式預註冊）。這比 `/goal` 強：`/goal` 的 checker 預設 **Haiku 且只看 transcript、不獨立讀檔跑工具**；我們的 P5 `/run-evaluator` 是 **fresh-context 唯讀 agent 真的去讀檔驗證**，再經 `pre_tier_upgrade_check.sh`（exit-2）擋在人類 P5→P6 之前。→ **C2 已覆蓋。**

2. **「自動自餵」被白紙黑字拒絕。**
   `research_direction.md` 啟動規則：**「所有 queue 項目需用戶手動批准才啟動」+「本檔案僅作為候選列表，不寫 execution trigger」**。文章最激進的「automation auto-pick → auto-run → auto-conclude」在這裡是*被設計掉的*，不是沒想到。

**為什麼這是對的——研究 loop ≠ 軟體 loop（文獻佐證，L3 subagent WebSearch）：**

- 軟體的「done」有便宜確定性 oracle（測試綠/build 過/exit 0），loop 可重跑可信任；**研究的「done」沒有可執行 oracle**——真偽取決於 falsification、confound 控制、out-of-distribution 泛化、最終 orthogonal/wet-lab 確認，*產生結果的 agent 無法自我認證*。
- **驗證者 agent 本身不可靠**：REPRO-Bench 最佳 agent 在「可重現性評估」僅 21.4% 準確（低於隨機）；Sakana AI-Scientist 42% 實驗因 coding error 失敗；「更強的 reasoning model 反而 hallucinate 更多」。→ 「加個 critic agent loop 就自我檢查」被*反證*。
- **agentic multiverse 放大 p-hacking**（「Many AI Analysts, One Dataset」）：不同 agent 選不同變數/轉換 → garden-of-forking-paths 自動化。
- **premature victory / progress-as-completion** 是頂級失敗模式（agent 做到 ~80% 就宣告 done 不驗證）——這正是本專案 2026-06-01 §13 捏造事故的根因類別。

**驚人收斂**：文獻推薦的科研 loop 護欄，harness 幾乎全已實作 ——

| 文獻推薦護欄 | harness 既有對應 |
|------|------|
| Pre-registration / frozen plan | `plan.json` OSF 4 段 |
| Falsification step（negative control / shuffled null）| `confound_sweep_plan` + `/scientific-rigor` + excess-over-null |
| Held-out / LOSO | LOSO（正是它把 Phase2 Cycle1 從 ⭐3 降為 NEGATIVE）|
| Permutation / shuffled null | cross-sample excess-over-null（已做）|
| Multiple-comparison correction | `/scientific-rigor`（部分）|
| Provenance grep-gate（數字必可 grep）| **§13 + `number_provenance_check.sh`(exit-2)** |
| compute 與 write 物理分離 | **§13.0 鐵則** |
| External judge（非自評）| **5 個 fresh-context 唯讀 evaluator** |
| 拒絕 reopen falsified hypothesis | `research_direction_guard.sh` + tombstone check |

→ harness 在*沒有*危險 auto-execute 的前提下，已獨立達到業界安全 loop 的標準。文章對它的作用是**驗證**，不是揭露缺口。

---

## §4 七個候選改進的對抗裁決（restraint-first）

> 對抗 reviewer 預設懷疑，目標是*防止過度採納*（延續 `feedback_harness_restraint_over_adoption`：21 研究 → 12 high-fit → 僅 4 存活）。本輪 **7 → 1 存活**，同模式。

| 候選 | 裁決 | 理由（一句）|
|------|:---:|------|
| **C1** 唯讀心跳（排程跑 health+triage 進 inbox）| **已覆蓋** | discovery+triage 內容已由 harness_health + /cycle-state + kb_freshness + /research-dashboard + SessionStart 建議完整提供；唯一 delta（timer）是低 ROI overlap，且 /loop+cron session-scoped 7 天過期不耐久 |
| **C2** plan.json machine-checkable stop condition | **已覆蓋** | schema 早有 `expected_effect.stop_criteria`，由 P5 fresh-context evaluator 檢查 + exit-2 tier gate |
| **C3** 全自動自餵研究迴圈（auto-pick→run→conclude→inject）| 🔴 **拒絕** | 違反 §1（NO-GO 是 Hard Gate 不可逆）+ §13.0/§13.7；無人值守科研捏造風險*四軸全中*；正是 research_direction.md 明文拒絕的形態 |
| **C4** 排程自動跑重 compute（ISM/BAM 過夜）| 🔴 **拒絕** | 違反 §1（模型自判長算=🔴）+ §8 長 job 規則（須主回合可見背景 Bash + Read 回真值）；過夜 errored/half-complete 落成「finding」= §13 捏造在機器尺度重演 |
| **C5** 外部 board（Linear/Notion）當 canonical memory | 🔴 **拒絕** | 本地 grep-able 檔案*是* §13/§8 機械閘賴以運作的承重層；搬到非 grep-able connector = 純退化 |
| **C6** 「/research-cycle」單鍵 orchestrator skill | **已覆蓋** | `/cycle-state` 已輸出 per-phase「下一步 invoke 哪個 skill」路由；包成單鍵會壓縮 Hard Gate 停點 + 增 comprehension debt，無新能力 |
| **C7** 把 harness_health 漂移偵測接成自動觸發 | 🟢 **採納（附護欄）** | 純唯讀穩定化，真實 drift 缺口，零捏造面；符合「修缺口非裝工具」 |

---

## §5 唯一建議：C7 — 變動觸發的唯讀健康心跳

**問題（跨多個要素的共通穩定化主題）**：盤點時 5/6 要素的 `stabilization_note` 都指向同一弱點——**invocation-dependence（只在手動觸發時更新）**。skills/checkers/memory/health 都是 event-driven，沒有 timer。後果：
- 一個 cycle 可停在 P3 不前進，跨 session 都沒人發現（`/cycle-state` 的 stale>7d 警告本身也要手動跑才看得到）。
- harness drift（skill/hook/agent 計數、Hard-Gate neutering、stale compile marker）只在有人想到跑 `harness_health.py` 時才被抓。

**解法（最小、唯讀、不碰研究）**：

```
新增 hook: scripts/hooks/health_drift_advisor.sh  (SessionStart)
邏輯:
  1. git diff --quiet 比對 .claude/{skills,agents,hooks,rules} + settings.local.json
     自上次 health snapshot 以來是否變動 (用 state/health_snapshots/ 的 mtime/commit)
  2. 若無 harness 變動 → exit 0 靜默 (絕大多數純 code-edit/單 doc session 不打擾, 對齊 §6)
  3. 若有變動 → 跑 harness_health.py 的「輕量燈號摘要」(非 686 行全 HTML)
     注入 1-3 行 advisory: "harness 自上次稽核後有變動, 燈號 X GREEN / Y YELLOW; 跑 /harness-health 看詳情"
  4. 順帶: 若有 active cycle last_advanced_at > 7d → 一行提醒
```

**護欄（對抗裁決要求，全部必含）**：
- ✅ **變動觸發，非每次 SessionStart**：git-detect harness 檔案有變才跑，避免對純 code-edit session 增啟動延遲/噪音（對齊 §6「純 code edits 不需詳細 context」）。
- ✅ **advisory-only，永遠 exit 0**：絕不 exit-2、絕不阻擋 session 啟動。
- ✅ **零研究數字**：只產 harness 維運遙測（計數/gate-truth drift），不寫 report/ledger/cycle state、不做 finding。
- ✅ **不是 `/loop`/`/goal`/cron**：是既有的 SessionStart event hook 範式（precedent: `session_start_inject_focus.sh`），非無人值守自動化。

**ROI**：把 harness_health 從「想到才跑」變成「一變動就自動提醒」，直接修掉盤點反覆出現的 invocation-dependence 漂移缺口。風險面 = 0（唯讀、advisory、change-gated）。

---

## §6 明確 DO-NOT 清單（tombstone，防再提案）

以下三項經對抗裁決*拒絕*，列此防未來重複提案（對齊 reopen-guard 精神）：

- 🔴 **不建全自動自餵研究迴圈**（C3）——`auto-conclude` 是不可逆科學決策，必須人類/機械 gate。
- 🔴 **不排程自動跑重 compute**（C4）——長 job 須主回合可見背景 Bash + Read 回真值（§8 + `feedback_workflow_lifetime_vs_long_jobs`）。
- 🔴 **不把 canonical memory 搬到外部 board**（C5）——本地 grep-able 檔案是 §13 反捏造承重層。

---

## §7 結論

| 問題 | 答案 |
|------|------|
| 能否建立研究流程循環？| **已存在**——P0→P6 七階段狀態機就是 loop body |
| 能否改進？| 能，但只有 **1 個**（C7 唯讀健康心跳）通過 restraint；其餘 6 個已覆蓋或危險 |
| 能否穩定整個架構？| C7 直接修掉跨要素共通的 invocation-dependence 漂移缺口；其餘穩定性（反捏造/Hard-Gate/external-judge/pre-registration）**已是業界安全 loop 標準且超前** |
| 文章對我們的價值 | **驗證 + 命名**現有設計，而非揭露缺口；最該內化的是它的*警告*（verification is still on you / stay the engineer），而非它的*工具清單* |

---

## §8 用戶決策與 resolution（2026-06-09 已定）

1. **C7 → 實作（已完成 ✓）**：`scripts/hooks/health_drift_advisor.sh` 寫好、wire 進 settings SessionStart、三情境測試通過、修掉註解觸發的 harness_health #2 RED false-positive、CLAUDE.md §4 計數同步 39→40 + 文件化 → harness_health **8 GREEN**。首次上線即抓到 G6/G1 7 天未推進（證明缺口真實）。
2. **「研究循環 skills 任務」澄清 → 雙層 SoT 漂移**：用戶意圖 = 確認「建立任務→分派→完成紀錄→下個任務」整個系統有無大問題。稽核結論 = 系統存在、迴圈閉合、**不丟資料**，但有 **雙 SoT 漂移**（機械 cycle 狀態機凍在 06-02 vs ledger/CURRENT_FOCUS 跑到 06-08；harness_health #4 有盲點未抓）。
   - **決議：接受雙層、明文化分工**（不物理刪 active.json，避免 cascade 進 user-managed 的 CURRENT_FOCUS）。規則落地 CLAUDE.md §6「🔑 雙層 SoT 分工」。
3. **延伸（2026-06-09 用戶追加）→ git 分支/commit 治理明文化**：補完 canonical 文件 §F（何時強制建議 commit：F1-F6）+ §G（衝突最小化+聚焦：shared-state 不入 commit / 完成即 commit / 一 branch 一主題 / worktree / WIP ≤2-3）+ 3 句決策流。檔：`InterSubMod/docs/references/migration/20260601_git_branch_commit_governance_01.md`；memory `feedback_git_branch_commit_governance` 同步。

**仍開放（未來，非本輪）**：是否把 maker/checker P5 驗證從 opt-in 改 default（盤點 stabilization_note 標出 Opus 4.8 預設不 spawn → P5 checker 分離目前 opt-in）；是否擴 harness_health #4 偵測 active-cycle phase staleness（C7 已部分覆蓋時間維度）。

---

> **Provenance footer**：harness 計數/機制/schema/rule = grep 驗證（L1，見 data_sources）。industry 統計與 arXiv ID（REPRO-Bench 21.4% / Sakana 42% / LangChain 60% / $47K loop / Codex automations 等）= workflow `wf_67aa09b0-29b` 17-agent WebSearch 回報（**L3，主 agent 未逐一重抓原始來源**）；本評估核心論證不依賴這些精確數字，僅依賴質性發現 + L1 harness 事實。`/loop`(v2.1.72+)、`/goal`(v2.1.139+, Haiku transcript-only checker, 無 spend cap, session-scoped 7-day expiry) = CC 官方文件，subagent 回報。
