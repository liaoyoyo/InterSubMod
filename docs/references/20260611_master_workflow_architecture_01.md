---
title: InterSubMod 主工作流架構（唯一 SoT）— Agentic Loop 四層 + Git 多工/分支/commit/合併決策
date: 2026-06-11
status: active / authoritative
type: governance / master-workflow
supersedes_as_index: 整併並指向 git governance（§A-§G）+ loop-engineering ADR + cycle 狀態機
owner: liaotzuyu000@gmail.com
framework: Diátaxis(reference) + 4-layer agentic loop
---

# 主工作流架構 — 之後「所有任務」都走這一張圖

> **目的**：把分散在多份 governance / ADR / skill 的規則收斂成**唯一一張清楚的工作流**。
> 任何任務啟動 → 照本檔走。git 的「何時開不衝突分支 / 何時同資料多工 / 何時 commit / 何時矯正合併」在 §3 做成可機械化的決策表 + §4 鉤子強制。
> **權威分工**：本檔 = 入口地圖 + git 決策表。細節仍指向 — git §A-§G `InterSubMod/docs/references/migration/20260601_git_branch_commit_governance_01.md`；loop-eng `InterSubMod/docs/plans/20260609_loop_engineering_research_cycle_architecture_review_01.md`；cycle 狀態機 `/cycle-state`。

---

## §0 L0 一眼結論

**每個任務 = 一次 agentic loop：收集脈絡 → 規劃 → 執行 → 驗證 → 紀錄。** git 紀律疊在「執行/紀錄」上：**預設一任務一分支（多工時各開 git worktree），完成一個已驗證單元就 commit，矯正/合併等並行安靜後做。** 鉤子在每個關卡把「宣稱」拉回「可機械驗證的真實」。

---

## §1 L1 唯一的端到端工作流（一張圖）

```
[L1 脈絡]            [L2 規劃]              [L2.5 git 佈署]         [L3 執行]            [L4 驗證]           [L4 紀錄]
任務分類(§0 6類)  →  研究先於實作       →  分支/worktree 決策   →  單回合 or            →  done=可執行      →  commit(F1-F6)
context 載入        plan mode(調查→         (§3 決策表 A/B)         subagent/workflow      gate 不可繞過        → 矯正/合併(§3 D)
(CLAUDE/AGENTS/      列範圍/風險/測試)                              (§8 路由)              (§4 hooks)          → ledger/memory/INDEX
 CURRENT_FOCUS)      → 用戶批准                                     maker≠checker分離
        │                                                                                      │
        └──────── 任務越大繞越多圈；單純問答只走「脈絡→答」 ◄───────────────────────────────────┘
```

### 四層對應（你的 agentic-loop 框架 → 本 harness 既有件，不是新發明）

| 層 | 你的定義 | 本 harness 對應件（已存在）| 啟動方式 |
|---|---|---|---|
| **L1 脈絡工程** | CLAUDE.md 把是什麼/為什麼/怎麼做寫清楚 | `.claude/CLAUDE.md` + `AGENTS.md`（自動載入）+ `CURRENT_FOCUS.md`（SessionStart hook 注入）+ memory + KB MCP | 自動 |
| **L2 工作流設計** | plan mode；研究/規劃 ≠ 實作 | plan mode + 7-phase cycle 狀態機（P0 cycle-init→P1 research-loop plan.json〔OSF 預註冊+stop_criteria〕→P2 check-staleness→…）+ `/pre-decision-audit` + `/problem-framing-ideation` | 規格缺項≥2 必回問；研究方向先 plan |
| **L3 編排** | 多開/平行/subagent；worktree 隔離；乾淨脈絡 review | **§3 git worktree 決策** + 18 agents（generator/evaluator 分離）+ Dynamic Workflow + 5 個 fresh-context 唯讀 verifier | §8 路由表 |
| **L4 驗證** | 測試/CI/pre-commit 不可繞過；done=可執行 | §4 hooks（6 真 exit-2 Hard Gate + §13 反捏造）+ harness_health 9 燈 + §13.7 完成宣稱 gate | 自動 + 收尾 |

---

## §2 任務啟動 5 步（每個任務固定走）

1. **分類**（CLAUDE.md §0 6 類）：A pilot / B validation / C production / D handoff / E hotfix / F demo → 決定 default scope（B/C/D = 全基因組全樣本）。模糊 → AskUserQuestion 必問 scope。
2. **載脈絡**（L1）：純 code edit/單 doc/簡答 → 跳過詳細；研究分析 → `/research-context-loader`。
3. **規劃**（L2）：非 trivial → plan mode 先「調查→列受影響範圍+風險+測試策略」→ 用戶批准 → 才執行。**直接跳寫 = 最常見失敗（解錯問題）**。
4. **git 佈署**（L2.5 §3）：先判斷分支/worktree（**動手前決定，不是事後**）。
5. **執行→驗證→紀錄**（L3-L4）：maker≠checker；done=可執行；commit 走 §3-C；evidence_ledger 每輪記。

---

## §3 Git 決策（核心 — 4 張可機械化的表）

> 一句總綱：**一任務一分支；多工 = 多 worktree；完成一驗證單元即 commit；矯正/合併等並行安靜。**

### A. 何時「開新不衝突分支 / worktree」

| 情境 | 動作 | 為什麼 |
|------|------|--------|
| **主題 ≠ 當前 branch 主題** | **開新 branch** | 防主題混入錯 branch（git §A）|
| **跨多 commit 獨立 feature / 重構 / 實驗可能放棄 / 需獨立 review** | **開新 branch** | 隔離、可乾淨 review/PR |
| **同時要跑「第二條工作線」** | **開新 git worktree**（非只 branch）| 兩線各自 working dir + HEAD，不互覆 |
| **⚠ 另有並行 Claude session 在主 dir** | **必開 worktree 起 session**（見 §3-B） | 共用主 dir+HEAD = commit 落錯 branch（2026-06-09 實際事故）|
| 同主題增量 / 小修 / 當前 branch hotfix | **持續 commit 當前 branch** | 同一件事不需新 branch |

**Base**：延續工作從當前 HEAD；全新獨立從 `main`/`develop`。

### B. 何時「同資料多工」→ 一律 git worktree 隔離（**這是 06-09 事故的根治**）

| 多工型態 | 機制 | 指令 |
|---------|------|------|
| **多個 Claude session 同 repo 並行** | **各開 worktree 起 session，主 dir 只留一個** | `git worktree add ../wt-<topic> <branch>` 後在該 dir 起 session；或 Claude Code 原生 `--worktree`（見 §4 機制⑤，啟用前確認版本旗標）|
| **一 session 內平行 subagent 改同檔** | Agent `isolation:'worktree'` / Workflow step `isolation:'worktree'` | 已內建；agent 各得 fresh checkout 自動清理 |
| **平行跑同 benchmark 多樣本（不改檔）** | `parallel-benchmark` agent；**先 resource preflight**（§8）| 不需 worktree（唯讀），但 fan-out 要 cap（§5）|
| **worktree 隔離的是 CODE 不是 RUNTIME** | TMPDIR/scratch/大 BAM 中間檔 + shared live-state（CURRENT_FOCUS/ledger/queue）仍共享 | 各 session 設自己 TMPDIR；shared-state 留用戶（§3-C）|

### C. 何時「清楚 commit」（F1-F6 觸發；advisory）

| # | 觸發 | 註 |
|---|------|----|
| F1 | **一個已驗證邏輯單元完成** | 立即 commit，別累積 |
| F2 | **切主題 / 開新 branch 前** | 先把當前主題 commit 乾淨 |
| F3 | **不可逆 git 操作前**（checkout 切 branch/reset/rebase/merge/大刪） | 先 commit 保護 |
| F4 | **「我的檔」累積 >~10 或跨 >3 邏輯單元** | 分批 commit（判準看「我改的」非總 dirty）|
| F5 | **C++ 編譯通過後** | 別讓編譯狀態漂移 |
| F6 | **session 結束 / 交接前** | 交接乾淨 |

**commit 紀律**：只 `git add <列檔>` 禁 `-A`；**shared live-state 檔不碰**（CURRENT_FOCUS/hypothesis_queue/evidence_ledger/active.json/*.log/state snapshots）留用戶；半成品/未驗證**不 commit**（§13.7）；message `type(scope): 摘要` + body + `Co-Authored-By:`。

> ⚠ **澄清「不碰」≠「不版控」（2026-06-15 audit D1-1 — 消除與 §6 的表面矛盾）**：shared live-state **是被 git 版控的**（CURRENT_FOCUS.md / evidence_ledger.jsonl 即 §6 雙層 SoT 的「敘述層權威 SoT」，必須版控才有歷史）。「不碰、留用戶」的精確意思 = **AI（尤其並行 session）不把 shared-state sweep 進自己的主題 commit**，由**檔案的擁有 session / 用戶**自行決定何時 commit 它。所以「`git log` 看到 CURRENT_FOCUS 常被 commit」是**正常**（擁有者在做事），**不是違規**；違規是「我這條 arc 的主題 commit 裡夾帶了我沒在改的 shared-state」。無機械 gate（git_branch_commit_guard 只擋 trunk）——靠此紀律 + 只 `git add <列檔>`。

### D. 何時「矯正 / 合併」

| 情境 | 動作 |
|------|------|
| **FF 可達**（`git merge-base --is-ancestor target feature` 為真）| `git branch -f target feature` → `git checkout target` → `git branch -d feature`（不動工作樹未提交變動）|
| **有 divergence** | `git merge` 產 merge commit；**衝突先停下討論** |
| **跨 session 污染**（commit 落錯 branch）| **不可在另一 session live 時 reset/rebase**；等持續安靜(others==0 + tip 跨區間不變 + 無 index.lock)→ 隔離 worktree cherry-pick stray 回原 branch → `git reset --mixed` 還原自己 branch（2026-06-10 實證流程）|
| **push / merge 回主線** | **永遠需用戶確認**，不自動 |

---

## §4 鉤子與機制（把 §3 從「建議」變「機械」）

### 既有（已落地）
| 機制 | 角色 |
|------|------|
| `pre_commit_compile_check`/`kb_schema_check`/`pipeline_block_check`/`no_binary_commit`/`search_scope_guard`/`pre_tier_upgrade_check`/`number_provenance_check` | 6+1 真 exit-2 Hard Gate（L4 驗證不可繞過）|
| `health_drift_advisor`（C7）| harness 變動 → nudge `/harness-health`（修 invocation-dependence 漂移）|
| `concurrent_session_advisor`（§G）| 主 repo >1 活躍 session → nudge 開 worktree（修跨 session 互撞）|
| `harness_health.py` 9 燈（含 Light#9 hook-wiring）| 監看「宣稱 vs 真實」漂移 |
| `precompact_autosave` | 壓縮前 milestone checkpoint 落檔 |

### 2026-06-11 新增（本檔落地）
| 機制 | 角色 |
|------|------|
| **③ `git_branch_commit_guard.sh`**（PreToolUse Bash，commit 前）| **exit-2 擋直接 commit 到 `main`/`master`/`develop`**（受保護主線；本 repo 一律 feature branch）+ **advisory warn 若偵測並行 session**。fail-OPEN（任何錯 → exit 0，絕不擋所有 commit）|
| **⑤ `.worktreeinclude`** + worktree-per-session 規則 | 新 session 用 worktree 時把 gitignored 的本地 config/data 指標帶進去，使隔離 session 真能跑（非空殼）。並把 §3-B 升為**首選原生 worktree**（啟用前確認 `claude --worktree` 旗標於當前版本）|
| **§8 concise-emit 一行** | return-contract 綁 output-token(500) 上限：單次回覆勿超量 → 寫檔回 path（修 friction #1 整 session 丟失）|

### 待用戶決定是否升級
- `git_branch_commit_guard` 的並行-session 部分目前 **advisory**（無法判斷意圖）；若要升 exit-2（並行時硬擋 commit）需用戶確認——風險是誤擋正當 worktree commit。

---

## §5 L4 驗證層 — 「done = 可執行」+ token 韌性

- **完成宣稱 gate（§13.7）**：宣稱 done/pass/修好前跑 IDENTIFY→RUN→READ→VERIFY→CLAIM；紅旗語言（should/probably/應該）= 停下驗證。
- **數據誠信（§13）**：報告每數字必能 grep 到來源；產數字與寫報告**永不同批**；`number_provenance` gate（已 wired）。
- **token / rate-limit 韌性**（friction #1 對策，2026-06-11）：① fan-out **cap 4-6**（≥5 樣本改循序）；② 長 job/多樣本**每單元落 disk checkpoint** → resume 只補缺；③ 單回覆 concise-emit 寫檔回 path；④ Dynamic Workflow resumable（resumeFromRunId）。

---

## §6 本次（2026-06-11）落地與待辦

- ✅ 本主工作流 SoT（本檔）。
- 套用中（others==0 窗口）：§8 concise-emit、`git_branch_commit_guard.sh` + wiring、`.worktreeinclude`、cross_sample_benchmark fan-out cap + checkpoint。
- 🔴 不做（restraint，web 研究確認過度工程）：DVC/git-annex 全裝、RO-Crate、Cookiecutter 重構、pre-commit framework 整套。
- 📌 最大「檔案整理」槓桿 = **commit 未 commit 的 script bucket**（318 untracked 中 169 是源碼）= §3-C 紀律問題，非新工具。

> **使用方式**：之後任何任務，先回到 §1 那張圖 + §2 五步；動 git 前查 §3 四表；不確定鉤子擋什麼看 §4。本檔有變動 → 同步 git governance 細節檔。
