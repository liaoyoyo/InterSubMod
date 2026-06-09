---
title: Git 分支 / Commit 治理規則
date: 2026-06-01
type: governance / reference
status: active
---

# Git 治理 — 何時開分支、何時持續 commit、如何分類、未來怎麼管

> **目的**：固定一套 git 決策規則，避免「主題混進錯 branch / pre-existing 混雜檔被誤 commit / commit 顆粒不清」。
> **背景**：本 repo 工作樹長期有大量 pre-existing 未 commit 研究產物（60+ ??）+ 多個 M 檔（CURRENT_FOCUS / queue / ledger / log），commit 時必須精準只動「自己這次的檔」。

---

## A. 開新 branch vs 持續在當前 branch commit — 決策表

| 情境 | 動作 | 範例 |
|------|------|------|
| **跨多 commit 的獨立主題 / feature / 重構** | **開新 branch** | 本次 harness 數據誠信（3 commits）→ `chore/harness-data-integrity-2026` |
| **工作主題 ≠ 當前 branch 主題** | **開新 branch** | 在 `refactor/phase1-safety`(C++) 上做 harness/doc → 該開 `chore/` 而非混入 |
| **實驗性 / 可能放棄的方向** | **開新 branch**（隔離，不污染主線）| 探索性 pilot、可能 revert 的大改 |
| **需獨立 review / 未來可能 PR 的單元** | **開新 branch** | 對外 handoff、跨人協作 |
| **同主題的連續推進（增量）** | **持續 commit 當前 branch** | phase1-safety 上的多個 C++ safety 修正 |
| **小修 / typo / 跟當前 branch 一致** | **持續 commit 當前 branch** | 修當前 branch 的文件 typo |
| **當前 branch 的 hotfix** | **持續 commit 當前 branch** | 修剛引入的 bug |

**一句話判準**：*「這批工作的主題，跟我現在站的 branch 是同一件事嗎？」* 同 → 持續 commit；不同 / 獨立 / 可能放棄 → 開新 branch。

**Base 選擇**：延續既有工作 → 從**當前 HEAD** 開；全新獨立 → 從 **main / develop** 開。

---

## B. Commit 分類原則（顆粒度）

1. **一個 commit = 一個邏輯變更**（atomic）；可獨立描述、獨立 revert。
2. **多主題 → 拆多 commit**（如本次：Batch A 防捏造系統 / Batch B 儀表板治理 / Batch C 鐵則強化）。
3. **只 stage 自己這次的檔**：`git add <明確列檔>`，**禁用 `git add -A` / `git add .`**（會掃進 pre-existing 混雜檔）。
4. **commit message 格式**：`type(scope): 摘要` + body（what/why + 驗收）+ `Co-Authored-By:` 行。
   - `type`：`feat` / `fix` / `docs` / `chore` / `refactor` / `test`。
5. **commit 前必 `git status` 分類**：明確區分「我這次改的」vs「session 開頭就在的 pre-existing」，後者**不碰**。

---

## C. 這 repo 特定治理（非通用，務必遵守）

| 規則 | 細節 |
|------|------|
| **pre-existing 混雜檔不動** | `CURRENT_FOCUS.md` / `hypothesis_queue.json` / `research_direction.md` / `active.json` / `evidence_ledger.jsonl` / `*.log` + 歷年 `??` 研究產物 → 除非任務明確要改，否則**留給用戶**，不混入 feature commit |
| **Hard Gate（不可繞過）** | C++ commit 必過 `pre_commit_compile_check`（須無 stale `/tmp/ism_cpp_pending_compile.txt`）；`no_binary_commit` 擋 binary；`kb_schema_check` |
| **敏感檔走 skill** | `state.json`/`active.json`（→ cycle skill）、`hypothesis_queue.json`/`evidence_ledger.jsonl`（→ skill + commit 標 `[manual-queue-edit]`）；**禁手改** |
| **生成物不入 git** | `state/health_snapshots/`（dashboard）等已 gitignore；勿強加 |
| **push / merge 需用戶確認** | `git push`（ask-gated）、merge 回主線 → **等用戶明示**，不自動 |
| **memory 在 repo 外** | `/bip7_disk/.../memory/*.md` 不入 git commit（獨立持久層）|
| **數字誠信（§13）** | 含 metric 的 .md/.html commit 前過 `number_provenance_check` gate（validated/pi 強制）|

---

## D. 標準工作流（branch → 分類 → commit → merge）

```bash
# 1. 開 branch（主題獨立時）
git checkout -b <type>/<topic>-<scope>           # 從正確 base

# 2. 分類：先看全狀態，分「我的」vs「pre-existing」
git status --porcelain=v1

# 3. 只 stage 我這次的檔（明確列，不用 -A）
git add path/a path/b ...
git diff --cached --name-only                    # 驗證：只有我的檔，無洩漏

# 4. commit（type(scope) + body + co-author）
git commit -F <msg-file>

# 5. merge 回主線（用戶確認後）
#    FF 情境（target 是 ancestor，無 divergence）— 不動工作樹未提交變動：
git branch -f <target> <feature>                 # 移 target ref 到 feature tip
git checkout <target>                            # 切 HEAD（target==current commit → 工作樹不變）
git branch -d <feature>                          # 刪已併入的 feature branch
#    非 FF（有 divergence）→ git merge 產 merge commit；衝突先停下討論
```

**FF vs merge commit 判斷**：`git merge-base --is-ancestor <target> <feature>` 為真 = 純快進（無 merge commit）；否則需 merge commit（先確認無衝突）。

---

## E. 本次（2026-06-01）作為範例

- 主題（harness 數據誠信）≠ 當前 branch（`refactor/phase1-safety` 是 C++ safety）→ **開 `chore/harness-data-integrity-2026`** ✓
- 3 個邏輯變更 → **3 commits**（防捏造系統 / 儀表板治理 / 鐵則強化）✓
- 只 stage 我的 9 檔，60+ pre-existing 全未動 ✓
- 用戶確認後 **FF merge 回 `refactor/phase1-safety`** + 刪 chore 分支 ✓

---

## F. 何時「強制建議 commit」— 觸發時機（2026-06-09 補；advisory 非 Hard Gate）

> 用戶 2026-06-09 要求補足「哪些時候強制建議 commit」。這是 advisory（提醒，不阻擋）——目的是**縮小未 commit 工作樹、減少衝突、保持聚焦**。判準一律看「**我這次新增/改的檔**」，不看總 dirty（本 repo 長期 200-300 dirty 多為 pre-existing）。

| # | 觸發時機 | 為什麼 |
|---|---------|--------|
| **F1** | **一個已驗證的邏輯單元完成** | 完成 + verified（例：C7 hook 測過、harness_health 8 GREEN）→ 立即 commit，別累積成大堆 |
| **F2** | **即將切換主題 / 開新 branch 前** | 先把當前主題 commit 乾淨，新 branch 才不會帶著舊主題半成品 → 反「主題混入錯 branch」 |
| **F3** | **即將做不可逆 / 高風險 git 操作前** | `git checkout`(切 branch) / `reset` / `rebase` / `merge` / 大刪除 前先 commit 保護當前工作 |
| **F4** | **「我的檔」累積 > ~10 changed 或跨 > 3 邏輯單元** | 未 commit 越多越難精準 stage（混雜風險↑）、越難回溯、conflict 面↑ → 分批 commit |
| **F5** | **C++ 改動完成並編譯通過後** | compile Hard Gate 一過即 commit，別讓編譯狀態與原始碼漂移 |
| **F6** | **session 結束 / 跨 session 交接前** | Stop hook 已提醒寫報告；同時把完成單元 commit，交接才乾淨 |

**不該為了湊 commit 而 commit**：半成品 / 未驗證 / 跑到一半的分析**不 commit**（對齊 §13 完成宣稱需 fresh 驗證）；留 `{{待填}}` 或續做完再 commit。

---

## G. 衝突最小化 + 聚焦機制（2026-06-09 補）

> 用戶要求「減少任務與修改衝突，同時更聚焦」。核心矛盾源 = **shared live-state 檔被多任務同時改** + **未 commit 工作樹無限膨脹**。

| 機制 | 做法 |
|------|------|
| **shared live-state 檔不入 feature commit** | `CURRENT_FOCUS.md` / `hypothesis_queue.json` / `evidence_ledger.jsonl` / `active.json` / `*.log` / `state/*snapshots/` 是**跨任務共享**，每個任務都改它 = 衝突主源 → 留用戶 / 走 skill，feature commit **不碰**（同 §C）|
| **完成即 commit（F1）縮小工作樹** | 下個任務從乾淨狀態開始 = 少衝突 + 主題單一 = 聚焦 |
| **一 branch 一主題** | branch 名即主題宣告；混入別主題 = 反模式（同 §A）|
| **平行第二線 / 大 refactor → git worktree 隔離** | 兩線改同檔不互覆（本 harness `parallel-benchmark` / `Workflow isolation:'worktree'` 已支援；Loop-Engineering 第 2 要素）|
| **WIP-limit** | active 主題線 ≤ 2-3（對齊 `active.json` ≤5 cycle + `CURRENT_FOCUS` ≤2 background slot）；超過先收斂再開新 → 聚焦 |

**3 句總決策流**（任何要動 git 時自問）：
1. **這批工作跟我站的 branch 同主題嗎？** 同 → 繼續 commit；不同/獨立/可能放棄/C++/需 review → **開新 branch**（§A）。
2. **現在該 commit 了嗎？** 命中 F1-F6 任一 → **commit**（只 stage 自己檔，§B/§C）。
3. **會不會跟別任務撞？** shared live-state 不碰、平行線用 worktree、WIP ≤ 2-3（§G）。
