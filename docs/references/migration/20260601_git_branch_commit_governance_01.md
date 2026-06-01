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
