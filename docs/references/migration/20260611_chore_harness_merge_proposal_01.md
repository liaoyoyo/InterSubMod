---
title: chore/harness-loop-governance-202606 — 5-commit Merge 提案
date: 2026-06-11
status: proposal / awaiting-user-decision
type: governance / merge-proposal
owner: liaotzuyu000@gmail.com
---

# Merge 提案：harness 治理 5 commit 如何整合

> **決策需用戶**：選 target + push/merge 一律需確認（git governance §3-D）。本檔只提案 + 給指令，不自動執行。

## §1 要整合什麼（5 commit，乾淨自包含的 11 檔，+722/−3）

| commit | 內容 |
|--------|------|
| `691a60b` | C7 `health_drift_advisor` + loop-engineering ADR |
| `12c5d54` | git 治理 §F 強制建議 commit 時機 + §G 衝突/聚焦 |
| `3101051` | `concurrent_session_advisor` + §G 並行 worktree 規則 |
| `ca1c5ba` | wire §13 `number_provenance` gate + 2 SKILL.md broken-ref + harness_health Light#9 |
| `c64be2b` | 主工作流架構 SoT + `git_branch_commit_guard` + worktree launcher + §8 concise-emit |

**觸碰檔（全屬 harness/docs/scripts，無 cis-asm 研究檔 → 邏輯自包含）**：
`.claude/CLAUDE.md`、2 個 SKILL.md、`docs/plans/20260609_loop_eng…`、`docs/references/20260611_master_workflow…`、`docs/references/migration/20260601_git_branch…`、`scripts/harness_health.py`、4 個 `scripts/hooks/*.sh` + launcher。

## §2 拓樸問題（為何不能直接 merge）

```
main e987fb6 ──┐ (三線皆分歧；main ⊄ develop)
develop 8d7cdc0 ┤
                │   065d9bc ─┬─ (cis-asm 線) ──→ feat/cis-asm-pipeline e0afc38 (065d9bc + 7 cherry-pick)
                │            └─ (harness 線) ──→ chore … c64be2b (065d9bc + 我的 5 commit)
```

- **chore base = `065d9bc`，NOT 在 main/develop 上**（chore 距 main **414 commit**，多為 cis-asm 研究）。→ **直接 `merge chore→main/develop` 會把 414 個 cis-asm commit 一起帶過去 = 錯。**
- **harness config 已跨線分歧**：`develop` 的 `.claude/CLAUDE.md` 對 chore-base 差 **284 ins / 205 del**；`scripts/harness_health.py` 在 develop 視角是 **686 行全新**（develop 沒這版）。→ **cherry-pick 5 commit 到 develop/main 會重度衝突**（trunk 的 harness 檔太舊/太不同）。
- 根因：harness 主要在**研究 feature 線上演進**，trunk（develop/main）相對 stale。

## §3 三個整合選項

| 選項 | 動作 | 衝突 | 取捨 |
|------|------|:---:|------|
| **A（推薦）merge chore → `feat/cis-asm-pipeline`** | 把 5 harness commit 併回它所基於的活躍研究線（同 harness-檔 lineage，cis-asm 的 7 commit 沒碰我的 harness 檔 → 低衝突；diverged 故產 merge commit） | 低 | harness 治理跟著活躍線走；之後 cis-asm→trunk 時一起上 trunk |
| B cherry-pick 5 → `develop`/`main` | 隔離出純 harness 整合到 trunk | **高**（CLAUDE.md 284/205、harness_health 全新）| 現在做要逐檔解大量衝突；不划算 → **建議延到 cis-asm 整體上 trunk 時一起處理** |
| C 維持 chore 獨立 | 不動，harness 工作留 branch | 無 | live 可用但未整合；之後再決定 |

**推薦 A**：5 commit 回到 cis-asm 活躍線（它們本就基於此線、harness 檔版本相符），低衝突、harness 治理與研究線同步前進。trunk 整合（B）是更大的「cis-asm 線 → trunk」工程的一部分，不適合現在硬拆。

## §4 推薦指令（選項 A；用戶確認後執行，建議並行 session 安靜時）

```bash
# 在 isolated worktree 做，避免動主 dir 共用 HEAD（§3-B）
git worktree add /tmp/merge_chore feat/cis-asm-pipeline
git -C /tmp/merge_chore merge --no-ff chore/harness-loop-governance-202606 \
    -m "merge(harness): C7+§F§G+concurrent-advisor+number_provenance+Light#9+主工作流SoT 治理 5 commit"
# 衝突？git -C /tmp/merge_chore merge --abort + 回報。乾淨 → 驗證:
git -C /tmp/merge_chore log --oneline -8
python3 scripts/harness_health.py   # 應 8 GREEN + 1 良性 YELLOW
git worktree remove /tmp/merge_chore
# 成功後 chore branch 可刪：git branch -d chore/harness-loop-governance-202606
```
> push 到遠端 / 進一步 merge 到 trunk → **另需用戶明示確認**。

## §5 風險與前置

- merge 前 SAFETY GATE：others==0 + 無 index.lock（同 untangle 紀律）。
- `settings.local.json` 的 hook wiring（number_provenance + git_branch_commit_guard）**仍 live-but-uncommitted**，依 §C 留用戶；本 merge 不含它（merge 後若要入版控，由用戶 commit settings）。
- merge 不可逆性低（merge commit 可 reset 回退）；但 push 後不可逆 → push 前確認。
