---
name: verification-loop
description: "**程式碼級綜合驗證迴圈**（不同於 /validate command 的研究假說 benchmark）— build、type check、lint、測試、安全掃描、diff review。USE WHEN：「驗證程式碼」「verify」「check quality」「PR 前 review」「程式碼品質」、建立 PR 前。涉及 .cpp/.hpp/.py 原始碼、build/ 產出。**職責邊界**：本 skill = 程式碼語法/build/test 級驗證；`/validate` command = 研究假說 benchmark + experiment_report；兩者互補不重疊。SKIP WHEN 研究假說 benchmark（用 /validate command）、tier 升級判定（用 run-evaluator）、純 docs / 純 .md 編輯、已通過 CI 不需重 verify。"
version: 1.1.0
user-invocable: true
---

# Verification Loop Skill

A comprehensive verification system for Claude Code sessions.

## When to Use

Invoke this skill:
- After completing a feature or significant code change
- Before creating a PR
- When you want to ensure quality gates pass
- After refactoring

## Verification Phases

Choose the commands adaptively for the current project instead of running every example blindly. Use the stack-appropriate command from `references/STACK-DETECTION.md` when the repo does not match the default examples below.


### Phase 1: Build Verification
```bash
# Python projects (uv)
uv build 2>&1 | tail -20
# OR
python -m build 2>&1 | tail -20

# Node.js projects
npm run build 2>&1 | tail -20
# OR
pnpm build 2>&1 | tail -20
```

If build fails, STOP and fix before continuing.

### Phase 2: Type Check
```bash
# TypeScript projects
npx tsc --noEmit 2>&1 | head -30

# Python projects
pyright . 2>&1 | head -30
```

Report all type errors. Fix critical ones before continuing.

### Phase 3: Lint Check
```bash
# JavaScript/TypeScript
npm run lint 2>&1 | head -30

# Python
ruff check . 2>&1 | head -30
```

### Phase 4: Test Suite
```bash
# Python projects
pytest --cov=src --cov-report=term-missing 2>&1 | tail -50

# Node.js projects
npm run test -- --coverage 2>&1 | tail -50
```

Report:
- Total tests: X
- Passed: X
- Failed: X
- Coverage: X%

### Phase 5: Security Scan
```bash
# Python: Check for secrets
grep -rn "sk-" --include="*.py" . 2>/dev/null | head -10
grep -rn "api_key" --include="*.py" . 2>/dev/null | head -10
pip-audit

# Node.js: Check for secrets
grep -rn "sk-" --include="*.ts" --include="*.js" . 2>/dev/null | head -10
grep -rn "api_key" --include="*.ts" --include="*.js" . 2>/dev/null | head -10

# Check for debug statements
grep -rn "print(" --include="*.py" src/ 2>/dev/null | head -10
grep -rn "console.log" --include="*.ts" --include="*.tsx" src/ 2>/dev/null | head -10
```

### Phase 6: Diff Review
```bash
# Show what changed
git diff --stat
git diff HEAD~1 --name-only
```

Review each changed file for:
- Unintended changes
- Missing error handling
- Potential edge cases

## Output Format

After running all phases, produce a verification report:

```
VERIFICATION REPORT
==================

Build:     [PASS/FAIL]
Types:     [PASS/FAIL] (X errors)
Lint:      [PASS/FAIL] (X warnings)
Tests:     [PASS/FAIL] (X/Y passed, Z% coverage)
Security:  [PASS/FAIL] (X issues)
Diff:      [X files changed]

Overall:   [READY/NOT READY] for PR

Issues to Fix:
1. ...
2. ...
```

## Continuous Mode

For long sessions, run verification every 15 minutes or after major changes:

```markdown
Set a mental checkpoint:
- After completing each function
- After finishing a component
- Before moving to next task

Run: /verify
```

## Integration with Hooks

This skill complements PostToolUse hooks but provides deeper verification.
Hooks catch issues immediately; this skill provides comprehensive review.


## Reference Files

Load only what is needed:
- `references/STACK-DETECTION.md` - how to choose the right verification command set for the current repo
- `references/REPORT-TEMPLATE.md` - report structure for final verification output
- `examples/example-verification-report.md` - example final report

---

## Phase & Chain Position

- **Phase**: **Governance / Cross-phase**（程式碼變更後通用，非單 phase）
- **Chain**: forward-link chain #6 (P5 失敗回溯) 的最後一環
  ```
  /run-evaluator (P5 失敗，risk>0.7) → pivot-direction → review-evidence
      ↓
  methodology-audit (C++ 修改前審查)
      ↓
  cpp-change (PDD 6 步驟)
      ↓
  verification-loop ← (本 skill: build/lint/test/security/diff)
      ↓
  返回 P3 PILOT 重跑
  ```
- **與 /validate command 互補**：本 skill = 程式碼語法層；`/validate` = 研究假說 benchmark 層
- **上游觸發**: cpp-change 完成 / PR 前 / 用戶手動「verify」「check quality」
- **下游 skill**: 通過後 → 回到 P3 PILOT or 進 PR；失敗 → methodology-audit 再評估

## Dependencies

| 類別 | 項目 |
|---|---|
| **Uses** | Bash(make / cmake / ctest / clang-format / clang-tidy / pylint / mypy)、Read（git diff stat 與原始碼）、Grep（搜 TODO/FIXME） |
| **Used by** | cpp-change Step 6（PDD 最後驗證）/ /run-evaluator P5 失敗後的 fallback / PR pre-flight / 用戶手動「verify」 |
| **Reads** | `src/**/*.cpp` `include/**/*.hpp` `tests/**/*.cpp` `scripts/**/*.py`、build/ 產出、`.clang-format`、git diff |
| **Writes** | 不寫永久檔案；輸出 verification report 到 stdout（依 `references/REPORT-TEMPLATE.md` 格式） |

## Failure Mode & Diagnostics

| # | 失敗症狀 | 先看哪 | 排查步驟 |
|---|---|---|---|
| 1 | build fail (compile error) | `cd build && make 2>&1 \| head -50` | 修 syntax → 重 build；若 CMake config 錯，刪 build/ 重 cmake |
| 2 | unit test fail | `ctest --output-on-failure` 輸出 | 找對應 `tests/test_*.cpp`；個人風格 anchor #1 強驗證 — 不允許 skip 失敗 test |
| 3 | clang-format diff | `clang-format -i <files>` | 直接套用；commit 前 hook 已強制 |
| 4 | lint warning | clang-tidy 或 pylint 輸出 | 個別評估；嚴重 (security / undefined behavior) 必修，stylistic 可暫忽略 |
| 5 | diff review 看到無關改動 | `git diff --stat` | 個人風格 anchor #5「One-turn freeze」— 拆 commit 不混 scope |
| 6 | security scan 警告（如 hardcoded credential） | scan output | 立即 Hard Gate；不可 commit 含敏感資料 |

**何時升級到別的 skill / agent / 人工審查**：
- build fail 連續 3 次 → 升級 `methodology-audit` 重新評估改動方向
- test fail 涉及 statistics 邏輯（非 syntax）→ 跳到 `auc-confound-guard` / `known-pitfalls` 確認方法學
- diff 跨 ≥3 modules → 拆 commit 並走 PR review

**個人風格適配**（依 `feedback_*` memory）：
- **Anchor #1 「L4 多層驗證必建」** → test fail 不允許 skip；必須查根因
- **Anchor #5 「One-turn mechanism freeze」** → 一次完整 verify，不在迴圈中逐步 fix（先收集所有 issue 再批次解）
- **/cpp-change Step 6 PDD 規範** → 本 skill 是 PDD 最後一步，必過才能 commit C++ 修改

---

## 與 /scientific-rigor 元方法論的關係

本 skill 為 `/scientific-rigor §7.2 可重現性 checklist` 的**程式碼層級具體驗證**:
- Phase 1-6 對應 `/scientific-rigor §7.2 可重現性 7 項 checklist` 中與程式碼相關的子集（build / type / lint / test / security / diff）
- **場景分流明示**: 研究實驗 reproducibility（seed / data version / commit hash）→ `/check-staleness`；程式碼 build/test 驗證 → 本 skill
- `/scientific-rigor §11 協作圖 step 8` 直接呼叫本 skill 作為「程式碼級驗證」節點
- `/scientific-rigor §8.3 Reflexion buffer` 的「下次避免方法」常含本 skill 的失敗診斷症狀

**級聯觸發**: `/scientific-rigor §6 消融` 完成後 → `/cpp-change` PDD step 5-6 → 本 skill 6 phase → 通過才進 git commit Hard Gate
