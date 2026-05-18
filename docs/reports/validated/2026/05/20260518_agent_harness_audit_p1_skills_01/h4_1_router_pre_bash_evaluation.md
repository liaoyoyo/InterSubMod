<!--
建立時間: 2026-05-18
parent: p4_to_p7_final.standalone.html H4-1 fix item
purpose: P3 audit H4-1 follow-up — evaluate 5 PreToolUse Bash hooks 是否該整合 router
status: validated, integration NOT recommended
-->

# H4-1 Router PreToolUse Bash Evaluation — 5 Hooks 整合分析

> **2026-05-18 P3 audit H4-1 結論**: **不整合**。5 hooks 各對 specific bash command pattern，整合反增複雜度且降可讀性。維持現狀。

---

## 5 Hooks 職責矩陣

| Hook | 觸發 Bash command pattern | 動作 | Hard Gate? |
|------|-------------------------|------|-----------|
| Inline rm-rf guard | `rm -rf` / `git reset --hard` / `git push --force` | echo warning | ❌ Advisory only |
| `pre_commit_compile_check.sh` | `git commit*` 含未編譯 C++ | exit 2 阻擋 | ✅ Hard Gate |
| `kb_schema_check.sh` | `git commit*` 含 staged `knowledge/**/*.md` | exit 2 阻擋（3 schema validators 任一失敗） | ✅ Hard Gate |
| `pipeline_block_check.sh` | 長 pipeline（`run_batch_vcf` / `run_vcf_all_snv` 等）+ disk critical | exit 2 阻擋 | ✅ Hard Gate |
| `no_binary_commit.sh` | `git commit*` 含 staged binary (.pptx / .pdf / .jpg / .png / large .html) | exit 2 阻擋 | ✅ Hard Gate |

---

## Trigger Specificity 分析

| Hook | Trigger 粗細度 | 重疊風險 |
|------|--------------|---------|
| rm-rf guard | 廣 (3 個危險 pattern) | — |
| compile_check | `git commit*` only | 與 kb_schema / no_binary 都對 git commit fire |
| kb_schema_check | `git commit*` + KB staged | 同上 |
| pipeline_block | 長 pipeline scripts | 不與 git commit 重疊 |
| no_binary_commit | `git commit*` + binary staged | 同 compile/schema |

→ **3 hooks 對 `git commit*` fire**（compile / schema / binary）— 但 **檢查不同條件**：
- compile: C++ 是否編譯
- schema: KB schema 是否合規
- binary: 是否誤含 binary

**整合 router 設計**：可寫 `git_commit_pre_router.sh` 統一處理。但：
1. 各 check 邏輯獨立可重用
2. 整合 router 必須維護「3-in-1」邏輯複雜度
3. 任一 check fail (exit 2) → 整 router 必須回傳 exit 2，與目前順序執行相同
4. **觀察當前運作正常**（無 false positive / false negative 紀錄）

---

## Integration 候選評估

### Option A: 全整合為 `pre_bash_router.sh`
- ✅ 統一入口
- ❌ 邏輯混雜 5 種完全不同檢查（rm-rf / compile / schema / pipeline / binary）
- ❌ 任一改動需重測全部
- ❌ Hard Gate 邏輯與 advisory 邏輯混在一起

### Option B: 部分整合 `git_commit_pre_router.sh`（合併 3 個 git commit hooks）
- ✅ 對 git commit 場景統一 entrypoint
- 🟡 整合需處理「3 個 hook 順序」+ early-exit
- ⚠ 與既有 `pipeline_block_check` / inline rm-rf 不在同 router
- 🟡 收益邊際（每 git commit 只多 ~30ms latency）

### Option C: 維持現狀（推薦）
- ✅ Each hook 單一職責（SRP）
- ✅ 可獨立測試 / 改 / deprecate
- ✅ 設計符合 Anthropic「PreToolUse 多 hook 平行執行」spec
- ✅ Cache hit rate 94.6% 證明 hook 不影響 prefix cache
- 🟡 唯一缺點: 5 hooks fire on every bash → ~150-250ms 額外 latency（可接受）

---

## 決策

**選 Option C（維持現狀）** — 理由：
1. 5 hooks 各偵測不同條件，無功能重複
2. 整合複雜度 > 收益（latency 邊際 + 維護成本上升）
3. Anthropic spec 已支援多 hook 平行 — 設計 idiomatic
4. **若未來新增第 6+ git commit pre-check**，再評估部分整合（Option B）

---

## 後續監控（P4+ 候選）

- `hook_latency_log.sh`: 每 hook 跑完寫 latency 到 monthly log
- 若任一 hook latency > 500ms → 該 hook 需重構

---

## 元層說明

本 evaluation 源自 2026-05-18 P3 audit H4-1。原始 fix proposal 為「整合 5 個 PreToolUse Bash hooks → 單一 router_pre_bash.sh」，實際 audit 後 **不採納整合動作**，維持單一職責設計。

與 H4-2 (`InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/h4_2_hook_overlap_evaluation.md`) 結論一致：hook architecture 已優化良好，不必過度整合。
