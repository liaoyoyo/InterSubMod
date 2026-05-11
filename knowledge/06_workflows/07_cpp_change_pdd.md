---
id: ism-kb-06-workflows-cpp-change-pdd
name: "C++ 修改 PDD 6 步驟協議"
description: "任何 C++ 修改前必經的 6 步驟協議：BASELINE → IMPLEMENT → UNIT TEST → FEATURE → VALIDATE → EVIDENCE LEDGER。前置 methodology-audit；commit 可追溯可回退。Hard rule，不因模型升級鬆動。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: howto
verified_scope: "PDD protocol against .claude/skills/cpp-change/SKILL.md + CLAUDE.md hard rules"
related_ids:
  - ism-kb-06-workflows-index
  - ism-kb-06-workflows-build-and-test
  - ism-kb-00-governance-query-protocol
  - ism-kb-10-research-status-evidence-ledger-format
  - ism-kb-00-governance-hooks-and-automation
  - ism-kb-00-governance-confirmation-protocol
  - ism-kb-00-governance-think-before-code
tags: [workflow, cpp, pdd, protocol, 6-steps, hard-rule]
canonical_paths: [06_workflows/07_cpp_change_pdd.md]
alias_paths: []
---

# C++ 修改 PDD 6 步驟協議

- 一句結論：任何 `src/*.cpp` / `include/*.hpp` 修改必經 6 commits（baseline → implement → test → feature → validate → evidence）；前置 methodology-audit；Hard rule
- 適用對象：修改 C++ 前、code review 檢查點
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 查完整 skill 定義
  cat /big7_disk/liaoyoyo2001/InterSubMod/.claude/skills/cpp-change/SKILL.md
  ```

---

## ⚠️ Hard Rule 聲明

本協議為**專案硬性要求**，`.claude/CLAUDE.md` 明列：
> 「C++ 修改的 6 步驟 PDD 協議（`/cpp-change`）不因模型升級而鬆動」

**違反後果**：
- `git commit` hook 會在 C++ 未編譯時 `exit 2` 阻擋（Hard Gate）
- 缺 evidence_ledger 紀錄的修改無法追溯成因
- 跳過 methodology-audit 可能導致方法學錯誤無人發現

---

## 🔀 前置：methodology-audit（若適用）

**觸發條件**：
- 用戶說「審查某方法」、「評估是否要改 XX」
- 即將修改任何 `.cpp` / `.hpp` 前
- 發現統計失效率過高、分類邊界模糊、特徵設計可疑時

**產出**：`docs/methodology/YYYYMMDD_{主題}_{流水號}.md`，含：
- 問題量化（失效率、AUROC、分佈）
- ≥3 個方案（含「不修改」選項）
- 驗收標準
- 用戶決策欄位（A/B/C）

**決策確認**後 → 進入 6 步驟實作協議

---

## 📋 6 步驟總覽

```
┌──────────────┐   ┌──────────────┐   ┌──────────────┐
│ Step 1       │   │ Step 2       │   │ Step 3       │
│ BASELINE     │ → │ IMPLEMENT    │ → │ UNIT TEST    │
│ chore:       │   │ (prep)       │   │ test:        │
└──────────────┘   └──────────────┘   └──────────────┘
                                              │
┌──────────────┐   ┌──────────────┐           ▼
│ Step 6       │   │ Step 5       │   ┌──────────────┐
│ EVIDENCE     │ ← │ VALIDATE     │ ← │ Step 4       │
│ chore:       │   │ docs:        │   │ FEATURE      │
│ ledger       │   │              │   │ feat:/fix:/  │
└──────────────┘   └──────────────┘   │ refactor:    │
                                      └──────────────┘
                                              │
                                        Step 4.5:
                                        CODE REVIEW
                                        (pr-review-toolkit)
```

---

## Step 1 — BASELINE（基線快照）

**目的**：記錄修改前狀態，失敗時可 `git revert` 回此點

```bash
# 確認當前 build 通過
cd build && make -j$(nproc) 2>&1 | tail -5

# /test-quick 確認基線可用
./scripts/run_vcf_all_snv.sh --mode chr19-verification

# 記錄基線 F1
python3 scripts/analysis/collect_baseline_metrics.py
```

**Commit**：
```
chore: snapshot baseline before [問題名稱] change

[基線 F1 摘要，如：HCC1395-5kHz F1=0.7167]
```

---

## Step 2 — IMPLEMENT（最小可測試修改）

**原則**：
- 一次只改一個變數/閾值/邏輯分支
- 先改最高影響的那一個
- `clang-format -i` 格式化後再 commit

```bash
clang-format -i src/core/XX.cpp

# 背景編譯（不阻塞前台準備 unit test）
cd build && make -j$(nproc) 2>&1 | tail -20  # run_in_background=True
```

**不使用背景編譯的情況**（串行更安全）：
- 修改 CMakeLists.txt
- 修改 header 檔（影響範圍大）
- 第一次修改該檔案

**注意**：Step 2 本身不 commit；與 Step 4 合併提交

---

## Step 3 — UNIT TEST（單元測試提交）

```bash
cd build && ./tests/test_<name> --gtest_filter="*<TestCase>*"

# 必要時新增測試案例
# tests/test_XX.cpp
```

**Commit**：
```
test: [問題名稱] unit test coverage

- 新增/修改測試: tests/test_XX.cpp
- 覆蓋邊界條件: [說明]
```

---

## Step 4 — FEATURE（功能提交）

```bash
cd build && make -j$(nproc)
ctest --output-on-failure
```

**Commit 類型**：

| 類型 | 場景 | 範例 |
|------|------|------|
| `feat:` | 新增功能/邏輯 | `feat: PERMANOVA fallback to subset when NaN>30%` |
| `fix:` | Bug 修正 | `fix: LOH flag propagated to VerificationClass` |
| `refactor:` | 重構無行為改變 | `refactor: extract quality_score penalty table` |

---

## Step 4.5 — CODE REVIEW（自動審查）

**目的**：攔截邏輯錯誤、thread safety、memory leak、OWASP

```
Agent(
  prompt="Review the C++ changes in the latest commit against project conventions:
    - Check: OWASP top-10, thread safety (OpenMP), memory leaks (HTSlib)
    - Check: clang-format compliance, naming conventions
    - Focus on: src/core/ changes only
    - Report: HIGH confidence issues only",
  subagent_type="pr-review-toolkit:code-reviewer"
)
```

**通過標準**：無 HIGH confidence issue → 進 Step 5；有則修復回 Step 4

**可選追加**（大改動時）：
- `pr-review-toolkit:silent-failure-hunter`（錯誤處理品質）
- `pr-review-toolkit:type-design-analyzer`（新增 struct/class 時）

---

## Step 5 — VALIDATE（驗證結果）

### 5a. 串行驗證（預設，小改動）
```bash
./scripts/run_vcf_all_snv.sh --mode chr19-verification    # <30s
# /test-data                                              # ~1min
# /test-full                                              # ~5-10min（可選）
```

### 5b. 平行驗證（大改動或跨樣本）
```
Agent(
  prompt="驗證 C++ 修改對多樣本的影響:
    TASKS:
    1. scripts/run_batch_vcf_analysis.sh --mode HCC1395_5kHz
    2. scripts/run_batch_vcf_analysis.sh --mode HCC1395_DORADO
    3. scripts/run_batch_vcf_analysis.sh --mode COLO829
    OUTPUT_DIR: research/autoresearch/cycles/cpp_change_validation/",
  subagent_type="parallel-benchmark"
)
```

### 5c. F1 delta 計算
```bash
python3 scripts/analysis/collect_baseline_metrics.py
```

**Commit**：
```
docs: [問題名稱] validation result

delta_f1:
  HCC1395-5kHz:  [+/-0.XXXX]
  COLO829:       [+/-0.XXXX]
  [其他已測樣本]

結論: [通過/回退/無顯著影響]
```

---

## Step 6 — EVIDENCE LEDGER（研究記錄）

**必做**：追加至 `research/autoresearch/evidence_ledger.jsonl`

```json
{
  "cycle_id": "YYYY-MM-DD_[問題名稱]",
  "hypothesis_id": "H0XX",
  "hypothesis": "[一句話說明修改了什麼]",
  "pipeline_track": "TO" | "paired" | "both",
  "type": "cpp_improvement",
  "tier": 3,
  "delta_f1": {
    "HCC1395-5kHz": 0.XXXX,
    "COLO829": 0.XXXX
  },
  "human_decision": "keep" | "revert",
  "key_observations": "[關鍵觀察]",
  "methodology_doc": "docs/methodology/YYYYMMDD_*.md",
  "operator": "[AI agent ID 或研究者名；v0.5+ 必填]",
  "reviewer": "[review verdict 的人；v0.5+ 必填]"
}
```

**v0.5 新增**：`operator`（誰執行 cycle）+ `reviewer`（誰 review verdict），對齊 ELN accountability。詳 `10_research_status/03_evidence_ledger_format.md`。

**Commit**：
```
chore: evidence_ledger record [H_ID] cpp_improvement [問題名稱]
```

---

## 🔄 回退檢查點

若 Step 5 驗證顯示 **F1 大幅下降（delta < -0.002）**：

1. **不要繼續 Step 6**
2. 更新 `docs/methodology/YYYYMMDD_*.md` 狀態為 `validated_failed`
3. `git revert` 回 Step 1 的 baseline commit
4. 回到 `methodology-audit` 討論其他方案

---

## Commit 類型完整定義

| Commit 類型 | 用途 | 強制場景 |
|-------------|------|----------|
| `chore: snapshot baseline` | 修改前快照 | Step 1 必用 |
| `test:` | 單元測試修改 | Step 3 必用 |
| `feat:` | 新功能/新邏輯 | Step 4，新增行為時 |
| `fix:` | Bug 修正 | Step 4，修正錯誤時 |
| `refactor:` | 重構無行為改變 | Step 4，只改結構時 |
| `docs:` | 驗證結果記錄 | Step 5 必用 |
| `chore: evidence_ledger` | 研究記錄 | Step 6 必用 |

---

## 🔗 相關協議

### Hard Gate 自動檢查
`.claude/settings.local.json` 設定 hook：
- `PreToolUse / git commit`：C++ 未編譯 → `exit 2` 阻擋
- `PostToolUse / Edit .cpp/.hpp`：建立待編譯標記
- `PostToolUse / make`：清除待編譯標記

### 前置 Skill
- `/methodology-audit` — 審查現有方法、量化問題、產出決策文件
- `/confirmation-protocol` — 確認時機判斷（Hard Gate vs Gate vs Review vs FYI）

### 後續 Skill
- `/test-quick`、`/test-data`、`/test-full`、`/validate` — 驗證工具
- `parallel-benchmark` subagent — 跨樣本平行驗證
- `pr-review-toolkit:code-reviewer` — Step 4.5 自動審查

---

## ✅ 最短 6 步驟 git log 範例

```
* chore: evidence_ledger record H017 cpp_improvement PERMANOVA fallback
* docs: PERMANOVA fallback validation result
|   delta_f1:
|     HCC1395-5kHz: +0.0023
|     COLO829:      -0.0003
|   結論: 通過
* feat: PERMANOVA fallback to subset when NaN>30%
* test: PERMANOVA fallback unit test coverage
* chore: snapshot baseline before PERMANOVA fallback change
|   HCC1395-5kHz F1=0.7167
```

---

## 🧭 何時該跑此協議

- ✅ 修改 `src/core/*.cpp`、`include/core/*.hpp`
- ✅ 修改 `include/utils/ArgParser.hpp`（CLI 行為變更）
- ✅ 修改 `include/core/Config.hpp` 預設值或常數
- ❌ **不需**：修改 `scripts/**/*.py`（Python 分析腳本走 evidence_ledger 但不需 6 steps）
- ❌ **不需**：修改 `docs/**/*.md`（文件變更走 `/doc-standards`）
- ❌ **不需**：修改 `tests/**/*.cpp`（測試新增本身是 Step 3，不獨立觸發）

---

## 相關

- Build/test workflow：[01_build_and_test.md](01_build_and_test.md)
- Evidence ledger 格式：[../10_research_status/03_evidence_ledger_format.md](../10_research_status/03_evidence_ledger_format.md)
- Query protocol：[../00_governance/05_query_protocol.md](../00_governance/05_query_protocol.md)
- Skill 定義：`.claude/skills/cpp-change/SKILL.md`、`.claude/skills/methodology-audit/SKILL.md`
- CLAUDE.md Hard rule 聲明：`.claude/CLAUDE.md` 「不可放寬的硬性規則」章節
