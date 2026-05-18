---
name: cpp-change
description: PDD C++ 修改協議（6步驟）。在 methodology-audit 完成且用戶選定方案後啟動。確保每步驟有對應 commit，修改可追溯、可回退。USE WHEN 「開始實作 [審查文件名]」「執行方案 B」「**改 src/core/ 或 include/core/ C++ 邏輯**」「**修 ReadParser / NormalBaseline / Haplotag / KDE / pipeline C++ 模組**」「**加 C++ filter/feature/邏輯**」。若無對應 methodology-audit 報告 → 先推 /methodology-audit + /problem-framing-ideation 走完規格 + 方案選擇 + Step→Verify 後再進本 skill。DO NOT USE WHEN：純 Python/R 分析、文檔修改、改 .md 或 schema 檔。
user-invocable: true
paths: ["src/**/*.cpp", "src/**/*.hpp", "src/**/*.h", "include/**/*.hpp", "include/**/*.h", "tests/**/*.cpp", "CMakeLists.txt"]
---

# CPP-Change Skill — PDD C++ 修改協議

## 前置條件

- `methodology-audit` skill 已完成，存在對應的 `docs/methodology/YYYYMMDD_*.md`
- 用戶已選定方案（B 或 C，非 A）
- `build` 目前可正常編譯（先執行 `/build` 確認）

---

## 6 步驟執行協議

### Step 1：BASELINE COMMIT（基線快照）

```bash
# 確認當前 build 通過
cd build && make -j$(nproc) 2>&1 | tail -5

# 執行快速測試確認基線可用
# /test-quick

# 記錄基線 F1（使用現有 metrics，不重跑）
python3 scripts/analysis/collect_baseline_metrics.py
```

commit 格式：
```
chore: snapshot baseline before [問題名稱] change

[基線 F1 摘要，如：HCC1395-5kHz F1=0.7167]
```

---

### Step 2：IMPLEMENT（最小可測試修改）

```
原則：
- 一次只改一個變數/閾值/邏輯分支
- 先改最高影響的那一個
- clang-format 格式化後再確認
```

```bash
# 格式化修改的檔案
clang-format -i src/core/XX.cpp
```

#### 背景編譯 + 前台繼續

修改完成並格式化後，使用背景編譯以節省等待時間：

```
# 背景編譯（不阻塞前台）
Bash("cd /big7_disk/liaoyoyo2001/InterSubMod/build && make -j$(nproc) 2>&1 | tail -20",
     run_in_background=True)

# 前台同時可進行：
#   - 撰寫 unit test（Step 3 的測試碼）
#   - 準備 commit message
#   - 閱讀相關程式碼確認修改完整性
#   - 更新 methodology audit 文件
```

**編譯完成通知到達後**：
- 若成功 → 直接進入 Step 3 執行測試
- 若失敗 → 查看錯誤，修復後重新背景編譯

**不使用背景編譯的情況**（串行更安全）：
- 修改了 CMakeLists.txt（需確認 build 結構正確）
- 修改了 header 檔（影響範圍大，需先確認編譯）
- 第一次修改本檔案（需確認語法正確）

---

### Step 3：UNIT TEST COMMIT（單元測試）

```bash
# 執行對應單元測試
cd build && ./tests/test_<name> --gtest_filter="*<TestCase>*"

# 若需補充測試案例，先在 tests/ 撰寫
```

commit 格式：
```
test: [問題名稱] unit test coverage

- 新增/修改測試: tests/test_XX.cpp
- 覆蓋邊界條件: [說明]
```

---

### Step 4：FEATURE COMMIT（功能提交）

```bash
cd build && make -j$(nproc)
ctest --output-on-failure
```

commit 格式（依修改類型選擇）：

| 類型 | 場景 | 範例 |
|------|------|------|
| `feat:` | 新增功能/邏輯 | `feat: PERMANOVA fallback to subset when NaN>30%` |
| `fix:` | Bug 修正 | `fix: LOH flag propagated to VerificationClass` |
| `refactor:` | 重構無行為改變 | `refactor: extract quality_score penalty table` |

---

### Step 4.5：CODE REVIEW（自動程式碼審查）

在功能提交後、驗證前，啟動自動審查以攔截邏輯錯誤：

```
# 使用 pr-review-toolkit 的 code-reviewer agent
Agent(prompt="Review the C++ changes in the latest commit against project conventions:
  - Check: OWASP top-10, thread safety (OpenMP), memory leaks (HTSlib)
  - Check: clang-format compliance, naming conventions
  - Focus on: src/core/ changes only
  - Report: HIGH confidence issues only",
  subagent_type="pr-review-toolkit:code-reviewer")
```

**審查通過標準**：
- 無 HIGH confidence issue → 直接進 Step 5
- 有 HIGH issue → 修復後重新 commit（回到 Step 4）

**可選追加審查**（大型修改時）：

```
# Silent failure 檢查（錯誤處理品質）
Agent(prompt="Hunt for silent failures in latest C++ changes",
  subagent_type="pr-review-toolkit:silent-failure-hunter",
  run_in_background=True)

# 型別設計檢查（新增 struct/class 時）
Agent(prompt="Analyze type design of new types in latest changes",
  subagent_type="pr-review-toolkit:type-design-analyzer",
  run_in_background=True)
```

---

### Step 5：VALIDATE COMMIT（驗證結果）

#### 串行驗證（預設，小改動）

```bash
# 快速驗證（< 30 秒）
# /test-quick

# 數據驗證（約 1 分鐘）
# /test-data

# 完整驗證（約 5-10 分鐘，可選）
# /test-full
```

#### 平行驗證（大改動或跨樣本確認）

當修改影響多個樣本時，使用 parallel-benchmark agent 同時驗證：

```
Agent(prompt="驗證 C++ 修改對多樣本的影響:
  TASKS:
  1. scripts/run_batch_vcf_analysis.sh --mode HCC1395_5kHz   # 主要 discovery
  2. scripts/run_batch_vcf_analysis.sh --mode HCC1395_DORADO  # cross-platform
  3. scripts/run_batch_vcf_analysis.sh --mode COLO829          # cross-sample
  OUTPUT_DIR: research/autoresearch/cycles/cpp_change_validation/
  ", subagent_type="parallel-benchmark")
```

#### F1 delta 計算

```python
# 使用 collect_baseline_metrics.py 對比修改前後
python3 scripts/analysis/collect_baseline_metrics.py
```

commit 格式：
```
docs: [問題名稱] validation result

delta_f1:
  HCC1395-5kHz:  [+/-0.XXXX]
  COLO829:       [+/-0.XXXX]
  [其他已測樣本]

結論: [通過/回退/無顯著影響]
```

---

### Step 6：EVIDENCE LEDGER RECORD（研究記錄）

追加至 `research/autoresearch/evidence_ledger.jsonl`：

```json
{
  "cycle_id": "YYYY-MM-DD_[問題名稱]",
  "hypothesis_id": "H0XX",
  "hypothesis": "[一句話說明修改了什麼]",
  "pipeline_track": "TO",
  "type": "cpp_improvement",
  "tier": 3,
  "delta_f1": {
    "HCC1395-5kHz": 0.XXXX,
    "COLO829": 0.XXXX
  },
  "human_decision": "keep",
  "key_observations": "[關鍵觀察]",
  "methodology_doc": "docs/methodology/YYYYMMDD_*.md"
}
```

commit 格式：
```
chore: evidence_ledger record [H_ID] cpp_improvement [問題名稱]
```

---

## Commit 類型完整定義

| 類型 | 用途 | 強制使用場景 |
|------|------|------------|
| `chore: snapshot baseline` | 修改前快照 | Step 1 必用 |
| `test:` | 單元測試修改 | Step 3 必用 |
| `feat:` | 新功能/新邏輯 | Step 4，新增行為時 |
| `fix:` | Bug 修正 | Step 4，修正錯誤時 |
| `refactor:` | 重構無行為改變 | Step 4，只改結構時 |
| `docs:` | 驗證結果記錄 | Step 5 必用 |
| `chore: evidence_ledger` | 研究記錄 | Step 6 必用 |

---

## 回退檢查點

若 Step 5 驗證顯示 F1 大幅下降（delta < -0.002）：

1. **不要繼續 Step 6**
2. 更新 `docs/methodology/YYYYMMDD_*.md` 的狀態為 `validated_failed`
3. `git revert` 回 Step 1 的 baseline commit
4. 回到 `methodology-audit` 討論其他方案

---

## 相關工具

- `/build` — 編譯
- `/test-quick` — 快速測試
- `/test-data` — 數據測試
- `/test-full` — 完整測試
- `methodology-audit` skill — 審查入口
- `collect_baseline_metrics.py` — 跨樣本基線對比
- `pr-review-toolkit:code-reviewer` — Step 4.5 自動程式碼審查
- `pr-review-toolkit:silent-failure-hunter` — 錯誤處理審查（可選）
- `pr-review-toolkit:type-design-analyzer` — 型別設計審查（可選）

---

## 與 /scientific-rigor 元方法論的關係

本 skill 為 `/scientific-rigor §6 消融實驗設計` 的 **C++ 程式碼層級具體實作**：
- §6 「最小單元改動 / 快速反饋 / 4 步紀錄」 → 本 skill Step 1 BASELINE + Step 2 IMPLEMENT 對應「一次只改一個變數」
- §9 PDCA cycle → 本 skill 6 steps 為 Do 階段；Step 5 VALIDATE 為 Check 階段；Step 6 EVIDENCE LEDGER 為 Act 階段
- §7.2 可重現性 7 項 → 本 skill Step 1 baseline F1 + Step 6 commit hash binding 滿足 reproducibility
- §8.3 Reflexion → Step 5 若 F1 大幅下降 → 走 §9.2 SRE Postmortem template 紀錄 lessons
