---
name: cpp-change
description: PDD C++ 修改協議（6步驟）。在 methodology-audit 完成且用戶選定方案後啟動。確保每步驟有對應 commit，修改可追溯、可回退。觸發條件：「開始實作 [審查文件名]」、「執行方案 B」。
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

# 確認編譯
cd build && make -j$(nproc) 2>&1 | tail -10
```

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

### Step 5：VALIDATE COMMIT（驗證結果）

```bash
# 快速驗證（< 30 秒）
# /test-quick

# 數據驗證（約 1 分鐘）
# /test-data

# 完整驗證（約 5-10 分鐘，可選）
# /test-full
```

執行後計算 F1 delta：

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
