# 核心測試紀律

每輪研究迴圈必須遵守以下四項原則。

## 原則 A：Git Commit 保護每一個測試點

每次修改都可能影響 F1，沒有版本控制就無法回溯。

```bash
# 1. 每輪開始前：確認乾淨 baseline
git status
git stash list

# 2. 執行測試前：commit baseline
git add research/autoresearch/ scripts/analysis/ docs/
git commit -m "research: [CYCLE_ID] pre-change baseline — [H_ID] 測試前快照"

# 3. 記錄結果後：commit result
git add research/autoresearch/ scripts/analysis/ docs/
git commit -m "research: [CYCLE_ID] result — [H_ID] delta_F1=[+/-X.XXXX] ([verdict])"
```

**Commit 命名規範**：
- `research: [CYCLE_ID] pre-change — [H_ID] [假設摘要]`
- `research: [CYCLE_ID] result — [H_ID] delta=[X] ([verdict])`
- `research: [CYCLE_ID] combination — [H_A]+[H_B] delta=[X]`

**回溯節點**：若 negative 或 destructive：
```bash
git log --oneline -10
git checkout [HASH] -- scripts/analysis/[修改的檔案]
```

---

## 原則 B：單獨測試優先（Isolated Testing）

每次只改動一個變數。

```
❌ 同時修改 QS threshold 和 hp_assign_rate filter
✓ 先測 QS>=50 alone → 記錄 → 再獨立測 hp_assign_rate
```

**流程**：
1. 選取單一假設 H_X
2. git commit baseline（pre-change）
3. 修改單一變數
4. 執行 benchmark → delta_F1
5. git commit result
6. 記錄 H_X 的獨立貢獻
7. 還原（negative）或保留（positive）

---

## 原則 C：組合分析（Combination Testing）

**觸發條件**：至少 2 個「單獨測試 positive」的假設。

| 組合類型 | 判斷 | 動作 |
|----------|------|------|
| 協同（Synergistic）| combo > A + B | 使用組合 |
| 加法（Additive）| combo ≈ A + B | 可組合 |
| 正交（Orthogonal）| 候選集無重疊 | 可組合 |
| 干涉（Interfering）| combo < max(A, B) | 選最佳單一 |

**組合測試流程**：
```bash
# 1. checkout 兩假設的共同 baseline
git checkout [BASELINE_HASH] -- scripts/analysis/[target_file]

# 2. 同時應用 H_A + H_B

# 3. benchmark → delta_combo

# 4. 分析
python3 -c "
deltaA, deltaB, delta_combo = [A], [B], [COMBO]
synergy = delta_combo - (deltaA + deltaB)
print('synergy =', synergy, '（> 0 協同, < 0 干涉）')
"
```

**組合結論**寫入 evidence_ledger：
```json
{
  "hypothesis_id": "[H_A]+[H_B]",
  "combination_type": "additive|synergistic|orthogonal|interfering",
  "delta_a_alone": 0.0,
  "delta_b_alone": 0.0,
  "delta_combined": 0.0,
  "recommendation": "use_combo | use_a_only | use_b_only | skip_both"
}
```

---

## 原則 D：結果分類與研究潛力標記

每個假設測試後標記：

> **verdict 層次說明**：此處 `verdict` 是**單輪測試結論**，不同於假說最終判定（positive/negative/NO-GO/conditional）。
> 單輪 verdict 匯入 evidence_ledger；最終判定寫在 Hypothesis Tracker。

| 維度 | 可選值 | 說明 |
|------|--------|------|
| `verdict` | positive_pilot / negative / dataset_specific / annotation_only | 單輪測試結論 |
| `orthogonality` | tested / untested | 是否測過組合 |
| `research_potential` | high / medium / low / exhausted | 後續潛力 |
| `mechanism_clarity` | clear / partial / unclear | 機制是否理解 |

**Research Potential 定義**：
- `high`：positive 且機制合理，有擴展空間
- `medium`：positive 但 delta 小或 dataset-specific
- `low`：negative 但有邊際訊號
- `exhausted`：完全 negative 且機制排除
