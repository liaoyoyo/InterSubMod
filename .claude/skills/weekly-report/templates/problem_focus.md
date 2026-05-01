# 問題型主線（Problem Focus）

## 用途

本週遇到明確 blocker、bug、或 anomaly，**未能產出預期結果**，需教授協助判斷或建議方向。

## 觸發 keyword

「問題、卡住、blocker、bug、anomaly、無法解釋、衝突、矛盾、failed、issue」

## narrative skeleton

```
方法（既有 baseline）→ 問題發現（具體現象）→ 目前判斷（[O][I][U] 為主）
→ 求建議（教授視角追問 §17 為主）
```

## 17 段填寫深度指引

| 段 | 重要度 | 問題型寫法 |
|:-:|:-:|---------|
| §1 主線（≤30 字）| ⭐⭐⭐ | 「X 現象與預期衝突，根因 [Y/Z] 待釐清」 |
| §2 一句話重點 | ⭐⭐⭐ | 突出衝突點 + 主要 hypothesis |
| §3 [F] 已確認 | ⭐⭐ | 限於「可重現的問題現象 + benchmark 數字」 |
| §4 [O][I] | ⭐⭐ | 列已排除的解釋 + 當前傾向假設 |
| §5 [U] 待確認 | ⭐⭐⭐ | **重點段**：列出無法解釋的觀察 + 各假設可否證條件 |
| §7-§8 排序 | ⭐⭐ | 把 root cause hypothesis 放最前 |
| §9-§10 PPT 模板 | ⭐⭐ | 通常推薦 `comparison_report`（對比預期 vs 實際）或 `improvement_report` |
| §11 補資料 | ⭐⭐⭐ | **重點段**：明確列「需要 X 資料才能下判斷」 |
| §16 下週 | ⭐⭐ | 短期 workaround + 長期 root cause investigation |
| §17 教授追問 | ⭐⭐⭐ | **最重點段**：問題型核心，明確列「我要問教授」 |

略寫段（佔比 < 15%）：§3 限於問題現象、§4 排除清單、§13、§14

## 範例：問題型一句話母稿主線

> ClairS-TO TO 模式 self-phasing 17.3:1 違反 1:1 預期，已排除 V1-V3 fix，根因待釐清需教授判斷下一步。

## 教授可能追問（5 個典型問題 + 預備回答骨架）

1. **「根因確認了嗎？」** → 列已排除清單 + 仍可能的 3 個假設
2. **「短期 workaround 可行嗎？」** → §16 列暫時 mitigation + 副作用
3. **「影響範圍多大？」** → 影響的下游分析 + 量化估計
4. **「是不是上游工具的問題？」** → 已查 ClairS-TO commit history + bug tracker
5. **「值得花多少時間 debug？」** → ROI 估算 + alternative paths 比較

## 與其他主線類型的差別

- **vs progress_focus**：問題型本週**沒有**達成預期；進展型有
- **vs advisor_consult**：問題型有明確 blocker；求協助型有多個 candidate path
- **vs new_direction_explore**：問題型是既有方向中的 stuck；探索型是新方向 pilot

## 反例（不該用 problem_focus 的場景）

- 本週其實達成 80%，只剩細節 → 應改 `progress_focus`，把問題列為附帶
- 本週發現 3 個 candidate fix 還沒選 → 應改 `advisor_consult`
- 本週是純粹的重跑驗證 → 應改 `progress_focus`
