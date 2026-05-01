# 探索型主線（New Direction Explore）

## 用途

本週啟動**新方向 pilot**（小規模 N≤3 樣本驗證），需向教授報告 **動機 + 假設 + pilot 結果 + 是否值得 scale up**。

## 觸發 keyword

「新方向、pilot、初步觀察、探索、試試看、preliminary、pilot study、feasibility」

## narrative skeleton

```
動機（為何冒出這方向）→ 假設（可否證條件）→ pilot 結果（[O] 為主）
→ scale up 風險評估 → 是否投入更多 N
```

## 17 段填寫深度指引

| 段 | 重要度 | 探索型寫法 |
|:-:|:-:|---------|
| §1 主線（≤30 字）| ⭐⭐⭐ | 「Pilot X，N=Y 結果 Z，是否 scale up 待議」 |
| §2 一句話重點 | ⭐⭐⭐ | 突出 motivation + pilot result preliminary nature |
| §3 [F] 已確認 | ⭐ | 限於「pilot setup 的 ground truth」（不是 pilot 結論）|
| §4 [O][I] | ⭐⭐⭐ | **重點段**：pilot 結果 + 推論可能 scale up 收益 |
| §5 [U] 待確認 | ⭐⭐⭐ | **重點段**：pilot 限制 + scale up 才能驗證的問題 |
| §7-§8 排序 | ⭐⭐ | 把「scale up vs 暫緩」決策放最前 |
| §9-§10 PPT 模板 | ⭐⭐ | 推薦 `data_showcase` 或 `concept_walkthrough` |
| §11 補資料 | ⭐⭐⭐ | **重點段**：scale up 需要哪些資料 + cost 估計 |
| §16 下週 | ⭐⭐ | 兩種 path：scale up vs 蹲守驗證更多 N |
| §17 教授追問 | ⭐⭐⭐ | 預測：「pilot 可信嗎？scale up 風險？資源？」 |

略寫段（佔比 < 15%）：§3 限於 setup、§13、§14

## 範例：探索型一句話母稿主線

> Pilot HPFineNGroups feature N=2 樣本顯示 AUC 0.62，需 7 樣本驗證 confound + scale up 風險評估。

## 探索型獨特區段：scale up risk assessment

新方向必填的「scale up 風險」5 項評估：

```markdown
### Scale up risk assessment

| 維度 | pilot (N=2) | 預期 scale (N=7) | 風險 |
|------|------------|----------------|------|
| 統計可信度 | weak (CI 寬) | strong | 低 |
| Confound 檢驗 | 未做 | within-group OLS | 中 |
| 計算成本 | 1 hr | 7 hr | 低 |
| 資料可用性 | OK | TODO 補 X 樣本 | 中 |
| Failure mode | unknown | spatial autocorr | 高 |
```

## 教授可能追問（5 個典型問題 + 預備回答骨架）

1. **「pilot 結果可信嗎？」** → CI / bootstrap / sanity check 結果
2. **「scale up 的風險？」** → 上方表格 + 列出最大未知
3. **「資源夠嗎？」** → 成本估算 + opportunity cost 對比現有方向
4. **「pilot 為什麼選這 2 樣本？」** → sampling rationale（典型 vs 邊緣）
5. **「失敗了怎麼辦？」** → fall-back plan + 對 main thesis 影響評估

## 與其他主線類型的差別

- **vs progress_focus**：探索型 N 不足、結論 preliminary；進展型 N 達標 + [F]
- **vs problem_focus**：探索型主動嘗試新方向；問題型被動卡在既有方向
- **vs advisor_consult**：探索型有具體 candidate（pilot 已跑）；求協助型 candidate 還在評估

## 反例（不該用 new_direction_explore 的場景）

- 已跑完 7 樣本 + confound check → 應改 `progress_focus`（已完成驗證）
- pilot 還沒跑，只是有 idea → 應改 `advisor_consult`（讓教授判斷要不要做 pilot）
- pilot 跑完且結果為 NEGATIVE 已棄用 → 應改 concluded report（不需 weekly-report）
