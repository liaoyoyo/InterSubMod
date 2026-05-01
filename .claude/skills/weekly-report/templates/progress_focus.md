# 進展型主線（Progress Focus）

## 用途

本週有具體研究進展、突破、或階段性完成，向教授報告 **what was achieved + impact + next**。

## 觸發 keyword

「進展、突破、完成、達成、驗證、確認、解決、shipped、milestone」

## narrative skeleton

```
背景（Layer 0.1 一段）→ 本週處理（Layer 2 Thread）→ 結果（[F] 主導）
→ 初步分析（Layer 3）→ 下週銜接（Layer 4 §16）
```

## 17 段填寫深度指引

| 段 | 重要度 | 進展型寫法 |
|:-:|:-:|---------|
| §1 主線（≤30 字）| ⭐⭐⭐ | 「本週驗證 X 假設，N 樣本得 ΔF1=+0.0Y」 |
| §2 一句話重點 | ⭐⭐⭐ | 強調 quantified delta + 對 main thesis 影響 |
| §3 [F] 已確認 | ⭐⭐⭐ | **重點段**：列出本週驗證的具體事實，含 source path |
| §4 [O][I] | ⭐⭐ | 限於支持進展的 supporting evidence |
| §5 [U] 待確認 | ⭐ | 略寫；通常進展型 [U] 不多 |
| §7-§8 排序 | ⭐⭐⭐ | **重點段**：把進展放最前 |
| §9-§10 PPT 模板 | ⭐⭐ | 通常推薦 `improvement_report` 或 `data_showcase` |
| §11 補資料 | ⭐ | 略寫 |
| §16 下週 | ⭐⭐⭐ | **重點段**：明確下週要做什麼接續本週發現 |
| §17 教授追問 | ⭐⭐ | 預測：「方法選擇？baseline 對比？scale up？」|

略寫段（佔比 < 15% 母稿）：§5、§6、§11、§12、§13、§15

## 範例：進展型一句話母稿主線

> 本週完成 7 樣本 paired-pure F1 benchmark，ΔF1=+0.0112（v.s. baseline 0.7042），跨樣本一致 7/7，下週進入 Phase 2 Normal Methylation。

## 教授可能追問（5 個典型問題 + 預備回答骨架）

1. **「為什麼是這個方法？」** → 回顧 Phase 1A 三方向比對表
2. **「跟 baseline 比 ΔF1 真的有意義嗎？統計顯著？」** → bootstrap CI / paired t-test 結果
3. **「7 樣本一致性怎麼定義？」** → /multi-sample-consistency skill 的 canonical 排序表
4. **「下一步邊界在哪？」** → §16 + 已知失敗條件（如 TO 模式不適用）
5. **「資源夠嗎（compute/storage/time）？」** → 預估時間 + 資源需求

## 與其他主線類型的差別

- **vs problem_focus**：進展型有具體成果可報告；問題型卡在 blocker 求建議
- **vs advisor_consult**：進展型用戶已決定方向；求協助型仍在多選項評估
- **vs new_direction_explore**：進展型是已驗證假設的 follow-through；探索型是 pilot 性質

## 反例（不該用 progress_focus 的場景）

- 本週只跑了 pilot，N=2 樣本 → 應改 `new_direction_explore`
- 本週遇到 blocker 無法產出結果 → 應改 `problem_focus`
- 本週主要是讀文獻 + 整理 → 應改 `concept_walkthrough`（pptx-build 的模板，非本 skill）
