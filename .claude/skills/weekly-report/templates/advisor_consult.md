# 求協助型主線（Advisor Consult）

## 用途

本週有 2-4 個 candidate research direction（或 candidate fix），自己無法決定優先序，**請教授協助判斷方向**。

## 觸發 keyword

「不確定、需要 advisor、方向選擇、candidate、option、不知道該、教授決定、優先序」

## narrative skeleton

```
情境（為何到這個交叉路口）→ 多選項（每個 candidate 各列利弊）
→ 我的傾向 + 不選的代價 → 待教授決策點 + 時程壓力
```

## 17 段填寫深度指引

| 段 | 重要度 | 求協助型寫法 |
|:-:|:-:|---------|
| §1 主線（≤30 字）| ⭐⭐⭐ | 「N 個 candidate direction，請教授判斷下一步」 |
| §2 一句話重點 | ⭐⭐⭐ | 強調「決策點 + 時程影響」|
| §3 [F] 已確認 | ⭐⭐ | 限於「為什麼到這個交叉路口」的背景事實 |
| §4 [O][I] | ⭐⭐⭐ | 列每個 candidate 的當前 evidence + 預期收益 |
| §5 [U] 待確認 | ⭐⭐⭐ | **重點段**：每個 candidate 的「需要驗證才能下判斷」清單 |
| §7-§8 排序 | ⭐⭐ | 列「我的傾向順序 + 理由」（但提示教授可推翻）|
| §9-§10 PPT 模板 | ⭐⭐ | 推薦 `comparison_report` 或 `decision_discussion` |
| §11 補資料 | ⭐⭐ | 每 candidate 缺哪些資料才能 fully evaluate |
| §16 下週 | ⭐⭐⭐ | **重點段**：依教授可能決策列 N 種 contingency plan |
| §17 教授追問 | ⭐⭐⭐ | **最重點段**：明確列「我要請教授判斷的問題」|

略寫段（佔比 < 15%）：§3 限於 setup、§13、§14

## 範例：求協助型一句話母稿主線

> Phase 1B Normal Methylation 有 3 個 candidate path（Sample ASM / LOH BED / Cross-region），各有利弊，請教授協助判斷下週投入優先序。

## 多 candidate 標準呈現格式（§4 + §5 推薦）

每個 candidate 必填欄位：

```markdown
### Candidate A: <名稱>

**evidence to date**: [F] 列出已驗證部分（如有）
**expected outcome**: [I] 推測：「若選 A，下週可達成 X」
**cost**: 預估時間 + 資源
**risk**: [U] 失敗條件（什麼情況下 A 會 NEGATIVE）
**dependency**: 需要哪些 prerequisite
**my preference**: 我傾向選 A 的理由（一段）
```

## 教授可能追問（5 個典型問題 + 預備回答骨架）

1. **「你傾向哪個？為什麼？」** → §4 my preference 段，evidence-based 理由
2. **「不選 A 的代價是什麼？」** → 預估錯失的收益 + opportunity cost
3. **「有沒有時程壓力？」** → §16 列「若 X 日前不決定，會錯過 Y」
4. **「能不能同時做？」** → 資源 / 時間 trade-off 估算
5. **「別人怎麼做？」** → 相關 paper / 同領域實踐參考

## 與其他主線類型的差別

- **vs progress_focus**：求協助型沒有單一 main thesis 可宣告；進展型有
- **vs problem_focus**：求協助型有多個可行 path；問題型卡在無法解
- **vs new_direction_explore**：求協助型 candidate 是 viable options；探索型是 pilot 階段

## 反例（不該用 advisor_consult 的場景）

- 只有 1 個方向但卡住 → 應改 `problem_focus`
- 候選方向已私下決定，只是想 sanity check → 應改 `progress_focus`（標方向選擇 rationale）
- pilot 結果未出，無法判斷 candidate 是否 viable → 應改 `new_direction_explore`
