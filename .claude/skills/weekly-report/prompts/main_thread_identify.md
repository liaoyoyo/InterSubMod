# main_thread_identify.md — W2 → C1 prompt（fast-track 必停）

## 使用時機

W1 raw data 確認後，識別本週主線類型（4 選 1）+ 鎖定 main statement（≤ 30 字）。

## 觸發前置條件

- C0 已通過
- W1 raw data 結構化清單已備齊

## AI 自動推斷（基於 W1 內容）

| W1 特徵 | 推斷主線 |
|---------|---------|
| commit 多 + 有具體量化 delta + evidence ledger 標 [F] | progress |
| commit 少 + experiments 標 [BLOCKER]/[ANOMALY] | problem |
| 多個 candidate path 都有部分 evidence | advisor_consult |
| pilot N≤3 + experiments 標 [PILOT] | new_direction_explore |

AI 顯示推斷結果 + 用戶確認 / 改選。

## AskUserQuestion 模板

```yaml
question: "本週主線類型？AI 推斷為 <X>，理由：<Y>。"
header: "C1 主線類型"
multiSelect: false
options:
  - label: "進展型 (progress_focus)"
    description: "本週有具體研究進展、突破、量化 delta；報告 what was achieved + impact + next。"
  - label: "問題型 (problem_focus)"
    description: "本週遇到 blocker / bug / anomaly；報告問題現象 + 目前判斷 + 求建議。"
  - label: "求協助型 (advisor_consult)"
    description: "本週有 2-4 個 candidate direction；請教授協助判斷優先序。"
  - label: "探索型 (new_direction_explore)"
    description: "本週啟動新方向 pilot (N≤3)；報告動機 + pilot 結果 + 是否值得 scale up。"
```

## 鎖定 main statement（≤ 30 字）

主線類型確認後，AI 起草 main statement，用戶可改：

```
AI 草稿 (進展型範例):
"本週驗證 X 假設，N=7 樣本得 ΔF1=+0.0Y，下週進入 Z"
                ↓
用戶確認 / 修改 / 重寫
                ↓
鎖定為 master_draft frontmatter `main_statement` 欄位
```

**規則**：
- ≤ 30 字（中文）/ ≤ 50 chars（英文）
- 必含具體數字或具體動詞（不可空泛）
- 必能對應到 §1 母稿主線

## 預期輸出格式

```markdown
### W2 主線確認

**Report Type**: progress | problem | consult | exploration
**Main Statement**: "<≤ 30 字>"
**Template Loaded**: templates/{type}_focus.md（後續 W3-W7 引用此模板）
```

## 用戶處置流程

| 用戶行動 | AI 回應 |
|---------|---------|
| 接受推斷 | 進入 W3 |
| 改選主線 | 重新載入對應 template，提示用戶補述為何不是 AI 推斷 |
| 主線混合（如進展+問題）| 取主敘事為主，附 sub-thread 標註 |
| 用戶 main statement 超 30 字 | 提示「過長，請壓縮」並建議刪減 |

## fast-track 模式

C1 為**必停** checkpoint（依 Q9 用戶決策）。即使 fast-track，仍須暫停 AskUserQuestion。
理由：主線類型決定整份母稿走向，AI 推斷可能誤判（如進展+問題混合場景），用戶必須確認。

---

## v2.1 升級：混合主線 + main statement follow-up（修正點 2 + 4）

### 混合主線 (修正點 4)

C1 主線可寫複合 `主線:子線`（如 `problem:progress`）。AskUserQuestion 加第 5 option（用「其他混合」option 觸發 free-text）：

```yaml
- label: "混合主線（注：主:sub）"
  description: "如 problem:progress（主線問題、sub 進展）。母稿 §1 後強制加 Sub-thread 段。詳見 references/LAYER_STRUCTURE.md §F"
```

混合主線時，AI 在 §1 後組「Sub-thread」段：
- 主線（依 main statement）
- Sub-thread（次要 thread）
- 教授視角優先序

### main statement follow-up (修正點 2)

C1 確認主線後，AI **加問第 2 個 AskUserQuestion**：

```yaml
question: "main statement '<AI 草稿 30 字>' 是否要修改？"
options:
  - 接受草稿
  - 微調（用戶提供改寫）
  - 完全重寫
  - 改主線類型（退回上一步）
```

防止用戶選對主線但希望調 statement 的場景。
