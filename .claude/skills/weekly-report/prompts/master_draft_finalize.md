# master_draft_finalize.md — W7 → C4 prompt（fast-track 必停）

## 使用時機

C3 邏輯+追問完成後，AI 組裝 17 段母稿（Layer 0-4 結構），用戶逐段批准 / 修改。

## 觸發前置條件

- C0/C1/C2/C3 全部通過
- internal state 已備齊：W1 raw data、W2 main_statement、W3 4 層分類、W4 4 桶分流、W5 修正後描述、W6 教授追問

## AI 自動處理

依 LAYER_STRUCTURE.md §B mapping，把 internal state 寫入 master_draft.md：

```
Layer 0.1 → §1 主線（C1 main_statement）+ §2 一句話重點
Layer 0.2 → §13 補充定義 / 背景知識
Layer 0.3 → 上週前情提要（從 git log 自動抽）
Layer 1   → 已建立知識參考（§3 [F] 中 cross-week 引用）
Layer 2   → §3 [F] + §4 [O][I] + §5 [U] 主體（含證據卡 Tier 1/2/3 + 因果鏈 Mermaid）
Layer 3   → §11 補資料 + §12 補圖表 + 整合更新（3-5 點本週新認知）
Layer 4   → §16 下週行動 + §17 教授追問 + 風險評估
4 桶分流  → §6 不放 PPT 暫存 + §7-§8 排序 + §15 暫存紀錄
handoff   → §9 PPT 模板建議 + §10 投影片架構建議
```

## 母稿輸出位置

```
InterSubMod/docs/reports/validated/YYYY/MM/YYYYMMDD_週報_<主題>/master_draft.md
```

frontmatter 寫入規範 → `references/HANDOFF_TO_PPTX_BUILD.md` §2

## AskUserQuestion 模板

```yaml
question: "母稿 N 段已組裝（總長 M 字）。逐段檢視結果如何？"
header: "C4 母稿確認"
multiSelect: false
options:
  - label: "批准全文"
    description: "母稿無需修改，進 handoff（4 選 A/B/C/D）。"
  - label: "修改特定段（指出 §N）"
    description: "用戶列「§3 改成 X / §17 增加 Y」之類的指令。"
  - label: "重跑某 W 階段"
    description: "如「§4 不對，重跑 W3 分類」→ 退回 W3 重跑後回到 W7。"
  - label: "棄用整份母稿"
    description: "重大方向錯誤（如選錯主線）→ 從 W2 重新開始。"
```

## 用戶處置流程

| 用戶選擇 | AI 行動 |
|---------|---------|
| 批准全文 | 寫入磁碟（status=ready_for_handoff）→ 進 handoff_to_pptx_build prompt |
| 修改 §N | 等用戶指令 → 修改該段 → 顯示 → 確認 |
| 重跑 W | 退回對應 W → 重跑 → 回 W7 → 再 C4 確認 |
| 棄用 | 確認後刪除 master_draft.md → 退回 W2 |

## 母稿示意（壓縮版，完整範例 → examples/master_draft_example.md）

```markdown
---
title: 週報 2026MMDD - <主題>
type: weekly_master_draft
report_type: progress
main_statement: "本週驗證 X 假設，N=7 樣本得 ΔF1=+0.0Y"
status: ready_for_handoff
material_classification:
  facts: 5
  observations: 3
  inferences: 2
  unconfirmed: 1
priority_buckets:
  ppt: 6
  speaker_note: 8
  appendix: 3
  shelf: 2
suggested_pptx_template: improvement_report
estimated_pptx_slides: 18
professor_qa_count: 6
---

# §1 本週主線
[C1 main_statement]

# §2 一句話重點
...

## Layer 0 研究脈絡

### Layer 0.1 宏觀問題定位（含 §1 §2）
### Layer 0.2 背景知識（含 §13）
### Layer 0.3 上週前情提要

## Layer 1 已建立知識參考

## Layer 2 本週調查

### Thread A: <主題>
#### 證據卡 Tier 1
- §3 [F] 已驗證項目
#### 證據卡 Tier 2
- §4 [O] 初步觀察
- §4 [I] 推論
#### 證據卡 Tier 3
- §5 [U] 待確認
#### 因果鏈圖 (Mermaid)

### §7 報告重點優先順序
### §8 建議報告順序

## Layer 3 整合更新
- §11 需補資料
- §12 需製圖
- 本週新增認知 3-5 點

## Layer 4 未來方向

### §16 下一步行動清單
### §17 教授可能提問 + 回答準備
### 風險評估

## 附錄

### §6 不建議放入 PPT
### §15 暫存紀錄
### §9 建議 PPT 模板：improvement_report
### §10 建議投影片架構：18 張，含 4 張 before-after
```

## §0 Highlights 強制檢查（v2.1 新增）★

組裝母稿時，**強制**在 §1 之前插入 §0：

```markdown
## §0 Highlights (TL;DR)

⭐⭐⭐ **This Week's Verdict**: [≤ 50 字一句決定性結論]

### Top Findings (3-5 條)
1. ⭐⭐⭐ [F] ...
2. ⭐⭐⭐ [F]/[O] ...
...

### Top Asks (≤ 3 條 [U])
1. ⭐⭐⭐ [U] ...
...

### ⭐⭐⭐ Decisive Next Step (1 條)
> Priority 1: ...
```

C4 確認時 AI 自動檢查：
- ✅ §0 存在 + Verdict ≤ 50 字
- ✅ Top Findings 3-5 條（多了壓縮）
- ✅ Top Asks ≤ 3 條
- ✅ Decisive Next Step 1 條
- ✅ §17 教授追問已分「⭐⭐⭐ 必問 ≤ 3」+「⭐⭐ 可能問」
- ✅ §16 priority 已分「Decisive / Operational / Maintenance」
- ✅ 每 Thread 結論首行為 black bold One-line Verdict
- ✅ Verdict Detail 段（optional，> 50 字場景必有）

**任一未通過 → 強制回 W7 重組。**

### v2.2 robust 驗證腳本（修正點 3 — regex 切割 bug fix）

舊版用單一 regex match `### Thread .*?(?=###)` 在 Thread B 後接 §7（不是 ###）時失敗。
新版用 multiline split：

```python
# robust 切割：以 H3/H2 之一作為終止
def find_threads(content):
    sections = re.split(r'^(?=##+ )', content, flags=re.MULTILINE)
    return [s for s in sections if s.startswith('### Thread ')]

# 對每個 Thread 檢查 One-line Verdict
for t in find_threads(content):
    has_verdict = '⭐⭐⭐ One-line Verdict' in t
    # ...
```

或更簡單的 grep 替代驗證：`grep -c "One-line Verdict" master_draft.md ≥ 2 * Thread_count`。

詳見 `references/LAYER_STRUCTURE.md` §E。

## 混合主線 §1 後 Sub-thread 段檢查（v2.1 新增）★

若 frontmatter `report_type` 含冒號（如 `problem:progress`），§1 後**強制**加 Sub-thread 段（含主線/Sub-thread/教授視角優先序 3 列）。詳見 `references/LAYER_STRUCTURE.md` §F。

## fast-track 模式

C4 為**必停** checkpoint（依 Q9 用戶決策）。即使 fast-track，仍須暫停 AskUserQuestion。
理由：母稿是最終產出，不能未經用戶確認就 handoff。
