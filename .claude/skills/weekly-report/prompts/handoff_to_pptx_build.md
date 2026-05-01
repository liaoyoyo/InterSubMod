# handoff_to_pptx_build.md — C4 後 4 選 handoff prompt

## 使用時機

C4 母稿確認完成後，銜接給 pptx-build sub-skill 或結束週報任務。

## 觸發前置條件

- C4 已通過
- master_draft.md 已寫入磁碟（`status: ready_for_handoff`）
- frontmatter 完整（含 report_type / main_statement / suggested_pptx_template / 等）

## AskUserQuestion 模板

```yaml
question: "母稿已存於 InterSubMod/docs/reports/validated/YYYY/MM/<dir>/master_draft.md（17 段，N 字）。下一步？"
header: "Handoff 4 選"
multiSelect: false
options:
  - label: "A. 立即觸發 pptx-build"
    description: "進入 pptx-build skill 跑 P1-P5（outline → section → slide → visual review → speaker script），讀母稿後跳過 main thesis 鎖定。預估 1-2 hr 完成 PPTX。"
  - label: "B. 母稿留檔，下次手動觸發"
    description: "之後以 /pptx-build --from-draft <path> 啟動。給用戶中間休息空間。"
  - label: "C. 母稿即終點（不產 PPT）"
    description: "本次任務僅需母稿（如 internal status check / weekly note 用途）。週報任務結束。"
  - label: "D. 母稿留檔 + 加寫下週計畫"
    description: "母稿保留，並萃取 §16 下一步行動 + §17 教授追問 + §11 補資料 + §15 暫存 → 加產 next_week_plan.md。"
```

## 各選項處置流程

### A. 立即觸發 pptx-build

```python
# 觸發邏輯
master_draft_path = "InterSubMod/docs/reports/validated/YYYY/MM/<dir>/master_draft.md"
trigger_skill(
    skill_name="pptx-build",
    args=f"--from-draft {master_draft_path}"
)
# pptx-build 自動：
# 1. 讀 frontmatter → 鎖定 main_thesis / report_type / audience
# 2. 跳過 P1 main thesis 鎖定 + P1 報告類型識別 + §20 Tier 評分
# 3. 進 P2 outline checkpoint（用戶仍須確認 5-7 段拆分）
```

### B. 母稿留檔，下次手動觸發

```
告知用戶：母稿留檔位置 = <path>
建議下次啟動指令：
  /pptx-build --from-draft InterSubMod/docs/reports/validated/YYYY/MM/<dir>/master_draft.md
週報任務本次結束。
```

### C. 母稿即終點

```
告知用戶：母稿存於 <path>
不產 PPT。
建議：將母稿加入 docs/experiments/INDEX.md 或 docs/reports/INDEX.md（依分類）。
週報任務結束。
```

### D. 加寫下週計畫

```
1. 自動萃取下列段到 next_week_plan.md：
   - 母稿 §16 下一步行動清單
   - 母稿 §17 教授可能提問 + 回答準備（轉成「下週要準備哪些 evidence」）
   - 母稿 §11 需補資料（轉成「下週要產的 artifacts」）
   - 母稿 §15 暫存紀錄（標 [SHELVED] 等下週判斷）
   - 母稿 §5 [U] 待確認項目（轉成「下週 blocker 追蹤」）

2. next_week_plan.md frontmatter:
   ---
   title: 下週計畫 YYYYMMDD（基於 master_draft.md）
   type: next_week_plan
   parent_master_draft: <path>
   date: YYYY-MM-DD
   ---

3. 完成後告知用戶：
   - 母稿位置：<path>/master_draft.md
   - 下週計畫：<path>/next_week_plan.md
   週報任務本次結束（不產 PPT，但有 actionable 下週清單）。
```

## fast-track 模式

C4 後 handoff 為**必停** AskUserQuestion（不論互動 / 全自動）。
理由：用戶必須明示是否繼續 PPT；自動觸發可能造成連續長時間工作疲勞。

## 失敗處置

| 情境 | 處置 |
|------|------|
| 用戶選 A 但 pptx-build skill 未載入 | 退回 B（留檔）+ 提示手動指令 |
| 用戶選 A 但 frontmatter status != ready_for_handoff | 警告 + 要求先回到 C4 確認 |
| 用戶選 A 但 master_draft.md 不存在 | 嚴重錯誤，要求重跑 W7 |
| 用戶選 D 但無法萃取（§16/§17 為空）| 提示「母稿缺必要段，無法產下週計畫」+ 退回 C 或 B |

## 完整 handoff 規範 → references/HANDOFF_TO_PPTX_BUILD.md
