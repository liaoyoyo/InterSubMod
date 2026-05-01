# raw_data_collect.md — W1 → C0 prompt

## 使用時機

skill 啟動後第一步。AI 自動掃描本週 raw data 來源後，向用戶確認完整性 / 補漏 / 移除誤抓。

## 觸發前置條件

- skill 觸發詞匹配（週報 / weekly report / 整理本週 等）
- 用戶模式：互動模式或全自動

## AI 自動掃描清單（W1 細則 → references/COLLECTION_PROTOCOL.md）

```bash
# 1. Git 活動
git log --since="7 days ago" --pretty=format:"%h %s (%ad)" --date=short

# 2. 新增 / 修改的 .md（experiments / reports / docs）
find docs/experiments/in_progress -newer .claude/last_weekly_report -name "*.md"

# 3. evidence_ledger 本週追加
tail -50 research/autoresearch/evidence_ledger.jsonl | jq 'select(.date >= "...")'

# 4. CURRENT_FOCUS 變動
git log --follow -p docs/CURRENT_FOCUS.md --since="7 days ago"

# 5. 上週 weekly report Layer 3-4
ls -t docs/reports/validated/YYYY/MM/ | head -2
```

## AskUserQuestion 模板

```yaml
question: "本週 raw data 自動掃描結果：N 個 commit、M 個 .md、K 個 evidence ledger 條目。是否完整？"
header: "C0 確認"
multiSelect: false
options:
  - label: "完整，繼續 W2"
    description: "AI 掃描已涵蓋本週主要活動。"
  - label: "補漏（用戶補述）"
    description: "AI 漏掉某些重要進展，用戶口述補充。"
  - label: "移除誤抓"
    description: "某些 commit / 檔案不該屬於本週主軸（如自動化雜項），請列明。"
  - label: "全部重來"
    description: "AI 掃描範圍錯誤（如時間區間不對），重新掃描。"
```

## 預期輸出格式

C0 確認後，輸出 W1 結構化清單（暫存於 in-conversation context）：

```markdown
### W1 Raw Data 結構化清單

**時間軸**：
- YYYY-MM-DD: <commit hash> <主題>
- ...

**產出 artifacts**：
- InterSubMod/docs/experiments/.../<file>.md
- ...

**evidence ledger 新增**：
- [hypothesis_id] <verdict>
- ...

**用戶補述**：
- ...
```

## 用戶處置流程

| 用戶選擇 | AI 行動 |
|---------|---------|
| 完整 | 進 W2 main_thread_identify |
| 補漏 | 接收用戶補述 → 加入清單 → 進 W2 |
| 移除誤抓 | 接收用戶指示 → 從清單移除 → 進 W2 |
| 全部重來 | 詢問正確時間區間或來源 → 重跑掃描 |

## fast-track 模式

C0 為非必停 checkpoint。fast-track 下：
- AI 掃描後直接列「以此清單推進，30 秒內無回覆視為確認」
- 不暫停 AskUserQuestion
- 結果存於 W1 lab notebook 風格清單，後續 C2/C3 用戶仍可調整
