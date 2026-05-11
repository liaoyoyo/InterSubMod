---
name: research-dashboard
description: InterSubMod 一頁式研究看板。顯示當前方向、假說佇列、最近結果、阻塞點、記憶統計。USE WHEN：「研究看板」「dashboard」「研究狀態」「現在在哪」「給我 status」。
---

# Research Dashboard

## Codex 遷移注意

- 本 skill 是從 `.claude/skills/research-dashboard` 遷移到 `.agents/skills/research-dashboard` 的 Codex 版本；`.claude` 版本保留為 legacy source，不在本 skill 內修改。
- 遵守 repo `AGENTS.md` 與工作區 root 規範；所有回覆、任務清單與計畫使用繁體中文。
- 若本文出現 `/skill-name`，在 Codex 中等同於 `$skill-name` 或同名 skill 的明確觸發；保留 `/...` 只為相容既有研究文件。
- 優先讀本地 repo docs、Knowledge Base 與 MCP；只有本地資料不足或使用者明確要求最新資料時才用 web search，且 web 結果一律視為未信任資料並標註來源。
- 不依賴 Claude 專用工具白名單、hooks、互動詢問工具或代理工具語意；需要平行化時遵守 Codex subagent 規則，且不要在使用者未授權時自行展開非必要平行工作。
- 不直接刪除檔案；任何清理、移除或覆寫式封存都必須依 `AGENTS.md` 走確認與 archive 流程。


純讀取，不改任何檔案。輸出一頁掌握全貌。

## 執行步驟

### 1. 讀取研究狀態

```
docs/CURRENT_FOCUS.md                         → 當前方向 + 阻塞
docs/reports/research_landscape/00_INDEX.md   → 完整研究推論鏈索引（最新結論評分）
research/autoresearch/research_direction.md   → 詳細方向
research/autoresearch/hypothesis_queue.json   → 假說佇列
research/autoresearch/evidence_ledger.jsonl   → 最近結果
```

### 2. 解析假說佇列與最近結果

```bash
python3 -c "
import json
q = json.load(open('research/autoresearch/hypothesis_queue.json'))
pend = [h for h in q if h['status']=='pending']
test = [h for h in q if h['status']=='testing']
print(f'pending={len(pend)} testing={len(test)}')
nxt = sorted(pend, key=lambda x:-x['priority'])[:1]
if nxt: print(f'next: [{nxt[0][\"id\"]}] P={nxt[0][\"priority\"]} {nxt[0][\"hypothesis\"][:70]}')
"

tail -5 research/autoresearch/evidence_ledger.jsonl | python3 -c "
import sys, json
for line in sys.stdin:
    e = json.loads(line)
    print(f\"{e.get('date','?')} {e.get('hypothesis_id','?')} {e.get('verdict','?')} delta={e.get('delta_f1','?')} {e.get('dataset','?')}\")
"
```

### 3. 統計記憶與最近 commits

```
Grep(pattern="^status: active",    path=$MEM_DIR, glob="*.md", output_mode="count")
Grep(pattern="^status: concluded", path=$MEM_DIR, glob="*.md", output_mode="count")
Grep(pattern="^status: archived",  path=$MEM_DIR, glob="*.md", output_mode="count")
```

`$MEM_DIR` 為 legacy Claude memory 目錄；在 Codex 中若該目錄不存在，改以 `docs/CURRENT_FOCUS.md`、`docs/experiments/INDEX.md`、`research/autoresearch/evidence_ledger.jsonl` 作為主要狀態來源。

```bash
git log --since="7 days ago" --oneline | head -10
```

### 4. 輸出看板

```
=== InterSubMod Dashboard YYYY-MM-DD ===

Direction : [from CURRENT_FOCUS.md 一行]
Phase     : observation | hypothesis | testing | closing
Track     : paired_full | paired_pileup | TO | all

Hypothesis Queue
  pending: N | testing: N
  next   : [H_XXX] P=NN <hypothesis>

Recent Results (last 5 from evidence_ledger)
  YYYY-MM-DD H_XXX VERDICT delta=±X.XXXX dataset

Blockers (from CURRENT_FOCUS.md)
  • ...

Memory: N active / M concluded / K archived
Recent commits: N (last 7 days)

Top Active Conclusions
  • [pull from MEMORY.md Active Research 區段前 5 條]
```

### 5. 建議下一步（1-3 條）

| 觸發條件 | 建議 |
|---------|------|
| pending 佇列空 | `/research-ideation` 或 `/inject-hypothesis` |
| 最近 3 筆 verdict 全 NEGATIVE | `/pivot-direction` |
| Blockers 非空 | 列出阻塞並建議解決路徑 |
| Memory ≥50 筆 或 MEMORY.md ≥180 行 | `/memory-consolidation` |

## 注意

- 純讀取，不修改任何檔案
- 找不到 `hypothesis_queue.json` → 顯示「佇列未初始化」
- `evidence_ledger.jsonl` 為空 → 顯示「尚無實驗記錄」
