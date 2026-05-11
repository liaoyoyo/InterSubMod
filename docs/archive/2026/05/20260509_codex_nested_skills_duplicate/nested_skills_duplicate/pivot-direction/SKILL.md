---
name: pivot-direction
description: 快速切換 InterSubMod 研究方向。注入用戶新觀察為高優先假設、降權連續失敗的方向、記錄轉向理由。USE WHEN：「換方向」「pivot」「我注意到 X 現象」「這方向不行」。
---

# 快速轉向 (Pivot Direction)

## Codex 遷移注意

- 本 skill 是從 `.claude/skills/pivot-direction` 遷移到 `.agents/skills/pivot-direction` 的 Codex 版本；`.claude` 版本保留為 legacy source，不在本 skill 內修改。
- 遵守 repo `AGENTS.md` 與工作區 root 規範；所有回覆、任務清單與計畫使用繁體中文。
- 若本文出現 `/skill-name`，在 Codex 中等同於 `$skill-name` 或同名 skill 的明確觸發；保留 `/...` 只為相容既有研究文件。
- 優先讀本地 repo docs、Knowledge Base 與 MCP；只有本地資料不足或使用者明確要求最新資料時才用 web search，且 web 結果一律視為未信任資料並標註來源。
- 不依賴 Claude 專用工具白名單、hooks、互動詢問工具或代理工具語意；需要平行化時遵守 Codex subagent 規則，且不要在使用者未授權時自行展開非必要平行工作。
- 不直接刪除檔案；任何清理、移除或覆寫式封存都必須依 `AGENTS.md` 走確認與 archive 流程。


## 觸發

- 「換方向」/「pivot」/「這個方向不行，改試 [X]」
- 「我注意到 X 現象」（用戶觀察 → 最高優先注入）
- research-loop 返回「建議切換」訊號

## 執行步驟

### 1. 讀取當前狀態

```bash
cat research/autoresearch/research_direction.md
python3 -c "
import json
q = json.load(open('research/autoresearch/hypothesis_queue.json'))
pending = sorted([h for h in q if h['status']=='pending'], key=lambda x:-x['priority'])
for h in pending[:5]:
    print(f\"[P={h['priority']:3d}] {h['id']} {h['hypothesis'][:60]}\")
"
```

### 2. 判斷轉向類型

| 類型 | 觸發 | 動作 |
|------|------|------|
| **focus_injection** | 用戶新觀察 | 用 `/inject-hypothesis` 注入 priority=90，其他降權無動作 |
| **direction_swap** | 用戶指定切換 | 升相關 +20、降被棄方向 -20 |
| **direction_close** | 連續 3 輪失敗 | 將該類所有 pending 假設 status 改 closed |

### 3. 執行（direction_swap 範本）

```bash
python3 -c "
import json
from datetime import datetime
queue = json.load(open('research/autoresearch/hypothesis_queue.json'))
TARGET_UP = '[要升權的 target_track 或 keyword]'
TARGET_DOWN = '[要降權的]'
for h in queue:
    if h['status'] != 'pending': continue
    text = (h.get('target_track','') + ' ' + h.get('hypothesis','')).lower()
    if TARGET_UP.lower() in text:
        h['priority'] = min(100, h['priority'] + 20)
    if TARGET_DOWN.lower() in text:
        h['priority'] = max(10, h['priority'] - 20)
json.dump(queue, open('research/autoresearch/hypothesis_queue.json','w'), ensure_ascii=False, indent=2)
print('已調整')
"
```

**focus_injection** → 直接呼叫 `/inject-hypothesis` skill，priority=90。

### 4. 記錄轉向

```bash
cat >> research/autoresearch/research_direction.md << EOF

---
## 轉向記錄 [$(date +%Y-%m-%d)]

**類型**: [focus_injection / direction_swap / direction_close]
**原因**: [一句]
**新焦點**: [新方向]
**棄置/降權**: [若有]
**前一輪結果**: [觸發轉向的最後 1-2 輪結果摘要]
EOF
```

### 5. 確認新佇列前 5 筆

重跑步驟 1 的 Python 命令，列印前 5 筆驗證優先序符合預期。

## 何時建議連續失敗 → direction_close

| 條件 | 動作 |
|------|------|
| 同一 target_track 連續 3 輪 \|delta F1\| < 0.001 | 該 track 全部 pending → closed，記錄至 evidence_ledger |
| 同一假設類型已耗盡所有合理變體 | 整類關閉，導向 inject 新方向 |

## 反模式

- ❌ 直接刪除舊假設（應改 status 為 closed 保留 provenance）
- ❌ 注入用戶觀察時忘記設 priority=90 → 排在後面
- ❌ 切換方向但未記錄到 research_direction.md → 後續無法追溯

## 與其他 skill 的關係

- 高優先注入 → 委派 `/inject-hypothesis`
- 連續失敗判定 → 對應 `/research-loop` 的 step 7 訊號
- 用戶決策需確認 → 對應 `/confirmation-protocol` 的 Gate 級別
