---
name: inject-hypothesis
description: 注入新假設到 InterSubMod 研究佇列（hypothesis_queue.json）。支援：用戶口述假設、從觀察自動生成、從文獻引導輸入。
allowed-tools: Read, Write, Bash
user-invocable: true
---

# 注入新假設 (Inject Hypothesis)

## 觸發時機

- 「加入新假設」/ 「inject hypothesis」/ 「新方向」
- 用戶描述一個想測試的想法
- 研究迴圈返回「佇列為空」
- 文獻閱讀後有新特徵要驗證

## 執行步驟

### 1. 讀取當前佇列

```bash
cat research/autoresearch/hypothesis_queue.json
```

### 2. 確認下一個 ID

```bash
python3 -c "
import json
with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)
if queue:
    ids = [int(h['id'][1:]) for h in queue]
    print('下一個 ID: H' + str(max(ids)+1).zfill(3))
else:
    print('下一個 ID: H001')
"
```

### 3. 依來源類型填寫假設

**用戶口述假設**：
依用戶說的內容，填入以下欄位：
```json
{
  "id": "H[XXX]",
  "type": "rule_change",
  "priority": 85,
  "hypothesis": "[用戶的假設內容，轉換為具體的可測試陳述]",
  "source": "用戶觀察 [YYYY-MM-DD]",
  "target_track": "TO",
  "target_datasets": ["HCC1395_5kHz_TO"],
  "scale": "pilot",
  "track": "fast",
  "tier": 1,
  "status": "pending",
  "added_at": "[今日日期]",
  "notes": "[用戶的原始說法]"
}
```

**來自文獻的特徵**：
```json
{
  "type": "literature_feature",
  "priority": 75,
  "source": "論文: [作者 et al., 年份] [期刊] [DOI或摘要說明]",
  "hypothesis": "[將文獻中的特徵轉換為具體的 InterSubMod 可測試假設]"
}
```

**來自設計弱點 W1/W2/W3**：
```json
{
  "type": "cpp_improvement",
  "priority": 90,
  "source": "設計弱點 W[1/2/3]：[弱點描述]",
  "hypothesis": "[具體的改進方向與預期效果]"
}
```

### 4. 寫入佇列

```bash
python3 -c "
import json
from datetime import datetime

new_hypothesis = {
    'id': '[H_ID]',
    'type': '[TYPE]',
    'priority': [PRIORITY],
    'hypothesis': '[HYPOTHESIS]',
    'source': '[SOURCE]',
    'target_track': '[TRACK]',
    'target_datasets': ['[DATASETS]'],
    'scale': 'pilot',
    'track': 'fast',
    'tier': [TIER],
    'status': 'pending',
    'added_at': datetime.now().strftime('%Y-%m-%d'),
    'notes': ''
}

with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)
queue.append(new_hypothesis)
queue.sort(key=lambda x: -x['priority'])
with open('research/autoresearch/hypothesis_queue.json', 'w') as f:
    json.dump(queue, f, ensure_ascii=False, indent=2)
print('已加入:', new_hypothesis['id'], '-', new_hypothesis['hypothesis'][:50])
print('目前佇列長度:', len([h for h in queue if h['status'] == 'pending']), '個 pending')
"
```

### 5. 確認結果

呈現更新後的佇列前 5 筆（按 priority 排序）：

```bash
python3 -c "
import json
with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)
pending = [h for h in queue if h['status'] == 'pending']
pending.sort(key=lambda x: -x['priority'])
print('=== 待測假設佇列（前5筆）===')
for i, h in enumerate(pending[:5]):
    print(f'{i+1}. [{h[\"id\"]}] P={h[\"priority\"]} [{h[\"type\"]}] {h[\"hypothesis\"][:60]}')
"
```

## 優先分數指引

| 情境 | 建議 priority |
|---|---|
| 用戶直接指示「優先測試」| 90-95 |
| 來自設計弱點 W1/W2/W3 | 85-90 |
| OBSERVE 中發現的 SIGNAL | 75-80 |
| 文獻引導特徵 | 70-80 |
| 參數組合探索 | 55-65 |
| 低風險嘗試性假設 | 40-55 |

## Type 類型選擇

- 修改篩選閾值 → `rule_change`（Tier 1）
- 新增特徵欄位（Python）→ `feature_hypothesis`（Tier 2）
- 修改 ISM 執行參數 → `param_combo`（Tier 1）
- 文獻中的新特徵 → `literature_feature`（Tier 1 or 2）
- 修改 C++ 核心邏輯 → `cpp_improvement`（Tier 3）
