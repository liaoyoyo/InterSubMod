---
name: pivot-direction
description: 快速切換 InterSubMod 研究方向。更新 research_direction.md、調整 hypothesis_queue 優先分數、記錄轉向理由。適用於：連續失敗後換方向、用戶觀察觸發新方向、track 切換（paired→TO）。
allowed-tools: Read, Write, Bash
user-invocable: true
---

# 快速轉向 (Pivot Direction)

## 觸發時機

- 「換方向」/ 「pivot」/ 「切換 track」
- 「換到 TO」/ 「改做 paired」
- 研究迴圈返回「建議切換 TO track」
- 用戶說「這個方向不行，改試 [X]」

## 執行步驟

### 1. 讀取當前方向

```bash
cat research/autoresearch/research_direction.md
```

### 2. 確認轉向類型

詢問或根據用戶輸入判斷：

| 轉向類型 | 說明 | 典型觸發 |
|---|---|---|
| `track_switch` | 切換 pipeline track（paired ↔ TO） | paired 連續無效 |
| `feature_pivot` | 從規則調整轉向特徵探索 | rule_change 方向耗盡 |
| `tier_escalation` | 從 Python 升到 C++ 改進 | Python 層驗證有效但需大規模 |
| `dataset_pivot` | 更換焦點資料集（5kHz ↔ DORADO ↔ 其他樣本）| 發現 dataset-specific 現象 |
| `focus_injection` | 注入用戶新觀察（優先度最高）| 用戶說「我注意到 X 現象」 |

### 3. 執行轉向

#### Track Switch（paired → TO）

```bash
python3 -c "
import json
from datetime import datetime

# 1. 讀取當前佇列
with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)

# 2. 更新 TO 相關假設的優先分數 +20
for h in queue:
    if h['status'] == 'pending' and h.get('target_track') in ['TO', 'all']:
        h['priority'] = min(100, h['priority'] + 20)

# 3. 降低 paired-only 假設優先分數 -20
for h in queue:
    if h['status'] == 'pending' and h.get('target_track') in ['paired_full', 'paired_pileup']:
        h['priority'] = max(10, h['priority'] - 20)

with open('research/autoresearch/hypothesis_queue.json', 'w') as f:
    json.dump(queue, f, ensure_ascii=False, indent=2)
print('TO track 假設優先分數已提升')
"
```

#### Focus Injection（用戶新觀察）

```bash
# 自動呼叫 inject-hypothesis 流程，priority = 90
python3 -c "
import json
from datetime import datetime

new_h = {
    'id': '[下一個 H_ID]',
    'type': 'rule_change',
    'priority': 90,
    'hypothesis': '[用戶觀察轉換為具體假設]',
    'source': '用戶觀察 [TODAY_DATE]',
    'target_track': '[用戶指定的 track]',
    'target_datasets': ['HCC1395_5kHz_TO'],
    'scale': 'pilot',
    'track': 'fast',
    'tier': 1,
    'status': 'pending',
    'added_at': datetime.now().strftime('%Y-%m-%d'),
    'notes': '[用戶原始說法]'
}

with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)
queue.append(new_h)
with open('research/autoresearch/hypothesis_queue.json', 'w') as f:
    json.dump(queue, f, ensure_ascii=False, indent=2)
print('已注入新方向:', new_h['hypothesis'][:60])
"
```

### 4. 更新 research_direction.md

在文件中加入轉向記錄：

```bash
cat >> research/autoresearch/research_direction.md << EOF

---
## 轉向記錄 [$(date +%Y-%m-%d)]

**類型**: [PIVOT_TYPE]
**原因**: [轉向理由]
**新焦點**: [新的研究方向]
**前一輪結果**: [觸發轉向的最後幾輪結果摘要]
**棄置方向**: [若有明確棄置的方向，記錄在此]
EOF
```

### 5. 確認轉向結果

```bash
python3 -c "
import json

with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)

pending = sorted([h for h in queue if h['status'] == 'pending'], key=lambda x: -x['priority'])
print('=== 轉向後的佇列前 5 筆 ===')
for h in pending[:5]:
    print(f'[P={h[\"priority\"]:3d}] {h[\"id\"]} [{h.get(\"target_track\",\"?\")}] {h[\"hypothesis\"][:60]}')
"
```

## TO Track 切換的判斷基準

若符合以下任一條件，建議切換到 TO track：

| 條件 | 說明 |
|---|---|
| paired 連續 3 輪 `abs(delta) < 0.001` | 空間太小 |
| 所有 FP 類型 A/B/C 都已分析過 | 方向耗盡 |
| 用戶說「想看更多 FP」 | 明確指示 |
| TO baseline F1 顯著低於 0.80 | 改善空間更大 |

**TO track 優勢**：
- TO FP 數量 ~300-800（vs paired FP ~50-100）
- TO 有更多 Type-B/C FP 可分析
- TO 三層架構改善空間更大（F1 起點更低）
- TO 測試時間不比 paired 長很多（~4 min）
