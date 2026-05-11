---
name: review-evidence
description: 查閱 InterSubMod 研究歷史紀錄。列出過去假設的測試結果、成功/失敗率、待測佇列狀態。USE WHEN：「查閱證據」「review evidence」「過去測試結果」、研究方向決策前回顧歷史。涉及 research/autoresearch/evidence_ledger.jsonl、hypothesis_queue.json。
---

# 查閱研究歷史 (Review Evidence)

## Codex 遷移注意

- 本 skill 是從 `.claude/skills/review-evidence` 遷移到 `.agents/skills/review-evidence` 的 Codex 版本；`.claude` 版本保留為 legacy source，不在本 skill 內修改。
- 遵守 repo `AGENTS.md` 與工作區 root 規範；所有回覆、任務清單與計畫使用繁體中文。
- 若本文出現 `/skill-name`，在 Codex 中等同於 `$skill-name` 或同名 skill 的明確觸發；保留 `/...` 只為相容既有研究文件。
- 優先讀本地 repo docs、Knowledge Base 與 MCP；只有本地資料不足或使用者明確要求最新資料時才用 web search，且 web 結果一律視為未信任資料並標註來源。
- 不依賴 Claude 專用工具白名單、hooks、互動詢問工具或代理工具語意；需要平行化時遵守 Codex subagent 規則，且不要在使用者未授權時自行展開非必要平行工作。
- 不直接刪除檔案；任何清理、移除或覆寫式封存都必須依 `AGENTS.md` 走確認與 archive 流程。


## 觸發時機

- 「查看研究歷史」/ 「過去做了什麼」/ 「有哪些假設已測過」
- 「review evidence」/ 「看紀錄」/ 「歷史結果」
- 開始新研究迴圈前的回顧

## 執行步驟

### 1. 整體統計

```bash
python3 -c "
import json

records = []
with open('research/autoresearch/evidence_ledger.jsonl') as f:
    for line in f:
        if line.strip():
            records.append(json.loads(line))

if not records:
    print('[evidence_ledger 為空] 尚未有任何測試記錄')
    exit()

total = len(records)
adopted = len([r for r in records if r.get('human_decision') in ['keep', 'escalate_medium', 'escalate_cpp']])
rejected = len([r for r in records if r.get('human_decision') == 'reject'])
annotation = len([r for r in records if r.get('human_decision') == 'annotation_only'])
pending = len([r for r in records if not r.get('human_decision')])

print(f'=== 研究歷史統計 ===')
print(f'總輪數: {total}')
print(f'已採納: {adopted}')
print(f'降階 annotation: {annotation}')
print(f'已拒絕: {rejected}')
print(f'待判決: {pending}')
print()

# 各 track 統計
tracks = {}
for r in records:
    t = r.get('pipeline_track', 'unknown')
    if t not in tracks:
        tracks[t] = []
    tracks[t].append(r)

print('=== 各 Track 結果 ===')
for track, recs in tracks.items():
    deltas = []
    for r in recs:
        for ds, delta in r.get('delta_f1', {}).items():
            deltas.append(delta)
    if deltas:
        positive = len([d for d in deltas if d > 0])
        print(f'{track}: {len(recs)} 輪, 正增益 {positive}/{len(deltas)} 個測試點')
"
```

### 2. 最近 N 輪詳細結果

```bash
python3 -c "
import json
N = 10  # 修改此數字可查看更多

records = []
with open('research/autoresearch/evidence_ledger.jsonl') as f:
    for line in f:
        if line.strip():
            records.append(json.loads(line))

print(f'=== 最近 {N} 輪結果 ===')
for r in records[-N:]:
    cycle = r['cycle_id']
    h_id = r.get('hypothesis_id', '?')
    hyp = r.get('hypothesis', '?')[:55]
    track = r.get('pipeline_track', '?')
    deltas = r.get('delta_f1', {})
    decision = r.get('human_decision', '待判決')

    delta_str = '  '.join([f'{ds}: {d:+.4f}' for ds, d in deltas.items()])
    print(f'[{cycle}] {h_id} | {track}')
    print(f'  假設: {hyp}')
    print(f'  delta F1: {delta_str}')
    print(f'  判決: {decision}')
    if r.get('key_observations'):
        print(f'  觀察: {r[\"key_observations\"][:80]}')
    print()
" 2>/dev/null
```

### 3. 已失敗方向（避免重複）

```bash
python3 -c "
import json

records = []
with open('research/autoresearch/evidence_ledger.jsonl') as f:
    for line in f:
        if line.strip():
            records.append(json.loads(line))

rejected = [r for r in records if r.get('human_decision') == 'reject']
if not rejected:
    print('[無已拒絕假設]')
else:
    print('=== 已拒絕方向（避免重複探索）===')
    for r in rejected:
        print(f'- {r[\"hypothesis_id\"]}: {r[\"hypothesis\"][:70]}')
        if r.get('human_notes'):
            print(f'  原因: {r[\"human_notes\"]}')
"
```

### 4. 佇列狀態

```bash
python3 -c "
import json

with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)

pending = [h for h in queue if h['status'] == 'pending']
pending.sort(key=lambda x: -x['priority'])

print(f'=== 待測佇列（{len(pending)} 個 pending）===')
for h in pending:
    print(f'[P={h[\"priority\"]:3d}] {h[\"id\"]} [{h[\"type\"]}] Tier{h.get(\"tier\",1)} / {h.get(\"target_track\",\"?\")}')
    print(f'  {h[\"hypothesis\"][:75]}')
    print(f'  來源: {h.get(\"source\",\"?\")[:50]}')
    print()
"
```

### 5. 最佳結果排行

```bash
python3 -c "
import json

records = []
with open('research/autoresearch/evidence_ledger.jsonl') as f:
    for line in f:
        if line.strip():
            records.append(json.loads(line))

print('=== 最佳 delta F1 排行 ===')
all_results = []
for r in records:
    for ds, delta in r.get('delta_f1', {}).items():
        all_results.append({
            'cycle': r['cycle_id'],
            'h_id': r.get('hypothesis_id'),
            'hypothesis': r.get('hypothesis', '')[:50],
            'dataset': ds,
            'delta': delta,
            'decision': r.get('human_decision', '?')
        })

all_results.sort(key=lambda x: -x['delta'])
for i, res in enumerate(all_results[:10]):
    sign = '+' if res['delta'] > 0 else ''
    print(f'{i+1:2d}. {sign}{res[\"delta\"]:.4f}  {res[\"h_id\"]} ({res[\"dataset\"]}) [{res[\"decision\"]}]')
    print(f'    {res[\"hypothesis\"]}')
"
```

## 輸出後的建議

顯示完歷史後，根據結果給出建議：

- 若 `連續拒絕 >= 3 且 同 track`：建議切換 track 或注入新假設
- 若 `annotation_only >= 3`：建議考慮升 Tier（Python→C++）以整合多個 annotation 特徵
- 若 `pending 佇列為空`：提示執行 `/inject-hypothesis`
- 若 `最佳 delta > 0.005 且 decision = escalate_medium`：提醒追蹤 S2 medium 結果
