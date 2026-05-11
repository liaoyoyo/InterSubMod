---
name: review-evidence
description: 查閱 InterSubMod 研究歷史紀錄。列出過去假設的測試結果、成功/失敗率、待測佇列狀態。USE WHEN：「查閱證據」「review evidence」「過去測試結果」、研究方向決策前回顧歷史。涉及 research/autoresearch/evidence_ledger.jsonl、hypothesis_queue.json。SKIP WHEN 全域系統審計（用 provenance-tier-audit）、cycle-level dashboard（用 cycle-state）、純 build / commit / docs、新假說注入（用 inject-hypothesis）。
allowed-tools: Read, Bash
user-invocable: true
---

# 查閱研究歷史 (Review Evidence)

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

## DECLARE 7-letter cognitive forcing — pivot decision 前必跑（v1.8 T2-5）

> 設計來源：DECLARE / Croskerry cognitive forcing strategies。**review-evidence 用於 pivot 或繼續決策的場景時必跑**；單純歷史查詢可跳過。

| 字母 | 提問 | 對 InterSubMod 怎麼答 |
|---|---|---|
| **D** ifferential | 「我接下來想驗證的假說，至少還有哪 2 種競爭解釋？」 | 列 ≥2 個 alternative；確認 review 結果不是 cherry-pick |
| **E** vidence | 「review 的 evidence_ledger 紀錄涵蓋 L1-L4 哪幾層？」 | 若只有 L1-L2 → 不足以 pivot；要拉到 L3/L4 |
| **C** onfounders | 「過去失敗結論的 root cause 與目前 pivot 計畫有重疊嗎？」 | 若有 → 修正 confound design 才 pivot；否則重蹈覆轍 |
| **L** ikelihood | 「最佳 delta 是真信號還是 noise？樣本數夠嗎？」 | 對照 cross-sample n / 雙 p-value；< 3 樣本不可作 tier 升級依據 |
| **A** lternatives | 「過去 NEGATIVE 結論裡，有沒有與我 pivot 方向衝突的？」 | 若有 → 必須在 plan.json 寫明為何此次能 escape 過去失敗模式 |
| **R** eassess | 「review 跨 cycle 結論彼此一致嗎？」 | 不一致 → 啟動 provenance-tier-audit 而非 pivot |
| **E** ngage | 「pivot 決策會更新 hypothesis_queue 與 CURRENT_FOCUS 嗎？」 | 必須同步更新；否則新方向無法被 dashboard 看見 |

## 輸出後的建議

顯示完歷史後，根據結果給出建議：

- 若 `連續拒絕 >= 3 且 同 track`：建議切換 track 或注入新假設
- 若 `annotation_only >= 3`：建議考慮升 Tier（Python→C++）以整合多個 annotation 特徵
- 若 `pending 佇列為空`：提示執行 `/inject-hypothesis`
- 若 `最佳 delta > 0.005 且 decision = escalate_medium`：提醒追蹤 S2 medium 結果
