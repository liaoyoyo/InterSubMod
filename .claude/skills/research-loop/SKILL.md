---
name: research-loop
description: InterSubMod 半自動研究迴圈。八步驟（觀察→定向→假設→設計→執行→記錄→呈現→回饋），驅動甲基化特徵探索與 F1 提升。USE WHEN：「開始研究迴圈」「research loop」「測試新假設」「下一輪假設」。涉及 research/autoresearch/ 下的 hypothesis_queue.json、evidence_ledger.jsonl、research_direction.md。適用 paired_full / paired_pileup / TO 三條 pipeline track。
allowed-tools: Read, Write, Bash, Glob, Grep
user-invocable: true
---

# InterSubMod 研究迴圈 (Research Loop)

協助 AI 與人類研究者共同：觀察 TP/FP/FN 模式差異 → 生成可測試假設 → benchmark 驗證 → 記錄推進。

## 觸發時機

「開始研究迴圈」「research loop」「測試新假設」「下一輪假設」

## 前置檢查

```
research/autoresearch/hypothesis_queue.json   ← 必存在（可為空 []）
research/autoresearch/evidence_ledger.jsonl   ← 可為空
research/autoresearch/research_direction.md   ← 必存在
```

佇列不存在 → 提示 `/inject-hypothesis`。

---

## 執行模式感知

本 skill 遵循「確認時機協議」（詳見 `/confirmation-protocol` skill）的模式切換。確認當前是**互動模式**（預設）或**全自動模式**。

---

## 核心紀律（摘要）

詳見 `references/TEST_DISCIPLINE.md`，每輪必遵守：

- **A. Git Commit**：改前 commit baseline、改後 commit result
- **B. 單獨測試**：每次只改一個變數
- **C. 組合測試**：≥2 個 positive 後測試協同/正交
- **D. 結果分類**：verdict + research_potential + mechanism_clarity

---

## 八步驟執行流程

### 步驟 0：OBSERVE — 觀察數據 `[FYI]`

詳見 `references/OBSERVE_PROTOCOL.md`。

從真實數據發現現象。確認 track → 讀特徵檔 → 產出 TP/FP/FN 特徵側寫 → 識別 SIGNAL。

全自動模式：靜默執行，結果記入 cycle artifact。

---

### 步驟 1：ORIENT — 定向 `[FYI]`

詳見 `references/OBSERVE_PROTOCOL.md` 步驟 1 區段。

讀 research_direction.md + evidence_ledger 最近 10 筆 + CURRENT_FOCUS.md → 產出 ORIENT 摘要。

全自動模式：靜默執行。

---

### 步驟 2：HYPOTHESIZE — 選取假設 `[Review]`

**步驟 2.0：Tombstone 檢查（自動）**

選取假設前，先掃描 evidence_ledger 中所有 NEGATIVE/NO-GO 記錄，阻擋與已失敗假說重複的新假設：

```bash
python3 -c "
import json
tombstones = []
try:
    with open('research/autoresearch/evidence_ledger.jsonl') as f:
        for line in f:
            if not line.strip(): continue
            rec = json.loads(line.strip())
            if rec.get('verdict','').upper() in ('NEGATIVE','NO-GO'):
                tombstones.append({
                    'id': rec.get('hypothesis_id','?'),
                    'hyp': rec.get('hypothesis','')[:60],
                    'track': rec.get('pipeline_track',''),
                    'conditions': rec.get('conditions',{})
                })
except FileNotFoundError:
    pass
if tombstones:
    print(f'[TOMBSTONE] {len(tombstones)} failed hypotheses:')
    for t in tombstones[-10:]:
        print(f'  x {t[\"id\"]}: {t[\"hyp\"]} ({t[\"track\"]})')
    print('New hypothesis must differ in premise, conditions, or method.')
else:
    print('[TOMBSTONE] No failed hypotheses on record.')
"
```

若候選假設與 tombstone 的**方向相同且前提條件未改變** → 標記 `BLOCKED_BY_TOMBSTONE` 並跳過。
復活條件：新數據、新方法、或前提條件明確改變（須標註差異）。

**步驟 2.1：選取假設**

從 `hypothesis_queue.json` 選 status=pending 且 priority 最高者。

```bash
python3 -c "
import json
with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)
pending = [h for h in queue if h['status'] == 'pending']
if not pending:
    print('[QUEUE EMPTY] 需注入新假設')
else:
    top = sorted(pending, key=lambda x: -x['priority'])[0]
    print(json.dumps(top, ensure_ascii=False, indent=2))
"
```

**互動模式**：展示假設 + tombstone 比對結果，等用戶確認或選擇其他。
**全自動模式**：自動選最高優先級（排除 BLOCKED_BY_TOMBSTONE），記錄選擇原因後繼續。

佇列為空 → 依 evidence_ledger 自動生成 3 候選，互動模式下等用戶選擇，全自動模式下選第一個。

---

### 步驟 3：DESIGN TEST — 設計測試 `[Gate(ML3) / Review(ML1-2)]`

詳見 `references/SCALE_LADDER.md`。
驗證閾值與路徑建議：參考 `/validation-protocol` skill（L1→L2→L3→L4 漸進驗證）。

依假設修改層級決定方式（ML1 Python 篩選 / ML2 Python 特徵 / ML3 C++）。

**互動模式**：展示修改計劃。ML3 必須等用戶「確認修改 C++」。
**全自動模式**：ML1-2 自動執行。**ML3 仍為 Hard Gate** — 必須暫停等確認。

修改前必須 Git commit baseline。

---

### 步驟 4：EXECUTE — 執行測試 `[FYI]`

```bash
CYCLE_ID=$(date +%Y%m%d_%H%M%S)
mkdir -p research/autoresearch/cycles/${CYCLE_ID}
```

依 track 執行 benchmark（S2/S3 自動啟用平行模式，詳見 `references/SCALE_LADDER.md`）。

ML3 → 先 `cd build && make -j$(nproc)` 再 benchmark。

全自動模式：靜默執行，完成後直接進 RECORD。

**Manifest 追蹤**：若 `research/{project}/manifest.yaml` 存在，每步完成後更新 artifacts 區段：
```bash
# 將新產出的 scripts/figures/data 路徑追加到 manifest.yaml artifacts
# 同時更新 project.status = "executing"
```

---

### 步驟 5：RECORD — 記錄結果 `[FYI]`

1. 從 result.txt 提取 F1/TP/FP/FN
2. 從 evidence_ledger 取 baseline
3. **自動追加** evidence_ledger.jsonl（不再需要手動）：

```bash
python3 -c "
import json
from datetime import datetime
record = {
    'cycle_id': '${CYCLE_ID}',
    'hypothesis_id': '[H_ID]',
    'hypothesis': '[假設摘要]',
    'pipeline_track': '[track]',
    'datasets_tested': ['[DATASET]'],
    'scale': '[scale]',
    'baseline_f1': {'[DATASET]': 0.0},
    'result_f1': {'[DATASET]': 0.0},
    'delta_f1': {'[DATASET]': 0.0},
    'delta_tp': {'[DATASET]': 0},
    'delta_fp': {'[DATASET]': 0},
    'delta_fn': {'[DATASET]': 0},
    'verdict': '[verdict]',
    'research_potential': '[potential]',
    'mechanism_clarity': '[clarity]',
    'conditions': {
        'mode': '[paired/TO]',
        'samples': 0,
        'metric': '[AUC/F1/effect_size]',
        'method': '[method_description]'
    },
    'conclusion_stability': 0.0,
    'superseded_by': None,
    'valid_until': None,
    'timestamp': datetime.now().isoformat(),
    'artifacts_path': 'research/autoresearch/cycles/${CYCLE_ID}/'
}
with open('research/autoresearch/evidence_ledger.jsonl', 'a') as f:
    f.write(json.dumps(record, ensure_ascii=False) + '\n')
"
```

4. 更新 hypothesis_queue.json 狀態 → pending_verdict
5. Git commit result

**Phase 2 Characterization Record（非 F1 型研究）**：

若研究不涉及 F1 delta（如 subclone characterization），使用以下格式：

```bash
python3 -c "
import json
from datetime import datetime
record = {
    'cycle_id': '${CYCLE_ID}',
    'type': 'characterization',
    'hypothesis': '[假說陳述]',
    'method': '[統計方法]',
    'result': {
        'effect_size': '[效應量]',
        'p_value': '[p值]',
        'cross_sample': '[N/7]',
        'confound_control': '[控制方法]'
    },
    'verdict': '[positive/negative/conditional]',
    'project': '[project_name]',
    'timestamp': datetime.now().isoformat(),
    'artifacts_path': 'research/${PROJECT}/'
}
with open('research/autoresearch/evidence_ledger.jsonl', 'a') as f:
    f.write(json.dumps(record, ensure_ascii=False) + '\n')
"
```

---

### 步驟 6：PRESENT — 呈現結果 `[Review]`

標準呈現格式：

```
╔═══════════════════════════════════════════════════════╗
║ 研究迴圈結果  CYCLE: [CYCLE_ID]                       ║
║ 假設 [H_ID]: [假設摘要]                               ║
╠═══════════════════════════════════════════════════════╣
║ 資料集        基線F1   本輪F1   delta   TP  FP  FN   ║
║ [DATASET]     [X.XX]  [X.XX]  [±X.XX] [N] [N] [N]   ║
╠═══════════════════════════════════════════════════════╣
║ 分析: [why positive/negative/mixed]                   ║
║ 建議: [keep/reject/annotation/escalate/pivot]         ║
╚═══════════════════════════════════════════════════════╝
```

**步驟 6.1：視覺化驗證**（若涉及 TP/FP 區分）：
- 前 5 個 AUC 最高 region 的 ISM 甲基化熱圖
- 前 3 個誤分類 case 的特徵分布圖
- 若可用，提供 IGV 截圖路徑

**互動模式**：展示後暫停，等用戶確認結論。
**全自動模式**：自動判定 verdict（delta > 0 → keep，delta ≤ 0 → reject），記錄後繼續。

---

### 步驟 7：FEEDBACK — 處理回饋 `[Gate]`

**互動模式**：等用戶輸入 `keep` / `reject` / `annotation` / `escalate` / `pivot [方向]`。

**全自動模式**：自動執行以下規則：
- delta > 0 → `keep`，同類假設 priority +15
- delta ≤ 0 → `reject`，同方向假設 priority -30
- delta > 0.002 且 scale=pilot → 自動 `escalate` 到 medium

動作執行後 Git commit。

### 自動推進觸發（步驟 7 後自動判斷）

| 條件 | 動作 |
|------|------|
| 連續 3 輪 negative 同 track | 建議 `/pivot-direction` |
| hypothesis_queue 為空 | 提示 `/inject-hypothesis` 或 `/research-ideation` |
| positive + scale=pilot | 建議升級 medium |
| verdict=concluded | 自動更新 MEMORY.md |
| ≥2 個 combination_ready 假設 | 觸發組合測試（原則 C） |

---

## 假說類型與優先分數

| 類型 | 基礎分 | 加分 |
|------|--------|------|
| `rule_change` | 70 | FP 匹配 +20；OBSERVE SIGNAL +15 |
| `param_combo` | 60 | 前輪同方向 delta > +0.001 +15 |
| `feature_hypothesis` | 65 | W1/W2/W3 設計弱點 +25 |
| `literature_feature` | 75 | IF > 10 +10 |
| `cpp_improvement` | 80 | 連 3 輪 rule_change 無效 +5；W1/W2/W3 +25 |

**新增假說類型**：在上表新增行（類型名 + 基礎分 + 加分條件），並同步在 `.claude/skills/validation-protocol/SKILL.md`「假說類型 → 驗證路徑建議」表新增對應的 L1-L4 路徑。

---

## ISM 使命追蹤（PRESENT 必含）

| 指標 | 計算 |
|------|------|
| 高關連甲基位點 | CramersV > 0.3 的 region 數 |
| Subclone 位點 | VerificationClass = Strong/Subclone 的 variant 數 |
| ISM 獨有貢獻 | ISM_call=keep 且 GQ < threshold 但 QS > 60 |

---

## 注意事項

1. 每輪必建 cycle artifact 目錄
2. 保留所有觀察與推論（含 negative）
3. TP/FP/FN 絕對數字必須記錄
4. 用戶回饋優先於自動判斷
5. ML3 修改之前，互動/全自動模式均必須等用戶確認
