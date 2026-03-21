---
name: research-loop
description: InterSubMod 半自動研究迴圈。每輪執行「觀察→定向→假設→設計→執行→記錄→呈現→回饋」八步驟，驅動甲基化特徵探索與 F1 提升。適用 paired_full / paired_pileup / TO 三條 pipeline track。
allowed-tools: Read, Write, Bash, Glob, Grep
user-invocable: true
---

<!--
InterSubMod AutoResearch Loop Skill
版本: 1.0.0 (2026-03-22)
設計原則:
  - 跨模型可用：所有步驟使用具體命令，不依賴特定模型推理
  - 觀察優先：數據分析在假設生成之前
  - 小步升規模：pilot → medium → full
  - TO 優先：paired 空間過小時自動建議切換
  - Python 先於 C++：修改層級明確
  - 全記錄：每輪建立 cycle artifact
-->

# InterSubMod 研究迴圈 (Research Loop)

## 此技能的目的

本技能實作 autoresearch 風格的半自動研究迴圈，協助 AI 模型與人類研究者共同：
1. 觀察甲基化特徵中的 TP/FP/FN 模式差異
2. 生成可測試的假設（filter 規則、新特徵、C++ 改進）
3. 執行 benchmark 並記錄 delta F1 與 TP/FP/FN 變化
4. 根據人類回饋調整探索方向
5. 追蹤 InterSubMod 核心使命：發現高關連甲基位點與 subclone 結構

## 觸發時機

當用戶說以下任何一種：
- 「開始研究迴圈」/ 「執行研究」/ 「下一輪假設」
- 「research loop」/ 「autoresearch」
- 「測試新假設」/ 「探索特徵」

## 前置檢查（每輪必做）

在開始八步驟之前，確認以下檔案存在：

```
research/autoresearch/hypothesis_queue.json   ← 必須存在，可為空陣列 []
research/autoresearch/evidence_ledger.jsonl   ← 可為空
research/autoresearch/research_direction.md   ← 必須存在
```

若 `hypothesis_queue.json` 不存在，提示用戶執行 `/inject-hypothesis` 先注入初始假設。

---

## 八步驟執行流程

### 步驟 0：OBSERVE — 觀察數據（每輪必做）

**目的**：先從真實數據中發現現象，避免盲目測試假設。

**確認當前研究 track 與資料集**：
先讀 `research/autoresearch/research_direction.md` 取得當前 track 與目標資料集。

**執行觀察指令**：

```bash
# 1. 確認最新的 per_region_features 位置
# 通常在 experiments/outputs/ 或 scripts/ 執行後的輸出目錄
ls scripts/outputs/features/ 2>/dev/null || ls /big8_disk/liaoyoyo2001/InterSubMod/experiments/outputs/ 2>/dev/null | head -20
```

若有 `per_region_features.tsv.gz` 或等效特徵檔：

```bash
# 2. 基礎 TP/FP/FN 特徵分佈統計（使用現有分析腳本）
# 先找可用的分析腳本
ls scripts/analyze_*.py 2>/dev/null

# 若有 research_common.py，嘗試：
python3 scripts/research_common.py --observe --dataset HCC1395_5kHz_TO \
  --output research/autoresearch/cycles/CYCLE_ID/feature_obs.txt 2>/dev/null \
  || echo "[OBSERVE] 需要手動分析或腳本尚未支援自動觀察"
```

**觀察記錄格式**（手動撰寫或腳本輸出）：

```
=== OBSERVE 結果 ===
資料集: [DATASET_NAME]  Track: [TO/paired_full/paired_pileup]
---
TP 特徵側寫（N=[TP_COUNT]）:
  CramersV:     中位數=[X]  Q1=[X]  Q3=[X]
  HPP:          中位數=[X]  Q1=[X]  Q3=[X]
  VAF:          中位數=[X]  Q1=[X]  Q3=[X]
  AlleleDelta:  中位數=[X]  Q1=[X]  Q3=[X]
  Quality_Score:中位數=[X]  Q1=[X]  Q3=[X]

FP 特徵側寫（N=[FP_COUNT]）:
  CramersV:     中位數=[X]  Q1=[X]  Q3=[X]
  HPP:          中位數=[X]  Q1=[X]  Q3=[X]  ← 若 >2x TP，標記 [SIGNAL]
  VAF:          中位數=[X]
  AlleleDelta:  中位數=[X]
  Quality_Score:中位數=[X]

FN 特徵側寫（N=[FN_COUNT]）:
  CramersV:     中位數=[X]
  HPP:          中位數=[X]
  VAF:          中位數=[X]

特殊關聯觀察:
  [若 HPP 與 AlleleDelta 相關性 > 0.5] → 標記「HPP-AlleleDelta 協同」
  [若 FP_median/TP_median > 2.0 for 某特徵] → 標記「[特徵] 具判別力」
  [若 FN_median 接近 FP_median] → 標記「FN/FP 混淆風險」

自動生成假設（若發現 SIGNAL）:
  → 將「[特徵] threshold 篩選 FP」加入 hypothesis_queue.json（priority=75，status=pending）
```

**若沒有現成特徵檔**：
記錄「本輪無 OBSERVE 數據，直接使用 evidence_ledger 的歷史觀察」並繼續。

---

### 步驟 1：ORIENT — 定向

**讀取以下三個檔案**（依序，不可省略）：

```bash
# 1. 研究方向（當前焦點）
cat research/autoresearch/research_direction.md

# 2. 最近 10 條歷史記錄
tail -10 research/autoresearch/evidence_ledger.jsonl 2>/dev/null | \
  python3 -c "import sys,json; [print(json.dumps(json.loads(l), ensure_ascii=False, indent=2)) for l in sys.stdin]" 2>/dev/null \
  || tail -10 research/autoresearch/evidence_ledger.jsonl

# 3. 當前焦點（快速瀏覽前 80 行）
head -80 docs/CURRENT_FOCUS.md 2>/dev/null
```

**產出 ORIENT 摘要**：

```
=== ORIENT 摘要 ===
當前研究 track: [TO / paired_full / paired_pileup]
當前焦點: [來自 research_direction.md]
最近 3 輪結果:
  輪次 [N-2]: H[ID] [假設摘要] → delta=[+/-X.XXXX] ([verdict])
  輪次 [N-1]: H[ID] [假設摘要] → delta=[+/-X.XXXX] ([verdict])
  輪次 [N]:   H[ID] [假設摘要] → delta=[+/-X.XXXX] ([verdict])

連續無效輪數（|delta| < 0.001）: [N] 輪
⚠ 若 N >= 3 且 track 為 paired_*: 建議切換到 TO track（FP 空間更大）
```

---

### 步驟 2：HYPOTHESIZE — 選取假設

**從 hypothesis_queue.json 選取**：

```bash
python3 -c "
import json
with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)
# 選 status=pending 且 priority 最高者
pending = [h for h in queue if h['status'] == 'pending']
if not pending:
    print('[QUEUE EMPTY] 佇列為空，需注入新假設')
else:
    top = sorted(pending, key=lambda x: -x['priority'])[0]
    print(json.dumps(top, ensure_ascii=False, indent=2))
"
```

**呈現選取的假設**：

```
=== 選取假設 ===
ID: H[XXX]
類型: [rule_change / param_combo / feature_hypothesis / literature_feature / cpp_improvement]
優先分數: [N]
假設: [假設內容]
來源: [來源]
目標 track: [track]
目標資料集: [datasets]
預計 scale: [pilot / medium / full]
修改層級: Tier [1/2/3]（[Python 腳本 / Python 特徵腳本 / C++ 核心]）
```

**若佇列為空**：
依 evidence_ledger 的空白區域自動生成 3 個候選假設，呈現給用戶選擇：

```
[佇列為空] 根據歷史紀錄，建議以下 3 個新方向：
1. [假設A] — 依據: [evidence_ledger 中的觀察]
2. [假設B] — 依據: [尚未測試的特徵組合]
3. [假設C] — 依據: [設計弱點 W1/W2/W3]
請確認後注入佇列（使用 /inject-hypothesis）
```

---

### 步驟 3：DESIGN TEST — 設計測試

依假設的 `tier` 欄位決定修改方式：

#### Tier 1：Python 篩選腳本（優先）

修改目標：`scripts/evaluate_rescue_with_methylation.py` 或對應的分析腳本

```
具體修改說明：
  [描述要修改的規則/參數/閾值]
  [說明修改的位置（函數名/行號）]
  [說明修改前後的邏輯差異]

執行命令：
  [具體的 bash 命令，可直接複製貼上執行]
```

**修改前必須**：
1. 讀取目標腳本的相關片段（使用 Read 工具）
2. 確認修改不會破壞其他功能
3. 先在沙盒中測試語法正確性

#### Tier 2：Python 特徵提取腳本

修改目標：`scripts/research_common.py` 或新建分析腳本

**規則**：必須先有 Tier 1 確認有效的特徵，才能升到 Tier 2 提取新欄位。

#### Tier 3：C++ 核心（最後手段）

**升級條件（全部滿足才執行）**：
- (a) Tier 1/2 在 S2 Medium scale 驗證 delta_f1 > 0.002
- (b) 需要全基因組或多樣本高效能執行
- (c) 人類明確確認：「請修改 C++」

若條件未全部滿足，呈現修改建議但等待人類確認後才動手：
```
[Tier 3 建議] 建議修改 src/core/[file].cpp:
  位置: [函數名:行號區間]
  改動: [說明]
  預期效益: [說明]
⚠ 需要您確認後才會執行修改（輸入「確認修改 C++」）
```

**Scale Ladder 選擇**：

| scale | 條件 | 資料集 | 預期時間 |
|---|---|---|---|
| S1 pilot | 任何新假設預設 | HCC1395 only（5kHz 或 DORADO） | ~2-5 min |
| S2 medium | pilot 通過（delta > 0，方向明確）| HCC1395 + 1-2 個不同 FP 類型的樣本 | ~10 min |
| S3 full | medium 一致（全樣本同方向）| HCC1395_5kHz + HCC1395_DORADO + COLO829 + H1437 + H2009 + HCC1937 | ~30-60 min |

---

### 步驟 4：EXECUTE — 執行測試

**生成 CYCLE_ID**（格式：YYYYMMDD_HHMMSS）：

```bash
CYCLE_ID=$(date +%Y%m%d_%H%M%S)
mkdir -p research/autoresearch/cycles/${CYCLE_ID}
echo "CYCLE_ID: ${CYCLE_ID}"
```

**建立 cycle 紀錄**：

```bash
# 儲存本輪假設說明
cat > research/autoresearch/cycles/${CYCLE_ID}/hypothesis.md << 'EOF'
假設 ID: [H_ID]
假設內容: [假設文字]
來源: [來源]
Pipeline Track: [track]
Scale: [S1/S2/S3]
修改層級: Tier [N]
預期方向: [正增益 / 移除 FP / 增加 TP]
EOF
```

**執行 benchmark**：

依 track 選擇對應命令：

```bash
# TO track（5kHz）
./scripts/run_batch_vcf_analysis.sh --mode HCC1395_5kHz_TO 2>&1 | \
  tee research/autoresearch/cycles/${CYCLE_ID}/result.txt

# TO track（DORADO）
./scripts/run_batch_vcf_analysis.sh --mode HCC1395_DORADO_TO 2>&1 | \
  tee research/autoresearch/cycles/${CYCLE_ID}/result.txt

# Paired track
./scripts/run_pure_research_round.sh 2>&1 | \
  tee research/autoresearch/cycles/${CYCLE_ID}/result.txt

# 快速單點驗證（chr19 only）
./scripts/run_vcf_all_snv.sh --mode chr19-verification 2>&1 | \
  tee research/autoresearch/cycles/${CYCLE_ID}/result.txt
```

若 Tier 3 C++ 修改：先編譯再執行：

```bash
cd build && make -j$(nproc) 2>&1 | tee research/autoresearch/cycles/${CYCLE_ID}/build.txt
cd .. && ./scripts/run_batch_vcf_analysis.sh ... 2>&1 | tee research/autoresearch/cycles/${CYCLE_ID}/result.txt
```

---

### 步驟 5：RECORD — 記錄結果

**從 result.txt 提取 F1/TP/FP/FN**：

```bash
# 嘗試自動提取（依輸出格式調整 grep 樣式）
grep -E "F1|precision|recall|TP:|FP:|FN:" research/autoresearch/cycles/${CYCLE_ID}/result.txt | head -30
```

**從 evidence_ledger.jsonl 取得 baseline**（同 dataset 最近一筆）：

```bash
python3 -c "
import json
dataset = 'HCC1395_5kHz_TO'  # 替換為本輪 dataset
results = []
with open('research/autoresearch/evidence_ledger.jsonl') as f:
    for line in f:
        r = json.loads(line)
        if dataset in r.get('datasets_tested', []):
            results.append(r)
if results:
    latest = results[-1]
    print('Baseline F1:', latest['result_f1'].get(dataset))
    print('Baseline TP/FP/FN from:', latest['cycle_id'])
else:
    print('No baseline found for', dataset)
    print('Use initial baseline from Pipeline Track table')
"
```

**寫入 evidence_ledger.jsonl**（append，不覆蓋）：

```bash
python3 -c "
import json
from datetime import datetime

record = {
    'cycle_id': '${CYCLE_ID}',
    'hypothesis_id': '[H_ID]',
    'hypothesis': '[假設摘要，60字以內]',
    'pipeline_track': '[TO/paired_full/paired_pileup]',
    'datasets_tested': ['[DATASET]'],
    'scale': '[pilot/medium/full]',
    'baseline_f1': {'[DATASET]': [BASELINE_F1]},
    'result_f1': {'[DATASET]': [RESULT_F1]},
    'delta_f1': {'[DATASET]': [DELTA_F1]},
    'delta_tp': {'[DATASET]': [DELTA_TP]},
    'delta_fp': {'[DATASET]': [DELTA_FP]},
    'delta_fn': {'[DATASET]': [DELTA_FN]},
    'key_observations': '[本輪重要觀察，來自 result.txt 分析]',
    'feature_correlations_noted': '[特徵間特殊關聯，若本輪有發現]',
    'verdict': '[positive_pilot/negative/dataset_specific/annotation_only]',
    'human_decision': '',
    'human_notes': '',
    'timestamp': datetime.now().isoformat(),
    'tier_used': [1/2/3],
    'artifacts_path': 'research/autoresearch/cycles/${CYCLE_ID}/'
}

with open('research/autoresearch/evidence_ledger.jsonl', 'a') as f:
    f.write(json.dumps(record, ensure_ascii=False) + '\n')
print('Recorded:', record['cycle_id'])
"
```

**更新 hypothesis_queue.json 狀態**（將 H_ID 的 status 改為 running→pending_verdict）：

```bash
python3 -c "
import json
with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)
for h in queue:
    if h['id'] == '[H_ID]':
        h['status'] = 'pending_verdict'
        break
with open('research/autoresearch/hypothesis_queue.json', 'w') as f:
    json.dump(queue, f, ensure_ascii=False, indent=2)
print('Updated H_ID status')
"
```

---

### 步驟 6：PRESENT — 呈現給人類

**標準呈現格式**（每輪固定輸出此格式）：

```
╔══════════════════════════════════════════════════════════════╗
║ 研究迴圈結果  CYCLE: [YYYYMMDD_HHMMSS]                       ║
║ 假設 [H_ID]: [假設摘要]                                       ║
╠══════════════════════════════════════════════════════════════╣
║ 資料集           基線 F1    本輪 F1    delta F1  TP  FP  FN ║
║ [DATASET_1]      [0.XXXX]  [0.XXXX]  [+/-0.XXXX] [N] [N] [N] ║
║ [DATASET_2]      [0.XXXX]  [0.XXXX]  [+/-0.XXXX] [N] [N] [N] ║
╠══════════════════════════════════════════════════════════════╣
║ 分析:                                                         ║
║   [說明為什麼結果是正/負/混合]                                 ║
║   [說明受影響的主要是哪類 TP/FP/FN]                           ║
║   [若結果跨樣本不一致，說明可能原因]                           ║
╠══════════════════════════════════════════════════════════════╣
║ 特徵關聯觀察:                                                  ║
║   [若本輪發現特徵間的特殊關聯，在此記錄]                       ║
║   [例如：HPP 與 AlleleDelta 協同出現在 FP-B 中]               ║
╠══════════════════════════════════════════════════════════════╣
║ InterSubMod 使命指標:                                          ║
║   高關連甲基位點（CramersV > 0.3）: [N] regions               ║
║   Subclone 位點（Strong/Subclone class）: [M] variants        ║
║   ISM 獨有貢獻（非 GQ/QUAL 可解釋）: [K] variants            ║
╠══════════════════════════════════════════════════════════════╣
║ 建議決策:                                                      ║
║   [keep / reject / annotation_only / escalate_medium]          ║
║   [理由]                                                       ║
╠══════════════════════════════════════════════════════════════╣
║ 等待您的指令:                                                  ║
║   keep         → 採納此規則，下輪升 scale 或繼續下一假設      ║
║   reject       → 放棄此假設，標記失敗方向                      ║
║   annotation   → 降階為 annotation only（不入核心規則）        ║
║   escalate     → 升 scale（pilot→medium 或 medium→full）      ║
║   pivot [說明] → 放棄當前方向，注入新假設                      ║
╚══════════════════════════════════════════════════════════════╝
```

---

### 步驟 7：FEEDBACK — 處理人類回饋

依用戶輸入執行對應動作：

**`keep`（採納）**：

```bash
python3 -c "
import json
with open('research/autoresearch/hypothesis_queue.json') as f:
    queue = json.load(f)
for h in queue:
    if h['id'] == '[H_ID]':
        h['status'] = 'adopted'
# 同類假設優先分數 +15
for h in queue:
    if h['status'] == 'pending' and h['type'] == '[ADOPTED_TYPE]':
        h['priority'] = min(100, h['priority'] + 15)
with open('research/autoresearch/hypothesis_queue.json', 'w') as f:
    json.dump(queue, f, ensure_ascii=False, indent=2)
"
# 更新 evidence_ledger 的 human_decision = 'keep'
```

**`reject`（拒絕）**：

```bash
# H_ID status = 'rejected'
# 同方向假設 priority -= 30（最低 10）
```

**`annotation`（降階 annotation）**：

```bash
# H_ID status = 'adopted_annotation'
# 在 research_direction.md 加入：[特徵] 已驗證為 dataset-specific annotation
```

**`escalate`（升 scale）**：

```bash
# H_ID scale 從 pilot → medium 或 medium → full
# 更新 status = 'pending'（重新測試）
# 自動補充 target_datasets（依 FP 類型選新樣本）
```

**`pivot [新方向說明]`**：

```bash
# 自動呼叫 inject-hypothesis 流程
# 將新方向設定 priority = 90（最高優先）
# 更新 research_direction.md 加入注記
```

**同步更新 evidence_ledger 的 human_decision 和 human_notes**：

```bash
python3 -c "
import json

# 讀取並更新最後一條記錄
records = []
with open('research/autoresearch/evidence_ledger.jsonl') as f:
    for line in f:
        records.append(json.loads(line.strip()))

# 更新最後一條
records[-1]['human_decision'] = '[keep/reject/annotation_only/escalate_medium/escalate_cpp]'
records[-1]['human_notes'] = '[用戶的補充說明]'

# 寫回
with open('research/autoresearch/evidence_ledger.jsonl', 'w') as f:
    for r in records:
        f.write(json.dumps(r, ensure_ascii=False) + '\n')
print('Updated human_decision')
"
```

---

## Pipeline Track 基線表

| Track | 資料集 | 基線 F1 | FP 數量 | FN 數量 | 測試命令關鍵字 |
|---|---|---|---|---|---|
| paired_full | HCC1395_5kHz_paired | 0.8532 | 少 | 少 | `--mode HCC1395_5kHz_paired` |
| paired_pileup | HCC1395_DORADO_paired | 0.8590 | 少 | 少 | `--mode HCC1395_DORADO_paired` |
| TO | HCC1395_5kHz_TO | 0.7127 | 多 | 多 | `--mode HCC1395_5kHz_TO` |
| TO | HCC1395_DORADO_TO | 0.7226 | 多 | 多 | `--mode HCC1395_DORADO_TO` |

**TO 優先規則**：若 paired track 連續 3 輪 `abs(delta_f1) < 0.001`，呈現建議：
```
⚠ 建議切換 TO Track
  - paired 連續 3 輪無效（|delta| < 0.001）
  - TO 有更多 FP 可觀察（paired FP ~數十個 vs TO FP ~數百個）
  - 輸入「切換 TO」或繼續輸入 keep/reject 維持當前 track
```

---

## 假設類型定義與優先分數算法

| 類型 | 說明 | 基礎優先分 | 加分條件 |
|---|---|---|---|
| `rule_change` | 修改 VCF 篩選規則或閾值 | 70 | FP 類型明確匹配 +20；OBSERVE 中有 SIGNAL +15 |
| `param_combo` | 修改 ISM 參數（distance_method, window_bp 等）| 60 | 前輪同方向 delta > +0.001 +15 |
| `feature_hypothesis` | 提取新特徵欄位（需 Tier 2 腳本修改）| 65 | W1/W2/W3 設計弱點直接對應 +25 |
| `literature_feature` | 文獻引導的特徵（附 source 標記）| 75 | 文獻 IF > 10 +10 |
| `cpp_improvement` | C++ 核心算法改進 | 80 | 連續 3 輪 rule_change 無效觸發 +5；W1/W2/W3 +25 |

**設計弱點分（來自 2026-03-17 分析）**：
- W1：HPP 未整合進分類邏輯 → cpp_improvement 假設 +25
- W2：HP 不均時 AlleleDelta 虛高 → feature_hypothesis（AlleleMethDelta 新欄位）+25
- W3：Phase block 邊界未標記 → cpp_improvement +20

---

## 程式碼修改層級規則

```
Tier 1：Python 篩選腳本（每次迭代的預設層）
  目標腳本: scripts/evaluate_rescue_with_methylation.py
           scripts/analyze_gq_methylation_rescue_matrix.py
           scripts/analyze_methylation_rescue_feature_space.py
           scripts/run_batch_vcf_analysis.sh（規則參數部分）
  升級不需條件：直接迭代

Tier 2：Python 特徵提取腳本
  目標腳本: scripts/research_common.py
           新建 scripts/analyze_feature_distributions.py
  升級條件: Tier 1 確認有效特徵，需從原始輸出提取新欄位

Tier 3：C++ 核心
  目標: src/core/*.cpp, include/core/*.hpp
  升級條件（全部需滿足）:
    (a) Tier 1/2 在 S2 Medium 驗證 delta_f1 > 0.002
    (b) 需全基因組或多樣本高效能執行
    (c) 人類輸入「確認修改 C++」
  注意: Tier 3 修改後必須執行 cd build && make -j$(nproc)
```

---

## InterSubMod 使命追蹤

每輪 PRESENT 必須包含以下使命指標（從 result.txt 或輸出統計提取）：

| 指標 | 計算方式 | 說明 |
|---|---|---|
| 高關連甲基位點 | `CramersV > 0.3` 的 region 數 | ISM 發現的強甲基訊號位點 |
| Subclone 位點 | `VerificationClass = Strong 或 Subclone` 的 variant 數 | ISM 主張有 subclone 結構的位點 |
| ISM 獨有貢獻 | `ISM_call = keep` 且 `GQ < threshold 但 Quality_Score > 60` | 非純 caller-quality 的貢獻 |
| Subclone class precision | `Subclone class 中 TP 數 / 全 Subclone class 數` | ISM subclone 判斷的精確度 |

**目標**：追蹤 ISM 是否真正在「發現 subclone 結構」，而非僅過濾 caller noise。

---

## 注意事項

1. **每輪必須建立 cycle artifact 目錄**，即使最後結果為 negative
2. **保留所有觀察與推論**，即使是「無效結果」也記錄原因
3. **TP/FP/FN 絕對數字必須記錄**，不得只記 F1 差值
4. **特徵間關聯若有發現，立即記錄**（可能比 F1 更有價值）
5. **用戶回饋優先於自動判斷**：若用戶說「換方向」，立即更新佇列
6. **Tier 3 修改之前**，必須等待用戶確認，不得自動執行 C++ 修改後編譯
