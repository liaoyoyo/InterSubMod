# OBSERVE 觀察協議

## 步驟 0：OBSERVE — 觀察數據（每輪必做）

**目的**：先從真實數據中發現現象，避免盲目測試假設。

**確認當前研究 track 與資料集**：
先讀 `research/autoresearch/research_direction.md` 取得當前 track 與目標資料集。

**執行觀察指令**：

```bash
# 1. 確認最新的 per_region_features 位置
ls scripts/outputs/features/ 2>/dev/null || ls /big8_disk/liaoyoyo2001/InterSubMod/experiments/outputs/ 2>/dev/null | head -20
```

若有 `per_region_features.tsv.gz` 或等效特徵檔：

```bash
# 2. 基礎 TP/FP/FN 特徵分佈統計
ls scripts/analyze_*.py 2>/dev/null

python3 scripts/research_common.py --observe --dataset HCC1395_5kHz_TO \
  --output research/autoresearch/cycles/CYCLE_ID/feature_obs.txt 2>/dev/null \
  || echo "[OBSERVE] 需要手動分析或腳本尚未支援自動觀察"
```

## 觀察記錄格式

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
  → 將「[特徵] threshold 篩選 FP」加入 hypothesis_queue.json（priority=75）
```

## 訊號強度分類

決定測試優先順序：

| 強度 | FP/TP 比值 | 動作 |
|------|-----------|------|
| `STRONG SIGNAL` | > 2x 或 < 0.5x | 優先測試 |
| `mild signal` | 1.5-2x 或 0.5-0.67x | 次要測試 |
| `noise` | 接近 1.0 | 跳過或降優先 |

## 步驟 1：ORIENT — 定向

**讀取以下三個檔案**（依序，不可省略）：

```bash
# 1. 研究方向
cat research/autoresearch/research_direction.md

# 2. 最近 10 條歷史記錄
tail -10 research/autoresearch/evidence_ledger.jsonl 2>/dev/null | \
  python3 -c "import sys,json; [print(json.dumps(json.loads(l), ensure_ascii=False, indent=2)) for l in sys.stdin]"

# 3. 當前焦點
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
⚠ 若 N >= 3 且 track 為 paired_*: 建議切換到 TO track
```
