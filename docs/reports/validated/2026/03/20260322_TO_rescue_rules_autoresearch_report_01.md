<!--
建立時間: 2026-03-22
目標: InterSubMod TO Track Rescue Rules 自動化研究循環完整報告
處理範圍: HCC1395_5kHz_TO 單獨測試 (S1 pilot)，12 個假設，8 個研究輪次
關聯檔案:
  - research/autoresearch/evidence_ledger.jsonl
  - research/autoresearch/hypothesis_queue.json
  - research/autoresearch/cycles/20260322_*/
-->

# InterSubMod TO Rescue Rules 自動化研究報告

> 生成時間：2026-03-22
> 資料集：HCC1395_5kHz_TO (pilot)
> 研究迴圈版本：v1.1.0（含 Git 紀律、單獨測試、組合分析三大原則）

---

## 一、研究背景

**目標**：提升 InterSubMod 在 Tumor-Only (TO) track 的 F1 分數，透過系統性探索 ISM 甲基化特徵作為 rescue 規則。

**基線**：
| 指標 | 值 |
|------|-----|
| 基線 F1 (HCC1395_5kHz_TO) | 0.712697 |
| TP | 28,396 |
| FP | 11,843 |
| FN | 11,051 |
| Truth Total | 39,447 |

**候選池**（ISM-analyzed rescue candidates）：
773 TP (currently FN) + 298 FP (filtered-out) = 1,071 candidates
ISM 覆蓋率：773/11,051 FN = **7.0%**

---

## 二、執行的假設清單與結果

### 2.1 採納的假設（Adopted）

| 假設 ID | 內容 | delta_F1 | 結論 |
|---------|------|----------|------|
| **H011** | QS>=50 alone (無 GQ gate) 作為 TO rescue | **+0.008556** | ✅ adopted |
| **H012** | GQ>=3 alone 作為 TO rescue（自動發現）| **+0.009365** | ✅ adopted — 新最佳 |

### 2.2 Annotation Only（不影響 F1，但作為 QC 標記）

| 假設 ID | 內容 | 關鍵指標 |
|---------|------|---------|
| H008 | PMD>=0.10 作為支援 annotation | PMD+QS 96% 相關，annotation only |
| H009 | CramersV>0.3 在 TO 的 precision | TO=85.4% vs paired=99.72%（跨 track 不穩定）|
| H010 | hp_assign_rate<0.2 作為 Low_Phase_Quality 標記 | hp<0.2 group 46.1% FP rate |

### 2.3 拒絕的假設（Rejected）

| 假設 ID | 拒絕原因 |
|---------|---------|
| H001 | HPP 在 TO 無判別力（TP 52.9% vs FP 52.1% HPP<0.05）|
| H002 | 依賴 H001（失敗）+ AlleleDelta（TO 近 0）|
| H003 | AlleleMethDelta AUROC=0.473（random），FP median > TP median |
| H004 | QS>=50+CramersV>0.15 = QS>=50 的嚴格子集，96%重疊 |
| H005 | VAF<0.08 在 TO 觸發率<1%；AlleleDelta≈0 in TO |
| H007 | 前提（H001 HPP 判別力）失效 |

### 2.4 延後的假設

| 假設 ID | 原因 |
|---------|------|
| H006 | 需要重跑 ISM (window_bp=1000)，資料準備中 |

---

## 三、核心發現

### 3.1 最佳 Rescue 規則（Tier 排名）

```
排名  規則                    delta_F1  rTP   rFP   precision
─────────────────────────────────────────────────────────────
[1]  GQ>=3 [H012]           +0.009365  728  255    74.1%  ← 最高 F1
[2]  GQ>=4                  +0.009146  695  221    75.9%
[3]  GQ>=5                  +0.008840  661  195    77.2%
[4]  QS>=50 [H011]          +0.008556  642  193    76.9%  ← 最高 precision
[5]  GQ>=10+QS>=60 [prev]   +0.005551  [from Phase 2]

理論最大值（rescue all ISM candidates）: +0.009692
```

### 3.2 OBSERVE 關鍵觀察

**TO track 的特殊性**（區別於 paired）：

| 特徵 | TO 中的表現 | Paired 中的表現 | 影響 |
|------|------------|----------------|------|
| HPP | TP=52.9% vs FP=52.1%（無判別力）| 有判別力 | H001/H002 失效 |
| AlleleDelta | 近 0（TP med=0.010, FP med=0.000）| 5.29x ratio | H005/H003 失效 |
| CramersV>0.3 precision | 85.4% | 99.72% | 跨 track 不穩定 |
| HP-based features | 普遍無效（haplotagging 不可靠）| 有效 | 根本性限制 |

### 3.3 組合測試結論

- **H011(QS>=50) ⊂ H012(GQ>=3)**：H011 是 H012 的嚴格子集（無正交性）
- **組合策略**：無正交組合可提升 F1 — GQ>=3 已是此候選池的最優
- **Annotation filters 均降低 F1**：hp>=0.2、VC!=Noise、PMD 過濾均損失 TP

### 3.4 TO HP-based 特徵的根本限制

```
根本原因：TO haplotagging 將大多數 reads 分配到 HP1，
         HP2 嚴重欠採樣，使所有 HP-based 差異特徵失效：
  - AlleleDelta ≈ 0
  - HPMergedDelta ≈ 0 (FP median > TP median)
  - AlleleMethDelta AUROC = 0.473
  - HPP ≈ 0.05 for both TP and FP

→ 在 TO track 中，甲基化特徵的有效信號僅來自 global 特徵
  （Quality_Score, PairwiseMedianDist），不來自 HP-specific 特徵
```

---

## 四、Annotation 體系（ISM 三層架構更新）

```
層級  名稱                      條件                     意義
────────────────────────────────────────────────────────────────────
L1   Caller-First              GQ / QUAL / FILTER       主過濾層（現有）
L2   Methylation-Rescue        GQ>=3 或 QS>=50          IS 在 TO 的 rescue 規則
     └─ High-F1 variant        GQ>=3                    最大 F1 增益 (+0.009365)
     └─ High-Precision variant QS>=50                   最高精確度 (+0.008556)
L3   Annotation / QC           hp_assign_rate, CramersV 品質標記（不影響主規則）
     └─ Low_Phase_Quality      hp_assign_rate < 0.2     46.1% FP rate 警告
     └─ High_Methylation       CramersV > 0.3           85.4% precision 強信號
     └─ Methylation_Support    PMD >= 0.10              corroborates rescue rule
```

---

## 五、研究迴圈執行記錄

| Cycle ID | 假設 | delta_F1 | verdict |
|----------|------|----------|---------|
| 20260322_041120 | H011 (QS>=50) | +0.008556 | positive_pilot → adopted |
| 20260322_044024 | H004 (QS+CV) | +0.000663 | negative → rejected |
| 20260322_044136 | H008 (PMD+GQ) | annotation | annotation_only |
| 20260322_044247 | H009 (CramersV>0.3) | annotation | annotation_only |
| 20260322_044500 | H005 (VAF+AD filter) | 0 | negative → rejected |
| 20260322_044612 | H010 (hp_assign QC) | annotation | annotation_only |
| 20260322_045305 | H003 (AlleleMethDelta) | AUROC=0.473 | negative → rejected |
| 20260322_045801 | **H012 (GQ>=3)** | **+0.009365** | **positive → adopted** |
| 20260322_045905 | COMBO test | max=+0.009365 | no orthogonal combination |

**Git Commits**: 18 commits（包含 pre-change baseline + result + feedback）

---

## 六、下一步行動建議

### 優先級 P1：驗證 H012 在 DORADO TO
```
目標: DORADO TO baseline F1=0.722573
預測: H012 在 DORADO TO 的 rescue opportunity 極少（+0.000476 max from Phase 2）
方法: 使用 DORADO TO rescue_joined_features（已知路徑）
```

### 優先級 P2：延後 H006（window_bp 測試）
```
需要: 重新對 HCC1395_5kHz_TO 跑 ISM，使用 window_bp=1000 vs 5000
BAM: /big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam (需確認大小)
注意: ISM 目前只覆蓋 7% FN pool — 擴大覆蓋率比調整 window 更有影響力
```

### 優先級 P3：擴大 ISM 覆蓋率
```
關鍵缺口: ISM 只分析了 773/11,051 FN 候選（7%）
若覆蓋率提升到 50% → rescue potential 可能達 +0.04 或更多
方法: 降低 ISM region selection 門檻，或對所有 FN 候選跑 ISM
```

### 優先級 P4：Medium Scale 驗證（HCC1395 DORADO + COLO829）
```
條件: H012 需在 S2 medium scale 確認 (DORADO TO + 1-2 其他樣本)
Note: DORADO TO 的 rescue 空間很小（Phase 2 確認），
      重點在確認 GQ>=3 不在 DORADO TO 引入大量 FP
```

---

## 七、研究潛力評估

| 方向 | 潛力 | 說明 |
|------|------|------|
| GQ>=3 escalation (DORADO TO) | medium | 確認跨樣本穩定性 |
| 擴大 ISM 覆蓋率 | **high** | 目前只覆蓋 7% FN，是最大改善空間 |
| window_bp 調整 (H006) | low | 機制不清，需 re-run |
| HP-based 特徵（TO）| exhausted | 所有 HP 特徵在 TO 無效，根本性限制 |
| CramersV annotation | medium | 有 annotation 價值，但不改 F1 |

---

## 八、方法論貢獻

本輪研究確立了 **研究迴圈 v1.1.0** 的四大核心紀律：

1. **Git Commit 紀律**：pre-change + result + feedback = 每輪 3 commits
2. **單獨測試**：一次改一個變數，確認 mechanism clarity
3. **組合分析**：有效結果後才組合，確認是否正交
4. **研究潛力分類**：high/medium/low/exhausted 精確追蹤

這些紀律使本輪研究能夠：
- 快速排除 8 個無效假設（H001/H002/H003/H004/H005/H007）
- 發現 H012（GQ>=3），比原假設 H011 多 +0.000809
- 確認 TO HP-based 特徵的根本性限制

---

*報告由 InterSubMod AutoResearch Loop v1.1.0 自動生成*
*研究者：Claude Sonnet 4.6 + 人類研究員*
