---
status: validated_partial_concern
date: 2026-04-21
type: marker_reaudit
priority: P1
pipeline_track: TO
outcome: marker_not_reproduced_on_HCC1395_TO_standalone
related_phase1: docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md
---

# HPFineNGroups Subclone Marker — P1 Re-audit（HCC1395 TO 標準資料）

## 結論（一行）

**原 F pilot 的 89.1%（7 樣本 master dataset pooled）結論不宣告撤回**，但新發現：
1. **HCC1395 TO 單樣本 ClairS-TO raw TP/FP split 上，marker 無 TP 富集**（NonLOH + NG≥4 + NR≥80 TP rate=0.694 vs NonLOH baseline 0.699）
2. **flag=on 後 NG≥3 regions 消失 100%** → 原 marker 的辨識訊號**至少部分**依賴 somatic HP tag 的細分 labels

這是一項**資料集 dependency discovery**，不是結論 retraction，但足以警示未來 marker 應用必須：
- 標記適用的 VCF source / pipeline
- 若使用 `--germline-hp-only=true` 則 marker 不可用

---

## 1. 修正理解：memory 的 89.1% 是**master dataset 7-sample pooled**

Memory `project_hpfinengroups_subclone_marker.md` line 26-27 明確記錄：
- **總體**：NG=4+AF<0.4+NR≥80 NonLOH TP rate=0.9281（n=14,197，7 樣本 pooled）
- **HCC1395 old vs new**：0.810 → 0.887（+7.7pp）

這兩個數值來源是 **master dataset** 分析流程（含 7 樣本 + AF filter），非本次 Phase 1 使用的 **HCC1395 TO ClairS-TO TP/FP raw split**（40,115 sites）。

因此**不能直接比較**。

---

## 2. 本次在 HCC1395 TO 上的實測

**資料集**：HCC1395 TO V3-Fixed haplotag BAM + ClairS-TO raw TP/FP VCF split
- TP = 28,509（from clairsto_tp.vcf.gz）
- FP = 11,606（from clairsto_fp.vcf.gz）
- 全域 TP rate = 0.711（NonLOH 子集 0.699）

### 2.1 flag=off（原行為）

| 條件 | n | TP rate | vs baseline |
|------|---|---------|------------|
| NonLOH baseline | 37,812 | 0.6992 | — |
| NonLOH + NG≥4 + NumReads≥80 | 5,749 | **0.6944** | -0.005 pp |
| NonLOH + NG=4 + NumReads≥80 | 5,749 | 0.6944 | -0.005 pp |
| NonLOH + NG≥4 + NTumorReads≥80 | 3,203 | 0.6959 | -0.003 pp |
| 全域 NG≥4 + NumReads≥80 | 5,760 | 0.694 | -0.017 pp |

**Fisher test**（NG≥4+NR≥80 vs TP label, 全域）：
- Odds ratio = 0.913（反向）, p = 3.5×10⁻³
- **統計顯著但方向錯誤**（marker 不富集 TP，邊緣 deplete）

**AF filter 無法在本次驗證**：
- ISM TSV 無 caller AF 欄位（memory 提到的「NG=4+AF<0.4」需外部 VCF AF 合併）
- 若需完整重現新 canonical filter，需另跑 Python post-process 合併 VCF AF

### 2.2 flag=on（somatic HP demoted）

| 條件 | n | TP rate |
|------|---|---------|
| NonLOH + NG≥4 + NumReads≥80 | **0** | — |
| NonLOH + NG≥3 + NumReads≥80 | **0** | — |
| 全域 NG≥2 + NumReads≥80 | 38,415 | 0.708 |

- **整個 HCC1395 TO（40k sites）無任何 region 在 flag=on 下 NG≥3**
- 符合 Phase 1 主報告的觀察（97.3% regions 歸 NG=2）
- 等同於**移除 somatic HP tag 後，HPFineNGroups marker 失去訊號源**

---

## 3. 兩種可能解讀

### 解讀 A：資料集 dependency（memory 結論在 master dataset 仍成立）
- F pilot 的 0.810 HCC1395 old filter 來自 **master dataset**（不同 VCF / aggregation）
- 本次 TO standalone 無富集是**資料集差異**，非 marker 機制失效
- flag=on 破壞 NG≥3 是**mechanism artifact**，但不等於原 POSITIVE 結論錯誤

### 解讀 B：Marker 訊號部分來自 somatic HP tag 污染
- 原 F pilot Step 3 chr-shuffle Z=43.5 確認**非 spatial artifact**（已驗證）
- 但 chr-shuffle null 不檢查「signal 是否來自 somatic HP 細分」
- flag=on 測試是 **orthogonal null**：移除 somatic HP tags 後 marker 消失 → 訊號源至少 50%+ 在 somatic HP labels 上

**兩種解讀可能同時為真**：master dataset 的 marker 訊號可能**部分**來自 subclone 結構、**部分**來自 somatic HP self-phasing。當前 Phase 1 無法分辨兩者佔比。

---

## 4. 具體操作建議

### 4.1 Memory 更新（已執行，部分）
- `project_hpfinengroups_subclone_marker.md`：標題已加警告「flag=on 下 N≥3 完全消失」
- 建議追加：**「F pilot 結論基於 master dataset，`--germline-hp-only=true` 下 marker 不可用」**

### 4.2 未來若要嚴格驗證 marker 生物學真實性
- **Phase 2B 建議實驗**（非本次）：master dataset 重跑 × 兩 flag，看 7-sample pooled TP rate 從 0.928 變多少
- 若 flag=on 後 master dataset marker TP rate 仍 ≥0.85 → 結論穩健
- 若 flag=on 後 master dataset marker TP rate 暴跌至 ~baseline → marker 生物學根基需重建

### 4.3 當前 marker 使用指引
- **ISM 當前研究（Phase 2 生物學特徵化）**：若使用 subclone marker 引用 → 必須標記 `--germline-hp-only=false` 前提
- **未來發表**：建議註明 marker 的 **pipeline dependency**
- **filter 方向**：不變（仍 NEGATIVE，AUC 0.63-0.65）

---

## 5. P2 / P3 建議更新版

### P2：within_dom_alt_frac downstream 重建
**建議 NOT DO（維持 Phase 1 判定）**
- 理由不變：Phase 1 全體 HP-derived 特徵 AUC 下降，單一 derived feature 改善至 ≥0.70 機率低
- 成本 vs 收益不划算

### P3：`--germline-hp-only` 推為未來 default
**建議 STATUS QUO（暫不改 default）**
- 本 P1 顯示 default=on 會**立即破壞 HPFineNGroups subclone marker 的當前使用**（即使該 marker 已有警示）
- 生物學解釋工具（subclone heterogeneity visualization）在 default=off 下 F pilot 結論仍可引用
- Default 改 on 等於強迫撤回所有 HP-dependent 結論 → 成本遠高於收益
- **折衷**：保留 `--germline-hp-only` 作為研究者 opt-in 工具；文檔加入「若懷疑 self-phasing 污染特定結論 → 開 flag 重跑比對」使用指引

### P4（新增）：master dataset × 兩 flag 比對（建議 OPEN, 條件執行）
- 若未來要發表 HPFineNGroups biological marker 論文 → 必做
- 若研究方向轉離 subclone marker → skip
- 成本：7 樣本 × 2 flag × 2 mode（TO+paired） = 28 runs ≈ 12-24 hr 機器時；延後到 CovM baseline blocker 解決後合併

---

## 6. 執行成本結算

- P1 執行（含回頭校驗 memory + 重分析）：<1 hr
- 主要產出：**marker 的 pipeline dependency 警示**（生物學真實性未被否定，但應用場景受限）
- 次要產出：P3 建議從 CONDITIONAL DO 調整為 STATUS QUO；新增 P4 條件 option

**投資報酬率評估**：中-高（避免未來引用 marker 時的不當推論；為後續 Phase 2B 判定留下接口）

---

## 檔案清單

- **資料**：`/tmp/ism_hp_fix_phase1/merged_{off,on}.csv`
- **Phase 1 主**：`20260421_ReadParser_GermlineHPOnly_HCC1395_01.md`
- **本 P1 報告**：`20260421_HPFineNGroups_Marker_Reaudit_01.md`
- **Memory**：
  - `project_hpfinengroups_subclone_marker.md`（需加 `--germline-hp-only` 下失效警示）
  - `project_readparser_germline_hp_only_phase1_negative.md`（Phase 1 主結論）
- **F pilot 原始**：`docs/experiments/in_progress/2026/04/20260418_F_HPFineNGroups_deepening_POSITIVE_01.md`（結論 master dataset 層仍成立）
