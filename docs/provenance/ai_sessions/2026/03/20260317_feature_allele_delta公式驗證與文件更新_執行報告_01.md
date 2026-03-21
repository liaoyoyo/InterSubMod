# feature_allele_delta 公式驗證與文件更新 執行報告

<!--
建立時間: 2026-03-17 00:00
目標: 驗證 InterSubMod feature_allele_delta 的真實計算公式，解釋 TP-4 符號差異，並更新相關文件
處理範圍: src/core/LabelTest.cpp 原始碼查驗 + 三份文件更新
關聯檔案:
  - src/core/LabelTest.cpp
  - docs/reports/validated/2026/03/20260317_甲基分類方法有效性觀察與InterSubMod改進重點_01.md
  - docs/references/manual/20260317_甲基位點觀察報告生成技能_01.md
  - docs/references/manual/20260317_case_obs_report_FP_B1_範例_01.md
-->

## 對話資訊

| 項目 | 內容 |
|------|------|
| 日期 | 2026-03-17 |
| 主要目的 | 查驗 feature_allele_delta 計算公式，解釋 TP-4 「符號相反」異常，更新文件中的錯誤假說 |
| AI 模型 | Claude Sonnet 4.6 |

## 對話背景

前次對話（同日）建立了甲基分類方法有效性觀察文件，其中標記 TP-4（chr16:35118902）的 `feature_allele_delta=+0.023`（TSV）與案例分析手算值 `ALT_meth−REF_meth=−0.198` 符號相反為「🔴 最嚴重異常」。本次對話查驗原始碼以確認真相。

## 關鍵決策

| 決策 | 原因 | 影響 |
|------|------|------|
| 查驗 `src/core/LabelTest.cpp` 的 `compute_delta()` 與 `compute_group_distances()` | 原始碼是唯一可信來源，不應憑假說修改 Pipeline | 確認公式定義，解除誤判 |
| 將 TP-4「符號相反」從「bug」重新分類為「兩個指標量不同面向」 | 驗證後確認不是 bug，是理解錯誤 | 移除誤導的 🔴 行動項目 |
| 建議新增 `AlleleMethDelta`（raw 甲基差）欄位與現有距離 delta 並列 | 兩個指標互補；距離 delta 在 ALT-HP 共線時可被污染，raw 甲基差提供獨立參考 | 改進建議調整為新增欄位，而非修正現有計算 |

## 產出成果

### 原始碼查驗結果

**`feature_allele_delta` 真實公式**（`src/core/LabelTest.cpp:compute_delta()`）：

```
feature_allele_delta = between_mean_distance(ALT, REF) − within_mean_distance(ALT+REF)
```

- Group 0 = ALT reads（AltSupport::ALT）
- Group 1 = REF reads（AltSupport::REF）
- 距離來源：read-read 甲基化距離矩陣（BERNOULLI 或其他 metric）
- **這是 PERMANOVA-style 距離空間分離指標，不是 raw 甲基差**

**TP-4 符號差異的解釋**：
- TSV `+0.023`：距離空間中 ALT/REF 群各自緊密、兩群互斥 → between > within → 正值 ✅
- 案例分析 `−0.198`：ALT reads 的 raw 平均甲基化低於 REF reads（ALT=低甲基位點）
- 兩者量的是不同面向，可以同時成立，無矛盾

**FP-B1 高 AlleleDelta（+0.058）虛高機制確認**：
- ALT reads 全在 HP1 側 → within_ALT 小（同質緊密）
- REF reads 跨 HP1/HP2 → within_REF 大（異質）
- between(ALT, REF) 高 → delta 虛高，但是 HP imbalance 造成，非 somatic allele 效應

### 修改檔案

- `docs/reports/validated/2026/03/20260317_甲基分類方法有效性觀察與InterSubMod改進重點_01.md`
  - 第 3.3 節：從「TP-4 符號相反是最嚴重異常（待驗證）」改為「公式已確認，兩個指標量不同面向（已驗證）」
  - 第 5.1 節：改進建議從「修正計算公式」改為「新增 AlleleMethDelta 欄位與現有 delta 並列」
  - 第 6 節優先級：第二優先從「修正 TP-4 符號錯誤」改為「新增 AlleleMethDelta 欄位」
  - 第 7 節 K3：更新知識點，說明 feature_allele_delta 是距離 delta，引用原始碼
  - 第 8 節行動清單：🔴 查驗行動改為 ✅ 已完成

- `docs/references/manual/20260317_甲基位點觀察報告生成技能_01.md`
  - 附錄特徵欄位說明：`AlleleDelta` 定義更新為 PERMANOVA 距離 delta，加入污染警示

- `docs/references/manual/20260317_case_obs_report_FP_B1_範例_01.md`
  - 數值摘要表：`AlleleDelta` 說明更新
  - Q2 解釋題 AI 初答：更新為正確的距離 delta 機制說明，並解釋 FP-B1 虛高原因

## 觀念更新

- 【舊】`feature_allele_delta = ALT_reads_mean_methylation − REF_reads_mean_methylation`（raw 甲基差）→ 【新】`feature_allele_delta = between_mean_distance(ALT,REF) − within_mean_distance(ALT+REF)`（PERMANOVA 距離 delta）
- 【舊】TP-4 的 TSV +0.023 vs 案例分析 −0.198 是「符號錯誤 bug」→ 【新】兩個指標量不同面向，不矛盾；正值代表距離空間分群良好，負值代表 ALT reads 的 raw 甲基低於 REF reads，可同時成立

## 後續行動

- [ ] 設計並實作 `AlleleMethDelta` 欄位（`ALT_reads_mean_meth − REF_reads_mean_meth`，per-read level）輸出到 feature matrix
- [ ] 在 InterSubMod 輸出中加入 `HP1MethMean` / `HP2MethMean` 欄位
- [ ] 設計 `HPDrivenFP_candidate` flag（HPP < 0.05 AND DominantLabel=hp AND AlleleHPCorrScore > 0.8）
- [ ] 對 FP-B2（chr9:75383880）補做完整觀察報告（v4.0 格式）
- [ ] 在更大批次驗證「HPP < 0.05 AND DominantLabel=hp」規則的 precision/recall

## 對話摘要

本次對話查驗 InterSubMod 原始碼（`src/core/LabelTest.cpp`），確認 `feature_allele_delta` 是 PERMANOVA 距離空間的 between-within delta，而非 raw 甲基差（ALT_meth−REF_meth）。TP-4 的「符號相反異常」解除為誤判——兩個指標量的是不同面向，可同時成立。FP-B1 的高 AlleleDelta（+0.058）虛高機制確認為 ALT-HP 完全共線導致 within_ALT 緊密、between(ALT,REF) 虛大。依此結果更新三份文件，移除錯誤假說，並將改進建議從「修正計算公式」調整為「新增 raw 甲基差 AlleleMethDelta 欄位以補充現有距離 delta」。
