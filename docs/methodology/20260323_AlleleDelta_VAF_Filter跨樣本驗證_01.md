<!--
建立時間: 2026-03-23
目標: 確認 AlleleDelta standalone filter 是否適合跨樣本部署
處理範圍: 7個樣本 AUROC 分析 + 決策記錄
關聯檔案:
  - scripts/pipeline/steps/03_filter_analysis.py (61-113行)
  - scripts/analysis/verify_class_decision_tree_audit.py
-->

# AlleleDelta Standalone Filter 跨樣本驗證

## 背景

前次分析（HCC1395 pileup 全量）發現：
- AlleleDelta > 0.25 作為 standalone filter（不需 VCF）
- 效果：FP removed=550 (11.4%), TP removed=214 (0.7%), ΔF1=+0.00264
- 原因：germline ASM FP 具有更高的 AlleleDelta

本文記錄跨樣本 AUROC 驗證結果與決策。

## 跨樣本 AUROC（FP-discriminating 方向）

| 樣本 | AlleleDelta AUROC | 結論 |
|------|------------------|------|
| HCC1395 5kHz paired_full | **0.763** | ✓ 有效 |
| HCC1395 DORADO paired_full | **0.412** | ✗ 無效（同樣本不同化學） |
| COLO829 paired_full | 0.545 | △ 邊際 |
| H1437 paired_full (FP=8) | 0.503 | — 樣本數不足 |
| H2009 paired_full | **0.379** | ✗ **TP AlleleDelta 反而更高** |
| HCC1937 paired_full | 0.212 | ✗ 反效果 |
| HCC1954 paired_full (FP=29) | 0.465 | — 樣本數不足 |

> AUROC > 0.5 表示 FP 的 AlleleDelta 更高（可用於過濾 FP）
> AUROC < 0.5 表示 TP 的 AlleleDelta 更高（過濾會誤傷 TP）

## 決策：**Rejected — 降為 annotation-only**

理由：
1. HCC1395 5kHz 有效，但 HCC1395 DORADO（同樣本、不同測序化學）無效
2. H2009 和 HCC1937 中，TP 的 AlleleDelta 反而更高，強制過濾會降低 F1
3. AlleleDelta 高的 FP 是 germline ASM，在 HCC1395 5kHz 特別明顯（可能是 5kHz 特有的偽影模式）

## 程式碼修改（v1 → v2）

**檔案**：`scripts/pipeline/steps/03_filter_analysis.py`

| 修改前 | 修改後 |
|--------|--------|
| `allele_delta_only_min: 0.25`（所有 purity bins） | `allele_delta_only_min: None`（所有 purity bins） |
| standalone filter 預設開啟 | standalone filter 停用 |
| 添加 `experimental_hcc1395_5khz_only: True` 標記 | — |

**影響**：
- HCC1395 5kHz：F1 從 +0.00264 回退至 baseline（可接受）
- COLO829/H2009/HCC1937：避免潛在 TP 誤傷
- VCF-based AD+VAF 主要規則不受影響（仍保留，需 VCF 輸入）

## AlleleDelta 特徵保留

AlleleDelta 欄位仍輸出至 significance_summary.csv，供未來分析使用：
- 跨樣本 AlleleDelta 分佈分析
- 配合其他特徵的組合 filter 設計

## 狀態

**Decision: Rejected (annotation-only)**

下一步：
- [ ] 確認 COLO829/H2009 實際 F1 無回退（預期：維持不變）
- [ ] 評估 LOH+Noise 組合 filter（H2009 Noise-FP LOH率=81.8%，跨樣本更一致）
