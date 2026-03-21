<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# HCC1395_DORADO Methylation Filter 條件分解與門檻搜尋報告
- 日期：2026-02-13
- 資料：`HCC1395_DORADO` run-tag `20260213_dorado_purity_full`
- 輸入 purity：`t7_n29, t19_n29, t30_n20, t40_n10, t50_n00`（`t00_n25` 已排除）
- Truth：SEQC2（`TRUTH_TOTAL=39447`）

## 1. 問題與目標
你要驗證：
1. 造成 F1 下降的主因是 `QUAL < 0.75` 還是 `AlleleDelta/CramersV/VAF` 條件。
2. 是否能透過調整數值讓 F1 轉為穩定正向。

## 2. 先確認目前流程實際邏輯
`03_filter_analysis.py` 現況：
- `QUAL < 0.75` 條件 **已停用（註解）**。
- 目前預設條件為：`AD > 0.15 AND CV < 0.03 AND VAF < 0.15`。

另外，`01_longphase_s.sh` 內建了 `GQ Integer -> Float` header 修正，避免 `bcftools` 在 pileup VCF 因 `GQ=.` 失敗。

## 3. 分解測試設計
對每個 purity 做四種條件比較（皆以 LongPhase baseline 為起點）：
1. `QUAL_only_q0.75`：只用 `QUAL < 0.75`
2. `ACV_only_a0.25_c0.05_v0.24`：只用 `(AD>0.25 & CV<0.05 & VAF<0.24)`
3. `QUAL_OR_ACV`：`(QUAL<0.75) OR ACV`
4. `CURRENT_a0.15_c0.03_v0.15`：目前程式預設

## 4. 分解結果（重點）
### 4.1 平均 F1 變化（5 組 purity）
- `QUAL_only_q0.75`：**avg ΔF1 = -0.0717**（最差）
- `QUAL_OR_ACV`：**avg ΔF1 = -0.0804**（更差）
- `ACV_only_a0.25_c0.05_v0.24`：**avg ΔF1 = -0.0121**
- `CURRENT_a0.15_c0.03_v0.15`：**avg ΔF1 = -0.0138**

### 4.2 關鍵 purity（低 purity t7_n29）
- Baseline F1：`0.4673`
- `QUAL_only_q0.75`：`0.2918`（Δ `-0.1755`，TP 移除 5299，FP 移除 25）
- `ACV_only_0.25/0.05/0.24`：`0.4200`（Δ `-0.0473`，TP 移除 1546，FP 移除 2）
- `QUAL_OR_ACV`：`0.2602`（Δ `-0.2071`，最嚴重）

結論：**主要傷害來源是 QUAL 條件**；`QUAL<0.75` 在此批 purity 會大量誤刪 TP，Recall 明顯下跌。

## 5. 門檻搜尋結果

### 5.1 OR 組合搜尋
搜尋空間：
- `q in [0.60,0.65,0.70,0.75,0.80]`
- `a in [0.15,0.20,0.25,0.30,0.35]`
- `c in [0.02,0.03,0.05,0.07]`
- `v in [0.10,0.15,0.20,0.24,0.30]`

最佳（依 `all_non_negative` 優先 + 平均 ΔF1）：
- `q=0.60, a=0.35, c=0.02, v=0.10`
- 但結果仍：`avg ΔF1 = -0.02017`，`min ΔF1 = -0.056736`

=> **OR 組合在這批資料找不到正向方案**。

### 5.2 ACV-only（不含 QUAL）搜尋
搜尋空間（擴大）：
- `a in [0.15..0.45]`
- `c in [0.005..0.07]`
- `v in [0.05..0.30]`

最佳全域：
- `a=0.15, c=0.005, v=0.05`
- `avg ΔF1 ≈ +0.000019`、`min ΔF1 ≈ -0.000021`
- 實際上幾乎等於「不過濾」（多數 purity `tp_removed=0, fp_removed=0`）

=> **要維持不傷 F1，只能把條件調到幾乎不動作**；有意義過濾時仍會整體負向。

## 6. 為何會這樣（資料特性）
1. LongPhase 後 FP 已低（尤其低 purity 的 FP 絕對數很小），可移除的 FP 空間有限。
2. QUAL 門檻在此資料分佈上對 TP 不友善（大量 TP 落在低 QUAL 區間）。
3. ACV 條件雖比 QUAL 好，但在低 purity 仍常是 TP 移除多於 FP 移除。

## 7. 可執行結論
1. 在 `HCC1395_DORADO` 這批 purity 上，`QUAL < 0.75` **不建議啟用**。
2. 固定單一門檻（含你指定的 `0.25/0.05/0.24`）**無法穩定提升全部 purity 的 F1**。
3. 若目標是「全 purity 不退步」，目前可行解接近「不過濾」或僅在高 purity 進行極弱過濾。
4. 下一步應改成 **purity-aware** 策略（低 purity 關閉或極保守，高 purity 才嘗試過濾）。

### 7.1 每個 purity 單獨最佳（ACV-only）補充
來自 `acv_only_best_global.json`：
- `t40_n10`：`a=0.25,c=0.005,v=0.15`，ΔF1 `+0.000057`
- `t50_n00`：`a=0.15,c=0.005,v=0.10`，ΔF1 `+0.000069`
- 其他 purity 最佳通常是「不過濾」（TP/FP 都不變）

上述提升量皆遠小於 0.001，實務上接近數值抖動等級，不能視為穩定可用提升。

## 8. 產出檔案
- 分解結果：
  - `.../component_analysis/components_eval.tsv`
- OR 搜尋前 20：
  - `.../component_analysis/grid_eval_top20.tsv`
- OR 最佳：
  - `.../component_analysis/best_global.json`
- ACV-only 搜尋前 20：
  - `.../component_analysis/acv_only_grid_top20.tsv`
- ACV-only 最佳（含每個 purity 單獨最佳）：
  - `.../component_analysis/acv_only_best_global.json`

## 9. 相關程式
- 分解/搜尋腳本（新增）：
  - `scripts/analysis/evaluate_dorado_purity_filter_components.py`
- GQ 修正（已在 pipeline）：
  - `scripts/pipeline/steps/01_longphase_s.sh`
