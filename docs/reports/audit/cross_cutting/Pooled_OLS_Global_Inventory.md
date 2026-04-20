# Cross-cutting: Pooled OLS Global Inventory（R-05 全面盤點）

> **建立日期**: 2026-04-19
>
> **範圍**：全專案 `scripts/analysis/` + `output/synthesis/` 所有使用 `residualize` / `residualise` / `residual` 的程式與報告清單
>
> **目的**：識別哪些結論可能受 **pooled OLS residualization trap**（見 `feedback_pooled_ols_residualization_trap.md`）影響；分類為「already within-group」/「pooled-only」兩類，供 R-05 決策與 P0-B 重算範圍界定。
>
> **來源決策**：CHECKLIST R-05 選 A（grep 盤點，不重跑）+ P0-B 選 B（C15/C16/C17 三者並行 within-group 重算）

---

## 1. 掃描方法

```bash
# Script 層級
grep -rn "residuali" scripts/analysis/*.py

# Report 層級
grep -rn "residuali" output/synthesis/

# 分類標記
grep -qi "within.group\|per[-_ ]sample\|groupby.*sample" <script>
# → HAS within-group/per-sample pattern → 標記「混合或 safe」
# → 否則 → 標記「POOLED」（需 P0-B 重算或 R-05 進一步檢驗）
```

---

## 2. Script 層級結果（8 個檔案）

| # | 檔案 | 類型 | 受影響結論 | 備註 |
|---|------|------|-----------|------|
| 1 | `verify_loh_cn_af_conclusions.py` | 🟢 **Mixed (per-sample + pooled)** | C04/C12/C15/C17 | Figure-level 驗證腳本；per-sample breakdown 存在 |
| 2 | `methylation_cn_validation.py` | 🔴 **POOLED-only** | C16, C20 | HPFineNGroups raw r=0.495 → residualized r=0.160（68% NumReads confound）為 pooled；P0-B 需重算 within-group |
| 3 | `run_m5_m6_m7_only.py` | 🔴 **POOLED-downstream** | M5/M6/M7 結論 | 消費 `m2_residualization_results.tsv`（pooled 來源）；隨 M2 重算後更新 |
| 4 | `build_beyond_auc_validation.py` | 🟡 **Mixed**（M2 pooled + per-sample stratum） | Beyond-AUC 耗盡 | M2 residualization 為 pooled；但 strata 含 per-sample；pooled 部分需 P0-B 重算 |
| 5 | `build_allele_deep_dive.py` | 🟢 **Mixed** | AlleleDelta 深度分析 | per-sample 驗證存在 |
| 6 | `build_allele_loh_auc_f1_analysis.py` | 🟢 **Within-group (per-sample partial corr)** | C17 LOH Subclone | per-sample partial correlation（AlleleDelta × NumReads × CovM）；R-04 無疑慮 |
| 7 | `build_observation_O13v2_cross_region_confound_control.py` | 🟢 **Within-group (per-sample median delta)** | C06 O13 | per-sample residualized median TP-TP vs FP-FP delta；典範方法 |
| 8 | `build_observation_O12_loh_methylation_scenarios.py` | 🟢 **Mixed** | C05 O12 | L1 residualize on reads+cpgs 為 pooled，L3 AF-bin 交叉驗證；已在 C05 R-04 補註說明 |

### 2.1 分類統計

| 分類 | 腳本數 | 風險等級 |
|------|--------|---------|
| 🔴 POOLED-only / POOLED-downstream | 2（`methylation_cn_validation.py`, `run_m5_m6_m7_only.py`） | **P0-B 需重算** |
| 🟡 Mixed（含 pooled 但有 strata） | 1（`build_beyond_auc_validation.py`） | pooled 段落需 P0-B 檢視 |
| 🟢 Mixed / Within-group | 5 | 已安全 |

### 2.2 受 P0-B 重算影響的已知結論

| 結論 | 相關 | 重算優先序 |
|------|------|-----------|
| **C16 HPFineNGroups residualized AUC=0.617** | `methylation_cn_validation.py` | **優先（P0-B 選項 B 並行 C15/C16/C17）** |
| **C17 LOH Subclone AF×Methylation** | `build_allele_loh_auc_f1_analysis.py`（per-sample 已 safe） | 已 within-group（確認現有結果穩固） |
| **C15 LOH Methylation Failure** | 對應 O15 `build_observation_O15*.py`（未列出） | 待補 grep（R-05 尾款） |
| **Beyond-AUC exhaustion** | `build_beyond_auc_validation.py` M2 pooled 段落 | 中優先（結論已 NEGATIVE 穩固，重算為補強） |

---

## 3. Report 層級結果（4 個檔案）

| # | 檔案 | 類型 | 內容 |
|---|------|------|------|
| 1 | `output/synthesis/observation_workspaces/20260401_O13v2_cross_region_confound_control/logic_review_report.md` | 🟢 基於 within-group 腳本 | O13 v2 邏輯審查 |
| 2 | 同目錄 `final_verdict_report.md` | 🟢 | O13 final verdict |
| 3 | 同目錄 `evidence_chain_report.md` | 🟢 | O13 證據鏈 |
| 4 | `output/synthesis/observation_workspaces/20260331_O12_loh_methylation_scenarios/observation_report.md` | 🟡 | O12 原始報告（含 L1 pooled + L3 AF-bin） |

4 個 reports 均為已完成觀察報告，無需重算，只需交叉連結本 inventory。

---

## 4. 未來新腳本的 residualize 規範（P0-B 後適用）

```python
# 禁止：pooled OLS residualize on group-linked confound
# residualized = y - OLS(y ~ n_reads).predict()  # ❌ 若 n_reads 與 sample/group 相關 → collider bias

# 建議：within-group OLS residualize
# for g in groupby(sample_id): residualized[g] = y[g] - OLS(y[g] ~ n_reads[g]).predict()
```

所有新增 residualize 腳本必須：
1. 明示 `# WITHIN-GROUP OLS` 或 `# POOLED OLS (acknowledged, with L3 AF-bin cross-check)` 註釋
2. 輸出同時包含 raw、pooled-residualized、per-sample-residualized 三版結果
3. 若 within-group 與 pooled 結果差異 >0.05 AUC，於 L3 AF-bin 交叉驗證

---

## 5. 與 audit cards 的交叉引用

| Card | 狀態 | 本 inventory 對應 |
|------|------|------------------|
| C04 O11 | ✅ 方法論典範（三方法交叉） | 非 residualize-only；R-04 已補註 |
| C05 O12 | ⚠️ L1 pooled + L3 AF-bin（安全） | 2.2 內已說明 |
| C06 O13 | ✅ 典範 within-group（per-sample median delta） | 🟢 段 #7 |
| **C15 LOH Methylation Failure** | 需 P0-B 重算（within-group） | 🔴 優先 |
| **C16 HPFineNGroups** | 需 P0-B 重算（within-group） | 🔴 優先 |
| C17 LOH Subclone AF | ✅ per-sample partial（safe） | 🟢 段 #6 |

---

## 6. 後續 TODO

- [ ] **P0-B 三者並行重算**：C15/C16/C17 within-group OLS AUC 與 pooled 差異 → 輸出 `p0b_within_group_vs_pooled_summary.tsv`
- [ ] **R-05 尾款**：對 C15 對應的 `build_observation_O15*.py` 腳本路徑 grep 補齊
- [ ] **建立 lint 規則**（低優先）：新 `.py` 檔 residualize 必須註解 WITHIN-GROUP 或 POOLED

---

## 關聯文件

- `docs/reports/audit/decisions/05_strategic_R_decisions.md#R-05`
- `docs/reports/audit/decisions/01_P0_critical_decisions.md#P0-B`
- `docs/reports/audit/cross_cutting/Pooled_OLS_Audit.md`
- `memory/feedback_pooled_ols_residualization_trap.md`
