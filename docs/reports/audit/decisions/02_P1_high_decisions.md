# 02 — P1 High Decisions（3 題）

> **建立日期**: 2026-04-19

> **優先等級**：P1（數值錯誤 / 統計瑕疵）
>
> **依賴**：P1-A 可獨立立即修；P1-B/P1-C 建議 P0 完成後執行

---

## 問題 P1-A: r=0.997 誤標修正範圍

**問題描述**：
`docs/reports/research_landscape/08_Zone_Aware.md:49,73` 記載 Coverage_Multiple 相關係數為 **r=0.997**，但實際數值應為 **r=0.831**（CovM vs SEQC2 CN proxy），或 **r=0.962 / Jaccard=1.0000**（LOH.bed SEQC2 concordance）。此誤標跨引用至 4 個 audit cards（C13, C17, C20, C22）。

**影響範圍**：
- 受影響結論：**C13**（LOH.bed SEQC2 Jaccard=0.847）, **C17**（LOH Subclone 引用）, **C20**（CovM vs CN 核心）, **C22**（Zone-Aware framework）
- 受影響功能：F3 Zone-Aware, F5 Variant Confidence（在 reviewer 眼中是誠信瑕疵）
- 若不處理的風險：論文引用時被 reviewer 抓到數值不符；內部文件一致性喪失

**現況證據**：
- Audit card: `cards/C13_LOH_Bed_SEQC2.md`, `cards/C17_LOH_Subclone_AF.md`, `cards/C20_CovM_Non_Independent.md`, `cards/C22_Zone_Aware.md`
- 錯標位置：`docs/reports/research_landscape/08_Zone_Aware.md:49` 與 `:73`
- 實際數值來源：
  - CovM vs SEQC2 CN Pearson r=0.831（Paired）, r=0.827（TO）— 見 `scripts/analysis/methylation_cn_validation.py`
  - LOH.bed SEQC2 Jaccard=0.928（HCC1395）— 見 `research/seqc2_cnv_stratification/`

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 同步修 `08_Zone_Aware.md:49,73` 主源 + 關聯 4 卡的引用（加註「修正自 r=0.997」）+ grep 全 docs 確認無殘留 | 30 分鐘 | 低 — 純文字編輯 |
| B | 僅修主要文件（`08_Zone_Aware.md`），audit card 引用保留（在 P1-A 狀態標示下視為舊 draft） | 10 分鐘 | 中 — audit card 仍錯 |
| C | 打補丁加註「原 draft 錯標」，保留原值供歷史追溯 | 5 分鐘 | 高 — 數值錯誤留在文件中，論文引用會出錯 |

**驗證標準（無論選 A/B/C）**：
- **必達**：`grep -r "r=0.997" docs/` 回傳 0 筆（或全部加註「舊錯值」）
- **必達**：`08_Zone_Aware.md:49,73` 顯示正確數值並附註來源腳本
- **驗收指令**：
  ```bash
  grep -rn "0.997" docs/ | grep -v "舊錯值\|原 draft"
  # 預期：無輸出
  ```

**依賴關係**：
- **前置依賴**：無（純文字修正）
- **被阻塞項**：無（但其他文件一致性依賴此修正）

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P1-B: FDR 校正範圍

**問題描述**：
60+ 特徵 AUC 掃描（含 NG × NR × AF 組合）無 BH-FDR 校正。「最佳組合」聲明可能是多重測試的 artifact。C16 HPFineNGroups 的 NG=4+AF<0.4+NR≥80 filter 是最大風險點（90 組合 screening）。

**影響範圍**：
- 受影響結論：**C03**（TO AUC ceiling）, **C04/C05/C06**（O11/O12/O13 NEGATIVE pre-stratify AUC）, **C07**（G1-G7 7 特徵）, **C11**（Phase 1A F1 3 模型）, **C16**（HPFineNGroups 90 組合）
- 受影響功能：F2 Subclone Marker（C16）, F5 Variant Confidence（C11）
- 若不處理的風險：C16 的「最佳 filter」聲明在 reviewer 下不成立

**現況證據**：
- Audit card: `cards/C16_HPFineNGroups.md`（無 FDR）, 其他受影響卡
- 相關 cross_cutting: `cross_cutting/Multiple_Testing_Correction.md`（含 FDR 方案）

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **C16 優先**（90 組合 BH-FDR）→ **C22 Zone 補算**（需 P0-A 完成後）→ **O11/O12/O13 追溯補算**（強化 NEGATIVE 敘事）→ C07 Bonferroni（n=7） | 1 週 | 低 — 每步可獨立驗收 |
| B | 全 7 卡同時 BH-FDR（一次性補完） | 5 天 | 中 — 平行但 C22 需等 P0-A |
| C | 僅 POSITIVE 結論（C11/C16）做 FDR；NEGATIVE 不追溯 | 2 天 | 中 — 錯失 NEGATIVE 穩固機會 |

**驗證標準（無論選 A/B/C）**：
- **必達**：C16 90 組合 BH-FDR adjusted p 表；NG=4+AF<0.4+NR≥80 survives（adjusted p<0.05）
- **觀察**：若該組合不 survive → C16 結論降為 CONDITIONAL，報告所有 surviving combos
- **觀察**：O11/O12/O13 pre-stratify AUC 多數仍 FDR-adjusted 顯著（強化「統計顯著 ≠ 真 signal」敘事）
- **驗收指令**：`scripts/analysis/fdr_correction.py`（新建或擴展）

**依賴關係**：
- **前置依賴**：部分需 P0-A 完成（C22 Zone 補算）
- **被阻塞項**：論文 Results 章節數據定稿

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P1-C: HPFineNGroups bootstrap CI

**問題描述**：
C16 HPFineNGroups NG=4+AF<0.4+NR≥80 的 residualized AUC 為單點聲明，無 CI。跨樣本一致性（7/7 或 5/7 ≥0.85）也無 CI。

**影響範圍**：
- 受影響結論：**C16**（HPFineNGroups 核心）
- 受影響功能：F2 Subclone Marker
- 若不處理的風險：「7/7 跨樣本一致」聲明在 reviewer 下不精確

**現況證據**：
- Audit card: `cards/C16_HPFineNGroups.md` D5 「CI 覆蓋率: ❌」
- 相關 cross_cutting: `cross_cutting/Multiple_Testing_Correction.md`（含 bootstrap 方案）

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 1000× bootstrap CI + per-AF-bin CI + **結合 P0-B within-group OLS**（一次完成） | 2-3 天 | 低 — 與 P0-B 合併執行最省時 |
| B | 僅整體 CI（不做 per-bin） | 1 天 | 中 — 無法揭露 per-bin 穩定性 |
| C | 省略（C16 跨樣本 7/7 穩定度已觀察為「壓倒性」） | 0 天 | 中-高 — 論文 reviewer 會要求 |

**驗證標準（無論選 A/B/C）**：
- **必達**：1000× bootstrap 樣本；每次重抽 regions（樣本內層級）
- **必達**：整體 AUC 95% CI；per-sample AUC 95% CI；per-AF-bin AUC 95% CI
- **觀察**：CI lower bound >0.58（POSITIVE 門檻）→ 結論穩固
- **驗收指令**：`scripts/analysis/bootstrap_ci_auc.py` 輸出 TSV + forest plot

**依賴關係**：
- **前置依賴**：P0-B within-group OLS（同時執行最佳）
- **被阻塞項**：論文 F2 圖表

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## P1 區段總結

**推薦總動作**（全選 A）：
1. **立即執行 P1-A**（30 分鐘）：修正 `08_Zone_Aware.md:49,73` + 全檔 grep
2. **P0-A 修正完成後執行 P1-B**：C16 BH-FDR 優先，其餘追溯
3. **P0-B 執行時同時做 P1-C**：bootstrap CI + within-group OLS 合併

**推薦理由**：
- P1-A 獨立可立即修（最易拿的分數）
- P1-B/P1-C 與 P0 耦合，合併執行最省時
- 三項補齊後論文聲明具備 reviewer-defensible 嚴謹度

**若全選 A 後的輸出**：
- `08_Zone_Aware.md:49,73` 修正後版本
- `output/analysis/c16_fdr_adjusted.tsv`（90 組合 BH-FDR）
- `output/analysis/c16_bootstrap_ci.tsv` + forest plot
- `cards/C16.md` D5 升級 ✅ FDR + bootstrap 完成

---

## 關聯文件

- EXECUTIVE_DECISION_BRIEF.md 第 4-5 節
- `cross_cutting/Multiple_Testing_Correction.md`
- `cross_cutting/Pooled_OLS_Audit.md`（P1-C 合併執行依據）
- CHECKLIST.md
