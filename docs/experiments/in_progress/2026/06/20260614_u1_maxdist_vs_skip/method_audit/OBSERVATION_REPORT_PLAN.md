---
title: "觀察報告計劃 + 給另一 session 的生成入口（要產哪些觀察報告、傳哪份說明）"
date: 2026-06-14
type: plan / handoff-index
purpose: ①規劃本 session 成果要整理成哪些「給人看的觀察報告」 ②明確指出哪份說明傳給另一 session 生成
binary_baseline: develop 0f9320f（SKIP 預設 + HP-AUC 欄）; 方法框架在 feat branch e472760
---

# 觀察報告計劃

## A. 規劃要生成的觀察報告（給人看，非技術盤點）

| # | 觀察報告 | 目標讀者 | 核心內容 | 數據來源 | 入口材料（傳給另一 session）|
|---|---|---|---|---|---|
| **R1** | **ISM 甲基-標籤驗證方法現況解釋報告** | 用戶/PI（理解 ISM 在做什麼、哪些有效、哪些是 gap）| 兩方向 A/B → 6 目標 → 各方法現況/有效性 → 6 gap → Δβ 模組改進方向 | METHODS_BY_GOAL + SIGNIFICANCE_INVENTORY + TWO_DIRECTIONS | **METHODS_BY_GOAL.md**（主）+ TWO_DIRECTIONS + SIGNIFICANCE_INVENTORY |
| **R2** | **tumor HP-AUC < normal 檢核觀察報告** | 研究（normal>tumor 是否 somatic）| HP-AUC 分布 + normal>tumor 配對 delta + **confound 分層**(CNV/LOH) + 案例 + 誠實判定 | wg_hpauc.json + tumor_landscape.json + significance_summary | **HANDOFF_tumor_vs_normal_hpauc.md**（已備，最完整）|
| **R3** | **SKIP 去污染 + germline-HP 主導觀察報告** | 研究/PI | U1 MAX_DIST vs SKIP + 案例熱圖 + gap2 germline 主導 + HP-AUC 地基 | U1 README + gap2 + figures/ + wg_hpauc.json | U1 `../README.md` + `gap2/README.md` + figures/ |

## B. 🎯 明確：哪份傳給另一 session（你問的核心）

**取決於你要哪種觀察報告**：

### ▶ 若要「ISM 方法現況解釋報告」（R1，最可能你要的「明確觀察報告與解釋說明」）
**傳這份為主入口**：`InterSubMod/docs/experiments/in_progress/2026/06/20260614_u1_maxdist_vs_skip/method_audit/METHODS_BY_GOAL.md`
- 附帶（讓另一 session 有完整脈絡）：`TWO_DIRECTIONS_FRAMEWORK.md` + `SIGNIFICANCE_INVENTORY.md`。
- 另一 session 任務：把「目標→方法→敘述→解釋→結果→有效→替代」轉成**給人看的解釋報告**（白話 + 為什麼這樣設計 + 哪些有效/待補），可加 HTML/圖。
- 數據引用：本 session 的 wg_hpauc.json / gap2 / U1（每數字標來源）。

### ▶ 若要「tumor HP-AUC<normal 檢核報告」（R2）
**傳這份**：`HANDOFF_tumor_vs_normal_hpauc.md`（已是完整 handoff：背景+3 競爭假說+confound 分層+案例選取+誠實鐵則+交付物）。**這份最現成、可直接接手。**

### ▶ 若要「SKIP/germline-HP 觀察報告」（R3）
傳：U1 `README.md` + `gap2/README.md` + `figures/`（案例熱圖已在）。

## C. 給另一 session 的通用生成指南（所有觀察報告共用）
1. **分層揭露**：L0 一句結論 → 表/圖 → 細節 → 誠實 caveat（[[feedback_reports_need_layered_disclosure]]）。
2. **每數字標來源**：grep 得到才寫（significance_summary 欄 / *.json）；[[feedback_no_fabricated_numbers_in_reports]]。
3. **claim 分級**：L1(源碼+數據) / L2(推論)；germline-HP 地基=L1，normal>tumor somatic=L2(待 confound)。
4. **誠實鐵則**：HP=germline 真值（非 subclone）；單樣本 HCC1395 partial；方法 gap 如實標。
5. **格式**：用戶看懂用 HTML、決策用 AskUserQuestion、AI 紀錄用 .md（[[feedback_review_format_html_ask_md]]）。

## D. 材料索引（另一 session 的完整資源清單）
- **方法框架**（feat e472760）：method_audit/ 的 5 份 .md（本計劃 §A 列）。
- **數據 json**（可重現，env U1_SKIP）：wg_hpauc.json / comparison_summary.json / interpretation.json / gap2_subclone_discrimination.json / tumor_landscape.json / hp_auc.json。
- **可重現 script**：wg_hpauc_analyze.py / compare_u1.py / hp_auc.py / plot_cases.py / gap2_analyze.py / tumor_landscape.py。
- **圖**：figures/case_*.png（4 案例距離熱圖+甲基圖）。
- **重生 summary**（/tmp 會清）：develop binary + filtered_snv_tp.vcf.gz + `-w5000 BERNOULLI`（SKIP 已預設）。

## E. 本 session 待實作（你消化框架對齊後）
- **Δβ 模組**（最高槓桿，cpp-change）：8 mean + 配對矩陣 + permutation 檢定 → 統一 #1+#3+subclone+對齊外部。
- 詳見 COMPUTE_FLOW_ARCHITECTURE §3.5 + METHODS_BY_GOAL 2a/3a。
