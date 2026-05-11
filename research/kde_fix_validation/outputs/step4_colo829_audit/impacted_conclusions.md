---
created: 2026-04-20
status: audit
scope: COLO829 conclusions impacted by stale baseline 75× (→ new KDE 29×)
---

# COLO829 既有結論審計：KDE baseline 修正後的衝擊

## 背景

- **Stale baseline**：`expected_coverage = 75.0`（全 7 樣本共用 default）
- **COLO829 新 KDE baseline**：29×（ratio 75/29 = **2.59×**）
- 這是**全 7 樣本最極端的 drift**，所有以 `Coverage_Multiple`、`Coverage_Category`、或受 CNV_Loss penalty 影響的 QS 結論都需重新審視
- **COLO829 新 CovM p50**：1.000（從 stale 0.387）— 代表 COLO829 **實際上是 CN-neutral dominant**，不是「定序深度偏低造成的 CNV_Loss artifact」

## 分類：影響嚴重度

### 🔴 CRITICAL — 結論直接錯誤，需撤回或重寫

#### C1. O1-O10 Fig 11: COLO829 Low QS Investigation
**檔案**：`docs/reports/validated/2026/04/20260401_comprehensive_observation_O1_O10_report_01.md:888`

**原結論**（引用）：
> COLO829 QS median ~60（Others ~85-100）。原因：coverage multiple 0.387（Others ~0.9）導致 CNV_Loss penalty (-30) 全面觸發。Low QS 由 low depth 驅動，非甲基化異常。

**問題**：`CovM=0.387` 是 stale baseline 75× 除以實際 29× 的**人為結果**，不是真實 CN 狀態。新 CovM=1.0，CNV_Loss 不會被觸發。「Low QS 由 low depth 驅動」這個解釋**可能不成立**。

**建議**：重跑 QS 計算於新 CovM 下，若 COLO829 QS 仍偏低 → 真實甲基化/定序問題；若 QS 回正 → 本段須撤回。

---

#### C2. TO LOH 遮罩與反轉：COLO829 91.7% 被排除
**檔案**：`docs/reports/validated/2026/04/20260408_TO_LOH額外研究_遮罩與反轉分析_01.md:379-387`

**原結論**：
> COLO829 | **91.7%** | 33.6% → 35.3% | 定序深度 ~33x → CovMul 系統性偏低（NumReads median=29, CovMul=0.39）
> M2 嚴重偏向 COLO829（91.7% 被排除），因原始定序深度僅 ~33x，而 Coverage_Multiple 分母固定為 75x，使 CovMul 系統性偏低至 0.387。

**問題**：93% 排除率是 M2 filter 依 CovMul ≥ 0.5 實作，stale 75× 下 COLO829 的 CovMul=0.39 幾乎全數被濾掉。新 baseline 29× 下 CovMul=1.0，**M2 不會再排除 COLO829 絕大多數**。此段已主動註記為「已知 bug」，本 cycle 確認修正後需重跑 M2/M4。

**建議**：新 cycle 重跑 M2/M4 masks；更新「M2 偏向 COLO829」為「M2 在 stale baseline 下偏向 COLO829（已修正）」。

---

#### C3. ISM 對黑色素瘤無效的推論
**檔案**：`docs/reports/validated/2026/03/20260322_cross_sample_TO_ISM_gradient_analysis_01.md:170, 202, 219, 231`

**原結論**：
> COLO829 的 QS 整體偏低（中位數 35），且 TP/FP 完全無法區分——黑色素瘤的 ISM Quality Score 比其他癌症類型低得多。
> ISM 在黑色素瘤中 QS 特別低（35 vs 75），且效果最弱——**ISM 可能不適合黑色素瘤**。

**問題**：若 C1 成立（QS 偏低由 CNV_Loss penalty 造成），則「ISM 不適合黑色素瘤」是**baseline artifact**，非真實生物學差異。

**建議**：**暫時標記為 CONDITIONAL**。重跑 QS 於新 CovM 後：
- 若 COLO829 QS 回升至 ≥75 → 撤回「ISM 不適合黑色素瘤」
- 若 QS 仍 <60 → 真實效應，保留

---

#### C4. COLO829 低 QS 跨樣本觀察
**檔案**：`docs/reports/validated/2026/03/20260317_跨樣本甲基特徵TP_FP分離觀察報告_01.md:91-499`

**原結論**：
> COLO829 的 TP median Quality_Score=60.0，其他 sample ≥75；TP/FP 箱型幾乎完全重疊；可能反映 PAO 平台校準問題。

**問題**：同 C1/C3，低 QS 可能 50%+ 來自 CNV_Loss penalty artifact。PAO 平台假說未被真正隔離。

**建議**：重跑 QS 於新 CovM 後分離兩效應（CNV penalty 貢獻 vs 平台貢獻）。

---

### 🟠 HIGH — 結論方向對但量化需更新

#### H1. TO LOH 內外 ISM 區分力推論鏈
**檔案**：`docs/reports/validated/2026/04/20260408_TO_LOH內外ISM特徵區分力完整推論鏈報告_01.md:252`

**原結論**：
> COLO829：黑色素瘤 | ~100% (cell line) | 12.0% | 22.5% | 0.546 | 定序深度最低 ~33x（PAO yield 102Gb）；高 TMB；**CovMul 系統性偏低**

**問題**："CovMul 系統性偏低" 是 stale baseline artifact。定序深度 33× 是事實，但 CovM 在新 baseline 下 **不再偏低**（p50=1.0）。

**建議**：修正為「定序深度 ~33x，新 baseline 29x 後 CovMul 中位數為 1.0（與其他樣本相當），**原 0.39 是 baseline 75× 的人為偏誤**」。

---

#### H2. O1-O10 Fig 10: PCA 跨樣本分離
**檔案**：`docs/reports/validated/2026/04/20260401_comprehensive_observation_O1_O10_report_01.md:824`

**原結論**：
> PC1 (34.6%) + PC2 (27.5%) = 62.1%。COLO829 孤立（low depth）。HCC1937 極端（high depth）。Depth 是跨 sample 分離第一驅動因素。

**問題**：COLO829 「孤立」可能部分來自 CovM 被壓縮到 0.39 造成 Coverage_Category 集中於 CNV_Loss。新 CovM 下 COLO829 的 Category 分佈與其他樣本更接近。

**建議**：重跑 PCA 於新 CovM 下，確認 COLO829 是否仍孤立；若不再孤立 → 「Depth 是第一驅動」需降級為「baseline artifact 造成的 Category 偏移」。

---

#### H3. LOH enrichment post-HP-fix 報告
**檔案**：`docs/reports/validated/2026/04/20260401_LOH_weekly_review/03_post_hp_fix_loh_enrichment.md`
**檔案**：`docs/reports/validated/2026/03/20260330_TO_LOH_enrichment_post_hp_fix_01.md`

**影響**：內容涉及 COLO829 coverage 相關觀察；具體行數未深挖但涉及 M2/M4 enrichment 計算。

**建議**：本 cycle 不深入修改；標記為「依賴 stale baseline，下一 cycle 重跑」。

---

### 🟡 LOW — 計劃性論述，不影響結論

#### L1. TO pipeline 多階段 characterization（planning doc, trash）
**檔案**：`docs/trash/to_pipeline_staging_v1/reports/20260412_TO_pipeline_multi_stage_characterization_01.md:292`

**內容**：「將本分析延伸至 COLO829, H1437 等其他 6 樣本，驗證 CNV 判別力是否一致」— planning 語句，非結論。

**建議**：無需修改。

---

#### L2. Z3 amplicon blacklist pilot（相關 cycle）
**檔案**：`docs/experiments/in_progress/2026/04/20260419_Z3_amplicon_blacklist_pilot_result_01.md:51`

**內容**：「COLO829 S2=−0.0109：chr5/8/17 上有 Z3 TP 被誤殺」— 與 Task B 重疊，本次 cycle 中 Task B 將重算 ΔF1。

**建議**：等 Task B 完成後連動更新。

---

## 總結

| 嚴重度 | 數量 | 建議 |
|--------|------|------|
| 🔴 CRITICAL | 4 | 重跑 QS 於新 CovM；若結論翻轉，撤回/重寫 |
| 🟠 HIGH | 3 | 修正量化描述；標記「baseline artifact」 |
| 🟡 LOW | 2 | 無需動作或等 Task B |

**下一 cycle 必要 follow-up**：
1. 重跑 QS 計算於新 CovM 下（QS 依賴 CNV_Loss penalty）
2. 重跑 M2/M4 masks（依賴 CovMul ≥ 0.5 閾值）
3. 重跑 跨樣本 PCA
4. 以新 QS + mask 重跑「ISM 對黑色素瘤無效」假說

**本 cycle 完成**：
- 盤點完成 9 筆衝擊（4 CRITICAL、3 HIGH、2 LOW）
- 不在本 cycle 重跑 QS（需重新跑 C++ pipeline with updated logic），先記錄
- 在結論文件加註 baseline 修正導致的 validity decay

## 關鍵數據對照

| 指標 | Stale (75×) | Fixed (KDE 29×) | 變化 |
|------|-------------|-----------------|------|
| COLO829 CovM p50 | 0.387 | 1.000 | +0.613 |
| COLO829 CNV_Loss % | ~60% (估) | 3.93% | -56 pp |
| COLO829 Normal % | ~20% (估) | 42.75% | +22 pp |
| QS median (舊) | 35-60 | 重跑 | 預期提升 |
| M2 mask 排除率 | 91.7% | 重跑 | 預期 <30% |
