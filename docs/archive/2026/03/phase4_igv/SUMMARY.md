<!--
建立時間: 2026-04-05 22:00
目標: Phase 4 Subclone Case Study / IGV 觀察 / FP 分型研究的封存摘要
處理範圍: 2026-03-15 ~ 2026-03-17 的 8 個研究文件
關聯檔案:
  - docs/experiments/INDEX.md（條目 25-28 相關）
  - docs/reports/research_landscape/01_TO_FP問題全貌.md
-->

# Phase 4 Subclone / IGV / FP 分型研究封存摘要

> **封存原因**：本方向已被 O1-O13 系統性觀察 + G1-G7 特徵搜索完全取代。
> 檔案中引用的圖片（casestudy heatmaps、IGV snapshots、FP typing plots）從未被 git 追蹤。

## 關鍵結論

### 1. Phase 4 Subclone Case Study（20260315）
- 對 HCC1395 5kHz paired 的 5 個 Strong Subclone 候選位點（CramersV>0.3, NumReads≥30）做甲基化視覺觀察
- 觀察到 cluster heatmap + distance heatmap + MethylArtist 的一致模式
- **後續判定**：O1-O13 證明這些視覺模式無法系統性區分 TP/FP（AUC<0.58）

### 2. IGV 自動截圖與觀察（20260316）
- 建立 full-template 9 loci（5 TP + 4 FP）的標準化 IGV PNG 觀察流程
- 驗證 AI 視覺初篩可辨識 ALT 成群、HP/allele 驅動與相鄰異常
- **後續判定**：視覺化輔助有價值，但無法自動化為定量 filter

### 3. FP 三類型分類（20260317）
- 建立 FP 分型框架：B1 共定位（ISM 可能有效）/ B2 弱定位 / B3 隨機（ISM 無效）
- H2009 FP triage 實測，跨樣本篩選規則精確度/召回率量化
- **後續判定**：G1-G7 證明即使 60+ 特徵組合，FP removal 在 TP loss ≤2% 下仍 = 0%

### 4. 週報整合（20260317）
- 兩份週報整合 Phase 4 + IGV + FP 分型 + 跨樣本甲基特徵分析
- 當週主要進展：docs/技能制度化、subclone case study 方法建立

## 檔案清單

| 檔案 | 行數 | 圖片引用 | 說明 |
|------|------|---------|------|
| `20260315_subclone_phase4_casestudy_01.md` | 854 | 27（缺失） | 5 位點甲基化視覺觀察 |
| `20260316_IGV_full_template忠實版案例導讀_01.md` | 383 | 19（缺失） | 9 loci IGV 忠實版導讀 |
| `20260316_IGV固定template_AI視覺初篩與正式驗證_01.md` | 300 | 9（缺失） | AI 視覺初篩驗證 |
| `20260316_IGV截圖與數據對照觀察_01.md` | 220 | 9（缺失） | IGV + 數據對照 |
| `20260317_FP分型與多樣本篩選規則觀察報告_01.md` | 367 | 3（缺失） | FP 三類型分類 |
| `20260317_研究主線週報_..._phase4_FP分型_01.md` | 585 | 11（缺失） | 週報：Phase 4 + FP 分型 |
| `20260317_研究主線週報_..._subclone_igv_governance整合_01.md` | 369 | 4（缺失） | 週報：整合版 |
| `20260317_case_obs_report_FP_B1_範例_01.md` | 283 | 8（缺失） | FP B1 案例報告範例 |
