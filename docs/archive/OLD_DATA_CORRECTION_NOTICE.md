# 舊數據矯正公告（P3-D）

> **建立日期**: 2026-04-19
>
> **⚠️ 本公告適用於本專案 2026-04-04 前所有報告的 TP/FP 數字**

---

## 矯正摘要

| 項目 | 舊數值（PRE-FIX） | 新數值（2026-04-04 後 canonical） | 矯正原因 |
|------|------------------|--------------------------------|---------|
| **TP** | 30,490 | **28,509**（paired） / **28,383**（TO） | VCF 來源錯誤矯正（2026-04-04）：原 pileup symlink 實為 ClairS-paired，非 TO |
| **FP** | 4,842 | **11,606**（paired） / **11,830**（TO） | 同上 |
| **F1** | ~0.89 | **0.7153 V3**（paired） / **0.7127 TO V5**（ClairS-TO） | 同上 + V3/V5 Haplotag 修正 |

**關鍵事件時間線**：

| 時間 | 事件 |
|------|------|
| 2026-03 及之前 | TP=30,490/FP=4,842 — 基於誤判的 pileup VCF |
| 2026-03-28 | Haplotag V3-Fixed：HP tag integer cast bug 修正 |
| 2026-04-04 | **VCF 來源錯誤矯正** — 重新跑 ClairS-paired/TO，得新 TP/FP |
| 2026-04-11 | Haplotag V5 Somatic Fallback：AMB% 17.5→8.0%，F1=0.7154 |

---

## 受影響文件列表（P3-D grep 盤點）

以下文件引用舊數值 30,490 / 4,842，**僅作歷史記錄，不可引用於新結論**：

### 高影響（核心報告）
- `docs/research/methylation_f1_optimization/README.md`
- `docs/research/methylation_f1_optimization/2026/01/20260115_*.md`（3 files）
- `docs/solutions/optimization/2026/01/20260107_F1_score影響評估_01.md`
- `docs/provenance/ai_sessions/2026/02/20260209_InterSubMod專案完整分析報告_01.md`
- `docs/plans/2026/01/20260113_甲基化F1研究計劃_01.md`
- `docs/plans/2026/04/20260403_LongPhase_TO_二次Phasing方案設計與驗證計畫_01.md`

### 中影響（細部分析）
- `docs/experiments/validated/2026/01/20260107_F1_Optimization_Deep_Analysis_01.md`
- `docs/experiments/validated/2026/01/20260107_F1_and_Data_Optimization_Report_01.md`
- `docs/experiments/validated/2026/01/20260107_report_round_06_F1_optimization_01.md`
- `docs/experiments/validated/2026/01/20260114_輸入與計算量效能分析_01.md`

### 低影響（架構/早期實作）
- `docs/concepts/2025/12/20251202_*.md`（2 files）
- `docs/archive/deep/2025-12_old_structure/reports/*.md`（2 files）

### 豁免（本文件本身）
- `docs/archive/OLD_DATA_CORRECTION_NOTICE.md`（本文件）
- `docs/reports/audit/decisions/04_P3_documentation_sync.md`（P3-D 決策定義）
- `docs/reports/audit/decisions/CHECKLIST.md`（引用本公告）

---

## 引用規範

### ❌ 不可寫

```markdown
本專案 TP=30,490, FP=4,842（HCC1395）
```

### ✅ 可寫（新結論）

```markdown
本專案 TP=28,509 / FP=11,606 (HCC1395 paired, ClairS-paired caller + V3 Haplotag, 2026-04-04 矯正後)
```

### ✅ 可寫（歷史追溯）

```markdown
2026-04-04 VCF 來源矯正前 TP=30,490（舊，見 `docs/archive/OLD_DATA_CORRECTION_NOTICE.md`）→ 矯正後 TP=28,509
```

---

## CHECKLIST P3-D 決策記錄

**採用選項**: 混合 A+B — 因 14+ 個受影響檔案分布廣、逐處加註成本高，採用「單一矯正公告 + 受影響檔案清單」取代逐處註記。未來引用時讀者可循 `docs/archive/OLD_DATA_CORRECTION_NOTICE.md` 了解矯正背景。

---

## 關聯文件

- `docs/reports/audit/decisions/04_P3_documentation_sync.md#P3-D`
- `docs/CURRENT_FOCUS.md`（當前 canonical 數據）
- `output/bip8_output_archive/README.md`（PRE-FIX 歷史數據警告）
