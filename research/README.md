<!--
建立時間: 2026-04-05 18:00
更新時間: 2026-04-14
目標: 研究工作區總索引
處理範圍: research/ 下所有研究主題目錄（10 個 + 1 範本）
關聯檔案:
  - docs/experiments/INDEX.md
  - docs/CURRENT_FOCUS.md
-->

# Research 工作區索引

本目錄存放 InterSubMod 各研究主題的工作區。每個子目錄為一個獨立研究，包含數據、腳本、圖表與報告。

## 研究目錄

| 目錄 | 研究主題 | 建立時間 | 狀態 | 結論 | 說明 | 對應 docs/ |
|------|---------|---------|------|------|------|-----------|
| `loh_investigation/` | LOH 深度調查與 Self-phasing 驗證 | 2026-04-02 | 進行中 | — | Round 1-4、PON-only、因果鏈驗證 | `reports/research_landscape/` |
| `loh_subclone_af/` | LOH Subclone AF × Methylation 雙重證據鏈 | 2026-04-14 | 已完成 | **POSITIVE** | H1-H4 全 supported，Inter AF→NGroups +0.705 | `experiments/INDEX.md` |
| `to_pipeline_staging/` | TO Pipeline 多階段特徵分析 (v2) | 2026-04-13 | 進行中 | — | ClairS-TO → LongPhase → ISM 三階段；v1 已棄用 | — |
| `literature_validation/` | 文獻假說 vs ISM 實證驗證 | 2026-04-12 | 已完成 | L1/L4 NEG, L3 FEASIBLE | 60+ 篇文獻 4 大假說 | `experiments/INDEX.md` |
| `seqc2_cnv_stratification/` | SEQC2 CNV 分層觀察 | 2026-04-10 | 已完成 | — | CNV + Coverage 交互分層 | `experiments/INDEX.md` |
| `independent_analyses_20260411/` | H2009 根因 + PON 跨樣本分析 | 2026-04-11 | 已完成 | — | 3 個獨立分析文件 | `experiments/INDEX.md` |
| `vaf_alleledelta_filter_study/` | VAF/AlleleDelta 過濾研究 | 2026-04-05 | 已完成 | — | CramersV + VAF 分析，過濾策略效果 | `experiments/INDEX.md` |
| `fp_provenance/` | TO FP 來源追蹤 | 2026-03-27 | 已完成 | NEGATIVE | 98.6% FP 為 raw_absent，ISM 無法過濾 | `experiments/INDEX.md` |
| `germline_asm_analysis/` | Germline ASM 方向性分析 | 2026-04-14 | 不完整 | — | 僅 figures/（6 張圖），L1 文獻驗證副產物 | — |
| `autoresearch/` | 自動化研究迴圈 | 2026-03-27 | 歸檔 | — | 自動觀察-假說-驗證循環歷史紀錄 | — |
| `_template/` | 新研究目錄範本 | — | 範本 | — | 建立新研究時複製此目錄 | — |

## 目錄結構規範

每個研究目錄應包含：
```
{study_name}/
├── README.md        # 研究概述與索引
├── figures/         # 圖表（gitignore，不追蹤）
├── data/            # 中間數據（gitignore，不追蹤）
├── scripts/         # 分析腳本（git 追蹤）
└── reports/         # 研究報告
```

## 導航

- **當前焦點**：`docs/CURRENT_FOCUS.md`
- **實驗歷史**：`docs/experiments/INDEX.md`
- **完整研究推論鏈**：`docs/reports/research_landscape/00_INDEX.md`
