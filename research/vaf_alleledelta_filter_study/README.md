<!--
建立時間: 2026-04-05 18:00
目標: VAF/AlleleDelta 過濾研究索引
處理範圍: CramersV + VAF 分析、過濾策略效果評估
關聯檔案:
  - research/vaf_alleledelta_filter_study/analysis_report.md
  - docs/experiments/INDEX.md
-->

# VAF/AlleleDelta Filter Study

分析 VAF 與 AlleleDelta 特徵在 TP/FP 區分上的效果，以及 CramersV 結合 VAF 的過濾策略。

## 內容

| 檔案 | 說明 |
|------|------|
| `analysis_report.md` | 主分析報告 |
| `cramersv_analysis_report.md` | CramersV 分析報告 |
| `cramersv_definitive_sites_report.md` | CramersV 確定性位點報告 |
| `cramersv_vaf_analysis_report.md` | CramersV+VAF 聯合分析報告 |
| `heatmap_case_comparison.md` | 熱圖案例比較 |
| `kept_sites_report.md` | 保留位點報告 |
| `upgma_cramersv_improvement_evaluation.md` | UPGMA CramersV 改進評估 |

## 腳本

| 腳本 | 用途 |
|------|------|
| `analyze_filter.py` | 過濾分析主腳本 |
| `analyze_kept_sites.py` | 保留位點分析 |
| `analyze_cramersv_definitive.py` | CramersV 確定性位點分析 |
| `analyze_cramersv_vaf_combined.py` | CramersV+VAF 聯合分析 |

## 結論

此研究屬於早期探索階段，後續被 O1-O10 系統性觀察和 G1-G7 更全面的特徵搜索所涵蓋。
