# data_showcase — 數據展示模板

## 觸發 keyword
「結果 / 量化 / 統計 / experiment / result / cross-sample」

## narrative skeleton
Hypothesis → Data → Stats → Caveats → Implications

## 常見頁面（6-8）
1. Hypothesis box (1)
2. 大圖 + caveat 小字 (1-3)
3. Stats table with CI (1)
4. Subgroup analysis (optional)
5. Caveats / limitations (1)
6. Implications (1)

## 範例
Cross-sample F1 benchmark 結果；HPFineNGroups 7 樣本一致性

## 對應 layouts
data_main_with_caveat / kpi_dashboard / before_after_split

## counter-example
教學報告（應改 concept_walkthrough）

---

## 業界框架對應（2026-05-08 Stage 1 補強）

### IMRAD 對應（學術論文標準結構）

| skeleton 段落 | IMRAD |
|---|---|
| Hypothesis box | **Introduction**（含研究問題與假說）|
| Data + Stats | **Methods + Results** 並陳 |
| Caveats + Limitations | **Discussion**（限制與 confound）|
| Implications | **Discussion + Conclusion**（推論與意義）|

IMRAD（Introduction / Methods / Results / And / Discussion）是生醫論文標準骨架；data_showcase 即其 slide 化版本。

### SCQA 對應

Hypothesis = Situation+Question；Data+Stats = Question's Answer；Caveats = 例外；Implications = 後續 Action。

### Quarto code-chunk 對應

Quarto 把 R/Python code chunk 與 figure 嵌入同 markdown 文件，符合 data_showcase「Hypothesis 緊鄰 Data」原則 — 這是業界 reproducible research 的標準流程。

### Reference

- IMRAD 結構: https://en.wikipedia.org/wiki/IMRAD
- Quarto reproducible reports: https://quarto.org/docs/presentations/revealjs/

