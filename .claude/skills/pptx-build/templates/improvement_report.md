# improvement_report — 改進模板

## 觸發 keyword
「修補/修復/升級/優化/fix/refactor/解決」

## narrative skeleton
Before → Problem → Root Cause → Solution → Verification → Impact

## 常見頁面（10）
1. Cover with thesis (1)
2. TL;DR 大字 delta (1)
3. Pipeline 4-stage (1)
4. Before-after split (1-2)
5. Root cause tree (1)
6. Solution mechanism flowchart (1-2)
7. Verification dashboard (1)
8. Sanity check 4-row table (1)
9. Impact metric (1)
10. Future steps (1)

## 範例
v3 Self-Phasing 24-slide（4-commit 修補：8b8c1fd PON-only → 41ff147 getVote → 380e8d2 INDEL → V5 Layer 1.5）

## 對應 layouts
before_after_split / pipeline_flowchart / root_cause_tree_with_trigger / data_main_with_caveat / 5card_grid

## counter-example
不該用 improvement_report：當本週只是 status report（無明顯 before-after）→ 改 executive_summary

---

## 業界框架對應（2026-05-08 Stage 1 補強）

### SCQA 對應（McKinsey/Minto）

| skeleton 段落 | SCQA |
|---|---|
| Before（修補前狀態）| **Situation** — 既有狀態 PI 已知 |
| Problem（問題定位）| **Complication** — 觸發改動的具體 issue |
| Root Cause（根因）| **Question** — 「為什麼 X 失敗？」隱含問題 |
| Solution + Verify + Impact | **Answer** — 具體解法 + 驗證 + 量化效益 |

### Toyota A3 Report 對應

A3 Report 是豐田品質改善方法，左欄寫 problem/cause、右欄寫 countermeasure/effect，共 1 頁 A3。我們的 improvement_report 7-8 slide 結構等同 A3 的多頁化版本：
- Slide 4-5 Before-after split = A3 左右雙欄擴展
- Slide 6 Root cause tree = A3 「5 Whys」分析
- Slide 8 Verification dashboard = A3 「Confirm Effectiveness」段

### Assertion-Evidence + PLOS Rules 強制

每張 slide 標題用 assertion sentence（≤30 字）；body ≤6 elements；中 ≤60 字 / 英 ≤30 word。

### Reference

- McKinsey SCQA / Pyramid Principle: https://www.theanalystacademy.com/powerpoint-storytelling/
- Toyota A3 Report: https://lean.org/lexicon-terms/a3-report/
- Penn State Assertion-Evidence: https://writing.engr.psu.edu/ae_comprehension.pdf

