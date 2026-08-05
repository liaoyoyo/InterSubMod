<!--
建立時間: 2026-07-25
目標: 導航布爾超立方體候選樹計算量教學 HTML、可重算 evidence 與 QA receipts
處理範圍: all7 chr1-22 technical snapshot；Task Type F（教學／示意）
關聯檔案: InterSubMod/research/20260725_hypercube_computation_complexity_explainer/20260725_布爾超立方體候選樹計算量教學_01.html
-->

# Boolean-hypercube 計算量教學：索引

> [!WARNING]
> **DEMO／方法教學。** all7 數字來自 2026-07-24 recurrence-allowed C++ technical run；
> 不是 production strict directed topology，也不能作 clone count 或 cellular lineage 證據。

## Canonical deliverables

- [互動式單檔 HTML](20260725_布爾超立方體候選樹計算量教學_01.html)
- [公式、數據與 HTML 驗證紀錄](20260725_計算量公式與HTML驗證紀錄_01.md)
- [Canonical evidence JSON](data/20260725_hypercube_computation_evidence_02.json)
- [數據重算 receipt](results/20260725_data_recalculation_receipt_02.json)
- [瀏覽器 QA receipt](results/20260725_html_qa_receipt_01.json)

## Reproducible scripts

- `InterSubMod/research/20260725_hypercube_computation_complexity_explainer/scripts/build_complexity_evidence.py`
- `InterSubMod/research/20260725_hypercube_computation_complexity_explainer/scripts/build_explainer_html.py`
- `InterSubMod/research/20260725_hypercube_computation_complexity_explainer/scripts/audit_complexity_explainer.py`

## Version note

`data/*_01.json` 與 `results/*_01.json` 是第一次重算保留件；其公式欄名
`full_cube_submask_visits` 過度暗示實際拜訪次數。Canonical `_02` 已改為
`connection_lb_submask_loose_bound`，數值不變，但語意更精確。

