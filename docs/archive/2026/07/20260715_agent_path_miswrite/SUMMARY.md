<!--
建立時間: 2026-07-15
目標: 封存誤寫於 repo top-level 的 matched-normal 與 strict methyl 實作/測試
處理範圍: scripts/analyze_matched_normal_candidate_controls.py 與兩個 top-level tests
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/analyze_matched_normal_candidate_controls.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_matched_normal_candidate_controls.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_strict_methyl_candidate_confirmation.py
-->

# Agent Path Miswrite Archive

2026-07-15 的修正曾誤寫至 repo top-level `scripts/` 與 `tests/`。實質變更已合併到正式 cooccurrence-validation topic；本目錄保留合併前版本供稽核，未刪除原始內容。

## 封存對照

| 原路徑 | 封存檔 | 正式維護路徑 |
|---|---|---|
| `InterSubMod/scripts/analyze_matched_normal_candidate_controls.py` | `archived_top_level_analyze_matched_normal_candidate_controls.py` | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/analyze_matched_normal_candidate_controls.py` |
| `InterSubMod/tests/test_matched_normal_candidate_controls.py` | `archived_top_level_matched_normal_candidate_controls_tests.py` | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_matched_normal_candidate_controls.py` |
| `InterSubMod/tests/test_strict_methyl_candidate_confirmation.py` | `archived_top_level_strict_methyl_candidate_confirmation_tests.py` | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_strict_methyl_candidate_confirmation.py` |

## Discovery 保護

- top-level analyzer 現為 deprecated CLI redirect，不保留第二份分析邏輯。
- 兩個 top-level test redirect 使用 `.ARCHIVED.md`，不符合 pytest Python test pattern。
- archive test 檔已改名，不以 `test_` 開頭，避免全 repo pytest 重複收集。
