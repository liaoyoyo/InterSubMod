<!--
建立時間: 2026-07-17T11:55:00+08:00
目標: 保存 M2 v5 方法與 final dataset/report release 接縫的獨立唯讀 agent 稽核及逐項處置
處理範圍: M2 asymmetric gate、denominator、planning levels、final dataset -> report contract、回歸測試
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/independent_m2_gate_recount.v3.json
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/build_all_ssnv_final_report_dataset.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/build_all_ssnv_report_artifact.py
-->

# M2 方法與 release 接縫獨立 agent 稽核

## 稽核身份與原始 verdict

- Agent：`019f6e03-b76d-73c0-8d93-2b07bec23c46`（Dewey）。
- 角色：唯讀方法 reviewer；未修改 repo。
- 原始 verdict：**NO-GO，2 High + 2 Medium**。此 verdict 針對稽核當時的 downstream release code，
  不否定已守恆的 M2 重算數字。
- 已獨立確認：`919 + 948 + 100,974 + 1 = 102,842`；121 是948個evaluable FAIL中的子集；
  `919/1,867=49.2234%`只能稱M2 screen-gate eligibility among evaluable M1 dataset-sites。

## Findings 與處置

| Severity | Finding | 處置 | 驗證狀態 |
|---|---|---|---|
| High | final dataset輸出6項`downstream_checks`，report validator只接受舊4項，必然觸發`M2 axis-statistic provenance contract drift` | report validator改為exact要求6項，包含assignment constant-axis proof與asymmetric positive/negative-power decision | CLOSED；完整final/report定向測試88 passed |
| High | final denominator與report methods仍使用「所有低power軸皆NOT_EVALUABLE」的舊對稱語意 | denominator改為aligned confound / powered non-alignment / assignment-proven constant三分法；methods明示陽性confound不因低power失效，並顯示121個實際位點 | CLOSED；專門n=140/最低n=152案例通過 |
| Medium | claim contract未明載categorical planning level ceilings 7/5/2，外部重作可能把observed levels誤當planning levels | 不修改已SHA鎖定的v5；改在final machine contract與report provenance結構化保存7/5/2，report validator強制exact，並明示observed levels只作constant-axis proof | CLOSED BY RELEASE SUPPLEMENT；planning drift故障注入通過 |
| Medium | 舊fixtures只提供4項checks；aligned-low-power測試讓其餘軸也低power，未覆蓋真實release接縫 | fixtures改為6項；新增缺asymmetric check fail-closed、planning ceiling drift fail-closed、單一HP-exact低power但其餘七軸可判的隔離案例 | CLOSED；相關完整測試305 passed |

## 修正後 claim ceiling

- M2支持：八個預定 measured axes 下，仍通過gate的 residual robust epigenetic partition候選。
- M2不支持：biological prevalence、排除未量測confound、cellular clone、clone數或lineage。
- Effect與499-permutation p仍來自frozen source-locked screen；terminal重算分類、n、p-grid、power與
  assignment category proof，但不是raw-read axis statistic的獨立再估計。

## 尚待完成的 release gates

此文件建立時，fresh cooccurrence producer、recovered tumor-REF、strict、matched-normal、CN/CCF、final dataset、
portable HTML及result-level external reviews尚未全部完成。因此本文件只關閉上述方法/code findings，不能單獨作
Task Type B最終發布證據。
