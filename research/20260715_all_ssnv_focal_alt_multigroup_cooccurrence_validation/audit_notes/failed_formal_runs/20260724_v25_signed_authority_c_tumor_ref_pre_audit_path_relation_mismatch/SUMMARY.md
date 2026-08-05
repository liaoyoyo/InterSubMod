<!--
建立時間: 2026-07-24
目標: 封存 v25 已簽 authority 的 C-stage tumor-REF audit lineage 失敗
處理範圍: v25 authority/reviews/V/R/C incident/downstream observations/key quarantine
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260724_v25_signed_authority_c_tumor_ref_pre_audit_path_relation_mismatch/failure_evidence.json
-->

# v25 Signed C Failure Archive

- v25 authority、V 與 R 均 PASS；C 在 final dataset integration fail-closed。
- 直接原因：tumor-REF manifest 綁定 v1 pre-audit，builder canonical input 為 v6 pre-audit。
- v1/v6 審查同一 `102,842` stable sites、`308,526` artifacts 與同一 artifact-set SHA，
  但 v25 沒有授權這個跨代 lineage relation，因此不可 bypass。
- 次要原因：v25 continuation 的 Python 預期槽與 rendered shell 實際槽不一致。
- 沒有 final dataset/report/supplemental receipt 或 signature；科學資料與 claim ceiling 未改變。
- authority key 已退役；未使用的 terminal-v15 key 也已封存並禁止重用。
