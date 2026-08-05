<!--
建立時間: 2026-07-22T10:20:00+08:00
目標: 封存 positional-singleton supplemental v4 test evidence 的失效嘗試
處理範圍: v28 JUnit、v4 signed test-evidence manifest、raw pytest cache
關聯檔案: InterSubMod/docs/provenance/source_authorities/20260722_all_ssnv_focal_alt_release_source_authority.v7.json
-->

# Supplemental v4 test evidence invalidated

## Verdict

`INVALIDATED / EXCLUDED / NOT CURRENT RELEASE EVIDENCE`。

v28 JUnit 與 signed v4 manifest 曾綁定一個臨時新增的
`test_build_positional_singleton_report.py`回歸測試。fresh singleton audit 隨後正確發現該test source
不再符合v7 signed test inventory。為恢復主Task-B authority，test source已精確還原為v7記錄的
mode `0444`、size `61,445`、SHA-256
`d1d1c1db7d1c13d3a7eaeecd792571f0aa7bac8871318ba2834f052f075698c0`。

還原後，v4 manifest綁定的臨時test SHA不再等於live source，因此即使其detached signature在歷史bytes上
仍可驗證，也不能作current supplemental release evidence。相關檔案全部保留於本目錄，未刪除、未覆寫。

## Preserved Evidence

- `pytest_full_supplemental_v28_v7_authority_alignment_canonical.xml`: 744-test historical run。
- `supplemental_report_test_evidence.v4.json`: invalidated signed manifest。
- `supplemental_report_test_evidence.v4.json.ed25519.sig`: historical detached signature。
- `pytest_full_supplemental_v28_v7_authority_alignment_raw_evidence/`: raw Python cache。

v4 one-time private key已按成功ceremony退役為mode `0000`，不得重用。修復鏈使用新v5 key、未修改的v7
test inventory、fresh v29 JUnit與`logs/supplemental_report_test_evidence.v5.json`。
