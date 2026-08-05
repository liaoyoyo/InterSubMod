<!--
建立時間: 2026-07-19
目標: 保留 post-reboot v24 canonical pytest 測試通過但 testcase inventory authority 過期的完整失敗紀錄
處理範圍: 7 datasets / 469,849 latest LongPhase-S recalibrated FILTER=PASS autosomal biallelic sSNV 的 Task Type B 正式證據鏈
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_attested_canonical_pytest_v2.py
-->

# v24 canonical pytest inventory authority 過期紀錄

## 判定

- 正式封裝：**FAIL-CLOSED；不可簽章、不可作 release authority**。
- pytest 實際執行：**738 passed、34 warnings、0 failed、0 error、0 skipped，return code 0**。
- 失敗原因：testcase 數量仍為 738，但凍結 inventory SHA 與實際 inventory SHA 不一致。

## 差異

- 原凍結 inventory SHA-256：`5206e2ed1a097cce773563cc3e23595631887d91ca2919e4df838a0c7747e5d6`
- 實際 v24 inventory SHA-256：`eb7afc3f13a61f14732d55e5156c602224137eb291d8766c278c09e71dbfbf98`
- v23 only：
  `test_attested_release_routes_exclude_v17_and_bind_v11_v22`
- v24 only：
  `test_attested_release_routes_bind_post_reboot_v24_and_v20_reviews`
- 結論：是一項 release-route test 的名稱與斷言語意同步更新；沒有新增、移除或跳過 testcase。

## 保留證據

- `pytest_full_pre_authority_v24_post_reboot_raw_evidence/pytest.xml`
  - SHA-256：`2f430a71b60a923eea58e2417a3f2b00a53e185de320d8d10ec360def5422ecd`
- `pytest_full_pre_authority_v24_post_reboot_attested_canonical.stdout.txt`
  - SHA-256：`c0f9d376425dd3e1e648aeec25193e3a5bede51ae0da30ef1248fbc1775b42a8`
- `pytest_full_pre_authority_v24_post_reboot_attested_canonical.stderr.txt`
  - SHA-256：`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`
- `.pytest_full_pre_authority_v24_post_reboot_bound_test_fds/`
  - 保留此次執行的 bound-FD testcase 路由記錄；程序結束後 `/proc/self/fd/*` symlink 預期失效。

## 修正與驗證要求

1. 將 canonical testcase inventory authority 旋轉為實際 SHA。
2. 同步更新所有受保護來源對 evidence validator 內容 SHA 的 pin。
3. 重新跑 targeted tests 與 738-test nonformal full suite。
4. 使用全新正式輸出槽重跑，只有新一輪通過封裝 guard 才可交給 one-time signer。
