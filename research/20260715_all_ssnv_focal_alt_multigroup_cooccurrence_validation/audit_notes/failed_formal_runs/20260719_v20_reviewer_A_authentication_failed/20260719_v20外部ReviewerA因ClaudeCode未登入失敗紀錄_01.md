<!--
建立時間: 2026-07-19
目標: 保留 Reviewer A post-reboot 外部 Claude Code 認證失敗的完整證據
處理範圍: Task Type B source/test/release-chain read-only external review
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_external_claude_review_v19_attested.py
-->

# Reviewer A v20 Claude Code 未登入失敗紀錄

## 判定

- 程序：**FAIL-CLOSED；未產生 attestation，不可簽章**。
- 科學／程式 verdict：**未進行**，不得解讀為 Reviewer A 拒絕。
- 根因：Claude Code transcript 回報
  `authentication_failed` 與 `Not logged in · Please run /login`。
- 獨立 `claude auth status` 回報：
  `loggedIn=false`、`authMethod=none`、`apiProvider=firstParty`。

## 歸檔證據

- command SHA-256：
  `f1ed5acd4af97d14fe8c5a9572a094babd2e99be8eea65cc89270b9aa23f9ea8`
- request SHA-256：
  `7d1aaea6e56e8bb778651dc9e540542c4bd5f216b9e23aa7061f1327b68553c5`
- stderr SHA-256：
  `222cc8ea4f254dae65e84ccc87b1f168a0ec967e310f2e7e92202ab8a9bc448a`
- stream transcript SHA-256：
  `c39aac48ab60d57c887eca46f0eead7d11aae24351aa857e310238676b954e2d`

## 復原條件

1. 使用同一 clean HOME 完成官方 Claude Code login。
2. `claude auth status` 必須回報 `loggedIn=true`。
3. 由全新 session ID 重跑 Reviewer A；不可重用或簽署本次失敗產物。
