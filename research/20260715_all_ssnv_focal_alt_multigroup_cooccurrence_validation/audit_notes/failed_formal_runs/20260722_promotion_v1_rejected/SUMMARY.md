<!--
建立時間: 2026-07-22T16:52:00+08:00
目標: 封存未執行且被獨立審查否決的 tumor-REF receipt promotion v1
處理範圍: promotion v1 source、Mendel/Nash review、未解 findings
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/completion_attempt2_archive_receipt.v1.json
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json
-->

# Tumor-REF receipt promotion v1 封存摘要

## 判定

`promote_tumor_ref_source_receipt_v1.py` 未曾執行，沒有建立 canonical receipt、promotion receipt 或簽章。來源 identity 為 size `24,505` bytes、SHA-256 `d04e9e7b5dd50cd9829bcd6d4a4351bd827f0c418b1d8d13d53e5473730f2cc9`。

兩位獨立 reviewer 均判定 `REJECT`。主要未解問題：

1. Review payload 未綁定執行中的精確 source identity 與固定 reviewer identity。
2. Canonical receipt 會早於 promotion receipt 與 detached signature 出現，可能形成可被下游消費的半完成狀態。
3. 關鍵 artifact 未維持 process-lifetime FD lease，也未使用 bound descriptor 驗章。
4. v3 approval、舊 review、attempt archive receipt 的 transitive binding 與 substantive replay 不完整。
5. 缺少 promotion 完成簽章後、任何 downstream command 前必跑的 fail-closed continuation gate。

## 後續設計要求

v2 必須採兩階段模型：先建立並簽署 promotion authorization，再以 O_EXCL promotion；promotion completion receipt 另用第二把一次性 Ed25519 key 簽署；最後由獨立 continuation gate 驗證兩組簽章、private-key retirement、canonical byte identity、v7 authority 與 runner-only gates，通過後才可執行 strict / matched-normal / CN / final dataset / report。

本封存不撤回有效的 mode-only chronology incident 判定，只否決 v1 的治理與 publication 實作。
