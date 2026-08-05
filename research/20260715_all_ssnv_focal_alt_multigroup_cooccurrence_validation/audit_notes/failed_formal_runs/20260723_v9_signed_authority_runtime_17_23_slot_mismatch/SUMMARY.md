<!--
建立時間: 2026-07-23T07:16:00+08:00
目標: 封存已簽但未通過第一個 runtime verifier gate 的 v9 recovery authority
處理範圍: v9 authority/reviews/runtime failure；科學資料與 downstream outputs 未變更
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260723_v9_signed_authority_runtime_17_23_slot_mismatch/failure_evidence.json
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/verify_tumor_ref_receipt_promotion_recovery_v9.py
-->

# v9 signed authority runtime rejection

v9 authority was atomically signed and published after three clean reviews, then the first recovery verifier rejected the existing signed authorization before creating a verification receipt. The exact mismatch was independently localized to `signed_authorization.downstream_output_absence`: the signed historical payload contains the exact 17 authorized legacy slots, while v9 incorrectly compared it with the expanded current 23-slot recovery inventory.

The two checks must remain separate: historical cross-links compare against 17 exact signed slots; live safety checks require all 23 current slots to be absent. The current inventory was absent, so this is a transition-adapter defect rather than evidence of a conflicting downstream run.

The v9 bundle and three reviews are preserved here with their hashes. Its private key is retired (`0000`), the original authority path is absent, and no verification/replay/continuation/downstream output was created. v9 cannot be reused or modified; a fresh v10 key and authority are required.
