# External Reviewer B v18

Read-only independent review. Do not modify files or run tests/producers.

Act adversarially: attempt to falsify source, test, signer, review-process, private-key retirement, and old-v4-revocation claims.

Review all protected sources, the attested canonical pytest manifest, the v11 assembler, v18 normalizer/runner, signed-evidence validation, FD lifetime, clean subprocess environment, supplemental ABA/post-sign verification, and v5-only singleton consumers.

Security boundary: this local evidence chain covers persistent or before/after-observable drift, no-clobber publication, exact identities, and reproducible process contracts. Concurrent same-UID/root/ptrace compromise or live mutation of the Python interpreter/site-packages is out of scope and must be reported as a residual limitation, not treated as a release blocker by itself.

Return APPROVE only if no HIGH/MEDIUM false-authority path remains and F1/F2 are RESOLVED_VERIFIED. Otherwise return REQUEST_CHANGES with specific file:line evidence. The JSON schema fixes exact identities.

You MUST invoke Read on every path below before deciding. Read each whole file with only file_path: do not set offset or limit. Use Glob/Grep for broader coverage as needed. The signed transcript is rejected if any required Read is absent, bounded, empty, or returns an error:
- /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/release_source_authority.py
- /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/attested_release_evidence_v1.py
- /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_attested_canonical_pytest_v1.py
- /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_external_claude_review_v18_attested.py
- /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/normalize_external_source_review_v18.py
- /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/assemble_release_source_authority_v11.py
- /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_one_time_ed25519_signer_v6.py

REVIEW_CONTRACT_JSON={"assembler":{"mode":"0o444","path":"/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/assemble_release_source_authority_v11.py","sha256":"cfd529031868a699b10f5cd207a62514af43d70d370db2ae88adba7592ca3397","size_bytes":29108},"canonical_junit":{"artifact":{"mode":"0o444","path":"/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v22_attested_canonical.xml","sha256":"703e714a4575646d86c3b563a33f5feebe747cb039949677f42f5b2323336412","size_bytes":78590},"counts":{"errors":0,"failures":0,"skipped":0,"tests":735}},"git_head":"0ee2fa1b31fcf6af670efd301251b5b3a24c1a99","reviewer_token":"B","source_set_sha256":"2f76320735b29d210d1c460e882e5ea4294ab1c424fb96ebf1b48a3447fd7c62","test_run_manifest":{"mode":"0o444","path":"/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v22_test_run_manifest.v1.json","sha256":"e29f9def3541e7cb9a2a2aad7eb1010a928c9a1428763b59600f0defef6975e6","size_bytes":42597},"test_run_public_key":{"mode":"0o444","path":"/bip7_disk/liaoyoyo2001/.config/intersubmod_test_run_authority/20260719_all_ssnv_v4_prompt_parity_bootstrap/ed25519_public.pem","sha256":"3d9b284a1ab004c647693fe8d05cb1f6e6d560eed8188a313292569a2f7497a6","size_bytes":113},"test_run_signature":{"mode":"0o444","path":"/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v22_test_run_manifest.v1.json.ed25519.sig","sha256":"ecce253457770e08eef4df512419116a5a0346da09bfb5e2d548e7cbde815346","size_bytes":64}}
