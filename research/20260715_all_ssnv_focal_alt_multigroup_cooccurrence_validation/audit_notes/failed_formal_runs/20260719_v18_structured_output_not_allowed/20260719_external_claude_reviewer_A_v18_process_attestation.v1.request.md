# External Reviewer A v18

Read-only independent review. Do not modify files or run tests/producers.

Trace every authorization claim to executable checks and look for coherent substitution, TOCTOU, environment injection, or false-pass paths.

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

REVIEW_CONTRACT_JSON={"assembler":{"mode":"0o444","path":"/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/assemble_release_source_authority_v11.py","sha256":"cfd529031868a699b10f5cd207a62514af43d70d370db2ae88adba7592ca3397","size_bytes":29108},"canonical_junit":{"artifact":{"mode":"0o444","path":"/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v22_attested_canonical.xml","sha256":"fb52e80f81e78f775ac4cec32bb0f8465ef1596b047dba734a525dee39dee7d5","size_bytes":78590},"counts":{"errors":0,"failures":0,"skipped":0,"tests":735}},"git_head":"0ee2fa1b31fcf6af670efd301251b5b3a24c1a99","reviewer_token":"A","source_set_sha256":"c8bd2f49932fd5ee165c1b7ff46cf0d0b410245d52880af312c75170e874980a","test_run_manifest":{"mode":"0o444","path":"/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v22_test_run_manifest.v1.json","sha256":"72e10fe9cb4ca04daf620035dfbbc73b42b6bc10b6ac1e11e3ef8d6ee35caa16","size_bytes":42581},"test_run_public_key":{"mode":"0o444","path":"/bip7_disk/liaoyoyo2001/.config/intersubmod_test_run_authority/20260719_all_ssnv_v2_bound_bootstrap/ed25519_public.pem","sha256":"f102c6c8f03021d6a6026bc03a15cf8920083b183e5d7f3ceb51ef13ae1e1e78","size_bytes":113},"test_run_signature":{"mode":"0o444","path":"/big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/logs/pytest_full_pre_authority_v22_test_run_manifest.v1.json.ed25519.sig","sha256":"ee676965c5679b4b317a4b93ac3d795ecc6d7d7537f07bb493a4c34d59fc4929","size_bytes":64}}
