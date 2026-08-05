<!--
建立時間: 2026-07-22
目標: 對 recovery v4 七個 frozen source 執行外部 Claude Opus 唯讀正式審查
處理範圍: bootstrap single-FD、sign-attempt retirement、atomic bundle terminal recheck
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/pre-decision-audit-recovery.md
-->

You are the external formal reviewer for InterSubMod Task-B transition-relation recovery v4.
This is READ-ONLY. Never edit, write, create, move, chmod, or delete files. Never run git
mutation commands. Use only Read, Grep, Glob, and read-only Bash.

Repository: `/big7_disk/liaoyoyo2001/InterSubMod`

Review exactly these mode-0444, nlink-1 frozen sources under
`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/`:

1. `validate_tumor_ref_schema_recovery_authority_v4.py`: size `52975`, SHA-256
   `e111b684146f7f263c810184c7b84cb5017585004a26c0e3088ea01e9a1e67ca`
2. `verify_tumor_ref_receipt_promotion_recovery_v4.py`: size `125930`, SHA-256
   `c576b7a80bab14d6411ed48317058bdccc9b38a93906b936b2f600b314e66182`
3. `replay_m2v5_runner_only_gates_recovery_v4.py`: size `131056`, SHA-256
   `64346e7b2c40b740ba5023436cd3cb4b0e960240e8bcb539f693e04ffcac28bd`
4. `continue_m2v5_after_tumor_ref_promotion_recovery_v4.py`: size `290849`, SHA-256
   `ed4510815188535c6400ded3f567e09f8c8a20519065b7f74017bbbcf46fd5ea`
5. `probe_tumor_ref_schema_recovery_sources_v4.py`: size `22692`, SHA-256
   `f24ae024ea65bab3b2522aa1089b7abe22bd6a63b469e495fbf427fe858cd7fa`
6. `schema_recovery_tests/test_tumor_ref_schema_recovery_v4.py`: size `30541`, SHA-256
   `aff345463e48af8d7c35a8970354c93a5571d370e092e78165df8b916b1c2f53`
7. `build_tumor_ref_schema_recovery_authority_v4.py`: size `29939`, SHA-256
   `411d474b43c65ae7ecbcbaeb1e368a6c562ad430d09a1ffc36a8359e54347b8a`

Run this exact probe and confirm zero output writes:

```bash
env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/probe_tumor_ref_schema_recovery_sources_v4.py
```

Required pre-review result: exit `0`, `pass=true`, `no_output_writes=true`, `39 passed`,
50 forbidden output slots absent, seven reviewed sources, ten prior inputs, and
`review_evidence_state=all_absent`.

History is append-only. Signed recovery v1 produced authority and Vrec, then Rrec v1
failed before any Rrec output because it rebound a signed historical eight-key
`during_execution` record as current live identity. Unsigned recovery v2 was formally
rejected because the builder performed a second unhashed validator pathname read, lacked
terminal rechecks, and accepted ctime inequality. Unsigned recovery v3 was formally
rejected because a post-sign exception could escape without retiring the key, the
validator source was reopened after bootstrap, and staged member/path identities were not
rechecked immediately before rename. Confirm v2/v3 authorities, reviews, runtime outputs,
and staging directories remain absent; both private keys are mode `0000`; rejection
evidence hashes are respectively
`f0c64e46ab34872149ca0cae830e79f48d818b591f8100160ac0cce44ed5ba84` and
`c31009328e2130449422e01a5fd766446f6673e6f2625be0f4f380e3f41e4ef5`.

Verify all v4 corrections and unchanged contracts:

1. Builder opens the validator exactly once with `O_NOFOLLOW`; the retained bootstrap FD
   supplies hash, compile/exec, and the authority source record. The validator is excluded
   from any later `_records` pathname reopen.
2. `signature_state.retirement_armed` is set before invoking OpenSSL. Any exception after
   the signing attempt begins retires the active v4 key and preserves staging evidence,
   even if the parent cannot observe a normal subprocess return.
3. Immediately before `renameat2(RENAME_NOREPLACE)`, the builder rechecks the staging
   directory path-to-FD inode, mode, nlink and exact inventory; each member path-to-FD
   inode, bytes, mode and nlink; authority signature; and commit binding. There is no
   sequential-file or cleanup fallback.
4. Recovery, legacy, original chain, signed v1 predecessor, both rejected-generation
   evidence records, reviews, public key, and retired keys remain descriptor leased with
   terminal stat/hash rechecks before signing and before atomic publication.
5. Preflight performs no signing or writes and leaves the v4 private key mode `0400`.
   Successful publication retires it to mode `0000`; all consumers require the exact
   three-member mode-0700 final bundle and runtime validator approval.
6. The historical exemption remains only at exact pointer
   `evidence.focal_source_identity_transition.during_execution`: exact eight historical
   keys, no inferred link count, and strict-forward historical `ctime_ns < current
   ctime_ns`. Current records remain exact nine keys, descriptor bound, mode `0444`, nlink
   `1`; equal/reversed ctime fails.
7. Vrec/Rrec/Crec use v4 paths and distinct `.recovery.v4` outputs. Signed v1 is the sole
   predecessor; v2/v3 are rejection evidence only. No old output is overwritten.
8. Tests preserve the v1 walker canary, v2/v3 exploit reproductions, strict JSON/type
   failures, sign-attempt retirement, bootstrap single-open, and staged directory/member
   mutation failures. Confirm all 39 tests cover the asserted contract without circular
   self-approval.
9. Builder cannot synthesize or transform reviews and accepts only three pre-existing,
   exact, mode-0444 APPROVE JSONs. Confirm no review fabrication, circular authority,
   exploitable TOCTOU, partial-as-complete publication, or post-sign unsafe retry remains.
10. Confirm sealed memfd, live-parent binding, exact-child normal exit zero, detached
    terminal signature, continuation-key retirement, incident fail-closed, runner lines
    1-358 only, downstream inventory, seven datasets/469,849 sites, and the biological
    claim ceiling are unchanged.

Canonical review bindings:

- `reviewed_source_set_sha256`: `cd6a9f55f9086bdc94a6f3fad9b094680ea56427ae188a0bb4ac318804ae4e4d`
- `legacy_source_set_sha256`: `92e6b11454fe04261dee8a281f1d20f23bdb131e400a5507814407bbb1ad82c3`
- `prior_recovery_chain_sha256`: `2800646ccbcedaf66aa9c6977c2ae48076c9a05d2e29f5255498352695133052`
- `rejected_generations_sha256`: `a2fccf38bb6591d707441d3c1e47c23244aa66e5509380723bf37abe9a7eecdf`
- `scope_sha256`: `366d32dd06327c1bb3cdce81efe05cac02f684468e3ddf74da9ce566af5fc32b`
- trusted public key SHA-256:
  `ee4172dae5b68fcecba93379e519a8a987f57cbde34fd6090588340afb1f0f6a`

Verdict rule: `APPROVE` only if probe exit `0`, `pass=true`, no writes, high findings `[]`,
medium findings `[]`, and unresolved conditions `[]`; otherwise `REJECT`. Low findings may
be recorded only when demonstrably non-blocking.

Return only the caller-supplied structured JSON schema. Use reviewer
`External Claude Opus`, model `claude-opus-4-8`, session ID
`a5b09aac-1f27-4eac-990d-9d7bd293e28e`.
