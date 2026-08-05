<!--
建立時間: 2026-07-22
目標: 外部 Claude Code 對 frozen v3 promotion/continuation chain 做唯讀正式審查
處理範圍: P/V/R/C、read-only probe、v4/v2 key lifecycle；不執行 producer/signer/downstream
關聯檔案: external_formal_readonly_probe_v1.py
-->

# External Claude Code Formal Review v3

Perform a formal, read-only release review for InterSubMod Task Type B (G4/G5).
Do not modify files and do not execute the producer, signer, promotion, replay, or
downstream pipeline. The threat model trusts the Linux kernel and the same research
account; malicious same-UID behavior is out of scope.

Review these frozen mode-0444 sources and independently verify their exact identity:

- `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/promote_tumor_ref_source_receipt_v2.py`
  - size 157669; SHA-256 `02fb9039b362fa261619b2236ddb764db278b23ae4467472fe2caa106770e06c`
- `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/verify_tumor_ref_receipt_promotion_v2.py`
  - size 120163; SHA-256 `03ff3b32368efffafa35491e04621508a46134d36407f060c3da12f90f2432a8`
- `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/replay_m2v5_runner_only_gates_v1.py`
  - size 118987; SHA-256 `10f1607aca3ef1a99da7fd77dcd6a207e0ba7003a6e3547a35b28926771fd694`
- `research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/continue_m2v5_after_tumor_ref_promotion_v1.py`
  - size 277381; SHA-256 `f7b77bd16bd86ae1cbd97e85eebb38a882998b09bc9228fa5b045abfc0ffcfbd`

Run this read-only probe and inspect its full JSON result:

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I -B \
  research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/external_formal_readonly_probe_v1.py
```

Probe identity: mode 0444, size 24050, SHA-256
`5d0c457d73547ab47810aabce647abc421ed73c2974a6257381094c45e907ed8`.

Review these failure surfaces:

1. Parent fork/no-exec producer session, exact `waitpid` raw-status decoding, and no
   parent attestation for exit 1, SIGTERM, or SIGKILL.
2. Prepare/promote exit-attestation exact schema and Ed25519 signature target; legacy
   payload signatures must remain absent.
3. P -> V -> R -> C exact schemas, strict bool/int comparison, full-stat/basic-stat
   boundaries, review/source/key cross-links, and 4+11 runtime closure.
4. FD-bound Python/finalizer commands and FD-bound signed-verifier source; confirm the
   parent checks all 22 result fields rather than a small boolean subset.
5. Key lifecycle: authorization/completion v4 keys are encrypted mode 0400 and each
   has a currently live one-time signer; continuation terminal v2 is an unencrypted
   resident mode-0400 key used only by C through bound FDs and retired to mode 0000
   immediately after signature publication.
6. Mode-000 private keys must not be passed to inotify file watches; metadata FDs and
   path/inode/mode checks must still be retained.
7. Formal review, promotion, legacy-signature, incident, and all 17 downstream slots
   must remain absent during this review.

Report only reproducible HIGH or MEDIUM findings. Return `APPROVE` only if there are
no unresolved findings and no conditions. Otherwise return `REJECT` with file/line,
reproduction path, impact, and required correction.

