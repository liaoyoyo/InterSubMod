<!--
created_at: 2026-07-23
goal: Preserve read-only syscall evidence for the v18 formal C child failure.
scope: Diagnostic reproduction only; no downstream release authority.
related: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/continue_m2v5_after_tumor_ref_promotion_recovery_v18.py
-->

# v18 C child failure diagnostic staging

This directory holds `strace -ff -e write,writev` evidence from a canonical,
clean-environment reproduction after the first v18 C supervisor reported child
exit code 1 and stderr SHA-256
`035cb2d90b4db9bdf0dd490a40f75d605bd9a404d1d838e5753b837a5c60b2da`.

The reproduction is diagnostic only. It grants no release authority and must be
archived with the failed v18 generation after the failure is classified.
