<!--
created: 2026-07-23T00:54:12Z
goal: Preserve the signed v10 failure and constrain the full v11 recovery.
scope: Recovery authority, V/R evidence, C pre-write failure; no scientific payload changes.
related: failure_evidence.json
-->

# v10 Signed Recovery Failure

v10 authority, V verification, and R runner-only replay completed successfully. C failed before its first downstream write because its verification consumer expected the full failed-attempt archive object where V correctly emitted the nested receipt artifact record.

The authority, three reviews, V receipt, R receipt/log, public key, and retired private key remain immutable evidence. v10 has no terminal release authority and produced no canonical downstream artifact.

The only allowed v11 semantic correction is:

`authorization.evidence.failed_attempt_archive -> authorization.evidence.failed_attempt_archive["receipt"]`

v11 must use fresh source paths, a fresh one-time authority key, fresh reviews, and fresh V/R/C recovery slots. It must cryptographically validate this archived v10 generation and add a real-receipt producer/consumer parity regression before signing.

