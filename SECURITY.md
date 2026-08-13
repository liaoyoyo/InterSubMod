<!--
建立時間: 2026-08-13
目標: 定義 InterSubMod 非 production research snapshot 的安全通報與支援邊界
處理範圍: 原始碼、workflow、依賴、憑證、路徑與資料外洩；不取代科學驗證
關聯檔案:
  - InterSubMod/CONTRIBUTING.md
  - InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md
狀態: RESEARCH_HANDOFF_SNAPSHOT / RELEASE_BLOCKED
-->

# Security Policy

## Support boundary

InterSubMod has **no production-supported or clinical version**. The current repository state
is a non-production, non-clinical research handoff snapshot and remains release-blocked.
Security fixes may be applied to the default branch on a best-effort basis, but there is no
response-time or backporting SLA.

| Source line | Security support | Intended use |
|---|---|---|
| Current default branch | Best-effort fixes after triage | Research and reproducibility only |
| Research handoff candidate | Immutable audit target once approved; currently blocked | Non-production handoff |
| Historical or superseded snapshots | No routine backports | Provenance/reference only |
| LongLineage preview | Separate security process; `P3/P4/P5/P7/P8` BLOCKED | NOT_READY / non-production |

This policy does not certify scientific validity. The frozen 2026-08-01 InterSubMod authority
supports **zero confirmed cellular subclones** and **zero confirmed linear ancestry
relationships**; methylation is association-only. Allele-specific copy number, copy-neutral
LOH, purity, and multiplicity are not integrated into the frozen candidate reconstruction;
optional `--loh-bed` input is stratification/annotation only and does not produce
CN/LOH-corrected cellular inference. The 88.2579% value is a model-conditional graph-shape
statistic rather than biological accuracy or prevalence.

## Reporting a vulnerability

Do not disclose an unpatched vulnerability, credential, sensitive path content, or exploit in
a public issue.

1. If GitHub shows a **Report a vulnerability** button on this repository's Security page, use
   that private reporting channel.
2. If private vulnerability reporting is unavailable, open a minimal public issue that says
   only that you need a private maintainer contact channel. Do not include technical exploit
   details, secrets, or sensitive data in that issue.
3. Include the affected immutable commit, operating system, compiler/runtime versions, minimal
   reproduction, expected/observed behavior, and impact once a private channel is established.

No project security email is asserted here because no repository-authoritative security
contact address has been designated. Maintainers will acknowledge and coordinate remediation
on a best-effort basis; silence must not be interpreted as an embargo deadline or permission
to expose credentials or protected data.

## In-scope security concerns

- memory-safety, command-injection, path-traversal, archive-extraction, or unsafe temporary-file
  behavior in tracked code and supported scripts;
- symlink or path handling that could overwrite data outside an explicitly selected output;
- credentials, tokens, local settings, or sensitive data accidentally committed to Git,
  generated documentation, CI logs, Pages, Wiki, or Release assets;
- dependency or CI workflow compromise, including unpinned executable actions;
- checksum, manifest, or provenance verification that can be bypassed while reporting PASS;
- a supported `--plan`, preflight, or validation command producing unexpected side effects.

Scientific disagreements, changed biological interpretations, missing benchmarks, and model
limitations are important but are not security vulnerabilities by themselves. Report them as
research or documentation issues with Claim–Evidence–Verdict records. A SHA-256 match proves
byte identity only; it does not prove code safety, scientific correctness, or clinical fitness.

## Data and credential incidents

This repository must not contain raw BAM/CRAM/FASTQ payloads, protected human data, private
keys, tokens, passwords, or machine-local profiles. If a real credential is exposed, revoke or
rotate it first. Preserve the incident evidence in a controlled archive, then coordinate any
necessary history cleanup as a separate reviewed operation; do not assume deleting the latest
file removes it from Git history.

Public research-source status does not automatically grant redistribution rights. Before
publishing a data or binary asset, verify its source, license, permitted redistribution,
intended public URI, size, and checksum. LongLineage assets remain outside this repository's
release boundary until its source-origin, license, dependency, and P8 gates pass.

The repository contains a GNU GPL version 3 license text, but its project-level SPDX choice
between `GPL-3.0-only` and `GPL-3.0-or-later` is not yet authoritatively resolved. This is a
release-governance blocker, not a security vulnerability; maintainers must reconcile source
headers, registries, citation metadata, and notices before approving a public handoff.

## Safe verification

Reproduce issues in a disposable checkout with synthetic data whenever possible. Do not run
proofs of concept against clinical systems or datasets you are not authorized to access. Use
out-of-tree builds and explicit temporary directories, and review commands before granting
network, write, or deletion privileges.
