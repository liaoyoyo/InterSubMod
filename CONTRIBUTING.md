<!--
建立時間: 2026-08-13
目標: 定義 InterSubMod research handoff snapshot 的可審核貢獻流程
處理範圍: 程式、測試、文件、公開 claims、資料治理與跨機可重現性
關聯檔案:
  - InterSubMod/AGENTS.md
  - InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md
  - InterSubMod/SECURITY.md
狀態: RESEARCH_HANDOFF_SNAPSHOT / RELEASE_BLOCKED
-->

# Contributing to InterSubMod

InterSubMod welcomes reproducible research and software contributions. The repository is
currently a **non-production, non-clinical research handoff snapshot**, not a release-ready
product. Do not use it for diagnosis, treatment, or patient-care decisions.

The scientific claim ceiling is fixed by the frozen 2026-08-01 authority bundle: **zero
confirmed cellular subclones**, **zero confirmed linear ancestry relationships**, methylation
is association-only, and allele-specific copy number, copy-neutral LOH, purity, and
multiplicity are not integrated into the frozen candidate reconstruction. Optional
`--loh-bed` input is stratification/annotation only; it does not produce CN/LOH-corrected
cellular inference. The reported 88.2579% is a model-conditional graph-shape statistic, not
biological accuracy or prevalence. LongLineage is a separate research-preview candidate;
`P3/P4/P5/P7/P8` remain **BLOCKED**. A contribution must not silently promote any of those
bounded statements.

Start with the [research handoff index](docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md),
[current focus](docs/CURRENT_FOCUS.md), and [repository governance](AGENTS.md). The handoff
index identifies the current authority, superseded material, machine profiles, and release
blockers.

## Before opening a change

1. State which research goal (`G1`–`G5`) and task type (`A`–`F` in `AGENTS.md`) the change
   serves.
2. For a scientific or high-impact change, write the assumptions, denominator, scope,
   confounders, stopping conditions, and planned verification before implementation.
3. Use a dedicated branch and an isolated clean worktree. Do not mix active CNV/drilldown
   research, LongLineage follow-up work, or unrelated local changes into a handoff PR.
4. Search the existing experiment index and authority/superseded crosswalk before creating a
   new result. A filename containing `final` is not evidence of finality.

## Development setup

Use an out-of-tree Release build so generated files do not contaminate the source checkout:

```bash
REPO_ROOT="$(pwd -P)"
BUILD_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/ism-build.XXXXXXXX")"
cmake -S "$REPO_ROOT" -B "$BUILD_ROOT" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_ROOT" -j"$(nproc)"
ctest --test-dir "$BUILD_ROOT" --output-on-failure --no-tests=error
```

For the acceptance Python environment, use Python 3.10 and the hash-locked requirements:

```bash
PYTHON_ENV="$(mktemp -d "${TMPDIR:-/tmp}/ism-python.XXXXXXXX")/venv"
python3.10 -m venv "$PYTHON_ENV"
"$PYTHON_ENV/bin/python" -m pip install --require-hashes \
  --requirement requirements-ci.lock
"$PYTHON_ENV/bin/python" -m pytest tests -q
```

The public tiny fixture is **DEMO/SMOKE only**:

```bash
E2E_PARENT="$(mktemp -d "${TMPDIR:-/tmp}/ism-e2e-parent.XXXXXXXX")"
scripts/handoff/run_tiny_public_e2e.sh \
  --source-repository "$REPO_ROOT" \
  --revision "$(git rev-parse HEAD)" \
  --work-dir "$E2E_PARENT/run" \
  --jobs 2
```

A passing tiny fixture proves plumbing and schema behavior; it is not biological validation.
Real-data workflows must use an untracked site profile derived from
[`config/site-profile.example.json`](config/site-profile.example.json), followed by
`scripts/site/doctor`. Never commit local settings, credentials, protected data, or raw
BAM/CRAM/FASTQ/VCF run payloads.

## Evidence and data contract

Scientific or documentation claims must use a Claim–Evidence–Verdict record and include:

- numerator, denominator, unit/grain, dataset/chromosome scope, and genome build;
- producer, exact source commit, inputs, command or verification receipt, schema, and SHA-256;
- `evidence_status`, `scope`, `finality`, availability, known limits, and supersession links;
- a distinction between technical datasets and biological IDs;
- a distinction between read/molecule-level observations, candidate mutation states, and
  cellular conclusions.

The current authority is the
[`authority_manifest.json`](docs/handoff/20260813_完整研究資料與軟體交接_01/evidence/authority_manifest.json)
plus the
[`denominator_registry.tsv`](docs/handoff/20260813_完整研究資料與軟體交接_01/evidence/denominator_registry.tsv).
Do not modify frozen authority bytes in place. Corrections must be append-only adjudications
or a new explicitly versioned authority bundle.

Regular Git must not gain files over 50 MiB. Final derived files above the repository policy
threshold belong in a governed GitHub Release or external data plane with a manifest, license,
size, and SHA-256. Runtime logs, caches, local profiles, failed scratch outputs, and secrets
remain ignored/local. Archive research evidence before replacing it; do not delete it silently
or rewrite repository history as routine cleanup.

## Pull request acceptance

Keep each PR single-purpose. In its description, record inputs, exact commands, output paths,
exit codes, representative output, and the bounded claim supported by each check. At minimum:

- run the relevant unit/regression tests, then the complete CTest and Python suites when the
  change can affect shared behavior;
- run the tiny synthetic E2E for supported workflow or schema changes;
- validate edited JSON against its schema and check Markdown/HTML links for documentation
  changes;
- regenerate derived Wiki/Pages/HTML artifacts only from their declared upstream source;
- run repository hygiene and secret scanning before proposing publication;
- leave a blocker visible when bip7, bip8, hosted CI, live GitHub surfaces, license, or source
  provenance has not been verified from the exact candidate commit.

Test and suite counts must come from the current receipt, not copied from an older document.
Do not describe a local or partial PASS as a release PASS. A tag or GitHub Release is forbidden
while any aggregate release gate is blocked.

### Release-agent manifest command

Only a release agent should run this after committing all source metadata and returning the
worktree to a clean state. Outputs must remain outside Git so the manifest can bind the exact
immutable source commit without self-reference. The command below intentionally produces a
**blocked candidate**; it does not authorize a tag or Release:

```bash
REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"
: "${PYTHON_ENV:?Set PYTHON_ENV to the hash-locked Python 3.10 venv described above}"
RELEASE_STAGE="$(mktemp -d "${TMPDIR:-/tmp}/ism-release-stage.XXXXXXXX")"
"$PYTHON_ENV/bin/python" scripts/handoff/build_release_manifest.py \
  --repo-root "$REPO_ROOT" \
  --revision "$(git rev-parse HEAD)" \
  --status BLOCKED \
  --blocker BIP7_DATA_PREFLIGHT_BLOCKED \
  --blocker BIP8_DATA_PREFLIGHT_BLOCKED \
  --blocker PUBLIC_CLAIMS_RELEASE_BLOCKED \
  --blocker GITHUB_LIVE_SURFACES_BLOCKED \
  --blocker LONGLINEAGE_PREVIEW_BLOCKED \
  --blocker INTERSUBMOD_PROJECT_LICENSE_CONFIRMATION_REQUIRED \
  --longlineage-preview-safety-status BLOCKED \
  --longlineage-candidate b9aaa12a11fa00606bd174dabd0f172a5d112359 \
  --output-json "$RELEASE_STAGE/research-handoff-manifest.json" \
  --output-checksums "$RELEASE_STAGE/SHA256SUMS"
```

Expected result while the current blockers remain:

```text
[RESULT] release_ready=false ... blockers=6 source_commit=<40-HEX-COMMIT>
```

`APPROVED_RESEARCH_HANDOFF` is fail-closed. It requires one schema-validated, closed
`aggregate_release_acceptance` receipt bound to the exact source commit, all 20 required
gates present and `PASS`, zero unknown gates, and zero blockers. Those gates include the
frozen 19/19 authority replay, registry/package and repository-hygiene checks, tracked-large
asset policy, clean build/CTest/Python/E2E, bip7 and bip8 acceptance, 84-case HTML QA,
claim/link validation, live GitHub publication alignment,
branch protection/required CI, full-history secret scan, Release-asset checksum round-trip,
six-question reader/AI acceptance, maintainer-confirmed project SPDX, and a separate
LongLineage preview-safety receipt. Authority and denominator bytes must match their pinned
2026-08-01 SHA-256 values. LongLineage must bind exactly
`b9aaa12a11fa00606bd174dabd0f172a5d112359`; its receipt may verify **preview safety only**,
while production readiness remains false and `P3/P4/P5/P7/P8` stay `BLOCKED`. Never edit
generated JSON or `SHA256SUMS` to simulate an approved state.

Every non-special aggregate gate must have a gate-specific semantic receipt plus its original
execution log/report supplied with `--gate-evidence-asset` and `--gate-log-asset`. The receipt,
aggregate gate, and declared log asset must all name the same `release-assets/...` URI and
SHA-256; a self-reported `exit_code: 0` or unbound log hash is rejected. Reader acceptance also
requires the tracked detailed receipt and re-executes the commit-pinned reader validator over
the exact 15-file package. LongLineage verification additionally requires
`--longlineage-repo-root <CLEAN_LONG_LINEAGE_CLONE>`: the builder resolves the exact safety
commit/tree, compares each copied support file with `git show` bytes, and requires the exact
commit to be advertised by the recorded `origin` branch. A local commit with only a matching
origin URL is insufficient.

The builder always includes its immutable core Git inventory: license and release metadata,
the main handoff/AI-reading path, frozen authority and denominator records, 19/19 replay,
full claim evidence, registries and schemas, portable-site tooling, CI lock/workflow, and the
tiny fixture/E2E validators. `--git-asset` only appends paths; it cannot replace that core.
PASS-only reader and LongLineage support receipts are conditional approved-candidate assets;
they are intentionally absent from the blocked inventory until independently produced. Static
algorithm-crosswalk `*_ASSET_READY` gates only establish bounded asset completeness; live
main/Wiki/Pages and repository release approval remain external aggregate gates.
Additional GitHub Release material uses
`--release-asset ARTIFACT_ID::[PUBLISHED_NAME=]PATH` and is accepted only when the artifact
registry marks the exact bytes `FINAL_FOR_SCOPE`, `AUTHORITY` or `VALIDATED_DERIVED`, licensed,
and available through `GITHUB_RELEASE`. Its size, SHA-256, and public URI must match, and the
aggregate receipt must bind the same bytes in both Release-asset round-trip and secret-scan
evidence. Raw sequencing payloads, indexes, runtime logs, databases, and secret-like files are
rejected even if a registry row exists.

## Code and documentation style

- C++17, `.hpp` headers, namespace `InterSubMod`, and the repository `.clang-format`.
- English source-code comments; reader-facing documentation may be English or Traditional
  Chinese, but scientific terms and units must remain unambiguous.
- Prefer portable paths and `<DATA_ROOT>`-style aliases in reusable commands. Absolute bip7/
  bip8 locations belong only in the governed InterSubMod machine-path registry or untracked
  site profiles.
- Treat `tree.nwk` as a read dendrogram, not a clone phylogeny, and index the current
  `significance_summary.csv` by column name.

## Issues, security, and licensing

Use GitHub issues for reproducible bugs and bounded research proposals. Do not place secrets,
exploit details, or sensitive data in a public issue; follow [SECURITY.md](SECURITY.md) instead.
The tracked [LICENSE](LICENSE) file contains the GNU GPL version 3 text, but this audit does
not infer a project-level `GPL-3.0-only` versus `GPL-3.0-or-later` grant from that text alone.
At this snapshot, three research solver sources declare `GPL-3.0-only`, while registry
metadata contains `GPL-3.0-or-later` declarations. The project-level SPDX choice therefore
has status **NEEDS_MAINTAINER_CONFIRMATION** and remains a release blocker. A maintainer must
choose the intended expression and synchronize `CITATION.cff`, source headers, registry
metadata, notices, and the aggregate license attestation before approval. Contributors must
not assume either expression while this conflict remains. Third-party code, data, and media
must retain their own license and source-origin notices; include only material you are
authorized to redistribute.
