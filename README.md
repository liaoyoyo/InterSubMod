# InterSubMod

**Local mutation-state candidate reconstruction from single-molecule somatic-mutation co-occurrence in ONT long reads.**

[繁體中文版 →](README.zh-TW.md) · [**Docs site →**](https://liaoyoyo.github.io/InterSubMod/) · [Wiki →](https://github.com/liaoyoyo/InterSubMod/wiki) · [How to run →](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

A tumour can contain several cellular populations carrying different mutation sets.
Learning their order and cellular co-membership matters for resistance and progression,
but those biological quantities are **not directly observed by an unbarcoded bulk long read**.

If the input is limited to per-locus marginal variant allele frequencies, with no linkage
or additional model assumptions, multiple joint structures can explain the same marginals;
the joint structure is therefore **not identifiable from those marginals alone**. ONT long
reads add a different observation: a single molecule can span
several called somatic mutations at once. Their co-occurrence on the **same physical
molecule** becomes a direct observation; cellular co-membership and lineage remain
model-dependent inferences because the originating cell is not identified.

This repository implements that idea end to end — and, just as importantly, it is explicit
about where the evidence runs out.

---

## Architecture

![System architecture](docs/images/architecture-overview.png)

Five layers, bottom to top. Only some of them are currently runnable — see
[Status](#status) below, which is deliberately honest about what does and does not work.

| Layer | Component | What it does |
|---|---|---|
| Raw data | tumour / normal BAM, reference FASTA | ONT long reads carrying methylation tags |
| Upstream | Dorado · ClairS · LongPhase-S · SAVANA | basecalling + methylation, somatic SNV calling, haplotagging, copy number |
| Interface | **HP/PS sidecar TSV** | one row per read; the only byte-identical contract between both engines |
| Engines | **`inter_sub_mod`** (C++17) · `longlineage` (C++17) | per-region methylation/statistics · version-scoped per-read artefacts |
| Presentation | Python analysis + standalone HTML | figures, funnels, interactive review workstations |

---

## The core methodological claim

The single most important design decision in this project is **what is allowed to drive
the reconstruction** — and methylation is deliberately *not*.

![Why methylation cannot drive reconstruction](docs/images/methylation-circularity.png)

When you observe two groups of reads with different methylation at a locus, there are at
least four possible causes: germline allele-specific methylation, unmasking by loss of
heterozygosity, copy-number dosage, and genuine lineage difference. **They are not
identifiable from the current single-bulk measurement set without orthogonal data or
additional assumptions.** Using methylation to independently "confirm" a cellular subclone
would therefore require an external cellular assignment; conditioning on mutation-defined
labels is concordant association, not independent confirmation.

So the backbone is **somatic mutation co-occurrence on the same physical molecule**, which
depends on no inferred label and is therefore non-circular. Methylation is retained as a
strictly *bounded auxiliary* signal: it is computed **after** the genetic candidate structure
is fixed, it annotates that structure, and it **cannot move a single edge**.

> Empirically this restraint is justified: of 811 evaluable methylation units, only
> **3 (0.37 %)** supported a robust pattern-conditioned association. This low yield is one
> reason not to use methylation to choose topology; it is not a test of every possible
> methylation signal.

---

## What the pipeline actually produces

![Seven-sample funnel](docs/images/funnel-7samples.png)

Canonical result across **7 datasets (6 biological samples), chr1–22**, frozen 2026-08-01.
Every layer is auditable and the arithmetic is self-consistent
(39,648 + 23,858 = 63,506; + 8,449 = 71,955; + 3,224 + 45 = 75,224; + 10,717 = 85,941;
+ 13,014 = 98,955).

| Quantity | Value |
|---|---|
| sSNV dataset records | 469,849 |
| `k=1` strict read-linkage components | 170,131 / 255,752 strict components (66.52 %) |
| Analysis units carrying mutations | 85,941 |
| Abandoned at the search-node guard | 10,717 (12.47 %) |
| Units rankable by read-AF | 71,955 |
| Resolving to one **rooted-unlabeled mathematical topology signature** | 63,506 (**88.26 %** of rankable) |
| **Confirmed cellular subclones** | **0** — see below |

> **How to read the 88.26 %.** It means: *of the 71,955 units that were already rankable,*
> 88.26 % converge to one mathematical shape under the frozen recurrence-allowed model.
> It is a **model-conditional graph statistic**, not "88 % of the tumour's evolutionary
> history is solved". Separately, 170,131 / 255,752 strict components (66.52 %) are `k=1`;
> relative to 469,849 sSNV dataset records, 170,131 is 36.21 %. These denominators are not
> interchangeable.

---

## Capability boundary

This project ships a machine-readable claim boundary. The canonical
`cohort_receipt.json` and `summary/all7_summary.json` (not the top level of
`authority_manifest.json`) record `technical_all_pass = true` alongside
`validation_evidence_eligible = false`:
every hash matches and all tests pass, yet the system declares the batch **not yet usable
as validation evidence**.

<table>
<tr><th width="50%">Permitted claims</th><th width="50%">Explicitly forbidden</th></tr>
<tr valign="top"><td>

- Strict read-linked **local** structure
- Complete minimal candidate families (when the family is complete)
- Local recurrence-allowed Hamming-1 candidate arborescences
- Deterministic read-AF ordering, **CN/LOH-uncorrected**
- Exact rooted-unlabeled topology census
- Pattern-conditioned methylation **association**

</td><td>

- **Confirmed cellular clone / subclone**
- CN/LOH-corrected CCF
- Calibrated likelihood or posterior
- Treating a read cluster as a cell clone
- Re-ranking topologies with methylation
- Treating a basecaller re-run as a biological replicate

</td></tr>
</table>

**A single bulk sample can characterise but not confirm.** Breaking that ceiling requires
single-cell or multi-region data — an information-theoretic limit, not an implementation gap.
Inferred internal nodes are labelled `inferred` rather than asserted.

---

## Status

The table is version-scoped. "Verified" means the listed artifact was run or inspected at
the stated audit version; it is not a blanket guarantee for every public command or file.

| Component | Status | Notes |
|---|:---:|---|
| `inter_sub_mod` | ✅ runnable | Fresh build/run audit at tracked core `73afaeac`; command receipt is linked below |
| C++ test suite | ✅ passing | **270 tests / 39 suites** and CTest **270/270**, fresh 2026-08-12 run |
| 7-sample HP/PS sidecars | ✅ complete | 7/7 PASS |
| Research exact-PS topology solver | ✅ runnable | Separate research solver path producing the exact-PS funnel above; not `inter_sub_mod` |
| `longlineage preflight` | ✅ runnable | Validates the 8-role manifest |
| `longlineage dataset-gate` | ⚠️ constrained | Only path yielding science output — but **hardcoded to one dataset** |
| `longlineage run` / `probe` | 🔴 blocked | `KernelBlocked` by design |
| `longlineage` topology output | 🔴 0 units | See [caveat](#a-note-on-longlineage) |
| Writing a tagged BAM | ✅ runnable | `inter_sub_mod` does not write BAM; LongLineage public main `583e03e` now builds `longlineage-tag-bam` (`CMakeLists.txt:221`), merged 2026-08-17 |
| Single script spanning both engines | ⚠️ per-engine only | LongLineage ships `scripts/run_sample.sh` (partition → tagged BAM). No single script yet spans LongLineage **and** `inter_sub_mod`; the run order is documented in [LongLineage's README](https://github.com/liaoyoyo/LongLineage#與-intersubmod-的關係先讀這節再決定要不要裝) |

### A note on LongLineage

`longlineage` is a clean-room industrial reimplementation (a *different git root* — not a
fork). Its contracts are excellent: every artefact is schema-locked with SHA-256 receipts,
and `topology_unit` splits "how far did we get" into four independent state fields with
named abstention reasons.

For the frozen **HCC1395 dataset-gate receipt** it emits 0 topology units; this does not
generalise to every possible LongLineage run:

![LongLineage funnel](docs/images/longlineage-funnel.png)

Its pipeline gates topology construction on methylation clustering stability, so 66,836 of
79,687 sites (83.9 %) never produce a co-occurrence row at all. This is a **direct conflict**
with the methodological position stated above, and any cross-engine comparison must say so
explicitly.

"Blocked" here means *the parity evidence does not exist yet* — **not** that the code is
missing. The kernels are implemented and have been executed.

---

## Quick start

Full walkthrough with expected output at every step:
[**How to run**](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

![Six steps](docs/images/howto-six-steps.png)

```bash
# 1. Build  (the checked-in binary is stale — always rebuild)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# 2. Verify the build   -> audited 2026-08-12: 270 tests from 39 suites, all passed
./build/bin/run_tests

# 3. Python dependencies
pip install -r requirements.txt

# 4. Supply your own indexed, licensed inputs. The public repository currently ships no
#    runnable BAM/FASTA/VCF fixture, so these are placeholders, not copy-paste data paths.
./build/bin/inter_sub_mod \
  --tumor-bam /path/to/tumor.mm_ml.bam \
  --reference /path/to/reference.fa \
  --vcf       /path/to/candidates.vcf \
  --output-dir out_min
```

<details>
<summary><b>Three things that will bite you</b> (all verified against source)</summary>

- **`--threads` defaults to 16, not 1.** The help text says `Default: 1`; the actual
  default in `Config.hpp` is `16`. Resource estimates will be off by 16×.
- **`--distance-metric` is silently forced to `NHD`.** `Config.hpp` declares `BERNOULLI`,
  but the argument parser unconditionally clears it. Always specify it explicitly.
- **`tree.nwk` leaves are *reads*, not clones.** It is a hierarchical clustering tree over
  read methylation similarity — **not** a subclonal phylogeny. This is the single most
  commonly misread output.

Two more version-sensitive details: `methylation.csv`'s first column is a row index rather
than a read name (the binding to reads is positional, with no key check). At audited core
`73afaeac`, `significance_summary.csv` has **199 columns**, including
`VerificationSchemaVersion=2` and `RegionStratificationSchemaVersion=1`. These are component
schema fields, not a single whole-file layout version; always index by column name and pin
the producing commit.

</details>

---

## Documentation

The **Explanation Center** (`docs/explain/`) is the entry point — 17 self-contained pages,
no external dependencies, openable offline.

| | Page | For |
|---|---|---|
| 🗺️ | [System overview](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) | **Start here** — architecture, honest status, funnel, capability boundary |
| ⚙️ | [InterSubMod I/O](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine) | 3 inputs, 8 stages, 17 output files with real headers, 9 pitfalls |
| 🧬 | [LongLineage I/O](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) | Subcommand status, M1→M2→topology chain, artefact contracts |
| 🔬 | [Upstream & data](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) | Dorado / ClairS / LongPhase-S / SAVANA, sidecar format, 7 samples |
| 📊 | [Analysis & presentation](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) | Which scripts to use, refuse-on-missing HTML generator |
| ▶️ | [How to run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run) | Six steps, each with an acceptance check |

> **Two ways to read the same material.** The table above links to the **Wiki**, which
> GitHub renders natively and is the quickest to skim. The same content is also served as
> fully self-contained HTML at **[liaoyoyo.github.io/InterSubMod](https://liaoyoyo.github.io/InterSubMod/)**.
> `docs/explain/` is the editorial upstream; the Wiki is a manually synchronized derivative,
> and publication is a separate step. At audited Pages deploy `fbdf7c7`, the 17 standalone
> pages contained **37 inline `<svg>` elements**, counted with the command documented in the
> correction receipt. Do not assume that this version-scoped element count represents 37
> distinct semantic figures or that Wiki and Pages are byte-identical.

Pages 01–10 cover the scientific method itself (glossary, ISM core, methylation read/filter,
worked case studies, statistical division of labour, capability vs. narrative).

---

## Design notes worth stealing

Two patterns in this repository generalise beyond genomics.

**1 · Streaming instead of materialising.** A tagged stream can be piped through a named FIFO
and reduced on the fly to the current 9-column sidecar contract. The seven audited sidecars
sum to **6,256,168,164 bytes (5.83 GiB)**. The previously displayed **1.67 TiB** tagged-BAM
total has no committed seven-file path/byte/hash/compression receipt, so it is **UNVERIFIED**
and no storage-reduction multiplier is claimed. The current sidecar omits `SEQ` and `QUAL`;
their byte share and downstream utility were not measured by a field-level census.

![Upstream chain](docs/images/upstream-toolchain.png)

**2 · Refuse missing required metrics.** The audited workstation generator takes a
declarative spec and **refuses to render (exit 3)** when a declared required metric is
missing. This prevents silent rendering of those declared fields; it does not validate
their truth, denominator or source, and it cannot detect omitted optional fields.

![Refuse-on-missing design](docs/images/workstation-refuse-design.png)

---

## Requirements

- **C++17** compiler, CMake ≥ 3.14, OpenMP
- **htslib** (BAM/VCF I/O), Eigen
- **Python ≥ 3.9** (`requirements.txt`); note that a few strict-partition scripts require **3.10**
- Reference FASTA **with `.fai`** (a missing index is a hard error)
- Input BAM must carry `MM`/`ML` methylation tags (reads without them are dropped)

## Repository layout

```
include/core/ src/core/     C++17 analysis engine
scripts/                    Python analysis entry points
docs/explain/               Explanation Center (17 standalone pages)
docs/images/                Figures used by this README and the wiki
tools/                      Rendering, QA and extraction utilities
state/                      Research cycle state machine
```

## Verification scope of this document

This README is bounded by the
[2026-08-12 public-claim audit](docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md).
The frozen exact-PS counts were checked against the authority manifest and denominator
registry; the tracked C++ core was freshly built and tested. The public quick-start data are
**not shipped**, GitHub About/Wiki/Pages have independent publication state, and LongLineage
capabilities are pinned to the commits named above. Therefore no blanket "every number and
command was verified" claim is made. Figures are extracted from `docs/explain/` by
`tools/extract_svg_for_github.py`; regenerate them after editing their source pages.

| Artifact / claim family | Verification identity and date | Reproducible check and result | Scope and known failure |
|---|---|---|---|
| Frozen exact-PS funnel | Frozen authority artifacts, rechecked 2026-08-12 | Manifest/hash census plus independent denominator recomputation; exact counts reproduced | Frozen 7-dataset analysis only; it is not the `inter_sub_mod` CLI and does not identify cellular clones |
| Tracked C++ core | `73afaeac-dirty`, GCC 11.4.0, htslib 1.18; run 2026-08-12 | `cmake -S . -B <build> -DCMAKE_BUILD_TYPE=Release`, build, direct GoogleTest and CTest: build exit 0; 270 tests / 39 suites; CTest 270/270 | Audited C++/CMake build inputs were byte-equivalent to remote feature `ddd8909`; the run was local-dirty, not a clean-checkout release certification |
| Public quick start | Current source, reviewed 2026-08-13 | Build/test path is reproducible; analysis invocation uses explicit user-supplied placeholders | No public tumor BAM/reference/VCF fixture is shipped, so no end-to-end biological result is claimed |
| LongLineage status | Public main `583e03e` (2026-08-17); frozen HCC1395 evidence | Merged the three `agent/public-preview-*` branches plus the tag-bam/solver commits into public `main`, then rebuilt and reran its suite in an isolated worktree: build exit 0, **ctest 49/49** | `scripts/ci/check_public_preview_gate.sh HEAD` still reports **FAIL with 5 open blockers** (license review pending, 4 `NO_GO` source rows, 21 unapproved source-license rows, 11 `NOASSERTION` dependencies, 4 historical blobs carrying developer absolute paths). Readable/buildable, **not** a license-cleared redistribution. No 7-sample runtime/memory or multi-dataset topology validation |
| GitHub surfaces | Local source corrected 2026-08-13; GitHub About is `RESOLVED_LIVE` | Live API check plus the P0/P1/P2 claim guards and source checks in this refresh cycle | Default `main`, Wiki and Pages still expose the earlier deployed bytes; they remain pending publication and must not be described as live-corrected |

The exact commands and captured output are preserved in the
[command receipts](research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md);
the per-claim correction status is tracked in the
[P0 correction cycle](research/20260813_public_docs_p0_correction/00_INDEX.md), and the
[2026-08-13 public-surfaces refresh cycle](research/20260813_intersubmod_public_surfaces_refresh/00_INDEX.md)
separates local corrections from the
[live remote-state receipt](research/20260813_intersubmod_public_surfaces_refresh/remote_state_receipt.md).

Known gaps, stated honestly: copy-number is `NOT_INTEGRATED`; LongLineage's 7-sample runtime
and memory ceiling have **never been measured** and its own documentation forbids
extrapolating from one sample.
