# InterSubMod

**Subclonal reconstruction from single-molecule somatic mutation co-occurrence in ONT long reads.**

[繁體中文版 →](README.zh-TW.md) · [**Docs site →**](https://liaoyoyo.github.io/InterSubMod/) · [Wiki →](https://github.com/liaoyoyo/InterSubMod/wiki) · [How to run →](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

A tumour is not one cell population — it is several subpopulations carrying different
mutation sets. Knowing *which mutation came first* and *which mutations live in the same
cells* matters for understanding resistance and progression.

With short reads you only observe each locus' marginal variant allele frequency, and
recovering the joint structure from those marginals is a **provably non-identifiable**
deconvolution problem. ONT long reads change the problem: a single molecule can span
several somatic mutations at once, so *"are these two mutations in the same cell lineage?"*
becomes a **direct observation** rather than an inference from frequencies.

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
| Engines | **`inter_sub_mod`** (C++17) · `longlineage` (C++17) | per-region statistics · per-read artefacts |
| Presentation | Python analysis + standalone HTML | figures, funnels, interactive review workstations |

---

## The core methodological claim

The single most important design decision in this project is **what is allowed to drive
the reconstruction** — and methylation is deliberately *not*.

![Why methylation cannot drive reconstruction](docs/images/methylation-circularity.png)

When you observe two groups of reads with different methylation at a locus, there are at
least four possible causes: germline allele-specific methylation, unmasking by loss of
heterozygosity, copy-number dosage, and genuine lineage difference. **A single bulk sample
cannot separate them.** Using methylation to "confirm" a subclone therefore requires
already knowing the subclone assignment — which is the very thing being proven.

So the backbone is **somatic mutation co-occurrence on the same physical molecule**, which
depends on no inferred label and is therefore non-circular. Methylation is retained as a
strictly *bounded auxiliary* signal: it is computed **after** the tree is fixed by genetic
evidence, it annotates, and it **cannot move a single edge**.

> Empirically this restraint is justified: of 811 evaluable methylation units, only
> **3 (0.37 %)** reached a robust association. Had it been used as the backbone, it would
> have been carrying almost no signal.

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
| Single sites with no co-occurrence partner (cannot form a tree) | 170,131 (66.52 %) |
| Analysis units carrying mutations | 85,941 |
| Abandoned at the search-node guard | 10,717 (12.47 %) |
| Units rankable by read-AF | 71,955 |
| Resolving to a **single rooted-unlabeled topology** | 63,506 (**88.26 %** of rankable) |
| **Confirmed cellular subclones** | **0** — see below |

> **How to read the 88.26 %.** It means: *of the 71,955 units that were already rankable,*
> 88.26 % converge to one tree shape. It is a **model-conditional graph statistic**, not
> "88 % of the tumour's evolutionary history is solved" — two thirds of all mutations were
> already lost upstream as isolated single sites.

---

## Capability boundary

This project ships a machine-readable claim boundary. The canonical output literally
records `technical_all_pass = true` alongside `validation_evidence_eligible = false`:
every hash matches and all tests pass, yet the system declares the batch **not yet usable
as validation evidence**.

<table>
<tr><th width="50%">Permitted claims</th><th width="50%">Explicitly forbidden</th></tr>
<tr valign="top"><td>

- Strict read-linked **local** structure
- Complete minimal candidate families (when the family is complete)
- Recurrence-allowed Hamming-1 parent trees
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

Everything below was verified by actually running it, not by reading documentation.

| Component | Status | Notes |
|---|:---:|---|
| `inter_sub_mod` | ✅ runnable | Minimal run completes in **2.9 s**, exit 0 |
| C++ test suite | ✅ passing | **265 tests / 38 suites**, 2.06 s |
| 7-sample HP/PS sidecars | ✅ complete | 7/7 PASS |
| Layered tree-enumeration solver | ✅ runnable | This is the path producing the numbers above |
| `longlineage preflight` | ✅ runnable | Validates the 8-role manifest |
| `longlineage dataset-gate` | ⚠️ constrained | Only path yielding science output — but **hardcoded to one dataset** |
| `longlineage run` / `probe` | 🔴 blocked | `KernelBlocked` by design |
| `longlineage` topology output | 🔴 0 units | See [caveat](#a-note-on-longlineage) |
| Writing a tagged BAM | 🔴 not supported | Neither engine can write BAM; forbidden by contract in LongLineage |
| Single script spanning both engines | 🔴 none | The two lines are currently independent |

### A note on LongLineage

`longlineage` is a clean-room industrial reimplementation (a *different git root* — not a
fork). Its contracts are excellent: every artefact is schema-locked with SHA-256 receipts,
and `topology_unit` splits "how far did we get" into four independent state fields with
named abstention reasons.

But **on real data it currently emits 0 topology units**, and this is *not* a bug:

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

# 2. Verify the build   -> expect "265 tests from 38 test suites ... PASSED"
./build/bin/run_tests

# 3. Python dependencies
pip install -r requirements.txt

# 4. Run on a single locus (~3 s)
./build/bin/inter_sub_mod \
  --tumor-bam data/bam/HCC1395/tumor.bam \
  --reference data/ref/hg38.fa \
  --vcf       one_snv.vcf \
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

Two more that fail *silently*: `methylation.csv`'s first column is a row index rather than
a read name (the binding to reads is positional, with no key check), and
`significance_summary.csv` has a **different column count across binary versions** with no
version field — always index by column name.

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
> fully self-contained HTML at **[liaoyoyo.github.io/InterSubMod](https://liaoyoyo.github.io/InterSubMod/)**
> — richer typography, interactive fold-out sections, and 29 hand-drawn diagrams rendered
> inline as SVG. Both are generated from `docs/explain/`, which is the single source.

Pages 01–10 cover the scientific method itself (glossary, ISM core, methylation read/filter,
worked case studies, statistical division of labour, capability vs. narrative).

---

## Design notes worth stealing

Two patterns in this repository generalise beyond genomics.

**1 · Streaming instead of materialising.** Haplotagged BAMs for 7 samples would occupy
**1.67 TiB**. The tagged stream is instead piped through a named FIFO and reduced on the fly
to a 9-column sidecar — **5.83 GiB, a 287× reduction** — because the analysis only ever needs
"which read, where, which haplotype". Sequence and quality strings were >99 % of the volume
and 0 % of the use.

![Upstream chain](docs/images/upstream-toolchain.png)

**2 · Make fabrication structurally impossible.** The report generator takes a declarative
spec and **refuses to render (exit 3)** when a required metric is missing — it does not
emit a dash or a blank. A missing number cannot be silently dressed up as a present one.

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

## Status of this document

Every number, command and file format in this README was verified on **2026-08-06** by
executing the commands and reading the source. Figures are generated from
`docs/explain/` by `tools/extract_svg_for_github.py` — regenerate them after editing the
source pages rather than editing the images.

Known gaps, stated honestly: copy-number is `NOT_INTEGRATED`; LongLineage's 7-sample runtime
and memory ceiling have **never been measured** and its own documentation forbids
extrapolating from one sample.
