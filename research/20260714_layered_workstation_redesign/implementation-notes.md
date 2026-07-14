<!--
建立時間: 2026-07-14 Asia/Taipei
目標: 逐步記錄 layered_workstation 全站重構的設計決定、偏離、折衷、未決與 gotcha
處理範圍: index + 7 sample pages；generator-first；canonical v5；桌機與手機
關聯檔案:
  - InterSubMod/research/20260714_layered_workstation_redesign/pre-decision-audit.md
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py
  - InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_workstation.py
cycle_id: 20260714-layered-workstation-redesign
spec_id: layered-workstation-redesign
status: complete
advisory: on
-->

# Implementation Notes: Layered Workstation Redesign

## 🔵 Design Decisions

### D-001 — Canonical v5 only

- **Status**: Accepted
- **Decision**: We will generate all 8 pages from 20260713_layered_reconstruction_v3_raw_all_lps_pass_v5; ClairS PASS v6 will be shown only as sensitivity context.
- **Rationale**: CURRENT_FOCUS and the 2026-07-14 machine summary explicitly supersede the 7/10–7/11 historical snapshot.
- **Evidence tier**: L1
- **Revisit if**: canonical run root changes.

### D-002 — Generator-first, all-page contract

- **Status**: Accepted
- **Decision**: We will modify the two tracked builders and regenerate index plus 7 ignored portable sample artifacts; we will not hand-edit numeric HTML or upstream JSON.
- **Rationale**: prevents index/sample drift and keeps every displayed number source-backed.
- **Evidence tier**: L1
- **Revisit if**: a packaged artifact runtime replaces the current offline generator.

### D-003 — Two-level reading architecture

- **Status**: Accepted
- **Decision**: We will make the index a cohort command center and every sample page a genome-to-region explorer: scope/status → whole-genome chromosome navigator → topology census → filters → selected region network → advanced evidence.
- **Rationale**: supports both a 30-second PI read and an auditable technical drill-down.
- **Evidence tier**: L3
- **Revisit if**: current-run screenshots show that a different task order is measurably clearer.

### D-004 — Topology language follows the candidate-set contract

- **Status**: Accepted
- **Decision**: We will distinguish observed read states, inferred hidden ancestors, forced versus variable edges, candidate combinations C, distinct shapes Topo, and complete versus incomplete enumeration. We will remove CCF/posterior and confirmed-clone wording.
- **Rationale**: bulk HP marginals/read AF cannot identify cellular clone K or confirm ancestry.
- **Evidence tier**: L1/L3 boundary
- **Revisit if**: orthogonal truth is joined.

### D-007 — Backbone sensitivity is part of the main interpretation layer

- **Status**: Accepted
- **Decision**: The cohort index and every dataset page visibly report the hash-bound sensitivity comparison; ClairS `FILTER=PASS` remains a sensitivity-only comparator and never changes the canonical denominator.
- **Rationale**: A `backbone_sensitive` result is interpretive context, not provenance-only metadata, so hiding it in raw JSON would overstate robustness.
- **Evidence tier**: L1

### D-008 — Candidate union and single-candidate networks have separate jobs

- **Status**: Accepted
- **Decision**: Candidate-union diagrams show edge presence/stability without visible `+S_i` labels; single-candidate diagrams show the per-edge `+S_i` acquisition labels. Repeated acquisition uses a magenta bold label plus explicit visible and accessible text.
- **Rationale**: A real H2009 eight-site tree demonstrated that combining union stability and per-edge acquisition labels produced collisions. Separation preserves both concepts without visual overplotting.
- **Evidence tier**: L2/L5

### D-009 — URL is durable workstation state

- **Status**: Accepted
- **Decision**: Discrete filters and region navigation use browser history; search typing updates the current entry; Back/Forward rehydrates filters, chromosome, selected region, detail, and focus.
- **Rationale**: Copyable deep links are only trustworthy if browser navigation restores the same evidence state.
- **Evidence tier**: L5

<!-- BEGIN USER-SPECIFIED -->
### D-005 — Preserve the whole-genome view

- **Status**: Accepted
- **Decision**: We will retain and strengthen the chr1-22 whole-genome view, with a direct path from chromosome summary to region detail.
- **DO NOT change**: the user explicitly required the whole-genome view to remain.
- **Rationale**: aggregate, canonical, outlier, and explained examples must coexist.
- **Evidence tier**: L1 user specification
<!-- END USER-SPECIFIED -->

<!-- BEGIN USER-SPECIFIED -->
### D-006 — Hide raw JSON links from the main reading flow

- **Status**: Accepted
- **Decision**: We will keep .json provenance accessible only inside collapsed evidence/source drawers; raw paths will not interrupt narrative numbers.
- **DO NOT change**: the user explicitly requested hidden JSON links in the preceding workstation/report refinement.
- **Evidence tier**: L1 user specification
<!-- END USER-SPECIFIED -->

## 🟠 Deviations

### DV-001 — Current L3 label is not pending

- **Status**: Accepted
- **Deviation**: The old UI says PENDING; canonical v5 payload says status=not_evaluated, role=bounded_auxiliary.
- **Action**: We will display the payload’s exact state and allowed/prohibited uses.
- **Evidence tier**: L1
- **Revisit if**: L3 data is actually joined.

### DV-002 — Region-level mixed-PS filtering is intentionally unavailable

- **Status**: Accepted upstream gap
- **Deviation**: Canonical v5 contains cohort/dataset mixed-PS totals but the region-view payload has no per-region PS membership.
- **Action**: Show the limitation beside filters, keep the aggregate sidecar, provide no fake PS filter, and create zero PS-derived topology edges.
- **Evidence tier**: L1
- **Revisit if**: upstream adds a hash-bound region-level PS field.

## 🟡 Trade-offs

### T-001 — Portable completeness versus page weight

- **Status**: Accepted
- **Decision**: We keep zero external dependencies and the complete region payload, split detail records into 22 embedded chromosome chunks, retain one lightweight all-genome index, and render at most 80 filtered rows into DOM initially.
- **Consequence**: portable files remain 13–52 MB, but the browser does not parse every network into DOM at startup and only decodes the selected chromosome chunk.
- **Evidence tier**: L2
- **Revisit if**: mobile load or interaction exceeds the agreed smoke thresholds.

### T-002 — Network beauty versus scientific honesty

- **Status**: Accepted
- **Decision**: We will prioritize legible state/edge semantics, direct labels, and non-color encodings over decorative network effects.
- **Consequence**: the network remains restrained and instrument-like.
- **Evidence tier**: L2
- **Revisit if**: users cannot distinguish candidate shapes in visual QA.

### T-003 — Dense union overview versus edge-level annotation

- **Status**: Accepted
- **Decision**: Union diagrams omit visible acquisition labels but retain complete parent→child/status/`+S_i` adjacency in their accessible description; exact candidate diagrams carry the visible labels.
- **Consequence**: The overview is legible while every edge remains auditable in the exact view and accessibility tree.
- **Evidence tier**: L2/L5

## 🔴 Open Questions

### Q-001 — Current-run performance ceiling

- **Status**: Resolved for browser acceptance
- **Question**: Can all 7 canonical v5 pages load and interact reliably at 390px with 11–46MB embedded payloads?
- **Falsifier**: console/page error, timeout, unresponsive filter, or greater-than-zero body overflow.
- **Resolution**: Final Chromium QA passed all 24 page/viewport contexts at 1440×1000, 390×844, and 320×720; 53/53 screenshots and 669/669 gates passed with source hashes unchanged. The independent 16-document desktop/mobile smoke also passed 365/365 checks.
- **Evidence tier**: L5 browser QA

### Q-002 — Best representative network example

- **Status**: Resolved
- **Question**: Which actual canonical region best explains C>1/Topo=1 versus C>1/Topo>1 without overloading the index?
- **Resolution**: Keep the cohort index at census/launcher level; place data-derived candidate-union and exact-tree networks in the selected region detail, where C and Topo are already bound to the correct candidate set. This avoids presenting one region as cohort-wide truth.
- **Evidence tier**: L5

## 🧪 Verification Log

### V-001 — Canonical build and freshness gate

- **Input**: `InterSubMod/research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json` and the seven canonical `layered_region_view_<sample>.json` files under `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/`.
- **Command**: `python3 InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_layered_per_sample.py` followed by the same command with `--index-only`.
- **Output**: `InterSubMod/docs/methodology/_assets/layered_workstation/index.html` plus seven portable sample HTML files in the same directory.
- **Observed result**: 7/7 pages reported `BUILT hash-bound page`, then 7/7 reported `VERIFIED hash-bound page`; summary SHA-256 was `71da78b69f8afe5fb8e618179ab7b38a6940fdb17be6282f99a6ec4b720e5de7`; both commands exited 0.

### V-002 — H2009 recurrence/incomplete overlap regression

- **Input**: `InterSubMod/docs/methodology/_assets/layered_workstation/H2009.html` served locally on Chromium.
- **Command**: Playwright selected `#ftopo=incomplete` and `#fsignal=recurrence` and read `.result-row[data-region]`.
- **Output**: browser-only regression observation; no source mutation.
- **Observed result**: exactly 4 regions — `chr8:79992384-80149786`, `chr9:275701-337149`, `chr13:93837736-93888639`, `chr15:31733893-31800487`; `expected_match=True`; page errors `[]`.

### V-003 — Generated artifact structure

- **Input**: all eight HTML files in `InterSubMod/docs/methodology/_assets/layered_workstation/`.
- **Command**: static Python/regex inspection plus `node --check` on the generated main runtime script.
- **Output**: terminal validation only.
- **Observed result**: each sample page carries the exact canonical meta SHA, one workstation payload, and 22 chromosome chunks; generated runtime syntax check exited 0. After runtime initialization, each sample has 4 JSON links and the index has 9; every link remains inside one default-collapsed evidence drawer and initial visible count is 0.

### V-004 — Pre-freeze full workstation smoke

- **Input**: layered and archived-topology workstations over `file://`, 16 documents at 1440×1000 and 390×844.
- **Command**: `python3 InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_workstation_ui.py --all --screenshot-mode off`.
- **Output**: terminal JSON report.
- **Observed result**: Chromium `147.0.7727.15`; 32 page runs; 365/365 checks PASS; 0 console/page/request errors; 0 body overflow. This was a pre-final-freeze regression and is repeated after the final visual freeze.

### V-005 — Claude Code Sonnet read-only cross-review

- **Input**: canonical machine summary, backbone comparison, renderer, driver, generated index and representative sample pages.
- **Command**: `claude -p --model sonnet --effort high --allowedTools Read Grep Glob ...` with a read-only semantic/IA/accessibility review prompt.
- **Output**: terminal review; no source mutation by Claude Code.
- **Observed result**: it found two material omissions—stored-shape scope could falsely imply unshown shapes, and backbone sensitivity was validated but not rendered. Both were corrected; the follow-up review reported no remaining P0/P1 under its review budget.

### V-006 — Iterative visual collision and hidden-control regressions

- **Input**: H2009 `chr16:65822713-65873626`, HCC1395 `chr7:57579795-57713483`, and before/after Chromium screenshots.
- **Command**: Playwright locator screenshots plus SVG `getBBox()` pairwise text-intersection checks.
- **Output**: `InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_after/comparisons/`.
- **Observed result**: candidate-union and exact network text overlap both reached 0 after task separation; one-result filtering now yields `load-more hidden=true`, computed `display=none`, remaining `0`, body overflow `0`, and page errors `[]`.

### V-007 — Final 24-context Chromium visual and interaction audit

- **Input**: `InterSubMod/docs/methodology/_assets/layered_workstation/index.html` plus all 7 dataset pages at the frozen HTML hashes.
- **Command**: `python3 InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_after/capture_layered_workstation_after.py`.
- **Output**: `InterSubMod/docs/methodology/_assets/workstation_visual_audit/20260714_layered_workstation_after/20260714_分層工作站最終視覺稽核_01.md`, `metrics.json`, and `screenshots/01_...png` through `53_...png`.
- **Observed result**: Chromium `147.0.7727.15`; 8/8 pages; 24/24 contexts; 53/53 screenshots; 669/669 checks PASS; 0 fatal error; source HTML hash changes `[]`. Scope includes 1440×1000, 390×844, and 320×720, keyboard/focus, local scrollers, URL history/deep links, candidate group semantics, acquisition-label collisions, special cases, hidden JSON, and no-network operation.

### V-008 — Final independent all-document smoke

- **Input**: current layered workstation plus archived topology workstation, 16 documents over desktop/mobile `file://` contexts.
- **Command**: `python3 InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_workstation_ui.py --all --screenshot-mode off`.
- **Output**: terminal JSON report.
- **Observed result**: 32/32 page runs; 365/365 checks PASS; 0 failed coordinates; 0 console/page/request error; exit code 0.

## 📚 Lore

- The sample HTML files are intentionally ignored because the portable artifacts total over 100MB; the tracked SoT is the builder plus index.html and README.
- region_determinacy now uses no_primary_lineage; legacy UI expects no_germline_lineage.
- H2009 has regions where recurrence and capped coexist; legacy precedence can show recurrence while candidate-level status is incomplete.
- HP3/H3?, HP4/H4?, none, and reference-only controls are auxiliary/control units, not primary parental lineages.
- PS is a QC axis only and must not become a topology edge.
