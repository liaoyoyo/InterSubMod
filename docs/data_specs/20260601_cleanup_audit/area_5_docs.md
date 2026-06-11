<!--
title: Cleanup Audit Area 5 — docs/
date: 2026-06-01
type: cleanup_audit_area_report
generated_by: AI scan (READ-ONLY; no deletion / move executed)
scope: /big7_disk/liaoyoyo2001/InterSubMod/docs/ (20 subdirs, ~110 MB)
-->

# Cleanup Audit — Area 5/6: `docs/`

> **READ-ONLY scan.** No file deleted or moved. Delete/move is Hard Gate, executed only by the parent agent after user confirmation.
> **最高原則**：結案 NEGATIVE ≠ 自動可刪。撐起已發佈結論/論文段落的原始證據 → 至少 ARCHIVE；任何不確定 → NEEDS_USER_DECISION。

## 0. Headline

- `docs/` total = **111,558,565 B (~106 MiB)**.
- Largest consumers: `presentations/` 52 MB (dominated by 教授報告 base+v2+v3 = 34 MB), `reports/` 32 MB, `experiments/` 12 MB, `trash/` 8 MB, `archive/` 3 MB.
- **Clearest reclaimable**: `docs/trash/to_pipeline_staging_v1/` (8.27 MB) — README self-declares "已棄用（待刪除）", wrong data source, corrected v2 confirmed present in `research/to_pipeline_staging/`.
- **Largest single-file reclaimables**: untracked `*.standalone.html` / `*.emailsafe.html` preview companions (total untracked = **9.05 MB**; most regenerable from their `.md`).
- **Anti-fragile note**: most `methyl_augmented_filter` references live in LIVE docs (CURRENT_FOCUS / INDEX) — the methyl-FP-filter **research output** is DEAD (⭐2 L4) but the **docs are G6 paper methods evidence** → KEEP docs, ARCHIVE (not delete) the research data.

---

## 1. (a) du --max-depth=1 整個 docs/ + stale/廢棄物標記

| subdir | size | verdict | note |
|---|---|---|---|
| `trash/` | 8.0 M | **SAFE_DELETE** (after user ack) | only child = `to_pipeline_staging_v1/`, README = "已棄用（待刪除）"; v2 corrected exists |
| `archive/` | 3.5 M | **KEEP** (it IS the archive) | 2025 + 2026 + deep/2025-12_old_structure;封存區，本就該留 |
| `drafts/` | 100 K | **NEEDS_USER_DECISION** | CLAUDE.md v2/v3 + AGENTS v3 drafts + `backup/` pre-v3deploy snapshots; superseded by live CLAUDE.md but tiny |
| `presentations/` | 52 M | mixed (see §1.1) | validated/2026/04 教授報告 = 34 M w/ cross-version dup figures |
| `reports/` | 32 M | mixed (see §2) | pi_reports 18 M + in_progress 8.2 M (large untracked HTML) |
| `experiments/` | 12 M | mixed | in_progress A0_assets 5.7 M (TSV matrices) + 2026/06 1.9 M HTML |
| `references/` | 1.3 M | KEEP | startup-context manuals |
| `concepts/ plans/ methodology/ ...` | <1 M each | KEEP | text-only research docs |

### 1.1 presentations/ detail

| path | size | bytes | purpose | trust_tier | conclusion | verdict | reclaimable | rationale |
|---|---|---|---|---|---|---|---|---|
| `presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告` | 34 M | 34,446,488 | 4/29 professor report (LOH/self-phasing) base + v2 (16 M latest) + v3 (4.5 M 精簡) | CURRENT (validated, V5-era) | LIVE (deliverable record) | **NEEDS_USER_DECISION** | ~7-11 M (dup figures across base/v2; base may be superseded by v2) | figures duplicated byte-identical across `base/figures` and `v2/figures` (fig01c/fig14/fig15/fig16 = 0.5 MB each ×2); validated → conservative, do not auto-delete base |
| `presentations/validated/2026/04/.../output/*.pptx` | 3.0+2.3+2.0 M | 7,319,847 | rendered PPTX (base/v2/v3) | CURRENT | LIVE | **KEEP** | 0 | final delivered artifacts |
| `presentations/validated/2026/04/.../__pycache__` | part of 110 K | 119,936 (all docs) | python bytecode cache | TRANSIENT | DEAD | **SAFE_DELETE** | 0.1 M | regenerable bytecode; 2 dirs (validated + draft 週報v2) |
| `presentations/validated/2026/05/self_phasing_synthesis_PI/preview_singlefile_archive.html` | 52 K | 53,487 | single-file archive of PI preview (user 21-commit deck) | CURRENT | LIVE | **KEEP** | 0 | memory note: this deck is user-hand-edited, do NOT rebuild |
| `presentations/draft/` | 164 K | — | early week-report drafts (週報v2_AE) | SUPERSEDED | DEAD | **NEEDS_USER_DECISION** | ~0.1 M | superseded by validated, but tiny + has __pycache__ |
| `presentations/in_progress/` | 4.2 M | — | weekly_PI 4versions + V6_signoff + BRCA2_ZAR1L decks (HTML slides) | PRE-FIX / in_progress | LIVE (active prep) | **KEEP** | 0 | active presentation prep |

---

## 2. (b) 重複/可重生的大 HTML 檔（*.standalone.html / *.emailsafe.html）

**Aggregate**: untracked standalone+emailsafe = **9,049,122 B (9.05 MB)**; tracked = **5,955,662 B (5.96 MB)**.

Rule applied: in_progress untracked `.standalone.html` **with .md companion** = regenerable → NEEDS_DECISION/SAFE_DELETE; validated/tracked終版 = KEEP; PI-report standalone **without .md companion** = primary deliverable, NOT regenerable → NEEDS_DECISION (conservative).

| path | size | bytes | git | .md companion | purpose | tier | verdict | reclaimable | rationale |
|---|---|---|---|---|---|---|---|---|---|
| `reports/in_progress/2026/05/20260523_HKU_collab_tag_proposal_01.emailsafe.html` | 4.9 M | 5,089,612 | **TRACKED** | yes (.md + .standalone) | HKU handoff email-embedded (base64 figs) deck; 5/24 delivered | CURRENT | LIVE (handoff record) | **KEEP** | 0 | tracked + handoff deliverable (HKU). do not delete |
| `experiments/in_progress/2026/06/20260601_ISM_aux_tag_figure_review_interactive_01.standalone.html` | 1.9 M | 1,948,887 | UNTRACKED | **yes** | interactive ISM aux-tag figure review (today) | in_progress | LIVE | **NEEDS_USER_DECISION** | 1.9 M | regenerable from .md but brand-new (today), likely in active use |
| `reports/in_progress/2026/05/20260527_HCC1395_TSG_promoter_ASM_reviewer_response_01.standalone.html` | 1.78 M | 1,781,913 | UNTRACKED | **yes** | TSG promoter ASM reviewer response (LIVE ASM ⭐3) | in_progress | LIVE | **NEEDS_USER_DECISION** | 1.78 M | .md companion exists → regenerable; LIVE characterization line |
| `reports/in_progress/2026/05/20260528_longread_ASM_method_workflow.standalone.html` | 1.13 M | 1,125,694 | UNTRACKED | **no** | longread ASM method workflow (PI explainer) | in_progress | LIVE | **NEEDS_USER_DECISION** | — | no .md → this HTML is the primary artifact, not regenerable; keep until .md exists |
| `reports/pi_reports/2026/05/20260531_ASM_H2H3H5_interactive.standalone.html` | 1.1 M | 1,100,219 | UNTRACKED | **no** | PI interactive ASM H2/H3/H5 (5/31) | CURRENT (pi) | LIVE | **NEEDS_USER_DECISION** | — | PI deliverable, no .md source → not regenerable; conservative KEEP-leaning |
| `reports/pi_reports/2026/05/20260531_ASM_mixture_decomposition.standalone.html` | 771 K | 771,338 | UNTRACKED | (check) | PI ASM mixture decomposition (5/31) | CURRENT (pi) | LIVE | **NEEDS_USER_DECISION** | — | PI deliverable; verify .md before any delete |
| `reports/pi_reports/2026/05/20260524_ClairS_TO_HKU_completion_dashboard_02.standalone.html` | 758 K | 758,257 | UNTRACKED | (check) | HKU completion dashboard (5/24 handoff) | CURRENT (pi) | LIVE (handoff) | **KEEP** | 0 | HKU handoff record; conservative KEEP |
| `experiments/in_progress/2026/05/20260531_ISM_aux_tag_observation_funnel_01.standalone.html` | 644 K | 643,892 | UNTRACKED | **yes** | ISM aux-tag observation funnel | in_progress | LIVE | **NEEDS_USER_DECISION** | 0.64 M | regenerable from .md |
| `experiments/in_progress/2026/05/20260531_ISM_aux_tag_critique_and_improvement_01.standalone.html` | 421 K | 420,960 | UNTRACKED | (check) | ISM aux-tag critique/improvement | in_progress | LIVE | **NEEDS_USER_DECISION** | 0.42 M | regenerable if .md present |
| `experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.standalone.html` | 99 K | 99,127 | UNTRACKED | yes (in git status ??) | V6 TPFP characterization (LIVE V6 production) | PRE-FIX | LIVE | **NEEDS_USER_DECISION** | 0.1 M | regenerable; LIVE V6 line |
| `reports/in_progress/2026/05/20260529_project_goal_landscape_dashboard_01.standalone.html` | 181 K | 180,626 | **TRACKED** | yes | goal-landscape dashboard (committed d300ba0/859c55a) | CURRENT | LIVE | **KEEP** | 0 | tracked deliverable |
| `architecture/20260519_研究路徑與五大目標路線圖_01.standalone.html` | (sm) | — | **TRACKED** | yes | research roadmap | CURRENT | LIVE | **KEEP** | 0 | tracked |
| `reports/validated/2026/05/*.standalone.html` (errata/HP-bug/harness/4Layer) | 44-100 K each | — | **TRACKED** | yes | validated 終版 PI/engineering reports | CURRENT (validated) | CONCLUDED_POSITIVE | **KEEP** | 0 | validated終版 = KEEP per rule |

> **Net reclaimable from untracked-with-.md regenerable HTML** (only if user opts to regenerate-on-demand): ISM_aux_tag_figure_review 1.9 M + TSG_promoter_ASM 1.78 M + observation_funnel 0.64 M + critique 0.42 M + V6_TPFP 0.1 M ≈ **4.84 MB**. All NEEDS_USER_DECISION (LIVE lines, conservative).

---

## 3. (c) docs → 即將移除路徑 對照表 + broken pointer 清單

### 3.1 docs → 即將被清理的 output/research 路徑（doc 側衝擊）

| to-be-cleaned path token | referencing docs (scoped grep docs/) | doc-side impact |
|---|---|---|
| `bip8_output_archive` / `big8_output_archive` | `data_specs/20260411_資料盤點快照_01.md`, `data_specs/20260414_output資料信任度與生命週期_01.md`, `standards/20260314_big7_canonical輸出與延續研究規範_01.md`, `archive/OLD_DATA_CORRECTION_NOTICE.md`, `experiments/.../20260510_external_dirs_inventory_proposal_01.md`, several drafts/backup | **Documentation OF the archive** (canonical path `big7_disk_output/{bip8,big8}_output_archive/`, SUPERSEDED). These docs DESCRIBE the data being cleaned → if external archive removed, update these specs (do not delete the specs — they are the provenance record). Impact = LOW (descriptive, not active pointer). |
| `research_rounds/archive` | `data_specs/20260414_...`, `experiments/.../20260423_KDE_Smoke_Test...`, `experiments/.../20260422_Archive_TO_Rerun_ISM_Requirement_01.md` | references to archived research_rounds output; `research/research_rounds/archive` path MISSING under InterSubMod/ (lives in big7_disk_output) → see broken pointers §3.2 |
| `methyl_augmented_filter` | **LIVE**: `CURRENT_FOCUS.md`, `experiments/INDEX.md`, `concepts/20260524_...`, `concepts/20260529_project_task_registry...`, multiple A2/A3/A4 experiment .md, V6_signoff presentations, `reports/validated/.../master_draft.md` | research output is DEAD (⭐2 L4 LOSO NEGATIVE) BUT docs are **G6 paper methods evidence**. `research/methyl_augmented_filter_phase2/` EXISTS. **Docs KEEP; research data ARCHIVE (not delete)** per最高原則. |
| `to_pipeline_staging` | `trash/to_pipeline_staging_v1/README.md`, `trash/.../reports/20260412_...md` | only the trash copy references it; v2 corrected lives in `research/to_pipeline_staging/` (EXISTS, verified). Trash v1 → SAFE_DELETE. |
| `batch_runs` | `data_specs/20260411_資料盤點快照_01.md`, `experiments/.../20260510_external_dirs_inventory_proposal_01.md` | inventory/snapshot descriptions only; LOW impact |
| `complete_matrix` | `plans/2026/04/20260426_S5_Erratum...`, `data_specs/20260411_工作區命名...`, several methodology/validated reports, `experiments/.../A0_assets/DATA_PROVENANCE_LEDGER.md` | mostly descriptive references to a matrix-build output; verify A0_assets ledger before touching that output |

### 3.2 Broken pointers (referenced data path 已不存在)

Probed from docs references (scoped):

| referenced path | exists? | referencing docs | status |
|---|---|---|---|
| `/big8_disk/output_archive` (i.e. big8-disk-side originals) | **MISSING** | data_specs/standards describe migration big8→big7 | EXPECTED — migration moved these to `big7_disk_output/{big8,bip8}_output_archive/`; docs already record the migration. NOT a true break (historical note). |
| `research/research_rounds/archive` (under InterSubMod/) | **MISSING** | `data_specs/20260414`, KDE_Smoke_Test, Archive_TO_Rerun_ISM_Requirement | **BROKEN within InterSubMod scope** — the archive lives in `big7_disk_output/`, not `research/`. Docs may imply a local path that does not exist. LOW impact (archive/rerun-requirement docs are historical). |
| `research/to_pipeline_staging` (v2) | **EXISTS** ✓ | trash v1 README points here as "正確版本" | OK — confirms trash v1 SAFE_DELETE is safe |
| `research/methyl_augmented_filter_phase2` | **EXISTS** ✓ | INDEX/CURRENT_FOCUS/A2-A4 | OK — not broken; DEAD-but-archive |

> No catastrophic broken pointers inside `docs/`. The `big8`/`research_rounds/archive` "misses" are migration/historical artifacts already documented; treat as advisory, not as blocking the cleanup.

---

## 4. Per-artifact verdict table (top reclaimable, full fields)

| path | size_human | bytes | purpose | trust_tier | conclusion_status | cleanup_verdict | reclaimable_bytes | referenced_by | rationale |
|---|---|---|---|---|---|---|---|---|---|
| `docs/trash/to_pipeline_staging_v1/` | 8.0 M | 8,273,343 | v1 multi-stage TO characterization built on WRONG VCF/ISM (F1=0.649 not 0.7127); moved to trash 2026-04-14 | DEPRECATED | DEAD | **SAFE_DELETE** | 8,273,343 | self (README) only; v2 corrected = research/to_pipeline_staging (EXISTS) | README explicitly "待刪除"; wrong data; corrected version verified present. Transient/superseded → safe to delete. |
| `docs/trash/.../data/hcc1395_pass_multimodal.csv` | 7.8 M | 8,175,945 | bulk of trash — merged table on wrong VCF | DEPRECATED | DEAD | **SAFE_DELETE** | 8,175,945 | none | regenerable-or-irrelevant (wrong source); 99% of trash size |
| `docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2 (+base dup figs)` | 16 M (v2) | 34,446,488 (whole) | 4/29 professor report; base+v2+v3 iterations w/ duplicated figures | CURRENT (validated) | LIVE | **NEEDS_USER_DECISION** | ~7-11 M (cross-version dup figures only; never the latest v2) | INDEX, weekly reports | validated deliverable — do NOT auto-delete; only dedup figures after user confirms which version is canonical |
| untracked `*.standalone.html` w/ .md companion (5 files) | ~4.8 M | 4,839,778 | in_progress preview companions, regenerable via /html-report-build | in_progress | LIVE | **NEEDS_USER_DECISION** | 4,839,778 | their own .md | regenerable from .md but LIVE research lines (ASM ⭐3, V6, ISM aux-tag) → conservative |
| `docs/**/__pycache__` (2 dirs) | 117 K | 119,936 | python bytecode cache from report-build scripts | TRANSIENT | DEAD | **SAFE_DELETE** | 119,936 | none | pure regenerable bytecode |
| PI-report standalone w/o .md (ASM_H2H3H5, ASM_mixture, longread_ASM) | ~3.0 M | 2,997,251 | PI deliverables, HTML is primary (no .md source) | CURRENT (pi) | LIVE | **NEEDS_USER_DECISION** | — | (PI consumption) | NOT regenerable (no .md) → must not delete; KEEP-leaning |
| `docs/drafts/` + `drafts/backup/` | 100 K | ~102,000 | CLAUDE/AGENTS md v2/v3 drafts + pre-v3deploy backups | SUPERSEDED | DEAD | **NEEDS_USER_DECISION** | ~102,000 | none (superseded by live CLAUDE.md) | superseded by deployed CLAUDE.md/AGENTS.md but tiny + may be wanted as governance history |
| `docs/archive/` | 3.5 M | 3,094,844 |封存區 (2025 + 2026 + deep old_structure) | SUPERSEDED (archived) | DEAD-but-archived | **KEEP** | 0 | various | it IS the archive; keep |

---

## 5. Totals (for parent aggregation)

- total_scanned_bytes (docs/) = **111,558,565**
- KEEP (validated/tracked deliverables, archive, references, active prep): bulk — see below
- SAFE_DELETE (trash v1 + __pycache__) = 8,273,343 + 119,936 = **8,393,279**
- ARCHIVE (none net-new in docs; archive/ already archived → counted in KEEP) = **0** *(note: the methyl-filter RESEARCH data — outside docs scope, area 1-4 — is ARCHIVE; docs themselves are KEEP)*
- NEEDS_USER_DECISION (cross-version dup figs in 教授報告 ~9 M + untracked regenerable standalone HTML ~4.84 M + PI standalone w/o .md ~3.0 M + drafts ~0.1 M) ≈ **~17 MB upper bound** (conservative; much is LIVE)

---

## 6. Anomalies vs prior先驗

1. **`methyl_augmented_filter` is LIVE in docs despite DEAD research verdict** — docs (CURRENT_FOCUS/INDEX/A2-A4/master_draft) cite it as G6 paper methods (LOSO NEGATIVE evidence). Confirms先驗: ARCHIVE research data, KEEP docs. Do NOT delete the experiment .md files.
2. **`research_rounds/archive` path MISSING under InterSubMod/research/** but referenced by 3 docs — the archive actually lives in `big7_disk_output/`; docs imply a local path that doesn't exist (minor broken pointer, historical docs, LOW impact).
3. **PI-report standalone HTMLs (1MB+) have no .md companion** (ASM_H2H3H5, ASM_mixture, longread_ASM) → they are NOT regenerable; the HTML IS the deliverable. Cannot treat as "regenerable preview". KEEP-leaning NEEDS_DECISION.
4. **教授報告 cross-version figure duplication** — fig01c/fig14/fig15/fig16 byte-identical in base/figures and v2/figures (0.5 MB each ×2). Dedup possible but it's a validated deliverable → user decides canonical version first.
5. **`/big8_disk/output_archive` MISSING** is EXPECTED (big8→big7 migration documented), not a real break.
