# {Group_ID} · {Feature_or_Group_Name} — Feature Layered Observation

**Date**: YYYY-MM-DD
**Group**: {G_ID} ({Group human name})
**Features**: f1, f2, f3, ...
**Data**: `research/feature_layered_observation/data/{master}.tsv.gz` (N rows × M cols, 7 samples × 2 modes)
**Script**: `research/feature_layered_observation/scripts/observe_{group_id}_{slug}.py`

---

## 1. Feature definitions and source

| Feature | Type | Source | Semantic |
|---------|------|--------|----------|
| {f1} | continuous | `src/core/X.cpp:NNN` / include-header | {one-line meaning} |
| {f2} | ordinal | `func_name()` | {one-line meaning} |

**Known bugs/caveats** (prior memory):
- {reference to memory id or landscape doc}
- {reference to pending infra bug if any}

## 2. Observation goals

{1-2 段 hypothesis + prior expectation + what we're trying to rule in/out}

## 3. Global distribution (Step 1)

`fig01_global_distribution.png`, table `data/{group}_global_stats.tsv`:

| feature | AUC [95% CI] | Cohen's d | mean TP | mean FP | MW-U p |
|---------|--------------|-----------|---------|---------|--------|
| {f1} | 0.XXX [0.XXX, 0.XXX] | +/-0.XX | X.X | X.X | 0.XX |

**TP/FP counts**: N_tp TP vs N_fp FP (overall TP rate X%).

Key observations:
- {bullet on where signal comes from}
- {bullet on assay-level vs region-level}
- {bullet on anti-predictive if relevant}

## 4. LOH × AF × CN 32-cell heatmap (Step 2)

Per-feature heatmaps (`{feat}_fig02_heatmap.png`), 4 panels: TP rate / feat_TP / feat_FP / Δ(TP−FP).

Stratification:
- LOH_Subtype (5): LOH_None, LOH_Weak, LOH_Noise, LOH_Strong, LOH_Subclone
- AF_class (3 from vcf_AF): Extreme, Near-half, Intermediate
- cn_tier_F (5 from Coverage_Multiple): CN_Loss, CN_Near1, CN_Diploid, CN_Gain, CN_HighGain

Highlights from `data/{group}_cell_delta.tsv`:

- {cell where |Δ| largest; n=X; Cohen's d=X}
- {cell with reversed sign if any}
- {coverage-related pattern}

## 5. Cross-sample consistency (Step 3)

`fig03_per_sample_consistency.png`. Spearman concordance matrix 7×7 on per-cell TP rate.

| Sample | Global AUC | Cohen's d | Direction | Notable cell |
|--------|-----------|-----------|-----------|--------------|
| HCC1395 | 0.XXX | +/- | ↑/↓ | — |
| HCC1395_DORADO | 0.XXX | | | |
| HCC1937 | 0.XXX | | | |
| HCC1954 | 0.XXX | | | |
| H2009 | 0.XXX | | | |
| H1437 | 0.XXX | | | |
| COLO829 | 0.XXX | | | |

**Concordance**: median ρ = X.XX; same-direction {n}/7 samples.

## 6. Stratified AUC (Step 4)

`fig04_stratified_auc.png`, table `data/{group}_auc_table.tsv`.

- Global AUC: 0.XXX
- Per-LOH range: [0.XX, 0.XX]
- Per-AF range: [0.XX, 0.XX] (Extreme=X, Near-half=X, Intermediate=X)
- Per-CN range: [0.XX, 0.XX]

{1-2 sentences on which stratum carries the signal}

## 7. Confound check (Step 5)

`fig05_confound_residualized.png`, table `data/{group}_confound.tsv`.

**Gate 1 (within-group OLS on NumReads + vcf_AF + Coverage_Multiple)**:
- Raw AUC: 0.XXX
- Residualized AUC: 0.XXX
- Δ = {raw - resid}

**Gate 2 (AF-bin stratification)**:
| AF bin | AUC | n |
|--------|-----|---|
| Extreme | 0.XXX | N |
| Near-half | 0.XXX | N |
| Intermediate | 0.XXX | N |

Range: X.XX (pass if <0.10)

**Gate 3 (permutation test, 1000 reps)**:
- Observed AUC: 0.XXX
- Null 95th percentile: 0.XXX
- p-value: 0.XXX (pass if <0.05)

**Verdict for confound guard**: PASS / FAIL (reason: ...)

## 8. Spatial autocorrelation (Step 6)

`fig06_spatial_autocorrelation.png`. Per 5Mb bin AUC genome-wide.

- N bins with n≥50: X
- Top-5 high-AUC bins: {chrX:Y_Z, ...}
- High-TP-rate overlap: {Yes/No} → {artifact warning / clean}

## 9. Papers and knowledge base context

1. {Citation 1}. {Relevance sentence.} {DOI or knowledge path}
2. {Citation 2}. ...
3. {Citation 3}. ...

Knowledge MCP:
- {knowledge/path/doc_id.md} — {relevance}

## 10. Conclusion and grilling

**Verdict**: POSITIVE / CONDITIONAL_POSITIVE / CONFOUND_COLLAPSED / NEGATIVE / SAMPLE_SPECIFIC / ANTI_SIGNAL / SPATIAL_ARTIFACT

**Logic chain**:
- Raw AUC {A} → {pass/fail 0.58 gate}
- Residualized AUC {B} → {pass/fail 0.55}
- Cross-sample {n}/7 → {pass/fail 5/7}
- Spatial → {pass/fail artifact check}
- ⟹ Verdict = {V}

**Three grilling questions**:

1. **Confound residual check**: 這個信號在 NumReads + vcf_AF + Coverage_Multiple 三個殘差後還剩 {X}? 答：{...}
2. **Cross-sample robustness**: 去掉 HCC1395 (TP 92%) 之後，AUC 還有 {X}? 答：{...}
3. **Pipeline dependency**: paired 與 TO 同向嗎？若只在 paired 有效，是否為 HP-tag 人工產物？答：{...}

**Follow-up suggestions**:
- {next step if POSITIVE: F1 filter pilot / per-sample deep}
- {next step if CONFOUND_COLLAPSED: what additional confound to test}
- {next step if SAMPLE_SPECIFIC: stratified modeling}

**Tier**: {1-5}. **evidence_ledger.jsonl entry**: {hypothesis_id, cycle_id}.
