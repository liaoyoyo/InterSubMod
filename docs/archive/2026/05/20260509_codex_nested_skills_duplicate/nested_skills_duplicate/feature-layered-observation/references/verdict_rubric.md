# Verdict Rubric — Feature Layered Observation

決策樹形式的 verdict 規則。每個 feature 觀察完成後，依下表得到單一 verdict。

## 主決策樹

```
Raw AUC (Step 1)
│
├─ < 0.50 ────────────────────────────────────────────────► ANTI_SIGNAL
│   （direction flipped；可觀察 1 - raw_auc 以查 FP > TP 方向）
│
├─ [0.50, 0.53) ──────────────────────────────────────────► NEGATIVE
│   （random-level；Step 2-6 可簡化，只做 fig02_tp_rate 確認無局部強）
│
├─ [0.53, 0.58) ──────────────────────────────────────────► NEGATIVE
│   （Beyond-AUC ceiling 以下；Step 2-6 跑完存檔但不升 tier）
│
└─ ≥ 0.58 ─► Confound Guard (Step 5)
              │
              ├─ Resid AUC < 0.55 ─────────────────────────► CONFOUND_COLLAPSED
              │   （NumReads + vcf_AF + Coverage_Multiple 解釋大部分 signal）
              │
              └─ Resid AUC ≥ 0.55 ─► AF-bin Range (Step 5 Gate 2)
                                      │
                                      ├─ range ≥ 0.10 ────► CONFOUND_COLLAPSED (AF-driven)
                                      │
                                      └─ range < 0.10 ────► Permutation p (Step 5 Gate 3)
                                                             │
                                                             ├─ p ≥ 0.05 ──► NEGATIVE (chance)
                                                             │
                                                             └─ p < 0.05 ──► Cross-sample (Step 3)
                                                                             │
                                                                             ├─ <5/7 同向 ──► SAMPLE_SPECIFIC
                                                                             │
                                                                             └─ ≥5/7 同向 ──► Spatial (Step 6)
                                                                                             │
                                                                                             ├─ artifact only ──► SPATIAL_ARTIFACT
                                                                                             │
                                                                                             └─ clean ──► Raw AUC magnitude
                                                                                                         │
                                                                                                         ├─ [0.58, 0.65) ──► CONDITIONAL_POSITIVE
                                                                                                         │
                                                                                                         └─ ≥ 0.65 ────────► POSITIVE
```

## Verdict 定義與 tier 對應

| Verdict | 定義 | Tier | evidence_ledger `tier_flags` |
|---------|------|:----:|------------------------------|
| **POSITIVE** | raw≥0.65 + resid≥0.55 + 5/7 跨樣本 + 非 spatial artifact | 3-4 | `within_group_ols, af_bin_stratified, permutation_tested, cross_sample_verified, spatial_clean` |
| **CONDITIONAL_POSITIVE** | raw [0.58, 0.65) + 其他 gate 全過 | 3 | 同上，但 `conditional_low_effect_size` |
| **CONFOUND_COLLAPSED** | raw≥0.58 但殘差 <0.55 或 AF-bin range ≥0.10 | 1-2 | `confound_suspected` |
| **SAMPLE_SPECIFIC** | 所有 gate 過但跨樣本 <5/7 | 2 | `sample_specific` |
| **NEGATIVE** | raw [0.50, 0.58) 或 permutation p≥0.05 | 1 | `not_discriminative` |
| **ANTI_SIGNAL** | raw <0.50；FP > TP 方向 | 1 | `direction_flipped`；可標註為 FP fingerprint |
| **SPATIAL_ARTIFACT** | 所有 gate 過但 AUC 只在 high-TP-rate bin | 1 | `spatial_artifact` |
| **CHARACTERIZATION_ONLY** | raw≥0.58 resid≥0.55 但 recall 太低（precision-oriented）| 2 | `annotation_only, recall_limited` |
| **DATA_GAP** | 樣本 coverage 不足或特徵大量 NaN | — | `data_gap`；不登 tier |

## 特殊情況

### 1. Assay-level shortcut detection

若特徵在同 sample×mode 的所有 rows 有相同值（變異 ≤1e-6），例 `Diploid_Coverage_Used`：

- 標註：`assay_level_shortcut`
- Verdict override：即使 raw AUC ≥ 0.58 也降為 NEGATIVE（"which sample are we in" shortcut）
- 參考：G1 `Diploid_Coverage_Used` 0.790 案例

### 2. Precision-recall tradeoff

對於 `same_hap_marker` 型 fingerprint（precision 極高但 recall 低）：

- 不計 AUC；直接標 CHARACTERIZATION_ONLY
- 報告 median precision + precision@recall=0.5

### 3. Pipeline-dependent

若 paired AUC 與 TO AUC 方向相反：

- 報告兩個 verdict（per-mode）
- 若 paired-only positive → 標註 `hp_tag_dependent`，警告可能為 self-phased artifact
- 參考：memory `project_readparser_germline_hp_only_phase1_negative`

### 4. LOH-gated verdict

若特徵僅在特定 LOH_Subtype 有效（per-LOH AUC range > 0.15）：

- 若只有 `LOH_Subclone` 或 `LOH_Strong` cells 有效 → `LOH_GATED_POSITIVE`
- 若只有 `LOH_None` cells 有效 → 可能是 non-LOH fingerprint，要小心

## Tier ceiling 規則

- 單 feature 最高 tier = 4（需 ≥3 pipeline 獨立驗證）
- Tier 5 = breakthrough + external literature replication；單 feature 觀察不可達此

## 寫入 evidence ledger 範本

```json
{
  "ts": "YYYY-MM-DDThh:mm:ss",
  "hypothesis_id": "FLO-Gn-{feature}",
  "cycle_id": "cycle_YYYYMMDD_flo_{feature}",
  "action": "observation",
  "verdict": "{VERDICT}",
  "tier": T,
  "tier_flags": ["within_group_ols", "af_bin_stratified", "permutation_tested", "cross_sample_verified", "spatial_clean"],
  "confidence_cap": C,
  "raw_auc": 0.XX,
  "resid_auc": 0.XX,
  "af_bin_range": 0.XX,
  "perm_p": 0.XXX,
  "cross_sample_n_same_direction": X,
  "cross_sample_spearman_median": 0.XX,
  "spatial_artifact": false,
  "confounds_controlled": ["NumReads", "vcf_AF", "Coverage_Multiple"],
  "report_path": "research/feature_layered_observation/features/{feature}.md"
}
```
