# Source Extraction — longphase_to_v5_somatic_04-29

## Fixture summary

- **expected dual intercept**:
  - **P2 PRECHECK**: warning on working tree dirty (binary_version contains "uncommitted" string → custom probe should flag)
  - **P5 EVALUATE**: warn/downgrade on multi_sample_consistency=1/7=0.143 (critical) + effect_size_stability low

## Source mapping

| field | value | source |
|---|---|---|
| binary_version | `380e8d2_uncommitted_layer1_5_alt_guard` | trigger §5.2 (working tree = 380e8d2 + Layer 1.5 + alt guard 兩塊未 commit) |
| dataset_id | `HCC1395_5kHz_pilot_only` | trigger pilot scope §4 |
| pilot.value | 0.083 (+8.3pp full-genome clean PS) | trigger §6.4 row |
| pilot.confound_guard.notes | SEQC2 F1 Δ=-0.0003 (noise) | trigger §3.4 PI report 4 §3.4 |
| generalize.samples × 7 | 1/7 validated; rest 0 (not run) | trigger §pilot scope; 7 sample expansion not done |

## Why dual intercept (P2 + P5)

- **P2 catch (warning level)**: binary_version 字串含 "uncommitted" — current `_staleness_check.py` 已加 working_tree_clean check (Day 4 增強)。因為 plan 標的 SHA 並非真實 git SHA (含 _uncommitted 後綴)，git_sha_exists() 會回 false → status="missing"。這也算是 BLOCKED 一種形式
- **P5 catch (主要)**: multi_sample_consistency=1/7 critical → forced pending_review

## Expected harness behavior

- **/check-staleness**: BLOCKED (binary SHA not in git history because "_uncommitted_" is not a valid SHA — the heuristic correctly flags reproducibility issue)
- **/run-evaluator** (if forced past P2):
  - multi_sample_consistency = 1/7 ≈ 0.143 (**critical** <0.2)
  - effect_size_stability = min(0.083, 0)/max ≈ 0 / 0.083 = 0 (**critical** <0.2)
  - 2 critical components → forced pending_review

## Fidelity note

- **fidelity: high** for HCC1395 single-sample data (Agent 2 §6.4 表)
- **fidelity: low** for 6/7 missing samples — they're filled with `metric_value=0.0` to represent "not run" (deliberately drives multi_sample_consistency low to stress evaluator)
- **caveat**: plan §10 #1「state.json schema 不可頻繁改版」— 這 fixture 的 binary_version 用人造字串而非真 SHA，是 deliberately 觸發 missing 條件的 P2 偵測；不違反 schema（schema 允許任意 string）但偏離正常 SHA 慣例
