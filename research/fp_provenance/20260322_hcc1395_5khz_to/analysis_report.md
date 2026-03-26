# TO FP provenance analysis workspace

## Sample summary

| sample | to_raw_fp_count | caller_pon_filtered | caller_nonpon_filtered | longphase_to_removed | to_postprocess_removed | to_residual_final_fp | paired_resolved_count | paired_persistent_count | to_final_f1 | paired_final_f1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HCC1395 | 2731541 | 2717339 | 2596 | 0 | 8 | 11598 | 2731016 | 525 | 0.716700 | 0.853200 |

## Discovery top rules

| sample | rule_id | rule_label | fp_removed | tp_removed | f1_after | delta_f1_vs_final |
| --- | --- | --- | --- | --- | --- | --- |
| HCC1395 | af_le_0.03_ad_ge_0.15 | AF<=0.03 and AlleleDelta>=0.15 | 0 | 0 | 0.716656 | -0.000044 |
| HCC1395 | af_le_0.03_ad_ge_0.20 | AF<=0.03 and AlleleDelta>=0.20 | 0 | 0 | 0.716656 | -0.000044 |
| HCC1395 | af_le_0.03_ad_ge_0.25 | AF<=0.03 and AlleleDelta>=0.25 | 0 | 0 | 0.716656 | -0.000044 |
| HCC1395 | af_le_0.03_ad_ge_0.30 | AF<=0.03 and AlleleDelta>=0.30 | 0 | 0 | 0.716656 | -0.000044 |
| HCC1395 | af_le_0.05_ad_ge_0.15 | AF<=0.05 and AlleleDelta>=0.15 | 0 | 0 | 0.716656 | -0.000044 |
| HCC1395 | af_le_0.05_ad_ge_0.20 | AF<=0.05 and AlleleDelta>=0.20 | 0 | 0 | 0.716656 | -0.000044 |
| HCC1395 | af_le_0.05_ad_ge_0.25 | AF<=0.05 and AlleleDelta>=0.25 | 0 | 0 | 0.716656 | -0.000044 |
| HCC1395 | af_le_0.05_ad_ge_0.30 | AF<=0.05 and AlleleDelta>=0.30 | 0 | 0 | 0.716656 | -0.000044 |
| HCC1395 | af_le_0.08_ad_ge_0.15 | AF<=0.08 and AlleleDelta>=0.15 | 0 | 0 | 0.716656 | -0.000044 |
| HCC1395 | af_le_0.08_ad_ge_0.20 | AF<=0.08 and AlleleDelta>=0.20 | 0 | 0 | 0.716656 | -0.000044 |
