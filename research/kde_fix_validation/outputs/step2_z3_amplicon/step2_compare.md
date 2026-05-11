# Step 2: Z3 amplicon blacklist — stale vs KDE-fixed ΔF1

## Method

- Stale CovM = NumReads / 75.0 (from master `Coverage_Multiple`)
- Fixed CovM = NumReads / per-sample KDE baseline
- Z3 flag (TO mode, loh_flag + af_extreme + HPFineNGroups<=1) unchanged
- S1 (literature arms) + S2 (whole-chr) coordinate-based → CovM-independent (sanity ∼0)
- S3 (CovM > 95%ile non-Z3) → per-sample cutoff changes with baseline

## Per-sample ΔF1

| Sample | Baseline new | F1 base | S1 Δ | S2 Δ | S3 stale Δ | S3 fixed Δ | S3 change |
|--------|-------------:|--------:|-----:|-----:|-----------:|-----------:|----------:|
| H2009 | 79x | 0.9545 | -0.0024 | -0.0025 | +0.0000 | +0.0000 | +0.0000 |
| H1437 | 69x | 0.8712 | -0.0072 | -0.0086 | -0.0002 | -0.0002 | +0.0000 |
| HCC1395 | 55x | 0.8309 | -0.0037 | -0.0041 | +0.0001 | +0.0001 | +0.0000 |
| HCC1395_DORADO | 53x | 0.8330 | -0.0041 | -0.0039 | +0.0003 | +0.0003 | +0.0000 |
| HCC1937 | 91x | 0.6772 | +0.0023 | +0.0038 | +0.0000 | +0.0000 | +0.0000 |
| HCC1954 | 61x | 0.4047 | +0.0058 | +0.0065 | +0.0002 | +0.0002 | +0.0000 |
| COLO829 | 29x | 0.7906 | -0.0108 | -0.0109 | -0.0004 | -0.0004 | +0.0000 |

## 95th percentile cutoffs

| Sample | Stale cutoff | Fixed cutoff | ratio stale/fixed |
|--------|-------------:|-------------:|------------------:|
| H2009 | 2.480 | 2.354 | 1.053 |
| H1437 | 1.587 | 1.725 | 0.920 |
| HCC1395 | 1.747 | 2.382 | 0.733 |
| HCC1395_DORADO | 1.853 | 2.623 | 0.707 |
| HCC1937 | 2.880 | 2.374 | 1.213 |
| HCC1954 | 1.640 | 2.016 | 0.813 |
| COLO829 | 0.640 | 1.655 | 0.387 |
