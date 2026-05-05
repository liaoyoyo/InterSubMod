# Step 1-6 SOP — Feature Layered Observation 詳細

## 共用設定

```python
SAMPLE_ORDER = ["HCC1395", "HCC1395_DORADO", "HCC1937", "HCC1954",
                "H2009", "H1437", "COLO829"]
LOH_ORDER = ["LOH_None", "LOH_Weak", "LOH_Noise", "LOH_Strong", "LOH_Subclone"]
AF_ORDER  = ["Extreme", "Near-half", "Intermediate"]
CN_BREAKS = [-np.inf, 0.65, 0.99, 1.33, 1.82, np.inf]
CN_LABELS = ["CN_Loss", "CN_Near1", "CN_Diploid", "CN_Gain", "CN_HighGain"]

TRUTH_PALETTE = {"TP": "#1f77b4", "FP": "#ff7f0e"}  # 藍/橘固定
MODE_PALETTE  = {"paired": "#2ca02c", "to": "#d62728"}  # 綠/紅
```

## Step 1 · 全域分佈

**目的**：判斷特徵在全樣本 pooled 下是否有 base-level 區分力。

**做法**：
1. 合併 7 samples × 2 modes → single pooled DataFrame
2. Groupby tp_label ∈ {0, 1}
3. 計算：
   - `mean_tp, mean_fp, std_tp, std_fp, median_tp, median_fp`
   - Cohen's d: `(mean_tp - mean_fp) / pooled_std`
   - `scipy.stats.mannwhitneyu(tp, fp, alternative='two-sided')`
   - `sklearn.metrics.roc_auc_score(label, feature)`
   - Wilson-style (Hanley-McNeil) 95% CI on AUC
4. 繪製 violin plot（tp_label=0/1）+ KDE overlay + n label + AUC annotation

**輸出**：
- 圖：`{group}/fig01_global_distribution.png`
- 表：`data/{group}_global_stats.tsv` (feature, auc, lo, hi, cohen_d, mwu_p, mean_tp, mean_fp, n_tp, n_fp)

**判定分流**：
- AUC ≥ 0.58 → 進 Step 5 confound guard
- 0.53 ≤ AUC < 0.58 → 跑完 Step 2-6 但 verdict 預設 NEGATIVE
- AUC < 0.53 → NEGATIVE，Step 2 做簡化 heatmap 後跳到結論

## Step 2 · LOH × AF × CN 32-cell heatmap

**目的**：確認信號在哪個生物學區塊最強；捕捉 stratum-level 效應。

**做法**：
1. 切分：`LOH_Subtype (5) × AF_class (3, 用 vcf_AF) × cn_tier_F (5, CovM)` = 75 cells
2. 丟掉 n < 20 的 cells；預期剩 ~32 有效 cells
3. 4 張 heatmap：
   - **A**: per-cell TP rate（色階 0-1，YlOrRd）
   - **B**: feature mean / median for TP only
   - **C**: feature mean / median for FP only
   - **D**: Δ(TP − FP) + Cohen's d（TwoSlopeNorm at 0, RdBu_r）
4. 每 cell 標註 n（n≥20 顯示整數，n<20 顯示 `·`）
5. row = LOH × AF（15 rows），col = cn_tier_F（5 cols），或轉置

**輸出**：
- 圖：`{group}/fig02_{tp_rate, feat_tp, feat_fp, delta}.png`
- 表：`data/{group}_cell_delta.tsv` (stratum_key, n_tp, n_fp, mean_tp, mean_fp, delta, cohen_d, tp_rate)

## Step 3 · 跨樣本一致性

**目的**：排除 sample-specific artifact；確認是 generic biological signal。

**做法**：
1. 7 samples × (LOH×AF 合併 cn_tier_F) 的 grid heatmap
2. 固定樣本順序
3. Per-sample 計算 AUC、Cohen's d、direction
4. Spearman concordance matrix 7×7 on per-cell TP rate 或 per-cell feature mean
5. Median ρ + 同向樣本數

**輸出**：
- 圖：`{group}/fig03_per_sample_consistency.png`
- 表：`data/{group}_per_sample_auc.tsv`

**門檻**：
- median ρ ≥ 0.5 → 一致
- 同向樣本 ≥ 5/7 → cross-sample pass
- 否則 → SAMPLE_SPECIFIC verdict

**呼叫 skill**：`/multi-sample-consistency`（結論依賴跨樣本時必觸發）

## Step 4 · 分層 AUC

**目的**：定位信號集中在哪個 stratum。

**做法**：
1. Compute Global AUC、per-LOH AUC（5 bars）、per-AF AUC（3 bars）、per-CN AUC（5 bars）
2. Bar chart 加虛線 0.50 (random)、0.58 (Beyond-AUC ceiling)
3. n < 50 的 stratum 標 `insufficient_n`

**輸出**：
- 圖：`{group}/fig04_stratified_auc.png`
- 表：`data/{group}_auc_table.tsv`

## Step 5 · Confound guard（呼叫 /auc-confound-guard）

**前置**：raw AUC ≥ 0.58 才跑；否則寫 skipped_reason=low_raw_auc。

**Gate 1 · Within-group OLS residualization**:

```python
# CORRECT: within-group OLS
from sklearn.linear_model import LinearRegression
confounds = ['NumReads', 'vcf_AF', 'Coverage_Multiple']  # G10 加 AlleleDelta
feat_resid = np.zeros(len(df))
for label in [0, 1]:
    mask = (df['tp_label'] == label).values
    X = df.loc[mask, confounds].values
    y = df.loc[mask, feature].values
    valid = ~np.isnan(y) & ~np.isnan(X).any(axis=1)
    lr = LinearRegression().fit(X[valid], y[valid])
    feat_resid_sub = np.full(mask.sum(), np.nan)
    feat_resid_sub[valid] = y[valid] - lr.predict(X[valid])
    feat_resid[mask] = feat_resid_sub
resid_auc = roc_auc_score(df['tp_label'], feat_resid)
```

**禁止**：pooled OLS（會保留 TP/FP 組間差 = 虛假 delta）。見 memory `feedback_pooled_ols_residualization_trap`。

**Gate 2 · AF-bin stratification**:

```python
bins = {'Extreme': df['AF_class']=='Extreme',
        'Near-half': df['AF_class']=='Near-half',
        'Intermediate': df['AF_class']=='Intermediate'}
af_aucs = {}
for name, mask in bins.items():
    if mask.sum() < 50: continue
    af_aucs[name] = roc_auc_score(df.loc[mask, 'tp_label'], df.loc[mask, feature])
# pass: range(af_aucs.values()) < 0.10
```

**Gate 3 · Permutation test**:

```python
observed = roc_auc_score(df['tp_label'], df[feature])
null_aucs = []
for _ in range(1000):
    shuffled = np.random.permutation(df['tp_label'].values)
    null_aucs.append(roc_auc_score(shuffled, df[feature]))
p = (np.array(null_aucs) >= observed).mean()
# pass: p < 0.05
```

**輸出**：
- 圖：`{group}/fig05_confound_residualized.png`（raw vs resid bar + AF-bin trace）
- 表：`data/{group}_confound.tsv` (feature, raw_auc, resid_auc, af_bin_aucs, perm_p)

## Step 6 · Spatial autocorrelation

**目的**：排除 genome-position-driven artifact。

**做法**：
1. 將 Chr + Pos 聚合為 5 Mb bins
2. Per-bin AUC（需 n ≥ 50）
3. 疊加 per-bin baseline TP rate
4. 若高 AUC 只出現在 baseline TP rate > 80% 的 bin → spatial artifact warning

**輸出**：
- 圖：`{group}/fig06_spatial_autocorrelation.png`
- 表：`data/{group}_spatial.tsv` (chr, bin_start, n, tp_rate, auc)

## Evidence ledger registration

完成後必寫入 `research/autoresearch/evidence_ledger.jsonl`：

```json
{"ts":"YYYY-MM-DDThh:mm:ss","hypothesis_id":"FLO-{Gn}-{feature}","cycle_id":"cycle_YYYYMMDD_{feature}","action":"observation","verdict":"{V}","tier":T,"raw_auc":0.XX,"resid_auc":0.XX,"cross_sample_n":X,"confounds_controlled":["NumReads","vcf_AF","Coverage_Multiple"],"report_path":"research/feature_layered_observation/features/{feature}.md"}
```
