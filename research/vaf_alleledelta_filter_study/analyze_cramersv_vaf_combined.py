#!/usr/bin/env python3
"""
CramersV AND VAF 組合過濾分析

分析不同 CramersV 和 VAF 門檻組合對 TP/FP 的影響

Author: Generated for InterSubMod project
Date: 2026-01-23
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Noto Sans CJK TC', 'Noto Sans CJK SC', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

import os
import json

# Data paths
DATA_DIR = "/big8_disk/liaoyoyo2001/InterSubMod/output/ml_feature_exploration/2026/01/data"
OUTPUT_DIR = "/big8_disk/liaoyoyo2001/InterSubMod/analysis/vaf_alleledelta_filter_study"
os.makedirs(f"{OUTPUT_DIR}/cramersv_vaf_analysis", exist_ok=True)

# Load data
print("=" * 70)
print("CramersV AND VAF 組合過濾分析")
print("=" * 70)
print("\n載入資料...")
tp_df = pd.read_csv(f"{DATA_DIR}/tp_with_qual_vaf.csv")
fp_df = pd.read_csv(f"{DATA_DIR}/fp_with_qual_vaf.csv")

print(f"TP 總數: {len(tp_df)}")
print(f"FP 總數: {len(fp_df)}")

SEQC2_TOTAL = 39447

#----------------------------------------------
# 1. CramersV 與 VAF 分布概覽
#----------------------------------------------
print("\n" + "=" * 70)
print("1. CramersV 與 VAF 分布概覽")
print("=" * 70)

# Define thresholds
vaf_thresholds = [0.10, 0.15, 0.20, 0.24, 0.30]
v_thresholds = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50]

print("\n基線統計:")
print(f"  TP VAF: mean={tp_df['VAF'].mean():.4f}, median={tp_df['VAF'].median():.4f}")
print(f"  FP VAF: mean={fp_df['VAF'].mean():.4f}, median={fp_df['VAF'].median():.4f}")
print(f"  TP CramersV: mean={tp_df['CramersV'].mean():.4f}, V>0: {(tp_df['CramersV']>0).sum()}")
print(f"  FP CramersV: mean={fp_df['CramersV'].mean():.4f}, V>0: {(fp_df['CramersV']>0).sum()}")

#----------------------------------------------
# 2. 過濾策略 1: 移除 (低 VAF AND 低 CramersV)
#----------------------------------------------
print("\n" + "=" * 70)
print("2. 過濾策略: 移除 (VAF < threshold AND CramersV < threshold)")
print("目的: 移除低支持度且無統計ASM證據的位點")
print("=" * 70)

results_remove_low = []

print("\n" + "-" * 85)
print(f"{'VAF <':<10} {'V <':<10} {'TP移除':<12} {'TP移除%':<10} {'FP移除':<12} {'FP移除%':<10} {'F1':<8}")
print("-" * 85)

for vaf_th in vaf_thresholds:
    for v_th in v_thresholds:
        # Filter: Remove if VAF < vaf_th AND CramersV < v_th
        tp_mask = (tp_df['VAF'] < vaf_th) & (tp_df['CramersV'] < v_th)
        fp_mask = (fp_df['VAF'] < vaf_th) & (fp_df['CramersV'] < v_th)
        
        tp_removed = tp_mask.sum()
        fp_removed = fp_mask.sum()
        
        tp_kept = len(tp_df) - tp_removed
        fp_kept = len(fp_df) - fp_removed
        
        precision = tp_kept / (tp_kept + fp_kept) if (tp_kept + fp_kept) > 0 else 0
        recall = tp_kept / SEQC2_TOTAL
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        
        results_remove_low.append({
            'vaf_threshold': vaf_th,
            'v_threshold': v_th,
            'tp_removed': tp_removed,
            'fp_removed': fp_removed,
            'tp_removed_pct': 100*tp_removed/len(tp_df),
            'fp_removed_pct': 100*fp_removed/len(fp_df),
            'tp_kept': tp_kept,
            'fp_kept': fp_kept,
            'precision': precision,
            'recall': recall,
            'f1_score': f1
        })
        
        tp_pct = 100*tp_removed/len(tp_df)
        fp_pct = 100*fp_removed/len(fp_df)
        print(f"{vaf_th:<10.2f} {v_th:<10.2f} {tp_removed:<12} {tp_pct:<10.2f} {fp_removed:<12} {fp_pct:<10.2f} {f1:<8.4f}")

results_low_df = pd.DataFrame(results_remove_low)
results_low_df.to_csv(f"{OUTPUT_DIR}/cramersv_vaf_analysis/filter_remove_low_vaf_and_low_v.csv", index=False)

# Find best F1
best_row = results_low_df.loc[results_low_df['f1_score'].idxmax()]
print(f"\n最佳 F1: {best_row['f1_score']:.4f} (VAF < {best_row['vaf_threshold']}, V < {best_row['v_threshold']})")

#----------------------------------------------
# 3. 過濾策略 2: 保留 (高 VAF OR 高 CramersV)
#----------------------------------------------
print("\n" + "=" * 70)
print("3. 過濾策略: 保留 (VAF >= threshold OR CramersV >= threshold)")
print("目的: 保留高支持度或有統計ASM證據的位點")
print("=" * 70)

results_keep_high = []

print("\n" + "-" * 85)
print(f"{'VAF >=':<10} {'V >=':<10} {'TP保留':<12} {'TP保留%':<10} {'FP保留':<12} {'FP保留%':<10} {'F1':<8}")
print("-" * 85)

for vaf_th in vaf_thresholds:
    for v_th in v_thresholds:
        # Keep if VAF >= vaf_th OR CramersV >= v_th
        tp_kept_mask = (tp_df['VAF'] >= vaf_th) | (tp_df['CramersV'] >= v_th)
        fp_kept_mask = (fp_df['VAF'] >= vaf_th) | (fp_df['CramersV'] >= v_th)
        
        tp_kept = tp_kept_mask.sum()
        fp_kept = fp_kept_mask.sum()
        tp_removed = len(tp_df) - tp_kept
        fp_removed = len(fp_df) - fp_kept
        
        precision = tp_kept / (tp_kept + fp_kept) if (tp_kept + fp_kept) > 0 else 0
        recall = tp_kept / SEQC2_TOTAL
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        
        results_keep_high.append({
            'vaf_threshold': vaf_th,
            'v_threshold': v_th,
            'tp_kept': tp_kept,
            'fp_kept': fp_kept,
            'tp_kept_pct': 100*tp_kept/len(tp_df),
            'fp_kept_pct': 100*fp_kept/len(fp_df),
            'tp_removed': tp_removed,
            'fp_removed': fp_removed,
            'precision': precision,
            'recall': recall,
            'f1_score': f1
        })
        
        tp_pct = 100*tp_kept/len(tp_df)
        fp_pct = 100*fp_kept/len(fp_df)
        print(f"{vaf_th:<10.2f} {v_th:<10.2f} {tp_kept:<12} {tp_pct:<10.2f} {fp_kept:<12} {fp_pct:<10.2f} {f1:<8.4f}")

results_high_df = pd.DataFrame(results_keep_high)
results_high_df.to_csv(f"{OUTPUT_DIR}/cramersv_vaf_analysis/filter_keep_high_vaf_or_high_v.csv", index=False)

#----------------------------------------------
# 4. 細粒度分析: VAF × CramersV 交叉表
#----------------------------------------------
print("\n" + "=" * 70)
print("4. VAF × CramersV 交叉表 (位點數量)")
print("=" * 70)

# Define bins
vaf_bins = [0, 0.10, 0.20, 0.30, 0.40, 0.50, 1.01]
vaf_labels = ['0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '>0.5']

v_bins = [0, 0.001, 0.10, 0.30, 0.50, 0.80, 1.01]
v_labels = ['0', '0-0.1', '0.1-0.3', '0.3-0.5', '0.5-0.8', '>0.8']

# Create binned columns
tp_df['VAF_bin'] = pd.cut(tp_df['VAF'], bins=vaf_bins, labels=vaf_labels, right=False)
tp_df['V_bin'] = pd.cut(tp_df['CramersV'], bins=v_bins, labels=v_labels, right=False)

fp_df['VAF_bin'] = pd.cut(fp_df['VAF'], bins=vaf_bins, labels=vaf_labels, right=False)
fp_df['V_bin'] = pd.cut(fp_df['CramersV'], bins=v_bins, labels=v_labels, right=False)

# Cross tables
tp_cross = pd.crosstab(tp_df['V_bin'], tp_df['VAF_bin'])
fp_cross = pd.crosstab(fp_df['V_bin'], fp_df['VAF_bin'])

print("\nTP 交叉表 (CramersV × VAF):")
print(tp_cross.to_string())

print("\nFP 交叉表 (CramersV × VAF):")
print(fp_cross.to_string())

# Calculate ratio table (FP/(TP+FP))
total_cross = tp_cross + fp_cross
ratio_cross = (fp_cross / total_cross * 100).round(1)
print("\nFP比例表 (FP/(TP+FP) × 100%):")
print(ratio_cross.fillna('-').to_string())

# Save cross tables
tp_cross.to_csv(f"{OUTPUT_DIR}/cramersv_vaf_analysis/crosstab_tp.csv")
fp_cross.to_csv(f"{OUTPUT_DIR}/cramersv_vaf_analysis/crosstab_fp.csv")
ratio_cross.to_csv(f"{OUTPUT_DIR}/cramersv_vaf_analysis/crosstab_fp_ratio.csv")

#----------------------------------------------
# 5. 特定門檻組合分析
#----------------------------------------------
print("\n" + "=" * 70)
print("5. 推薦門檻組合分析")
print("=" * 70)

# Baseline
baseline_precision = len(tp_df) / (len(tp_df) + len(fp_df))
baseline_recall = len(tp_df) / SEQC2_TOTAL
baseline_f1 = 2 * (baseline_precision * baseline_recall) / (baseline_precision + baseline_recall)

print(f"\n基線 F1: {baseline_f1:.4f} (Precision={baseline_precision:.4f}, Recall={baseline_recall:.4f})")

recommended_filters = [
    ("VAF < 0.15 AND V < 0.05", lambda df: (df['VAF'] < 0.15) & (df['CramersV'] < 0.05)),
    ("VAF < 0.20 AND V < 0.05", lambda df: (df['VAF'] < 0.20) & (df['CramersV'] < 0.05)),
    ("VAF < 0.24 AND V < 0.05", lambda df: (df['VAF'] < 0.24) & (df['CramersV'] < 0.05)),
    ("VAF < 0.20 AND V < 0.10", lambda df: (df['VAF'] < 0.20) & (df['CramersV'] < 0.10)),
    ("VAF < 0.15 AND V = 0", lambda df: (df['VAF'] < 0.15) & (df['CramersV'] == 0)),
    ("VAF < 0.20 AND V = 0", lambda df: (df['VAF'] < 0.20) & (df['CramersV'] == 0)),
    ("VAF < 0.24 AND V = 0", lambda df: (df['VAF'] < 0.24) & (df['CramersV'] == 0)),
]

print("\n" + "-" * 95)
print(f"{'過濾條件 (移除)':<30} {'TP移除':<10} {'FP移除':<10} {'Precision':<12} {'Recall':<12} {'F1':<10} {'ΔF1':<8}")
print("-" * 95)

recommended_results = []

for name, condition in recommended_filters:
    tp_mask = condition(tp_df)
    fp_mask = condition(fp_df)
    
    tp_removed = tp_mask.sum()
    fp_removed = fp_mask.sum()
    
    tp_kept = len(tp_df) - tp_removed
    fp_kept = len(fp_df) - fp_removed
    
    precision = tp_kept / (tp_kept + fp_kept) if (tp_kept + fp_kept) > 0 else 0
    recall = tp_kept / SEQC2_TOTAL
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    delta_f1 = f1 - baseline_f1
    
    print(f"{name:<30} {tp_removed:<10} {fp_removed:<10} {precision:<12.4f} {recall:<12.4f} {f1:<10.4f} {delta_f1:+8.4f}")
    
    recommended_results.append({
        'filter': name,
        'tp_removed': tp_removed,
        'fp_removed': fp_removed,
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'delta_f1': delta_f1
    })

rec_df = pd.DataFrame(recommended_results)
rec_df.to_csv(f"{OUTPUT_DIR}/cramersv_vaf_analysis/recommended_filters.csv", index=False)

#----------------------------------------------
# 6. 低 VAF 位點的詳細分析
#----------------------------------------------
print("\n" + "=" * 70)
print("6. 低 VAF 位點的詳細分析 (VAF < 0.24)")
print("=" * 70)

tp_low_vaf = tp_df[tp_df['VAF'] < 0.24]
fp_low_vaf = fp_df[fp_df['VAF'] < 0.24]

print(f"\n低 VAF 位點數量:")
print(f"  TP: {len(tp_low_vaf)} ({100*len(tp_low_vaf)/len(tp_df):.2f}%)")
print(f"  FP: {len(fp_low_vaf)} ({100*len(fp_low_vaf)/len(fp_df):.2f}%)")

print(f"\n低 VAF 位點中 CramersV 分布:")
print(f"  TP: V=0: {(tp_low_vaf['CramersV']==0).sum()}, V>0: {(tp_low_vaf['CramersV']>0).sum()}")
print(f"  FP: V=0: {(fp_low_vaf['CramersV']==0).sum()}, V>0: {(fp_low_vaf['CramersV']>0).sum()}")

# Calculate FP ratio in low VAF
fp_ratio_v0 = (fp_low_vaf['CramersV']==0).sum() / ((tp_low_vaf['CramersV']==0).sum() + (fp_low_vaf['CramersV']==0).sum()) * 100
fp_ratio_v_pos = (fp_low_vaf['CramersV']>0).sum() / ((tp_low_vaf['CramersV']>0).sum() + (fp_low_vaf['CramersV']>0).sum()) * 100 if ((tp_low_vaf['CramersV']>0).sum() + (fp_low_vaf['CramersV']>0).sum()) > 0 else 0

print(f"\n低 VAF 區 FP 比例:")
print(f"  V=0 位點的 FP 比例: {fp_ratio_v0:.1f}%")
print(f"  V>0 位點的 FP 比例: {fp_ratio_v_pos:.1f}%")

#----------------------------------------------
# 7. 生成 Heatmap 可視化
#----------------------------------------------
print("\n" + "=" * 70)
print("7. 生成可視化圖表")
print("=" * 70)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# 7.1 TP Heatmap
ax = axes[0, 0]
im = ax.imshow(tp_cross.values, cmap='Blues', aspect='auto')
ax.set_xticks(range(len(vaf_labels)))
ax.set_xticklabels(vaf_labels)
ax.set_yticks(range(len(v_labels)))
ax.set_yticklabels(v_labels)
ax.set_xlabel('VAF')
ax.set_ylabel('CramersV')
ax.set_title(f'TP Count (n={len(tp_df)})')
for i in range(len(v_labels)):
    for j in range(len(vaf_labels)):
        val = tp_cross.values[i, j]
        ax.text(j, i, f'{val}', ha='center', va='center', fontsize=9,
               color='white' if val > tp_cross.values.max()/2 else 'black')

# 7.2 FP Heatmap
ax = axes[0, 1]
im = ax.imshow(fp_cross.values, cmap='Reds', aspect='auto')
ax.set_xticks(range(len(vaf_labels)))
ax.set_xticklabels(vaf_labels)
ax.set_yticks(range(len(v_labels)))
ax.set_yticklabels(v_labels)
ax.set_xlabel('VAF')
ax.set_ylabel('CramersV')
ax.set_title(f'FP Count (n={len(fp_df)})')
for i in range(len(v_labels)):
    for j in range(len(vaf_labels)):
        val = fp_cross.values[i, j]
        ax.text(j, i, f'{val}', ha='center', va='center', fontsize=9,
               color='white' if val > fp_cross.values.max()/2 else 'black')

# 7.3 FP Ratio Heatmap
ax = axes[0, 2]
ratio_vals = ratio_cross.fillna(0).values
im = ax.imshow(ratio_vals, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=100)
ax.set_xticks(range(len(vaf_labels)))
ax.set_xticklabels(vaf_labels)
ax.set_yticks(range(len(v_labels)))
ax.set_yticklabels(v_labels)
ax.set_xlabel('VAF')
ax.set_ylabel('CramersV')
ax.set_title('FP Ratio (%)')
plt.colorbar(im, ax=ax, shrink=0.8)
for i in range(len(v_labels)):
    for j in range(len(vaf_labels)):
        val = ratio_vals[i, j]
        ax.text(j, i, f'{val:.0f}%' if val > 0 else '-', ha='center', va='center', fontsize=9,
               color='white' if val > 50 else 'black')

# 7.4 Scatter: VAF vs CramersV
ax = axes[1, 0]
ax.scatter(tp_df['VAF'], tp_df['CramersV'], alpha=0.3, s=5, c='blue', label='TP')
ax.scatter(fp_df['VAF'], fp_df['CramersV'], alpha=0.3, s=5, c='red', label='FP')
ax.axhline(y=0.05, color='gray', linestyle='--', alpha=0.7)
ax.axvline(x=0.24, color='gray', linestyle='--', alpha=0.7)
ax.set_xlabel('VAF')
ax.set_ylabel('CramersV')
ax.set_title('VAF vs CramersV')
ax.legend()
ax.set_xlim(0, 1)
ax.set_ylim(-0.02, 1.02)

# 7.5 F1 Score Heatmap
# Pivot results_low_df for heatmap
f1_pivot = results_low_df.pivot(index='v_threshold', columns='vaf_threshold', values='f1_score')
ax = axes[1, 1]
im = ax.imshow(f1_pivot.values, cmap='viridis', aspect='auto')
ax.set_xticks(range(len(vaf_thresholds)))
ax.set_xticklabels([f'{v:.2f}' for v in vaf_thresholds])
ax.set_yticks(range(len(v_thresholds)))
ax.set_yticklabels([f'{v:.2f}' for v in v_thresholds])
ax.set_xlabel('VAF Threshold')
ax.set_ylabel('CramersV Threshold')
ax.set_title('F1 Score (Remove if VAF<T AND V<T)')
plt.colorbar(im, ax=ax, shrink=0.8)
for i, v_th in enumerate(v_thresholds):
    for j, vaf_th in enumerate(vaf_thresholds):
        val = f1_pivot.values[i, j]
        ax.text(j, i, f'{val:.4f}', ha='center', va='center', fontsize=8,
               color='white' if val < f1_pivot.values.mean() else 'black')

# 7.6 Recommended filters comparison
ax = axes[1, 2]
rec_subset = rec_df.copy()
colors = ['red' if d < 0 else 'green' for d in rec_subset['delta_f1']]
bars = ax.barh(range(len(rec_subset)), rec_subset['delta_f1'], color=colors, edgecolor='black')
ax.set_yticks(range(len(rec_subset)))
ax.set_yticklabels(rec_subset['filter'], fontsize=9)
ax.set_xlabel('ΔF1 Score')
ax.set_title('Filter Impact on F1 Score')
ax.axvline(x=0, color='black', linewidth=1)
ax.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/cramersv_vaf_analysis/cramersv_vaf_analysis.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/cramersv_vaf_analysis/cramersv_vaf_analysis.png")

#----------------------------------------------
# 8. 最終總結
#----------------------------------------------
print("\n" + "=" * 70)
print("8. 分析總結")
print("=" * 70)

# Find best filter
best_filter = rec_df.loc[rec_df['f1_score'].idxmax()]
print(f"\n推薦過濾條件: {best_filter['filter']}")
print(f"  F1 Score: {best_filter['f1_score']:.4f}")
print(f"  ΔF1: {best_filter['delta_f1']:+.4f}")
print(f"  TP 移除: {best_filter['tp_removed']}")
print(f"  FP 移除: {best_filter['fp_removed']}")

# Save summary
summary = {
    "analysis_date": "2026-01-23",
    "description": "CramersV AND VAF 組合過濾分析",
    "baseline_f1": baseline_f1,
    "best_filter": {
        "condition": best_filter['filter'],
        "f1_score": best_filter['f1_score'],
        "delta_f1": best_filter['delta_f1'],
        "tp_removed": int(best_filter['tp_removed']),
        "fp_removed": int(best_filter['fp_removed'])
    },
    "low_vaf_analysis": {
        "tp_count": len(tp_low_vaf),
        "fp_count": len(fp_low_vaf),
        "fp_ratio_v0": fp_ratio_v0,
        "fp_ratio_v_positive": fp_ratio_v_pos
    }
}

with open(f"{OUTPUT_DIR}/cramersv_vaf_analysis/analysis_summary.json", "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print(f"\n所有結果已儲存至: {OUTPUT_DIR}/cramersv_vaf_analysis/")
print("分析完成!")
