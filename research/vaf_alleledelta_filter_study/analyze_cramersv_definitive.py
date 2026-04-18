#!/usr/bin/env python3
"""
CramersV 明確位點 (V > 0) 的 TP/FP 分析

目的:
1. 找出 CramersV > 0 (即 is_reliable=True) 的 TP 和 FP 位點
2. 評估這些位點的特殊狀況或數值差異
3. 生成觀察列表供深入分析

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
from scipy import stats

# Data paths
DATA_DIR = "/big8_disk/liaoyoyo2001/InterSubMod/output/ml_feature_exploration/2026/01/data"
OUTPUT_DIR = "/big8_disk/liaoyoyo2001/InterSubMod/analysis/vaf_alleledelta_filter_study"
os.makedirs(f"{OUTPUT_DIR}/cramersv_sites", exist_ok=True)

# Load data
print("=" * 60)
print("CramersV 明確位點 (V > 0) 分析")
print("=" * 60)
print("\n載入資料...")
tp_df = pd.read_csv(f"{DATA_DIR}/tp_with_qual_vaf.csv")
fp_df = pd.read_csv(f"{DATA_DIR}/fp_with_qual_vaf.csv")

print(f"TP 總數: {len(tp_df)}")
print(f"FP 總數: {len(fp_df)}")

#----------------------------------------------
# 1. 找出 CramersV > 0 的位點
#----------------------------------------------
print("\n" + "=" * 60)
print("1. CramersV 分布概覽")
print("=" * 60)

# CramersV thresholds for analysis
V_DEFINITIVE = 0.0  # V > 0 表示 is_reliable=True

tp_v_definitive = tp_df[tp_df['CramersV'] > V_DEFINITIVE]
fp_v_definitive = fp_df[fp_df['CramersV'] > V_DEFINITIVE]

print(f"\nCramersV > 0 (明確計算出的 V 值):")
print(f"  TP: {len(tp_v_definitive)} / {len(tp_df)} ({100*len(tp_v_definitive)/len(tp_df):.2f}%)")
print(f"  FP: {len(fp_v_definitive)} / {len(fp_df)} ({100*len(fp_v_definitive)/len(fp_df):.2f}%)")

# More detailed breakdown
v_ranges = [
    (0.0, 0.05, "極低 (0 ~ 0.05)"),
    (0.05, 0.20, "低 (0.05 ~ 0.20)"),
    (0.20, 0.50, "中 (0.20 ~ 0.50)"),
    (0.50, 0.80, "高 (0.50 ~ 0.80)"),
    (0.80, 1.01, "極高 (0.80 ~ 1.00)")
]

print("\n\nCramersV 區間分布 (僅 V > 0 的位點):")
print("-" * 60)
print(f"{'區間':<20} {'TP 數量':<12} {'TP %':<10} {'FP 數量':<12} {'FP %':<10}")
print("-" * 60)

tp_v_pos = tp_df[tp_df['CramersV'] > 0]
fp_v_pos = fp_df[fp_df['CramersV'] > 0]

for low, high, label in v_ranges:
    tp_in_range = tp_v_pos[(tp_v_pos['CramersV'] > low) & (tp_v_pos['CramersV'] <= high)]
    fp_in_range = fp_v_pos[(fp_v_pos['CramersV'] > low) & (fp_v_pos['CramersV'] <= high)]
    
    tp_pct = 100 * len(tp_in_range) / len(tp_v_pos) if len(tp_v_pos) > 0 else 0
    fp_pct = 100 * len(fp_in_range) / len(fp_v_pos) if len(fp_v_pos) > 0 else 0
    
    print(f"{label:<20} {len(tp_in_range):<12} {tp_pct:<10.2f} {len(fp_in_range):<12} {fp_pct:<10.2f}")

#----------------------------------------------
# 2. 統計比較分析
#----------------------------------------------
print("\n" + "=" * 60)
print("2. TP vs FP 的 CramersV 統計比較 (僅 V > 0)")
print("=" * 60)

def stat_summary(data, name):
    """Print statistical summary"""
    if len(data) == 0:
        print(f"  {name}: 無數據")
        return None
    return {
        'n': len(data),
        'mean': data.mean(),
        'std': data.std(),
        'median': data.median(),
        'q25': data.quantile(0.25),
        'q75': data.quantile(0.75),
        'min': data.min(),
        'max': data.max()
    }

tp_v_stats = stat_summary(tp_v_pos['CramersV'], 'TP')
fp_v_stats = stat_summary(fp_v_pos['CramersV'], 'FP')

print("\nTP (CramersV > 0):")
if tp_v_stats:
    print(f"  n={tp_v_stats['n']}, mean={tp_v_stats['mean']:.4f}, std={tp_v_stats['std']:.4f}")
    print(f"  median={tp_v_stats['median']:.4f}, IQR=[{tp_v_stats['q25']:.4f}, {tp_v_stats['q75']:.4f}]")
    print(f"  range=[{tp_v_stats['min']:.4f}, {tp_v_stats['max']:.4f}]")

print("\nFP (CramersV > 0):")
if fp_v_stats:
    print(f"  n={fp_v_stats['n']}, mean={fp_v_stats['mean']:.4f}, std={fp_v_stats['std']:.4f}")
    print(f"  median={fp_v_stats['median']:.4f}, IQR=[{fp_v_stats['q25']:.4f}, {fp_v_stats['q75']:.4f}]")
    print(f"  range=[{fp_v_stats['min']:.4f}, {fp_v_stats['max']:.4f}]")

# Mann-Whitney U test
if len(tp_v_pos) > 0 and len(fp_v_pos) > 0:
    u_stat, u_pval = stats.mannwhitneyu(tp_v_pos['CramersV'], fp_v_pos['CramersV'], alternative='two-sided')
    print(f"\nMann-Whitney U 檢定: U={u_stat:.0f}, p={u_pval:.4e}")
    
    # Effect size (rank-biserial correlation)
    n1, n2 = len(tp_v_pos), len(fp_v_pos)
    r_biserial = 1 - (2*u_stat) / (n1 * n2)
    print(f"Rank-biserial correlation: {r_biserial:.4f}")

#----------------------------------------------
# 3. 多特徵比較
#----------------------------------------------
print("\n" + "=" * 60)
print("3. CramersV > 0 位點的多特徵比較")
print("=" * 60)

features_to_compare = ['AlleleDelta', 'VAF', 'QUAL', 'NumReads', 'HPMergedDelta', 'HeuristicScore']

comparison_results = []
print("\n" + "-" * 80)
print(f"{'特徵':<18} {'TP mean':<12} {'FP mean':<12} {'TP median':<12} {'FP median':<12} {'差異顯著'}")
print("-" * 80)

for feat in features_to_compare:
    if feat in tp_v_pos.columns and feat in fp_v_pos.columns:
        tp_vals = tp_v_pos[feat].dropna()
        fp_vals = fp_v_pos[feat].dropna()
        
        if len(tp_vals) > 0 and len(fp_vals) > 0:
            u_stat, p_val = stats.mannwhitneyu(tp_vals, fp_vals, alternative='two-sided')
            sig = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else ""))
            
            print(f"{feat:<18} {tp_vals.mean():<12.4f} {fp_vals.mean():<12.4f} {tp_vals.median():<12.4f} {fp_vals.median():<12.4f} {sig} (p={p_val:.4e})")
            
            comparison_results.append({
                'feature': feat,
                'tp_mean': tp_vals.mean(),
                'fp_mean': fp_vals.mean(),
                'tp_median': tp_vals.median(),
                'fp_median': fp_vals.median(),
                'p_value': p_val,
                'significant': p_val < 0.05
            })

#----------------------------------------------
# 4. VerificationClass 分布
#----------------------------------------------
print("\n" + "=" * 60)
print("4. VerificationClass 分布 (CramersV > 0)")
print("=" * 60)

if 'VerificationClass' in tp_v_pos.columns:
    print("\nTP (CramersV > 0) VerificationClass 分布:")
    for vc, cnt in tp_v_pos['VerificationClass'].value_counts().items():
        print(f"  {vc}: {cnt} ({100*cnt/len(tp_v_pos):.1f}%)")

if 'VerificationClass' in fp_v_pos.columns:
    print("\nFP (CramersV > 0) VerificationClass 分布:")
    for vc, cnt in fp_v_pos['VerificationClass'].value_counts().items():
        print(f"  {vc}: {cnt} ({100*cnt/len(fp_v_pos):.1f}%)")

#----------------------------------------------
# 5. 生成觀察列表
#----------------------------------------------
print("\n" + "=" * 60)
print("5. 生成觀察列表")
print("=" * 60)

# Sort by CramersV descending to see highest V sites first
tp_v_sorted = tp_v_pos.sort_values('CramersV', ascending=False)
fp_v_sorted = fp_v_pos.sort_values('CramersV', ascending=False)

# Select key columns for observation
obs_cols = ['RegionID', 'Chr', 'Pos', 'Ref', 'Alt', 'NumReads', 'CramersV', 
            'AlleleDelta', 'VAF', 'QUAL', 'HPMergedDelta', 'HPMergedSig', 
            'HeuristicScore', 'VerificationClass', 'DominantLabel', 'PassedGating']

# Ensure columns exist
obs_cols = [c for c in obs_cols if c in tp_v_pos.columns]

# Save full lists
tp_v_sorted[obs_cols].to_csv(f"{OUTPUT_DIR}/cramersv_sites/tp_cramersv_definitive.csv", index=False)
fp_v_sorted[obs_cols].to_csv(f"{OUTPUT_DIR}/cramersv_sites/fp_cramersv_definitive.csv", index=False)

print(f"已儲存 TP 觀察列表: {OUTPUT_DIR}/cramersv_sites/tp_cramersv_definitive.csv ({len(tp_v_sorted)} 筆)")
print(f"已儲存 FP 觀察列表: {OUTPUT_DIR}/cramersv_sites/fp_cramersv_definitive.csv ({len(fp_v_sorted)} 筆)")

#----------------------------------------------
# 6. 高 V 位點的詳細分析 (V >= 0.5)
#----------------------------------------------
print("\n" + "=" * 60)
print("6. 高 CramersV 位點分析 (V >= 0.5)")
print("=" * 60)

tp_high_v = tp_v_pos[tp_v_pos['CramersV'] >= 0.5]
fp_high_v = fp_v_pos[fp_v_pos['CramersV'] >= 0.5]

print(f"\nTP with V >= 0.5: {len(tp_high_v)}")
print(f"FP with V >= 0.5: {len(fp_high_v)}")

if len(tp_high_v) > 0:
    print(f"\nTP 高 V 位點特徵:")
    print(f"  AlleleDelta: mean={tp_high_v['AlleleDelta'].mean():.4f}, median={tp_high_v['AlleleDelta'].median():.4f}")
    print(f"  VAF: mean={tp_high_v['VAF'].mean():.4f}, median={tp_high_v['VAF'].median():.4f}")
    print(f"  QUAL: mean={tp_high_v['QUAL'].mean():.4f}, median={tp_high_v['QUAL'].median():.4f}")
    print(f"  NumReads: mean={tp_high_v['NumReads'].mean():.1f}, median={tp_high_v['NumReads'].median():.1f}")
    if 'PassedGating' in tp_high_v.columns:
        print(f"  PassedGating=True: {tp_high_v['PassedGating'].sum()} ({100*tp_high_v['PassedGating'].sum()/len(tp_high_v):.1f}%)") 

if len(fp_high_v) > 0:
    print(f"\nFP 高 V 位點特徵:")
    print(f"  AlleleDelta: mean={fp_high_v['AlleleDelta'].mean():.4f}, median={fp_high_v['AlleleDelta'].median():.4f}")
    print(f"  VAF: mean={fp_high_v['VAF'].mean():.4f}, median={fp_high_v['VAF'].median():.4f}")
    print(f"  QUAL: mean={fp_high_v['QUAL'].mean():.4f}, median={fp_high_v['QUAL'].median():.4f}")
    print(f"  NumReads: mean={fp_high_v['NumReads'].mean():.1f}, median={fp_high_v['NumReads'].median():.1f}")
    if 'PassedGating' in fp_high_v.columns:
        print(f"  PassedGating=True: {fp_high_v['PassedGating'].sum()} ({100*fp_high_v['PassedGating'].sum()/len(fp_high_v):.1f}%)")

# Save high V sites
tp_high_v[obs_cols].to_csv(f"{OUTPUT_DIR}/cramersv_sites/tp_high_cramersv.csv", index=False)
fp_high_v[obs_cols].to_csv(f"{OUTPUT_DIR}/cramersv_sites/fp_high_cramersv.csv", index=False)
print(f"\n已儲存高 V 位點列表:")
print(f"  TP: {OUTPUT_DIR}/cramersv_sites/tp_high_cramersv.csv ({len(tp_high_v)} 筆)")
print(f"  FP: {OUTPUT_DIR}/cramersv_sites/fp_high_cramersv.csv ({len(fp_high_v)} 筆)")

#----------------------------------------------
# 7. 特殊情況分析
#----------------------------------------------
print("\n" + "=" * 60)
print("7. 特殊情況分析")
print("=" * 60)

# Case 1: High V + High AlleleDelta (strong ASM signal)
tp_strong_asm = tp_v_pos[(tp_v_pos['CramersV'] >= 0.5) & (tp_v_pos['AlleleDelta'] > 0.25)]
fp_strong_asm = fp_v_pos[(fp_v_pos['CramersV'] >= 0.5) & (fp_v_pos['AlleleDelta'] > 0.25)]
print(f"\n強 ASM 信號 (V >= 0.5 AND AD > 0.25):")
print(f"  TP: {len(tp_strong_asm)}")
print(f"  FP: {len(fp_strong_asm)}")

# Case 2: High V but Low AlleleDelta (paradox)
tp_paradox = tp_v_pos[(tp_v_pos['CramersV'] >= 0.5) & (tp_v_pos['AlleleDelta'] <= 0.1)]
fp_paradox = fp_v_pos[(fp_v_pos['CramersV'] >= 0.5) & (fp_v_pos['AlleleDelta'] <= 0.1)]
print(f"\n矛盾情況 (V >= 0.5 但 AD <= 0.1):")
print(f"  TP: {len(tp_paradox)}")
print(f"  FP: {len(fp_paradox)}")

# Case 3: Low VAF with high V (interesting for filtering)
tp_low_vaf_high_v = tp_v_pos[(tp_v_pos['VAF'] < 0.24) & (tp_v_pos['CramersV'] >= 0.5)]
fp_low_vaf_high_v = fp_v_pos[(fp_v_pos['VAF'] < 0.24) & (fp_v_pos['CramersV'] >= 0.5)]
print(f"\n低 VAF + 高 V (VAF < 0.24 AND V >= 0.5):")
print(f"  TP: {len(tp_low_vaf_high_v)}")
print(f"  FP: {len(fp_low_vaf_high_v)}")

# Case 4: HPMergedSig correlation
if 'HPMergedSig' in tp_v_pos.columns:
    tp_hpsig_true = tp_v_pos[tp_v_pos['HPMergedSig'] == True]
    fp_hpsig_true = fp_v_pos[fp_v_pos['HPMergedSig'] == True]
    print(f"\nHPMergedSig = True 且 CramersV > 0:")
    print(f"  TP: {len(tp_hpsig_true)} ({100*len(tp_hpsig_true)/len(tp_v_pos):.1f}% of V>0 TP)")
    print(f"  FP: {len(fp_hpsig_true)} ({100*len(fp_hpsig_true)/len(fp_v_pos):.1f}% of V>0 FP)")

#----------------------------------------------
# 8. 生成可視化
#----------------------------------------------
print("\n" + "=" * 60)
print("8. 生成可視化圖表")
print("=" * 60)

# Figure 1: CramersV distribution for V > 0
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 8.1 CramersV histogram (V > 0 only)
ax = axes[0, 0]
ax.hist(tp_v_pos['CramersV'], bins=30, alpha=0.6, label=f'TP (n={len(tp_v_pos)})', color='blue', edgecolor='black')
ax.hist(fp_v_pos['CramersV'], bins=30, alpha=0.6, label=f'FP (n={len(fp_v_pos)})', color='red', edgecolor='black')
ax.set_xlabel('CramersV')
ax.set_ylabel('Count')
ax.set_title('CramersV Distribution (V > 0 only)')
ax.legend()
ax.grid(True, alpha=0.3)

# 8.2 CramersV vs AlleleDelta scatter
ax = axes[0, 1]
ax.scatter(tp_v_pos['AlleleDelta'], tp_v_pos['CramersV'], alpha=0.4, s=20, label=f'TP', c='blue')
ax.scatter(fp_v_pos['AlleleDelta'], fp_v_pos['CramersV'], alpha=0.4, s=20, label=f'FP', c='red')
ax.set_xlabel('AlleleDelta')
ax.set_ylabel('CramersV')
ax.set_title('CramersV vs AlleleDelta (V > 0)')
ax.legend()
ax.grid(True, alpha=0.3)

# 8.3 CramersV vs VAF scatter
ax = axes[1, 0]
ax.scatter(tp_v_pos['VAF'], tp_v_pos['CramersV'], alpha=0.4, s=20, label=f'TP', c='blue')
ax.scatter(fp_v_pos['VAF'], fp_v_pos['CramersV'], alpha=0.4, s=20, label=f'FP', c='red')
ax.set_xlabel('VAF')
ax.set_ylabel('CramersV')
ax.set_title('CramersV vs VAF (V > 0)')
ax.axvline(x=0.24, color='green', linestyle='--', alpha=0.7, label='VAF=0.24')
ax.legend()
ax.grid(True, alpha=0.3)

# 8.4 CramersV vs NumReads
ax = axes[1, 1]
ax.scatter(tp_v_pos['NumReads'], tp_v_pos['CramersV'], alpha=0.4, s=20, label=f'TP', c='blue')
ax.scatter(fp_v_pos['NumReads'], fp_v_pos['CramersV'], alpha=0.4, s=20, label=f'FP', c='red')
ax.set_xlabel('NumReads')
ax.set_ylabel('CramersV')
ax.set_title('CramersV vs NumReads (V > 0)')
ax.set_xlim(0, 200)  # Limit for better visualization
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/cramersv_sites/cramersv_definitive_analysis.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/cramersv_sites/cramersv_definitive_analysis.png")

# Figure 2: Box plots comparison
fig, axes = plt.subplots(1, 4, figsize=(16, 5))

for i, feat in enumerate(['CramersV', 'AlleleDelta', 'VAF', 'NumReads']):
    ax = axes[i]
    data_tp = tp_v_pos[feat].dropna()
    data_fp = fp_v_pos[feat].dropna()
    
    bp = ax.boxplot([data_tp, data_fp], labels=['TP', 'FP'], patch_artist=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    ax.set_ylabel(feat)
    ax.set_title(f'{feat} (V > 0)')
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/cramersv_sites/cramersv_boxplots.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/cramersv_sites/cramersv_boxplots.png")

#----------------------------------------------
# 9. 儲存完整分析結果
#----------------------------------------------
print("\n" + "=" * 60)
print("9. 儲存分析結果摘要")
print("=" * 60)

analysis_summary = {
    "analysis_date": "2026-01-23",
    "description": "CramersV 明確位點 (V > 0) 的 TP 和 FP 分析",
    
    "overview": {
        "tp_total": len(tp_df),
        "fp_total": len(fp_df),
        "tp_v_positive": len(tp_v_pos),
        "fp_v_positive": len(fp_v_pos),
        "tp_v_positive_pct": round(100*len(tp_v_pos)/len(tp_df), 2),
        "fp_v_positive_pct": round(100*len(fp_v_pos)/len(fp_df), 2)
    },
    
    "cramersv_stats": {
        "tp": tp_v_stats if tp_v_stats else None,
        "fp": fp_v_stats if fp_v_stats else None
    },
    
    "high_v_sites": {
        "threshold": 0.5,
        "tp_count": len(tp_high_v),
        "fp_count": len(fp_high_v)
    },
    
    "special_cases": {
        "strong_asm": {
            "description": "V >= 0.5 AND AD > 0.25",
            "tp": len(tp_strong_asm),
            "fp": len(fp_strong_asm)
        },
        "paradox": {
            "description": "V >= 0.5 AND AD <= 0.1",
            "tp": len(tp_paradox),
            "fp": len(fp_paradox)
        },
        "low_vaf_high_v": {
            "description": "VAF < 0.24 AND V >= 0.5",
            "tp": len(tp_low_vaf_high_v),
            "fp": len(fp_low_vaf_high_v)
        }
    },
    
    "output_files": {
        "tp_all_v_positive": "cramersv_sites/tp_cramersv_definitive.csv",
        "fp_all_v_positive": "cramersv_sites/fp_cramersv_definitive.csv",
        "tp_high_v": "cramersv_sites/tp_high_cramersv.csv",
        "fp_high_v": "cramersv_sites/fp_high_cramersv.csv"
    }
}

# Convert numpy types for JSON serialization
def convert_numpy(obj):
    if isinstance(obj, dict):
        return {k: convert_numpy(v) for k, v in obj.items()}
    elif isinstance(obj, (np.integer, np.int64)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

analysis_summary = convert_numpy(analysis_summary)

with open(f"{OUTPUT_DIR}/cramersv_sites/analysis_summary.json", "w") as f:
    json.dump(analysis_summary, f, indent=2, ensure_ascii=False)

print(f"已儲存分析摘要: {OUTPUT_DIR}/cramersv_sites/analysis_summary.json")

print("\n" + "=" * 60)
print("分析完成!")
print("=" * 60)
print(f"\n所有輸出檔案位於: {OUTPUT_DIR}/cramersv_sites/")
print("\n主要發現:")
print(f"  1. TP 中有 {100*len(tp_v_pos)/len(tp_df):.2f}% 位點具有明確的 CramersV (V > 0)")
print(f"  2. FP 中有 {100*len(fp_v_pos)/len(fp_df):.2f}% 位點具有明確的 CramersV (V > 0)")
print(f"  3. 高 V (>= 0.5) 位點: TP={len(tp_high_v)}, FP={len(fp_high_v)}")
