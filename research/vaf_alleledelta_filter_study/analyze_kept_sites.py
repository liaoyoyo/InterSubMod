#!/usr/bin/env python3
"""
過濾後剩餘 TP/FP 狀況分析

分析過濾條件：AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24 刪除後
剩餘的 TP、FP 統計、特徵分布、VerificationClass 變化與案例分析

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
os.makedirs(f"{OUTPUT_DIR}/figures", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/data", exist_ok=True)

# Load data
print("載入資料...")
tp_df = pd.read_csv(f"{DATA_DIR}/tp_with_qual_vaf.csv")
fp_df = pd.read_csv(f"{DATA_DIR}/fp_with_qual_vaf.csv")

print(f"TP 總數: {len(tp_df)}")
print(f"FP 總數: {len(fp_df)}")

#----------------------------------------------
# 1. 應用過濾條件
#----------------------------------------------
def apply_filter(df, ad_thresh=0.25, v_thresh=0.05, vaf_thresh=0.24):
    """Return mask for rows that should be REMOVED"""
    return (df['AlleleDelta'] > ad_thresh) & (df['CramersV'] < v_thresh) & (df['VAF'] < vaf_thresh)

tp_removed_mask = apply_filter(tp_df)
fp_removed_mask = apply_filter(fp_df)

tp_kept = tp_df[~tp_removed_mask].copy()
fp_kept = fp_df[~fp_removed_mask].copy()
tp_removed = tp_df[tp_removed_mask].copy()
fp_removed = fp_df[fp_removed_mask].copy()

print(f"\n過濾條件: AD > 0.25 AND V < 0.05 AND VAF < 0.24")
print(f"TP 保留: {len(tp_kept)} ({100*len(tp_kept)/len(tp_df):.2f}%)")
print(f"FP 保留: {len(fp_kept)} ({100*len(fp_kept)/len(fp_df):.2f}%)")

#----------------------------------------------
# 2. 剩餘位點統計分析
#----------------------------------------------
print("\n" + "="*60)
print("剩餘 TP 統計")
print("="*60)

def detailed_stats(df, label):
    print(f"\n=== {label} ===")
    print(f"總數: {len(df)}")
    
    # 數值特徵
    print(f"\n數值特徵:")
    for col in ['AlleleDelta', 'CramersV', 'VAF', 'QUAL', 'HeuristicScore', 'NumReads', 'NumCpGs']:
        if col in df.columns:
            print(f"  {col}: mean={df[col].mean():.4f}, median={df[col].median():.4f}, std={df[col].std():.4f}")
    
    # VerificationClass
    if 'VerificationClass' in df.columns:
        print(f"\nVerificationClass 分布:")
        for vc, cnt in df['VerificationClass'].value_counts().items():
            print(f"  {vc}: {cnt} ({100*cnt/len(df):.2f}%)")
    
    # Significant
    if 'Significant' in df.columns:
        sig_count = df['Significant'].sum() if df['Significant'].dtype == bool else (df['Significant'] == True).sum()
        print(f"\nSignificant: {sig_count} ({100*sig_count/len(df):.2f}%)")
    
    # PassedGating
    if 'PassedGating' in df.columns:
        gated = df['PassedGating'].sum() if df['PassedGating'].dtype == bool else (df['PassedGating'] == True).sum()
        print(f"PassedGating: {gated} ({100*gated/len(df):.2f}%)")
    
    # DominantLabel
    if 'DominantLabel' in df.columns:
        print(f"\nDominantLabel 分布:")
        for dl, cnt in df['DominantLabel'].value_counts().items():
            print(f"  {dl}: {cnt} ({100*cnt/len(df):.2f}%)")
    
    return df.describe()

tp_kept_stats = detailed_stats(tp_kept, "剩餘 TP")
fp_kept_stats = detailed_stats(fp_kept, "剩餘 FP")

#----------------------------------------------
# 3. 過濾前後比較
#----------------------------------------------
print("\n" + "="*60)
print("過濾前後 VerificationClass 比較")
print("="*60)

def vc_comparison(before_df, after_df, label):
    print(f"\n{label}:")
    vc_before = before_df['VerificationClass'].value_counts()
    vc_after = after_df['VerificationClass'].value_counts()
    
    for vc in ['Strong', 'Weak', 'Subclone', 'Noise']:
        before_n = vc_before.get(vc, 0)
        after_n = vc_after.get(vc, 0)
        removed = before_n - after_n
        print(f"  {vc}: {before_n} -> {after_n} (移除 {removed}, {100*removed/before_n if before_n > 0 else 0:.1f}%)")

vc_comparison(tp_df, tp_kept, "TP VerificationClass 變化")
vc_comparison(fp_df, fp_kept, "FP VerificationClass 變化")

#----------------------------------------------
# 4. 剩餘 FP 分析（仍需關注的假陽性）
#----------------------------------------------
print("\n" + "="*60)
print("剩餘 FP 深入分析")
print("="*60)

# FP by VerificationClass
print("\n剩餘 FP 中的 Strong 位點（最需關注）:")
fp_kept_strong = fp_kept[fp_kept['VerificationClass'] == 'Strong']
print(f"  數量: {len(fp_kept_strong)} ({100*len(fp_kept_strong)/len(fp_kept):.2f}% of kept FP)")

if len(fp_kept_strong) > 0:
    print(f"  AlleleDelta: mean={fp_kept_strong['AlleleDelta'].mean():.4f}")
    print(f"  VAF: mean={fp_kept_strong['VAF'].mean():.4f}")
    print(f"  QUAL: mean={fp_kept_strong['QUAL'].mean():.4f}")

# High QUAL but still FP
fp_high_qual = fp_kept[fp_kept['QUAL'] > 0.9]
print(f"\n高品質 FP (QUAL > 0.9): {len(fp_high_qual)}")

# FP with high AlleleDelta but high VAF (not filtered)
fp_high_ad_high_vaf = fp_kept[(fp_kept['AlleleDelta'] > 0.25) & (fp_kept['VAF'] >= 0.24)]
print(f"高 AD + 高 VAF FP (AD>0.25, VAF>=0.24): {len(fp_high_ad_high_vaf)}")

#----------------------------------------------
# 5. 統計顯著性分析
#----------------------------------------------
print("\n" + "="*60)
print("統計顯著性分析（過濾後）")
print("="*60)

# Significance distribution
print("\nSignificant 位點分布:")
tp_sig = tp_kept['Significant'].sum() if tp_kept['Significant'].dtype == bool else (tp_kept['Significant'] == True).sum()
fp_sig = fp_kept['Significant'].sum() if fp_kept['Significant'].dtype == bool else (fp_kept['Significant'] == True).sum()
print(f"  TP Significant: {tp_sig}/{len(tp_kept)} ({100*tp_sig/len(tp_kept):.2f}%)")
print(f"  FP Significant: {fp_sig}/{len(fp_kept)} ({100*fp_sig/len(fp_kept):.2f}%)")

# PassedGating distribution
print("\nPassedGating 位點分布:")
tp_gated = tp_kept['PassedGating'].sum() if tp_kept['PassedGating'].dtype == bool else (tp_kept['PassedGating'] == True).sum()
fp_gated = fp_kept['PassedGating'].sum() if fp_kept['PassedGating'].dtype == bool else (fp_kept['PassedGating'] == True).sum()
print(f"  TP PassedGating: {tp_gated}/{len(tp_kept)} ({100*tp_gated/len(tp_kept):.2f}%)")
print(f"  FP PassedGating: {fp_gated}/{len(fp_kept)} ({100*fp_gated/len(fp_kept):.2f}%)")

#----------------------------------------------
# 6. 案例分析 - 剩餘的問題 FP
#----------------------------------------------
print("\n" + "="*60)
print("案例分析 - 剩餘的高風險 FP")
print("="*60)

# Top 10 highest AD FP that were kept
fp_kept_sorted = fp_kept.sort_values('AlleleDelta', ascending=False)
print("\n剩餘 FP 中 AlleleDelta 最高的 10 個位點:")
print(fp_kept_sorted[['Chr', 'Pos', 'AlleleDelta', 'CramersV', 'VAF', 'QUAL', 'VerificationClass', 'NumReads']].head(10).to_string())

# Top 10 highest QUAL FP that were kept  
fp_kept_qual_sorted = fp_kept.sort_values('QUAL', ascending=False)
print("\n剩餘 FP 中 QUAL 最高的 10 個位點:")
print(fp_kept_qual_sorted[['Chr', 'Pos', 'AlleleDelta', 'CramersV', 'VAF', 'QUAL', 'VerificationClass', 'NumReads']].head(10).to_string())

#----------------------------------------------
# 7. 生成圖表
#----------------------------------------------
print("\n" + "="*60)
print("生成可視化圖表")
print("="*60)

# Figure: Kept sites comparison
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Row 1: TP kept
axes[0, 0].hist(tp_kept['AlleleDelta'], bins=50, alpha=0.7, color='blue', edgecolor='black')
axes[0, 0].axvline(x=tp_kept['AlleleDelta'].mean(), color='red', linestyle='--', label=f'Mean={tp_kept["AlleleDelta"].mean():.3f}')
axes[0, 0].set_xlabel('AlleleDelta')
axes[0, 0].set_ylabel('Count')
axes[0, 0].set_title(f'TP Kept: AlleleDelta Distribution (n={len(tp_kept)})')
axes[0, 0].legend()

axes[0, 1].hist(tp_kept['VAF'], bins=50, alpha=0.7, color='blue', edgecolor='black')
axes[0, 1].axvline(x=tp_kept['VAF'].mean(), color='red', linestyle='--', label=f'Mean={tp_kept["VAF"].mean():.3f}')
axes[0, 1].set_xlabel('VAF')
axes[0, 1].set_title(f'TP Kept: VAF Distribution')
axes[0, 1].legend()

# TP VerificationClass pie
tp_vc = tp_kept['VerificationClass'].value_counts()
axes[0, 2].pie(tp_vc.values, labels=tp_vc.index, autopct='%1.1f%%', colors=['#2ecc71', '#f1c40f', '#e74c3c', '#95a5a6'])
axes[0, 2].set_title(f'TP Kept: VerificationClass')

# Row 2: FP kept
axes[1, 0].hist(fp_kept['AlleleDelta'], bins=50, alpha=0.7, color='red', edgecolor='black')
axes[1, 0].axvline(x=fp_kept['AlleleDelta'].mean(), color='blue', linestyle='--', label=f'Mean={fp_kept["AlleleDelta"].mean():.3f}')
axes[1, 0].set_xlabel('AlleleDelta')
axes[1, 0].set_ylabel('Count')
axes[1, 0].set_title(f'FP Kept: AlleleDelta Distribution (n={len(fp_kept)})')
axes[1, 0].legend()

axes[1, 1].hist(fp_kept['VAF'], bins=50, alpha=0.7, color='red', edgecolor='black')
axes[1, 1].axvline(x=fp_kept['VAF'].mean(), color='blue', linestyle='--', label=f'Mean={fp_kept["VAF"].mean():.3f}')
axes[1, 1].set_xlabel('VAF')
axes[1, 1].set_title(f'FP Kept: VAF Distribution')
axes[1, 1].legend()

# FP VerificationClass pie
fp_vc = fp_kept['VerificationClass'].value_counts()
axes[1, 2].pie(fp_vc.values, labels=fp_vc.index, autopct='%1.1f%%', colors=['#2ecc71', '#f1c40f', '#e74c3c', '#95a5a6'])
axes[1, 2].set_title(f'FP Kept: VerificationClass')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/kept_sites_overview.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/figures/kept_sites_overview.png")

# Figure: Before/After comparison
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# VerificationClass comparison
vc_order = ['Strong', 'Weak', 'Subclone', 'Noise']
x = np.arange(len(vc_order))
width = 0.35

# TP
tp_before_vc = [tp_df['VerificationClass'].value_counts().get(v, 0) for v in vc_order]
tp_after_vc = [tp_kept['VerificationClass'].value_counts().get(v, 0) for v in vc_order]
axes[0, 0].bar(x - width/2, tp_before_vc, width, label='Before', color='lightblue', edgecolor='black')
axes[0, 0].bar(x + width/2, tp_after_vc, width, label='After', color='blue', edgecolor='black')
axes[0, 0].set_xlabel('VerificationClass')
axes[0, 0].set_ylabel('Count')
axes[0, 0].set_title('TP: VerificationClass Before/After Filter')
axes[0, 0].set_xticks(x)
axes[0, 0].set_xticklabels(vc_order)
axes[0, 0].legend()

# FP
fp_before_vc = [fp_df['VerificationClass'].value_counts().get(v, 0) for v in vc_order]
fp_after_vc = [fp_kept['VerificationClass'].value_counts().get(v, 0) for v in vc_order]
axes[0, 1].bar(x - width/2, fp_before_vc, width, label='Before', color='lightcoral', edgecolor='black')
axes[0, 1].bar(x + width/2, fp_after_vc, width, label='After', color='red', edgecolor='black')
axes[0, 1].set_xlabel('VerificationClass')
axes[0, 1].set_ylabel('Count')
axes[0, 1].set_title('FP: VerificationClass Before/After Filter')
axes[0, 1].set_xticks(x)
axes[0, 1].set_xticklabels(vc_order)
axes[0, 1].legend()

# Summary metrics
metrics = ['Total', 'Strong', 'Significant', 'PassedGating']
tp_before = [len(tp_df), 
             tp_df['VerificationClass'].value_counts().get('Strong', 0),
             tp_df['Significant'].sum() if tp_df['Significant'].dtype == bool else (tp_df['Significant'] == True).sum(),
             tp_df['PassedGating'].sum() if tp_df['PassedGating'].dtype == bool else (tp_df['PassedGating'] == True).sum()]
tp_after = [len(tp_kept),
            tp_kept['VerificationClass'].value_counts().get('Strong', 0),
            tp_kept['Significant'].sum() if tp_kept['Significant'].dtype == bool else (tp_kept['Significant'] == True).sum(),
            tp_kept['PassedGating'].sum() if tp_kept['PassedGating'].dtype == bool else (tp_kept['PassedGating'] == True).sum()]

x2 = np.arange(len(metrics))
axes[1, 0].bar(x2 - width/2, tp_before, width, label='Before', color='lightblue', edgecolor='black')
axes[1, 0].bar(x2 + width/2, tp_after, width, label='After', color='blue', edgecolor='black')
axes[1, 0].set_xlabel('Metric')
axes[1, 0].set_ylabel('Count')
axes[1, 0].set_title('TP: Key Metrics Before/After Filter')
axes[1, 0].set_xticks(x2)
axes[1, 0].set_xticklabels(metrics)
axes[1, 0].legend()

# FP summary
fp_before = [len(fp_df),
             fp_df['VerificationClass'].value_counts().get('Strong', 0),
             fp_df['Significant'].sum() if fp_df['Significant'].dtype == bool else (fp_df['Significant'] == True).sum(),
             fp_df['PassedGating'].sum() if fp_df['PassedGating'].dtype == bool else (fp_df['PassedGating'] == True).sum()]
fp_after = [len(fp_kept),
            fp_kept['VerificationClass'].value_counts().get('Strong', 0),
            fp_kept['Significant'].sum() if fp_kept['Significant'].dtype == bool else (fp_kept['Significant'] == True).sum(),
            fp_kept['PassedGating'].sum() if fp_kept['PassedGating'].dtype == bool else (fp_kept['PassedGating'] == True).sum()]

axes[1, 1].bar(x2 - width/2, fp_before, width, label='Before', color='lightcoral', edgecolor='black')
axes[1, 1].bar(x2 + width/2, fp_after, width, label='After', color='red', edgecolor='black')
axes[1, 1].set_xlabel('Metric')
axes[1, 1].set_ylabel('Count')
axes[1, 1].set_title('FP: Key Metrics Before/After Filter')
axes[1, 1].set_xticks(x2)
axes[1, 1].set_xticklabels(metrics)
axes[1, 1].legend()

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/before_after_comparison.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/figures/before_after_comparison.png")

# Figure: Scatter of kept FP (highlighting risky zones)
fig, ax = plt.subplots(figsize=(12, 8))

# Plot kept TP
ax.scatter(tp_kept['VAF'], tp_kept['AlleleDelta'], c='blue', alpha=0.2, s=5, label='TP Kept')
# Plot kept FP
ax.scatter(fp_kept['VAF'], fp_kept['AlleleDelta'], c='red', alpha=0.4, s=10, label='FP Kept')

# Filter zone (removed)
ax.axhline(y=0.25, color='green', linestyle='--', linewidth=1, label='AD=0.25 threshold')
ax.axvline(x=0.24, color='green', linestyle='--', linewidth=1, label='VAF=0.24 threshold')
ax.fill_between([0, 0.24], 0.25, 0.8, alpha=0.1, color='green', label='Filtered Zone')

ax.set_xlabel('VAF')
ax.set_ylabel('AlleleDelta')
ax.set_title(f'Kept Sites: TP (n={len(tp_kept)}) vs FP (n={len(fp_kept)})')
ax.set_xlim(0, 1)
ax.set_ylim(0, 0.8)
ax.legend(loc='upper right')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/kept_sites_scatter.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/figures/kept_sites_scatter.png")

#----------------------------------------------
# 8. 儲存資料
#----------------------------------------------
print("\n" + "="*60)
print("儲存資料")
print("="*60)

# Save kept sites
tp_kept.to_csv(f"{OUTPUT_DIR}/data/tp_kept_after_filter.csv", index=False)
fp_kept.to_csv(f"{OUTPUT_DIR}/data/fp_kept_after_filter.csv", index=False)
print(f"已儲存: {OUTPUT_DIR}/data/tp_kept_after_filter.csv")
print(f"已儲存: {OUTPUT_DIR}/data/fp_kept_after_filter.csv")

# Save summary statistics
summary = {
    "filter_condition": "AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24",
    "tp_before": len(tp_df),
    "tp_after": len(tp_kept),
    "tp_removed": len(tp_removed),
    "fp_before": len(fp_df),
    "fp_after": len(fp_kept),
    "fp_removed": len(fp_removed),
    "tp_kept_stats": {
        "AlleleDelta_mean": round(tp_kept['AlleleDelta'].mean(), 4),
        "VAF_mean": round(tp_kept['VAF'].mean(), 4),
        "QUAL_mean": round(tp_kept['QUAL'].mean(), 4),
        "Strong_count": int(tp_kept['VerificationClass'].value_counts().get('Strong', 0)),
        "Strong_pct": round(100*tp_kept['VerificationClass'].value_counts().get('Strong', 0)/len(tp_kept), 2),
        "Significant_count": int(tp_kept['Significant'].sum()) if tp_kept['Significant'].dtype == bool else int((tp_kept['Significant'] == True).sum()),
        "PassedGating_count": int(tp_kept['PassedGating'].sum()) if tp_kept['PassedGating'].dtype == bool else int((tp_kept['PassedGating'] == True).sum())
    },
    "fp_kept_stats": {
        "AlleleDelta_mean": round(fp_kept['AlleleDelta'].mean(), 4),
        "VAF_mean": round(fp_kept['VAF'].mean(), 4),
        "QUAL_mean": round(fp_kept['QUAL'].mean(), 4),
        "Strong_count": int(fp_kept['VerificationClass'].value_counts().get('Strong', 0)),
        "Strong_pct": round(100*fp_kept['VerificationClass'].value_counts().get('Strong', 0)/len(fp_kept), 2),
        "Significant_count": int(fp_kept['Significant'].sum()) if fp_kept['Significant'].dtype == bool else int((fp_kept['Significant'] == True).sum()),
        "PassedGating_count": int(fp_kept['PassedGating'].sum()) if fp_kept['PassedGating'].dtype == bool else int((fp_kept['PassedGating'] == True).sum())
    }
}

with open(f"{OUTPUT_DIR}/data/kept_sites_summary.json", "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
print(f"已儲存: {OUTPUT_DIR}/data/kept_sites_summary.json")

# Sample cases for report
print("\n" + "="*60)
print("案例樣本 (用於報告)")
print("="*60)

# High risk FP samples
fp_high_risk = fp_kept[
    (fp_kept['VerificationClass'] == 'Strong') & 
    (fp_kept['QUAL'] > 0.8)
].sort_values('AlleleDelta', ascending=False).head(10)

print("\n高風險 FP 案例 (Strong + QUAL > 0.8):")
print(fp_high_risk[['Chr', 'Pos', 'Ref', 'Alt', 'AlleleDelta', 'CramersV', 'VAF', 'QUAL', 'NumReads', 'VerificationClass']].to_string())

fp_high_risk.to_csv(f"{OUTPUT_DIR}/data/high_risk_fp_cases.csv", index=False)

# Good TP samples
tp_good = tp_kept[
    (tp_kept['VerificationClass'] == 'Strong') & 
    (tp_kept['QUAL'] > 0.95) &
    (tp_kept['VAF'] > 0.4)
].sort_values('AlleleDelta').head(10)

print("\n良好 TP 案例 (Strong + QUAL > 0.95 + VAF > 0.4):")
print(tp_good[['Chr', 'Pos', 'Ref', 'Alt', 'AlleleDelta', 'CramersV', 'VAF', 'QUAL', 'NumReads', 'VerificationClass']].to_string())

tp_good.to_csv(f"{OUTPUT_DIR}/data/good_tp_cases.csv", index=False)

print("\n=== 分析完成 ===")
