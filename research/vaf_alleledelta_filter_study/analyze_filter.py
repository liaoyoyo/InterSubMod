#!/usr/bin/env python3
"""
VAF + AlleleDelta 組合篩選策略分析

分析過濾條件：AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24
不採用 QUAL 過濾，專注於甲基化特徵與 VAF 組合效果

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
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/figures", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/data", exist_ok=True)

# Load data
print("載入資料...")
tp_df = pd.read_csv(f"{DATA_DIR}/tp_with_qual_vaf.csv")
fp_df = pd.read_csv(f"{DATA_DIR}/fp_with_qual_vaf.csv")

print(f"TP 總數: {len(tp_df)}")
print(f"FP 總數: {len(fp_df)}")

#----------------------------------------------
# 1. 核心過濾條件分析
#----------------------------------------------
print("\n=== 1. 核心過濾條件分析 ===")

# Filter condition: AlleleDelta > 0.25 AND CramersV < 0.05 AND VAF < 0.24
def apply_filter(df, ad_thresh=0.25, v_thresh=0.05, vaf_thresh=0.24):
    """Return mask for rows that should be REMOVED"""
    return (df['AlleleDelta'] > ad_thresh) & (df['CramersV'] < v_thresh) & (df['VAF'] < vaf_thresh)

# Apply filter
tp_removed_mask = apply_filter(tp_df)
fp_removed_mask = apply_filter(fp_df)

tp_removed = tp_df[tp_removed_mask]
tp_kept = tp_df[~tp_removed_mask]
fp_removed = fp_df[fp_removed_mask]
fp_kept = fp_df[~fp_removed_mask]

print(f"\n過濾條件: AD > 0.25 AND V < 0.05 AND VAF < 0.24 (不含 QUAL)")
print(f"TP 被移除: {len(tp_removed)} ({100*len(tp_removed)/len(tp_df):.4f}%)")
print(f"FP 被移除: {len(fp_removed)} ({100*len(fp_removed)/len(fp_df):.2f}%)")
print(f"TP 保留: {len(tp_kept)}")
print(f"FP 保留: {len(fp_kept)}")

# Calculate F1
SEQC2_TOTAL = 39447
tp_kept_n = len(tp_kept)
fp_kept_n = len(fp_kept)
fn = SEQC2_TOTAL - tp_kept_n

precision = tp_kept_n / (tp_kept_n + fp_kept_n)
recall = tp_kept_n / SEQC2_TOTAL
f1 = 2 * (precision * recall) / (precision + recall)

print(f"\n=== F1-Score 計算 (SEQC2 標準: {SEQC2_TOTAL}) ===")
print(f"Precision: {precision:.4f}")
print(f"Recall: {recall:.4f}")
print(f"F1-Score: {f1:.4f}")

# Baseline (no filter)
baseline_precision = len(tp_df) / (len(tp_df) + len(fp_df))
baseline_recall = len(tp_df) / SEQC2_TOTAL
baseline_f1 = 2 * (baseline_precision * baseline_recall) / (baseline_precision + baseline_recall)

print(f"\n=== 基線 (無過濾) ===")
print(f"Baseline F1: {baseline_f1:.4f}")
print(f"F1 提升: +{f1 - baseline_f1:.4f} ({100*(f1 - baseline_f1)/baseline_f1:.2f}%)")

results_summary = {
    "filter_condition": "AD > 0.25 AND V < 0.05 AND VAF < 0.24",
    "qual_filter": "未採用",
    "tp_total": len(tp_df),
    "fp_total": len(fp_df),
    "tp_removed": len(tp_removed),
    "fp_removed": len(fp_removed),
    "tp_kept": int(tp_kept_n),
    "fp_kept": int(fp_kept_n),
    "seqc2_total": SEQC2_TOTAL,
    "precision": round(precision, 4),
    "recall": round(recall, 4),
    "f1_score": round(f1, 4),
    "baseline_f1": round(baseline_f1, 4),
    "f1_improvement": round(f1 - baseline_f1, 4)
}

with open(f"{OUTPUT_DIR}/data/filter_results_summary.json", "w") as f:
    json.dump(results_summary, f, indent=2)

#----------------------------------------------
# 2. 被移除位點特徵分析
#----------------------------------------------
print("\n=== 2. 被移除位點特徵分析 ===")

def analyze_removed(df, label):
    print(f"\n{label} 被移除位點統計:")
    if len(df) == 0:
        print("  無被移除位點")
        return
    
    print(f"  總數: {len(df)}")
    print(f"  AlleleDelta: mean={df['AlleleDelta'].mean():.4f}, median={df['AlleleDelta'].median():.4f}, max={df['AlleleDelta'].max():.4f}")
    print(f"  CramersV: mean={df['CramersV'].mean():.4f}, median={df['CramersV'].median():.4f}, max={df['CramersV'].max():.4f}")
    print(f"  VAF: mean={df['VAF'].mean():.4f}, median={df['VAF'].median():.4f}, max={df['VAF'].max():.4f}")
    print(f"  QUAL: mean={df['QUAL'].mean():.4f}, median={df['QUAL'].median():.4f}")
    print(f"  NumReads: mean={df['NumReads'].mean():.2f}, median={df['NumReads'].median():.2f}")
    
    if 'VerificationClass' in df.columns:
        print(f"  VerificationClass 分布:")
        for vc, cnt in df['VerificationClass'].value_counts().items():
            print(f"    {vc}: {cnt} ({100*cnt/len(df):.1f}%)")

analyze_removed(tp_removed, "TP")
analyze_removed(fp_removed, "FP")

# Save removed sites
tp_removed.to_csv(f"{OUTPUT_DIR}/data/tp_removed_by_filter.csv", index=False)
fp_removed.to_csv(f"{OUTPUT_DIR}/data/fp_removed_by_filter.csv", index=False)

#----------------------------------------------
# 3. 染色體分布分析
#----------------------------------------------
print("\n=== 3. 染色體分布分析 ===")

tp_removed_chr = tp_removed.groupby('Chr').size().sort_index()
fp_removed_chr = fp_removed.groupby('Chr').size().sort_index()

print("\nTP 被移除 - 染色體分布:")
for chr_name, cnt in tp_removed_chr.items():
    print(f"  {chr_name}: {cnt}")

print("\nFP 被移除 - 染色體分布:")
for chr_name, cnt in fp_removed_chr.items():
    print(f"  {chr_name}: {cnt}")

#----------------------------------------------
# 4. VAF 門檻敏感度分析
#----------------------------------------------
print("\n=== 4. VAF 門檻敏感度分析 ===")

vaf_thresholds = np.arange(0.10, 0.35, 0.02)
sensitivity_results = []

for vaf_th in vaf_thresholds:
    tp_rem = apply_filter(tp_df, vaf_thresh=vaf_th).sum()
    fp_rem = apply_filter(fp_df, vaf_thresh=vaf_th).sum()
    tp_kpt = len(tp_df) - tp_rem
    fp_kpt = len(fp_df) - fp_rem
    
    prec = tp_kpt / (tp_kpt + fp_kpt) if (tp_kpt + fp_kpt) > 0 else 0
    rec = tp_kpt / SEQC2_TOTAL
    f1_val = 2 * (prec * rec) / (prec + rec) if (prec + rec) > 0 else 0
    
    sensitivity_results.append({
        "vaf_threshold": round(vaf_th, 2),
        "tp_removed": tp_rem,
        "fp_removed": fp_rem,
        "tp_kept": tp_kpt,
        "fp_kept": fp_kpt,
        "precision": round(prec, 4),
        "recall": round(rec, 4),
        "f1_score": round(f1_val, 4)
    })
    print(f"VAF < {vaf_th:.2f}: TP移除={tp_rem}, FP移除={fp_rem}, F1={f1_val:.4f}")

sensitivity_df = pd.DataFrame(sensitivity_results)
sensitivity_df.to_csv(f"{OUTPUT_DIR}/data/vaf_sensitivity_analysis.csv", index=False)

#----------------------------------------------
# 5. AlleleDelta 門檻敏感度分析
#----------------------------------------------
print("\n=== 5. AlleleDelta 門檻敏感度分析 ===")

ad_thresholds = np.arange(0.15, 0.35, 0.02)
ad_sensitivity_results = []

for ad_th in ad_thresholds:
    tp_rem = apply_filter(tp_df, ad_thresh=ad_th).sum()
    fp_rem = apply_filter(fp_df, ad_thresh=ad_th).sum()
    tp_kpt = len(tp_df) - tp_rem
    fp_kpt = len(fp_df) - fp_rem
    
    prec = tp_kpt / (tp_kpt + fp_kpt) if (tp_kpt + fp_kpt) > 0 else 0
    rec = tp_kpt / SEQC2_TOTAL
    f1_val = 2 * (prec * rec) / (prec + rec) if (prec + rec) > 0 else 0
    
    ad_sensitivity_results.append({
        "ad_threshold": round(ad_th, 2),
        "tp_removed": tp_rem,
        "fp_removed": fp_rem,
        "tp_kept": tp_kpt,
        "fp_kept": fp_kpt,
        "precision": round(prec, 4),
        "recall": round(rec, 4),
        "f1_score": round(f1_val, 4)
    })
    print(f"AD > {ad_th:.2f}: TP移除={tp_rem}, FP移除={fp_rem}, F1={f1_val:.4f}")

ad_sensitivity_df = pd.DataFrame(ad_sensitivity_results)
ad_sensitivity_df.to_csv(f"{OUTPUT_DIR}/data/ad_sensitivity_analysis.csv", index=False)

#----------------------------------------------
# 6. 生成可視化圖表
#----------------------------------------------
print("\n=== 6. 生成可視化圖表 ===")

# Figure 1: Distribution comparison
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# AlleleDelta distribution
axes[0, 0].hist(tp_df['AlleleDelta'], bins=50, alpha=0.6, label='TP', color='blue', density=True)
axes[0, 0].hist(fp_df['AlleleDelta'], bins=50, alpha=0.6, label='FP', color='red', density=True)
axes[0, 0].axvline(x=0.25, color='black', linestyle='--', label='Threshold=0.25')
axes[0, 0].set_xlabel('AlleleDelta')
axes[0, 0].set_ylabel('Density')
axes[0, 0].set_title('AlleleDelta Distribution: TP vs FP')
axes[0, 0].legend()

# VAF distribution  
axes[0, 1].hist(tp_df['VAF'], bins=50, alpha=0.6, label='TP', color='blue', density=True)
axes[0, 1].hist(fp_df['VAF'], bins=50, alpha=0.6, label='FP', color='red', density=True)
axes[0, 1].axvline(x=0.24, color='black', linestyle='--', label='Threshold=0.24')
axes[0, 1].set_xlabel('VAF')
axes[0, 1].set_ylabel('Density')
axes[0, 1].set_title('VAF Distribution: TP vs FP')
axes[0, 1].legend()

# CramersV distribution
axes[0, 2].hist(tp_df['CramersV'], bins=50, alpha=0.6, label='TP', color='blue', density=True)
axes[0, 2].hist(fp_df['CramersV'], bins=50, alpha=0.6, label='FP', color='red', density=True)
axes[0, 2].axvline(x=0.05, color='black', linestyle='--', label='Threshold=0.05')
axes[0, 2].set_xlabel('CramersV')
axes[0, 2].set_ylabel('Density')
axes[0, 2].set_title('CramersV Distribution: TP vs FP')
axes[0, 2].legend()

# VAF sensitivity
axes[1, 0].plot(sensitivity_df['vaf_threshold'], sensitivity_df['f1_score'], 'b-o', label='F1-Score')
axes[1, 0].axvline(x=0.24, color='red', linestyle='--', label='Selected=0.24')
axes[1, 0].set_xlabel('VAF Threshold')
axes[1, 0].set_ylabel('F1-Score')
axes[1, 0].set_title('VAF Threshold Sensitivity Analysis')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# AD sensitivity
axes[1, 1].plot(ad_sensitivity_df['ad_threshold'], ad_sensitivity_df['f1_score'], 'g-o', label='F1-Score')
axes[1, 1].axvline(x=0.25, color='red', linestyle='--', label='Selected=0.25')
axes[1, 1].set_xlabel('AlleleDelta Threshold')
axes[1, 1].set_ylabel('F1-Score')
axes[1, 1].set_title('AlleleDelta Threshold Sensitivity Analysis')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

# TP/FP removed by VAF threshold
ax2 = axes[1, 2]
ax2.plot(sensitivity_df['vaf_threshold'], sensitivity_df['tp_removed'], 'b-o', label='TP Removed')
ax2.plot(sensitivity_df['vaf_threshold'], sensitivity_df['fp_removed'], 'r-s', label='FP Removed')
ax2.axvline(x=0.24, color='green', linestyle='--', label='Selected=0.24')
ax2.set_xlabel('VAF Threshold')
ax2.set_ylabel('Count Removed')
ax2.set_title('Removal Count by VAF Threshold')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/distribution_and_sensitivity.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/figures/distribution_and_sensitivity.png")

# Figure 2: Scatter plot - VAF vs AlleleDelta
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# TP scatter
ax1 = axes[0]
scatter1 = ax1.scatter(tp_df['VAF'], tp_df['AlleleDelta'], c='blue', alpha=0.3, s=5, label='TP')
ax1.axhline(y=0.25, color='red', linestyle='--', linewidth=1)
ax1.axvline(x=0.24, color='red', linestyle='--', linewidth=1)
ax1.fill_between([0, 0.24], 0.25, 1.0, alpha=0.2, color='red', label='Filter Zone')
ax1.set_xlabel('VAF')
ax1.set_ylabel('AlleleDelta')
ax1.set_title(f'TP Sites (n={len(tp_df)}, removed={len(tp_removed)})')
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 0.8)
ax1.legend()

# FP scatter
ax2 = axes[1]
scatter2 = ax2.scatter(fp_df['VAF'], fp_df['AlleleDelta'], c='red', alpha=0.3, s=5, label='FP')
ax2.axhline(y=0.25, color='black', linestyle='--', linewidth=1)
ax2.axvline(x=0.24, color='black', linestyle='--', linewidth=1)
ax2.fill_between([0, 0.24], 0.25, 1.0, alpha=0.2, color='green', label='Filter Zone')
ax2.set_xlabel('VAF')
ax2.set_ylabel('AlleleDelta')
ax2.set_title(f'FP Sites (n={len(fp_df)}, removed={len(fp_removed)})')
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 0.8)
ax2.legend()

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/vaf_vs_alleledelta_scatter.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/figures/vaf_vs_alleledelta_scatter.png")

# Figure 3: Filter effectiveness comparison
fig, ax = plt.subplots(figsize=(10, 6))

categories = ['No Filter', 'Meth+VAF Filter\n(AD>0.25 & V<0.05 & VAF<0.24)']
f1_scores = [baseline_f1, f1]
colors = ['lightgray', 'steelblue']

bars = ax.bar(categories, f1_scores, color=colors, edgecolor='black')
ax.set_ylabel('F1-Score')
ax.set_title('F1-Score Comparison: Filter Effectiveness\n(Without QUAL Filter)')
ax.set_ylim(0.80, 0.85)

# Add value labels
for bar, val in zip(bars, f1_scores):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.002, 
            f'{val:.4f}', ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.axhline(y=baseline_f1, color='red', linestyle='--', alpha=0.5, label=f'Baseline={baseline_f1:.4f}')
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/f1_comparison.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/figures/f1_comparison.png")

# Figure 4: Removed sites analysis
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# TP removed characteristics
if len(tp_removed) > 0:
    axes[0, 0].hist(tp_removed['VAF'], bins=20, alpha=0.7, color='blue', edgecolor='black')
    axes[0, 0].set_xlabel('VAF')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title(f'TP Removed: VAF Distribution (n={len(tp_removed)})')
else:
    axes[0, 0].text(0.5, 0.5, 'No TP removed', ha='center', va='center', transform=axes[0, 0].transAxes)

if len(tp_removed) > 0:
    axes[0, 1].hist(tp_removed['AlleleDelta'], bins=20, alpha=0.7, color='blue', edgecolor='black')
    axes[0, 1].set_xlabel('AlleleDelta')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title(f'TP Removed: AlleleDelta Distribution')
else:
    axes[0, 1].text(0.5, 0.5, 'No TP removed', ha='center', va='center', transform=axes[0, 1].transAxes)

# FP removed characteristics
if len(fp_removed) > 0:
    axes[1, 0].hist(fp_removed['VAF'], bins=20, alpha=0.7, color='red', edgecolor='black')
    axes[1, 0].set_xlabel('VAF')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title(f'FP Removed: VAF Distribution (n={len(fp_removed)})')

if len(fp_removed) > 0:
    axes[1, 1].hist(fp_removed['AlleleDelta'], bins=20, alpha=0.7, color='red', edgecolor='black')
    axes[1, 1].set_xlabel('AlleleDelta')
    axes[1, 1].set_ylabel('Count')
    axes[1, 1].set_title(f'FP Removed: AlleleDelta Distribution')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/removed_sites_analysis.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"已儲存: {OUTPUT_DIR}/figures/removed_sites_analysis.png")

#----------------------------------------------
# 7. 生物學意義分析
#----------------------------------------------
print("\n=== 7. 生物學意義分析 ===")

# Analyze high AD + low VAF pattern (mapping error hypothesis)
print("\n假設驗證: 高 AlleleDelta + 低 VAF = Mapping 錯誤假 ASM")
print("-" * 50)

# FP with high AD + low VAF (our target)
fp_target = fp_df[(fp_df['AlleleDelta'] > 0.25) & (fp_df['VAF'] < 0.24)]
print(f"\nFP 中符合 AD>0.25 AND VAF<0.24 的位點: {len(fp_target)}")
if len(fp_target) > 0:
    print(f"  這些位點特徵:")
    print(f"  - AlleleDelta 平均: {fp_target['AlleleDelta'].mean():.4f}")
    print(f"  - VAF 平均: {fp_target['VAF'].mean():.4f}")
    print(f"  - CramersV 分布: V=0 佔 {(fp_target['CramersV'] < 0.05).sum()/len(fp_target)*100:.1f}%")
    print(f"  - NumReads 平均: {fp_target['NumReads'].mean():.1f}")

# TP with same pattern (potential misclassification)
tp_target = tp_df[(tp_df['AlleleDelta'] > 0.25) & (tp_df['VAF'] < 0.24)]
print(f"\nTP 中符合 AD>0.25 AND VAF<0.24 的位點: {len(tp_target)}")
if len(tp_target) > 0:
    print(f"  這些位點特徵:")
    print(f"  - AlleleDelta 平均: {tp_target['AlleleDelta'].mean():.4f}")
    print(f"  - VAF 平均: {tp_target['VAF'].mean():.4f}")
    print(f"  - CramersV 分布: V=0 佔 {(tp_target['CramersV'] < 0.05).sum()/len(tp_target)*100:.1f}%")
    print(f"  - QUAL 平均: {tp_target['QUAL'].mean():.4f}")

# High VAF + low AlleleDelta (good quality TP)
tp_highvaf = tp_df[(tp_df['VAF'] > 0.4) & (tp_df['AlleleDelta'] < 0.15)]
print(f"\n高品質 TP (VAF>0.4 AND AD<0.15): {len(tp_highvaf)}")

# Summary statistics
analysis_summary = {
    "hypothesis": "高 AlleleDelta + 低 VAF = Mapping 錯誤導致的假 ASM",
    "fp_target_count": len(fp_target),
    "fp_target_pct": f"{100*len(fp_target)/len(fp_df):.2f}%",
    "tp_target_count": len(tp_target),
    "tp_target_pct": f"{100*len(tp_target)/len(tp_df):.4f}%",
    "filter_specificity": f"{100*len(fp_target)/(len(fp_target)+len(tp_target)):.2f}%" if (len(fp_target)+len(tp_target)) > 0 else "N/A",
    "biological_interpretation": {
        "high_ad_low_vaf_fp": "這些 FP 很可能是 mapping 錯誤造成的，少量錯誤 reads 被貼到錯誤位置，帶來不同的甲基化模式",
        "high_ad_low_vaf_tp": "極少數 TP 也符合此模式，可能是真實的低頻亞克隆突變",
        "purity_consideration": "在 purity=1.0 的情況下，VAF<0.24 且 Ref>Alt 很可能代表 Alt reads 是錯誤貼過來的"
    }
}

with open(f"{OUTPUT_DIR}/data/biological_analysis.json", "w") as f:
    json.dump(analysis_summary, f, indent=2, ensure_ascii=False)

print("\n=== 分析完成 ===")
print(f"所有結果已儲存至: {OUTPUT_DIR}")
print(f"圖表: {OUTPUT_DIR}/figures/")
print(f"數據: {OUTPUT_DIR}/data/")
