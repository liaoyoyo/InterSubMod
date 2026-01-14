#!/usr/bin/env python3
"""
Phase 1: 快速驗證實驗 - 甲基化 F1 研究計劃

包含五個實驗:
1.1 VCF QUAL 分布分析
1.2 VCF AF 分布分析
1.3 Strand Bias 分析
1.4 Regional Clustering 分析
1.5 現有特徵組合 AUC

輸入:
- TP/FP VCF: data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz, filtered_snv_fp.vcf.gz
- TP/FP 分析結果: output/bip8_disk_output/20260113_all-with-w5000_3/

輸出:
- 圖表: docs/research/methylation_f1_optimization/05_results/phase1_plots/
- 報告: docs/research/methylation_f1_optimization/05_results/phase1_report.md
"""

import subprocess
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set Chinese font for matplotlib
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei', 'Taipei Sans TC Beta']
plt.rcParams['axes.unicode_minus'] = False

# --- Paths Configuration ---
PROJECT_ROOT = Path("/big8_disk/liaoyoyo2001/InterSubMod")
VCF_DIR = PROJECT_ROOT / "data/vcf/HCC1395/pileup"
TP_VCF = VCF_DIR / "filtered_snv_tp.vcf.gz"
FP_VCF = VCF_DIR / "filtered_snv_fp.vcf.gz"

# Analysis output directory
ANALYSIS_DIR = PROJECT_ROOT / "output/bip8_disk_output/20260113_all-with-w5000_3"
TP_SUMMARY = ANALYSIS_DIR / "filtered_snv_tp/significance_summary.csv"
FP_SUMMARY = ANALYSIS_DIR / "filtered_snv_fp/significance_summary.csv"

# Output paths
RESULTS_DIR = PROJECT_ROOT / "docs/research/methylation_f1_optimization/05_results"
PLOT_DIR = RESULTS_DIR / "phase1_plots"
REPORT_PATH = RESULTS_DIR / "phase1_report.md"


def setup_directories():
    """Create output directories"""
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output directory: {PLOT_DIR}")


def extract_vcf_qual(vcf_path: Path) -> pd.DataFrame:
    """Extract QUAL from VCF using bcftools"""
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t%QUAL\\n' {vcf_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    lines = result.stdout.strip().split('\n')
    data = []
    for line in lines:
        if line:
            parts = line.split('\t')
            if len(parts) >= 3:
                try:
                    qual = float(parts[2]) if parts[2] != '.' else np.nan
                    data.append({'Chr': parts[0], 'Pos': int(parts[1]), 'QUAL': qual})
                except ValueError:
                    continue
    return pd.DataFrame(data)


def extract_vcf_af(vcf_path: Path) -> pd.DataFrame:
    """Extract AF (Allele Frequency) from VCF"""
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t[%AF]\\n' {vcf_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    lines = result.stdout.strip().split('\n')
    data = []
    for line in lines:
        if line:
            parts = line.split('\t')
            if len(parts) >= 3:
                try:
                    af = float(parts[2]) if parts[2] != '.' else np.nan
                    data.append({'Chr': parts[0], 'Pos': int(parts[1]), 'AF': af})
                except ValueError:
                    continue
    return pd.DataFrame(data)


def extract_vcf_strand_bias(vcf_path: Path) -> pd.DataFrame:
    """Extract strand counts (FAU, RAU, etc.) from VCF for strand bias calculation"""
    # Extract forward and reverse strand counts
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t[%FAU]\\t[%RAU]\\n' {vcf_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    lines = result.stdout.strip().split('\n')
    data = []
    for line in lines:
        if line:
            parts = line.split('\t')
            if len(parts) >= 4:
                try:
                    # FAU and RAU are comma-separated values for each allele
                    fau_vals = [int(x) for x in parts[2].split(',') if x.isdigit()]
                    rau_vals = [int(x) for x in parts[3].split(',') if x.isdigit()]
                    
                    fau_sum = sum(fau_vals) if fau_vals else 0
                    rau_sum = sum(rau_vals) if rau_vals else 0
                    
                    # Calculate strand bias ratio
                    total = fau_sum + rau_sum
                    if total > 0:
                        strand_ratio = fau_sum / total
                        strand_bias = abs(strand_ratio - 0.5) * 2  # 0 = balanced, 1 = complete bias
                    else:
                        strand_bias = np.nan
                    
                    data.append({
                        'Chr': parts[0], 
                        'Pos': int(parts[1]), 
                        'FAU': fau_sum,
                        'RAU': rau_sum,
                        'StrandBias': strand_bias
                    })
                except (ValueError, IndexError):
                    continue
    return pd.DataFrame(data)


def calculate_regional_clustering(df: pd.DataFrame, window_size: int = 50000) -> pd.DataFrame:
    """Calculate regional clustering - number of variants within window"""
    result = df.copy()
    result['ClusterCount'] = 0
    
    for chrom in df['Chr'].unique():
        chrom_mask = result['Chr'] == chrom
        chrom_df = result[chrom_mask].sort_values('Pos')
        
        for idx, row in chrom_df.iterrows():
            pos = row['Pos']
            # Count variants within window
            count = ((chrom_df['Pos'] >= pos - window_size) & 
                     (chrom_df['Pos'] <= pos + window_size)).sum() - 1
            result.loc[idx, 'ClusterCount'] = count
    
    return result


def calculate_auc(y_true, y_score, feature_name: str):
    """Calculate AUC and return statistics"""
    # Filter out NaN values
    mask = ~np.isnan(y_score)
    y_true_clean = y_true[mask]
    y_score_clean = y_score[mask]
    
    if len(y_true_clean) < 10 or len(np.unique(y_true_clean)) < 2:
        return None
    
    try:
        auc = roc_auc_score(y_true_clean, y_score_clean)
        return {
            'Feature': feature_name,
            'AUC': auc,
            'N': len(y_true_clean),
            'N_TP': sum(y_true_clean),
            'N_FP': len(y_true_clean) - sum(y_true_clean)
        }
    except Exception as e:
        print(f"[WARN] Cannot calculate AUC for {feature_name}: {e}")
        return None


def experiment_1_qual_analysis():
    """Experiment 1.1: VCF QUAL Distribution Analysis"""
    print("\n" + "="*60)
    print("實驗 1.1: VCF QUAL 分布分析")
    print("="*60)
    
    # Extract QUAL from VCFs
    print("[INFO] 提取 VCF QUAL 值...")
    tp_qual = extract_vcf_qual(TP_VCF)
    fp_qual = extract_vcf_qual(FP_VCF)
    
    tp_qual['Label'] = 1
    fp_qual['Label'] = 0
    
    combined = pd.concat([tp_qual, fp_qual], ignore_index=True)
    combined = combined.dropna(subset=['QUAL'])
    
    print(f"[INFO] TP: {len(tp_qual)} sites, FP: {len(fp_qual)} sites")
    print(f"[INFO] TP QUAL: mean={tp_qual['QUAL'].mean():.2f}, median={tp_qual['QUAL'].median():.2f}")
    print(f"[INFO] FP QUAL: mean={fp_qual['QUAL'].mean():.2f}, median={fp_qual['QUAL'].median():.2f}")
    
    # Calculate AUC
    auc_result = calculate_auc(combined['Label'].values, combined['QUAL'].values, 'QUAL')
    if auc_result:
        print(f"[RESULT] QUAL AUC: {auc_result['AUC']:.4f}")
    
    # Plot distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    axes[0].hist(tp_qual['QUAL'].dropna(), bins=50, alpha=0.6, label='TP', density=True, color='green')
    axes[0].hist(fp_qual['QUAL'].dropna(), bins=50, alpha=0.6, label='FP', density=True, color='red')
    axes[0].set_xlabel('QUAL Score')
    axes[0].set_ylabel('Density')
    axes[0].set_title('QUAL Distribution: TP vs FP')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Box plot
    combined['Type'] = combined['Label'].map({1: 'TP', 0: 'FP'})
    combined.boxplot(column='QUAL', by='Type', ax=axes[1])
    axes[1].set_xlabel('Type')
    axes[1].set_ylabel('QUAL Score')
    axes[1].set_title(f'QUAL Score by Type (AUC={auc_result["AUC"]:.4f})' if auc_result else 'QUAL Score by Type')
    
    plt.suptitle('')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'exp1_1_qual_distribution.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Mann-Whitney U test
    stat, pval = stats.mannwhitneyu(
        tp_qual['QUAL'].dropna(), 
        fp_qual['QUAL'].dropna(), 
        alternative='greater'
    )
    print(f"[STATS] Mann-Whitney U test: stat={stat:.2f}, p-value={pval:.4e}")
    
    return {
        'experiment': '1.1 QUAL Analysis',
        'auc': auc_result['AUC'] if auc_result else None,
        'tp_mean': tp_qual['QUAL'].mean(),
        'fp_mean': fp_qual['QUAL'].mean(),
        'pvalue': pval,
        'conclusion': 'Effective' if auc_result and auc_result['AUC'] > 0.55 else 'Weak/Ineffective'
    }


def experiment_2_af_analysis():
    """Experiment 1.2: VCF AF Distribution Analysis"""
    print("\n" + "="*60)
    print("實驗 1.2: VCF AF 分布分析")
    print("="*60)
    
    print("[INFO] 提取 VCF AF 值...")
    tp_af = extract_vcf_af(TP_VCF)
    fp_af = extract_vcf_af(FP_VCF)
    
    tp_af['Label'] = 1
    fp_af['Label'] = 0
    
    combined = pd.concat([tp_af, fp_af], ignore_index=True)
    combined = combined.dropna(subset=['AF'])
    
    print(f"[INFO] TP AF: mean={tp_af['AF'].mean():.4f}, median={tp_af['AF'].median():.4f}")
    print(f"[INFO] FP AF: mean={fp_af['AF'].mean():.4f}, median={fp_af['AF'].median():.4f}")
    
    # Calculate AUC
    auc_result = calculate_auc(combined['Label'].values, combined['AF'].values, 'AF')
    if auc_result:
        print(f"[RESULT] AF AUC: {auc_result['AUC']:.4f}")
    
    # Plot distribution
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Histogram
    axes[0].hist(tp_af['AF'].dropna(), bins=50, alpha=0.6, label='TP', density=True, color='green')
    axes[0].hist(fp_af['AF'].dropna(), bins=50, alpha=0.6, label='FP', density=True, color='red')
    axes[0].set_xlabel('Allele Frequency')
    axes[0].set_ylabel('Density')
    axes[0].set_title('AF Distribution: TP vs FP')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # KDE plot
    if len(tp_af['AF'].dropna()) > 10 and len(fp_af['AF'].dropna()) > 10:
        sns.kdeplot(data=tp_af['AF'].dropna(), ax=axes[1], label='TP', color='green')
        sns.kdeplot(data=fp_af['AF'].dropna(), ax=axes[1], label='FP', color='red')
    axes[1].set_xlabel('Allele Frequency')
    axes[1].set_ylabel('Density')
    axes[1].set_title('AF KDE: TP vs FP')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # AF range analysis
    af_ranges = [(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 1.0)]
    tp_ratios = []
    labels = []
    
    for low, high in af_ranges:
        tp_in_range = ((tp_af['AF'] >= low) & (tp_af['AF'] < high)).sum()
        fp_in_range = ((fp_af['AF'] >= low) & (fp_af['AF'] < high)).sum()
        total = tp_in_range + fp_in_range
        if total > 0:
            tp_ratio = tp_in_range / total
        else:
            tp_ratio = 0
        tp_ratios.append(tp_ratio)
        labels.append(f'{low:.1f}-{high:.1f}')
    
    axes[2].bar(labels, tp_ratios, color='skyblue', edgecolor='navy')
    axes[2].axhline(y=0.5, color='red', linestyle='--', label='50% baseline')
    axes[2].set_xlabel('AF Range')
    axes[2].set_ylabel('TP Ratio')
    axes[2].set_title('TP Ratio by AF Range')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'exp1_2_af_distribution.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Stats
    stat, pval = stats.mannwhitneyu(
        tp_af['AF'].dropna(), 
        fp_af['AF'].dropna()
    )
    print(f"[STATS] Mann-Whitney U test: stat={stat:.2f}, p-value={pval:.4e}")
    
    return {
        'experiment': '1.2 AF Analysis',
        'auc': auc_result['AUC'] if auc_result else None,
        'tp_mean': tp_af['AF'].mean(),
        'fp_mean': fp_af['AF'].mean(),
        'pvalue': pval,
        'conclusion': 'Effective' if auc_result and auc_result['AUC'] > 0.55 else 'Weak/Ineffective'
    }


def experiment_3_strand_bias():
    """Experiment 1.3: Strand Bias Analysis"""
    print("\n" + "="*60)
    print("實驗 1.3: Strand Bias 分析")
    print("="*60)
    
    print("[INFO] 提取 Strand 資訊...")
    tp_strand = extract_vcf_strand_bias(TP_VCF)
    fp_strand = extract_vcf_strand_bias(FP_VCF)
    
    tp_strand['Label'] = 1
    fp_strand['Label'] = 0
    
    # Check if we have valid data
    if tp_strand.empty or fp_strand.empty or 'StrandBias' not in tp_strand.columns:
        print("[WARN] FAU/RAU 欄位在 VCF 中不存在或無法解析，跳過 Strand Bias 分析")
        return {
            'experiment': '1.3 Strand Bias',
            'auc': None,
            'conclusion': 'SkippedNoData'
        }
    
    combined = pd.concat([tp_strand, fp_strand], ignore_index=True)
    combined = combined.dropna(subset=['StrandBias'])
    
    if len(combined) < 10:
        print("[WARN] 有效樣本數過少，跳過分析")
        return {
            'experiment': '1.3 Strand Bias',
            'auc': None,
            'conclusion': 'SkippedInsufficientData'
        }
    
    print(f"[INFO] TP StrandBias: mean={tp_strand['StrandBias'].mean():.4f}")
    print(f"[INFO] FP StrandBias: mean={fp_strand['StrandBias'].mean():.4f}")
    
    # Calculate AUC
    auc_result = calculate_auc(combined['Label'].values, combined['StrandBias'].values, 'StrandBias')
    if auc_result:
        print(f"[RESULT] StrandBias AUC: {auc_result['AUC']:.4f}")
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    combined['Type'] = combined['Label'].map({1: 'TP', 0: 'FP'})
    combined.boxplot(column='StrandBias', by='Type', ax=ax)
    ax.set_xlabel('Type')
    ax.set_ylabel('Strand Bias (0=balanced, 1=biased)')
    ax.set_title(f'Strand Bias by Type (AUC={auc_result["AUC"]:.4f})' if auc_result else 'Strand Bias by Type')
    plt.suptitle('')
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'exp1_3_strand_bias.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return {
        'experiment': '1.3 Strand Bias',
        'auc': auc_result['AUC'] if auc_result else None,
        'conclusion': 'Effective' if auc_result and auc_result['AUC'] > 0.55 else 'Weak/Ineffective'
    }


def experiment_4_regional_clustering():
    """Experiment 1.4: Regional Clustering Analysis"""
    print("\n" + "="*60)
    print("實驗 1.4: Regional Clustering 分析 (50kb window)")
    print("="*60)
    
    print("[INFO] 計算區域聚集指數...")
    
    # Load positions from VCF
    tp_pos = extract_vcf_qual(TP_VCF)[['Chr', 'Pos']]
    fp_pos = extract_vcf_qual(FP_VCF)[['Chr', 'Pos']]
    
    tp_pos['Label'] = 1
    fp_pos['Label'] = 0
    
    combined = pd.concat([tp_pos, fp_pos], ignore_index=True)
    combined = calculate_regional_clustering(combined, window_size=50000)
    
    tp_cluster = combined[combined['Label'] == 1]['ClusterCount']
    fp_cluster = combined[combined['Label'] == 0]['ClusterCount']
    
    print(f"[INFO] TP ClusterCount: mean={tp_cluster.mean():.2f}, median={tp_cluster.median():.0f}")
    print(f"[INFO] FP ClusterCount: mean={fp_cluster.mean():.2f}, median={fp_cluster.median():.0f}")
    
    # Calculate AUC (inverted - higher cluster might indicate FP)
    auc_result = calculate_auc(combined['Label'].values, -combined['ClusterCount'].values, 'ClusterCount')
    if auc_result:
        auc_result['AUC'] = 1 - auc_result['AUC']  # Invert back
        print(f"[RESULT] ClusterCount AUC: {auc_result['AUC']:.4f} (lower = more likely FP)")
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    axes[0].hist(tp_cluster, bins=30, alpha=0.6, label='TP', density=True, color='green')
    axes[0].hist(fp_cluster, bins=30, alpha=0.6, label='FP', density=True, color='red')
    axes[0].set_xlabel('Cluster Count (50kb window)')
    axes[0].set_ylabel('Density')
    axes[0].set_title('Regional Clustering Distribution')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # TP ratio by cluster count
    cluster_bins = [0, 5, 10, 20, 50, 100, combined['ClusterCount'].max()]
    tp_ratios = []
    labels = []
    
    for i in range(len(cluster_bins)-1):
        low, high = cluster_bins[i], cluster_bins[i+1]
        tp_in_range = ((combined['Label'] == 1) & 
                       (combined['ClusterCount'] >= low) & 
                       (combined['ClusterCount'] < high)).sum()
        fp_in_range = ((combined['Label'] == 0) & 
                       (combined['ClusterCount'] >= low) & 
                       (combined['ClusterCount'] < high)).sum()
        total = tp_in_range + fp_in_range
        if total > 0:
            tp_ratio = tp_in_range / total
        else:
            tp_ratio = 0
        tp_ratios.append(tp_ratio)
        labels.append(f'{low}-{high}')
    
    axes[1].bar(labels, tp_ratios, color='skyblue', edgecolor='navy')
    axes[1].axhline(y=0.5, color='red', linestyle='--', label='50% baseline')
    axes[1].set_xlabel('Cluster Count Range')
    axes[1].set_ylabel('TP Ratio')
    axes[1].set_title('TP Ratio by Clustering Level')
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'exp1_4_regional_clustering.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return {
        'experiment': '1.4 Regional Clustering',
        'auc': auc_result['AUC'] if auc_result else None,
        'tp_mean': tp_cluster.mean(),
        'fp_mean': fp_cluster.mean(),
        'conclusion': 'Effective' if auc_result and abs(auc_result['AUC'] - 0.5) > 0.05 else 'Weak/Ineffective'
    }


def experiment_5_combined_features():
    """Experiment 1.5: Combined Feature AUC with Random Forest"""
    print("\n" + "="*60)
    print("實驗 1.5: 現有特徵組合 AUC (Random Forest)")
    print("="*60)
    
    print("[INFO] 載入甲基化分析結果...")
    tp_summary = pd.read_csv(TP_SUMMARY)
    fp_summary = pd.read_csv(FP_SUMMARY)
    
    tp_summary['Label'] = 1
    fp_summary['Label'] = 0
    
    combined = pd.concat([tp_summary, fp_summary], ignore_index=True)
    
    # Select features for analysis
    feature_cols = ['GlobalP', 'CramersV', 'HeuristicScore', 'NumReads', 'NumCpGs']
    available_features = [c for c in feature_cols if c in combined.columns]
    
    print(f"[INFO] Available features: {available_features}")
    print(f"[INFO] TP samples: {len(tp_summary)}, FP samples: {len(fp_summary)}")
    
    # Clean data
    combined_clean = combined[available_features + ['Label']].dropna()
    X = combined_clean[available_features].values
    y = combined_clean['Label'].values
    
    print(f"[INFO] Clean samples: {len(combined_clean)} (after removing NaN)")
    
    # Individual feature AUC
    individual_aucs = {}
    for feat in available_features:
        auc_result = calculate_auc(combined_clean['Label'].values, combined_clean[feat].values, feat)
        if auc_result:
            individual_aucs[feat] = auc_result['AUC']
            print(f"[RESULT] {feat} AUC: {auc_result['AUC']:.4f}")
    
    # Random Forest with cross-validation
    print("\n[INFO] 訓練 Random Forest 模型...")
    rf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42, n_jobs=-1)
    cv_scores = cross_val_score(rf, X, y, cv=5, scoring='roc_auc')
    
    print(f"[RESULT] 5-Fold CV AUC: {cv_scores.mean():.4f} (+/- {cv_scores.std():.4f})")
    
    # Feature importance
    rf.fit(X, y)
    feature_importance = dict(zip(available_features, rf.feature_importances_))
    
    print("\n[INFO] Feature Importance:")
    for feat, imp in sorted(feature_importance.items(), key=lambda x: x[1], reverse=True):
        print(f"  {feat}: {imp:.4f}")
    
    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Individual AUC
    features = list(individual_aucs.keys())
    aucs = list(individual_aucs.values())
    colors = ['green' if a > 0.55 else 'orange' if a > 0.52 else 'red' for a in aucs]
    axes[0].barh(features, aucs, color=colors, edgecolor='black')
    axes[0].axvline(x=0.5, color='gray', linestyle='--', label='No discrimination')
    axes[0].axvline(x=0.55, color='orange', linestyle='--', label='Weak threshold')
    axes[0].set_xlabel('AUC')
    axes[0].set_title('Individual Feature AUC')
    axes[0].legend(loc='lower right')
    axes[0].set_xlim(0.4, 0.7)
    
    # Feature importance
    axes[1].barh(available_features, [feature_importance[f] for f in available_features], 
                 color='steelblue', edgecolor='black')
    axes[1].set_xlabel('Importance')
    axes[1].set_title('Random Forest Feature Importance')
    
    # CV scores
    axes[2].bar(range(1, 6), cv_scores, color='lightgreen', edgecolor='green')
    axes[2].axhline(y=cv_scores.mean(), color='red', linestyle='--', label=f'Mean: {cv_scores.mean():.4f}')
    axes[2].set_xlabel('Fold')
    axes[2].set_ylabel('AUC')
    axes[2].set_title('5-Fold Cross-Validation AUC')
    axes[2].legend()
    axes[2].set_ylim(0.4, 0.7)
    
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'exp1_5_combined_features.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return {
        'experiment': '1.5 Combined Features',
        'cv_auc_mean': cv_scores.mean(),
        'cv_auc_std': cv_scores.std(),
        'individual_aucs': individual_aucs,
        'feature_importance': feature_importance,
        'conclusion': 'Effective' if cv_scores.mean() > 0.6 else 'Weak' if cv_scores.mean() > 0.55 else 'Ineffective'
    }


def generate_report(results: list):
    """Generate final report in Markdown"""
    print("\n" + "="*60)
    print("生成分析報告...")
    print("="*60)
    
    report = """# Phase 1: 快速驗證實驗報告

生成時間: 2026-01-13

## 實驗總覽

| 實驗 | AUC | 結論 |
|:---|:---:|:---|
"""
    
    for r in results:
        exp_name = r['experiment']
        auc = r.get('auc') or r.get('cv_auc_mean')
        conclusion = r['conclusion']
        
        if auc is not None:
            report += f"| {exp_name} | {auc:.4f} | {conclusion} |\n"
        else:
            report += f"| {exp_name} | N/A | {conclusion} |\n"
    
    report += """
## 判定標準

- **有效 (Effective)**: AUC > 0.6
- **弱有效 (Weak)**: 0.55 < AUC < 0.6
- **無效 (Ineffective)**: AUC < 0.55

## 詳細結果

"""
    
    for i, r in enumerate(results, 1):
        report += f"### {r['experiment']}\n\n"
        
        for key, val in r.items():
            if key not in ['experiment', 'individual_aucs', 'feature_importance']:
                if isinstance(val, float):
                    report += f"- **{key}**: {val:.4f}\n"
                else:
                    report += f"- **{key}**: {val}\n"
        
        if 'individual_aucs' in r:
            report += "\n**Individual Feature AUC:**\n"
            for feat, auc in r['individual_aucs'].items():
                report += f"- {feat}: {auc:.4f}\n"
        
        if 'feature_importance' in r:
            report += "\n**Feature Importance:**\n"
            for feat, imp in sorted(r['feature_importance'].items(), key=lambda x: x[1], reverse=True):
                report += f"- {feat}: {imp:.4f}\n"
        
        report += f"\n![{r['experiment']}](phase1_plots/exp1_{i}_*.png)\n\n"
    
    report += """
## 結論

根據 Phase 1 實驗結果，評估各特徵對 TP/FP 區分的有效性。

### 下一步建議

1. 對有效特徵進行深度分析
2. 嘗試特徵組合過濾策略
3. 評估對 F1-score 的實際影響
"""
    
    with open(REPORT_PATH, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"[INFO] 報告已儲存至: {REPORT_PATH}")


def main():
    print("="*60)
    print(" Phase 1: 快速驗證實驗")
    print(" 甲基化 F1 研究計劃")
    print("="*60)
    
    setup_directories()
    
    results = []
    
    # Run all experiments
    results.append(experiment_1_qual_analysis())
    results.append(experiment_2_af_analysis())
    results.append(experiment_3_strand_bias())
    results.append(experiment_4_regional_clustering())
    results.append(experiment_5_combined_features())
    
    # Generate report
    generate_report(results)
    
    print("\n" + "="*60)
    print(" 所有實驗完成！")
    print(f" 圖表目錄: {PLOT_DIR}")
    print(f" 報告: {REPORT_PATH}")
    print("="*60)
    
    return results


if __name__ == "__main__":
    main()
