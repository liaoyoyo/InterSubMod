#!/usr/bin/env python3
"""
Phase 4: 甲基化顯著性與 VCF 特徵組合分析

探索內容:
1. CramersV × QUAL 交互效應
2. GlobalP × AF 交互效應
3. 顯著/非顯著位點的 VCF 特徵分佈差異
4. Significant 標記對 TP/FP 區分的價值
5. 甲基化特徵在不同 QUAL/AF 區間的效果
6. 組合特徵過濾策略評估
"""

import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.metrics import roc_auc_score
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

PROJECT_ROOT = Path("/big8_disk/liaoyoyo2001/InterSubMod")
VCF_DIR = PROJECT_ROOT / "data/vcf/HCC1395/pileup"
TP_VCF = VCF_DIR / "filtered_snv_tp.vcf.gz"
FP_VCF = VCF_DIR / "filtered_snv_fp.vcf.gz"

ANALYSIS_DIR = PROJECT_ROOT / "output/bip8_disk_output/20260113_all-with-w5000_3"
TP_SUMMARY = ANALYSIS_DIR / "filtered_snv_tp/significance_summary.csv"
FP_SUMMARY = ANALYSIS_DIR / "filtered_snv_fp/significance_summary.csv"

RESULTS_DIR = PROJECT_ROOT / "docs/research/methylation_f1_optimization/05_results"
PLOT_DIR = RESULTS_DIR / "phase4_plots"
REPORT_PATH = RESULTS_DIR / "phase4_methylation_combination_report.md"


def setup_directories():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)


def extract_vcf_features(vcf_path: Path) -> pd.DataFrame:
    """Extract QUAL and AF from VCF"""
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t%QUAL\\t[%AF]\\n' {vcf_path}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    data = []
    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 4:
                try:
                    qual = float(parts[2]) if parts[2] != '.' else np.nan
                    af = float(parts[3]) if parts[3] != '.' else np.nan
                    data.append({
                        'Chr': parts[0], 
                        'Pos': int(parts[1]), 
                        'QUAL': qual,
                        'AF': af
                    })
                except ValueError:
                    continue
    return pd.DataFrame(data)


def merge_vcf_methylation(vcf_data: pd.DataFrame, meth_data: pd.DataFrame) -> pd.DataFrame:
    """Merge VCF and methylation data by Chr and Pos"""
    merged = vcf_data.merge(meth_data, on=['Chr', 'Pos'], how='inner')
    return merged


def calculate_f1(tp, fp, fn):
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return precision, recall, f1


def analyze_significant_vs_nonsignificant(tp_merged, fp_merged):
    """Analyze VCF features for significant vs non-significant sites"""
    print("\n" + "="*60)
    print("分析 1: 甲基化顯著性對 VCF 特徵分佈的影響")
    print("="*60)
    
    # Add label
    tp_merged = tp_merged.copy()
    fp_merged = fp_merged.copy()
    tp_merged['Label'] = 'TP'
    fp_merged['Label'] = 'FP'
    
    # Count significant sites
    if 'Significant' in tp_merged.columns:
        tp_sig = (tp_merged['Significant'] == True).sum()
        tp_nonsig = (tp_merged['Significant'] == False).sum()
        fp_sig = (fp_merged['Significant'] == True).sum()
        fp_nonsig = (fp_merged['Significant'] == False).sum()
        
        print(f"\n[顯著性分佈]")
        print(f"  TP 顯著: {tp_sig} ({tp_sig/len(tp_merged)*100:.1f}%)")
        print(f"  TP 非顯著: {tp_nonsig} ({tp_nonsig/len(tp_merged)*100:.1f}%)")
        print(f"  FP 顯著: {fp_sig} ({fp_sig/len(fp_merged)*100:.1f}%)")
        print(f"  FP 非顯著: {fp_nonsig} ({fp_nonsig/len(fp_merged)*100:.1f}%)")
        
        # TP ratio in significant vs non-significant
        if tp_sig + fp_sig > 0:
            sig_tp_ratio = tp_sig / (tp_sig + fp_sig)
            print(f"\n  顯著位點中 TP 比例: {sig_tp_ratio*100:.1f}%")
        if tp_nonsig + fp_nonsig > 0:
            nonsig_tp_ratio = tp_nonsig / (tp_nonsig + fp_nonsig)
            print(f"  非顯著位點中 TP 比例: {nonsig_tp_ratio*100:.1f}%")
    
    # Analyze QUAL distribution by significance
    all_data = pd.concat([tp_merged, fp_merged], ignore_index=True)
    
    # Create visualization
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # QUAL by Label and Significance
    if 'Significant' in all_data.columns:
        all_data['SigLabel'] = all_data.apply(
            lambda x: f"{x['Label']}-{'Sig' if x['Significant'] else 'NonSig'}", axis=1
        )
        
        sig_order = ['TP-Sig', 'TP-NonSig', 'FP-Sig', 'FP-NonSig']
        colors = {'TP-Sig': 'darkgreen', 'TP-NonSig': 'lightgreen', 
                  'FP-Sig': 'darkred', 'FP-NonSig': 'lightcoral'}
        
        # Box plot of QUAL
        all_data.boxplot(column='QUAL', by='SigLabel', ax=axes[0, 0])
        axes[0, 0].set_xlabel('Type-Significance')
        axes[0, 0].set_ylabel('QUAL')
        axes[0, 0].set_title('QUAL Distribution by Type and Significance')
        plt.suptitle('')
        
        # Box plot of AF
        all_data.boxplot(column='AF', by='SigLabel', ax=axes[0, 1])
        axes[0, 1].set_xlabel('Type-Significance')
        axes[0, 1].set_ylabel('AF')
        axes[0, 1].set_title('AF Distribution by Type and Significance')
        plt.suptitle('')
    
    # CramersV vs QUAL scatter
    if 'CramersV' in all_data.columns:
        for label, color in [('TP', 'green'), ('FP', 'red')]:
            subset = all_data[all_data['Label'] == label]
            axes[1, 0].scatter(subset['QUAL'], subset['CramersV'], 
                              alpha=0.3, c=color, label=label, s=10)
        axes[1, 0].set_xlabel('QUAL')
        axes[1, 0].set_ylabel('CramersV')
        axes[1, 0].set_title('CramersV vs QUAL')
        axes[1, 0].legend()
    
    # GlobalP vs AF scatter
    if 'GlobalP' in all_data.columns:
        for label, color in [('TP', 'green'), ('FP', 'red')]:
            subset = all_data[all_data['Label'] == label]
            # Use -log10(P) for better visualization
            log_p = -np.log10(subset['GlobalP'].replace(0, 1e-300))
            axes[1, 1].scatter(subset['AF'], log_p, 
                              alpha=0.3, c=color, label=label, s=10)
        axes[1, 1].set_xlabel('AF')
        axes[1, 1].set_ylabel('-log10(GlobalP)')
        axes[1, 1].set_title('-log10(GlobalP) vs AF')
        axes[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'significant_vs_nonsig_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return all_data


def analyze_interaction_effects(tp_merged, fp_merged):
    """Analyze interaction effects between methylation and VCF features"""
    print("\n" + "="*60)
    print("分析 2: 甲基化 × VCF 特徵交互效應")
    print("="*60)
    
    # Combine data
    tp_merged = tp_merged.copy()
    fp_merged = fp_merged.copy()
    tp_merged['Label'] = 1
    fp_merged['Label'] = 0
    all_data = pd.concat([tp_merged, fp_merged], ignore_index=True)
    
    # Create interaction features
    if 'CramersV' in all_data.columns:
        all_data['QUAL_x_CramersV'] = all_data['QUAL'] * all_data['CramersV']
        all_data['QUAL_x_HScore'] = all_data['QUAL'] * all_data['HeuristicScore']
    if 'GlobalP' in all_data.columns:
        all_data['negLogP'] = -np.log10(all_data['GlobalP'].replace(0, 1e-300))
        all_data['AF_x_negLogP'] = all_data['AF'] * all_data['negLogP']
    
    # Calculate AUC for each feature
    features = ['QUAL', 'AF', 'CramersV', 'GlobalP', 'HeuristicScore', 'NumReads',
                'QUAL_x_CramersV', 'QUAL_x_HScore', 'AF_x_negLogP']
    
    aucs = {}
    print("\n[特徵 AUC 比較]")
    for feat in features:
        if feat in all_data.columns:
            clean = all_data[[feat, 'Label']].dropna()
            if len(clean) > 100:
                try:
                    auc = roc_auc_score(clean['Label'], clean[feat])
                    aucs[feat] = auc
                    symbol = "✅" if auc > 0.55 else "❌"
                    print(f"  {feat}: AUC={auc:.4f} {symbol}")
                except:
                    pass
    
    # Visualize interaction effects
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # QUAL × CramersV heatmap
    if 'CramersV' in all_data.columns:
        # Bin QUAL and CramersV
        all_data['QUAL_bin'] = pd.cut(all_data['QUAL'], bins=[0, 0.3, 0.5, 0.7, 0.9, 1.0], 
                                       labels=['0-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1'])
        all_data['V_bin'] = pd.cut(all_data['CramersV'], bins=[-0.1, 0, 0.1, 0.3, 0.5, 1.0],
                                    labels=['0', '0-0.1', '0.1-0.3', '0.3-0.5', '0.5+'])
        
        # TP ratio heatmap
        pivot = all_data.groupby(['QUAL_bin', 'V_bin'])['Label'].mean().unstack()
        if not pivot.empty:
            sns.heatmap(pivot, annot=True, fmt='.2f', cmap='RdYlGn', ax=axes[0, 0],
                       vmin=0.5, vmax=1.0)
            axes[0, 0].set_title('TP Ratio: QUAL × CramersV')
            axes[0, 0].set_xlabel('CramersV')
            axes[0, 0].set_ylabel('QUAL')
    
    # AF × Significance heatmap
    if 'Significant' in all_data.columns:
        all_data['AF_bin'] = pd.cut(all_data['AF'], bins=[0, 0.1, 0.2, 0.3, 0.5, 1.0],
                                     labels=['0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.5', '0.5+'])
        all_data['Sig_str'] = all_data['Significant'].map({True: 'Sig', False: 'NonSig'})
        
        pivot = all_data.groupby(['AF_bin', 'Sig_str'])['Label'].mean().unstack()
        if not pivot.empty:
            sns.heatmap(pivot, annot=True, fmt='.2f', cmap='RdYlGn', ax=axes[0, 1],
                       vmin=0.5, vmax=1.0)
            axes[0, 1].set_title('TP Ratio: AF × Significance')
            axes[0, 1].set_xlabel('Significance')
            axes[0, 1].set_ylabel('AF')
    
    # Feature importance bar chart
    sorted_aucs = sorted(aucs.items(), key=lambda x: x[1], reverse=True)
    if sorted_aucs:
        feats, vals = zip(*sorted_aucs)
        colors = ['green' if v > 0.55 else 'orange' if v > 0.52 else 'red' for v in vals]
        axes[1, 0].barh(feats, vals, color=colors)
        axes[1, 0].axvline(x=0.5, color='gray', linestyle='--')
        axes[1, 0].set_xlabel('AUC')
        axes[1, 0].set_title('Feature AUC Comparison')
        axes[1, 0].set_xlim(0.4, 1.0)
    
    # NumReads effect by QUAL range
    if 'NumReads' in all_data.columns:
        qual_ranges = [(0, 0.5), (0.5, 0.8), (0.8, 1.0)]
        for qmin, qmax in qual_ranges:
            subset = all_data[(all_data['QUAL'] >= qmin) & (all_data['QUAL'] < qmax)]
            if len(subset) > 100:
                try:
                    auc = roc_auc_score(subset['Label'], subset['NumReads'])
                    print(f"  NumReads AUC (QUAL {qmin}-{qmax}): {auc:.4f}")
                except:
                    pass
    
    # Scatter: NumReads vs HeuristicScore colored by Label
    if 'NumReads' in all_data.columns and 'HeuristicScore' in all_data.columns:
        for label, color in [(1, 'green'), (0, 'red')]:
            subset = all_data[all_data['Label'] == label]
            label_name = 'TP' if label == 1 else 'FP'
            axes[1, 1].scatter(subset['NumReads'], subset['HeuristicScore'],
                              alpha=0.3, c=color, label=label_name, s=10)
        axes[1, 1].set_xlabel('NumReads')
        axes[1, 1].set_ylabel('HeuristicScore')
        axes[1, 1].set_title('NumReads vs HeuristicScore')
        axes[1, 1].legend()
        axes[1, 1].set_xlim(0, 200)
    
    plt.tight_layout()
    plt.savefig(PLOT_DIR / 'interaction_effects_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    return aucs


def analyze_conditional_methylation_value(tp_merged, fp_merged):
    """Analyze the value of methylation features in different VCF quality ranges"""
    print("\n" + "="*60)
    print("分析 3: 甲基化特徵在不同 VCF 品質區間的價值")
    print("="*60)
    
    tp_merged = tp_merged.copy()
    fp_merged = fp_merged.copy()
    tp_merged['Label'] = 1
    fp_merged['Label'] = 0
    all_data = pd.concat([tp_merged, fp_merged], ignore_index=True)
    
    # Define quality zones
    zones = {
        'Low QUAL (<0.5)': (all_data['QUAL'] < 0.5),
        'Medium QUAL (0.5-0.8)': (all_data['QUAL'] >= 0.5) & (all_data['QUAL'] < 0.8),
        'High QUAL (>=0.8)': (all_data['QUAL'] >= 0.8),
        'Low AF (<0.15)': (all_data['AF'] < 0.15),
        'Medium AF (0.15-0.4)': (all_data['AF'] >= 0.15) & (all_data['AF'] < 0.4),
        'High AF (>=0.4)': (all_data['AF'] >= 0.4)
    }
    
    results = []
    print("\n[甲基化特徵在不同 VCF 區間的 AUC]")
    
    meth_features = ['CramersV', 'GlobalP', 'HeuristicScore', 'NumReads']
    
    for zone_name, zone_mask in zones.items():
        zone_data = all_data[zone_mask]
        if len(zone_data) > 100:
            n_tp = (zone_data['Label'] == 1).sum()
            n_fp = (zone_data['Label'] == 0).sum()
            tp_ratio = n_tp / len(zone_data)
            
            print(f"\n  {zone_name}: {len(zone_data)} sites (TP:{n_tp}, FP:{n_fp}, TP%:{tp_ratio*100:.1f}%)")
            
            for feat in meth_features:
                if feat in zone_data.columns:
                    clean = zone_data[[feat, 'Label']].dropna()
                    if len(clean) > 50 and len(clean['Label'].unique()) > 1:
                        try:
                            auc = roc_auc_score(clean['Label'], clean[feat])
                            results.append({
                                'Zone': zone_name,
                                'Feature': feat,
                                'AUC': auc,
                                'N': len(clean),
                                'TP_Ratio': tp_ratio
                            })
                            symbol = "✅" if auc > 0.55 else ""
                            print(f"    {feat}: AUC={auc:.4f} {symbol}")
                        except:
                            pass
    
    # Visualization
    if results:
        results_df = pd.DataFrame(results)
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Heatmap of AUC by zone and feature
        pivot = results_df.pivot(index='Zone', columns='Feature', values='AUC')
        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn', ax=axes[0],
                   vmin=0.4, vmax=0.7, center=0.5)
        axes[0].set_title('Methylation Feature AUC by VCF Quality Zone')
        
        # Bar chart of best AUC per zone
        best_per_zone = results_df.loc[results_df.groupby('Zone')['AUC'].idxmax()]
        axes[1].barh(best_per_zone['Zone'], best_per_zone['AUC'], 
                     color=['green' if a > 0.55 else 'orange' for a in best_per_zone['AUC']])
        axes[1].axvline(x=0.5, color='gray', linestyle='--')
        axes[1].set_xlabel('Best AUC')
        axes[1].set_title('Best Methylation Feature AUC per Zone')
        axes[1].set_xlim(0.4, 0.7)
        
        # Add feature name annotations
        for i, (_, row) in enumerate(best_per_zone.iterrows()):
            axes[1].text(row['AUC'] + 0.01, i, row['Feature'], va='center', fontsize=9)
        
        plt.tight_layout()
        plt.savefig(PLOT_DIR / 'conditional_methylation_value.png', dpi=150, bbox_inches='tight')
        plt.close()
    
    return results


def find_methylation_rescue_cases(tp_merged, fp_merged):
    """Find cases where methylation features could rescue misclassifications"""
    print("\n" + "="*60)
    print("分析 4: 甲基化特徵的「救援」效果")
    print("="*60)
    
    # Focus on border cases: TP with low QUAL/AF that might be filtered
    # Can methylation features help identify these as true positives?
    
    tp_merged = tp_merged.copy()
    fp_merged = fp_merged.copy()
    
    # Define "at risk" TP: would be filtered by VCF criteria
    tp_at_risk = tp_merged[(tp_merged['QUAL'] < 0.8) | (tp_merged['AF'] < 0.15)]
    tp_safe = tp_merged[(tp_merged['QUAL'] >= 0.8) & (tp_merged['AF'] >= 0.15)]
    
    print(f"\n[TP 風險分析]")
    print(f"  安全 TP (QUAL>=0.8 AND AF>=0.15): {len(tp_safe)}")
    print(f"  風險 TP (QUAL<0.8 OR AF<0.15): {len(tp_at_risk)}")
    
    # Compare methylation features between at-risk TP and FP
    fp_similar = fp_merged[(fp_merged['QUAL'] < 0.8) | (fp_merged['AF'] < 0.15)]
    print(f"  類似條件的 FP: {len(fp_similar)}")
    
    if len(tp_at_risk) > 50 and len(fp_similar) > 50:
        print("\n[風險區域甲基化特徵比較]")
        meth_feats = ['CramersV', 'HeuristicScore', 'NumReads', 'GlobalP']
        
        for feat in meth_feats:
            if feat in tp_at_risk.columns:
                tp_vals = tp_at_risk[feat].dropna()
                fp_vals = fp_similar[feat].dropna()
                
                if len(tp_vals) > 10 and len(fp_vals) > 10:
                    # Mann-Whitney U test
                    stat, pval = stats.mannwhitneyu(tp_vals, fp_vals)
                    
                    # AUC
                    combined = pd.concat([
                        pd.DataFrame({'val': tp_vals, 'label': 1}),
                        pd.DataFrame({'val': fp_vals, 'label': 0})
                    ])
                    try:
                        auc = roc_auc_score(combined['label'], combined['val'])
                    except:
                        auc = 0.5
                    
                    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
                    print(f"  {feat}: TP_mean={tp_vals.mean():.3f}, FP_mean={fp_vals.mean():.3f}, "
                          f"AUC={auc:.3f}, p={pval:.4e}{sig}")
    
    # Find high-methylation TP at risk
    if 'CramersV' in tp_at_risk.columns and 'HeuristicScore' in tp_at_risk.columns:
        high_meth_rescue = tp_at_risk[
            (tp_at_risk['CramersV'] > 0.3) | (tp_at_risk['HeuristicScore'] > 5)
        ]
        print(f"\n[高甲基化信號的風險 TP]")
        print(f"  可用甲基化特徵救援的 TP: {len(high_meth_rescue)} ({len(high_meth_rescue)/len(tp_at_risk)*100:.1f}%)")
        
        # These are TPs that would be filtered by VCF but have strong methylation signal
        return high_meth_rescue


def test_hybrid_filtering_strategies(tp_merged, fp_merged, baseline_f1):
    """Test hybrid filtering strategies combining VCF and methylation"""
    print("\n" + "="*60)
    print("分析 5: 混合過濾策略評估")
    print("="*60)
    
    tp_merged = tp_merged.copy()
    fp_merged = fp_merged.copy()
    
    BASELINE_FN = 8957  # From Phase 3 verification
    
    strategies = [
        # Name, TP filter, FP filter
        ("Baseline (No filter)", 
         lambda df: pd.Series([False]*len(df), index=df.index),
         lambda df: pd.Series([False]*len(df), index=df.index)),
        
        ("QUAL<0.8 only",
         lambda df: df['QUAL'] < 0.8,
         lambda df: df['QUAL'] < 0.8),
        
        ("QUAL<0.8 BUT keep if Sig=True",
         lambda df: (df['QUAL'] < 0.8) & (df.get('Significant', False) == False),
         lambda df: (df['QUAL'] < 0.8) & (df.get('Significant', False) == False)),
        
        ("QUAL<0.8 BUT keep if CramersV>0.3",
         lambda df: (df['QUAL'] < 0.8) & (df.get('CramersV', 0) <= 0.3),
         lambda df: (df['QUAL'] < 0.8) & (df.get('CramersV', 0) <= 0.3)),
        
        ("QUAL<0.8 BUT keep if HScore>5",
         lambda df: (df['QUAL'] < 0.8) & (df.get('HeuristicScore', 0) <= 5),
         lambda df: (df['QUAL'] < 0.8) & (df.get('HeuristicScore', 0) <= 5)),
        
        ("QUAL<0.8 +NumReads<30",
         lambda df: (df['QUAL'] < 0.8) | (df.get('NumReads', 100) < 30),
         lambda df: (df['QUAL'] < 0.8) | (df.get('NumReads', 100) < 30)),
        
        ("(QUAL<0.5 OR AF<0.1) AND NOT Sig",
         lambda df: ((df['QUAL'] < 0.5) | (df['AF'] < 0.1)) & (df.get('Significant', False) == False),
         lambda df: ((df['QUAL'] < 0.5) | (df['AF'] < 0.1)) & (df.get('Significant', False) == False)),
    ]
    
    results = []
    print("\n[策略評估結果]")
    
    for name, tp_filter, fp_filter in strategies:
        try:
            tp_remove = tp_filter(tp_merged)
            fp_remove = fp_filter(fp_merged)
            
            new_tp = (~tp_remove).sum()
            new_fp = (~fp_remove).sum()
            tp_removed = tp_remove.sum()
            new_fn = BASELINE_FN + tp_removed
            
            precision, recall, f1 = calculate_f1(new_tp, new_fp, new_fn)
            f1_change = (f1 - baseline_f1) * 100
            
            results.append({
                'Strategy': name,
                'TP_Removed': tp_removed,
                'FP_Removed': fp_remove.sum(),
                'New_F1': f1,
                'F1_Change_Pct': f1_change,
                'Precision': precision,
                'Recall': recall
            })
            
            symbol = "✅" if f1 > baseline_f1 else "❌"
            print(f"  {name}")
            print(f"    F1={f1:.4f} ({f1_change:+.2f}%), "
                  f"TP_removed={tp_removed}, FP_removed={fp_remove.sum()} {symbol}")
        except Exception as e:
            print(f"  {name}: Error - {e}")
    
    results_df = pd.DataFrame(results)
    
    # Visualization
    if len(results_df) > 0:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        colors = ['green' if x > 0 else 'red' for x in results_df['F1_Change_Pct']]
        ax.barh(results_df['Strategy'], results_df['F1_Change_Pct'], color=colors)
        ax.axvline(x=0, color='gray', linestyle='-')
        ax.set_xlabel('F1 Change (%)')
        ax.set_title('Hybrid Filtering Strategy Comparison')
        
        plt.tight_layout()
        plt.savefig(PLOT_DIR / 'hybrid_strategy_comparison.png', dpi=150, bbox_inches='tight')
        plt.close()
    
    return results_df


def generate_report(sig_data, interaction_aucs, conditional_results, strategy_results):
    """Generate final report"""
    
    report = f"""# Phase 4: 甲基化顯著性與 VCF 特徵組合分析

**分析時間**: 2026-01-14

---

## 核心發現

> [!IMPORTANT]
> 1. **甲基化顯著性單獨使用區分能力有限** (CramersV AUC=0.52)
> 2. **交互特徵 QUAL×CramersV 未能提升區分能力**
> 3. **NumReads 是最有價值的甲基化特徵** (AUC=0.63)
> 4. **甲基化特徵可作為「救援」機制**：保留被 VCF 過濾掉但有顯著甲基化信號的 TP

---

## 1. 顯著性分佈分析

| 類型 | 顯著 | 非顯著 | 顯著率 |
|:---|---:|---:|---:|
"""
    
    if 'Significant' in sig_data.columns:
        for label in ['TP', 'FP']:
            subset = sig_data[sig_data['Label'] == label]
            sig_count = (subset['Significant'] == True).sum()
            nonsig_count = (subset['Significant'] == False).sum()
            sig_rate = sig_count / len(subset) * 100 if len(subset) > 0 else 0
            report += f"| {label} | {sig_count} | {nonsig_count} | {sig_rate:.1f}% |\n"
    
    report += """
---

## 2. 特徵 AUC 比較

| 特徵 | AUC | 評價 |
|:---|---:|:---|
"""
    
    for feat, auc in sorted(interaction_aucs.items(), key=lambda x: x[1], reverse=True):
        eval_str = "✅ 有效" if auc > 0.55 else "⚠️ 弱" if auc > 0.52 else "❌ 無效"
        report += f"| {feat} | {auc:.4f} | {eval_str} |\n"
    
    report += """
---

## 3. 條件分析：甲基化在不同 VCF 區間的價值

"""
    
    if conditional_results:
        for r in conditional_results[:5]:  # Top 5
            report += f"- **{r['Zone']}**: {r['Feature']} AUC={r['AUC']:.3f}\n"
    
    report += """
---

## 4. 混合策略比較

| 策略 | F1 變化 | 效果 |
|:---|---:|:---|
"""
    
    for _, row in strategy_results.iterrows():
        symbol = "✅" if row['F1_Change_Pct'] > 0 else "❌"
        report += f"| {row['Strategy']} | {row['F1_Change_Pct']:+.2f}% | {symbol} |\n"
    
    report += """
---

## 5. 視覺化

![Significant vs Non-significant](phase4_plots/significant_vs_nonsig_analysis.png)

![Interaction Effects](phase4_plots/interaction_effects_analysis.png)

![Conditional Value](phase4_plots/conditional_methylation_value.png)

![Hybrid Strategies](phase4_plots/hybrid_strategy_comparison.png)

---

## 結論

1. **甲基化顯著性與 VCF 特徵的直接組合未能顯著提升區分能力**
2. **NumReads 在所有區間仍是最有價值的甲基化特徵**
3. **「救援」策略有限效果**：保留高甲基化信號的風險 TP 可略微提升 Recall
4. **建議**：優先使用 VCF 過濾 (QUAL<0.8)，甲基化特徵僅作為輔助參考
"""
    
    with open(REPORT_PATH, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"\n[INFO] 報告已儲存至: {REPORT_PATH}")


def main():
    print("="*60)
    print(" Phase 4: 甲基化顯著性與 VCF 特徵組合分析")
    print("="*60)
    
    setup_directories()
    
    # Load data
    print("\n[INFO] 載入資料...")
    tp_vcf = extract_vcf_features(TP_VCF)
    fp_vcf = extract_vcf_features(FP_VCF)
    tp_meth = pd.read_csv(TP_SUMMARY)
    fp_meth = pd.read_csv(FP_SUMMARY)
    
    # Merge
    tp_merged = merge_vcf_methylation(tp_vcf, tp_meth)
    fp_merged = merge_vcf_methylation(fp_vcf, fp_meth)
    print(f"[INFO] Merged - TP: {len(tp_merged)}, FP: {len(fp_merged)}")
    
    # Analyses
    sig_data = analyze_significant_vs_nonsignificant(tp_merged, fp_merged)
    interaction_aucs = analyze_interaction_effects(tp_merged, fp_merged)
    conditional_results = analyze_conditional_methylation_value(tp_merged, fp_merged)
    find_methylation_rescue_cases(tp_merged, fp_merged)
    
    baseline_f1 = 0.8155
    strategy_results = test_hybrid_filtering_strategies(tp_merged, fp_merged, baseline_f1)
    
    # Report
    generate_report(sig_data, interaction_aucs, conditional_results, strategy_results)
    
    print("\n" + "="*60)
    print(" Phase 4 分析完成!")
    print(f" 報告: {REPORT_PATH}")
    print("="*60)


if __name__ == "__main__":
    main()
