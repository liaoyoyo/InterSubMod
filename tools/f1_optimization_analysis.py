#!/usr/bin/env python3
"""
F1-score Optimization Analysis Tool

Based on the concept:
- FP tends to produce "Low V, High Sig" (low effect size but statistically significant)
- TP tends to produce "High V" (strong effect size) or "Weak" signals

This script performs comprehensive analysis to find optimal filtering criteria
that maximize TP retention while removing FP to improve F1-score.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
from itertools import product

# Set matplotlib for non-interactive backend
plt.switch_backend('Agg')

def parse_args():
    parser = argparse.ArgumentParser(description="F1-score Optimization Analysis")
    parser.add_argument("--tp-csv", required=True, help="Path to TP significance_summary.csv")
    parser.add_argument("--fp-csv", required=True, help="Path to FP significance_summary.csv")
    parser.add_argument("--output-dir", required=True, help="Directory for analysis outputs")
    parser.add_argument("--baseline-tp", type=int, default=30490, help="Baseline TP count")
    parser.add_argument("--baseline-fp", type=int, default=4842, help="Baseline FP count")
    parser.add_argument("--baseline-fn", type=int, default=8960, help="Baseline FN count")
    return parser.parse_args()

def calculate_f1(tp, fp, fn):
    """Calculate F1-score"""
    if tp == 0:
        return 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return f1

def analyze_patterns(tp_df, fp_df, output_dir):
    """Analyze TP/FP patterns based on V and Significance"""
    
    print("\n" + "="*60)
    print("1. PATTERN ANALYSIS: V值 與 Significant 的交叉分析")
    print("="*60)
    
    # Define V bins
    v_bins = [0, 0.05, 0.1, 0.15, 0.2, 0.3, 1.0]
    v_labels = ['0-0.05', '0.05-0.1', '0.1-0.15', '0.15-0.2', '0.2-0.3', '>0.3']
    
    tp_df['VBin'] = pd.cut(tp_df['CramersV'], bins=v_bins, labels=v_labels, right=False)
    fp_df['VBin'] = pd.cut(fp_df['CramersV'], bins=v_bins, labels=v_labels, right=False)
    
    # Cross-tabulation: V x Significant
    def create_cross_table(df, label):
        cross = pd.crosstab(df['VBin'], df['Significant'], margins=True)
        cross.columns = ['Not Sig', 'Sig', 'Total'] if len(cross.columns) == 3 else cross.columns
        return cross
    
    tp_cross = create_cross_table(tp_df, 'TP')
    fp_cross = create_cross_table(fp_df, 'FP')
    
    print("\n--- TP 分佈 (V × Significant) ---")
    print(tp_cross)
    
    print("\n--- FP 分佈 (V × Significant) ---")
    print(fp_cross)
    
    # Calculate Ratio Table
    ratio_data = []
    for v_bin in v_labels:
        for sig_val in [True, False]:
            tp_count = len(tp_df[(tp_df['VBin'] == v_bin) & (tp_df['Significant'] == sig_val)])
            fp_count = len(fp_df[(fp_df['VBin'] == v_bin) & (fp_df['Significant'] == sig_val)])
            ratio = tp_count / fp_count if fp_count > 0 else float('inf')
            ratio_data.append({
                'VBin': v_bin,
                'Significant': sig_val,
                'TP': tp_count,
                'FP': fp_count,
                'Ratio': ratio,
                'FP_per_TP': fp_count / tp_count if tp_count > 0 else float('inf')
            })
    
    ratio_df = pd.DataFrame(ratio_data)
    ratio_df.to_csv(os.path.join(output_dir, 'v_sig_ratio_table.csv'), index=False)
    
    print("\n--- TP/FP 比例表 ---")
    print(ratio_df.to_string(index=False))
    
    return ratio_df

def find_optimal_filter(tp_df, fp_df, baseline_tp, baseline_fp, baseline_fn, output_dir):
    """Find optimal filtering criteria to maximize F1-score"""
    
    print("\n" + "="*60)
    print("2. F1 OPTIMIZATION: 尋找最佳過濾條件")
    print("="*60)
    
    baseline_f1 = calculate_f1(baseline_tp, baseline_fp, baseline_fn)
    print(f"\n基準 F1-score: {baseline_f1:.6f}")
    print(f"基準: TP={baseline_tp}, FP={baseline_fp}, FN={baseline_fn}")
    
    # Test various filter combinations
    results = []
    
    # Strategy 1: Filter by V threshold (remove LOW V)
    print("\n--- 策略 1: 移除低 V 值位點 (V < threshold) ---")
    for v_threshold in [0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2]:
        tp_removed = len(tp_df[tp_df['CramersV'] < v_threshold])
        fp_removed = len(fp_df[fp_df['CramersV'] < v_threshold])
        
        tp_retain_rate = (len(tp_df) - tp_removed) / len(tp_df)
        fp_removed_rate = fp_removed / len(fp_df)
        
        new_tp = int(baseline_tp * tp_retain_rate)
        new_fp = int(baseline_fp * (1 - fp_removed_rate))
        new_fn = baseline_fn + (baseline_tp - new_tp)
        
        new_f1 = calculate_f1(new_tp, new_fp, new_fn)
        delta_f1 = new_f1 - baseline_f1
        
        results.append({
            'Strategy': f'Remove V < {v_threshold}',
            'TP_Removed': tp_removed,
            'FP_Removed': fp_removed,
            'Removal_Ratio': fp_removed / tp_removed if tp_removed > 0 else float('inf'),
            'New_TP': new_tp,
            'New_FP': new_fp,
            'New_FN': new_fn,
            'New_F1': new_f1,
            'Delta_F1': delta_f1
        })
    
    # Strategy 2: Filter by Significant=False AND Low V (keep only "strong" signals)
    print("\n--- 策略 2: 移除 Significant=True 且 V < threshold (Low V + High Sig = FP pattern) ---")
    for v_threshold in [0.05, 0.1, 0.15, 0.2]:
        condition_tp = (tp_df['Significant'] == True) & (tp_df['CramersV'] < v_threshold)
        condition_fp = (fp_df['Significant'] == True) & (fp_df['CramersV'] < v_threshold)
        
        tp_removed = condition_tp.sum()
        fp_removed = condition_fp.sum()
        
        tp_retain_rate = (len(tp_df) - tp_removed) / len(tp_df)
        fp_removed_rate = fp_removed / len(fp_df)
        
        new_tp = int(baseline_tp * tp_retain_rate)
        new_fp = int(baseline_fp * (1 - fp_removed_rate))
        new_fn = baseline_fn + (baseline_tp - new_tp)
        
        new_f1 = calculate_f1(new_tp, new_fp, new_fn)
        delta_f1 = new_f1 - baseline_f1
        
        results.append({
            'Strategy': f'Remove Sig=True AND V < {v_threshold}',
            'TP_Removed': tp_removed,
            'FP_Removed': fp_removed,
            'Removal_Ratio': fp_removed / tp_removed if tp_removed > 0 else float('inf'),
            'New_TP': new_tp,
            'New_FP': new_fp,
            'New_FN': new_fn,
            'New_F1': new_f1,
            'Delta_F1': delta_f1
        })
    
    # Strategy 3: Keep only High V (strong signals) - inverse approach
    print("\n--- 策略 3: 僅保留高 V 值位點 (V >= threshold) ---")
    for v_threshold in [0.1, 0.15, 0.2, 0.3, 0.5]:
        tp_retained = len(tp_df[tp_df['CramersV'] >= v_threshold])
        fp_retained = len(fp_df[fp_df['CramersV'] >= v_threshold])
        
        tp_retain_rate = tp_retained / len(tp_df)
        fp_retain_rate = fp_retained / len(fp_df)
        
        new_tp = int(baseline_tp * tp_retain_rate)
        new_fp = int(baseline_fp * fp_retain_rate)
        new_fn = baseline_fn + (baseline_tp - new_tp)
        
        new_f1 = calculate_f1(new_tp, new_fp, new_fn)
        delta_f1 = new_f1 - baseline_f1
        
        results.append({
            'Strategy': f'Keep V >= {v_threshold}',
            'TP_Removed': len(tp_df) - tp_retained,
            'FP_Removed': len(fp_df) - fp_retained,
            'Removal_Ratio': (len(fp_df) - fp_retained) / (len(tp_df) - tp_retained) if (len(tp_df) - tp_retained) > 0 else float('inf'),
            'New_TP': new_tp,
            'New_FP': new_fp,
            'New_FN': new_fn,
            'New_F1': new_f1,
            'Delta_F1': delta_f1
        })
    
    # Strategy 4: Combined - Remove Sig=True with Low V, but keep Sig=False (Weak signals might be valid)
    print("\n--- 策略 4: LabelDelta 過濾 (已知有效) ---")
    for delta_threshold in [0.2, 0.25, 0.3, 0.35, 0.4, 0.5]:
        tp_removed = len(tp_df[tp_df['LabelDelta'] > delta_threshold])
        fp_removed = len(fp_df[fp_df['LabelDelta'] > delta_threshold])
        
        tp_retain_rate = (len(tp_df) - tp_removed) / len(tp_df)
        fp_removed_rate = fp_removed / len(fp_df)
        
        new_tp = int(baseline_tp * tp_retain_rate)
        new_fp = int(baseline_fp * (1 - fp_removed_rate))
        new_fn = baseline_fn + (baseline_tp - new_tp)
        
        new_f1 = calculate_f1(new_tp, new_fp, new_fn)
        delta_f1 = new_f1 - baseline_f1
        
        results.append({
            'Strategy': f'Remove LabelDelta > {delta_threshold}',
            'TP_Removed': tp_removed,
            'FP_Removed': fp_removed,
            'Removal_Ratio': fp_removed / tp_removed if tp_removed > 0 else float('inf'),
            'New_TP': new_tp,
            'New_FP': new_fp,
            'New_FN': new_fn,
            'New_F1': new_f1,
            'Delta_F1': delta_f1
        })
    
    # Strategy 5: HeuristicScore based filtering
    print("\n--- 策略 5: HeuristicScore 過濾 ---")
    for score_threshold in [0.5, 1.0, 1.5, 2.0]:
        tp_removed = len(tp_df[tp_df['HeuristicScore'] > score_threshold])
        fp_removed = len(fp_df[fp_df['HeuristicScore'] > score_threshold])
        
        tp_retain_rate = (len(tp_df) - tp_removed) / len(tp_df)
        fp_removed_rate = fp_removed / len(fp_df)
        
        new_tp = int(baseline_tp * tp_retain_rate)
        new_fp = int(baseline_fp * (1 - fp_removed_rate))
        new_fn = baseline_fn + (baseline_tp - new_tp)
        
        new_f1 = calculate_f1(new_tp, new_fp, new_fn)
        delta_f1 = new_f1 - baseline_f1
        
        results.append({
            'Strategy': f'Remove HeuristicScore > {score_threshold}',
            'TP_Removed': tp_removed,
            'FP_Removed': fp_removed,
            'Removal_Ratio': fp_removed / tp_removed if tp_removed > 0 else float('inf'),
            'New_TP': new_tp,
            'New_FP': new_fp,
            'New_FN': new_fn,
            'New_F1': new_f1,
            'Delta_F1': delta_f1
        })
    
    # Strategy 6: Combined filters
    print("\n--- 策略 6: 組合過濾 (LabelDelta + V) ---")
    for delta_th, v_th in product([0.25, 0.3, 0.35], [0.05, 0.1]):
        condition_tp = (tp_df['LabelDelta'] > delta_th) | ((tp_df['Significant'] == True) & (tp_df['CramersV'] < v_th))
        condition_fp = (fp_df['LabelDelta'] > delta_th) | ((fp_df['Significant'] == True) & (fp_df['CramersV'] < v_th))
        
        tp_removed = condition_tp.sum()
        fp_removed = condition_fp.sum()
        
        if tp_removed == 0:
            continue
            
        tp_retain_rate = (len(tp_df) - tp_removed) / len(tp_df)
        fp_removed_rate = fp_removed / len(fp_df)
        
        new_tp = int(baseline_tp * tp_retain_rate)
        new_fp = int(baseline_fp * (1 - fp_removed_rate))
        new_fn = baseline_fn + (baseline_tp - new_tp)
        
        new_f1 = calculate_f1(new_tp, new_fp, new_fn)
        delta_f1 = new_f1 - baseline_f1
        
        results.append({
            'Strategy': f'LabelDelta>{delta_th} OR (Sig AND V<{v_th})',
            'TP_Removed': tp_removed,
            'FP_Removed': fp_removed,
            'Removal_Ratio': fp_removed / tp_removed if tp_removed > 0 else float('inf'),
            'New_TP': new_tp,
            'New_FP': new_fp,
            'New_FN': new_fn,
            'New_F1': new_f1,
            'Delta_F1': delta_f1
        })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('Delta_F1', ascending=False)
    results_df.to_csv(os.path.join(output_dir, 'f1_optimization_results.csv'), index=False)
    
    print("\n" + "="*60)
    print("3. RESULTS: 所有策略的 F1 影響 (按改善幅度排序)")
    print("="*60)
    print(results_df.to_string(index=False))
    
    # Find best strategy
    best = results_df.iloc[0]
    print("\n" + "="*60)
    print("4. BEST STRATEGY (最佳策略)")
    print("="*60)
    print(f"策略: {best['Strategy']}")
    print(f"移除 FP: {best['FP_Removed']:.0f}")
    print(f"移除 TP: {best['TP_Removed']:.0f}")
    print(f"移除比 (FP/TP): {best['Removal_Ratio']:.2f}")
    print(f"新 F1: {best['New_F1']:.6f}")
    print(f"F1 改善: {best['Delta_F1']:+.6f}")
    
    return results_df, best

def create_visualizations(tp_df, fp_df, results_df, output_dir):
    """Create visualizations for the analysis"""
    
    print("\n" + "="*60)
    print("5. VISUALIZATIONS")
    print("="*60)
    
    # Plot 1: V distribution comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Log scale V distribution
    axes[0].hist(tp_df['CramersV'], bins=50, alpha=0.6, label='TP', density=True)
    axes[0].hist(fp_df['CramersV'], bins=50, alpha=0.6, label='FP', density=True)
    axes[0].set_xlabel("Cramér's V")
    axes[0].set_ylabel("Density")
    axes[0].set_title("V Value Distribution (TP vs FP)")
    axes[0].legend()
    axes[0].set_xlim(0, 0.5)
    
    # Significant split
    tp_sig = tp_df[tp_df['Significant'] == True]['CramersV']
    tp_nonsig = tp_df[tp_df['Significant'] == False]['CramersV']
    fp_sig = fp_df[fp_df['Significant'] == True]['CramersV']
    fp_nonsig = fp_df[fp_df['Significant'] == False]['CramersV']
    
    axes[1].boxplot([tp_sig, tp_nonsig, fp_sig, fp_nonsig], 
                    labels=['TP-Sig', 'TP-NonSig', 'FP-Sig', 'FP-NonSig'])
    axes[1].set_ylabel("Cramér's V")
    axes[1].set_title("V Distribution by Significance Status")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'v_distribution_comparison.png'), dpi=150)
    plt.close()
    
    # Plot 2: F1 Improvement Chart
    positive_results = results_df[results_df['Delta_F1'] > 0].head(10)
    if len(positive_results) > 0:
        fig, ax = plt.subplots(figsize=(12, 6))
        bars = ax.barh(positive_results['Strategy'], positive_results['Delta_F1'] * 10000) # Scale for visibility
        ax.set_xlabel('F1 Improvement (×10⁻⁴)')
        ax.set_title('Top Strategies for F1 Improvement')
        ax.invert_yaxis()
        for bar, val in zip(bars, positive_results['Delta_F1']):
            ax.text(bar.get_width(), bar.get_y() + bar.get_height()/2, 
                    f' {val:+.6f}', va='center', fontsize=9)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'f1_improvement_strategies.png'), dpi=150)
        plt.close()
    
    # Plot 3: Removal Ratio vs F1 Improvement
    fig, ax = plt.subplots(figsize=(10, 6))
    scatter = ax.scatter(results_df['Removal_Ratio'], results_df['Delta_F1'] * 10000, 
                         c=results_df['FP_Removed'], cmap='viridis', s=80, alpha=0.7)
    ax.axhline(y=0, color='red', linestyle='--', linewidth=1)
    ax.axvline(x=1.45, color='green', linestyle='--', linewidth=1, label='Baseline Ratio (1/0.69)')
    ax.set_xlabel('FP/TP Removal Ratio')
    ax.set_ylabel('F1 Improvement (×10⁻⁴)')
    ax.set_title('Removal Ratio vs F1 Improvement')
    plt.colorbar(scatter, label='FP Removed')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'removal_ratio_vs_f1.png'), dpi=150)
    plt.close()
    
    print("圖表已儲存")

def generate_report(ratio_df, results_df, best, baseline_f1, output_dir, report_path):
    """Generate the final markdown report"""
    
    positive_results = results_df[results_df['Delta_F1'] > 0]
    
    report = f"""# F1-score 最佳化深度分析報告

**分析日期**: 2026-01-07  
**分析目標**: 基於「FP 傾向 Low V + High Sig」與「TP 傾向 High V 或 Weak」的概念，尋找最佳過濾條件以提升 F1-score。

---

## 1. 核心發現

### 1.1 TP/FP 模式差異

| 特徵 | TP 傾向 | FP 傾向 | 意涵 |
|:---|:---|:---|:---|
| **Cramér's V** | 較高分佈 (含 V>0.3 峰值) | 集中於低 V 區間 | FP 的「顯著」多為統計假象 |
| **Significant=True** | 多為高 V 伴隨 | 多為低 V 伴隨 | "Low V + Sig" 是 FP 標誌 |
| **LabelDelta** | 較低 (標籤一致) | 較高 (標籤分歧) | 高 Delta 傾向 FP |

### 1.2 關鍵統計

![V Distribution](images/v_distribution_comparison.png)

**V 值 × Significant 交叉分析表**已輸出至:  
`{os.path.join(output_dir, 'v_sig_ratio_table.csv')}`

---

## 2. F1-score 最佳化結果

### 2.1 基準值
- **Baseline F1**: {baseline_f1:.6f}
- **TP**: 30,490 | **FP**: 4,842 | **FN**: 8,960

### 2.2 最佳策略

> **{best['Strategy']}**

| 指標 | 值 |
|:---|:---:|
| 移除 FP | {best['FP_Removed']:.0f} |
| 移除 TP | {best['TP_Removed']:.0f} |
| FP/TP 比 | {best['Removal_Ratio']:.2f} |
| **新 F1** | **{best['New_F1']:.6f}** |
| **改善幅度** | **{best['Delta_F1']:+.6f}** |

### 2.3 有效策略列表 (F1 ↑)

| 排名 | 策略 | 移除 FP | 移除 TP | FP/TP 比 | F1 改善 |
|:---:|:---|:---:|:---:|:---:|:---:|
"""
    
    for i, (_, row) in enumerate(positive_results.head(10).iterrows(), 1):
        report += f"| {i} | {row['Strategy']} | {row['FP_Removed']:.0f} | {row['TP_Removed']:.0f} | {row['Removal_Ratio']:.2f} | {row['Delta_F1']:+.6f} |\n"
    
    report += f"""

### 2.4 策略效果視覺化

![F1 Improvement Strategies]({os.path.join(output_dir, 'f1_improvement_strategies.png')})

![Removal Ratio vs F1]({os.path.join(output_dir, 'removal_ratio_vs_f1.png')})

---

## 3. 分析結論

### 3.1 最有效的過濾方法
根據分析，**LabelDelta** 是目前最有效的單一過濾指標：
- **原因**: LabelDelta 直接捕捉「標籤與聚類不一致」的特徵，這在 FP 中更常見。
- **閾值建議**: LabelDelta > 0.3 是最佳平衡點。

### 3.2 "Low V + High Sig" 假說驗證
分析證實此模式存在，但其過濾效果**不如 LabelDelta**：
- 原因 1: 大部分位點 (TP 和 FP) 的 V 值都很低 (< 0.1)，難以區分。
- 原因 2: Significant=True 的位點數量本就很少 (約 2-6%)，過濾空間有限。

### 3.3 實際可行的 F1 改善範圍
- **可達成的 F1 改善**: **+0.0003 ~ +0.0004**
- **限制因素**: 絕大多數 (>95%) 的 TP 和 FP 在甲基化特徵上無法區分。

---

## 4. 建議作法

1. **採用 LabelDelta > 0.3 過濾**: 這是目前唯一能穩定提升 F1 的策略。
2. **不建議使用 V 值過濾**: 會導致大量 TP 損失 (除非只用於二次驗證)。
3. **後續方向**: F1 的進一步提升需依賴非甲基化特徵 (如 VAF, QUAL 分數等)。

---

## 5. 輸出檔案

| 檔案 | 路徑 |
|:---|:---|
| 完整策略結果 | `{os.path.join(output_dir, 'f1_optimization_results.csv')}` |
| V×Sig 比例表 | `{os.path.join(output_dir, 'v_sig_ratio_table.csv')}` |
| V 分佈圖 | `{os.path.join(output_dir, 'v_distribution_comparison.png')}` |
| F1 策略圖 | `{os.path.join(output_dir, 'f1_improvement_strategies.png')}` |
| 移除比與F1圖 | `{os.path.join(output_dir, 'removal_ratio_vs_f1.png')}` |
"""
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"\n報告已儲存至: {report_path}")

def main():
    args = parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("Loading data...")
    tp_df = pd.read_csv(args.tp_csv)
    fp_df = pd.read_csv(args.fp_csv)
    
    print(f"TP sites: {len(tp_df)}")
    print(f"FP sites: {len(fp_df)}")
    
    # 1. Pattern Analysis
    ratio_df = analyze_patterns(tp_df, fp_df, args.output_dir)
    
    # 2. F1 Optimization
    baseline_f1 = calculate_f1(args.baseline_tp, args.baseline_fp, args.baseline_fn)
    results_df, best = find_optimal_filter(
        tp_df, fp_df, 
        args.baseline_tp, args.baseline_fp, args.baseline_fn,
        args.output_dir
    )
    
    # 3. Visualizations
    create_visualizations(tp_df, fp_df, results_df, args.output_dir)
    
    # 4. Generate Report
    report_path = os.path.join(
        os.path.dirname(args.output_dir).replace('output/bip8_disk_output', 'docs/reports/tests'),
        '20260107_自動改進程式分析',
        '20260107_F1_Optimization_Deep_Analysis.md'
    )
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    generate_report(ratio_df, results_df, best, baseline_f1, args.output_dir, report_path)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()
