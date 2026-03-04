#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Depth Effect Analysis Tool for InterSubMod
"""

import pandas as pd
import numpy as np
import argparse
import os
import sys

# Add tools directory to path for font_config import
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Configure CJK font support (must be before importing pyplot)
try:
    from font_config import configure_matplotlib_fonts
    configure_matplotlib_fonts()
except ImportError:
    pass  # Fallback: font_config not available

import matplotlib.pyplot as plt
import seaborn as sns

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze relationship between Depth, Cramers V and TP/FP ratio.")
    parser.add_argument("--tp-csv", required=True, help="Path to TP significance_summary.csv")
    parser.add_argument("--fp-csv", required=True, help="Path to FP significance_summary.csv")
    parser.add_argument("--output-dir", required=True, help="Directory to save plots and reports")
    return parser.parse_args()

def analyze_depth_effect(tp_csv, fp_csv, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Loading TP: {tp_csv}")
    tp_df = pd.read_csv(tp_csv)
    tp_df['Label'] = 'TP'
    
    print(f"Loading FP: {fp_csv}")
    fp_df = pd.read_csv(fp_csv)
    fp_df['Label'] = 'FP'
    
    # Combine
    df = pd.concat([tp_df, fp_df], ignore_index=True)
    
    # Filter for reads >= 3 as requested
    df = df[df['NumReads'] >= 3]
    
    print("-" * 60)
    print(f"Total sites with Reads >= 3: {len(df)}")
    print(f"TP: {len(df[df['Label']=='TP'])}")
    print(f"FP: {len(df[df['Label']=='FP'])}")
    print("-" * 60)

    # ==========================================
    # Granular Analysis (Per Depth Unit)
    # ==========================================
    
    # Group by NumReads
    depth_stats = df.groupby('NumReads')['Label'].value_counts().unstack(fill_value=0)
    
    # Ensure both TP and FP columns exist
    if 'TP' not in depth_stats.columns: depth_stats['TP'] = 0
    if 'FP' not in depth_stats.columns: depth_stats['FP'] = 0
    
    depth_stats['Total'] = depth_stats['TP'] + depth_stats['FP']
    depth_stats['Ratio_TP_FP'] = depth_stats['TP'] / (depth_stats['FP'].replace(0, 0.001))
    depth_stats['Precision'] = depth_stats['TP'] / depth_stats['Total']
    
    # Reset index to make NumReads a column
    depth_stats = depth_stats.reset_index()
    
    # Separate low/mid depth for detailed view, and maybe bin high depths if needed
    # But user asked for "per unit", so let's output everything to CSV, 
    # and maybe truncate plot x-axis or use log scale if it goes too high.
    
    # Save Detailed Table
    table_path = os.path.join(output_dir, "depth_performance_stats.csv")
    depth_stats.to_csv(table_path, index=False)
    print(f"Saved detailed stats to {table_path}")
    
    # ==========================================
    # Visualization
    # ==========================================
    
    # 1. TP/FP Ratio vs Depth (Line Plot)
    plt.figure(figsize=(15, 6))
    
    # Limit plot to reasonable range to see details (e.g., up to 100 or 200)
    # or just plot all if not too sparse.
    plot_df = depth_stats[depth_stats['NumReads'] <= 200]
    
    sns.lineplot(data=plot_df, x='NumReads', y='Ratio_TP_FP', marker='o', markersize=4, label='TP/FP Ratio')
    plt.axhline(1.0, color='red', linestyle='--', label='Ratio = 1.0 (Neutral)')
    plt.title("TP/FP Ratio vs Read Depth (up to 200x)")
    plt.xlabel("Read Depth")
    plt.ylabel("TP/FP Ratio")
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "plot_ratio_vs_depth.png"), dpi=150)
    plt.close()
    
    # 2. Counts Distribution (Stacked Bar or Area)
    plt.figure(figsize=(15, 6))
    plt.bar(plot_df['NumReads'], plot_df['TP'], label='TP', alpha=0.7, color='green')
    plt.bar(plot_df['NumReads'], plot_df['FP'], bottom=plot_df['TP'], label='FP', alpha=0.7, color='red')
    plt.title("TP and FP Counts vs Read Depth (up to 200x)")
    plt.xlabel("Read Depth")
    plt.ylabel("Count")
    plt.legend()
    plt.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "plot_counts_vs_depth.png"), dpi=150)
    plt.close()
    
    # 3. Precision vs Depth
    plt.figure(figsize=(15, 6))
    sns.lineplot(data=plot_df, x='NumReads', y='Precision', marker='o', markersize=4, color='purple')
    plt.axhline(0.5, color='gray', linestyle='--', label='50% Precision')
    plt.title("Precision (TP / Total) vs Read Depth")
    plt.xlabel("Read Depth")
    plt.ylabel("Precision")
    plt.ylim(0, 1.05)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "plot_precision_vs_depth.png"), dpi=150)
    plt.close()

    # ==========================================
    # Legacy Heatmap (V vs Depth) - Kept but simplified or adjusted?
    # User asked for "analyze depth effect... output table and images"
    # The V-value interaction is still useful, but maybe we just bin V and keep Depth granular?
    # Or keep the old heatmap as a supplementary view. Let's keep a simplified version.
    # ==========================================
    
    # We will still generate the V-Heatmap but maybe keep the bins for V, 
    # but for Depth logic, the user specifically wanted granular depth stats.
    # The heatmap requires bins to be readable. Let's make a "Detailed Heatmap" if possible,
    # or just stick to the line plots which are better for "per unit".
    # I'll add a section in the report for the "Sweet Spots" based on the granular data.

    # ==========================================
    # Report Generation
    # ==========================================
    
    report_path = os.path.join(output_dir, "depth_detailed_report.md")
    with open(report_path, 'w') as f:
        f.write("# 深度效應詳細分析報告 (Granular Depth Analysis)\n\n")
        f.write(f"**來源**: TP={tp_csv}, FP={fp_csv}\n")
        f.write(f"**總位點數**: {len(df)}\n\n")
        
        f.write("## 1. 關鍵發現 (基於前 50 個深度單位)\n\n")
        f.write("下表列出深度 3 到 50 的詳細表現：\n\n")
        
        # Markdown Table for low depths
        f.write("| Depth | TP | FP | Total | Ratio (TP/FP) | Precision (TP/Tot) |\n")
        f.write("|---|---|---|---|---|---|\n")
        
        for _, row in depth_stats[depth_stats['NumReads'] <= 50].iterrows():
            d = int(row['NumReads'])
            tp = int(row['TP'])
            fp = int(row['FP'])
            tot = int(row['Total'])
            ratio = row['Ratio_TP_FP']
            prec = row['Precision']
            
            # Highlight good rows
            formatting = ""
            if ratio > 2.0: formatting = "**"
            
            f.write(f"| {d} | {tp} | {fp} | {tot} | {formatting}{ratio:.2f}{formatting} | {prec:.2%}{formatting} |\n")
            
        f.write("\n\n")
        
        f.write("## 2. 說明\n")
        f.write("- **Ratio > 1.0**: TP 多於 FP (好的信號)\n")
        f.write("- **Ratio < 1.0**: FP 多於 TP (噪音主導)\n")
        f.write("- 請參考 `depth_performance_stats.csv` 查看完整數據。\n")
        f.write("- 產生的圖表:\n")
        f.write("  - `plot_ratio_vs_depth.png`: TP/FP 比例隨深度變化\n")
        f.write("  - `plot_counts_vs_depth.png`: TP/FP 數量堆疊圖\n")
        f.write("  - `plot_precision_vs_depth.png`: 準確率隨深度變化\n")
        
    print(f"Analysis complete. Report saved to {report_path}")

if __name__ == "__main__":
    args = parse_args()
    analyze_depth_effect(args.tp_csv, args.fp_csv, args.output_dir)
