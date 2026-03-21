#!/usr/bin/env python3
"""
analyze_feature_distribution.py - Analyze AD/CV distribution for TP vs FP.

Target Samples: HCC1395 (Reference), COLO829 (High FP)
"""

import os
import csv
import gzip
import numpy as np

SAMPLES = {
    "HCC1395": "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395/20260212_2",
    "COLO829": "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/COLO829/20260212_1",
}

def read_values(csv_path, col_name):
    vals = []
    if not os.path.exists(csv_path):
        return vals
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                val = float(row.get(col_name, 0))
                if col_name == "AlleleDelta":
                    val = abs(val)
                vals.append(val)
            except:
                pass
    return vals

def print_histo(name, vals, bins):
    hist, bin_edges = np.histogram(vals, bins=bins)
    total = len(vals)
    print(f"--- {name} (N={total}) ---")
    for i in range(len(hist)):
        pct = 100 * hist[i] / total if total > 0 else 0
        bar = "#" * int(pct / 2)
        range_str = f"{bin_edges[i]:.2f}-{bin_edges[i+1]:.2f}"
        print(f"{range_str:<12} | {hist[i]:<6} ({pct:5.1f}%) | {bar}")

def analyze(sample, base_dir):
    print(f"\n================ {sample} ================")
    tp_csv = os.path.join(base_dir, "intersubmod_tp", "significance_summary.csv")
    fp_csv = os.path.join(base_dir, "intersubmod_fp", "significance_summary.csv")
    
    # 1. AlleleDelta
    print("\n[AlleleDelta Distribution]")
    bins_ad = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0]
    tp_ad = read_values(tp_csv, "AlleleDelta")
    fp_ad = read_values(fp_csv, "AlleleDelta")
    print_histo("TP", tp_ad, bins_ad)
    print_histo("FP", fp_ad, bins_ad)

    # 2. CramersV
    print("\n[CramersV Distribution]")
    bins_cv = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
    tp_cv = read_values(tp_csv, "CramersV")
    fp_cv = read_values(fp_csv, "CramersV")
    print_histo("TP", tp_cv, bins_cv)
    print_histo("FP", fp_cv, bins_cv)

if __name__ == "__main__":
    for s, p in SAMPLES.items():
        analyze(s, p)
