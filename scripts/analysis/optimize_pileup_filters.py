#!/usr/bin/env python3
"""
optimize_pileup_filters.py - Find optimal filter parameters for Pileup mode samples.

Goal: Maximize F1-score by adjusting thresholds for:
- AlleleDelta (AD)
- CramersV (CV)
- VAF

Features:
- Loads data from 7 samples (TP/FP significance summaries + VCF features)
- Performs grid search on filter parameters
- Outputs optimal parameters per sample and global best
"""

import os
import json
import csv
import gzip
import glob
import numpy as np
from collections import defaultdict
import multiprocessing

# Configuration
SAMPLES = {
    "HCC1395": "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395/20260212_2",
    "HCC1395_DORADO": "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1395_DORADO/20260212",
    "COLO829": "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/COLO829/20260212_1",
    "H1437": "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/H1437/20260212",
    "H2009": "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/H2009/20260212",
    "HCC1937": "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1937/20260212",
    "HCC1954": "/big8_disk/liaoyoyo2001/InterSubMod/output/bip8_disk_output/s-pure-pileup/HCC1954/20260212",
}

# Grid Search Ranges
# Current: AD > 0.25, CV < 0.05, VAF < 0.24
AD_RANGE = [0.10, 0.15, 0.20, 0.25, 0.30]  # Higher = More strict (less FP removed, more TP kept)
CV_RANGE = [0.03, 0.05, 0.07, 0.10]        # Lower = More strict (less FP removed) -- Wait, current is CV < 0.05 (remove low association). 
                                           # If we augment CV threshold (e.g. < 0.10), we remove MORE variants (including potentially TPs).
VAF_RANGE = [0.15, 0.20, 0.24, 0.30]       # Higher = More strict (remove more variants with higher VAF). 

# Special mode: "No VAF check" (always remove if AD/CV met) -> simulated by VAF < 1.0

def parse_vcf_features(vcf_path):
    """Extract VAF. Returns dict: (chrom, pos) -> vaf"""
    features = {}
    if not os.path.exists(vcf_path):
        return features
    
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            
            # Extract VAF
            vaf = 0.0
            fmt_cols = parts[8].split(":")
            s_cols = parts[9].split(":")
            data = dict(zip(fmt_cols, s_cols))
            
            if "VAF" in data:
                try: vaf = float(data["VAF"])
                except: pass
            elif "AF" in data:
                try: vaf = float(data["AF"])
                except: pass
            
            features[(chrom, pos)] = vaf
    return features

def read_summary(csv_path):
    rows = []
    if not os.path.exists(csv_path):
        return rows
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Pre-parse needed fields for speed
            try:
                row["_AD"] = abs(float(row.get("AlleleDelta", 0)))
                row["_CV"] = float(row.get("CramersV", 0))
                row["_Chr"] = row.get("Chr", "")
                row["_Pos"] = int(row.get("Pos", 0))
                rows.append(row)
            except:
                pass
    return rows

def load_sample_data(sample_name, base_dir):
    print(f"Loading data for {sample_name}...")
    metrics_path = os.path.join(base_dir, "metrics.json")
    with open(metrics_path, "r") as f:
        metrics = json.load(f)
    
    truth_total = metrics["truth_total"]
    baseline_tp = metrics["baseline"]["tp"] # This is AFTER longphase-s, BEFORE any intersubmod filter. 
    # Wait, metrics.json baseline is derived from input VCF. For pileup, it's correct.
    baseline_fp = metrics["baseline"]["fp"]

    # Paths
    tp_csv = os.path.join(base_dir, "intersubmod_tp", "significance_summary.csv")
    fp_csv = os.path.join(base_dir, "intersubmod_fp", "significance_summary.csv")
    tp_vcf = os.path.join(base_dir, "longphase_s", "filtered_snv_tp.vcf.gz")
    fp_vcf = os.path.join(base_dir, "longphase_s", "filtered_snv_fp.vcf.gz")
    
    tp_rows = read_summary(tp_csv)
    fp_rows = read_summary(fp_csv)
    tp_vaf = parse_vcf_features(tp_vcf)
    fp_vaf = parse_vcf_features(fp_vcf)

    # Attach VAF to rows
    for row in tp_rows:
        row["_VAF"] = tp_vaf.get((row["_Chr"], row["_Pos"]), 0.0)
    for row in fp_rows:
        row["_VAF"] = fp_vaf.get((row["_Chr"], row["_Pos"]), 0.0)
        
    return {
        "name": sample_name,
        "truth": truth_total,
        "base_tp": baseline_tp,
        "base_fp": baseline_fp,
        "tp_rows": tp_rows,
        "fp_rows": fp_rows
    }

def evaluate_filter(data, ad_thresh, cv_thresh, vaf_thresh):
    # Filter Logic: Remove if (AD > thresh) AND (CV < thresh) AND (VAF < thresh)
    # Rationale: High AD (methylation diff) but Low CV (weak association) -> suspicious, likely artifacts. 
    # Low VAF check reinforces that these artifacts usually appear at lower VAFs.
    
    tp_rem = 0
    fp_rem = 0
    
    for row in data["tp_rows"]:
        if (row["_AD"] > ad_thresh) and (row["_CV"] < cv_thresh) and (row["_VAF"] < vaf_thresh):
            tp_rem += 1
            
    for row in data["fp_rows"]:
        if (row["_AD"] > ad_thresh) and (row["_CV"] < cv_thresh) and (row["_VAF"] < vaf_thresh):
            fp_rem += 1
            
    final_tp = data["base_tp"] - tp_rem
    final_fp = data["base_fp"] - fp_rem
    
    precision = final_tp / (final_tp + final_fp) if (final_tp + final_fp) > 0 else 0.0
    recall = final_tp / data["truth"] if data["truth"] > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    
    return f1, tp_rem, fp_rem

def optimize():
    # Load all data first
    all_data = []
    for s, p in SAMPLES.items():
        all_data.append(load_sample_data(s, p))
        
    print(f"Data loaded. Starting grid search...")
    
    best_global_params = None
    best_global_avg_f1 = -1
    
    # Store results per sample to find per-sample best
    sample_best = {d["name"]: {"f1": -1, "params": None} for d in all_data}
    
    iterations = 0
    
    print(f"{'AD':<6} {'CV':<6} {'VAF':<6} | {'Avg F1':<8} | {'Details (Diff from Base F1)'}")
    print("-" * 80)

    # Baseline F1s
    base_f1s = {}
    for d in all_data:
        f1, _, _ = evaluate_filter(d, 1.0, 0.0, 0.0) # Always false filter
        base_f1s[d["name"]] = f1

    for ad in AD_RANGE:
        for cv in CV_RANGE:
            for vaf in VAF_RANGE:
                iterations += 1
                
                current_f1s = []
                f1_deltas = []
                
                for d in all_data:
                    f1, tp_rem, fp_rem = evaluate_filter(d, ad, cv, vaf)
                    current_f1s.append(f1)
                    f1_deltas.append(f1 - base_f1s[d["name"]])
                    
                    if f1 > sample_best[d["name"]]["f1"]:
                        sample_best[d["name"]]["f1"] = f1
                        sample_best[d["name"]]["params"] = (ad, cv, vaf)
                
                avg_f1 = np.mean(current_f1s)
                
                if avg_f1 > best_global_avg_f1:
                    best_global_avg_f1 = avg_f1
                    best_global_params = (ad, cv, vaf)
                    
                # Print progress for some interesting combinations or periodically
                # print(f"{ad:<6} {cv:<6} {vaf:<6} | {avg_f1:.4f}   | {[round(x, 4) for x in f1_deltas]}")

    print("-" * 80)
    print(f"Optimization Complete.")
    print(f"Best Limit (Avg F1): AD>{best_global_params[0]}, CV<{best_global_params[1]}, VAF<{best_global_params[2]} -> F1={best_global_avg_f1:.4f}")
    
    print("\nPer-Sample Best:")
    for name, res in sample_best.items():
        best_p = res["params"]
        base = base_f1s[name]
        impr = res["f1"] - base
        print(f"{name:<15}: AD>{best_p[0]}, CV<{best_p[1]}, VAF<{best_p[2]} -> F1={res['f1']:.4f} (Base {base:.4f}, Delta {impr:+.4f})")
        
    print("\nComparison with Default (AD>0.25, CV<0.05, VAF<0.24):")
    default_f1s = []
    for d in all_data:
        f1, _, _ = evaluate_filter(d, 0.25, 0.05, 0.24)
        default_f1s.append(f1)
    
    print(f"Default Avg F1: {np.mean(default_f1s):.4f}")
    print(f"Optimized Avg F1: {best_global_avg_f1:.4f}")

if __name__ == "__main__":
    optimize()
