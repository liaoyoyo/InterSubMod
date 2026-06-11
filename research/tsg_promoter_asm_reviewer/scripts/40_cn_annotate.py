#!/usr/bin/env python3
"""
Stage A — CN annotation for the HCC1395 ASM x SEQC2 CN confound pilot.

Builds a per-site master table from the O1 dual-axis ASM data (TP + FP) annotated
with SEQC2 continuous median copy number (Masood 2024 short-read NGS ground truth).

CN annotation logic (matches 18_dual_axis_pivot.py:47-52 closed-interval style):
  - Build per-chrom interval maps from:
      Additional_file_3_cnv_gain_cn_median.txt   (gains, CN>=3)
      Additional_file_4_cnv_loss_cn_median.txt   (losses, CN 0-1)
  - For each somatic_pos (closed interval s <= pos <= e):
      median_cn = gain CN if in a gain interval
                  elif in a loss interval -> loss CN
                  else 2  (diploid / cnLOH default)
  - cn_class   = 'gain' if median_cn>2, 'loss' if median_cn<2, else 'neutral'
  - cnloh_flag = 1 if (loh_status=='LOH' and median_cn==2) else 0  (copy-neutral LOH)

Output: genome_survey_v2/cn_confound/master_o1_cn.tsv
  (all original O1 cols + is_tp, abs_delta, median_cn, cn_class, cnloh_flag)

Usage: python3 40_cn_annotate.py
"""
import os
from collections import defaultdict
import pandas as pd

BASE = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
TP_TSV = os.path.join(BASE, "genome_survey_v2", "asm_dualaxis_tp.tsv")
FP_TSV = os.path.join(BASE, "genome_survey_v2", "asm_dualaxis_fp.tsv")
GAIN_TXT = "/big8_disk/data/HCC1395/SEQC2/CNV/Additional_file_3_cnv_gain_cn_median.txt"
LOSS_TXT = "/big8_disk/data/HCC1395/SEQC2/CNV/Additional_file_4_cnv_loss_cn_median.txt"
OUT_TSV = os.path.join(BASE, "genome_survey_v2", "cn_confound", "master_o1_cn.tsv")


def load_intervals(path):
    """chrom -> list of (start, end, median_CN). Closed-interval semantics."""
    m = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            p = line.split("\t")
            if len(p) < 4:
                continue
            # skip a header row if present (non-int start)
            try:
                start, end, cn = int(p[1]), int(p[2]), int(float(p[3]))
            except ValueError:
                continue
            m[p[0]].append((start, end, cn))
    # sort intervals per chrom for deterministic first-match
    for chrom in m:
        m[chrom].sort()
    return m


def lookup_cn(chrom, pos, gain_map, loss_map):
    """median_cn: gain CN if in gain interval; elif loss CN; else 2 (diploid/cnLOH)."""
    pos = int(pos)
    for s, e, cn in gain_map.get(chrom, []):
        if s <= pos <= e:
            return cn
    for s, e, cn in loss_map.get(chrom, []):
        if s <= pos <= e:
            return cn
    return 2


def main():
    # 1. Load TP/FP dual-axis ASM, tag is_tp, concat, add abs_delta
    tp = pd.read_csv(TP_TSV, sep="\t")
    fp = pd.read_csv(FP_TSV, sep="\t")
    tp["is_tp"] = 1
    fp["is_tp"] = 0
    df = pd.concat([tp, fp], ignore_index=True)
    df["abs_delta"] = df["mean_delta"].abs()
    print(f"[load] TP rows={len(tp)}  FP rows={len(fp)}  concat={len(df)}")

    # 2. Build CN interval maps
    gain_map = load_intervals(GAIN_TXT)
    loss_map = load_intervals(LOSS_TXT)
    n_gain_iv = sum(len(v) for v in gain_map.values())
    n_loss_iv = sum(len(v) for v in loss_map.values())
    print(f"[cn-maps] gain intervals={n_gain_iv} (chroms={len(gain_map)})  "
          f"loss intervals={n_loss_iv} (chroms={len(loss_map)})")

    # 3. Annotate median_cn per row, then cn_class + cnloh_flag
    df["median_cn"] = [
        lookup_cn(c, p, gain_map, loss_map)
        for c, p in zip(df["chrom"], df["somatic_pos"])
    ]
    df["cn_class"] = df["median_cn"].apply(
        lambda cn: "gain" if cn > 2 else ("loss" if cn < 2 else "neutral")
    )
    df["cnloh_flag"] = (
        (df["loh_status"] == "LOH") & (df["median_cn"] == 2)
    ).astype(int)

    # 4. Save master table (original cols + new cols)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[save] {OUT_TSV}  ({len(df)} rows, {df.shape[1]} cols)")

    # 5. Sanity
    print("\n==== SANITY ====")
    print(f"total rows: {len(df)}")
    print("by is_tp:")
    print(df["is_tp"].value_counts().sort_index().to_string())

    print("\ncn_class value_counts (overall):")
    print(df["cn_class"].value_counts().to_string())

    # HP-axis significant subset (axis_type==HP & wilcoxon_p<0.05)
    hp_sig = df[(df["axis_type"] == "HP") & (df["wilcoxon_p"] < 0.05)]
    print(f"\nHP-axis sig (axis_type==HP & p<0.05) rows: {len(hp_sig)}")
    print("cn_class value_counts (HP-axis sig):")
    print(hp_sig["cn_class"].value_counts().to_string())

    n_cnloh = int(df["cnloh_flag"].sum())
    print(f"\nn cnloh (copy-neutral LOH): {n_cnloh}")

    print("\nmedian_cn distribution (overall):")
    print(df["median_cn"].value_counts().sort_index().to_string())

    n_gt2 = int((df["median_cn"] > 2).sum())
    n_lt2 = int((df["median_cn"] < 2).sum())
    print(f"\nn median_cn>2 (gain): {n_gt2}")
    print(f"n median_cn<2 (loss): {n_lt2}")
    print(f"n median_cn==2 (neutral): {int((df['median_cn'] == 2).sum())}")

    # Return-contract compact summary (printed for parent capture)
    overall_counts = df["cn_class"].value_counts().to_dict()
    hpsig_counts = hp_sig["cn_class"].value_counts().to_dict()
    print("\n==== COMPACT ====")
    print({
        "n_total": len(df),
        "n_tp": int(tp.shape[0]),
        "n_fp": int(fp.shape[0]),
        "cn_class_overall": overall_counts,
        "cn_class_hp_sig": hpsig_counts,
        "n_cn_gt2": n_gt2,
        "n_cn_lt2": n_lt2,
        "n_cnloh": n_cnloh,
        "output_path": OUT_TSV,
    })


if __name__ == "__main__":
    main()
