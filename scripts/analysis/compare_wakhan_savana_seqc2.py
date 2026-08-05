#!/usr/bin/env python3
"""
compare_wakhan_savana_seqc2.py — 三方 CN 交叉比對 (10kb bin, chr1-22)
  Wakhan rank-1 (HP1/HP2 integer CN + loh_regions) vs SAVANA cna_normalhet vs SEQC2 truth
LOH: Wakhan=loh_regions.bed ; SAVANA=minorAlleleCopyNumber<0.5 ; SEQC2=loh segment
total CN: Wakhan=HP1_CN+HP2_CN ; SAVANA=copyNumber
"""
import pandas as pd, numpy as np
from collections import defaultdict
import scipy.stats as stats
from pathlib import Path

WO = Path("/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/wakhan_out/2.77_0.92_0.8/bed_output")
HP1 = WO / "HCC1395_2.77_0.92_0.8_copynumbers_segments_HP_1.bed"
HP2 = WO / "HCC1395_2.77_0.92_0.8_copynumbers_segments_HP_2.bed"
WLOH = WO / "loh_regions.bed"
SAV = Path("/big7_disk/liaoyoyo2001/cnv_sv_work/HCC1395/savana_wgs/cna_normalhet/HCC1395_segmented_absolute_copy_number.tsv")
SEQC2 = Path("/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed")
CH = [f"chr{i}" for i in range(1, 23)]
BIN = 10000

def ff(x):
    try: return float(x)
    except Exception: return None

def read_wakhan_cn(path):
    seg = defaultdict(list)
    for ln in open(path):
        if ln.startswith("#"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 5 or f[0] not in CH: continue
        cn = ff(f[4])
        seg[f[0]].append((int(f[1]), int(f[2]), cn))
    return seg

def read_bed_intervals(path, chrom_col=0, s=1, e=2):
    iv = defaultdict(list)
    for ln in open(path):
        if ln.startswith("#"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) <= e or f[chrom_col] not in CH: continue
        iv[f[chrom_col]].append((int(f[s]), int(f[e])))
    return iv

def binmap_cn(seg):
    """chrom -> {bin: cn}"""
    out = defaultdict(dict)
    for c, lst in seg.items():
        for s, e, cn in lst:
            if cn is None: continue
            for b in range(s // BIN, e // BIN + 1):
                out[c][b] = cn
    return out

def binset(iv):
    out = defaultdict(set)
    for c, lst in iv.items():
        for s, e in lst:
            out[c].update(range(s // BIN, e // BIN + 1))
    return out

# --- load ---
hp1 = binmap_cn(read_wakhan_cn(HP1)); hp2 = binmap_cn(read_wakhan_cn(HP2))
wloh = binset(read_bed_intervals(WLOH))
# SAVANA
sav_cn = defaultdict(dict); sav_loh = defaultdict(set); sav_cov = defaultdict(set)
for _, r in pd.read_csv(SAV, sep="\t").iterrows():
    if r["chromosome"] not in CH: continue
    cn = ff(r["copyNumber"]); mn = ff(r["minorAlleleCopyNumber"])
    for b in range(int(r["start"]) // BIN, int(r["end"]) // BIN + 1):
        sav_cov[r["chromosome"]].add(b)
        if cn is not None: sav_cn[r["chromosome"]][b] = cn
        if mn is not None and mn < 0.5: sav_loh[r["chromosome"]].add(b)
# SEQC2
seq_loh = defaultdict(set); seq_gain = defaultdict(set)
for ln in open(SEQC2):
    f = ln.rstrip().split("\t")
    if len(f) < 4 or f[0] not in CH: continue
    for b in range(int(f[1]) // BIN, int(f[2]) // BIN + 1):
        if f[3] == "loh": seq_loh[f[0]].add(b)
        elif f[3] == "gain": seq_gain[f[0]].add(b)

# --- aggregate over bins covered by all 3 (SAVANA cov ∩ Wakhan cov) ---
def jac(a, b): u = len(a | b); return len(a & b) / u if u else float("nan")
W_LOH = set(); S_LOH = set(); Q_LOH = set(); COV = set()
wtot_list = []; stot_list = []
for ci, c in enumerate(CH):
    off = ci * 10_000_000
    wcov = set(hp1.get(c, {}).keys()) & set(hp2.get(c, {}).keys())
    cov = wcov & sav_cov.get(c, set())   # bins where Wakhan(HP1&HP2) and SAVANA both定義
    for b in cov:
        COV.add(off + b)
        if b in wloh.get(c, set()): W_LOH.add(off + b)
        if b in sav_loh.get(c, set()): S_LOH.add(off + b)
        if b in seq_loh.get(c, set()): Q_LOH.add(off + b)
        wt = hp1[c][b] + hp2[c][b]; st = sav_cn[c].get(b)
        if st is not None:
            wtot_list.append(wt); stot_list.append(st)

print(f"=== 三方 CN 交叉比對 (10kb bin, chr1-22, 共同覆蓋 {len(COV)*BIN/1e6:.0f}Mb) ===")
print(f"purity/ploidy: Wakhan rank-1 = 0.92/2.77 ; SAVANA cna_normalhet = 0.96/2.79  (近乎一致)")
print(f"\nLOH 覆蓋率(共同 bin): Wakhan {len(W_LOH)/len(COV):.1%}  SAVANA {len(S_LOH)/len(COV):.1%}  SEQC2 {len(Q_LOH)/len(COV):.1%}")
print(f"\nLOH 三方 Jaccard:")
print(f"  Wakhan vs SEQC2 : {jac(W_LOH, Q_LOH):.3f}")
print(f"  SAVANA vs SEQC2 : {jac(S_LOH, Q_LOH):.3f}  (前報 0.962)")
print(f"  Wakhan vs SAVANA: {jac(W_LOH, S_LOH):.3f}")
tri = len(W_LOH & S_LOH & Q_LOH)
print(f"  三方同時 LOH bin: {tri} ; 三方 union: {len(W_LOH|S_LOH|Q_LOH)} ; 三方一致(all-agree LOH or all-not)= "
      f"{(len(COV)-len((W_LOH|S_LOH|Q_LOH) - (W_LOH&S_LOH&Q_LOH)))/len(COV):.1%}")
wt = np.array(wtot_list); st = np.array(stot_list)
rho, _ = stats.spearmanr(wt, st); r, _ = stats.pearsonr(wt, st)
print(f"\ntotal CN (Wakhan HP1+HP2 vs SAVANA copyNumber) n={len(wt)}: Spearman {rho:.3f} / Pearson {r:.3f}")
print(f"  Wakhan total CN range {wt.min():.1f}-{wt.max():.1f} median {np.median(wt):.2f} ; SAVANA {st.min():.1f}-{st.max():.1f} median {np.median(st):.2f}")
