"""補強 B.1-4 + 調查 COLO829/H2009 異常"""
import pandas as pd
import numpy as np
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod"
TSV = f"{ROOT}/output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz"

df = pd.read_csv(TSV, sep="\t",
    usecols=["sample","mode","truth_label","NumReads","HPFineNGroups","to_loh_bed_hit","caller_af"])
df["is_tp"] = (df["truth_label"]=="TP").astype(int)
df["is_loh"] = df["to_loh_bed_hit"].astype(str).str.lower().isin(["true","1"])

# B.1-4 補強: LOH × NGroups × all NR (no NR filter) TP rate
print("=== B.1-4 LOH × NGroups TP rate (TO mode, all NR) ===")
to = df[df["mode"]=="to"]
for loh in [False, True]:
    for ng in [1,2,3,4]:
        s = to[(to["is_loh"]==loh) & (to["HPFineNGroups"]==ng)]
        rate = s["is_tp"].mean() if len(s) else np.nan
        print(f"LOH={loh} NGroups={ng}: n={len(s):>7}  TPrate={rate:.4f}")
print()

# 調查 COLO829 vs H2009 為何特殊
print("=== COLO829 TO NonLOH NumReads & AF distribution ===")
c = to[(to["sample"]=="COLO829") & (to["is_loh"]==False)]
print(f"Total: {len(c)}")
print(f"NR>=80: {(c['NumReads']>=80).sum()}")
print("NR histogram:", c["NumReads"].describe().round(1).to_dict())
# 各 NGroups × TP/FP
print("\nCOLO829 NonLOH TP/FP × NGroups:")
xt = pd.crosstab([c["HPFineNGroups"]], c["truth_label"], margins=True)
print(xt)

print("\n=== H2009 TO NonLOH TP/FP × NGroups ===")
h = to[(to["sample"]=="H2009") & (to["is_loh"]==False)]
print(f"Total: {len(h)}, TP rate overall: {h['is_tp'].mean():.4f}")
xt = pd.crosstab([h["HPFineNGroups"]], h["truth_label"], margins=True)
print(xt)
# H2009 baseline TP rate 已極高 -> ceiling effect
print("\nH2009 NR≥80 NonLOH baseline TP rate:", h[h["NumReads"]>=80]["is_tp"].mean())

# residualized on NR bin 後 effect size
print("\n=== residualized (NR 20pp bin) per-sample delta N4-N3 (TO NonLOH) ===")
nr_edges = [20,40,60,80,100,150,500]
nr_labels = [f"{a}-{b}" for a,b in zip(nr_edges[:-1], nr_edges[1:])]
for sample in sorted(df["sample"].unique()):
    s = to[(to["sample"]==sample) & (to["is_loh"]==False)]
    weighted_delta = []
    for a, b, label in zip(nr_edges[:-1], nr_edges[1:], nr_labels):
        sub = s[(s["NumReads"]>=a) & (s["NumReads"]<b)]
        n3 = sub[sub["HPFineNGroups"]==3]
        n4 = sub[sub["HPFineNGroups"]==4]
        if len(n3)>=20 and len(n4)>=20:
            w = len(n4)
            d = n4["is_tp"].mean() - n3["is_tp"].mean()
            weighted_delta.append((w, d))
    if weighted_delta:
        total_w = sum(w for w,_ in weighted_delta)
        wd = sum(w*d for w,d in weighted_delta)/total_w
        print(f"{sample}: weighted delta = {wd*100:+.2f}pp (across {len(weighted_delta)} bins)")
    else:
        print(f"{sample}: insufficient")
