"""想法2 (子任務 #6-#8): modkit/DSS 逐 CpG 率差 vs 我方 read-distance, 同一 ISM 矩陣 apples-to-apples.
在 subclone 軸 (HP1 germline vs HP1-1 carrier) 上, 四法吃同一 56-read 矩陣:
  (1) modkit-style: per-CpG 二值化(β>0.5) Fisher 率差 (= modkit/PerCpgAsm 邊際率差精神)
  (2) DSS: beta-binomial (counts 表另跑 audit_B_dss.R)
  (3) 我方連續: per-CpG MWU (想法1 同口徑)
  (4) read-distance: region 級 PERMANOVA (significance.json) — haplotype-coherence
輸出 per-CpG sig 集合 Venn + 位點級 2×2。read-set 一致(同 56 reads), 隔離「統計量」差異。
真 modkit 工具(pooled 全 reads)的 region 級結果引既有 07 (BRCA2 somatic −0.159 vs ISM −0.122)。
"""
import json
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests

import lib_region as L

OUT = "docs/experiments/in_progress/2026/06/20260620_subclone_label_methods_pilot/_assets/data"
MIN_G = 3
os.makedirs(OUT, exist_ok=True)

PILOT = {22305: ("chr13", 32315128), 22306: ("chr13", 32317522)}  # canonical + strong


def per_cpg_compare(rid, chrom, pos):
    rdir = L.find_region_dir(chrom, pos)
    R = L.load_region(rdir)
    meth, reads = R["meth"], R["reads"]
    t = reads[reads.is_tumor == 1]
    g_germ = t[t.hp == "1"].read_id.tolist()       # HP1 germline
    g_carr = t[t.hp == "1-1"].read_id.tolist()      # HP1-1 carrier
    rows = []
    for cpg in meth.columns:
        a = meth[cpg].reindex(g_germ).dropna()
        b = meth[cpg].reindex(g_carr).dropna()
        if len(a) < MIN_G or len(b) < MIN_G:
            continue
        # (1) modkit-style 二值化 Fisher 率差
        a_m, b_m = int((a > 0.5).sum()), int((b > 0.5).sum())
        a_n, b_n = len(a), len(b)
        try:
            fish_p = fisher_exact([[a_m, a_n - a_m], [b_m, b_n - b_m]])[1]
        except ValueError:
            fish_p = np.nan
        rate_delta = a_m / a_n - b_m / b_n
        # (3) 我方連續 MWU
        try:
            mwu_p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError:
            mwu_p = np.nan
        rows.append(dict(region_id=rid, cpg_pos=cpg, germ_n=a_n, carr_n=b_n,
                         germ_meth=a_m, carr_meth=b_m, rate_delta=rate_delta,
                         beta_delta=float(a.mean() - b.mean()),
                         modkit_fisher_p=fish_p, our_mwu_p=mwu_p))
    df = pd.DataFrame(rows)
    for col, q in [("modkit_fisher_p", "modkit_fisher_q"), ("our_mwu_p", "our_mwu_q")]:
        m = df[col].notna()
        df[q] = np.nan
        df.loc[m, q] = multipletests(df.loc[m, col], method="fdr_bh")[1]
    return df, rdir


# read-distance structure (region 級 PERMANOVA)
def region_structure(rdir):
    sig = json.load(open(f"{rdir}/clustering/significance.json"))
    ls = sig.get("label_structure", {})
    return dict(hp_permanova_p=ls.get("hp_permanova_p"), allele_permanova_p=ls.get("allele_permanova_p"),
                cluster_permanova_p=sig.get("cluster_structure", {}).get("permanova_p"))


all_cpg, dss_tables, summ = [], [], []
for rid, (chrom, pos) in PILOT.items():
    df, rdir = per_cpg_compare(rid, chrom, pos)
    struct = region_structure(rdir)
    all_cpg.append(df)
    # DSS counts table (chr pos hp1_N hp1_X hp2_N hp2_X)
    dss = df[["cpg_pos", "germ_n", "germ_meth", "carr_n", "carr_meth"]].copy()
    dss.columns = ["pos", "hp1_N", "hp1_X", "hp2_N", "hp2_X"]
    dss.insert(0, "chr", chrom)
    dss.to_csv(f"{OUT}/idea2_dss_counts_{rid}.tsv", sep="\t", index=False)
    n_fish = int((df.modkit_fisher_q < 0.05).sum())
    n_mwu = int((df.our_mwu_q < 0.05).sum())
    both = int(((df.modkit_fisher_q < 0.05) & (df.our_mwu_q < 0.05)).sum())
    summ.append(dict(region_id=rid, snv=f"{chrom}:{pos}", n_cpg_tested=len(df),
                     n_modkit_fisher_sig=n_fish, n_our_mwu_sig=n_mwu, n_both_sig=both,
                     structure=struct))
    print(f"[想法2] {rid} {chrom}:{pos}: {len(df)} CpG tested | modkit-Fisher sig={n_fish} "
          f"our-MWU sig={n_mwu} both={both} | read-dist HP-PERMANOVA p={struct['hp_permanova_p']}")

pd.concat(all_cpg, ignore_index=True).to_csv(f"{OUT}/idea2_percpg_compare.tsv", sep="\t", index=False, float_format="%.4g")
json.dump(summ, open(f"{OUT}/idea2_summary.json", "w"), indent=1)
print("\nWROTE idea2_percpg_compare.tsv + dss_counts tables + summary")
