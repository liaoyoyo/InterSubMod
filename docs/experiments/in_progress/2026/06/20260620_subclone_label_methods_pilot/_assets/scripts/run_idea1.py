"""想法1 (子任務 #2-#5): per-CpG × 三軸標籤關聯 + 歸屬定位表 + normal-ASM-control.
對每個 BRCA2 region 窗內每個 CpG，分別對 HP-family / HP-fine / subclone / allele 軸
跑 MWU(2 群) 或 Kruskal(多群) + 窗內各軸獨立 BH-FDR -> per_cpg_axis_attribution.tsv。
characterization 非 TP/FP 判別。§13.0 落檔→Read。
"""
import json
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests

import lib_region as L

OUT = "docs/experiments/in_progress/2026/06/20260620_subclone_label_methods_pilot/_assets/data"
MIN_G = 3  # 每組最小 read 數
os.makedirs(OUT, exist_ok=True)


def group_betas(meth_col, ids):
    v = meth_col.reindex(ids).dropna().values
    return v


def test_2group(meth_col, g0, g1):
    a, b = group_betas(meth_col, g0), group_betas(meth_col, g1)
    if len(a) < MIN_G or len(b) < MIN_G:
        return np.nan, np.nan, len(a), len(b)
    try:
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
    except ValueError:
        return np.nan, np.nan, len(a), len(b)
    return p, float(np.mean(a) - np.mean(b)), len(a), len(b)


def test_kw(meth_col, groups):
    arrs = [group_betas(meth_col, ids) for _, ids in groups]
    arrs = [a for a in arrs if len(a) >= MIN_G]
    if len(arrs) < 2:
        return np.nan, np.nan
    try:
        p = kruskal(*arrs).pvalue
    except ValueError:
        return np.nan, np.nan
    means = [np.mean(a) for a in arrs]
    return p, float(max(means) - min(means))  # max pairwise |delta| proxy


def run_region(rid, chrom, pos):
    rdir = L.find_region_dir(chrom, pos)
    R = L.load_region(rdir)
    meth, reads, snv = R["meth"], R["reads"], R["snv_pos"]
    G = L.tumor_axis_groups(reads)
    nhp1, nhp2 = L.normal_hp_groups(reads)
    normal_ok = len(nhp1) >= MIN_G and len(nhp2) >= MIN_G

    rows = []
    for cpg in meth.columns:
        col = meth[cpg]
        hpfam_p, hpfam_d, n0, n1 = test_2group(col, G["hpfam"][0][1], G["hpfam"][1][1])
        sub_p, sub_d, _, _ = test_2group(col, G["subclone"][0][1], G["subclone"][1][1])
        al_p, al_d, _, _ = test_2group(col, G["allele"][0][1], G["allele"][1][1])
        hpfine_p, hpfine_d = test_kw(col, G["hpfine"])
        # normal-anchored ASM (HP1 vs HP2 on normal reads)
        if normal_ok:
            norm_p, norm_d, _, _ = test_2group(col, nhp1, nhp2)
        else:
            norm_p, norm_d = np.nan, np.nan
        rows.append(dict(region_id=rid, snv_anchor=f"{chrom}:{snv}", cpg_pos=cpg,
                         dist_to_snv=cpg - snv,
                         hpfam_p=hpfam_p, hpfam_delta=hpfam_d,
                         hpfine_p=hpfine_p, hpfine_maxdelta=hpfine_d,
                         sub_p=sub_p, sub_delta=sub_d,
                         allele_p=al_p, allele_delta=al_d,
                         normal_hp_p=norm_p, normal_hp_delta=norm_d))
    df = pd.DataFrame(rows)
    # 窗內各軸獨立 BH-FDR
    for ax in ["hpfam", "hpfine", "sub", "allele", "normal_hp"]:
        pcol = f"{ax}_p"
        qcol = f"{ax}_q"
        mask = df[pcol].notna()
        df[qcol] = np.nan
        if mask.sum() > 0:
            df.loc[mask, qcol] = multipletests(df.loc[mask, pcol], method="fdr_bh")[1]
    # dominant axis 在 {hpfam, sub, allele} 之中 (2 群可解釋軸); 過 q<0.05
    def dom(row):
        cand = {a: row[f"{a}_q"] for a in ["hpfam", "sub", "allele"] if pd.notna(row[f"{a}_q"]) and row[f"{a}_q"] < 0.05}
        if not cand:
            return "none", 0
        best = min(cand, key=cand.get)
        return best, len(cand)
    dd = df.apply(dom, axis=1, result_type="expand")
    df["dominant_axis"] = dd[0]
    df["multi_axis_n"] = dd[1]
    # germline_ASM_confounded: normal HP1-vs-HP2 顯著 -> True; normal 不可評估 -> NA
    if normal_ok:
        df["germline_ASM_confounded"] = (df["normal_hp_q"] < 0.05)
    else:
        df["germline_ASM_confounded"] = pd.NA
    df["normal_control_evaluable"] = normal_ok
    return df, dict(region_id=rid, snv=f"{chrom}:{snv}", n_cpg=len(df),
                    n_tumor=int((reads.is_tumor == 1).sum()),
                    n_normal=int((reads.is_tumor == 0).sum()),
                    normal_hp1=len(nhp1), normal_hp2=len(nhp2), normal_control_evaluable=normal_ok,
                    sig_hpfam=int((df.hpfam_q < 0.05).sum()),
                    sig_hpfine=int((df.hpfine_q < 0.05).sum()),
                    sig_sub=int((df.sub_q < 0.05).sum()),
                    sig_allele=int((df.allele_q < 0.05).sum()),
                    dom_counts=df.dominant_axis.value_counts().to_dict())


all_df, summ = [], []
for rid, (chrom, pos) in L.BRCA2_REGIONS.items():
    df, s = run_region(rid, chrom, pos)
    all_df.append(df)
    summ.append(s)
    print(f"region {rid} {s['snv']}: {s['n_cpg']} CpG | sig hpfam={s['sig_hpfam']} "
          f"sub={s['sig_sub']} allele={s['sig_allele']} | normal_ctrl={s['normal_control_evaluable']} "
          f"| dom={s['dom_counts']}")

full = pd.concat(all_df, ignore_index=True)
full.to_csv(f"{OUT}/idea1_per_cpg_axis_attribution.tsv", sep="\t", index=False, float_format="%.4g")
json.dump(summ, open(f"{OUT}/idea1_summary.json", "w"), indent=1)
print("\nWROTE", f"{OUT}/idea1_per_cpg_axis_attribution.tsv", "rows=", len(full))

# === cross-check vs BRCA2 bespoke 真值 (region 22305, HP-axis Δβ=-0.122) ===
r0 = full[full.region_id == 22305]
hp_sig = r0[r0.hpfam_q < 0.05]
print(f"\n[cross-check 22305] HP-family sig CpG = {len(hp_sig)}/{len(r0)}; "
      f"mean hpfam_delta(sig) = {hp_sig.hpfam_delta.mean():.4f} (bespoke 全域 Δβ=-0.122, 預期負)")
# promoter -500..0bp 區 (snv=32315128) 的 HP delta 方向
prox = r0[(r0.dist_to_snv >= -800) & (r0.dist_to_snv <= 200)]
print(f"[cross-check] proximal(-800..+200bp) mean hpfam_delta = {prox.hpfam_delta.mean():.4f}, "
      f"n={len(prox)}, sig={int((prox.hpfam_q<0.05).sum())}")
json.dump(dict(hp_sig_cpg=len(hp_sig), n_cpg_22305=len(r0),
               mean_hpfam_delta_sig=float(hp_sig.hpfam_delta.mean()) if len(hp_sig) else None,
               proximal_mean_delta=float(prox.hpfam_delta.mean()),
               proximal_n=len(prox)),
          open(f"{OUT}/idea1_crosscheck_brca2.json", "w"), indent=1)
