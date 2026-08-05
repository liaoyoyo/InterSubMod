"""想法3 (子任務 #9-#11): per-label 平均 Δβ + 閾值 sweep + 與「距離依標籤顯著」集合 Venn.
A) pilot (5 BRCA2 region): per-read mean β -> per-label mean -> Δβ + permutation p (HP/allele/subclone 軸),
   cross-check vs summary CSV 既有欄。
B) full-scope (29,754 region): 直接從 summary CSV 算 set A(Δβ 法命中) vs set B(label-PERMANOVA 顯著) Venn
   + 閾值 sweep + 差集案例 + LEVEL-shift confound。read-only，非重跑。
characterization 非判別；與 tumor-only 非監督 NEGATIVE 區隔。
"""
import json
import os
import numpy as np
import pandas as pd

import lib_region as L

OUT = "docs/experiments/in_progress/2026/06/20260620_subclone_label_methods_pilot/_assets/data"
SUMMARY = "output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/intersubmod_tp/significance_summary.csv"
MIN_G = 3
NPERM = 999
rng = np.random.default_rng(42)
os.makedirs(OUT, exist_ok=True)


# ---------- A) pilot: Python per-label Δβ (per-read mean) ----------
def read_means(meth):
    return meth.mean(axis=1, skipna=True)  # per-read mean β over region CpGs


def perm_dbeta(rmeans, g0, g1):
    a = rmeans.reindex(g0).dropna()
    b = rmeans.reindex(g1).dropna()
    if len(a) < MIN_G or len(b) < MIN_G:
        return np.nan, np.nan, len(a), len(b)
    obs = a.mean() - b.mean()
    pool = np.concatenate([a.values, b.values])
    na = len(a)
    cnt = 0
    for _ in range(NPERM):
        rng.shuffle(pool)
        d = pool[:na].mean() - pool[na:].mean()
        if abs(d) >= abs(obs):
            cnt += 1
    p = (cnt + 1) / (NPERM + 1)
    return float(obs), float(p), len(a), len(b)


pilot = []
for rid, (chrom, pos) in L.BRCA2_REGIONS.items():
    rdir = L.find_region_dir(chrom, pos)
    R = L.load_region(rdir)
    rmeans = read_means(R["meth"])
    G = L.tumor_axis_groups(R["reads"])
    row = {"region_id": rid, "snv": f"{chrom}:{R['snv_pos']}"}
    for ax in ["hpfam", "subclone", "allele"]:
        (g0n, g0), (g1n, g1) = G[ax]
        d, p, n0, n1 = perm_dbeta(rmeans, g0, g1)
        row[f"{ax}_dbeta"] = d
        row[f"{ax}_p"] = p
        row[f"{ax}_n0"], row[f"{ax}_n1"] = n0, n1
    pilot.append(row)
    print(f"[pilot 想法3] {rid} {row['snv']}: HPΔβ={row['hpfam_dbeta']:.3f}(p={row['hpfam_p']:.3f}) "
          f"subΔβ={row['subclone_dbeta']:.3f}(p={row['subclone_p']:.3f}) "
          f"alleleΔβ={row['allele_dbeta']:.3f}(p={row['allele_p']:.3f})")
pilot_df = pd.DataFrame(pilot)
pilot_df.to_csv(f"{OUT}/idea3_pilot_dbeta.tsv", sep="\t", index=False, float_format="%.4g")

# cross-check vs summary CSV (HPMergedDelta etc.)
S = pd.read_csv(SUMMARY)
for c in ["LabelHPPermanovaValid", "LabelAllelePermanovaValid", "HPMergedSig", "AlleleSig", "HPFineSig"]:
    if c in S.columns:
        S[c] = S[c].astype(str).str.lower().isin(["true", "1", "1.0"])
cc = S[S.RegionID.isin(L.BRCA2_REGIONS.keys())][
    ["RegionID", "HPMergedDelta", "HPMergedP", "AlleleDelta", "AlleleP", "HPFineD_HP1_HP1S",
     "LabelHPPermanovaP", "LabelAllelePermanovaP"]].set_index("RegionID")
print("\n[cross-check 想法3 pilot vs summary CSV]")
for rid in L.BRCA2_REGIONS:
    py = pilot_df[pilot_df.region_id == rid].iloc[0]
    cs = cc.loc[rid]
    print(f"  {rid}: Python HPΔβ={py['hpfam_dbeta']:.3f} vs CSV HPMergedDelta={cs['HPMergedDelta']:.3f} | "
          f"Python subΔβ={py['subclone_dbeta']:.3f} vs CSV HPFineD_HP1_HP1S={cs['HPFineD_HP1_HP1S']:.3f}")
cc.to_csv(f"{OUT}/idea3_pilot_crosscheck_csv.tsv", sep="\t", float_format="%.4g")


# ---------- B) full-scope Venn (29,754 region, 從 summary CSV) ----------
def to_bool(s):
    return s.astype(str).str.lower().isin(["true", "1", "1.0"])


setB_hp = to_bool(S["LabelHPPermanovaValid"]) & (S["LabelHPPermanovaP"] <= 0.05)
setB_al = to_bool(S["LabelAllelePermanovaValid"]) & (S["LabelAllelePermanovaP"] <= 0.05)
setB = setB_hp | setB_al  # 距離依標籤顯著 (label-PERMANOVA)

# set B 另一定義 (對齊 / GlobalP) 供並排
setB_global = (S["GlobalP"] <= 0.05) | (S["GlobalP_HPFamily"] <= 0.05) | (S["GlobalP_HPFine"] <= 0.05)

THRESH = [0.10, 0.15, 0.20, 0.25]
N = len(S)


def venn(A, B):
    inter = int((A & B).sum())
    return dict(A=int(A.sum()), B=int(B.sum()), inter=inter,
                A_minus_B=int((A & ~B).sum()), B_minus_A=int((~A & B).sum()),
                jaccard=round(inter / max(1, int((A | B).sum())), 4),
                A_subset_B=bool((A & ~B).sum() == 0), B_subset_A=bool((~A & B).sum() == 0))


sweep = {}
for tau in THRESH:
    hp_hit = (S["HPMergedDelta"].abs() > tau) & (S["HPMergedP"] <= 0.05)
    al_hit = (S["AlleleDelta"].abs() > tau) & (S["AlleleP"] <= 0.05)
    sub_hit = (S["HPFineD_HP1_HP1S"].abs() > tau) & (S["HPFineP"] <= 0.05)
    setA = hp_hit | al_hit  # 主 Venn 用 HP+allele (各有自己 P)
    sweep[f"tau_{tau}"] = {
        "setA_HP_allele_vs_setB_labelPERMANOVA": venn(setA, setB),
        "setA_vs_setB_GlobalP": venn(setA, setB_global),
        "n_HP_hit": int(hp_hit.sum()), "n_allele_hit": int(al_hit.sum()),
        "n_subclone_hit": int(sub_hit.sum()),
    }
    print(f"[full-scope τ={tau}] |A|={int(setA.sum())} |B|={int(setB.sum())} "
          f"A∩B={int((setA&setB).sum())} A−B={int((setA&~setB).sum())} B−A={int((~setA&setB).sum())} "
          f"Jaccard={venn(setA,setB)['jaccard']}")

# 差集案例 (τ=0.15 為代表)
tau = 0.15
hp_hit = (S["HPMergedDelta"].abs() > tau) & (S["HPMergedP"] <= 0.05)
al_hit = (S["AlleleDelta"].abs() > tau) & (S["AlleleP"] <= 0.05)
setA = hp_hit | al_hit
S["_maxdelta"] = S[["HPMergedDelta", "AlleleDelta"]].abs().max(axis=1)
S["_maxF"] = S[["LabelHPPermanovaF", "LabelAllelePermanovaF"]].max(axis=1)
AmB = S[setA & ~setB].nlargest(15, "_maxdelta")[
    ["RegionID", "Chr", "Pos", "HPMergedDelta", "AlleleDelta", "Fisher_N_Sig", "Potential_LOH"]]
BmA = S[~setA & setB].nlargest(15, "_maxF")[
    ["RegionID", "Chr", "Pos", "LabelHPPermanovaF", "LabelAllelePermanovaF", "HPMergedDelta", "Fisher_N_Sig"]]
AmB.to_csv(f"{OUT}/idea3_diffcase_A_minus_B.tsv", sep="\t", index=False, float_format="%.4g")
BmA.to_csv(f"{OUT}/idea3_diffcase_B_minus_A.tsv", sep="\t", index=False, float_format="%.4g")

# confound: LEVEL-shift = HP Δβ 顯著 (|Δβ|>0.15 & p<=0.05) 但 per-CpG driver=0 (Fisher_N_Sig==0)
hp_dbeta_sig = (S["HPMergedDelta"].abs() > 0.15) & (S["HPMergedP"] <= 0.05)
level_shift = hp_dbeta_sig & (S["Fisher_N_Sig"] == 0)
confound = dict(
    n_HP_dbeta_sig=int(hp_dbeta_sig.sum()),
    n_with_percpg_driver=int((hp_dbeta_sig & (S["Fisher_N_Sig"] > 0)).sum()),
    n_LEVEL_shift_no_driver=int(level_shift.sum()),
    level_shift_frac=round(int(level_shift.sum()) / max(1, int(hp_dbeta_sig.sum())), 4),
)
print(f"\n[confound LEVEL-shift] HP-Δβ-sig={confound['n_HP_dbeta_sig']} "
      f"有per-CpG driver={confound['n_with_percpg_driver']} "
      f"LEVEL-shift(無driver)={confound['n_LEVEL_shift_no_driver']} ({confound['level_shift_frac']*100:.1f}%)")

out = dict(N_regions=N, setB_labelPERMANOVA_size=int(setB.sum()),
           setB_global_size=int(setB_global.sum()),
           sweep=sweep, confound_level_shift=confound)
json.dump(out, open(f"{OUT}/idea3_fullscope_venn.json", "w"), indent=1)
print("\nWROTE idea3_fullscope_venn.json + pilot/diffcase tsv")
