"""想法2 #8: join modkit-Fisher + DSS + 我方 MWU per-CpG -> 三法 Venn + 位點級 2x2."""
import json
import os
import pandas as pd

OUT = "docs/experiments/in_progress/2026/06/20260620_subclone_label_methods_pilot/_assets/data"
cmp = pd.read_csv(f"{OUT}/idea2_percpg_compare.tsv", sep="\t")
summ = json.load(open(f"{OUT}/idea2_summary.json"))
struct = {s["region_id"]: s["structure"] for s in summ}

out = {}
for rid in [22305, 22306]:
    d = cmp[cmp.region_id == rid].copy()
    dss = pd.read_csv(f"{OUT}/idea2_dss_out_{rid}.tsv", sep="\t")  # chr pos mu1 mu2 diff stat pval fdr
    d = d.merge(dss[["pos", "fdr"]].rename(columns={"pos": "cpg_pos", "fdr": "dss_fdr"}), on="cpg_pos", how="left")
    sig_mod = set(d.loc[d.modkit_fisher_q < 0.05, "cpg_pos"])
    sig_dss = set(d.loc[d.dss_fdr < 0.05, "cpg_pos"])
    sig_our = set(d.loc[d.our_mwu_q < 0.05, "cpg_pos"])
    allsig = sig_mod | sig_dss | sig_our
    hp_p = struct[rid]["hp_permanova_p"]
    out[rid] = dict(
        n_cpg=len(d),
        n_modkit_fisher=len(sig_mod), n_dss=len(sig_dss), n_our_mwu=len(sig_our),
        all3_overlap=len(sig_mod & sig_dss & sig_our),
        modkit_only=len(sig_mod - sig_dss - sig_our),
        dss_only=len(sig_dss - sig_mod - sig_our),
        our_only=len(sig_our - sig_mod - sig_dss),
        union_any=len(allsig),
        jaccard_modkit_our=round(len(sig_mod & sig_our) / max(1, len(sig_mod | sig_our)), 3),
        jaccard_dss_our=round(len(sig_dss & sig_our) / max(1, len(sig_dss | sig_our)), 3),
        # 位點級 2x2: read-distance 結構顯著 (PERMANOVA p<=0.05) × 有 >=1 sig CpG
        locus_structure_sig=bool(hp_p is not None and hp_p <= 0.05),
        locus_has_sig_cpg=len(allsig) > 0,
        hp_permanova_p=hp_p,
    )
    print(f"[想法2 join] {rid}: modkit={len(sig_mod)} DSS={len(sig_dss)} our={len(sig_our)} "
          f"| all3={out[rid]['all3_overlap']} dss_only={out[rid]['dss_only']} our_only={out[rid]['our_only']} "
          f"| Jaccard(modkit,our)={out[rid]['jaccard_modkit_our']} (dss,our)={out[rid]['jaccard_dss_our']} "
          f"| struct_sig={out[rid]['locus_structure_sig']} has_sig_cpg={out[rid]['locus_has_sig_cpg']}")

json.dump(out, open(f"{OUT}/idea2_concordance.json", "w"), indent=1)
print("\nWROTE idea2_concordance.json")
