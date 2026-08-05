"""後處理驗證: 用既有 per-region 距離矩陣 + a-priori 標籤，以 skbio 算真 PERMANOVA + PERMDISP。
目的: (1) 重現 pipeline 的 label-PERMANOVA p (驗證對齊正確) (2) 算真 PERMDISP p
      (3) 比對 CSV 的 stub 1.0，確認散度未被計算 + 量化 location-clean vs dispersion-confounded。
路 A = 在既有距離矩陣後處理，與 C++ 用同一矩陣同一標籤 → 等價且免重跑。
"""
import glob
import json
import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.stats.distance import permanova, permdisp

RUN_TP = "output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/intersubmod_tp"
REGIONS = {22305: ("chr13", 32315128), 22306: ("chr13", 32317522), 22307: ("chr13", 32324831),
           22308: ("chr13", 32339132), 22309: ("chr13", 32345630)}
NPERM = 199


def load(region_dir):
    M = pd.read_csv(f"{region_dir}/distance/BERNOULLI/matrix.csv", index_col=0)
    M.index = M.index.astype(int)
    M.columns = M.columns.astype(int)
    reads = pd.read_csv(f"{region_dir}/reads/reads.tsv", sep="\t")
    reads["hp"] = reads["hp"].astype(str)
    return M, reads


def run_axis(M, ids, labels):
    """ids: list of read_id; labels: matching group labels. 回 (permanova_p, permanova_F, permdisp_p, permdisp_F, n0, n1)."""
    keep = [i for i in ids if i in M.index]
    sub = M.loc[keep, keep].values
    sub = (sub + sub.T) / 2.0
    np.fill_diagonal(sub, 0.0)
    lab = [labels[i] for i in keep]
    from collections import Counter
    cnt = Counter(lab)
    if len(cnt) < 2 or min(cnt.values()) < 3:
        return None
    dm = DistanceMatrix(sub, ids=[str(i) for i in keep])
    pa = permanova(dm, lab, permutations=NPERM)
    pd_ = permdisp(dm, lab, permutations=NPERM)
    return dict(permanova_p=float(pa["p-value"]), permanova_F=float(pa["test statistic"]),
                permdisp_p=float(pd_["p-value"]), permdisp_F=float(pd_["test statistic"]),
                groups=dict(cnt))


# CSV stub 對照
S = pd.read_csv(f"{RUN_TP}/significance_summary.csv")

rows = []
for rid, (chrom, pos) in REGIONS.items():
    rdir = glob.glob(f"{RUN_TP}/filtered_snv_tp/{chrom}/{chrom}_{pos}/{chrom}_*")[0]
    M, reads = load(rdir)
    t = reads[reads.is_tumor == 1]
    # HP-family 軸: HP1{1,1-1} vs HP2{2}
    hp_ids = t[t.hp.isin(["1", "1-1", "2"])].read_id.tolist()
    hp_lab = {r.read_id: ("HP1" if r.hp in ("1", "1-1") else "HP2") for r in t.itertuples()}
    hp = run_axis(M, hp_ids, hp_lab)
    # allele 軸: ALT vs REF
    al_ids = t[t.alt_support.isin(["ALT", "REF"])].read_id.tolist()
    al_lab = {r.read_id: r.alt_support for r in t.itertuples()}
    al = run_axis(M, al_ids, al_lab)
    csv = S[S.RegionID == rid].iloc[0]
    rec = dict(region_id=rid, snv=f"{chrom}:{pos}",
               csv_LabelHPPermanovaP=float(csv["LabelHPPermanovaP"]),
               csv_LabelHPDispersionP=float(csv["LabelHPDispersionP"]),
               hp=hp, allele=al)
    rows.append(rec)
    h = hp or {}
    print(f"[{rid} {chrom}:{pos}] HP軸: 我算 PERMANOVA p={h.get('permanova_p')} (CSV={csv['LabelHPPermanovaP']}) "
          f"| 我算 PERMDISP p={h.get('permdisp_p')} (CSV stub={csv['LabelHPDispersionP']}) groups={h.get('groups')}")
    a = al or {}
    print(f"        allele軸: PERMANOVA p={a.get('permanova_p')} | PERMDISP p={a.get('permdisp_p')}")

OUT = "docs/experiments/in_progress/2026/06/20260620_subclone_label_methods_pilot/_assets/data"
json.dump(rows, open(f"{OUT}/permdisp_check_brca2.json", "w"), indent=1, default=str)
print(f"\nWROTE {OUT}/permdisp_check_brca2.json")
print("\n=== 判讀 ===")
for r in rows:
    h = r["hp"]
    if not h:
        continue
    loc_clean = h["permanova_p"] <= 0.05 and h["permdisp_p"] > 0.05
    tag = "location-clean(真差異)" if loc_clean else ("dispersion-confounded" if h["permanova_p"] <= 0.05 and h["permdisp_p"] <= 0.05 else "n.s.")
    print(f"  {r['region_id']} HP軸: PERMANOVA p={h['permanova_p']:.3f} / PERMDISP p={h['permdisp_p']:.3f} → {tag}")
