"""(b) 全 scope 分軸 PERMDISP (analytic betadisper, vegan 式 + 非歐負特徵值校正).
對 label-significant(loc PERMANOVA p≤0.05) 的 region，分 HP 軸/allele 軸算 dispersion analytic ANOVA-F p，
與 CSV 的 location p 合併分類: location-clean(loc≤0.05 & disp>0.05) vs dispersion-confounded(loc≤0.05 & disp≤0.05)。
免 permutation 迴圈 → 可全 scope。先對 BRCA2 5 region 交叉驗證 analytic vs skbio。
"""
import glob
import json
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import f_oneway

RUN_TP = "output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/intersubmod_tp"
OUT = "docs/experiments/in_progress/2026/06/20260620_subclone_label_methods_pilot/_assets/data"
MIN_G = 3


def analytic_betadisper(D, labels):
    """vegan betadisper analytic (ANOVA-F on distance-to-centroid; 負特徵值校正)。回 (F, p)。"""
    n = len(D)
    A = -0.5 * (D ** 2)
    J = np.eye(n) - np.ones((n, n)) / n
    G = J @ A @ J
    G = (G + G.T) / 2
    ev, U = np.linalg.eigh(G)
    tol = 1e-8 * (np.abs(ev).max() if n else 1)
    pos = ev > tol
    neg = ev < -tol
    xr = U[:, pos] * np.sqrt(ev[pos]) if pos.any() else np.zeros((n, 0))
    xi = U[:, neg] * np.sqrt(-ev[neg]) if neg.any() else np.zeros((n, 0))
    labs = np.array(labels)
    z = np.zeros(n)
    for g in np.unique(labs):
        m = labs == g
        cr = xr[m].mean(0) if xr.shape[1] else np.zeros(0)
        ci = xi[m].mean(0) if xi.shape[1] else np.zeros(0)
        dr = ((xr[m] - cr) ** 2).sum(1) if xr.shape[1] else np.zeros(m.sum())
        di = ((xi[m] - ci) ** 2).sum(1) if xi.shape[1] else np.zeros(m.sum())
        z[m] = np.sqrt(np.abs(dr - di))
    groups = [z[labs == g] for g in np.unique(labs)]
    if any(len(g) < 2 for g in groups):
        return np.nan, np.nan
    F, p = f_oneway(*groups)
    return float(F), float(p)


def region_dir(chrom, pos):
    h = glob.glob(f"{RUN_TP}/filtered_snv_tp/{chrom}/{chrom}_{pos}/{chrom}_*")
    return h[0] if h else None


def load_dist_reads(rdir):
    M = pd.read_csv(f"{rdir}/distance/BERNOULLI/matrix.csv", index_col=0)
    M.index = M.index.astype(int); M.columns = M.columns.astype(int)
    reads = pd.read_csv(f"{rdir}/reads/reads.tsv", sep="\t")
    reads["hp"] = reads["hp"].astype(str)
    return M, reads


def axis_disp(M, reads, axis):
    t = reads[reads.is_tumor == 1]
    if axis == "hp":
        ids = t[t.hp.isin(["1", "1-1", "2"])].read_id.tolist()
        lab = ["HP1" if h in ("1", "1-1") else "HP2" for h in t[t.hp.isin(["1", "1-1", "2"])].hp]
    else:
        sel = t[t.alt_support.isin(["ALT", "REF"])]
        ids = sel.read_id.tolist(); lab = sel.alt_support.tolist()
    ids = [i for i in ids if i in M.index]
    lab = lab[:len(ids)]
    from collections import Counter
    c = Counter(lab)
    if len(c) < 2 or min(c.values()) < MIN_G:
        return np.nan
    sub = M.loc[ids, ids].values; sub = (sub + sub.T) / 2; np.fill_diagonal(sub, 0)
    return analytic_betadisper(sub, lab)[1]


S = pd.read_csv(f"{RUN_TP}/significance_summary.csv")
hp_loc = pd.to_numeric(S["LabelHPPermanovaP"], errors="coerce")
al_loc = pd.to_numeric(S["LabelAllelePermanovaP"], errors="coerce")

# --- 交叉驗證 analytic vs skbio (BRCA2 5 region) ---
print("=== 交叉驗證 analytic betadisper vs skbio (BRCA2) ===")
sk = {r["region_id"]: r for r in json.load(open(f"{OUT}/permdisp_check_brca2.json"))}
for rid in [22305, 22306, 22307]:
    row = S[S.RegionID == rid].iloc[0]
    M, reads = load_dist_reads(region_dir(row["Chr"], int(row["Pos"])))
    ap = axis_disp(M, reads, "hp")
    skp = sk[rid]["hp"]["permdisp_p"] if sk[rid]["hp"] else None
    print(f"  {rid} HP: analytic disp p={ap:.3f} vs skbio permdisp p={skp}")

# --- 全 scope: 只算 label-sig 的 region ---
atrisk = S[(hp_loc <= 0.05) | (al_loc <= 0.05)].copy()
print(f"\n=== 全 scope 分軸 PERMDISP: label-sig region = {len(atrisk)} ===", flush=True)
recs = []
t0 = time.time()
for n, (_, row) in enumerate(atrisk.iterrows()):
    rdir = region_dir(row["Chr"], int(row["Pos"]))
    if not rdir:
        continue
    try:
        M, reads = load_dist_reads(rdir)
        hp_d = axis_disp(M, reads, "hp") if hp_loc[row.name] <= 0.05 else np.nan
        al_d = axis_disp(M, reads, "allele") if al_loc[row.name] <= 0.05 else np.nan
    except Exception:
        hp_d = al_d = np.nan
    recs.append((row["RegionID"], hp_loc[row.name], hp_d, al_loc[row.name], al_d))
    if (n + 1) % 3000 == 0:
        print(f"  ...{n+1}/{len(atrisk)} ({time.time()-t0:.0f}s)", flush=True)

R = pd.DataFrame(recs, columns=["RegionID", "hp_loc_p", "hp_disp_p", "al_loc_p", "al_disp_p"])
R.to_csv(f"{OUT}/fullscope_permdisp.tsv", sep="\t", index=False, float_format="%.4g")


def classify(loc, disp):
    sig = loc <= 0.05
    testable = sig & disp.notna()
    clean = testable & (disp > 0.05)
    conf = testable & (disp <= 0.05)
    return dict(loc_sig=int(sig.sum()), testable=int(testable.sum()),
                clean=int(clean.sum()), confounded=int(conf.sum()),
                pct_clean=round(100 * clean.sum() / max(1, testable.sum()), 1),
                pct_confounded=round(100 * conf.sum() / max(1, testable.sum()), 1))


summ = {"N_atrisk": len(R),
        "HP_axis": classify(R["hp_loc_p"], R["hp_disp_p"]),
        "allele_axis": classify(R["al_loc_p"], R["al_disp_p"])}
json.dump(summ, open(f"{OUT}/fullscope_permdisp_summary.json", "w"), indent=1)
print("\n=== 分軸 location-clean vs dispersion-confounded (全 scope) ===")
for ax in ["HP_axis", "allele_axis"]:
    s = summ[ax]
    print(f"{ax}: loc-sig={s['loc_sig']} testable={s['testable']} | "
          f"location-clean={s['clean']}({s['pct_clean']}%) dispersion-confounded={s['confounded']}({s['pct_confounded']}%)")
print(f"\n總耗時 {time.time()-t0:.0f}s; WROTE fullscope_permdisp_summary.json")
