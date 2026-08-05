#!/usr/bin/env python3
"""Per-sample subcluster 切割總帳 — 從單次 ISM scan 輸出(distance matrix + significance_summary)
複製 wg_contingency.py 的逐字分群邏輯 + permanova_clean_and_4group.py 的 Q1 Venn,
產出 A/C/Jaccard 主表 + 切割品質分帶(clear/mid/weak/vweak/degen)。

A (cansplit)  = all-read tumor 分群 best_k>=2  (wg_contingency cluster_bestk 逐字)
C (clean)     = significance_summary LabelHP 或 LabelAllele PERMANOVA p<0.05 & !DispersionWarn
分帶          = cansplit 內 cluster x label Cramer's V (依 build_spectrum_html.py 顯示定義重建)

用法: compute_split_accounting.py <SAMPLE> <ISM_OUT_DIR> <RESULT_DIR> [NPROC]
  ISM_OUT_DIR 下含 significance_summary.csv + <run>/chr*/chrN_POS/chrN_S_E/{distance/BERNOULLI/matrix.csv, reads/reads.tsv}
輸出: <RESULT_DIR>/records_wg2_<SAMPLE>.json + result_<SAMPLE>.json
"""
import os, sys, csv, glob, json, time
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency, fisher_exact
from sklearn.metrics import silhouette_score
from collections import Counter, defaultdict
from multiprocessing import Pool

# ===== 逐字複製自 wg_contingency.py =====
MIN_SZ = 3
LABMAP = {"1": "1", "HP1": "1", "1-1": "1-1", "2": "2", "HP2": "2", "2-1": "2-1"}

def _f(x):
    x = x.strip()
    return np.nan if x in ("", "NA", "nan", "NaN") else float(x)

def load_matrix(path):
    rows = open(path).read().strip().split("\n"); ids = []; M = []
    for ln in rows[1:]:
        p = ln.split(","); ids.append(p[0]); M.append([_f(x) for x in p[1:]])
    return ids, np.array(M)

def peel(sub):
    idx = list(range(sub.shape[0]))
    while True:
        S = sub[np.ix_(idx, idx)]; bad = (S < 0) | np.isnan(S); np.fill_diagonal(bad, False)
        if not bad.any():
            return idx
        idx.remove(idx[int(np.argmax(bad.sum(1)))])
        if len(idx) < MIN_SZ * 2:
            return idx

def cluster_bestk(sub):
    n = sub.shape[0]
    if n < MIN_SZ * 2:
        return None
    np.fill_diagonal(sub, 0.0); sub = np.maximum(sub, sub.T)
    try:
        Z = linkage(squareform(sub, checks=False), method="average")
    except Exception:
        return None
    best = None
    for k in range(2, min(6, n // MIN_SZ) + 1):
        lab = fcluster(Z, k, criterion="maxclust"); sz = Counter(lab)
        if len(sz) < k or min(sz.values()) < MIN_SZ:
            continue
        try:
            s = silhouette_score(sub, lab, metric="precomputed")
        except Exception:
            continue
        if best is None or s > best[2]:
            best = (k, lab, float(s))
    return best

def istum(t):
    return str(t) in ("1", "true", "True")

# ===== per-region 處理(逐字 wg_contingency 的 all + within 結構) =====
def process_region(mp):
    """回 record dict 或 None。mp = .../chrN_POS/chrN_S_E/distance/BERNOULLI/matrix.csv"""
    rdir = os.path.dirname(os.path.dirname(os.path.dirname(mp)))
    rt = os.path.join(rdir, "reads", "reads.tsv")
    if not os.path.exists(rt):
        return None
    chrom = None; pos = None
    for part in rdir.split("/"):
        if part.startswith("chr") and part.count("_") == 1:
            chrom = part.split("_")[0]; pos = part.split("_")[1]; break
    try:
        reads = list(csv.DictReader(open(rt), delimiter="\t"))
    except Exception:
        return None
    rid2 = {x["read_id"]: x for x in reads}
    try:
        ids, M = load_matrix(mp)
    except Exception:
        return None
    if M.size == 0:
        return None
    id2row = {rid: i for i, rid in enumerate(ids)}
    # ----- all-read tumor clustering + contingency -----
    all_rows = []; all_lab = []
    for rid in ids:
        if rid in rid2 and istum(rid2[rid]["is_tumor"]) and rid2[rid]["hp"] in LABMAP:
            all_rows.append(id2row[rid]); all_lab.append(LABMAP[rid2[rid]["hp"]])
    allrec = {"n": len(all_rows), "n_complete": None, "best_k": None, "best_sil": None,
              "contingency": None, "labels": dict(Counter(all_lab))}
    if len(all_rows) >= MIN_SZ * 2:
        sub = M[np.ix_(all_rows, all_rows)]; keep = peel(sub); allrec["n_complete"] = len(keep)
        if len(keep) >= MIN_SZ * 2:
            cb = cluster_bestk(sub[np.ix_(keep, keep)].copy())
            if cb:
                k, labs, sil = cb; allrec["best_k"] = int(k); allrec["best_sil"] = round(sil, 3)
                kept_lab = [all_lab[i] for i in keep]
                cont = defaultdict(lambda: defaultdict(int))
                for cl, lb in zip(labs, kept_lab):
                    cont[lb][int(cl)] += 1
                allrec["contingency"] = {lb: dict(d) for lb, d in cont.items()}
    # ----- within-label subcluster -----
    within = {}
    for lab in ("1-1", "2-1", "1", "2"):
        rows = [id2row[rid] for rid in ids
                if rid in rid2 and istum(rid2[rid]["is_tumor"]) and LABMAP.get(rid2[rid]["hp"]) == lab]
        w = {"n": len(rows), "best_k": None, "best_sil": None, "min_frac": None}
        if len(rows) >= MIN_SZ * 2:
            sub = M[np.ix_(rows, rows)]; keep = peel(sub)
            if len(keep) >= MIN_SZ * 2:
                cb = cluster_bestk(sub[np.ix_(keep, keep)].copy())
                if cb:
                    k, labs, sil = cb; sz = Counter(labs)
                    w["best_k"] = int(k); w["best_sil"] = round(sil, 3)
                    w["min_frac"] = round(min(sz.values()) / len(keep), 3)
        within[lab] = w
    return {"chrom": chrom, "pos": pos, "all": allrec, "within": within}

# ===== 切割品質分帶: cluster x label Cramer's V (依 build_spectrum_html.py 定義重建) =====
def cramers_v_and_sig(contingency):
    """contingency = {label:{cluster:count}}。回 (V, fisher_or_chi2_sig_bool, n_label, n_cluster)。
    degen = 單一 label 或 單一 cluster。"""
    labels = sorted(contingency.keys())
    clusters = sorted({c for d in contingency.values() for c in d.keys()})
    if len(labels) < 2 or len(clusters) < 2:
        return None, False, len(labels), len(clusters)  # degenerate
    T = np.zeros((len(labels), len(clusters)))
    li = {l: i for i, l in enumerate(labels)}; ci = {c: j for j, c in enumerate(clusters)}
    for l, d in contingency.items():
        for c, n in d.items():
            T[li[l], ci[c]] += n
    # 移除全零行列
    T = T[T.sum(1) > 0][:, T.sum(0) > 0]
    if min(T.shape) < 2:
        return None, False, len(labels), len(clusters)
    nT = T.sum()
    try:
        chi2, p, _, _ = chi2_contingency(T)
    except Exception:
        return None, False, len(labels), len(clusters)
    V = float(np.sqrt(chi2 / (nT * (min(T.shape) - 1))))
    # 顯著性: 2x2 用 fisher, 否則 chi2 p
    if T.shape == (2, 2):
        try:
            _, p = fisher_exact(T.astype(int))
        except Exception:
            pass
    return V, (p < 0.05), len(labels), len(clusters)

def band_of(V, sig):
    if V is None:
        return "degen"
    if V >= 0.7 and sig:
        return "clear"
    if V >= 0.5 and sig:
        return "mid"
    if V >= 0.3:
        return "weak"
    return "vweak"

# ===== C: clean PERMANOVA (逐字 permanova_clean_and_4group.py) =====
def tru(v):
    return str(v).strip().lower() in ("true", "1", "yes")

def ff(v):
    try:
        return float(v)
    except Exception:
        return None

def main():
    sample = sys.argv[1]; ism_out = sys.argv[2]; result_dir = sys.argv[3]
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else max(1, os.cpu_count() - 2)
    os.makedirs(result_dir, exist_ok=True)
    t0 = time.time()

    mats = glob.glob(f"{ism_out}/**/distance/BERNOULLI/matrix.csv", recursive=True)
    print(f"[{sample}] found {len(mats)} distance matrices, clustering with {nproc} procs...", flush=True)
    with Pool(nproc) as pool:
        records = [r for r in pool.map(process_region, mats, chunksize=64) if r is not None]
    json.dump(records, open(f"{result_dir}/records_wg2_{sample}.json", "w"))
    print(f"[{sample}] records={len(records)} elapsed={int(time.time()-t0)}s", flush=True)

    N = len(records)
    recs = {f"{r['chrom']}_{r['pos']}": r for r in records}

    # A = cansplit (best_k>=2)
    cansplit = {k for k, r in recs.items() if r["all"]["best_k"] is not None and r["all"]["best_k"] >= 2}

    # C = clean PERMANOVA (任一軸) — 讀 significance_summary.csv
    sig_path = f"{ism_out}/significance_summary.csv"
    clean = set(); hp_sig = hp_clean = al_sig = al_clean = 0; matched = 0
    if os.path.exists(sig_path):
        for row in csv.DictReader(open(sig_path)):
            key = f"{row['Chr']}_{row['Pos']}"
            if key not in recs:
                continue
            matched += 1
            hv = tru(row.get("LabelHPPermanovaValid")); hp = ff(row.get("LabelHPPermanovaP")); hw = tru(row.get("LabelHPDispersionWarn"))
            av = tru(row.get("LabelAllelePermanovaValid")); ap = ff(row.get("LabelAllelePermanovaP")); aw = tru(row.get("LabelAlleleDispersionWarn"))
            hp_c = hv and hp is not None and hp < 0.05 and not hw
            al_c = av and ap is not None and ap < 0.05 and not aw
            if hv and hp is not None and hp < 0.05:
                hp_sig += 1
            if al_c:
                al_clean += 1
            if hp_c:
                hp_clean += 1
            if hp_c or al_c:
                clean.add(key)

    A_ = cansplit; C = clean
    inter = len(A_ & C); uni = len(A_ | C)
    jaccard = round(inter / uni, 3) if uni else None

    # 切割總帳分帶: insuff (n<6 不能分群) / cant (n>=6 但 best_k None or <2) / cansplit 內分帶
    insuff = cant = 0
    bands = Counter()
    for k, r in recs.items():
        a = r["all"]
        if a["n"] < MIN_SZ * 2:
            insuff += 1; continue
        if a["best_k"] is None or a["best_k"] < 2:
            cant += 1; continue
        # cansplit -> 分帶
        V, sig, nlab, nclu = cramers_v_and_sig(a["contingency"]) if a["contingency"] else (None, False, 0, 0)
        bands[band_of(V, sig)] += 1

    cant_total = insuff + cant

    result = {
        "sample": sample,
        "N": N,
        "matched_significance": matched,
        # 主表 (A/C Venn)
        "A_cansplit": len(A_),
        "C_clean_permanova": len(C),
        "AnC": inter,
        "A_minus_C": len(A_ - C),
        "C_minus_A": len(C - A_),
        "jaccard": jaccard,
        "A_pct": round(100 * len(A_) / N, 1) if N else None,
        "C_pct": round(100 * len(C) / N, 1) if N else None,
        "AnC_pct": round(100 * inter / N, 1) if N else None,
        "A_minus_C_pct": round(100 * len(A_ - C) / N, 1) if N else None,
        "C_minus_A_pct": round(100 * len(C - A_) / N, 1) if N else None,
        # C 拆解
        "C_hp_clean": hp_clean,
        "C_al_clean": al_clean,
        # 切割總帳分帶
        "insuff": insuff,
        "cant": cant,
        "cant_total": cant_total,
        "cansplit": len(A_),
        "clear": bands["clear"],
        "mid": bands["mid"],
        "weak": bands["weak"],
        "vweak": bands["vweak"],
        "degen": bands["degen"],
        "elapsed_s": int(time.time() - t0),
    }
    json.dump(result, open(f"{result_dir}/result_{sample}.json", "w"), indent=1)
    print(json.dumps(result, ensure_ascii=False, indent=1), flush=True)

if __name__ == "__main__":
    main()
