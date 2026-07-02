#!/usr/bin/env python3
"""全基因組 phylo-v3.1 — TP+FP, big7 本機, 平行。對照 4-gate 31.3% 的『有結構』比例。
每 chr: binary 跑 TP+FP VCF → 各 locus 矩陣 → v3.1 標籤 + 對齊分類 → 收 {ngroups, maxdepth, aligned} → 刪暫存。
v3.1=descend quarantine(物有所值)+per-subgroup null(RNULL=25 提速)。對齊=top split CramerV(hp|allele)≥0.3。
輸出 phylo_v3_wg_records.json + summary。用法: python3 phylo_v3_wg.py [CHRS]  e.g. 22(smoke) / 留空=1-22。"""
import os, csv, glob, json, sys, subprocess, shutil, time
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
import numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; BIN = f"{WT}/build/bin/inter_sub_mod"
TUMOR = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"
NORMAL = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF = "/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"
VCFDIR = "/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
os.environ["TMPDIR"] = "/big7_disk/liaoyoyo2001/tmp"; os.makedirs("/big7_disk/liaoyoyo2001/tmp", exist_ok=True)
NPROC = 24; MINSZ = 3; SEP_MIN = 1.3; RNULL = 25
_raw = sys.argv[1].split(",") if len(sys.argv) > 1 else [str(c) for c in range(1, 23)]
CHRS = [c if str(c).startswith("chr") else f"chr{c}" for c in _raw]


def _bw(s, S1, S2):
    bet = s[np.ix_(S1, S2)]; bet = bet[bet >= 0]
    w1 = s[np.ix_(S1, S1)][np.triu_indices(len(S1), 1)]; w1 = w1[w1 >= 0]
    w2 = s[np.ix_(S2, S2)][np.triu_indices(len(S2), 1)]; w2 = w2[w2 >= 0]
    wm = np.concatenate([w1, w2]) if (w1.size or w2.size) else np.array([])
    if bet.size == 0 or wm.size == 0 or wm.mean() <= 1e-6: return None
    return float(bet.mean()) / float(wm.mean())


def _tree(D):
    n = D.shape[0]; Z, s = CR.linkZ(D)
    desc = {i: [i] for i in range(n)}; ch = {}
    for i in range(len(Z)):
        a, b = int(Z[i, 0]), int(Z[i, 1]); desc[n + i] = desc[a] + desc[b]; ch[n + i] = (a, b)
    return Z, s, desc, ch, n


def cramerV(a, b):
    ca = sorted(set(a)); cb = sorted(set(b)); ia = {c: i for i, c in enumerate(ca)}; ib = {c: i for i, c in enumerate(cb)}
    tab = np.zeros((len(ca), len(cb)))
    for x, y in zip(a, b): tab[ia[x], ib[y]] += 1
    nn = tab.sum()
    if nn == 0 or min(len(ca), len(cb)) < 2: return 0.0
    row = tab.sum(1, keepdims=True); col = tab.sum(0, keepdims=True); exp = row * col / nn
    chi2 = np.nansum((tab - exp) ** 2 / np.where(exp > 0, exp, 1))
    return float(np.sqrt(chi2 / (nn * (min(len(ca), len(cb)) - 1))))


def v3_label(sub, P, rng):
    Z, s, desc, ch, n = _tree(sub)

    def subnull(leaves):
        S = np.array(leaves); m = len(S); Psub = P[S]; ns = []
        for _ in range(RNULL):
            Pn = Psub.copy()
            for cc in range(Pn.shape[1]):
                col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
                if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
            Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); Dn = np.maximum(Dn, Dn.T)
            _, sn, dn, cn, _ = _tree(Dn); nc1, nc2 = cn[2 * m - 2]
            ns.append(_bw(sn, np.array(dn[nc1]), np.array(dn[nc2])))
        ns = [x for x in ns if x is not None]
        return float(np.percentile(ns, 95)) if ns else 0

    def descend(node):
        cur = node; quar = []
        while cur in ch:
            c1, c2 = ch[cur]; s1, s2 = len(desc[c1]), len(desc[c2])
            if min(s1, s2) >= MINSZ: return cur, quar
            small, big = (c1, c2) if s1 < s2 else (c2, c1); quar += desc[small]; cur = big
        return None, quar

    def test(bnode):
        c1, c2 = ch[bnode]; r = _bw(s, np.array(desc[c1]), np.array(desc[c2]))
        if r is None or r < SEP_MIN: return False
        return r > subnull(desc[bnode])

    lab = [None] * n

    def rec(node, label):
        leaves = desc[node]
        if len(leaves) < 2 * MINSZ:
            for i in leaves: lab[i] = label; return
        bnode, quar = descend(node)
        if bnode is not None and test(bnode):
            for i in quar: lab[i] = "outlier"
            c1, c2 = ch[bnode]; big, small = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
            rec(big, label + "-1"); rec(small, label + "-2")
        else:
            for i in leaves: lab[i] = label
    rec(2 * n - 2, "1")
    lab = [l if l else "outlier" for l in lab]
    sm = {L for L, c in Counter(l for l in lab if l != "outlier").items() if c < MINSZ}
    return [("outlier" if l in sm else l) for l in lab]


def proc_locus(arg):
    rd, setlab = arg
    try:
        reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
        dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
        rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); nC = len(rows[0].split(",")) - 1
        mi = {}; M = []
        for j, ln in enumerate(rows[1:]):
            q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
        M = np.array(M); it = lambda t: str(t) in ("1", "true", "True")
        ids = [x for x in dids if x in reads and it(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
        if len(ids) < MINSZ * 2: return None
        sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
        if len(kp) < MINSZ * 2: return None
        ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
        lab = v3_label(sub, P, np.random.default_rng(20260622))
        terms = [l for l in lab if l != "outlier"]; ng = len(set(terms)); nout = sum(1 for l in lab if l == "outlier")
        maxdepth = max((l.count("-") for l in set(terms)), default=0)
        aligned = False; Vhp = Val = 0.0
        if ng >= 2:
            # 用頂層 2 群(標籤第一字元 1 vs 其餘? 用完整標籤對齊)
            hp = [reads[ids[i]]["hp"] for i in range(n) if lab[i] != "outlier"]
            al = [reads[ids[i]]["alt_support"] for i in range(n) if lab[i] != "outlier"]
            tl = [lab[i] for i in range(n) if lab[i] != "outlier"]
            Vhp = cramerV(tl, hp); Val = cramerV(tl, al); aligned = max(Vhp, Val) >= 0.3
        return {"pos": rd.split("/")[-1].split("_")[1] if "_" in rd.split("/")[-1] else None, "set": setlab, "n": n, "C": nC,
                "ngroups": ng, "maxdepth": maxdepth, "n_outlier": nout, "aligned": aligned,
                "V_hp": round(Vhp, 2), "V_allele": round(Val, 2)}
    except Exception as e:
        return {"err": str(e)[:80], "set": setlab}


def run():
    out = []; t0 = time.time()
    for chrom in CHRS:
        for setlab, pref in (("TP", "filtered_snv_tp"), ("FP", "filtered_snv_fp")):
            vcf = f"{VCFDIR}/{pref}_{chrom}.vcf.gz"
            if not os.path.exists(vcf): continue
            od = f"{WT}/output/_pv3wg_{chrom}_{setlab}"; shutil.rmtree(od, ignore_errors=True); os.makedirs(od, exist_ok=True)
            subprocess.run([BIN, "-t", TUMOR, "-n", NORMAL, "-r", REF, "-v", vcf, "-w", "5000", "-j", "16",
                            "--distance-metric", "BERNOULLI", "--nan-distance-strategy", "SKIP", "-o", od],
                           stdout=open(f"{od}/run.log", "w"), stderr=subprocess.STDOUT)
            mats = glob.glob(f"{od}/**/distance/BERNOULLI/matrix.csv", recursive=True)
            args = [(m.rsplit("/distance/", 1)[0], setlab) for m in mats]
            with Pool(NPROC) as pool:
                res = [r for r in pool.map(proc_locus, args) if r]
            for r in res:
                if r and "pos" in r: r["chrom"] = chrom
            out += [r for r in res if r and "err" not in r]
            shutil.rmtree(od, ignore_errors=True)
            print(f"[{chrom}/{setlab}] loci={len(mats)} kept={len([r for r in res if 'err' not in r])} elapsed={int(time.time()-t0)}s", flush=True)
            json.dump(out, open(f"{A}/phylo_v3_wg_records.json", "w"))
    TP = [r for r in out if r["set"] == "TP"]; FP = [r for r in out if r["set"] == "FP"]
    def summ(g):
        if not g: return {}
        multi = [r for r in g if r["ngroups"] >= 2]
        return {"n": len(g), "structure_pct": round(100 * len(multi) / len(g), 2),
                "aligned_pct": round(100 * sum(1 for r in multi if r["aligned"]) / len(g), 2),
                "unaligned_pct": round(100 * sum(1 for r in multi if not r["aligned"]) / len(g), 2),
                "no_structure_pct": round(100 * sum(1 for r in g if r["ngroups"] < 2) / len(g), 2),
                "ngroups_dist": dict(Counter(r["ngroups"] for r in g)),
                "maxdepth_dist": dict(Counter(r["maxdepth"] for r in multi))}
    S = {"TP": summ(TP), "FP": summ(FP), "elapsed_s": int(time.time() - t0), "params": {"RNULL": RNULL, "SEP_MIN": SEP_MIN, "chrs": CHRS}}
    json.dump(S, open(f"{A}/phylo_v3_wg_summary.json", "w"), indent=2)
    print("DONE", json.dumps(S, ensure_ascii=False))


if __name__ == "__main__":
    run()
