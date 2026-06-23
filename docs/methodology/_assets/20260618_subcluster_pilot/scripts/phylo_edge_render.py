#!/usr/bin/env python3
"""[邊界渲染] 抽 genome-wide 邊界代表位點 → mini-VCF → binary → v3 render，供肉眼確認邊界合理。
渲染: E1 0群退化 / E2 高群數 / E4 FP unaligned / E6 大n。每張=樹(v3群色)+甲基+HP/ALT+距離。
重用 phylo_v3_render_all 的上色+render；C++ binary 在 worktree build。"""
import os, csv, glob, json, sys, subprocess, shutil, gzip
import numpy as np
from collections import Counter
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; BIN = f"{WT}/build/bin/inter_sub_mod"
TUMOR = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"; NORMAL = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF = "/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"; VCFDIR = "/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
os.environ["TMPDIR"] = "/big7_disk/liaoyoyo2001/tmp"
FIGS = f"{A}/figs_edge"; os.makedirs(FIGS, exist_ok=True)
MINSZ = 3; SEP_MIN = 1.3; RNULL = 40
PALETTE = H.CLUSTER_COL

# 選渲染的邊界位點: (chrom, pos, set, 類別標)
samples = json.load(open(f"{A}/phylo_edge_samples.json"))
SEL = []
for cat, key in [("E1_degenerate", "E1 0群退化"), ("E2_highk", "E2 高群數"), ("E4_fp_unaligned", "E4 FP-unaligned"), ("E6_bign", "E6 大n")]:
    for s in samples[cat][:2]:
        SEL.append((s[0], str(s[1]), s[2], key))


def _bw(s, S1, S2):
    bet = s[np.ix_(S1, S2)]; bet = bet[bet >= 0]
    w1 = s[np.ix_(S1, S1)][np.triu_indices(len(S1), 1)]; w1 = w1[w1 >= 0]
    w2 = s[np.ix_(S2, S2)][np.triu_indices(len(S2), 1)]; w2 = w2[w2 >= 0]
    wm = np.concatenate([w1, w2]) if (w1.size or w2.size) else np.array([])
    if bet.size == 0 or wm.size == 0 or wm.mean() <= 1e-6: return None
    return float(bet.mean()) / float(wm.mean())


def _tree(D):
    n = D.shape[0]; Z, s = CR.linkZ(D)
    dc = {i: [i] for i in range(n)}; ch = {}
    for i in range(len(Z)):
        a, b = int(Z[i, 0]), int(Z[i, 1]); dc[n + i] = dc[a] + dc[b]; ch[n + i] = (a, b)
    return Z, s, dc, ch, n


def v3_label(sub, P, rng):
    Z, s, desc, ch, n = _tree(sub)
    def subnull(leaves):
        S = np.array(leaves); m = len(S); Ps = P[S]; ns = []
        for _ in range(RNULL):
            Pn = Ps.copy()
            for cc in range(Pn.shape[1]):
                col = Pn[:, cc]; vi = np.where(~np.isnan(col))[0]
                if vi.size > 1: Pn[vi, cc] = col[rng.permutation(vi)]
            Dn = CR.bernoulli_dist(Pn); np.fill_diagonal(Dn, 0); Dn = np.maximum(Dn, Dn.T)
            _, sn, dn, cn, _ = _tree(Dn); nc1, nc2 = cn[2 * m - 2]
            ns.append(_bw(sn, np.array(dn[nc1]), np.array(dn[nc2])))
        ns = [x for x in ns if x is not None]; return float(np.percentile(ns, 95)) if ns else 0
    def descend(node):
        cur = node; q = []
        while cur in ch:
            c1, c2 = ch[cur]; s1, s2 = len(desc[c1]), len(desc[c2])
            if min(s1, s2) >= MINSZ: return cur, q
            sm, bg = (c1, c2) if s1 < s2 else (c2, c1); q += desc[sm]; cur = bg
        return None, q
    def test(bn):
        c1, c2 = ch[bn]; r = _bw(s, np.array(desc[c1]), np.array(desc[c2]))
        if r is None or r < SEP_MIN: return False
        return r > subnull(desc[bn])
    lab = [None] * n
    def rec(node, label):
        lv = desc[node]
        if len(lv) < 2 * MINSZ:
            for i in lv: lab[i] = label; return
        bn, q = descend(node)
        if bn is not None and test(bn):
            for i in q: lab[i] = "outlier"
            c1, c2 = ch[bn]; big, sm = (c1, c2) if len(desc[c1]) >= len(desc[c2]) else (c2, c1)
            rec(big, label + "-1"); rec(sm, label + "-2")
        else:
            for i in lv: lab[i] = label
    rec(2 * n - 2, "1")
    lab = [l if l else "outlier" for l in lab]
    sms = {L for L, c in Counter(l for l in lab if l != "outlier").items() if c < MINSZ}
    return [("outlier" if l in sms else l) for l in lab], Z, desc


# --- 建 mini-VCF (per set) ---
def vcf_lines(chrom, pos, setlab):
    pref = "filtered_snv_tp" if setlab == "TP" else "filtered_snv_fp"
    fp = f"{VCFDIR}/{pref}_{chrom}.vcf.gz"
    target = int(pos) + 5000  # region 目錄 pos = SNV − window(5000) → 還原 SNV pos
    hdr = []; rec = None
    with gzip.open(fp, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"): hdr.append(ln)
            else:
                f = ln.split("\t")
                if f[0] == chrom and f[1] == str(target): rec = ln
    return hdr, rec


od = f"{WT}/output/_edge_render"; shutil.rmtree(od, ignore_errors=True); os.makedirs(od, exist_ok=True)
by_set = {}
for chrom, pos, setlab, cat in SEL:
    by_set.setdefault(setlab, []).append((chrom, pos, cat))
locus_meta = {}
for setlab, items in by_set.items():
    hdr = None; recs = []
    for chrom, pos, cat in items:
        h, r = vcf_lines(chrom, pos, setlab)
        if r: hdr = hdr or h; recs.append(r); locus_meta[f"{chrom}_{pos}"] = cat
    if not recs: continue
    def _ck(ln):
        f = ln.split("\t"); cn = f[0].replace("chr", "").replace("X", "23").replace("Y", "24")
        return (int(cn) if cn.isdigit() else 99, int(f[1]))
    recs.sort(key=_ck)  # tabix 需 chrom+pos 排序
    mvcf = f"{od}/mini_{setlab}.vcf"
    with open(mvcf, "w") as o:
        o.writelines(hdr); o.writelines(recs)
    subprocess.run(["bgzip", "-f", mvcf]); subprocess.run(["tabix", "-p", "vcf", f"{mvcf}.gz"])
    oo = f"{od}/{setlab}"; os.makedirs(oo, exist_ok=True)
    subprocess.run([BIN, "-t", TUMOR, "-n", NORMAL, "-r", REF, "-v", f"{mvcf}.gz", "-w", "5000", "-j", "16",
                    "--distance-metric", "BERNOULLI", "--nan-distance-strategy", "SKIP", "-o", oo],
                   stdout=open(f"{oo}/run.log", "w"), stderr=subprocess.STDOUT)

# --- render ---
HPC = H.HP_COL
alc = H.ALT_COL
out = []
for mp in glob.glob(f"{od}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    base = rd.split("/")[-1]; parts = base.split("_")  # chrN_start_end (窗座標)
    if not base.startswith("chr") or len(parts) < 2: continue
    key = f"{parts[0]}_{parts[1]}"  # = records key (窗起點 = SNV-5000)
    snv_pos = (int(parts[1]) + int(parts[2])) // 2 if len(parts) >= 3 else int(parts[1]) + 5000
    if key not in locus_meta: continue
    cat = locus_meta[key]
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(mp); di = {x: i for i, x in enumerate(dids)}
    rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs = [int(c) for c in rows[0].split(",")[1:]]
    mi = {}; M = []
    for j, ln in enumerate(rows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    if len(ids) < MINSZ * 2: continue
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]; n = len(ids); P = np.array([M[mi[x]] for x in ids])
    lab, Z, desc = v3_label(sub, P, np.random.default_rng(20260622))
    terms = sorted(set(l for l in lab if l != "outlier")); cmap = {t: PALETTE[i % len(PALETTE)] for i, t in enumerate(terms)}
    lc = lambda l: "#cfcfcf" if l in (None, "outlier", "-") else cmap.get(l, "#888")
    dn = dendrogram(Z, orientation="left", no_plot=True); order = dn["leaves"][::-1]
    meth = np.array([P[i] for i in order]); dist = sub[np.ix_(order, order)].copy(); np.fill_diagonal(dist, 0); dist[dist < 0] = np.nan
    lab_o = [lab[i] for i in order]; hp_o = [reads[ids[i]]["hp"] for i in order]; al_o = [reads[ids[i]]["alt_support"] for i in order]
    nodelab = {nd: (list({lab[i] for i in desc[nd]})[0] if len({lab[i] for i in desc[nd]}) == 1 else None) for nd in range(2 * n - 1)}
    sb = [("v3群", [lc(l) for l in lab_o]), ("HP", [HPC.get(h, "#ddd") for h in hp_o]), ("ALT", [alc.get(a, "#ddd") for a in al_o])]
    mc, dc = H.mpl_cmaps(); nsb = len(sb); wr = [1.1] + [0.045] * nsb + [1.7, 0.14] + [0.045] * nsb + [1.4]
    fig = plt.figure(figsize=(13.5, 5.6)); gs = fig.add_gridspec(1, len(wr), width_ratios=wr, wspace=0.04)
    c = 0; axdn = fig.add_subplot(gs[0, c]); c += 1
    dendrogram(Z, orientation="left", link_color_func=lambda k: lc(nodelab.get(k)) if nodelab.get(k) else "#ccc", ax=axdn, no_labels=True)
    axdn.set_xticks([]); axdn.set_yticks([]); axdn.set_title("UPGMA(v3群)", fontsize=8); [axdn.spines[x].set_visible(False) for x in axdn.spines]
    for lb, hx in sb: H._sb(fig.add_subplot(gs[0, c]), hx, lb); c += 1
    axm = fig.add_subplot(gs[0, c]); c += 1; axm.imshow(meth, aspect="auto", cmap=mc, vmin=0, vmax=1, interpolation="nearest")
    sx = H.snv_fractional_x(cpgs, snv_pos)
    if sx is not None: axm.axvline(sx, color=H.SNV_COL, lw=2)
    axm.set_xticks([]); axm.set_yticks([]); axm.set_title("甲基 read×CpG", fontsize=8.5)
    fig.add_subplot(gs[0, c]).axis("off"); c += 1
    for lb, hx in sb: H._sb(fig.add_subplot(gs[0, c]), hx, lb); c += 1
    axd = fig.add_subplot(gs[0, c]); vmax = max(0.5, float(np.nanmax(dist)) if np.isfinite(np.nanmax(dist)) else 0.5)
    axd.imshow(dist, aspect="auto", cmap=dc, vmin=0, vmax=vmax, interpolation="nearest"); axd.set_xticks([]); axd.set_yticks([]); axd.set_title("距離", fontsize=8.5)
    ng = len(terms); nout = sum(1 for l in lab if l == "outlier")
    fig.suptitle(f"[{cat}] {key}  n={n} CpG={len(cpgs)} → v3 {ng}群 離群{nout}  標籤{terms[:6]}", fontsize=9.5, y=1.02)
    fn = f"{FIGS}/edge_{key}.png"; fig.savefig(fn, dpi=155, bbox_inches="tight"); plt.close(fig)
    out.append({"key": key, "cat": cat, "n": n, "C": len(cpgs), "ngroups": ng, "n_outlier": nout, "png": f"figs_edge/edge_{key}.png"})
    print(f"  [{cat}] {key} n={n} C={len(cpgs)} → v3 {ng}群 離群{nout}", flush=True)
json.dump(out, open(f"{A}/phylo_edge_render.json", "w"), ensure_ascii=False, indent=1)
shutil.rmtree(od, ignore_errors=True)
print(f"\nDONE {len(out)} 邊界圖 → figs_edge/")
