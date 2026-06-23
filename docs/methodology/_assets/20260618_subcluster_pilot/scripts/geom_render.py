#!/usr/bin/env python3
"""[geometry-vs-null 觀察渲染] 每位點畫: UPGMA 樹 + 並列雙色欄(C++ null95 coarse vs 幾何 0.7×max) + 甲基熱圖 + 標準側欄 + 距離。
divergence(coarse 單色但 geom 多色)= 肉眼一眼可辨。Pool 平行。讀 geom_records.json + 重跑保留的矩陣/tsv。
⚠ geometry 欄=純幾何無 null, 本就過切; 是「閘保守度」診斷非 subclone 宣稱。"""
import os, csv, json, sys
import numpy as np
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0, os.path.dirname(__file__)); import cluster_redesign as CR
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, fcluster
from collections import Counter
from multiprocessing import Pool

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
FIGS = f"{A}/figs_geomdiv"; os.makedirs(FIGS, exist_ok=True)
ODS = {"TP": f"{WT}/output/_geomdiv_TP", "FP": f"{WT}/output/_geomdiv_FP"}
# canonical 配色一律取自 ism_heatmap_std(H),禁止本地硬編以防 drift:
#   HP1 族藍(1淺#60a5fa/1-1深#1e3a8a)、HP2 族紫(2淺#c084fc/2-1深#6b21a8); ALT 紅/REF 琥珀; T/N 橘/綠; Strand 板岩;
#   群 id(C++ coarse + 幾何)用 CLUSTER_COL 專用 teal/pink 盤(不撞 HP/ALT/T/N 語意軸)
PAL = H.CLUSTER_COL
HPC = H.HP_COL; ALC = H.ALT_COL; TNC = H.TN_COL; STR = H.STRAND_COL
NA = H.NA_GREY
MINSZ = 3


def find_region(setlab, chrom, rec_pos):
    od = ODS[setlab]
    # region dir = {chrom}_{start}_{end}, start = rec_pos
    import glob
    hits = glob.glob(f"{od}/**/{chrom}_{rec_pos}_*/clustering/phylo_groups.tsv", recursive=True)
    return os.path.dirname(os.path.dirname(hits[0])) if hits else None


def load_reads_meta(region):
    """reads/reads.tsv → read_id(str) -> (is_tumor, strand). join on read_id(同 phylo_groups.tsv 索引)。"""
    p = f"{region}/reads/reads.tsv"
    if not os.path.exists(p):
        return {}
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    if not rows:
        return {}
    cols = rows[0].keys()
    tcol = next((c for c in cols if c.lower() in ("is_tumor", "tumor", "sample_is_tumor")), None)
    scol = next((c for c in cols if c.lower() in ("strand", "read_strand")), None)
    out = {}
    for r in rows:
        out[r.get("read_id", "")] = (r.get(tcol, "") if tcol else "", r.get(scol, "") if scol else "")
    return out


def render_one(rec):
    try:
        setlab = rec["set"]; chrom = rec["chrom"]; rp = rec["rec_pos"]
        region = find_region(setlab, chrom, rp)
        if not region:
            return None
        tsv = f"{region}/clustering/phylo_groups.tsv"
        rows = list(csv.DictReader(open(tsv), delimiter="\t"))
        if len(rows) < 2 * MINSZ:
            return None
        rid = [int(r["read_id"]) for r in rows]; rname = [r["read_name"] for r in rows]
        clab = [r["coarse_label"] for r in rows]; hp = [r["hp"] for r in rows]; al = [r["alt_support"] for r in rows]
        dids, D = CR.loadm(f"{region}/distance/BERNOULLI/matrix.csv")
        idx = [dids.index(str(x)) for x in rid if str(x) in dids]
        sub = D[np.ix_(idx, idx)]
        Z, _ = CR.linkZ(sub)
        # geometry 0.7×max flat cut → per-read 群 id(size>=MINSZ 著色, 否則 gray)
        gcut = 0.7 * float(Z[:, 2].max())
        glab = fcluster(Z, t=gcut, criterion="distance")
        gsz = Counter(glab); gbig = [k for k, c in gsz.items() if c >= MINSZ]
        gcolor = {k: PAL[i % len(PAL)] for i, k in enumerate(sorted(gbig))}
        gcol = [gcolor.get(g, "#dcdcdc") for g in glab]
        # C++ coarse 著色
        terms = sorted(set(l for l in clab if l not in ("other", "outlier", "-")))
        cmap = {t: PAL[i % len(PAL)] for i, t in enumerate(terms)}
        ccol = ["#555" if l == "other" else "#dcdcdc" if l in ("outlier", "-") else cmap.get(l, "#888") for l in clab]
        # 甲基矩陣
        mr = open(f"{region}/methylation/methylation.csv").read().strip().split("\n")
        cpgs = [int(c) for c in mr[0].split(",")[1:]]
        Mall = {}
        for ln in mr[1:]:
            q = ln.split(","); Mall[q[0]] = [np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]]
        meth = np.array([Mall[str(x)] for x in rid])
        # T/N + strand
        rm = load_reads_meta(region)
        tn = [rm.get(str(x), ("", ""))[0] for x in rid]; st = [rm.get(str(x), ("", ""))[1] for x in rid]
        # 樹葉序
        dn = dendrogram(Z, orientation="left", no_plot=True); order = dn["leaves"][::-1]
        mo = meth[order]; do = sub[np.ix_(order, order)].copy(); np.fill_diagonal(do, 0); do[do < 0] = np.nan
        cc = [ccol[i] for i in order]; gg = [gcol[i] for i in order]
        hpo = [HPC.get(hp[i], NA) for i in order]; alo = [ALC.get(al[i], NA) for i in order]
        tno = [TNC.get(str(tn[i]), NA) for i in order]; sto = [STR.get(str(st[i]), NA) for i in order]
        sb = [("C++ null95", cc), ("幾何0.7cut", gg), ("HP", hpo), ("ALT", alo), ("T/N", tno), ("Strand", sto)]
        # 樹枝著色: 依 C++ coarse(子樹單一標籤=該色)
        n = len(rid); descn = {i: [i] for i in range(n)}
        for i in range(len(Z)):
            aa, bb = int(Z[i, 0]), int(Z[i, 1]); descn[n + i] = descn[aa] + descn[bb]
        nodelab = {nd: (list({clab[i] for i in descn[nd]})[0] if len({clab[i] for i in descn[nd]}) == 1 else None) for nd in range(2 * n - 1)}
        def lc(l): return "#555" if l == "other" else "#dcdcdc" if l in ("outlier", "-") else cmap.get(l, "#888")
        mc, dc = H.mpl_cmaps()
        nsb = len(sb); wr = [1.0] + [0.05] * nsb + [1.5, 0.12] + [1.2]
        fig = plt.figure(figsize=(13.5, 5.4)); gs = fig.add_gridspec(1, len(wr), width_ratios=wr, wspace=0.04)
        c = 0; ax = fig.add_subplot(gs[0, c]); c += 1
        dendrogram(Z, orientation="left", no_labels=True, ax=ax,
                   link_color_func=lambda k: lc(nodelab.get(k)) if nodelab.get(k) else "#888")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title("UPGMA(枝=C++ null95)", fontsize=8); [ax.spines[s].set_visible(False) for s in ax.spines]
        for lb, hx in sb: H._sb(fig.add_subplot(gs[0, c]), hx, lb); c += 1
        am = fig.add_subplot(gs[0, c]); c += 1
        am.imshow(mo, aspect="auto", cmap=mc, vmin=0, vmax=1, interpolation="nearest")
        try:
            sx = H.snv_fractional_x(cpgs, rp + 5000)
            if sx is not None: am.axvline(sx, color=H.SNV_COL, lw=2)
        except Exception: pass
        am.set_xticks([]); am.set_yticks([]); am.set_title("甲基 read×CpG", fontsize=8.5)
        fig.add_subplot(gs[0, c]).axis("off"); c += 1
        ad = fig.add_subplot(gs[0, c]); vmax = max(0.5, float(np.nanmax(do)) if np.isfinite(np.nanmax(do)) else 0.5)
        ad.imshow(do, aspect="auto", cmap=dc, vmin=0, vmax=vmax, interpolation="nearest"); ad.set_xticks([]); ad.set_yticks([]); ad.set_title("距離", fontsize=8.5)
        dv = "★DIVERGENCE " if rec["divergence"] else ""
        fig.suptitle(f"{dv}{chrom}:{rp+5000} [{setlab}/{rec['wstate']}] n={rec['n']}  null95={rec['coarse_ng']} null90={rec['fine_ng']} 幾何={rec['geometry_ng']}  V_hp={rec['V_hp']} V_al={rec['V_allele']}",
                     fontsize=9.5, y=1.02)
        fn = f"{FIGS}/geom_{setlab}_{chrom}_{rp}.png"; fig.savefig(fn, dpi=125, bbox_inches="tight"); plt.close(fig)
        return (f"{chrom}:{rp}", rec["divergence"], f"figs_geomdiv/geom_{setlab}_{chrom}_{rp}.png")
    except Exception as e:
        return ("ERR " + rec.get("chrom", "?") + ":" + str(rec.get("rec_pos", "?")) + " " + str(e)[:70], False, None)


if __name__ == "__main__":
    recs = json.load(open(f"{A}/geom_records.json"))
    print(f"rendering {len(recs)} loci (Pool)...", flush=True)
    with Pool(24) as p:
        out = p.map(render_one, recs)
    ok = [r for r in out if r and r[2]]; err = [r for r in out if r and not r[2]]
    json.dump({f"{r[0]}": r[2] for r in ok}, open(f"{A}/geom_render_index.json", "w"), ensure_ascii=False, indent=1)
    print(f"DONE rendered={len(ok)} divergence_figs={sum(1 for r in ok if r[1])} errors={len(err)}", flush=True)
    for e in err[:10]: print("  ", e[0])
