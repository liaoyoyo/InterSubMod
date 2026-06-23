#!/usr/bin/env python3
"""[v4 重驗] 用 v4(coarse/fine + other + instability) 重渲染用戶回饋的邊界位點。
genome-wide(chr8/chr7/chr4): mini-VCF→binary; pilot(chr20/chr22): 快取。
圖: 樹+甲基+HP/ALT+距離; coarse 群上色 + 'other' 深灰 + caption 報 coarse/fine/other/unstable。"""
import os, csv, glob, json, sys, subprocess, shutil, gzip
import numpy as np
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/scripts"); import ism_heatmap_std as H
sys.path.insert(0, os.path.dirname(__file__))
import cluster_redesign as CR
from phylo_v4 import v4_label, analyze_v4, ng_of  # guarded import (no main run)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

WT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
A = f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"; BIN = f"{WT}/build/bin/inter_sub_mod"
TUMOR = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam"; NORMAL = "/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam"
REF = "/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa"; VCFDIR = "/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup"
CACHE = f"{WT}/output/_ws_render"; os.environ["TMPDIR"] = "/big7_disk/liaoyoyo2001/tmp"
FIGS = f"{A}/figs_v4edge"; os.makedirs(FIGS, exist_ok=True)
MINSZ = 3
PALETTE = H.CLUSTER_COL
HPC = H.HP_COL
alc = H.ALT_COL

# genome-wide(set, chrom, records_pos, 標) ; pilot(chrom_pos 快取)
GW = [("FP", "chr8", "144597485", "E4 other(殘離群)"), ("FP", "chr7", "155610406", "E4 fine候選"),
      ("TP", "chr4", "190105905", "E6 fine候選"), ("TP", "chr4", "190103740", "E6 fine候選")]
PIL = [("chr20_42981498", "F1 instability"), ("chr22_30454004", "other已生效")]


def load_mat(rd):
    reads = {x["read_id"]: x for x in csv.DictReader(open(f"{rd}/reads/reads.tsv"), delimiter="\t")}
    dids, D = CR.loadm(f"{rd}/distance/BERNOULLI/matrix.csv"); di = {x: i for i, x in enumerate(dids)}
    rows = open(f"{rd}/methylation/methylation.csv").read().strip().split("\n"); cpgs = [int(c) for c in rows[0].split(",")[1:]]
    mi = {}; M = []
    for j, ln in enumerate(rows[1:]):
        q = ln.split(","); mi[q[0]] = j; M.append([np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]])
    M = np.array(M); itf = lambda t: str(t) in ("1", "true", "True")
    ids = [x for x in dids if x in reads and itf(reads[x]["is_tumor"]) and reads[x]["hp"] in CR.LABMAP and x in mi]
    if len(ids) < MINSZ * 2: return None
    sub = D[np.ix_([di[x] for x in ids], [di[x] for x in ids])]; kp = CR.peel(sub)
    ids = [ids[i] for i in kp]; sub = sub[np.ix_(kp, kp)]
    return ids, sub, np.array([M[mi[x]] for x in ids]), reads, cpgs


def render(key, snv_pos, sub, P, ids, reads, cpgs, tag):
    a = analyze_v4(sub, P); lab = a["coarse_lab"]; n = len(ids)
    terms = sorted(set(l for l in lab if l not in ("outlier", "other"))); cmap = {t: PALETTE[i % len(PALETTE)] for i, t in enumerate(terms)}
    def lc(l): return "#555555" if l == "other" else "#dcdcdc" if l in (None, "outlier", "-") else cmap.get(l, "#888")
    Z, _, desc, _, _ = (lambda D: (CR.linkZ(D)[0], None, None, None, None))(sub)  # placeholder
    Z, _ = CR.linkZ(sub)
    desc = {i: [i] for i in range(n)}
    for i in range(len(Z)):
        aa, bb = int(Z[i, 0]), int(Z[i, 1]); desc[n + i] = desc[aa] + desc[bb]
    dn = dendrogram(Z, orientation="left", no_plot=True); order = dn["leaves"][::-1]
    meth = np.array([P[i] for i in order]); dist = sub[np.ix_(order, order)].copy(); np.fill_diagonal(dist, 0); dist[dist < 0] = np.nan
    lab_o = [lab[i] for i in order]; hp_o = [reads[ids[i]]["hp"] for i in order]; al_o = [reads[ids[i]]["alt_support"] for i in order]
    nodelab = {nd: (list({lab[i] for i in desc[nd]})[0] if len({lab[i] for i in desc[nd]}) == 1 else None) for nd in range(2 * n - 1)}
    sb = [("coarse群", [lc(l) for l in lab_o]), ("HP", [HPC.get(h, "#ddd") for h in hp_o]), ("ALT", [alc.get(x, "#ddd") for x in al_o])]
    mc, dc = H.mpl_cmaps(); nsb = len(sb); wr = [1.1] + [0.045] * nsb + [1.7, 0.14] + [0.045] * nsb + [1.4]
    fig = plt.figure(figsize=(13.5, 5.6)); gs = fig.add_gridspec(1, len(wr), width_ratios=wr, wspace=0.04)
    c = 0; axdn = fig.add_subplot(gs[0, c]); c += 1
    dendrogram(Z, orientation="left", link_color_func=lambda k: lc(nodelab.get(k)) if nodelab.get(k) else "#ccc", ax=axdn, no_labels=True)
    axdn.set_xticks([]); axdn.set_yticks([]); axdn.set_title("UPGMA(coarse群)", fontsize=8); [axdn.spines[x].set_visible(False) for x in axdn.spines]
    for lb, hx in sb: H._sb(fig.add_subplot(gs[0, c]), hx, lb); c += 1
    axm = fig.add_subplot(gs[0, c]); c += 1; axm.imshow(meth, aspect="auto", cmap=mc, vmin=0, vmax=1, interpolation="nearest")
    sx = H.snv_fractional_x(cpgs, snv_pos)
    if sx is not None: axm.axvline(sx, color=H.SNV_COL, lw=2)
    axm.set_xticks([]); axm.set_yticks([]); axm.set_title("甲基 read×CpG (灰=other殘群)", fontsize=8.5)
    fig.add_subplot(gs[0, c]).axis("off"); c += 1
    for lb, hx in sb: H._sb(fig.add_subplot(gs[0, c]), hx, lb); c += 1
    axd = fig.add_subplot(gs[0, c]); vmax = max(0.5, float(np.nanmax(dist)) if np.isfinite(np.nanmax(dist)) else 0.5)
    axd.imshow(dist, aspect="auto", cmap=dc, vmin=0, vmax=vmax, interpolation="nearest"); axd.set_xticks([]); axd.set_yticks([]); axd.set_title("距離", fontsize=8.5)
    uns = f" ⚠unstable{a['ng_range']}" if a["unstable"] else ""
    fine = f" · fine候選 {a['fine_ng']}群" if a["fine_ng"] > a["coarse_ng"] else ""
    oth = f" · other殘群 {a['n_other']}" if a["n_other"] > 0 else ""
    fig.suptitle(f"[{tag}] {key}  n={n} CpG={len(cpgs)} → coarse {a['coarse_ng']}群{fine}{oth}{uns}", fontsize=10, y=1.02)
    fn = f"{FIGS}/v4_{key}.png"; fig.savefig(fn, dpi=155, bbox_inches="tight"); plt.close(fig)
    return {"key": key, "tag": tag, "n": n, "C": len(cpgs), "coarse_ng": a["coarse_ng"], "fine_ng": a["fine_ng"],
            "n_other": a["n_other"], "unstable": a["unstable"], "ng_range": list(a["ng_range"]), "png": f"figs_v4edge/v4_{key}.png"}


out = []
# pilot 快取
cdir = {}
for mp in glob.glob(f"{CACHE}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]
    for part in rd.split("/"):
        if part.count("_") == 1 and part.startswith("chr"): cdir[part] = rd
for key, tag in PIL:
    rd = cdir.get(key)
    if not rd: print("skip pilot", key); continue
    r = load_mat(rd)
    if not r: continue
    ids, sub, P, reads, cpgs = r
    out.append(render(key, int(key.split("_")[1]), sub, P, ids, reads, cpgs, tag))
    print(f"  [{tag}] {key} done", flush=True)

# genome-wide binary
od = f"{WT}/output/_v4edge"; shutil.rmtree(od, ignore_errors=True); os.makedirs(od, exist_ok=True)
bs = {}
for setlab, chrom, pos, tag in GW: bs.setdefault(setlab, []).append((chrom, pos, tag))
meta = {}
for setlab, items in bs.items():
    pref = "filtered_snv_tp" if setlab == "TP" else "filtered_snv_fp"
    recs = []; hdr = None
    for chrom, pos, tag in items:
        target = int(pos) + 5000
        with gzip.open(f"{VCFDIR}/{pref}_{chrom}.vcf.gz", "rt") as fh:
            for ln in fh:
                if ln.startswith("#"):
                    if hdr is None: pass
                    hdr = hdr or []
                    if hdr is not None and ln.startswith("#"): pass
                else:
                    f = ln.split("\t")
                    if f[0] == chrom and f[1] == str(target): recs.append(ln); meta[f"{chrom}_{pos}"] = tag
        # 收 header 一次
    # 重收 header(完整)
    with gzip.open(f"{VCFDIR}/{pref}_{items[0][0]}.vcf.gz", "rt") as fh:
        hdr = [ln for ln in fh if ln.startswith("#")]
    if not recs: continue
    recs.sort(key=lambda ln: (int(ln.split("\t")[0].replace("chr", "").replace("X", "23").replace("Y", "24")) if ln.split("\t")[0].replace("chr", "").replace("X", "23").replace("Y", "24").isdigit() else 99, int(ln.split("\t")[1])))
    mvcf = f"{od}/m_{setlab}.vcf"; open(mvcf, "w").writelines(hdr + recs)
    subprocess.run(["bgzip", "-f", mvcf]); subprocess.run(["tabix", "-p", "vcf", f"{mvcf}.gz"])
    oo = f"{od}/{setlab}"; os.makedirs(oo, exist_ok=True)
    subprocess.run([BIN, "-t", TUMOR, "-n", NORMAL, "-r", REF, "-v", f"{mvcf}.gz", "-w", "5000", "-j", "16",
                    "--distance-metric", "BERNOULLI", "--nan-distance-strategy", "SKIP", "-o", oo],
                   stdout=open(f"{oo}/run.log", "w"), stderr=subprocess.STDOUT)
for mp in glob.glob(f"{od}/**/distance/BERNOULLI/matrix.csv", recursive=True):
    rd = mp.rsplit("/distance/", 1)[0]; base = rd.split("/")[-1]; parts = base.split("_")
    if not base.startswith("chr") or len(parts) < 3: continue
    key = f"{parts[0]}_{parts[1]}"; snv = (int(parts[1]) + int(parts[2])) // 2
    if key not in meta: continue
    r = load_mat(rd)
    if not r: continue
    ids, sub, P, reads, cpgs = r
    out.append(render(key, snv, sub, P, ids, reads, cpgs, meta[key]))
    print(f"  [{meta[key]}] {key} done", flush=True)
shutil.rmtree(od, ignore_errors=True)
json.dump(out, open(f"{A}/phylo_v4_edge.json", "w"), ensure_ascii=False, indent=1)
print(f"\nDONE {len(out)} v4 邊界圖 → figs_v4edge/")
for o in out: print(f"  {o['key']}: coarse {o['coarse_ng']} fine {o['fine_ng']} other {o['n_other']} unstable {o['unstable']}{o['ng_range'] if o['unstable'] else ''}")
