#!/usr/bin/env python3
"""[多位點完整觀察 gallery] 對每個 multi-sSNV 連鎖位點: 自動收 sSNV(anchor + 窗內 TP+FP, normal 確認真 somatic)
→ per-read 基因型 + 甲基 + HP → genotype-group 分群 → BERNOULLI 距離(精確) + UPGMA + fcluster + ISM coarse_label
交叉 → 渲染 距離熱圖/UPGMA樹(按基因型上色)/甲基熱圖(按基因型分組)/sSNV barcode SVG + 表。組成 gallery HTML。
輸出 ../../20260625_multilocus_subclone_observation_01.standalone.html。"""
import json, csv, io, base64, gzip, sys
from collections import Counter, defaultdict
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster, dendrogram
from scipy.spatial.distance import squareform
sys.path.insert(0, "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot/scripts")
import p2_snv_linkage as P2

A = P2.A; VD = P2.VD
LOCI = [("chr17", "48360161", "TP"), ("chr11", "91134191", "TP"), ("chr1", "57478998", "TP"),
        ("chr1", "50031328", "TP"), ("chr22", "33135662", "FP"), ("chr19", "17533317", "TP"), ("chr17", "39668348", "TP")]
RD = {f"{c['chrom']}:{c['pos']}:{c['set']}": c["region_dir"] for c in json.load(open(f"{A}/cis_candidates_resolved.json"))}
TX = "#141413"; BD = "#E3DACC"; MUT = "#6B6862"; AC = "#D97757"
PAL = ["#4A6E8A", "#D97757", "#8A6FA0", "#5B8A5B", "#C0563F", "#C98A5B", "#6B6862", "#9A958C"]
plt.rcParams.update({"font.size": 9, "text.color": TX, "axes.edgecolor": BD, "figure.facecolor": "white", "axes.facecolor": "white"})


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=108, bbox_inches="tight"); plt.close(fig); return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def vcf_window(chrom, lo, hi):
    out = []
    for tag in ("tp", "fp"):
        try:
            with gzip.open(f"{VD}/filtered_snv_{tag}_{chrom}.vcf.gz", "rt") as fh:
                for ln in fh:
                    if ln.startswith("#"):
                        continue
                    p = ln.split("\t")
                    if lo <= int(p[1]) <= hi:
                        out.append((int(p[1]), p[3], p[4], tag))
        except FileNotFoundError:
            pass
    return sorted(set(out))


def process(chrom, pos, sset):
    region = RD.get(f"{chrom}:{pos}:{sset}")
    if not region:
        return None
    anchor = int(pos) + 5000
    cand = vcf_window(chrom, anchor - 8000, anchor + 8000)
    tb = pysam.AlignmentFile(P2.TBAM, "rb"); nb = pysam.AlignmentFile(P2.NBAM, "rb")
    tcalls, thp = P2.per_read_allele(tb, chrom, [(p, r, a) for p, r, a, _ in cand])
    ncalls, _ = P2.per_read_allele(nb, chrom, [(p, r, a) for p, r, a, _ in cand])
    tb.close(); nb.close()
    # normal 確認真 somatic (normal_ALT<=1) + tumor 有變異
    snvs = []
    for p, r, a, tag in cand:
        nc = Counter(c.get(p) for c in ncalls.values() if p in c); tc = Counter(c.get(p) for c in tcalls.values() if p in c)
        if nc["ALT"] <= 1 and nc.get("REF", 0) >= 4 and tc["ALT"] >= 4:
            snvs.append({"pos": p, "ref": r, "alt": a, "tag": tag, "nREF": nc["REF"], "nALT": nc["ALT"], "tALT": tc["ALT"]})
    if len(snvs) < 2:
        return {"locus": f"{chrom}:{pos}", "set": sset, "skip": f"only {len(snvs)} confirmed somatic"}
    snvs.sort(key=lambda s: s["pos"]); SP = [s["pos"] for s in snvs]
    # reads.tsv + methylation join
    rows = open(f"{region}/reads/reads.tsv").read().splitlines(); hdr = rows[0].split("\t")
    ix = {k: hdr.index(k) for k in ("read_id", "read_name", "is_tumor", "hp")}
    name2rid = {}; rid_meta = {}
    for rr in rows[1:]:
        c = rr.split("\t"); name2rid[c[ix["read_name"]]] = c[ix["read_id"]]; rid_meta[c[ix["read_id"]]] = {"is_tumor": c[ix["is_tumor"]], "hp": c[ix["hp"]]}
    mr = open(f"{region}/methylation/methylation.csv").read().strip().split("\n"); cpgs = [int(x) for x in mr[0].split(",")[1:]]
    meth = {}
    for ln in mr[1:]:
        q = ln.split(","); meth[q[0]] = [None if v in ("", "NA", "nan", "NaN") else float(v) for v in q[1:]]
    coarse = {row["read_id"]: row.get("coarse_label", "") for row in csv.DictReader(open(f"{region}/clustering/phylo_groups.tsv"), delimiter="\t")}
    # per-read 基因型 + 甲基 (tumor, 覆蓋全 sSNV)
    reads = []
    for nm, g in tcalls.items():
        rid = name2rid.get(nm)
        if rid is None or rid not in meth or rid_meta.get(rid, {}).get("is_tumor") != "1":
            continue
        gv = tuple(g.get(p, "?") for p in SP)
        if gv.count("?") > len(SP) // 2:
            continue
        reads.append({"rid": rid, "gv": gv, "nALT": sum(1 for x in gv if x == "ALT"), "hp": rid_meta[rid]["hp"], "coarse": coarse.get(rid, ""), "meth": meth[rid]})
    if len(reads) < 8:
        return {"locus": f"{chrom}:{pos}", "set": sset, "skip": "n_reads<8"}
    # genotype-group (依 ALT 數 + vector 排序)
    groups = sorted(set(r["gv"] for r in reads), key=lambda v: (v.count("ALT"), v))
    gid = {v: i for i, v in enumerate(groups)}
    for r in reads:
        r["grp"] = gid[r["gv"]]
    # BERNOULLI 距離
    M = np.array([[np.nan if v is None else v for v in r["meth"]] for r in reads]); n = len(reads)

    def bern(i, j):
        m = ~np.isnan(M[i]) & ~np.isnan(M[j])
        if m.sum() < 3:
            return 0.5
        p, q = M[i][m], M[j][m]; w = (2 * np.abs(p - 0.5)) * (2 * np.abs(q - 0.5)); sw = w.sum()
        return float((w * (p * (1 - q) + (1 - p) * q)).sum() / sw) if sw > 1e-9 else 0.0
    Dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            Dm[i, j] = Dm[j, i] = bern(i, j)
    Z = linkage(squareform(Dm, checks=False), method="average"); order = [int(x) for x in leaves_list(Z)]
    nclu = min(len(groups), 4); clu = [int(x) for x in fcluster(Z, nclu, "maxclust")]
    # 交叉: ISM coarse × genotype-group ; UPGMA cluster × genotype-group
    cross_c = defaultdict(Counter); cross_u = defaultdict(Counter)
    for r, cl in zip(reads, clu):
        cross_c[r["coarse"]][r["grp"]] += 1; cross_u[str(cl)][r["grp"]] += 1
    # ---- 渲染 ----
    gcol = {i: PAL[i % len(PAL)] for i in range(len(groups))}

    def fig_dist():
        Do = Dm[np.ix_(order, order)]
        fig, ax = plt.subplots(figsize=(4.4, 4.0)); im = ax.imshow(Do, cmap="magma_r", vmin=0, vmax=max(0.3, Do.max()))
        for i, oi in enumerate(order):
            ax.add_patch(plt.Rectangle((-1.6, i - 0.5), 1.3, 1, color=gcol[reads[oi]["grp"]], clip_on=False))
        ax.set_xlim(-1.8, n - 0.5); ax.set_ylim(n - 0.5, -0.5); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
        ax.set_title("甲基距離(UPGMA序;邊色=基因型群)", fontsize=8.5)
        return b64(fig)

    def fig_dendro():
        fig, ax = plt.subplots(figsize=(5.4, 2.7)); dd = dendrogram(Z, ax=ax, no_labels=True, color_threshold=0, above_threshold_color=MUT)
        for i, leaf in enumerate(dd["leaves"]):
            ax.add_patch(plt.Rectangle((i * 10 + 5, -0.02), 9, 0.016, color=gcol[reads[leaf]["grp"]], clip_on=False, transform=ax.get_xaxis_transform()))
        ax.set_xticks([]); [ax.spines[s].set_visible(False) for s in ("top", "right")]
        ax.set_title("UPGMA(葉色=基因型群)", fontsize=8.5)
        return b64(fig)

    def fig_meth():
        ordr = sorted(range(n), key=lambda i: (reads[i]["grp"], np.nanmean(M[i]) if np.isfinite(M[i]).any() else 0))
        fig, ax = plt.subplots(figsize=(8, 3.4)); cm = plt.cm.RdBu_r.copy(); cm.set_bad("#e8e4dc")
        ax.imshow(M[ordr], aspect="auto", cmap=cm, vmin=0, vmax=1)
        y = 0
        for gi in range(len(groups)):
            c = sum(1 for i in ordr if reads[i]["grp"] == gi)
            if c == 0:
                continue
            if y > 0:
                ax.axhline(y - 0.5, color=TX, lw=1)
            y += c
        ax.set_xlabel(f"CpG (n={len(cpgs)})"); ax.set_yticks([]); ax.set_title("甲基熱圖(按基因型群分組)", fontsize=8.5)
        return b64(fig)

    def svg_bar():
        ordr = sorted(range(n), key=lambda i: (reads[i]["grp"], -reads[i]["nALT"]))
        rh = 7; x0 = 8; cw = max(22, min(40, 360 // len(SP)))
        H = n * rh + 26 + len(groups) * 4
        s = [f'<svg viewBox="0 0 {x0+len(SP)*cw+120} {H}" xmlns="http://www.w3.org/2000/svg" font-size="8" font-family="system-ui">']
        for k, sn in enumerate(snvs):
            s.append(f'<text x="{x0+k*cw}" y="9" font-size="7">{sn["pos"]%100000}{"*" if sn["tag"]=="fp" else ""}</text>')
        y = 16; pg = -1
        for i in ordr:
            if reads[i]["grp"] != pg:
                pg = reads[i]["grp"]; y += 3
            for k in range(len(SP)):
                al = reads[i]["gv"][k]; col = "#cfcabf" if al == "?" else ("#4A6E8A" if al == "REF" else "#C0563F")
                s.append(f'<rect x="{x0+k*cw}" y="{y}" width="{cw-2}" height="{rh-1}" fill="{col}"/>')
            mv = np.nanmean(M[i]); mv = 0 if np.isnan(mv) else mv
            s.append(f'<rect x="{x0+len(SP)*cw+6}" y="{y}" width="{mv*90:.0f}" height="{rh-1}" fill="#7a9a7a"/>')
            y += rh
        s.append("</svg>")
        return "".join(s)
    return {"locus": f"{chrom}:{pos}", "set": sset, "n_reads": n, "n_snv": len(snvs), "n_groups": len(groups),
            "snvs": snvs, "groups": [list(g) for g in groups], "fig_dist": fig_dist(), "fig_dendro": fig_dendro(),
            "fig_meth": fig_meth(), "svg": svg_bar(),
            "cross_coarse": {k: dict(v) for k, v in cross_c.items()}, "cross_upgma": {k: dict(v) for k, v in cross_u.items()},
            "lv": json.load(open(f"{A}/p2_linkage.json")).get(f"{chrom}:{pos}", {}).get("lineage_verdict", "?")}


results = [process(*l) for l in LOCI]
results = [r for r in results if r]
json.dump([{k: v for k, v in r.items() if not k.startswith("fig") and k != "svg"} for r in results],
          open(f"{A}/multilocus_observe_summary.json", "w"), ensure_ascii=False, indent=1)

# ---- gallery HTML ----
def gv_str(g, snvs):
    return " ".join(("·" if a == "?" else ("R" if a == "REF" else "A")) for a in g)


secs = []
nav = []
for r in results:
    lid = r["locus"].replace(":", "_")
    nav.append(f'<a href="#{lid}">{r["locus"]}</a>')
    if r.get("skip"):
        secs.append(f'<section id="{lid}"><h2>{r["locus"]} [{r["set"]}]</h2><p class="skip">跳過: {r["skip"]}</p></section>')
        continue
    snv_tbl = "".join(f"<tr><td>{s['pos']}</td><td>{s['ref']}→{s['alt']}</td><td>{'FP*' if s['tag']=='fp' else 'TP'}</td><td>{s['tALT']}</td><td>{s['nREF']}/{s['nALT']}</td></tr>" for s in r["snvs"])
    grp_tbl = "".join(f"<tr><td>G{i}</td><td>{gv_str(g, r['snvs'])}</td><td>{sum(1 for x in g if x=='ALT')}</td></tr>" for i, g in enumerate(r["groups"]))
    cc = "".join(f"<tr><td>{k or 'NA'}</td>" + "".join(f"<td>{r['cross_coarse'][k].get(str(gi),r['cross_coarse'][k].get(gi,0))}</td>" for gi in range(r['n_groups'])) + "</tr>" for k in r["cross_coarse"])
    secs.append(f"""<section id="{lid}"><h2>{r['locus']} [{r['set']}] <span class="badge">{r['lv']}</span></h2>
<p>{r['n_reads']} reads · {r['n_snv']} somatic SNV · {r['n_groups']} 基因型群</p>
<div class="grid"><div>{r['svg']}</div><div><img src="{r['fig_dist']}"></div></div>
<img src="{r['fig_dendro']}"><img src="{r['fig_meth']}">
<div class="tbls"><table><caption>somatic SNV (normal 確認;* =SEQC2-FP但真)</caption><tr><th>pos</th><th>ref→alt</th><th>類</th><th>tALT</th><th>nR/nA</th></tr>{snv_tbl}</table>
<table><caption>基因型群 (R=REF A=ALT)</caption><tr><th>群</th><th>基因型</th><th>#ALT</th></tr>{grp_tbl}</table>
<table><caption>ISM 甲基分類 coarse × 基因型群</caption><tr><th>coarse\\群</th>{''.join(f'<th>G{i}</th>' for i in range(r['n_groups']))}</tr>{cc}</table></div></section>""")

html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>多位點 subclone 觀察 gallery</title><style>
*{{box-sizing:border-box}} body{{margin:0;font-family:system-ui,"Noto Sans CJK TC",sans-serif;color:{TX};background:#FAF9F5;line-height:1.55}}
.wrap{{max-width:1080px;margin:0 auto;padding:0 20px 80px}}
header{{background:{TX};color:#FAF9F5;padding:24px 20px}} header h1{{margin:0;font-size:20px}} header .sub{{color:#cfc9bf;font-size:13px;margin-top:4px}}
nav{{position:sticky;top:0;background:#FAF9F5;border-bottom:1px solid {BD};padding:8px 20px;font-size:12px;z-index:9}} nav a{{color:{MUT};text-decoration:none;margin-right:11px}} nav a:hover{{color:{AC}}}
h2{{font-size:16px;margin:24px 0 8px;padding-bottom:5px;border-bottom:2px solid {AC}}}
img,svg{{max-width:100%;border:1px solid {BD};border-radius:5px;margin:5px 0;background:white}} svg{{padding:6px}}
.grid{{display:flex;gap:12px;flex-wrap:wrap}} .grid>div{{flex:1;min-width:280px}}
.tbls{{display:flex;gap:14px;flex-wrap:wrap}} table{{border-collapse:collapse;font-size:11.5px;margin:6px 0}} th,td{{border:1px solid {BD};padding:3px 7px;text-align:center}} th{{background:#efe9df}} caption{{font-size:11px;color:{MUT};padding:2px}}
.badge{{font-size:11px;background:{AC};color:white;padding:1px 8px;border-radius:10px}}
.skip{{color:{MUT};font-style:italic}} .intro{{background:white;border:1px solid {BD};border-radius:6px;padding:11px 14px;margin:12px 0;font-size:13px}}
footer{{margin-top:30px;padding-top:12px;border-top:1px solid {BD};font-size:11.5px;color:{MUT}}}
</style></head><body>
<header><div class="wrap" style="padding:0"><h1>多位點 subclone 完整觀察 gallery</h1><div class="sub">每位點：甲基距離熱圖 · UPGMA 樹 · ISM 分群×基因型交叉 · 甲基熱圖 · SVG read 關聯 · sSNV 數據（含 SEQC2-FP 但真 somatic 補回）· HCC1395 ⭐3</div></div></header>
<nav class="wrap" style="max-width:1080px">{''.join(nav)}</nav>
<div class="wrap">
<div class="intro">每個位點自動收集窗內 somatic SNV（TP + <b>SEQC2-FP 但 normal=REF 的真 somatic，標 *</b>）→ per-read 基因型分群 → BERNOULLI 距離(精確公式) + UPGMA(葉色=基因型群)。🔑 觀察重點：<b>甲基樹/ISM 分類是否 recover 基因型定義的克隆群</b>（多數只分粗、細分需遺傳錨）。</div>
{''.join(secs)}
<footer>HCC1395 ⭐3 · §13-A JSON 注入 · 腳本 p3_multilocus_observe.py · BERNOULLI src/core/DistanceMatrix.cpp:243-246 · 數據 multilocus_observe_summary.json</footer>
</div></body></html>"""
open("/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260625_multilocus_subclone_observation_01.standalone.html", "w").write(html)
print(f"loci processed: {len(results)}; " + "; ".join(f"{r['locus']}:{r.get('skip') or str(r['n_snv'])+'snv/'+str(r['n_groups'])+'grp'}" for r in results))
print("[-> 20260625_multilocus_subclone_observation_01.standalone.html]")
