#!/usr/bin/env python3
"""[chr17:48360161 IGV 式 SVG 集大成解釋] 真實基因組座標的 IGV pileup: SKAP1 gene track + 4 sSNV 標記 +
每條 read(按 lineage 分組) 的等位 tick + 甲基 CpG 點 + 克隆樹 + 方法解釋(偵測/建樹/甲基驗證) + 基因/癌症/subclone 意義。
§13-A 由 chr17_complete_data.json 注入。輸出 ../../20260625_chr17_igv_subclone_explainer_01.standalone.html。"""
import json, io, base64, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260625_chr17_igv_subclone_explainer_01.standalone.html"
D = json.load(open(f"{A}/chr17_complete_data.json"))
reads = D["reads"]; cpgs = D["cpgs"]; sn = D["snv_norm"]
TX = "#141413"; BD = "#E3DACC"; MUT = "#6B6862"; AC = "#D97757"; BLU = "#4A6E8A"; RED = "#C0563F"; PUR = "#8A6FA0"; GRN = "#5B8A5B"; TEAL = "#2E8B8B"
# read start/end from reads.tsv
region = [c["region_dir"] for c in json.load(open(f"{A}/cis_candidates_resolved.json")) if c["chrom"] == "chr17" and c["pos"] == "48360161"][0]
rmeta = {}
for row in csv.DictReader(open(f"{region}/reads/reads.tsv"), delimiter="\t"):
    rmeta[row["read_id"]] = (int(row["start"]), int(row["end"]))
for r in reads:
    r["span"] = rmeta.get(r["rid"], (None, None))

# ---- per-CpG 對齊軸 attribution + per-sSNV silhouette（已落檔驗證: chr17_tree_data.json / chr17_persnv.json）----
TD = json.load(open(f"{A}/chr17_tree_data.json"))
PSNV = json.load(open(f"{A}/chr17_persnv.json"))["snv_res"]
AXCNT = TD["best_axis_counts"]; CROSS = TD["cross_coarse_x_lineage"]
CPGAX = {p["cpg"]: p["best_axis"] for p in TD["percpg_attrib"]}  # cpg pos → best_axis|None
N_AX = sum(1 for p in TD["percpg_attrib"] if p["best_axis"]); N_CPG_ATTR = len(TD["percpg_attrib"])

# ---- stored ISM BERNOULLI 距離 (region/distance/BERNOULLI/matrix.csv) → 50 lineage reads 子集 ----
def _load_stored_dist():
    dl = open(f"{region}/distance/BERNOULLI/matrix.csv").read().strip().split("\n")
    colids = dl[0].split(",")[1:]; colidx = {c: k for k, c in enumerate(colids)}
    rowmap = {q[0]: q[1:] for q in (ln.split(",") for ln in dl[1:])}
    rids = [r["rid"] for r in reads]; n = len(rids)
    Dm = np.full((n, n), np.nan)
    for i, ri in enumerate(rids):
        rv = rowmap.get(ri)
        if rv is None:
            continue
        for j, rj in enumerate(rids):
            if rj in colidx:
                v = float(rv[colidx[rj]]); Dm[i, j] = v if v >= 0 else np.nan
    mx = np.nanmax(Dm) if np.isfinite(Dm).any() else 1.0
    Dm = np.where(np.isnan(Dm), mx, Dm); np.fill_diagonal(Dm, 0.0)
    return (Dm + Dm.T) / 2

DM = _load_stored_dist(); DLIN = [r["lineage"] for r in reads]
ZLINK = linkage(squareform(DM, checks=False), method="average")

SNVS = [("γ", 48357368, PUR), ("β1", 48362515, GRN), ("α", 48365089, BLU), ("β2", 48365161, AC)]
SNVPOS = {nm: p for nm, p, _ in SNVS}
LABS = ["ancestral(no_somatic)", "γ_subclone(sibling)", "L1_α_only(ancestor)", "L2_αβ(descendant)"]
LS = {"ancestral(no_somatic)": "ancestral", "γ_subclone(sibling)": "γ subclone", "L1_α_only(ancestor)": "L1 α-only", "L2_αβ(descendant)": "L2 α+β"}
LC = {"ancestral(no_somatic)": MUT, "γ_subclone(sibling)": PUR, "L1_α_only(ancestor)": BLU, "L2_αβ(descendant)": AC}
# best_axis 字串 → (顯示名, 顏色)：CpG 甲基對齊到哪個 sSNV REF/ALT 軸
AXMAP = {"ALT@48365089(α)": ("α", BLU), "ALT@48362515(β1)": ("β1", GRN),
         "ALT@48365161(β2)": ("β2", AC), "lineage_L1_vs_L2": ("L1→L2", TEAL)}

VIEW_LO, VIEW_HI = 48356800, 48370400
W = 1080; XPAD = 150; PW = W - XPAD - 20


def X(pos):
    return round(XPAD + (pos - VIEW_LO) / (VIEW_HI - VIEW_LO) * PW, 1)  # 1-decimal: 縮小 inline SVG 體積


def betacol(b):  # hex（比 rgb() 短 ~半）：縮小 inline SVG 體積
    if b is None:
        return "#e3ded4"
    if b < 0.5:
        t = b / 0.5
        r, g, bl = int(74 + t * 150), int(110 + t * 120), int(170 + t * 60)  # blue→light
    else:
        t = (b - 0.5) / 0.5
        r, g, bl = int(224 - t * 32), int(220 - t * 130), int(200 - t * 150)  # light→red
    return f"#{r:02x}{g:02x}{bl:02x}"


def svg_igv():
    rh = 7; lin_gap = 14; top = 92
    grp = {g: sorted([r for r in reads if r["lineage"] == g], key=lambda r: -r["nALT"] if "nALT" in r else 0) for g in LABS}
    nread = sum(len(grp[g]) for g in LABS)
    H = top + nread * rh + len([g for g in LABS if grp[g]]) * lin_gap + 30
    s = [f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" font-family="system-ui" font-size="11">']
    # 座標軸
    s.append(f'<text x="{XPAD}" y="14" font-size="11" fill="{TX}">chr17 (hg38) · SKAP1 內含子 · 17q21.32</text>')
    for t in range(48357000, 48371000, 2000):
        s.append(f'<line x1="{X(t)}" y1="20" x2="{X(t)}" y2="24" stroke="{MUT}"/><text x="{X(t)}" y="34" font-size="8.5" fill="{MUT}" text-anchor="middle">{t/1e6:.3f}M</text>')
    # 基因 track (SKAP1 內含子線 + 26bp exon)
    gy = 44
    s.append(f'<line x1="{X(VIEW_LO)}" y1="{gy}" x2="{X(VIEW_HI)}" y2="{gy}" stroke="{MUT}" stroke-width="2"/>')
    s.append(f'<rect x="{X(48363789)}" y="{gy-5}" width="{max(3,X(48363814)-X(48363789))}" height="10" fill="{MUT}"/><text x="{X(48363801)}" y="{gy-8}" font-size="8" fill="{MUT}" text-anchor="middle">exon</text>')
    s.append(f'<text x="{XPAD-6}" y="{gy+4}" font-size="10" fill="{MUT}" text-anchor="end">SKAP1 (−)</text>')
    # sSNV 標記：黃色下三角（base 上 / tip 下）+ 全高虛線 + 名稱/座標
    tri_top, tri_tip = 50, 61
    for nm, p, col in SNVS:
        x = X(p); tag = "*" if nm == "γ" else ""
        s.append(f'<text x="{x}" y="{tri_top-2}" font-size="10" fill="{col}" text-anchor="middle" font-weight="bold">{nm}{tag}</text>')
        s.append(f'<polygon points="{x-6},{tri_top} {x+6},{tri_top} {x},{tri_tip}" fill="#E8C24A" stroke="{col}" stroke-width="1.2"/>')
        s.append(f'<line x1="{x}" y1="{tri_tip}" x2="{x}" y2="{H-20}" stroke="{col}" stroke-width="0.8" stroke-dasharray="2,2" opacity="0.45"/>')
        s.append(f'<text x="{x}" y="{H-8}" font-size="7.5" fill="{col}" text-anchor="middle">{p}</text>')
    # CpG 對齊軸 track：每個 CpG 標其甲基(β)最對齊的 sSNV REF/ALT 軸（best_axis = FDR-sig 最大 |Δβ|；資料 chr17_tree_data.json）
    ax_y, ax_h = 74, 9
    s.append(f'<text x="6" y="{ax_y+ax_h-1}" font-size="9" fill="{MUT}">CpG 對齊軸</text>')
    s.append(f'<line x1="{X(VIEW_LO)}" y1="{ax_y+ax_h+1.5}" x2="{X(VIEW_HI)}" y2="{ax_y+ax_h+1.5}" stroke="{BD}" stroke-width="0.6"/>')
    for cp in cpgs:
        if not (VIEW_LO <= cp <= VIEW_HI):
            continue
        ba = CPGAX.get(cp)
        if ba and ba in AXMAP:
            s.append(f'<rect x="{X(cp)-1.5}" y="{ax_y}" width="3" height="{ax_h}" fill="{AXMAP[ba][1]}"><title>CpG {cp} → {AXMAP[ba][0]} 軸</title></rect>')
        else:
            s.append(f'<rect x="{X(cp)-0.7}" y="{ax_y+ax_h-2}" width="1.4" height="2" fill="#cfc8ba"><title>CpG {cp} 無顯著對齊軸</title></rect>')
    # read pileup (按 lineage)
    y = top
    for g in LABS:
        if not grp[g]:
            continue
        s.append(f'<text x="6" y="{y+10}" font-size="10.5" fill="{LC[g]}" font-weight="bold">{LS[g]} ({len(grp[g])})</text>')
        for r in grp[g]:
            st, en = r["span"]
            if st is None:
                continue
            x1 = X(max(st, VIEW_LO)); x2 = X(min(en, VIEW_HI))
            s.append(f'<line x1="{x1}" y1="{y+rh/2}" x2="{x2}" y2="{y+rh/2}" stroke="#d8d2c6" stroke-width="{rh-2}"/>')
            # 甲基 CpG 點
            for j, cp in enumerate(cpgs):
                b = r["meth"][j]
                if b is not None and VIEW_LO <= cp <= VIEW_HI:
                    s.append(f'<rect x="{X(cp)-1}" y="{y+1}" width="2.2" height="{rh-2}" fill="{betacol(b)}"/>')
            # sSNV 等位 tick (覆蓋 + REF/ALT)
            for nm, p, _ in SNVS:
                al = r["geno"].get(nm)
                if al in ("REF", "ALT") and st <= p <= en:
                    c = BLU if al == "REF" else RED
                    s.append(f'<rect x="{X(p)-2.5}" y="{y}" width="5" height="{rh}" fill="{c}" stroke="white" stroke-width="0.4"/>')
            y += rh
        y += lin_gap
    # 圖例
    ly = H - 4
    s.append(f'<rect x="{XPAD}" y="{ly-12}" width="9" height="9" fill="{BLU}"/><text x="{XPAD+12}" y="{ly-4}" font-size="9">sSNV REF</text>'
             f'<rect x="{XPAD+70}" y="{ly-12}" width="9" height="9" fill="{RED}"/><text x="{XPAD+82}" y="{ly-4}" font-size="9">sSNV ALT</text>'
             f'<rect x="{XPAD+150}" y="{ly-12}" width="9" height="9" fill="{betacol(0.1)}"/><rect x="{XPAD+160}" y="{ly-12}" width="9" height="9" fill="{betacol(0.9)}"/><text x="{XPAD+172}" y="{ly-4}" font-size="9">CpG 甲基 低→高</text>')
    s.append("</svg>")
    return "".join(s)


# 克隆樹 (matplotlib)
def fig_tree():
    lc = D["lineage_counts"]
    fig, ax = plt.subplots(figsize=(7.5, 3.4)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 8)
    N = {"anc": (6, 7, MUT, f"ancestral\n{lc.get('ancestral(no_somatic)',0)}"), "γ": (2.5, 4, PUR, f"γ subclone\n+48357368*\n{lc.get('γ_subclone(sibling)',0)}"),
         "α": (8.5, 4.5, BLU, f"α subclone\n+48365089\n{lc.get('L1_α_only(ancestor)',0)+lc.get('L2_αβ(descendant)',0)}"),
         "L1": (7, 1.5, BLU, f"L1 α-only\n{lc.get('L1_α_only(ancestor)',0)}"), "L2": (10.5, 1.5, AC, f"L2 α+β\n{lc.get('L2_αβ(descendant)',0)}")}
    for a, b in [("anc", "γ"), ("anc", "α"), ("α", "L1"), ("α", "L2")]:
        ax.annotate("", xy=(N[b][0], N[b][1] + 0.6), xytext=(N[a][0], N[a][1] - 0.6), arrowprops=dict(arrowstyle="-|>", color=TX, lw=1.4))
    for k, (x, y, c, l) in N.items():
        ax.add_patch(plt.Circle((x, y), 0.9, color=c, alpha=0.85)); ax.text(x, y, l, ha="center", va="center", fontsize=7, color="white")
    ax.text(3.6, 5.8, "sibling 互斥", fontsize=8, color=PUR); ax.text(10.6, 3.1, "nested", fontsize=8, color=AC)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=100, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


# 甲基距離熱圖（stored ISM BERNOULLI，UPGMA 序，上條=lineage）
def fig_dist():
    from matplotlib.gridspec import GridSpec
    order = leaves_list(ZLINK)
    fig = plt.figure(figsize=(6.0, 5.4)); gs = GridSpec(2, 1, height_ratios=[1, 24], hspace=0.015)
    axs = fig.add_subplot(gs[0]); axh = fig.add_subplot(gs[1])
    axh.imshow(DM[np.ix_(order, order)], cmap="magma_r", aspect="auto")
    axh.set_xticks([]); axh.set_yticks([]); axh.set_xlabel("reads (UPGMA 序)", fontsize=9); axh.set_ylabel("reads", fontsize=9)
    for x, i in enumerate(order):
        axs.add_patch(plt.Rectangle((x - 0.5, 0), 1, 1, color=LC[DLIN[i]]))
    axs.set_xlim(-0.5, len(order) - 0.5); axs.set_ylim(0, 1); axs.axis("off")
    axs.set_title("BERNOULLI 距離熱圖（stored ISM）· 上條=lineage 色", fontsize=9)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=100, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


# UPGMA dendrogram（葉按 lineage 上色）
def fig_dendro():
    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    dn = dendrogram(ZLINK, ax=ax, no_labels=True, color_threshold=0, above_threshold_color=MUT)
    for x, idx in enumerate(dn["leaves"]):
        ax.plot(5 + 10 * x, 0, marker="s", color=LC[DLIN[idx]], ms=5, clip_on=False)
    ax.set_ylabel("BERNOULLI dist", fontsize=9); ax.set_xticks([]); ax.set_title("UPGMA 樹（葉=lineage 色）", fontsize=9)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=100, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


igv = svg_igv(); tree = fig_tree(); dist = fig_dist(); dendro = fig_dendro()
n_methA = sum(CROSS["1-2"].values()); n_methB = sum(CROSS["1-1"].values())
sil_rows = "".join(f"<tr><td>{nm}</td><td>{PSNV[nm]['silhouette']}</td><td>{PSNV[nm]['n_sig_cpg']}</td><td>{PSNV[nm]['median_abs_db']}</td><td>{'✅' if PSNV[nm]['corroborated'] else '—'}</td></tr>" for nm in ['γ', 'α', 'β1', 'β2'])
ax_legend = "".join(f'<span style="display:inline-block;width:11px;height:11px;background:{c};border-radius:2px;margin:0 3px 0 12px;vertical-align:middle"></span>{disp} 軸 ({AXCNT.get(k,0)})' for k, (disp, c) in AXMAP.items())
lin_order = ["L0", "L1", "L2", "other"]
cross_rows = "".join(f'<tr><td>{("1-2 (8條,L0為主)" if cl=="1-2" else "1-1 (42條,L1+L2混)")}</td>' + "".join(f"<td>{CROSS[cl].get(l,0)}</td>" for l in lin_order) + "</tr>" for cl in ["1-2", "1-1"])
snv_rows = "".join(f"<tr><td>{nm}{'*' if nm=='γ' else ''}</td><td>{v['pos']}</td><td>{v['ref']}→{v['alt']}</td><td>{v['tumor_REF']}/{v['tumor_ALT']}</td><td>{v['normal_REF']}/{v['normal_ALT']}</td><td>{'⚠ ONT候選' if nm=='γ' else '✅ HighConf'}</td></tr>" for nm, v in zip(['γ','β1','α','β2'], [sn['γ'], sn['β1'], sn['α'], sn['β2']]))
pairs = json.load(open(f"{A}/p2_linkage.json"))["chr17:48360161"]["pairs"]
html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>chr17:48360161 IGV-style subclone explainer</title><style>
:root{{--ac:{AC};--tx:{TX};--bd:{BD};--mut:{MUT};--grn:{GRN};--red:{RED};--blu:{BLU}}}
*{{box-sizing:border-box}} body{{margin:0;font-family:system-ui,"Noto Sans CJK TC",sans-serif;color:{TX};background:#FAF9F5;line-height:1.7}}
.wrap{{max-width:1140px;margin:0 auto;padding:0 22px 80px}}
header{{background:{TX};color:#FAF9F5;padding:28px 22px}} header h1{{margin:0;font-size:21px}} header .sub{{color:#cfc9bf;font-size:13.5px;margin-top:5px}}
nav{{position:sticky;top:0;background:#FAF9F5;border-bottom:1px solid {BD};padding:9px 22px;font-size:12.5px;z-index:9}} nav a{{color:{MUT};text-decoration:none;margin-right:12px}} nav a:hover{{color:{AC}}}
h2{{font-size:18px;margin:30px 0 10px;padding-bottom:5px;border-bottom:2px solid {AC}}} h2 .n{{color:{AC};margin-right:8px}}
.igv{{overflow-x:auto;border:1px solid {BD};border-radius:8px;background:white;padding:8px;margin:12px 0}}
img{{max-width:100%;border:1px solid {BD};border-radius:6px;margin:8px 0;background:white}}
table{{border-collapse:collapse;font-size:13px;margin:8px 0}} th,td{{border:1px solid {BD};padding:5px 9px;text-align:center}} th{{background:#efe9df}} td.mini{{display:inline-table}}
.mini{{display:inline-table;margin:5px 12px 5px 0;font-size:12px}}
.box{{background:white;border:1px solid {BD};border-left:4px solid {AC};border-radius:0 6px 6px 0;padding:12px 16px;margin:13px 0}} .red{{border-left-color:{RED};background:#fbf0ec}} .ok{{border-left-color:{GRN};background:#eef4ee}}
.step{{display:flex;gap:12px;margin:10px 0}} .step .num{{flex:0 0 30px;height:30px;background:{AC};color:white;border-radius:50%;display:flex;align-items:center;justify-content:center;font-weight:bold}} .step .body{{flex:1}}
footer{{margin-top:32px;padding-top:14px;border-top:1px solid {BD};font-size:12px;color:{MUT}}} .tag{{display:inline-block;background:{AC};color:white;font-size:11px;padding:1px 7px;border-radius:10px;margin-left:6px}}
b.k{{color:{AC}}}
</style></head><body>
<header><div class="wrap" style="padding:0"><h1>chr17:48360161 — 多 sSNV 連鎖 + 甲基 + subclone 重建（IGV 式解釋）<span class="tag">⭐3 HCC1395</span></h1>
<div class="sub">SKAP1 內含子 · 真實基因組座標 IGV pileup · 偵測→建樹→甲基驗證 全流程 · 基因/癌症/subclone 意義</div></div></header>
<nav class="wrap" style="max-width:1140px"><a href="#igv">①IGV 圖</a><a href="#detect">②如何偵測連鎖</a><a href="#tree">③如何建分類樹</a><a href="#meth">④甲基輔助驗證</a><a href="#dist">④B 距離聚類/對齊軸</a><a href="#gene">⑤基因/癌症/意義</a></nav>
<div class="wrap">

<section id="igv"><h2><span class="n">①</span>IGV 式 pileup（真實基因組座標）</h2>
<p>橫軸 = chr17 基因組座標（SKAP1 內含子）。每列一條 tumor read（按基因型 lineage 分組）：<b class="k">方塊 = 該 read 在 4 個 sSNV 的等位</b>（藍 REF / 紅 ALT），<b class="k">細條 = read 上每個 CpG 的甲基</b>（藍 unmeth → 紅 meth）。一眼看出：γ-ALT 那群（紫）在 α/β 全 REF；α-ALT 群在 β 分 REF/ALT。</p>
<div class="igv">{igv}</div>
<table><tr><th>sSNV</th><th>pos (hg38)</th><th>ref→alt</th><th>tumor R/A</th><th>normal R/A</th><th>SEQC2 狀態</th></tr>{snv_rows}</table>
<p style="font-size:12px;color:{MUT}">* γ=48357368 不在 SEQC2 HighConf/superSet（在 callable HC region 內），ONT normal=REF + VAF 0.18 → <b>ONT 偵測候選，非金標準確認</b>；α/β1/β2 皆 SEQC2 HighConf PASS。</p></section>

<section id="detect"><h2><span class="n">②</span>如何偵測 sSNV 並建立關聯</h2>
<div class="step"><div class="num">1</div><div class="body"><b>取 per-read 等位</b>：pysam 對每個 somatic SNV 位置，逐 read 讀出 base → REF/ALT（MAPQ≥20）。同一條 long read 覆蓋多個 sSNV → 得該 read 的多位點基因型向量。</div></div>
<div class="step"><div class="num">2</div><div class="body"><b>normal 確認 somatic</b>：normal BAM 同位置必須 REF（≤1 ALT 容錯）→ 排除 germline。此處 α/β/γ normal 皆 ~0 ALT。</div></div>
<div class="step"><div class="num">3</div><div class="body"><b>建共現 2×2 表</b>：對每對 sSNV，數覆蓋兩者的 read 的 (REF/ALT)×(REF/ALT) 四格。例 α×β2：</div></div>
<div>{''.join(f'<table class="mini"><caption>{k}</caption><tr><th></th><th>R</th><th>A</th></tr><tr><th>R</th><td>{v["REF_REF"]}</td><td>{v["REF_ALT"]}</td></tr><tr><th>A</th><td>{v["ALT_REF"]}</td><td>{v["ALT_ALT"]}</td></tr></table>' for k, v in pairs.items())}</div>
<div class="box"><b>關聯判讀（從 2×2）</b>：<b class="k">互斥</b>(ALT_ALT=0,離對角>0)=兩突變在不同細胞=sibling；<b class="k">嵌套</b>(一方向=0,ALT_ALT>0)=一突變 ⊂ 另一=祖先-後代；<b class="k">共連</b>(兩離對角=0)=同 lineage。🔴 須<b>同 germline HP</b> 才算 subclone（異 HP=等位）；此處 17q LOH → 只剩單一 haplotype → 互斥必為細胞層級。</div></section>

<section id="tree"><h2><span class="n">③</span>如何建立分類（克隆）樹</h2>
<img src="{tree}" alt="tree" style="max-width:680px">
<p>把 2×2 關係串起來：γ 與 α <b>互斥</b> → 兩 sibling 分支；β1/β2 <b>嵌套</b>在 α 內且彼此<b>共連</b> → α 的後代。得 <b class="k">4 細胞群</b>：ancestral / γ subclone / L1(α) / L2(α+β)。VAF 階層 α 0.82 > β 0.47 > γ 0.18 與拓樸一致（祖先廣、後代次、minor sibling 最小）。</p></section>

<section id="meth"><h2><span class="n">④</span>如何用甲基建立輔助驗證資訊</h2>
<p>克隆群一旦由<b>遺傳（sSNV）定義好</b>，甲基就能<b>非循環</b>地 characterize：</p>
<div class="box ok"><b>甲基輔助驗證的兩個用途</b>：<br>
1. <b>佐證分群</b>：每個 sSNV 的 REF/ALT 分群若被甲基距離聚類支持（silhouette），加強信心——此處 γ silhouette 0.504 最強、α 0.372，全 4 sSNV corroborated。<br>
2. <b>刻畫 lineage</b>：祖先 L1 vs 後代 L2 有 16 個 CpG 甲基顯著不同；per-CpG 可歸因（23 CpG→α、6→L1-L2 轉變）。IGV 圖上每條 read 的甲基細條即此資訊。</div>
<div class="red">🔴 <b>誠實邊界</b>：甲基<b>不能單獨</b>定 subclone（genotype-同質群內甲基切群=double-dip 循環）；ISM 無監督甲基 clustering 在此只切得出「有無 somatic」，切不出 L1↔L2。甲基是<b>輔助/characterize</b>，遺傳 sSNV 連鎖才是非循環錨。</div></section>

<section id="dist"><h2><span class="n">④B</span>甲基距離聚類 + UPGMA + CpG 對齊軸</h2>
<p>用 <b>stored ISM BERNOULLI 距離矩陣</b>（<code>region/distance/BERNOULLI/matrix.csv</code>，50 條有 lineage 的 read）做 UPGMA。葉/上條按<b>遺傳 lineage</b> 上色，直接看「甲基距離能否分出克隆群」：</p>
<div style="display:flex;gap:16px;flex-wrap:wrap;align-items:flex-start">
<div><img src="{dist}" alt="distance heatmap" style="max-width:360px"></div>
<div style="flex:1;min-width:330px"><img src="{dendro}" alt="UPGMA dendrogram" style="max-width:100%">
<table style="margin-top:8px"><caption style="font-size:12px;color:{MUT};text-align:left">無監督甲基 coarse cluster × 遺傳 lineage（reads 數）</caption>
<tr><th>甲基 cluster</th><th>L0</th><th>L1</th><th>L2</th><th>other</th></tr>{cross_rows}</table></div></div>
<div class="red">🔴 <b>關鍵觀察</b>：甲基樹只把 <b>L0（無 somatic）</b>分出去；<b>L1 祖先（藍）與 L2 後代（橘）在甲基空間混在一起</b>。無監督甲基 clustering 只切出 cluster 1-2（{n_methA} 條，L0 為主）vs 1-1（{n_methB} 條，L1+L2 混）——<b>切得出「有沒有 somatic」，切不出 L1↔L2 subclone</b>。subclone 細分必須靠遺傳錨（sSNV 連鎖）。</div>

<h3 style="font-size:15px;margin:20px 0 8px">CpG 對齊軸：哪些甲基位點佐證哪個 sSNV 的 REF/ALT 分群</h3>
<p>IGV 圖頂的 <b>CpG 對齊軸</b> track（基因track 與 read pileup 之間）：每個 CpG 標其甲基（β）<b>最對齊</b>的 sSNV REF/ALT 軸 = <code>best_axis</code>（FDR-顯著、最大 |Δβ|）。色：{ax_legend}。<br>共 <b>{N_AX}/{N_CPG_ATTR}</b> 個 CpG 有顯著對齊軸 → 即「能以某個 sSNV 的 REF/ALT 佐證聚類結構」的甲基位點；其餘灰點 = 無顯著對齊（不偏向任何 sSNV 軸）。<b>α 軸吸走最多（{AXCNT.get('ALT@48365089(α)',0)} 個）</b>＝祖先 somatic 的 cis-ASM 訊號最強；<b>{AXCNT.get('lineage_L1_vs_L2',0)} 個對齊 L1→L2 lineage</b>＝後代 subclone 特有甲基。</p>
<p>每個 sSNV 的 REF/ALT 分群被甲基距離 silhouette 佐證的程度（全 4 sSNV corroborated）：</p>
<table><tr><th>sSNV</th><th>silhouette</th><th>顯著 CpG 數</th><th>median |Δβ|</th><th>corroborated</th></tr>{sil_rows}</table>
<p style="font-size:12px;color:{MUT}">🔴 <b>對齊軸 ≠ 偵測 subclone</b>：CpG 對齊到 α 軸多半是 cis-ASM（甲基隨 somatic 共分離），是 characterize / 佐證分群信心，<b>非</b>用甲基偵測 subclone；非循環錨仍是 sSNV 連鎖。</p></section>

<section id="gene"><h2><span class="n">⑤</span>基因資訊 + 癌症影響 + subclone 意義</h2>
<div class="box"><b>基因 = SKAP1</b>（Src kinase-associated phosphoprotein 1，ENSG00000141293，17q21.32，負股）。4 個 sSNV <b>全在內含子</b>（唯一 26bp 外顯子 48363789-814 無命中）→ <b>passenger 突變</b>（lineage tracing 標準標記，非功能）。</div>
<div class="box"><b>癌症影響（SKAP1 基因層）</b>：文獻報 oncogene-like / 預後不良（胃癌、大腸癌、<b>乳癌</b>），機制 JAK1/PI3K/AKT；但<b>非 TNBC/HCC1395-specific、非 census driver</b>。🔴 <b>SKAP1 DNA 甲基化 × 癌症 = 文獻空白</b>（無人研究）。<b>內含子 passenger 突變 ≠ SKAP1 促癌功能</b>，兩條獨立軸不可相乘。</div>
<div class="box"><b>subclone 意義與應用</b>：HCC1395 有 10 subclone 分支演化（Fang2021），S1=MRCA 帶 clonal driver（TP53/PIK3CA/17q NLOH）。此位點在 <b>17q clonal LOH</b> 之上——LOH 是主幹、α/β/γ somatic 定義其上的<b>局部 subclonal 分支</b>。<b>重要性</b>：示範「somatic haplotagging（sSNV 連鎖）重建局部克隆樹 + methylation profiles characterize」——即論文標題。<b>應用</b>：long-read 直接讀出 read 級克隆結構（短讀共識可能漏低-CCF subclone，如 γ）；但⭐3 單樣本只能 characterize，confirmation 需 single-cell/multi-region。</div></section>

<footer>chr17:48360161 · SKAP1 內含子 · HCC1395 ⭐3 · §13-A 由 chr17_complete_data.json + p2_linkage.json + SEQC2 本地驗證注入<br>
γ=48357368(ONT候選,SEQC2未收) · α=48365089/β1=48362515/β2=48365161(SEQC2 HighConf) · 腳本 build_chr17_igv.py · 來源 Fang2021 Nat Biotechnol / SKAP1 PMC10380396</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes, svg + tree")
