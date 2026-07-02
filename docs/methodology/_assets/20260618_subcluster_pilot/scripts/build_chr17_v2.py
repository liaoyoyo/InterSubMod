#!/usr/bin/env python3
"""[chr17:48360161 完整觀察工作站 v2] SVG per-read sSNV 關聯 + BERNOULLI 距離熱圖 + UPGMA 樹(按lineage上色)
+ 甲基分類×基因型lineage 交叉 + per-CpG 基因軸歸因 + 甲基熱圖 + 完整數據。§13-A 由 JSON 注入。
輸出 ../../20260625_chr17_48360161_subclone_workstation_01.standalone.html。"""
import json, io, base64
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import numpy as np

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260625_chr17_48360161_subclone_workstation_01.standalone.html"
D = json.load(open(f"{A}/chr17_subclone_data.json"))
T = json.load(open(f"{A}/chr17_tree_data.json"))
LK = json.load(open(f"{A}/p2_linkage.json"))["chr17:48360161"]
reads = D["reads"]; cpgs = D["cpgs"]; n = len(reads)
AC = "#D97757"; TX = "#141413"; BD = "#E3DACC"; MUT = "#6B6862"; GRN = "#5B8A5B"; RED = "#C0563F"; BLU = "#4A6E8A"
LCOL = {"L0": MUT, "L1": BLU, "L2": AC, "other": "#cfcabf"}
LINS = {"L0_ancestral_root": "L0", "L1_alpha_only(ancestor)": "L1", "L2_alpha_beta(descendant)": "L2", "other": "other"}
lin = [LINS[r["lineage"]] for r in reads]
plt.rcParams.update({"font.size": 10, "text.color": TX, "axes.edgecolor": BD, "figure.facecolor": "white", "axes.facecolor": "white"})


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=115, bbox_inches="tight"); plt.close(fig); return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


# ---- SVG: per-read sSNV 關聯 barcode ----
def svg_linkage():
    order = ["L0", "L1", "L2", "other"]
    grp = {g: [r for r, l in zip(reads, lin) if l == g] for g in order}
    rh = 9; gap = 18; x0 = 95; cw = 46
    snvs = [("48365089", "α"), ("48365161", "β2"), ("48362515", "β1")]
    H = sum(len(grp[g]) * rh + gap for g in order if grp[g]) + 40
    s = [f'<svg viewBox="0 0 560 {H}" xmlns="http://www.w3.org/2000/svg" font-family="system-ui" font-size="10">']
    s.append(f'<text x="{x0}" y="14" font-size="11">α 48365089</text><text x="{x0+cw}" y="14" font-size="11">β2 48365161</text><text x="{x0+2*cw}" y="14" font-size="11">β1 48362515</text><text x="{x0+3*cw+14}" y="14" font-size="11">methylation mean</text>')
    y = 24
    for g in order:
        if not grp[g]:
            continue
        s.append(f'<text x="6" y="{y+12}" fill="{LCOL[g]}" font-weight="bold">{g} (n={len(grp[g])})</text>')
        for r in sorted(grp[g], key=lambda r: np.nanmean([np.nan if v is None else v for v in r["meth"]]) if any(v is not None for v in r["meth"]) else 0):
            for k, (pos, _) in enumerate(snvs):
                al = r["geno"].get(pos)
                col = "#cfcabf" if al not in ("REF", "ALT") else (BLU if al == "REF" else RED)
                s.append(f'<rect x="{x0+k*cw}" y="{y}" width="{cw-3}" height="{rh-1.5}" fill="{col}"/>')
            mv = np.nanmean([np.nan if v is None else v for v in r["meth"]])
            mv = 0 if np.isnan(mv) else mv
            s.append(f'<rect x="{x0+3*cw+14}" y="{y}" width="{mv*120:.1f}" height="{rh-1.5}" fill="#7a9a7a"/>')
            y += rh
        y += gap
    s.append(f'<rect x="{x0}" y="{H-16}" width="11" height="11" fill="{BLU}"/><text x="{x0+15}" y="{H-7}">REF</text><rect x="{x0+55}" y="{H-16}" width="11" height="11" fill="{RED}"/><text x="{x0+70}" y="{H-7}">ALT</text><rect x="{x0+110}" y="{H-16}" width="11" height="11" fill="#cfcabf"/><text x="{x0+125}" y="{H-7}">NA</text>')
    s.append("</svg>")
    return "".join(s)


# ---- BERNOULLI 距離熱圖 (UPGMA 序) ----
def fig_dist():
    Dm = np.array(T["distance_matrix"]); order = T["read_order_upgma"]
    Dord = Dm[np.ix_(order, order)]
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    im = ax.imshow(Dord, cmap="magma_r", vmin=0, vmax=max(0.3, Dord.max()))
    # lineage 色條
    for i, oi in enumerate(order):
        ax.add_patch(plt.Rectangle((-1.8, i - 0.5), 1.5, 1, color=LCOL[lin[oi]], clip_on=False))
        ax.add_patch(plt.Rectangle((i - 0.5, -1.8), 1, 1.5, color=LCOL[lin[oi]], clip_on=False))
    ax.set_xlim(-2, n - 0.5); ax.set_ylim(n - 0.5, -2); ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02).set_label("BERNOULLI distance")
    ax.set_title("read×read 甲基距離熱圖 (UPGMA 序; 邊色=基因型 lineage)", fontsize=10, color=TX)
    return b64(fig)


# ---- UPGMA 樹 (葉按 lineage 上色) ----
def fig_dendro():
    Dm = np.array(T["distance_matrix"]); Z = linkage(squareform(Dm, checks=False), method="average")
    fig, ax = plt.subplots(figsize=(11, 3.4))
    dd = dendrogram(Z, ax=ax, no_labels=True, color_threshold=0, above_threshold_color=MUT)
    for i, leaf in enumerate(dd["leaves"]):
        ax.add_patch(plt.Rectangle((i * 10 + 5, -0.015), 9, 0.012, color=LCOL[lin[leaf]], clip_on=False, transform=ax.get_xaxis_transform()))
    ax.set_ylabel("BERNOULLI dist"); ax.set_xticks([])
    ax.set_title("UPGMA 樹（葉色=基因型 lineage L0/L1/L2）→ 甲基樹只分出 L0，不分 L1/L2", fontsize=10, color=TX)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return b64(fig)


# ---- 甲基熱圖 (按 lineage 分組) ----
def fig_methheat():
    order = ["L0", "L1", "L2", "other"]; ordered = []
    for g in order:
        grp = [r for r, l in zip(reads, lin) if l == g]
        grp.sort(key=lambda r: np.nanmean([np.nan if v is None else v for v in r["meth"]]) if any(v is not None for v in r["meth"]) else 0)
        ordered += [(g, r) for r in grp]
    mat = np.array([[np.nan if v is None else v for v in r["meth"]] for _, r in ordered])
    fig, ax = plt.subplots(figsize=(11, 5))
    cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#e8e4dc")
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    y = 0
    for g in order:
        c = sum(1 for l, _ in ordered if l == g)
        if c == 0:
            continue
        if y > 0:
            ax.axhline(y - 0.5, color=TX, lw=1.2)
        ax.text(-2.5, y + c / 2, g, va="center", ha="right", fontsize=10, color=LCOL[g], fontweight="bold")
        y += c
    ax.set_xlabel(f"CpG sites (n={len(cpgs)})"); ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01).set_label("β")
    ax.set_title("per-read × CpG 甲基熱圖（reads 按基因型 lineage 分組）", fontsize=10, color=TX)
    return b64(fig)


# ---- per-CpG 基因軸歸因 ----
def fig_attrib():
    pc = T["percpg_attrib"]; bc = T["best_axis_counts"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.2), gridspec_kw={"width_ratios": [1, 2.4]})
    axn = list(bc.keys()); vals = [bc[k] for k in axn]
    cmap = {"ALT@48365089(α)": BLU, "lineage_L1_vs_L2": AC, "ALT@48365161(β2)": GRN, "ALT@48362515(β1)": "#C98A5B"}
    a1.barh(range(len(axn)), vals, color=[cmap.get(k, MUT) for k in axn])
    a1.set_yticks(range(len(axn))); a1.set_yticklabels([k.replace("ALT@", "").replace("lineage_", "") for k in axn], fontsize=8)
    for i, v in enumerate(vals):
        a1.text(v + 0.2, i, str(v), va="center", fontsize=8)
    a1.set_title("每 CpG 最佳歸因軸 (n)", fontsize=9.5); a1.spines[["top", "right"]].set_visible(False)
    # per-CpG strip along genome
    for p in pc:
        ba = p["best_axis"]
        col = cmap.get(ba, "#e8e4dc") if ba else "#e8e4dc"
        a2.add_patch(plt.Rectangle((p["cpg"], 0), 60, 1, color=col))
    a2.set_xlim(min(c["cpg"] for c in pc) - 100, max(c["cpg"] for c in pc) + 100); a2.set_ylim(0, 1)
    a2.set_yticks([]); a2.set_xlabel("CpG 基因座標")
    a2.set_title("每 CpG 沿基因組的最佳歸因軸（藍=α / 橘=lineage / 綠=β2 / 棕=β1 / 灰=無）", fontsize=9.5)
    return b64(fig)


charts = {"dist": fig_dist(), "dendro": fig_dendro(), "methheat": fig_methheat(), "attrib": fig_attrib()}
svg = svg_linkage()

# tables
pairs = LK["pairs"]
def tbl(pk, pv):
    a, b = pk.split("_")
    return (f'<table class="mini"><caption>{a}×{b} <b>{pv["rel"]}</b></caption><tr><th></th><th>{b}R</th><th>{b}A</th></tr>'
            f'<tr><th>{a}R</th><td>{pv["REF_REF"]}</td><td>{pv["REF_ALT"]}</td></tr><tr><th>{a}A</th><td>{pv["ALT_REF"]}</td><td>{pv["ALT_ALT"]}</td></tr></table>')
pair_html = "".join(tbl(k, v) for k, v in pairs.items())


def cross_html(ct, title):
    lins = ["L0", "L1", "L2", "other"]
    rows = "".join(f"<tr><th>{k}</th>" + "".join(f"<td>{ct[k].get(l,0)}</td>" for l in lins) + "</tr>" for k in ct)
    return f'<table><caption>{title}</caption><tr><th>cluster\\lineage</th>' + "".join(f"<th>{l}</th>" for l in lins) + f"</tr>{rows}</table>"


cc = cross_html(T["cross_coarse_x_lineage"], "ISM 甲基分類 coarse_label × 基因型 lineage")
asig = T["axis_sig_cpg_count"]
asig_rows = "".join(f"<tr><td>{k.replace('ALT@','').replace('lineage_','')}</td><td>{v}</td></tr>" for k, v in asig.items())
html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>chr17:48360161 subclone 觀察工作站</title><style>
:root{{--ac:{AC};--tx:{TX};--bd:{BD};--mut:{MUT};--grn:{GRN};--red:{RED}}}
*{{box-sizing:border-box}} body{{margin:0;font-family:system-ui,"Noto Sans CJK TC","PingFang TC",sans-serif;color:var(--tx);background:#FAF9F5;line-height:1.6}}
.wrap{{max-width:1120px;margin:0 auto;padding:0 22px 80px}}
header{{background:var(--tx);color:#FAF9F5;padding:28px 22px;margin-bottom:8px}} header h1{{margin:0 0 5px;font-size:22px}} header .sub{{color:#cfc9bf;font-size:13.5px}}
nav{{position:sticky;top:0;background:#FAF9F5;border-bottom:1px solid var(--bd);padding:9px 22px;font-size:12.5px;z-index:9;margin-bottom:14px}} nav a{{color:var(--mut);text-decoration:none;margin-right:12px}} nav a:hover{{color:var(--ac)}}
h2{{font-size:17px;margin:28px 0 9px;padding-bottom:5px;border-bottom:2px solid var(--ac)}}
img,svg{{max-width:100%;border:1px solid var(--bd);border-radius:6px;margin:8px 0;background:white}}
svg{{padding:8px}}
table{{border-collapse:collapse;font-size:12.5px;margin:8px 0}} th,td{{border:1px solid var(--bd);padding:4px 8px;text-align:center}} th{{background:#efe9df}} caption{{font-size:12px;color:var(--mut);padding:3px;caption-side:top}}
table.mini{{display:inline-table;margin:6px 12px 6px 0}}
.ok{{border-left:4px solid var(--grn);background:#eef4ee;padding:11px 14px;border-radius:0 6px 6px 0;margin:12px 0}}
.red{{border-left:4px solid var(--red);background:#fbf0ec;padding:11px 14px;border-radius:0 6px 6px 0;margin:12px 0}}
.note{{background:white;border:1px solid var(--bd);border-radius:6px;padding:10px 14px;margin:10px 0;font-size:13.5px}}
details{{margin:10px 0}} summary{{cursor:pointer;color:var(--ac);font-weight:bold}}
footer{{margin-top:34px;padding-top:14px;border-top:1px solid var(--bd);font-size:12px;color:var(--mut)}}
.tag{{display:inline-block;background:var(--ac);color:white;font-size:11px;padding:1px 7px;border-radius:10px;margin-left:6px}}
</style></head><body>
<header><div class="wrap" style="padding:0"><h1>chr17:48360161 — subclone 觀察工作站 <span class="tag">⭐3 HCC1395</span></h1>
<div class="sub">read 關聯(SVG) · 甲基距離熱圖 · UPGMA 樹 · 甲基分類×基因型交叉 · per-CpG 基因軸歸因 · 完整數據</div></div></header>
<nav class="wrap" style="max-width:1120px"><a href="#link">①read關聯</a><a href="#tree">②距離/UPGMA</a><a href="#cls">③分類交叉</a><a href="#heat">④甲基熱圖</a><a href="#attr">⑤per-CpG歸因</a><a href="#data">⑥完整數據</a></nav>
<div class="wrap">

<section id="link"><h2>① read 關聯（SVG，每條 read 的 sSNV 基因型）</h2>
<p>每列一條 tumor read（按基因型 lineage 分組）；3 格 = α(48365089)/β2(48365161)/β1(48362515)，藍=REF 紅=ALT；右綠條=該 read 甲基均值。直接看出分支：<b>L0 全 REF（無 somatic）→ L1 只 α-ALT（祖先）→ L2 α+β 全 ALT（後代）</b>。</p>
{svg}
<div class="note">🔑 連鎖關係：α-ALT 的 read 分兩支——一支 β=REF（L1 祖先，{LK['pairs'].get('48365089_48365161',{}).get('ALT_REF','?')} 條）、一支 β=ALT（L2 後代，{LK['pairs'].get('48365089_48365161',{}).get('ALT_ALT','?')} 條）；β 從不單獨出現（REF_ALT=0）→ nested 階層。</div></section>

<section id="tree"><h2>② 甲基距離熱圖 + UPGMA 樹</h2>
<img src="{charts['dist']}" alt="dist"><img src="{charts['dendro']}" alt="dendro">
<div class="red">🔴 <b>關鍵觀察</b>：BERNOULLI 距離精確重算（DistanceMatrix.cpp 公式）；UPGMA 樹的葉按基因型 lineage 上色 → <b>甲基樹只把 L0（無 somatic）分出去，L1 祖先與 L2 後代在甲基空間混在一起</b>。即<b>無監督甲基 clustering 切得出「有沒有 somatic」，切不出 L1↔L2 subclone</b>——subclone 的細分必須靠遺傳錨（sSNV 連鎖）。</div></section>

<section id="cls"><h2>③ ISM 甲基分類 × 基因型 lineage 交叉</h2>
{cc}
<p>ISM phylo coarse_label 把 reads 分成 <b>1-2(8 條，L0 為主)</b> 與 <b>1-1(42 條，L1+L2 混)</b> → 證實 ②：甲基分類對齊 L0-vs-somatic，不對齊 L1-vs-L2。</p></section>

<section id="heat"><h2>④ per-read × CpG 甲基熱圖（按 lineage 分組）</h2>
<img src="{charts['methheat']}" alt="methheat">
<p>reads 用<b>基因型</b>分組後，可見甲基模式隨 lineage（尤其 L0 vs L1/L2）變化 = 甲基 characterize 基因型定義的 lineage。</p></section>

<section id="attr"><h2>⑤ per-CpG 基因軸歸因（哪些 CpG 關連哪個 sSNV/HP）</h2>
<img src="{charts['attrib']}" alt="attrib">
<div class="ok"><b>定位結果</b>：差異 CpG 可逐一歸因到基因軸 — <b>{T['best_axis_counts'].get('ALT@48365089(α)',0)} 個關連 α(48365089 REF/ALT)</b>（祖先 somatic）、<b>{T['best_axis_counts'].get('lineage_L1_vs_L2',0)} 個關連 L1→L2 lineage 轉變</b>（後代 subclone）、少數關連 β2/β1。</div>
<table><caption>每軸顯著 CpG 數（MWU FDR<0.05 且 |Δβ|≥0.2）</caption><tr><th>基因軸</th><th>顯著 CpG 數</th></tr>{asig_rows}</table>
<p>🔴 HP1 vs HP1-1 = 0（此位點 HP1 僅 4 條 read，無力測；結構主要由 sSNV lineage 解釋）。</p></section>

<section id="data"><h2>⑥ 完整數據觀察</h2>
{pair_html}
<details><summary>展開：每條 read 的基因型 + lineage + HP + 甲基均值（{n} reads）</summary>
<table><tr><th>read_id</th><th>α 89</th><th>β2 161</th><th>β1 2515</th><th>lineage</th><th>hp</th><th>甲基cluster</th><th>β均值</th></tr>
{''.join(f"<tr><td>{r['rid']}</td><td>{r['geno'].get('48365089')}</td><td>{r['geno'].get('48365161')}</td><td>{r['geno'].get('48362515')}</td><td>{LINS[r['lineage']]}</td><td>{r['hp']}</td><td>{r.get('coarse','')}</td><td>{round(float(np.nanmean([np.nan if v is None else v for v in r['meth']])),2) if any(v is not None for v in r['meth']) else 'NA'}</td></tr>" for r in reads)}
</table></details></section>

<footer>chr17:48360161 · HCC1395 ⭐3 · §13-A 由 chr17_subclone_data.json + chr17_tree_data.json + p2_linkage.json 注入<br>
腳本 p3_chr17_extract.py + p3_chr17_full.py + build_chr17_v2.py · BERNOULLI 公式 src/core/DistanceMatrix.cpp:243-246</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes")
