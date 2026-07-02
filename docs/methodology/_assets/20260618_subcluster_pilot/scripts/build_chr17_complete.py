#!/usr/bin/env python3
"""[chr17:48360161 完整 4-sSNV subclone 定版 HTML] 完整克隆樹(γ sibling + α→α+β nested) + SVG per-read 4-sSNV 關聯
+ 6 對 2×2 + per-lineage 甲基熱圖 + normal 確認 + SEQC2 未收 γ（ONT 候選）方法學洞見。§13-A JSON 注入。
輸出 ../../20260625_chr17_complete_subclone_example_01.standalone.html。"""
import json, io, base64
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260625_chr17_complete_subclone_example_01.standalone.html"
D = json.load(open(f"{A}/chr17_complete_data.json"))
reads = D["reads"]; cpgs = D["cpgs"]; lc = D["lineage_counts"]
AC = "#D97757"; TX = "#141413"; BD = "#E3DACC"; MUT = "#6B6862"; GRN = "#5B8A5B"; RED = "#C0563F"; BLU = "#4A6E8A"; PUR = "#8A6FA0"
LABS = ["ancestral(no_somatic)", "γ_subclone(sibling)", "L1_α_only(ancestor)", "L2_αβ(descendant)"]
LSHORT = {"ancestral(no_somatic)": "ancestral", "γ_subclone(sibling)": "γ subclone", "L1_α_only(ancestor)": "L1 α-only", "L2_αβ(descendant)": "L2 α+β"}
LCOL = {"ancestral(no_somatic)": MUT, "γ_subclone(sibling)": PUR, "L1_α_only(ancestor)": BLU, "L2_αβ(descendant)": AC, "complex": "#ccc"}
plt.rcParams.update({"font.size": 10, "text.color": TX, "axes.edgecolor": BD, "figure.facecolor": "white", "axes.facecolor": "white"})


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=115, bbox_inches="tight"); plt.close(fig); return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def fig_tree():
    fig, ax = plt.subplots(figsize=(8.5, 4.2)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 8)
    N = {"anc": (6, 7, MUT, f"germline hap1\nancestral\n{lc.get('ancestral(no_somatic)',0)} reads"),
         "γ": (2.5, 4, PUR, f"γ subclone\n+48357368\n{lc.get('γ_subclone(sibling)',0)} reads"),
         "α": (8.5, 4.5, BLU, f"α subclone\n+48365089\n{lc.get('L1_α_only(ancestor)',0)+lc.get('L2_αβ(descendant)',0)} reads"),
         "L1": (7, 1.5, BLU, f"L1 α-only\n{lc.get('L1_α_only(ancestor)',0)} reads"),
         "L2": (10.5, 1.5, AC, f"L2 α+β\n+48365161,48362515\n{lc.get('L2_αβ(descendant)',0)} reads")}
    edges = [("anc", "γ"), ("anc", "α"), ("α", "L1"), ("α", "L2")]
    for a, b in edges:
        ax.annotate("", xy=(N[b][0], N[b][1] + 0.6), xytext=(N[a][0], N[a][1] - 0.6), arrowprops=dict(arrowstyle="-|>", color=TX, lw=1.5))
    for k, (x, y, col, lab) in N.items():
        ax.add_patch(plt.Circle((x, y), 0.92, color=col, alpha=0.85)); ax.text(x, y, lab, ha="center", va="center", fontsize=7, color="white")
    ax.text(3.7, 5.8, "sibling", fontsize=8, color=PUR); ax.text(7.6, 6, "sibling", fontsize=8, color=BLU)
    ax.text(10.7, 3.1, "nested", fontsize=8, color=AC)
    ax.set_title("完整局部克隆樹（4 somatic SNV，read 連鎖重建，全 normal=REF）", fontsize=11, color=TX)
    return b64(fig)


def svg_linkage():
    grp = {g: [r for r in reads if r["lineage"] == g] for g in LABS}
    rh = 9; gap = 16; x0 = 88; cw = 40
    cols = [("γ", "γ"), ("α", "α"), ("β1", "β1"), ("β2", "β2")]
    H = sum(len(grp[g]) * rh + gap for g in LABS if grp[g]) + 42
    s = [f'<svg viewBox="0 0 560 {H}" xmlns="http://www.w3.org/2000/svg" font-family="system-ui" font-size="10">']
    for k, (nm, _) in enumerate(cols):
        s.append(f'<text x="{x0+k*cw}" y="14" font-size="10">{nm}</text>')
    s.append(f'<text x="{x0+4*cw+10}" y="14" font-size="10">methyl mean</text>')
    y = 24
    for g in LABS:
        if not grp[g]:
            continue
        s.append(f'<text x="6" y="{y+11}" fill="{LCOL[g]}" font-weight="bold" font-size="10">{LSHORT[g]} ({len(grp[g])})</text>')
        for r in sorted(grp[g], key=lambda r: np.nanmean([np.nan if v is None else v for v in r["meth"]]) if any(v is not None for v in r["meth"]) else 0):
            for k, (nm, _) in enumerate(cols):
                al = r["geno"].get(nm); col = "#cfcabf" if al not in ("REF", "ALT") else (BLU if al == "REF" else RED)
                s.append(f'<rect x="{x0+k*cw}" y="{y}" width="{cw-3}" height="{rh-1.5}" fill="{col}"/>')
            mv = np.nanmean([np.nan if v is None else v for v in r["meth"]]); mv = 0 if np.isnan(mv) else mv
            s.append(f'<rect x="{x0+4*cw+10}" y="{y}" width="{mv*110:.1f}" height="{rh-1.5}" fill="#7a9a7a"/>')
            y += rh
        y += gap
    s.append(f'<rect x="{x0}" y="{H-15}" width="10" height="10" fill="{BLU}"/><text x="{x0+13}" y="{H-6}">REF</text><rect x="{x0+50}" y="{H-15}" width="10" height="10" fill="{RED}"/><text x="{x0+63}" y="{H-6}">ALT</text>')
    s.append("</svg>")
    return "".join(s)


def fig_methheat():
    ordered = []
    for g in LABS:
        grp = [r for r in reads if r["lineage"] == g]
        grp.sort(key=lambda r: np.nanmean([np.nan if v is None else v for v in r["meth"]]) if any(v is not None for v in r["meth"]) else 0)
        ordered += [(g, r) for r in grp]
    mat = np.array([[np.nan if v is None else v for v in r["meth"]] for _, r in ordered])
    fig, ax = plt.subplots(figsize=(11, 4.6)); cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#e8e4dc")
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    y = 0
    for g in LABS:
        c = sum(1 for l, _ in ordered if l == g)
        if c == 0:
            continue
        if y > 0:
            ax.axhline(y - 0.5, color=TX, lw=1.2)
        ax.text(-2.5, y + c / 2, LSHORT[g], va="center", ha="right", fontsize=9.5, color=LCOL[g], fontweight="bold"); y += c
    ax.set_xlabel(f"CpG sites (n={len(cpgs)})"); ax.set_yticks([]); fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01).set_label("β")
    ax.set_title("per-read × CpG 甲基熱圖（按完整 4-lineage 分組）", fontsize=10, color=TX)
    return b64(fig)


charts = {"tree": fig_tree(), "methheat": fig_methheat()}
svg = svg_linkage()
pairs = D["pairs_2x2"]
def t2(pk, pv):
    a, b = pk.split("×")
    excl = "互斥(sibling)" if pv["AA"] == 0 and (pv["RA"] + pv["AR"]) >= 2 else ("nested" if (pv["RA"] == 0) != (pv["AR"] == 0) else ("co-link" if pv["RA"] == 0 and pv["AR"] == 0 else "—"))
    return (f'<table class="mini"><caption>{a}×{b} <b>{excl}</b></caption><tr><th></th><th>{b}R</th><th>{b}A</th></tr>'
            f'<tr><th>{a}R</th><td>{pv["RR"]}</td><td>{pv["RA"]}</td></tr><tr><th>{a}A</th><td>{pv["AR"]}</td><td>{pv["AA"]}</td></tr></table>')
pair_html = "".join(t2(k, v) for k, v in pairs.items())
snv_rows = "".join(f"<tr><td>{nm}</td><td>{v['pos']}</td><td>{v['ref']}→{v['alt']}</td><td>{v['tumor_REF']}/{v['tumor_ALT']}</td><td>{v['normal_REF']}/{v['normal_ALT']}</td><td>{'✓' if v['somatic'] else '⚠'}</td></tr>" for nm, v in D["snv_norm"].items())
html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>chr17:48360161 完整 subclone 例</title><style>
:root{{--ac:{AC};--tx:{TX};--bd:{BD};--mut:{MUT};--grn:{GRN};--red:{RED}}}
*{{box-sizing:border-box}} body{{margin:0;font-family:system-ui,"Noto Sans CJK TC","PingFang TC",sans-serif;color:var(--tx);background:#FAF9F5;line-height:1.62}}
.wrap{{max-width:1100px;margin:0 auto;padding:0 22px 80px}}
header{{background:var(--tx);color:#FAF9F5;padding:28px 22px;margin-bottom:8px}} header h1{{margin:0 0 5px;font-size:22px}} header .sub{{color:#cfc9bf;font-size:13.5px}}
nav{{position:sticky;top:0;background:#FAF9F5;border-bottom:1px solid var(--bd);padding:9px 22px;font-size:12.5px;z-index:9;margin-bottom:14px}} nav a{{color:var(--mut);text-decoration:none;margin-right:12px}}
h2{{font-size:17px;margin:26px 0 9px;padding-bottom:5px;border-bottom:2px solid var(--ac)}}
img,svg{{max-width:100%;border:1px solid var(--bd);border-radius:6px;margin:8px 0;background:white}} svg{{padding:8px}}
table{{border-collapse:collapse;font-size:12.5px;margin:8px 0}} th,td{{border:1px solid var(--bd);padding:4px 8px;text-align:center}} th{{background:#efe9df}} caption{{font-size:11.5px;color:var(--mut);padding:2px}}
table.mini{{display:inline-table;margin:5px 10px 5px 0}}
.kpi{{display:flex;gap:10px;flex-wrap:wrap;margin:12px 0}} .kpi div{{background:white;border:1px solid var(--bd);border-radius:8px;padding:10px 14px;min-width:115px}} .kpi b{{display:block;font-size:20px;color:var(--ac)}} .kpi span{{font-size:11.5px;color:var(--mut)}}
.ok{{border-left:4px solid var(--grn);background:#eef4ee;padding:11px 14px;border-radius:0 6px 6px 0;margin:12px 0}}
.red{{border-left:4px solid var(--red);background:#fbf0ec;padding:11px 14px;border-radius:0 6px 6px 0;margin:12px 0}}
details{{margin:10px 0}} summary{{cursor:pointer;color:var(--ac);font-weight:bold}}
footer{{margin-top:32px;padding-top:13px;border-top:1px solid var(--bd);font-size:12px;color:var(--mut)}}
.tag{{display:inline-block;background:var(--ac);color:white;font-size:11px;padding:1px 7px;border-radius:10px;margin-left:6px}}
</style></head><body>
<header><div class="wrap" style="padding:0"><h1>chr17:48360161 — 完整 subclone 例（4-sSNV 克隆樹）<span class="tag">⭐3 HCC1395</span></h1>
<div class="sub">2 sibling 分支(γ,α) + α 內 nested 後代。含 SEQC2 未收 γ（ONT 候選）的方法學洞見。</div></div></header>
<nav class="wrap" style="max-width:1100px"><a href="#tldr">結論</a><a href="#tree">完整克隆樹</a><a href="#link">read 關聯</a><a href="#pairs">6對2×2</a><a href="#meth">甲基</a><a href="#insight">方法學洞見</a></nav>
<div class="wrap">

<section id="tldr"><h2>結論</h2>
<div class="kpi">
<div><b>4</b><span>位點（α/β SEQC2＋γ ONT候選；全 normal=REF）</span></div>
<div><b>4</b><span>細胞群（ancestral/γ/L1/L2）</span></div>
<div><b>2+1</b><span>sibling 分支 + nested 後代</span></div>
</div>
<div class="ok"><b>✅ 完整局部克隆樹</b>：germline hap1 → <b>γ subclone(5)</b> 與 <b>α subclone(39)</b> 兩 sibling 分支（γ×α 互斥 AA=0）；α 內再分 <b>L1 α-only(20)</b> 與 <b>L2 α+β(19)</b> nested 後代（β 嵌套在 α）。全 4 位點 normal=REF（α/β1/β2 SEQC2 HighConf 確認；γ 為 ONT 候選，SEQC2 未收）。這是 read-level 連鎖直接重建的乾淨局部 subclonal phylogeny。</div></section>

<section id="tree"><h2>① 完整克隆樹</h2><img src="{charts['tree']}" alt="tree">
<table><tr><th>somatic</th><th>pos</th><th>ref→alt</th><th>tumor R/A</th><th>normal R/A</th><th>somatic</th></tr>{snv_rows}</table></section>

<section id="link"><h2>② read 關聯（SVG，4-sSNV 基因型）</h2>
<p>每列一條 read（按 lineage 分組），4 格=γ/α/β1/β2（藍REF 紅ALT）；右綠=甲基均值。<b>γ subclone 只 γ-ALT；α 分支 α-ALT；L2 再加 β</b>。</p>{svg}</section>

<section id="pairs"><h2>③ 6 對 pairwise 共現 2×2</h2>{pair_html}
<p>🔑 γ×α/γ×β1/γ×β2 全 <b>AA=0</b>（γ 與 α 分支互斥=sibling）；α×β2/β1×α <b>AR=0</b>（β 嵌套 α）；β1×β2 只 RR+AA（共連）。</p></section>

<section id="meth"><h2>④ per-lineage 甲基熱圖</h2><img src="{charts['methheat']}" alt="meth">
<p>reads 按完整 4-lineage 分組；甲基模式可 characterize 各遺傳定義的克隆群（非循環）。</p></section>

<section id="insight"><h2>⑤ 🔴 方法學洞見：SEQC2 未收的 clonally-informative ONT 候選</h2>
<div class="red">γ=<b>48357368</b> <b>不在 SEQC2 HighConf/superSet</b>（落在 callable HC region 內＝缺金標準確認，非主動判 FP）→ 落到我們 pipeline 的 filtered-FP 集；最初只用 filtered TP 集時<b>漏掉它</b>，把 γ subclone 誤標成「ancestral root」。但 γ <b>normal 35/0 REF + tumor VAF 0.18 = ONT 偵測候選</b>（somatic-like，惟 SEQC2 未確認），且 longphase 用它 phase（故那些 read 被標 HP1-1）。<br><br>
<b>教訓</b>：① <b>「HP1-1 但所有已知 sSNV 皆 REF」是發現漏掉 somatic 的訊號</b>（longphase HP tag 比 filtered VCF 知道更多）；② 完整克隆重建應用 <b>longphase 輸入集（methyl_PASS）+ normal 確認</b>，非只 SEQC2-filtered TP；③ 這也解釋為何 phase-block（PS=46496608，~1.9Mb）能標 tag 但 read 只連局部——HP tag 是 Mb 級、read 共現是 ≤read 長級。</div></section>

<footer>chr17:48360161 · HCC1395 ⭐3 · §13-A 由 chr17_complete_data.json 注入 · 腳本 p3_chr17_complete.py + build_chr17_complete.py<br>
γ=48357368(longphase methyl_PASS, SEQC2 未收/HC region 內 ONT 候選, normal=REF) · α=48365089 · β1=48362515 · β2=48365161(SEQC2 HighConf)</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes")
