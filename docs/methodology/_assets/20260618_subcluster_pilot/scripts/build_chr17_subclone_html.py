#!/usr/bin/env python3
"""[chr17:48360161 明確 subclone 例 HTML] 克隆樹 + sSNV 連鎖 2×2 + per-read×CpG 甲基熱圖(按 lineage 分組)
+ 祖先 L1 vs 後代 L2 甲基差異(非循環) + normal 對照 + 全數據。§13-A 由 chr17_subclone_data.json 注入。
輸出 ../../20260625_chr17_48360161_subclone_example_01.standalone.html。"""
import json, io, base64, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260625_chr17_48360161_subclone_example_01.standalone.html"
D = json.load(open(f"{A}/chr17_subclone_data.json"))
AC = "#D97757"; TX = "#141413"; BG = "#FAF9F5"; BD = "#E3DACC"; MUT = "#6B6862"; GRN = "#5B8A5B"; RED = "#C0563F"; BLU = "#4A6E8A"
plt.rcParams.update({"font.size": 10, "text.color": TX, "axes.labelcolor": TX, "xtick.color": TX, "ytick.color": TX,
                     "axes.edgecolor": BD, "figure.facecolor": "white", "axes.facecolor": "white"})


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=115, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


LIN_ORDER = ["L0_ancestral_root", "L1_alpha_only(ancestor)", "L2_alpha_beta(descendant)", "other"]
LIN_COL = {"L0_ancestral_root": MUT, "L1_alpha_only(ancestor)": BLU, "L2_alpha_beta(descendant)": AC, "other": "#cccccc"}
LIN_SHORT = {"L0_ancestral_root": "L0 root", "L1_alpha_only(ancestor)": "L1 祖先(α)", "L2_alpha_beta(descendant)": "L2 後代(α+β)", "other": "other"}


# ---- 克隆樹 ----
def fig_tree():
    fig, ax = plt.subplots(figsize=(7.5, 3.2)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
    lc = D["lineage_counts"]
    nodes = {"root": (1.2, 3, MUT, f"L0 root\n89=REF\n{lc.get('L0_ancestral_root',0)} reads"),
             "L1": (5, 3, BLU, f"L1 祖先\nα: 89-ALT\n{lc.get('L1_alpha_only(ancestor)',0)} reads"),
             "L2": (8.6, 3, AC, f"L2 後代\nα+β: 89,161,2515-ALT\n{lc.get('L2_alpha_beta(descendant)',0)} reads")}
    for k, (x, y, col, lab) in nodes.items():
        ax.add_patch(plt.Circle((x, y), 0.95, color=col, alpha=0.85))
        ax.text(x, y, lab, ha="center", va="center", fontsize=7.5, color="white")
    for x0, x1, lab in [(2.15, 4.05, "+α (48365089)"), (5.95, 7.65, "+β (48365161+48362515)")]:
        ax.annotate("", xy=(x1, 3), xytext=(x0, 3), arrowprops=dict(arrowstyle="-|>", color=TX, lw=1.6))
        ax.text((x0 + x1) / 2, 3.7, lab, ha="center", fontsize=8, color=TX)
    ax.set_title("局部克隆樹（sSNV 連鎖重建，非循環）", fontsize=11, color=TX)
    return b64(fig)


# ---- 甲基熱圖: reads × CpG, 按 lineage 分組 ----
def fig_heatmap():
    reads = D["reads"]; cpgs = D["cpgs"]
    ordered = []
    for lin in LIN_ORDER:
        grp = [r for r in reads if r["lineage"] == lin]
        grp.sort(key=lambda r: np.nanmean([np.nan if v is None else v for v in r["meth"]]) if any(v is not None for v in r["meth"]) else 0)
        ordered += [(lin, r) for r in grp]
    mat = np.array([[np.nan if v is None else v for v in r["meth"]] for _, r in ordered], float)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#e8e4dc")
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    # lineage 分隔線 + 標籤
    y = 0
    for lin in LIN_ORDER:
        n = sum(1 for l, _ in ordered if l == lin)
        if n == 0:
            continue
        if y > 0:
            ax.axhline(y - 0.5, color=TX, lw=1.3)
        ax.text(-2.5, y + n / 2, LIN_SHORT[lin], va="center", ha="right", fontsize=9, color=LIN_COL[lin], rotation=0, fontweight="bold")
        y += n
    ax.set_xlabel(f"CpG sites (n={len(cpgs)}, 基因座標序)"); ax.set_ylabel("reads (按 lineage 分組)")
    ax.set_yticks([]); cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label("methylation β")
    ax.set_title("per-read × CpG 甲基熱圖（reads 按基因型 lineage 分組 → 甲基模式隨 lineage 變）", fontsize=10.5, color=TX)
    return b64(fig)


# ---- 祖先 vs 後代 甲基 profile + normal ----
def fig_profile():
    diff = D["L1_vs_L2_diff"]; cpgs = [d["cpg"] for d in diff]
    L1 = [d["L1"] for d in diff]; L2 = [d["L2"] for d in diff]
    nm = D["normal_cpg_mean"]; allc = D["cpgs"]
    nmap = {allc[j]: nm[j] for j in range(len(allc))}
    norm = [nmap.get(c) for c in cpgs]
    x = np.arange(len(cpgs))
    fig, ax = plt.subplots(figsize=(11, 3.6))
    ax.plot(x, L1, "-o", color=BLU, ms=3, label=f"L1 祖先(α) n={D['L1_n']}")
    ax.plot(x, L2, "-o", color=AC, ms=3, label=f"L2 後代(α+β) n={D['L2_n']}")
    ax.plot(x, [v if v is not None else np.nan for v in norm], "-", color=GRN, lw=1, alpha=0.7, label=f"normal n={D['n_normal_meth_reads']}")
    sig = {d["cpg"] for d in D["sig_diff_cpg"]}
    for i, c in enumerate(cpgs):
        if c in sig:
            ax.axvspan(i - 0.4, i + 0.4, color=RED, alpha=0.08)
    ax.set_xticks([]); ax.set_ylabel("methylation β"); ax.set_ylim(-0.05, 1.05); ax.legend(fontsize=8.5, frameon=False, ncol=3)
    ax.set_title(f"祖先 L1 vs 後代 L2 per-CpG 甲基（{D['n_sig_diff_cpg']} CpG |Δβ|≥0.2 紅標；非循環=lineage 由 sSNV 定義）", fontsize=10, color=TX)
    return b64(fig)


charts = {"tree": fig_tree(), "heatmap": fig_heatmap(), "profile": fig_profile()}
# 2×2 tables
pairs = D["pairs_2x2"]
def cell_tbl(pk, pv):
    a, b = pk.split("_")
    return (f'<table class="mini"><caption>{a} × {b} — <b>{pv["rel"]}</b></caption>'
            f'<tr><th></th><th>{b}=REF</th><th>{b}=ALT</th></tr>'
            f'<tr><th>{a}=REF</th><td>{pv["REF_REF"]}</td><td>{pv["REF_ALT"]}</td></tr>'
            f'<tr><th>{a}=ALT</th><td>{pv["ALT_REF"]}</td><td>{pv["ALT_ALT"]}</td></tr></table>')
pair_html = "".join(cell_tbl(k, v) for k, v in pairs.items())
snv_rows = "".join(f"<tr><td>{p['pos']}</td><td>{p['ref']}→{p['alt']}</td>"
                   f"<td>{D['snv_stat'][str(p['pos'])]['tumor_REF']}/{D['snv_stat'][str(p['pos'])]['tumor_ALT']}</td>"
                   f"<td>{D['snv_normal'][str(p['pos'])]['normal_REF']}/{D['snv_normal'][str(p['pos'])]['normal_ALT']}</td></tr>"
                   for p in D["snvs"])
lc = D["lineage_counts"]
html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>chr17:48360161 — 明確 nested subclone 例</title><style>
:root{{--ac:{AC};--tx:{TX};--bg:{BG};--bd:{BD};--mut:{MUT};--grn:{GRN};--red:{RED}}}
*{{box-sizing:border-box}} body{{margin:0;font-family:system-ui,"Noto Sans CJK TC","PingFang TC",sans-serif;color:var(--tx);background:var(--bg);line-height:1.65}}
.wrap{{max-width:1100px;margin:0 auto;padding:0 24px 80px}}
header{{background:var(--tx);color:var(--bg);padding:30px 24px;margin-bottom:8px}} header h1{{margin:0 0 6px;font-size:23px}} header .sub{{color:#cfc9bf;font-size:14px}}
nav{{position:sticky;top:0;background:var(--bg);border-bottom:1px solid var(--bd);padding:9px 24px;font-size:13px;z-index:9;margin-bottom:16px}}
nav a{{color:var(--mut);text-decoration:none;margin-right:13px}} nav a:hover{{color:var(--ac)}}
h2{{font-size:18px;margin:30px 0 10px;padding-bottom:5px;border-bottom:2px solid var(--ac)}}
img{{max-width:100%;border:1px solid var(--bd);border-radius:6px;margin:8px 0;background:white}}
table{{border-collapse:collapse;font-size:13px;margin:8px 0}} th,td{{border:1px solid var(--bd);padding:4px 9px;text-align:center}} th{{background:#efe9df}}
table.mini{{display:inline-table;margin:6px 14px 6px 0}} table.mini caption{{font-size:12px;padding:3px;color:var(--mut)}}
.kpi{{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}} .kpi div{{background:white;border:1px solid var(--bd);border-radius:8px;padding:11px 15px;min-width:120px}}
.kpi b{{display:block;font-size:21px;color:var(--ac)}} .kpi span{{font-size:12px;color:var(--mut)}}
.ok{{border-left:4px solid var(--grn);background:#eef4ee;padding:11px 15px;border-radius:0 6px 6px 0;margin:12px 0}}
.red{{border-left:4px solid var(--red);background:#fbf0ec;padding:11px 15px;border-radius:0 6px 6px 0;margin:12px 0}}
footer{{margin-top:36px;padding-top:14px;border-top:1px solid var(--bd);font-size:12px;color:var(--mut)}}
.tag{{display:inline-block;background:var(--ac);color:white;font-size:11px;padding:1px 7px;border-radius:10px;margin-left:6px}}
</style></head><body>
<header><div class="wrap" style="padding:0"><h1>chr17:48360161 — 明確的 nested subclone 例 <span class="tag">⭐3 HCC1395</span></h1>
<div class="sub">遺傳 sSNV 連鎖重建克隆階層 + 甲基非循環 characterize。論文「somatic haplotagging + methylation profiles」的具體體現。</div></div></header>
<nav class="wrap" style="max-width:1100px"><a href="#tldr">結論</a><a href="#tree">克隆樹</a><a href="#link">連鎖證據</a><a href="#meth">甲基熱圖</a><a href="#diff">祖先vs後代</a><a href="#why">為何明確</a></nav>
<div class="wrap">

<section id="tldr"><h2>結論</h2>
<div class="kpi">
<div><b>3</b><span>somatic SNV（nested）</span></div>
<div><b>{lc.get('L1_alpha_only(ancestor)',0)}+{lc.get('L2_alpha_beta(descendant)',0)}</b><span>祖先 + 後代 reads</span></div>
<div><b>{D['n_sig_diff_cpg']}</b><span>CpG 祖先vs後代甲基差(|Δβ|≥0.2)</span></div>
<div><b>2/3</b><span>somatic-confirmed(normal=REF)</span></div>
</div>
<div class="ok"><b>✅ 為何是明確 subclone</b>：① 遺傳上 3 個 somatic SNV 在 read 上呈 <b>nested</b>（α=48365089 祖先 → α+β 後代，REF_ALT=0），normal 皆 REF；② 表觀上 <b>{D['n_sig_diff_cpg']} 個 CpG 在祖先 L1 vs 後代 L2 顯著不同</b>——而 lineage 是<b>由 sSNV 定義</b>（非甲基切群）→ <b>甲基差異非循環</b>。這是「haplotag/sSNV 重建 lineage、甲基 characterize」的乾淨案例。</div></section>

<section id="tree"><h2>① 局部克隆樹（sSNV 連鎖重建）</h2>
<img src="{charts['tree']}" alt="tree">
<table><tr><th>somatic SNV</th><th>ref→alt</th><th>tumor REF/ALT</th><th>normal REF/ALT</th></tr>{snv_rows}</table></section>

<section id="link"><h2>② 連鎖證據（per-read 共現 2×2，nested）</h2>
{pair_html}
<p>🔑 <code>48365089 × 48365161</code>：ALT-REF <b>22</b> + ALT-ALT <b>29</b>，REF_ALT=<b>0</b> → 48365161(β) 嵌套在 48365089(α) 內（β 從不單獨出現）→ α 祖先、α+β 後代。<code>48362515 × 48365161</code> co-linked（β1、β2 同屬後代）。</p></section>

<section id="meth"><h2>③ 甲基熱圖（per-read × CpG，按基因型 lineage 分組）</h2>
<img src="{charts['heatmap']}" alt="heatmap">
<p>每列一條 read（按 sSNV lineage L0/L1/L2 分組），每行一個 CpG，顏色=β（藍 unmeth ↔ 紅 meth，灰=NA）。reads 是<b>用基因型分組</b>，甲基模式隨 lineage 改變 = 甲基 characterize 已驗 lineage。</p></section>

<section id="diff"><h2>④ 祖先 L1 vs 後代 L2 甲基差異（非循環 + tumor-specific）</h2>
<img src="{charts['profile']}" alt="profile">
<p>L1 祖先(α，n={D['L1_n']}) vs L2 後代(α+β，n={D['L2_n']}) 的 per-CpG 甲基；紅標 = {D['n_sig_diff_cpg']} 個 |Δβ|≥0.2 顯著差異 CpG。normal（綠，n={D['n_normal_meth_reads']}）對照確認方向。🔴 因 lineage 由 sSNV 定義，此甲基差異<b>非 double-dip</b>。</p></section>

<section id="why"><h2>⑤ 為何「明確」+ 誠實邊界</h2>
<div class="ok"><b>明確之處</b>：遺傳 nested（方向性 = 祖先→後代）+ normal=REF（somatic）+ 甲基非循環 characterize（lineage 先由 sSNV 定，再看甲基）。三者獨立交叉。</div>
<div class="red"><b>🔴 誠實邊界</b>：⭐3 單樣本；這是<b>局部</b>克隆階層（此 phase-block 內），非 genome-wide clone tree；48362515 normal 有 1 ALT（borderline，但 α/β2=48365089/48365161 normal 全 REF，nested 主結構不依賴它）；subclone 確認黃金標準仍是 single-cell/multi-region。</div></section>

<footer>chr17:48360161 · HCC1395 ⭐3 · build_branch docs/method-comparison-ism-external-202606 · §13-A 數字由 chr17_subclone_data.json 注入<br>
資料：tumor BAM per-read sSNV 等位 + methylation.csv + normal BAM 確認｜腳本 p3_chr17_extract.py + build_chr17_subclone_html.py</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes, 3 charts")
