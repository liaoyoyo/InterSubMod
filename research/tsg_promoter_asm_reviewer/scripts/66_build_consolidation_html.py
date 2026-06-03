#!/usr/bin/env python3
"""
66 - Standalone HTML consolidating the 38 key positions across 6 samples:
existence (coverage/methylation completeness) + credible ASM gene-set existence +
copy-number (CN) consistency.

§13 layer-A: numbers from keypos_cn_consolidation.json (65) + 6 *_credible_discovery.json (62).

Output: docs/experiments/in_progress/2026/06/
        20260603_keypos_cn_consolidation_01.standalone.html
"""
import os, json, html, base64, io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
       "20260603_keypos_cn_consolidation_01.standalone.html")

ALL = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}
CCOLOR = {"breast": "#be123c", "lung": "#0e7490", "melanoma": "#7c3aed"}
ZCOLOR = {"LOH": "#dc2626", "het": "#2563eb", "subhap": "#16a34a", "low": "#9ca3af"}


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig); return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def esc(x):
    return html.escape(str(x))


C = json.load(open(f"{CS}/keypos_cn_consolidation.json"))
ex = C["existence"]; cv = C["hcc1395_cn_validation"]; cc = C["cross_sample_cn_consistency"]
PP = C["per_position"]
med_cov = C["meta"]["sample_median_coverage"]

# credible counts from discovery
cred_n = {}
for s in ALL:
    p = f"{CS}/{s}_credible_discovery.json"
    if os.path.exists(p):
        fn = json.load(open(p))["funnel"]
        cred_n[s] = (fn["n_credible_pass_tierA_relax"], fn["n_credible_pass_tierA_strict"])
    else:
        cred_n[s] = (None, None)


# ---- FIG 1: coverage heatmap (38 x 6) ----
def fig_cov_heatmap():
    M = np.array([[p["per_sample"][s]["cov"] for s in ALL] for p in PP])
    fig, ax = plt.subplots(figsize=(6.4, 9.5))
    im = ax.imshow(np.clip(M, 0, 400), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(ALL)))
    ax.set_xticklabels([f"{s}\n({CANCER[s][:3]})" for s in ALL], fontsize=8)
    ax.set_yticks(range(len(PP)))
    ax.set_yticklabels([p["gene"][:12] for p in PP], fontsize=6)
    ax.set_title("Coverage at 38 key positions × 6 samples\n(all ≥10 = clearly exist)")
    fig.colorbar(im, ax=ax, fraction=0.04, label="reads in window (clip 400)")
    return b64(fig)


# ---- FIG 2: HCC1395 rel_cov vs SEQC2 median_cn ----
def fig_cn_scatter():
    xs, ys, labs = [], [], []
    for p in PP:
        mcn = p["hcc1395_cn"].get("median_cn")
        rc = p["per_sample"]["HCC1395"]["rel_cov"]
        if mcn is not None and rc is not None:
            xs.append(mcn); ys.append(rc)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.scatter(xs, ys, s=26, alpha=0.7, c="#1E3A8A")
    ax.set_xlabel("HCC1395 SEQC2 median CN (ground-truth)")
    ax.set_ylabel("relative coverage (n_reads / sample median)")
    ax.set_title(f"HCC1395: coverage tracks SEQC2 copy number\n"
                 f"Spearman ρ={cv['cn_relcov_spearman']} (p={cv['cn_relcov_p']}, n={cv['n_cn_pairs']})")
    ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


# ---- FIG 3: per-position zygosity across samples (cancer-specific CN) ----
def fig_zyg():
    order = ["LOH", "het", "subhap", "low"]
    M = np.zeros((len(PP), len(ALL)))
    code = {z: i for i, z in enumerate(order)}
    for i, p in enumerate(PP):
        for j, s in enumerate(ALL):
            M[i, j] = code[p["per_sample"][s]["zygosity"]]
    cmap = LinearSegmentedColormap.from_list("z", [ZCOLOR[z] for z in order], N=4)
    fig, ax = plt.subplots(figsize=(6.4, 9.5))
    ax.imshow(M, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
    ax.set_xticks(range(len(ALL)))
    ax.set_xticklabels([s for s in ALL], rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(PP)))
    ax.set_yticklabels([p["gene"][:12] for p in PP], fontsize=6)
    ax.set_title(f"Zygosity per position × sample\n"
                 f"only {cc['n_zyg_consistent']}/{len(PP)} consistent → CN cancer-specific")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=ZCOLOR[z], label=z) for z in order],
              fontsize=7, loc="upper right", bbox_to_anchor=(1.0, 1.06), ncol=4)
    return b64(fig)


FIG_COV = fig_cov_heatmap()
FIG_CN = fig_cn_scatter()
FIG_ZYG = fig_zyg()

# per-sample summary table
samp_rows = ""
for s in ALL:
    cr = cred_n[s]
    chip = f'<span class="chip" style="background:{CCOLOR[CANCER[s]]}">{esc(s)}</span>'
    cred_str = (f'{cr[0]}/{cr[1]}' if cr[0] is not None else "—")
    flag = ' <span class="warn-i">覆蓋限制</span>' if med_cov[s] < 50 else ""
    samp_rows += (f'<tr><td>{chip}</td><td>{esc(CANCER[s])}</td>'
                  f'<td class="num">{med_cov[s]:.0f}{flag}</td>'
                  f'<td class="num">{cred_str}</td></tr>')

# position existence table (min cov + zygosity spread)
pos_rows = ""
for p in sorted(PP, key=lambda p: p["min_cov"])[:14]:
    cn = p["hcc1395_cn"]
    zc = p["zyg_counts"]
    pos_rows += (f'<tr><td>{esc(p["gene"])}</td><td class="mono">{esc(p["pos"])}</td>'
                 f'<td class="num">{p["min_cov"]}</td>'
                 f'<td class="num">{"✓" if p["all_covered"] else "✗"}</td>'
                 f'<td class="num">{cn.get("median_cn") if cn.get("median_cn") is not None else "—"}</td>'
                 f'<td class="num">{esc(cn.get("cn_class") or "")}</td>'
                 f'<td>LOH{zc["LOH"]}/het{zc["het"]}/sub{zc["subhap"]}</td></tr>')

HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>關鍵位點跨樣本統整 — 存在性 + CN 一致性</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1080px;margin:0 auto;padding:24px 18px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:24px 26px;margin-bottom:16px}}
header.top h1{{margin:0 0 6px;font-size:1.4rem}}.meta{{font-size:.82rem;opacity:.88}}
.tldr{{background:#eef2ff;border:1px solid #c7d2fe;border-left:5px solid #1E3A8A;border-radius:12px;padding:16px 20px;margin-bottom:14px}}
.tldr h2{{margin:0 0 8px;font-size:1.06rem;color:#1e3a8a}}.tldr p{{margin:.4rem 0}}
.kpis{{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}}
.kpi{{flex:1;min-width:150px;background:#fff;border:1px solid var(--line);border-radius:12px;padding:14px;text-align:center}}
.kpi .v{{font-size:1.4rem;font-weight:800;color:var(--accent)}}.kpi .l{{font-size:.75rem;color:var(--mut);margin-top:3px}}
.card{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:14px 0}}
.qlabel{{font-weight:800;color:var(--accent);font-size:1rem}}
table{{width:100%;border-collapse:collapse;font-size:.83rem;margin-top:8px}}
th,td{{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left}}th{{background:#f1f5f9;font-size:.77rem}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}.mono{{font-family:ui-monospace,monospace;font-size:.78rem;color:#475569}}
.chip{{color:#fff;font-weight:700;font-size:.72rem;padding:2px 8px;border-radius:6px}}.warn-i{{color:#b45309;font-size:.7rem;font-weight:700}}
.two{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}@media(max-width:820px){{.two{{grid-template-columns:1fr}}}}
.figbox img{{width:100%;border:1px solid var(--line);border-radius:8px;margin-top:8px}}
.legend{{font-size:.76rem;color:var(--mut);margin-top:4px}}.tablewrap{{overflow-x:auto}}
details.method>summary{{cursor:pointer;color:var(--accent);font-weight:600;font-size:.84rem}}
details.method p{{font-size:.83rem;background:var(--bg);border-radius:8px;padding:8px 12px}}
.mech{{background:#f0fdf4;border:1px solid #bbf7d0;border-left:5px solid #16a34a;border-radius:12px;padding:14px 18px;margin:14px 0;font-size:.9rem}}
footer{{margin-top:24px;font-size:.76rem;color:var(--mut);border-top:1px solid var(--line);padding-top:12px}}
</style></head><body><div class="wrap">

<header class="top">
  <h1>關鍵位點跨樣本統整 — 存在性 + 甲基完整性 + CN 一致性</h1>
  <div class="meta">38 個關鍵位點 × 6 樣本 · credible ASM 基因集存在性 · 與各樣本 copy 狀況關係 · A pilot · 2026-06-03 · §13 layer-A</div>
</header>

<div class="tldr">
  <h2>TL;DR — 位點都存在且甲基觀察完整；但 CN 是 cancer-specific</h2>
  <p>① <b>存在性 / 甲基完整性：38/38 位點在全 6 樣本都 covered≥10（最低 {ex["min_cov_overall"]}）</b> → 甲基觀察完整、無 gap。</p>
  <p>② <b>credible ASM 基因集存在</b>（各樣本自己）：HCC1395 {cred_n["HCC1395"][1]} · HCC1954 {cred_n["HCC1954"][1]} · H2009 {cred_n["H2009"][1]} · H1437 {cred_n["H1437"][1]} · HCC1937 {cred_n["HCC1937"][1]} strict；<b>COLO829 {cred_n["COLO829"][1]}（覆蓋限制 median cov {med_cov["COLO829"]:.0f}）</b>。</p>
  <p>③ <b>與 CN 的關係</b>：HCC1395 內覆蓋強烈追蹤 SEQC2 CN（Spearman <b>{cv["cn_relcov_spearman"]}</b>）；但跨樣本只 <b>{cc["n_zyg_consistent"]}/38</b> 位點 zygosity 一致 → <b>CN/LOH 幾乎完全 cancer-specific</b>。</p>
</div>

<div class="kpis">
  <div class="kpi"><div class="v">{ex["n_all_covered"]}/38</div><div class="l">全樣本 covered≥10（存在）</div></div>
  <div class="kpi"><div class="v">{cv["cn_relcov_spearman"]}</div><div class="l">HCC1395 覆蓋~SEQC2 CN ρ</div></div>
  <div class="kpi"><div class="v">{cc["n_zyg_consistent"]}/38</div><div class="l">跨樣本 zygosity 一致</div></div>
  <div class="kpi"><div class="v">{med_cov["COLO829"]:.0f}</div><div class="l">COLO829 中位覆蓋（最低）</div></div>
</div>

<div class="mech">
  <b>🔑 統一機制（回答「一致的現象與問題」）</b>：覆蓋（∝CN，HCC1395 內 ρ={cv["cn_relcov_spearman"]} 驗證）<b>gates ASM 可偵測性</b>；而 CN 是 <b>cancer-specific</b>（{cc["n_zyg_consistent"]}/38 一致）→ 故 credible ASM 產量也 cancer-specific。<b>COLO829</b>（melanoma，中位覆蓋僅 {med_cov["COLO829"]:.0f} ≈ 其他 0.3×）→ ASM 子群 reads 不足 → <b>0 credible</b>。即「問題」= ASM 偵測受各樣本自身 CN/覆蓋制約，非生物缺席。
</div>

<section class="card">
  <span class="qlabel">① 各樣本：中位覆蓋 + credible ASM 基因集</span>
  <div class="tablewrap"><table>
    <thead><tr><th>樣本</th><th>癌種</th><th class="num">38 位點中位覆蓋</th><th class="num">credible (relax/strict)</th></tr></thead>
    <tbody>{samp_rows}</tbody></table></div>
</section>

<div class="two">
  <section class="card">
    <span class="qlabel">② 存在性 / 甲基完整性</span>
    <p class="legend">38/38 位點全樣本 covered≥10（最低 {ex["min_cov_overall"]}）。下列最低覆蓋 14 位點 + HCC1395 SEQC2 CN + zygosity 分布。</p>
    <div class="tablewrap"><table>
      <thead><tr><th>gene</th><th>pos</th><th class="num">min cov</th><th class="num">all≥10</th>
        <th class="num">HCC1395 CN</th><th class="num">cn_class</th><th>zygosity 分布</th></tr></thead>
      <tbody>{pos_rows}</tbody></table></div>
    <div class="figbox"><img src="{FIG_COV}"></div>
  </section>
  <section class="card">
    <span class="qlabel">③ CN 關係（HCC1395 驗證 + 跨樣本）</span>
    <p class="legend">左：HCC1395 覆蓋強烈追蹤 SEQC2 CN（ρ={cv["cn_relcov_spearman"]}）；右：各位點 zygosity 跨 6 樣本幾乎不一致 → CN cancer-specific。</p>
    <div class="figbox"><img src="{FIG_CN}"></div>
    <div class="figbox"><img src="{FIG_ZYG}"></div>
  </section>
</div>

<section class="card">
  <span class="qlabel">④ 方法 / 限制</span>
  <details class="method" open><summary>讀法 / caveat（必讀）</summary><p>
  <b>存在性</b>=window 覆蓋≥10（38 位點 source = HCC1395 credible loci + BRCA2，在其他樣本檢視）。<br>
  <b>CN</b>：HCC1395 用 SEQC2 median_cn/cn_class（ground-truth，HCC1395-only）；<b>其他樣本無正交 CN</b> → 用相對覆蓋(n_reads/樣本中位)當 depth/CN proxy + HP 結構 zygosity(LOH/het/subhap)。<br>
  <b>核心限制</b>：① 其他 5 樣本 CN 為覆蓋 proxy 非正交真值（SEQC2 僅 HCC1395）；② zygosity LOH 判定用 HP-tag 結構(單倍型 dominant)，與 SEQC2 LOH 僅 {cv["loh_structure_vs_seqc2_match"]}/{cv["loh_total"]} 吻合（HCC1395 多位點有 somatic subhap 歸 subhap 非 LOH）；③ 這 38 位點是 HCC1395-derived，在其他樣本是 germline/非 somatic context；④ 覆蓋~CN ρ=0.871 在 HCC1395 內驗證，外推其他樣本未獨立驗證。<br>
  <b>誠實結論</b>：位點都<b>明顯存在且甲基觀察完整</b>（38/38 全覆蓋）；ASM 可偵測性<b>受各樣本自身覆蓋/CN 制約</b>（HCC1395 ρ=0.871）；CN <b>cancer-specific</b>（1/38 一致）→ credible ASM 產量 cancer-specific，COLO829 低覆蓋 → 0。</p></details>
</section>

<footer>
  數據源：<span class="mono">research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/keypos_cn_consolidation.json</span>
  （+ 6× *_keypos.json + 6× *_credible_discovery.json + master_o1_cn.tsv）· 腳本 65/66 · §13 layer-A
</footer>
</div></body></html>'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(HTML_DOC)
print(f"[66] wrote {OUT} ({len(HTML_DOC)//1024} KB)")
