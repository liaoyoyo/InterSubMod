#!/usr/bin/env python3
"""
64 - Cross-sample ASM landscape: compare each sample's credible ASM gene set.

§13 layer-A: numbers from the 6 <sample>_credible_discovery.json (script 62).
Builds a deterministic landscape comparison + standalone HTML. NO heavy compute.

ANALYSES
  1. Per-sample credible (pass_tierA) gene set + the gene x sample matrix.
  2. Cross-sample RECURRENCE (genes credible in >=2 samples).
  3. UNDERPOWER null: expected pairwise gene-set overlap under random draws from
     a ~20k-gene universe -> contextualises whether "0 recurrence" is meaningful
     privacy or just small/coverage-limited sets (CRITICAL honesty).
  4. Jaccard similarity matrix between samples (and whether same-cancer pairs
     cluster).
  5. Positional character of credible loci (CpG context, |TSS dist|) — pooled.
     NOTE: CpG-island enrichment is partly induced by the discovery's CpG-island
     PRE-FILTER, so only the within-island-proximal TSS-distance is informative.
  6. Direction (hypo/hyper) + axis (HP1/HP2) of credible deltas.

Output:
  genome_survey_v2/cn_confound/cross_sample/asm_landscape.json
  docs/experiments/in_progress/2026/06/20260603_cross_sample_ASM_landscape_01.standalone.html
"""
import os, json, html, base64, io, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
OUTJSON = f"{CS}/asm_landscape.json"
OUTHTML = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
           "20260603_cross_sample_ASM_landscape_01.standalone.html")

ALL = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}
CCOLOR = {"breast": "#be123c", "lung": "#0e7490", "melanoma": "#7c3aed"}
GENE_UNIVERSE = 20000   # ~ protein-coding+lncRNA for null overlap estimate


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig); return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def esc(x):
    return html.escape(str(x))


# ---- load credible loci ----
cred = {}   # sample -> list of pass_tierA loci dicts
for s in ALL:
    d = json.load(open(f"{CS}/{s}_credible_discovery.json"))
    cred[s] = [r for r in d["credible_loci"] if r["pass_tierA"]]

gene_set = {s: sorted(set(r["nearest_gene"] for r in cred[s])) for s in ALL}
locus_set = {s: set(r["locus"] for r in cred[s]) for s in ALL}

# recurrence (gene-level)
from collections import Counter, defaultdict
gc = Counter()
gsamp = defaultdict(list)
for s in ALL:
    for g in gene_set[s]:
        gc[g] += 1
        gsamp[g].append(s)
union_genes = sorted(gc)
recurrent = {g: gsamp[g] for g in union_genes if gc[g] >= 2}

# underpower null: expected pairwise overlap E = n1*n2/G; prob>=1 ~ 1-exp(-E)
pair_expected = {}
for a, b in itertools.combinations(ALL, 2):
    n1, n2 = len(gene_set[a]), len(gene_set[b])
    E = n1 * n2 / GENE_UNIVERSE
    pair_expected[f"{a}|{b}"] = dict(n1=n1, n2=n2, expected_overlap=round(E, 4),
                                     observed_overlap=len(set(gene_set[a]) & set(gene_set[b])))
sum_expected = round(sum(v["expected_overlap"] for v in pair_expected.values()), 4)
sum_observed = sum(v["observed_overlap"] for v in pair_expected.values())

# Jaccard
def jacc(a, b):
    A, B = set(gene_set[a]), set(gene_set[b])
    return 0.0 if not (A | B) else round(len(A & B) / len(A | B), 4)

jaccard = {f"{a}|{b}": jacc(a, b) for a, b in itertools.combinations(ALL, 2)}

# positional + direction (pooled credible)
all_cred = [r for s in ALL for r in cred[s]]
ctx_counts = Counter(r.get("cpg_context") for r in all_cred)
tss_abs = [abs(r["tss_dist"]) for r in all_cred if r.get("tss_dist") is not None]
axis_counts = Counter(r["axis"] for r in all_cred)
n_hypo = sum(1 for r in all_cred if r["delta"] < 0)
n_hyper = sum(1 for r in all_cred if r["delta"] > 0)

# per-cancer-type union
ct_union = defaultdict(set)
for s in ALL:
    ct_union[CANCER[s]].update(gene_set[s])
ct_union = {k: sorted(v) for k, v in ct_union.items()}

landscape = dict(
    meta=dict(script="64_asm_landscape.py", samples=ALL, cancer_types=CANCER,
              gene_universe_for_null=GENE_UNIVERSE,
              note=("credible = pass_tierA (ARI>=0.30 & placebo<0.10) from per-sample "
                    "discovery (62). CpG-island enrichment is partly induced by the "
                    "discovery CpG-island PRE-FILTER. '0 recurrence' must be read against "
                    "the random-null expected overlap (sets are small + coverage-limited).")),
    per_sample=dict((s, dict(cancer=CANCER[s], n_credible_loci=len(cred[s]),
                             n_credible_genes=len(gene_set[s]), genes=gene_set[s]))
                    for s in ALL),
    recurrence=dict(n_union_genes=len(union_genes), n_recurrent_genes=len(recurrent),
                    recurrent=recurrent),
    underpower_null=dict(sum_expected_overlap=sum_expected, sum_observed_overlap=sum_observed,
                         pairwise=pair_expected,
                         interpretation=("0 observed cross-sample gene recurrence; under "
                                         f"random draws the TOTAL expected pairwise overlap "
                                         f"is only {sum_expected} genes -> 0 observed is "
                                         "CONSISTENT WITH CHANCE and does NOT establish "
                                         "biological privacy (UNDERPOWERED).")),
    jaccard=jaccard,
    positional=dict(cpg_context=dict(ctx_counts),
                    tss_abs_median=(round(float(np.median(tss_abs)), 1) if tss_abs else None),
                    tss_within_2kb_frac=(round(float(np.mean([t <= 2000 for t in tss_abs])), 3)
                                         if tss_abs else None),
                    n_pooled_credible=len(all_cred)),
    direction=dict(axis=dict(axis_counts), n_hypo=n_hypo, n_hyper=n_hyper),
    cancer_type_union={k: len(v) for k, v in ct_union.items()},
)
with open(OUTJSON, "w") as f:
    json.dump(landscape, f, indent=2,
              default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
print(f"[64] wrote {OUTJSON}")
print(f"[64] union_genes={len(union_genes)} recurrent={len(recurrent)} "
      f"sum_expected_overlap={sum_expected} (observed={sum_observed})")
print(f"[64] positional: ctx={dict(ctx_counts)} tss_median={landscape['positional']['tss_abs_median']} "
      f"within2kb={landscape['positional']['tss_within_2kb_frac']}")
print(f"[64] direction: hypo/hyper={n_hypo}/{n_hyper} axis={dict(axis_counts)}")


# ---------------- FIGURES ----------------
def fig_genes_bar():
    fig, ax = plt.subplots(figsize=(8.6, 3.4))
    xs = np.arange(len(ALL))
    ng = [len(gene_set[s]) for s in ALL]
    cols = [CCOLOR[CANCER[s]] for s in ALL]
    ax.bar(xs, ng, color=cols)
    for i, v in enumerate(ng):
        ax.text(i, v + 0.4, str(v), ha="center", fontsize=9, fontweight="bold")
    ax.set_xticks(xs); ax.set_xticklabels([f"{s}\n({CANCER[s]})" for s in ALL], fontsize=8)
    ax.set_ylabel("# credible ASM genes (pass tierA)")
    ax.set_title("Per-sample credible ASM gene set size — 0/94 recur across samples\n"
                 "(all private; but see underpower null: expected random overlap ~0 too)")
    ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


def fig_tss():
    if not tss_abs:
        return None
    fig, ax = plt.subplots(figsize=(8.6, 3.0))
    ax.hist(np.clip(tss_abs, 0, 50000), bins=30, color="#1E3A8A", alpha=0.8, edgecolor="white")
    ax.axvline(2000, color="#b45309", ls="--", lw=1, label="2kb (promoter)")
    ax.set_xlabel("|distance to nearest TSS| (bp, clipped 50kb)")
    ax.set_ylabel("# credible loci")
    ax.set_title(f"Credible ASM loci promoter proximity (pooled n={len(tss_abs)}; "
                 f"{landscape['positional']['tss_within_2kb_frac']:.0%} within 2kb)")
    ax.legend(fontsize=8); ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


FIG_GENES = fig_genes_bar()
FIG_TSS = fig_tss()

# gene matrix (each gene -> its single sample, grouped) — show top per sample
matrix_rows = ""
for s in ALL:
    chip = f'<span class="chip" style="background:{CCOLOR[CANCER[s]]}">{esc(s)}</span>'
    # top credible by |delta|
    top = sorted(cred[s], key=lambda r: -abs(r["delta"]))[:10]
    genes = ", ".join(f'{esc(r["nearest_gene"])}(<span class="mono">Δ{r["delta"]:+.2f}</span>'
                      f'/ari{r["ari"]:.2f})' for r in top)
    matrix_rows += (f'<tr><td>{chip}</td><td>{esc(CANCER[s])}</td>'
                    f'<td class="num">{len(gene_set[s])}</td>'
                    f'<td style="font-size:.78rem">{genes or "—（0 credible）"}</td></tr>')

# jaccard table
jrows = ""
for a, b in itertools.combinations(ALL, 2):
    v = jaccard[f"{a}|{b}"]
    pe = pair_expected[f"{a}|{b}"]
    same = "同癌種" if CANCER[a] == CANCER[b] else ""
    jrows += (f'<tr><td>{esc(a)} × {esc(b)}</td><td>{esc(CANCER[a])}/{esc(CANCER[b])} {same}</td>'
              f'<td class="num">{pe["observed_overlap"]}</td>'
              f'<td class="num">{pe["expected_overlap"]}</td>'
              f'<td class="num">{v}</td></tr>')

ctx_str = " · ".join(f"{k}:{v}" for k, v in ctx_counts.items())

HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>跨樣本 ASM landscape — credible 基因集比對</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1060px;margin:0 auto;padding:24px 18px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:24px 26px;margin-bottom:16px}}
header.top h1{{margin:0 0 6px;font-size:1.4rem}}.meta{{font-size:.82rem;opacity:.88}}
.tldr{{background:#eef2ff;border:1px solid #c7d2fe;border-left:5px solid #1E3A8A;border-radius:12px;padding:16px 20px;margin-bottom:14px}}
.tldr h2{{margin:0 0 8px;font-size:1.06rem;color:#1e3a8a}}.tldr p{{margin:.4rem 0}}
.warn{{background:#fffbeb;border:1px solid #fde68a;border-left:5px solid #b45309;border-radius:12px;padding:14px 18px;margin:14px 0}}
.kpis{{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}}
.kpi{{flex:1;min-width:140px;background:#fff;border:1px solid var(--line);border-radius:12px;padding:14px;text-align:center}}
.kpi .v{{font-size:1.5rem;font-weight:800;color:var(--accent)}}.kpi .l{{font-size:.76rem;color:var(--mut);margin-top:3px}}
.card{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:14px 0}}
.qlabel{{font-weight:800;color:var(--accent);font-size:1rem}}
table{{width:100%;border-collapse:collapse;font-size:.83rem;margin-top:8px}}
th,td{{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}th{{background:#f1f5f9;font-size:.77rem}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}.mono{{font-family:ui-monospace,monospace;font-size:.76rem;color:#475569}}
.chip{{color:#fff;font-weight:700;font-size:.72rem;padding:2px 8px;border-radius:6px}}
.figbox img{{width:100%;border:1px solid var(--line);border-radius:8px;margin-top:8px}}
.legend{{font-size:.76rem;color:var(--mut);margin-top:4px}}.tablewrap{{overflow-x:auto}}
details.method>summary{{cursor:pointer;color:var(--accent);font-weight:600;font-size:.84rem}}
details.method p{{font-size:.83rem;background:var(--bg);border-radius:8px;padding:8px 12px}}
footer{{margin-top:24px;font-size:.76rem;color:var(--mut);border-top:1px solid var(--line);padding-top:12px}}
</style></head><body><div class="wrap">

<header class="top">
  <h1>跨樣本 ASM landscape — credible 基因集比對</h1>
  <div class="meta">6 癌症樣本各自的 credible somatic-ASM 基因集（pass tierA: ARI≥0.30 & placebo&lt;0.10）· A pilot · 2026-06-03 · §13 layer-A</div>
</header>

<div class="tldr">
  <h2>TL;DR — 各樣本 credible ASM 基因集完全不重疊（但統計 underpowered）</h2>
  <p>① <b>union {len(union_genes)} 個 credible 基因，跨樣本復現 = {len(recurrent)}</b>：沒有任何基因在 ≥2 個癌症同時 credible。</p>
  <p>② 各樣本基因數：HCC1395 {len(gene_set["HCC1395"])} · HCC1954 {len(gene_set["HCC1954"])} · H2009 {len(gene_set["H2009"])} · H1437 {len(gene_set["H1437"])} · HCC1937 {len(gene_set["HCC1937"])} · COLO829 {len(gene_set["COLO829"])}（覆蓋限制=0）。</p>
  <p>③ credible loci 位置：{ctx_str}；{(f"{landscape['positional']['tss_within_2kb_frac']:.0%} 在 TSS ±2kb（promoter）" if landscape['positional']['tss_within_2kb_frac'] is not None else "")} —— 但 island 富集部分來自 discovery 的 CpG-island 預篩（非純生物）。</p>
</div>

<div class="warn">
  <b>⚠ 統計誠實（關鍵）</b>：「0 跨樣本復現」<b>不能</b>證明 ASM 基因 private —— 各樣本只 {len(gene_set["HCC1937"])}–{len(gene_set["HCC1395"])} 個 credible 基因，從 ~{GENE_UNIVERSE//1000}k 基因中抽，<b>隨機期望的總成對重疊也僅 {sum_expected} 個基因</b>（觀測 {sum_observed}）。故 0 觀測<b>與純隨機一致</b>、是 <b>UNDERPOWERED</b>（基因集太小 + 覆蓋限制），非生物 privacy 的證據。真正的 private 證據在 targeted 層（同位點 somatic 復發 0/38，見前一份報告）。
</div>

<div class="kpis">
  <div class="kpi"><div class="v">{len(recurrent)} / {len(union_genes)}</div><div class="l">跨樣本復現 credible 基因</div></div>
  <div class="kpi"><div class="v">{sum_expected}</div><div class="l">隨機 null 期望總重疊</div></div>
  <div class="kpi"><div class="v">{landscape['positional']['tss_within_2kb_frac']:.0%}</div><div class="l">credible 在 promoter ±2kb</div></div>
  <div class="kpi"><div class="v">{n_hypo}/{n_hyper}</div><div class="l">hypo/hyper credible</div></div>
</div>

<section class="card">
  <span class="qlabel">① 各樣本 credible 基因數 + top 基因</span>
  <div class="tablewrap"><table>
    <thead><tr><th>樣本</th><th>癌種</th><th class="num"># genes</th><th>top credible 基因（|Δβ| 排序）</th></tr></thead>
    <tbody>{matrix_rows}</tbody></table></div>
  {(f'<div class="figbox"><img src="{FIG_GENES}"></div>') if FIG_GENES else ""}
</section>

<section class="card">
  <span class="qlabel">② 成對重疊 / Jaccard（含隨機期望對照）</span>
  <div class="tablewrap"><table>
    <thead><tr><th>樣本對</th><th>癌種</th><th class="num">觀測重疊</th><th class="num">隨機期望</th><th class="num">Jaccard</th></tr></thead>
    <tbody>{jrows}</tbody></table></div>
  <p class="legend">所有對觀測重疊=0、Jaccard=0；隨機期望本就接近 0 → 無法區分同癌種 vs 跨癌種相似度（underpowered）。</p>
</section>

<section class="card">
  <span class="qlabel">③ credible loci 位置/方向特性（pooled n={len(all_cred)}）</span>
  <p style="font-size:.85rem">CpG context: {ctx_str}；|TSS dist| median = {landscape['positional']['tss_abs_median']} bp；hypo/hyper = {n_hypo}/{n_hyper}；axis {dict(axis_counts)}。</p>
  {(f'<div class="figbox"><img src="{FIG_TSS}"></div>') if FIG_TSS else ""}
  <p class="legend">⚠ island/shore 富集部分由 CpG-island 預篩誘導；promoter(±2kb) proximity 在預篩內仍 informative。</p>
</section>

<section class="card">
  <span class="qlabel">④ 方法 / 限制</span>
  <details class="method" open><summary>讀法 / caveat（必讀）</summary><p>
  <b>資料</b>：各樣本 credible = per-sample discovery（script 62）pass_tierA loci 的 nearest_gene。<br>
  <b>核心限制</b>：① <b>「0 復現」underpowered</b> —— credible 集太小（3–30 基因）+ 覆蓋限制，隨機期望重疊本就 ~{sum_expected}，0 觀測不證明 private（真 private 證據在 targeted 同位點 0/38 層）；② <b>COLO829=0</b>（melanoma 覆蓋最低）對 landscape 無貢獻；③ <b>island 富集是預篩誘導</b>非純生物；④ credible 基因偏 LOC/lncRNA，本 pilot 不宣稱功能/TSG 富集；⑤ 同管線同 caller（現象層級非獨立驗證）。<br>
  <b>誠實結論</b>：各樣本<b>都有</b>自己的 credible ASM 基因集（discovery 有效），credible loci promoter-proximal；但<b>跨樣本基因復現無法在此覆蓋/集大小下檢定</b>。</p></details>
</section>

<footer>
  數據源：<span class="mono">research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/asm_landscape.json</span>（+ 6× *_credible_discovery.json）· 腳本 64 · §13 layer-A
</footer>
</div></body></html>'''

os.makedirs(os.path.dirname(OUTHTML), exist_ok=True)
with open(OUTHTML, "w") as f:
    f.write(HTML_DOC)
print(f"[64] wrote {OUTHTML} ({len(HTML_DOC)//1024} KB)")
