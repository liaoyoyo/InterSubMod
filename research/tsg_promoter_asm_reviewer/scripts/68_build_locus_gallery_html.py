#!/usr/bin/env python3
"""
68 - Credible ASM locus GALLERY HTML: per-sample, per-locus methylation heatmaps
organized by tier + ARI, for complete visual confirmation across samples.

§13 layer-A: figures + numbers from <sample>_locus_figs.json (67) +
<sample>_credible_discovery.json funnels (62).

Output: docs/experiments/in_progress/2026/06/
        20260603_credible_ASM_locus_gallery_01.standalone.html
"""
import os, json, html
import numpy as np

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
       "20260603_credible_ASM_locus_gallery_01.standalone.html")

ALL = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CANCER_FULL = {"HCC1395": "breast (SEQC2 ref)", "HCC1937": "breast (BRCA1-mut)",
               "HCC1954": "breast", "H1437": "lung adeno", "H2009": "lung adeno",
               "COLO829": "melanoma"}
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}
CCOLOR = {"breast": "#be123c", "lung": "#0e7490", "melanoma": "#7c3aed"}


def esc(x):
    return html.escape(str(x))


figs = {}
for s in ALL:
    p = f"{CS}/{s}_locus_figs.json"
    figs[s] = json.load(open(p)) if os.path.exists(p) else None

cred_n = {}
for s in ALL:
    p = f"{CS}/{s}_credible_discovery.json"
    if os.path.exists(p):
        fn = json.load(open(p))["funnel"]
        cred_n[s] = (fn["n_credible_pass_tierA_relax"], fn["n_credible_pass_tierA_strict"],
                     fn["n_ari_evaluable"], fn["n_hp_axis_survey"])
    else:
        cred_n[s] = (0, 0, 0, 0)

total_rendered = sum(figs[s]["n_rendered"] for s in ALL if figs[s])
n_samples_with = sum(1 for s in ALL if figs[s] and figs[s]["n_rendered"] > 0)

# summary table
sum_rows = ""
for s in ALL:
    cr = cred_n[s]
    nrend = figs[s]["n_rendered"] if figs[s] else 0
    chip = f'<span class="chip" style="background:{CCOLOR[CANCER[s]]}">{esc(s)}</span>'
    note = "覆蓋限制 → 0 credible" if cr[0] == 0 else f"{nrend} 圖"
    sum_rows += (f'<tr><td>{chip}</td><td>{esc(CANCER_FULL[s])}</td>'
                 f'<td class="num">{cr[3]}</td><td class="num">{cr[2]}</td>'
                 f'<td class="num"><b>{cr[0]}/{cr[1]}</b></td><td class="num">{note}</td></tr>')

# per-sample galleries
galleries = ""
for s in ALL:
    chip = f'<span class="chip" style="background:{CCOLOR[CANCER[s]]}">{esc(s)}</span>'
    if not figs[s] or figs[s]["n_rendered"] == 0:
        galleries += (f'<section class="card"><div class="gh">{chip} '
                      f'<span class="gl">{esc(CANCER_FULL[s])}</span> '
                      f'<span class="badge warn">0 credible（覆蓋限制，中位覆蓋最低）</span></div>'
                      f'<p class="legend">此樣本無 pass_tierA credible locus → 無圖。'
                      f'非生物缺席，是 somatic 子單倍型 reads &lt; 8（覆蓋限制）使 ARI 不可評估。</p></section>')
        continue
    F = figs[s]["figs"]
    # split strict / relax-only
    strict = [x for x in F if x.get("tier") == "strict_ge100"]
    relax = [x for x in F if x.get("tier") != "strict_ge100"]

    def cards(lst):
        out = ""
        for x in lst:
            plcb = x["placebo_ari"] if x["placebo_ari"] is not None else "—"
            out += (f'<figure class="locus"><img src="{x["fig"]}" loading="lazy">'
                    f'<figcaption><b>{esc(x["gene"])}</b> '
                    f'<span class="mono">{esc(x["locus"])}</span><br>'
                    f'Δβ={x["delta"]:+.2f} · ARI={x["ari"]:.2f} · plcb={plcb} · '
                    f'{x["n_paired_cpg"]}CpG · {esc(x.get("cpg_context"))}</figcaption></figure>')
        return out

    block = (f'<div class="gh">{chip} <span class="gl">{esc(CANCER_FULL[s])}</span> '
             f'<span class="badge">{len(F)} credible（strict {len(strict)} / relax-only {len(relax)}）</span></div>')
    if strict:
        block += f'<h4 class="tierhdr">strict tier（n_paired_cpg≥100，依 ARI 排序）</h4><div class="grid">{cards(strict)}</div>'
    if relax:
        block += f'<h4 class="tierhdr">relaxed tier（n_paired_cpg 30–99）</h4><div class="grid">{cards(relax)}</div>'
    galleries += f'<section class="card">{block}</section>'

HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Credible ASM locus gallery — 逐位點甲基觀察</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1180px;margin:0 auto;padding:24px 18px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:24px 26px;margin-bottom:16px}}
header.top h1{{margin:0 0 6px;font-size:1.4rem}}.meta{{font-size:.82rem;opacity:.88}}
.tldr{{background:#eef2ff;border:1px solid #c7d2fe;border-left:5px solid #1E3A8A;border-radius:12px;padding:16px 20px;margin-bottom:14px}}
.tldr h2{{margin:0 0 8px;font-size:1.05rem;color:#1e3a8a}}.tldr p{{margin:.35rem 0}}
.legendbox{{background:#fff;border:1px solid var(--line);border-radius:10px;padding:10px 14px;margin:12px 0;font-size:.82rem;display:flex;gap:16px;flex-wrap:wrap;align-items:center}}
.sw{{display:inline-block;width:14px;height:14px;border-radius:3px;vertical-align:middle;margin-right:4px;border:1px solid #ccc}}
.card{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:16px 18px;margin:14px 0}}
.gh{{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:6px}}
.gl{{font-weight:800;color:var(--accent);font-size:1rem;flex:1}}
.badge{{color:#fff;font-weight:700;font-size:.74rem;padding:3px 11px;border-radius:999px;background:#15803d}}
.badge.warn{{background:#b45309}}
.chip{{color:#fff;font-weight:700;font-size:.72rem;padding:2px 8px;border-radius:6px}}
.tierhdr{{margin:12px 0 6px;font-size:.86rem;color:#374151;border-bottom:1px solid var(--line);padding-bottom:3px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(240px,1fr));gap:12px}}
figure.locus{{margin:0;background:var(--bg);border:1px solid var(--line);border-radius:8px;padding:6px}}
figure.locus img{{width:100%;border-radius:4px;display:block}}
figcaption{{font-size:.72rem;color:var(--mut);margin-top:4px;line-height:1.4}}
.mono{{font-family:ui-monospace,monospace;font-size:.7rem}}
table{{width:100%;border-collapse:collapse;font-size:.83rem;margin-top:8px}}
th,td{{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left}}th{{background:#f1f5f9;font-size:.77rem}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}.legend{{font-size:.78rem;color:var(--mut)}}
details.method>summary{{cursor:pointer;color:var(--accent);font-weight:600;font-size:.84rem}}
details.method p{{font-size:.83rem;background:var(--bg);border-radius:8px;padding:8px 12px}}
footer{{margin-top:24px;font-size:.76rem;color:var(--mut);border-top:1px solid var(--line);padding-top:12px}}
</style></head><body><div class="wrap">

<header class="top">
  <h1>Credible ASM locus gallery — 逐位點甲基觀察</h1>
  <div class="meta">各樣本每個 credible（pass tierA）位點的 read×CpG 甲基矩陣，依 HP 分群供肉眼確認 ASM clustering ·
  {total_rendered} 圖 / {n_samples_with} 樣本 · A pilot · 2026-06-03 · §13 layer-A</div>
</header>

<div class="tldr">
  <h2>用途 — 完整視覺確認各樣本是否真有 ASM clustering</h2>
  <p>每張圖 = 一個 credible locus 的<b>實際甲基觀察</b>：rows=reads（上=germline 主單倍型、下=somatic 子單倍型 HPx-1、再下=另一 HP），cols=CpG，組內依平均甲基排序。</p>
  <p>ASM 表現 = somatic 區塊與 germline 區塊<b>整體甲基水平差異</b>（Δβ）+ read-level <b>clustering 分離</b>（ARI≥0.30，placebo&lt;0.10 排除長度 artifact）。橘虛線=somatic 變異位置。</p>
  <p>共 <b>{total_rendered}</b> 個 credible 位點圖（{n_samples_with}/6 樣本有 credible；<b>COLO829 覆蓋限制 0 credible</b>）。</p>
</div>

<div class="legendbox">
  <b>圖例：</b>
  <span><span class="sw" style="background:#dc2626"></span>甲基化(5mC≥0.78)</span>
  <span><span class="sw" style="background:#2563eb"></span>非甲基(≤0.20)</span>
  <span><span class="sw" style="background:#d1d5db"></span>中間</span>
  <span><span class="sw" style="background:#ffffff"></span>無 call</span>
  <span><span class="sw" style="background:#f59e0b"></span>變異位置</span>
</div>

<section class="card">
  <div class="gh"><span class="gl">各樣本 credible 摘要</span></div>
  <table>
    <thead><tr><th>樣本</th><th>癌種</th><th class="num">HP-axis survey</th><th class="num">ARI-eval</th>
      <th class="num">credible relax/strict</th><th class="num">圖</th></tr></thead>
    <tbody>{sum_rows}</tbody></table>
</section>

{galleries}

<section class="card">
  <div class="gh"><span class="gl">方法 / 限制</span></div>
  <details class="method" open><summary>讀法 / caveat</summary><p>
  <b>圖</b>：window ±600bp，每 read 用 pysam r.modified_bases（5mC，ML≥200 甲基/≤50 非甲基/中間），依 HP:Z tag 分群（germline 主 HP / somatic HPx-1 / 另一 HP），每組最多抽 45 reads，組內依平均甲基排序。<br>
  <b>判讀</b>：① somatic 區塊整體偏紅(或偏藍)、germline 區塊相反 = ASM（Δβ 符號）；② 小 Δβ 但高 ARI = 甲基 pattern（非僅水平）clustering 分離；③ somatic 區塊常很小（低 AF 子單倍型 reads 少）—— 這是真實覆蓋現象。<br>
  <b>限制</b>：① 圖為 binarized 顯示（連續 β 在 discovery 計算）；② 子群抽樣 ≤45 reads 僅供視覺；③ 這些是各樣本<b>自己</b>的 credible 位點（非同位點跨樣本，各癌 private）；④ COLO829 覆蓋限制 0 credible 無圖；⑤ strict/relax = n_paired_cpg≥100 / 30–99。</p></details>
</section>

<footer>
  數據源：<span class="mono">research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/*_locus_figs.json</span>
  （+ *_credible_discovery.json）· 腳本 67（圖）/68（gallery）· §13 layer-A
</footer>
</div></body></html>'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(HTML_DOC)
print(f"[68] wrote {OUT} ({len(HTML_DOC)//1024} KB, {total_rendered} locus figures)")
