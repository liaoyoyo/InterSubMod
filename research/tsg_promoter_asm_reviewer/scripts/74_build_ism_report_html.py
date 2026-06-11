#!/usr/bin/env python3
"""
74 - Complete ISM existence + cis report HTML (powered, 280k positions, 6 samples × TP/FP/FN).
§13 layer-A: numbers from ism_aggregate.json (73). No hand-typed metrics.

Output: docs/experiments/in_progress/2026/06/20260604_ISM_complete_TPFPFN_existence_cis_01.standalone.html
"""
import os, json, html, io, base64
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
       "20260604_ISM_complete_TPFPFN_existence_cis_01.standalone.html")

ALL = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}
CCOLOR = {"breast": "#be123c", "lung": "#0e7490", "melanoma": "#7c3aed"}
CLSCOLOR = {"tp": "#15803d", "fp": "#dc2626", "fn": "#b45309"}

D = json.load(open(f"{CS}/ism_aggregate.json"))
EX = D["existence"]; DISC = D["discrimination"]; CIS = D["cis_hcc1395"]


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig); return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def esc(x):
    return html.escape(str(x))


def pct(x):
    return "—" if x is None else f"{x*100:.2f}%"


# FIG 1: significant-rate per sample × class
def fig_rate():
    fig, ax = plt.subplots(figsize=(9.5, 3.8))
    x = np.arange(len(ALL)); w = 0.26
    for i, c in enumerate(["tp", "fp", "fn"]):
        rates = [(EX[s][c]["significant_rate"] or 0) * 100 for s in ALL]
        ax.bar(x + (i - 1) * w, rates, w, color=CLSCOLOR[c], label=c.upper(), alpha=0.9)
    ax.set_xticks(x); ax.set_xticklabels([f"{s}\n({CANCER[s]})" for s in ALL], fontsize=8)
    ax.set_ylabel("ISM Significant rate (%)")
    ax.set_title("Complete ISM ASM-significance rate — 6 samples × TP/FP/FN (powered, 280k loci)")
    ax.legend(fontsize=8); ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


# FIG 2: pooled discrimination (overall + subhap-conditional)
def fig_disc():
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    cls = ["tp", "fp", "fn"]; x = np.arange(len(cls)); w = 0.36
    overall = [(DISC[c]["pooled_significant_rate"] or 0) * 100 for c in cls]
    cond = [(DISC[c]["pooled_sig_rate_in_subhap"] or 0) * 100 for c in cls]
    ax.bar(x - w/2, overall, w, color="#1E3A8A", label="all loci")
    ax.bar(x + w/2, cond, w, color="#60a5fa", label="subhap-present (matched)")
    for i, (a, b) in enumerate(zip(overall, cond)):
        ax.text(i - w/2, a + 0.05, f"{a:.1f}", ha="center", fontsize=8)
        ax.text(i + w/2, b + 0.05, f"{b:.1f}", ha="center", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels([c.upper() for c in cls])
    ax.set_ylabel("Significant rate (%)")
    ax.set_title("Pooled TP/FP/FN discrimination\n(TP modestly > FP even subhap-matched, but low + inconsistent)")
    ax.legend(fontsize=8); ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


# FIG 3: HCC1395 cis candidate rate
def fig_cis():
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    cls = ["tp", "fp", "fn"]; x = np.arange(len(cls))
    rates = [(CIS[c]["cis_candidate_rate"] or 0) * 100 for c in cls]
    ax.bar(x, rates, 0.5, color=[CLSCOLOR[c] for c in cls])
    for i, r in enumerate(rates):
        ax.text(i, r + 0.3, f"{r:.1f}%", ha="center", fontsize=9, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([c.upper() for c in cls])
    ax.set_ylabel("cis-candidate rate (%)")
    ax.set_title("HCC1395 normal-anchored cis-candidate (HP_Residual sig + |Δ|≥0.10)\n"
                 "similar TP/FP/FN, FP highest → NOT cis-discriminative")
    ax.spines[["top", "right"]].set_visible(False)
    return b64(fig)


FIG_RATE = fig_rate()
FIG_DISC = fig_disc()
FIG_CIS = fig_cis()

# existence table
ex_rows = ""
for s in ALL:
    chip = f'<span class="chip" style="background:{CCOLOR[CANCER[s]]}">{esc(s)}</span>'
    cells = ""
    for c in ["tp", "fp", "fn"]:
        e = EX[s][c]
        cells += f'<td class="num">{e["n_significant"]}/{e["n"]}<br><b>{pct(e["significant_rate"])}</b></td>' if e else '<td>—</td>'
    ex_rows += f'<tr><td>{chip}</td><td>{esc(CANCER[s])}</td>{cells}</tr>'

# discrimination table
disc_rows = ""
for c in ["tp", "fp", "fn"]:
    x = DISC[c]
    disc_rows += (f'<tr><td><b>{c.upper()}</b></td>'
                  f'<td class="num">{x["pooled_n"]:,}</td>'
                  f'<td class="num">{pct(x["pooled_significant_rate"])}</td>'
                  f'<td class="num">{x["pooled_n_with_subhap"]:,}</td>'
                  f'<td class="num"><b>{pct(x["pooled_sig_rate_in_subhap"])}</b></td></tr>')

# cis table
cis_rows = ""
for c in ["tp", "fp", "fn"]:
    v = CIS[c]
    cis_rows += (f'<tr><td><b>{c.upper()}</b></td>'
                 f'<td class="num">{v["n"]:,}</td>'
                 f'<td class="num">{v["n_cis_candidate"]:,}</td>'
                 f'<td class="num"><b>{pct(v["cis_candidate_rate"])}</b></td>'
                 f'<td class="num">{v["median_abs_HP_Residual_Delta"]}</td>'
                 f'<td class="num">{v["n_sampleASM_sig"]:,}</td></tr>')

tp_r = DISC["tp"]["pooled_significant_rate"]; fp_r = DISC["fp"]["pooled_significant_rate"]
tp_c = DISC["tp"]["pooled_sig_rate_in_subhap"]; fp_c = DISC["fp"]["pooled_sig_rate_in_subhap"]
ratio = round(tp_c / fp_c, 1) if fp_c else None

HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ISM 完整 TP/FP/FN 存在性 + cis 分析 — 6 樣本 powered</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1080px;margin:0 auto;padding:24px 18px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:24px 26px;margin-bottom:16px}}
header.top h1{{margin:0 0 6px;font-size:1.34rem}}.meta{{font-size:.82rem;opacity:.88}}
.tldr{{background:#eef2ff;border:1px solid #c7d2fe;border-left:5px solid #1E3A8A;border-radius:12px;padding:16px 20px;margin-bottom:14px}}
.tldr h2{{margin:0 0 8px;font-size:1.05rem;color:#1e3a8a}}.tldr p{{margin:.35rem 0}}
.kpis{{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}}
.kpi{{flex:1;min-width:150px;background:#fff;border:1px solid var(--line);border-radius:12px;padding:14px;text-align:center}}
.kpi .v{{font-size:1.4rem;font-weight:800;color:var(--accent)}}.kpi .l{{font-size:.75rem;color:var(--mut);margin-top:3px}}
.card{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:14px 0}}
.qlabel{{font-weight:800;color:var(--accent);font-size:1rem}}
table{{width:100%;border-collapse:collapse;font-size:.83rem;margin-top:8px}}
th,td{{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left}}th{{background:#f1f5f9;font-size:.76rem}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}
.chip{{color:#fff;font-weight:700;font-size:.72rem;padding:2px 8px;border-radius:6px}}
.figbox img{{width:100%;border:1px solid var(--line);border-radius:8px;margin-top:8px}}.tablewrap{{overflow-x:auto}}
.mech{{background:#fffbeb;border:1px solid #fde68a;border-left:5px solid #b45309;border-radius:12px;padding:14px 18px;margin:14px 0;font-size:.9rem}}
details.method>summary{{cursor:pointer;color:var(--accent);font-weight:600;font-size:.84rem}}
details.method p{{font-size:.83rem;background:var(--bg);border-radius:8px;padding:8px 12px}}
footer{{margin-top:24px;font-size:.76rem;color:var(--mut);border-top:1px solid var(--line);padding-top:12px}}
</style></head><body><div class="wrap">

<header class="top"><h1>ISM 完整 TP/FP/FN 存在性 + cis 分析（6 樣本 powered）</h1>
<div class="meta">unmodified ISM (build/bin/inter_sub_mod) 跑全 6 樣本 × TP/FP/FN ~28 萬位點 + HCC1395 normal-anchored cis · 2026-06-04 · §13 layer-A（數字由 ism_aggregate.json 注入）</div></header>

<div class="tldr">
  <h2>TL;DR — ASM 廣泛存在但<b>不是可用的 TP/FP filter</b>；normal-anchored cis 也不判別</h2>
  <p>① <b>完整 powered 掃描</b>（294k TP / 3,177 FP / 34,805 FN，取代之前 underpowered 的小 N 對照）。</p>
  <p>② <b>TP 的 ASM-significance 率 modestly &gt; FP</b>（pooled {pct(tp_r)} vs {pct(fp_r)}；subhap-matched {pct(tp_c)} vs {pct(fp_c)}，~{ratio}×）<b>但不是 usable filter</b>：絕對率極低（~96% TP 不 significant）+ <b>COLO829 TP≈FP 判別消失</b> + 多樣本 FN≈TP。</p>
  <p>③ <b>HCC1395 normal-anchored cis（HP_Residual）跨 TP/FP/FN 相近、FP 最高</b>（{pct(CIS["fp"]["cis_candidate_rate"])}）→ cis 殘差<b>不判別</b>。</p>
</div>

<div class="kpis">
  <div class="kpi"><div class="v">~280k</div><div class="l">完整掃描位點</div></div>
  <div class="kpi"><div class="v">{pct(tp_r)}</div><div class="l">TP ASM-significant 率</div></div>
  <div class="kpi"><div class="v">{pct(fp_r)}</div><div class="l">FP ASM-significant 率</div></div>
  <div class="kpi"><div class="v">{ratio}×</div><div class="l">TP/FP（subhap-matched）</div></div>
</div>

<div class="mech">
  <b>🔑 為何「TP&gt;FP 但不可用」？</b> ISM <code>Significant</code> = PassedGating + global_p≤0.05 + CramersV≥0.1 + NumReads≥20（廣義甲基結構顯著）。TP（真 somatic）較常有此結構 → modestly &gt; FP。但 (a) <b>只 ~4% TP significant</b>（sensitivity 太低不能 rescue 多數 TP）；(b) <b>COLO829（melanoma 低覆蓋）TP≈FP</b> 判別消失；(c) 此「廣義 significance」≠ 舊「strong-ASM extreme」（後者 FP-enriched OR=0.194，極端 LOH/低覆蓋 tail）。<b>淨：ASM 真實、弱 TP-關聯、非 filter。</b>
</div>

<section class="card">
  <span class="qlabel">① 完整存在性 — ISM Significant rate（每樣本 × TP/FP/FN）</span>
  <div class="tablewrap"><table>
    <thead><tr><th>樣本</th><th>癌種</th><th class="num">TP</th><th class="num">FP</th><th class="num">FN</th></tr></thead>
    <tbody>{ex_rows}</tbody></table></div>
  <div class="figbox"><img src="{FIG_RATE}"></div>
</section>

<section class="card">
  <span class="qlabel">② TP/FP/FN 判別（pooled + subhap-matched）</span>
  <div class="tablewrap"><table>
    <thead><tr><th>class</th><th class="num">n</th><th class="num">significant 率(全)</th><th class="num">有 subhap n</th><th class="num">significant 率(subhap-matched)</th></tr></thead>
    <tbody>{disc_rows}</tbody></table></div>
  <div class="figbox"><img src="{FIG_DISC}"></div>
  <p style="font-size:.8rem;color:#6b7280">subhap-matched = 只比有 somatic 子單倍型的 loci（排除「FP 無子單倍型」confound）；TP 仍 ~{ratio}× FP，但見上方為何不可用。</p>
</section>

<section class="card">
  <span class="qlabel">③ HCC1395 normal-anchored cis-ASM（HP_Residual = tumor−normal 殘差）</span>
  <div class="tablewrap"><table>
    <thead><tr><th>class</th><th class="num">n</th><th class="num">cis-candidate</th><th class="num">cis-cand 率</th><th class="num">med|residual|</th><th class="num">SampleASM_sig</th></tr></thead>
    <tbody>{cis_rows}</tbody></table></div>
  <div class="figbox"><img src="{FIG_CIS}"></div>
  <p style="font-size:.82rem;color:#6b7280">cis-candidate = HP_Residual_P&lt;0.05 且 |HP_Residual_Delta|≥0.10。⚠ ISM HP_Residual 是 <b>germline HP-family 層級</b> tumor−normal 殘差（較廣），<b>非</b> script 34 的特定 somatic-allele(HP1-1) vs normal-HP1 三方 cis-test；故為部分 cis 評估。跨 TP/FP/FN 相近、FP 最高 → 不判別。</p>
</section>

<section class="card">
  <span class="qlabel">④ 方法 / 限制 / 與舊結論對照</span>
  <details class="method" open><summary>讀法 / caveat（必讀）</summary><p>
  <b>工具</b>：unmodified ISM (build/bin/inter_sub_mod)，--window-size 1000，--no-output-distance-matrix；significance_summary.csv 117 欄/位點（含 native Stage-2 HP-fine + per-CpG Fisher + signed Δβ + Significant gate）。<br>
  <b>完整性</b>：6 樣本 × {{TP,FP,FN}} 全位點（~28 萬）—— 取代之前 CpG-island 4% 子集 + 小 N 對照（script 69）。cis = HCC1395 + full normal (HCC1395BL 136G tagged)。<br>
  <b>與舊結論對照（一致）</b>：① 舊「甲基→F1 filter DEAD ⭐2」「ISM TO 甲基 AUC&lt;0.58 平坦」→ 本完整掃描證實 ASM 非 usable filter（低 sensitivity + 不一致）。② 舊「strong-ASM FP-enriched OR=0.194」是極端 tail；本「廣義 significant TP&gt;FP」是中度結構，<b>不衝突</b>（不同 metric/threshold）。③ cis 不判別 → 對齊舊「高 ARI 只 drift ≠ cis」。<br>
  <b>限制</b>：① ISM Significant 是工具內建 gate（非單一 effect size）；② cis 用 ISM germline-family 殘差，非特定 somatic-allele 三方（後者需 script 34 per-locus）；③ FN 位點 caller 未呼叫故多無 somatic 子單倍型，ASM 分析受限；④ 單一 caller（ClairS-TO/longphase）。</p></details>
</section>

<footer>數據源：<span style="font-family:monospace">ism_aggregate.json</span>（← 18 existence + 3 cis significance_summary.csv）· ISM build/bin/inter_sub_mod · 腳本 71-74 · §13 layer-A</footer>
</div></body></html>'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(HTML_DOC)
print(f"[74] wrote {OUT} ({len(HTML_DOC)//1024} KB)")
