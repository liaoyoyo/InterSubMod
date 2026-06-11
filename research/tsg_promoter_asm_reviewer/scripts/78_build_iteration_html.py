#!/usr/bin/env python3
"""
78 - Final explanatory HTML for the meaningful-ASM filter iteration experiment.
§13 layer-A: numbers from iteration_history.json (77). Explains definition, each
iteration's rationale/hypothesis/result/diff/verification, the key finding +
mechanism, best pick + caveat + rollback.

Output: docs/experiments/in_progress/2026/06/20260607_meaningful_ASM_filter_iteration_01.standalone.html
"""
import os, json, html, io, base64
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CS = ("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/"
      "genome_survey_v2/cn_confound/cross_sample")
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
       "20260607_meaningful_ASM_filter_iteration_01.standalone.html")
D = json.load(open(f"{CS}/iteration_history.json"))
ITS = D["iterations"]
meta = D["meta"]
best = D["best_pick"]


def esc(x):
    return html.escape(str(x))


def b64(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig); return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def fig_journey():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.8))
    ids = [r["id"] for r in ITS]
    x = np.arange(len(ids))
    tp = [(r["pooled_tp_rate"] or 0) * 100 for r in ITS]
    fp = [(r["pooled_fp_rate"] or 0) * 100 for r in ITS]
    a1.bar(x - 0.2, tp, 0.38, color="#15803d", label="TP %")
    a1.bar(x + 0.2, fp, 0.38, color="#dc2626", label="FP %")
    for i, (t, f) in enumerate(zip(tp, fp)):
        a1.text(i - 0.2, t + 0.1, f"{t:.1f}", ha="center", fontsize=7)
        a1.text(i + 0.2, f + 0.1, f"{f:.1f}", ha="center", fontsize=7)
    a1.set_xticks(x); a1.set_xticklabels(ids); a1.set_ylabel("pass rate (%)")
    a1.set_title("TP vs FP 通過率（每迭代）"); a1.legend(fontsize=8)
    a1.spines[["top", "right"]].set_visible(False)
    ratio = [r["tp_fp_ratio"] or 0 for r in ITS]
    a2.plot(x, ratio, "o-", color="#1E3A8A", markersize=7)
    a2.axhline(1.0, color="#9ca3af", ls="--", lw=1, label="比=1（無判別）")
    for i, rr in enumerate(ratio):
        a2.text(i, rr + 0.08, f"{rr}", ha="center", fontsize=8, fontweight="bold")
    a2.set_xticks(x); a2.set_xticklabels(ids); a2.set_ylabel("TP/FP 分離比")
    a2.set_title("加 Δβ branch → 判別比崩潰"); a2.legend(fontsize=8)
    a2.spines[["top", "right"]].set_visible(False)
    return b64(fig)

FIG = fig_journey()

# iteration table
rows = ""
for r in ITS:
    mark = ' <span class="best">★ best</span>' if r["id"] == best["id"] else (' <span class="froze">凍結</span>' if r["id"] == "I0" else "")
    rows += (f'<tr><td><b>{esc(r["id"])}</b>{mark}<br><span class="sub">{esc(r["name"])}</span></td>'
             f'<td class="num">{(r["pooled_tp_rate"] or 0)*100:.2f}%<br><span class="sub">{r["pooled_tp_pass"]}</span></td>'
             f'<td class="num">{(r["pooled_fp_rate"] or 0)*100:.2f}%<br><span class="sub">{r["pooled_fp_pass"]}</span></td>'
             f'<td class="num"><b>{r["tp_fp_ratio"]}</b></td>'
             f'<td class="num">+{r["added_vs_I0"]}<br><span class="sub">-{r["removed_vs_I0"]}</span></td>'
             f'<td class="sub">cv:{r["branch_cramers_only"]}<br>Δβ:{r["branch_dbeta_only"]}<br>both:{r["branch_both"]}</td></tr>')

# per-iteration cards
cards = ""
for r in ITS:
    ver = "".join(f'<li>{esc(v)}</li>' for v in r.get("verification", []))
    cards += f'''<section class="card">
      <div class="ch"><span class="cid">{esc(r["id"])}</span> <span class="cn">{esc(r["name"])}</span></div>
      <p><b>理由</b>：{esc(r["rationale"])}</p>
      <p><b>假設</b>：{esc(r["hypothesis"])}</p>
      <p class="mono2"><b>filter</b>：{esc(r["filter_desc"])}</p>
      <p><b>結果</b>：TP {r["pooled_tp_pass"]}/{r["pooled_tp_n"]} ({(r["pooled_tp_rate"] or 0)*100:.2f}%) · FP {r["pooled_fp_pass"]}/{r["pooled_fp_n"]} ({(r["pooled_fp_rate"] or 0)*100:.2f}%) · <b>TP/FP={r["tp_fp_ratio"]}</b> · vs I0 +{r["added_vs_I0"]}/-{r["removed_vs_I0"]}</p>
      <details open><summary>驗證 / 假設是否成立 / 例外</summary><ul>{ver}</ul></details>
    </section>'''

HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Meaningful-ASM filter 迭代實驗 — 講解 + 驗證</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1060px;margin:0 auto;padding:24px 18px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:22px 26px;margin-bottom:14px}}
header.top h1{{margin:0 0 5px;font-size:1.32rem}}.meta{{font-size:.82rem;opacity:.88}}
.def{{background:#eef2ff;border:1px solid #c7d2fe;border-left:5px solid #1E3A8A;border-radius:12px;padding:14px 18px;margin-bottom:14px}}
.def h2{{margin:0 0 6px;font-size:1.02rem;color:#1e3a8a}}
.find{{background:#fef2f2;border:1px solid #fecaca;border-left:5px solid #dc2626;border-radius:12px;padding:14px 18px;margin:14px 0}}
.find h2{{margin:0 0 6px;font-size:1.02rem;color:#b91c1c}}
.keep{{background:#f0fdf4;border:1px solid #bbf7d0;border-left:5px solid #16a34a;border-radius:12px;padding:14px 18px;margin:14px 0}}
.card{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:14px 18px;margin:12px 0}}
.ch{{margin-bottom:6px}}.cid{{background:#1E3A8A;color:#fff;font-weight:800;padding:2px 9px;border-radius:6px;font-size:.85rem}}.cn{{font-weight:700;color:var(--accent)}}
.card p{{margin:.35rem 0;font-size:.88rem}}.mono2{{font-family:ui-monospace,monospace;font-size:.78rem;background:var(--bg);padding:4px 8px;border-radius:6px}}
table{{width:100%;border-collapse:collapse;font-size:.84rem;margin-top:8px}}
th,td{{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}th{{background:#f1f5f9;font-size:.77rem}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}.sub{{color:var(--mut);font-size:.74rem}}
.best{{background:#15803d;color:#fff;font-size:.68rem;padding:1px 6px;border-radius:5px}}.froze{{background:#6b7280;color:#fff;font-size:.68rem;padding:1px 6px;border-radius:5px}}
.figbox img{{width:100%;border:1px solid var(--line);border-radius:8px;margin-top:8px}}
details summary{{cursor:pointer;color:var(--accent);font-weight:600;font-size:.82rem}}details ul{{font-size:.82rem;margin:6px 0}}
footer{{margin-top:22px;font-size:.76rem;color:var(--mut);border-top:1px solid var(--line);padding-top:12px}}
</style></head><body><div class="wrap">

<header class="top"><h1>「有甲基意義」位點 filter 迭代實驗 — 講解 + 驗證</h1>
<div class="meta">目標 (c) 聯集（CramersV-clustering OR |Δβ|-meanshift）· 全 6 樣本 × TP/FP（~30 萬位點）· post-hoc 不動源頭 · 2026-06-07 · §13 layer-A</div></header>

<div class="def">
  <h2>「有甲基意義」操作定義（4 要件，迭代調整內部門檻）</h2>
  <p>① <b>真實</b>：≥2/4 獨立 permutation 檢定 p≤0.05　② <b>夠強</b>：CramersV≥0.1 <b>OR</b> |Δβ|≥0.10（← 你選的聯集 (c)）　③ <b>乾淨</b>：無 dispersion-warn artifact + valid　④ <b>有檢力</b>：reads≥20</p>
  <p class="sub">Baseline I0（ISM 原生 Significant）<b>凍結</b>：源頭 significance_summary.csv 不動 → 隨時可回溯。</p>
</div>

<section class="card">
  <div class="ch"><span class="cn">迭代旅程（每步 diff + 判別比）</span></div>
  <table>
    <thead><tr><th>迭代</th><th class="num">TP 率<br>(數)</th><th class="num">FP 率<br>(數)</th><th class="num">TP/FP 比</th><th class="num">vs I0<br>+/-</th><th>branch<br>組成</th></tr></thead>
    <tbody>{rows}</tbody></table>
  <div class="figbox"><img src="{FIG}"></div>
</section>

<div class="find">
  <h2>🔴 關鍵發現 + 機制解釋（驗證假設 → 找到例外）</h2>
  <p>{esc(D["key_finding"])}</p>
  <p style="font-size:.85rem"><b>白話</b>：你假設「兩訊號都能找到有意義資訊」—— 機械上對（聯集多抓 ~85% 位點），但<b>加 Δβ branch 後 TP/FP 判別比從 3.7 崩到 ~1.0</b>（FP 升得比 TP 還多）→ Δβ 平均差多抓的是「基因組困難區」（LOH/低覆蓋，FP 集中處）非 somatic 真實訊號。**這是假設的例外，且與舊機制 (strong-ASM FP-enriched) 一致。**</p>
</div>

<h3 style="color:var(--accent)">每個迭代：理由 → 假設 → 結果 → 驗證/例外</h3>
{cards}

<div class="keep">
  <h2>★ 保留的最佳結果（best pick: {esc(best["id"])}）</h2>
  <p><b>用途</b>：{esc(best["purpose"])}</p>
  <p><b>理由</b>：{esc(best["reason"])}</p>
  <p><b>⚠ caveat</b>：{esc(best["caveat"])}</p>
</div>

<section class="card">
  <div class="ch"><span class="cn">結論（讓你理解 + 驗證）</span></div>
  <ul style="font-size:.88rem">
    <li><b>定義清楚</b>：meaningful = 真實+夠強(聯集)+乾淨+有檢力，門檻見上。</li>
    <li><b>過程可回溯</b>：I0 凍結；I1→I3 全 post-hoc，iteration_history.json 記錄每步 diff/假設/驗證/例外。</li>
    <li><b>保留最佳</b>：I2（聯集+corroboration+FDR）= 15,391 個乾淨 meaningful 位點當<b>表徵集</b>。</li>
    <li><b>誠實邊界</b>：此集<b>不是 TP/FP filter</b>（比≈1.0）；要判別用 I0 但仍弱。Δβ branch 的 FP-balanced 是真實例外，已機制解釋。</li>
    <li><b>每步例外檢查</b>：dispersion-warn 漏入 = 0（4 迭代皆過）。</li>
  </ul>
</section>

<footer>數據源：<span class="mono2">iteration_history.json</span>（← 6 樣本 × TP/FP significance_summary，源頭未改）· 腳本 76（驗證象限）→ 77（迭代）→ 78（HTML）· FDR=BH-q on GlobalP（min-p anti-conservative caveat）· §13 layer-A</footer>
</div></body></html>'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(HTML_DOC)
print(f"[78] wrote {OUT} ({len(HTML_DOC)//1024} KB)")
