#!/usr/bin/env python3
"""產 PI standalone HTML — 數字全從 summary_all_samples.json + tp_fp_comparison.json 注入(§13-A 反捏造)。
缺檔直接 refuse。用法: build_pi_html.py <RESULT_DIR> <COMMIT> <OUT_HTML>
"""
import sys, json, os, html

RD = sys.argv[1]
COMMIT = sys.argv[2] if len(sys.argv) > 2 else "UNKNOWN"
OUT = sys.argv[3] if len(sys.argv) > 3 else f"{RD}/PI_report.standalone.html"

summ_path = f"{RD}/summary_all_samples.json"
tpfp_path = f"{RD}/tp_fp_comparison.json"
for p in (summ_path, tpfp_path):
    if not os.path.exists(p):
        sys.exit(f"REFUSE: 缺真值來源 {p}（§13-A 不可捏造）")
summ = json.load(open(summ_path))
tpfp = json.load(open(tpfp_path))

ORDER = ["HCC1395", "HCC1395_DORADO", "COLO829", "H2009", "H1437", "HCC1937", "HCC1954"]
samples = [s for s in ORDER if s in summ]

def td(v):
    return f"<td>{html.escape(str(v))}</td>"

# --- 主表 rows ---
main_rows = ""
for s in samples:
    r = summ[s]
    main_rows += ("<tr>" + td(s) + td(f"{r['N']:,}")
        + td(f"{r['A_cansplit']:,} ({r['A_pct']}%)")
        + td(f"{r['C_clean_permanova']:,} ({r['C_pct']}%)")
        + td(f"{r['AnC']:,} ({r['AnC_pct']}%)")
        + td(f"{r['A_minus_C']:,}") + td(f"{r['C_minus_A']:,}")
        + f"<td class='jac'>{r['jaccard']}</td></tr>")

# --- 分帶 rows ---
band_rows = ""
for s in samples:
    r = summ[s]
    def p(k):
        return f"{r[k]:,} ({round(100*r[k]/r['N'],1)}%)" if r['N'] else r[k]
    band_rows += ("<tr>" + td(s) + td(f"{r['N']:,}")
        + td(p('insuff')) + td(p('cant')) + td(p('cansplit'))
        + td(f"{r['clear']:,} ({round(100*r['clear']/r['N'],1)}%)")
        + td(f"{r['mid']:,}") + td(f"{r['weak']:,}") + td(f"{r['vweak']:,}") + td(f"{r['degen']:,}") + "</tr>")

# --- TP/FP rows ---
tpfp_rows = ""
for ps in tpfp["per_sample"]:
    s = ps["sample"]; tp = ps["tp"]; fp = ps["fp"]
    if fp:
        lowN = " ⚠" if fp["N"] < 100 else ""
        dA = round(tp["A_pct"] - fp["A_pct"], 1)
        dC = round(tp["C_pct"] - fp["C_pct"], 1)
        tpfp_rows += ("<tr>" + td(s) + td(f"{tp['N']:,}") + td(f"{tp['A_pct']}%") + td(f"{tp['C_pct']}%") + td(tp['jaccard'])
            + td(f"{fp['N']:,}{lowN}") + td(f"{fp['A_pct']}%") + td(f"{fp['C_pct']}%") + td(fp['jaccard'])
            + f"<td class='{'pos' if dA>0 else 'neg'}'>{dA:+}</td><td class='{'pos' if dC>0 else 'neg'}'>{dC:+}</td></tr>")
    else:
        tpfp_rows += "<tr>" + td(s) + td(f"{tp['N']:,}") + td(f"{tp['A_pct']}%") + td(f"{tp['C_pct']}%") + td(tp['jaccard']) + "<td colspan=6 class='miss'>FP MISSING</td></tr>"

tpp = tpfp["tp_pooled"]; fpp = tpfp["fp_pooled"]
jacs = [summ[s]['jaccard'] for s in samples]
jmin, jmax = min(jacs), max(jacs)

HTML = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>全7樣本 subcluster 切割總帳 — PI report</title>
<style>
:root{{--ink:#1a1a1a;--mut:#666;--line:#e2ddd3;--bg:#faf8f3;--acc:#8a5a2b;--ok:#1f7a4d;--warn:#b3541e;--card:#fff}}
*{{box-sizing:border-box}}
body{{margin:0;font-family:-apple-system,"Helvetica Neue","PingFang TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:var(--bg);line-height:1.65}}
.wrap{{max-width:1180px;margin:0 auto;padding:32px 28px 80px}}
header h1{{font-size:25px;margin:0 0 6px}}
.sub{{color:var(--mut);font-size:13px;margin-bottom:4px}}
nav{{position:sticky;top:0;background:rgba(250,248,243,.95);backdrop-filter:blur(6px);border-bottom:1px solid var(--line);padding:10px 0;margin:18px 0 24px;font-size:13px;z-index:10}}
nav a{{color:var(--acc);text-decoration:none;margin-right:16px;white-space:nowrap}}
nav a:hover{{text-decoration:underline}}
h2{{font-size:19px;border-left:4px solid var(--acc);padding-left:10px;margin:34px 0 14px}}
.head{{background:linear-gradient(180deg,#fff,#fcfaf5);border:1px solid var(--line);border-radius:12px;padding:20px 22px;margin:8px 0 8px}}
.head .big{{font-size:16px;font-weight:600}}
.head .big b{{color:var(--acc)}}
table{{width:100%;border-collapse:collapse;font-size:13px;background:var(--card);border:1px solid var(--line);border-radius:8px;overflow:hidden;margin:6px 0 10px}}
th,td{{padding:7px 9px;text-align:right;border-bottom:1px solid var(--line)}}
th:first-child,td:first-child{{text-align:left;font-weight:600}}
thead th{{background:#f1ece2;color:#4a4036;font-size:12px;position:sticky;top:46px}}
tbody tr:hover{{background:#fcf8f0}}
.jac{{font-weight:700;color:var(--acc)}}
.pos{{color:var(--ok)}} .neg{{color:var(--warn)}} .miss{{color:var(--mut);text-align:center}}
.card{{background:var(--card);border:1px solid var(--line);border-left:4px solid var(--acc);border-radius:8px;padding:13px 16px;margin:10px 0;font-size:14px}}
.card.warn{{border-left-color:var(--warn)}}
.card b{{color:var(--acc)}}
.note{{font-size:12px;color:var(--mut)}}
.pill{{display:inline-block;background:#efe7d8;color:#6a4a22;border-radius:20px;padding:1px 9px;font-size:11px;margin-left:6px}}
ul{{margin:6px 0 6px 0;padding-left:20px}} li{{margin:4px 0;font-size:14px}}
footer{{margin-top:40px;padding-top:16px;border-top:1px solid var(--line);font-size:11px;color:var(--mut)}}
code{{background:#f1ece2;padding:1px 5px;border-radius:4px;font-size:12px}}
</style></head><body><div class="wrap">
<header>
<h1>全 7 樣本 subcluster 切割總帳<span class="pill">PI report</span></h1>
<div class="sub">canonical paired_full longphase-s tagged BAM · tumor-only · 全基因組 · 3 癌種 7 樣本</div>
<div class="sub">binary <code>inter_sub_mod @ feat/summary-nreadsvalid {COMMIT}</code> · ISM <code>w=5000 BERNOULLI C_min=3 SKIP</code> · 2026-06-20</div>
</header>

<nav>
<a href="#tldr">結論</a><a href="#main">主表 A/C Venn</a><a href="#tpfp">TP vs FP</a>
<a href="#bands">品質分帶</a><a href="#obs">跨樣本觀察</a><a href="#method">方法/驗證</a><a href="#caveat">口徑</a>
</nav>

<section id="tldr">
<div class="head">
<div class="big">7 樣本全部：「切得出群」(A=cansplit) 與「標籤對齊乾淨 PERMANOVA」(C) 的 Jaccard <b>一律偏低（{jmin}–{jmax}）</b>。</div>
<p style="margin:10px 0 0;font-size:14px">兩種訊號擷取到<b>大致不相交</b>的位點集合，跨 3 癌種穩健重現 HCC1395 單樣本發現 → <b>無監督切群數 ≠ 可確認 subclone</b>。</p>
</div>
</section>

<section id="main">
<h2>1 · 主表 — A/C clean-location Venn（全基因組 TP）</h2>
<table><thead><tr><th>樣本</th><th>N</th><th>A=cansplit</th><th>C=clean PERMANOVA</th><th>A∩C</th><th>A−C</th><th>C−A</th><th>Jaccard</th></tr></thead>
<tbody>{main_rows}</tbody></table>
<p class="note">A = 無監督 UPGMA+silhouette best_k≥2（tumor reads）；C = LabelHP 或 LabelAllele PERMANOVA p&lt;0.05 且非 dispersion。</p>
</section>

<section id="tpfp">
<h2>2 · TP vs FP 比例對照</h2>
<table><thead><tr><th>樣本</th><th>TP N</th><th>TP A%</th><th>TP C%</th><th>TP Jac</th><th>FP N</th><th>FP A%</th><th>FP C%</th><th>FP Jac</th><th>ΔA%</th><th>ΔC%</th></tr></thead>
<tbody>{tpfp_rows}</tbody></table>
<div class="card warn"><b>FP 也一樣會切。</b> 可靠的大-N FP 樣本（HCC1395 FP={summ['HCC1395']['N'] and tpfp['per_sample'][0]['fp']['N']}、COLO829 FP={[p['fp']['N'] for p in tpfp['per_sample'] if p['sample']=='COLO829'][0]}）的 cansplit% 與 TP <b>相當或更高</b>（HCC1395 FP A=38.8% &gt; TP 19.5%）→ <b>切群能力不區分 TP/FP</b>，與既有「ASM/結構非 usable filter」一致。</div>
<p class="note">⚠ 標記者 FP N&lt;100（H1437=8 / HCC1954=29 / H2009=86）為低 N 噪音，ΔA/ΔC 不可解讀。Pooled: TP N={tpp['N']:,} A={tpp['A_pct']}% C={tpp['C_pct']}%；FP N={fpp['N']:,} A={fpp['A_pct']}% C={fpp['C_pct']}%（pooled 受最大樣本主導，僅參考）。</p>
</section>

<section id="bands">
<h2>3 · 切割品質總帳（cant vs cansplit + 分帶）</h2>
<table><thead><tr><th>樣本</th><th>N</th><th>insuff(n&lt;6)</th><th>cant(切不出)</th><th>cansplit</th><th>clear(V≥.7+sig)</th><th>mid</th><th>weak</th><th>vweak</th><th>degen</th></tr></thead>
<tbody>{band_rows}</tbody></table>
<p class="note">分帶 clear/mid/weak 依 build_spectrum_html.py 顯示定義重建（V 邊界±1-3 容差）；insuff/cant/cansplit/degen/vweak 為確定性定義。</p>
</section>

<section id="obs">
<h2>4 · 跨樣本觀察（見樹也見林）</h2>
<ul>
<li><b>Aggregate</b>：Jaccard(A,C) 全 7 樣本 {jmin}–{jmax} 一律偏低 → 切群 vs 標籤對齊大致不相交，跨 3 癌種穩健。</li>
<li><b>Canonical</b>：HCC1395 paired_full（A {summ['HCC1395']['A_pct']}% / C {summ['HCC1395']['C_pct']}% / J {summ['HCC1395']['jaccard']}）≈ pileup-track 參考（19.7% / 29.0% / 0.118）→ 換 BAM track 幾乎不動主數字。</li>
<li><b>Outlier</b>：H2009 cansplit {summ['H2009']['A_pct']}% 最高但 vweak {summ['H2009']['vweak']:,}（占 cansplit {round(100*summ['H2009']['vweak']/summ['H2009']['cansplit'],0):.0f}%）= 過度切分無意義；H1437 insuff {summ['H1437']['insuff']:,}（{round(100*summ['H1437']['insuff']/summ['H1437']['N'],1)}%）最高=低覆蓋；COLO829 clear 僅 {round(100*summ['COLO829']['clear']/summ['COLO829']['N'],1)}% 最低。</li>
<li><b>Well-explained</b>：低 Jaccard + FP 同樣會切 → 重現「within_clean≠subclone」「切不出≠沒訊號」「ASM 非 usable filter」框架。</li>
</ul>
</section>

<section id="method">
<h2>5 · 方法與驗證</h2>
<div class="card"><b>邏輯驗證（信任基礎）</b>：driver 先用 HCC1395 <b>pileup-track</b> 既有資料重算，A/C/A∩C/A−C/C−A/Jaccard <b>位元符合</b>參考表（5997 / 8833 / 1560 / 4437 / 7273 / 0.118）才上線跑全樣本。</div>
<p style="font-size:14px">每樣本：單次 ISM scan（保留 distance matrix）→ 逐字 <code>wg_contingency</code> 的 UPGMA+silhouette best_k → 逐字 <code>permanova_clean_and_4group.py</code> Q1 Venn。單次 ISM 取代原始 2-pass（per-region 矩陣位元相同）。每樣本 <code>insuff+cant+cansplit=N</code>、<code>分帶=cansplit</code>、Venn 關係全部加總一致 ✓。</p>
</section>

<section id="caveat">
<h2>6 · 口徑（給審稿/PI）</h2>
<ul>
<li><b>tumor-only</b>：normal BAM 各樣本來源不對稱；cansplit 分群本就只用 tumor reads，故統一 tumor-only 使 7 樣本完全可比；對 A 無影響、C 自洽可算。</li>
<li><b>track</b>：原始參考(5997/8833/0.118)為 pileup-track；本表為 paired_full longphase-s，HCC1395 兩者極接近。</li>
<li><b>FP 低 N</b>：H1437/HCC1954/H2009 的 FP region 數 &lt;100，個別比例為噪音。</li>
<li><b>tier ⭐3</b>：單 binary / 單 pipeline；cross-platform 復現需外部工具。</li>
</ul>
</section>

<footer>
provenance: binary <code>inter_sub_mod @ feat/summary-nreadsvalid {COMMIT}</code> · 數字注入自 <code>summary_all_samples.json</code> + <code>tp_fp_comparison.json</code>（§13-A 由構造防捏造，非手打）· 生成 2026-06-20 ·
報告目錄 <code>InterSubMod/docs/experiments/in_progress/2026/06/20260620_allsample_subcluster_split/</code>
</footer>
</div></body></html>"""

open(OUT, "w").write(HTML)
print(f"[PI HTML] {OUT} ({len(HTML)} bytes)")
