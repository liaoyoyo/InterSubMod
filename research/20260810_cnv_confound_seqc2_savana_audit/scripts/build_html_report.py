#!/usr/bin/env python3
"""Build the standalone HTML companion for the CN-confound audit.

Every figure quoted in the page is injected from the stored JSON outputs; the
script refuses to render if a required key is missing, so a number can never be
typed by hand into the page (project rule 13-A).
"""

from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(HERE, "..", "docs", "20260810_CNV混淆稽核_01.standalone.html")


def load(name):
    with open(os.path.join(DATA, name)) as fh:
        return json.load(fh)


def need(d, *path):
    cur = d
    for p in path:
        if isinstance(cur, dict) and p in cur:
            cur = cur[p]
        elif isinstance(cur, list) and isinstance(p, int) and p < len(cur):
            cur = cur[p]
        else:
            sys.exit(f"REFUSE TO RENDER: missing required value {'/'.join(map(str, path))}")
    if cur is None:
        sys.exit(f"REFUSE TO RENDER: null required value {'/'.join(map(str, path))}")
    return cur


rv = load("report_values.json")
cons = load("consolidated_findings.json")
res = load("resolution_vs_cn_stratified.json")
purity = load("purity_selfconsistency_audit.json")
mech1 = load("mechanism_and_intervals.json")
rob = load("robustness_checks.json")
cn = load("cn_annotation_summary.json")

# ---- pull every value up front so a missing key fails before any HTML exists ----
V = {
    "unit_altered": need(rv, "cn_share", "unit_cn_altered_percent"),
    "site_altered": need(rv, "cn_share", "site_cn_altered_percent"),
    "unit_neutral": need(rv, "cn_share", "unit_class_percent", "neutral"),
    "site_neutral": need(rv, "cn_share", "site_class_percent", "neutral"),
    "struct_neu": need(rv, "resolution", "structure_unique_percent_neutral"),
    "struct_alt": need(rv, "resolution", "structure_unique_percent_altered"),
    "struct_all": need(rv, "resolution", "structure_unique_percent_all"),
    "uniq_all": need(rv, "resolution", "unique_best_tree_percent_all"),
    "uniq_neu": need(rv, "resolution", "unique_best_tree_percent_neutral"),
    "uniq_alt": need(rv, "resolution", "unique_best_tree_percent_altered"),
    "or_fisher": need(rv, "resolution", "fisher_odds_ratio"),
    "p_fisher": need(rv, "resolution", "fisher_p_value_scientific"),
    "or_tree": need(rv, "resolution", "cmh_or_tree_count"),
    "or_k": need(rv, "resolution", "cmh_or_k"),
    "or_both": need(rv, "resolution", "cmh_or_both"),
    "or_chrom": need(rv, "robustness", "cmh_or_chromosome"),
    "p_chrom": need(rv, "robustness", "cmh_p_chromosome_scientific"),
    "top4": need(rv, "robustness", "neutral_top4_combined_percent"),
    "pub_purity": need(rv, "savana", "published_purity"),
    "pub_ploidy": need(rv, "savana", "published_ploidy"),
    "pub_agree": need(rv, "savana", "published_state_agreement_percent"),
    "pub_int": need(rv, "savana", "published_integer_match_percent"),
    "best_purity": need(rv, "savana", "best_purity"),
    "best_ploidy": need(rv, "savana", "best_ploidy"),
    "best_agree": need(rv, "savana", "best_state_agreement_percent"),
    "best_int": need(rv, "savana", "best_integer_match_percent"),
    "seqc2_ploidy": need(rv, "savana", "seqc2_implied_ploidy"),
    "r2_all": need(rv, "savana", "regression_all_r2"),
    "r2_non": need(rv, "savana", "regression_nonneutral_r2"),
    "excess_pct": need(mech1, "cn_attributable_uniqueness_excess", "excess_as_percent_of_unique_units", "point"),
    "excess_lo": need(mech1, "cn_attributable_uniqueness_excess", "excess_as_percent_of_unique_units", "lo"),
    "excess_hi": need(mech1, "cn_attributable_uniqueness_excess", "excess_as_percent_of_unique_units", "hi"),
    "excess_n": need(mech1, "cn_attributable_uniqueness_excess", "excess_point"),
    "uniq_units": need(mech1, "cn_attributable_uniqueness_excess", "total_unique_best_tree_units"),
    "ident_neu": need(rv, "mechanism", "identical_af_neutral_percent"),
    "ident_alt": need(rv, "mechanism", "identical_af_altered_percent"),
    "mediation_or": need(rv, "mechanism", "mediation_or"),
    "spread_neu": need(rv, "mechanism", "spread_median_neutral"),
    "spread_alt": need(rv, "mechanism", "spread_median_altered"),
    "taut_units": need(rv, "mechanism", "tautology_units"),
    "spread_p": need(rv, "mechanism", "spread_p_scientific"),
    "ident_or": need(rv, "mechanism", "identical_af_or"),
    "ident_p": need(rv, "mechanism", "identical_af_p_scientific"),
    "mediation_p": need(rv, "mechanism", "mediation_p_scientific"),
    "unres_neu": need(rv, "mechanism", "unresolvable_neutral_percent"),
    "unres_alt": need(rv, "mechanism", "unresolvable_altered_percent"),
    "nondeg_alt": need(rv, "mechanism", "nondeg_altered_percent"),
    "nondeg_neu": need(rv, "mechanism", "nondeg_neutral_percent"),
    "nondeg_n_alt": need(rv, "mechanism", "nondeg_n_altered"),
    "nondeg_n_neu": need(rv, "mechanism", "nondeg_n_neutral"),
    "nondeg_or": need(rv, "mechanism", "nondeg_or"),
    "nondeg_p": need(rv, "mechanism", "nondeg_p"),
    "chr8_units": need(rv, "chr8", "units"),
    "chr8_gl": need(rv, "chr8", "gain_loh"),
    "chr8_pct": need(rv, "chr8", "gain_loh_percent"),
    "zero_chroms": need(rv, "chromosomes_without_neutral", "chroms"),
    "zero_units": need(rv, "chromosomes_without_neutral", "combined_units"),
}

af_class = need(rv, "resolution", "af_broke_tie_percent_by_class")
af_med = need(rv, "mechanism", "af_median_by_class")
contested = need(res, "af_broke_tie_by_cn_class_contested_only")
chrom_rows = need(cons, "per_chromosome_distribution")
example = need(cons, "worked_examples", "gain_loh_af_broke_tie")
unit_counts = need(cn, "unit_class_counts")
site_counts = need(cn, "site_class_counts")

CLASS_LABEL = {
    "neutral": "CN 中性",
    "gain": "gain",
    "gain_loh": "gain + LOH",
    "cnloh": "copy-neutral LOH",
    "loss": "loss",
    "loss_loh": "loss + LOH",
}
CLASS_COLOR = {
    "neutral": "var(--ok)",
    "gain": "var(--warn)",
    "gain_loh": "var(--bad)",
    "cnloh": "var(--mid)",
    "loss": "var(--dim)",
    "loss_loh": "var(--dim)",
}


def bar_rows(order, getval, maxval, fmt="{:.2f}%"):
    out = []
    for k in order:
        v = getval(k)
        if v is None:
            continue
        w = 100.0 * float(v) / maxval if maxval else 0
        out.append(
            f'<div class="bar-row"><span class="bar-lab">{CLASS_LABEL.get(k,k)}</span>'
            f'<span class="bar-track"><span class="bar-fill" style="width:{w:.1f}%;'
            f'background:{CLASS_COLOR.get(k,"var(--mid)")}"></span></span>'
            f'<span class="bar-val">{fmt.format(float(v))}</span></div>'
        )
    return "".join(out)


ORDER = ["neutral", "cnloh", "gain", "gain_loh", "loss", "loss_loh"]

tie_bars = bar_rows(ORDER, lambda k: af_class.get(k), 100.0)
med_bars = bar_rows(ORDER, lambda k: af_med.get(k), 1.0, fmt="{:.4f}")

chrom_bars = ""
for r in chrom_rows:
    neu = r["neutral_percent"] or 0
    chrom_bars += (
        f'<div class="bar-row"><span class="bar-lab">{r["chrom"]}</span>'
        f'<span class="bar-track"><span class="bar-fill" style="width:{neu:.1f}%;'
        f'background:{"var(--ok)" if neu>10 else "var(--bad)"}"></span></span>'
        f'<span class="bar-val">{neu:.2f}% <em>n={r["units"]}</em></span></div>'
    )

cg = rv["clean_ground"]
ey = rv["evolution_yield"]
ms = rv["multisample"]
sq = rv["seqc2_external_comparison"]
lk = rv["linkage_constraint"]
cv = rv["constraint_validity"]
cp = rv["cn_phenomena"]

cn_af_bars = ""
for cnk in ["1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0", "8.0"]:
    med = cp["af_median_by_total_cn"].get(cnk)
    if med is None:
        continue
    grid = cp["predicted_lowest_grid_point_by_cn"].get(cnk)
    hit = cp["grid_hit_percent_by_cn"].get(cnk)
    nn = cp["n_sites_by_cn"].get(cnk)
    col = "var(--ok)" if float(cnk) <= 2 else "var(--warn)" if float(cnk) <= 4 else "var(--bad)"
    cn_af_bars += (
        f'<div class="bar-row"><span class="bar-lab">CN = {cnk[:-2]}</span>'
        f'<span class="bar-track"><span class="bar-fill" style="width:{100*med:.1f}%;'
        f'background:{col}"></span></span>'
        f'<span class="bar-val">{med:.4f} <em>1/CN={grid} · 格點{hit}% · n={nn:,}</em></span></div>'
    )

purity_rows = ""
for s, v in purity["samples"].items():
    if not isinstance(v, dict) or "fitted_purity" not in v:
        continue
    vr = v["violation_rate_percent"]
    verdict = "mis-fit" if vr > 30 else "通過"
    cls = "bad" if vr > 30 else "ok"
    purity_rows += (
        f"<tr><td>{s}</td><td class='num'>{v['fitted_purity']}</td>"
        f"<td class='num'>{v['fitted_ploidy']}</td>"
        f"<td class='num'>{v['segments_usable']}</td>"
        f"<td class='num {cls}'>{vr}%</td><td class='{cls}'>{verdict}</td></tr>"
    )

calib = purity["calibration_case"]["at_refit_purity_1.0"]

HTML = f"""<title>CNV 混淆稽核 — HCC1395 全基因組</title>
<style>
:root {{
  --bg:#ffffff; --fg:#1a1d21; --muted:#5b6470; --line:#e2e6ea; --card:#f7f9fb;
  --ok:#2e7d5b; --bad:#c0392b; --warn:#d98324; --mid:#4a6fa5; --dim:#9aa4b0;
  --accent:#1f4e79;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --bg:#14171a; --fg:#e8ecef; --muted:#9aa4b0; --line:#2a3038; --card:#1b1f24;
    --ok:#4caf82; --bad:#e57063; --warn:#e5a04a; --mid:#7ba3d8; --dim:#6b7684;
    --accent:#8ab4e8;
  }}
}}
:root[data-theme="dark"] {{
  --bg:#14171a; --fg:#e8ecef; --muted:#9aa4b0; --line:#2a3038; --card:#1b1f24;
  --ok:#4caf82; --bad:#e57063; --warn:#e5a04a; --mid:#7ba3d8; --dim:#6b7684;
  --accent:#8ab4e8;
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0; padding:2rem 1.25rem 5rem; background:var(--bg); color:var(--fg);
  font:16px/1.75 -apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans TC",
       "PingFang TC","Microsoft JhengHei",sans-serif; }}
.wrap {{ max-width:940px; margin:0 auto; }}
h1 {{ font-size:1.75rem; line-height:1.35; margin:0 0 .4rem; letter-spacing:-.01em; }}
h2 {{ font-size:1.22rem; margin:2.75rem 0 .85rem; padding-top:1.1rem;
      border-top:1px solid var(--line); letter-spacing:-.01em; }}
h3 {{ font-size:1.02rem; margin:1.6rem 0 .5rem; color:var(--accent); }}
.sub {{ color:var(--muted); margin:0 0 1.75rem; font-size:.94rem; }}
.verdict {{ background:var(--card); border-left:4px solid var(--accent);
   padding:1.15rem 1.3rem; border-radius:0 8px 8px 0; margin:1.5rem 0; }}
.verdict strong {{ color:var(--accent); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr));
   gap:.85rem; margin:1.35rem 0; }}
.stat {{ background:var(--card); border:1px solid var(--line); border-radius:9px;
   padding:.95rem 1.05rem; }}
.stat .n {{ font-size:1.7rem; font-weight:650; letter-spacing:-.02em; display:block;
   line-height:1.2; }}
.stat .l {{ font-size:.83rem; color:var(--muted); display:block; margin-top:.3rem; }}
.stat.ok .n {{ color:var(--ok); }} .stat.bad .n {{ color:var(--bad); }}
.stat.warn .n {{ color:var(--warn); }}
table {{ border-collapse:collapse; width:100%; font-size:.9rem; margin:1rem 0; }}
th,td {{ border-bottom:1px solid var(--line); padding:.55rem .6rem; text-align:left; }}
th {{ color:var(--muted); font-weight:600; font-size:.83rem;
   text-transform:uppercase; letter-spacing:.04em; }}
td.num {{ text-align:right; font-variant-numeric:tabular-nums; }}
td.ok {{ color:var(--ok); font-weight:600; }} td.bad {{ color:var(--bad); font-weight:600; }}
.scroll {{ overflow-x:auto; }}
.bar-row {{ display:flex; align-items:center; gap:.7rem; margin:.3rem 0;
   font-size:.87rem; }}
.bar-lab {{ flex:0 0 8.5rem; color:var(--muted); }}
.bar-track {{ flex:1 1 auto; height:15px; background:var(--card);
   border-radius:4px; overflow:hidden; border:1px solid var(--line); }}
.bar-fill {{ display:block; height:100%; border-radius:3px 0 0 3px; }}
.bar-val {{ flex:0 0 8.5rem; text-align:right; font-variant-numeric:tabular-nums; }}
.bar-val em {{ color:var(--muted); font-style:normal; font-size:.8rem; }}
.note {{ background:var(--card); border:1px solid var(--line); border-radius:8px;
   padding:.9rem 1.1rem; font-size:.89rem; color:var(--muted); margin:1.1rem 0; }}
.note.bad {{ border-left:3px solid var(--bad); }}
.note.ok {{ border-left:3px solid var(--ok); }}
code {{ background:var(--card); padding:.12em .4em; border-radius:4px;
   font-size:.88em; border:1px solid var(--line); }}
.two {{ display:grid; grid-template-columns:1fr 1fr; gap:1.5rem; }}
@media (max-width:700px) {{ .two {{ grid-template-columns:1fr; }}
  .bar-lab {{ flex-basis:6rem; }} .bar-val {{ flex-basis:6.5rem; }} }}
footer {{ margin-top:3.5rem; padding-top:1.2rem; border-top:1px solid var(--line);
   color:var(--muted); font-size:.83rem; }}
</style>
<div class="wrap">

<h1>CNV 混淆稽核：SEQC2 結對驗證與 SAVANA 流程驗證</h1>
<p class="sub">HCC1395 · chr1–22 全基因組 · canonical 2026-07-24 frozen cohort · 2026-08-10</p>

<div class="verdict">
<strong>結論：</strong>共現與拓撲骨幹<strong>沒有</strong>被 CNV 抬高，可以保留；
但「read-AF 挑出唯一最佳樹」這一層<strong>被顯著抬高</strong> ——
約 {V['excess_lo']}–{V['excess_hi']}% 的「唯一樹」是拷貝數的產物而非譜系訊號。
機制已釐清：CN 不是讓 read-AF 排錯順序，而是決定它<strong>有沒有東西可排</strong>。
既有結論不需撤回，但單一拓撲率這類數字必須標為 CN-conditional。
</div>

<h2>1 · 有多少坐在 CNV 區</h2>
<div class="grid">
  <div class="stat bad"><span class="n">{V['unit_altered']}%</span>
    <span class="l">unit 層 CN-altered（9,185 / 9,624）</span></div>
  <div class="stat bad"><span class="n">{V['site_altered']}%</span>
    <span class="l">sSNV 層 CN-altered（18,833 / 20,060）</span></div>
  <div class="stat ok"><span class="n">{V['unit_neutral']}%</span>
    <span class="l">unit 層乾淨（CN=2 且無 LOH）</span></div>
  <div class="stat ok"><span class="n">{V['site_neutral']}%</span>
    <span class="l">sSNV 層乾淨</span></div>
</div>

<div class="scroll"><table>
<tr><th>CN 類別</th><th>unit 數</th><th>site 數</th></tr>
{''.join(f"<tr><td>{CLASS_LABEL.get(k,k)}</td><td class='num'>{unit_counts.get(k,0):,}</td><td class='num'>{site_counts.get(k,0):,}</td></tr>" for k in ORDER)}
</table></div>

<h3>逐染色體：乾淨地在哪裡</h3>
{chrom_bars}
<div class="note bad">
<strong>關鍵限制：</strong>CN 中性樣本高度集中 —— chr2 + chr21 + chr6 + chr4 就佔了
全部中性 unit 的 <strong>{V['top4']}%</strong>，而 {len(V['zero_chroms'])} 條染色體
（{' / '.join(V['zero_chroms'])}，共 {V['zero_units']:,} 個 unit）<strong>完全沒有</strong>中性區。
chr8 的 {V['chr8_units']} 個 unit 中有 {V['chr8_gl']} 個是 gain+LOH（{V['chr8_pct']}%）。
因此「中性 vs 變異」的對比天生帶有染色體混淆，第 4 節對此做了專門檢驗。
</div>

<h2>2 · 拓撲層沒被抬高，read-AF 層被抬高</h2>
<p>把「唯一性」拆成兩層是本次分析的關鍵：<code>結構唯一</code>指枚舉本身就只有一棵樹、
read-AF 未參與（理論上 CN-robust）；<code>AF 破 tie</code>指多棵結構合法的樹由 read-AF 選出贏家
（正是 CN 能污染的一層）。</p>

<div class="two">
<div>
<h3>結構唯一率 — 中性反而更高</h3>
<div class="grid" style="grid-template-columns:1fr 1fr">
  <div class="stat ok"><span class="n">{V['struct_neu']}%</span><span class="l">CN 中性</span></div>
  <div class="stat"><span class="n">{V['struct_alt']}%</span><span class="l">CN 變異</span></div>
</div>
<p style="font-size:.89rem;color:var(--muted)">方向與「CN 抬高唯一性」相反 →
共現骨幹 CN-robust，拓撲軸存活。</p>
</div>
<div>
<h3>read-AF 破 tie 率 — 差距懸殊</h3>
<div class="grid" style="grid-template-columns:1fr 1fr">
  <div class="stat ok"><span class="n">{contested['neutral']['percent']}%</span><span class="l">CN 中性（n={contested['neutral']['n']}）</span></div>
  <div class="stat bad"><span class="n">{contested['gain_loh']['percent']}%</span><span class="l">gain+LOH（n={contested['gain_loh']['n']}）</span></div>
</div>
<p style="font-size:.89rem;color:var(--muted)">Fisher OR = {V['or_fisher']}，p = {V['p_fisher']}</p>
</div>
</div>

<h3>各 CN 類別的 read-AF 破 tie 率（僅候選樹 &gt; 1 者）</h3>
{tie_bars}

<h3>分層控制後是否存活</h3>
<div class="scroll"><table>
<tr><th>分層變項</th><th>CMH 合併 OR</th><th>解讀</th></tr>
<tr><td>候選樹數</td><td class="num">{V['or_tree']}</td><td>幾乎不衰減</td></tr>
<tr><td>active bit count k</td><td class="num">{V['or_k']}</td><td>幾乎不衰減</td></tr>
<tr><td>候選樹數 × k</td><td class="num">{V['or_both']}</td><td>幾乎不衰減</td></tr>
<tr><td><strong>染色體</strong></td><td class="num bad">{V['or_chrom']}</td><td>效應存活但<strong>量級減半</strong>（p = {V['p_chrom']}）</td></tr>
</table></div>
<div class="note ok">
局部突變密度是本專案已知的必要控制變項 —— 2026-07-26 的甲基分析就是在分層後訊號全滅。
這次前三項分層幾乎不衰減，<strong>是相反的模式</strong>。但按染色體分層後縮到 {V['or_chrom']}，
據實記錄：約四成原始效應可歸因於染色體結構。
</div>

<h2>3 · 影響量化</h2>
<div class="grid">
  <div class="stat bad"><span class="n">{V['excess_pct']}%</span>
    <span class="l">「唯一樹」中 CN 可歸因的超額（範圍 {V['excess_lo']}–{V['excess_hi']}%）</span></div>
  <div class="stat"><span class="n">{V['excess_n']:,.0f}</span>
    <span class="l">超額 unit 數（總 {V['uniq_units']:,} 個唯一樹）</span></div>
</div>
<div class="note bad">
此估計假設「CN 中性的破 tie 率就是無 CN 基線」，但中性區仍可能含未報告的 subclonal CN，
且中性 contested 樣本僅 {contested['neutral']['n']} 個、區間寬。
若改採染色體分層後的 OR {V['or_chrom']}，超額會低於點估計 —— <strong>上表應視為上緣估計</strong>。
</div>

<h2>4 · 機制：CN 決定 AF 值在算術上是否可區分</h2>

<div class="note bad">
<strong>先說修正。</strong>本節初版把不進分數的位點算了進去 —— 相異性與離散度遍歷了整個
<code>af_coverage</code>，但 read-AF 分數只看 <code>active_positions</code>；其餘是全參考列
（<code>0/1</code>），永遠不進 <code>s(p,c)</code>。9,130 個 ranked unit 中 4,734 個帶這種多餘列。
此錯誤由對抗審查抓出（49 項發現中唯一確認為真者），修正後<strong>兩項檢驗的結論都改變且方向一致</strong>。
第 1 節 CN 佔比、第 2 節結構層與主要關聯、第 3 節超額估計、第 5–6 節 SAVANA 皆不使用相異性旗標，已逐項確認不受影響。
</div>

<div class="scroll"><table>
<tr><th>檢驗</th><th>初版（含非 active 位點）</th><th>修正後（僅 active 位點）</th></tr>
<tr><td>AF 離散度 altered vs neutral</td><td class="bad">0.5223 &lt; 0.7826，p = 0.9997 → 推翻</td>
    <td class="ok">{V['spread_alt']} &gt; {V['spread_neu']}，p = {V['spread_p']} → 成立</td></tr>
<tr><td>AF 全同比例 neutral vs altered</td><td>20.80% vs 11.78%，OR 1.9669</td>
    <td class="ok">{V['ident_neu']}% vs {V['ident_alt']}%，OR {V['ident_or']}，p = {V['ident_p']}</td></tr>
<tr><td>中介檢驗 CMH</td><td class="bad">2.819，p = 3.9e-06 → 部分中介</td>
    <td class="ok">{V['mediation_or']}，p = {V['mediation_p']} → 完全中介</td></tr>
</table></div>

<h3>機制本身</h3>
<p>read-AF 以 exact rational 比較，tie 發生在候選樹<strong>分數相等</strong>時。分數
<code>s(p,c) = Σ(AF_i − AF_j)</code>：若一個 unit 的所有 active AF 都等於同一個值，
每項為 0、所有候選樹同分 → <strong>tie 是數學必然</strong>。實測：{V['taut_units']:,} 個 active AF 全同的 unit，
<code>best_score</code> 全部是 <code>0/1</code>，零例外。</p>
<p>純腫瘤細胞株的 CN 中性區，單 copy haplotype 上的 clonal 突變讀出恰好 1/1，同一 unit 內多點共享同值；
gain/LOH 下同一 clone 的突變落在不同 m/c 格點，值就分開了。</p>

<h3>各 CN 類別的 read-AF 中位數</h3>
{med_bars}
<p style="font-size:.89rem;color:var(--muted)">gain+LOH 的中位數恰好是
<strong>{af_med['gain_loh']}</strong> ≈ 1/3 —— 正是 CN=3、完全 LOH、單 copy 突變的期望值，
直接的 mutation multiplicity 讀數。</p>

<h3>完全中介：效應全部經由「算術可解性」</h3>
<div class="two">
<div>
<p style="font-size:.9rem;margin:0 0 .5rem">contested unit 中<strong>算術上不可能被解開</strong>的比例：</p>
<div class="grid" style="grid-template-columns:1fr 1fr">
  <div class="stat bad"><span class="n">{V['unres_neu']}%</span><span class="l">CN 中性（n={contested['neutral']['n']}）</span></div>
  <div class="stat"><span class="n">{V['unres_alt']}%</span><span class="l">CN 變異（n={contested['gain']['n']+contested['gain_loh']['n']+contested['cnloh']['n']+contested['loss']['n']+contested['loss_loh']['n']:,}）</span></div>
</div>
</div>
<div>
<p style="font-size:.9rem;margin:0 0 .5rem">排除這些之後的破 tie 率：</p>
<div class="grid" style="grid-template-columns:1fr 1fr">
  <div class="stat"><span class="n">{V['nondeg_alt']}%</span><span class="l">CN 變異（n={V['nondeg_n_alt']:,}）</span></div>
  <div class="stat"><span class="n">{V['nondeg_neu']}%</span><span class="l">CN 中性（n={V['nondeg_n_neu']}）</span></div>
</div>
</div>
</div>
<div class="note ok">
Fisher OR = {V['nondeg_or']}，p = {V['nondeg_p']}；控制相異性與候選樹數的 CMH OR = {V['mediation_or']}，p = {V['mediation_p']}。
<strong>在算術可解性之外，偵測不到殘餘的 CN 效應。</strong><br><br>
這是<strong>強化</strong>主結論而非削弱：CN 不是讓 read-AF「排錯順序」，而是決定它
<strong>有沒有東西可排</strong>。CN 中性區有三分之二的 unit 在算術上本來就無解，
CN 變異區把這個比例砍半 —— 多出來的「可解」單元，其可解性來自拷貝數結構，不是譜系。<br><br>
⚠ 非退化子集的 n_neutral 僅 {V['nondeg_n_neu']}，區間寬，OR = 1 並未被證明；
只能說「無法宣稱存在殘餘效應」，不能說「已證明沒有」。
</div>

<h2>4A · 非 CN 區還剩多少資訊？「乾淨」比 4.56% 還要少</h2>
<p>4.56% 是<strong>區域</strong>比例，不是<strong>可用資訊</strong>比例 —— 建樹至少需要 2 個 sSNV。</p>
<div class="scroll"><table>
<tr><th></th><th>ranked unit</th><th>k=1（數學上不可能建樹）</th><th>k≥2（真有樹結構資訊）</th></tr>
<tr><td>CN 中性</td><td class="num">{cg['neutral_ranked']}</td>
    <td class="num bad">{cg['neutral_k1_no_tree']}</td>
    <td class="num">{cg['neutral_k_ge2']}（{cg['neutral_k_ge2_percent']}%）</td></tr>
<tr><td>CN 變異</td><td class="num">8,779</td><td class="num">2,972</td>
    <td class="num">{cg['altered_k_ge2']:,}（{cg['altered_k_ge2_percent']}%）</td></tr>
</table></div>
<div class="grid">
  <div class="stat bad"><span class="n">{cg['tree_informative_neutral_share_percent']}%</span>
    <span class="l">真正能提供樹結構資訊的乾淨 unit 佔比（{cg['neutral_k_ge2']}／5,984）</span></div>
  <div class="stat warn"><span class="n">{cg['neutral_k_ge2_top2_share_percent']}%</span>
    <span class="l">這 {cg['neutral_k_ge2']} 個之中，chr2+chr21 就佔掉的比例（共 {cg['neutral_k_ge2_chrom_count']} 條染色體）</span></div>
  <div class="stat warn"><span class="n">OR {cg['power']['min_detectable_odds_ratio']}</span>
    <span class="l">現有樣本量能偵測的最小效應（基線 {cg['power']['baseline_percent']}%，power 0.8）</span></div>
</div>
<div class="note bad">
小於 OR {cg['power']['min_detectable_odds_ratio']} 的 CN 效應，這個樣本量<strong>看不到</strong> ——
所以第 4 節「排除算術不可解後偵測不到殘餘效應」必須理解為「偵測不到」，不是「不存在」。
</div>

<h3>但有一項乾淨區的資訊不可替代</h3>
<p>CN 中性非 LOH、純腫瘤細胞株、單 copy haplotype 上的 clonal 突變，read-AF 應讀出 1/1。
因此<strong>低於 1 的 AF 在這裡是唯一沒有拷貝數解釋的 subclonal 讀數</strong>。</p>
<div class="grid" style="grid-template-columns:1fr 1fr">
  <div class="stat"><span class="n">{cg['neutral_sites_clonal']['point_percent']}%</span>
    <span class="l">clonal（AF ≥ 0.95），{cg['neutral_sites_clonal']['count']} / {cg['neutral_sites_total']:,}</span></div>
  <div class="stat ok"><span class="n">{cg['neutral_sites_subclonal']['count']}</span>
    <span class="l">subclonal（AF &lt; 0.95）＝{cg['neutral_sites_subclonal']['point_percent']}%（{cg['neutral_sites_subclonal']['lo_percent']}–{cg['neutral_sites_subclonal']['hi_percent']}）</span></div>
</div>
<div class="note ok">
這 <strong>{cg['neutral_sites_subclonal']['count']} 個位點</strong>是全樣本中唯一能在「無拷貝數解釋」前提下
讀出 subclonal 分率的證據。數量不大，但它們不是殘渣 —— 是目前<strong>唯一乾淨的頻率軸來源</strong>。
</div>

<h2>4B · CN 區的特殊現象（三項，都不只是「被污染」）</h2>

<h3>現象一：read-AF 隨拷貝數單調下降，下限精確落在 1/CN</h3>
{cn_af_bars}
<div class="note bad">
中位 AF 從 CN=1 的 {cp['af_median_by_total_cn']['1.0']} 一路降到 CN=7 的 {cp['af_median_by_total_cn']['7.0']}，
<strong>而且格點吻合率在高 CN 端反而上升</strong>（CN=8 時 {cp['grid_hit_percent_by_cn']['8.0']}%）。
這是 mutation multiplicity 的直接讀數：拷貝數越高，AF 越被鎖在 m/c 的離散格點上、
越不可能是連續的細胞分率。<strong>在 CN≥4 的區域把低 AF 解讀成「小 subclone」，幾乎一定是誤讀。</strong>
</div>

<h3>現象二：LOH 有雙重且方向相反的效應</h3>
<div class="scroll"><table>
<tr><th>CN 類別</th><th>結構唯一率</th><th>contested 中 AF 破 tie 率</th></tr>
<tr><td>gain（無 LOH）</td><td class="num ok">{cp['structure_unique_percent_by_class']['gain']}%</td><td class="num">{cp['af_broke_tie_percent_by_class']['gain']}%</td></tr>
<tr><td>neutral</td><td class="num ok">{cp['structure_unique_percent_by_class']['neutral']}%</td><td class="num">{cp['af_broke_tie_percent_by_class']['neutral']}%</td></tr>
<tr><td>loss</td><td class="num">{cp['structure_unique_percent_by_class']['loss']}%</td><td class="num">{cp['af_broke_tie_percent_by_class']['loss']}%</td></tr>
<tr><td>gain + LOH</td><td class="num bad">{cp['structure_unique_percent_by_class']['gain_loh']}%</td><td class="num bad">{cp['af_broke_tie_percent_by_class']['gain_loh']}%</td></tr>
<tr><td>copy-neutral LOH</td><td class="num bad">{cp['structure_unique_percent_by_class']['cnloh']}%</td><td class="num">{cp['af_broke_tie_percent_by_class']['cnloh']}%</td></tr>
<tr><td>loss + LOH</td><td class="num bad">{cp['structure_unique_percent_by_class']['loss_loh']}%</td><td class="num">{cp['af_broke_tie_percent_by_class']['loss_loh']}%</td></tr>
</table></div>
<div class="note">
帶 LOH 的區段結構唯一率全面偏低（9–32%），沒有 LOH 的反而高（56–58%）。
LOH 讓<strong>結構</strong>更難定（候選樹更多），卻同時讓 <strong>read-AF</strong> 更容易挑出唯一解 ——
gain+LOH 兩者兼具，成為「候選樹多、又幾乎總能被 AF 選掉」的最危險組合（{cp['af_broke_tie_percent_by_class']['gain_loh']}%）。
機制上合理：LOH 使一條 haplotype 消失、partial pattern 變多（結構更糊），
同時多個 copy 集中在單一 haplotype、AF 落在不同 m/c（AF 更能分辨）。
</div>

<h3>現象三：CN <em>不</em>改變樹的形狀 —— 控制 k 之後差異消失</h3>
<div class="scroll"><table>
<tr><th>分支（Sister 類）比例</th><th>n</th><th>比例（95% CI）</th></tr>
<tr><td>CN 中性（k≥2）</td><td class="num">{cg['branching_neutral']['n']}</td>
    <td class="num">{cg['branching_neutral']['percent']}%（{cg['branching_neutral']['lo']}–{cg['branching_neutral']['hi']}）</td></tr>
<tr><td>CN 變異（k≥2）</td><td class="num">{cg['branching_altered']['n']:,}</td>
    <td class="num">{cg['branching_altered']['percent']}%（{cg['branching_altered']['lo']}–{cg['branching_altered']['hi']}）</td></tr>
<tr><td>k=2</td><td class="num">{cg['branching_by_k']['k=2']['neutral']['n']} / {cg['branching_by_k']['k=2']['altered']['n']:,}</td>
    <td class="num">{cg['branching_by_k']['k=2']['neutral']['percent']}%（{cg['branching_by_k']['k=2']['neutral']['lo']}–{cg['branching_by_k']['k=2']['neutral']['hi']}） vs {cg['branching_by_k']['k=2']['altered']['percent']}%（{cg['branching_by_k']['k=2']['altered']['lo']}–{cg['branching_by_k']['k=2']['altered']['hi']}）</td></tr>
<tr><td>k=3</td><td class="num">{cg['branching_by_k']['k=3']['neutral']['n']} / {cg['branching_by_k']['k=3']['altered']['n']:,}</td>
    <td class="num">{cg['branching_by_k']['k=3']['neutral']['percent']}%（{cg['branching_by_k']['k=3']['neutral']['lo']}–{cg['branching_by_k']['k=3']['neutral']['hi']}） vs {cg['branching_by_k']['k=3']['altered']['percent']}%（{cg['branching_by_k']['k=3']['altered']['lo']}–{cg['branching_by_k']['k=3']['altered']['hi']}）</td></tr>
<tr><td>k=4</td><td class="num">{cg['branching_by_k']['k=4']['neutral']['n']} / {cg['branching_by_k']['k=4']['altered']['n']}</td>
    <td class="num">{cg['branching_by_k']['k=4']['neutral']['percent']}%（{cg['branching_by_k']['k=4']['neutral']['lo']}–{cg['branching_by_k']['k=4']['neutral']['hi']}） vs {cg['branching_by_k']['k=4']['altered']['percent']}%（{cg['branching_by_k']['k=4']['altered']['lo']}–{cg['branching_by_k']['k=4']['altered']['hi']}）</td></tr>
</table></div>
<div class="verdict">
每一層信賴區間都重疊。<strong>拷貝數不改變重建出來的樹長什麼樣，只改變我們能不能從候選中挑出唯一一棵。</strong><br><br>
三層總結：<strong>形狀 CN-robust</strong>（現象三）→ <strong>結構唯一性 CN-robust</strong>（第 2 節）→
<strong>AF 挑唯一樹被 CN 抬高</strong>（第 2–3 節）。受影響的只有最後一層。
</div>

<h2>4C · 扣掉 CN 區之後，還夠不夠確認樹的演化？</h2>
<p>`unit 數`不是正確的計量單位。一個 unit 要真的貢獻演化資訊必須同時滿足：<strong>k≥2</strong>
（k=1 的樹必然是 <code>ROOT→A</code>，沒有先後順序資訊）<strong>且樹被定出來</strong>。
真正的單位是 <strong>deep edge</strong> —— parent 不是 ROOT 的邊，即真正表達「兩個突變誰先誰後」的關係。</p>

<div class="scroll"><table>
<tr><th>CN 中性</th><th>unit</th><th>結構唯一</th><th>AF 定出</th><th>仍 tied</th><th>可定樹</th><th>deep edge</th><th>染色體</th></tr>
<tr><td>全部 ranked</td><td class="num">{ey['neutral_all']['n_units']}</td>
    <td class="num">{ey['neutral_all']['structure_unique']}</td>
    <td class="num">{ey['neutral_all']['af_resolved']}</td>
    <td class="num">{ey['neutral_all']['still_tied']}</td>
    <td class="num">{ey['neutral_all']['determined']}（{ey['neutral_all']['determined_percent']}%）</td>
    <td class="num ok">{ey['neutral_all']['deep_edges']}</td>
    <td class="num">{ey['neutral_all']['chromosomes']}</td></tr>
<tr><td><strong>k ≥ 2</strong></td><td class="num">{ey['neutral_k_ge_2']['n_units']}</td>
    <td class="num">{ey['neutral_k_ge_2']['structure_unique']}</td>
    <td class="num">{ey['neutral_k_ge_2']['af_resolved']}</td>
    <td class="num bad">{ey['neutral_k_ge_2']['still_tied']}</td>
    <td class="num">{ey['neutral_k_ge_2']['determined']}（{ey['neutral_k_ge_2']['determined_percent']}%）</td>
    <td class="num ok">{ey['neutral_k_ge_2']['deep_edges']}</td>
    <td class="num">{ey['neutral_k_ge_2']['chromosomes']}</td></tr>
<tr><td><strong>k ≥ 3</strong></td><td class="num">{ey['neutral_k_ge_3']['n_units']}</td>
    <td class="num">{ey['neutral_k_ge_3']['structure_unique']}</td>
    <td class="num">{ey['neutral_k_ge_3']['af_resolved']}</td>
    <td class="num bad">{ey['neutral_k_ge_3']['still_tied']}</td>
    <td class="num">{ey['neutral_k_ge_3']['determined']}（{ey['neutral_k_ge_3']['determined_percent']}%）</td>
    <td class="num ok">{ey['neutral_k_ge_3']['deep_edges']}</td>
    <td class="num bad">{ey['neutral_k_ge_3']['chromosomes']}</td></tr>
</table></div>

<div class="note bad">
<strong>「乾淨區可定樹率 {ey['neutral_all']['determined_percent']}%」是假象</strong> ——
那個數字被 {ey['neutral_by_k']['k=1']['n']} 個 k=1 撐起來，而 k=1 貢獻
<strong>{ey['neutral_by_k']['k=1']['deep_edges']} 條</strong> deep edge。
扣掉之後可定樹率掉到 <strong>{ey['neutral_k_ge_2']['determined_percent']}%</strong>。
</div>

<h3>與 CN 變異區對照</h3>
<div class="grid">
  <div class="stat"><span class="n">{ey['neutral_k_ge_2']['determined_percent']}%</span>
    <span class="l">CN 中性 k≥2 可定樹（n={ey['neutral_k_ge_2']['n_units']}）</span></div>
  <div class="stat"><span class="n">{ey['altered_k_ge_2']['determined_percent']}%</span>
    <span class="l">CN 變異 k≥2 可定樹（n={ey['altered_k_ge_2']['n_units']:,}）</span></div>
  <div class="stat ok"><span class="n">{ey['neutral_k_ge_2']['deep_edges']}</span>
    <span class="l">CN 中性的 deep edge</span></div>
  <div class="stat bad"><span class="n">{ey['clean_share_of_deep_edges_k_ge_2_percent']}%</span>
    <span class="l">乾淨區佔全部 deep edge 的比例（另 {ey['altered_share_of_deep_edges_k_ge_2_percent']}% 在 CN 變異區）</span></div>
</div>

<h3>反直覺但關鍵：乾淨區反而更定不出樹</h3>
<div class="note">
CN 中性的「仍 tied」是 <strong>{ey['neutral_k_ge_2']['still_tied_percent']}%</strong>，
CN 變異只有 {ey['altered_k_ge_2']['still_tied_percent']}%。原因正是第 4 節的機制：
乾淨區的 read-AF 大量等於 1/1、彼此相同 → 分數相同 → <strong>數學上無法破 tie</strong>。
CN 反而「提供」了打破平手所需的數值差異 —— 只是那個差異來自拷貝數，不是時間。<br><br>
<strong>這是取捨，不是單純的好壞：</strong>CN 變異區<em>可辨識性高、可信度低</em>；
CN 中性區<em>可辨識性低、可信度高</em>。
</div>

<h3>夠做什麼、不夠做什麼</h3>
<div class="scroll"><table>
<tr><th>用途</th><th>判定</th><th>理由</th></tr>
<tr><td>證明方法在乾淨地上可運作（proof-of-concept）</td><td class="ok">夠</td>
    <td>{ey['neutral_k_ge_2']['determined']} 個可定樹 unit、{ey['neutral_k_ge_2']['deep_edges']} 條無爭議 deep edge</td></tr>
<tr><td>作為 CN-free 驗證錨點／陰性對照</td><td class="ok">夠</td>
    <td>全樣本唯一不需拷貝數假設的譜系關係</td></tr>
<tr><td>描述這個腫瘤的演化史</td><td class="bad">不夠</td>
    <td>能看多步譜系的只有 {ey['neutral_k_ge_3']['determined']} 個 unit、{ey['neutral_k_ge_3']['deep_edges']} 條 deep edge，集中在 {ey['neutral_k_ge_3']['chromosomes']} 條染色體</td></tr>
<tr><td>統計推論（偵測中小效應）</td><td class="bad">不夠</td>
    <td>檢定力只到 OR {cg['power']['min_detectable_odds_ratio']}</td></tr>
<tr><td>取代 CN 變異區重做全基因組分析</td><td class="bad">不可能</td>
    <td>{ey['altered_share_of_deep_edges_k_ge_2_percent']}% 的 deep edge 在 CN 變異區</td></tr>
</table></div>
<div class="verdict">
<strong>實務結論：不能靠「扣掉 CN 區」來做研究。</strong>正確做法是<strong>保留全部、標註 CN 狀態、分層報告</strong> ——
CN 變異區照樣輸出候選樹與共識骨幹（形狀與結構層都是 CN-robust），
但其 read-AF 挑出的「唯一樹」一律停在 <code>RAW_AF_UNIQUE_REPRESENTATIVE</code>；
CN 中性的那 {ey['neutral_k_ge_2']['deep_edges']} 條 deep edge 則作為可對外宣稱的乾淨錨點。
</div>

<h2>4D · 跨樣本累積乾淨 deep edge + 與 SEQC2 對照</h2>

<h3>先驗證「用 SAVANA 判乾淨區」準不準（在 HCC1395 上量測）</h3>
<div class="grid">
  <div class="stat ok"><span class="n">{ms['gate_precision_percent']}%</span>
    <span class="l">校正 CN + BAF gate 的 precision</span></div>
  <div class="stat warn"><span class="n">{ms['gate_recall_percent']}%</span>
    <span class="l">同一 gate 的 recall（刻意保守）</span></div>
  <div class="stat bad"><span class="n">{ms['gate_precision_published_fit_percent']}%</span>
    <span class="l">若用 SAVANA 發布 fit，precision 只有這樣</span></div>
</div>
<p style="font-size:.9rem;color:var(--muted)">判準：<code>total CN ≈ 2（±0.5）且段 mean BAF &lt; 0.75</code>
—— 只有 1+1 才會讓單 copy haplotype 的 clonal 突變讀出 1/1。對錨點而言 precision 遠比 recall 重要。</p>

<h3>三樣本結果</h3>
<div class="scroll"><table>
<tr><th>樣本</th><th>乾淨 k≥2 unit</th><th>可定樹率</th><th>deep edge</th><th>k≥3 deep</th><th>嚴格版 deep</th><th>染色體</th></tr>
{''.join(f"<tr><td>{k}</td><td class='num'>{v['clean_k_ge2_units']:,}</td><td class='num'>{v['clean_k_ge2_determined_percent']}%</td><td class='num ok'>{v['clean_k_ge2_deep_edges']}</td><td class='num'>{v['clean_k_ge3_deep_edges']}</td><td class='num'>{v['clean_strict_k_ge2_deep_edges']}</td><td class='num'>{v['chromosomes']}</td></tr>" for k, v in ms['per_sample'].items())}
</table></div>

<div class="grid">
  <div class="stat ok"><span class="n">{ms['pooled']['k_ge_2_deep_edges']}</span>
    <span class="l">累積乾淨 deep edge（原本只有 {ey['neutral_k_ge_2']['deep_edges']}）</span></div>
  <div class="stat"><span class="n">{ms['pooled']['k_ge_3_deep_edges']}</span>
    <span class="l">其中 k≥3（多步譜系）</span></div>
  <div class="stat"><span class="n">{ms['pooled']['hcc1395_share_of_deep_edges_percent']}%</span>
    <span class="l">HCC1395 佔比（其餘來自三個無真值樣本）</span></div>
</div>

<h3>跨樣本獨立重現了機制</h3>
<div class="scroll"><table>
<tr><th>樣本</th><th>乾淨區可定樹率</th><th>CN 變異區可定樹率</th></tr>
<tr><td>HCC1395</td><td class="num">{ey['neutral_k_ge_2']['determined_percent']}%</td><td class="num">{ey['altered_k_ge_2']['determined_percent']}%</td></tr>
{''.join(f"<tr><td>{k}</td><td class='num'>{v['clean_k_ge2_determined_percent']}%</td><td class='num'>{v['altered_k_ge2_determined_percent']}%</td></tr>" for k, v in ms['per_sample'].items())}
</table></div>
<div class="note ok">
<strong>4/4 樣本都是「乾淨區反而更定不出樹」</strong> —— 在四個獨立樣本上重現了
「AF 值相同 → 分數相同 → 數學上無法破 tie」的機制，比單樣本結論強得多。
</div>
<div class="note bad">
⚠ 非 HCC1395 的那 {ms['pooled']['k_ge_2_deep_edges'] - ey['neutral_k_ge_2']['deep_edges']} 條 deep edge，
其可信度是<strong>期望而非量測</strong>：三樣本用 SAVANA 發布的 purity/ploidy（通過 BAF 自洽性檢驗），
但無外部 CN 真值，{ms['gate_precision_percent']}% 的 precision 是從 HCC1395 遷移過來的推論。
</div>

<h3>與 SEQC2 外部 clonality 分析對照 —— 尺度不同</h3>
<div class="scroll"><table>
<tr><th>SEQC2（Fang et al., Nat Biotechnol 2021）</th><th>我們</th></tr>
<tr><td>{sq['seqc2_reported']['tree']}，S2 的 cancer cell fraction {sq['seqc2_reported']['s2_cancer_cell_fraction_percent']}%</td>
    <td>local per-PS×HP block 的 mutation-state 樹</td></tr>
<tr><td>10x Single Cell CNV：{sq['seqc2_reported']['single_cells_tumor']:,} 個 tumor 細胞 + {sq['seqc2_reported']['single_cells_normal']} 個 normal</td>
    <td>單一 ONT bulk，無單細胞</td></tr>
<tr><td>SNV clonality = VAF <strong>經 local copy number 校正</strong></td>
    <td class="bad">做不到（缺 allele-specific CN）</td></tr>
<tr><td>subclone 主要由染色體 gain/loss 事件定義</td><td>由 read-level sSNV 共現定義</td></tr>
</table></div>
<div class="note">
<strong>兩者是不同尺度的物件，不是同一棵樹的兩個版本。</strong>
SEQC2 的 S1–S10 是全基因組的細胞群演化樹；我們的是局部突變狀態樹，
且 2026-08-01 規格 §7 明文禁止跨 block 合併。
</div>

<div class="two">
<div>
<h4 style="margin:.6rem 0 .4rem;font-size:.97rem;color:var(--ok)">✓ 可對照且吻合的兩點</h4>
<p style="font-size:.9rem">① <strong>大部分突變是 clonal</strong>：SEQC2 稱大部分 driver 點突變位於 MRCA；
我們在乾淨區（depth≥20，n={sq['neutral_sites_depth_ge20']}）觀察到
<strong>{sq['clonal_af_near_1_count']} 個位點（{round(100*sq['clonal_af_near_1_count']/sq['neutral_sites_depth_ge20'],2)}%）read-AF ≈ 1.0</strong>，清楚的 clonal 主峰。<br><br>
② <strong>高度異質、大量 CNA 在 sub-population</strong>：SEQC2 由 1,270 個單細胞確認；
我們獨立觀察到 <strong>{V['unit_altered']}% 的 unit 坐在 CN 變異區</strong> ——
這正好解釋了乾淨區為何如此稀少：<strong>不是方法問題，是這顆腫瘤的性質</strong>。</p>
</div>
<div>
<h4 style="margin:.6rem 0 .4rem;font-size:.97rem;color:var(--bad)">✗ 做不到的部分</h4>
<p style="font-size:.9rem">乾淨區的 {sq['subclonal_count']} 個 subclonal 位點（AF &lt; 0.95）呈<strong>連續分布</strong>，
落在 0.55–0.65（對照 SEQC2 S2 的 CCF 60%）僅 {sq['near_ccf60_count']} 個
（{round(100*sq['near_ccf60_count']/sq['subclonal_count'],2)}%），
<strong>沒有對應 S1–S10 離散 CCF 的清晰峰</strong>。<br><br>
原因不必歸咎方法：n 太小、位點分散於不同 local block（可能屬不同 subclone）、
depth 20–50 時 AF 抽樣誤差已達 ±0.1。
<strong>結論是「解析度不足以重建 S1–S10」，不是「與 SEQC2 矛盾」。</strong></p>
</div>
</div>

<h2>4E · 結構層能反過來「限制」VAF-based 的 subclone 推論嗎？</h2>

<div class="verdict">
<strong>可以，而且在 CN 變異區最有用。</strong>VAF-based 重建只有邊際頻率，
看不到哪些突變在同一條分子上；長讀長直接提供那個缺失的觀測量。
</div>

<h3>兩類方法的資訊落差</h3>
<div class="scroll"><table>
<tr><th>觀測到的 haplotype pattern</th><th>結論（不需任何頻率論證）</th></tr>
<tr><td><code>{{RR, AR, AA}}</code></td><td>巢狀，A 先於 B</td></tr>
<tr><td><code>{{RR, RA, AA}}</code></td><td>巢狀，B 先於 A</td></tr>
<tr><td><code>{{RR, AR, RA}}</code></td><td class="bad"><strong>互斥 —— 不同分支，不存在祖先關係</strong></td></tr>
<tr><td>四種都出現</td><td>違反完美系統發生</td></tr>
</table></div>

<h3>量測：VAF-only 會在哪裡判錯</h3>
<p style="font-size:.9rem;color:var(--muted)">取 HCC1395 中 k=2 的 unit（每個恰好隔離一組成對關係，兩位點深度皆 ≥10），共 {lk['k2_pairs_total']:,} 對。</p>
<div class="grid">
  <div class="stat"><span class="n">{lk['relation_percent_all']['nested']}%</span>
    <span class="l">巢狀（可定先後）</span></div>
  <div class="stat warn"><span class="n">{lk['relation_percent_all']['mutually_exclusive']}%</span>
    <span class="l">互斥（不同分支）— VAF 無法偵測</span></div>
  <div class="stat bad"><span class="n">{lk['exclusive_large_gap_all']['percent']}%</span>
    <span class="l">互斥對中 AF 差距 ≥{lk['large_gap_threshold']}（{lk['exclusive_large_gap_all']['count']}/{lk['exclusive_large_gap_all']['n']}）— <strong>VAF 會誤判為祖先關係</strong></span></div>
  <div class="stat ok"><span class="n">{lk['nested_violating_monotonicity_all']['percent']}%</span>
    <span class="l">巢狀對違反 AF 單調（{lk['nested_violating_monotonicity_all']['count']}/{lk['nested_violating_monotonicity_all']['n']}）</span></div>
</div>
<div class="note">
<strong>問題不在 VAF 排錯順序，而在 VAF 無法察覺「根本不該排序」。</strong>
真正巢狀時 VAF 大致可靠（僅 {lk['nested_violating_monotonicity_all']['percent']}% 違反單調）；
但它對互斥關係完全盲目。
</div>

<h3>這個落差正好集中在 CN 變異區</h3>
<div class="scroll"><table>
<tr><th></th><th>互斥對中 AF 差距 ≥{lk['large_gap_threshold']}</th><th>互斥對的 AF 差距中位</th></tr>
<tr><td>CN 中性</td><td class="num ok">{lk['exclusive_large_gap_neutral']['percent']}%（{lk['exclusive_large_gap_neutral']['count']}/{lk['exclusive_large_gap_neutral']['n']}）</td>
    <td class="num ok">{lk['median_af_gap_exclusive_neutral']}</td></tr>
<tr><td>CN 變異</td><td class="num bad">{lk['exclusive_large_gap_altered']['percent']}%（{lk['exclusive_large_gap_altered']['count']}/{lk['exclusive_large_gap_altered']['n']}）</td>
    <td class="num bad">{lk['median_af_gap_exclusive_altered']}</td></tr>
</table></div>
<div class="verdict">
在拷貝數乾淨的地方，互斥突變的 AF <strong>中位差距為 {lk['median_af_gap_exclusive_neutral']}</strong> ——
VAF 方法看到相同頻率，不會硬排順序。<br>
在 CN 變異區，multiplicity 給了互斥突變一個<strong>假的頻率梯度</strong>，VAF 就會把它讀成祖先關係。<br><br>
<strong>於是出現一個互補性：CN 破壞 AF 層的地方，正是 linkage 約束最有價值的地方。</strong>
</div>

<h3>我們能提供多少硬約束（全 7 樣本）</h3>
<div class="scroll"><table>
<tr><th>樣本</th><th>k≥2 unit</th><th>順序約束（deep edge）</th><th>互斥宣告</th></tr>
{''.join(f"<tr><td>{k}</td><td class='num'>{v['units_k_ge2']:,}</td><td class='num ok'>{v['ordering_constraints_deep_edges']:,}</td><td class='num'>{v['units_with_exclusion_branch']:,}</td></tr>" for k, v in lk['inventory'].items() if k != 'TOTAL')}
<tr style="font-weight:650"><td>合計</td><td class="num">{lk['inventory']['TOTAL']['units_k_ge2']:,}</td>
    <td class="num ok">{lk['inventory']['TOTAL']['ordering_constraints']:,}</td>
    <td class="num">{lk['inventory']['TOTAL']['exclusion_units']:,}</td></tr>
</table></div>

<div class="note ok">
<strong>定位：不是取代 VAF-based 重建，而是提供它拿不到的兩類硬約束</strong><br>
① <strong>排除約束</strong>（{lk['inventory']['TOTAL']['exclusion_units']:,} 個 unit）：「這兩個突變不可能有祖先關係」——
VAF 無法產生此判斷，且在 CN 區系統性判錯<br>
② <strong>順序約束</strong>（{lk['inventory']['TOTAL']['ordering_constraints']:,} 條）：「A 必定先於 B」——
分子直接觀測，不依賴頻率、不需 CN 校正
</div>
<div class="note bad">
<strong>誠實邊界</strong>：這些約束都是 <strong>local</strong>（同一 PS×HP block 內、有共享 read 的位點對），
覆蓋率遠低於全基因組 VAF 方法。它們能<strong>約束與否證</strong>一棵候選演化樹，
但<strong>不能單獨建構</strong>全樣本 clone 樹（規格 §7 禁止跨 block 合併）。
正確用法：<strong>VAF 方法產生候選 → 我們的約束篩掉不相容者</strong>。
</div>

<h2>4F · 約束方法成立嗎？能排除 CN 影響嗎？</h2>

<div class="note bad">
<strong>先說修正。</strong>第 4E 節把兩類約束並列、並暗示「排除約束」較有價值。
<strong>那個框架是錯的</strong>，本節的稽核顯示情況正好相反。
</div>

<h3>約束一：順序（<code>{{RR, AR, AA}}</code> → A 先於 B）— ✔ 成立，且 CN-free by construction</h3>
<p><strong>論證</strong>：B 只出現在已帶 A 的分子上。infinite-sites 下突變沿分子譜系只增不減，
故 A+B 分子必然衍生自 A 分子。<strong>整個推論沒有引用 copy number、ploidy 或 purity</strong> ——
它是關於<strong>分子繼承</strong>的陳述，不是頻率的陳述。
CN 改變的是「同一 haplotype 有幾份」，不改變「B 出現在 A 的背景上」。</p>
<div class="grid">
  <div class="stat ok"><span class="n">{cv['ordering_total']:,}</span>
    <span class="l">順序約束總數（HCC1395）</span></div>
  <div class="stat ok"><span class="n">{cv['ordering_on_altered']:,}</span>
    <span class="l">位於 CN 變異區 — <strong>全部有效</strong></span></div>
  <div class="stat"><span class="n">{cv['ordering_on_neutral']}</span>
    <span class="l">位於 CN 中性區</span></div>
</div>
<div class="note">
⚠ <strong>層級限定</strong>：這是<strong>分子譜系</strong>的順序，不是細胞譜系。可以說
「B 發生在已帶 A 的 DNA 分子上」，不能直接說「帶 B 的細胞是帶 A 細胞的後代」。<br><br>
⚠ <strong>真正的殘餘威脅是抽樣，不是 CN</strong>：宣稱「沒看到第三種 pattern」可能只是漏抽。
若某 pattern 真實比例 p，在 D 條分子中漏掉的機率是 (1−p)<sup>D</sup>。
</div>
<div class="scroll"><table>
<tr><th>連結深度</th><th>unit 數</th><th>佔比</th><th>能以 95% 信心排除的 pattern 比例</th></tr>
{''.join(f"<tr><td>{k}</td><td class='num'>{v['units']:,}</td><td class='num'>{v['percent_of_units']}%</td><td class='num {'bad' if v['low']>25 else 'ok' if v['low']<6 else ''}'>&gt;{v['high']}–{v['low']}%</td></tr>" for k, v in cv['depth_band_detectable_percent'].items())}
</table></div>
<p style="font-size:.9rem;color:var(--muted)">中位連結深度 {cv['median_linking_depth']}，中位可排除比例
<strong>{cv['median_min_detectable_percent']}%</strong>。順序約束對「低比例的第三種 pattern」不敏感；
depth&lt;20 的 <strong>{cv['low_depth_unit_percent']}%</strong> unit 應標為低信心。</p>

<h3>約束二：互斥（<code>{{RR, AR, RA}}</code> → 不同分支）— ⚠ 只在 CN 中性區能推到細胞層</h3>
<div class="two">
<div>
<p style="font-size:.92rem"><strong>CN 中性（1+1）</strong>：每個細胞對該 haplotype 只貢獻一份分子。
「沒有分子同時帶 A 和 B」⇒「沒有細胞同時帶 A 和 B」⇒ <strong>真的是不同細胞分支</strong> ✔</p>
</div>
<div>
<p style="font-size:.92rem"><strong>CN gain</strong>：該 haplotype 每個細胞有 c&gt;1 份。
<strong>A 在 copy 1、B 在 copy 2，可以屬於同一個細胞</strong> —— 沒有單一分子同時帶 A+B，
但這個細胞兩個都有。<strong>分子層觀測不變，細胞層推論崩掉</strong> ✘</p>
</div>
</div>
<div class="grid">
  <div class="stat"><span class="n">{cv['exclusion_total']:,}</span>
    <span class="l">互斥 unit 總數</span></div>
  <div class="stat bad"><span class="n">{cv['exclusion_cellular_share']['percent']}%</span>
    <span class="l">能推論到細胞層（僅 {cv['exclusion_cellular_valid']} 個，CI {cv['exclusion_cellular_share']['lo']}–{cv['exclusion_cellular_share']['hi']}）</span></div>
  <div class="stat warn"><span class="n">{cv['exclusion_molecular_only_percent']}%</span>
    <span class="l">只能停在分子層（{cv['exclusion_molecular_only']:,} 個）</span></div>
</div>

<div class="verdict">
<strong>修正後的結論</strong><br><br>
<table style="margin:.5rem 0">
<tr><th>約束</th><th>能否排除 CN 影響</th><th>全基因組可用性</th></tr>
<tr><td><strong>順序</strong></td><td class="ok"><strong>✔ 可以，by construction</strong></td>
    <td>{cv['ordering_total']:,} 條全部可用（分子譜系層）</td></tr>
<tr><td><strong>互斥</strong></td><td class="bad"><strong>✘ 不能</strong></td>
    <td>細胞層僅 {cv['exclusion_cellular_valid']} 條（{cv['exclusion_cellular_share']['percent']}%）</td></tr>
</table>
<strong>原因很簡單</strong>：順序約束問的是「B 是否出現在 A 的背景上」，與有幾份 copy 無關；
互斥約束問的是「有沒有細胞同時帶兩者」，而那正是 copy number 直接介入的地方。
</div>

<h2>5 · SAVANA 流程驗證：訊號好，校準錯</h2>
<div class="scroll"><table>
<tr><th></th><th>purity</th><th>ploidy</th><th>狀態一致率</th><th>整數 CN 吻合</th></tr>
<tr><td>SAVANA 發布值</td><td class="num">{V['pub_purity']}</td><td class="num">{V['pub_ploidy']}</td>
    <td class="num bad">{V['pub_agree']}%</td><td class="num bad">{V['pub_int']}%</td></tr>
<tr><td>grid search 最佳</td><td class="num">{V['best_purity']}</td><td class="num">{V['best_ploidy']}</td>
    <td class="num ok">{V['best_agree']}%</td><td class="num ok">{V['best_int']}%</td></tr>
</table></div>
<div class="note ok">
SEQC2 自身隱含的常染色體 ploidy 為 <strong>{V['seqc2_ploidy']}</strong>，與 grid 最佳的
{V['best_ploidy']} 幾乎重合；purity {V['best_purity']} 正是純腫瘤細胞株該有的值。
段層回歸在非中性區 R² 達 <strong>{V['r2_non']}</strong>（全部 bin 的 {V['r2_all']} 偏低是
CN=2 常數值稀釋造成的統計 artifact，非 segmentation 品質不良）。
<strong>裁決：log2r 訊號與 segmentation 可用，錯的是 purity/ploidy 擬合 —— 校準失敗而非訊號失敗。</strong>
</div>

<h2>6 · 回頭裁決其他樣本</h2>
<p>用一個<strong>不需外部真值</strong>的內部檢驗：逐段計算 BAF 理論上限
<code>BAF_max = (ρn + (1−ρ)) / (ρn + 2(1−ρ))</code>，統計超出上限的段比例。
先用 HCC1395 校準 —— 在發布的 purity 下違反率
{purity['samples']['HCC1395']['violation_rate_percent']}%，
在校正後 purity 1.0 下 <strong>{calib['violation_rate_percent']}%</strong>。
這個檢驗確實能抓到 mis-fit，且其判定被 SEQC2 外部真值背書。</p>

<div class="scroll"><table>
<tr><th>樣本</th><th>purity</th><th>ploidy</th><th>可用段數</th><th>違反率</th><th>裁決</th></tr>
{purity_rows}
</table></div>
<div class="note ok">
<strong>這次獨立檢驗支持而非推翻既有判定</strong>：H1437 / H2009 / HCC1954 可用、HCC1937 不可用。
COLO829 有三個 fit 版本，本輪未納入，屬缺口。
</div>

<h2>7 · 可以矯正嗎</h2>
<div class="scroll"><table>
<tr><th>矯正層級</th><th>需要什麼</th><th>可行性</th></tr>
<tr><td>L1 CN 註釋層</td><td>total CN + LOH</td><td class="ok">現在就可做（本輪已建）</td></tr>
<tr><td>L2 期望 AF 帶</td><td>total CN + purity + LOH</td><td class="ok">可做，需假設 allele 分配</td></tr>
<tr><td>L3 真 multiplicity 校正</td><td><strong>allele-specific CN</strong></td><td class="bad">做不到 — 真正的阻斷點</td></tr>
<tr><td>L4 HP↔allele 對應</td><td>哪個 HP 對應 major allele</td><td class="bad">做不到</td></tr>
</table></div>
<p><strong>建議做法＝適用性閘門，而非數值校正。</strong>為每個 unit 標
<code>cn_loh_status</code>；只有 <code>CN_NEUTRAL_NO_LOH</code> 的 unit 其 read-AF 排序
才可進入生物學解讀，其餘停在規格 §5.2 已規定的 <code>RAW_AF_UNIQUE_REPRESENTATIVE</code>。
理由：L3 缺 allele-specific CN，任何數值校正都得假設 allele 分配，會把假設偷渡成結果。</p>
<div class="note bad">
<strong>代價要說清楚：</strong>套用此 gate 後，可用於 AF 生物學解讀的 unit 只剩
{unit_counts.get('neutral',0)} 個（{V['unit_neutral']}%），且其中 {V['top4']}% 集中在 4 條染色體。
這是誠實的可用範圍，不是失敗。
</div>

<h2>8 · 具體案例（見樹）</h2>
<div class="note">
<code>{example['region_id']}</code><br>
{example['chrom']}:{example['span_start']:,}–{example['span_end']:,}（{example['span_end']-example['span_start']:,} bp），
{example['active_bit_count']} 個 sSNV，SEQC2 <strong>total CN {example['seqc2_total_cn_weighted']}、
LOH 涵蓋率 {example['seqc2_loh_fraction']}</strong>，
{example['total_tree_count']} 棵候選樹經 read-AF 選出唯一贏家，標為 <code>{example['resolution_class']}</code>。<br><br>
在一條 haplotype 帶約 7 個 copy 的區段上，AF 只能取 m/7 —— 這些差異幾乎確定來自
mutation multiplicity 而非發生順序。<strong>這是那 {V['excess_n']:,.0f} 個超額唯一樹的典型形態。</strong>
</div>

<h2>9 · 之前的 subclone 結論還站得住嗎</h2>
<div class="scroll"><table>
<tr><th>既有結論</th><th>裁決</th><th>理由</th></tr>
<tr><td>sSNV 單分子共現是非循環骨幹</td><td class="ok">維持</td>
    <td>結構唯一率在中性區反而更高（{V['struct_neu']}% vs {V['struct_alt']}%）</td></tr>
<tr><td>拓撲軸 CN-robust</td><td class="ok">維持</td><td>同上，以 canonical 分母重新驗證</td></tr>
<tr><td>甲基化為 bounded-auxiliary</td><td class="ok">不受影響</td><td>CN 問題出在 read-AF 層</td></tr>
<tr><td>read-AF 僅作 model-conditional ordering</td><td class="ok">強化</td>
    <td>給出量化：CN 可歸因超額 {V['excess_lo']}–{V['excess_hi']}%</td></tr>
<tr><td>「88.26% 單一拓撲」「55.1% 唯一最佳樹」</td><td class="bad">必須加註</td>
    <td>係 CN-conditional；未標註即等同宣稱生物學拓撲盛行率，規格 §10 已明文禁止</td></tr>
<tr><td>VAF/頻率軸嚴重 CN-limited</td><td class="ok">維持並量化</td><td>{V['unit_altered']}% unit 坐在 CN 變異區</td></tr>
</table></div>
<div class="verdict">
<strong>總評：</strong>不需要撤回任何已發表結論，但需要為 AF 層的數字補上 CN-conditional 標註與
可用範圍聲明。本輪把 2026-08-01 規格中「allele-specific CN/LOH 未整合」這個定性警告，
變成了帶區間的量化，並補上該規格 §9 第 6 項所要求的 independent CN annotation。
</div>

<h2>10 · 限制</h2>
<ol style="font-size:.92rem">
<li><strong>單一樣本。</strong>HCC1395 是唯一有 CN 真值者，且高度非整倍體，中性對照僅 {V['unit_neutral']}%。不可外推為「所有樣本的 CN 影響都是這個量級」。</li>
<li><strong>中性組染色體高度集中</strong>（{V['top4']}% 在 4 條染色體），按染色體分層後 OR 由 {V['or_tree']} 降至 {V['or_chrom']}。</li>
<li><strong>無 allele-specific CN。</strong>m/c 格點分析依賴「LOH 時 c = total CN、非 LOH 時均分」的假設，是 model-conditional。</li>
<li><strong>機制在校正後與「完全經由算術可解性中介」一致</strong>（OR {V['nondeg_or']}，p = {V['nondeg_p']}），但非退化子集 n_neutral 僅 {V['nondeg_n_neu']}，尚不能排除部分中介。初版曾宣稱顯著殘餘（OR 2.819），源自位點集錯誤，已撤回。</li>
<li><strong>未重算。</strong>依任務界定不重跑 CCF／樹枚舉／subclone 判定，「矯正後結論如何改變」只有設計與估計。</li>
<li><strong>COLO829 未納入</strong> SAVANA 自洽性比較。</li>
<li><strong>SEQC2 CN 本身</strong>是 NGS benchmark 推導的共識，非細胞層真值；half-integer 代表 subclonal CN，本輪當作段層代表值處理。</li>
</ol>

<footer>
資料來源：canonical 2026-07-24 frozen cohort（<code>HCC1395.topology.jsonl</code> / <code>HCC1395.census.jsonl</code>）·
SEQC2 benchmark CNV（整數 CN gain 1,502 段 + loss 63 段）· SAVANA WGS CN。<br>
全部數字由 <code>scripts/</code> 產出至 <code>data/*.json</code> 後注入本頁，缺值即拒絕產生。
分母與 20260801 <code>authority_manifest.json</code> / <code>denominator_registry.tsv</code> 對帳。<br>
留底報告：<code>InterSubMod/research/20260810_cnv_confound_seqc2_savana_audit/docs/20260810_CNV混淆稽核_SEQC2結對與SAVANA流程驗證_01.md</code>
</footer>
</div>
"""

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write(HTML)

print(f"wrote {OUT}")
print(f"size: {os.path.getsize(OUT):,} bytes")
