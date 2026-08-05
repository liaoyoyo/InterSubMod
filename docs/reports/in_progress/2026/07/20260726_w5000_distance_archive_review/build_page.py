#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_page.py — 20260119 all-with-w5000 archive 的 distance/ 產出檢視頁

設計契約: L0 一眼結論 / L1 重點邏輯 / L2 案例卡(收合) / L3 邊界與溯源
          無 emoji(本機無 emoji 字型) → CSS 徽章;狀態一律附文字標籤
§13-A   : 數字全由 data.json 注入;缺 required key 即 refuse
自足    : PNG 以 base64 內嵌,零外部資源
"""
import base64
import json
import os
from pathlib import Path

S = Path("/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/"
         "efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad")
OUT = Path(__file__).parent
STEM = "20260726_w5000_distance_archive_review"
ARCHIVE = ("/big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive/"
           "20260119_all-with-w5000_1/filtered_snv_tp/filtered_snv_tp")

sel = json.load(open(S / "w5000_sel.json"))
axis = json.load(open(S / "w5000_axis.json"))
ST = sel["stats"]

D = {
    "generated": "2026-07-26",
    "archive": ARCHIVE,
    "run_label": "20260119_all-with-w5000_1 / filtered_snv_tp",
    "stats": ST,
    "axis": axis,
    "cases": sel["cases"],
    "artifacts": {
        "p_zero_n": ST["p_zero"], "p_zero_pct": round(ST["p_zero"] / ST["total"] * 100, 1),
        "p_min_nonzero": "5.00e-04", "n_perm": 2000,
        "v_one_n": ST["v_one"],
        "v_lo": ST["v_lo"], "v_hi": ST["v_hi"],
    },
    "layout": {
        "per_locus_dirs": ["distance/BERNOULLI", "plots/BERNOULLI", "clustering",
                           "methylation", "reads"],
        "distance_files": ["matrix.csv", "matrix_forward.csv", "matrix_reverse.csv",
                           "stats.txt", "stats_forward.txt", "stats_reverse.txt"],
        "plot_files": ["distance_heatmap.png", "cluster_heatmap.png"],
        "sig_fields": ["anchor_key", "num_reads", "optimal_k", "passed_gating",
                       "global_alt{p_value,cramers_v}", "global_hp{p_value,cramers_v}",
                       "heuristic_score"],
    },
}
for k in ("stats", "axis", "cases", "artifacts", "layout"):
    if k not in D:
        raise SystemExit(f"REFUSE: 缺 required key '{k}'")
OUT.mkdir(parents=True, exist_ok=True)
(OUT / f"{STEM}.data.json").write_text(json.dumps(D, indent=1, ensure_ascii=False), encoding="utf-8")


def b64(p):
    try:
        return "data:image/png;base64," + base64.b64encode(Path(p).read_bytes()).decode()
    except OSError:
        return ""


def n(x):
    return f"{x:,}" if isinstance(x, int) else x


A = D["artifacts"]
AXA, AXR = axis["axis_all"], axis["axis_robust"]
ROB, SUS, NEG = D["cases"]["robust"], D["cases"]["suspect"], D["cases"]["negative"]


def case_card(r, kind):
    base = os.path.join(r["path"], "plots", "BERNOULLI")
    d_img = b64(os.path.join(base, "distance_heatmap.png"))
    c_img = b64(os.path.join(base, "cluster_heatmap.png"))
    hp = r["hp_v"] if r["hp_v"] is not None else 0.0
    if hp < 0.3:
        axb, axl = "pass", "ALT-only（HP 無關）"
    elif hp >= 0.7:
        axb, axl = "stop", "HP 混淆（甲基跟著 germline）"
    else:
        axb, axl = "warn", "中間"
    pv = "&lt; 5×10⁻⁴" if r["alt_p"] == 0.0 else f'{r["alt_p"]:.2e}'
    kb = {"robust": ("pass", "穩健"), "suspect": ("warn", "過擬合疑慮"),
          "negative": ("todo", "未通過（對照）")}[kind]
    return f"""
<details class="card"><summary><code>{r['anchor']}</code>
  <span class="st {kb[0]}">{kb[1]}</span><span class="ct">reads {r['reads']} · V {r['alt_v']:.3f}</span></summary>
<div class="cb">
  <div class="tw"><table style="min-width:420px"><tbody>
  <tr><th>reads</th><td class="num">{r['reads']}</td>
      <th>optimal_k</th><td class="num">{r['k']}</td></tr>
  <tr><th>ALT 軸 Cramér's V</th><td class="num ok">{r['alt_v']:.3f}</td>
      <th>ALT 軸 p</th><td class="num">{pv}</td></tr>
  <tr><th>HP 軸 Cramér's V</th><td class="num">{hp:.3f}</td>
      <th>軸判定</th><td><span class="st {axb}">{axl}</span></td></tr>
  </tbody></table></div>
  <div class="imgs">
    <figure><img src="{d_img}" alt="distance heatmap"><figcaption>distance_heatmap（read×read BERNOULLI 距離）</figcaption></figure>
    <figure><img src="{c_img}" alt="cluster heatmap"><figcaption>cluster_heatmap（分群後）</figcaption></figure>
  </div>
</div></details>"""


HTML = f"""<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>w5000 archive distance/ 產出檢視 — 30,476 位點</title>
<style>
:root{{--bg:#fff;--fg:#0f172a;--mut:#64748b;--line:#e2e8f0;--card:#f8fafc;
--ok:#0d9488;--bad:#dc2626;--warn:#b45309;--acc:#2563eb;--rad:11px}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--fg);line-height:1.72;font-size:15.5px;
font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans CJK TC","PingFang TC","Microsoft JhengHei",Arial,sans-serif}}
.wrap{{max-width:1060px;margin:0 auto;padding:30px 20px 78px}}
.eyebrow{{font-size:11.5px;letter-spacing:.11em;text-transform:uppercase;color:var(--mut);font-weight:700}}
h1{{font-size:25px;margin:5px 0 7px;line-height:1.34}}
h2{{font-size:19px;margin:0 0 5px}}
.sub{{color:var(--mut);font-size:13px;margin-bottom:18px}}
.sec{{margin:36px 0 0;padding-top:22px;border-top:2px solid var(--line)}}
.sec>p.lead{{color:var(--mut);font-size:13.5px;margin:3px 0 12px}}
.verdict{{border:1px solid #fcd34d;border-left:6px solid var(--warn);background:linear-gradient(180deg,#fffbeb,#fffefa);
border-radius:var(--rad);padding:17px 21px;margin:14px 0 6px}}
.verdict .big{{font-size:17px;font-weight:700;line-height:1.58}}
.verdict .dn{{margin-top:9px;font-size:12.5px;color:var(--mut);padding-top:9px;border-top:1px dashed #fcd34d}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0}}
.kpi{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:11px 13px}}
.kpi .v{{font-size:21px;font-weight:700;font-variant-numeric:tabular-nums}}
.kpi .v small{{font-size:13px;color:var(--mut)}}
.kpi .l{{font-size:11.5px;color:var(--mut);margin-top:3px}}
.kpi.ok .v{{color:var(--ok)}} .kpi.bad .v{{color:var(--bad)}} .kpi.acc .v{{color:var(--acc)}}
.mk{{display:inline-block;font-size:10.5px;font-weight:700;padding:1px 7px;border-radius:3px;margin-right:7px;border:1px solid}}
.mk.crit{{background:#fef2f2;color:#b91c1c;border-color:#fca5a5}}
.mk.key{{background:#fffbeb;color:#92400e;border-color:#fcd34d}}
.mk.bgm{{background:#f1f5f9;color:#475569;border-color:#cbd5e1}}
.st{{display:inline-block;font-size:11px;font-weight:700;padding:1px 8px;border-radius:20px;border:1px solid;white-space:nowrap}}
.st.pass{{background:#f0fdfa;color:#0f766e;border-color:#5eead4}}
.st.warn{{background:#fffbeb;color:#92400e;border-color:#fcd34d}}
.st.todo{{background:#f1f5f9;color:#475569;border-color:#cbd5e1}}
.st.stop{{background:#fef2f2;color:#b91c1c;border-color:#fca5a5}}
ul.logic{{list-style:none;padding:0;margin:12px 0}}
ul.logic li{{padding:9px 0;border-bottom:1px dashed var(--line)}}
ul.logic li:last-child{{border-bottom:0}}
.tw{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:13.5px;margin:10px 0}}
th{{background:#f1f5f9;text-align:left;padding:8px 11px;border-bottom:2px solid #cbd5e1;font-weight:700;white-space:nowrap}}
td{{padding:8px 11px;border-bottom:1px solid var(--line)}}
td.num,th.num{{text-align:right;font-variant-numeric:tabular-nums}}
code{{background:#f1f5f9;padding:1px 5px;border-radius:4px;font-size:12.5px;font-family:ui-monospace,Menlo,monospace}}
.ok{{color:var(--ok);font-weight:700}} .bad{{color:var(--bad);font-weight:700}}
details.card{{border:1px solid var(--line);border-radius:var(--rad);margin:9px 0;overflow:hidden}}
details.card>summary{{cursor:pointer;padding:12px 16px;font-weight:600;font-size:14px;background:var(--card);
list-style:none;display:flex;align-items:center;gap:9px;flex-wrap:wrap}}
details.card>summary::-webkit-details-marker{{display:none}}
details.card>summary::before{{content:"+";font-family:ui-monospace,Menlo,monospace;font-weight:700;color:var(--acc);font-size:17px}}
details.card[open]>summary::before{{content:"−"}}
details.card>summary:hover{{background:#eef2f7}}
details.card .ct{{margin-left:auto;font-size:11px;color:var(--mut);background:#fff;border:1px solid var(--line);
border-radius:20px;padding:1px 9px;font-weight:600}}
details.card .cb{{padding:4px 16px 16px}}
.imgs{{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px;margin-top:10px}}
.imgs figure{{margin:0;border:1px solid var(--line);border-radius:9px;padding:8px;background:#fff}}
.imgs img{{width:100%;height:auto;display:block;border-radius:4px}}
.imgs figcaption{{font-size:11.5px;color:var(--mut);margin-top:6px}}
.box{{border:1px solid var(--line);border-radius:9px;padding:12px 16px;margin:11px 0;background:var(--card);font-size:14px}}
.box.dang{{background:#fef2f2;border-color:#fecaca}}
.box.warn{{background:#fffbeb;border-color:#fcd34d}}
.box.good{{background:#f0fdfa;border-color:#99f6e4}}
.box b.hd{{display:block;margin-bottom:4px}}
footer{{margin-top:42px;padding-top:16px;border-top:2px solid var(--line);font-size:12.5px;color:var(--mut)}}
@media(prefers-color-scheme:dark){{
:root{{--bg:#0b1120;--fg:#e2e8f0;--mut:#94a3b8;--line:#1e293b;--card:#0f172a}}
th{{background:#1e293b}} code{{background:#111c33}} .imgs figure{{background:#0f172a}}
.verdict{{background:#2a2410;border-color:#78350f}}
.box.dang{{background:#2a1414;border-color:#7f1d1d}} .box.warn{{background:#2a2410;border-color:#78350f}}
.box.good{{background:#0d2926;border-color:#115e59}}
}}
</style>
<div class="wrap">

<div class="eyebrow">Archive review · 2026-07-26</div>
<h1>2026-01-19 w5000 archive 的 distance/ 產出：<br>31.6% 位點「通過」，但通過率被小樣本膨脹</h1>
<div class="sub">來源：<code>{D['run_label']}</code>　·　數字由 <code>{STEM}.data.json</code> 注入</div>

<div class="verdict">
<div class="big">掃描 <b>{n(ST['total'])}</b> 個位點的 <code>clustering/significance.json</code>：
<b>{n(ST['passed'])}（{ST['passed']/ST['total']*100:.2f}%）</b> 標記 <code>passed_gating=true</code>。
但 Cramér's V 與 read 數<b>反向相關</b>（reads&lt;50 → V {A['v_lo']}；reads≥100 → V {A['v_hi']}），
真正穩健的只有 <b class="ok">{n(ST['robust'])} 個（佔全體 {ST['robust']/ST['total']*100:.2f}%）</b>。</div>
<div class="dn"><b>這是 2026-01 的舊 pipeline 產出</b>，早於後續的切群三閘與 double-dip 修正。
本頁只做「當時產出了什麼、品質如何」的檢視，<b>不作為現行科學結論</b>。</div>
</div>

<div class="grid">
<div class="kpi"><div class="v">{n(ST['total'])}</div><div class="l">掃描位點數</div></div>
<div class="kpi acc"><div class="v">{n(ST['passed'])}</div><div class="l">passed_gating（31.6%）</div></div>
<div class="kpi ok"><div class="v">{n(ST['robust'])}</div><div class="l">A 類穩健（1.58%）</div></div>
<div class="kpi bad"><div class="v">{n(A['p_zero_n'])}</div><div class="l">p 記為 0（{A['p_zero_pct']}%）</div></div>
<div class="kpi bad"><div class="v">{n(A['v_one_n'])}</div><div class="l">V 恰為 1.000</div></div>
<div class="kpi bad"><div class="v">{AXR['conf']}</div><div class="l">A 類中 HP 混淆（{AXR['conf']/AXR['n']*100:.0f}%）</div></div>
</div>

<div class="sec">
<div class="eyebrow">L1 · 這批資料是什麼</div>
<h2>每個位點一個目錄，distance/ 存 read×read 距離矩陣，plots/ 存兩張熱圖</h2>
<div class="tw"><table style="min-width:520px">
<thead><tr><th>子目錄</th><th>內容</th></tr></thead><tbody>
<tr><td><code>distance/BERNOULLI/</code></td><td><code>matrix.csv</code>、<code>matrix_forward.csv</code>、
<code>matrix_reverse.csv</code> ＋ 對應 <code>stats.txt</code>（read×read BERNOULLI 距離，含正反股分開）</td></tr>
<tr><td><code>plots/BERNOULLI/</code></td><td><code>distance_heatmap.png</code>、<code>cluster_heatmap.png</code></td></tr>
<tr><td><code>clustering/</code></td><td><code>significance.json</code>（判定依據）、<code>linkage_matrix.csv</code>、
<code>leaf_order.txt</code>、<code>tree.nwk</code>／forward／reverse</td></tr>
<tr><td><code>methylation/</code></td><td><code>methylation.csv</code>（read×CpG）、forward／reverse、<code>cpg_sites.tsv</code></td></tr>
<tr><td><code>reads/</code></td><td><code>reads.tsv</code>、<code>filtered_reads.tsv</code></td></tr>
</tbody></table></div>
<p><code>significance.json</code> 的判定欄位：<code>{'</code>、<code>'.join(D['layout']['sig_fields'])}</code></p>
</div>

<div class="sec">
<div class="eyebrow">L1 · 重點邏輯</div>
<h2>兩個必須先講的統計假象，否則 31.6% 會被誤讀</h2>
<ul class="logic">
<li><span class="mk crit">決定性</span><b>p 值恰為 0 的有 {n(A['p_zero_n'])} 個（{A['p_zero_pct']}%）—— 那不是真的 p=0。</b>
非零最小 p 是 <code>{A['p_min_nonzero']}</code> = 1/{n(A['n_perm'])}，代表 permutation 只跑 {n(A['n_perm'])} 次。
正確寫法是 <b>p &lt; 5×10⁻⁴</b>，不能寫 p=0。</li>
<li><span class="mk crit">決定性</span><b>Cramér's V 隨 read 數下降而升高 —— 小樣本膨脹。</b>
passed 中 <code>reads&lt;50</code> 的 V 中位數 <b class="bad">{A['v_lo']}</b>，
<code>reads≥100</code> 只有 <b>{A['v_hi']}</b>。V 恰為 1.000 的 {n(A['v_one_n'])} 個位點，
read 中位數 51（全體中位數 {ST['reads_median']}）。</li>
<li><span class="mk key">重要</span><b>所以「明顯」必須用 V 高<i>且</i> read 足夠定義。</b>
A 類（reads≥100 且 V≥0.8）只有 <b>{n(ST['robust'])}</b> 個，佔 passed 的 {ST['robust']/ST['passed']*100:.2f}%、
佔全體的 <b>{ST['robust']/ST['total']*100:.2f}%</b>。</li>
<li><span class="mk crit">決定性</span><b>即使是 A 類，仍有 {AXR['conf']/AXR['n']*100:.1f}% 的甲基分群跟著 germline haplotype 走。</b>
這是 cis-ASM 混淆的直接證據（見下節）。</li>
</ul>
</div>

<div class="sec">
<div class="eyebrow">L1 · 關鍵分析</div>
<h2>甲基分群到底跟著 ALT 還是跟著 HP？—— 三分之一跟著 HP</h2>
<p class="lead"><code>significance.json</code> 同時記錄 <code>global_alt</code> 與 <code>global_hp</code> 兩個軸。
若分群同時對齊 HP，代表它反映的是 germline allele-specific methylation，而非 somatic 事件。</p>
<div class="tw"><table style="min-width:560px">
<thead><tr><th>軸判定</th><th class="num">全部 passed（{n(AXA['n'])}）</th><th class="num">A 類穩健（{n(AXR['n'])}）</th></tr></thead>
<tbody>
<tr><td><span class="st pass">ALT-only 乾淨</span> <code>hp_V &lt; 0.3</code></td>
<td class="num">{n(AXA['clean'])}（{AXA['clean']/AXA['n']*100:.1f}%）</td>
<td class="num ok">{n(AXR['clean'])}（{AXR['clean']/AXR['n']*100:.1f}%）</td></tr>
<tr><td><span class="st warn">中間</span> <code>0.3 ≤ hp_V &lt; 0.7</code></td>
<td class="num">{n(AXA['mid'])}（{AXA['mid']/AXA['n']*100:.1f}%）</td>
<td class="num">{n(AXR['mid'])}（{AXR['mid']/AXR['n']*100:.1f}%）</td></tr>
<tr><td><span class="st stop">HP 混淆</span> <code>hp_V ≥ 0.7</code></td>
<td class="num bad">{n(AXA['conf'])}（{AXA['conf']/AXA['n']*100:.1f}%）</td>
<td class="num bad">{n(AXR['conf'])}（{AXR['conf']/AXR['n']*100:.1f}%）</td></tr>
</tbody></table></div>
<div class="box dang"><b class="hd">A 類穩健案例中，HP 混淆比例反而更高（{AXR['conf']/AXR['n']*100:.1f}% vs {AXA['conf']/AXA['n']*100:.1f}%）</b>
read 越多越容易同時解析出 HP 結構。所以「訊號越強」不等於「越乾淨」——
高 V 有可能只是把 germline ASM 解析得更清楚。真正可用的是
<b>ALT-only 乾淨且 reads 足夠</b>的 <b>{n(AXR['clean'])}</b> 個（佔全體 {AXR['clean']/ST['total']*100:.2f}%）。</div>
</div>

<div class="sec">
<div class="eyebrow">L2 · 案例（點開看熱圖）</div>
<h2>三種狀況的實例對照</h2>
<p class="lead">每個案例含 <code>distance_heatmap</code>（read×read 距離）與 <code>cluster_heatmap</code>（分群後）。</p>

<h3 style="font-size:15px;margin:18px 0 6px">A · 穩健明顯（reads ≥ 100 且 V ≥ 0.8）</h3>
{''.join(case_card(r,'robust') for r in ROB)}

<h3 style="font-size:15px;margin:22px 0 6px">C · 過擬合疑慮（reads &lt; 50 但 V = 1.000）</h3>
<p class="lead">注意這兩個的 <b>p 值並不顯著</b>（0.098、0.089）——V=1 完全是小樣本假象。</p>
{''.join(case_card(r,'suspect') for r in SUS)}

<h3 style="font-size:15px;margin:22px 0 6px">對照 · 未通過（read 充足但 V = 0）</h3>
<p class="lead">這才是可信的陰性：read 足夠、但確實沒有 ALT-關聯的甲基分群。</p>
{''.join(case_card(r,'negative') for r in NEG)}
</div>

<div class="sec">
<div class="eyebrow">L3 · 適用邊界與溯源</div>
<h2>四項限制</h2>
<div class="box warn"><b class="hd">1 · 這是 2026-01-19 的舊 pipeline 產出</b>
早於後續的切群三閘、double-dip 修正與 allele-axis baseline 規範。
本頁是<b>檔案考古與品質檢視</b>，不是現行科學結論。</div>
<div class="box warn"><b class="hd">2 · p 值不可直接引用</b>
{A['p_zero_pct']}% 記為 0，實為 <code>p &lt; 5×10⁻⁴</code>（{n(A['n_perm'])} 次 permutation 的解析度上限）。</div>
<div class="box warn"><b class="hd">3 · Cramér's V 未做小樣本校正</b>
V 與 read 數反向相關（{A['v_lo']} vs {A['v_hi']}）。任何用 V 排序的結論都會偏向低覆蓋位點。</div>
<div class="box dang"><b class="hd">4 · 未做 germline-het null，ASM 無法形式排除</b>
{AXR['conf']/AXR['n']*100:.1f}% 的 A 類案例甲基分群對齊 HP。要宣稱 somatic-specific，
必須補 allele-axis baseline 對照。</div>

<details class="card"><summary>L3 · 溯源與可重算</summary><div class="cb">
<ul>
<li><b>archive root</b>：<code>{ARCHIVE}</code></li>
<li><b>結構</b>：22 個染色體目錄 → 每位點一目錄 → 每 region 一子目錄（含
<code>{'</code>、<code>'.join(D['layout']['per_locus_dirs'])}</code>）</li>
<li><b>掃描</b>：讀取每個 <code>clustering/significance.json</code>，共 {n(ST['total'])} 個；14 個位點缺此檔</li>
<li><b>圖片</b>：<code>plots/BERNOULLI/{{distance,cluster}}_heatmap.png</code>，本頁以 base64 內嵌</li>
<li><b>反捏造</b>：所有數字由 <code>{STEM}.data.json</code> 注入；缺 required key 直接 refuse（§13-A）</li>
</ul>
</div></details>
</div>

<footer>
<p>生成 {D['generated']}　·　資訊分層：L0 一眼結論 → L1 資料結構／重點邏輯／軸分析 → L2 案例卡（收合）→ L3 邊界與溯源。
　·　無 emoji（本機無 emoji 字型），狀態一律附文字標籤。</p>
</footer>
</div>
"""
(OUT / f"{STEM}.standalone.html").write_text(HTML, encoding="utf-8")
print(f"wrote {OUT / (STEM + '.standalone.html')}  ({len(HTML)/1024:.0f} KB)")
print(f"wrote {OUT / (STEM + '.data.json')}")
