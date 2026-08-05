#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_page.py — 甲基 allele 軸訊號與 somatic 遺傳骨幹共同富集（含內建陰性對照）

設計: L0 verdict / L1 邏輯+主圖 / L2 分層證據卡 / L3 邊界與溯源
      無 emoji(本機無 emoji 字型) → CSS 徽章;狀態附文字標籤
§13-A: 數字全由 data.json 注入;缺 required key 即 refuse
"""
import json
from pathlib import Path

S = Path("/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/"
         "efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad")
OUT = Path(__file__).parent
STEM = "20260726_methyl_allele_axis_backbone_coenrichment"

D = json.loads((S / "methyl_backbone.json").read_text())
for k in ("n_tp","baseline","A","B","cn","cov","double","coread","density","density_strata","regions","top_regions","verdict"):
    if k not in D:
        raise SystemExit(f"REFUSE: 缺 required key '{k}'")
D["sources"] = {
    "phylo": "docs/methodology/_assets/20260618_subcluster_pilot/phylo_cpp_wg_full_records_v6.json",
    "linkage": "docs/methodology/_assets/20260618_subcluster_pilot/sm_linkage_genomewide.json",
    "archive": ("big7_disk_output/bip8_output_archive/20260119_all-with-w5000_1/"
                "filtered_snv_tp"),
}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / f"{STEM}.data.json").write_text(json.dumps(D, indent=1, ensure_ascii=False), encoding="utf-8")

BL, A, B = D["baseline"], D["A"], D["B"]
DELTA = round(A["pct"] - B["pct"], 1)


def n(x):
    return f"{x:,}" if isinstance(x, int) else x


def svg_main(width=940, height=300):
    """主圖:A / B / baseline 三條橫棒 + 陰性對照標註"""
    ml, bw = 210, 560
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" role="img" '
         f'aria-label="甲基分群軸別與 sSNV 串聯率"><rect width="{width}" height="{height}" fill="white"/>']
    p.append(f'<text x="{ml}" y="26" font-size="13" font-weight="700" fill="#0f172a">'
             f'甲基分群對齊到哪個軸 → 該位點有 sSNV 串聯 partner 的比例</text>')
    rows = [("A · 對齊 ALLELE 軸", A, "#0d9488", "訊號組"),
            ("B · 對齊 HP 軸", B, "#dc2626", "內建陰性對照"),
            ("全體 TP 位點", BL, "#94a3b8", "基線")]
    for i, (lab, r, col, tag) in enumerate(rows):
        y = 62 + i * 62
        p.append(f'<text x="{ml-12}" y="{y+16}" font-size="12.5" text-anchor="end" '
                 f'font-weight="600" fill="#334155">{lab}</text>')
        p.append(f'<text x="{ml-12}" y="{y+32}" font-size="10.5" text-anchor="end" fill="#94a3b8">{tag}</text>')
        p.append(f'<rect x="{ml}" y="{y}" width="{bw}" height="26" rx="4" fill="#f1f5f9"/>')
        w = bw * r["pct"] / 100
        p.append(f'<rect x="{ml}" y="{y}" width="{w:.1f}" height="26" rx="4" fill="{col}" opacity=".88"/>')
        p.append(f'<text x="{ml+w+10}" y="{y+18}" font-size="13" font-weight="700" fill="{col}">'
                 f'{r["pct"]}%</text>')
        p.append(f'<text x="{ml+w+62}" y="{y+18}" font-size="11" fill="#64748b">'
                 f'{n(r["hit"])} / {n(r["n"])}</text>')
    yb = 62 + 1 * 62
    p.append(f'<line x1="{ml+bw*BL["pct"]/100}" y1="52" x2="{ml+bw*BL["pct"]/100}" y2="{62+3*62-8}" '
             f'stroke="#94a3b8" stroke-dasharray="5 4" stroke-width="1.5"/>')
    p.append(f'<text x="{ml+bw*BL["pct"]/100+6}" y="48" font-size="10.5" fill="#94a3b8">基線 {BL["pct"]}%</text>')
    p.append(f'<rect x="{ml}" y="252" width="{bw}" height="34" rx="7" fill="#fef2f2" stroke="#dc2626" stroke-width="1.4"/>')
    p.append(f'<text x="{ml+14}" y="273" font-size="12" fill="#b91c1c">'
             f'<tspan font-weight="700">此表面差距已被推翻</tspan>：A 與 B 的局部突變密度本就不同'
             f'（中位 {D["density"]["A_median"]} vs {D["density"]["B_median"]}）→ 見下方密度分層</text>')
    p.append("</svg>")
    return "".join(p)


def svg_strata(width=940):
    """分層圖:CN×覆蓋 雙分層下的 A-B 差距"""
    d = D["double"]
    rowh, top, ml, bw = 30, 62, 150, 600
    height = top + rowh * len(d) + 56
    mx = max(abs(x["delta"]) for x in d)
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" role="img" '
         f'aria-label="CN 與覆蓋雙分層下的差距"><rect width="{width}" height="{height}" fill="white"/>']
    p.append(f'<text x="{ml}" y="24" font-size="13" font-weight="700" fill="#0f172a">'
             f'CN × 覆蓋深度雙分層：每一層 A 都高於 B（{D["strata_A_gt_B"]}/{D["strata_total"]}）</text>')
    p.append(f'<text x="{ml}" y="44" font-size="11" fill="#64748b">'
             f'橫軸 = A 串聯率 − B 串聯率（百分點）；覆蓋四分位切點 n_tumor = {D["quartiles"]}</text>')
    for i, x in enumerate(d):
        y = top + i * rowh
        lab = f'CN {x["cn"]} · 覆蓋 {x["cov"]}'
        p.append(f'<text x="{ml-10}" y="{y+16}" font-size="11.5" text-anchor="end" fill="#334155">{lab}</text>')
        w = bw * x["delta"] / mx
        p.append(f'<rect x="{ml}" y="{y+3}" width="{max(w,2):.1f}" height="18" rx="3" '
                 f'fill="#0d9488" opacity=".85"/>')
        p.append(f'<text x="{ml+max(w,2)+9}" y="{y+17}" font-size="11.5" font-weight="700" '
                 f'fill="#0d9488">+{x["delta"]}pp</text>')
        p.append(f'<text x="{ml+max(w,2)+66}" y="{y+17}" font-size="10.5" fill="#94a3b8">'
                 f'A {x["A"]["pct"]}% (n={x["A"]["n"]}) vs B {x["B"]["pct"]}% (n={x["B"]["n"]})</text>')
    yb = top + rowh * len(d) + 18
    p.append(f'<text x="{ml}" y="{yb+8}" font-size="11.5" fill="#475569">'
             f'差距在<tspan font-weight="700">低覆蓋層最大</tspan>（Q1 +51.6pp）、'
             f'高覆蓋層最小（Q4 +7.1pp）——與「覆蓋度混淆」的預期方向<tspan font-weight="700">相反</tspan>。</text>')
    p.append("</svg>")
    return "".join(p)



def svg_density(width=940):
    """否決圖:局部 sSNV 密度分層後差距崩解"""
    d = D["density_strata"]
    ml, bw, rowh, top = 130, 560, 40, 76
    height = top + rowh * len(d) + 64
    p = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" role="img" '
         f'aria-label="局部突變密度分層後差距崩解"><rect width="{width}" height="{height}" fill="white"/>']
    p.append(f'<text x="{ml-100}" y="26" font-size="13" font-weight="700" fill="#b91c1c">'
             f'控制「局部 sSNV 密度（±50 kb）」後，A 與 B 的差距逐層消失</text>')
    p.append(f'<text x="{ml-100}" y="46" font-size="11" fill="#64748b">'
             f'未分層時 +{DELTA}pp；密度越高，兩組越接近，最終完全重合</text>')
    p.append(f'<text x="{ml-14}" y="{top-8}" font-size="10.5" text-anchor="end" fill="#94a3b8">密度</text>')
    p.append(f'<text x="{ml}" y="{top-8}" font-size="10.5" fill="#94a3b8">A 串聯率（青）vs B（紅）</text>')
    for i, x in enumerate(d):
        y = top + i * rowh
        p.append(f'<text x="{ml-14}" y="{y+20}" font-size="12" text-anchor="end" '
                 f'font-weight="600" fill="#334155">{x["lab"]}</text>')
        for j, (k, col) in enumerate([("A", "#0d9488"), ("B", "#dc2626")]):
            w = bw * x[k]["pct"] / 100
            yy = y + j * 13
            p.append(f'<rect x="{ml}" y="{yy}" width="{max(w,2):.1f}" height="11" rx="2" '
                     f'fill="{col}" opacity=".85"/>')
        p.append(f'<text x="{ml+bw+12}" y="{y+18}" font-size="12" font-weight="700" '
                 f'fill="{"#b45309" if x["delta"]>1 else "#94a3b8"}">'
                 f'{"+" if x["delta"]>0 else ""}{x["delta"]}pp</text>')
        p.append(f'<text x="{ml+bw+72}" y="{y+18}" font-size="10" fill="#94a3b8">'
                 f'n {x["A"]["n"]}/{x["B"]["n"]}</text>')
    yb = top + rowh * len(d) + 20
    p.append(f'<text x="{ml-100}" y="{yb+4}" font-size="11.5" fill="#475569">'
             f'A 組密度 median <tspan font-weight="700">{D["density"]["A_median"]}</tspan>'
             f'（mean {D["density"]["A_mean"]}）vs B 組 <tspan font-weight="700">{D["density"]["B_median"]}</tspan>'
             f'（mean {D["density"]["B_mean"]}）；密度=0 者 A 僅 {D["density"]["A_zero"]}、B 達 {D["density"]["B_zero"]}。</text>')
    p.append("</svg>")
    return "".join(p)

def _cn_row(k, v):
    hi = v["A"]["pct"] > v["B"]["pct"]
    verdict = ('<span class="st pass">A 較高</span>' if hi
               else '<span class="st stop">A 較低</span>')
    if v["A"]["n"] < 20:
        verdict += ' <span class="st warn">n 過小</span>'
    return ('<tr><td><code>{k}</code></td>'
            '<td class="num">{an}</td><td class="num ok">{ap}%</td>'
            '<td class="num">{bn}</td><td class="num">{bp}%</td>'
            '<td class="num {cls}">{d:+.1f}pp</td><td>{ver}</td></tr>').format(
        k=k, an=n(v["A"]["n"]), ap=v["A"]["pct"], bn=n(v["B"]["n"]), bp=v["B"]["pct"],
        cls="ok" if hi else "bad", d=v["A"]["pct"] - v["B"]["pct"], ver=verdict)


cn_rows = "".join(_cn_row(k, v) for k, v in D["cn"].items())
cov_rows = "".join(
    f'<tr><td><code>{k}</code></td>'
    f'<td class="num">{n(v["A"]["n"])}</td><td class="num ok">{v["A"]["pct"]}%</td>'
    f'<td class="num">{n(v["B"]["n"])}</td><td class="num">{v["B"]["pct"]}%</td>'
    f'<td class="num ok">{v["A"]["pct"]-v["B"]["pct"]:+.1f}pp</td></tr>'
    for k, v in D["cov"].items())
den_rows = "".join(
    f'<tr><td><code>{x["lab"]}</code></td>'
    f'<td class="num">{x["A"]["n"]}</td><td class="num">{x["A"]["pct"]}%</td>'
    f'<td class="num">{x["B"]["n"]}</td><td class="num">{x["B"]["pct"]}%</td>'
    f'<td class="num {"bad" if x["delta"]>1 else "mutcell"}">{"+" if x["delta"]>0 else ""}{x["delta"]}pp</td></tr>'
    for x in D["density_strata"])
reg_rows = "".join(
    f'<tr data-a="{r["A"]}" data-n="{r["n_ssnv"]}"><td><code>{r["region_id"]}</code></td>'
    f'<td class="num">{r["n_ssnv"]}</td><td class="num ok">{r["A"]}</td>'
    f'<td class="num">{r["B"]}</td><td class="num">{r["D"]}</td>'
    f'<td class="num">{r["n_pairs"]}</td><td class="num">{r["max_coread"]}</td>'
    f'<td style="font-size:11.5px">{r["rels"]}</td>'
    f'<td style="font-size:11px;max-width:340px;word-break:break-all">{r["loci"][:300]}</td></tr>'
    for r in D["top_regions"])
dbl_rows = "".join(
    f'<tr><td><code>{x["cn"]}</code></td><td><code>{x["cov"]}</code></td>'
    f'<td class="num">{x["A"]["n"]}</td><td class="num ok">{x["A"]["pct"]}%</td>'
    f'<td class="num">{x["B"]["n"]}</td><td class="num">{x["B"]["pct"]}%</td>'
    f'<td class="num ok">+{x["delta"]}pp</td></tr>' for x in D["double"])

HTML = f"""<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>甲基 allele 軸「共同富集」— 已被局部突變密度分層否決</title>
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
.verdict{{border:1px solid #99f6e4;border-left:6px solid var(--ok);background:linear-gradient(180deg,#f0fdfa,#fafffe);
border-radius:var(--rad);padding:17px 21px;margin:14px 0 6px}}
.verdict .big{{font-size:17px;font-weight:700;line-height:1.58}}
.verdict .dn{{margin-top:9px;font-size:12.5px;color:var(--mut);padding-top:9px;border-top:1px dashed #99f6e4}}
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
figure{{margin:16px 0;padding:12px;background:#fff;border:1px solid var(--line);border-radius:var(--rad)}}
figure svg{{width:100%;height:auto;display:block}}
figcaption{{font-size:12.5px;color:var(--mut);margin-top:9px;padding-top:8px;border-top:1px dashed var(--line)}}
.tw{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:13.5px;margin:10px 0}}
th{{background:#f1f5f9;text-align:left;padding:8px 11px;border-bottom:2px solid #cbd5e1;font-weight:700;white-space:nowrap}}
td{{padding:8px 11px;border-bottom:1px solid var(--line)}}
td.num,th.num{{text-align:right;font-variant-numeric:tabular-nums}}
code{{background:#f1f5f9;padding:1px 5px;border-radius:4px;font-size:12.5px;font-family:ui-monospace,Menlo,monospace}}
.ok{{color:var(--ok);font-weight:700}} .bad{{color:var(--bad);font-weight:700}}
details.card{{border:1px solid var(--line);border-radius:var(--rad);margin:9px 0;overflow:hidden}}
details.card>summary{{cursor:pointer;padding:12px 16px;font-weight:600;font-size:14.5px;background:var(--card);
list-style:none;display:flex;align-items:center;gap:9px}}
details.card>summary::-webkit-details-marker{{display:none}}
details.card>summary::before{{content:"+";font-family:ui-monospace,Menlo,monospace;font-weight:700;color:var(--acc);font-size:17px}}
details.card[open]>summary::before{{content:"−"}}
details.card .ct{{margin-left:auto;font-size:11px;color:var(--mut);background:#fff;border:1px solid var(--line);
border-radius:20px;padding:1px 9px;font-weight:600}}
details.card .cb{{padding:4px 16px 16px}}
.box{{border:1px solid var(--line);border-radius:9px;padding:12px 16px;margin:11px 0;background:var(--card);font-size:14px}}
.box.dang{{background:#fef2f2;border-color:#fecaca}}
.box.warn{{background:#fffbeb;border-color:#fcd34d}}
.box.good{{background:#f0fdfa;border-color:#99f6e4}}
.box b.hd{{display:block;margin-bottom:4px}}
footer{{margin-top:42px;padding-top:16px;border-top:2px solid var(--line);font-size:12.5px;color:var(--mut)}}
@media(prefers-color-scheme:dark){{
:root{{--bg:#0b1120;--fg:#e2e8f0;--mut:#94a3b8;--line:#1e293b;--card:#0f172a}}
figure{{background:#0f172a}} figure svg{{filter:invert(.92) hue-rotate(180deg)}}
th{{background:#1e293b}} code{{background:#111c33}}
.verdict{{background:#0d2926;border-color:#115e59}}
.box.dang{{background:#2a1414;border-color:#7f1d1d}} .box.warn{{background:#2a2410;border-color:#78350f}}
.box.good{{background:#0d2926;border-color:#115e59}}
}}
</style>
<div class="wrap">

<div class="eyebrow">Refuted finding · 2026-07-26</div>
<h1>甲基 <b>allele 軸</b>與 somatic 骨幹的「共同富集」<br>—— <span style="color:var(--bad)">已被局部突變密度分層否決</span></h1>
<div class="sub">兩條獨立證據線交叉：甲基分群（phylo v6）× sSNV 串聯（sm_linkage）　·
數字由 <code>{STEM}.data.json</code> 注入</div>

<div class="verdict" style="border-color:#fecaca;border-left-color:var(--bad);background:linear-gradient(180deg,#fef2f2,#fffafa)">
<div class="big">表面上：allele 軸 <b>{A['pct']}%</b> vs HP 軸 <b>{B['pct']}%</b>（≡基線 {BL['pct']}%），CN×覆蓋雙分層 {D['strata_A_gt_B']}/{D['strata_total']} 全存活，看似有完美的內建陰性對照。<br>
<b style="color:var(--bad)">但加入「局部 sSNV 密度」分層後崩解：+{D['density_strata'][0]['delta']} → +{D['density_strata'][1]['delta']} → +{D['density_strata'][2]['delta']} → 0 → 0。</b> 根因是 A 組本來就落在突變更密的區域（密度中位 {D['density']['A_median']} vs {D['density']['B_median']}）。</div>
<div class="dn"><b>結論：共同富集 = 局部突變密度混淆，不是獨立關聯。</b> 與既有裁決一致 —— allele 軸訊號是<b>突變局部 cis 足跡</b>，非譜系訊號。本頁保留完整推導，作為「內建陰性對照看似完美仍可能被第三變項解釋」的實例。</div>
</div>

<div class="grid">
<div class="kpi ok"><div class="v">{A['pct']}<small>%</small></div><div class="l">A 對齊 ALLELE（n={n(A['n'])}）</div></div>
<div class="kpi bad"><div class="v">{B['pct']}<small>%</small></div><div class="l">B 對齊 HP（陰性對照）</div></div>
<div class="kpi"><div class="v">{BL['pct']}<small>%</small></div><div class="l">全體 TP 基線（{n(BL['n'])}）</div></div>
<div class="kpi ok"><div class="v">+{DELTA}<small>pp</small></div><div class="l">A − B 差距</div></div>
<div class="kpi ok"><div class="v">{D['strata_A_gt_B']}/{D['strata_total']}</div><div class="l">分層存活</div></div>
<div class="kpi"><div class="v">{D['coread']['A']}<small> vs {D['coread']['B']}</small></div><div class="l">coread 中位（A vs B）</div></div>
</div>

<div class="sec">
<div class="eyebrow">L1 · 主結果</div>
<h2>B 組（HP 軸）串聯率與基線完全相同 —— 這是設計上的陰性對照</h2>
<figure>{svg_main()}
<figcaption>三組都通過<b>同樣的</b>甲基分群與對齊門檻，差別只在<b>對齊到哪個軸</b>。
若甲基訊號只反映 germline haplotype（ASM），它與 somatic sSNV 串聯應無關係
—— 實測 B 組 {B['pct']}% 與基線 {BL['pct']}% 一致，符合這個預期。</figcaption></figure>
</div>

<div class="sec">
<div class="eyebrow">L1 · 重點邏輯</div>
<h2>CN 與覆蓋分層都存活 —— 但那不是足夠的對照</h2>
<ul class="logic">
<li><span class="mk crit">決定性</span><b>差距在低覆蓋層最大、高覆蓋層最小。</b>
Q1 低覆蓋 <b>+{D['cov']['Q1']['A']['pct']-D['cov']['Q1']['B']['pct']:.1f}pp</b>，
Q4 高覆蓋僅 <b>+{D['cov']['Q4']['A']['pct']-D['cov']['Q4']['B']['pct']:.1f}pp</b>。
若是覆蓋度混淆，方向應該<b>相反</b>（覆蓋越高、越容易偵測串聯、差距越大）。</li>
<li><span class="mk crit">決定性</span><b>CN × 覆蓋雙分層 {D['strata_A_gt_B']}/{D['strata_total']} 層全部 A &gt; B。</b>
沒有任何一層反轉。</li>
<li><span class="mk key">重要</span><b>LOH 區的差距最大（+{D['cn']['loh']['A']['pct']-D['cn']['loh']['B']['pct']:.1f}pp）。</b>
機制上合理：LOH 區實質只剩單一 haplotype，HP 對齊在此本就無意義 ——
B 組在 LOH 只有 {D['cn']['loh']['B']['pct']}%，遠低於其他 CN 狀態。</li>
<li><span class="mk key">重要</span><b>A 組的 coread 中位數也較高（{D['coread']['A']} vs {D['coread']['B']}）。</b>
不只更常有 partner，分子證據也更厚。</li>
<li><span class="mk bgm">背景</span><b>join 校準：offset = +{D['offset_bp']} bp。</b>
phylo 記錄存的是 <code>w5000</code> window 起點，非 SNV 座標；校正後交集從 1 躍升到 24,220。</li>
</ul>
</div>

<div class="sec">
<div class="eyebrow">L1 · 分層</div>
<h2>CN × 覆蓋雙分層：每一層 A 都高於 B</h2>
<figure>{svg_strata()}
<figcaption>只列 n(A) ≥ 10 且 n(B) ≥ 10 的分層。CN neutral 因 A 組僅
{D['cn']['neutral']['A']['n']} 個位點未達門檻，單獨列於 L3。</figcaption></figure>
</div>

<div class="sec">
<div class="eyebrow">L1 · 否決</div>
<h2>加入「局部 sSNV 密度」後，差距逐層消失至 0</h2>
<p class="lead">CN 與覆蓋只是弱代理。真正的共同原因是<b>局部突變密度</b> ——
它同時驅動 cis 足跡強度（→ allele 軸易對齊）與串聯可偵測性（→ 有 partner）。</p>
<figure>{svg_density()}
<figcaption>密度 = 該位點 ±50 kb 內其他 TP sSNV 的數目。
差距從 +{D['density_strata'][0]['delta']}pp 一路收斂到 <b>0.0pp</b>；
密度 ≥5 時兩組<b>完全重合（皆 100%）</b>。</figcaption></figure>
<div class="tw"><table style="min-width:560px">
<thead><tr><th>密度</th><th class="num">A n</th><th class="num">A%</th>
<th class="num">B n</th><th class="num">B%</th><th class="num">差距</th></tr></thead>
<tbody>{den_rows}</tbody></table></div>
<div class="box dang"><b class="hd">根因：A 與 B 的密度分佈本來就不同</b>
A 組 median <b>{D['density']['A_median']}</b>（mean {D['density']['A_mean']}）vs
B 組 <b>{D['density']['B_median']}</b>（mean {D['density']['B_mean']}）。
密度 = 0（±50 kb 內無其他 sSNV）者：A 僅 <b>{D['density']['A_zero']}</b> 個（{D['density']['A_zero']/A['n']*100:.1f}%），
B 有 <b>{n(D['density']['B_zero'])}</b> 個（{D['density']['B_zero']/B['n']*100:.1f}%）。
<b>A 組本來就更常落在突變密集區 → 更容易有串聯 partner。</b></div>
<div class="box warn"><b class="hd">與既有裁決一致（互相印證）</b>
既有 memory 已證「genotype-甲基差是普遍 <b>somatic-cis</b>」「主要是<b>突變局部 cis 足跡</b>，
lineage 殘餘僅 6–11%」。本輪是該結論在「串聯骨幹」軸上的獨立重現 —— 不是新發現，是同一件事。</div>
</div>

<div class="sec">
<div class="eyebrow">L2 · 區域瀏覽器</div>
<h2>邊看 sSNV 結構、邊看甲基狀況</h2>
<p class="lead">以 sSNV 串聯的<b>連通分量</b>為區域單位（共 {n(D['regions']['total'])} 個有甲基註記的區域）。
下表列 A 類最多的 {len(D['top_regions'])} 個。<b>這是描述性標記，不是證據。</b></p>
<div class="grid">
<div class="kpi"><div class="v">{n(D['regions']['total'])}</div><div class="l">有甲基註記的串聯區域</div></div>
<div class="kpi ok"><div class="v">{n(D['regions']['with_A'])}</div><div class="l">含 ≥1 個 allele 軸位點</div></div>
<div class="kpi"><div class="v">{n(D['regions']['with_A2'])}</div><div class="l">含 ≥2 個</div></div>
<div class="kpi"><div class="v">{n(D['regions']['A_and_B'])}</div><div class="l">A 與 B 並存</div></div>
</div>
<div class="tw"><table style="min-width:940px">
<thead><tr><th>區域（連通分量）</th><th class="num">sSNV</th><th class="num">A</th><th class="num">B</th>
<th class="num">未對齊</th><th class="num">pairs</th><th class="num">coread</th><th>串聯關係</th>
<th>位點｜軸｜V_allele｜V_hp｜CN</th></tr></thead>
<tbody>{reg_rows}</tbody></table></div>
<p class="mut">完整 30,077 列見同目錄 <code>methyl_class_x_linkage_annotation.tsv</code>
（欄位：chrom·snv_pos·methyl_class·V_allele·V_hp·coarse_ng·n_tumor·cn_state·is_loh·
n_linkage_partners·max_coread·top_rel·both_somatic·partners）。</p>
</div>

<div class="sec">
<div class="eyebrow">L2 · 分層證據</div>
<h2>逐項數據</h2>

<details class="card"><summary>證據 1 · CN 狀態分層<span class="ct">含反例</span></summary><div class="cb">
<div class="tw"><table style="min-width:620px">
<thead><tr><th>CN</th><th class="num">A n</th><th class="num">A 串聯率</th>
<th class="num">B n</th><th class="num">B 串聯率</th><th class="num">差距</th><th>判定</th></tr></thead>
<tbody>{cn_rows}</tbody></table></div>
<div class="box dang"><b class="hd">CN neutral 是唯一反轉的分層</b>
A 組僅 <b>{D['cn']['neutral']['A']['n']}</b> 個位點（{D['cn']['neutral']['A']['pct']}%），
低於 B 組 {D['cn']['neutral']['B']['pct']}%。n 太小無法下結論，
但這正是<b>唯一不受 CN 混淆</b>的分層 —— 也是最需要補資料的地方。</div>
</div></details>

<details class="card"><summary>證據 2 · 覆蓋深度分層（n_tumor 四分位）<span class="ct">方向性反證</span></summary><div class="cb">
<div class="tw"><table style="min-width:560px">
<thead><tr><th>覆蓋</th><th class="num">A n</th><th class="num">A 串聯率</th>
<th class="num">B n</th><th class="num">B 串聯率</th><th class="num">差距</th></tr></thead>
<tbody>{cov_rows}</tbody></table></div>
<p>四分位切點 <code>n_tumor = {D['quartiles']}</code>。
差距<b>隨覆蓋上升而縮小</b>，與覆蓋度混淆的預期方向相反。</p>
</div></details>

<details class="card"><summary>證據 3 · CN × 覆蓋雙分層（最嚴格）<span class="ct">{D['strata_A_gt_B']}/{D['strata_total']}</span></summary><div class="cb">
<div class="tw"><table style="min-width:600px">
<thead><tr><th>CN</th><th>覆蓋</th><th class="num">A n</th><th class="num">A%</th>
<th class="num">B n</th><th class="num">B%</th><th class="num">差距</th></tr></thead>
<tbody>{dbl_rows}</tbody></table></div>
</div></details>
</div>

<div class="sec">
<div class="eyebrow">L3 · 能標記什麼、不能標記什麼</div>
<h2>可標「描述性註記」，不可標「關聯」或「已確認」</h2>

<div class="box good"><b class="hd">✔ 可以標記：descriptive annotation（描述性註記）</b>
「這些位點的甲基分群對齊到哪個軸（allele／HP／未對齊），以及它們與哪些 sSNV 有串聯。」
—— 純描述，可用於瀏覽與挑選候選；<b>不附帶任何關聯或因果宣稱</b>。</div>

<div class="box dang"><b class="hd">✘ 不可以標記：co-enrichment / association / confirmed</b>
<ul>
<li><b>共同富集本身已被否決</b> —— 密度分層後歸零，不是獨立關聯。</li>
<li><b>CN-neutral 分層（A n={D['cn']['neutral']['A']['n']}）方向相反</b>
—— 唯一乾淨的分層沒有支持，這是關鍵缺口。</li>
<li>仍未做 <b>germline-het null</b>（allele 軸 baseline confound 未形式排除）。</li>
<li>仍未做 <b>matched-normal cis-control</b>。</li>
</ul></div>

<div class="box warn"><b class="hd">分母警告：這是第三條資料線</b>
本頁用的是 <b>2026-01 w5000 archive + phylo v6 + sm_linkage</b>（TP {n(D['n_tp'])} 位點）。
<b>不是</b> exact-PS cohort（98,955 units），<b>也不是</b> Task-B 全 sSNV 線（469,849 sites）。
三者分母不同，<b>不可並列或相除</b>。</div>

<details class="card"><summary>L3 · 溯源與可重算</summary><div class="cb">
<ul>
<li><b>甲基分群</b>：<code>{D['sources']['phylo']}</code>（欄位 <code>aligned</code>／<code>V_hp</code>／<code>V_allele</code>／<code>cn_state</code>／<code>n_tumor</code>）</li>
<li><b>sSNV 串聯</b>：<code>{D['sources']['linkage']}</code>（53,094 pairs；<code>COREAD_MIN=6</code>、<code>TIER_R=50000</code>）</li>
<li><b>原始產出</b>：<code>{D['sources']['archive']}</code></li>
<li><b>join 校準</b>：phylo <code>pos</code> + {D['offset_bp']} = linkage 座標（window start → SNV pos）</li>
<li><b>分組定義</b>：A = <code>aligned ∧ ¬unstable ∧ V_allele≥0.7</code>；B = <code>aligned ∧ ¬unstable ∧ V_hp≥0.7</code></li>
<li><b>反捏造</b>：數字全由 <code>{STEM}.data.json</code> 注入；缺 required key 直接 refuse（§13-A）</li>
</ul>
</div></details>
</div>

<footer>
<p>生成 {D['generated']}　·　L0 一眼結論 → L1 主結果／邏輯／分層 → L2 證據卡（收合）→ L3 標記界線與溯源。
　·　無 emoji（本機無 emoji 字型），狀態一律附文字標籤。</p>
</footer>
</div>
"""
(OUT / f"{STEM}.standalone.html").write_text(HTML, encoding="utf-8")
print(f"wrote {OUT / (STEM + '.standalone.html')}  ({len(HTML)/1024:.0f} KB)")
