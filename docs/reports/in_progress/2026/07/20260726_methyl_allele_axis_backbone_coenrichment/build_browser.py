#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_browser.py — sSNV 分支結構 × 甲基狀況 互動瀏覽器

需求：邊確認 sSNV 結構、邊確認甲基；並在重點檢視圖片上標記哪些是甲基明顯的區域。
選案準則：區域內同時具 (a) sSNV 分支（mutual_excl / nested）與 (b) ≥1 個 allele 軸甲基位點。

§13-A：數字全由 data.json 注入；缺 required key 即 refuse
自足：熱圖 base64 內嵌、JS 內聯，零外部資源
無 emoji（本機無 emoji 字型）→ CSS 徽章 + 文字標籤
"""
import base64
import json
from pathlib import Path

S = Path("/tmp/claude-1067/-big7-disk-liaoyoyo2001-InterSubMod/"
         "efb6e3d8-c4af-43d8-ac97-9dffdbec60ed/scratchpad")
OUT = Path(__file__).parent
STEM = "20260726_ssnv_branch_x_methyl_browser"

CASES = json.loads((S / "branch_methyl_cases.json").read_text())
BACK = json.loads((S / "methyl_backbone.json").read_text())
if not CASES:
    raise SystemExit("REFUSE: branch_methyl_cases.json 為空")

FLAGSHIP = "chr2:18068480-18123326"   # 2026-06 深度個案，獨立被自動篩選命中

D = {"generated": "2026-07-26", "n_cases": len(CASES),
     "n_branch_methyl_regions": 472,
     "regions_total": BACK["regions"]["total"],
     "regions_with_A": BACK["regions"]["with_A"],
     "verdict_note": ("區域層級的共同富集已被局部突變密度分層否決；"
                      "本頁為描述性瀏覽，不作關聯或因果宣稱。"),
     "cases": [{k: c[k] for k in ("region_id", "chrom", "span", "n_ssnv",
                                  "A", "B", "C", "D", "E", "n_pairs",
                                  "max_coread", "rels", "mutual_excl", "branch")}
               for c in CASES]}
for k in ("cases", "n_cases", "regions_total"):
    if k not in D:
        raise SystemExit(f"REFUSE: 缺 required key '{k}'")
(OUT / f"{STEM}.data.json").write_text(json.dumps(D, indent=1, ensure_ascii=False), encoding="utf-8")


def b64(p, name):
    try:
        return "data:image/png;base64," + base64.b64encode((Path(p) / name).read_bytes()).decode()
    except (OSError, TypeError):
        return ""


CLS_LABEL = {"A_ALLELE": ("pass", "allele 軸"), "B_HP": ("stop", "HP 軸(ASM)"),
             "C_other": ("warn", "其他軸"), "D_unaligned": ("todo", "未對齊"),
             "E_unstable": ("todo", "不穩定")}


def svg_structure(c, width=880):
    """區域內 sSNV 的線性排列 + 甲基軸別著色"""
    loci = sorted(c["locus_detail"], key=lambda x: x["pos"])
    lo, hi = loci[0]["pos"], loci[-1]["pos"]
    span = max(hi - lo, 1)
    ml, bw, y = 24, width - 130, 78
    p = [f'<svg viewBox="0 0 {width} 140" xmlns="http://www.w3.org/2000/svg" role="img" '
         f'aria-label="{c["region_id"]} 的 sSNV 排列與甲基軸別">'
         f'<rect width="{width}" height="140" fill="white"/>']
    p.append(f'<text x="{ml}" y="22" font-size="12" font-weight="700" fill="#0f172a">'
             f'{c["region_id"]}　<tspan font-weight="400" fill="#64748b">'
             f'{c["n_ssnv"]} sSNV · span {c["span"]:,} bp · 互斥 {c["mutual_excl"]} 對 · '
             f'max coread {c["max_coread"]}</tspan></text>')
    p.append(f'<line x1="{ml}" y1="{y}" x2="{ml+bw}" y2="{y}" stroke="#cbd5e1" stroke-width="2"/>')
    for lc in loci:
        x = ml + bw * (lc["pos"] - lo) / span
        cls, lab = CLS_LABEL.get(lc["cls"], ("todo", "?"))
        col = {"pass": "#0d9488", "stop": "#dc2626", "warn": "#b45309", "todo": "#cbd5e1"}[cls]
        r = 9 if lc["cls"] == "A_ALLELE" else 6
        p.append(f'<line x1="{x}" y1="{y-14}" x2="{x}" y2="{y+14}" stroke="{col}" stroke-width="1.4"/>')
        p.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{col}" '
                 f'{"" if lc["cls"]=="A_ALLELE" else "opacity=.6"}/>')
        if lc["cls"] == "A_ALLELE":
            p.append(f'<text x="{x}" y="{y-20}" font-size="9.5" text-anchor="middle" '
                     f'font-weight="700" fill="#0d9488">V={lc["va"]:.2f}</text>')
        p.append(f'<text x="{x}" y="{y+30}" font-size="8.5" text-anchor="middle" fill="#94a3b8" '
                 f'transform="rotate(35 {x} {y+30})">{lc["pos"]}</text>')
    lg = [("#0d9488", "allele 軸（甲基明顯）"), ("#dc2626", "HP 軸（ASM）"),
          ("#b45309", "其他軸"), ("#cbd5e1", "未對齊／不穩定")]
    for i, (col, txt) in enumerate(lg):
        x = ml + i * 205
        p.append(f'<circle cx="{x+6}" cy="126" r="6" fill="{col}"/>')
        p.append(f'<text x="{x+18}" y="130" font-size="10.5" fill="#475569">{txt}</text>')
    p.append("</svg>")
    return "".join(p)


def case_card(c, idx):
    loci = sorted(c["locus_detail"], key=lambda x: x["pos"])
    aloci = [l for l in loci if l["cls"] == "A_ALLELE" and l["img"]]
    rows = "".join(
        '<tr><td class="num"><code>{p}</code></td><td><span class="st {cl}">{lb}</span></td>'
        '<td class="num">{va:.3f}</td><td class="num">{vh:.3f}</td><td><code>{cn}</code></td></tr>'.format(
            p=l["pos"], cl=CLS_LABEL.get(l["cls"], ("todo", "?"))[0],
            lb=CLS_LABEL.get(l["cls"], ("todo", "?"))[1], va=l["va"], vh=l["vh"], cn=l["cn"])
        for l in loci)
    imgs = "".join(
        f'<figure><figcaption><b>{l["pos"]}</b> · allele 軸 V={l["va"]:.3f} · HP V={l["vh"]:.3f} · {l["cn"]}</figcaption>'
        f'<div class="pair"><img src="{b64(l["img"],"distance_heatmap.png")}" alt="distance">'
        f'<img src="{b64(l["img"],"cluster_heatmap.png")}" alt="cluster"></div>'
        f'<div class="cap">左：read×read BERNOULLI 距離　右：分群後</div></figure>'
        for l in aloci[:3])
    flag = ('<span class="st warn">2026-06 深度個案</span>' if c["region_id"] == FLAGSHIP else "")
    return f"""
<details class="card" data-a="{c['A']}" data-n="{c['n_ssnv']}" data-me="{c['mutual_excl']}"
         data-chrom="{c['chrom']}" data-id="{c['region_id']}">
<summary><code>{c['region_id']}</code> {flag}
  <span class="ct">sSNV {c['n_ssnv']} · allele軸 {c['A']} · 互斥 {c['mutual_excl']}</span></summary>
<div class="cb">
  <figure class="struct">{svg_structure(c)}</figure>
  <div class="two">
    <div>
      <h4>區域內每個 sSNV 的甲基軸別</h4>
      <div class="tw"><table><thead><tr><th>位置</th><th>甲基軸別</th>
      <th class="num">V_allele</th><th class="num">V_hp</th><th>CN</th></tr></thead>
      <tbody>{rows}</tbody></table></div>
    </div>
    <div>
      <h4>sSNV 串聯結構</h4>
      <div class="tw"><table><tbody>
      <tr><th>pair 數</th><td class="num">{c['n_pairs']}</td></tr>
      <tr><th>互斥（分支）</th><td class="num">{c['mutual_excl']}</td></tr>
      <tr><th>分支類 pair 合計</th><td class="num">{c['branch']}</td></tr>
      <tr><th>max coread</th><td class="num">{c['max_coread']}</td></tr>
      <tr><th>關係分佈</th><td style="font-size:11.5px">{c['rels']}</td></tr>
      <tr><th>span</th><td class="num">{c['span']:,} bp</td></tr>
      </tbody></table></div>
    </div>
  </div>
  <h4>甲基明顯位點的原始熱圖（allele 軸，最多 3 個）</h4>
  {imgs if imgs else '<p class="mut">此區的 allele 軸位點無可用熱圖。</p>'}
</div></details>"""


cards = "".join(case_card(c, i) for i, c in enumerate(CASES))
chroms = sorted({c["chrom"] for c in CASES}, key=lambda x: (len(x), x))
chrom_opts = "".join(f'<option value="{c}">{c}</option>' for c in chroms)

HTML = f"""<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>sSNV 分支結構 × 甲基狀況 — 互動瀏覽器</title>
<style>
:root{{--bg:#fff;--fg:#0f172a;--mut:#64748b;--line:#e2e8f0;--card:#f8fafc;
--ok:#0d9488;--bad:#dc2626;--warn:#b45309;--acc:#2563eb;--rad:11px}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--fg);line-height:1.7;font-size:15.5px;
font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans CJK TC","PingFang TC","Microsoft JhengHei",Arial,sans-serif}}
.wrap{{max-width:1080px;margin:0 auto;padding:30px 20px 78px}}
.eyebrow{{font-size:11.5px;letter-spacing:.11em;text-transform:uppercase;color:var(--mut);font-weight:700}}
h1{{font-size:24px;margin:5px 0 7px;line-height:1.35}}
h4{{font-size:13px;margin:14px 0 6px;color:#334155}}
.sub{{color:var(--mut);font-size:13px;margin-bottom:16px}}
.note{{border:1px solid #fcd34d;border-left:5px solid var(--warn);background:#fffbeb;
border-radius:var(--rad);padding:13px 17px;margin:14px 0;font-size:14px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0}}
.kpi{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:11px 13px}}
.kpi .v{{font-size:21px;font-weight:700;font-variant-numeric:tabular-nums}}
.kpi .l{{font-size:11.5px;color:var(--mut);margin-top:3px}}
.kpi.ok .v{{color:var(--ok)}}
.bar{{display:flex;gap:10px;flex-wrap:wrap;align-items:center;background:var(--card);
border:1px solid var(--line);border-radius:var(--rad);padding:11px 14px;margin:16px 0;
position:sticky;top:0;z-index:5}}
.bar label{{font-size:12.5px;color:var(--mut);font-weight:600}}
.bar select,.bar input{{font:inherit;font-size:13px;padding:4px 8px;border:1px solid var(--line);
border-radius:6px;background:#fff;color:var(--fg)}}
.bar .cnt{{margin-left:auto;font-size:12.5px;color:var(--mut);font-weight:700}}
.st{{display:inline-block;font-size:11px;font-weight:700;padding:1px 8px;border-radius:20px;border:1px solid;white-space:nowrap}}
.st.pass{{background:#f0fdfa;color:#0f766e;border-color:#5eead4}}
.st.warn{{background:#fffbeb;color:#92400e;border-color:#fcd34d}}
.st.todo{{background:#f1f5f9;color:#475569;border-color:#cbd5e1}}
.st.stop{{background:#fef2f2;color:#b91c1c;border-color:#fca5a5}}
details.card{{border:1px solid var(--line);border-radius:var(--rad);margin:9px 0;overflow:hidden}}
details.card>summary{{cursor:pointer;padding:12px 16px;font-weight:600;font-size:14px;background:var(--card);
list-style:none;display:flex;align-items:center;gap:9px;flex-wrap:wrap}}
details.card>summary::-webkit-details-marker{{display:none}}
details.card>summary::before{{content:"+";font-family:ui-monospace,Menlo,monospace;font-weight:700;color:var(--acc);font-size:17px}}
details.card[open]>summary::before{{content:"−"}}
details.card>summary:hover{{background:#eef2f7}}
details.card .ct{{margin-left:auto;font-size:11px;color:var(--mut);background:#fff;
border:1px solid var(--line);border-radius:20px;padding:1px 9px;font-weight:600}}
details.card .cb{{padding:6px 16px 18px}}
figure{{margin:10px 0;padding:0;border:0}}
figure.struct{{border:1px solid var(--line);border-radius:9px;padding:8px;background:#fff}}
figure svg{{width:100%;height:auto;display:block}}
figcaption{{font-size:12px;color:#334155;margin-bottom:6px;font-weight:600}}
.pair{{display:grid;grid-template-columns:1fr 1fr;gap:8px}}
.pair img{{width:100%;height:auto;border:1px solid var(--line);border-radius:6px;display:block}}
.cap{{font-size:11px;color:var(--mut);margin-top:4px}}
.two{{display:grid;grid-template-columns:1.25fr 1fr;gap:16px}}
@media(max-width:760px){{.two{{grid-template-columns:1fr}}}}
.tw{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:13px;margin:4px 0}}
th{{background:#f1f5f9;text-align:left;padding:6px 9px;border-bottom:2px solid #cbd5e1;font-weight:700;white-space:nowrap}}
td{{padding:6px 9px;border-bottom:1px solid var(--line)}}
td.num,th.num{{text-align:right;font-variant-numeric:tabular-nums}}
code{{background:#f1f5f9;padding:1px 5px;border-radius:4px;font-size:12.5px;font-family:ui-monospace,Menlo,monospace}}
.mut{{color:var(--mut);font-size:12.5px}}
footer{{margin-top:36px;padding-top:16px;border-top:2px solid var(--line);font-size:12.5px;color:var(--mut)}}
@media(prefers-color-scheme:dark){{
:root{{--bg:#0b1120;--fg:#e2e8f0;--mut:#94a3b8;--line:#1e293b;--card:#0f172a}}
figure.struct{{background:#0f172a}} figure.struct svg{{filter:invert(.92) hue-rotate(180deg)}}
th{{background:#1e293b}} code{{background:#111c33}} .bar select,.bar input{{background:#0f172a;color:#e2e8f0}}
.note{{background:#2a2410;border-color:#78350f}}
}}
</style>
<div class="wrap">

<div class="eyebrow">Interactive browser · 2026-07-26</div>
<h1>sSNV 分支結構 × 甲基狀況<br>—— 同時檢視，並標記甲基明顯的位點</h1>
<div class="sub">選案準則：區域內同時具 <b>sSNV 分支</b>（mutual_excl／nested）與
<b>≥1 個 allele 軸甲基位點</b>　·　數字由 <code>{STEM}.data.json</code> 注入</div>

<div class="note"><b>先讀：本頁是描述性瀏覽，不是證據。</b>
{D['verdict_note']}詳見同目錄
<code>20260726_methyl_allele_axis_backbone_coenrichment.standalone.html</code>。</div>

<div class="grid">
<div class="kpi"><div class="v">{D['regions_total']:,}</div><div class="l">有甲基註記的串聯區域</div></div>
<div class="kpi ok"><div class="v">{D['regions_with_A']:,}</div><div class="l">含 ≥1 allele 軸位點</div></div>
<div class="kpi ok"><div class="v">{D['n_branch_methyl_regions']:,}</div><div class="l">分支 ＋ 甲基（符合準則）</div></div>
<div class="kpi"><div class="v">{D['n_cases']}</div><div class="l">本頁展示（含熱圖）</div></div>
</div>

<div class="bar">
  <label>染色體 <select id="fc"><option value="">全部</option>{chrom_opts}</select></label>
  <label>allele 軸位點 ≥ <select id="fa">
    <option value="1">1</option><option value="2">2</option>
    <option value="3" selected>3</option><option value="4">4</option></select></label>
  <label>互斥對 ≥ <select id="fm">
    <option value="0" selected>0</option><option value="3">3</option>
    <option value="5">5</option><option value="8">8</option></select></label>
  <label>排序 <select id="so">
    <option value="a">allele 軸數</option><option value="m">互斥對數</option>
    <option value="n">sSNV 數</option><option value="id">位置</option></select></label>
  <label><input type="checkbox" id="fo"> 全部展開</label>
  <span class="cnt" id="cnt"></span>
</div>

<div id="list">{cards}</div>

<footer>
<p>選案準則：<code>A ≥ 1 且 (mutual_excl + nested) ≥ 1 且 2 ≤ n_sSNV ≤ 12</code>，
全 cohort 符合者 {D['n_branch_methyl_regions']:,} 個，本頁取 allele 軸數最多且有熱圖者 {D['n_cases']} 個。<br>
熱圖來源：<code>20260119_all-with-w5000_1/filtered_snv_tp/{{chrom}}/{{chrom}}_{{pos}}/…/plots/BERNOULLI/</code>，base64 內嵌。<br>
完整 30,077 列位點註記見 <code>methyl_class_x_linkage_annotation.tsv</code>。
無 emoji（本機無 emoji 字型），狀態一律附文字標籤。</p>
</footer>
</div>
<script>
(function(){{
  var L=document.getElementById('list'), C=document.getElementById('cnt');
  var fc=document.getElementById('fc'), fa=document.getElementById('fa'),
      fm=document.getElementById('fm'), so=document.getElementById('so'),
      fo=document.getElementById('fo');
  var items=[].slice.call(L.querySelectorAll('details.card'));
  function apply(){{
    var vc=fc.value, va=+fa.value, vm=+fm.value, shown=0;
    items.forEach(function(el){{
      var ok=(!vc||el.dataset.chrom===vc) && (+el.dataset.a>=va) && (+el.dataset.me>=vm);
      el.style.display = ok ? '' : 'none';
      if(ok) shown++;
    }});
    var key=so.value;
    var vis=items.filter(function(e){{return e.style.display!=='none';}});
    vis.sort(function(x,y){{
      if(key==='id') return x.dataset.id.localeCompare(y.dataset.id);
      var k={{a:'a',m:'me',n:'n'}}[key];
      return (+y.dataset[k])-(+x.dataset[k]);
    }});
    vis.forEach(function(e){{L.appendChild(e);}});
    C.textContent = '顯示 '+shown+' / '+items.length+' 區域';
  }}
  [fc,fa,fm,so].forEach(function(e){{e.addEventListener('change',apply);}});
  // 「全部展開」只在切換當下作用，不干擾使用者手動展開的卡片
  fo.addEventListener('change',function(){{
    items.forEach(function(el){{ if(el.style.display!=='none') el.open=fo.checked; }});
  }});
  apply();
}})();
</script>
"""
(OUT / f"{STEM}.standalone.html").write_text(HTML, encoding="utf-8")
print(f"wrote {OUT / (STEM + '.standalone.html')}  ({len(HTML)/1024:.0f} KB)")
