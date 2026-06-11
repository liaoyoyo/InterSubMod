#!/usr/bin/env python3
"""
75 - Interactive credible-locus EXPLORER HTML: per-locus methylation heatmaps +
full data, with filter (sample/cancer/tier/ARI/|Δβ|/gene-search) + sort + click-to-zoom.
Vanilla JS (no CDN) => fully standalone/offline. §13 layer-A (data from JSONs/CSVs).

Sources:
  - 5x <sample>_locus_figs.json (script 67): per-locus base64 heatmap + ARI/Δβ/...
  - ism_existence_scan/<sample>_tp/significance_summary.csv: enrich with ISM native
    stats (Significant, CramersV_HPFine, HPMergedDelta, ClusterPermanovaP, NumReads)
  - ism_cis_scan/HCC1395_tp_cis/significance_summary.csv: HP_Residual (HCC1395 cis)

Output: docs/experiments/in_progress/2026/06/20260604_credible_locus_interactive_explorer_01.standalone.html
"""
import os, csv, json, html

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
CISD = "/big7_disk/liaoyoyo2001/ism_cis_scan"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
       "20260604_credible_locus_interactive_explorer_01.standalone.html")

SAMP_WITH_FIGS = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009"]
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}


def fnum(v):
    try:
        f = float(v)
        return None if f != f else round(f, 4)
    except (ValueError, TypeError):
        return None


def load_ism(path, want):
    """key (chrom,pos) -> {col:val} for wanted cols."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for r in csv.DictReader(f):
            ch = r.get("Chr") or r.get("Chrom") or r.get("chrom") or r.get("CHROM")
            ps = r.get("Pos") or r.get("pos") or r.get("POS")
            if ch is None or ps is None:
                continue
            out[(ch, str(ps))] = {w: r.get(w) for w in want}
    return out


# build locus records
records = []
for s in SAMP_WITH_FIGS:
    fp = f"{CS}/{s}_locus_figs.json"
    if not os.path.exists(fp):
        continue
    figs = json.load(open(fp))["figs"]
    ism = load_ism(f"{EX}/{s}_tp/significance_summary.csv",
                   ["Significant", "CramersV_HPFine", "HPMergedDelta", "ClusterPermanovaP", "NumReads"])
    cis = load_ism(f"{CISD}/HCC1395_tp_cis/significance_summary.csv",
                   ["HP_Residual_Delta", "HP_Residual_P"]) if s == "HCC1395" else {}
    for x in figs:
        ch, ps = x["locus"].split(":")
        im = ism.get((ch, ps), {})
        ci = cis.get((ch, ps), {})
        records.append(dict(
            sample=s, cancer=CANCER[s], gene=x["gene"], locus=x["locus"],
            axis=x["axis"].replace("HP", "").replace("_vs_", "v"),
            delta=x["delta"], ari=x["ari"], placebo=x.get("placebo_ari"),
            ncpg=x["n_paired_cpg"], ctx=x.get("cpg_context"),
            tier="strict" if x.get("tier") == "strict_ge100" else "relax",
            ism_sig=str(im.get("Significant", "")).lower() == "true",
            cramers=fnum(im.get("CramersV_HPFine")),
            hpmerged=fnum(im.get("HPMergedDelta")),
            cis_resid=fnum(ci.get("HP_Residual_Delta")),
            fig=x["fig"],
        ))

# sort by |delta| desc default
records.sort(key=lambda r: -abs(r["delta"]))
DATA = json.dumps(records)
n = len(records)
genes = sorted(set(r["gene"] for r in records))

HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Credible ASM locus 互動 explorer</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC",sans-serif;color:var(--ink);background:var(--bg);line-height:1.5}}
.wrap{{max-width:1280px;margin:0 auto;padding:18px 14px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:12px;padding:18px 22px;margin-bottom:12px}}
header.top h1{{margin:0 0 4px;font-size:1.24rem}}.meta{{font-size:.8rem;opacity:.88}}
.bar{{position:sticky;top:0;z-index:5;background:#fff;border:1px solid var(--line);border-radius:10px;padding:10px 12px;margin-bottom:10px;display:flex;gap:10px;flex-wrap:wrap;align-items:center;font-size:.82rem;box-shadow:0 2px 8px rgba(0,0,0,.04)}}
.bar label{{font-size:.74rem;color:var(--mut);display:flex;flex-direction:column;gap:2px}}
.bar select,.bar input{{font-size:.8rem;padding:3px 6px;border:1px solid var(--line);border-radius:6px}}
.bar .rng{{width:120px}}
.count{{margin-left:auto;font-weight:700;color:var(--accent)}}
table{{width:100%;border-collapse:collapse;font-size:.78rem;background:#fff}}
th,td{{padding:5px 7px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}}
th{{background:#f1f5f9;font-size:.72rem;cursor:pointer;position:sticky;top:62px;user-select:none}}
th:hover{{background:#e2e8f0}} th.sorted::after{{content:" ▾"}} th.sortedup::after{{content:" ▴"}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}.mono{{font-family:ui-monospace,monospace;font-size:.74rem;color:#475569}}
.chip{{color:#fff;font-weight:700;font-size:.68rem;padding:1px 6px;border-radius:5px}}
.thumb{{height:46px;border:1px solid var(--line);border-radius:4px;cursor:zoom-in;display:block}}
.yes{{color:#15803d;font-weight:700}}.no{{color:#9ca3af}}
.tablewrap{{overflow-x:auto;border:1px solid var(--line);border-radius:10px;max-height:74vh;overflow-y:auto}}
#modal{{display:none;position:fixed;inset:0;background:rgba(0,0,0,.82);z-index:50;align-items:center;justify-content:center;cursor:zoom-out}}
#modal img{{max-width:94vw;max-height:88vh;border:3px solid #fff;border-radius:6px}}
#modal .cap{{position:absolute;top:12px;left:0;right:0;text-align:center;color:#fff;font-size:.9rem}}
.legend{{font-size:.74rem;color:var(--mut);margin:6px 2px}}
</style></head><body><div class="wrap">

<header class="top"><h1>Credible ASM locus 互動 explorer — 圖 + 數據 + 篩選 + 縮放</h1>
<div class="meta">{n} 個 credible（pass tierA）位點，每個含 read×CpG 甲基熱圖 + discovery 數據 + ISM 完整統計 enrich · 2026-06-04 · vanilla JS standalone（離線可用）· §13 layer-A</div></header>

<div class="bar">
  <label>樣本<select id="f_sample"><option value="">全部</option>{"".join(f'<option>{esc}</option>' for esc in SAMP_WITH_FIGS)}</select></label>
  <label>癌種<select id="f_cancer"><option value="">全部</option><option>breast</option><option>lung</option></select></label>
  <label>tier<select id="f_tier"><option value="">全部</option><option>strict</option><option>relax</option></select></label>
  <label>ARI ≥ <span id="v_ari">0</span><input class="rng" id="f_ari" type="range" min="0" max="1" step="0.05" value="0"></label>
  <label>|Δβ| ≥ <span id="v_delta">0</span><input class="rng" id="f_delta" type="range" min="0" max="0.5" step="0.02" value="0"></label>
  <label>ISM<select id="f_ism"><option value="">全部</option><option value="1">Significant</option><option value="0">非</option></select></label>
  <label>基因搜尋<input id="f_gene" type="text" placeholder="gene..." size="10"></label>
  <span class="count" id="count"></span>
</div>
<div class="legend">點熱圖縮放放大。熱圖：上=germline 主單倍型、下=somatic HPx-1（紅=甲基/藍=非甲基/橘=變異）。點欄位標題排序。</div>

<div class="tablewrap"><table id="tbl">
<thead><tr>
  <th data-k="fig">圖</th><th data-k="sample">樣本</th><th data-k="cancer">癌種</th>
  <th data-k="gene">gene</th><th data-k="locus">locus</th><th data-k="axis">axis</th>
  <th data-k="delta" class="num">Δβ</th><th data-k="ari" class="num">ARI</th><th data-k="placebo" class="num">placebo</th>
  <th data-k="ncpg" class="num">nCpG</th><th data-k="ctx">CpG</th><th data-k="tier">tier</th>
  <th data-k="ism_sig">ISM-Sig</th><th data-k="cramers" class="num">CramersV</th>
  <th data-k="hpmerged" class="num">ISM Δ</th><th data-k="cis_resid" class="num">cis resid</th>
</tr></thead><tbody id="tbody"></tbody></table></div>

<div id="modal"><div class="cap" id="mcap"></div><img id="mimg"></div>

<script>
const DATA = {DATA};
const CC = {{breast:"#be123c",lung:"#0e7490",melanoma:"#7c3aed"}};
let sortK="delta", sortDir=-1;
const $=id=>document.getElementById(id);
function fmt(v,d=2){{ if(v===null||v===undefined||v==="")return "—"; if(typeof v==="number")return (v>=0?"+":"")+v.toFixed(d); return v; }}
function filtered(){{
  const fs=$("f_sample").value, fc=$("f_cancer").value, ft=$("f_tier").value,
        fa=parseFloat($("f_ari").value), fd=parseFloat($("f_delta").value),
        fi=$("f_ism").value, fg=$("f_gene").value.toLowerCase();
  return DATA.filter(r=>
    (!fs||r.sample===fs)&&(!fc||r.cancer===fc)&&(!ft||r.tier===ft)&&
    (r.ari>=fa)&&(Math.abs(r.delta)>=fd)&&
    (fi===""||(fi==="1")===r.ism_sig)&&
    (!fg||(r.gene||"").toLowerCase().includes(fg))
  );
}}
function render(){{
  let rows=filtered();
  rows.sort((a,b)=>{{
    let x=a[sortK],y=b[sortK];
    if(sortK==="delta"){{x=Math.abs(x);y=Math.abs(y);}}
    if(typeof x==="boolean"){{x=x?1:0;y=y?1:0;}}
    if(x===null)x=-1e9; if(y===null)y=-1e9;
    if(typeof x==="string")return sortDir*x.localeCompare(y);
    return sortDir*(x-y);
  }});
  const tb=$("tbody"); tb.innerHTML="";
  const fr=document.createDocumentFragment();
  for(const r of rows){{
    const tr=document.createElement("tr");
    tr.innerHTML=
      `<td><img class="thumb" src="${{r.fig}}" data-f="${{r.fig}}" data-c="${{r.gene}} ${{r.locus}} (${{r.sample}})"></td>`+
      `<td><span class="chip" style="background:${{CC[r.cancer]}}">${{r.sample}}</span></td>`+
      `<td>${{r.cancer}}</td><td>${{r.gene||""}}</td><td class="mono">${{r.locus}}</td><td>${{r.axis}}</td>`+
      `<td class="num">${{fmt(r.delta)}}</td><td class="num">${{fmt(r.ari)}}</td><td class="num">${{fmt(r.placebo)}}</td>`+
      `<td class="num">${{r.ncpg}}</td><td>${{r.ctx||""}}</td><td>${{r.tier}}</td>`+
      `<td class="${{r.ism_sig?'yes':'no'}}">${{r.ism_sig?'✓':'·'}}</td>`+
      `<td class="num">${{fmt(r.cramers,3)}}</td><td class="num">${{fmt(r.hpmerged)}}</td>`+
      `<td class="num">${{r.cis_resid===null?'—':fmt(r.cis_resid)}}</td>`;
    fr.appendChild(tr);
  }}
  tb.appendChild(fr);
  $("count").textContent=rows.length+" / "+DATA.length+" 位點";
  document.querySelectorAll("#tbl th").forEach(th=>{{th.classList.remove("sorted","sortedup");
    if(th.dataset.k===sortK)th.classList.add(sortDir<0?"sorted":"sortedup");}});
}}
document.querySelectorAll("#tbl th").forEach(th=>th.onclick=()=>{{
  const k=th.dataset.k; if(k==="fig")return;
  if(sortK===k)sortDir*=-1; else {{sortK=k;sortDir=(k==="gene"||k==="locus"||k==="sample"||k==="cancer"||k==="axis"||k==="ctx"||k==="tier")?1:-1;}}
  render();
}});
["f_sample","f_cancer","f_tier","f_ism","f_gene"].forEach(id=>$(id).oninput=render);
$("f_ari").oninput=()=>{{$("v_ari").textContent=$("f_ari").value;render();}};
$("f_delta").oninput=()=>{{$("v_delta").textContent=$("f_delta").value;render();}};
document.addEventListener("click",e=>{{
  if(e.target.classList.contains("thumb")){{$("mimg").src=e.target.dataset.f;$("mcap").textContent=e.target.dataset.c;$("modal").style.display="flex";}}
}});
$("modal").onclick=()=>$("modal").style.display="none";
render();
</script>
</div></body></html>'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(HTML_DOC)
print(f"[75] wrote {OUT} ({len(HTML_DOC)//1024} KB, {n} loci, {len(genes)} genes)")
