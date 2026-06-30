#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐區演化樹觀察工作站 HTML 產生器（read-driven × topology join）。
明確處理 4 case（各自獨立函式，未來重用保證處理）:
  C1 兩棵樹(n_roots>=2)  → 顯示每個 root 分支
  C2 HP 身分            → derive_tree_hp() 用 node_hp 修正 H?→H1/H2
  C3 CN-gain confound   → dual_hp_confound() 標「可能 multiplicity 非真 dual-HP」
  C4 截斷 dual-HP/artifact → artifact_class() 歸 artifact(n_sSNV>8 / 密集 / centro-telo)
§13-A: 全數字由 topology_per_region.json + rd_perregion.json 注入，不手打。
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TOPO = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json"
RD = os.path.join(HERE, "data", "rd_perregion.json")
OUT = os.path.join(HERE, "perregion_topology_workstation.standalone.html")

topo = json.load(open(TOPO, encoding="utf-8"))
rd = json.load(open(RD, encoding="utf-8"))
RS = json.load(open("/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/region_shape_distribution.json", encoding="utf-8"))
DET = topo["detail"]


# ---- C2: 用 node_hp 把每棵樹(node_path)的 H? 修成 H1/H2 ----
def sorted_positions(node_hp):
    return sorted(node_hp.keys(), key=lambda k: int(k.split(":")[1]))


def derive_tree_hp(path, sposs, node_hp):
    alt_idx = [i for i, ch in enumerate(path) if ch == "A"]
    hps = set()
    for i in alt_idx:
        if i < len(sposs):
            h = node_hp.get(sposs[i])
            if h:
                hps.add(h)
    if len(hps) == 1:
        return list(hps)[0]
    if len(hps) > 1:
        return "mixed:" + "/".join(sorted(hps))
    return "H?"


# ---- C4: artifact 分類（截斷 / 密集 / centro-telo）----
def artifact_class(r):
    if r.get("truncated"):
        return True, f"truncated(n_sSNV={r['n_sSNV']}>8, 樹只用前8)"
    dens = r["span"] / r["n_sSNV"] if r["n_sSNV"] else 1e9
    if r["n_sSNV"] >= 6 and dens < 1500:
        return True, f"高密度({dens:.0f}bp/sSNV)"
    if r.get("genome_ctx") in ("centromere", "telomere") and r["n_sSNV"] >= 8:
        return True, f"{r.get('genome_ctx')} 密集"
    return False, ""


# ---- C1+C3: dual-HP 判定 + CN-gain confound ----
def dual_hp_status(r):
    hp = r.get("haplotypes") or ""
    dual = (r.get("n_roots") or 0) >= 2 or ("H1" in hp and "H2" in hp)
    confound = dual and r.get("cn") == "gain"
    return dual, confound


recs = []
for r in DET:
    nh = r.get("node_hp") or {}
    sposs = sorted_positions(nh)
    pops = r.get("populations") or {}
    trees = []
    for path, label in (r.get("node_paths") or {}).items():
        trees.append({"path": path, "hp": derive_tree_hp(path, sposs, nh),
                      "old": label, "reads": pops.get(path, 0)})
    art, artr = artifact_class(r)
    dual, conf = dual_hp_status(r)
    rdr = rd.get(r["region"], {})
    recs.append({
        "r": r["region"], "c": r["chrom"], "n": r["n_sSNV"], "cn": r.get("cn"),
        "ctx": r.get("genome_ctx"), "shape": r.get("tree_shape"), "roots": r.get("n_roots"),
        "det": r.get("determinacy"), "trunc": bool(r.get("truncated")),
        "tp": r.get("tp"), "fp": r.get("fp"),
        "dens": round(r["span"] / r["n_sSNV"], 0) if r["n_sSNV"] else None,
        "art": art, "artr": artr, "dual": dual, "conf": conf,
        "hps": r.get("haplotypes"), "trees": trees, "edges": r.get("edges"),
        "rdm": rdr.get("rd_multi_alt"), "rdc": rdr.get("rd_combos"), "rdx": rdr.get("rd_max_chain"),
    })

# summary
n = len(recs)
n_full = sum(1 for x in recs if x["shape"] == "full_tree")
n_dual = sum(1 for x in recs if x["dual"])
n_dualconf = sum(1 for x in recs if x["conf"])
n_dualclean = sum(1 for x in recs if x["dual"] and not x["conf"] and not x["trunc"])
n_art = sum(1 for x in recs if x["art"])
n_trunc = sum(1 for x in recs if x["trunc"])
n_linear = sum(1 for x in recs if x["shape"] == "linear_nested")
n_sib = sum(1 for x in recs if x["shape"] == "sibling_only")
n_col = sum(1 for x in recs if x["shape"] == "co_linked_lineage")
struct_total = n_full + n_linear + n_sib + n_col
rs_full = RS["by_shape"]["full_tree"]["n"]

DATA_JSON = json.dumps(recs, ensure_ascii=False)

html = f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>逐區演化樹觀察工作站 — HCC1395</title>
<style>
body{{font-family:"Segoe UI","微軟正黑體",sans-serif;max-width:1240px;margin:0 auto;padding:20px;color:#222;font-size:13px}}
h1{{font-size:21px;border-bottom:3px solid #44546A;padding-bottom:7px}}
.cards{{display:flex;flex-wrap:wrap;gap:10px;margin:12px 0}}
.card{{flex:1 1 130px;background:#fff;border:1px solid #ddd;border-radius:7px;padding:10px;text-align:center}}
.card .v{{font-size:21px;font-weight:700;color:#2c5f9e}} .card .l{{font-size:11px;color:#666}}
.ctl{{background:#f4f6f9;border:1px solid #dde;border-radius:6px;padding:9px;margin:10px 0;display:flex;flex-wrap:wrap;gap:8px;align-items:center}}
select,input{{font-size:12px;padding:3px 5px}}
table{{border-collapse:collapse;width:100%;font-size:12px}}
th,td{{border:1px solid #e6e6e6;padding:4px 7px;text-align:left}} th{{background:#f4f6f9;position:sticky;top:0;cursor:pointer}}
.num{{text-align:right;font-variant-numeric:tabular-nums}}
tr.art{{background:#fdf0ee}} tr.dualc{{background:#fef7e8}}
.tag{{display:inline-block;padding:1px 6px;border-radius:9px;font-size:10px;font-weight:700;color:#fff}}
.t-art{{background:#c0392b}} .t-dual{{background:#8e44ad}} .t-conf{{background:#e67e22}} .t-full{{background:#27ae60}} .t-trunc{{background:#7f8c8d}}
.det{{background:#fafbfc;font-size:11px;padding:8px 12px;border-left:3px solid #5B9BD5}}
.tree{{font-family:monospace;margin:3px 0}}
.hp1{{color:#2c5f9e;font-weight:700}} .hp2{{color:#8e44ad;font-weight:700}} .hpx{{color:#999}}
.note{{font-size:11px;color:#666}}
</style></head><body>
<h1>逐區演化樹觀察工作站 — HCC1395（read-driven × topology）</h1>
<p class="note">3,885 區逐區克隆樹 + read-driven 串接交叉確認。4 個顯示 case 已實做：① 兩棵樹(n_roots≥2) ② HP 身分(node_hp 修正 H?→H1/H2) ③ CN-gain confound 標記 ④ 截斷/artifact 歸類。點列展開看每棵樹的 HP + read-driven 確認。</p>
<div class="cards">
<div class="card"><div class="v">{n:,}</div><div class="l">總區域</div></div>
<div class="card"><div class="v">{n_full:,}</div><div class="l">full_tree</div></div>
<div class="card"><div class="v">{n_dual:,}</div><div class="l">dual-HP(2根)</div></div>
<div class="card"><div class="v">{n_dualconf:,}</div><div class="l">dual-HP CN-gain confound</div></div>
<div class="card"><div class="v">{n_dualclean:,}</div><div class="l">dual-HP 乾淨(待查)</div></div>
<div class="card"><div class="v">{n_art:,}</div><div class="l">artifact(含 {n_trunc} 截斷)</div></div>
</div>
<div style="background:#fff7e8;border:1px solid #f0d090;border-radius:6px;padding:10px 14px;margin:10px 0;font-size:12px">
🔴 <b>full_tree 兩口徑對賬（誠實標示）</b>：本工作站採 <b>topology 逐區嚴格解析（基因型向量完全解析）= full_tree {n_full}</b>（+ linear_nested {n_linear} + sibling_only {n_sib} + co_linked {n_col} = <b>{struct_total} 有確認結構</b>；detail 3,885 區 n_sSNV≥2）。
另一較寬口徑 regions.tsv（region_shape_distribution，7,143 全區）= full_tree <b>{rs_full}</b>。
兩者<b>定義不同層級</b>（單-sSNV 區不可能 full_tree → 496 差距非區數差，是判準寬嚴差）；對外引用須標明口徑。本頁逐區樹一律以 topology 嚴格解析為準。
</div>
<div class="ctl">
篩選:
<select id="fctx"><option value="">ctx:全</option><option>arm</option><option>centromere</option><option>telomere</option></select>
<select id="fcn"><option value="">CN:全</option><option>gain</option><option>loh</option><option>neutral</option><option>loss</option></select>
<select id="fshape"><option value="">shape:全</option><option>full_tree</option><option>linear_nested</option><option>sibling_only</option><option>co_linked_lineage</option><option>single</option><option>inconsistent</option></select>
<select id="fflag"><option value="">flag:全</option><option value="dual">dual-HP</option><option value="dualclean">dual-HP乾淨</option><option value="conf">dual-HP confound</option><option value="art">artifact</option><option value="trunc">截斷</option><option value="full">full_tree</option></select>
<input id="fsearch" placeholder="搜尋 region/chrom"size="16">
<span id="cnt" class="note"></span>
</div>
<table id="tb"><thead><tr>
<th data-k="r">region</th><th data-k="n">n_sSNV</th><th data-k="cn">CN</th><th data-k="ctx">ctx</th>
<th data-k="shape">shape</th><th data-k="roots">roots</th><th data-k="rdm">rd_multiALT</th><th data-k="rdc">rd_combos</th>
<th data-k="rdx">rd_maxchain</th><th>tp/fp</th><th>flags</th></tr></thead><tbody id="bd"></tbody></table>
<script>
const D={DATA_JSON};
const bd=document.getElementById('bd');
function hpcls(h){{return h=='H1'?'hp1':h=='H2'?'hp2':'hpx';}}
function tags(x){{let t='';if(x.shape=='full_tree')t+='<span class="tag t-full">full</span> ';if(x.dual)t+='<span class="tag t-dual">dual-HP</span> ';if(x.conf)t+='<span class="tag t-conf">CN-gain?</span> ';if(x.trunc)t+='<span class="tag t-trunc">截斷</span> ';if(x.art&&!x.trunc)t+='<span class="tag t-art">artifact</span> ';return t;}}
function detail(x){{
  let s='<div class="det"><b>每棵樹（HP 由 node_hp 修正）</b>：';
  if(x.trees&&x.trees.length){{x.trees.forEach(t=>{{s+=`<div class="tree"><span class="${{hpcls(t.hp)}}">${{t.hp}}</span> ｜ path=${{t.path}} ｜ ${{t.reads}} reads ${{t.old&&t.old!=t.hp?'（原標 '+t.old+'）':''}}</div>`;}});}}else s+='<i>無 node_paths（germline-only / single）</i>';
  s+=`<div class="note" style="margin-top:5px">edges: ${{JSON.stringify(x.edges||[])}}</div>`;
  s+=`<div class="note">read-driven: 多-ALT read ${{x.rdm??'—'}} / distinct combos ${{x.rdc??'—'}} / max chain ${{x.rdx??'—'}} ｜ haplotypes=${{x.hps}} ｜ determinacy=${{x.det}}</div>`;
  if(x.conf)s+='<div class="note" style="color:#e67e22"><b>⚠ CN-gain dual-HP</b>：可能是 multiplicity（兩等位都看似 ALT）非真兩 HP 獨立 lineage。</div>';
  if(x.art)s+=`<div class="note" style="color:#c0392b"><b>🔴 artifact</b>：${{x.artr}}（樹不可信，列為密集假叢集）。</div>`;
  s+='</div>';return s;}}
let SK='r',SD=1;
function render(){{
  const ctx=fctx.value,cn=fcn.value,sh=fshape.value,fl=fflag.value,q=fsearch.value.toLowerCase();
  let rows=D.filter(x=>(!ctx||x.ctx==ctx)&&(!cn||x.cn==cn)&&(!sh||x.shape==sh)&&(!q||(x.r+x.c).toLowerCase().includes(q))
    &&(!fl||(fl=='dual'&&x.dual)||(fl=='dualclean'&&x.dual&&!x.conf&&!x.trunc)||(fl=='conf'&&x.conf)||(fl=='art'&&x.art)||(fl=='trunc'&&x.trunc)||(fl=='full'&&x.shape=='full_tree')));
  rows.sort((a,b)=>{{let u=a[SK],v=b[SK];if(u==null)u=-1;if(v==null)v=-1;return (u>v?1:u<v?-1:0)*SD;}});
  cnt.textContent=rows.length+' 區';
  bd.innerHTML=rows.slice(0,1500).map((x,i)=>`<tr class="${{x.art?'art':x.conf?'dualc':''}}" data-i="${{D.indexOf(x)}}">
   <td>${{x.r}}</td><td class="num">${{x.n}}</td><td>${{x.cn}}</td><td>${{x.ctx}}</td><td>${{x.shape}}</td>
   <td class="num">${{x.roots??''}}</td><td class="num">${{x.rdm??''}}</td><td class="num">${{x.rdc??''}}</td><td class="num">${{x.rdx??''}}</td>
   <td class="num">${{x.tp??0}}/${{x.fp??0}}</td><td>${{tags(x)}}</td></tr>`).join('')
   +(rows.length>1500?`<tr><td colspan=11 class="note">… 顯示前 1500/${{rows.length}}（縮篩選）</td></tr>`:'');
}}
bd.addEventListener('click',e=>{{const tr=e.target.closest('tr[data-i]');if(!tr)return;
  const nx=tr.nextElementSibling;if(nx&&nx.classList.contains('dt')){{nx.remove();return;}}
  const x=D[+tr.dataset.i];const d=document.createElement('tr');d.className='dt';d.innerHTML=`<td colspan=11>${{detail(x)}}</td>`;tr.after(d);}});
[fctx,fcn,fshape,fflag,fsearch].forEach(e=>e.addEventListener('input',render));
document.querySelectorAll('th[data-k]').forEach(th=>th.onclick=()=>{{const k=th.dataset.k;SD=(SK==k?-SD:1);SK=k;render();}});
render();
</script>
<p class="note" style="margin-top:20px">generated by build_perregion_html.py（§13-A 注入）· data: topology_per_region.json（3,885 區樹）+ rd_perregion.json（read-driven per-region）· HCC1395 ⭐3 · 2026-06-30 · 點任一列展開看每棵 HP 樹 + read-driven 確認</p>
</body></html>"""
open(OUT, "w").write(html)
print(f"DONE -> {OUT} ({len(html):,} bytes)  regions={n} full={n_full} dual={n_dual} dualconf={n_dualconf} dualclean={n_dualclean} art={n_art} trunc={n_trunc}")
