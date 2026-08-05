#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""逐區演化樹觀察工作站 HTML 產生器 v2（畫真 clone 樹 SVG + read-driven join）。
4 case（各獨立處理，未來重用保證）:
  C1 兩棵樹(n_roots>=2)   → 樹 SVG 顯示多 root 分支
  C2 HP 身分             → tree_hp() 以 node_paths 標籤為主(pipeline 已算)、H? 且 node_hp 完整才補；節點依 HP 上色
  C3 CN-gain confound    → dual_hp_status() 標「可能 multiplicity」
  C4 截斷/artifact        → artifact_class() 不畫假樹、改 artifact 註記 + read-driven 原始串接
  + HP 未確認(H?/HP3/unphased) 誠實標示與計數（單樣本相位真實限制）
§13-A: 全數字由 topology_per_region.json + rd_perregion.json 注入。
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
TOPO = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/topology_per_region.json"
RD = os.path.join(HERE, "data", "rd_perregion.json")
RSP = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/region_shape_distribution.json"
OUT = os.path.join(HERE, "perregion_topology_workstation.standalone.html")

topo = json.load(open(TOPO, encoding="utf-8"))
rd = json.load(open(RD, encoding="utf-8"))
RS = json.load(open(RSP, encoding="utf-8"))
DET = topo["detail"]


# ---- C2: HP 判定（標籤為主，H?+node_hp 完整才從 node_hp 補；修正先前 index bug）----
def tree_hp(label, path, node_hp, n_sSNV):
    for h in ("H1", "H2"):
        if label.startswith(h):
            return h
    if label.startswith("H3"):
        return "H3?"
    # label = H? → dual-HP fallback：node_hp 完整(==n_sSNV)時才用 index 推
    if len(node_hp) == n_sSNV and n_sSNV > 0:
        sp = sorted(node_hp.keys(), key=lambda k: int(k.split(":")[1]))
        hs = {node_hp.get(sp[i]) for i, c in enumerate(path) if c == "A" and i < len(sp)}
        hs = {h for h in hs if h}
        if len(hs) == 1:
            return list(hs)[0]
        if len(hs) > 1:
            return "mixed"
    return "H?"


# ---- C4: artifact 分類 ----
def artifact_class(r):
    if r.get("truncated"):
        return True, f"截斷 n_sSNV={r['n_sSNV']}>8（樹只用前8）"
    dens = r["span"] / r["n_sSNV"] if r["n_sSNV"] else 1e9
    if r["n_sSNV"] >= 6 and dens < 1500:
        return True, f"高密度 {dens:.0f}bp/sSNV"
    if r.get("genome_ctx") in ("centromere", "telomere") and r["n_sSNV"] >= 8:
        return True, f"{r.get('genome_ctx')} 密集"
    return False, ""


# ---- C1+C3: dual-HP + CN-gain confound ----
def dual_hp_status(r):
    hp = r.get("haplotypes") or ""
    dual = (r.get("n_roots") or 0) >= 2 or ("H1" in hp and "H2" in hp)
    return dual, (dual and r.get("cn") == "gain")


recs = []
for r in DET:
    nh = r.get("node_hp") or {}
    pops = r.get("populations") or {}
    trees = []
    hpbad = 0
    for path, label in (r.get("node_paths") or {}).items():
        h = tree_hp(label, path, nh, r["n_sSNV"])
        if h in ("H?", "mixed"):
            hpbad += 1
        trees.append({"path": path, "hp": h, "reads": pops.get(path, 0)})
    art, artr = artifact_class(r)
    dual, conf = dual_hp_status(r)
    rdr = rd.get(r["region"], {})
    recs.append({
        "r": r["region"], "c": r["chrom"], "n": r["n_sSNV"], "cn": r.get("cn"),
        "ctx": r.get("genome_ctx"), "shape": r.get("tree_shape"), "roots": r.get("n_roots") or 0,
        "det": r.get("determinacy"), "trunc": bool(r.get("truncated")),
        "tp": r.get("tp") or 0, "fp": r.get("fp") or 0,
        "art": art, "artr": artr, "dual": dual, "conf": conf, "hpbad": hpbad,
        "hps": r.get("haplotypes"), "germR": r.get("germline_reads") or 0,
        "trees": trees, "edges": r.get("edges") or [],
        "rdm": rdr.get("rd_multi_alt"), "rdc": rdr.get("rd_combos"), "rdx": rdr.get("rd_max_chain"),
    })

n = len(recs)
n_full = sum(1 for x in recs if x["shape"] == "full_tree")
n_lin = sum(1 for x in recs if x["shape"] == "linear_nested")
n_sib = sum(1 for x in recs if x["shape"] == "sibling_only")
n_col = sum(1 for x in recs if x["shape"] == "co_linked_lineage")
struct = n_full + n_lin + n_sib + n_col
n_dual = sum(1 for x in recs if x["dual"])
n_conf = sum(1 for x in recs if x["conf"])
n_art = sum(1 for x in recs if x["art"])
n_trunc = sum(1 for x in recs if x["trunc"])
n_hpbad = sum(1 for x in recs if x["hpbad"] > 0)
rs_full = RS["by_shape"]["full_tree"]["n"]
DATA_JSON = json.dumps(recs, ensure_ascii=False)

HEAD = f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>逐區演化樹工作站 — HCC1395</title>
<style>
body{{font-family:"Segoe UI","微軟正黑體",sans-serif;max-width:1240px;margin:0 auto;padding:18px;color:#222;font-size:13px}}
h1{{font-size:20px;border-bottom:3px solid #44546A;padding-bottom:6px}}
.cards{{display:flex;flex-wrap:wrap;gap:9px;margin:10px 0}}
.card{{flex:1 1 120px;background:#fff;border:1px solid #ddd;border-radius:7px;padding:9px;text-align:center}}
.card .v{{font-size:20px;font-weight:700;color:#2c5f9e}} .card .l{{font-size:10.5px;color:#666}}
.recon{{background:#fff7e8;border:1px solid #f0d090;border-radius:6px;padding:9px 13px;margin:9px 0;font-size:12px}}
.ctl{{background:#f4f6f9;border:1px solid #dde;border-radius:6px;padding:8px;margin:9px 0;display:flex;flex-wrap:wrap;gap:7px;align-items:center}}
select,input{{font-size:12px;padding:3px 5px}}
table{{border-collapse:collapse;width:100%;font-size:12px}}
th,td{{border:1px solid #e6e6e6;padding:4px 7px;text-align:left}} th{{background:#f4f6f9;cursor:pointer;position:sticky;top:0}}
.num{{text-align:right;font-variant-numeric:tabular-nums}}
tr.art{{background:#fdf0ee}} tr.dualc{{background:#fef7e8}} tr[data-i]{{cursor:pointer}}
.tag{{display:inline-block;padding:1px 6px;border-radius:9px;font-size:10px;font-weight:700;color:#fff;margin-right:2px}}
.t-art{{background:#c0392b}}.t-dual{{background:#8e44ad}}.t-conf{{background:#e67e22}}.t-full{{background:#27ae60}}.t-trunc{{background:#7f8c8d}}.t-hp{{background:#d39e00}}
.det{{background:#fafbfc;padding:8px 12px;border-left:3px solid #5B9BD5}}
.note{{font-size:11px;color:#666;margin:3px 0}}
</style></head><body>
<h1>逐區演化樹觀察工作站 — HCC1395（clone 樹 SVG × read-driven）</h1>
<p class="note">3,885 區逐區克隆樹（點列展開看樹 SVG）。節點依 HP 上色：<b style="color:#2c5f9e">H1</b>／<b style="color:#8e44ad">H2</b>／<b style="color:#2f9e44">H3?</b>／<span style="color:#999">HP未確認</span>。4 case：①兩棵樹(n_roots≥2) ②HP標籤上色 ③CN-gain confound ④截斷→artifact 不畫假樹。</p>
<div class="cards">
<div class="card"><div class="v">{n:,}</div><div class="l">總區域</div></div>
<div class="card"><div class="v">{n_full:,}</div><div class="l">full_tree(嚴格)</div></div>
<div class="card"><div class="v">{struct:,}</div><div class="l">有確認結構</div></div>
<div class="card"><div class="v">{n_dual:,}</div><div class="l">dual-HP(2根)</div></div>
<div class="card"><div class="v">{n_conf:,}</div><div class="l">dual-HP CN-gain</div></div>
<div class="card"><div class="v">{n_hpbad:,}</div><div class="l">含 HP未確認節點</div></div>
<div class="card"><div class="v">{n_art:,}</div><div class="l">artifact(截斷{n_trunc})</div></div>
</div>
<div class="recon">🔴 <b>full_tree 兩口徑</b>：topology 嚴格解析 = <b>{n_full}</b>（+linear {n_lin}+sibling {n_sib}+co_linked {n_col} = <b>{struct} 有確認結構</b>，3,885 區 n_sSNV≥2）；regions.tsv 寬鬆 = <b>{rs_full}</b>。本頁用嚴格口徑。｜
🔴 <b>HP 未確認</b>：{n_hpbad} 區含 H? 節點 = somatic reads 無 confident germline-HP（HP3/unphased，單樣本相位真實限制，非 bug）。</div>
<div class="ctl">篩選:
<select id="fcn"><option value="">CN:全</option><option>gain</option><option>loh</option><option>neutral</option><option>loss</option></select>
<select id="fctx"><option value="">ctx:全</option><option>arm</option><option>centromere</option><option>telomere</option></select>
<select id="fshape"><option value="">shape:全</option><option>full_tree</option><option>linear_nested</option><option>sibling_only</option><option>co_linked_lineage</option><option>single</option><option>no_confirmed_structure</option><option>inconsistent</option></select>
<select id="fflag"><option value="">flag:全</option><option value="full">full_tree</option><option value="dual">dual-HP</option><option value="conf">dual-HP CN-gain</option><option value="hp">含HP未確認</option><option value="art">artifact</option><option value="trunc">截斷</option></select>
<input id="fsearch" placeholder="搜尋 region/chrom" size="15"><span id="cnt" class="note"></span></div>
<table><thead><tr><th data-k="r">region</th><th data-k="n">n_sSNV</th><th data-k="cn">CN</th><th data-k="ctx">ctx</th><th data-k="shape">shape</th><th data-k="roots">roots</th><th data-k="rdx">rd_maxchain</th><th data-k="hpbad">HP未確認</th><th>tp/fp</th><th>flags</th></tr></thead><tbody id="bd"></tbody></table>
"""

JS = r'''<script>
var D=__DATA__;
function col(h){return h=='H1'?'#dbe9f7':h=='H2'?'#ece0f5':h=='H3?'?'#e3f4e1':(h=='mixed'?'#fde8d0':'#eeeeee');}
function stk(h){return h=='H1'?'#2c5f9e':h=='H2'?'#8e44ad':h=='H3?'?'#2f9e44':(h=='mixed'?'#e67e22':'#999999');}
function nz(v){return v==null?'—':v;}
function noteArt(msg,x){return '<div class="note" style="color:#c0392b">🔴 '+msg+'<br>read-driven 原始串接：max chain '+nz(x.rdx)+' / combos '+nz(x.rdc)+' / 多-ALT reads '+nz(x.rdm)+'（artifact 不建假樹）</div>';}
function treeSVG(x){
  if(x.trunc) return noteArt('截斷 n_sSNV='+x.n+'>8'+(x.cn=='gain'?'（CN-gain 密集）':''),x);
  if(x.art) return noteArt(x.artr,x);
  var edges=x.edges||[];
  if(!edges.length) return '<div class="note">germline-only / 單 sSNV：無樹結構</div>';
  var hpOf={},rdOf={};(x.trees||[]).forEach(function(t){hpOf[t.path]=t.hp;rdOf[t.path]=t.reads;});
  var ch={},nodes={};edges.forEach(function(e){(ch[e[0]]=ch[e[0]]||[]).push(e[1]);nodes[e[1]]=1;if(e[0]!='ROOT')nodes[e[0]]=1;});
  var depth={ROOT:0};(function dfs(nn,d){(ch[nn]||[]).forEach(function(c){if(depth[c]==null){depth[c]=d+1;dfs(c,d+1);}});})('ROOT',0);
  Object.keys(nodes).forEach(function(nn){if(depth[nn]==null)depth[nn]=1;});
  var byD={};Object.keys(nodes).forEach(function(nn){var d=depth[nn];(byD[d]=byD[d]||[]).push(nn);});
  var ds=Object.keys(depth).map(function(k){return depth[k];});var maxd=Math.max.apply(null,[1].concat(ds));
  var ws=Object.keys(byD).map(function(d){return byD[d].length;});var maxw=Math.max.apply(null,[1].concat(ws));
  var W=Math.max(380,150*maxw),VG=74,H=(maxd+1)*VG+14,gx=W/2;
  var pos={};Object.keys(byD).forEach(function(d){var a=byD[d];a.forEach(function(nn,i){pos[nn]=W*(i+1)/(a.length+1);});});
  function Y(d){return 26+d*VG;}
  var s='<svg viewBox="0 0 '+W+' '+H+'" width="100%" height="'+H+'" style="max-width:740px;font-family:monospace">';
  s+='<rect x="'+(gx-48)+'" y="'+(Y(0)-15)+'" width="96" height="30" rx="6" fill="#fff" stroke="#495057"/><text x="'+gx+'" y="'+(Y(0)+4)+'" text-anchor="middle" font-size="10">germline '+(x.germR||0)+'r</text>';
  edges.forEach(function(e){var p=e[0],c=e[1];var px=(p=='ROOT'?gx:pos[p]),py=(p=='ROOT'?Y(0)+15:Y(depth[p])+15);s+='<line x1="'+px+'" y1="'+py+'" x2="'+pos[c]+'" y2="'+(Y(depth[c])-15)+'" stroke="#ccc" stroke-width="1.6"/>';});
  Object.keys(nodes).forEach(function(nn){var h=hpOf[nn]||'H?';var rd=rdOf[nn]||0;var y=Y(depth[nn]);
    s+='<rect x="'+(pos[nn]-54)+'" y="'+(y-15)+'" width="108" height="30" rx="6" fill="'+col(h)+'" stroke="'+stk(h)+'" stroke-width="1.8"/>';
    s+='<text x="'+pos[nn]+'" y="'+(y-2)+'" text-anchor="middle" font-size="10" font-weight="700" fill="'+stk(h)+'">'+(h=='H?'?'HP未確認':h)+'</text>';
    s+='<text x="'+pos[nn]+'" y="'+(y+10)+'" text-anchor="middle" font-size="8" fill="#666">'+nn+' · '+rd+'r</text>';});
  return s+'</svg>';
}
function detail(x){
  var s='<div class="det">'+treeSVG(x);
  if(x.hpbad>0)s+='<div class="note" style="color:#d39e00">⚠ '+x.hpbad+' 個節點 HP 未確認/mixed（H?=somatic reads 無 confident germline-HP=HP3/unphased；mixed=ALT 跨 H1/H2=CN-gain multiplicity）</div>';
  if(x.conf)s+='<div class="note" style="color:#e67e22"><b>⚠ CN-gain dual-HP</b>：可能 multiplicity（兩等位都看似 ALT）非真兩 HP 獨立 lineage</div>';
  s+='<div class="note">determinacy='+x.det+' ｜ haplotypes='+x.hps+' ｜ tree_shape='+x.shape+' ｜ read-driven: 多-ALT '+nz(x.rdm)+' / combos '+nz(x.rdc)+' / maxchain '+nz(x.rdx)+'</div>';
  return s+'</div>';
}
var SK='r',SD=1,bd=document.getElementById('bd');
function tags(x){var t='';if(x.shape=='full_tree')t+='<span class="tag t-full">full</span>';if(x.dual)t+='<span class="tag t-dual">dual-HP</span>';if(x.conf)t+='<span class="tag t-conf">CN-gain?</span>';if(x.hpbad>0)t+='<span class="tag t-hp">HP?</span>';if(x.trunc)t+='<span class="tag t-trunc">截斷</span>';if(x.art&&!x.trunc)t+='<span class="tag t-art">artifact</span>';return t;}
function render(){
  var cn=fcn.value,ctx=fctx.value,sh=fshape.value,fl=fflag.value,q=fsearch.value.toLowerCase();
  var rows=D.filter(function(x){return (!cn||x.cn==cn)&&(!ctx||x.ctx==ctx)&&(!sh||x.shape==sh)&&(!q||(x.r+x.c).toLowerCase().indexOf(q)>=0)&&(!fl||(fl=='full'&&x.shape=='full_tree')||(fl=='dual'&&x.dual)||(fl=='conf'&&x.conf)||(fl=='hp'&&x.hpbad>0)||(fl=='art'&&x.art)||(fl=='trunc'&&x.trunc));});
  rows.sort(function(a,b){var u=a[SK],v=b[SK];if(u==null)u=-1;if(v==null)v=-1;return (u>v?1:u<v?-1:0)*SD;});
  document.getElementById('cnt').textContent=rows.length+' 區';
  var html=rows.slice(0,1500).map(function(x){return '<tr class="'+(x.art?'art':x.conf?'dualc':'')+'" data-i="'+D.indexOf(x)+'"><td>'+x.r+'</td><td class="num">'+x.n+'</td><td>'+x.cn+'</td><td>'+(x.ctx||'')+'</td><td>'+x.shape+'</td><td class="num">'+x.roots+'</td><td class="num">'+nz(x.rdx)+'</td><td class="num">'+(x.hpbad||'')+'</td><td class="num">'+x.tp+'/'+x.fp+'</td><td>'+tags(x)+'</td></tr>';}).join('');
  if(rows.length>1500)html+='<tr><td colspan="10" class="note">… 顯示前 1500/'+rows.length+'（縮篩選）</td></tr>';
  bd.innerHTML=html;
}
bd.addEventListener('click',function(e){var tr=e.target.closest('tr[data-i]');if(!tr)return;var nx=tr.nextElementSibling;if(nx&&nx.className=='dt'){nx.remove();return;}var d=document.createElement('tr');d.className='dt';d.innerHTML='<td colspan="10">'+detail(D[+tr.dataset.i])+'</td>';tr.after(d);});
['fcn','fctx','fshape','fflag','fsearch'].forEach(function(id){document.getElementById(id).addEventListener('input',render);});
document.querySelectorAll('th[data-k]').forEach(function(th){th.onclick=function(){var k=th.dataset.k;SD=(SK==k?-SD:1);SK=k;render();};});
render();
</script>'''

FOOT = '<p class="note" style="margin-top:18px">build_perregion_html.py v2（§13-A）· data: topology_per_region.json（3,885 區）+ rd_perregion.json · HCC1395 ⭐3 · 點列展開看 clone 樹 SVG（HP 上色）+ read-driven 確認</p></body></html>'

open(OUT, "w").write(HEAD + JS.replace("__DATA__", DATA_JSON) + FOOT)
print(f"DONE -> {OUT}  regions={n} full={n_full} struct={struct} dual={n_dual} conf={n_conf} hpbad_regions={n_hpbad} art={n_art} trunc={n_trunc}")
