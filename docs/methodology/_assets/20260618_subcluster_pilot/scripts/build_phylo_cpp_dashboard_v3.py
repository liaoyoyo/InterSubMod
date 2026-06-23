#!/usr/bin/env python3
"""[全位點觀察分類儀表板 v3] 對齊參考 HCC1395 FULL 的 5 層摺疊 + SVG 比例圖/直方/funnel + 最終判別 verdict/tier 篩選。
數字全由 records_v3 / dashboard_stats_v3 注入(SVG 由 Python 從 stats 生成, generator 不手打數字)。圖名 JS 端決定式(無 FIGMAP)。"""
import json, os

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
OUTD = f"{A}/obs_ws/cpp_wg"; os.makedirs(OUTD, exist_ok=True)
rec = json.load(open(f"{A}/phylo_cpp_wg_full_records_v3.json"))
S = json.load(open(f"{A}/dashboard_stats_v3.json"))
CHRS = [f"chr{c}" for c in range(1, 23)]; ck = {c: i for i, c in enumerate(CHRS)}
CNS = ["gain", "loss", "loh", "neutral"]; cni = {s: i for i, s in enumerate(CNS)}
WST = ["S1", "S2", "S3", "S4", "S6"]; wi = {s: i for i, s in enumerate(WST)}
ALS = ["aligned", "unaligned", "NA"]; ali = {s: i for i, s in enumerate(ALS)}

# figs symlink
figlink = f"{OUTD}/figs"
if not (os.path.islink(figlink) and os.path.realpath(figlink) == f"{A}/figs_cpp_wg_full"):
    try:
        if os.path.islink(figlink) or os.path.exists(figlink): os.remove(figlink)
        os.symlink(f"{A}/figs_cpp_wg_full", figlink)
    except Exception: pass
figset = set(os.listdir(f"{A}/figs_cpp_wg_full")) if os.path.isdir(f"{A}/figs_cpp_wg_full") else set()

# ---- 緊湊 D (v2 + wstate/tier/gdiv/als) ----
D = []
for r in rec:
    p = int(r["pos"]); snv = p + 5000
    fn = f"cpp_{r['chrom']}_{p}_{p+10000}.png"
    cv = r.get("cn_value")
    D.append([ck.get(r["chrom"], 0), snv, 0 if r["set"] == "TP" else 1, r["n"], r["coarse_ng"], r["fine_ng"],
              r["n_other"], 1 if r["unstable"] else 0, 1 if r["hidden_het"] else 0, r.get("V_hp", 0), r.get("V_allele", 0),
              cni.get(r.get("cn_state", "neutral"), 3), (cv if isinstance(cv, (int, float)) else -1), 1 if r.get("is_loh") else 0,
              1 if fn in figset else 0, r.get("n_tumor", -1) if r.get("n_tumor") is not None else -1,
              r.get("n_normal", -1) if r.get("n_normal") is not None else -1, r.get("geometry_ng", -1) if r.get("geometry_ng") is not None else -1,
              wi.get(r.get("wstate", "S6"), 4), r.get("structure_tier", 3), r.get("geom_divergence", 0), ali.get(r.get("align_status", "NA"), 2)])


# ---- SVG 生成器(從 stats 注入, 不手打數字) ----
def svg_stacked(segs, w=440, h=26):
    tot = sum(v for _, v, _ in segs) or 1; x = 0.0; p = []
    for lab, v, c in segs:
        ww = w * v / tot
        p.append(f'<rect x="{x:.1f}" y="0" width="{ww:.1f}" height="{h}" fill="{c}"/>')
        if ww > 40: p.append(f'<text x="{x+ww/2:.1f}" y="{h/2+4:.0f}" font-size="10" fill="#fff" text-anchor="middle">{lab} {100*v/tot:.0f}%</text>')
        x += ww
    return f'<svg viewBox="0 0 {w} {h}" width="100%" height="{h}" role="img">{"".join(p)}</svg>'


def svg_hist(dist, color="#0891b2", w=440, h=72):
    keys = sorted(dist, key=lambda k: float(k)); mx = max(dist.values()) or 1; bw = w / max(1, len(keys)); p = []
    for i, k in enumerate(keys):
        bh = h * 0.78 * dist[k] / mx; x = i * bw
        p.append(f'<rect x="{x+1:.1f}" y="{h-bh-12:.1f}" width="{bw-2:.1f}" height="{bh:.1f}" fill="{color}"/>')
        p.append(f'<text x="{x+bw/2:.1f}" y="{h-2:.0f}" font-size="8.5" fill="#9aa3b2" text-anchor="middle">{k}</text>')
    return f'<svg viewBox="0 0 {w} {h}" width="100%" height="{h}">{"".join(p)}</svg>'


def svg_funnel(stages, w=440):
    mx = stages[0]["n"] or 1; rh = 26; gap = 6; p = []
    for i, s in enumerate(stages):
        ww = w * s["n"] / mx; y = i * (rh + gap); x = (w - ww) / 2
        p.append(f'<rect x="{x:.1f}" y="{y}" width="{ww:.1f}" height="{rh}" rx="3" fill="#0891b2" opacity="{1-i*0.13:.2f}"/>')
        p.append(f'<text x="{w/2:.0f}" y="{y+rh/2+4:.0f}" font-size="10.5" fill="#fff" text-anchor="middle">{s["stage"]} {s["n"]:,}</text>')
    h = len(stages) * (rh + gap); return f'<svg viewBox="0 0 {w} {h}" width="100%" height="{h}">{"".join(p)}</svg>'


def svg_compare(pairs, w=300, h=92):
    mx = max(v for _, v, _ in pairs) or 1; slot = w / len(pairs); bw = slot * 0.5; p = []
    for i, (lab, v, c) in enumerate(pairs):
        bh = h * 0.66 * v / mx; x = i * slot + (slot - bw) / 2
        p.append(f'<rect x="{x:.1f}" y="{h-bh-16:.1f}" width="{bw:.1f}" height="{bh:.1f}" fill="{c}"/>')
        p.append(f'<text x="{x+bw/2:.1f}" y="{h-bh-20:.0f}" font-size="11" fill="#e6e8ec" text-anchor="middle">{v}%</text>')
        p.append(f'<text x="{x+bw/2:.1f}" y="{h-2:.0f}" font-size="9.5" fill="#9aa3b2" text-anchor="middle">{lab}</text>')
    return f'<svg viewBox="0 0 {w} {h}" width="100%" height="{h}">{"".join(p)}</svg>'


CNC = {"gain": "#ef4444", "loss": "#3b82f6", "loh": "#a855f7", "neutral": "#64748b"}
WC = {"S1": "#db2777", "S2": "#9333ea", "S3": "#d97706", "S4": "#64748b", "S6": "#334155"}


def stat_card(title, svg, note=""):
    return f'<div class="stcard"><div class="stt">{title}</div>{svg}<div class="sx">{note}</div></div>'


# L2 SVG 統計卡
cn_bar = svg_stacked([(s, S["cn_state"][s], CNC[s]) for s in CNS])
ws_bar = svg_stacked([(s, S["wstate"][s], WC[s]) for s in WST])
tn_bar = svg_stacked([("tumor", S["tn_total"]["tumor"], "#f97316"), ("normal", S["tn_total"]["normal"], "#22c55e")])
tier_bar = svg_stacked([("確認", S["tier"]["1"], "#16a34a"), ("邊緣", S["tier"]["2"], "#d97706"), ("無結構", S["tier"]["3"], "#475569")])
L2 = '<div class="statgrid">' + "".join([
    stat_card("SEQC2 CN 狀態", cn_bar, f"gain {S['cn_state']['gain']:,} / loh {S['cn_state']['loh']:,} / loss {S['cn_state']['loss']:,} / neutral {S['cn_state']['neutral']:,}"),
    stat_card("structure 5 態 (S1–S6)", ws_bar, "S1 確認對齊 / S2 不對齊 / S3 不穩定 / S4 次閾值 / S6 乾淨單群"),
    stat_card("structure tier", tier_bar, f"確認 {S['tier']['1']:,} / 邊緣 {S['tier']['2']:,} / 無 {S['tier']['3']:,}"),
    stat_card("tumor/normal read 組成（總）", tn_bar, f"tumor {S['tn_total']['tumor']:,} / normal {S['tn_total']['normal']:,}；中位 {S['median_nt']}/{S['median_nn']}"),
    stat_card("coarse 切群數 (null95)", svg_hist(S["coarse_dist"], "#db2777")),
    stat_card("fine 切群數 (null90)", svg_hist(S["fine_dist"], "#9333ea")),
    stat_card("幾何切群數 (0.7cut)", svg_hist(S["geom_dist"], "#0891b2"), "純幾何過切, 僅肉眼/null 對照, 非 subclone"),
    stat_card("判別 funnel", svg_funnel(S["funnel"])),
    stat_card("confident-multi: TP vs FP（反判別）", svg_compare([("TP", S["confident_multi_tp_pct"], "#60a5fa"), ("FP", S["confident_multi_fp_pct"], "#f87171")]), "FP&gt;TP = 反判別, confident-multi=cis-ASM 非 subclone"),
], ) + '</div>'

# L3 cross-tab tables
def xtab_cn_struct():
    h = '<table><tr><th>CN狀態</th><th>多群</th><th>單群</th><th>多群%</th></tr>'
    for s in CNS:
        m = S["cn_x_structure"][s]["multi"]; sg = S["cn_x_structure"][s]["single"]; t = m + sg
        h += f'<tr><td class="l"><span style="color:{CNC[s]}">●</span>{s}</td><td>{m:,}</td><td>{sg:,}</td><td>{round(100*m/t,1) if t else 0}</td></tr>'
    return h + '</table>'


def xtab_reads_wstate():
    h = '<table><tr><th>state</th>' + "".join(f'<th>{w}</th>' for w in WST) + '</tr><tr><td class="l">中位 n reads</td>' + "".join(f'<td>{S["reads_by_wstate"][w]}</td>' for w in WST) + '</tr></table>'
    return h


L3 = f'<div style="display:flex;gap:18px;flex-wrap:wrap"><div style="flex:1;min-width:300px"><b style="font-size:12px">CN 狀態 × 結構</b>{xtab_cn_struct()}</div><div style="flex:1;min-width:300px"><b style="font-size:12px">reads × structure 態</b>{xtab_reads_wstate()}</div></div>'

SECTIONS = f"""
<details open><summary>① 摘要 + 紅線</summary><div class="note">⚠ <b>觀察/刻畫層,非 subclone caller、非 TP/FP filter</b>。① 對齊(V_hp/V_allele)只敘述=cis-ASM 跡象<b>非 subclone</b>；② confident-multi 在 FP({S['confident_multi_fp_pct']}%)比 TP({S['confident_multi_tp_pct']}%)多=<b>反判別</b>；③ 幾何切群=純過切診斷非結構宣稱。單樣本 HCC1395 探索 ⭐2-3, 無外部真值。</div></details>
<details open><summary>② 整體統計觀察（SVG 比例圖）</summary>{L2}</details>
<details><summary>③ 分類交叉表</summary>{L3}</details>
"""

TPL = r"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>全位點觀察分類儀表板 v3</title><style>
:root{--bg:#0a0e16;--card:#111826;--fg:#e7edf5;--mut:#8b9bb4;--line:#2a3344;--ac:#D97757}
*{box-sizing:border-box}body{font-family:system-ui,"Noto Sans CJK TC",sans-serif;background:var(--bg);color:var(--fg);margin:0;line-height:1.5}
.wrap{max-width:1480px;margin:0 auto;padding:14px}
h1{font-size:19px;margin:2px 0}.sub{color:var(--mut);font-size:12.5px}
details{background:var(--card);border:1px solid var(--line);border-radius:8px;margin:7px 0;padding:3px 12px}
summary{cursor:pointer;font-weight:600;font-size:13.5px;padding:6px 0}
.note{background:#2a2410;border:1px solid #d97706;border-radius:7px;padding:8px 12px;font-size:12px;margin:6px 0}.note b{color:#fbbf24}
.statgrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(280px,1fr));gap:10px;margin:8px 0}
.stcard{background:#0b1222;border:1px solid var(--line);border-radius:7px;padding:8px 10px}.stt{font-size:12px;font-weight:600;margin-bottom:5px}.sx{color:var(--mut);font-size:10.5px;margin-top:4px}
table{border-collapse:collapse;font-size:12px;margin:4px 0;width:100%}td,th{border:1px solid var(--line);padding:3px 8px;text-align:right}th{background:#1a2130}td.l{text-align:left}
.bar{display:flex;flex-wrap:wrap;gap:7px;align-items:center;background:var(--card);border:1px solid var(--line);border-radius:8px;padding:9px;margin:8px 0;position:sticky;top:0;z-index:9}
.bar label{font-size:11.5px;color:var(--mut)}input,select{background:#0b1222;color:var(--fg);border:1px solid var(--line);border-radius:5px;padding:4px 6px;font-size:12px}
input[type=search]{width:150px}.btn{background:#1a2130;border:1px solid var(--line);border-radius:5px;color:var(--fg);padding:5px 9px;cursor:pointer;font-size:12px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(330px,1fr));gap:10px;margin:9px 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:9px;overflow:hidden}.card.j-in{border-color:#16a34a}.card.j-mb{border-color:#d97706}.card.j-ex{border-color:#dc2626}
.hd{padding:8px 10px}.ttl{font-weight:700;font-size:13px;cursor:pointer}.badges{margin-top:3px;display:flex;flex-wrap:wrap;gap:4px}
.badge{font-size:10.5px;padding:1px 6px;border-radius:9px;border:1px solid var(--line)}
.b-tp{background:#10243a;color:#7cc4ff}.b-fp{background:#3a1a1a;color:#ffa3a3}.b-gain{background:#3a1414;color:#ffb3b3;border-color:#ef4444}.b-loss{background:#10233a;color:#9fc6ff;border-color:#3b82f6}.b-loh{background:#241038;color:#d6b3ff;border-color:#a855f7}.b-neu{background:#1a2130;color:var(--mut)}
.b-al{background:#13261c;color:#9ff0c4}.b-tn{background:#0f2a24;color:#7fe3c4}.b-div{background:#3a2418;color:#ffba80;border-color:#d97706}
.b-t1{background:#13261c;color:#86efac;border-color:#16a34a}.b-t2{background:#2a2410;color:#fcd34d;border-color:#d97706}.b-t3{background:#1a2130;color:var(--mut)}
.desc{color:var(--mut);font-size:11px;margin-top:3px}
.figs{cursor:pointer;background:#0b1222;border-top:1px solid var(--line)}.figs img{width:100%;display:block}.nofig{padding:18px;text-align:center;color:var(--mut);font-size:11px}
.exp{padding:7px 10px;border-top:1px solid var(--line);font-size:11.5px}.exp .k{color:var(--mut)}
.jrow{display:flex;gap:5px;padding:7px 10px;border-top:1px solid var(--line)}.jb{flex:1;padding:4px;border:1px solid var(--line);border-radius:5px;background:#0b1222;color:var(--fg);cursor:pointer;font-size:11px}.jb.on-in{background:#14532d}.jb.on-mb{background:#3a2f10}.jb.on-ex{background:#3a1414}
.modal{display:none;position:fixed;inset:0;background:#000a;z-index:50;overflow:auto;padding:20px}.mc{max-width:1180px;margin:0 auto;background:var(--card);border:1px solid var(--ac);border-radius:10px;padding:14px}.mc img{width:100%;border-radius:6px;background:#fff}
.pager{text-align:center;margin:10px 0}.pager button{padding:5px 12px}
</style></head><body><div class="wrap">
<h1>全位點觀察分類儀表板 v3 — 完整數據 + SEQC2 CN + T/N + 切群數</h1>
<div class="sub">本 session full run（演化樹 tree-aware）+ SEQC2 CN + tumor/normal + 切群數(null95/null90/幾何)。路徑式, 需與同層 figs/ 一起開。__STAMP__</div>
__SECTIONS__
<div class="bar">
<label>🔍 <input type="search" id="q" placeholder="chr1:39078224 或 chr8"></label>
<label>判別 <select id="fws"><option value="">全</option><option value="0">S1 確認對齊</option><option value="1">S2 不對齊</option><option value="2">S3 不穩定</option><option value="3">S4 次閾值</option><option value="4">S6 乾淨單群</option></select></label>
<label>tier <select id="ftier"><option value="">全</option><option value="1">1 確認</option><option value="2">2 邊緣</option><option value="3">3 無</option></select></label>
<label>set <select id="fset"><option value="">全</option><option value="0">TP</option><option value="1">FP</option></select></label>
<label>CN <select id="fcn"><option value="">全</option><option value="0">gain</option><option value="1">loss</option><option value="2">loh</option><option value="3">neutral</option></select></label>
<label>CN值 <input type="number" id="cnmin" placeholder="min" style="width:54px" step="0.5">-<input type="number" id="cnmax" placeholder="max" style="width:54px" step="0.5"></label>
<label><input type="checkbox" id="floh"> LOH</label>
<label>coarse <input type="number" id="gmin" value="1" style="width:42px">-<input type="number" id="gmax" value="9" style="width:42px"></label>
<label>n≥ <input type="number" id="nmin" placeholder="0" style="width:54px"></label>
<label><input type="checkbox" id="funs"> unstable</label>
<label><input type="checkbox" id="ffc"> fine&gt;coarse</label>
<label><input type="checkbox" id="fdiv"> 可疑漏切</label>
<label>排序 <select id="sort"><option value="1">位置</option><option value="3">n</option><option value="4">coarse</option><option value="17">幾何</option><option value="12">CN值</option><option value="vmax">對齊V</option></select><button class="btn" id="sdir">▼</button></label>
<span class="sub" id="cnt"></span><button class="btn" id="ej">匯出判讀</button>
</div>
<div id="grid" class="grid"></div><div class="pager" id="pager"></div>
<details><summary>⑤ 參考 / 方法 / provenance</summary><div class="sub" style="padding:6px 0">
卡片圖：左 UPGMA 樹(tree-aware 群色) ｜ 中 甲基 read×CpG(RdBu_r) ｜ 右 read×read 距離(magma, 暗=近)。側欄 HP|ALT|T/N|Strand。<br>
判別 verdict = weak-structure 5 態(S1 確認多群對齊=cis-ASM 候選 / S2 不對齊 / S3 不穩定 / S4 次閾值候選 chance-level / S6 乾淨單群)。對齊=cis-ASM 跡象只敘述不歸類。切群數三法: null95(嚴主判)/null90(寬候選)/幾何0.7cut(肉眼過切)。<br>
來源: phylo_cpp_wg_full_records_v3.json + dashboard_stats_v3.json + SEQC2 CNV ngs_benchmark(hg38) + reads.tsv(is_tumor)。§13 數字注入不手打。單樣本 HCC1395 ⭐2-3。
</div></details>
</div>
<div class="modal" id="modal" onclick="if(event.target===this)closeM()"><div class="mc" id="mc"></div></div>
<script>
const D=__D__, CHRS=__CHRS__, CNS=['gain','loss','loh','neutral'], WST=['S1','S2','S3','S4','S6'];
const C={ck:0,ps:1,set:2,n:3,cg:4,fg:5,oth:6,uns:7,hh:8,vhp:9,val:10,cns:11,cnv:12,loh:13,fig:14,nt:15,nn:16,geom:17,wst:18,tier:19,gdiv:20,als:21};
const $=id=>document.getElementById(id);
const figName=r=>'cpp_'+CHRS[r[C.ck]]+'_'+(r[C.ps]-5000)+'_'+(r[C.ps]+5000)+'.png';
const key=r=>CHRS[r[C.ck]]+':'+r[C.ps];
let J=JSON.parse(localStorage.getItem('phylo_v3_judge')||'{}'),page=0,sdir=-1,PER=24;
function filt(){let a=D.slice();
 const q=$('q').value.trim().toLowerCase();
 if(q)a=a.filter(r=>(CHRS[r[C.ck]]+':'+r[C.ps]).toLowerCase().includes(q)||CHRS[r[C.ck]].toLowerCase()===q);
 const fw=$('fws').value;if(fw!=='')a=a.filter(r=>r[C.wst]==fw);
 const ft=$('ftier').value;if(ft!=='')a=a.filter(r=>r[C.tier]==ft);
 const fs=$('fset').value;if(fs!=='')a=a.filter(r=>r[C.set]==fs);
 const fc=$('fcn').value;if(fc!=='')a=a.filter(r=>r[C.cns]==fc);
 if($('cnmin').value!=='')a=a.filter(r=>r[C.cnv]>=0&&r[C.cnv]>=+$('cnmin').value);if($('cnmax').value!=='')a=a.filter(r=>r[C.cnv]>=0&&r[C.cnv]<=+$('cnmax').value);
 if($('floh').checked)a=a.filter(r=>r[C.loh]);
 a=a.filter(r=>r[C.cg]>=+$('gmin').value&&r[C.cg]<=+$('gmax').value);
 if($('nmin').value!=='')a=a.filter(r=>r[C.n]>=+$('nmin').value);
 if($('funs').checked)a=a.filter(r=>r[C.uns]);if($('ffc').checked)a=a.filter(r=>r[C.fg]>r[C.cg]);if($('fdiv').checked)a=a.filter(r=>r[C.gdiv]);
 const sk=$('sort').value;a.sort((x,y)=>{const f=v=>sk==='vmax'?Math.max(v[C.vhp],v[C.val]):(sk==1?v[C.ck]*1e10+v[C.ps]:v[+sk]);return sdir*(f(x)-f(y));});
 return a;}
function cnBadge(r){const st=CNS[r[C.cns]];const cls={0:'b-gain',1:'b-loss',2:'b-loh',3:'b-neu'}[r[C.cns]];return `<span class="badge ${cls}">${st==='loh'?'LOH':st+(r[C.cnv]>=0?' CN'+r[C.cnv]:'')}</span>`;}
function tierBadge(r){const t=r[C.tier];return `<span class="badge b-t${t}">${WST[r[C.wst]]}·tier${t}</span>`;}
function alignDesc(r){const vh=r[C.vhp],va=r[C.val];return Math.max(vh,va)>=0.3?`對齊 ${va>=vh?'allele':'hp'} 軸 (V_hp ${vh}/V_al ${va}) — cis-ASM 跡象，非 subclone`:`未對齊 germline 軸 (V_hp ${vh}/V_al ${va})`;}
function structText(r){return {0:'確認多群·對齊 germline（cis-ASM 候選，非 subclone）',1:'確認多群·不對齊（結構但無標籤對應）',2:'不穩定多群（seed 分歧 borderline）',3:'次閾值候選（僅 null90 過，chance-level）',4:'乾淨單群（無監督切不出≥2）'}[r[C.wst]];}
function clusterText(r){const cg=r[C.cg],fg=r[C.fg],g=r[C.geom];const v=(g>=2&&cg<2)?`<span style="color:#ffba80">幾何 ${g} 群但嚴閘判單群=可疑漏切</span>`:(fg>cg?`寬閘多切（${cg}→${fg}）`:'各法一致');return `null95 <b>${cg}</b> ｜ null90 <b>${fg}</b> ｜ 幾何 <b>${g>=0?g:'?'}</b> → ${v}`;}
function card(r){const k=key(r),fn=r[C.fig]?figName(r):null;const j=(J[k]||{}).c;
 return `<div class="card${j?' j-'+j:''}"><div class="hd"><div class="ttl" onclick="openM('${k}')">${k}</div>
  <div class="badges">${r[C.set]===0?'<span class="badge b-tp">TP</span>':'<span class="badge b-fp">FP</span>'}${tierBadge(r)}${cnBadge(r)}<span class="badge b-tn">T${r[C.nt]>=0?r[C.nt]:'?'}/N${r[C.nn]>=0?r[C.nn]:'?'}</span><span class="badge ${r[C.gdiv]?'b-div':'b-al'}">切群 ${r[C.cg]}/${r[C.fg]}/${r[C.geom]>=0?r[C.geom]:'?'}${r[C.gdiv]?' ⚠':''}</span>${r[C.loh]?'<span class="badge b-loh">LOH</span>':''}</div>
  <div class="desc">${alignDesc(r)}</div></div>
  ${fn?`<div class="figs" onclick="openM('${k}')"><img loading="lazy" src="figs/${fn}"></div>`:'<div class="nofig">無圖（reads 過少）</div>'}
  <details class="exp"><summary>歸類 + 切群數 + tumor/normal</summary>
   <div><span class="k">歸類:</span> <b>${structText(r)}</b></div>
   <div><span class="k">切群數:</span> ${clusterText(r)}</div>
   <div><span class="k">tumor/normal:</span> tumor <b>${r[C.nt]>=0?r[C.nt]:'?'}</b>/normal <b>${r[C.nn]>=0?r[C.nn]:'?'}</b>（phylo 用 tumor ${r[C.n]} 條）</div>
   <div><span class="k">SEQC2:</span> ${CNS[r[C.cns]]}${r[C.cnv]>=0?' CN='+r[C.cnv]:''}${r[C.loh]?' LOH':''} ｜ other ${r[C.oth]} ｜ unstable ${r[C.uns]?'是':'否'} ｜ hidden-het ${r[C.hh]?'是':'否'}</div></details>
  <div class="jrow"><button class="jb${j==='in'?' on-in':''}" onclick="setJ('${k}','in')">應含</button><button class="jb${j==='mb'?' on-mb':''}" onclick="setJ('${k}','mb')">可能漏</button><button class="jb${j==='ex'?' on-ex':''}" onclick="setJ('${k}','ex')">排除</button></div></div>`;}
function draw(){const a=filt();const pages=Math.max(1,Math.ceil(a.length/PER));page=Math.min(page,pages-1);
 $('cnt').textContent=a.length.toLocaleString()+' 個';
 $('grid').innerHTML=a.slice(page*PER,page*PER+PER).map(card).join('')||'<div class="sub">無符合</div>';
 $('pager').innerHTML=`<button class="btn" ${page<=0?'disabled':''} onclick="page--;draw()">‹</button> 第 ${page+1}/${pages} 頁 <button class="btn" ${page>=pages-1?'disabled':''} onclick="page++;draw()">›</button>`;}
['q','fws','ftier','fset','fcn','cnmin','cnmax','floh','gmin','gmax','nmin','funs','ffc','fdiv','sort'].forEach(id=>$(id).addEventListener('input',()=>{page=0;draw();}));
$('sdir').onclick=()=>{sdir*=-1;$('sdir').textContent=sdir===-1?'▼':'▲';page=0;draw();};
window.setJ=(k,c)=>{J[k]=J[k]||{};J[k].c=(J[k].c===c?null:c);localStorage.setItem('phylo_v3_judge',JSON.stringify(J));draw();};
window.openM=k=>{const r=D.find(x=>key(x)===k);if(!r)return;const fn=r[C.fig]?figName(r):null;
 $('mc').innerHTML=`<div style="display:flex;justify-content:space-between"><h2 style="margin:0">${k} ${tierBadge(r)}</h2><button class="btn" onclick="closeM()">關閉 ✕</button></div>
  <div class="desc" style="margin:6px 0">${alignDesc(r)}</div>${fn?`<img src="figs/${fn}">`:'<div class="nofig">無圖</div>'}
  <div class="exp" style="border:none"><div><span class="k">歸類:</span> <b>${structText(r)}</b></div><div><span class="k">切群數:</span> ${clusterText(r)}</div>
   <div><span class="k">tumor/normal:</span> tumor <b>${r[C.nt]>=0?r[C.nt]:'?'}</b>/normal <b>${r[C.nn]>=0?r[C.nn]:'?'}</b>（phylo 用 tumor ${r[C.n]} 條）</div>
   <div><span class="k">SEQC2 CN:</span> ${CNS[r[C.cns]]}${r[C.cnv]>=0?' CN='+r[C.cnv]:''}${r[C.loh]?' LOH':''} ｜ set ${r[C.set]===0?'TP':'FP'} ｜ 對齊 V_hp ${r[C.vhp]}/V_al ${r[C.val]} ｜ unstable ${r[C.uns]?'是':'否'} ｜ hidden-het ${r[C.hh]?'是':'否'}</div></div>
  <div class="sub" style="margin-top:8px">距離對角塊乾淨+coarse 側欄對齊 HP/ALT=真分群。對齊=cis-ASM 非 subclone。</div>`;$('modal').style.display='block';};
window.closeM=()=>$('modal').style.display='none';document.addEventListener('keydown',e=>{if(e.key==='Escape')closeM();});
$('ej').onclick=()=>{const o=Object.entries(J).filter(([k,v])=>v.c).map(([k,v])=>({locus:k,choice:v.c}));const b=new Blob([JSON.stringify({n:o.length,judgments:o},null,1)],{type:'application/json'});const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download='phylo_v3_judgments.json';a.click();};
draw();
</script></body></html>"""

out = (TPL.replace("__D__", json.dumps(D, separators=(",", ":")))
          .replace("__CHRS__", json.dumps(CHRS))
          .replace("__SECTIONS__", SECTIONS)
          .replace("__STAMP__", "2026-06-23 v3 · SEQC2 CNV(hg38) + reads.tsv T/N · single-sample ⭐2-3"))
outp = f"{OUTD}/20260623_phylo_cpp_observation_dashboard_v3.html"
open(outp, "w").write(out)
print(f"WROTE {outp} ({len(out)//1024} KB) | {len(D)} loci")
