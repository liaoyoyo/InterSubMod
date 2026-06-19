#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_workstation.py — reusable generator for an INTERACTIVE per-item
judgment + verification HTML workstation (/verify-workstation skill).

It turns a single declarative spec (items + metrics + figures + changelog +
provenance) into a self-contained HTML where a HUMAN records per-item verdicts
(agree / doubt / reject) + reason, persisted to localStorage, with JSON/CSV
export-import, ⌖ provenance src badges, an authoritative-vs-draft / correction
changelog, and §13-A construction-level anti-fabrication.

§13-A (CLAUDE.md): every metric is INJECTED from the spec by key. A required
metric whose value is null/missing makes the build REFUSE (exit 3) instead of
silently rendering a dash — the exact failure mode that froze the page-04 Fig4.
Classification must be COMPUTED in the data-prep step (never hand-assigned) and
arrive in the spec already resolved.

Usage:
    python3 build_workstation.py <workstation_spec.json> [-o out.html]

Spec schema: see ../templates/workstation_spec.schema.json
Worked examples: ../examples/README.md (chr2:18M inline-SVG ; HCC1395 ASM external-PNG)

Scale note: this generator targets the COMMON case (<= a few hundred rich cards,
inline SVG or external PNG). For genome-scale (>~1000 items) use the compact
column-array pattern documented in ../references/section_modules.md (e.g. the
HCC1395 ASM workstation script 102) — same UI core, denser data payload.
"""
import json, sys, os, html, argparse, datetime

EXIT_REFUSE = 3


class Refuse(RuntimeError):
    pass


def esc(s):
    return html.escape(str(s))


def req(value, name):
    """§13-A: REFUSE the build on a missing REQUIRED value (no silent dash)."""
    if value is None or (isinstance(value, float) and value != value):
        raise Refuse(f"REQUIRED value missing -> refuse to render: {name}")
    return value


def get(d, k, default=None):
    return d.get(k, default) if isinstance(d, dict) else default


# ---------------------------------------------------------------- CSS / JS (raw, no f-string)
CSS = """
:root{--bg:#0f172a;--card:#1e293b;--mut:#94a3b8;--fg:#e2e8f0;--ac:#38bdf8;
--ok:#22c55e;--warn:#f59e0b;--no:#ef4444;--line:#334155;--gold:#f2c94c}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:14px/1.6 system-ui,-apple-system,'Segoe UI','Noto Sans CJK TC','Microsoft JhengHei',sans-serif}
header{padding:16px 24px;border-bottom:1px solid var(--line);background:#0b1222;position:sticky;top:0;z-index:40}
h1{margin:0 0 3px;font-size:19px}h2{font-size:16px;margin:24px 0 8px;border-left:3px solid var(--ac);padding-left:9px}
.sub{color:var(--mut);font-size:12px}.wrap{max-width:1320px;margin:0 auto;padding:0 20px 90px}
.panel{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px;margin:12px 0}
.banner{margin-top:8px;background:#3a1d05;border:1px solid #92400e;border-radius:7px;padding:7px 12px;font-size:12px;color:#fcd34d}
.row{display:flex;gap:14px;flex-wrap:wrap;align-items:center}
button,select,input[type=file]{background:#0b1222;color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:6px 10px;font-size:12.5px;cursor:pointer}
.stat4{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:10px}
.sb{background:#0b1222;border:1px solid var(--line);border-radius:8px;padding:9px 12px}.sb .n{font-size:20px;font-weight:800}.sb .l{font-size:11px;color:var(--mut)}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(310px,1fr));gap:12px}
.card{background:var(--card);border:1px solid var(--line);border-radius:9px;overflow:hidden}
.card.j-in{border-color:var(--ok)}.card.j-ex{border-color:var(--no)}.card.j-doubt{border-color:var(--warn)}
.card .hd{padding:8px 11px;border-bottom:1px solid var(--line);cursor:pointer}
.card .ttl{font-weight:700;font-size:13px}.card .meta{font-size:10.5px;color:var(--mut);margin-top:3px}
.kv{display:flex;justify-content:space-between;gap:8px;font-size:11.5px;padding:1px 11px}
.kv .k{color:var(--mut)}.kv .v{font-variant-numeric:tabular-nums}
.figs{display:flex;gap:2px;background:#0b1222;cursor:pointer}.figs img{width:50%;display:block;background:#fff;min-height:40px}
.figsvg{background:#11161f;padding:4px;cursor:pointer}.figsvg svg{width:100%;display:block}
.nofig{padding:12px;color:#64748b;font-size:11px;text-align:center}
.badge{display:inline-block;font-size:9.5px;padding:1px 7px;border-radius:9px;font-weight:700;margin-left:3px}
.b-ok{background:#14532d;color:#86efac}.b-warn{background:#78350f;color:#fcd34d}.b-no{background:#450a0a;color:#fca5a5}
.b-info{background:#0c2a4d;color:#60a5fa}.b-mut{background:#1e293b;color:#cbd5e1}.b-gold{background:#3a300a;color:#f2c94c}
.b-src{background:#0b1222;color:#7c8aa0;font-weight:500;font-size:10px;border:1px solid var(--line);cursor:help}
.jrow{display:flex;gap:4px;padding:6px 8px;border-top:1px solid var(--line)}
.jrow button{flex:1;font-size:10.5px;padding:4px 2px}
.jb.on-in{background:#14532d;color:#86efac;border-color:#22c55e}.jb.on-doubt{background:#78350f;color:#fcd34d;border-color:#f59e0b}.jb.on-ex{background:#450a0a;color:#fca5a5;border-color:#ef4444}
.rsel{width:100%;padding:5px 8px;font-size:10.5px;border-top:1px solid var(--line);border-radius:0}
table{border-collapse:collapse;width:100%;font-size:11.5px;margin:4px 0}
th,td{border:1px solid var(--line);padding:4px 8px;text-align:left}th{background:#0b1222;color:var(--mut)}
.note{background:#1c1407;border:1px solid #78350f;border-radius:8px;padding:9px 13px;font-size:12px;color:#fcd34d;margin:10px 0}
code{background:#0b1222;padding:1px 5px;border-radius:4px;font-size:11.5px;color:#fcd34d}
#modal{position:fixed;inset:0;background:rgba(0,0,0,.86);display:none;z-index:100;overflow:auto;padding:22px}
#modal .mc{max-width:1040px;margin:0 auto;background:var(--card);border:1px solid var(--ac);border-radius:12px;padding:18px}
#modal img{width:100%;background:#fff;border-radius:6px}.mfigs{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin:12px 0}
.statg{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:7px}
.st{background:#0b1222;border:1px solid var(--line);border-radius:6px;padding:6px 9px}.st .l{font-size:10px;color:var(--mut)}.st .n{font-size:14px;font-weight:700}
footer{margin-top:30px;border-top:1px solid var(--line);padding-top:12px;color:var(--mut);font-size:12px}
.pager{display:flex;gap:8px;align-items:center;justify-content:center;margin:14px 0}
"""

JS = r"""
const SPEC=__SPEC__;
const $=id=>document.getElementById(id);
const ITEMS=SPEC.items, IC=SPEC.item_config||{};
const LSKEY=SPEC.meta.lskey||'verify_ws_v1';
let J=JSON.parse(localStorage.getItem(LSKEY)||'{}');
function saveJ(){localStorage.setItem(LSKEY,JSON.stringify(J));}
const VSTATES=IC.verdict_states||[{key:'in',label:'同意'},{key:'doubt',label:'存疑'},{key:'ex',label:'否定'}];
const REASONS=IC.reason_options||[{key:'',label:'— 原因（可選）—'},{key:'criteria',label:'標準不好'},{key:'method',label:'方法問題'},{key:'visual',label:'沒看清楚'},{key:'other',label:'其他'}];
function setJ(id,c){const cur=J[id]||{};cur.c=(cur.c===c?'un':c);cur.t=new Date().toISOString().slice(0,19);J[id]=cur;saveJ();draw();prog();}
function setR(id,r){const cur=J[id]||{c:'un'};cur.r=r;cur.t=new Date().toISOString().slice(0,19);J[id]=cur;saveJ();}
window.setJ=setJ;window.setR=setR;

function prog(){const v=Object.values(J);const n={};for(const s of VSTATES)n[s.key]=v.filter(x=>x.c===s.key).length;
 const un=ITEMS.length-Object.values(n).reduce((a,b)=>a+b,0);
 const cells=VSTATES.map(s=>`<div class="sb"><div class="n">${n[s.key]||0}</div><div class="l">${s.label}</div></div>`).join('')+
  `<div class="sb"><div class="n">${un}</div><div class="l">未判讀</div></div>`;
 $('prog').innerHTML=cells;}

function badge(b){const cls={ok:'b-ok',warn:'b-warn',no:'b-no',info:'b-info',mut:'b-mut',gold:'b-gold'}[b.kind||'mut']||'b-mut';
 return `<span class="badge ${cls}"${b.title?` title="${b.title}"`:''}>${b.label}</span>`;}
function srcBadge(src){return src?`<span class="b-src badge" title="provenance: ${src}">⌖</span>`:'';}
function figHTML(it,big){
 const f=it.figures;if(!f)return '<div class="nofig">無圖</div>';
 if(f.mode==='svg')return `<div class="figsvg" onclick="openM('${it.id}')">${f.svg}</div>`;
 if(f.mode==='png'){const imgs=f.imgs||[];if(!imgs.length)return '<div class="nofig">無圖</div>';
  if(big)return `<div class="mfigs">`+imgs.map(g=>`<div><div class="sub" style="text-align:center">${g.caption||''}</div><img src="${g.path}"></div>`).join('')+`</div>`;
  return `<div class="figs" onclick="openM('${it.id}')">`+imgs.map(g=>`<img loading="lazy" src="${g.path}">`).join('')+`</div>`;}
 return '<div class="nofig">無圖</div>';}
function metaLine(it){return (it.meta_metrics||[]).map(k=>{const m=(it.metrics||[]).find(x=>x.k===k);return m?`${m.k}=${m.v}`:'';}).filter(Boolean).join(' ｜ ');}
function jbtn(id,c,lab){const on=((J[id]||{}).c===c)?(' on-'+c):'';return `<button class="jb${on}" onclick="event.stopPropagation();setJ('${id}','${c}')">${lab}</button>`;}
function card(it){const j=J[it.id]||{};const cc=j.c==='in'?' j-in':j.c==='ex'?' j-ex':j.c==='doubt'?' j-doubt':'';
 const badges=(it.badges||[]).map(badge).join('');
 const kv=(it.metrics||[]).filter(m=>m.show!==false).map(m=>`<div class="kv"><span class="k">${m.k} ${srcBadge(m.src)}</span><span class="v">${m.v}</span></div>`).join('');
 const reasons=REASONS.map(r=>`<option value="${r.key}"${(j.r||'')===r.key?' selected':''}>${r.label}</option>`).join('');
 return `<div class="card${cc}" data-id="${it.id}">
   <div class="hd" onclick="openM('${it.id}')"><div class="ttl">${it.title} ${badges}</div>${it.subtitle?`<div class="meta">${it.subtitle}</div>`:''}${metaLine(it)?`<div class="meta">${metaLine(it)}</div>`:''}</div>
   ${figHTML(it,false)}
   <div style="padding:4px 0">${kv}</div>
   <div class="jrow">${VSTATES.map(s=>jbtn(it.id,s.key,s.label)).join('')}<button class="jb${(J[it.id]||{}).c==='un'?' on-mut':''}" onclick="event.stopPropagation();setJ('${it.id}','un')">清除</button></div>
   <select class="rsel" onchange="event.stopPropagation();setR('${it.id}',this.value)">${reasons}</select>
 </div>`;}

let page=0;const PER=IC.per_page||24;
function items(){let a=ITEMS.slice();const jf=$('jf')?$('jf').value:'all';
 if(jf!=='all')a=a.filter(it=>{const c=(J[it.id]||{}).c||'un';return c===jf;});
 return a;}
function draw(){const a=items();const pages=Math.max(1,Math.ceil(a.length/PER));page=Math.min(page,pages-1);
 $('grid').innerHTML=a.slice(page*PER,page*PER+PER).map(card).join('')||'<div class="sub">無符合</div>';
 $('gcount').textContent=a.length+' 個';
 $('pager').innerHTML=`<button ${page<=0?'disabled':''} onclick="page--;draw()">‹</button><span>第 ${page+1}/${pages} 頁</span><button ${page>=pages-1?'disabled':''} onclick="page++;draw()">›</button>`;}
if($('jf'))$('jf').addEventListener('input',()=>{page=0;draw();});

window.openM=id=>{const it=ITEMS.find(x=>x.id===id);if(!it)return;
 const stats=(it.modal_stats||it.metrics||[]).map(m=>`<div class="st"><div class="l">${m.k} ${srcBadge(m.src)}</div><div class="n">${m.v}</div></div>`).join('');
 $('mc').innerHTML=`<div class="row" style="justify-content:space-between"><h2 style="border:none;margin:0">${it.title} ${(it.badges||[]).map(badge).join('')}</h2><button onclick="closeM()">關閉 ✕</button></div>
  ${figHTML(it,true)}<div class="statg">${stats}</div>
  <div class="jrow" style="margin-top:10px">${VSTATES.map(s=>jbtn(it.id,s.key,s.label)).join('')}</div>
  ${it.reading_guide?`<div class="sub" style="margin-top:8px">${it.reading_guide}</div>`:''}`;
 $('modal').style.display='block';};
window.closeM=()=>$('modal').style.display='none';
document.addEventListener('keydown',e=>{if(e.key==='Escape')closeM();});

function rows4export(){return ITEMS.filter(it=>{const c=(J[it.id]||{}).c;return c&&c!=='un';}).map(it=>{const j=J[it.id];
 return {id:it.id,title:it.title,verdict:j.c,reason:j.r||'',ts:j.t||''};});}
window.exportJ=()=>{const o={spec:SPEC.meta.title,build:SPEC.meta.build_commit,generated:new Date().toISOString(),judgments:rows4export()};
 dl((SPEC.meta.lskey||'verdicts')+'_'+new Date().toISOString().slice(0,10)+'.json',JSON.stringify(o,null,1),'application/json');};
window.exportCSV=()=>{const r=rows4export(),h=['id','title','verdict','reason','ts'];
 const csv=[h.join(',')].concat(r.map(x=>h.map(k=>JSON.stringify(x[k]??'')).join(','))).join('\n');
 dl((SPEC.meta.lskey||'verdicts')+'_'+new Date().toISOString().slice(0,10)+'.csv',csv,'text/csv');};
function dl(n,t,m){const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([t],{type:m}));a.download=n;a.click();}
window.importJ=e=>{const f=e.target.files[0];if(!f)return;const rd=new FileReader();
 rd.onload=()=>{try{const o=JSON.parse(rd.result),arr=o.judgments||o;let n=0;
  for(const x of arr){if(x.id&&x.verdict){J[x.id]={c:x.verdict,r:x.reason||'',t:x.ts||''};n++;}}
  saveJ();draw();prog();alert('匯入 '+n+' 筆');}catch(err){alert('匯入失敗：'+err);}};rd.readAsText(f);};
window.clearJ=()=>{if(confirm('清空所有判讀？')){J={};saveJ();draw();prog();}};
prog();draw();
"""


def render_changelog(cl):
    if not cl:
        return ""
    sb = {"in-HEAD": ("#052e16", "#4ade80", "✓ 已落地"),
          "planned-phase2": ("#3a1d05", "#fcd34d", "◔ 計畫中"),
          "not-done": ("#450a0a", "#fca5a5", "✗ 尚未做"),
          "authoritative": ("#052e16", "#86efac", "✓ 權威"),
          "superseded": ("#450a0a", "#fca5a5", "✗ 已過時")}
    rows = ""
    for c in cl.get("corrections", []):
        s = sb.get(c.get("status"), ("#1e293b", "#cbd5e1", esc(c.get("status", ""))))
        rows += (f'<tr><td><b>{esc(c.get("id",""))}</b></td><td>{esc(c.get("what",""))}</td>'
                 f'<td><span class="badge" style="background:{s[0]};color:{s[1]}">{s[2]}</span></td>'
                 f'<td>{esc(c.get("effect",""))}</td><td><code>{esc(c.get("src",""))}</code></td></tr>')
    head = ""
    for key, lab, col in (("data_status", "資料狀態", "--warn"), ("binary_status", "程式/版本", "--ac"),
                          ("phase", "階段", "--fg")):
        v = cl.get(key)
        if isinstance(v, dict):
            txt = v.get("summary") or " → ".join(str(x) for x in v.values() if isinstance(x, str))
            head += f'<div><b style="color:var({col})">{lab}：</b>{esc(txt)}</div>'
    concl = cl.get("audit_conclusion", {})
    note = f'<div class="note"><b>結論：</b>{esc(concl.get("summary",""))}</div>' if concl else ""
    return (f'<h2>⓪b 修正過程紀錄 / changelog</h2><div class="panel">'
            f'<div class="sub" style="line-height:1.9;margin-bottom:8px">{head}</div>'
            f'<table><tr><th>ID</th><th>內容</th><th>狀態</th><th>影響</th><th>來源(已驗證)</th></tr>{rows}</table>{note}</div>')


def render_sections(sections):
    out = ""
    for s in sections or []:
        out += f'<h2>{esc(s.get("title",""))}</h2><div class="panel">{s.get("html","")}</div>'
    return out


def build(spec):
    meta = spec["meta"]
    req(meta.get("title"), "meta.title")
    req(spec.get("items"), "items")
    # validate required item fields + required metrics (§13-A)
    reqd = (spec.get("item_config") or {}).get("required_metrics", [])
    for i, it in enumerate(spec["items"]):
        req(it.get("id"), f"items[{i}].id")
        req(it.get("title"), f"items[{i}].title")
        present = {m["k"] for m in it.get("metrics", []) if "v" in m and m["v"] is not None}
        for rk in reqd:
            if rk not in present:
                raise Refuse(f"items[{it.get('id')}] missing required metric '{rk}' (§13-A refuse-on-missing)")

    banner = ""
    if spec.get("banner"):
        banner = f'<div class="banner">{spec["banner"].get("text","")}</div>'
    elif (spec.get("changelog") or {}).get("data_status"):
        ds = spec["changelog"]["data_status"]
        dstxt = ds.get("summary") if isinstance(ds, dict) else ds
        banner = f'<div class="banner">⚠ 資料狀態：{esc(dstxt)}　build <code>{esc(meta.get("build_commit","?"))}</code></div>'

    jf = ('<div class="row"><span class="sub">篩選判讀</span>'
          '<select id="jf"><option value="all">全部</option>'
          + "".join(f'<option value="{s["key"]}">{esc(s["label"])}</option>'
                    for s in (spec.get("item_config") or {}).get("verdict_states",
                              [{"key": "in", "label": "同意"}, {"key": "doubt", "label": "存疑"}, {"key": "ex", "label": "否定"}]))
          + '<option value="un">未判讀</option></select><span id="gcount" class="sub"></span></div>')

    prov = spec.get("provenance", {})
    prov_items = "".join(f"<li><code>{esc(s)}</code></li>" for s in prov.get("sources", []))
    footer = (f'<footer><b>Provenance</b> · build <code>{esc(meta.get("build_commit","?"))}</code>'
              f' · branch <code>{esc(meta.get("build_branch","?"))}</code> · 生成 {datetime.date.today().isoformat()}'
              f' · tier {esc(meta.get("tier","?"))} · scope {esc(meta.get("scope","?"))}<br>'
              f'{esc(prov.get("note",""))}<ul>{prov_items}</ul>'
              f'所有 metric 由 spec 注入（§13-A req()-refuse），人工判讀存 localStorage 並可匯出。</footer>')

    body = f"""<header>
 <h1>{esc(meta["title"])}</h1>
 <div class="sub">{esc(meta.get("subtitle",""))}</div>
 {banner}
</header>
<div class="wrap">
{render_changelog(spec.get("changelog"))}
{render_sections(spec.get("sections"))}
<h2>判讀紀錄 — 匯出 / 匯入 / 進度</h2>
<div class="panel">
 <div class="row">
  <button onclick="exportJ()">⬇ 匯出 JSON</button>
  <button onclick="exportCSV()">⬇ 匯出 CSV</button>
  <label class="row" style="gap:5px"><span class="sub">⬆ 匯入再校正</span><input type="file" accept=".json" onchange="importJ(event)"></label>
  <button onclick="clearJ()">✕ 清空</button>
 </div>
 <div class="stat4" id="prog" style="margin-top:12px"></div>
</div>
<h2>逐項判讀（每項可記 {esc("/".join(s["label"] for s in (spec.get("item_config") or {}).get("verdict_states", [{"label":"同意"},{"label":"存疑"},{"label":"否定"}])))} + 原因）</h2>
<div class="panel">
 {jf}
 <div class="grid" id="grid"></div>
 <div class="pager" id="pager"></div>
</div>
{footer}
</div>
<div id="modal" onclick="if(event.target.id==='modal')closeM()"><div class="mc" id="mc"></div></div>"""

    # Optional caller-supplied JS, injected ONCE at end of body and POPPED before
    # serialization so it never lands inside the <script type="application/json"> spec
    # block (a nested </script> there would close that block early). Additive: specs
    # without extra_js are byte-identical to before. Must not itself contain "</script>".
    extra_js = spec.pop("extra_js", None)
    extra = f"<script>{extra_js}</script>" if extra_js else ""
    spec_json = json.dumps(spec, ensure_ascii=False)
    return f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(meta["title"])}</title>
<style>{CSS}</style></head><body>
{body}
<script id="wsspec" type="application/json">{spec_json}</script>
<script>{JS.replace("__SPEC__", "JSON.parse(document.getElementById('wsspec').textContent)")}</script>
{extra}
</body></html>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("-o", "--out")
    a = ap.parse_args()
    with open(a.spec, encoding="utf-8") as fh:
        spec = json.load(fh)
    try:
        doc = build(spec)
    except Refuse as e:
        sys.stderr.write(f"[REFUSE §13-A] {e}\n")
        sys.exit(EXIT_REFUSE)
    out = a.out or os.path.splitext(a.spec)[0] + ".html"
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print(f"[ok] wrote {out} ({os.path.getsize(out)//1024} KB, {len(spec['items'])} items)")


if __name__ == "__main__":
    main()
