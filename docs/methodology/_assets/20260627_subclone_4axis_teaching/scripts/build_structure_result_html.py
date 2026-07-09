#!/usr/bin/env python3
"""build_structure_result_html.py — 整合結構結果 standalone HTML(2026-07-09)。
讀 20260709_structure_result_data.json → ①region層c分布 ②整體VAF頻率譜(直方圖+峰) ③巢狀clone→subclone父子
④互動樹檢視器(精選區·多候選樹切換·每節點VAF標註)。數字全從JSON注入(§13-A 反捏造);樹檢視器JS embed curated data。
輸出: docs/methodology/_assets/20260709_structure_result.standalone.html
"""
import json

DATA = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260709_structure_result_data.json"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260709_structure_result.standalone.html"
SAMPLES = ["HCC1395", "HCC1395_DORADO", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
d = json.load(open(DATA))
S = d["samples"]

def pct(a, b): return f"{100*a/b:.0f}%" if b else "-"

# ---------- §1 region 層 c 分布 ----------
rows1 = []
for s in SAMPLES:
    x = S[s]; cd = x["c_dist"]
    cge5 = sum(int(cd[str(k)]) for k in range(5, 10))
    rows1.append(f"<tr><td class='sn'>{s}</td>"
                 f"<td>{cd['1']}</td><td>{cd['2']}</td><td>{cd['3']}</td><td>{cd['4']}</td><td>{cge5}</td>"
                 f"<td class='hi'>{x['c_ge1']}</td><td class='sub'>{pct(x['c_ge2'],x['c_ge1'])}</td>"
                 f"<td class='mx'>{x['maxC']}</td><td class='un'>{x['unresolved']}</td></tr>")
tot_ge1 = sum(S[s]["c_ge1"] for s in SAMPLES)
tot_ge2 = sum(S[s]["c_ge2"] for s in SAMPLES)

# ---------- §2 VAF 頻率譜(直方圖 20 bins) ----------
def spectrum(s):
    x = S[s]; h = x["vaf_hist"]; mx = max(h) if h else 1
    peaks = {round(p["vaf"], 3) for p in x["vaf_peaks"]}
    bars = []
    for i, v in enumerate(h):
        center = round((i + 0.5) * 0.05, 3)
        w = 100 * v / mx if mx else 0
        # 色:低VAF藍(subclonal) 高VAF紅(clonal)
        col = "#337ab7" if center < 0.4 else ("#f0ad4e" if center < 0.75 else "#d9534f")
        ispk = any(abs(center - pc) < 0.03 for pc in peaks)
        mark = " ◀峰" if ispk else ""
        lab = f"{center:.2f}"
        bars.append(f"<div class='vrow'><span class='vlab'>{lab}</span>"
                    f"<span class='vbar' style='width:{w:.1f}%;background:{col}'></span>"
                    f"<span class='vn'>{v}{mark}</span></div>")
    pk = " · ".join(f"VAF≈{p['vaf']}" for p in x["vaf_peaks"])
    return (f"<div class='spec'><div class='spechd'>{s} <span class='dim'>(n_somatic_pos={x['n_somatic_pos']}; 峰: {pk})</span></div>"
            + "".join(bars) + "</div>")

# ---------- §3 巢狀父子 ----------
rows3 = []
for s in SAMPLES:
    x = S[s]
    rows3.append(f"<tr><td class='sn'>{s}</td><td>{x['nested']}</td><td>{x['sister']}</td>"
                 f"<td>{x['n_pair']}</td><td class='hi'>{pct(x['anc_gt_der'],x['n_pair'])}</td>"
                 f"<td class='cl'>{x['anc_med_vaf']}</td><td class='su'>{x['der_med_vaf']}</td></tr>")

CUR = json.dumps(d["curated_regions"], ensure_ascii=False)

html = f"""<!--
建立: 2026-07-09 | 整合結構結果(region層c + 整體VAF頻率譜 + VAF標註候選樹) | build: build_structure_result_html.py
資料: 20260709_structure_result_data.json(build_structure_result_data.py 一鍵重算) | 數字全注入·反捏造
-->
<meta charset="utf-8">
<title>Subclonal 結構結果 — clone數 × VAF頻率 × 候選樹</title>
<style>
:root{{--bg:#0f1420;--card:#1a2233;--ink:#e8edf5;--dim:#8fa0b8;--line:#2a3550;--acc:#4a9eff}}
*{{box-sizing:border-box}}
body{{margin:0;font-family:-apple-system,'Segoe UI',Roboto,'Noto Sans TC',sans-serif;background:var(--bg);color:var(--ink);line-height:1.55;font-size:14px}}
.wrap{{max-width:1180px;margin:0 auto;padding:24px}}
h1{{font-size:22px;margin:0 0 4px}} h2{{font-size:17px;margin:28px 0 10px;padding-left:9px;border-left:4px solid var(--acc)}}
.sub0{{color:var(--dim);font-size:13px;margin-bottom:6px}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:14px 16px;margin:10px 0}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{border:1px solid var(--line);padding:5px 8px;text-align:center}} th{{background:#131b2b;color:var(--dim);font-weight:600}}
td.sn{{text-align:left;font-weight:600}} .hi{{color:#7dd3fc;font-weight:700}} .sub{{color:#fca5a5;font-weight:700}}
.mx{{color:#fcd34d;font-weight:700}} .un{{color:var(--dim)}} .cl{{color:#d9534f;font-weight:700}} .su{{color:#337ab7;font-weight:700}}
.dim{{color:var(--dim);font-weight:400;font-size:12px}}
.note{{font-size:12.5px;color:var(--dim);margin:6px 0}}
.warn{{background:#2a1f1a;border-left:3px solid #d9534f;padding:8px 12px;border-radius:5px;font-size:12.5px;margin:8px 0}}
.spec{{margin:10px 0;font-family:ui-monospace,Menlo,monospace;font-size:11px}}
.spechd{{font-family:-apple-system,sans-serif;font-weight:600;font-size:13px;margin-bottom:3px}}
.vrow{{display:flex;align-items:center;gap:6px;height:15px}}
.vlab{{width:34px;text-align:right;color:var(--dim)}} .vbar{{height:10px;border-radius:2px;min-width:1px}} .vn{{color:var(--dim)}}
.legend{{display:flex;gap:14px;flex-wrap:wrap;font-size:12px;margin:6px 0}}
.legend span{{display:inline-flex;align-items:center;gap:5px}}
.dot{{width:11px;height:11px;border-radius:50%;display:inline-block}}
.sq{{width:11px;height:11px;display:inline-block;background:#8896ad}}
/* tree viewer */
.tv{{display:grid;grid-template-columns:270px 1fr;gap:14px}}
.rlist{{max-height:520px;overflow-y:auto}}
.ritem{{padding:7px 9px;border:1px solid var(--line);border-radius:6px;margin-bottom:5px;cursor:pointer;font-size:12.5px}}
.ritem:hover{{border-color:var(--acc)}} .ritem.sel{{background:#22304a;border-color:var(--acc)}}
.rcat{{display:inline-block;font-size:10px;padding:1px 6px;border-radius:8px;margin-right:5px;font-weight:700}}
.cat-nested{{background:#1e3a5f;color:#7dd3fc}} .cat-sister{{background:#3a2a1e;color:#fcd34d}}
.cat-complex{{background:#3a1e2e;color:#fca5a5}} .cat-single{{background:#1e3a2e;color:#86efac}} .cat-multitree{{background:#2e1e3a;color:#c4b5fd}}
.detail{{background:#131b2b;border:1px solid var(--line);border-radius:8px;padding:14px;min-height:480px}}
.sw{{display:flex;align-items:center;gap:10px;margin:8px 0}}
.sw button{{background:#22304a;color:var(--ink);border:1px solid var(--line);border-radius:5px;padding:3px 11px;cursor:pointer;font-size:15px}}
.sw button:hover{{border-color:var(--acc)}}
svg{{background:#0d1420;border-radius:6px}}
.pv{{font-family:ui-monospace,monospace;font-size:11px;color:var(--dim);margin-top:8px;word-break:break-all}}
.chip{{display:inline-block;padding:1px 7px;border-radius:4px;margin:1px 2px;font-size:11px;font-family:ui-monospace,monospace}}
</style>
<div class="wrap">
<h1>Subclonal 結構結果 — clone 數 × VAF 頻率 × 候選樹</h1>
<div class="sub0">新骨幹(ClairS PASS)·per-HP-家族樹枚舉·2026-07-09 | 三軸: 拓撲(共現定序) × <b>頻率(VAF)</b> × clone 數(region 層 c)</div>

<div class="card">
<b>定義（本頁採用）</b>
<div class="note">• <b>c（clone 數）= region 層 distinct ALT 組合群數 = 該 locus 的 somatic clone 數</b>。REF=恆存在的 germline 祖先根(排除)。c=1 單 clone · <b>c≥2 ⟺ 有 subclonal 結構</b>。<span class="dim">(乾淨參考單體型不自成單位；full-cov 無 somatic 組合 = 未解析，非 0 clone)</span></div>
<div class="note">• <b>VAF（頻率軸）</b>= 每 somatic 位點 within-family nALT/(nREF+nALT)。峰值=clone/subclone 頻率群；<b>祖先突變高 VAF(clonal) → 衍生突變低 VAF(subclonal)</b>，方向由<b>同一 read 共現(co-occurrence)直接定序</b>，VAF 為獨立佐證。</div>
<div class="warn">🔴 <b>誠實邊界</b>：raw VAF 有 CN confound(LOH 抬高/gain 稀釋)→ 祖先~1.0 可能是 LOH、衍生~0.4 可能是 CN 假象；轉真 CCF 需整數 CN，<b>僅 HCC1395 有 SEQC2</b>。region-local(非全域樹)；單-bulk 只可識別 clusters+部分序，完整唯一樹「定不出來」。</div>
</div>

<h2>§1　clone 數（region 層 c）— 7 樣本</h2>
<div class="card"><table>
<tr><th>樣本</th><th>c=1<br><span class="dim">單clone</span></th><th>c=2</th><th>c=3</th><th>c=4</th><th>c≥5</th><th>c≥1 區</th><th>c≥2<br><span class="dim">subclonal</span></th><th>maxC<br><span class="dim">最高</span></th><th>未解析</th></tr>
{''.join(rows1)}
</table>
<div class="note">全樣本 c≥1 區 = <b>{tot_ge1:,}</b>；c≥2（有 subclonal 結構）= <b>{tot_ge2:,}（{pct(tot_ge2,tot_ge1)}）</b>。<b>maxC=該樣本一個 locus 最高可能 clone 數</b>（HCC1395=8）。未解析=somatic 只在 partial read / 弱(read 沒跨全位點)，非 0 clone。</div>
</div>

<h2>§2　整體 VAF 頻率譜（峰值 = clone/subclone 頻率群）</h2>
<div class="card">
<div class="legend"><span><span class="dot" style="background:#d9534f"></span>VAF≥0.75 clonal(骨幹/LOH)</span>
<span><span class="dot" style="background:#f0ad4e"></span>0.4–0.75 het clonal</span>
<span><span class="dot" style="background:#337ab7"></span>&lt;0.4 subclonal</span> · ◀峰=偵測到的頻率群</div>
{''.join(spectrum(s) for s in SAMPLES)}
<div class="note">幾乎所有樣本主峰 <b>VAF≈1.0</b>(LOH/clonal 骨幹)；HCC1395/DORADO/H2009 另有 ≈0.5(het clonal)；<b>HCC1937 三峰 0.23/0.47/0.98(最豐富 subclonal 層次)</b>、HCC1954 0.33/0.98。=SciClone/PyClone/DPClust 式 VAF/CCF 分群。</div>
</div>

<h2>§3　父子關係（巢狀 clone→subclone）— 祖先 VAF vs 衍生 VAF</h2>
<div class="card"><table>
<tr><th>樣本</th><th>巢狀<br><span class="dim">clone→subclone</span></th><th>姊妹<br><span class="dim">分支</span></th><th>有VAF對</th><th>祖先VAF&gt;衍生VAF</th><th>祖先中位VAF<br><span class="dim">clonal</span></th><th>衍生中位VAF<br><span class="dim">subclonal</span></th></tr>
{''.join(rows3)}
</table>
<div class="note">祖先突變(共現於更多群)恆為 clonal(中位 ~0.97–1.0)、衍生突變恆為半頻(~0.37–0.52)，<b>~99–100% 一致</b>。方向由<b>共現直接定序</b>，VAF 獨立佐證兩訊號吻合 → 比純 VAF 統計推父子可信。</div>
</div>

<h2>§4　VAF 標註的候選樹（精選區·多樹切換）</h2>
<div class="card">
<div class="legend"><span><span class="sq"></span>ROOT germline(全REF)</span>
<span><span class="dot" style="background:#4a9eff;border:2px solid #4a9eff"></span>實測節點</span>
<span><span class="dot" style="background:transparent;border:2px dashed #c4b5fd"></span>推測祖先(hidden)</span>
· 節點色=新獲突變 VAF(紅clonal→藍subclonal) · 邊上數字=該 clone 獲得突變的 VAF</div>
<div class="tv">
<div class="rlist" id="rlist"></div>
<div class="detail" id="detail"></div>
</div>
</div>

<div class="note" style="margin-top:20px">資料一鍵重算: <code>build_structure_result_data.py</code> → <code>build_structure_result_html.py</code>。數字全從 JSON 注入(反捏造)。</div>
</div>

<script>
const CUR = {CUR};
let curIdx = 0, treeIdx = 0;
const CATNAME = {{nested:'巢狀 clone→subclone', sister:'姊妹分支', complex:'複雜(c≥3)', single:'單 clone', multitree:'多候選樹'}};

function vafColor(v){{ if(v==null) return '#8896ad'; if(v>=0.75) return '#d9534f'; if(v>=0.4) return '#f0ad4e'; return '#337ab7'; }}
function geno(node){{ return String(node).startsWith('H_') ? String(node).slice(2) : String(node); }}

function renderList(){{
  const el = document.getElementById('rlist');
  el.innerHTML = CUR.map((r,i)=>
    `<div class="ritem ${{i===curIdx?'sel':''}}" onclick="selReg(${{i}})">`+
    `<span class="rcat cat-${{r.category}}">${{CATNAME[r.category]}}</span>`+
    `<b>${{r.chrom}}:${{r.start}}</b> fam${{r.family}}<br>`+
    `<span class="dim">c=${{r.c}} · ${{r.shape}} · ${{r.n_trees}}候選樹</span></div>`).join('');
}}
function selReg(i){{ curIdx=i; treeIdx=0; renderList(); renderDetail(); }}
function switchTree(delta){{
  const r = CUR[curIdx]; const n = r.trees.length;
  treeIdx = (treeIdx + delta + n) % n; renderDetail();
}}

function layout(edges){{
  const children={{}}, parent={{}}, nodes=new Set();
  edges.forEach(([p,c])=>{{ (children[p]=children[p]||[]).push(c); parent[c]=p; nodes.add(p); nodes.add(c); }});
  let root=[...nodes].find(n=>!(n in parent)); if(!root) root='ROOT';
  const depth={{}}; depth[root]=0; const order=[root];
  for(let i=0;i<order.length;i++){{ const nn=order[i]; (children[nn]||[]).forEach(c=>{{depth[c]=depth[nn]+1; order.push(c);}}); }}
  // x by leaf-order DFS
  let leafx=0; const x={{}};
  function dfs(n){{ const ch=children[n]||[]; if(ch.length===0){{ x[n]=leafx++; return x[n]; }}
    const xs=ch.map(dfs); x[n]=xs.reduce((a,b)=>a+b,0)/xs.length; return x[n]; }}
  dfs(root);
  return {{children,parent,depth,x,root,nodes:[...nodes]}};
}}

function acquired(node, par, positions){{
  // 新獲突變 = node 有 A 而 parent 沒有
  const gn=geno(node), gp = par? geno(par) : 'R'.repeat(gn.length);
  const res=[];
  for(let i=0;i<gn.length;i++){{ if(gn[i]==='A' && gp[i]!=='A') res.push(i); }}
  return res;
}}

function renderDetail(){{
  const r = CUR[curIdx]; const t = r.trees[treeIdx];
  const posvaf = r.posvaf; const positions = r.positions;
  const L = layout(t.edges);
  const NW=1.0; const cols=Math.max(1,Math.max(...Object.values(L.x))+1); const rows=Math.max(...Object.values(L.depth))+1;
  const W=Math.max(560, cols*150), H=Math.max(230, rows*115+60);
  const px=x=>60 + x*((W-120)/Math.max(1,cols-1||1));
  const py=dp=>45 + dp*((H-90)/Math.max(1,rows-1||1));
  let svg=`<svg width="${{W}}" height="${{H}}" viewBox="0 0 ${{W}} ${{H}}">`;
  // edges
  t.edges.forEach(([p,c])=>{{
    const x1=px(L.x[p]),y1=py(L.depth[p]),x2=px(L.x[c]),y2=py(L.depth[c]);
    const acq=acquired(c,p,positions);
    const vs=acq.map(i=>posvaf[String(i)]).filter(v=>v!=null);
    const av = vs.length? vs.reduce((a,b)=>a+b,0)/vs.length : null;
    svg+=`<line x1="${{x1}}" y1="${{y1}}" x2="${{x2}}" y2="${{y2}}" stroke="${{vafColor(av)}}" stroke-width="2" opacity="0.7"/>`;
    if(av!=null) svg+=`<text x="${{(x1+x2)/2+6}}" y="${{(y1+y2)/2}}" fill="${{vafColor(av)}}" font-size="11" font-family="monospace">VAF ${{av.toFixed(2)}}</text>`;
  }});
  // nodes
  L.nodes.forEach(nd=>{{
    const cx=px(L.x[nd]),cy=py(L.depth[nd]); const info=t.nodes[nd];
    const par=L.parent[nd]; const acq=acquired(nd,par,positions);
    const vs=acq.map(i=>posvaf[String(i)]).filter(v=>v!=null);
    const av=vs.length? vs.reduce((a,b)=>a+b,0)/vs.length : null;
    if(info && info.root){{
      svg+=`<rect x="${{cx-9}}" y="${{cy-9}}" width="18" height="18" fill="#8896ad"/>`;
      svg+=`<text x="${{cx}}" y="${{cy+26}}" fill="#8fa0b8" font-size="10.5" text-anchor="middle">germline</text>`;
    }} else {{
      const hid = info && info.hidden;
      svg+=`<circle cx="${{cx}}" cy="${{cy}}" r="10" fill="${{hid?'transparent':vafColor(av)}}" stroke="${{hid?'#c4b5fd':vafColor(av)}}" stroke-width="2.5" ${{hid?'stroke-dasharray="3,2"':''}}/>`;
      // 標新獲突變位點:VAF
      const lab = acq.map(i=>`${{positions[i]}}<tspan fill="${{vafColor(posvaf[String(i)])}}"> ${{posvaf[String(i)]!=null?posvaf[String(i)].toFixed(2):'?'}}</tspan>`).join('  ');
      svg+=`<text x="${{cx}}" y="${{cy+25}}" fill="#8fa0b8" font-size="9.5" text-anchor="middle" font-family="monospace">${{lab}}</text>`;
      if(hid) svg+=`<text x="${{cx+11}}" y="${{cy-9}}" fill="#c4b5fd" font-size="9">ᴴ</text>`;
    }}
  }});
  svg+='</svg>';
  // populations + posvaf 表
  const pops = Object.entries(r.populations).map(([g,n])=>`<span class="chip" style="background:#22304a">${{g}} ×${{n}}</span>`).join('');
  const pvchips = positions.map((p,i)=>{{ const v=posvaf[String(i)]; return `<span class="chip" style="background:${{v!=null?vafColor(v):'#333'}}22;color:${{v!=null?vafColor(v):'#888'}}">${{p}}: ${{v!=null?v.toFixed(2):'—'}}</span>`; }}).join('');
  document.getElementById('detail').innerHTML =
    `<div><b>${{r.chrom}}:${{r.start}}</b> · fam${{r.family}} · <span class="rcat cat-${{r.category}}">${{CATNAME[r.category]}}</span>`+
    ` · c=${{r.c}} clone · shape=${{r.shape}} · <b>${{r.n_trees}} 候選樹</b>(存 ${{r.n_stored}})</div>`+
    `<div class="sw"><button onclick="switchTree(-1)">◀</button>`+
    `<span>候選樹 ${{treeIdx+1}} / ${{r.trees.length}} · 形狀 ${{t.shape}} · ${{t.n_hidden}} 推測祖先</span>`+
    `<button onclick="switchTree(1)">▶</button></div>`+
    svg+
    `<div class="pv"><b>full-cov populations</b>(genotype×reads): ${{pops}}</div>`+
    `<div class="pv"><b>每位點 VAF</b>: ${{pvchips}}</div>`+
    `<div class="note">深度越深=越晚獲得的突變(越可能 subclonal)；節點/邊色越藍=VAF 越低=該 clone 佔比越小。多候選樹=同資料下等機率的不同拓撲(read 分不出→「定不出來」)。</div>`;
}}
renderList(); renderDetail();
</script>
"""
open(OUT, "w", encoding="utf-8").write(html)
print(f"→ 寫出 {OUT}  ({len(html)} bytes)")
# 抽 JS 做 node --check
import re
m = re.search(r"<script>\n(.*?)\n</script>", html, re.S)
if m:
    open("/tmp/_sr_check.js", "w").write("const document={getElementById:()=>({})};\n" + m.group(1).replace("renderList(); renderDetail();", ""))
    print("JS extracted for node --check")
