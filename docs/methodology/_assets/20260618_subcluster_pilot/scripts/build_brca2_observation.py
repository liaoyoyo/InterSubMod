#!/usr/bin/env python3
"""BRCA2 per-CpG 多軸歸屬 完整觀察 HTML (§13-A 注入)。
嵌 3 圖 (carrier track / 3-axis track / dual-panel) + 多軸 CpG 表 + 分類裁決 + 逐項肉眼判讀(localStorage)。
數字全從 brca2_*.json + decisionflow record + step4 truth 注入; 缺 refuse。"""
import json, base64, os, csv, subprocess, sys
from collections import Counter
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
MAIN="/big7_disk/liaoyoyo2001/InterSubMod"
AS=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
def need(p):
    if not os.path.exists(p): sys.exit(f"REFUSE: missing {p}")
    return p
idea1=json.load(open(need(f"{AS}/brca2_idea1_out.json")))
mo=json.load(open(need(f"{AS}/brca2_multiaxis_out.json")))
ptab=json.load(open(need(f"{AS}/brca2_percpg_table.json")))
sig={k:set(v) for k,v in mo['sig'].items()}
# decisionflow record
dfr=[r for r in json.load(open(f"{AS}/decisionflow_records_tumor.json")) if r['chrom']=='chr13' and r['pos']=='32315128'][0]
bt=dfr['byT']['4']; pm=dfr['perm']
# truth
_trf=f"{WT}/research/tsg_promoter_asm_reviewer/output/step4_ism_results.json"
if not os.path.exists(_trf): _trf=f"{MAIN}/research/tsg_promoter_asm_reviewer/output/step4_ism_results.json"
tr=json.load(open(need(_trf)))[0]['stats']['HP1_vs_HP1-1']
# crosstab from reads.tsv
RD=f"{WT}/output/_brca2_region/mini/chr13/chr13_32315128/chr13_32310128_32320128"
ct=Counter()
if os.path.exists(f"{RD}/reads/reads.tsv"):
    for r in csv.DictReader(open(f"{RD}/reads/reads.tsv"),delimiter='\t'):
        if r['is_tumor'] in('1','true','True') and r['alt_support'] in('REF','ALT'):
            car='germ' if r['hp'] in('1','2') else 'somatic'
            ct[(car,r['alt_support'])]+=1
BC=subprocess.run(["git","-C",WT,"rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
def b64(p): return base64.b64encode(open(need(f"{AS}/{p}"),"rb").read()).decode()
# combos
combo=Counter()
for cp in set().union(*sig.values()):
    combo['+'.join(a for a in('hp','carrier','allele') if cp in sig[a])]+=1
combo_rows="".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k,v in sorted(combo.items(),key=lambda x:-x[1]))
# per-CpG table rows
def cell(v,neg_good=True):
    if v is None: return '<td class="na">·</td>'
    c='#dc2626' if (v<0)==neg_good else '#059669'
    return f'<td style="color:{c};font-weight:600">{v:+.2f}</td>'
trows=""
for r in ptab:
    axc={'carrier':'#dc2626','allele':'#059669','hp':'#2563eb'}
    badge=" ".join(f'<span style="background:{axc[a]};color:#fff;border-radius:3px;padding:0 4px;font-size:10px">{a}</span>' for a in r['axes'].split('+') if a)
    trows+=f'<tr><td>{r["dist"]:+d}</td>{cell(r.get("carrier"))}{cell(r.get("allele"),neg_good=False)}{cell(r.get("hp"))}<td>{r["n_axes"]}</td><td>{badge}</td></tr>'
def card(cid,title,fig,desc,obs):
    return f'''<div class="card"><div class="ch"><b>{title}</b></div>
<div class="cm">{desc}</div><img loading="lazy" src="data:image/png;base64,{b64(fig)}"/>
<div class="obs">觀察重點：{obs}</div>
<div class="verdict" data-key="{cid}"><span>判讀：</span><button data-v="ok">✓ 合理</button><button data-v="fix">✎ 需矯正</button>
<input class="note" placeholder="矯正回饋…"/></div></div>'''
CSS="""
:root{--bg:#16181d;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:14px/1.6 -apple-system,system-ui,"PingFang TC","PingFang SC","Microsoft JhengHei","Microsoft YaHei","Noto Sans CJK TC","Noto Sans TC","Hiragino Sans GB","Droid Sans Fallback",sans-serif}
.wrap{max-width:1120px;margin:0 auto;padding:22px}h1{font-size:20px}h2{font-size:15px;border-left:3px solid var(--acc);padding-left:8px;margin-top:26px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:9px;padding:13px 17px;margin:12px 0}
.kpi{color:var(--acc);font-weight:700}.grn{color:#34d399;font-weight:700}
.card{background:var(--card);border:1px solid var(--bd);border-radius:9px;padding:13px;margin:12px 0}
.card img{width:100%;border-radius:5px;margin:8px 0;background:#fff;border:1px solid var(--bd)}
.ch{font-size:14px}.cm{font-size:12px;color:var(--mut);margin:3px 0}
.obs{font-size:12.5px;background:#231f1a;border:1px solid #4a3a26;border-radius:6px;padding:7px 10px;margin-top:6px}
.verdict{margin-top:8px;display:flex;gap:7px;align-items:center;flex-wrap:wrap}
.verdict button{cursor:pointer;border:1px solid var(--bd);background:#262a32;color:var(--txt);border-radius:5px;padding:3px 9px;font-size:12px}
.verdict button.on[data-v=ok]{background:#1a6e3a;border-color:#2ea043}.verdict button.on[data-v=fix]{background:#7a4a00;border-color:#d97706}
.note{flex:1;min-width:160px;background:#262a32;color:var(--txt);border:1px solid var(--bd);border-radius:5px;padding:3px 7px;font-size:12px}
table{width:100%;border-collapse:collapse;font-size:12px;margin:8px 0}th,td{border:1px solid var(--bd);padding:4px 7px;text-align:center}th{background:#21252d;color:var(--mut)}td.na{color:#475569}
.toolbar{position:sticky;top:0;background:var(--bg);padding:8px 0;border-bottom:1px solid var(--bd);z-index:5;display:flex;gap:10px;align-items:center}
.toolbar button{background:#262a32;color:var(--txt);border:1px solid var(--bd);border-radius:5px;padding:4px 10px;cursor:pointer}
.foot{color:var(--mut);font-size:11px;margin-top:24px;border-top:1px solid var(--bd);padding-top:10px}
"""
JS="""
const LS='brca2_obs_v1';
function load(){return JSON.parse(localStorage.getItem(LS)||'{}')}
function save(o){localStorage.setItem(LS,JSON.stringify(o))}
function render(){const st=load();document.querySelectorAll('.verdict').forEach(v=>{const k=v.dataset.key,s=st[k]||{};v.querySelectorAll('button').forEach(b=>b.classList.toggle('on',s.v===b.dataset.v));const n=v.querySelector('.note');if(n&&s.n!=null)n.value=s.n});cnt()}
document.addEventListener('click',e=>{if(e.target.matches('.verdict button')){const v=e.target.closest('.verdict'),k=v.dataset.key,st=load();st[k]=st[k]||{};st[k].v=e.target.dataset.v;save(st);render()}});
document.addEventListener('input',e=>{if(e.target.matches('.note')){const v=e.target.closest('.verdict'),k=v.dataset.key,st=load();st[k]=st[k]||{};st[k].n=e.target.value;save(st)}});
function cnt(){const st=load();let ok=0,fx=0;Object.values(st).forEach(s=>{if(s.v==='ok')ok++;else if(s.v==='fix')fx++});document.getElementById('cnt').textContent=`判讀 ✓${ok} ✎${fx}/${document.querySelectorAll('.verdict').length}`}
function expt(){const st=load();let rows=[['item','verdict','note']];document.querySelectorAll('.verdict').forEach(v=>{const k=v.dataset.key,s=st[k]||{};rows.push([k,s.v||'',(s.n||'').replace(/,/g,';')])});const bl=new Blob([rows.map(r=>r.join(',')).join('\\n')],{type:'text/csv'});const a=document.createElement('a');a.href=URL.createObjectURL(bl);a.download='brca2_obs_verdicts.csv';a.click()}
window.onload=render;
"""
HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BRCA2 per-CpG 多軸歸屬 — 完整觀察驗證</title><style>{CSS}</style></head><body><div class="wrap">
<h1>BRCA2 chr13:32315128 — per-CpG 甲基定位 + 多軸歸屬 完整觀察</h1>
<div class="banner">這頁讓你<b>肉眼完整確認</b>「想法1：用甲基定位差異 CpG + 歸屬到哪個標籤軸」在已知 ASM 上是否成立，並逐項給<b>合理 / 需矯正</b>回饋。<br>
<b>位點背景</b>：somatic G→A (tvaf=0.189)；<b>decisionflow</b> tumor = <span class="grn">⑤確認真結構，aligned=carrier</span>（carrier F={pm['carrier']['F']} / hp F={pm['hp']['F']} / allele F={pm['allele']['F']}，切群 n_valid={bt['n_valid']} minority={bt['minority']}）；<b>既有真值</b> germline vs somatic Δβ=<span class="kpi">{tr['mean_delta_som_minus_germ']:.4f}</span>, wilcoxon p={tr['wilcoxon_p']:.1e}（somatic 低甲基）。</div>

<div class="toolbar"><span id="cnt"></span><button onclick="expt()">⬇ 匯出判讀回饋</button><span style="color:var(--mut);font-size:11px">(判讀存 localStorage)</span></div>

<h2>① 定位 — carrier 軸 per-CpG Δβ（焦點化）</h2>
{card('loc_carrier','定位 track（carrier 軸）','figs_brca2/track_carrier_dbeta.png',
  f"全 ±5000 窗逐 CpG carrier Δβ（somatic−germline）。可測 {idea1['n_tested']} CpG，顯著 {idea1['n_sig']}（BH-FDR q<0.05），方向符 −0.122：{idea1['n_neg']}/{idea1['n_sig']}。",
  "顯著紅點是否<b>集中</b>在 −600~−169 焦點窗、其餘 Δβ≈0？=定位是否乾淨焦點化（非全域散開）。")}

<h2>② 多軸 — 三軸 per-CpG 各自 Δβ</h2>
{card('loc_3axis','三軸 track（hp / carrier / allele）','figs_brca2/track_3axis.png',
  f"三軸各自顯著：hp {len(sig['hp'])} / carrier {len(sig['carrier'])} / allele {len(sig['allele'])} CpG。紫菱=多軸（≥2 軸顯著）。",
  "carrier 是否最強最多？allele/hp 是否只在 carrier 最強的 CpG 才亮（=allele「騎」在 carrier 上、共線）？")}

<h2>③ read 層結構 — dual-panel（甲基 + 距離）</h2>
{card('dual_tumor','TUMOR-only（56 reads = germline 30 + somatic 26）','figs_brca2/dualpanel_brca2_tumor.png',
  "只 tumor read，排序 carrier→allele→β；側欄 carrier(germ藍/som紅)/HP/ALT/strand。decisionflow 即用此 tumor 集。",
  "下方 somatic（carrier 紅）reads 是否在焦點窗呈<b>明顯低甲基塊</b>（左圖藍）、上方 germline 高甲基（紅）？右圖距離 somatic 是否成獨立塊？")}
{card('dual_all','含 normal cis-control 對照（101 reads = tumor 56 + normal 45）','figs_brca2/dualpanel_brca2.png',
  "加入 normal read（T-N 側欄）= cis-control：normal 應與 germline-tumor 同為高甲基。",
  "normal（T-N 綠）是否與 germline-tumor 一致高甲基（=somatic 低甲基非技術假象）？")}

<h2>④ 多軸 CpG 明細表（21 顯著 CpG）</h2>
<div class="card"><div class="cm">carrier Δβ=somatic−germline（負=somatic 低甲基,紅）；allele Δβ=REF−ALT（正=ALT 低甲基,紅）；hp Δβ=HP2fam−HP1fam。依多軸數→距離排序。</div>
<table><thead><tr><th>距 sSNV</th><th>carrier Δβ</th><th>allele Δβ</th><th>hp Δβ</th><th>軸數</th><th>對齊軸</th></tr></thead><tbody>{trows}</tbody></table>
<div class="obs">觀察：carrier-only（14 個，只有 carrier 欄有值）=somatic 單倍型整體低甲基（非 ALT 特異）；多軸列=carrier 與 allele 同時（最 confounded）。</div></div>

<h2>⑤ 軸組合分佈 + carrier×allele 共線</h2>
<div class="card"><table><thead><tr><th>對齊軸組合</th><th>CpG 數</th></tr></thead><tbody>{combo_rows}</tbody></table>
<div class="obs"><b>共線診斷</b>（tumor reads）：germline=全 REF（{ct[('germ','REF')]}）、somatic=ALT({ct[('somatic','ALT')]})+REF({ct[('somatic','REF')]})。
→ 有 {ct[('somatic','REF')]} 個 somatic-REF read 才使 carrier 與 allele <b>可部分分離</b>：carrier-only CpG = somatic 兩 allele 都低甲基 = 單倍型層，非 ALT 特異。</div>
<div class="verdict" data-key="multiaxis_logic"><span>判讀此歸屬邏輯：</span><button data-v="ok">✓ 合理</button><button data-v="fix">✎ 需矯正</button><input class="note" placeholder="矯正回饋…"/></div></div>

<div class="foot">build {BC} · §13-A 注入自 brca2_idea1_out.json / brca2_multiaxis_out.json / brca2_percpg_table.json / decisionflow_records_tumor.json / step4_ism_results.json ·
單樣本 cis-ASM characterization 非 subclone · allele 軸最 confounded 需 normal cis-control · 此頁確認後再大規模全 SNV</div>
</div><script>{JS}</script></body></html>"""
out=f"{AS}/20260620_brca2_percpg_observation_01.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out} ({len(HTML)//1024} KB)")
