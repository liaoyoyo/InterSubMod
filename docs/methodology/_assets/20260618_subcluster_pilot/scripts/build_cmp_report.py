#!/usr/bin/env python3
"""ISM joint vs modkit-equiv per-CpG 全資料對照 HTML (§13-A 注入 percpg_compare_summary.json)。"""
import json, base64, os, subprocess, sys
WT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
AS=f"{WT}/docs/methodology/_assets/20260618_subcluster_pilot"
S=json.load(open(f"{AS}/percpg_compare_summary.json"))
BC=subprocess.run(["git","-C",WT,"rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
def b64(p):
    fp=f"{AS}/{p}";
    if not os.path.exists(fp): sys.exit(f"REFUSE missing {fp}")
    return base64.b64encode(open(fp,"rb").read()).decode()
N=S['n_merged']; o1=S['overall_K1']; o3=S['overall_K3']; v=S['venn_K1']; bs=S['by_state']
AX=['hp','carrier','allele']
def vrow(a):
    x=v[a]; return f"<tr><td>{a}</td><td>{x['n_test']:,}</td><td>{x['both']:,}</td><td class=g>{x['ism_only']:,} ({x['ism_only_pct']}%)</td><td class=o>{x['mod_only']:,} ({x['mod_only_pct']}%)</td><td>{x['neither']:,}</td></tr>"
def srow(s,lab):
    d=bs[s]; return f"<tr><td>{lab}</td><td>{d['n']:,}</td><td>{d['median_nsig']:.0f}</td><td class=hl>{d['percpg_zero']:,} ({d['percpg_zero_pct']}%)</td></tr>"
CSS="""
:root{--bg:#16181d;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:14px/1.6 -apple-system,system-ui,"PingFang TC","Microsoft JhengHei","Noto Sans CJK TC","Droid Sans Fallback",sans-serif}
.wrap{max-width:1080px;margin:0 auto;padding:22px}h1{font-size:20px}h2{font-size:15px;border-left:3px solid var(--acc);padding-left:8px;margin-top:26px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:9px;padding:13px 17px;margin:12px 0}
.kpi{color:var(--acc);font-weight:700}.g{color:#34d399;font-weight:600}.o{color:#fbbf24;font-weight:600}.hl{color:#fbbf24;font-weight:600}
img{width:100%;border-radius:6px;margin:10px 0;background:#fff;border:1px solid var(--bd)}
table{width:100%;border-collapse:collapse;font-size:12.5px;margin:8px 0}th,td{border:1px solid var(--bd);padding:5px 8px;text-align:center}th{background:#21252d;color:var(--mut)}td:first-child{text-align:left}
.note{background:#231f1a;border:1px solid #4a3a26;border-radius:8px;padding:11px 15px;margin:12px 0;font-size:12.5px}
.foot{color:var(--mut);font-size:11px;margin-top:24px;border-top:1px solid var(--bd);padding-top:10px}
"""
H=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ISM joint vs modkit-equiv per-CpG 全資料對照</title><style>{CSS}</style></head><body><div class="wrap">
<h1>ISM(read 層 joint/cluster)vs modkit-equiv(per-CpG marginal)— 全資料對照</h1>
<div class="banner">回答：<b>modkit 式逐-CpG 能涵蓋多少、ISM 的 read 層 joint 多抓多少？</b> N={N:,} somatic 位點（全 22 autosome，100% full（chr20-22 因 CPU 競爭較慢但已完成））。
modkit 用 <b>per-CpG-by-label equiv（Fisher，= modkit marginal 檢定本質）</b>；ISM 用 <b>joint PERMANOVA（read×read 距離）</b> + 無監督 cluster。</div>

<img src="data:image/png;base64,{b64('figs_cmp/cmp_venn_state.png')}"/>

<h2>① per-軸 4-cell Venn（K=1：≥1 顯著 CpG = modkit 偵測到）</h2>
<table><thead><tr><th>軸</th><th>可比 n</th><th>both（modkit也抓）</th><th>ISM-only（joint顯著,per-CpG貧乏）</th><th>modkit-only（per-CpG有,joint不顯著）</th><th>neither</th></tr></thead>
<tbody>{vrow('hp')}{vrow('carrier')}{vrow('allele')}</tbody></table>
<div class="note"><b>整體 any-軸</b>（K=1）：both <b>{o1['both']:,}</b> · <span class=g>ISM-only <b>{o1['ism_only']:,}（{o1['ism_only_pct']}%）</b></span> · <span class=o>modkit-only {o1['mod_only']:,}（{o1['mod_only_pct']}%）</span> · neither {o1['neither']:,}。
嚴格門檻 K=3：ISM-only {o3['ism_only']:,}（{o3['ism_only_pct']}%）/ modkit-only {o3['mod_only']:,}（{o3['mod_only_pct']}%）。</div>

<h2>② 各 ISM state 中「modkit per-CpG 找不到」（=0）比例</h2>
<table><thead><tr><th>ISM state</th><th>位點數</th><th>per-CpG 顯著數 中位</th><th>per-CpG=0（modkit 一無所獲）</th></tr></thead>
<tbody>{srow('S5','⑤對齊真結構')}{srow('S4','④可分未對齊')}{srow('S3','③切不出但 joint 顯著')}{srow('S2','②無訊號')}</tbody></table>

<h2>解讀（誠實）</h2>
<div class="note">
<b>1. modkit 式 per-CpG 涵蓋多數。</b> both = {o1['both']:,}（佔有訊號位點絕大部分）→ 強 ASM 焦點訊號 modkit 也抓得到（如 BRCA2）。⑤ 位點 per-CpG=0 僅 {bs['S5']['percpg_zero_pct']}%（中位 {bs['S5']['median_nsig']:.0f} 顯著 CpG）。<br>
<b>2. ISM 的 read 層 joint 淨多抓 {o1['ism_only_pct']}%（{o1['ism_only']:,} 位點）</b>，> modkit-only {o1['mod_only_pct']}%（{o1['mod_only']:,}）≈ <b>{o1['ism_only']/o1['mod_only']:.1f}×</b>。這些是「分散弱訊號跨多 CpG」— 單點各自不顯著、聯合距離抓得到。<br>
<b>3. 最清楚的 ISM 優勢證據：③（切不出但 joint 顯著）有 {bs['S3']['percpg_zero_pct']}% 位點 per-CpG 完全=0</b> → ISM 說「有結構」、modkit 逐點「一無所獲」。<br>
<b>4. ⚠ ④ 未對齊有 {bs['S4']['percpg_zero_pct']}% per-CpG=0</b> → 三分之一的「可分未對齊」切群無 per-CpG 支持 = 雜訊嫌疑（呼應 ④ 桶異質）。<br>
<b>5. 此 Venn 只比「supervised per-CpG vs joint」；ISM 的無監督 clustering（⑤/④ 把 read 切群）modkit 完全沒有對應功能 → 是另一層 ISM-only 能力,未計入此表。</b></div>

<div class="foot">build {BC} · §13-A 注入 percpg_compare_summary.json · per-CpG-by-label = modkit marginal equiv（Fisher）· joint = ISM PERMANOVA · 100% full (30490) · 單樣本 ⭐2-3 characterization</div>
</div></body></html>"""
out=f"{AS}/20260620_ism_vs_modkit_percpg_comparison_01.standalone.html"
open(out,"w").write(H); print(f"WROTE {out} ({len(H)//1024} KB)")
