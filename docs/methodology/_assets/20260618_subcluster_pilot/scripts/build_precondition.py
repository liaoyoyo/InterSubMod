#!/usr/bin/env python3
"""build 驗證前提覆蓋率 HTML: 標籤種類分佈 + 各軸可驗證覆蓋 bar。注入 precondition_coverage.json(§13-A)。"""
import json, subprocess
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
P=json.load(open(f"{A}/precondition_coverage.json"))
BC=subprocess.run(["git","-C","/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra","rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
N=P["N"]
def bar(label,val,pct,color,note=""):
    return (f'<div class="brow"><span class="bl">{label}</span>'
            f'<span class="bt"><span class="bf" style="width:{pct}%;background:{color}">{pct}%</span></span>'
            f'<span class="bv">{val}{note}</span></div>')
nd=P["n_distinct_labels"]; ndp=P["n_distinct_pct"]
dist_bars="".join(bar(f"{k} 種標籤"+(" (無法驗證)" if k=="1" else ""),nd[k],ndp[k],
    "#6a6f78" if k=="1" else ["#4C9EE0","#5FB85F","#D97757"][int(k)-2]) for k in ("1","2","3","4"))
ax=[("HP 軸 (HP1-fam vs HP2-fam, tumor)",P["HP_axis"]["both_ge3"],P["HP_axis"]["both_ge3_pct"],"#4C9EE0"," ·45.8%只單側HP測不了"),
    ("CARRIER 軸 (germline vs carrier, tumor) ★subclone",P["carrier_axis"]["both_ge3"],P["carrier_axis"]["both_ge3_pct"],"#5FB85F",""),
    ("ALLELE 軸 (REF vs ALT)",P["allele_axis"]["csv_valid"],P["allele_axis"]["csv_valid_pct"],"#B05FA8",""),
    ("HP-fine 全4組各≥3",P["hpfine"]["all4"],P["hpfine"]["all4_pct"],"#caa05a"," ·子標籤稀疏")]
axis_bars="".join(bar(l,v,p,c,n) for l,v,p,c,n in ax)
CSS="""
:root{--bg:#15171c;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038;--grn:#7ee787}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:15px/1.6 system-ui,'Noto Sans CJK TC',sans-serif}
.wrap{max-width:920px;margin:0 auto;padding:24px}h1{font-size:20px}h2{font-size:16px;border-left:3px solid var(--acc);padding-left:9px;margin-top:26px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:13px 17px;margin:10px 0}
.key{background:#1a2520;border:1px solid #2c4a38;border-radius:8px;padding:12px 16px;margin:12px 0}.key b{color:var(--grn)}
.warn{background:#231f12;border:1px solid #4a432c;border-left:3px solid #cc5;border-radius:6px;padding:10px 14px;margin:10px 0}
.brow{display:flex;align-items:center;margin:6px 0;font-size:13px}.bl{width:330px;color:var(--mut)}
.bt{flex:1;background:#262a32;border-radius:4px;height:22px;overflow:hidden}.bf{display:block;height:22px;line-height:22px;color:#111;font-weight:700;text-align:right;padding-right:6px;font-size:11px;min-width:30px}
.bv{width:150px;text-align:right;font-family:ui-monospace,monospace;font-size:12px;color:#cdd}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}th,td{border:1px solid var(--bd);padding:5px 9px;text-align:left}th{background:#22262e;color:var(--mut)}
.foot{color:var(--mut);font-size:11px;margin-top:24px;border-top:1px solid var(--bd);padding-top:10px}code{color:#9ecbff;font-size:12px}
"""
HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>驗證前提覆蓋率</title><style>{CSS}</style></head><body><div class="wrap">
<h1>驗證前提覆蓋率 — 多少位點有 ≥2 種標籤可驗證</h1>
<div class="banner"><b>HCC1395 全基因組 N={N} · tumor 標籤 · ⭐2</b><br>
標籤檢定(Fisher+V / PERMANOVA / Δβ)要能跑，位點必須有 <b>≥2 種不同標籤(各≥3 reads)</b>。本頁量化前提覆蓋。</div>

<h2>① 標籤種類分佈（會不會很多全是 1 種？）</h2>
{dist_bars}
<div class="key"><b>不會 —— 單一標籤只 {P['single_pct']}%</b>（{P['single_label']} 個）。<b>{P['multi_pct']}% 有 ≥2 種標籤</b>、前提多數成立。</div>

<h2>② 各軸可驗證覆蓋（both ≥3 reads）</h2>
{axis_bars}
<table><tr><th>軸</th><th>可驗證</th><th>佔比</th><th>意義</th></tr>
<tr><td>HP 軸 (tumor)</td><td>{P['HP_axis']['both_ge3']}</td><td>{P['HP_axis']['both_ge3_pct']}%</td><td>45.8% 只單側 HP(LOH/主導)→測不了</td></tr>
<tr><td><b>CARRIER 軸 (tumor)★</b></td><td>{P['carrier_axis']['both_ge3']}</td><td><b>{P['carrier_axis']['both_ge3_pct']}%</b></td><td>subclone 前提(germline+carrier 都在)</td></tr>
<tr><td>ALLELE 軸</td><td>{P['allele_axis']['csv_valid']}</td><td>{P['allele_axis']['csv_valid_pct']}%</td><td>somatic 位點幾乎都兩等位</td></tr>
<tr><td>HP-fine 全4組</td><td>{P['hpfine']['all4']}</td><td>{P['hpfine']['all4_pct']}%</td><td>子標籤稀疏(≥3組僅 {P['hpfine']['ge3_pct']}%)</td></tr></table>

<h2>🔴 ③ 重要釐清：tumor-only vs paired</h2>
<div class="warn">C++ 的 <code>LabelHPPermanovaValid</code> = <b>{P['HP_axis']['csv_valid_pct']}%</b>，但 tumor-only 算 HP 軸只 <b>{P['HP_axis']['both_ge3_pct']}%</b>。差別在 <b>C++ PERMANOVA 是 paired(含 normal reads，normal 補滿 HP1/HP2)</b>。<br>
→ <b>「85% PERMANOVA sig」的分母是 paired 全集；只看 tumor-only，HP 軸前提只 {P['HP_axis']['both_ge3_pct']}%。</b> subclone(CARRIER 軸 tumor-only)前提 {P['carrier_axis']['both_ge3_pct']}%。</div>

<h2>④ 整合觀察</h2>
<div class="key">
• 前提覆蓋**不是問題**：90.8% 有 ≥2 標籤、CARRIER(subclone)軸 79.5% tumor-only 可測。<br>
• **HP 軸 tumor-only 只 54%**(45.8% LOH/單側)→ 報 HP-ASM 時分母要講清楚。<br>
• **4 組 HP-fine 全有的只 7.5%**→ 多組細分受限子標籤稀疏(呼應前面 PERMANOVA 4 組瓶頸)。<br>
• ALLELE 軸幾乎全覆蓋(99.9%)→ 等位軸前提最寬。</div>

<div class="foot">build_commit {BC} · 數字注入自 precondition_coverage.json(§13-A) · tumor 標籤計數自 records_wg2.json·all.labels / allele valid 自 significance_csv(paired) · 單樣本 ⭐2</div>
</div></body></html>"""
out=f"{A}/20260620_precondition_coverage_01.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out} ({len(HTML)//1024} KB)")
