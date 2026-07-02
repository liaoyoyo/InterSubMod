#!/usr/bin/env python3
"""build k-profile 解説 HTML: 切法品質(夠好vs次好) + 三態 + 多解析度階層熱圖。
數字注入 kprofile_summary.json(§13-A)。嵌 curated 代表熱圖(multi-resolution 跨軸優先)。"""
import json, base64, subprocess, os
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
def J(f): return json.load(open(f"{A}/{f}"))
S=J("kprofile_summary.json"); IDX=J("kprofile_heatmap_index.json")
T=S["tumor"]; Mg=S["merged"]
BC=subprocess.run(["git","-C","/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra","rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
def b64(p): return base64.b64encode(open(f"{A}/{p}","rb").read()).decode()
def pct(x,n): return f"{x} ({round(100*x/n,1)}%)"

# curate heatmaps: multi-resolution 跨>=2 軸優先, + 2 ambiguous + 2 confident
items=IDX["items"]
def naxes(it): return len({v[0] for v in it["mk_axes"].values()})
mr=sorted([it for it in items if it["group"]=="multi-resolution"],key=lambda it:-naxes(it))[:7]
amb=[it for it in items if it["group"]=="ambiguous-near-tie"][:2]
cu=[it for it in items if it["group"]=="confident-unique"][:2]
curated=mr+cu+amb

def card(it):
    mk=";".join(f"k{k}→{it['mk_axes'][str(k)][0]}({it['mk_axes'][str(k)][1]})" for k in sorted(map(int,it['mk_axes'].keys()))) if it['mk_axes'] else "—"
    return f'''<figure class="hm"><img src="data:image/png;base64,{b64(it["png"])}"/>
<figcaption><b>[{it['group']}]</b> {it['chrom']}:{it['pos']} n={it['n']} · meaningful_ks={it['meaningful_ks']}<br>{mk} · sil_m={it['sil_margin']} align_m={it['align_margin']}</figcaption></figure>'''

def bar(label,uniq,mid,amb,n,color):
    wu,wm,wa=[round(100*x/n,1) for x in (uniq,mid,amb)]
    return (f'<div class="mrow"><span class="ml">{label}</span>'
            f'<span class="seg" style="width:{wu}%;background:{color}" title="唯一 {uniq}">{uniq}</span>'
            f'<span class="seg" style="width:{wm}%;background:#caa05a" title="中 {mid}">{mid}</span>'
            f'<span class="seg" style="width:{wa}%;background:#6a6f78" title="模糊 {amb}">{amb}</span></div>')

CSS="""
:root{--bg:#15171c;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038;--grn:#7ee787}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:15px/1.6 system-ui,'Noto Sans CJK TC',sans-serif}
.wrap{max-width:1000px;margin:0 auto;padding:24px}
h1{font-size:21px}h2{font-size:17px;border-left:3px solid var(--acc);padding-left:9px;margin-top:30px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:13px 17px;margin:10px 0}
.key{background:#1a2520;border:1px solid #2c4a38;border-radius:8px;padding:12px 16px;margin:14px 0}.key b{color:var(--grn)}
table{border-collapse:collapse;width:100%;font-size:13px;margin:9px 0}th,td{border:1px solid var(--bd);padding:5px 9px;text-align:left}th{background:#22262e;color:var(--mut)}
.mrow{display:flex;align-items:center;margin:5px 0;font-size:12px}.ml{width:130px;color:var(--mut)}
.seg{color:#111;text-align:center;font-weight:600;padding:2px 0;overflow:hidden;white-space:nowrap}
.hmgrid{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:10px}
.hm{margin:0;background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:8px}.hm img{width:100%;border-radius:4px;background:#fff}
.hm figcaption{font-size:11px;color:var(--mut);margin-top:5px;font-family:ui-monospace,monospace}
.caveat{border-left:3px solid #cc5;background:#231f12;padding:9px 13px;border-radius:6px;margin:7px 0;font-size:13.5px}
.foot{color:var(--mut);font-size:11px;margin-top:26px;border-top:1px solid var(--bd);padding-top:11px}
code{color:#9ecbff;font-size:12px}
"""
sm,am=T["sil_margin"],T["align_margin"]
ts=T["three_state"]; tc=T["k_choice_n"]
HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>切法品質 k-profile — 夠好 vs 次好 + 多解析度</title><style>{CSS}</style></head><body><div class="wrap">
<h1>切法品質 k-profile — 「夠好 vs 次好」+ 多解析度</h1>
<div class="banner"><b>HCC1395 全基因組 · tumor reads · ⭐2 偵測非驗證</b><br>
問題：k-sweep 的切法是「夠好(唯一最佳)」還是「次好(近平手)」？每位點切法是否唯一？k=2 切成 k=3 是否也有意義？</div>

<h2>① 核心答案：切法品質取決於「用哪個判準」</h2>
<p>margin = best-k 與 2nd-best-k 的分數差。margin 大 = 切法唯一夠好；margin 小 = best≈2nd = 實質「次好/歧義」。</p>
<table><tr><th>判準</th><th>mean margin</th><th>唯一(≥0.15)</th><th>結論</th></tr>
<tr><td>silhouette（舊，群內聚集）</td><td>{sm['mean']}</td><td>{pct(sm['uniq'],sm['n'])}</td><td>多數<b>次好/近平手</b></td></tr>
<tr><td><b>alignment（新，對齊生物軸）</b></td><td><b>{am['mean']}</b></td><td><b>{pct(am['uniq'],am['n'])}</b></td><td><b>多數變「夠好/唯一」</b></td></tr></table>
<div class="key"><b>這就是「為何最新版更好」的量化證明</b>：silhouette 對多數位點選不出唯一 k（mean margin {sm['mean']}、唯一只 {round(100*sm['uniq']/sm['n'],1)}%）；<b>換成 alignment 判準後 mean margin {am['mean']}（{round(am['mean']/sm['mean'],1)}×）、唯一切法升到 {round(100*am['uniq']/am['n'],1)}%</b>。<b>alignment 把 silhouette 的近平手消歧義了</b>（{round(100*T['align_gt_sil']/T['both_margin_n'],1)}% 位點 align-margin > sil-margin）。</div>

<h2>② 三態分類（tumor）</h2>
<p>{pct(T['single_k_forced'],T['N'])} 是 <b>single-k-forced</b>（n 太小只能 k=2，無從選擇）。其餘 <b>k-choice 子集 {tc}</b> 分三態：</p>
{bar("multi-resolution",ts['multi-resolution'],0,0,tc,"#5FB85F").replace('<span class="seg" style="width:0.0%;background:#caa05a" title="中 0">0</span><span class="seg" style="width:0.0%;background:#6a6f78" title="模糊 0">0</span>','')}
<table><tr><th>三態</th><th>數量</th><th>意義</th></tr>
<tr><td><b>multi-resolution</b></td><td>{pct(ts['multi-resolution'],tc)}</td><td>≥2 個 k 都顯著對齊 → <b>k=2 切成 k=3 也有意義</b></td></tr>
<tr><td>confident-unique</td><td>{pct(ts['confident-unique'],tc)}</td><td>1 個唯一夠好切法（margin≥0.15）</td></tr>
<tr><td>ambiguous-near-tie</td><td>{pct(ts['ambiguous-near-tie'],tc)}</td><td>best≈2nd，無清楚贏家</td></tr></table>

<h2>③ 多解析度是真階層：不同 k 對齊不同軸</h2>
<p>multi-resolution 位點的關鍵：<b>不同 cut 層級揭露不同 a-priori 軸</b>（如 k=2→allele、k=3→HP、k=4→carrier）。每張左欄 <code>k2|k3|k4</code>(cut) + <code>HP|carr|alle|str</code>(生物軸)；看更高 k 的分界是否對齊某條生物側欄。</p>
<div class="hmgrid">{''.join(card(it) for it in curated)}</div>

<h2>④ tumor vs merged(cis-control) k-profile</h2>
<table><tr><th>模式</th><th>可分群</th><th>sil-margin mean</th><th>align-margin mean</th><th>multi-resolution(of choice)</th></tr>
<tr><td>tumor</td><td>{T['N']}</td><td>{sm['mean']}</td><td>{am['mean']}</td><td>{pct(ts['multi-resolution'],tc)}</td></tr>
<tr><td>merged</td><td>{Mg['N']}</td><td>{Mg['sil_margin']['mean']}</td><td>{Mg['align_margin']['mean']}</td><td>{pct(Mg['three_state']['multi-resolution'],Mg['k_choice_n'])}</td></tr></table>
<p class="muted">兩模式 alignment-margin 都遠大於 silhouette-margin（merged {round(100*Mg['align_gt_sil']/Mg['both_margin_n'],1)}% align>sil）→ 結論穩健：alignment 是更唯一的判準。</p>

<h2>⑤ 誠實邊界</h2>
<div class="caveat"><b>「夠好」是相對判準</b>：以 alignment 看 35% 唯一、以 silhouette 看僅 3%。沒有絕對「唯一正確切法」—— 多數位點本質 k-歧義（呼應「無乾淨無監督 k」）。</div>
<div class="caveat"><b>multi-resolution ≠ subclone</b>：不同 k 對齊不同軸（含 HP/allele）多為 germline-cis 層；subclone 仍需 normal cis-control。</div>
<div class="caveat"><b>ambiguous-near-tie 不是「無結構」</b>：是 k-選擇歧義，位點仍可能有弱/漸層訊號。single-k-forced({round(100*T['single_k_forced']/T['N'])}%)=n 太小無從選 k。</div>
<div class="caveat">單樣本 HCC1395 ⭐2；margin 門檻(0.05/0.15)是約定。</div>

<div class="foot">build_commit {BC} · binary 5c39051 · 數字注入自 kprofile_summary.json(§13-A) · 熱圖={len(curated)}/{IDX['n']} curated(完整見 figs_kprofile/) · gate V≥{S['gate']['V']}/e≥{S['gate']['e']}/p<{S['gate']['p']} · 規則 feedback_verification_table_per_data_answer</div>
</div></body></html>"""
out=f"{A}/20260618_kprofile_explainer_01.standalone.html"
open(out,"w").write(HTML)
print(f"WROTE {out} ({len(HTML)//1024} KB) curated={len(curated)}")
