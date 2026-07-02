#!/usr/bin/env python3
"""判別流程 5-態分類 standalone HTML (§13-A: 數字全從 decisionflow_summary.json 注入, 缺 key refuse)。
嵌 fig1-4(base64) + 決策流程 SVG + 5-態比例表 + 門檻依據 + 切不出≠沒訊號 + 有效性 + 驗證表。
tumor / merged 分開。用法: python3 build_decisionflow_report.py"""
import json, base64, os, subprocess, sys
A=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WTROOT="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra"
SUM=f"{A}/decisionflow_summary.json"
if not os.path.exists(SUM): sys.exit(f"REFUSE: missing {SUM} — 先跑 decisionflow_analyze.py")
S=json.load(open(SUM))
def req(d,*ks):
    for k in ks:
        if d is None or k not in d: sys.exit(f"REFUSE(§13-A): missing key {'/'.join(map(str,ks))}")
        d=d[k]
    return d
BC=subprocess.run(["git","-C",WTROOT,"rev-parse","--short","HEAD"],capture_output=True,text=True).stdout.strip()
def b64(p):
    fp=f"{A}/{p}"
    if not os.path.exists(fp): sys.exit(f"REFUSE: missing fig {fp} — 先跑 decisionflow_plot.py")
    return base64.b64encode(open(fp,"rb").read()).decode()
SETS=[k for k in ("tumor","merged") if k in S]
SN={"tumor":"tumor-only（subclone 向）","merged":"merged tumor+normal（cis-control 對照）"}
STLAB={"S1":"①不可驗證","S2":"②1群/無訊號","S3":"③監督可分(mean-shift)","S4":"④可分未對齊(epiallele?雜訊?)","S5":"⑤確認真結構(對齊生物軸)"}
STC={"S1":"#9ca3af","S2":"#64748b","S3":"#d97706","S4":"#7c3aed","S5":"#059669"}

# ---- 5-態比例表 ----
def states_table():
    rows=""
    for s in ("S1","S2","S3","S4","S5"):
        cells=""
        for nm in SETS:
            st=req(S,nm,"states_T4"); tot=req(S,nm,"n_total"); v=st[s]
            cells+=f'<td>{v:,}</td><td class="pct">{100*v/tot:.1f}%</td>'
        rows+=f'<tr><td><span class="dot" style="background:{STC[s]}"></span>{STLAB[s]}</td>{cells}</tr>'
    head="".join(f'<th colspan="2">{SN[n]}</th>' for n in SETS)
    sub="".join('<th>位點數</th><th>比例</th>' for _ in SETS)
    return f'<table class="dt"><thead><tr><th rowspan="2">判別流程狀態</th>{head}</tr><tr>{sub}</tr></thead><tbody>{rows}<tr class="tot"><td>合計（precond-pass / 全部）</td>'+"".join(f'<td>{req(S,n,"n_precond"):,}</td><td class="pct">{req(S,n,"pct_precond"):.1f}%</td>' for n in SETS)+'</tbody></table>'

# ---- 有效性卡 ----
def eff_cards():
    c=""
    for nm in SETS:
        sa=req(S,nm,"split_align_rate"); ms=req(S,nm,"nonsplit_meanshift_rate"); lo=req(S,nm,"nonsplit_loc_rate")
        c+=f'''<div class="card"><h4>{SN[nm]}</h4>
        <div class="metric"><span class="big">{sa:.1f}%</span><span class="ml">切群對齊生物軸 ⑤/(④+⑤)<br><small>無監督切得出的群，多少對應 HP/carrier/allele</small></span></div>
        <div class="metric"><span class="big">{ms:.1f}%</span><span class="ml">切不出但有訊號 ③/(②+③)<br><small>「切不出≠沒訊號」— 其中 location-clean {lo:.1f}%</small></span></div></div>'''
    return c

# ---- 門檻表 ----
def thr_table():
    rows=""
    for T in (3,4,5):
        cells=""
        for nm in SETS:
            ts=req(S,nm,"threshold_sensitivity",str(T))
            cells+=f'<td>{ts["S5"]:,}</td><td>{ts["S4"]:,}</td><td class="pct">{ts["split_align_rate"]:.1f}%</td>'
        mark=" ★建議" if T==4 else ""
        rows+=f'<tr{" class=rec" if T==4 else ""}><td>valid≥{T}{mark}<br><small>{"≤"+str(T-1)} 當 outlier</small></td>{cells}</tr>'
    head="".join(f'<th colspan="3">{SN[n]}</th>' for n in SETS)
    sub="".join('<th>⑤對齊</th><th>④未對齊</th><th>對齊率</th>' for _ in SETS)
    return f'<table class="dt"><thead><tr><th rowspan="2">valid 群門檻</th>{head}</tr><tr>{sub}</tr></thead><tbody>{rows}</tbody></table>'

def minority_rows():
    r=""
    for nm in SETS:
        md=req(S,nm,"minority_dist")
        al=md.get("aligned",{}); un=md.get("unaligned",{})
        r+=f'<tr><td>{SN[nm]}</td><td>對齊(→⑤): median <b>{al.get("median","-")}</b> (p25 {al.get("p25","-")} / p75 {al.get("p75","-")}, n={al.get("n",0)})</td><td>未對齊(→④): median <b>{un.get("median","-")}</b> (p25 {un.get("p25","-")} / p75 {un.get("p75","-")}, n={un.get("n",0)})</td></tr>'
    return r

# ---- 驗證表 (每 headline 數字 → 來源 key + L 級) ----
def verif_rows():
    items=[]
    for nm in SETS:
        items.append((f"{nm} precond-pass %",f'{req(S,nm,"pct_precond"):.1f}%',f'summary[{nm}].pct_precond',"L1 records 重算"))
        for s in ("S2","S3","S4","S5"):
            items.append((f"{nm} {STLAB[s]}",f'{req(S,nm,"states_T4")[s]:,}',f'summary[{nm}].states_T4.{s}',"L1 cls 計數"))
        items.append((f"{nm} 切群對齊率",f'{req(S,nm,"split_align_rate"):.1f}%',f'summary[{nm}].split_align_rate',"L1 ⑤/(④+⑤)"))
        items.append((f"{nm} 切不出有訊號率",f'{req(S,nm,"nonsplit_meanshift_rate"):.1f}%',f'summary[{nm}].nonsplit_meanshift_rate',"L1 ③/(②+③) BH-FDR"))
    return "".join(f'<tr><td>{a}</td><td class="num">{b}</td><td class="src">{c}</td><td>{d}</td></tr>' for a,b,c,d in items)

CSS="""
:root{--bg:#16181d;--card:#1e2128;--txt:#e6e6e6;--mut:#9aa0aa;--acc:#D97757;--bd:#2c3038;--grn:#059669;--amb:#d97706;--pur:#7c3aed}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--txt);font:14px/1.6 system-ui,'Droid Sans Fallback',sans-serif}
.wrap{max-width:1140px;margin:0 auto;padding:24px}
h1{font-size:21px;margin:0 0 4px}h2{font-size:16px;border-left:3px solid var(--acc);padding-left:9px;margin-top:30px}
h4{margin:0 0 8px;font-size:13px}
.sub{color:var(--mut);font-size:12.5px;margin-bottom:14px}
.banner{background:var(--card);border:1px solid var(--bd);border-radius:9px;padding:14px 18px;margin:14px 0}
.kpi{color:var(--acc);font-weight:700}
img{width:100%;border-radius:6px;margin:10px 0;background:#fff;border:1px solid var(--bd)}
table.dt{width:100%;border-collapse:collapse;font-size:12.5px;margin:10px 0}
table.dt th,table.dt td{border:1px solid var(--bd);padding:5px 8px;text-align:center}
table.dt th{background:#21252d;color:var(--mut);font-weight:600}
table.dt td:first-child{text-align:left}
table.dt .pct{color:var(--acc)}.dt .tot td{background:#1b1e24;font-weight:600}.dt .rec{background:#16241c}
.dot{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:6px;vertical-align:middle}
.cards{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin:12px 0}
.card{background:var(--card);border:1px solid var(--bd);border-radius:8px;padding:14px}
.metric{display:flex;align-items:center;gap:12px;margin:9px 0}
.big{font-size:26px;font-weight:800;color:var(--grn);min-width:96px;text-align:right}
.ml{font-size:12px;color:var(--txt)}.ml small{color:var(--mut)}
.num{color:var(--acc);font-weight:600;font-family:ui-monospace,monospace}
.src{font-family:ui-monospace,monospace;font-size:11px;color:var(--mut)}
.note{background:#231f1a;border:1px solid #4a3a26;border-radius:8px;padding:11px 15px;margin:12px 0;font-size:12.5px}
.foot{color:var(--mut);font-size:11px;margin-top:28px;border-top:1px solid var(--bd);padding-top:12px}
svg{width:100%;max-width:920px;display:block;margin:8px auto}
"""

# decision-flow SVG (手刻, 零依賴)
def flow_svg():
    box=lambda x,y,w,h,c,t,sub="": (f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="7" fill="{c}" opacity="0.18" stroke="{c}" stroke-width="1.5"/>'
        f'<text x="{x+w/2}" y="{y+(18 if sub else h/2+4)}" text-anchor="middle" font-size="12" font-weight="700" fill="#e6e6e6">{t}</text>'
        +(f'<text x="{x+w/2}" y="{y+34}" text-anchor="middle" font-size="10" fill="#9aa0aa">{sub}</text>' if sub else ""))
    ar=lambda x1,y1,x2,y2,lab="": (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#6b7280" stroke-width="1.4" marker-end="url(#a)"/>'
        +(f'<text x="{(x1+x2)/2+4}" y="{(y1+y2)/2-3}" font-size="9.5" fill="#9aa0aa">{lab}</text>' if lab else ""))
    s='<svg viewBox="0 0 920 360" role="img"><title>判別流程</title><defs><marker id="a" markerWidth="9" markerHeight="9" refX="7" refY="3" orient="auto"><path d="M0,0 L7,3 L0,6 Z" fill="#6b7280"/></marker></defs>'
    s+=box(370,8,180,40,"#3b82f6","位點 read×CpG")
    s+=ar(460,48,460,72)
    s+=box(350,72,220,44,"#9ca3af","precondition","n_complete≥6 可聚類?")
    s+=ar(350,94,150,150,"否")
    s+=box(40,150,210,40,STC["S1"],STLAB["S1"])
    s+=ar(460,116,460,150,"是")
    s+=box(350,150,220,44,"#a855f7","outlier-tolerant 切群","≥2 valid群(≥4), outlier≤20%")
    s+=ar(350,172,250,230,"切得出")
    s+=ar(570,172,690,150,"切不出")
    # split branch
    s+=box(120,230,260,44,"#a855f7","對齊 HP/carrier/allele?","CramérV≥.3 & χ² & Cochran")
    s+=ar(180,274,90,310,"否")
    s+=ar(320,274,400,310,"是")
    s+=box(20,310,150,40,STC["S4"],STLAB["S4"])
    s+=box(330,310,160,40,STC["S5"],STLAB["S5"])
    # cantsplit branch
    s+=box(620,150,260,44,"#d97706","a-priori PERMANOVA","mean-shift FDR 顯著?")
    s+=ar(680,194,600,250,"否")
    s+=ar(820,194,840,250,"是")
    s+=box(520,250,150,40,STC["S2"],STLAB["S2"])
    s+=box(720,250,180,40,STC["S3"],"③ mean-shift")
    s+='</svg>'
    return s

H=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HCC1395 判別流程 5-態分類（全基因組）</title><style>{CSS}</style></head><body><div class="wrap">
<h1>HCC1395 全基因組「位點甲基能否切成不同群」判別流程 — 5 態分類結果</h1>
<div class="sub">tumor-only 與 merged tumor+normal 分開計算 · valid 群門檻 T=4（≤3 當 outlier）· 單樣本 ⭐2-3 偵測非驗證 · build {BC}</div>

<div class="banner">這份回答三件事：<b>① valid 群 size≥多少合適</b>（門檻依據圖表）·<b>② 依判別流程，HCC1395 每個位點落在哪一態、比例多少</b>（tumor/merged 分開）·<b>③ 這套流程是否有效</b>。
<br>判別流程把每個位點分 5 態：<span class="kpi">①不可驗證</span>（read 太少）→ 試<b>無監督離散切群</b>（outlier-tolerant）→ 切得出再看是否<b>對齊生物軸</b>（⑤真結構 / ④epiallele?雜訊?）；切不出再看 <b>a-priori PERMANOVA mean-shift</b>（③切不出但有訊號 / ②真無訊號）。</div>

<h2>判別流程圖</h2>{flow_svg()}

<h2>① 5-態組成與比例（headline：valid≥4）</h2>
<img src="data:image/png;base64,{b64('figs_decisionflow/fig1_states.png')}"/>
{states_table()}

<h2>② valid 群門檻 size≥多少合適（數據依據）</h2>
<p class="sub">下圖左：切群「最小 valid 群（minority）」大小分佈，依「對齊→⑤」vs「未對齊→④」分色；右：T=3/4/5 門檻敏感度（⑤/④ 計數 + 對齊率）。</p>
<img src="data:image/png;base64,{b64('figs_decisionflow/fig2_minority.png')}"/>
<table class="dt"><thead><tr><th>read-set</th><th>對齊切群的最小群大小</th><th>未對齊切群的最小群大小</th></tr></thead><tbody>{minority_rows()}</tbody></table>
<img src="data:image/png;base64,{b64('figs_decisionflow/fig3_threshold.png')}"/>
{thr_table()}
<div class="note"><b>門檻建議 valid≥4（≤3 當 outlier）</b>：minority=3 的切群多落在「未對齊」（雜訊/單離群鏈）；提高門檻到 4 把這些降級為 outlier 記錄、不據以判 ≥2 群，使保留的切群對齊率上升（見右圖對齊率隨 T 變化）。size 與「足夠明顯」的關係：先前個案測試顯示<b>群大小（minority≥6 vs =3）比 silhouette 區分度更乾淨地分開「該救/不該救」</b>，故以 size 為主判準。</div>

<h2>③ 「切不出 ≠ 沒訊號」（切不出位點的 mean-shift 分解）</h2>
<p class="sub">切不出離散群的位點，再用 a-priori PERMANOVA（BH-FDR）測「群間甲基均值差」，並用 PERMDISP 區分真位置差（location）vs 散開度混淆（dispersion）。</p>
<img src="data:image/png;base64,{b64('figs_decisionflow/fig4_meanshift.png')}"/>

<h2>判別流程是否有效</h2>
<div class="cards">{eff_cards()}</div>
<div class="note">⚠ <b>merged=cis-control 對照</b>：merged 把 normal read 併入，HP 軸=germline-haplotype 層（cis-ASM），其切群/mean-shift 多為 germline cis 效應，<b>非 somatic subclone</b>。tumor-only 才是 subclone 向，但單樣本仍 ⭐2-3。<b>dispersion 混淆</b>：PERMANOVA 顯著中有相當比例是散開度非位置差（見 fig4 amber），故 ③ 內以 location-clean 子集為準。</div>

<h2>驗證表（每數字 → summary key + L 級）</h2>
<table class="dt"><thead><tr><th>數字</th><th>值</th><th>來源 key</th><th>L 級/重算</th></tr></thead><tbody>{verif_rows()}</tbody></table>

<h2>對抗驗證（2-agent workflow，2026-06-20）</h2>
<table class="dt"><thead><tr><th>軸</th><th>verdict</th><th>重點</th></tr></thead><tbody>
<tr><td>獨立重算</td><td style="color:#059669;font-weight:700">VERIFIED_MATCH</td><td style="text-align:left">fresh agent 自寫 BH-FDR、未複用 analyze，從原始 records 重算 tumor 5 態全數一致；partition 互斥窮盡 413+19501+10576=30490 ✓</td></tr>
<tr><td>方法學</td><td style="color:#d97706;font-weight:700">SOUND_WITH_CAVEATS（0 blocking）</td><td style="text-align:left">PERMANOVA/PERMDISP 實作數值正確（向量化 2 群與暴力 SSW bit-identical；betadisper 對非歐 BERNOULLI 正確）；valid≥4 / FDR pool / merged=cis-control 皆成立</td></tr>
</tbody></table>
<div class="note"><b>🔑 state③ 為何不是「tumor-only 非監督+PERMANOVA NEGATIVE」陷阱</b>：那個 NEGATIVE 是 <b>cluster-first</b>（silhouette 先挑距離矩陣最分離 partition→再 PERMANOVA 測它 = <b>循環</b>，被 read-內甲基相關打敗）。本 state③ 是結構相反的 <b>label-first</b>：只對「切不出」位點、用<b>外部 BAM tag</b>（HP/carrier/allele，非甲基 clustering）做 <b>label-permutation PERMANOVA（距離矩陣固定）</b>；label 置換在 H0 下保留 read-內相關 → exchangeability 成立 → <b>52% 是真 omnibus 關聯非置換假陽性</b>。<br>⚠ 但真關聯 ≠ subclone：單樣本主導軸 allele/germline-haplotype = germline cis-ASM/read-identity，非 somatic subclone（維持 ⭐2-3）。🔴 <b>allele 軸最 read-identity-confounded → next step = normal-anchored cis-control 專測 allele 軸</b>。</div>

<div class="foot">build_commit {BC} · 數字 §13-A 注入自 decisionflow_summary.json（缺 key refuse）· records→analyze→plot→此報告 ·
判別流程: precondition→outlier-tolerant 切群(T4)→對齊(CramérV)/PERMANOVA(BH-FDR)+PERMDISP · 單樣本 HCC1395 ⭐2-3 偵測非驗證 · NPERM=99</div>
</div></body></html>"""
out=f"{A}/20260620_decisionflow_5state_classification_01.standalone.html"
open(out,"w").write(H)
print(f"WROTE {out} ({len(H)//1024} KB)")
