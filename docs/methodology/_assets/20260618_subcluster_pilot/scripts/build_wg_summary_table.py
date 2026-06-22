#!/usr/bin/env python3
"""全基因組分類觀察表(單頁 HTML)— 所有比例 + 甲基訊號層級 + TP/FP 特異性。§13-A 從 records 注入。"""
import json, os, html, base64
from collections import Counter
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
rec=[r for r in json.load(open(f"{A}/cluster_redesign_wg_records.json")) if "err" not in r]
TP=[r for r in rec if r["set"]=="TP"]; FP=[r for r in rec if r["set"]=="FP"]
nTP,nFP=len(TP),len(FP); rawTP,rawFP=30490,4842   # decisionflow 全集(含 n<6 不可驗證)
def cnt(g,k): return sum(1 for r in g if r["fine_conf"]==k)
def bb(g): return sum(1 for r in g if r["coarse_conf"] in("CONFIRMED","NEAR_CONFIRMED"))
order=[("CONFIRMED","對齊 germline cis-ASM","#0d9488"),("NEAR_CONFIRMED","近確認(邊界)","#0891b2"),
       ("REAL_NOVEL","真實+大跳不對齊=subclone候選","#7c3aed"),("REAL_DIFFUSE","真實但diffuse/無法歸因","#d97706"),
       ("NO_CLEAR_SPLIT","真無結構(無甲基切群訊號)","#6b7280")]
def b64(p):
    fp=os.path.join(A,p); return "data:image/png;base64,"+base64.b64encode(open(fp,"rb").read()).decode() if os.path.exists(fp) else None

# 表1 rows
t1=""
for k,desc,col in order:
    t,f=cnt(TP,k),cnt(FP,k); tp,fp=100*t/nTP,100*f/nFP; en=tp/fp if fp>0 else 0
    enc="#0d9488" if en>1.3 else "#db2777" if en<0.8 else "#888"
    t1+=(f"<tr><td><span class='dot' style='background:{col}'></span><b>{k.replace('_SPLIT','')}</b><br><span class='mut'>{desc}</span></td>"
         f"<td>{t}</td><td><b>{tp:.1f}%</b></td><td>{f}</td><td>{fp:.1f}%</td>"
         f"<td style='color:{enc};font-weight:700'>{en:.2f}×</td></tr>")
# 層級(TP)
sig=nTP-cnt(TP,"NO_CLEAR_SPLIT"); gm=cnt(TP,"CONFIRMED")+cnt(TP,"NEAR_CONFIRMED"); un=cnt(TP,"REAL_NOVEL")+cnt(TP,"REAL_DIFFUSE")
def kd(g): c=Counter(r["fine_k"] for r in g if r["fine_k"]>=2); return c
kT,kF=kd(TP),kd(FP)
HTML=f"""<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>全基因組切群分類觀察表</title><style>
body{{font-family:system-ui,"PingFang TC","Microsoft JhengHei","Noto Sans CJK TC","Droid Sans Fallback",sans-serif;background:#FAF9F5;color:#141413;margin:0;line-height:1.6}}
.wrap{{max-width:1000px;margin:0 auto;padding:26px 20px}} h1{{font-size:21px}} h2{{font-size:16px;border-left:4px solid #D97757;padding-left:9px;margin-top:26px}}
table{{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}} th,td{{border:1px solid #E3DACC;padding:5px 9px;text-align:center}} th{{background:#f0ebe2}} td:first-child{{text-align:left}}
.mut{{color:#888;font-size:11px}} .dot{{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px;vertical-align:middle}}
.box{{border:1px solid #E3DACC;border-radius:8px;padding:11px 15px;background:#fff;margin:10px 0}} .box.key{{border-color:#D97757;background:#fdf6f1}}
.bar{{height:18px;border-radius:3px;display:inline-block}}</style></head><body><div class="wrap">
<h1>全基因組切群分類觀察表 — HCC1395 tumor-only (TP+FP)</h1>
<div class="mut">big7 本機 binary · 處理 TP {nTP} + FP {nFP} = {nTP+nFP} 可分群位點（decisionflow 全集 TP {rawTP}/FP {rawFP}，差額 = n&lt;6 不可驗證）· 單樣本 ⭐2-3</div>

<div class="box key"><b>🔴 一句話</b>：唯一 somatic-相關訊號 = <b>CONFIRMED(germline 對齊 cis-ASM)TP 富集 3.3×</b>；<b>subclone 候選(REAL_NOVEL)非 somatic 特異(0.65×,FP 反更多)</b>=不能宣稱 subclone。</div>

<h2>判別定義 · 三閘 → 五類</h2>
<div class="box"><b>三閘(每解析度算)</b>：
<b>① real(真實)</b> = clusterboot Jaccard(重抽 80% read 重切)− within-1-group null(打散 read 間結構的 95 百分位)<b>≥0.10</b>(切群重現且超過「無結構假穩」基線) ·
<b>② big_gap(結構)</b> = UPGMA 樹分支跳躍 <b>≥ max(8×中位, 0.4×最大 gap)</b>(樹上明顯斷層) ·
<b>③ aligned(歸因)</b> = vs germline 軸 <b>CramérV≥0.3 &amp; p&lt;0.05 &amp; Cochran e≥5</b>(對齊且樣本量可靠)。
</div>
<table>
<tr><th>類別</th><th>① real</th><th>② big_gap</th><th>③ aligned</th><th>判別定義（白話）</th></tr>
<tr><td><span class='dot' style='background:#0d9488'></span><b>CONFIRMED</b></td><td>✓</td><td>–</td><td>✓</td><td>真實 + 可靠對齊 germline = <b>cis-ASM</b>(somatic-相關）</td></tr>
<tr><td><span class='dot' style='background:#0891b2'></span><b>NEAR_CONFIRMED</b></td><td>近 0.08–0.10</td><td>–</td><td>✓</td><td>excess 差一點但對齊可靠 = <b>邊界 cis-ASM</b></td></tr>
<tr><td><span class='dot' style='background:#7c3aed'></span><b>REAL_NOVEL</b></td><td>✓</td><td>✓</td><td>✗</td><td>真實 + 樹大斷層 <b>但不跟 germline 走</b> = <b>subclone 候選</b></td></tr>
<tr><td><span class='dot' style='background:#d97706'></span><b>REAL_DIFFUSE</b></td><td>✓</td><td>✗</td><td>✗</td><td>真實但 <b>無大跳 + 無可靠對齊</b>（散/單樣本無法歸因）</td></tr>
<tr><td><span class='dot' style='background:#6b7280'></span><b>NO_CLEAR</b></td><td>✗(全 k)</td><td>–</td><td>–</td><td>甲基均勻,<b>切不出超過 null 的群</b> = 真無結構</td></tr>
</table>
<div class="mut">優先序：real&amp;aligned→CONFIRMED；否則 near&amp;aligned→NEAR；否則 real&amp;big_gap→REAL_NOVEL；否則有 real→REAL_DIFFUSE；否則→NO_CLEAR。位點取 fine=最細 supported(real 且(big_gap 或 aligned))解析度。</div>

<h2>表 1 · 五類分類 × TP/FP + 特異性富集</h2>
<table><tr><th>類別</th><th>TP n</th><th>TP %</th><th>FP n</th><th>FP %</th><th>富集 TP/FP</th></tr>{t1}</table>
<div class="mut">富集 &gt;1.3 = somatic 特異(綠) · &lt;0.8 = 非特異/FP偏多(紅) · 中間灰</div>

<h2>表 2 · 甲基訊號層級（TP {nTP}）</h2>
<table>
<tr><th>層級</th><th>n</th><th>% of TP</th><th>說明</th></tr>
<tr><td><b>有甲基結構訊號</b>(非 NO_CLEAR)</td><td>{sig}</td><td><b>{100*sig/nTP:.1f}%</b></td><td>距離矩陣可切出真實群(excess≥0.10)</td></tr>
<tr><td>├ germline 對齊(CONFIRMED+NEAR)</td><td>{gm}</td><td>{100*gm/nTP:.1f}%</td><td>✅ cis-ASM,somatic-相關(可 characterize)</td></tr>
<tr><td>├ 真實但未對齊(NOVEL+DIFFUSE)</td><td>{un}</td><td>{100*un/nTP:.1f}%</td><td>🔴 subclone 候選,但非 somatic 特異(FP 一樣多)</td></tr>
<tr><td>└ <b>無結構</b>(NO_CLEAR)</td><td>{cnt(TP,'NO_CLEAR_SPLIT')}</td><td>{100*cnt(TP,'NO_CLEAR_SPLIT')/nTP:.1f}%</td><td>多 read 但甲基均勻/無離散群</td></tr>
</table>
<div class="mut">「多少有甲基訊號」= <b>{100*sig/nTP:.1f}% 的 TP 有真實甲基切群結構</b>;其中只 <b>{100*gm/nTP:.1f}% 對齊 germline(可信 cis-ASM)</b>,{100*un/nTP:.1f}% 真實但無法歸因。</div>

<h2>表 3 · germline 骨幹 + 切幾群分布</h2>
<table>
<tr><th>set</th><th>有 germline 骨幹(coarse)</th><th>fine k=2</th><th>k=3</th><th>k=4</th><th>k=5</th><th>k=6</th></tr>
<tr><td>TP</td><td><b>{bb(TP)} ({100*bb(TP)/nTP:.1f}%)</b></td>{''.join(f'<td>{kT.get(k,0)}</td>' for k in range(2,7))}</tr>
<tr><td>FP</td><td>{bb(FP)} ({100*bb(FP)/nFP:.1f}%)</td>{''.join(f'<td>{kF.get(k,0)}</td>' for k in range(2,7))}</tr>
</table>
<div class="mut">germline 骨幹 TP 22.9% vs FP 8.0% = 2.86× → 對齊結構是 somatic-相關;切群數多落 k=2-3。</div>

{f'<h2>對比圖</h2><img src="{b64("figs_dashboard/wg_tpfp_contrast.png")}" style="max-width:100%;border:1px solid #E3DACC;border-radius:6px">' if b64("figs_dashboard/wg_tpfp_contrast.png") else ''}

<h2>結論</h2>
<div class="box">
<b>① 有甲基切群訊號:TP 73.1%</b>(多數 TP 位點附近甲基可切出真實群)。<br>
<b>② 但只 21.0% 對齊 germline(可信 cis-ASM,somatic 富集 3.3×)</b> = 唯一可宣稱 somatic-相關 characterization 訊號。<br>
<b>③ 52.1% 真實但未對齊(subclone 候選)非 somatic 特異</b>(FP 一樣多甚至更多)→ <b>單樣本無法當 subclone</b>,需 normal cis-control / 多樣本。<br>
<b>④ 26.9% 真無結構(NO_CLEAR)</b>,方法正確不誤切。
</div>
<div class="mut">資料源 cluster_redesign_wg_records.json(build f429313/a678f0e)。單樣本 ⭐2-3 characterization 非 subclone 確認。</div>
</div></body></html>"""
out=f"{A}/20260622_wg_classification_table_01.standalone.html"; open(out,"w").write(HTML)
print(f"WROTE {out}")
print(f"有訊號 TP {100*sig/nTP:.1f}% | germline對齊 {100*gm/nTP:.1f}% | 未對齊 {100*un/nTP:.1f}% | NO_CLEAR {100*cnt(TP,'NO_CLEAR_SPLIT')/nTP:.1f}%")
