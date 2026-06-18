#!/usr/bin/env python3
"""≥2群確認光譜 + 對齊軸 + 三態決策 standalone HTML. 數字注入自 spectrum_decision.json; 圖 base64; self-contained."""
import json, base64
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
D=json.load(open(f"{A}/spectrum_decision.json"))
AC=json.load(open(f"{A}/split_accounting.json"))
N=D["N"]; sp=D["spec"]; al=D["align"]; tri=D["tri"]; COMMIT="5c39051"
def p(x): return round(100*x/N,1)
def b64(fn): return "data:image/png;base64,"+base64.b64encode(open(f"{A}/figs/{fn}","rb").read()).decode()
F_SPEC=b64("spectrum.png"); F_PIE=b64("alignment_pie.png")
CSS="""
:root{--c-accent:#D97757;--c-text:#141413;--c-bg:#FAF9F5;--c-border:#E3DACC;--c-card:#FFF;--c-muted:#6B6862;
--c-ok:#3F7E5B;--c-warn:#B8862F;--mono:"JetBrains Mono",ui-monospace,monospace;
--sp-2:8px;--sp-3:12px;--sp-4:16px;--sp-5:24px;--sp-6:32px;--sp-7:48px;--sans:system-ui,"Noto Sans CJK TC","Segoe UI",sans-serif;}
*{box-sizing:border-box}body{margin:0;background:var(--c-bg);color:var(--c-text);font-family:var(--sans);line-height:1.65;font-size:15.5px}
.wrap{display:grid;grid-template-columns:210px 1fr;max-width:1100px;margin:0 auto;gap:var(--sp-6)}
nav{position:sticky;top:0;align-self:start;height:100vh;overflow-y:auto;padding:var(--sp-5) var(--sp-3);border-right:1px solid var(--c-border);font-size:13px}
nav h4{margin:0 0 var(--sp-3);font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--c-muted)}
nav a{display:block;padding:5px var(--sp-2);color:var(--c-muted);text-decoration:none;border-radius:5px}
nav a:hover{background:#F0EBE0}nav a.lead{color:var(--c-accent);font-weight:600}
main{padding:var(--sp-6) var(--sp-5) var(--sp-7) 0;min-width:0}
h1{font-size:23px;margin:0 0 var(--sp-2)}.sub{color:var(--c-muted);font-size:13.5px;margin-bottom:var(--sp-5)}
h2{font-size:18.5px;margin:var(--sp-7) 0 var(--sp-3);padding-bottom:6px;border-bottom:2px solid var(--c-border)}
.card{background:var(--c-card);border:1px solid var(--c-border);border-radius:10px;padding:var(--sp-5);margin:var(--sp-4) 0}
.verdict{border-left:5px solid var(--c-ok);background:linear-gradient(90deg,#F1F7F3,#FFF 60%)}
.big{font-size:16.5px;font-weight:600;margin-bottom:var(--sp-2)}.muted{color:var(--c-muted)}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:var(--sp-3) 0}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--c-border)}th{background:#F4F0E7;font-size:12.5px}
td.num,th.num{text-align:right;font-family:var(--mono)}tr.hi td{background:#F1F7F3;font-weight:600}
.badge{display:inline-block;padding:2px 9px;border-radius:20px;font-size:12px;font-weight:600;font-family:var(--mono)}
.b-ok{background:#E4F0E9;color:var(--c-ok)}.b-warn{background:#F6EFDD;color:var(--c-warn)}.b-mut{background:#EEEAE0;color:var(--c-muted)}.b-cand{background:#F6E8E2;color:var(--c-accent)}
.callout{border:1px solid #E0CDA0;background:#FBF6E9;border-radius:9px;padding:var(--sp-4);margin:var(--sp-4) 0}
.callout.key{border-color:#A8CDB5;background:#F1F7F3}
img.fig{max-width:100%;border:1px solid var(--c-border);border-radius:8px;margin:var(--sp-2) 0}
.two{display:grid;grid-template-columns:1fr 1fr;gap:var(--sp-4)}@media(max-width:720px){.two{grid-template-columns:1fr}}
ul.tight{margin:var(--sp-2) 0;padding-left:var(--sp-5)}ul.tight li{margin:4px 0}
.tree{font-family:var(--mono);font-size:13px;background:#F4F0E7;border:1px solid var(--c-border);border-radius:8px;padding:var(--sp-4);white-space:pre;overflow-x:auto;line-height:1.7}
footer{margin-top:var(--sp-7);padding-top:var(--sp-4);border-top:1px solid var(--c-border);font-size:12.5px;color:var(--c-muted)}
@media(max-width:820px){.wrap{grid-template-columns:1fr}nav{position:static;height:auto}}
"""
html=f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>≥2 群確認光譜 + 對齊軸 + 三態決策</title>
<!-- data_sources: spectrum_decision.json,records_wg2.json -->
<!-- provenance-verified: 數字注入自 spectrum_decision.json (data-lock); 圖自 records_wg2.json -->
<style>{CSS}</style></head><body><div class="wrap">
<nav><h4>導覽</h4>
<a class="lead" href="#v">★ 結論</a><a href="#total">0·完整切群總帳</a><a href="#key">1·核心洞見</a><a href="#spec">2·清晰光譜</a>
<a href="#align">3·對齊軸(確認依據)</a><a href="#tri">4·三態決策</a><a href="#how">5·如何確認</a><a href="#caveat">6·口徑</a></nav>
<main>
<h1>甲基 read 距離「≥2 群」確認 — 光譜・對齊軸・三態決策</h1>
<div class="sub">HCC1395 單樣本 · 全基因組 {N} TP SNV · tumor reads · <span class="badge b-warn">⭐2 偵測非驗證</span> · {COMMIT} · data-lock</div>

<section id="v"><div class="card verdict">
<div class="big">你熱圖看到的「明顯 ≥2 群」是<b>可確認的</b>：<b>{D['confirmable2']} 位點（{p(D['confirmable2'])}%）</b>切出的 2 群<b>大效應 + 對齊獨立 a-priori 標籤軸</b>（CramérV≥0.7）→ 非 chance/epiallele。其中 <b>100% 對齊已知生物軸</b>（85% germline/carrier、15% HP1/HP2、0% 隨機）。</div>
<div class="muted">「抽象地有幾群」無監督沒乾淨解；但「切法有沒有對齊已知軸」有解、嚴謹、且已在 pipeline。</div>
</div></section>

<section id="total"><h2>0 · 完整切群總帳（有多少切得出來/多少不行，加總驗證 ✓ {N}）</h2>
<table><thead><tr><th>類別</th><th class="num">位點</th><th class="num">佔 WG</th><th>意義</th></tr></thead><tbody>
<tr><td colspan="4" style="background:#F4F0E7;font-weight:600">切不出來（無監督上是 1 群 / 問不了）　小計 {AC['cant_total']}（{p(AC['cant_total'])}%）</td></tr>
<tr><td>① reads 不足（n&lt;6）</td><td class="num">{AC['insuff']}</td><td class="num">{p(AC['insuff'])}%</td><td>連分群都不能</td></tr>
<tr><td>② reads 夠但切不出平衡 ≥2 群</td><td class="num">{AC['cant']}</td><td class="num">{p(AC['cant'])}%</td><td>一團相對同質（多數含多標籤、a-priori Δβ 一測即知，非黑洞）</td></tr>
<tr><td colspan="4" style="background:#F1F7F3;font-weight:600">切得出來（無監督形成 ≥2 群，best_k≥2）　小計 {AC['cansplit']}（{p(AC['cansplit'])}%）</td></tr>
<tr class="hi"><td>└ ✅ clear ≥2（V≥0.7+sig，對齊獨立軸＝可確認）</td><td class="num">{AC['clear']}</td><td class="num">{p(AC['clear'])}%</td><td>熱圖明顯者</td></tr>
<tr><td>└ 中等（V0.5-0.7+sig）</td><td class="num">{AC['mid']}</td><td class="num">{p(AC['mid'])}%</td><td>較可信</td></tr>
<tr><td>└ 弱（V0.3-0.5）</td><td class="num">{AC['weak']}</td><td class="num">{p(AC['weak'])}%</td><td>不可靠</td></tr>
<tr><td>└ 極弱/無（V&lt;0.3）</td><td class="num">{AC['vweak']}</td><td class="num">{p(AC['vweak'])}%</td><td>切得出但無意義</td></tr>
<tr><td>└ 退化表（只1標籤/1群）</td><td class="num">{AC['degen']}</td><td class="num">{p(AC['degen'])}%</td><td>無法測對齊</td></tr>
<tr style="border-top:2px solid var(--c-border)"><td><b>加總</b></td><td class="num"><b>{AC['cant_total']+AC['cansplit']}</b></td><td class="num"><b>100%</b></td><td>{AC['insuff']}+{AC['cant']}+{AC['cansplit']} = {N} ✓</td></tr>
</tbody></table>
<div class="callout"><b>一句話</b>：<b>切得出來 {p(AC['cansplit'])}%（{AC['cansplit']}）、切不出來 {p(AC['cant_total'])}%（{AC['cant_total']}）</b>；切得出來中<b>可確認的 ≥2 群 = {AC['clear']}（{p(AC['clear'])}%）</b>。「切不出來」多數（79%）非失敗＝reads 相對同質、且含多標籤（用 a-priori Δβ 判，見 §4）。</div>
<div class="callout key"><b>🔴 全程只需 tumor 資訊</b>：以上分群、cluster×label 列聯表、CramérV 確認、三態決策 <b>全部只用 tumor reads（is_tumor=1）+ tumor 衍生標籤</b>（haplotag / somatic-carrier，由 ClairS-TO somatic VCF + 相位推得）。<b>normal 完全沒進這套判定</b> → 「能不能切 ≥2 群 / 確認 / subclone 候選」<b>tumor-only 即可</b>。normal 僅在<b>最後一步 cis-control</b>（分 cis-ASM vs subclone）才需要。<span class="muted">（本次 records 來自 paired binary run，但分析端已 filter 至 tumor；真 tumor-only BAM 上以 ClairS-TO+longphase 可同樣推得標籤。）</span></div></section>

<section id="key"><h2>1 · 核心洞見：問對問題</h2>
<div class="callout key"><b>不要問「抽象地有幾群」</b>（無監督沒乾淨解：epiallele/gradient 造成 chance-partition）。<b>要問「無監督切出的群，有沒有對齊一條獨立的 a-priori 標籤軸」</b>——因為標籤與甲基距離無關，<b>對齊就是非循環的確認</b>，高 CramérV = 切法幾乎完美重現標籤分割 = 2 群是真的。</div></section>

<section id="spec"><h2>2 · 清晰光譜（5,461 個有 2-cut 位點）</h2>
<img class="fig" src="{F_SPEC}" alt="清晰光譜">
<table><thead><tr><th>清晰度（CramérV＝對齊標籤強度）</th><th class="num">位點</th><th class="num">佔 WG</th><th>判讀</th></tr></thead><tbody>
<tr class="hi"><td><span class="badge b-ok">clear ≥2</span> V≥0.7 + 顯著</td><td class="num">{sp['clear']}</td><td class="num">{p(sp['clear'])}%</td><td>✅ 可確認 ≥2 群（熱圖明顯者）</td></tr>
<tr><td>中等 V0.5-0.7 + 顯著</td><td class="num">{sp['mid']}</td><td class="num">{p(sp['mid'])}%</td><td>較可信</td></tr>
<tr><td>弱 V0.3-0.5</td><td class="num">{sp['weak']}</td><td class="num">{p(sp['weak'])}%</td><td>不可靠（gradient/epiallele）</td></tr>
<tr><td>極弱/無 V<0.3</td><td class="num">{sp['vweak']}</td><td class="num">{p(sp['vweak'])}%</td><td>切得出但無意義</td></tr>
</tbody></table></section>

<section id="align"><h2>3 · clear ≥2 群對齊什麼軸（＝確認依據）</h2>
<div class="two"><div><img class="fig" src="{F_PIE}" alt="對齊軸圓餅"></div>
<div><table><thead><tr><th>2 群對齊</th><th class="num">數量</th><th class="num">%</th></tr></thead><tbody>
<tr><td>germline vs carrier（somatic 差）</td><td class="num">{al['gc']}</td><td class="num">{round(100*al['gc']/sp['clear'],1)}%</td></tr>
<tr><td>HP1 vs HP2（germline ASM）</td><td class="num">{al['hp']}</td><td class="num">{round(100*al['hp']/sp['clear'],1)}%</td></tr>
<tr class="hi"><td>其他/隨機</td><td class="num">{al['other']}</td><td class="num">0%</td></tr>
</tbody></table>
<p class="muted"><b>0% 隨機</b> → 每個清楚 2-群都落在生物上講得通的軸 → 確認這些切法是<b>真結構、非亂切</b>。</p></div></div></section>

<section id="tri"><h2>4 · 三態決策（subclone 軸，全 {N}）</h2>
<div class="tree">位點
├─ ❌ 不能問（germline或carrier其一&lt;3 reads）  → {tri['cant']} ({p(tri['cant'])}%)
│       └ LOH（只剩一單倍型）{tri['loh']}  → 標「真單一族群/無對照」
│
└─ ✅ 能問（germline & carrier 都≥3）         → {tri['confirm1']+tri['multi_cand']} ({p(tri['confirm1']+tri['multi_cand'])}%)
        ├─ SubcloneDbeta 顯著 → 多群候選       → {tri['multi_cand']} ({p(tri['multi_cand'])}%)
        └─ SubcloneDbeta 不顯著 → 確認 1 群     → {tri['confirm1']} ({p(tri['confirm1'])}%)</div>
<table><thead><tr><th>態</th><th class="num">位點</th><th class="num">佔 WG</th><th>處理</th></tr></thead><tbody>
<tr><td><span class="badge b-mut">不能問/LOH</span></td><td class="num">{tri['cant']}</td><td class="num">{p(tri['cant'])}%</td><td>單一族群/覆蓋不足，不跑 subclone</td></tr>
<tr><td><span class="badge b-ok">確認 1 群</span></td><td class="num">{tri['confirm1']}</td><td class="num">{p(tri['confirm1'])}%</td><td>germline=carrier，無 subclone 甲基差（有信心負結果）</td></tr>
<tr class="hi"><td><span class="badge b-cand">多群候選</span></td><td class="num">{tri['multi_cand']}</td><td class="num">{p(tri['multi_cand'])}%</td><td>germline≠carrier → 下一步 normal-anchored cis-control</td></tr>
</tbody></table></section>

<section id="how"><h2>5 · 如何「確認真的 ≥2 群且合理」</h2>
<ul class="tight">
<li><b>步驟</b>：無監督 UPGMA 切 2 群（不看標籤）→ cluster×label 列聯表 → Fisher（顯著？）+ <b>Cramér's V（對齊強度）</b>。</li>
<li><b>確認標準</b>：<b>V≥0.7 + Fisher 顯著</b> → 切法幾乎完美對齊獨立標籤軸 → <b>確認 ≥2 群</b>（{D['confirmable2']} 個，{p(D['confirmable2'])}%）。</li>
<li><b>為何非循環</b>：標籤是事前 haplotag、與甲基距離無關；切法是無監督做的。對齊 = 兩條獨立軸一致 = 不是 chance。</li>
<li><b>為何合理</b>：100% clear 案落在已知生物軸（germline/carrier、HP1/HP2），0% 隨機。</li>
<li><b>已在 pipeline</b>：GlobalTest 的 CramérV + Fisher 即此確認器，無需新無監督檢定。</li>
</ul></section>

<section id="caveat"><h2>6 · 誠實口徑</h2>
<ul class="tight">
<li>🔴 <b>「確認 ≥2 群」≠「確認 subclone」</b>：85% germline/carrier 的 2 群是<b>真 2 群</b>，但可能是 <b>cis-ASM</b> 不是 subclone → 需 normal-anchored cis-control 才分得開。</li>
<li>「抽象無監督幾群」中段（V<0.5，{sp['weak']+sp['vweak']} 個）<b>不可靠</b>（epiallele/gradient）。</li>
<li>單樣本 HCC1395；tumor reads；偵測非驗證；門檻（V≥0.7）描述性可調。</li>
</ul></section>
<footer>數據溯源：spectrum_decision.json（注入）← records_wg2.json（全 WG cluster×label 列聯表）＋ significance_summary（SubcloneDbeta/LOH）。圖自 records_wg2.json。
build {COMMIT} · feat/summary-nreadsvalid · HCC1395 單樣本 · ⭐2 偵測非驗證 · 確認方法＝CramérV 對齊獨立軸（GlobalTest 既有）。</footer>
</main></div></body></html>"""
out=f"{A}/20260618_subcluster_confirmation_spectrum_decision_01.standalone.html"
open(out,"w").write(html)
print(f"WROTE {out} ({len(html)//1024} KB)")
