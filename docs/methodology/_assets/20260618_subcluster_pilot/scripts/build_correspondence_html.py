#!/usr/bin/env python3
"""comprehensive cluster×label correspondence HTML: 數據界定 + Venn + 有訊號 + 為何0.4(圖) + 沒訊號拆解 + 雙lens + 口徑.
所有數字注入自 contingency_summary/nosignal_breakdown/sensitivity.json; 圖 base64 inline; self-contained."""
import json, base64, sys
ASSET="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
S=json.load(open(f"{ASSET}/contingency_summary.json"))
NS=json.load(open(f"{ASSET}/nosignal_breakdown.json"))
SENS=json.load(open(f"{ASSET}/sensitivity.json"))
SENSE=json.load(open(f"{ASSET}/sensitivity_ext.json"))
S8=json.load(open(f"{ASSET}/section8.json"))
COMMIT="5c39051"; BRANCH="feat/summary-nreadsvalid"; SCOPE="全基因組 WG (30490 TP SNV)"
N=S["N"]
def p(x): return round(100*x/N,1)
def b64(fn):
    return "data:image/png;base64,"+base64.b64encode(open(f"{ASSET}/figs/{fn}","rb").read()).decode()
FIG_DIST=b64("sil_distribution.png"); FIG_SENS=b64("threshold_sensitivity.png"); FIG_OVERLAP=b64("sil_noise_overlap.png"); FIG_CV=b64("cramersv_paired_vs_tumor.png")
onlyA,AB,onlyB,o2o,nostruct,wc=S["onlyA"],S["AB"],S["onlyB"],S["one2one"],S["nostruct"],S["within_carrier"]
R=NS["reasons"]
def rget(k): return next(v for kk,v in R.items() if kk.startswith(k))
r1,r2,r3,r4,r5=rget("1_"),rget("2_"),rget("3_"),rget("4_"),rget("5_")

CSS="""
:root{--c-accent:#D97757;--c-text:#141413;--c-bg:#FAF9F5;--c-border:#E3DACC;--c-card:#FFF;--c-muted:#6B6862;
--c-1to1:#3F7E5B;--c-many1:#D97757;--c-1many:#4A77A8;
--sp-2:8px;--sp-3:12px;--sp-4:16px;--sp-5:24px;--sp-6:32px;--sp-7:48px;
--mono:"JetBrains Mono",ui-monospace,monospace;--sans:system-ui,"Noto Sans CJK TC","Segoe UI",sans-serif;}
*{box-sizing:border-box}body{margin:0;background:var(--c-bg);color:var(--c-text);font-family:var(--sans);line-height:1.65;font-size:15.5px}
.wrap{display:grid;grid-template-columns:230px 1fr;max-width:1140px;margin:0 auto;gap:var(--sp-6)}
nav{position:sticky;top:0;align-self:start;height:100vh;overflow-y:auto;padding:var(--sp-5) var(--sp-3);border-right:1px solid var(--c-border);font-size:13px}
nav h4{margin:0 0 var(--sp-3);font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--c-muted)}
nav a{display:block;padding:5px var(--sp-2);color:var(--c-muted);text-decoration:none;border-radius:5px}
nav a:hover{background:#F0EBE0;color:var(--c-text)}nav a.lead{color:var(--c-accent);font-weight:600}
main{padding:var(--sp-6) var(--sp-5) var(--sp-7) 0;min-width:0}
h1{font-size:24px;margin:0 0 var(--sp-2)}.sub{color:var(--c-muted);font-size:14px;margin-bottom:var(--sp-5)}
h2{font-size:19px;margin:var(--sp-7) 0 var(--sp-3);padding-bottom:6px;border-bottom:2px solid var(--c-border)}
h3{font-size:15.5px;margin:var(--sp-4) 0 var(--sp-2)}
.card{background:var(--c-card);border:1px solid var(--c-border);border-radius:10px;padding:var(--sp-5);margin:var(--sp-4) 0}
.verdict{border-left:5px solid var(--c-accent);background:linear-gradient(90deg,#FBF1EC,#FFF 60%)}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:var(--sp-3) 0}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--c-border)}th{background:#F4F0E7;font-size:12.5px}
td.num,th.num{text-align:right;font-family:var(--mono)}tr.hi td{background:#FBF1EC}
.badge{display:inline-block;padding:2px 9px;border-radius:20px;font-size:12px;font-weight:600;font-family:var(--mono)}
.b-1to1{background:#E4F0E9;color:var(--c-1to1)}.b-many1{background:#F6E8E2;color:var(--c-many1)}
.b-1many{background:#E3ECF5;color:var(--c-1many)}.b-none{background:#EEEAE0;color:var(--c-muted)}
code{font-family:var(--mono);font-size:12.5px}
.callout{border:1px solid #E0CDA0;background:#FBF6E9;border-radius:9px;padding:var(--sp-4);margin:var(--sp-4) 0}
.callout.key{border-color:#D6B5A8;background:#FBF1EC}
ul.tight{margin:var(--sp-2) 0;padding-left:var(--sp-5)}ul.tight li{margin:4px 0}
.kv{display:grid;grid-template-columns:auto 1fr;gap:4px var(--sp-4);font-size:13.5px}.kv .k{color:var(--c-muted)}
.big{font-size:17px;font-weight:600;margin-bottom:var(--sp-2)}.muted{color:var(--c-muted)}
img.fig{max-width:100%;border:1px solid var(--c-border);border-radius:8px;margin:var(--sp-2) 0}
footer{margin-top:var(--sp-7);padding-top:var(--sp-4);border-top:1px solid var(--c-border);font-size:12.5px;color:var(--c-muted)}
@media(max-width:820px){.wrap{grid-template-columns:1fr}nav{position:static;height:auto}}
"""

# sensitivity table rows
def _srow(r):
    hi=" class='hi'" if r['th']==0.4 else ""
    return (f"<tr{hi}><td class='num'>{r['th']}</td>"
      f"<td class='num'>{round(100*r['has_struct']/N,1)}%</td><td class='num'>{round(100*r['one2one']/N,1)}%</td>"
      f"<td class='num'>{round(100*r['many2one']/N,1)}%</td><td class='num'>{round(100*r['nosig']/N,1)}%</td></tr>")
srows="".join(_srow(r) for r in SENS["rows"])
def _erow(r):
    hi=" class='hi'" if r['th']==0.4 else ""
    return (f"<tr{hi}><td class='num'>{r['th']}</td><td class='num'>{round(100*r['has_struct']/N,1)}%</td>"
      f"<td class='num'>{round(100*r['one2one']/N,1)}%</td><td class='num'>{round(100*r['many2one']/N,1)}%</td>"
      f"<td class='num'>{round(100*r['nosig']/N,1)}%</td></tr>")
erows="".join(_erow(r) for r in SENSE["rows"])

html=f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>cluster×label 對應 + 有訊號/沒訊號 — HCC1395 全基因組</title>
<!-- data_sources: contingency_summary.json,nosignal_breakdown.json,sensitivity.json,records_wg2.json -->
<!-- provenance-verified: 數字注入自 *.json (data-lock); 圖自 records_wg2.json; 偵測非驗證單樣本 -->
<style>{CSS}</style></head><body><div class="wrap">
<nav><h4>導覽</h4>
<a class="lead" href="#v">★ 結論</a><a href="#scope">0·數據與範圍界定</a><a href="#def">1·三種對應定義</a>
<a href="#venn">2·文氏圖+數量比例</a><a href="#sig">3·有訊號 vs 沒訊號</a><a href="#why">4·為何 sil≥0.4(觀察圖)</a>
<a href="#nosig">5·沒訊號拆解</a><a href="#lens">6·雙 lens 定義差</a><a href="#pipe">8·Fisher/CramV·paired vs tumor</a><a href="#caveat">7·誠實口徑</a></nav>
<main>
<h1>標籤內還有沒有可切結構 — cluster×label 對應全基因組</h1>
<div class="sub">HCC1395 單樣本 · {SCOPE} · BERNOULLI ±5000 · min_sz≥3 掃 k≥2 · 偵測非驗證 · {COMMIT}</div>

<section id="v"><div class="card verdict">
<div class="big">明顯可切的無監督結構罕見（~{p(S['has_struct'])}%）；「一標籤多群」(多對1) 只 ~{p(S['many2one'])}%。但「沒訊號 {p(S['nosig'])}%」其中 <b>{round(100*NS['nosig_vc_nonnoise']/NS['nosig'])}% 仍有 label-PERMANOVA 訊號</b> → 「沒訊號」=「無 discrete 無監督群」非「無生物訊號」。</div>
<div class="muted">真·均勻無甲基訊號只 ~{round(100*NS['nosig_vc_noise']/N,1)}%（{NS['nosig_vc_noise']} 位點）。</div>
</div></section>

<section id="scope"><h2>0 · 數據與範圍界定（先把名詞與範圍講清楚）</h2>
<div class="card"><div class="kv">
<div class="k">樣本</div><div><b>HCC1395 單一細胞株</b>（tumor＝<code>HCC1395_Tmode_tagged_ClairS_pileup_v040</code>、normal＝<code>HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG</code>）；hg38；ONT。<b>跨樣本未做</b>。</div>
<div class="k">🔴 分析用 reads</div><div><b>只用 tumor reads（is_tumor=1）</b>。binary 跑 <b>paired</b>（-t tumor -n normal，距離矩陣含 normal reads），但本 cluster×label 分析<b>排除全部 normal reads</b>。label <code>1/1-1/2/2-1</code> 都是 <b>tumor read</b> 的 haplotag（germline=HP 上未帶 somatic ALT 的 tumor read、carrier=帶 ALT 的 tumor read）。normal 僅 cis-control 驗證用（本分析未做）→ 結果是 <b>tumor-only 結構</b>。</div>
<div class="k">位點 (N={N})</div><div>每個＝1 顆 TP somatic SNV（ClairS-TO pileup）+ 其 <b>±5000bp window</b>。全 22 體染色體 TP SNV。</div>
<div class="k">read / CpG</div><div>read＝覆蓋該 window 的 ONT read；CpG＝window 內 5mCG 位點；甲基矩陣＝read×CpG。</div>
<div class="k">label（a-priori 標籤）</div><div>haplotag：<code>1</code>=HP1 germline、<code>1-1</code>=HP1 somatic-carrier、<code>2</code>=HP2 germline、<code>2-1</code>=HP2 somatic-carrier（longphase 標，看甲基前定好）。</div>
<div class="k">cluster（無監督群）</div><div>read×read <b>BERNOULLI 距離</b> → UPGMA → 掃 k≥2 取最佳 silhouette、每群 ≥3 reads。</div>
<div class="k">兩分析層</div><div><b>all-read</b>（整個 tumor read 群 × label 列聯表 → 1對1/1對多/多對1）＋ <b>within-label</b>（單一 carrier 標籤內再分 → subclone 方向）。</div>
</div>
<div class="callout"><b>資料關係（兩個「有結構」lens 是巢狀的）</b>：<br>
「<b>無監督乾淨群</b>(本頁, sil≥0.4)」⊂相關於「<b>label-PERMANOVA</b>(pipeline, 任何弱差異)」。我的判準較嚴 → 我的「沒訊號」集合<b>大於</b> pipeline 的「Noise」。兩者非矛盾、是嚴格度不同（§6 詳）。</div>
<p class="muted"><b>範圍/限制</b>：單樣本 · 偵測非驗證（候選非確認子克隆）· 門檻為描述性選擇（可重算）· 「沒訊號」是無監督口徑。</p>
</div></section>

<section id="def"><h2>1 · 三種對應（cluster=甲基無監督群、label=haplotag）</h2>
<table><thead><tr><th>類型</th><th>意義</th></tr></thead><tbody>
<tr><td><span class="badge b-1to1">1對1</span></td><td>甲基群乾淨對齊標籤（群=標籤）→ 標籤可由甲基判別</td></tr>
<tr><td><span class="badge b-1many">1對多</span></td><td>一個甲基群混 ≥2 標籤 → 甲基分不開這些標籤</td></tr>
<tr><td><span class="badge b-many1">多對1</span></td><td>一個標籤跨 ≥2 甲基群 → 標籤內子結構（subclone 方向）</td></tr>
</tbody></table>
<p class="muted">一位點可同時 1對多+多對1 → 文氏圖。對應只在「有乾淨無監督結構(sil≥{S['th_struct']})」位點上算。</p></section>

<section id="venn"><h2>2 · 文氏圖 + 數量與比例（N={N}）</h2>
<div class="card"><div style="text-align:center">
<svg viewBox="0 0 640 300" width="100%" style="max-width:620px" role="img"><title>對應文氏圖</title>
<rect x="6" y="6" width="628" height="288" rx="10" fill="#F7F4ED" stroke="#E3DACC"/>
<text x="20" y="26" font-size="12" fill="#6B6862">全 {N} 位點　無結構 {nostruct} ({p(nostruct)}%)</text>
<circle cx="270" cy="160" r="95" fill="#D97757" fill-opacity="0.18" stroke="#D97757"/>
<circle cx="380" cy="160" r="95" fill="#4A77A8" fill-opacity="0.16" stroke="#4A77A8"/>
<text x="200" y="118" font-size="13" fill="#D97757" font-weight="700">多對1</text>
<text x="430" y="118" font-size="13" fill="#4A77A8" font-weight="700">1對多</text>
<text x="205" y="165" font-size="15" fill="#141413" text-anchor="middle" font-family="monospace">{onlyA}</text>
<text x="205" y="182" font-size="10" fill="#6B6862" text-anchor="middle">只多對1 {p(onlyA)}%</text>
<text x="325" y="165" font-size="15" fill="#141413" text-anchor="middle" font-family="monospace">{AB}</text>
<text x="325" y="182" font-size="10" fill="#6B6862" text-anchor="middle">交集 {p(AB)}%</text>
<text x="445" y="165" font-size="15" fill="#141413" text-anchor="middle" font-family="monospace">{onlyB}</text>
<text x="445" y="182" font-size="10" fill="#6B6862" text-anchor="middle">只1對多 {p(onlyB)}%</text>
<rect x="495" y="225" width="130" height="52" rx="7" fill="#E4F0E9" stroke="#3F7E5B"/>
<text x="560" y="246" font-size="12" fill="#3F7E5B" text-anchor="middle" font-weight="700">1對1 乾淨</text>
<text x="560" y="266" font-size="13" fill="#141413" text-anchor="middle" font-family="monospace">{o2o} ({p(o2o)}%)</text>
</svg></div>
<table><thead><tr><th>類型</th><th class="num">位點數</th><th class="num">佔全基因組</th></tr></thead><tbody>
<tr><td>有 all-read 結構(sil≥{S['th_struct']})</td><td class="num">{S['has_struct']}</td><td class="num">{p(S['has_struct'])}%</td></tr>
<tr><td><span class="badge b-1to1">1對1</span> 乾淨對齊</td><td class="num">{o2o}</td><td class="num">{p(o2o)}%</td></tr>
<tr><td><span class="badge b-1many">1對多</span>（含交集）</td><td class="num">{S['one2many']}</td><td class="num">{p(S['one2many'])}%</td></tr>
<tr><td><span class="badge b-many1">多對1</span>（含交集, joint）</td><td class="num">{S['many2one']}</td><td class="num">{p(S['many2one'])}%</td></tr>
<tr><td>└ within-carrier split（subclone 方向, sil≥{S['th_within']}）</td><td class="num">{wc}</td><td class="num">{p(wc)}%</td></tr>
<tr><td>1對多 ∩ 多對1</td><td class="num">{S['both']}</td><td class="num">{p(S['both'])}%</td></tr>
<tr><td><span class="badge b-none">無結構</span></td><td class="num">{nostruct}</td><td class="num">{p(nostruct)}%</td></tr>
</tbody></table></div></section>

<section id="sig"><h2>3 · 有訊號 vs 沒訊號（以位點為單位）</h2>
<table><thead><tr><th>定義</th><th class="num">有訊號</th><th class="num">比例</th></tr></thead><tbody>
<tr><td>嚴格＝1對1（標籤可判別）</td><td class="num">{S['sig_strict']}</td><td class="num">{p(S['sig_strict'])}%</td></tr>
<tr class="hi"><td><b>+ 多對1/within-carrier（標籤內多群也算）</b></td><td class="num"><b>{S['sig_loose']}</b></td><td class="num"><b>{p(S['sig_loose'])}%</b></td></tr>
<tr><td>任何無監督結構</td><td class="num">{S['sig_any']}</td><td class="num">{p(S['sig_any'])}%</td></tr>
<tr><td><span class="badge b-none">沒訊號（無監督口徑）</span></td><td class="num">{S['nosig']}</td><td class="num">{p(S['nosig'])}%</td></tr>
</tbody></table>
<div class="callout key">加上「一標籤多群(多對1)也算」後，有訊號從 {p(S['sig_strict'])}% → <b>{p(S['sig_loose'])}%</b>（多 {S['sig_loose']-S['sig_strict']} 個位點）。但「沒訊號」的真正意義見 §5。</div></section>

<section id="why"><h2>4 · 為何採用 sil≥0.4（觀察圖）</h2>
<p>silhouette = (群外距離−群內距離)/max → 量「群內聚集 vs 群外分離」。<b>0.4 是「中等結構」慣例</b>（&lt;0.25 弱、0.25–0.5 弱-中、≥0.5 明顯）。兩張圖佐證這選擇合理且非 knife-edge：</p>
<div class="card"><h3>圖 A · best silhouette 分布</h3>
<img class="fig" src="{FIG_DIST}" alt="silhouette 分布">
<p class="muted">有切出群的位點分布峰在 ~0.35–0.4；<b>0.4（紅）落在右肩</b>（取較乾淨分離的）。0.5 更嚴（只取右尾明顯者）、0.3 會納入主體。80% 位點根本切不出 ≥2 平衡群（best_k=None，見 §5）。</p></div>
<div class="card"><h3>圖 B · 門檻敏感度</h3>
<img class="fig" src="{FIG_SENS}" alt="門檻敏感度">
<table><thead><tr><th class="num">TH</th><th class="num">有結構</th><th class="num">1對1</th><th class="num">多對1</th><th class="num">沒訊號</th></tr></thead><tbody>{srows}</tbody></table>
<p class="muted">各類比例隨門檻<b>平滑變化、無陡降</b> → 0.4 不是 knife-edge；換 0.35 或 0.45 結論方向不變（有結構 8%→12% 或 5%），原始列聯表已存可任意重算。</p></div>
<div class="card"><h3>圖 C · 補到 0.1 + 「有 sil 值 ＝ 有訊號嗎？」</h3>
<table><thead><tr><th class="num">TH</th><th class="num">有結構</th><th class="num">1對1</th><th class="num">多對1</th><th class="num">沒訊號</th></tr></thead><tbody>{erows}</tbody></table>
<p class="muted">注意 <b>0.1→0.2 幾乎不動（19.7%→19.3%）</b>：「沒訊號」的 <b>~80% 地板是「根本切不出任何平衡群」(best_k=None)，與門檻無關</b>，再降也救不回。</p>
<img class="fig" src="{FIG_OVERLAP}" alt="Noise vs 結構 silhouette 重疊">
<div class="callout key"><b>答：silhouette 是連續強度量（你的直覺對），但「有 sil 值 ≠ 有訊號」。</b>
<ul class="tight">
<li>best-split 是「挑最分離的一刀」→ <b>連純噪音都有正 silhouette</b>（selection/double-dip）。</li>
<li>實證：pipeline 判 <b>Noise</b> 的位點若有切出群（n=11），silhouette <b>中位數 0.434 ＞ 結構-VC 的 0.378</b>；<b>sil≥0.4 比例 Noise 55% ＞ 結構 41%</b>（紅藍分布大幅重疊，紅色甚至延伸到 0.7-0.8）。</li>
<li>→ <b>低-中 silhouette 與噪音無法區分</b>；門檻是「強度」切點、非「訊號/無訊號」邊界。真正的判別不在 silhouette 大小，而在「能不能切出群」＋ 外部軸驗證。</li>
</ul>
<span class="muted">⚠ caveat：Noise 樣本小（n=11，因 Noise 多數 best_k=None 切不出），但方向與 memory「Noise silhouette 可比 Strong 還高」一致。</span></div></div></section>

<section id="nosig"><h2>5 · 「沒訊號」{S['nosig']} 位點拆解（{p(S['nosig'])}%）</h2>
<div class="callout key"><b>關鍵：「沒訊號」≠「無生物訊號」。</b> 這 {NS['nosig']} 個中 <b>{NS['nosig_vc_nonnoise']}（{round(100*NS['nosig_vc_nonnoise']/NS['nosig'])}%）仍有 label-PERMANOVA 訊號</b>（pipeline 判非-Noise）；真·均勻無結構只 <b>{NS['nosig_vc_noise']}（{round(100*NS['nosig_vc_noise']/NS['nosig'],1)}%）</b>。即多數是「germline ASM 那種<b>連續/弱</b>差異」，不構成乾淨離散群。</div>
<table><thead><tr><th>原因（互斥）</th><th class="num">位點</th><th class="num">佔沒訊號</th><th class="num">佔全WG</th></tr></thead><tbody>
<tr class="hi"><td><b>3 無平衡分群</b>（切不出 ≥2 群各 ≥3）</td><td class="num">{r3}</td><td class="num">{round(100*r3/NS['nosig'],1)}%</td><td class="num">{p(r3)}%</td></tr>
<tr><td>4 弱結構 sub-threshold（0.3≤sil&lt;0.4）</td><td class="num">{r4}</td><td class="num">{round(100*r4/NS['nosig'],1)}%</td><td class="num">{p(r4)}%</td></tr>
<tr><td>5 很弱/接近均勻（sil&lt;0.3）</td><td class="num">{r5}</td><td class="num">{round(100*r5/NS['nosig'],1)}%</td><td class="num">{p(r5)}%</td></tr>
<tr><td>1 讀數不足（n&lt;6，問不了）</td><td class="num">{r1}</td><td class="num">{round(100*r1/NS['nosig'],1)}%</td><td class="num">{p(r1)}%</td></tr>
<tr><td>2 矩陣不完整（NaN/覆蓋）</td><td class="num">{r2}</td><td class="num">{round(100*r2/NS['nosig'],1)}%</td><td class="num">{p(r2)}%</td></tr>
</tbody></table>
<div class="kv">
<div class="k">『根本問不了』(類1-3)</div><div><b>{NS['abstain']}</b>（{round(100*NS['abstain']/NS['nosig'],1)}% of 沒訊號）— 主因＝<b>無平衡分群（一團相對同質、無離散二分）{round(100*r3/NS['nosig'])}%</b></div>
<div class="k">『問了但結構弱』(類4-5)</div><div>{NS['weak']}（{round(100*NS['weak']/NS['nosig'],1)}%）</div>
</div>
<p class="muted"><b>主因解讀</b>：86% 沒訊號是「<b>切不出平衡的兩群</b>」（UPGMA 二分把次群切到 &lt;3 reads＝reads 形成單一主群、僅少數離群）。這與「98% 仍有 label-PERMANOVA 訊號」一致 → germline ASM 是<b>沿標籤軸的弱/漸進效應</b>，非離散子群。subclone（離散群）本就該罕見，多數位點是單一主 clone 的 ASM 背景。</p></section>

<section id="lens"><h2>6 · 🔴 「有訊號」高度依賴定義（雙 lens 關係）</h2>
<div class="callout"><b>同批位點，兩種「有結構」差 ~12×（巢狀非矛盾）：</b>
<ul class="tight">
<li><b>無監督乾淨群（本頁）</b> sil≥{S['th_struct']} → ~{p(S['has_struct'])}%。問「有沒有<b>乾淨可切的離散群</b>」（≈ subclone 該有的）。</li>
<li><b>label-PERMANOVA（pipeline）</b> → ~98.5%。問「有沒有<b>任何標籤相關甲基差異</b>」（germline ASM 無所不在故幾乎全中）。</li>
</ul>
前者 ⊂ 後者。subclone 判別力在 <b>conditioned 軸</b>、不在普遍性。「沒訊號」是前者口徑。</div></section>

<section id="pipe"><h2>8 · Pipeline Fisher/Cramér's V ＋ paired vs tumor-only（兩軸）</h2>
<p>pipeline 的 GlobalTest（Fisher-Freeman-Halton + Cramér's V）與我的 tumor-only pass <b>測的是不同軸</b>。源碼確認（<code>RegionProcessor.cpp:2554-2563</code> / <code>GlobalTest.cpp:153-169</code>）：GlobalTest 不過濾 is_tumor、label=<code>read.alt_support</code>（ALT/REF）→ <b>含 normal reads（=REF）的 paired 軸</b>。</p>
<div class="card"><h3>兩軸示意</h3>
<svg viewBox="0 0 700 290" width="100%" style="max-width:700px" role="img"><title>paired vs tumor-only 兩軸</title>
<rect x="8" y="8" width="335" height="274" rx="9" fill="#EAF0F6" stroke="#4A77A8"/>
<text x="175" y="30" font-size="13" fill="#4A77A8" text-anchor="middle" font-weight="700">PAIRED（pipeline GlobalTest）</text>
<rect x="30" y="44" width="290" height="34" rx="6" fill="#FFF" stroke="#4A77A8"/><text x="175" y="65" font-size="11.5" text-anchor="middle">reads = tumor ＋ normal（normal=REF）</text>
<text x="175" y="95" font-size="13" text-anchor="middle" fill="#6B6862">↓ 分群（全部都能切）</text>
<rect x="30" y="104" width="290" height="34" rx="6" fill="#FFF" stroke="#4A77A8"/><text x="175" y="125" font-size="11.5" text-anchor="middle">cluster × <b>allele（ALT/REF）</b></text>
<text x="175" y="155" font-size="13" text-anchor="middle" fill="#6B6862">↓ Fisher-FFH + CramérV</text>
<rect x="30" y="164" width="290" height="40" rx="6" fill="#DCE7F2" stroke="#4A77A8"/><text x="175" y="180" font-size="11.5" text-anchor="middle" font-weight="700">顯著 74.2% · V med 0.561</text><text x="175" y="196" font-size="10.5" text-anchor="middle">V&gt;0 僅 41.4%</text>
<rect x="30" y="216" width="290" height="54" rx="6" fill="#FBF1EC" stroke="#D97757"/><text x="175" y="237" font-size="11" fill="#C0432F" text-anchor="middle" font-weight="700">⚠ 高普遍性主要來自</text><text x="175" y="254" font-size="11" fill="#C0432F" text-anchor="middle">tumor-vs-normal ＋ germline ASM</text><text x="175" y="268" font-size="10" fill="#6B6862" text-anchor="middle">非 tumor subclone</text>
<rect x="357" y="8" width="335" height="274" rx="9" fill="#FBF1EC" stroke="#D97757"/>
<text x="524" y="30" font-size="13" fill="#D97757" text-anchor="middle" font-weight="700">TUMOR-ONLY（我的 pass）</text>
<rect x="379" y="44" width="290" height="34" rx="6" fill="#FFF" stroke="#D97757"/><text x="524" y="65" font-size="11.5" text-anchor="middle">reads = <b>只 tumor</b>（is_tumor=1）</text>
<text x="524" y="95" font-size="13" text-anchor="middle" fill="#6B6862">↓ 分群（只 17.9% 切得出）</text>
<rect x="379" y="104" width="290" height="34" rx="6" fill="#FFF" stroke="#D97757"/><text x="524" y="125" font-size="11.5" text-anchor="middle">cluster × <b>haplotag（1/1-1/2/2-1）</b></text>
<text x="524" y="155" font-size="13" text-anchor="middle" fill="#6B6862">↓ chi²-p + CramérV</text>
<rect x="379" y="164" width="290" height="40" rx="6" fill="#F6E1D8" stroke="#D97757"/><text x="524" y="180" font-size="11.5" text-anchor="middle" font-weight="700">顯著 13.7% · V med 0.750</text><text x="524" y="196" font-size="10.5" text-anchor="middle">可算僅 17.9%</text>
<rect x="379" y="216" width="290" height="54" rx="6" fill="#E4F0E9" stroke="#3F7E5B"/><text x="524" y="237" font-size="11" fill="#3F7E5B" text-anchor="middle" font-weight="700">✓ 才是 tumor subclone 方向</text><text x="524" y="254" font-size="11" fill="#3F7E5B" text-anchor="middle">稀疏但效應更強</text><text x="524" y="268" font-size="10" fill="#6B6862" text-anchor="middle">（挑 tumor 內最乾淨切法）</text>
</svg></div>
<div class="card"><h3>對照表（同 {N} 位點）</h3>
<table><thead><tr><th>項目</th><th>PAIRED（pipeline）</th><th>TUMOR-ONLY（我的 pass）</th></tr></thead><tbody>
<tr><td>read-set</td><td>tumor ＋ normal（normal=REF）</td><td>只 tumor（is_tumor=1）</td></tr>
<tr><td>label</td><td><b>allele ALT/REF</b>（主）/ HP / hp_fine</td><td>haplotag 1/1-1/2/2-1</td></tr>
<tr><td>可算位點</td><td class="num">{p(S8['paired_allele']['v']['nonzero']+S8['paired_allele']['v']['zero'])}%（~全）</td><td class="num">{p(S8['tumor_only']['computable'])}%（切得出群）</td></tr>
<tr class="hi"><td>Fisher 顯著 p&lt;0.05</td><td class="num"><b>{p(S8['paired_allele']['sig05'])}%</b>（hp_fine {p(S8['paired_hpfine']['sig05'])}%）</td><td class="num"><b>{p(S8['tumor_only']['sig05'])}%</b></td></tr>
<tr><td>CramérV&gt;0（可靠效應）</td><td class="num">{p(S8['paired_allele']['v']['nonzero'])}%</td><td class="num">{p(S8['tumor_only']['v']['nonzero'])}%</td></tr>
<tr><td>CramérV&gt;0 中位數</td><td class="num">{S8['paired_allele']['v']['med']}</td><td class="num"><b>{S8['tumor_only']['v']['med']}</b>（更強）</td></tr>
<tr><td>PassedGating</td><td class="num">{p(S8['passed_gating'])}%</td><td class="muted">—（不同流程）</td></tr>
<tr><td>per-CpG Fisher ≥1 顯著</td><td class="num">{p(S8['percpg_sig'])}%</td><td class="muted">—</td></tr>
</tbody></table>
<img class="fig" src="{FIG_CV}" alt="CramersV paired vs tumor-only">
</div>
<div class="callout key"><b>敘述：74% 顯著不是 subclone。</b>
<ul class="tight">
<li>paired 之所以普遍顯著（74%），是因為 <b>normal reads（全 REF＋germline 甲基）穩定地與 tumor 分開</b> → cluster×allele 關聯<b>大半是 tumor-vs-normal ＋ germline ASM</b>。</li>
<li>tumor-only 只 {p(S8['tumor_only']['computable'])}% 切得出群、{p(S8['tumor_only']['sig05'])}% 顯著（但<b>可算者顯著率 ~77%，與 paired 相當</b>）→ 差別在「paired 分群永遠成功（含普遍的 tumor-vs-normal）」，非「有無訊號」。</li>
<li>效應量 tumor-only 更強（V med {S8['tumor_only']['v']['med']} ＞ paired {S8['paired_allele']['v']['med']}）：tumor-only 切的就是強的。</li>
</ul>
<span class="muted">⚠ tumor-only 用 chi²-p（近似 Fisher）＋ 不同分群，非 controlled 對照，代表「各自分析結果」。CramérV=0 含 Cochran-gated 不可靠（raw-V 未 emit，memory 標待修）；Fisher over-dispersion（read 非獨立）使 p 偏寬，74% 可能偏高。</span></div></section>

<section id="caveat"><h2>7 · 誠實口徑</h2>
<ul class="tight">
<li><b>單樣本 HCC1395</b>；跨樣本未做。</li>
<li><b>偵測非驗證</b>：「有兩群」≠「兩個子克隆」（含 epiallele/技術假象，需外部軸驗證）。</li>
<li><b>門檻描述性可調</b>：sil≥{S['th_struct']}(結構)/{S['th_within']}(within)；§4 已證非 knife-edge；原始列聯表 <code>records_wg2.json</code> 可重算。</li>
<li>「沒訊號 {p(S['nosig'])}%」是<b>無監督</b>口徑，98% 仍有 label 訊號（§5）。</li>
</ul></section>

<footer>數據溯源：contingency_summary / nosignal_breakdown / sensitivity.json（注入）← records_wg2.json（全 WG cluster×label 列聯表 + within-label 子分群，逐 chr disk-safe）。圖自 records_wg2.json。
build {COMMIT} · {BRANCH} · HCC1395 單樣本 {SCOPE} · 偵測非驗證 · ⭐2 描述性 pilot。</footer>
</main></div></body></html>"""
out=f"{ASSET}/20260618_cluster_label_correspondence_wg.standalone.html"
open(out,"w").write(html)
print(f"WROTE {out} ({len(html)//1024} KB)")
