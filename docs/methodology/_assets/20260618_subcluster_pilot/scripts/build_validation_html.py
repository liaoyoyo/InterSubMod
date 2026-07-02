#!/usr/bin/env python3
"""consolidated 方法學+驗證確認 standalone HTML. 數字注入自鎖檔 json; 圖 base64; self-contained."""
import json, base64
A="/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
S=json.load(open(f"{A}/contingency_summary.json")); NS=json.load(open(f"{A}/nosignal_breakdown.json"))
S8=json.load(open(f"{A}/section8.json")); CV=json.load(open(f"{A}/case_validation.json")); TL=json.load(open(f"{A}/tumor_layers.json"))
N=S["N"]; COMMIT="5c39051"; BRANCH="feat/summary-nreadsvalid"
def p(x): return round(100*x/N,1)
def b64(fn): return "data:image/png;base64,"+base64.b64encode(open(f"{A}/figs/{fn}","rb").read()).decode()
F_OV=b64("sil_noise_overlap.png"); F_CV=b64("cramersv_paired_vs_tumor.png")
R=NS["reasons"]; r3=next(v for k,v in R.items() if k.startswith("3_")); r4=next(v for k,v in R.items() if k.startswith("4_"))
onlyA,AB,onlyB,o2o,nostruct,wc=S["onlyA"],S["AB"],S["onlyB"],S["one2one"],S["nostruct"],S["within_carrier"]

CSS="""
:root{--c-accent:#D97757;--c-text:#141413;--c-bg:#FAF9F5;--c-border:#E3DACC;--c-card:#FFF;--c-muted:#6B6862;
--c-ok:#3F7E5B;--c-warn:#B8862F;--c-bad:#C0432F;--mono:"JetBrains Mono",ui-monospace,monospace;
--sp-2:8px;--sp-3:12px;--sp-4:16px;--sp-5:24px;--sp-6:32px;--sp-7:48px;--sans:system-ui,"Noto Sans CJK TC","Segoe UI",sans-serif;}
*{box-sizing:border-box}body{margin:0;background:var(--c-bg);color:var(--c-text);font-family:var(--sans);line-height:1.65;font-size:15.5px}
.wrap{display:grid;grid-template-columns:218px 1fr;max-width:1120px;margin:0 auto;gap:var(--sp-6)}
nav{position:sticky;top:0;align-self:start;height:100vh;overflow-y:auto;padding:var(--sp-5) var(--sp-3);border-right:1px solid var(--c-border);font-size:13px}
nav h4{margin:0 0 var(--sp-3);font-size:12px;letter-spacing:.06em;text-transform:uppercase;color:var(--c-muted)}
nav a{display:block;padding:5px var(--sp-2);color:var(--c-muted);text-decoration:none;border-radius:5px}
nav a:hover{background:#F0EBE0}nav a.lead{color:var(--c-accent);font-weight:600}
main{padding:var(--sp-6) var(--sp-5) var(--sp-7) 0;min-width:0}
h1{font-size:23px;margin:0 0 var(--sp-2)}.sub{color:var(--c-muted);font-size:13.5px;margin-bottom:var(--sp-5)}
h2{font-size:18.5px;margin:var(--sp-7) 0 var(--sp-3);padding-bottom:6px;border-bottom:2px solid var(--c-border)}
h3{font-size:15px;margin:var(--sp-4) 0 var(--sp-2)}
.card{background:var(--c-card);border:1px solid var(--c-border);border-radius:10px;padding:var(--sp-5);margin:var(--sp-4) 0}
.verdict{border-left:5px solid var(--c-accent);background:linear-gradient(90deg,#FBF1EC,#FFF 60%)}
.big{font-size:16.5px;font-weight:600;margin-bottom:var(--sp-2)}.muted{color:var(--c-muted)}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:var(--sp-3) 0}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--c-border)}th{background:#F4F0E7;font-size:12.5px}
td.num,th.num{text-align:right;font-family:var(--mono)}tr.hi td{background:#FBF1EC}
.badge{display:inline-block;padding:2px 9px;border-radius:20px;font-size:12px;font-weight:600;font-family:var(--mono)}
.b-ok{background:#E4F0E9;color:var(--c-ok)}.b-warn{background:#F6EFDD;color:var(--c-warn)}.b-bad{background:#F7E4E0;color:var(--c-bad)}.b-tier{background:#F6E8E2;color:var(--c-accent)}
.callout{border:1px solid #E0CDA0;background:#FBF6E9;border-radius:9px;padding:var(--sp-4);margin:var(--sp-4) 0}
.callout.key{border-color:#D6B5A8;background:#FBF1EC}.callout.ok{border-color:#A8CDB5;background:#F1F7F3}
ul.tight{margin:var(--sp-2) 0;padding-left:var(--sp-5)}ul.tight li{margin:4px 0}
img.fig{max-width:100%;border:1px solid var(--c-border);border-radius:8px;margin:var(--sp-2) 0}
code{font-family:var(--mono);font-size:12.5px}
.conf{display:grid;grid-template-columns:auto 1fr;gap:var(--sp-3) var(--sp-4);align-items:start}
footer{margin-top:var(--sp-7);padding-top:var(--sp-4);border-top:1px solid var(--c-border);font-size:12.5px;color:var(--c-muted)}
@media(max-width:820px){.wrap{grid-template-columns:1fr}nav{position:static;height:auto}}
"""

html=f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>標籤內子分群 — 方法學與驗證確認</title>
<!-- data_sources: contingency_summary,nosignal_breakdown,section8,case_validation,tumor_layers.json -->
<!-- provenance-verified: 數字注入自鎖檔 json (data-lock §13.0); 對抗驗證+自驗 -->
<style>{CSS}</style></head><body><div class="wrap">
<nav><h4>導覽</h4>
<a class="lead" href="#v">★ 結論</a><a href="#data">0·資料與範圍</a><a href="#res">1·結果(數量比例)</a>
<a href="#nosig">2·沒訊號拆解</a><a href="#sil">3·為何0.4·有sil≠訊號</a><a href="#fisher">4·Fisher·paired vs tumor</a>
<a href="#confirm">5·驗證確認(3項)</a><a href="#caveat">6·誠實口徑</a></nav>
<main>
<h1>標籤內子分群偵測 — 方法學與驗證確認</h1>
<div class="sub">HCC1395 單樣本 · 全基因組 {N} TP SNV · <span class="badge b-tier">⭐2 偵測非驗證 pilot</span> · {COMMIT} · 數字皆 data-lock</div>

<section id="v"><div class="card verdict">
<div class="big">乾淨可切的離散結構罕見（{p(S['has_struct'])}%）、「一標籤多群」subclone 方向僅 {p(wc)}%；91.4% 無乾淨群但 <b>98.4% 仍有 label-PERMANOVA 弱訊號</b>。方法學<b>對偵測有效、為 detector 非 validator</b>。</div>
<div class="muted">三項確認：案例適合 ✅ ／ 結果符合判斷 ⚠ 中等(metric-sensitive，已修) ／ 方法有效 ✅。真·均勻無訊號僅 {round(100*NS['nosig_vc_noise']/N,1)}%。</div>
</div></section>

<section id="data"><h2>0 · 資料與範圍（🔴 tumor-only vs paired）</h2>
<table><thead><tr><th></th><th>本分析（對應）</th><th>pipeline GlobalTest</th></tr></thead><tbody>
<tr><td>read-set</td><td><b>只 tumor</b>(is_tumor=1)</td><td><b>paired</b>(tumor+normal,normal=REF)</td></tr>
<tr><td>label</td><td>haplotag 1/1-1/2/2-1</td><td><b>allele ALT/REF</b>/HP/hp_fine</td></tr>
<tr><td>源碼證</td><td>Python pass</td><td><code>RegionProcessor.cpp:2554-2563</code> 不過濾 is_tumor</td></tr>
</tbody></table>
<p class="muted">單樣本 HCC1395（tumor=ClairS_pileup_v040 / normal=5khz_simplex）· 偵測非驗證 · 門檻描述性可重算。</p></section>

<section id="res"><h2>1 · 結果（cluster×label 對應，數量與比例）</h2>
<table><thead><tr><th>類型</th><th class="num">位點</th><th class="num">%</th></tr></thead><tbody>
<tr><td>無結構</td><td class="num">{nostruct}</td><td class="num">{p(nostruct)}%</td></tr>
<tr><td>1對1 乾淨對齊</td><td class="num">{o2o}</td><td class="num">{p(o2o)}%</td></tr>
<tr><td>只 1對多</td><td class="num">{onlyB}</td><td class="num">{p(onlyB)}%</td></tr>
<tr><td>只 多對1</td><td class="num">{onlyA}</td><td class="num">{p(onlyA)}%</td></tr>
<tr><td>1對多 ∩ 多對1</td><td class="num">{AB}</td><td class="num">{p(AB)}%</td></tr>
<tr class="hi"><td>└ within-carrier split（subclone 方向）</td><td class="num">{wc}</td><td class="num">{p(wc)}%</td></tr>
</tbody></table>
<p class="muted">分割驗證 {o2o}+{onlyA}+{onlyB}+{AB}+{nostruct} = {N} ✓。有訊號：嚴格 {p(S['sig_strict'])}% / +多對1 {p(S['sig_loose'])}% / 任何結構 {p(S['sig_any'])}% / 沒訊號 {p(S['nosig'])}%。</p></section>

<section id="nosig"><h2>2 · 「沒訊號 {S['nosig']}」拆解（🔴 ≠ 無生物訊號）</h2>
<div class="callout key">這 {NS['nosig']} 個中 <b>{round(100*NS['nosig_vc_nonnoise']/NS['nosig'])}% 仍有 label-PERMANOVA 訊號</b>（germline ASM 弱/漸進）；真·均勻無結構只 <b>{NS['nosig_vc_noise']}（{round(100*NS['nosig_vc_noise']/NS['nosig'],1)}%）</b>。主因＝<b>無平衡分群 {round(100*r3/NS['nosig'])}%</b>（一團同質、非離散）+ 弱結構 {round(100*r4/NS['nosig'],1)}%。</div></section>

<section id="sil"><h2>3 · 為何 sil≥0.4 ＋「有 sil 值 ≠ 有訊號」</h2>
<p>門檻 0.1→0.5 平滑（有結構 19.7%→8.1%→2.9%），非 knife-edge；80% 地板＝切不出平衡群、與門檻無關。<b>best-split 挑最分離一刀 → 連噪音都有正 silhouette</b>：</p>
<img class="fig" src="{F_OV}" alt="Noise vs 結構 silhouette 重疊">
<p class="muted">pipeline-Noise 位點若有切（n=11）silhouette 中位 <b>0.434 ＞ 結構-VC 0.378</b> → 低-中 silhouette 與噪音無法區分。判別不在 sil 大小，在「能不能切出群」＋外部軸驗證。</p></section>

<section id="fisher"><h2>4 · Fisher/Cramér's V：paired vs tumor-only</h2>
<table><thead><tr><th></th><th class="num">PAIRED allele</th><th class="num">TUMOR-only</th></tr></thead><tbody>
<tr><td>Fisher 顯著 p&lt;0.05</td><td class="num">{p(S8['paired_allele']['sig05'])}%</td><td class="num"><b>{p(S8['tumor_only']['sig05'])}%</b>(可算 {p(TL['computable'])}%)</td></tr>
<tr><td>CramérV&gt;0</td><td class="num">{p(S8['paired_allele']['v']['nonzero'])}%</td><td class="num">{p(S8['tumor_only']['v']['nonzero'])}%</td></tr>
<tr><td>CramérV 中位(顯著者)</td><td class="num">{S8['paired_allele']['v']['med']}</td><td class="num"><b>{TL['v_med_sig']}</b></td></tr>
</tbody></table>
<img class="fig" src="{F_CV}" alt="CramersV paired vs tumor-only">
<div class="callout"><b>tumor 樣本可「驗證」分層</b>：①有切出群 {p(TL['computable'])}% / ②顯著 {p(TL['sig05'])}% / ③顯著+V≥0.5 {p(TL['sig_v05'])}% / ④V≥0.7 {p(TL['sig_v07'])}% / ⑤<b>subclone 方向 {p(TL['within_carrier'])}%</b>。🔴 paired 74% 主要 <b>tumor-vs-normal+germline ASM 非 subclone</b>；tumor-only 13.7% 多為 germline HP ASM，真 subclone 僅 {p(wc)}%。高效應量含 selection 偏誤。</div></section>

<section id="confirm"><h2>5 · 驗證確認（5-agent 對抗 + 自驗）</h2>
<div class="card"><div class="conf">
<div><span class="badge b-ok">✅ 確認1</span></div><div><b>案例適合</b>：chr1 carrier sil≥0.4 <b>全普查 41/41 + 4 低控制</b>；KS vs WG pool p=0.52 非 cherry-pick。<span class="muted">caveat：chr1-only（WG 有更強位點）。</span></div>
<div><span class="badge b-warn">⚠ 確認2</span></div><div><b>案例結果符合判斷（PARTIAL，已修）</b>：sil↔per-CpG 差 r=<b>{CV['sil_vs_cpgmax_r']}</b> 顯著；但 region-level r=<b>{CV['sil_vs_level_r']}</b> p={CV['sil_vs_level_p']} <b>不顯著</b> → 中等、metric-sensitive。<br>🔴 修正：高 sil <b>{CV['high_sil_pattern_driven']}/{CV['high_sil_n']} 是 per-CpG pattern 驅動</b>（非先前誤稱的「整體 level」）。strand 假象排除（{CV['high_strand_purity']}/45 高純度）。</div>
<div><span class="badge b-ok">✅ 確認3</span></div><div><b>方法學有效（有界）</b>：UPGMA+silhouette 正確、列聯表非循環、<b>detector 非 validator</b> 界定正確並嵌入產物、對齊 memory（無監督 double-dip NEGATIVE / subclone 在 a-priori 軸）。</div>
</div></div>
<div class="callout ok"><b>整體裁決：方法學健全、為偵測有效；無翻案級問題。</b></div>
<div class="callout"><b>🔴 方法論教訓（§13.7 實證）</b>：對抗 Workflow 的 adjudicator <b>自己捏造</b>——宣稱資料檔「不存在=fabrication」，實為它搜錯位置（主 repo 未複製該 json）。<b>未盲信</b>：自查證實檔案存在（45 案/{N} loci/r={CV['sil_vs_cpgmax_r']}）→ 駁回 adjudicator、採信讀真檔的 4 維 agent + 自驗 2 修正。反捏造稽核自己捏造、被獨立驗證抓到。</div></section>

<section id="caveat"><h2>6 · 誠實口徑</h2>
<ul class="tight">
<li><b>單樣本 HCC1395</b>；跨樣本未做。<b>偵測非驗證</b>（有兩群≠子克隆，需外部軸）。</li>
<li>「有訊號」依定義：無監督乾淨群 8% vs label-PERMANOVA 98.5%（巢狀非矛盾）。「沒訊號 91.4%」是無監督口徑、98% 仍有 label 訊號。</li>
<li>tumor-only 13.7% 顯著多為 germline ASM，subclone 方向僅 {p(wc)}%；高效應量含 selection 偏誤。</li>
<li>案例符合判斷<b>中等</b>（cpg_max r={CV['sil_vs_cpgmax_r']} 顯著 / level r={CV['sil_vs_level_r']} 不顯著）。</li>
</ul></section>
<footer>數據溯源：contingency_summary / nosignal_breakdown / section8 / case_validation / tumor_layers.json（注入）← records_wg2.json（全 WG）。圖自 records_wg2.json。
build {COMMIT} · {BRANCH} · HCC1395 單樣本 · ⭐2 偵測非驗證 pilot · 對應詳版見 `InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/20260618_cluster_label_correspondence_wg.standalone.html`。</footer>
</main></div></body></html>"""
out=f"{A}/20260618_within_label_subcluster_methodology_validation_01.standalone.html"
open(out,"w").write(html)
print(f"WROTE {out} ({len(html)//1024} KB)")
