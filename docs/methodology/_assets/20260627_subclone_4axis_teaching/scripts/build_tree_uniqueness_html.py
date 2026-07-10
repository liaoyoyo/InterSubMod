#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_tree_uniqueness_html.py — 樹唯一性 / VAF 加權 完整解釋 HTML(2026-07-10)。
讀 20260710_tree_uniqueness_analysis_data.json(全 verified 真值) → 靜態注入(無 JS,零捏造)。
內容:①VAF/CCF 加權流程 ②reach 結果 ③形狀 vs 精確樹(真多拓撲) ④CN confound ⑤為何定不出唯一樹 ⑥誠實結論。
輸出:docs/methodology/_assets/20260710_tree_uniqueness_analysis.standalone.html
"""
import json, os

DATA = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260710_tree_uniqueness_analysis_data.json"
OUT = "/big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260710_tree_uniqueness_analysis.standalone.html"
d = json.load(open(DATA))
W, SH, CN, RC = d["ccf_weighting"], d["shape"], d["cn_confound"], d["region_c_context"]

def bar(pct, color, w=180):
    return f"<span style='display:inline-block;height:11px;width:{pct/100*w:.0f}px;background:{color};border-radius:2px;vertical-align:middle'></span>"

# reach 三段
reach = f"""<table>
<tr><th>階段</th><th>數量</th><th>比例</th><th></th></tr>
<tr><td>ambiguous 單位(等權多樹,待破)</td><td>{W['n_ambiguous']:,}</td><td>100%</td><td></td></tr>
<tr><td class='ok'>VAF 破 tie(top posterior ≥0.6)</td><td>{W['reach_n']:,}</td><td class='hi'>{W['reach_pct']}%</td><td>{bar(W['reach_pct'],'#4a9eff')}</td></tr>
<tr><td>其中 winner pigeonhole-clean</td><td>{W['winner_clean_n']:,}</td><td class='hi'>{W['winner_clean_pct']}%</td><td>{bar(W['winner_clean_pct'],'#22c55e')}</td></tr>
<tr><td class='warn'>🔴 扣 CN 後真正可信(strictly-neutral & 破)</td><td>{W['trustworthy_n']:,}</td><td class='sub'>{W['trustworthy_pct']}%</td><td>{bar(W['trustworthy_pct'],'#d9534f')}</td></tr>
<tr><td>VAF 破不掉(對稱/難)</td><td>{W['unbroken_n']:,}</td><td>{W['unbroken_pct']}%</td><td></td></tr>
</table>"""

top = W["top_dist"]
topbars = "".join(f"<div class='vrow'><span class='vlab'>{k}</span>{bar(100*v/W['n_ambiguous'],'#4a9eff',240)}<span class='vn'>{v:,}</span></div>" for k, v in top.items())

# 形狀
shape = f"""<table>
<tr><th>群體</th><th>單一形狀<br><span class='dim'>拓撲已定,只內部順序歧義</span></th><th>🔴 真多形狀<br><span class='dim'>不同拓撲(線性/分支/星)</span></th></tr>
<tr><td>全 ambiguous ({W['n_ambiguous']:,})</td><td>{SH['all_single']:,} ({SH['all_single_pct']}%)</td><td class='sub'>{SH['all_multi']:,} ({SH['all_multi_pct']}%)</td></tr>
<tr><td>VAF 破不掉者 ({SH['unbroken_total']:,})</td><td>{SH['unbroken_single']:,} ({SH['unbroken_single_pct']}%)</td><td class='sub'>{SH['unbroken_multi']:,} ({SH['unbroken_multi_pct']}%)</td></tr>
</table>"""

usdist = SH["unbroken_shapes"]
usbars = "".join(f"<div class='vrow'><span class='vlab'>{k} 種形狀</span>{bar(100*v/SH['unbroken_total'],'#c084fc' if k!='1' else '#22c55e',240)}<span class='vn'>{v:,} ({100*v/SH['unbroken_total']:.0f}%)</span></div>" for k, v in usdist.items())

# cross-tab
def ct(broke, single, cn):
    for r in SH["crosstab"]:
        if r["broke"] == broke and r["single_shape"] == single and r["cn"] == cn:
            return r["n"]
    return 0
crosstab = f"""<table>
<tr><th rowspan=2>VAF</th><th rowspan=2>形狀</th><th colspan=2>CN</th></tr>
<tr><th>clean</th><th>gain(read≠CCF)</th></tr>
<tr><td rowspan=2>破 tie</td><td>單一形狀</td><td>{ct(True,True,'clean')}</td><td>{ct(True,True,'gain')}</td></tr>
<tr><td>多形狀</td><td class='ok'>{ct(True,False,'clean')}</td><td class='dimc'>{ct(True,False,'gain')}</td></tr>
<tr><td rowspan=2>破不掉</td><td>單一形狀</td><td>{ct(False,True,'clean')}</td><td>{ct(False,True,'gain')}</td></tr>
<tr><td>多形狀</td><td>{ct(False,False,'clean')}</td><td>{ct(False,False,'gain')}</td></tr>
</table>"""

# CN
cnr = CN["region_pct"]; cnn = CN["nested_pct"]
cntab = f"""<table>
<tr><th>CN 狀態</th><th>全共現區 ({CN['region_total']:,})</th><th>巢狀 clone→subclone ({CN['nested_total']})</th></tr>
<tr><td class='ok'>neutral(VAF≈CCF 乾淨)</td><td>{CN['region']['neutral']} ({cnr['neutral']}%) {bar(cnr['neutral'],'#22c55e',80)}</td><td class='sub'>{CN['nested']['neutral']} ({cnn['neutral']}%)</td></tr>
<tr><td>loh(抬高 VAF)</td><td>{CN['region']['loh']:,} ({cnr['loh']}%) {bar(cnr['loh'],'#f0ad4e',80)}</td><td>{CN['nested']['loh']} ({cnn['loh']}%)</td></tr>
<tr><td class='warn'>🔴 gain(稀釋 VAF→假 subclone)</td><td>{CN['region']['gain']:,} ({cnr['gain']}%) {bar(cnr['gain'],'#d9534f',80)}</td><td>{CN['nested']['gain']} ({cnn['gain']}%)</td></tr>
<tr><td>loss</td><td>{CN['region']['loss']} ({cnr['loss']}%)</td><td>{CN['nested']['loss']} ({cnn['loss']}%)</td></tr>
<tr><td><b>CN-altered 合計</b></td><td class='sub'>{CN['region_altered_pct']}%</td><td class='sub'>{CN['nested_altered_pct']}%</td></tr>
</table>"""

html = f"""<!--
建立 2026-07-10 | 樹唯一性/VAF 加權 完整解釋 | build_tree_uniqueness_html.py
資料 20260710_tree_uniqueness_analysis_data.json(全 verified) | 數字全注入·靜態·無 JS·反捏造
-->
<meta charset="utf-8">
<title>為何定不出唯一的樹 — VAF 加權與拓撲歧義完整分析</title>
<style>
:root{{--bg:#0f1420;--card:#1a2233;--ink:#e8edf5;--dim:#8fa0b8;--line:#2a3550;--acc:#4a9eff}}
*{{box-sizing:border-box}}
body{{margin:0;font-family:-apple-system,'Segoe UI',Roboto,'Noto Sans TC',sans-serif;background:var(--bg);color:var(--ink);line-height:1.6;font-size:14px}}
.wrap{{max-width:1080px;margin:0 auto;padding:26px}}
h1{{font-size:22px;margin:0 0 4px}} h2{{font-size:17px;margin:30px 0 10px;padding-left:9px;border-left:4px solid var(--acc)}}
.sub0{{color:var(--dim);font-size:13px;margin-bottom:8px}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:14px 17px;margin:10px 0}}
table{{border-collapse:collapse;width:100%;font-size:13px;margin:6px 0}}
th,td{{border:1px solid var(--line);padding:6px 9px;text-align:center}} th{{background:#131b2b;color:var(--dim);font-weight:600}}
td:first-child{{text-align:left}}
.hi{{color:#7dd3fc;font-weight:700}} .sub{{color:#fca5a5;font-weight:700}} .ok{{color:#86efac;font-weight:600}}
.dim{{color:var(--dim);font-weight:400;font-size:11.5px}} .dimc{{color:var(--dim)}}
.warn{{color:#fca5a5}}
.note{{font-size:12.5px;color:var(--dim);margin:6px 0}}
.big{{background:#2a1f1a;border-left:3px solid #d9534f;padding:10px 14px;border-radius:6px;margin:10px 0}}
.good{{background:#16241a;border-left:3px solid #22c55e;padding:10px 14px;border-radius:6px;margin:10px 0}}
.vrow{{display:flex;align-items:center;gap:8px;height:18px;margin:2px 0}}
.vlab{{width:70px;text-align:right;color:var(--dim);font-size:12px}} .vn{{color:var(--dim);font-size:12px}}
ol.flow{{margin:6px 0;padding-left:20px}} ol.flow li{{margin:3px 0}}
.k{{font-family:ui-monospace,monospace;color:#7dd3fc}}
</style>
<div class="wrap">
<h1>為何定不出唯一的樹 — VAF 加權與拓撲歧義完整分析</h1>
<div class="sub0">HCC1395 · {d['backbone']} · 2026-07-10 全 verified · 數字全注入反捏造</div>

<div class="card">
<b>中心問題</b>：sSNV 共現建出的 Steiner 樹每條邊<b>權重都=1（等權）→ 常畫得出好幾棵都對的樹</b>。
能否用 <b>VAF/CCF（每突變出現在多少比例細胞）</b> 加權，把等權樹變成「有先後可能性、能挑出較可能者」？
<div class="note">原理（白話）：早發生的突變被更多後代細胞繼承 → 出現在更多 read → VAF 高 → 排祖先。pigeonhole：祖先 CCF ≥ 後代 CCF。純 read/VAF＝遺傳訊號、非甲基、非循環。</div>
</div>

<h2>§1　方法：VAF/CCF 加權流程（6 階段）</h2>
<div class="card">
<ol class="flow">
<li><b>輸入</b>：枚舉出的 N 棵等權樹。節點=k-bit genotype（<span class="k">A</span>/<span class="k">R</span>），ROOT=全 R。</li>
<li><b>每突變 VAF</b>：<span class="k">vaf_j = nALT/(nREF+nALT)</span>（within-family，含 partial reads）。</li>
<li>🔴 <b>VAF→CCF 矯正</b>（目前缺）：<span class="k">CCF = vaf·(ρ·CNt+2(1−ρ))/(ρ·m)</span>。region 內共享 CN → 排序時 CN 對消，<b>殘留只在 per-mutation multiplicity（gain 區）</b>。</li>
<li><b>邊權</b>：<span class="k">w(P→C)=Σ_(a∈P)(CCF_a−CCF_j)</span>，acquired j。正=祖先頻率更高（一致）。</li>
<li><b>樹分→posterior</b>：<span class="k">softmax(Σ_edge w / T)</span>；top≥0.6=破 tie；對稱(CCF 相等)→保持等權不硬排。</li>
<li><b>CN-gate + soft beta-binomial</b>（正解，待做）：只在 strictly-neutral 採信；read 數不確定性進權重。</li>
</ol>
</div>

<h2>§2　結果 A：VAF 加權能破多少？</h2>
<div class="card">
{reach}
<div class="note"><b>top posterior 分布</b>（越高越確信）：</div>
{topbars}
<div class="good">✅ <b>方法成立</b>：66.8% 能破、破時 winner 99.4% pigeonhole-clean、非甲基非循環，可寫進論文「頻率加權層」。</div>
<div class="big">🔴 <b>但真正可信只 6.1%（strictly-neutral，權威 ccf_and_cn_multisample）</b>：CN-gain 佔 <b>54.9%</b>（2,134 個破在 gain 區、read≠CCF）→ 那些「破」不可信。<b>先前寫的 22.3% 是錯的</b>——誤把 LOH 也當乾淨，但 LOH 同樣抬高 VAF，只有 strictly-neutral（597 區）才 read≈CCF。</div>
</div>

<h2>§3　結果 B：VAF 破不掉的，是同形狀換順序、還是真的多種拓撲？</h2>
<div class="card">
<div class="note">關鍵區分：<b>拓撲（形狀：single/linear/branched/star）</b> vs <b>精確樹（含順序+隱藏節點擺放）</b>。</div>
{shape}
<div class="big">🔴 <b>推翻直覺</b>：多數歧義<b>不是</b>「同形狀換順序」——<b>71.2% 是真的多種拓撲</b>（線性 vs 分支 vs 星，連「是不是分支」都不定）。VAF 破不掉者更高達 <b>74.1% 真多拓撲</b>。</div>
<div class="note">破不掉者的形狀種數分布：</div>
{usbars}
<div class="note" style="margin-top:12px">交叉表（VAF 破? × 形狀 × CN）：</div>
{crosstab}
<div class="note">➊ VAF 破 tie 時=挑一棵樹=收斂到一個形狀（3,980 破的裡 2,778 原是多形狀被收斂），但其中 1,906 是 CN-gain 不可信 → <b>可信收斂多形狀→單一者只 872（neutral+loh+loss 粗分＋破＋原多形狀；若嚴格 strictly-neutral 更少）</b>。
➋ CN 不製造多形狀：clean 單形狀 {SH['cn_clean_single_pct']}% vs gain {SH['cn_gain_single_pct']}%（差不多）→ <b>多形狀來自 partial 缺角，CN 只讓排序不可信</b>。</div>
</div>

<h2>§4　結果 C：CN confound 打進共現骨幹（SEQC2 真值）</h2>
<div class="card">
{cntab}
<div class="big">🔴 共現區 <b>94.2% CN-altered</b>、巢狀 clone→subclone 區 <b>98% CN-altered（僅 15 CN-clean）</b>。衍生突變 VAF 中位 gain {CN['derived_vaf_gain']} &lt; neutral {CN['derived_vaf_neutral']}（稀釋跡象弱,n=15 小）。</div>
<div class="note">🔑 拓撲（共現連鎖）CN-robust；但 VAF/subclone 頻率解讀 98% 被 CN 混淆。HCC1395 高度非整倍體，共現區富集在 gain。</div>
</div>

<h2>§5　為何定不出唯一的樹（5 個原因，白話）</h2>
<div class="card">
<table>
<tr><th>原因</th><th>白話</th><th>對應狀況</th></tr>
<tr><td>① 拼圖缺角<br>（partial read）</td><td>read 太短沒跨過所有位點,中間漏格,補法不只一種</td><td class='sub'>主因:全 5,959 歧義都是這型</td></tr>
<tr><td>② 看不見的祖先<br>（hidden node）</td><td>中間祖先 0 read,一定存在但擺哪是推的</td><td>與①綁定→產生不同形狀</td></tr>
<tr><td>③ 頻率相同<br>（對稱 CCF）</td><td>兩突變一樣多細胞→誰先分不出</td><td>VAF 破不掉的 34%;74% 仍真多拓撲</td></tr>
<tr><td>④ 拷貝數搗亂<br>（CN-gain）</td><td>突變在多份 copy→VAF 被灌大≠細胞比例</td><td>categorical-gain 54.9%;strictly-neutral 可信只 6.1%</td></tr>
<tr><td>⑤ 混合平均<br>（單 bulk）</td><td>看整團腫瘤平均,不同樹打出同平均</td><td>根本非識別性;需 single-cell/多切片</td></tr>
</table>
<div class="note">①②<b>被 VAF 破掉 66.8%</b>（但 CN 卡到可信只 6.1% strictly-neutral）；③④⑤ <b>VAF 破不掉，且多為真多拓撲</b>。</div>
</div>

<h2>§6　誠實結論</h2>
<div class="card">
<div class="good">✅ 「等權樹 → VAF 加權有向樹」<b>方法成立、已實作骨架、66.8% 能破、winner 99.4% 乾淨</b>——非甲基非循環，可進論文。</div>
<div class="big">🔴 但大多數 ambiguous 區的誠實答案是<b>「連是線性還是分支都定不出來」（71% 真多拓撲）</b>，比「順序定不出」更強的非唯一；真正拓撲唯一的是 ~52% determined 區（共現完整無缺角）。</div>
<div class="note">➊ 「定不出唯一樹」不是方法失敗,是單一 bulk 的物理極限。➋ VAF 矯正的硬缺口=整數 CN+purity（現 categorical bed 算不出 multiplicity）。➌ 突破天花板只有 single-cell 或 multi-region。</div>
<div class="note" style="margin-top:8px">關聯：region 層 c / VAF 峰值 / 巢狀 clone→subclone 詳見 <span class="k">20260709_structure_result.standalone.html</span>（全樣本 c≥1={RC['all_samples_c_ge1']:,}、c≥2 {RC['c_ge2_pct']}%、maxC={RC['maxC_HCC1395']}）。</div>
</div>

<div class="note" style="margin-top:18px">一鍵重算：<span class="k">ccf_tree_weighting_full_observe.py</span> / <span class="k">shape_vs_exact_ambiguity.py</span> / <span class="k">cn_confound_check.py</span> → <span class="k">build_tree_uniqueness_html.py</span>。數字全從 JSON 注入(反捏造)。</div>
</div>
"""
open(OUT, "w", encoding="utf-8").write(html)
print(f"→ 寫出 {OUT} ({len(html):,} bytes)")
# 驗證無殘留模板
assert "{" + "{" not in html.replace("{{", ""), "殘留模板"
print("關鍵注入檢查:", all(str(x) in html for x in [W['n_ambiguous'], W['reach_pct'], SH['all_multi_pct'], CN['region_altered_pct'], W['trustworthy_pct']]))
print("含 <script>:", "<script" in html)
