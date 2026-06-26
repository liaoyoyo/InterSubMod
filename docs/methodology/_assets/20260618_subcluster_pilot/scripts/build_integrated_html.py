"""[L7 整合 HTML] §13-A 由各層 JSON 注入 + base64 嵌圖 → index.standalone.html。
整合 clone/subclone 完整敘述: precedence 階層 + L0-L4 各層 data-supported 結果 + 整合 narrative + 限制 + 檔案索引。
"""
import json
import base64
import os

RPT = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/20260627_clone_subclone_integrated_report"
OUT = f"{RPT}/index.standalone.html"
TX = "#141413"; BD = "#E3DACC"; MUT = "#6B6862"; AC = "#D97757"; BLU = "#4A6E8A"; GRN = "#5B8A5B"; RED = "#C0563F"; TEAL = "#2E8B8B"


def J(name):
    p = f"{RPT}/data/{name}"
    return json.load(open(p)) if os.path.exists(p) else {}


def img(name):
    p = f"{RPT}/figures/{name}"
    if not os.path.exists(p):
        return ""
    b = base64.b64encode(open(p, "rb").read()).decode()
    return f'<img src="data:image/png;base64,{b}" style="max-width:100%;border:1px solid {BD};border-radius:6px;margin:6px 0">'


m = J("sm_locus_master_summary.json")
hp = J("sm_hp_contribution.json")
ccf = J("sm_ccf_tiers.json")
psx = J("sm_phaseset_extension.json")
me = J("sm_methyl_corroboration.json")
mr = J("sm_methyl_reextract_merged.json")
rst = J("sm_region_stats.json")
br = J("sm_branch_analysis.json")


def g(d, *keys, default="—"):
    for k in keys:
        if isinstance(d, dict) and k in d:
            d = d[k]
        else:
            return default
    return d


hp_ab = g(hp, "mutual_excl_ablation", default={})
ccf_grad = g(ccf, "ancestor_ge_descendant_VAF_gradient", default={})
ps_rel = g(psx, "ps_reliability_per_region", default={})

html = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HCC1395 clone/subclone 整合分析</title><style>
body{{margin:0;font-family:system-ui,"Noto Sans CJK TC",sans-serif;color:{TX};background:#FAF9F5;line-height:1.7}}
.wrap{{max-width:1080px;margin:0 auto;padding:0 22px 80px}}
header{{background:{TX};color:#FAF9F5;padding:26px 22px}} header h1{{margin:0;font-size:21px}} header .s{{color:#cfc9bf;font-size:13px;margin-top:5px}}
nav{{position:sticky;top:0;background:#FAF9F5;border-bottom:1px solid {BD};padding:9px 22px;font-size:12px;z-index:9}} nav a{{color:{MUT};text-decoration:none;margin-right:11px}} nav a:hover{{color:{AC}}}
h2{{font-size:18px;margin:28px 0 9px;padding-bottom:5px;border-bottom:2px solid {AC}}} h3{{font-size:15px;margin:16px 0 6px}}
.cards{{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}} .card{{flex:1;min-width:135px;background:white;border:1px solid {BD};border-radius:8px;padding:11px 13px}}
.card b{{font-size:21px;color:{AC}}} .card span{{display:block;font-size:11.5px;color:{MUT};margin-top:3px}}
table{{border-collapse:collapse;font-size:12.5px;margin:8px 0;width:100%}} th,td{{border:1px solid {BD};padding:4px 8px;text-align:center}} th{{background:#efe9df}}
.box{{background:white;border:1px solid {BD};border-left:4px solid {AC};border-radius:0 6px 6px 0;padding:11px 15px;margin:11px 0;font-size:13.5px}} .red{{border-left-color:{RED};background:#fbf0ec}} .ok{{border-left-color:{GRN};background:#eef4ee}}
code{{background:#efe9df;padding:1px 5px;border-radius:3px;font-size:11px}}
footer{{margin-top:28px;padding-top:12px;border-top:1px solid {BD};font-size:12px;color:{MUT}}}
</style></head><body>
<header><div class="wrap" style="padding:0"><h1>HCC1395 — 用 sSNV + ONT read 建構 clone/subclone：多層整合分析 <span style="background:{AC};color:#fff;font-size:11px;padding:1px 7px;border-radius:10px">⭐3</span></h1>
<div class="s">sSNV 連鎖骨幹 + HP/CCF/PS/甲基 單變量驗證整合 · Tier-R · 數字 §13-A 注入 · 已過對抗稽核</div></div></header>
<nav class="wrap" style="max-width:1080px"><a href="#prec">可信度階層</a><a href="#l0">L0 主表</a><a href="#l1">L1 HP</a><a href="#l2">L2 CCF</a><a href="#l3">L3 PS</a><a href="#l4">L4 甲基</a><a href="#narr">整合敘述</a><a href="#lim">限制</a><a href="#idx">索引</a></nav>
<div class="wrap">

<h2 id="prec">🔴 可信度先後（precedence）— 整合鐵則</h2>
<div class="box">subclone 判定鏈（順序不可顛倒）：<b>① sSNV read 連鎖（非循環·ground truth）→ ② HP（genetic·sibling vs allelic）→ ③ 甲基（epigenetic·corroborate 非 detect）</b>。衝突時 <b>genetic 勝</b>（甲基異質=同克隆內可塑性）。基礎見 <code>00b_methods_grounding_HP_ISM.md</code>。</div>

<h2 id="l0">L0 — Per-locus 主表（基底）</h2>
<div class="cards">
<div class="card"><b>{g(m,'n_loci'):,}</b><span>sSNV 主表（sum-check {('✓' if g(m,'sum_check_ok') else '✗')}）</span></div>
<div class="card"><b>{g(m,'hp_dist_linked_somatic','1-1',default=0):,}/{g(m,'hp_dist_linked_somatic','2-1',default=0):,}</b><span>HP1-1/HP2-1 linked-somatic（平衡）</span></div>
<div class="card"><b>{g(m,'tree_role_somatic','ancestor',default=0):,}/{g(m,'tree_role_somatic','descendant',default=0):,}</b><span>ancestor/descendant</span></div>
<div class="card"><b>{g(m,'tree_role_somatic','sibling',default=0):,}</b><span>sibling</span></div>
</div>
{img('01_locus_master_overview.png')}

<h2 id="l1">L1 — HP tag 貢獻（單變量 ablation）</h2>
<div class="box ok"><b>HP = sibling 判別不可或缺</b>：ablation 無 HP gate {g(hp_ab,'without_HP_gate_called_subclone',default='—'):,} 對全當 sibling；有 HP 只留 {g(hp_ab,'with_HP_gate_true_subclone(same-HP)',default='—'):,} → <b>移除 {g(hp_ab,'removed_as_allelic(diff-HP)',default='—'):,} allelic（{g(hp_ab,'false_subclone_removed_pct')}%）</b>。nested/co_linked same-HP 85–91%（~2× chance {g(hp,'chance_same_hp_baseline')}）= 真同單倍型克隆連鎖。</div>
{img('02_hp_contribution.png')}

<h2 id="l2">L2 — CCF tier + 克隆階層梯度（單變量）</h2>
<div class="box ok"><b>祖先≥後代 VAF 梯度 = {ccf_grad.get('data_support_rate','—')}（CN-clean {g(ccf_grad,'clean_only(loh+neutral)','rate')}）</b> vs null 0.5 → read-樹方向被 VAF <b>獨立驗證</b>。離散 CCF 峰 ~0.9 clonal / ~0.45 major / ~0.05 rare = 一致 subclone 層級。</div>
{img('03_ccf_tiers.png')}

<h2 id="l3">L3 — phase-set（單變量 PS）</h2>
<div class="box"><b>PS-reliable 區域 = {g(ps_rel,'single_ps_reliable',default='—'):,} / uncertain {g(ps_rel,'multi_ps_uncertain',default='—'):,}（rate {g(ps_rel,'reliable_rate')}）</b>。多 PS region = phase-switch 風險，HP 判別不可信（最乾淨 sibling 結論排除之）。Tier-PS（同 PS >50kb）= germline 單倍型 context，<b>非克隆連鎖</b>（同 PS ≠ 同克隆）。{img('04_phaseset.png')}</div>

<h2 id="l4">L4 — 甲基 corroboration（單變量；已從 BAM 重抽補完 genome-wide）</h2>
<div class="box red">既有 ISM 甲基輸出只覆蓋 0.19%（9/4678）→ <b>已直接從 tumor BAM MM/ML 重抽（驗證 corr=1.000 vs ISM）</b>。genome-wide CN-clean：<b>{g(mr,'n_tested',default='—')} 區測試，{g(mr,'n_corroborated',default='—')}（{g(mr,'corroboration_rate')}）甲基 corroborate 遺傳群，全部 subclone-specific（cis-ASM-explained {g(mr,'cis_explained_rate')}）</b>。median sig CpG=0 → <b>93% 區域甲基不區分遺傳 subclone</b>。🔴 甲基 = <b>弱 corroborator 非 detector</b>（6.6% 區域有獨立表觀支持；其餘無）。</div>
{img('05_methylation_corroboration.png')}

<h2 id="narr">整合敘述（sSNV + read → clone/subclone）</h2>
<div class="box"><b>骨幹</b>：ONT long read 同分子讀多 sSNV → 共現 2×2 → 局部克隆樹（7,143 區域，full_tree {g(rst,'tree_shape','full_tree',default='—')}）。<br>
<b>HP 加值</b>：把 read 級互斥轉成細胞層 sibling（移除 {g(hp_ab,'false_subclone_removed_pct')}% allelic）。<br>
<b>CCF 加值</b>：祖先≥後代 VAF {ccf_grad.get('data_support_rate','—')} 獨立驗證樹方向 + 離散 subclone CCF 層級。<br>
<b>PS</b>：reliability 旗標（排除 phase-uncertain 區）；不延伸克隆連鎖。<br>
<b>甲基</b>：已從 BAM 重抽 genome-wide（{g(mr,'n_tested',default='—')} 區）→ 僅 {g(mr,'corroboration_rate')} 區域 corroborate（弱 corroborator，全 subclone-specific）。<br>
→ <b>結論</b>：sSNV+HP+CCF 三軸一致支撐局部克隆階層（單樣本 ⭐3 分子證據）；甲基為**弱**佐證（6.6% 區域有獨立表觀支持）。</div>

<h2 id="lim">🔴 限制（誠實，無可被質疑的未佐證推論）</h2>
<div class="box red">⭐3 單樣本；regional（≤read-span）<b>非 genome-wide tree</b>；HP 有 ~85-90% phasing 誤差 + problem PS block（已標記排除）；甲基 corroborate 非 detect（既有輸出覆蓋 0.19% → 已重抽 740 區，僅 6.6% 弱 corroborate）；CCF 僅 CN-clean 可估（gain multiplicity 歧義）；64% sSNV 在 CN-gain 混淆（乾淨集=LOH/neutral）；分子證據非 single-cell confirmation。對外勿稱「甲基偵測 subclone / genome-wide tree / 對手缺檢定」。</div>

<h2 id="idx">檔案索引（可驗證 + 可查詢）</h2>
<table><tr><th>檔</th><th>內容</th></tr>
<tr><td><code>00_INDEX.md</code></td><td>定義/格式/索引</td></tr>
<tr><td><code>00b_methods_grounding_HP_ISM.md</code></td><td>HP/ISM 定義+precedence+錯誤模式</td></tr>
<tr><td><code>01..07_*.md</code></td><td>L0-L6 各層觀察+解釋+分析</td></tr>
<tr><td><code>06_integrated_narrative.md</code></td><td>整合敘述（對抗稽核後）</td></tr>
<tr><td><code>data/*.json,*.tsv</code></td><td>各層數據（可 grep 重算）</td></tr>
<tr><td><code>figures/*.png</code></td><td>各層圖</td></tr></table>

<footer>HCC1395 ⭐3 · 多層單變量整合 · §13-A 由 data/*.json 注入 · 腳本 build_integrated_html.py · 對抗稽核見 06_*.md</footer>
</div></body></html>"""
open(OUT, "w").write(html)
print(f"[-> {OUT}] {len(html)} bytes")
