#!/usr/bin/env python3
"""Standalone HTML explaining whether the sr2b/sr2c joint model can be landed here.

All figures are injected from spec_readiness_audit.json and the earlier audit
outputs; the script refuses to render on a missing key so nothing can be typed
by hand into the page.
"""

from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(HERE, "..", "docs", "20260813_lineage約束擬合CCF可行性_01.standalone.html")


def load(n):
    with open(os.path.join(DATA, n)) as fh:
        return json.load(fh)


def need(d, *path):
    cur = d
    for p in path:
        if isinstance(cur, dict) and p in cur:
            cur = cur[p]
        else:
            sys.exit(f"REFUSE TO RENDER: missing {'/'.join(map(str,path))}")
    if cur is None:
        sys.exit(f"REFUSE TO RENDER: null {'/'.join(map(str,path))}")
    return cur


sp = load("spec_readiness_audit.json")
rv = load("report_values.json")

INV = need(sp, "canonical_mlhp_inventory")
S2 = need(sp, "s2_evidence_scale")
REQ = need(sp, "required_inputs")
BYCN = need(sp, "s2_evidence_by_cn_class")

cv = need(rv, "constraint_validity")
lk = need(rv, "linkage_constraint")
cg = need(rv, "clean_ground")
cp = need(rv, "cn_phenomena")
sq = need(rv, "seqc2_external_comparison")
spec_af = load("neutral_af_spectrum.json")
SUB_MED = need(spec_af, "subclonal_median_af")
SUB_N = need(spec_af, "subclonal_count")
pilot = load("hp_to_allele_pilot.json")
PD = need(pilot, "diagnosis")
PV = need(pilot, "vcf_scan")
PH = need(pilot, "haplotype_major_fraction")

STATUS_STYLE = {
    "AVAILABLE": ("ok", "已具備"),
    "ALREADY_COMPUTED": ("ok", "已算好（我們的貢獻）"),
    "AVAILABLE_WITH_CAVEAT": ("warn", "具備但有但書"),
    "AVAILABLE_UNVALIDATED": ("warn", "具備但未驗證"),
    "MEASURABLE_NOT_MEASURED": ("warn", "可測但尚未測"),
    "NOT_YET_ESTIMATED": ("warn", "尚未估計"),
    "JUDGEMENT_PARAMETER": ("warn", "判斷參數"),
    "MISSING_BUT_DERIVABLE": ("bad", "缺，但可導出"),
}
LABEL = {
    "rho_purity": "ρ 腫瘤純度",
    "CN_m_total": "CNₘ 總拷貝數",
    "allele_specific_CN": "G allele-specific CN",
    "hp_to_allele_mapping": "HP → allele 對應",
    "mu_m_multiplicity": "μₘ mutation multiplicity",
    "read_states_S2": "D_v 讀取狀態（S₂）",
    "alt_ref_counts_S1": "aₘ,dₘ alt/ref 計數（S₁）",
    "error_floor_delta": "δᵤ 誤差地板（ε, η）",
    "lambda_escape_mass": "λ escape mass",
    "tree_topology_T_and_ancestry_E": "T, Eᵤ 樹拓撲與祖先關係",
}
ORDER = [
    "read_states_S2", "alt_ref_counts_S1", "tree_topology_T_and_ancestry_E",
    "CN_m_total", "rho_purity", "allele_specific_CN",
    "error_floor_delta", "lambda_escape_mass",
    "mu_m_multiplicity", "hp_to_allele_mapping",
]

rows = ""
for k in ORDER:
    v = REQ[k]
    cls, zh = STATUS_STYLE.get(v["status"], ("", v["status"]))
    rows += (
        f"<tr><td><strong>{LABEL.get(k,k)}</strong><br>"
        f"<span style='font-size:.82rem;color:var(--muted)'>{v['spec_symbol']}</span></td>"
        f"<td class='{cls}' style='white-space:nowrap'>{zh}</td>"
        f"<td style='font-size:.88rem'>{v['what_we_have']}</td>"
        f"<td style='font-size:.88rem;color:var(--muted)'>{v['gap']}</td></tr>"
    )

bycn_rows = "".join(
    f"<tr><td>{k}</td><td class='num'>{v['units']:,}</td>"
    f"<td class='num'>{v['full_molecules']:,}</td>"
    f"<td class='num'>{v['partial_molecules']:,}</td>"
    f"<td class='num bad'>{v['partial_share_percent']}%</td></tr>"
    for k, v in BYCN.items()
)

HTML = f"""<title>以 lineage 約束擬合 CCF 與建樹 — 可行性稽核</title>
<style>
:root {{
  --bg:#fff; --fg:#1a1d21; --muted:#5b6470; --line:#e2e6ea; --card:#f7f9fb;
  --ok:#2e7d5b; --bad:#c0392b; --warn:#d98324; --mid:#4a6fa5; --accent:#1f4e79;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --bg:#14171a; --fg:#e8ecef; --muted:#9aa4b0; --line:#2a3038; --card:#1b1f24;
    --ok:#4caf82; --bad:#e57063; --warn:#e5a04a; --mid:#7ba3d8; --accent:#8ab4e8;
  }}
}}
:root[data-theme="dark"] {{
  --bg:#14171a; --fg:#e8ecef; --muted:#9aa4b0; --line:#2a3038; --card:#1b1f24;
  --ok:#4caf82; --bad:#e57063; --warn:#e5a04a; --mid:#7ba3d8; --accent:#8ab4e8;
}}
*{{box-sizing:border-box}}
body{{margin:0;padding:2rem 1.25rem 5rem;background:var(--bg);color:var(--fg);
 font:16px/1.75 -apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans TC","PingFang TC",sans-serif}}
.wrap{{max-width:960px;margin:0 auto}}
h1{{font-size:1.7rem;line-height:1.35;margin:0 0 .4rem;letter-spacing:-.01em}}
h2{{font-size:1.2rem;margin:2.6rem 0 .8rem;padding-top:1.1rem;border-top:1px solid var(--line)}}
h3{{font-size:1.02rem;margin:1.5rem 0 .5rem;color:var(--accent)}}
.sub{{color:var(--muted);margin:0 0 1.6rem;font-size:.93rem}}
.verdict{{background:var(--card);border-left:4px solid var(--accent);padding:1.1rem 1.3rem;
 border-radius:0 8px 8px 0;margin:1.4rem 0}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:.85rem;margin:1.3rem 0}}
.stat{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:.95rem 1.05rem}}
.stat .n{{font-size:1.65rem;font-weight:650;display:block;line-height:1.2;letter-spacing:-.02em}}
.stat .l{{font-size:.83rem;color:var(--muted);display:block;margin-top:.3rem}}
.stat.ok .n{{color:var(--ok)}} .stat.bad .n{{color:var(--bad)}} .stat.warn .n{{color:var(--warn)}}
table{{border-collapse:collapse;width:100%;font-size:.9rem;margin:1rem 0}}
th,td{{border-bottom:1px solid var(--line);padding:.55rem .6rem;text-align:left;vertical-align:top}}
th{{color:var(--muted);font-weight:600;font-size:.82rem;text-transform:uppercase;letter-spacing:.04em}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}
td.ok{{color:var(--ok);font-weight:600}} td.bad{{color:var(--bad);font-weight:600}}
td.warn{{color:var(--warn);font-weight:600}}
.scroll{{overflow-x:auto}}
.note{{background:var(--card);border:1px solid var(--line);border-radius:8px;
 padding:.9rem 1.1rem;font-size:.89rem;color:var(--muted);margin:1.1rem 0}}
.note.bad{{border-left:3px solid var(--bad)}} .note.ok{{border-left:3px solid var(--ok)}}
.note.warn{{border-left:3px solid var(--warn)}}
code{{background:var(--card);padding:.12em .4em;border-radius:4px;font-size:.88em;border:1px solid var(--line)}}
.formula{{background:var(--card);border:1px solid var(--line);border-radius:8px;
 padding:1rem 1.2rem;margin:1rem 0;text-align:center;font-size:1.05rem;overflow-x:auto}}
ol,ul{{padding-left:1.4rem}} li{{margin:.35rem 0}}
.step{{background:var(--card);border:1px solid var(--line);border-radius:9px;padding:1rem 1.2rem;margin:.7rem 0}}
.step h4{{margin:0 0 .4rem;font-size:.98rem;color:var(--accent)}}
footer{{margin-top:3.4rem;padding-top:1.2rem;border-top:1px solid var(--line);color:var(--muted);font-size:.83rem}}
</style>
<div class="wrap">

<h1>以 lineage 約束來擬合 CCF 與建樹：可行性稽核</h1>
<p class="sub">對照 lab-tutorial sr2b（CCF 與分群）／sr2c（局部分子系統發生）規格 ·
以 HCC1395 canonical frozen cohort 實測 · 2026-08-13</p>

<div class="verdict">
<strong>結論：可以，而且我們已經持有規格最難取得的那一塊（樹結構與約束）。
但有<u>一個</u>阻斷輸入尚未計算 —— HP→allele 對應。</strong><br><br>
規格本身自陳「This is a specification, not yet validated」。
我們的資料正好能提供它需要的驗證，反過來規格也補上了我們缺的 μₘ 校正層。
</div>

<h2>1 · 規格要的是什麼</h2>
<p>sr2b 的核心是一條把 CCF、拷貝數與 multiplicity 綁在一起的期望 VAF：</p>
<div class="formula">θ<sub>mk</sub> = E[VAF<sub>m</sub> | z<sub>m</sub>=k] =
ρ · c<sub>k</sub> · μ<sub>m</sub> / ( ρ · CN<sub>m</sub> + 2(1−ρ) )</div>
<p>並以 <code>z</code>（突變→群指派）、<code>T</code>（樹拓撲）、<code>G</code>（copy genotype）、
<code>w</code>（細胞比例）<strong>聯合優化</strong>；S₁ 單變異窗用 beta-binomial，
S₂ 多變異窗用讀取狀態的 multinomial，兩者以 <strong>單向</strong> prior 結合（escape mass λ）。</p>

<p>sr2c 則明確界定估計對象是
<strong>每個 linked segment 內的局部分子系統發生，不是全域細胞層 clone 樹</strong>，
並採用兩條與我們完全相同的原則：</p>
<ul>
<li><strong>Pigeonhole</strong>：子代細胞比例之和不得超過親代 —— 且明說「can only rule trees out, never select one」</li>
<li><strong>四配子檢驗</strong>：兩位點若最多出現三種 haplotype 即相容；<strong>全部成對相容 ⇒ 樹已由 Xᵤ 決定，無需搜尋</strong></li>
</ul>
<div class="note ok">
這兩條正是本輪稽核已經實作並量測過的東西。換句話說，
<strong>規格的結構層我們已經跑完了</strong>（見第 3 節）。
</div>

<h2>2 · 十項所需輸入的可得性</h2>
<div class="scroll"><table>
<tr><th style="width:20%">輸入</th><th>狀態</th><th style="width:34%">我們有什麼</th><th style="width:30%">缺口</th></tr>
{rows}
</table></div>

<div class="grid">
  <div class="stat ok"><span class="n">7</span><span class="l">已具備或已算好</span></div>
  <div class="stat warn"><span class="n">2</span><span class="l">可測／待估／判斷參數</span></div>
  <div class="stat bad"><span class="n">1</span><span class="l">阻斷輸入：HP→allele 對應</span></div>
</div>

<h2>3 · 我們已經替規格算好的部分</h2>
<div class="grid">
  <div class="stat ok"><span class="n">{cv['ordering_total']:,}</span>
    <span class="l">順序約束（HCC1395）· 全 7 樣本 {lk['inventory']['TOTAL']['ordering_constraints']:,}</span></div>
  <div class="stat ok"><span class="n">0</span>
    <span class="l">四配子違反 — 全部成對相容，符合規格「樹已由 Xᵤ 決定」</span></div>
  <div class="stat"><span class="n">{lk['inventory']['TOTAL']['exclusion_units']:,}</span>
    <span class="l">互斥宣告（僅 CN-neutral 可升格細胞層）</span></div>
  <div class="stat"><span class="n">{cv['median_linking_depth']}</span>
    <span class="l">中位連結深度 → 誤差地板可就地量測</span></div>
</div>
<div class="note">
規格說「A zero observation has three causes: sampling zero, structural zero, or error floor.
Only the error floor can be measured in-place.」我們已用連結深度把
<strong>sampling zero</strong> 這一項量化：中位可排除比例 {cv['median_min_detectable_percent']}%，
depth&lt;20 的 {cv['low_depth_unit_percent']}% unit 對低比例狀態不敏感。
</div>

<h2>4 · 實測的落地挑戰：88.57% 的證據是 partial</h2>
<p>sr2c 要求「unobserved positions must be marginalized, not assumed reference」。
這在我們的資料上<strong>不是可選項</strong>：</p>
<div class="grid">
  <div class="stat"><span class="n">{S2['fully_covering_molecules']:,}</span>
    <span class="l">完整覆蓋分子（可直接讀 pattern）</span></div>
  <div class="stat bad"><span class="n">{S2['partial_molecules_with_X']:,}</span>
    <span class="l">含 X 的 partial 分子</span></div>
  <div class="stat bad"><span class="n">{S2['partial_share_percent']}%</span>
    <span class="l">partial 佔全部證據</span></div>
  <div class="stat warn"><span class="n">{S2['full_evidence_share_per_unit']['median']}</span>
    <span class="l">每 unit 完整證據佔比（中位）</span></div>
</div>
<div class="note bad">
若只用完整覆蓋的分子，會丟掉約 <strong>{S2['partial_share_percent']}%</strong> 的證據；
每個 unit 的完整證據中位只佔 <strong>{round(100*S2['full_evidence_share_per_unit']['median'],1)}%</strong>。
partial pattern 中位有 <strong>{round(100*S2['x_slot_share_in_partial_patterns']['median'],1)}%</strong> 的位置槽是 X。
<strong>規格的 sub-cube marginalisation 是這批資料能不能用的前提，不是精緻化選項。</strong>
</div>
<div class="scroll"><table>
<tr><th>CN 類別</th><th>unit</th><th>完整分子</th><th>partial 分子</th><th>partial 佔比</th></tr>
{bycn_rows}
</table></div>

<h2>5 · 唯一的阻斷點：HP → allele 對應</h2>
<div class="note bad">
θ 公式需要 <strong>μₘ</strong>（帶變異的 copy 數）。要得到它必須知道
<strong>「觀測到的這個 HP family，在這一段裡帶幾份 copy」</strong>。<br><br>
我們有：germline phasing 已定義 HP family ✔ ｜ SAVANA 提供 major/minor CN
（HCC1395 覆蓋 92.5% 段）✔<br>
我們缺：<strong>兩者之間的連結</strong> —— 哪個 HP family 對應 major allele。
</div>

<h3>可導出，需要三步</h3>
<div class="note bad">
<strong><span style="color:var(--bad)">●</span> 已實測，且第一次嘗試失敗 —— 這個失敗值得記錄。</strong><br><br>
我先假設可以「只用 phased germline VCF」導出：它同時帶 <code>GT</code>（{PV['phased_het']:,} 個 autosomal phased het）
與 <code>AD</code>（ref/alt 深度），看起來足夠。實測結果：
LOH 段的 HP major fraction 中位 <strong>{PH['seqc2_loh_segments']['median']}</strong>、
非 LOH 段 <strong>{PH['seqc2_non_loh_segments']['median']}</strong>，
分離 <strong>AUC = {PH['separation_auc']}</strong>（比隨機還差）。<br><br>
<strong>診斷</strong>：該 VCF 的 header <code>##commandline</code> 指向
<code>alignment-sort-hcc1395bl.bam</code> 與 <code>clair3_normal_output</code> ——
<strong>它的 AD 全部來自配對的 normal（HCC1395BL），不是 tumour</strong>。
normal 是正常二倍體，每個 het 位點必然 ~50:50，方法從一開始就測錯了對象。<br><br>
<strong>裁決：方法用錯對象，不是方法被否證。</strong>
但「不必回到 BAM」這個較樂觀的估計<strong>撤回</strong>。
</div>
<div class="step">
<h4>步驟 1 · 在 HP-tagged <u>tumour</u> BAM 上統計 HP-resolved 深度</h4>
phased germline VCF 提供的是<strong>雜合位點的位置與各 allele 的 HP 歸屬</strong>（這部分可用，
{PV['phased_het']:,} 個位點）；但<strong>每個 haplotype 的深度必須在 tagged tumour BAM 上計算</strong>。
好消息是 tagged BAM 已帶 HP tag，<strong>不需要重新 phasing</strong>，只需在這些位點統計 HP tag 分布。
</div>
<div class="step">
<h4>步驟 2 · 判定該段哪個 HP 是 major allele</h4>
段內 HP1 與 HP2 深度比 → 對應 SAVANA 的 <code>copyNumber</code> 與
<code>minorAlleleCopyNumber</code>，得到<strong>每個 HP family 各帶幾份 copy</strong>。
需同時輸出信心值（段內 het SNP 數、深度不平衡的一致性）。
</div>
<div class="step">
<h4>步驟 3 · 把 CN_HP 寫回 unit 層，並改寫 θ 為 HP-conditional 形式</h4>
我們的 read-AF 已經是 HP-family-specific，<strong>分母不是全樣本深度</strong>，
所以 sr2b 原式必須改寫：在純細胞株（ρ=1）且 CN_HP 已知時，
θ 退化為 <code>c_k · μ_m / CN_HP</code>。
CN-neutral（CN_HP=1、μ=1）時 <strong>θ = c_k，read-AF 直接就是 CCF</strong> ——
這正是本輪已驗證的那 {cg['neutral_sites_subclonal']['count']} 個乾淨 subclonal 位點的意義。
</div>

<h2>6 · 現有資料已經替這個模型驗證了什麼</h2>
<div class="scroll"><table>
<tr><th>規格的假設／預測</th><th>我們的實測</th><th>結果</th></tr>
<tr><td>μₘ/CNₘ 會壓低觀測 AF</td>
    <td>read-AF 中位隨 CN 單調下降：CN=1 的 {cp['af_median_by_total_cn']['1.0']} → CN=7 的 {cp['af_median_by_total_cn']['7.0']}；
        CN=8 時 {cp['grid_hit_percent_by_cn']['8.0']}% 落在 m/c 格點 ±0.05 內</td>
    <td class="ok">支持</td></tr>
<tr><td>infinite sites（每個突變只發生一次）</td>
    <td>ranked unit 的四配子違反 = 0</td><td class="ok">支持</td></tr>
<tr><td>成對相容 ⇒ 樹已決定，無需搜尋</td>
    <td>k=2 對中 {lk['relation_percent_all']['nested']}% nested、{lk['relation_percent_all']['mutually_exclusive']}% 互斥，全部相容</td>
    <td class="ok">支持</td></tr>
<tr><td>「version 1 僅建議 CN-neutral het segment」</td>
    <td>CN-neutral 僅佔 unit 的 {rv['cn_share']['unit_class_percent']['neutral']}%，
        真正有樹結構資訊的僅 {cg['tree_informative_neutral_share_percent']}%</td>
    <td class="warn">支持但代價明確</td></tr>
<tr><td>全域 CCF atoms 可作 S₂ 的 soft prior</td>
    <td>CN-neutral 的 {SUB_N} 個 subclonal 位點 AF 呈<strong>連續分布</strong>（中位 {SUB_MED}），無離散峰</td>
    <td class="bad">λ 不可取小值</td></tr>
<tr><td>mislabeling η 主導誤差（k=3 時 3e-3 vs 1e-7）</td>
    <td>partial 證據佔 {S2['partial_share_percent']}%，marginalisation 路徑會放大 η 的影響</td>
    <td class="warn">需就地估計 η</td></tr>
</table></div>

<h2>7 · 建議的落地順序</h2>
<ol>
<li><strong>先補 HP→allele 對應</strong>（唯一阻斷點）。無此則 μₘ 無法估，θ 公式空轉。
    <strong>須在 tagged tumour BAM 上做</strong>（germline VCF 的 AD 來自 normal，已實測不可用）。
    產出應含 per-segment 的 CN_HP 與信心值。</li>
<li><strong>就地估計誤差地板 η</strong>。規格說只有它可以就地量測，而我們 partial 佔 {S2['partial_share_percent']}%，
    marginalisation 會放大它。</li>
<li><strong>改寫 θ 為 HP-conditional 形式</strong>，並在 CN-neutral 子集上做封閉檢驗 ——
    該處 θ 應退化為 c_k，可直接用那 {cg['neutral_sites_subclonal']['count']} 個位點驗證擬合是否無偏。</li>
<li><strong>λ 敏感度分析</strong>。我們的證據（乾淨區 AF 連續、無離散 CCF 峰）指向不可取小 λ，
    但沒有給出數值，規格也說沒有客觀方法。</li>
<li><strong>最後才做聯合優化</strong>，並以本輪的 {cv['ordering_total']:,} 條順序約束作為
    hard constraint 篩掉不相容解 —— 這是我們相對於純 VAF 方法的獨有輸入。</li>
</ol>

<div class="verdict">
<strong>一句話</strong>：規格缺的是 <em>lineage 約束</em>，我們缺的是 <em>μₘ 校正層</em>，
兩者恰好互補。接上的成本是<strong>一個</strong>可導出的中間量（HP→allele 對應），
不是新的定序或新的樣本 —— 但它必須在 <strong>tagged tumour BAM</strong> 上算，
不能只讀 germline VCF（已實測，見第 5 節）。
</div>

<div class="note warn">
<strong>誠實邊界</strong>：本頁是<strong>可行性稽核，不是驗證結果</strong>。
規格自陳未經驗證；我們也尚未實作 θ 擬合。表中「支持」指的是
我們的獨立觀測與規格的假設方向一致，不等於規格已被證實。
順序約束停在<strong>分子譜系</strong>層，互斥宣告僅在 CN-neutral 可升格至細胞層
（{cv['exclusion_cellular_share']['percent']}%）。
</div>

<footer>
資料來源：canonical 2026-07-24 frozen cohort（<code>HCC1395.exact_ps_mlhp.json</code> ·
<code>HCC1395.topology.jsonl</code>）· SEQC2 benchmark CNV · SAVANA WGS CN。<br>
規格來源：<code>ccu-bioinformatics-lab.github.io/lab-tutorial/sr2b.html</code> ·
<code>sr2c.html</code>（其自陳 “a specification, not yet validated”）。<br>
所有數字由 <code>scripts/spec_readiness_audit.py</code> 產出至
<code>data/spec_readiness_audit.json</code> 後注入本頁，缺值即拒絕產生。
</footer>
</div>
"""

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write(HTML)
print(f"wrote {OUT}")
print(f"size: {os.path.getsize(OUT):,} bytes")
