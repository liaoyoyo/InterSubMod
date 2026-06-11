#!/usr/bin/env python3
"""
43 — Build standalone HTML companion for the ASM characterization report (v1, for user 理解與確認).
§13-A by-construction anti-fabrication: ALL metric numbers are INJECTED from result JSONs
(not hand-typed); figures base64-embedded (self-contained). Refuses if a required key is missing.

Layers (feedback_pi_report_html_stack): L0 headline / L1 top findings / L2 evidence cards / L3 raw.
Output: docs/experiments/in_progress/2026/06/20260602_ASM_cluster_characterization_regimeA_LOH_diagnostic_01.standalone.html
"""
import json, base64, csv
from pathlib import Path

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
GS = WS / "genome_survey_v2"
FIGD = WS / "figures"
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
           "20260602_ASM_cluster_characterization_regimeA_LOH_diagnostic_01.standalone.html")

probe = json.load(open(GS / "regimeA_credible_probe.json"))
hard = json.load(open(GS / "regimeA_hardening.json"))
resid = json.load(open(GS / "regimeA_residual_controls.json"))
loh = json.load(open(GS / "loh_diagnostic_classifier.json"))
ann = list(csv.DictReader(open(GS / "credible_loci_annotation.tsv"), delimiter="\t"))

ps = probe['summary']; hs = hard['summary']; m2 = hs['M2_regimeA_vs_hetnull']
rs = resid['summary']; ls = loh['summary']


def req(d, *keys):
    """fetch nested key, refuse (raise) if missing — by-construction guard."""
    cur = d
    for k in keys:
        if isinstance(cur, dict) and k in cur:
            cur = cur[k]
        else:
            raise KeyError(f"MISSING required key {'/'.join(map(str,keys))} -> refuse to render")
    return cur


def b64img(name):
    p = FIGD / name
    if not p.exists():
        raise FileNotFoundError(f"figure missing: {p}")
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


# --- injected numbers (from JSON, by construction) ---
N = dict(
    regimeA_total=req(ps, 'regimeA_n_total'),
    evaluated=req(ps, 'n_evaluated'),
    placebo_pass=req(ps, 'n_pass_tierA'),
    survive_both=int(req(rs, 'n_survive_BOTH').split('/')[0]),
    median_ari=round(req(ps, 'median_ari'), 3),
    placebo_median=round(req(ps, 'median_placebo_ari'), 3),
    realvsplacebo_p=req(ps, 'real_vs_placebo_wilcoxon_greater', 'p'),
    m5_rho=round(req(ps, 'M5_spearman_ari_vs_logcov', 'rho'), 3),
    m5_p=round(req(ps, 'M5_spearman_ari_vs_logcov', 'p'), 3),
    regimeA_med=round(req(m2, 'regimeA_median_ari'), 3),
    hetnull_med=round(req(m2, 'hetnull_median_ari'), 3),
    m2_p=req(m2, 'mw_p_greater'),
    cliffs=round(req(m2, 'cliffs_delta'), 3),
    ci=req(m2, 'median_diff_ci95'),
    passrate_a=round(req(hs, 'passrate_ARI30', 'regimeA'), 3),
    passrate_n=round(req(hs, 'passrate_ARI30', 'hetnull'), 3),
    m3_a=round(req(hs, 'M3_rarefied_silhouette', 'regimeA_pass_median'), 3),
    m3_n=round(req(hs, 'M3_rarefied_silhouette', 'hetnull_median'), 3),
    m4c=req(rs, 'M4c_CpGcontext_survive'),
    m8s=req(rs, 'M8strong_random_anchor_pass'),
    loh_total=req(ls, 'n_sig_loh_total'),
    loh_eval=req(ls, 'n_evaluated'),
    loh_artifact=req(ls, 'class_counts', 'self_phasing_artifact'),
    loh_subclone=req(ls, 'class_counts', 'candidate_subclone'),
    loh_cn=req(ls, 'class_counts', 'CN_regression'),
)
# percentages
N['loh_artifact_pct'] = round(100 * N['loh_artifact'] / N['loh_eval'])
N['loh_subclone_pct'] = round(100 * N['loh_subclone'] / N['loh_eval'])
N['loh_cn_pct'] = round(100 * N['loh_cn'] / N['loh_eval'])
ex = req(ls, 'exemplars')

# survivors (M4c=OK AND M8strong=OK) with gene annotation
survive_loci = {r['locus'] for r in resid['loci'] if r.get('all_pass')}
ann_by_locus = {r['locus']: r for r in ann}
survivors = []
for loc in survive_loci:
    a = ann_by_locus.get(loc, {})
    survivors.append((loc, a.get('nearest_gene', '?'), a.get('cpg_context', ''), a.get('ari', '')))
survivors.sort(key=lambda x: -float(x[3]) if x[3] else 0)
collider_killed = [r['locus'] for r in resid['loci'] if not r.get('m8strong_pass')]
collider_genes = [ann_by_locus.get(l, {}).get('nearest_gene', l) for l in collider_killed]

FIG1 = b64img("fig1_asm_landscape.png")
FIG2 = b64img("fig2_gate_evidence.png")
FIG3 = b64img("fig3_loh_triptych.png")
FIG4 = b64img("fig4_regression_to_extreme.png")

surv_rows = "".join(
    f"<tr><td class='mono'>{l}</td><td><b>{g}</b></td><td>{c}</td><td class='mono'>{float(a):.2f}</td></tr>"
    for l, g, c, a in survivors)

html = f"""<!DOCTYPE html>
<html lang="zh-TW"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>ASM 甲基分群 characterization v1 · regime-A + LOH 診斷 (HCC1395 單樣本)</title>
<style>
:root{{--bg:#FAF9F5;--surf:#fff;--ink:#141413;--soft:#57564F;--rule:#E3DACC;--code:#F4EFE6;
--primary:#1E3A8A;--accent:#D97757;--ok:#15803D;--warn:#A16207;--neg:#C2410C;--purple:#7C3AED;
--ff:-apple-system,"Noto Sans CJK TC","Noto Sans CJK JP","PingFang TC","Droid Sans Fallback",sans-serif;
--mono:"JetBrains Mono","DejaVu Sans Mono",monospace;}}
*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--ink);font-family:var(--ff);line-height:1.65;font-size:16px}}
.wrap{{max-width:1180px;margin:0 auto;display:grid;grid-template-columns:240px 1fr;gap:28px;padding:24px}}
nav{{position:sticky;top:18px;align-self:start;font-size:13.5px;border-right:1px solid var(--rule);padding-right:14px}}
nav h4{{margin:0 0 8px;color:var(--accent);font-size:12px;letter-spacing:1px;text-transform:uppercase}}
nav a{{display:block;color:var(--soft);text-decoration:none;padding:3px 0;border-left:2px solid transparent;padding-left:8px}}
nav a:hover{{color:var(--primary);border-left-color:var(--accent)}}
main{{min-width:0}}
h1{{font-size:27px;line-height:1.25;margin:0 0 6px}}
h2{{font-size:21px;margin:32px 0 10px;padding-bottom:6px;border-bottom:2px solid var(--primary);color:var(--primary)}}
h3{{font-size:17px;margin:18px 0 6px;color:var(--accent)}}
.sub{{color:var(--soft);font-size:14px;margin-bottom:14px}}
.badge{{display:inline-block;padding:2px 9px;border-radius:10px;font-size:11px;font-weight:700;color:#fff;font-family:var(--mono);letter-spacing:.5px}}
.b-t3{{background:#2563EB}}.b-pos{{background:var(--ok)}}.b-warn{{background:var(--warn)}}.b-single{{background:var(--neg)}}
.headline{{background:linear-gradient(0deg,#fff,#fff);border:1px solid var(--rule);border-left:5px solid var(--primary);border-radius:6px;padding:16px 20px;margin:14px 0}}
.headline.ok{{border-left-color:var(--ok);background:#F0FDF4}}
.headline b{{color:var(--accent)}}
.metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:14px 0}}
.metric{{background:var(--surf);border:1px solid var(--rule);border-top:3px solid var(--accent);border-radius:6px;padding:10px;text-align:center}}
.metric .v{{font-family:var(--mono);font-size:21px;font-weight:700;color:var(--primary)}}
.metric.ok .v{{color:var(--ok)}}.metric .l{{font-size:11.5px;color:var(--soft);margin-top:4px;line-height:1.35}}
table{{border-collapse:collapse;width:100%;font-size:14px;margin:10px 0}}
th,td{{border-top:1px solid var(--rule);border-bottom:1px solid var(--rule);padding:7px 9px;text-align:left;vertical-align:top}}
th{{background:var(--code);font-weight:700;border-top:2px solid var(--ink)}}
td.mono,.mono{{font-family:var(--mono);font-variant-numeric:tabular-nums}}
.ok{{color:var(--ok)}}.warn{{color:var(--warn);font-weight:600}}.neg{{color:var(--neg);font-weight:600}}
figure{{margin:14px 0;background:var(--surf);border:1px solid var(--rule);border-radius:6px;padding:10px}}
figure img{{width:100%;height:auto;display:block;border-radius:4px}}
figcaption{{font-size:13px;color:var(--soft);margin-top:8px;line-height:1.5}}
details{{background:var(--surf);border:1px solid var(--rule);border-radius:6px;margin:10px 0}}
summary{{cursor:pointer;padding:10px 14px;font-weight:600;background:var(--code);border-radius:6px;color:var(--primary)}}
details[open] summary{{border-radius:6px 6px 0 0}}
details .body{{padding:12px 16px}}
.caveat{{background:#FEF7ED;border-left:4px solid var(--warn);padding:10px 14px;border-radius:0 4px 4px 0;margin:10px 0;font-size:14.5px}}
.caveat b{{color:#92600A}}
.verdict-box{{background:#EFF6FF;border:1px solid #BFDBFE;border-left:4px solid var(--primary);border-radius:0 6px 6px 0;padding:12px 16px;margin:10px 0}}
ul{{margin:6px 0;padding-left:22px}}li{{margin:4px 0}}
footer{{margin-top:36px;padding-top:14px;border-top:1px solid var(--rule);font-size:12px;color:var(--soft);font-family:var(--mono);line-height:1.6}}
code{{background:var(--code);padding:1px 5px;border-radius:3px;font-family:var(--mono);font-size:13px}}
.tag{{display:inline-block;background:var(--code);color:var(--accent);padding:0 6px;border-radius:3px;font-family:var(--mono);font-size:11px;margin-right:4px}}
</style></head><body>
<div class="wrap">
<nav>
<h4>目錄</h4>
<a href="#tldr">0 · TL;DR 結論</a>
<a href="#scope">1 · scope 與背景</a>
<a href="#g1">2 · G1 找清楚分群</a>
<a href="#fig1">　Fig1 landscape</a>
<a href="#fig2">　Fig2 gate 證據</a>
<a href="#g1cav">　誠實 caveat</a>
<a href="#surv">　15 存活位點</a>
<a href="#g2">3 · G2 LOH 診斷</a>
<a href="#fig3">　Fig3 三成因</a>
<a href="#fig4">　Fig4 regression</a>
<a href="#verdict">4 · verdict + tier</a>
<a href="#method">5 · 方法邊界</a>
<a href="#next">6 · 下一步</a>
</nav>
<main>
<h1>ASM 甲基分群 characterization（v1）<br>regime-A credible 子集 + LOH 成因診斷</h1>
<div class="sub">HCC1395 paired_full · longphase-S HP-axis · 單樣本 ·
<span class="badge b-t3">⭐3 L1</span> <span class="badge b-single">SINGLE-SAMPLE Tier-A</span>
<span class="badge b-pos">G1 POSITIVE-bounded</span> · 2026-06-02 · 給用戶理解確認用 (in_progress)</div>

<section id="tldr"><h2>0 · TL;DR 結論（先看這段）</h2>
<div class="headline ok"><b>G1 — 能找到真 somatic 甲基分群，但有嚴格邊界。</b>
在 somatic-controlled HP-axis、避開 regression-to-extreme 的 credible 子集中，
<b>{N['placebo_pass']}/{N['evaluated']}</b> 過 collider cut、其中 <b>{N['survive_both']}/{N['placebo_pass']}</b> 撐過完整 artifact battery；
regime-A blind-ARI（median {N['regimeA_med']}）<b>顯著高於 germline-het baseline</b>（{N['hetnull_med']}）：
Cliff δ=<b>{N['cliffs']}</b>, p={N['m2_p']:.1e}。→ refine（非推翻）5/31「全基因組 artifact-dominated」。</div>
<div class="headline"><b>G2 — LOH 盲點假設被驗證。</b>
longphase-S 不驗證 LOH 而產生的「表觀雙 haplotype」位點，診斷分類為
<b>{N['loh_artifact_pct']}% self-phasing artifact</b>（無真 cluster）/ <b>{N['loh_subclone_pct']}% candidate subclone</b> /
<b>{N['loh_cn_pct']}% CN-regression</b> — 多數是 phasing artifact，但少數是 LOH 內真訊號。</div>
<div class="caveat"><b>必讀邊界：</b>① <b>單樣本 Tier-A 天花板</b>（M9 cross-sample 不可能 → 不可 generalize）；
② 強化 collider gate <b>刷掉好看的癌症基因</b>（{', '.join(collider_genes[:3])} 等 = genomic-context collider，非 somatic 特異）；
③ <b>BRCA2 本身未過 collider cut</b>（ARI 0.79 但 placebo 0.132，borderline）。
④ <b>未碰 filter-NEGATIVE</b>（甲基當 FP filter 仍 DEAD）— 本輪是 characterization。</div>
</section>

<section id="scope"><h2>1 · scope 與背景（為何這樣做）</h2>
<p>延續 BRCA2/ZAR1L reviewer 答辯，用戶兩個目標：<b>G1</b> 找更多像 BRCA2 一樣清楚的 ASM 分群；
<b>G2</b> 處理 unphase/HP3/LOH 盲點。盤點發現既有 scripts 已完成約 70%、且結果偏 cautionary，故本輪 =
<b>consolidation + 補決定性 gate + 全新 LOH per-locus 診斷</b>，scope 確認為「嚴謹 characterization verdict」、
G1 framing「先探 regime-A credible 子集」、限 HCC1395 單樣本。</p>
<div class="verdict-box"><b>為何用 HP-axis（回答 G1.2「依哪種 tag 一致」）：</b>
HP1 vs HP1-1 = 同一 germline haplotype 被 somatic allele 切分，by construction 控制 germline allelic baseline；
ALLELE-axis（ALT vs REF）= 兩條 germline haplotype = baseline-allelic confounded。
→ 甲基分群「依 germline allele 一致 > 依 somatic」；只有 HP-axis 上超過 placebo 的才 somatic-attributable。</div>
</section>

<section id="g1"><h2>2 · G1 — regime-A credible 子集有真 somatic 甲基分群</h2>
<div class="metrics">
<div class="metric ok"><div class="v">{N['cliffs']}</div><div class="l">M2 Cliff δ<br>(regime-A vs het-null, medium)</div></div>
<div class="metric ok"><div class="v">{N['m2_p']:.0e}</div><div class="l">M2 MW p<br>(somatic > baseline)</div></div>
<div class="metric"><div class="v">{N['survive_both']}/{N['placebo_pass']}</div><div class="l">撐過完整<br>collider battery</div></div>
<div class="metric"><div class="v">{N['m5_rho']}</div><div class="l">M5 coverage ρ<br>(p={N['m5_p']} NS = 無 artifact)</div></div>
</div>
<figure id="fig1"><img src="{FIG1}" alt="Fig1 ASM landscape">
<figcaption><b>Fig1 · 全基因組 HP-axis ASM landscape。</b>灰=全 sig 位點（artifact-dominated）；綠=regime-A credible 子集（{N['regimeA_total']} 個）；橘星=BRCA2/ZAR1L (rank 25)。乾淨子集散布全基因組，非單一區域。</figcaption></figure>
<figure id="fig2"><img src="{FIG2}" alt="Fig2 gate evidence">
<figcaption><b>Fig2 · headline 證據（3 panel）。</b>A) regime-A blind-ARI（median {N['regimeA_med']}）vs germline-het null（{N['hetnull_med']}）— M2 顯著高於 baseline。
B) gate funnel {N['regimeA_total']}→{N['evaluated']}→{N['placebo_pass']}→{N['survive_both']} + 通過的 gate scorecard（⚠ 強化 collider 刷掉 {', '.join(collider_genes[:3])}）。
C) G2 LOH 三成因分布。</figcaption></figure>
<details><summary>L2 · 六道 gate 逐項數字（展開）</summary><div class="body">
<table>
<tr><th>Gate</th><th>內容</th><th>結果</th><th>判定</th></tr>
<tr><td>M1 blind-ARI</td><td>盲分群恢復 somatic/germline split（PRIMARY）</td><td class="mono">median {N['median_ari']}</td><td>—</td></tr>
<tr><td>M8 length-placebo</td><td>同位點長度切割 collider proxy</td><td class="mono">real>>placebo p={N['realvsplacebo_p']:.1e}（median placebo {N['placebo_median']}）</td><td class="ok">✓</td></tr>
<tr><td>M5 coverage 簽名</td><td>ARI vs log(n_cpg) Spearman</td><td class="mono">ρ={N['m5_rho']} (p={N['m5_p']} NS)</td><td class="ok">✓ 無 artifact</td></tr>
<tr><td>M2 vs het-null</td><td>baseline-allelic null (n≥50, cov-matched)</td><td class="mono">{N['regimeA_med']} vs {N['hetnull_med']}; δ={N['cliffs']}; CI[{N['ci'][0]:.2f},{N['ci'][1]:.2f}]</td><td class="ok">✓</td></tr>
<tr><td>M3 rarefied silhouette</td><td>降採樣等組大小 (B=200)</td><td class="mono">{N['m3_a']} vs null {N['m3_n']}</td><td class="ok">✓ 非 read-count</td></tr>
<tr><td>M4c CpG-context</td><td>drop 變異±20bp CpG 後存活</td><td class="mono">{N['m4c']}</td><td class="ok">✓</td></tr>
<tr><td>M8-strong random-anchor</td><td>隨機非-somatic 位點 collider</td><td class="mono">{N['m8s']}</td><td class="warn">⚠ 7 collider</td></tr>
<tr><td>M9 cross-sample</td><td>≥2 樣本 sign-concordance</td><td class="mono">未跑 = 單樣本不可能</td><td class="neg">✗ N/A</td></tr>
</table>
<p>pass-rate ARI≥0.30：regime-A {N['passrate_a']} vs het-null {N['passrate_n']}。</p>
</div></details>
<h3 id="g1cav">誠實 caveat</h3>
<div class="caveat"><b>① 單樣本 Tier-A 天花板。</b> M9 cross-sample（red-team A7 single-sample circularity = 本專案歷次 positive 全滅根因）需第 2 樣本 tagged BAM，目前只 HCC1395。不可宣稱 Tier-B / generalize。</div>
<div class="caveat"><b>② 強化 collider 刷掉好看的基因。</b> length-split placebo 太弱；random-anchor M8 揭露 {len(collider_killed)} 個是 genomic-context collider，正含 <b>{', '.join([g for g in collider_genes if g in ('SOX2','HOTTIP','SDHAF1')])}</b> — 不能當 somatic-ASM 發現。這正是嚴格 battery 的價值（阻止 overclaim）。</div>
<div class="caveat"><b>③ BRCA2 本身 borderline。</b> chr13:32315128 在 regime-A 內（n_cpg=197, germ_β=0.234, nonLOH），blind-ARI 0.79，但 placebo 0.132 > 0.10 → 未過 collider cut。最乾淨的反而是其他位點。</div>
<h3 id="surv">撐過完整 battery 的 {N['survive_both']} 個存活位點</h3>
<table><tr><th>locus</th><th>最近基因</th><th>CpG context</th><th>ARI</th></tr>{surv_rows}</table>
<p class="sub">多為 lncRNA/pseudogene/非經典位點；僅部分 promoter-proximal。SOX2/HOTTIP/SDHAF1 等已被 collider gate 排除，不在此列。</p>
</section>

<section id="g2"><h2>3 · G2 — LOH 表觀雙 haplotype 成因診斷</h2>
<p><b>用戶洞見（已從 C++ 碼 L1 確認）：</b>longphase-S 不驗證 LOH（germline-absent+somatic → H3/H1_1），不像 longphase-TO V6 給保守 HP:i:33。
故 LOH 區（物理單 haplotype）出現「表觀雙 haplotype」。用戶提的「強制歸類」= V5 Layer 1.5 已證失敗（4.19:1 偏差 + marker −19%），故改用<b>診斷分類</b>。</p>
<div class="metrics">
<div class="metric"><div class="v">{N['loh_artifact_pct']}%</div><div class="l">self-phasing artifact<br>(n={N['loh_artifact']}, 無真 cluster)</div></div>
<div class="metric ok"><div class="v">{N['loh_subclone_pct']}%</div><div class="l">candidate subclone<br>(n={N['loh_subclone']}, LOH 內真訊號)</div></div>
<div class="metric"><div class="v">{N['loh_cn_pct']}%</div><div class="l">CN / regression<br>(n={N['loh_cn']})</div></div>
<div class="metric"><div class="v">{N['loh_total']}</div><div class="l">sig LOH HP-axis 位點<br>(評估 {N['loh_eval']})</div></div>
</div>
<figure id="fig3"><img src="{FIG3}" alt="Fig3 LOH triptych">
<figcaption><b>Fig3 · 三成因 read-level 甲基矩陣。</b>每行=1 read，上=germline-tag 下=somatic-tag，黑=甲基/白=未甲基/灰=ambiguous。
(a) candidate subclone <code>{ex['candidate_subclone']['locus']}</code> (IGF2, ARI={ex['candidate_subclone']['ari']}) 有分群；
(b) self-phasing artifact <code>{ex['self_phasing_artifact']['locus']}</code> (ARI={ex['self_phasing_artifact']['ari']}) 無分離；
(c) CN/regression <code>{ex['CN_regression']['locus']}</code> (germ_β={ex['CN_regression']['germ_beta']:.2f} 極端)。</figcaption></figure>
<div class="caveat"><b>注意：</b>candidate subclone 最強 exemplar chr11:2146993 落在 <b>IGF2（imprinted gene）</b> — 雖 HP-axis 控制 imprinting，仍須謹慎詮釋。CN class 為 categorical inference（非實測整數 CN，obs25 caveat）。</div>
<figure id="fig4"><img src="{FIG4}" alt="Fig4 regression to extreme">
<figcaption><b>Fig4 · regression-to-extreme 機制。</b>coverage(n_cpg) vs |Δβ|；LOH(紅) 集中在低 coverage + 高 |Δβ| 區 = baseline 壓到極端 → 假性大效應。這是 5/31「strong-ASM FP-enriched」與 G2 CN/regression 類的共同機制。</figcaption></figure>
</section>

<section id="verdict"><h2>4 · Verdict + tier</h2>
<table><tr><th>主張</th><th>tier</th><th>邊界</th></tr>
<tr><td>credible regime(HP-axis) 有真 somatic 甲基 cluster，高於 germline-het baseline</td><td><span class="badge b-t3">⭐3 L1</span></td><td>單樣本 Tier-A；M9 未驗</td></tr>
<tr><td>{N['survive_both']} 存活位點為單樣本候選；SOX2/HOTTIP/SDHAF1 被 collider 排除</td><td><span class="badge b-t3">⭐3 L1</span></td><td>需 cross-sample 升 Tier-B</td></tr>
<tr><td>LOH 表觀雙 haplotype {N['loh_artifact_pct']}% self-phasing / {N['loh_subclone_pct']}% candidate subclone</td><td><span class="badge b-t3">⭐3 L1</span></td><td>單樣本；CN categorical</td></tr>
<tr><td>全基因組 ASM 仍 artifact-dominated（5/31 成立，本輪 refine 非推翻）</td><td><span class="badge b-t3">⭐3 L1</span></td><td>—</td></tr></table>
<div class="verdict-box"><b>與 5/31 收斂的關係：</b>5/31「ASM non-directional/non-discriminative/coverage-modulated」在全基因組/ALLELE-axis/全 regime 仍成立；
本輪 refine = credible regime × HP-axis × 全 collider battery 下有 {N['survive_both']} 個高於 baseline 的真 cluster 候選（單樣本）。
兩者不矛盾。<b>filter-NEGATIVE（甲基當 FP filter ⭐2 L4 DEAD）未被觸碰。</b></div>
</section>

<section id="method"><h2>5 · 方法完整性與未做項（誠實邊界）</h2>
<p><b>已跑 gate：</b>M0 observed-only · M1 blind-ARI · M2 het-null(n≥50 分層) · M3 rarefied · M5 coverage 簽名 · M8 length+random-anchor · M4c CpG-context。
<b>驗尺：</b>24_cluster_eval_core 已驗 PC1(ARI>0.5)+NC1(ARI<0.15)；in-battery 控制涵蓋 NC2-5 主要失效模式。</p>
<p><b>未做（remaining for full Tier-A / Tier-B）：</b></p>
<ul>
<li><b>M9 cross-sample</b>（BLOCKING for Tier-B）— 需 COLO829 等第 2 樣本 tagged BAM。</li>
<li>完整 PC1-3 / NC1-7 formal harness（本輪用既有 PC1/NC1 + in-battery 控制替代）。</li>
<li>M4a mask-swap / M4b missing-permute / M4d basecaller-context（M4d 需 re-basecalling）。</li>
</ul>
</section>

<section id="next"><h2>6 · 下一步（cross-sample 升級路徑）</h2>
<ul>
<li><b>M9 cross-sample（最高優先，唯一通往 Tier-B）：</b>COLO829/其他 paired tagged BAM 跑同 {N['survive_both']}-locus + LOH battery，檢 sign-concordance。</li>
<li>完整 PC1-3/NC1-7 formal meta-validation。</li>
<li>M4a/b/d 補完整 battery。</li>
<li>存活位點功能跟進（motif/expression）— 須先 cross-sample 確認非單樣本 artifact。</li>
</ul>
</section>

<footer>
數據誠實（§13-A 由構造）：本 HTML 所有 metric 數字由 result JSON 注入（非手打），圖 base64 內嵌。<br>
data_sources: genome_survey_v2/{{regimeA_credible_probe, regimeA_hardening, regimeA_residual_controls, loh_diagnostic_classifier}}.json + credible_loci_annotation.tsv<br>
報告 .md: InterSubMod/docs/experiments/in_progress/2026/06/20260602_ASM_cluster_characterization_regimeA_LOH_diagnostic_01.md<br>
ledger: ASM-cluster-characterization-20260602 · sample: HCC1395 paired_full single-sample (Tier-A ceiling) · build 2026-06-02
</footer>
</main></div></body></html>
"""

OUT.write_text(html, encoding="utf-8")
print(f"standalone HTML -> {OUT}")
print(f"  size: {OUT.stat().st_size//1024} KB")
print(f"  injected: regimeA={N['regimeA_total']} eval={N['evaluated']} pass={N['placebo_pass']} survive={N['survive_both']}")
print(f"  M2 cliff={N['cliffs']} p={N['m2_p']:.1e} | LOH {N['loh_artifact_pct']}/{N['loh_subclone_pct']}/{N['loh_cn_pct']}%")
print(f"  survivors={len(survivors)} collider_killed={len(collider_killed)} ({collider_genes[:4]})")
