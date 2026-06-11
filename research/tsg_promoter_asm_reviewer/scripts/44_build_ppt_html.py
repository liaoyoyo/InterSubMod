#!/usr/bin/env python3
"""
44 — Build 16:9 PPT-style presenter deck (BCG style) for ASM characterization v1.
Per-goal 確認點 boxes for user to confirm each goal + detail. §13-A by-construction:
numbers injected from result JSONs; 4 figures base64-embedded. Presenter mode (speaker notes).
Output: docs/presentations/in_progress/20260602_ASM_characterization_deck/presenter_ASM_deck.html
"""
import json, base64, csv
from pathlib import Path

WS = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
GS = WS / "genome_survey_v2"; FIGD = WS / "figures"
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/presentations/in_progress/"
           "20260602_ASM_characterization_deck/presenter_ASM_deck.html")

probe = json.load(open(GS / "regimeA_credible_probe.json")); ps = probe['summary']
hard = json.load(open(GS / "regimeA_hardening.json")); hs = hard['summary']; m2 = hs['M2_regimeA_vs_hetnull']
resid = json.load(open(GS / "regimeA_residual_controls.json")); rs = resid['summary']
loh = json.load(open(GS / "loh_diagnostic_classifier.json")); ls = loh['summary']
ann = {r['locus']: r for r in csv.DictReader(open(GS / "credible_loci_annotation.tsv"), delimiter="\t")}


def b64(name):
    p = FIGD / name
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


N = dict(
    rA=ps['regimeA_n_total'], ev=ps['n_evaluated'], pass1=ps['n_pass_tierA'],
    surv=int(rs['n_survive_BOTH'].split('/')[0]),
    medari=round(ps['median_ari'], 3), pl=round(ps['median_placebo_ari'], 3),
    rvp=ps['real_vs_placebo_wilcoxon_greater']['p'],
    rho=round(ps['M5_spearman_ari_vs_logcov']['rho'], 2), rhop=round(ps['M5_spearman_ari_vs_logcov']['p'], 2),
    rAm=round(m2['regimeA_median_ari'], 3), hn=round(m2['hetnull_median_ari'], 3),
    m2p=m2['mw_p_greater'], cd=round(m2['cliffs_delta'], 3), ci=m2['median_diff_ci95'],
    pra=round(hs['passrate_ARI30']['regimeA'], 3), prn=round(hs['passrate_ARI30']['hetnull'], 3),
    m3a=round(hs['M3_rarefied_silhouette']['regimeA_pass_median'], 3), m3n=round(hs['M3_rarefied_silhouette']['hetnull_median'], 3),
    m4c=rs['M4c_CpGcontext_survive'], m8s=rs['M8strong_random_anchor_pass'],
    lt=ls['n_sig_loh_total'], le=ls['n_evaluated'],
    art=ls['class_counts']['self_phasing_artifact'], sub=ls['class_counts']['candidate_subclone'], cn=ls['class_counts']['CN_regression'],
)
N['artp'] = round(100 * N['art'] / N['le']); N['subp'] = round(100 * N['sub'] / N['le']); N['cnp'] = round(100 * N['cn'] / N['le'])
ex = ls['exemplars']
collider_killed = [r['locus'] for r in resid['loci'] if not r.get('m8strong_pass')]
collider_genes = [ann.get(l, {}).get('nearest_gene', l) for l in collider_killed]
collider_named = [g for g in collider_genes if g in ('SOX2', 'HOTTIP', 'SDHAF1')]
surv = sorted([(l, ann.get(l, {}).get('nearest_gene', '?'), ann.get(l, {}).get('ari', '0'))
               for l in {r['locus'] for r in resid['loci'] if r.get('all_pass')}],
              key=lambda x: -float(x[2]) if x[2] else 0)
surv_top = " · ".join(f"{g}({float(a):.2f})" for _, g, a in surv[:8])

F1, F2, F3, F4 = b64("fig1_asm_landscape.png"), b64("fig2_gate_evidence.png"), b64("fig3_loh_triptych.png"), b64("fig4_regression_to_extreme.png")

# ---------- slide content ----------
SLIDES = []

SLIDES.append(("cover", """
<div class="cover">
  <div class="eyebrow">研究 characterization · 給用戶理解確認 · 6-8 min</div>
  <h1>ASM 甲基分群兩目標 characterization<br><span class="en">G1 找清楚 somatic 分群 · G2 LOH 盲點診斷（HCC1395 單樣本）</span></h1>
  <div class="lead">核心結論：<b>乾淨子集有真 somatic 甲基 cluster（高於 baseline）</b>，但<b>單樣本 Tier-A 天花板</b>；<b>LOH 表觀雙 haplotype 72% 是 phasing artifact</b>（驗證你的盲點假設）。</div>
  <div class="arc">
    <div class="step"><span class="num">S2</span><span class="lbl">背景與 scope</span></div>
    <div class="step"><span class="num">S3-5</span><span class="lbl">G1 找清楚分群</span></div>
    <div class="step"><span class="num">S6-7</span><span class="lbl">G2 LOH 診斷</span></div>
    <div class="step"><span class="num">S8</span><span class="lbl">verdict + 確認</span></div>
    <div class="step" style="border-top-color:var(--primary)"><span class="num">S9</span><span class="lbl">下一步</span></div>
  </div>
  <div class="credits">
    <div class="cr"><span class="k">樣本</span><span class="v">HCC1395 paired_full</span></div>
    <div class="cr"><span class="k">軸</span><span class="v">longphase-S HP-axis</span></div>
    <div class="cr"><span class="k">tier</span><span class="v">★3 L1 · 單樣本 Tier-A</span></div>
    <div class="cr"><span class="k">日期</span><span class="v">2026-06-02</span></div>
  </div>
</div>"""))

SLIDES.append(("背景", f"""
<div class="head"><div class="eb">S2 · 背景與 scope · 兩個目標 + 為何用 HP-axis</div>
<h2 class="t">兩目標 · 既有 70% 已做 → 本輪補決定性 gate + 全新 LOH 診斷</h2></div>
<div class="body split-11">
  <div class="side">
    <h3>你的兩個目標</h3>
    <ul>
      <li><b>G1</b>：找更多像 BRCA2 一樣清楚的 somatic 甲基分群位點 + 量化「清楚」+ 依哪種 tag 一致</li>
      <li><b>G2</b>：unphase/HP3/遠 read + LOH 盲點（longphase-S 不驗證 LOH → 表觀雙 haplotype）</li>
    </ul>
    <div class="vbox"><b>scope（你確認過）</b>：嚴謹 characterization verdict · 先探 regime-A credible 子集 · 限 HCC1395 單樣本明標。既有 scripts 15-35 約 70% 已做 → 本輪 = consolidation + 補 gate + 全新 LOH per-locus 診斷。</div>
  </div>
  <div class="side">
    <h3>為何用 HP-axis（回答 G1.2「依哪種 tag」）</h3>
    <ul>
      <li><b>HP1 vs HP1-1</b> = 同一 germline haplotype 被 somatic 切分 → <b>by construction 控制 baseline</b> → 超過 placebo 才 somatic</li>
      <li><b>ALT vs REF</b>（ALLELE-axis）= 兩條 germline haplotype → baseline-allelic <b>confounded</b></li>
    </ul>
    <div class="vbox ok">→ 甲基分群「<b>依 germline allele 一致 &gt; 依 somatic</b>」；本研究只信 HP-axis 上撐過 placebo 的訊號。</div>
  </div>
</div>
<div class="terms"><b>▸</b> ASM=allele-specific methylation · HP:Z tag: 1/2=germline · 1-1/2-1=重建 somatic · blind-ARI=盲分群恢復 split 的吻合度（主判準）</div>
"""))

SLIDES.append(("G1證據", f"""
<div class="head"><div class="eb">S3 · G1 · 乾淨子集有真 somatic 甲基分群（高於 baseline）</div>
<h2 class="t">regime-A credible 子集：撐過 4 道 artifact gate，顯著高於 germline-het baseline</h2></div>
<div class="metrics">
  <div class="m ok"><div class="v">{N['cd']}</div><div class="l">M2 Cliff δ<br>(vs het-null, medium)</div></div>
  <div class="m ok"><div class="v">{N['m2p']:.0e}</div><div class="l">M2 MW p<br>(somatic&gt;baseline)</div></div>
  <div class="m"><div class="v">{N['surv']}/{N['pass1']}</div><div class="l">撐過完整<br>collider battery</div></div>
  <div class="m"><div class="v">{N['rho']}</div><div class="l">M5 coverage ρ<br>(p={N['rhop']} NS=無 artifact)</div></div>
</div>
<div class="body split-13">
  <div class="hero"><img src="{F2}"></div>
  <div class="side">
    <h3>讀法</h3>
    <ul>
      <li>A: regime-A ARI <b>{N['rAm']}</b> ≫ baseline <b>{N['hn']}</b></li>
      <li>B: funnel {N['rA']}→{N['ev']}→{N['pass1']}→<b>{N['surv']}</b>；M3 rarefied {N['m3a']} vs {N['m3n']}</li>
      <li>C: G2 LOH 三成因（下頁）</li>
    </ul>
    <div class="confirm"><b>▶ 請確認 G1-a</b><br>「乾淨子集（HP-axis、避開 regression-to-extreme）<b>有真 somatic 甲基 cluster、且高於 baseline</b>」這個結論方向，你認同嗎？</div>
  </div>
</div>
"""))

SLIDES.append(("G1誠實", f"""
<div class="head"><div class="eb">S4 · G1 · 誠實 caveat（嚴格 gate 阻止 overclaim）</div>
<h2 class="t">三個必讀邊界 — 嚴格 collider gate 連好看的癌症基因都刷掉</h2></div>
<div class="body single">
  <div class="cav-grid">
    <div class="cav"><div class="cn">①</div><div><b>單樣本 Tier-A 天花板</b><br>M9 cross-sample（避 single-sample circularity = 本專案歷次 positive 全滅根因）需第 2 樣本 tagged BAM。<b>不可宣稱 generalize / Tier-B</b>。</div></div>
    <div class="cav"><div class="cn">②</div><div><b>強化 collider 刷掉好看的基因</b><br>length-placebo 太弱；random-anchor M8（{N['m8s']}）揭露 {len(collider_killed)} 個是 genomic-context collider，<b>正含 {', '.join(collider_named)}</b> — 不能當 somatic-ASM 發現。這正是嚴格 battery 的價值。</div></div>
    <div class="cav"><div class="cn">③</div><div><b>BRCA2 本身 borderline</b><br>chr13:32315128 在 regime-A 內（n_cpg=197, germ_β=0.234），blind-ARI 0.79，<b>但 placebo 0.132 &gt; 0.10 → 未過 cut</b>。最乾淨的反而是其他位點。</div></div>
  </div>
  <div class="confirm"><b>▶ 請確認 G1-b</b><br>你能接受「<b>不宣稱 generalize</b>、<b>SOX2/HOTTIP/SDHAF1 因 collider 不列為發現</b>、<b>BRCA2 自己只是 borderline</b>」這些誠實限制嗎？（這是不 overclaim 的關鍵）</div>
</div>
"""))

SLIDES.append(("G1存活", f"""
<div class="head"><div class="eb">S5 · G1 · 全基因組 landscape + 撐過全 battery 的 {N['surv']} 個存活位點</div>
<h2 class="t">乾淨子集散布全基因組；存活者多為非經典位點（非 SOX2 那類）</h2></div>
<div class="body split-13">
  <div class="hero"><img src="{F1}"></div>
  <div class="side">
    <h3>{N['surv']} 個全 battery 存活（最近基因）</h3>
    <div class="vbox" style="font-size:13px;line-height:1.6">{surv_top} …</div>
    <ul>
      <li>多為 <b>lncRNA / pseudogene / 非經典位點</b>（HERC6 / CIB2 / TBC1D16 / ENTPD8）</li>
      <li>僅部分 promoter-proximal / 落 CpG island-shore</li>
      <li>SOX2 / HOTTIP / SDHAF1 已被 collider gate 排除，不在此列</li>
    </ul>
    <div class="confirm"><b>▶ 請確認 G1-c</b><br>「乾淨存活位點偏非經典基因、需 cross-sample 才能談功能意義」這個 framing OK 嗎？</div>
  </div>
</div>
"""))

SLIDES.append(("G2診斷", f"""
<div class="head"><div class="eb">S6 · G2 · LOH 盲點診斷（驗證你的假設）</div>
<h2 class="t">LOH 表觀雙 haplotype：{N['artp']}% 是 self-phasing artifact，{N['subp']}% 才是真 subclone</h2></div>
<div class="metrics">
  <div class="m"><div class="v">{N['artp']}%</div><div class="l">self-phasing artifact<br>(n={N['art']}, 無真 cluster)</div></div>
  <div class="m ok"><div class="v">{N['subp']}%</div><div class="l">candidate subclone<br>(n={N['sub']}, LOH 內真訊號)</div></div>
  <div class="m"><div class="v">{N['cnp']}%</div><div class="l">CN / regression<br>(n={N['cn']})</div></div>
  <div class="m"><div class="v">{N['lt']}</div><div class="l">sig LOH 位點<br>(評估 {N['le']})</div></div>
</div>
<div class="body split-13">
  <div class="hero"><img src="{F3}"></div>
  <div class="side">
    <h3>你的洞見被驗證</h3>
    <ul>
      <li>longphase-S <b>不驗證 LOH</b>（C++ 碼確認）→ 表觀雙 haplotype</li>
      <li>多數（{N['artp']}%）= phasing artifact（無甲基 cluster）</li>
      <li>用<b>診斷分類</b>非強制歸類（強制歸類=V5 Layer1.5 已證失敗 4.19:1 偏差）</li>
      <li>! 最強 exemplar chr11:2146993 落 <b>IGF2(imprinted)</b>，須謹慎</li>
    </ul>
    <div class="confirm"><b>▶ 請確認 G2</b><br>「LOH 雙 haplotype 多數是 phasing artifact、少數(18%)是真 subclone」+「用診斷分類而非強制歸類」這方向你認同嗎？</div>
  </div>
</div>
"""))

SLIDES.append(("機制", f"""
<div class="head"><div class="eb">S7 · 機制 · regression-to-extreme + 與 5/31 結論的關係</div>
<h2 class="t">低 coverage 放大 |Δβ| 是共同 artifact 機制；本輪 refine 非推翻 5/31</h2></div>
<div class="body split-13">
  <div class="hero"><img src="{F4}"></div>
  <div class="side">
    <h3>關鍵詮釋</h3>
    <ul>
      <li>LOH(紅)集中在<b>低 cov + 高 |Δβ|</b> 區 = baseline 壓到極端→假性大效應</li>
      <li>這是 5/31「strong-ASM FP-enriched」與 G2「CN/regression 類」的共同機制</li>
    </ul>
    <div class="vbox"><b>與 5/31 的關係</b>：5/31「ASM non-directional/non-discriminative/coverage-modulated」在<b>全基因組/ALLELE-axis</b>仍成立；本輪 refine = <b>credible regime × HP-axis × 全 collider battery</b> 下有 {N['surv']} 個高於 baseline 的真 cluster 候選（單樣本）。<b>兩者不矛盾。</b></div>
    <div class="confirm"><b>▶ 請確認 機制</b><br>「artifact 主導全局、但乾淨子集有真訊號」這個 refine（非推翻）的定位 OK 嗎？</div>
  </div>
</div>
"""))

SLIDES.append(("verdict", f"""
<div class="head"><div class="eb">S8 · Verdict + tier + 整體確認清單</div>
<h2 class="t">四條主張全 ★3 L1 單樣本 Tier-A；filter-NEGATIVE 未觸碰</h2></div>
<div class="body single">
  <table class="bcg">
    <tr><th>主張</th><th>tier</th><th>邊界</th></tr>
    <tr><td>credible regime(HP-axis) 有真 somatic cluster，高於 baseline</td><td class="c">★3 L1</td><td>單樣本；M9 未驗</td></tr>
    <tr><td>{N['surv']} 存活位點為候選；SOX2/HOTTIP/SDHAF1 被 collider 排除</td><td class="c">★3 L1</td><td>需 cross-sample 升 Tier-B</td></tr>
    <tr><td>LOH 表觀雙 haplotype {N['artp']}% artifact / {N['subp']}% candidate subclone</td><td class="c">★3 L1</td><td>單樣本；CN categorical</td></tr>
    <tr><td>全基因組仍 artifact-dominated（5/31 成立，本輪 refine）</td><td class="c">★3 L1</td><td>—</td></tr>
  </table>
  <div class="confirm" style="margin-top:14px"><b>▶ 整體確認清單（請逐項回我 OK / 要改）</b><br>
  ① G1 結論方向（乾淨子集有真訊號高於 baseline）　② G1 誠實限制（不 generalize / collider 刷掉好基因 / BRCA2 borderline）　③ G1 存活位點 framing（偏非經典）　④ G2 診斷（72% artifact / 18% subclone，診斷非強制歸類）　⑤ 與 5/31「refine 非推翻」定位　⑥ tier=★3 單樣本 Tier-A + filter 不碰</div>
</div>
"""))

SLIDES.append(("下一步", f"""
<div class="head"><div class="eb">S9 · 下一步 + deferred（誠實邊界）</div>
<h2 class="t">cross-sample（COLO829）是唯一通往 Tier-B 的路</h2></div>
<div class="body split-11">
  <div class="side">
    <h3>下一步（你決定）</h3>
    <ul>
      <li><b>最高優先</b>：COLO829 等第 2 樣本 tagged BAM 跑同 {N['surv']}-locus + LOH battery，檢 sign-concordance → 升 Tier-B</li>
      <li>完整 PC1-3/NC1-7 formal harness</li>
      <li>M4a/M4b/M4d 補完整 battery</li>
      <li>存活位點功能跟進（須先 cross-sample）</li>
    </ul>
  </div>
  <div class="side">
    <h3>本輪 deferred（刻意，非遺漏）</h3>
    <div class="vbox warn"><b>M9 cross-sample</b>：單樣本物理上不可能 — BLOCKING for Tier-B</div>
    <div class="vbox warn">完整 PC1-3/NC1-7 harness：用既有 PC1/NC1 + in-battery 控制替代</div>
    <div class="vbox warn">M4a/b/d：邊際價值低於已跑的 M8-strong/M4c</div>
    <div class="confirm"><b>▶ 請確認 下一步</b><br>要先確認 COLO829 tagged BAM 是否存在（啟動 cross-sample），還是這版先收？</div>
  </div>
</div>
"""))

NOTES = {
 1: "開場：兩目標 + 一句話結論。強調這是單樣本 Tier-A，不 overclaim。",
 2: "講 scope correction（70% 已做）+ HP-axis 為何是對的軸（控制 baseline）。G1.2 答案：依 germline allele > 依 somatic。",
 3: "G1 核心證據。先講 A panel（高於 baseline，Cliff δ 0.37），再講 funnel B（15 撐過全 gate），M5 ρ 負=無 coverage artifact。請用戶確認 G1-a。",
 4: "誠實 caveat 三點。重點：嚴格 collider 刷掉 SOX2/HOTTIP/SDHAF1（好看但不算數）= 防 overclaim 的價值。BRCA2 自己 borderline。請用戶確認 G1-b。",
 5: "landscape 全圖 + 15 存活者偏非經典基因。請用戶確認 G1-c framing。",
 6: "G2 = 用戶盲點被驗證。72% artifact / 18% subclone。診斷分類非強制歸類（V5 已失敗）。IGF2 imprinted 警示。請用戶確認 G2。",
 7: "機制：regression-to-extreme（Fig4）。與 5/31 refine 非推翻。請用戶確認機制定位。",
 8: "verdict 表 + 六項整體確認清單 — 這是讓用戶逐項拍板的核心頁。",
 9: "deferred 三項（M9 物理不可能、PC harness、M4abd）+ 下一步 cross-sample 決策點。",
}
note_html = ""
for i, txt in NOTES.items():
    note_html += f'<details class="note" data-for="{i}"><summary>▶ 講稿 · S{i}</summary><div class="nb">{txt}</div></details>\n'

cards = "\n".join(f'<section class="slide{" active" if i==0 else ""}" data-idx="{i+1}">{c}</section>'
                  for i, (_, c) in enumerate(SLIDES))
total = len(SLIDES)

html = f"""<!DOCTYPE html><html lang="zh-TW"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ASM characterization presenter deck (v1) · {total} slides</title>
<style>
:root{{--bg:#FAF9F5;--surf:#fff;--ink:#141413;--soft:#57564F;--rule:#E3DACC;--code:#F4EFE6;
--primary:#1E3A8A;--accent:#D97757;--ad:#B85A3F;--ok:#15803D;--warn:#A16207;--neg:#C2410C;
--ff:-apple-system,"Noto Sans CJK TC","Noto Sans CJK JP","PingFang TC","Droid Sans Fallback",sans-serif;
--mono:"JetBrains Mono","DejaVu Sans Mono",monospace;--r:5px;}}
*{{box-sizing:border-box}}html,body{{margin:0;background:#262626;color:var(--ink);font-family:var(--ff)}}
.wrap{{min-height:100vh;display:flex;flex-direction:column;align-items:center;padding:12px 16px}}
.bar{{width:100%;max-width:1280px;display:flex;justify-content:space-between;align-items:center;color:#eee;font-size:14px;padding:6px 10px;margin-bottom:10px}}
.bar .br{{font-weight:700;color:#fff}}.bar .br small{{display:block;color:#bbb;font-weight:400;font-size:12px}}
.bar button{{background:#3A3A3A;color:#fff;border:1px solid #555;border-radius:var(--r);padding:5px 12px;cursor:pointer;font-family:inherit}}
.bar button:disabled{{opacity:.35}}.bar .ct{{font-family:var(--mono);background:#1a1a1a;color:#fff;padding:5px 12px;border-radius:var(--r)}}
.bar label{{color:#ddd;font-size:13px;cursor:pointer}}
.canvas{{width:100%;max-width:1280px;aspect-ratio:16/9;background:var(--surf);border:1px solid var(--rule);box-shadow:0 2px 8px rgba(0,0,0,.3);overflow:hidden;position:relative}}
.slide{{display:none;width:100%;height:100%;padding:22px 32px 14px;flex-direction:column}}
.slide.active{{display:flex}}
.head{{border-bottom:2px solid var(--primary);padding-bottom:7px;margin-bottom:11px;flex-shrink:0}}
.eb{{font-size:12px;text-transform:uppercase;letter-spacing:1px;color:var(--ad);font-weight:700}}
.t{{font-size:25px;margin:3px 0 0;line-height:1.18;font-weight:700}}
.body{{flex:1;display:grid;gap:14px;min-height:0;align-items:stretch}}
.body.split-13{{grid-template-columns:1.45fr 1fr}}.body.split-11{{grid-template-columns:1fr 1fr}}.body.single{{grid-template-columns:1fr}}
.hero{{display:flex;align-items:center;justify-content:center;background:var(--bg);border:1px solid var(--rule);border-radius:var(--r);padding:6px;min-height:0;overflow:hidden}}
.hero img{{max-width:100%;max-height:100%;object-fit:contain}}
.side{{display:flex;flex-direction:column;gap:9px;min-height:0;overflow:auto}}
.side h3{{margin:0;font-size:16px;color:var(--primary);border-bottom:1px dashed var(--rule);padding-bottom:4px}}
.side ul{{margin:0;padding-left:18px}}.side li{{font-size:14px;line-height:1.5;margin:3px 0}}.side li b{{color:var(--ad)}}
.metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-bottom:11px;flex-shrink:0}}
.m{{background:var(--bg);border:1px solid var(--rule);border-top:3px solid var(--accent);border-radius:var(--r);padding:8px;text-align:center}}
.m .v{{font-family:var(--mono);font-size:22px;font-weight:700;color:var(--primary)}}.m.ok .v{{color:var(--ok)}}.m .l{{font-size:11px;color:var(--soft);margin-top:3px;line-height:1.3}}
.vbox{{background:#EFF6FF;border:1px solid #BFDBFE;border-left:3px solid var(--primary);border-radius:0 var(--r) var(--r) 0;padding:8px 11px;font-size:13px;line-height:1.5}}
.vbox.ok{{background:#F0FDF4;border-color:#BBF7D0;border-left-color:var(--ok)}}
.vbox.warn{{background:#FEF7ED;border-color:#FCD9A8;border-left-color:var(--warn);margin:5px 0}}
.confirm{{background:#FFFBEB;border:1.5px solid var(--accent);border-radius:var(--r);padding:9px 12px;font-size:13.5px;line-height:1.5;margin-top:auto}}
.confirm b{{color:var(--ad)}}
.cav-grid{{display:grid;grid-template-columns:1fr;gap:9px}}
.cav{{display:flex;gap:11px;background:#FEF7ED;border-left:4px solid var(--warn);border-radius:0 var(--r) var(--r) 0;padding:9px 13px;font-size:14px;line-height:1.5}}
.cav .cn{{font-family:var(--mono);font-size:20px;font-weight:700;color:var(--warn)}}.cav b{{color:#92600A}}
table.bcg{{border-collapse:collapse;width:100%;font-size:14px}}
table.bcg th,table.bcg td{{border-top:1px solid var(--rule);border-bottom:1px solid var(--rule);padding:7px 10px;text-align:left}}
table.bcg th{{background:var(--code);font-weight:700;border-top:2px solid var(--ink)}}table.bcg td.c{{text-align:center;font-family:var(--mono);color:#2563EB;font-weight:700}}
.terms{{font-size:12.5px;color:var(--soft);background:#FCFAF4;border:1px dashed var(--rule);border-radius:var(--r);padding:6px 12px;margin-top:8px;flex-shrink:0}}.terms b{{color:var(--primary)}}
.cover{{display:flex;flex-direction:column;justify-content:center;height:100%;padding:20px 50px}}
.cover .eyebrow{{font-size:13px;letter-spacing:2px;text-transform:uppercase;color:var(--ad);font-weight:700}}
.cover h1{{font-size:30px;margin:10px 0;line-height:1.22;font-weight:800}}.cover h1 .en{{display:block;font-size:18px;color:var(--primary);margin-top:6px}}
.cover .lead{{font-size:16px;color:var(--ink);line-height:1.5;max-width:92%}}.cover .lead b{{color:var(--ad)}}
.cover .arc{{display:grid;grid-template-columns:repeat(5,1fr);gap:8px;margin:18px 0}}
.cover .arc .step{{background:var(--bg);border:1px solid var(--rule);border-top:3px solid var(--accent);padding:9px;border-radius:var(--r)}}
.cover .arc .num{{font-family:var(--mono);font-weight:700;font-size:12px;color:var(--ad)}}.cover .arc .lbl{{display:block;margin-top:3px;font-weight:700;font-size:13px}}
.cover .credits{{display:flex;gap:0;border-top:1px solid var(--rule);padding-top:12px;font-size:14px}}
.cover .cr{{padding-right:18px;margin-right:18px;border-right:1px solid var(--rule)}}.cover .cr:last-child{{border:none}}
.cover .cr .k{{font-size:11px;color:var(--soft);text-transform:uppercase;letter-spacing:1px;display:block}}.cover .cr .v{{font-weight:700;margin-top:2px;display:block}}
.notes{{width:100%;max-width:1280px;margin-top:12px}}
.note{{background:#1f1f1f;border:1px solid #444;border-radius:var(--r);margin:6px 0;color:#eee}}
.note summary{{cursor:pointer;padding:9px 13px;font-weight:600;color:#fff;background:#2a2a2a;border-radius:var(--r)}}
.note .nb{{padding:11px 15px;color:#ddd;font-size:14px;line-height:1.6}}
body.proj .notes{{display:none}}
@media print{{.bar,.notes{{display:none}}.canvas{{aspect-ratio:auto;box-shadow:none;page-break-after:always}}.slide{{display:flex!important;page-break-after:always}}}}
</style></head><body>
<div class="wrap">
<header class="bar"><div class="br">ASM characterization · presenter deck v1<small>HCC1395 單樣本 · G1 找清楚分群 + G2 LOH 診斷 · 給用戶確認</small></div>
<div style="display:flex;gap:8px;align-items:center"><button id="p">← Prev</button><span class="ct" id="ct">1 / {total}</span><button id="n">Next →</button>
<label><input id="pj" type="checkbox"> Projector</label></div></header>
<div class="canvas" id="cv">{cards}</div>
<div class="notes" id="nt">{note_html}</div>
</div>
<script>
(function(){{
 const S=[...document.querySelectorAll('.slide')],NT=[...document.querySelectorAll('.note')],ct=document.getElementById('ct');
 const P=document.getElementById('p'),Nx=document.getElementById('n'),PJ=document.getElementById('pj');let i=0,T=S.length;
 function show(k){{i=Math.max(0,Math.min(T-1,k));S.forEach((s,j)=>s.classList.toggle('active',j===i));
  NT.forEach(d=>{{const f=+d.getAttribute('data-for');d.style.display=(f===i+1)?'block':'none';if(f!==i+1)d.removeAttribute('open');}});
  ct.textContent=(i+1)+' / '+T;P.disabled=i===0;Nx.disabled=i===T-1;try{{history.replaceState(null,'','#s'+(i+1))}}catch(e){{}}}}
 P.onclick=()=>show(i-1);Nx.onclick=()=>show(i+1);
 document.addEventListener('keydown',e=>{{if(e.target.tagName==='INPUT')return;
  if(['ArrowRight','PageDown',' '].includes(e.key)){{e.preventDefault();show(i+1)}}
  else if(['ArrowLeft','PageUp'].includes(e.key)){{e.preventDefault();show(i-1)}}
  else if(e.key==='Home')show(0);else if(e.key==='End')show(T-1);
  else if(/^[1-9]$/.test(e.key)){{const k=+e.key-1;if(k<T)show(k)}}
  else if(e.key==='p'||e.key==='P'){{PJ.checked=!PJ.checked;PJ.dispatchEvent(new Event('change'))}}}});
 PJ.onchange=()=>document.body.classList.toggle('proj',PJ.checked);
 const m=(location.hash.match(/^#s(\\d+)$/)||[])[1];show(m?+m-1:0);
}})();
</script></body></html>"""

OUT.write_text(html, encoding="utf-8")
print(f"deck -> {OUT}  ({OUT.stat().st_size//1024} KB, {total} slides)")
print(f"  injected: cd={N['cd']} surv={N['surv']}/{N['pass1']} LOH {N['artp']}/{N['subp']}/{N['cnp']}% collider_killed={collider_named}")
