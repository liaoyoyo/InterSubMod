#!/usr/bin/env python3
"""Build standalone interactive HTML for the ASM x SEQC2 CN confound pilot.
LLM-authored content + base64-inline 4 figures. No external deps at view time."""
import base64, html
from pathlib import Path

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
FIG = ROOT / "figures/cn_confound"
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
           "20260602_ASM_CN_confound_disentanglement_pilot_01.standalone.html")


def b64(p):
    p = Path(p)
    if not p.exists():
        return None
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


FIGS = {k: b64(FIG / f"{k}.png") for k in
        ["q2_dose_response", "q3_reverse_predict", "q4_two_axis", "q5_error_align"]}

# ---- verdict cards (data from p4_synthesis.json / workflow results) ----
HYP = [
    dict(id="H-A", q="Q2 · CN dose-response", badge="CONFIRM", cls="ok",
         claim="HP-axis ASM |Δβ| 獨立於 ground-truth median CN（非 CN-driven）",
         nums=[("HP-axis partial ρ(|Δβ|,CN | cov)", "−0.055", "p=7.7e-5, n=5142 · 反向且 |ρ|≪0.2"),
               ("4 個 coverage strata", "ρ∈[−0.073,−0.053]", "全一致負向 → 非 Simpson's paradox"),
               ("integer-CN median |Δβ|", "0.082→0.062 (CN1→8)", "無單調上升")],
         rule="Falsifier = cov-controlled partial ρ > 0.2 且正向 p<0.05。觀測 ρ=−0.055：方向相反 + 遠低於 0.2 → falsifier 未觸發，H-A supported。",
         note="顯著只因 n 大（power），effect 可忽略。方向與 dose-response artifact 相反（artifact 會在高 CN 灌大 |Δβ|=正 ρ），強化『非 CN gain/loss artifact』。",
         fig="q2_dose_response"),
    dict(id="H-B", q="Q3 · 反向預測倍體", badge="INCONCLUSIVE", cls="warn",
         claim="甲基特徵無法預測 CN（超 coverage）— 落灰帶，唯一未乾淨收斂",
         nums=[("dAUC max (RF / logistic)", "0.046 / 0.009", "0.02 < 0.046 < 0.05 = 灰帶"),
               ("全 absolute AUC", "< 0.58", "RF full 0.566 / cov-only 0.520 / meth-only 0.547"),
               ("meth-only permutation", "AUC 0.549, p=0.001", "弱但非全 null")],
         rule="Falsifier = dAUC>0.05（未觸發）；strict confirm = dAUC≤0.02（RF 未過）→ 灰帶。o2-augmented dAUC=0.144 會 REFUTE 但棄用（非 pre-reg 特徵 + selection bias + overfit）。",
         note="甲基帶弱-但-真的 CN trace（perm p=0.001），實務無用（AUC<0.58）。需 COLO829 複製 + 測 CN→coverage/region vs CN→methylation 路徑才能推出灰帶。",
         fig="q3_reverse_predict"),
    dict(id="H-C", q="Q4 · 兩軸診斷", badge="CONFIRM", cls="ok",
         claim="HP 與 ALLELE 軸一致；ALLELE excess 非 allele-dosage artifact",
         nums=[("OLS excess~CN slope", "−0.00196", "p=4.9e-8 但反向；r²=0.0015 ≈ 平"),
               ("MW gain-vs-neutral excess", "p=0.967", "無 class 差"),
               ("HP-vs-ALLELE |Δβ| 相關", "0.59 / 0.78", "overall / within cnLOH（CN=2 平衡區更一致）")],
         rule="Falsifier = excess 隨 CN 上升（正向 slope p<0.05）。觀測 slope 反向 → falsifier 未觸發。p<0.05 不誤判為 falsification（方向不符）。",
         note="cnLOH 內一致性更高（0.78）是 dosage artifact 的反面，強化 H-C。n=19,937 well-powered。",
         fig="q4_two_axis"),
    dict(id="H-D", q="Q5 · error / alignment", badge="CONFIRM-partial", cls="ok",
         claim="NM/supplementary error 非 ASM 主因；MAPQ 通道無法測（blind spot）",
         nums=[("resid(|Δβ| | cov+CN) vs error", "全 |ρ|≤0.05", "median_nm −0.008 / nm_per_kb 0.052 / supp −0.042, p≥0.21"),
               ("resid model r²(cov+CN)", "0.129", "87% |Δβ| variance 殘餘且與 error 無關"),
               ("median_mapq", "恆=60 (零變異)", "⚠ low-MAPQ confound 無法測，非 refute")],
         rule="Falsifier = 殘差 vs error proxy ρ>0.2 p<0.05。可測 proxy 全 ρ≤0.05 → falsifier 未觸發。但 MAPQ 飽和零變異 → 限縮為『NM/supp error 非主因』。",
         note="結論不可寫成『無 error confound』；MAPQ/alignment-difficulty 通道是已知 blind spot，需有變異的 per-read MAPQ 才能關閉。",
         fig="q5_error_align"),
    dict(id="H-E", q="O2 · 聚類 × CN/error", badge="CONFIRM", cls="ok",
         claim="read-level 甲基聚類=真結構（≫placebo）；對 CN 有弱殘餘相關",
         nums=[("blind-ARI vs placebo (length null)", "0.267 vs 0.015", "paired Wilcoxon p=1.8e-43 · collider_flag=false"),
               ("blind-ARI vs CN / nm_per_kb", "ρ=0.116 / 0.119", "顯著 p<0.005 但 <0.2 falsifier bar"),
               ("canonical BRCA2 chr13:32,315,128", "blind_ari=0.790", "placebo 0.132, n0/n1=34/21 · positive control 過")],
         rule="Falsifier = ARI~CN 或 ~error |ρ|>0.2 p<0.05。無 proxy 越過 0.2 → falsifier 未觸發。Load-bearing = blind≫placebo（p=1.8e-43）。",
         note="非乾淨 null：ARI 對 CN/error 有弱殘餘（ρ≈0.12）。應 foreground placebo 分離（p=1.8e-43）而非裸 null 宣稱。",
         fig=None),
]

NEEDS = [
    "H-B 未解：dAUC=0.046 在灰帶，perm p=0.001 弱非 null CN leakage。需 COLO829 複製 + 區分 CN→coverage/region vs CN→methylation。",
    "H-D MAPQ blind spot：median_mapq 恆 60（零變異）→ low-MAPQ/alignment confound 無法測。結論限縮為『NM/supp error 非主因』。",
    "單樣本天花板：所有 verdict 僅 HCC1395（唯一 SEQC2-CN 樣本）。downstream 全標 partial_flag，不宣稱跨樣本。",
    "Q1 reproducibility scope 窄：raw Level1 disk-safe-deleted；re-audit 只重生 20 chr19 HP1 位點（~0.04% rows）。可考慮多 chrom/ALLELE 軸 spot re-audit。",
    "H-E ARI 非完全正交：ari~CN ρ=0.116、ari~nm/kb ρ=0.119（p<0.005，<0.2 bar）。Foreground blind-vs-placebo Wilcoxon。",
    "跨平台 caveat：SEQC2 median CN 是 short-read NGS，非 ONT gold-standard。附於所有 CN-based verdict。",
]

badge_css = {"ok": "#15803D", "warn": "#B45309", "bad": "#B91C1C"}


def card(h):
    rows = "".join(
        f'<tr><td class="k">{html.escape(k)}</td><td class="v">{html.escape(v)}</td>'
        f'<td class="d">{html.escape(d)}</td></tr>' for k, v, d in h["nums"])
    fig = ""
    if h["fig"] and FIGS.get(h["fig"]):
        fig = (f'<details class="figbox"><summary>圖：{html.escape(h["q"])}</summary>'
               f'<img src="{FIGS[h["fig"]]}" alt="{h["fig"]}"></details>')
    bc = badge_css[h["cls"]]
    return f'''
    <section class="card" id="{h["id"]}">
      <div class="card-h">
        <span class="hid">{h["id"]}</span>
        <span class="qlabel">{html.escape(h["q"])}</span>
        <span class="badge" style="background:{bc}">{h["badge"]}</span>
      </div>
      <p class="claim">{html.escape(h["claim"])}</p>
      <table class="nums"><tbody>{rows}</tbody></table>
      <details class="method"><summary>方法判準 + caveat（點開審查）</summary>
        <p><b>Decision rule：</b>{html.escape(h["rule"])}</p>
        <p><b>Caveat：</b>{html.escape(h["note"])}</p>
      </details>
      {fig}
    </section>'''


cards = "\n".join(card(h) for h in HYP)
needs = "\n".join(f"<li>{html.escape(n)}</li>" for n in NEEDS)
strip = "".join(
    f'<a href="#{h["id"]}" class="chip" style="border-color:{badge_css[h["cls"]]};color:{badge_css[h["cls"]]}">'
    f'{h["id"]} · {h["badge"]}</a>' for h in HYP)

HTML = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ASM × SEQC2 CN confound pilot — HCC1395</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--card:#fff;--accent:#1E3A8A}}
*{{box-sizing:border-box}}
body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;
color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:980px;margin:0 auto;padding:24px 18px 80px}}
header.top{{background:linear-gradient(135deg,#1E3A8A,#0f1f4d);color:#fff;border-radius:14px;padding:26px 28px;margin-bottom:18px}}
header.top h1{{margin:0 0 6px;font-size:1.5rem}}
header.top .meta{{font-size:.84rem;opacity:.85}}
.tldr{{background:#ecfdf5;border:1px solid #a7f3d0;border-left:5px solid #15803D;border-radius:12px;padding:18px 20px;margin-bottom:16px}}
.tldr h2{{margin:0 0 8px;font-size:1.06rem;color:#065f46}}
.tldr p{{margin:.3rem 0}}
.strip{{display:flex;flex-wrap:wrap;gap:8px;margin:14px 0 4px}}
.chip{{text-decoration:none;font-weight:600;font-size:.82rem;border:1.5px solid;border-radius:999px;padding:4px 12px;background:#fff}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:14px 0;scroll-margin-top:14px}}
.card-h{{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:6px}}
.hid{{font-weight:800;font-size:1.05rem;color:var(--accent)}}
.qlabel{{color:var(--mut);font-size:.9rem;flex:1}}
.badge{{color:#fff;font-weight:700;font-size:.78rem;padding:3px 11px;border-radius:999px;letter-spacing:.3px}}
.claim{{font-weight:600;margin:.2rem 0 .7rem}}
table.nums{{width:100%;border-collapse:collapse;font-size:.86rem}}
table.nums td{{padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}}
table.nums td.k{{color:var(--mut);width:34%}}
table.nums td.v{{font-weight:700;font-variant-numeric:tabular-nums;width:22%;color:#0f172a}}
table.nums td.d{{color:var(--mut);font-size:.82rem}}
details{{margin-top:10px}}
details.method>summary{{cursor:pointer;color:var(--accent);font-weight:600;font-size:.85rem}}
details.method p{{font-size:.85rem;background:var(--bg);border-radius:8px;padding:8px 12px;margin:.5rem 0}}
details.figbox>summary{{cursor:pointer;color:#374151;font-weight:600;font-size:.85rem;margin-top:8px}}
details.figbox img{{width:100%;border:1px solid var(--line);border-radius:8px;margin-top:8px}}
.sec-h{{font-size:1.15rem;margin:26px 0 6px;border-bottom:2px solid var(--accent);padding-bottom:4px}}
.forest{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:10px 0}}
.fc{{background:#fff;border:1px solid var(--line);border-radius:10px;padding:12px 14px;text-align:center}}
.fc .n{{font-size:1.3rem;font-weight:800;color:var(--accent)}}
.fc .l{{font-size:.78rem;color:var(--mut)}}
.dag{{background:#fff;border:1px solid var(--line);border-radius:10px;padding:14px 16px;font-size:.86rem}}
.warnbox{{background:#fffbeb;border:1px solid #fde68a;border-left:5px solid #B45309;border-radius:10px;padding:14px 18px;margin:12px 0}}
.warnbox ul{{margin:.4rem 0;padding-left:1.2rem}}.warnbox li{{margin:.3rem 0;font-size:.86rem}}
footer{{margin-top:30px;font-size:.78rem;color:var(--mut);border-top:1px solid var(--line);padding-top:14px}}
code{{background:#eef2ff;padding:1px 6px;border-radius:5px;font-size:.82em}}
</style></head><body><div class="wrap">

<header class="top">
  <h1>甲基分群差異是 mutation-linked，還是 CN / ploidy / LOH / error / alignment confound？</h1>
  <div class="meta">ASM × SEQC2 連續 median CN ground-truth · confound disentanglement pilot
  &nbsp;|&nbsp; 樣本 <b>HCC1395</b>（唯一有 SEQC2 CN truth；single-sample A pilot）
  &nbsp;|&nbsp; 2026-06-02 &nbsp;|&nbsp; anchor <code>8a5faa8</code> · workflow <code>wf_8fe183fa-b09</code></div>
</header>

<div class="tldr">
  <h2>TL;DR — 在可測通道上，confound 被大致排除</h2>
  <p><b>甲基（ASM Δβ）與 read-level 聚類的分群差異主要是 mutation-linked</b>，而非 copy-number / ploidy / LOH / error-read / alignment 的產物。</p>
  <p>機制鑰匙：<b>HP-axis 設計上把 CN/ploidy/alignment/region held constant</b>（同位點同單倍型 germline-read vs somatic-read）→ de-confounded estimator。Q1 方法 gate 通過（β 重算 20/20 完全重現；聚類驗尺 PC1=0.557/NC1=−0.005）。</p>
  <p><b>3 個誠實限制</b>：① 單樣本天花板（不可跨樣本一般化）② H-B 反向預測落灰帶 + MAPQ 通道無法測 ③ SEQC2 CN 是 short-read NGS（非 ONT gold-standard）。對齊既有 ZAR1L/BRCA2 ⭐3（真實但需 COLO829）。</p>
  <div class="strip">{strip}</div>
</div>

<h2 class="sec-h">見林 — 資料層（HCC1395 高度擴增基因組）</h2>
<div class="forest">
  <div class="fc"><div class="n">56,320</div><div class="l">per-site 列（TP 51,171 / FP 5,149）</div></div>
  <div class="fc"><div class="n">77.8%</div><div class="l">gain (CN&gt;2) · neutral 21.3% · loss 0.9%</div></div>
  <div class="fc"><div class="n">10,342</div><div class="l">cnLOH（copy-neutral, CN=2）</div></div>
  <div class="fc"><div class="n">2.85</div><div class="l">genome ploidy（SEQC2/KB）</div></div>
  <div class="fc"><div class="n">596</div><div class="l">BAM-pass loci（O2 聚類 + Q5 error）</div></div>
</div>
<div class="dag"><b>Confound DAG：</b> <code>MUT → METH</code>（目標因果邊）；confounders <code>{{CN/ploidy, LOH, coverage, error-read, alignment}} → METH</code>；結構鏈 <code>CN → alignment 難度 → error-read</code>。HP-axis 為唯一同時 held-constant CN/ploidy/alignment 的乾淨估計子，4 個 CONFIRM 全建於此軸。</div>

<h2 class="sec-h">見樹 — 五假說逐條 verdict（對照 pre-reg falsifier）</h2>
{cards}

<h2 class="sec-h">誠實限制 / needs-work（封鎖 tier 升級）</h2>
<div class="warnbox"><ul>{needs}</ul></div>

<h2 class="sec-h">Next step（不自動啟動）</h2>
<div class="dag">
(a) <b>COLO829 ground-truth 複製</b> → H-A/B/C 跨樣本（需查 COLO829 是否有 SEQC2-style CN truth）；
(b) <b>H-B 收尾</b> — 測 CN→coverage/region vs CN→methylation 路徑把灰帶推出；
(c) <b>MAPQ blind spot</b> — 需 source 有變異的 per-read MAPQ；
(d) 回 <b>phasing 主軸</b>，本 pilot 當 ASM 真實性的 methods 支撐。
</div>

<footer>
Provenance（scientific-rigor §8.4）：anchor commit <code>8a5faa8</code> · workflow <code>wf_8fe183fa-b09</code>（8 agents）·
pre-reg + 結果 <code>InterSubMod/docs/experiments/in_progress/2026/06/20260602_ASM_CN_confound_disentanglement_pilot_01.md</code> ·
scripts <code>40_cn_annotate.py</code>–<code>46_q5_error_align.py</code> ·
evaluator synthesis <code>genome_survey_v2/cn_confound/p4_synthesis.json</code> ·
master tables <code>master_o1_cn.tsv</code> / <code>master_o2_error.tsv</code>。
Single-sample exploratory ⭐3；所有數值可追溯至上述 JSON / TSV。圖為 workflow agent matplotlib 產出（base64 內嵌）。
</footer>
</div></body></html>'''

OUT.write_text(HTML, encoding="utf-8")
sz = OUT.stat().st_size // 1024
nfig = sum(1 for v in FIGS.values() if v)
print(f"[html] {OUT} ({sz} KB, {nfig}/4 figs inlined)")
