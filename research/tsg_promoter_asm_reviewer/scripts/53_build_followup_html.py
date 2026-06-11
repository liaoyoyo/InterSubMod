#!/usr/bin/env python3
"""Standalone HTML for the ASM follow-up (count fix + context TP/FP + measurement + subclone)."""
import base64, html
from pathlib import Path

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
FIG = ROOT / "figures/cn_confound"
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
           "20260602_ASM_CN_confound_disentanglement_pilot_01_followup.standalone.html")


def b64(p):
    p = Path(p)
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode() if p.exists() else None


FIGS = {k: b64(FIG / f"{k}.png") for k in
        ["m1_perpos_reanalysis", "m2_meth_measurement", "c1_context_tpfp_full",
         "c2_context_tpfp_clustering", "s_subclone_assessment"]}

SECTIONS = [
    dict(id="count", q="Q1 · count 校正 + pseudoreplication（M1）", badge="已校正 · H-A robust", cls="ok",
         claim="56,320 是 records 非位點；pseudoreplication 不影響 H-A（反而更穩）",
         nums=[("56,320 records → 唯一位點", "34,154", "TP 30,511 + FP 3,643"),
               ("SEQC2 truth / 分析 TP", "39,447 / 30,490", "ClairS-TO recall 77%；先前『全39,447 TP』是 mislabel"),
               ("pseudoreplication: per-record→per-position ρ", "−0.055 → −0.039", "inflation 1.017× · 只 86/5056 位點雙軸顯著"),
               ("H-A 結論", "HOLDS（effect 縮小）", "ASM 非 CN-driven 對 pseudoreplication robust")],
         note="只有 86/5056 HP-sig 位點同時 HP1+HP2 顯著 → 先前 per-record 統計幾乎沒膨脹；H-A verdict 穩固。",
         fig="m1_perpos_reanalysis"),
    dict(id="context", q="Q2 · 倍體是否影響甲基分群（C1）", badge="否（分兩層）", cls="ok",
         claim="倍體不驅動甲基分群結構；但 CN 強烈影響 caller TP/FP rate（兩回事，不可混）",
         nums=[("甲基 |Δβ| / clustering vs CN", "CN-independent", "H-A partial ρ=−0.039；H-E blind≫placebo p=1.8e-43"),
               ("caller TP-rate vs CN（logistic）", "強相關", "LOH OR=0.36 / loss OR=0.12 / gain OR=1.57 / cov OR=1.004"),
               ("最差 / 最佳 cell", "LOH.loss 0.36 / nonLOH.gain 0.951", "n=239 / n=18,321"),
               ("LOH 內 abnormal-CN vs cnLOH FP", "OR=0.892 (p=0.009)", "abnormal NOT more FP（反向小，1.7pp）")],
         note="關鍵區分：『倍體不驅動甲基分群』成立（clustering 結構 CN-independent），但 CN-class 確實移動 caller 的 TP/FP 精度 —— 兩個 claim 必須分開講。",
         fig="c1_context_tpfp_full"),
    dict(id="hyp", q="Q2 · 你的兩個具體假設（C2）", badge="都不成立 · 其一反轉", cls="warn",
         claim="H-ctx1 REFUTED+反轉（LOH 內高 clustering 反而 TP）；H-ctx2 UNDERPOWERED",
         nums=[("H-ctx1: LOH∩異常CN∩clust→FP", "REFUTED", "cell n=84(FP=9) TP-rate 89.3% vs 93.4%, OR=0.59 p=0.18 ns"),
               ("★ 反轉發現", "LOH 內 high-clust = TP-enriched", "OR=2.385 p=0.019；TP blind_ari 0.258 > FP 0.168 (MW p=0.0075)"),
               ("→ FP 實際落在", "LOW-clustering LOH-gain", "與直覺相反"),
               ("H-ctx2: 非LOH∩高cov∩clust→TP", "UNDERPOWERED", "nonLOH 全域僅 2 FP → 結構上無法證")],
         note="⚠ 不可據此建『high blind_ari = FP filter』—— 機制與你的直覺相反，且 FP 僅 43/596 嚴重 underpower。auc-confound-guard：這是 context 效應非 classifier。",
         fig="c2_context_tpfp_clustering"),
    dict(id="meas", q="Q3 · 甲基測量法 scrutiny（M2）", badge="targeted 升級", cls="warn",
         claim="非全域偏差（r=0.998），但特定位點 lossy：5hmC 被丟 + intermediate 被壓平",
         nums=[("binarized vs 連續 β", "r=0.998", "全域不偏；dead-zone 僅 12.5%"),
               ("MAX-collapse 與 5hmC 相關", "corr=0.033（正交）", "5hmC 通道被完全丟棄；5hmC>5mC 佔 12%"),
               ("binarization 壓平 intermediate", "BRCA2 連續 0.162 vs bin 0.072", "15% intermediate reads 被壓 0/1"),
               ("建議", "分離連續 5mC + 5hmC", "5mC-only canonical；高 intermediate 位點用連續 β")],
         note="全域 max-collapse binarization 安全（average r=0.998），但兩機制在特定位點藏訊號 → targeted（非全域）升級即可。",
         fig="m2_meth_measurement"),
    dict(id="sub", q="Q3 · Subclone（S）", badge="候選未確認", cls="warn",
         claim="19 個 subclone-like 候選，0 個 blind_ari 佐證；BRCA2 unimodal=clean ASM",
         nums=[("per-read 雙峰 loci", "43/247 (17.4%)", "但 24/43 = HP-split（germ vs som）"),
               ("subclone-like / blind_ari 佐證", "19 / 0", "M2 雙峰與 C2 clustering 抽的 loci 幾乎不重疊"),
               ("BRCA2 chr13:32,315,128", "UNIMODAL (dBIC=−6.6)", "均勻 somatic hypomethylation = 乾淨 ASM 非 subclone"),
               ("measurement_hides_subclone", "YES（特定位點）", "subclone claim 卡在第二樣本 COLO829")],
         note="dominant per-read 結構 = somatic HP split，非新 minor clone。subclone 需 COLO829 複製 + M2/C2 共用 loci 重測才能宣稱。",
         fig="s_subclone_assessment"),
]

badge_css = {"ok": "#15803D", "warn": "#B45309", "bad": "#B91C1C"}


def card(h):
    rows = "".join(f'<tr><td class="k">{html.escape(k)}</td><td class="v">{html.escape(v)}</td>'
                   f'<td class="d">{html.escape(d)}</td></tr>' for k, v, d in h["nums"])
    fig = (f'<details class="figbox" open><summary>圖：{html.escape(h["q"])}</summary>'
           f'<img src="{FIGS[h["fig"]]}"></details>') if FIGS.get(h["fig"]) else ""
    bc = badge_css[h["cls"]]
    return f'''<section class="card" id="{h["id"]}"><div class="card-h">
      <span class="qlabel">{html.escape(h["q"])}</span>
      <span class="badge" style="background:{bc}">{html.escape(h["badge"])}</span></div>
      <p class="claim">{html.escape(h["claim"])}</p>
      <table class="nums"><tbody>{rows}</tbody></table>
      <details class="method"><summary>caveat / 讀法</summary><p>{html.escape(h["note"])}</p></details>
      {fig}</section>'''


cards = "\n".join(card(h) for h in SECTIONS)
strip = "".join(f'<a href="#{h["id"]}" class="chip" style="border-color:{badge_css[h["cls"]]};'
                f'color:{badge_css[h["cls"]]}">{html.escape(h["q"].split(" · ")[1])}</a>' for h in SECTIONS)

HTML = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ASM follow-up — count + context + measurement + subclone</title><style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1000px;margin:0 auto;padding:24px 18px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:24px 26px;margin-bottom:16px}}
header.top h1{{margin:0 0 6px;font-size:1.4rem}}.meta{{font-size:.82rem;opacity:.85}}
.tldr{{background:#eef2ff;border:1px solid #c7d2fe;border-left:5px solid #1E3A8A;border-radius:12px;padding:16px 20px;margin-bottom:14px}}
.tldr h2{{margin:0 0 8px;font-size:1.05rem;color:#1e3a8a}}.tldr p{{margin:.35rem 0}}
.strip{{display:flex;flex-wrap:wrap;gap:8px;margin-top:12px}}
.chip{{text-decoration:none;font-weight:600;font-size:.78rem;border:1.5px solid;border-radius:999px;padding:4px 11px;background:#fff}}
.card{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:14px 0;scroll-margin-top:14px}}
.card-h{{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:6px}}
.qlabel{{font-weight:800;color:var(--accent);flex:1;font-size:1rem}}
.badge{{color:#fff;font-weight:700;font-size:.76rem;padding:3px 11px;border-radius:999px}}
.claim{{font-weight:600;margin:.2rem 0 .7rem}}
table.nums{{width:100%;border-collapse:collapse;font-size:.85rem}}
table.nums td{{padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}}
table.nums td.k{{color:var(--mut);width:32%}}td.v{{font-weight:700;width:26%;color:#0f172a;font-variant-numeric:tabular-nums}}td.d{{color:var(--mut);font-size:.8rem}}
details{{margin-top:10px}}details.method>summary{{cursor:pointer;color:var(--accent);font-weight:600;font-size:.84rem}}
details.method p{{font-size:.84rem;background:var(--bg);border-radius:8px;padding:8px 12px}}
details.figbox>summary{{cursor:pointer;color:#374151;font-weight:600;font-size:.84rem;margin-top:8px}}
details.figbox img{{width:100%;border:1px solid var(--line);border-radius:8px;margin-top:8px}}
.warnbox{{background:#fffbeb;border:1px solid #fde68a;border-left:5px solid #B45309;border-radius:10px;padding:14px 18px;margin:14px 0}}
.warnbox ul{{margin:.3rem 0;padding-left:1.2rem}}.warnbox li{{margin:.3rem 0;font-size:.85rem}}
footer{{margin-top:28px;font-size:.77rem;color:var(--mut);border-top:1px solid var(--line);padding-top:12px}}
code{{background:#eef2ff;padding:1px 6px;border-radius:5px;font-size:.82em}}
</style></head><body><div class="wrap">
<header class="top"><h1>ASM follow-up — count 校正 · 倍體×甲基分群 · 測量法 · subclone</h1>
<div class="meta">HCC1395 single-sample A pilot · 2026-06-02 · workflow <code>wf_97d18c3b-42a</code>（6 agents）· 接續 <code>wf_8fe183fa-b09</code></div></header>

<div class="tldr"><h2>TL;DR — 三問答案</h2>
<p>① <b>Count</b>：56,320 是 <b>records 非位點</b>（唯一位點 34,154）；先前「全 39,447 TP」是把 SEQC2 truth 誤當分析 TP（實際分析 30,490 ClairS-TO calls）。pseudoreplication 修正後 <b>H-A 仍成立且更穩</b>。</p>
<p>② <b>倍體是否影響甲基分群</b>：<b>不驅動甲基分群結構</b>（clustering CN-independent）—— 但 CN 確實影響 caller TP/FP rate（兩層分開）。你的兩個假設都不成立，且 <b>H-ctx1 反轉</b>：LOH 內<b>高</b> clustering 反而是 <b>TP</b>，FP 落在<b>低</b> clustering。</p>
<p>③ <b>甲基測量法</b>：非全域偏差（r=0.998），但<b>特定位點 lossy</b>（5hmC 被 max-collapse 丟、intermediate 被 binarization 壓平）→ 建議 <b>targeted 升級為分離連續 5mC+5hmC</b>。<b>subclone 有 19 候選但 0 佐證</b>（BRCA2 是 unimodal 乾淨 ASM），需 COLO829 複製。</p>
<div class="strip">{strip}</div></div>

{cards}

<h2 style="font-size:1.1rem;border-bottom:2px solid var(--accent);padding-bottom:4px;margin-top:24px">誠實限制 / needs-work</h2>
<div class="warnbox"><ul>
<li>「倍體不驅動甲基分群」≠「CN 無影響」—— CN 確實移動 caller TP/FP rate（gain OR=1.57 / loss OR=0.12）；CN-independence 專指甲基 |Δβ|/blind_ari 結構。</li>
<li>H-ctx1 是 <b>非顯著 trend</b>（p=0.18；ARI≥0.3 時 0.088），9 個 in-cell FP 故是真 ns 非 underpower；並注意<b>反轉</b>（高 clustering=TP-enriched）。</li>
<li>H-ctx2 <b>結構性 underpower</b>（cell FP=0；nonLOH 全域僅 2 FP）；100% TP-rate 不可當證據。</li>
<li>Subclone <b>卡第二樣本</b>：19 候選 0 個 blind_ari 佐證（M2 雙峰與 C2 clustering 抽的 loci 幾乎不重疊）；HCC1395 單樣本不可升 hypothesis tier 以上。</li>
<li>30,490（caller TP）vs 30,511（ASM-eligible 位點）差 21 —— 已標明，不靜默平均。</li>
<li>單樣本天花板：三問皆 HCC1395-only；tier ≤3 直到 COLO829 複製。</li>
</ul></div>

<h2 style="font-size:1.1rem;border-bottom:2px solid var(--accent);padding-bottom:4px;margin-top:24px">Next step（不自動啟動）</h2>
<div class="warnbox" style="background:#f0fdf4;border-color:#bbf7d0;border-left-color:#15803D"><ul>
<li>(a) <b>targeted 測量升級</b>：~10 top-divergent 位點 + 19 subclone 候選改用<b>分離連續 5mC+5hmC β</b>，per-read GMM 與 blind_ari <b>共用 loci</b> 重測。</li>
<li>(b) <b>COLO829 複製</b>：subclone 候選 + H-A/context 跨樣本 → 才可升 ⭐4 / 宣稱 subclone。</li>
<li>(c) 回 <b>phasing 主軸</b>，本 pilot（ASM 真實 + CN/error confound 排除 + 測量法已 scrutinize）當 methods 支撐。</li>
</ul></div>

<footer>Provenance（§8.4）：workflow <code>wf_97d18c3b-42a</code>（6 agents）· scripts <code>48_perpos_reanalysis.py</code>–<code>52_subclone_assessment.py</code> · evaluator <code>p4b_synthesis.json</code> · per-position master <code>master_perpos.tsv</code> · 主 doc <code>InterSubMod/docs/experiments/in_progress/2026/06/20260602_ASM_CN_confound_disentanglement_pilot_01.md §8</code>。Single-sample exploratory ⭐3；所有數值可追溯至 JSON/TSV。</footer>
</div></body></html>'''

OUT.write_text(HTML, encoding="utf-8")
print(f"[html] {OUT} ({OUT.stat().st_size//1024} KB, {sum(1 for v in FIGS.values() if v)}/5 figs)")
