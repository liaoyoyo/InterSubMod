#!/usr/bin/env python3
"""Standalone HTML for the modkit trial (install+crossval gate, GMM re-test, tool comparison)."""
import base64, html
from pathlib import Path

ROOT = Path("/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer")
FIG = ROOT / "figures/cn_confound"
OUT = Path("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
           "20260602_ASM_CN_confound_disentanglement_pilot_01_modkit_trial.standalone.html")


def b64(p):
    p = Path(p)
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode() if p.exists() else None


FIGS = {k: b64(FIG / f"{k}.png") for k in ["modkit_crossval", "modkit_gmm_retest"]}

CARDS = [
    dict(id="gate", q="① 安裝 + 對照驗證 gate", badge="PASS · AGREE", cls="ok",
         claim="modkit 跑我們 modBAM 與 pysam 近乎完美一致；?-mode 正確",
         nums=[("BRCA2 per-read 5mC Pearson", "1.0", "0.99999998；aggregate Pearson=Spearman=1.0"),
               ("call-set 一致", "shared 13,805 / 0 / 0", "modkit_only=0, pysam_only=0, match_frac=1.0"),
               ("?-mode（關鍵風險）", "0 / 45,332 inferred", "modkit 不把 skipped/unknown 當 canonical → 驗證關閉"),
               ("操作要點", "mod_qual∈[0,1] · BED6", "非 0-255；modkit 拒 BED4；m/h 預設分開 row")],
         note="max diff 0.00195 = ONT 8-bit ML 量化，非真分歧。modkit = 可信賴、可互換的 per-read 5mC extractor（單樣本）。",
         fig="modkit_crossval"),
    dict(id="gmm", q="② GMM 重測（modkit 連續分離 β）", badge="NEGATIVE/NULL", cls="warn",
         claim="更乾淨的測量沒改變 subclone 結論 —— 否證「MAX-collapse 藏了 5hmC 亞克隆」",
         nums=[("13 loci bimodal", "9 confirmed", "3 LOST（全小 n 邊界）、0 GAINED"),
               ("5hmC 新結構", "0", "5hmC 均勻低 0.06-0.10，與 M2 corr=0.033 一致"),
               ("3 個 lost", "chr2/chr7/chr9", "n=28-68，dBIC 9.30/8.25/7.40 剛跌破 ≥10"),
               ("BRCA2", "兩法皆 UNIMODAL", "modkit 連續 minor weight 0.023 = outlier tail = 乾淨 ASM")],
         note="更乾淨測量反而更保守（MAX-collapse 在邊界輕微膨脹 BIC）。pysam MAX-collapse 本來就沒扭曲 subclone 圖像（M2 binarized-vs-continuous r=0.998 已預示）；subclone 卡跨樣本不是測量法。",
         fig="modkit_gmm_retest"),
]

badge_css = {"ok": "#15803D", "warn": "#B45309"}

COMPARE_ROWS = [
    ("角色", "modBAM → per-read 連續抽取 / pileup", "somatic ASM 軸 (HP+ALLELE, tumor+normal)", "per-read subclone 聚類"),
    ("5mC/5hmC", "✅ 預設分開 (m/h rows)", "都帶，但下游 MAX-collapse 有 bug", "5mCG-pattern 為主"),
    ("somatic-aware", "❌ (dmr 僅 region-level)", "✅ 核心", "partial (TP/FP + HP context)"),
    ("subclone-aware", "❌", "partial", "✅ 唯一專做"),
    ("?-mode modBAM", "✅ spec-compliant", "依 modBAM", "依 modBAM"),
    ("可引用 / reviewer 信任", "✅✅ ONT 官方", "internal", "internal"),
    ("速度", "快 (Rust 多執行緒)", "C++ 多執行緒", "C++ O(reads²) 較重"),
    ("我們任務 fit", "extraction 前端 (分離連續 m/h)", "somatic 軸定義來源 (+ collapse bug)", "subclone 聚類本家"),
]


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
      <details class="method"><summary>讀法 / caveat</summary><p>{html.escape(h["note"])}</p></details>
      {fig}</section>'''


cards = "\n".join(card(h) for h in CARDS)
cmp_rows = "".join(
    f'<tr><td class="cdim">{html.escape(d)}</td><td>{html.escape(mk)}</td><td>{html.escape(ms)}</td><td>{html.escape(im)}</td></tr>'
    for d, mk, ms, im in COMPARE_ROWS)

HTML = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>modkit trial — HCC1395</title><style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1000px;margin:0 auto;padding:24px 18px 70px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:24px 26px;margin-bottom:16px}}
header.top h1{{margin:0 0 6px;font-size:1.38rem}}.meta{{font-size:.82rem;opacity:.85}}
.tldr{{background:#eef2ff;border:1px solid #c7d2fe;border-left:5px solid #1E3A8A;border-radius:12px;padding:16px 20px;margin-bottom:14px}}
.tldr h2{{margin:0 0 8px;font-size:1.05rem;color:#1e3a8a}}.tldr p{{margin:.35rem 0}}
.card{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:14px 0}}
.card-h{{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:6px}}
.qlabel{{font-weight:800;color:var(--accent);flex:1;font-size:1rem}}
.badge{{color:#fff;font-weight:700;font-size:.76rem;padding:3px 11px;border-radius:999px}}
.claim{{font-weight:600;margin:.2rem 0 .7rem}}
table.nums{{width:100%;border-collapse:collapse;font-size:.85rem}}table.nums td{{padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}}
table.nums td.k{{color:var(--mut);width:30%}}td.v{{font-weight:700;width:26%;color:#0f172a;font-variant-numeric:tabular-nums}}td.d{{color:var(--mut);font-size:.8rem}}
details{{margin-top:10px}}details.method>summary{{cursor:pointer;color:var(--accent);font-weight:600;font-size:.84rem}}
details.method p{{font-size:.84rem;background:var(--bg);border-radius:8px;padding:8px 12px}}
details.figbox>summary{{cursor:pointer;color:#374151;font-weight:600;font-size:.84rem;margin-top:8px}}details.figbox img{{width:100%;border:1px solid var(--line);border-radius:8px;margin-top:8px}}
.sec-h{{font-size:1.12rem;margin:24px 0 6px;border-bottom:2px solid var(--accent);padding-bottom:4px}}
table.cmp{{width:100%;border-collapse:collapse;font-size:.84rem;background:#fff;border:1px solid var(--line);border-radius:10px;overflow:hidden}}
table.cmp th{{background:#1E3A8A;color:#fff;padding:8px 10px;text-align:left;font-size:.82rem}}
table.cmp td{{padding:7px 10px;border-bottom:1px solid var(--line)}}td.cdim{{font-weight:700;color:#0f172a;width:18%}}
.arch{{background:#f0fdf4;border:1px solid #bbf7d0;border-left:5px solid #15803D;border-radius:10px;padding:14px 18px;margin:12px 0;font-size:.9rem}}
.warnbox{{background:#fffbeb;border:1px solid #fde68a;border-left:5px solid #B45309;border-radius:10px;padding:12px 18px;margin:12px 0}}
.warnbox ul{{margin:.3rem 0;padding-left:1.2rem}}.warnbox li{{margin:.3rem 0;font-size:.85rem}}
footer{{margin-top:26px;font-size:.77rem;color:var(--mut);border-top:1px solid var(--line);padding-top:12px}}
code{{background:#eef2ff;padding:1px 6px;border-radius:5px;font-size:.82em}}
</style></head><body><div class="wrap">
<header class="top"><h1>modkit trial — 測量法升級可行性 + 工具定位</h1>
<div class="meta">ONT modkit v0.6.3 · HCC1395 single-sample A pilot · 2026-06-02 · workflow <code>wf_5f8faad4-a08</code>（4 agents）· 裝於 <code>/big7_disk/liaoyoyo2001/modkit/</code></div></header>

<div class="tldr"><h2>TL;DR</h2>
<p>① <b>裝得起 + 對照過</b>：modkit 跑我們 modBAM 與 pysam <b>per-read 5mC Pearson=1.0</b>、call 完全一致、<b>?-mode 正確</b>（0/45332 inferred）→ 可信賴可互換的 extractor。</p>
<p>② <b>vs MSA/ISM</b>：modkit 是<b>互補的 extraction 前端</b>（官方/可引用/分離 m·h/快），<b>不取代</b> —— 它無 somatic 軸、無 subclone 聚類；那是 MSA / ISM 的領域。</p>
<p>③ <b>能否解決</b>：modkit 能供分離連續抽取，<b>但 GMM 重測沒改變 subclone 結論</b>（9 confirmed / 3 lost 小n / 0 gained / <b>0 個 5hmC 新結構</b>）→ pysam MAX-collapse 本來就沒扭曲；subclone 卡<b>跨樣本</b>不是測量法。modkit 真正 ROI = reviewer-grade provenance + <b>退役 MSA collapse bug</b>。</p></div>

{cards}

<h2 class="sec-h">③ modkit vs MSA vs ISM</h2>
<table class="cmp"><thead><tr><th>維度</th><th>modkit v0.6.3</th><th>MSA</th><th>ISM (本專案)</th></tr></thead>
<tbody>{cmp_rows}</tbody></table>
<div class="arch"><b>建議架構（evaluator 認可）</b>：<code>modkit extract（前端）→ MSA 定 somatic 軸 → ISM 做 subclone 聚類</code>。modkit 只取代「原始 MM/ML 測量步」+ <b>順便退役 MSA 的 5mC+5hmC double-row MAX-collapse bug</b>（−0.054 砍半 artifact 根因）。</div>

<h2 class="sec-h">誠實限制 / next</h2>
<div class="warnbox"><ul>
<li>單樣本（HCC1395）；subclone / 升 ⭐4 需 <b>COLO829 跨樣本</b>複製。</li>
<li>GMM 3 個 flip 全小 n（n=28-68）剛跌破 dBIC≥10 → 報 <b>bootstrap/CI</b> 而非硬 binary，避免 n-driven flip 被當生物訊號。</li>
<li>GMM 重測當 <b>NEGATIVE/NULL 確認</b> —— 不可僅憑測量法 reopen/re-grade subclone 候選。</li>
<li><b>count provenance</b>：19 = genome-wide subclone-like（M2，43 bimodal 中）；12 = s_subclone q3 array top 候選；13 = 12 + BRCA2 實測。勿混。</li>
<li>modkit 尚未進 KB（/big8 shared）；建議 pipeline-manifest 接好後補 KB tool doc。</li>
</ul></div>
<div class="arch"><b>next（不自動啟動）</b>：(a) 把 modkit-extract 寫進 pipeline manifest（附 cross-val log）；(b) COLO829 跨樣本 + dBIC bootstrap/CI。</div>

<footer>Provenance（§8.4）：workflow <code>wf_5f8faad4-a08</code>（4 agents）· modkit v0.6.3 <code>/big7_disk/liaoyoyo2001/modkit/</code> · scripts <code>54_modkit_extract_crossval.py</code> + <code>55_modkit_gmm_retest.py</code> · <code>modkit_trial/{{crossval,gmm_retest,tool_comparison,modkit_trial_synthesis}}.json</code> · 主 doc <code>InterSubMod/docs/experiments/in_progress/2026/06/20260602_ASM_CN_confound_disentanglement_pilot_01.md §9</code>。Single-sample A pilot；數值可追溯至 JSON/TSV。</footer>
</div></body></html>'''

OUT.write_text(HTML, encoding="utf-8")
print(f"[html] {OUT} ({OUT.stat().st_size//1024} KB, {sum(1 for v in FIGS.values() if v)}/2 figs)")
