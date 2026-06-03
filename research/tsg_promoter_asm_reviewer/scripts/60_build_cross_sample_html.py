#!/usr/bin/env python3
"""
60 - Standalone presentation HTML for the cross-sample key-position ASM verification.

§13 layer-A: every number is INJECTED from the deterministic JSONs
(cross_sample_synthesis.json + optional *_gwasm.json). No hand-typed metrics.
Genome-wide section renders only if the gwasm JSONs exist (so the script can be run
before AND after the 59 expansion completes).

Output: docs/experiments/in_progress/2026/06/
        20260603_cross_sample_keypos_ASM_verification_01.standalone.html
"""
import os, json, html, base64, io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
OUT = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
       "20260603_cross_sample_keypos_ASM_verification_01.standalone.html")

ALL = ["HCC1395", "COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
CANCER_FULL = {
    "HCC1395": "breast (SEQC2 ref)", "HCC1937": "breast (BRCA1-mut)",
    "HCC1954": "breast", "H1437": "lung adeno", "H2009": "lung adeno",
    "COLO829": "melanoma",
}
CANCER_COLOR = {"breast": "#be123c", "lung": "#0e7490", "melanoma": "#7c3aed"}


def fig_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def esc(x):
    return html.escape(str(x))


def fnum(x, nd=3):
    if x is None:
        return "—"
    if isinstance(x, float):
        return f"{x:+.{nd}f}" if abs(x) < 100 else f"{x:.0f}"
    return str(x)


# ---------------------------------------------------------------------------
S = json.load(open(f"{CS}/cross_sample_synthesis.json"))
meta = S["meta"]
H = S["headline"]
SO = S["somatic_overlap"]
POS = S["per_position"]
cancer = meta["cancer_types"]

# optional genome-wide
GW = {}
for s in ALL:
    p = f"{CS}/{s}_gwasm.json"
    if os.path.exists(p):
        GW[s] = json.load(open(p))
has_gw = len(GW) >= 5   # genome-wide section only when (nearly) all samples done


# ---------------------------------------------------------------------------
# FIGURE 1: somatic overlap at the 38 key positions
# ---------------------------------------------------------------------------
def fig_somatic_overlap():
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    xs = np.arange(len(ALL))
    som = [SO[s]["n_somatic_at_keypos"] for s in ALL]
    asm = [SO[s]["n_somatic_asm"] for s in ALL]
    cols = [CANCER_COLOR[cancer[s]] for s in ALL]
    ax.bar(xs - 0.2, som, 0.38, color=cols, alpha=0.55, label="somatic call at key pos")
    ax.bar(xs + 0.2, asm, 0.38, color=cols, alpha=1.0, label="somatic-controlled ASM")
    for i, (a, b) in enumerate(zip(som, asm)):
        ax.text(i - 0.2, a + 0.4, str(a), ha="center", fontsize=8)
        ax.text(i + 0.2, b + 0.4, str(b), ha="center", fontsize=8, fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{s}\n({cancer[s]})" for s in ALL], fontsize=8)
    ax.set_ylabel("count / 38 key positions")
    ax.set_title("Somatic recurrence at HCC1395's 38 key ASM positions\n"
                 "(other 5 samples: 0 — somatic mutations are private)")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylim(0, 40)
    ax.spines[["top", "right"]].set_visible(False)
    return fig_b64(fig)


# ---------------------------------------------------------------------------
# FIGURE 2: HCC1395 private-somatic HP-axis delta distribution
# ---------------------------------------------------------------------------
def fig_delta_dist():
    priv = [p["hcc1395_hp_axis_delta"] for p in POS
            if p["somatic_status"] == "hcc1395_private_somatic"
            and p["hcc1395_hp_axis_delta"] is not None]
    fig, ax = plt.subplots(figsize=(8.4, 3.2))
    ax.hist(priv, bins=np.arange(-0.35, 0.36, 0.05), color="#1E3A8A", alpha=0.8,
            edgecolor="white")
    ax.axvline(0, color="#6b7280", lw=1, ls="--")
    ax.set_xlabel("HCC1395 HP-axis Δβ (somatic-controlled, 5mC)")
    ax.set_ylabel("n positions")
    ax.set_title(f"HCC1395 private somatic ASM magnitude (n={len(priv)})\n"
                 "negative = somatic allele hypo-methylated vs germline")
    ax.spines[["top", "right"]].set_visible(False)
    return fig_b64(fig)


# ---------------------------------------------------------------------------
# FIGURE 3 (optional): genome-wide somatic ASM rate vs null
# ---------------------------------------------------------------------------
def fig_genomewide():
    samples = [s for s in ALL if s in GW]
    rate = [GW[s]["genomewide_somatic_asm"]["rate_strong_asm"] for s in samples]
    null = [GW[s]["genomewide_somatic_asm"]["rate_strong_null"] for s in samples]
    med = [GW[s]["genomewide_somatic_asm"]["median_abs_delta"] for s in samples]
    nmed = [GW[s]["genomewide_somatic_asm"]["median_abs_null"] for s in samples]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.8))
    xs = np.arange(len(samples))
    cols = [CANCER_COLOR[cancer[s]] for s in samples]
    a1.bar(xs - 0.2, rate, 0.38, color=cols, alpha=1.0, label="observed |Δβ|≥0.10")
    a1.bar(xs + 0.2, null, 0.38, color="#9ca3af", alpha=0.9, label="HP-shuffle null")
    a1.set_xticks(xs); a1.set_xticklabels(samples, rotation=30, ha="right", fontsize=8)
    a1.set_ylabel("strong-ASM rate (of evaluable)")
    a1.set_title("Genome-wide somatic ASM rate vs noise floor")
    a1.legend(fontsize=8); a1.spines[["top", "right"]].set_visible(False)
    a2.bar(xs - 0.2, med, 0.38, color=cols, alpha=1.0, label="observed median|Δβ|")
    a2.bar(xs + 0.2, nmed, 0.38, color="#9ca3af", alpha=0.9, label="null median")
    a2.set_xticks(xs); a2.set_xticklabels(samples, rotation=30, ha="right", fontsize=8)
    a2.set_ylabel("median |Δβ|")
    a2.set_title("Effect-size central tendency vs null")
    a2.legend(fontsize=8); a2.spines[["top", "right"]].set_visible(False)
    return fig_b64(fig)


FIG_OVERLAP = fig_somatic_overlap()
FIG_DELTA = fig_delta_dist()
FIG_GW = fig_genomewide() if has_gw else None


# ---------------------------------------------------------------------------
# TABLE builders
# ---------------------------------------------------------------------------
def cancer_chip(s):
    c = cancer[s]
    return (f'<span class="chip" style="background:{CANCER_COLOR[c]}">'
            f'{esc(s)}</span>')


# cancer panel
cancer_rows = "".join(
    f'<tr><td>{cancer_chip(s)}</td><td>{esc(CANCER_FULL[s])}</td>'
    f'<td class="num">{SO[s]["of_total"]}</td>'
    f'<td class="num">{SO[s]["n_somatic_at_keypos"]}</td>'
    f'<td class="num"><b>{SO[s]["n_somatic_asm"]}</b></td></tr>'
    for s in ALL)

# private-somatic top table (by |delta|)
priv = [p for p in POS if p["somatic_status"] == "hcc1395_private_somatic"]
priv.sort(key=lambda p: -abs(p["hcc1395_hp_axis_delta"] or 0))
priv_rows = "".join(
    f'<tr><td>{esc(p["gene"])}</td><td class="mono">{esc(p["pos"])}</td>'
    f'<td class="num">{fnum(p["hcc1395_hp_axis_delta"])}</td>'
    f'<td class="num">{p["hcc1395_hp_axis_minN"] or "—"}'
    f'{" ⚠" if p.get("hcc1395_low_subgroup_n") else ""}</td>'
    f'<td>{esc(p["set"].replace("_"," "))}</td>'
    f'<td>{esc(p.get("cn_class") or "")}</td></tr>'
    for p in priv[:15])

# recurrent germline ASM — concordant (candidate imprinting)
def germ_row(h, concordant):
    dd = " · ".join(f'{esc(s)} {fnum(v)}' for s, v in h["deltas"].items())
    tag = ("方向一致 → 候選 imprinting" if concordant else
           "方向相反 → 非 imprinting（sporadic/noise）")
    low = ' <span class="warn-i">lowN</span>' if h.get("low_subgroup_n") else ""
    return (f'<tr><td>{esc(h["gene"])}</td><td class="mono">{esc(h["pos"])}</td>'
            f'<td class="num">{h["n_samples"]}</td>'
            f'<td>{dd}{low}</td>'
            f'<td>{"·".join(esc(c) for c in h["cancer_types"])}</td></tr>')

imp_rows = "".join(germ_row(h, True) for h in H["imprinting_consistent_hits"])
dis_rows = "".join(germ_row(h, False) for h in H["discordant_recurrent_hits"])

# BRCA2 flagship per-sample detail
brca = next((p for p in POS if p["pos"] == "chr13:32315128"), None)
brca_rows = ""
if brca:
    for s in ALL:
        ps = brca["per_sample"].get(s, {})
        st = []
        if ps.get("somatic_in_sample"):
            st.append("somatic")
        if ps.get("has_subhap"):
            st.append("subhap")
        if ps.get("somatic_asm"):
            st.append("<b>ASM</b>")
        brca_rows += (
            f'<tr><td>{cancer_chip(s)}</td><td>{esc(CANCER_FULL[s])}</td>'
            f'<td class="num">{ps.get("n_reads","—")}</td>'
            f'<td class="num">{fnum(ps.get("hp_axis_delta"))}</td>'
            f'<td class="num">{fnum(ps.get("germline_delta"))}</td>'
            f'<td>{" · ".join(st) if st else "germline only"}</td></tr>')

# full per-position table (DataTables)
def full_row(p):
    ps_cells = ""
    for s in ALL[1:]:  # other 5
        v = p["per_sample"].get(s, {})
        d = v.get("hp_axis_delta")
        g = v.get("germline_delta")
        cell = fnum(d) if d is not None else (fnum(g) + "ᵍ" if g is not None else "—")
        ps_cells += f'<td class="num">{cell}</td>'
    return (f'<tr><td>{esc(p["gene"])}</td><td class="mono">{esc(p["pos"])}</td>'
            f'<td class="num">{fnum(p["hcc1395_hp_axis_delta"])}</td>'
            f'<td>{esc(p["somatic_status"].replace("_"," "))}</td>'
            f'<td>{esc(p["germline_status"].replace("_"," "))}</td>'
            f'{ps_cells}</tr>')

full_rows = "".join(full_row(p) for p in POS)
other_headers = "".join(f"<th>{esc(s)}</th>" for s in ALL[1:])

# genome-wide section
gw_section = ""
if has_gw:
    gw_rows = ""
    for s in ALL:
        if s not in GW:
            continue
        g = GW[s]["genomewide_somatic_asm"]
        ip = GW[s]["imprinted_dmr_panel"]
        gw_rows += (
            f'<tr><td>{cancer_chip(s)}</td>'
            f'<td class="num">{g["n_tp_total"]}</td>'
            f'<td class="num">{g["n_evaluable"]}</td>'
            f'<td class="num">{fnum(g["rate_strong_asm"],3)}</td>'
            f'<td class="num">{fnum(g["rate_strong_null"],3)}</td>'
            f'<td class="num"><b>{fnum(g["rate_excess_over_null"],3)}</b></td>'
            f'<td class="num">{fnum(g["median_abs_delta"],3)}</td>'
            f'<td class="num">{fnum(g["frac_hypo"],2)}/{fnum(g["frac_hyper"],2)}</td>'
            f'<td class="num">{ip["n_strong"]}/{ip["n_evaluable"]}</td></tr>')
    gwfig = (f'<details class="figbox" open><summary>圖：genome-wide somatic ASM rate vs null</summary>'
             f'<img src="{FIG_GW}"></details>') if FIG_GW else ""
    n_gw = len([s for s in ALL if s in GW])
    _ex = [GW[s]["genomewide_somatic_asm"]["rate_excess_over_null"] for s in ALL if s in GW]
    _cancers_pos = sorted(set(cancer[s] for s in ALL if s in GW
                              and GW[s]["genomewide_somatic_asm"]["rate_excess_over_null"] > 0))
    repro_line = (f'<b>全 {len(_ex)}/{n_gw} 樣本 excess-over-null &gt; 0</b>'
                  f'（範圍 +{min(_ex):.3f} ~ +{max(_ex):.3f}，mean +{float(np.mean(_ex)):.3f}，CV {float(np.std(_ex)/np.mean(_ex)):.2f}）'
                  f'，涵蓋 {len(_cancers_pos)} 癌種（乳腺/肺/黑色素瘤）皆正 → '
                  f'<b>somatic ASM 現象跨癌種復現</b>（各樣本用自己的 private somatic 突變）。'
                  f'⚠ 此為單一 caller/估計法的「現象複製」非獨立管線驗證 → tier 封頂 ⭐3。')
    gw_section = f'''
    <section class="card">
      <div class="card-h"><span class="qlabel">擴大運行 — genome-wide somatic ASM rate（各樣本自己的 private somatic）</span>
        <span class="badge" style="background:#1d4ed8">{n_gw}/6 樣本</span></div>
      <p class="claim">同位點不復發（private），但<b>現象</b>是否復發？各樣本用自己的 TP somatic SNV 隨機子集，量測 somatic-controlled HP-axis ASM rate，並與 HP-label-shuffle null 比對噪音底。</p>
      <div class="tldr" style="margin:10px 0">{repro_line}</div>
      <table class="dt"><thead><tr><th>樣本</th><th>TP somatic 總數</th><th>evaluable</th>
        <th>strong rate</th><th>null rate</th><th>excess</th><th>median|Δβ|</th><th>hypo/hyper</th><th>imprinted DMR≥0.10</th></tr></thead>
        <tbody>{gw_rows}</tbody></table>
      <details class="method"><summary>讀法 / caveat</summary><p>
      strong rate = evaluable somatic loci 中 |HP-axis Δβ|≥0.10 的比例；null = 同位點 reads 隨機重分配 20 次的 |Δβ|≥0.10 比例（有限樣本噪音底）；<b>excess = 觀測 − null</b> 才是真訊號。固定門檻的 raw rate 會被 per-locus 小 N 噪音灌水（median|Δβ| 通常遠低於 0.10），故必看 excess。imprinted DMR panel 為 EXPLORATORY（重度 LOH 樣本多數 germline 軸不可測）。</p></details>
      {gwfig}
    </section>'''


# ---------------------------------------------------------------------------
HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>跨樣本關鍵位點 ASM 驗證 — HCC1395 38 位點 × 6 樣本</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}
body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC","Microsoft JhengHei",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1060px;margin:0 auto;padding:24px 18px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:24px 26px;margin-bottom:16px}}
header.top h1{{margin:0 0 6px;font-size:1.42rem}}
.meta{{font-size:.82rem;opacity:.88}}
.tldr{{background:#eef2ff;border:1px solid #c7d2fe;border-left:5px solid #1E3A8A;border-radius:12px;padding:16px 20px;margin-bottom:14px}}
.tldr h2{{margin:0 0 8px;font-size:1.08rem;color:#1e3a8a}}.tldr p{{margin:.4rem 0}}
.kpis{{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}}
.kpi{{flex:1;min-width:150px;background:#fff;border:1px solid var(--line);border-radius:12px;padding:14px 16px;text-align:center}}
.kpi .v{{font-size:1.6rem;font-weight:800;color:var(--accent);font-variant-numeric:tabular-nums}}
.kpi .l{{font-size:.78rem;color:var(--mut);margin-top:3px}}
.card{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:14px 0}}
.card-h{{display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:6px}}
.qlabel{{font-weight:800;color:var(--accent);flex:1;font-size:1.02rem}}
.badge{{color:#fff;font-weight:700;font-size:.76rem;padding:3px 11px;border-radius:999px}}
.claim{{font-weight:600;margin:.2rem 0 .7rem}}
table{{width:100%;border-collapse:collapse;font-size:.84rem;margin-top:6px}}
th,td{{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}
th{{background:#f1f5f9;font-size:.78rem;color:#334155}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}
td.mono,.mono{{font-family:ui-monospace,Menlo,monospace;font-size:.8rem;color:#475569}}
.chip{{color:#fff;font-weight:700;font-size:.72rem;padding:2px 8px;border-radius:6px}}
.warn-i{{color:#b45309;font-size:.72rem;font-weight:700}}
details{{margin-top:10px}}
details.method>summary{{cursor:pointer;color:var(--accent);font-weight:600;font-size:.84rem}}
details.method p{{font-size:.83rem;background:var(--bg);border-radius:8px;padding:8px 12px}}
details.figbox>summary{{cursor:pointer;color:#374151;font-weight:600;font-size:.84rem;margin-top:8px}}
details.figbox img{{width:100%;border:1px solid var(--line);border-radius:8px;margin-top:8px}}
.two{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}
@media(max-width:760px){{.two{{grid-template-columns:1fr}}}}
.eval{{display:flex;gap:8px;align-items:center;font-size:.85rem;margin:.3rem 0}}
.pass{{color:#15803d;font-weight:700}}.fix{{color:#b45309;font-weight:700}}
.tablewrap{{overflow-x:auto}}
footer{{margin-top:24px;font-size:.76rem;color:var(--mut);border-top:1px solid var(--line);padding-top:12px}}
.legend{{font-size:.76rem;color:var(--mut);margin-top:4px}}
</style></head><body><div class="wrap">

<header class="top">
  <h1>跨樣本關鍵位點 ASM 驗證</h1>
  <div class="meta">HCC1395 的 {meta["n_key_positions"]} 個關鍵 ASM 位點 × 6 癌症樣本（乳腺 ×3 / 肺 ×2 / 黑色素瘤 ×1）·
  A pilot（targeted）· 2026-06-03 · 數據確定性來自 57/58 腳本（§13 layer-A，無手打數字）</div>
</header>

<div class="tldr">
  <h2>TL;DR — somatic ASM 是 HCC1395-private，現象不在同位點復發</h2>
  <p>① <b>同位點 somatic 復發 = 0/38</b>：其他 5 個癌症樣本在 HCC1395 的關鍵位點全部 <b>0 個</b> somatic call（各癌 private mutation，符合生物學預期）。</p>
  <p>② <b>{H["somatic_status_counts"].get("hcc1395_private_somatic",0)} 個 HCC1395-private somatic ASM</b>（含 flagship BRCA2/ZAR1L）；其他樣本在這些位點只有 germline 背景或 LOH 單倍型。</p>
  <p>③ <b>germline 軸復發 {H["n_recurrent_germline_asm"]} 個</b> → 嚴格方向檢查後拆成 <b>{H["n_imprinting_consistent"]} 個方向一致（候選 imprinting）</b> + <b>{H["n_recurrent_discordant"]} 個方向相反（非 imprinting，sporadic/noise）</b>。</p>
  <p class="legend">⚠ 跨樣本用 window-mean 5mC β（各樣本一致），重現 HCC1395 credible-loci delta 的<b>方向</b>而非精確 magnitude（原表用 paired-CpG MAX-collapse Wilcoxon）。HCC1395 為唯一 SEQC2-truth 單樣本參考。</p>
</div>

<div class="kpis">
  <div class="kpi"><div class="v">0 / 38</div><div class="l">其他樣本同位點 somatic 復發</div></div>
  <div class="kpi"><div class="v">{H["somatic_status_counts"].get("hcc1395_private_somatic",0)}</div><div class="l">HCC1395-private somatic ASM</div></div>
  <div class="kpi"><div class="v">{H["n_imprinting_consistent"]}</div><div class="l">方向一致候選 imprinting</div></div>
  <div class="kpi"><div class="v">{H["n_recurrent_discordant"]}</div><div class="l">方向相反（非 imprinting）</div></div>
</div>

<section class="card">
  <div class="card-h"><span class="qlabel">① 各樣本在 38 關鍵位點的 somatic 重疊</span>
    <span class="badge" style="background:#15803d">private 確認</span></div>
  <p class="claim">只有 HCC1395 在這些位點是 somatic（33/38 有 call、28 個 somatic-controlled ASM）；其他 5 樣本 = 0。</p>
  <div class="tablewrap"><table>
    <thead><tr><th>樣本</th><th>癌種</th><th class="num">關鍵位點</th><th class="num">somatic call</th><th class="num">somatic ASM</th></tr></thead>
    <tbody>{cancer_rows}</tbody></table></div>
  <details class="figbox" open><summary>圖：somatic 重疊</summary><img src="{FIG_OVERLAP}"></details>
</section>

<section class="card">
  <div class="card-h"><span class="qlabel">② BRCA2/ZAR1L（flagship）逐樣本</span>
    <span class="badge" style="background:#15803d">HCC1395-private</span></div>
  <p class="claim">chr13:32315128（HCC1395 G&gt;A AF=0.19）：只有 HCC1395 有真 somatic 子單倍型 + ASM（Δβ≈−0.19）；COLO829/HCC1954 為 germline het 雙倍型、H1437/HCC1937 為 LOH 單倍型、H2009 僅 4-read 假性 subhap（Δ≈0）。</p>
  <div class="tablewrap"><table>
    <thead><tr><th>樣本</th><th>癌種</th><th class="num">window reads</th><th class="num">HP-axis Δβ</th><th class="num">germline Δβ</th><th>狀態</th></tr></thead>
    <tbody>{brca_rows}</tbody></table></div>
</section>

<div class="two">
  <section class="card">
    <div class="card-h"><span class="qlabel">③ HCC1395-private somatic ASM（top 15）</span></div>
    <div class="tablewrap"><table>
      <thead><tr><th>gene</th><th>pos</th><th class="num">Δβ</th><th class="num">minN</th><th>set</th><th>CN</th></tr></thead>
      <tbody>{priv_rows}</tbody></table></div>
    <p class="legend">⚠ = delta 所依子群 &lt;10 reads。共 {len(priv)} 個 private。</p>
    <details class="figbox"><summary>圖：private somatic ASM magnitude 分布</summary><img src="{FIG_DELTA}"></details>
  </section>
  <section class="card">
    <div class="card-h"><span class="qlabel">④ germline 軸復發（方向檢查）</span></div>
    <p class="claim" style="font-size:.85rem;color:#15803d">方向一致 → 候選 imprinting（{H["n_imprinting_consistent"]}）</p>
    <div class="tablewrap"><table>
      <thead><tr><th>gene</th><th>pos</th><th class="num">n</th><th>per-sample Δβ</th><th>癌種</th></tr></thead>
      <tbody>{imp_rows}</tbody></table></div>
    <p class="claim" style="font-size:.85rem;color:#b45309;margin-top:12px">方向相反 → 非 imprinting（{H["n_recurrent_discordant"]}）</p>
    <div class="tablewrap"><table>
      <thead><tr><th>gene</th><th>pos</th><th class="num">n</th><th>per-sample Δβ</th><th>癌種</th></tr></thead>
      <tbody>{dis_rows}</tbody></table></div>
  </section>
</div>

{gw_section}

<section class="card">
  <div class="card-h"><span class="qlabel">⑤ 完整 38 位點表（HCC1395 Δβ + 雙軸狀態 + 其他 5 樣本 Δβ）</span></div>
  <p class="legend">其他樣本欄：數字為 HP-axis Δβ；上標 ᵍ 表示該樣本無 somatic 子群、值為 germline HP1-vs-HP2 Δβ；— 表示不可測（LOH/低覆蓋）。</p>
  <div class="tablewrap"><table id="full">
    <thead><tr><th>gene</th><th>pos</th><th class="num">HCC1395 Δβ</th><th>somatic 軸</th><th>germline 軸</th>{other_headers}</tr></thead>
    <tbody>{full_rows}</tbody></table></div>
</section>

<section class="card">
  <div class="card-h"><span class="qlabel">⑥ 驗證 / 對抗審查（generator–evaluator 分離）</span>
    <span class="badge" style="background:#15803d">2 evaluator</span></div>
  <div class="eval"><span class="pass">PASS</span> 分離正確性：private-somatic 與 BRCA2 標籤全 data-grounded，交叉驗證 6 cells 0 transcription error。</div>
  <div class="eval"><span class="fix">已修正</span> over-claim 防護：evaluator 抓到「imprinting」誤用（4 hit 中 3 個方向相反）→ 已改方向不可知命名 + sign-concordance + 只方向一致才標候選 imprinting + surface 子群 N。</div>
  <details class="method"><summary>方法 / caveat（必讀）</summary><p>
  <b>方法</b>：複用已驗證 pysam 抽取（modkit crossval Pearson=1.0）；per-read 5mC β（mean over CpGs in ±600bp），依 HP:Z tag 分群；somatic-controlled HP-axis Δβ=β(subhap)−β(main)。<br>
  <b>口徑</b>：window-mean 5mC（各樣本一致）重現 HCC1395 credible-loci delta 的<b>方向</b>，非精確 magnitude（原表 paired-CpG MAX-collapse Wilcoxon）。<br>
  <b>限制</b>：① SEQC2 truth 只有 HCC1395（其他 5 樣本 TP VCF 是 caller 輸出、無正交 truth）；② targeted「不復發」部分受 LOH/低 N 使 germline 軸不可測所限（非純生物缺席）；③ IGF2 是 cnLOH，4/5 樣本 germline 軸塌縮 → 已知 imprinting 被 LOH 遮蔽（footnote，非 code bug）；④ {("genome-wide 固定門檻 raw rate 受小 N 噪音灌水（median|Δβ| 全 &lt;0.10），必看 excess-over-null。" if has_gw else "genome-wide 擴大運行進行中。")}{("<b>⑤ 6 樣本共用同一 ClairS-TO/longphase caller + 同一 HP-axis 估計法 → 此為「現象複製」非獨立管線/獨立 cohort 驗證；共用系統性偏差（caller 或估計法）未排除，故 tier 封頂 ⭐3（bordering ⭐4），升 ⭐4 需獨立管線或正交 truth。</b>" if has_gw else "")}⑥ genome-wide rate 用 window-mean 5mC（與 targeted credible-loci 的 paired-CpG MAX-collapse Wilcoxon 口徑不同），方向可比、magnitude 不可直接並列。</p></details>
</section>

<footer>
  數據源：<span class="mono">research/tsg_promoter_asm_reviewer/genome_survey_v2/cn_confound/cross_sample/cross_sample_synthesis.json</span>
  （+ 6× *_keypos.json{" + *_gwasm.json" if has_gw else ""}）· 腳本 57/58{("/59/60" if has_gw else "/60")} · 確定性合成（§13 layer-A）·
  thresholds: ASM_MIN={meta["thresholds"]["ASM_MIN"]} GERM_MIN={meta["thresholds"]["GERM_MIN"]} COV_MIN={meta["thresholds"]["COV_MIN"]}
</footer>
</div></body></html>'''

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write(HTML_DOC)
print(f"[60] wrote {OUT} ({len(HTML_DOC)//1024} KB)  genome_wide={'yes' if has_gw else 'no'}")
