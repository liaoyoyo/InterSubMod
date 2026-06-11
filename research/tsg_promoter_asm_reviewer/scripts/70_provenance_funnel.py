#!/usr/bin/env python3
"""
70 - Complete provenance funnel: from FULL longphase-S output down to credible ASM loci.

Answers: "are the credible loci a COMPLETE statistic+filtering of ALL longphase-S
output positions?" -> NO. This makes the full funnel explicit, per sample, so the
exact subset (and everything excluded at each stage) is visible.

LAYERS (per sample, longphase_s/):
  somatic_pass.vcf.gz   = caller's FULL somatic PASS output (all types, genome-wide)
  truth_scoped.vcf.gz   = benchmark truth set (in scored regions)
  filtered_snv_tp/fp/fn = SNV benchmarked subset vs truth (TP/FP/FN)
then the discovery funnel (62) on TP only:
  TP -> CpG-island-proximal -> HP-axis survey -> regimeA(relax/strict) -> ARI-eval -> credible

§13 layer-A: VCF counts re-counted here via pysam; funnel from *_credible_discovery.json.

Output: genome_survey_v2/cn_confound/cross_sample/provenance_funnel.json
        docs/experiments/in_progress/2026/06/20260603_provenance_funnel_01.standalone.html
"""
import os, json, html, io, base64
import pysam
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"
OUTJSON = f"{CS}/provenance_funnel.json"
OUTHTML = ("/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/06/"
           "20260603_provenance_funnel_01.standalone.html")

ALL = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}
CCOLOR = {"breast": "#be123c", "lung": "#0e7490", "melanoma": "#7c3aed"}
RUNDIR = {s: ("20260314_HCC1395_paired_full_full_complete_matrix" if s == "HCC1395"
              else f"20260315_{s}_paired_full_full_complete_matrix") for s in ALL}


def count(vcf):
    try:
        vf = pysam.VariantFile(vcf)
        return sum(1 for _ in vf)
    except Exception:
        return None


def esc(x):
    return html.escape(str(x))


rows = {}
for s in ALL:
    base = (f"/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{s}/paired_full/"
            f"{RUNDIR[s]}/longphase_s")
    fn = json.load(open(f"{CS}/{s}_credible_discovery.json"))["funnel"]
    rows[s] = dict(
        cancer=CANCER[s],
        somatic_pass=count(f"{base}/somatic_pass.vcf.gz"),
        truth=count(f"{base}/truth_scoped.vcf.gz"),
        tp=count(f"{base}/filtered_snv_tp.vcf.gz"),
        fp=count(f"{base}/filtered_snv_fp.vcf.gz"),
        fn=count(f"{base}/filtered_snv_fn.vcf.gz"),
        tp_used=fn["n_tp_somatic"],
        island_prox=fn["n_cpg_island_proximal"],
        hp_survey=fn["n_hp_axis_survey"],
        regimeA_relax=fn["n_regimeA_relax"], regimeA_strict=fn["n_regimeA_strict"],
        ari_eval=fn["n_ari_evaluable"],
        credible_relax=fn["n_credible_pass_tierA_relax"],
        credible_strict=fn["n_credible_pass_tierA_strict"],
    )
    print(f"[70] {s}: somatic_pass={rows[s]['somatic_pass']} TP={rows[s]['tp']} "
          f"island={rows[s]['island_prox']} survey={rows[s]['hp_survey']} "
          f"regimeA={rows[s]['regimeA_relax']}/{rows[s]['regimeA_strict']} "
          f"ARIeval={rows[s]['ari_eval']} credible={rows[s]['credible_relax']}/{rows[s]['credible_strict']}")

out = dict(meta=dict(script="70_provenance_funnel.py", samples=ALL,
                     note=("credible loci are NOT a complete statistic of all longphase-S "
                           "positions: TP-only (excludes FP/FN/all non-SNV/all non-benchmarked "
                           "somatic_pass), then CpG-island-proximal only, then require somatic "
                           "subhaplotype + >=100/30 paired CpG + ARI>=0.30 + placebo<0.10.")),
           per_sample=rows)
with open(OUTJSON, "w") as f:
    json.dump(out, f, indent=2, default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)
print(f"[70] wrote {OUTJSON}")

# ---- figure: HCC1395 funnel (log) ----
def fig_funnel():
    s = "HCC1395"; r = rows[s]
    stages = ["somatic_pass\n(全輸出)", "TP\n(benchmarked)", "CpG-island\nproximal",
              "HP-axis\nsurvey", "regimeA\nrelax", "ARI\nevaluable", "credible\nstrict"]
    vals = [r["somatic_pass"], r["tp"], r["island_prox"], r["hp_survey"],
            r["regimeA_relax"], r["ari_eval"], r["credible_strict"]]
    fig, ax = plt.subplots(figsize=(9.2, 3.8))
    x = np.arange(len(stages))
    ax.bar(x, vals, color="#1E3A8A", alpha=0.85)
    for i, v in enumerate(vals):
        ax.text(i, v * 1.4 + 1, f"{v:,}", ha="center", fontsize=8, fontweight="bold")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(stages, fontsize=7.5)
    ax.set_ylabel("位點數 (log)")
    pct = 100 * r["credible_strict"] / r["somatic_pass"]
    ax.set_title(f"HCC1395 provenance funnel — credible {r['credible_strict']} = "
                 f"{pct:.3f}% of全 somatic_pass 輸出 ({r['somatic_pass']:,})")
    ax.spines[["top", "right"]].set_visible(False)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=130, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()

FIG = fig_funnel()

trows = ""
for s in ALL:
    r = rows[s]
    chip = f'<span class="chip" style="background:{CCOLOR[r["cancer"]]}">{esc(s)}</span>'
    pct_tp = (100 * r["credible_strict"] / r["tp"]) if r["tp"] else 0
    pct_all = (100 * r["credible_strict"] / r["somatic_pass"]) if r["somatic_pass"] else 0
    trows += (f'<tr><td>{chip}</td>'
              f'<td class="num">{r["somatic_pass"]:,}</td>'
              f'<td class="num">{r["tp"]:,}</td><td class="num">{r["fp"]:,}</td><td class="num">{r["fn"]:,}</td>'
              f'<td class="num">{r["island_prox"]:,}</td>'
              f'<td class="num">{r["hp_survey"]:,}</td>'
              f'<td class="num">{r["regimeA_relax"]}/{r["regimeA_strict"]}</td>'
              f'<td class="num">{r["ari_eval"]}</td>'
              f'<td class="num"><b>{r["credible_relax"]}/{r["credible_strict"]}</b></td>'
              f'<td class="num">{pct_tp:.2f}% / {pct_all:.3f}%</td></tr>')

HTML_DOC = f'''<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Provenance funnel — longphase-S → credible</title>
<style>
:root{{--ink:#1f2937;--mut:#6b7280;--line:#e5e7eb;--bg:#f8fafc;--accent:#1E3A8A}}
*{{box-sizing:border-box}}body{{margin:0;font-family:-apple-system,"Noto Sans CJK TC",sans-serif;color:var(--ink);background:var(--bg);line-height:1.6}}
.wrap{{max-width:1120px;margin:0 auto;padding:24px 18px 80px}}
header.top{{background:linear-gradient(135deg,#0f1f4d,#1E3A8A);color:#fff;border-radius:14px;padding:24px 26px;margin-bottom:16px}}
header.top h1{{margin:0 0 6px;font-size:1.36rem}}.meta{{font-size:.82rem;opacity:.88}}
.tldr{{background:#fef2f2;border:1px solid #fecaca;border-left:5px solid #dc2626;border-radius:12px;padding:16px 20px;margin-bottom:14px}}
.tldr h2{{margin:0 0 8px;font-size:1.05rem;color:#b91c1c}}.tldr p{{margin:.35rem 0}}
.card{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:14px 0}}
.qlabel{{font-weight:800;color:var(--accent);font-size:1rem}}
table{{width:100%;border-collapse:collapse;font-size:.82rem;margin-top:8px}}
th,td{{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left}}th{{background:#f1f5f9;font-size:.74rem}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}
.chip{{color:#fff;font-weight:700;font-size:.72rem;padding:2px 8px;border-radius:6px}}
.figbox img{{width:100%;border:1px solid var(--line);border-radius:8px;margin-top:8px}}.tablewrap{{overflow-x:auto}}
.excl{{font-size:.85rem}}.excl li{{margin:.2rem 0}}
footer{{margin-top:24px;font-size:.76rem;color:var(--mut);border-top:1px solid var(--line);padding-top:12px}}
</style></head><body><div class="wrap">
<header class="top"><h1>Provenance funnel — longphase-S 全輸出 → credible ASM 位點</h1>
<div class="meta">回答「credible 是否=全 longphase-S 輸出位點的完整統計」· A pilot · 2026-06-03 · §13 layer-A</div></header>

<div class="tldr">
  <h2>❌ 不是完整統計 — credible 是極小的過濾子集（~0.02–0.08% of 全 somatic 輸出）</h2>
  <p>credible 位點<b>只取 TP somatic SNV</b>（排除 FP/FN/全部非 SNV/未 benchmark 的 somatic_pass），再<b>只取 CpG-island 鄰近</b>，再要求<b>有 somatic 子單倍型 + ≥100/30 paired CpG + ARI≥0.30 + placebo&lt;0.10</b>。</p>
  <p>HCC1395：全 somatic_pass <b>{rows["HCC1395"]["somatic_pass"]:,}</b> → credible strict <b>{rows["HCC1395"]["credible_strict"]}</b>（{100*rows["HCC1395"]["credible_strict"]/rows["HCC1395"]["somatic_pass"]:.3f}%）。</p>
</div>

<section class="card">
  <span class="qlabel">完整 funnel（每樣本，每層位點數）</span>
  <div class="tablewrap"><table>
    <thead><tr><th>樣本</th><th class="num">somatic_pass<br>(全輸出)</th><th class="num">TP</th><th class="num">FP</th><th class="num">FN</th>
      <th class="num">island-prox<br>(of TP)</th><th class="num">HP-axis<br>survey</th><th class="num">regimeA<br>relax/strict</th>
      <th class="num">ARI-eval</th><th class="num">credible<br>relax/strict</th><th class="num">credible<br>%TP / %全輸出</th></tr></thead>
    <tbody>{trows}</tbody></table></div>
  <div class="figbox"><img src="{FIG}"></div>
</section>

<section class="card">
  <span class="qlabel">每層排除了什麼（為何不是完整統計）</span>
  <ul class="excl">
    <li><b>somatic_pass → TP</b>：排除 <b>FP</b>（caller 誤判，HCC1395 627）、<b>未在 truth-scoped SNV benchmark 內</b>的所有 PASS 呼叫（含 indel/multiallelic/truth 區外）。注意 somatic_pass({rows["HCC1395"]["somatic_pass"]:,}) ≫ TP+FP，故大量 somatic 呼叫根本沒進此分析。</li>
    <li><b>(另排除 FN {rows["HCC1395"]["fn"]:,})</b>：truth 有但 caller 漏判的位點，完全沒看。</li>
    <li><b>TP → CpG-island-proximal</b>：為 tractability + credible 集中於 island，只掃 ±2kb CpG-island 鄰近（HCC1395 1248/29754 = 4%）→ <b>非 island 區的 TP（96%）未掃</b>。</li>
    <li><b>island → HP-axis survey</b>：需有 somatic 子單倍型(HP:Z:N-M) + ≥5 paired CpG。</li>
    <li><b>survey → regimeA</b>：wilcoxon_p&lt;0.05 + n_paired_cpg≥100(strict)/30(relax) + extremity≤0.3。</li>
    <li><b>regimeA → ARI-eval</b>：兩 HP 組各 ≥8 reads（覆蓋 gate，COLO829 在此塌縮）。</li>
    <li><b>ARI-eval → credible</b>：ARI≥0.30 且 placebo&lt;0.10。</li>
  </ul>
  <p style="font-size:.84rem;color:#6b7280">⚠ 故 gallery/credible = <b>「TP somatic 中、CpG-island 鄰近、有強 ASM clustering」的子集表徵</b>，<b>不是</b> longphase-S 全輸出的完整篩選統計；也<b>不含</b> FP/FN（TP/FP 判別力另見 <span style="font-family:monospace">69_tp_vs_fp_ari</span>：ASM 對 TP/FP 無判別力）。</p>
</section>
<footer>數據源：<span style="font-family:monospace">cross_sample/provenance_funnel.json</span>（VCF 計數 pysam re-count + discovery funnel）· 腳本 70 · §13 layer-A</footer>
</div></body></html>'''

os.makedirs(os.path.dirname(OUTHTML), exist_ok=True)
with open(OUTHTML, "w") as f:
    f.write(HTML_DOC)
print(f"[70] wrote {OUTHTML} ({len(HTML_DOC)//1024} KB)")
